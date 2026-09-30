"""A common objective for greedy selection and bounded exhaustive comparison."""
from __future__ import annotations

import itertools
import math

from .models import PanelBudgetDecision


def panel_metrics(selected, hypotheses):
    ids = [h.hypothesis_id for h in hypotheses]
    signatures = {hid: [(c.hypothesis_perturbation or c.hypothesis_response).get(hid, 0.0) for c in selected] for hid in ids}
    coverage = {hid: max(signatures[hid], default=0.0) for hid in ids}
    prior = {h.hypothesis_id: max(0.05, min(1.0, h.discotope_prior if h.discotope_prior is not None else 0.5)) for h in hypotheses}
    pairs = [{"left": a, "right": b,
              "distance": math.sqrt(sum((x-y)**2 for x, y in zip(signatures[a], signatures[b])) / 8.0)}
             for a, b in itertools.combinations(ids, 2)]
    distances = [p["distance"] for p in pairs]
    uniform_coverage = sum(coverage.values()) / max(1, len(ids))
    weighted_coverage = sum(coverage[hid]*prior[hid] for hid in ids) / max(1e-12, sum(prior.values()))
    mean = sum(distances) / max(1, len(pairs))
    pair_weight = sum(prior[p["left"]]*prior[p["right"]] for p in pairs)
    weighted_mean = sum(p["distance"]*prior[p["left"]]*prior[p["right"]] for p in pairs) / max(1e-12, pair_weight)
    near = [p for p in pairs if p["distance"] < 0.1]
    resolved_fraction = (len(pairs)-len(near))/max(1, len(pairs))
    information = .5*(uniform_coverage+weighted_coverage) + .75*(mean+weighted_mean) + .5*min(distances, default=0) + .25*resolved_fraction
    adjacency = {hid: set() for hid in ids}
    for p in near:
        adjacency[p["left"]].add(p["right"])
        adjacency[p["right"]].add(p["left"])
    groups, remaining = [], set(ids)
    while remaining:
        start = min(remaining)
        group, stack = {start}, [start]
        remaining.remove(start)
        while stack:
            for neighbor in sorted(adjacency[stack.pop()] & remaining):
                remaining.remove(neighbor)
                group.add(neighbor)
                stack.append(neighbor)
        if len(group) > 1:
            groups.append(sorted(group))
    return dict(information=information, coverage=uniform_coverage,
                weighted_coverage=weighted_coverage, mean_separation=mean,
                weighted_separation=weighted_mean, minimum_distance=min(distances, default=0),
                near_pairs=near, pairwise_distances=pairs, groups=groups,
                uncovered=[hid for hid in ids if coverage[hid] < .1], signatures=signatures)


def construct_signature(candidate):
    return tuple(sorted((m.residue.uniprot_pos, m.mutant_aa) for m in candidate.mutations))


def objective(selected, hypotheses, risk_function, risk_weights=None):
    metrics = panel_metrics(selected, hypotheses)
    risk = sum(risk_function(c, risk_weights) + (.01 if c.order > 2 else 0) for c in selected)
    redundancy = 0.0
    ids = [h.hypothesis_id for h in hypotheses]
    for a, b in itertools.combinations(selected, 2):
        av, bv = a.hypothesis_perturbation or a.hypothesis_response, b.hypothesis_perturbation or b.hypothesis_response
        distance = math.sqrt(sum((av.get(hid, 0)-bv.get(hid, 0))**2 for hid in ids)/max(1, len(ids)))
        redundancy += max(0., float((risk_weights or {}).get("redundancy", .03))) * max(0., 1-distance/.25)
    return metrics["information"] - risk - redundancy, metrics, risk, redundancy


def greedy_panel(pool, hypotheses, budget, higher_order_max, risk_function, risk_weights=None, stop_parameters=None):
    decision = PanelBudgetDecision(**budget)
    stop_parameters = stop_parameters or {}
    selected = []
    old_value, old, old_risk, old_redundancy = objective([], hypotheses, risk_function, risk_weights)
    for _ in range(decision.maximum_constructs):
        choices = []
        for candidate in pool:
            if candidate.qualification == "REVIEW_REQUIRED" or any(m.review_required for m in candidate.mutations):
                continue
            if candidate.order >= 4 or candidate.search_role == "EXPLORATORY" or not candidate.connected:
                continue
            if candidate.order > 2 and sum(c.order > 2 for c in selected) >= higher_order_max:
                continue
            if any(construct_signature(candidate) == construct_signature(c) for c in selected):
                continue
            response = candidate.hypothesis_perturbation or candidate.hypothesis_response
            if any(all(abs(response.get(h.hypothesis_id, 0)-(c.hypothesis_perturbation or c.hypothesis_response).get(h.hypothesis_id, 0)) < 1e-6 for h in hypotheses) for c in selected):
                continue
            value, metrics, risk, redundancy = objective(selected+[candidate], hypotheses, risk_function, risk_weights)
            choices.append((value-old_value, candidate, metrics, value, risk, redundancy))
        if not choices:
            decision.stop_reason = "NO_FEASIBLE_NONREDUNDANT_CANDIDATES"
            break
        gain, winner, metrics, value, risk, redundancy = min(choices, key=lambda c: (-c[0], c[1].candidate_id))
        information_gain = metrics["information"] - old["information"]
        relative = information_gain / max(old["information"], 1e-12)
        pair_reduction = len(old["near_pairs"]) - len(metrics["near_pairs"])
        separation_improved = metrics["mean_separation"] > old["mean_separation"] + 1e-12 or pair_reduction > 0
        if gain <= 1e-12:
            decision.stop_reason = "NO_POSITIVE_NET_INFORMATION_GAIN"
            break
        if decision.adaptive_budget and len(selected) >= decision.minimum_constructs and (
            information_gain < float(stop_parameters.get("minimum_marginal_information_gain", .02))
            or relative < float(stop_parameters.get("minimum_relative_information_gain", .05))
            or not separation_improved
        ):
            decision.stop_reason = "INFORMATION_GAIN_BELOW_STOPPING_POLICY"
            break
        selected.append(winner)
        decision.budget_curve.append(dict(construct_count=len(selected), candidate_id=winner.candidate_id,
            objective=value, net_gain=gain, information_gain=information_gain, relative_information_gain=relative,
            coverage_gain=metrics["coverage"]-old["coverage"],
            mean_separation_gain=metrics["mean_separation"]-old["mean_separation"],
            minimum_separation_gain=metrics["minimum_distance"]-old["minimum_distance"],
            unresolved_pair_reduction=pair_reduction, risk_delta=risk-old_risk,
            redundancy_delta=redundancy-old_redundancy))
        old_value, old, old_risk, old_redundancy = value, metrics, risk, redundancy
    else:
        decision.stop_reason = "MAXIMUM_CONSTRUCT_BUDGET_REACHED"
    decision.selected_constructs = len(selected)
    decision.minimum_target_met = len(selected) >= decision.minimum_constructs
    return selected, old, decision


def exhaustive_comparison(pool, hypotheses, budget, risk_function, risk_weights=None, higher_order_max=2):
    """Offline check only; refuse searches larger than 300,000 panels."""
    if len(pool) > 20 or not 1 <= budget <= 8:
        raise ValueError("Exhaustive comparison requires <=20 candidates and budget <=8")
    if sum(math.comb(len(pool), k) for k in range(1, min(budget, len(pool))+1)) > 300000:
        raise ValueError("Exhaustive comparison exceeds 300000-panel work limit")
    best_value, best_ids = 0.0, []
    for k in range(1, min(budget, len(pool))+1):
        for subset in itertools.combinations(pool, k):
            if any(c.order >= 4 or not c.connected or c.search_role == "EXPLORATORY" or c.qualification == "REVIEW_REQUIRED" or any(m.review_required for m in c.mutations) for c in subset):
                continue
            if sum(c.order > 2 for c in subset) > higher_order_max or len({construct_signature(c) for c in subset}) != k:
                continue
            if any(all(abs((a.hypothesis_perturbation or a.hypothesis_response).get(h.hypothesis_id, 0)-(b.hypothesis_perturbation or b.hypothesis_response).get(h.hypothesis_id, 0)) < 1e-6 for h in hypotheses) for a,b in itertools.combinations(subset,2)):
                continue
            value = objective(subset, hypotheses, risk_function, risk_weights)[0]
            ids = sorted(c.candidate_id for c in subset)
            if value > best_value+1e-12 or (abs(value-best_value)<1e-12 and ids < best_ids):
                best_value, best_ids = value, ids
    greedy, _, _ = greedy_panel(pool, hypotheses, dict(minimum_constructs=min(4,budget), target_constructs=budget, maximum_constructs=budget, adaptive_budget=False), higher_order_max, risk_function, risk_weights)
    greedy_value = objective(greedy, hypotheses, risk_function, risk_weights)[0]
    return dict(optimal_candidate_ids=best_ids, optimal_objective=best_value,
                greedy_objective=greedy_value, absolute_gap=max(0., best_value-greedy_value),
                relative_gap=max(0., best_value-greedy_value)/max(best_value,1e-12))
