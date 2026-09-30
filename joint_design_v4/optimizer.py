from __future__ import annotations

import math
from collections import defaultdict

from .models import V4Candidate, V4Panel, V4Sector, V4SurfaceHypothesis


def _candidate_rank(candidate: V4Candidate) -> tuple:
    stability = candidate.stability
    spurs = (
        stability.spurs_multi_ddg_kcal_mol
        if stability and stability.spurs_multi_ddg_kcal_mol is not None
        else float("inf")
    )
    foldx = (
        stability.foldx_ddg_kcal_mol
        if stability and stability.foldx_ddg_kcal_mol is not None
        else float("inf")
    )
    return (
        candidate.order > 2,
        max(spurs, foldx),
        candidate.structural_risk,
        -min(candidate.coverage_by_scenario.values(), default=0.0),
        -sum(candidate.coverage_by_scenario.values()),
        candidate.order,
        -candidate.epitope_score,
        candidate.candidate_id,
    )


def _position_signature(candidate: V4Candidate) -> tuple[int, ...]:
    return tuple(sorted(option.residue.uniprot_pos or -1 for option in candidate.mutations))


def _diverse_take(
    pool: list[V4Candidate], count: int, already_selected: list[V4Candidate]
) -> list[V4Candidate]:
    if count <= 0:
        return []
    ranked = sorted(pool, key=_candidate_rank)
    used_candidates = {item.candidate_id for item in already_selected}
    used_positions = {_position_signature(item) for item in already_selected}
    chosen: list[V4Candidate] = []
    for candidate in ranked:
        signature = _position_signature(candidate)
        if candidate.candidate_id in used_candidates or signature in used_positions:
            continue
        chosen.append(candidate)
        used_candidates.add(candidate.candidate_id)
        used_positions.add(signature)
        if len(chosen) >= count:
            return chosen
    for candidate in ranked:
        if candidate.candidate_id in used_candidates:
            continue
        chosen.append(candidate)
        used_candidates.add(candidate.candidate_id)
        if len(chosen) >= count:
            break
    return chosen


def _information_take(
    pool: list[V4Candidate],
    hypotheses: list[V4SurfaceHypothesis],
    budget: int,
    higher_order_max: int,
) -> tuple[list[V4Candidate], float, float, int, float, float, list[list[str]], list[str], dict[str, list[float]]]:
    """Greedy panel selection by hypothesis coverage and pair separation."""
    ids = [item.hypothesis_id for item in hypotheses]
    pairs = [(left, right) for left in range(len(ids)) for right in range(left + 1, len(ids))]
    responses = {
        item.candidate_id: [item.hypothesis_response.get(hypothesis_id, 0.0) for hypothesis_id in ids]
        for item in pool
    }
    selected: list[V4Candidate] = []
    coverage = [0.0] * len(ids)
    separation_sq = [0.0] * len(pairs)
    signatures: set[tuple[int, ...]] = set()
    for _ in range(max(0, budget)):
        choices = []
        for candidate in pool:
            signature = _position_signature(candidate)
            if candidate.candidate_id in {item.candidate_id for item in selected} or signature in signatures:
                continue
            if candidate.order > 2 and sum(item.order > 2 for item in selected) >= higher_order_max:
                continue
            row = responses[candidate.candidate_id]
            if any(all(abs(left - right) < 1e-6 for left, right in zip(row, responses[item.candidate_id])) for item in selected):
                continue
            coverage_gain = sum(max(0.0, value - coverage[index]) for index, value in enumerate(row)) / max(1, len(ids))
            separation_gain = sum(
                math.sqrt(old + (row[left] - row[right]) ** 2) - math.sqrt(old)
                for index, (left, right) in enumerate(pairs)
                for old in [separation_sq[index]]
            ) / max(1, len(pairs))
            minimum_gain = (
                min(math.sqrt(separation_sq[index] + (row[left] - row[right]) ** 2)
                for index, (left, right) in enumerate(pairs)
                ) - min(math.sqrt(value) for value in separation_sq)
            ) if pairs else 0.0
            # Stability is already a formal gate. This small penalty breaks
            # near ties without turning ddG into an epitope probability.
            penalty = 0.005 * max(0.0, candidate.structural_risk) + (0.01 if candidate.order > 2 else 0.0)
            value = coverage_gain + 1.5 * separation_gain + 0.5 * minimum_gain - penalty
            choices.append((value, coverage_gain, separation_gain, candidate))
        if not choices:
            break
        value, coverage_gain, separation_gain, winner = min(
            choices, key=lambda item: (-item[0], -item[2], -item[1], item[3].candidate_id)
        )
        if value <= 0 and selected:
            break
        selected.append(winner)
        signatures.add(_position_signature(winner))
        row = responses[winner.candidate_id]
        coverage = [max(old, new) for old, new in zip(coverage, row)]
        separation_sq = [
            old + (row[left] - row[right]) ** 2
            for old, (left, right) in zip(separation_sq, pairs)
        ]
    separation = [math.sqrt(value / max(1, len(selected))) for value in separation_sq]
    prior = {
        item.hypothesis_id: max(0.05, min(1.0, item.discotope_prior if item.discotope_prior is not None else 0.5))
        for item in hypotheses
    }
    pair_weights = [prior[ids[left]] * prior[ids[right]] for left, right in pairs]
    weighted_separation = (
        sum(value * weight for value, weight in zip(separation, pair_weights)) / sum(pair_weights)
        if pair_weights and sum(pair_weights) else 0.0
    )
    signatures_by_hypothesis = {
        hypothesis_id: [responses[item.candidate_id][index] for item in selected]
        for index, hypothesis_id in enumerate(ids)
    }
    adjacency = {hypothesis_id: set() for hypothesis_id in ids}
    for distance, (left, right) in zip(separation, pairs):
        if distance < 0.1:
            adjacency[ids[left]].add(ids[right])
            adjacency[ids[right]].add(ids[left])
    unresolved = []
    remaining = set(ids)
    while remaining:
        start = min(remaining)
        group = {start}
        queue = [start]
        remaining.remove(start)
        while queue:
            for neighbor in sorted(adjacency[queue.pop()] & remaining):
                remaining.remove(neighbor)
                group.add(neighbor)
                queue.append(neighbor)
        if len(group) > 1:
            unresolved.append(sorted(group))
    return (
        selected,
        sum(coverage) / max(1, len(coverage)),
        sum(separation) / max(1, len(separation)),
        sum(value < 0.1 for value in separation),
        min(separation, default=0.0),
        weighted_separation,
        unresolved,
        [ids[index] for index, value in enumerate(coverage) if value < 0.1],
        signatures_by_hypothesis,
    )


def select_panel(
    candidates: list[V4Candidate],
    max_primary: int = 7,
    higher_order_max: int = 2,
    recommendations_per_sector: int = 6,
    epitope_priority_per_sector: int = 3,
    stability_fallback_per_sector: int = 3,
    sectors: list[V4Sector] | None = None,
    surface_areas: dict[int, float] | None = None,
    hypotheses: list[V4SurfaceHypothesis] | None = None,
) -> tuple[V4Panel, list[V4Candidate]]:
    """Select a panel; use footprint identifiability when available."""
    by_sector: dict[str, list[V4Candidate]] = defaultdict(list)
    for candidate in candidates:
        if candidate.connected:
            by_sector[candidate.sector_id].append(candidate)

    sector_series: dict[str, list[V4Candidate]] = {}
    per_sector = max(1, int(recommendations_per_sector))
    for sector_id, sector_candidates in sorted(by_sector.items()):
        prioritized = [
            item
            for item in sector_candidates
            if item.recommendation_class == "EPITOPE_PRIORITIZED"
        ]
        fallback = [
            item
            for item in sector_candidates
            if item.recommendation_class == "STABILITY_FALLBACK"
        ]
        series: list[V4Candidate] = []
        series.extend(
            _diverse_take(prioritized, min(per_sector, epitope_priority_per_sector), series)
        )
        series.extend(
            _diverse_take(
                fallback,
                min(per_sector - len(series), stability_fallback_per_sector),
                series,
            )
        )
        if len(series) < per_sector:
            series.extend(_diverse_take(sector_candidates, per_sector - len(series), series))
        sector_series[sector_id] = series[:per_sector]

    selected: list[V4Candidate] = []
    mean_response = mean_separation = None
    indistinguishable = None
    minimum_distance = None
    prior_weighted_separation = None
    unresolved_groups: list[list[str]] = []
    uncovered_hypotheses: list[str] = []
    hypothesis_signatures: dict[str, list[float]] = {}
    if hypotheses:
        (selected, mean_response, mean_separation, indistinguishable,
         minimum_distance, prior_weighted_separation, unresolved_groups, uncovered_hypotheses,
         hypothesis_signatures) = _information_take(
            [item for item in candidates if item.connected],
            hypotheses,
            max(1, int(max_primary)),
            max(0, int(higher_order_max)),
        )
    high_order_count = 0
    offset = 0
    while not hypotheses and len(selected) < max(1, int(max_primary)):
        added = False
        for sector_id in sorted(sector_series):
            series = sector_series[sector_id]
            if offset >= len(series):
                continue
            candidate = series[offset]
            if candidate.order > 2:
                if high_order_count >= max(0, int(higher_order_max)):
                    continue
                high_order_count += 1
            selected.append(candidate)
            added = True
            if len(selected) >= max(1, int(max_primary)):
                break
        if not added:
            break
        offset += 1

    roles: dict[str, list[str]] = {"PRIMARY": []}
    for candidate in selected:
        candidate.role = candidate.recommendation_class
        roles["PRIMARY"].append(candidate.candidate_id)
        roles.setdefault(candidate.recommendation_class, []).append(candidate.candidate_id)
    denominator = (
        {position for sector in sectors for position in sector.residue_positions}
        if sectors is not None
        else {
            position
            for candidate in candidates
            for positions in candidate.covered_positions_by_scenario.values()
            for position in positions
        }
    )
    coverage_positions: dict[str, set[int]] = defaultdict(set)
    for candidate in selected:
        for scenario, positions in candidate.covered_positions_by_scenario.items():
            coverage_positions[scenario].update(set(positions) & denominator)
    coverage = {
        scenario: round(len(positions) / len(denominator), 6) if denominator else 0.0
        for scenario, positions in sorted(coverage_positions.items())
    }
    areas = surface_areas or {}
    total_area = sum(max(0.0, areas.get(position, 0.0)) for position in denominator)
    area_coverage = {
        scenario: round(
            sum(max(0.0, areas.get(position, 0.0)) for position in positions) / total_area,
            6,
        ) if total_area else 0.0
        for scenario, positions in sorted(coverage_positions.items())
    }
    selected_sectors = {item.sector_id for item in selected}
    uncovered = [
        {"sector_id": sector_id, "reason": "NO_SELECTED_FEASIBLE_CANDIDATE"}
        for sector_id in sorted(by_sector)
        if sector_id not in selected_sectors
    ]
    panel = V4Panel(
        selected_candidate_ids=[item.candidate_id for item in selected],
        candidates_by_role=roles,
        achieved_coverage=coverage,
        uncovered_components=uncovered,
        solver_status="FOOTPRINT_IDENTIFIABILITY_GREEDY" if hypotheses else "SECTOR_SERIES_DETERMINISTIC",
        coverage_denominator_residues=len(denominator),
        covered_residues_by_scenario={
            scenario: len(positions) for scenario, positions in sorted(coverage_positions.items())
        },
        coverage_denominator_area_a2=round(total_area, 3),
        achieved_area_coverage=area_coverage,
        mean_hypothesis_response=round(mean_response, 6) if mean_response is not None else None,
        mean_pairwise_separation=round(mean_separation, 6) if mean_separation is not None else None,
        indistinguishable_hypothesis_pairs=indistinguishable,
        minimum_signature_distance=round(minimum_distance, 6) if minimum_distance is not None else None,
        unresolved_hypothesis_groups=unresolved_groups,
        uncovered_hypotheses=uncovered_hypotheses,
        hypothesis_signatures=hypothesis_signatures,
        prior_weighted_mean_pairwise_separation=(
            round(prior_weighted_separation, 6) if prior_weighted_separation is not None else None
        ),
    )
    return panel, selected
