"""Lazy, site-first queues for ordinary candidate evaluation."""
from __future__ import annotations

import hashlib
import itertools
import math
import time
from collections import defaultdict
from copy import copy, deepcopy
from dataclasses import replace

from .models import V4Candidate
from .candidate_search import _coverage_from_distances, _distances
from .regional_perturbation import (_candidate, _option, _position, _span_ok, candidate_impact,
                                    complete_candidate_impact, dominates, rank_impacts)
from .spatial_support import annotate_spatial_support, spatial_priority

POLICY = "ordinary_site_first_v2"
PRINCIPLES = ("epitope_support", "mechanism_pareto", "sector_coverage")


def _site_id(sector_id, positions):
    return f"{sector_id}|" + ",".join(map(str, positions))


def _identity(options, sector_id):
    sig = tuple(sorted((_position(option), option.mutant_aa) for option in options))
    return f"{sector_id}|V41_" + hashlib.sha256(repr(sig).encode()).hexdigest()[:16]


def _representatives(sector, options):
    """Rank actual single substitutions by the existing mechanism channel utilities."""
    roles = ("charge_polarity", "size_aromatic", "combined")
    candidates = [_candidate(sector, [option], 1) for option in options]
    ranked = rank_impacts(candidates, assess_sensitivity=False)
    channels = {candidate.mutations[0].mutant_aa: candidate.impact_evidence["ranking"]
                for candidate in ranked}
    by_role = {}
    for role in roles:
        def key(option):
            row = channels[option.mutant_aa]
            complete = -int(row.get("complete", False))
            if role == "combined":
                values = row.get("leximin_utility", ())
                unknown = not values
                utility = tuple(-float(value) for value in values)
            else:
                value = row.get("mechanism_group_utilities", {}).get(role)
                unknown = value is None
                utility = -float(value) if value is not None else math.inf
            risk = float(option.spurs_ddg_kcal_mol) if option.spurs_ddg_kcal_mol is not None else math.inf
            return (complete, unknown, utility, risk, option.mutant_aa)
        by_role[role] = sorted(options, key=key)
    selected, used = [], set()
    # Three distinct role witnesses when the evidence provides them.
    for role in roles:
        option = next((item for item in by_role[role] if item.mutant_aa not in used), None)
        if option is not None:
            selected.append((role, option))
            used.add(option.mutant_aa)
    initial = selected[:]
    # Fill to six by cycling each role's next available observed substitution.
    while len(selected) < 6:
        progressed = False
        for role in roles:
            option = next((item for item in by_role[role] if item.mutant_aa not in used), None)
            if option is not None:
                selected.append((role, option))
                used.add(option.mutant_aa)
                progressed = True
                if len(selected) == 6:
                    break
        if not progressed:
            break
    # A later seven-way product adds the next observed substitution only after
    # the first two physical tiers fail.
    expanded = selected[:6]
    while len(selected) < min(len(options), 7):
        progressed = False
        for role in roles:
            option = next((item for item in by_role[role] if item.mutant_aa not in used), None)
            if option is not None:
                selected.append((role, option))
                used.add(option.mutant_aa)
                progressed = True
                break
        if not progressed:
            break
    return initial, expanded, selected


def _principle_key(candidate, principle):
    if principle == "epitope_support":
        return (-candidate.epitope_supported_sites / max(1, candidate.order),
                -candidate.epitope_score / max(1, candidate.order), candidate.candidate_id)
    if principle == "mechanism_pareto":
        row = candidate.impact_evidence.get("ranking", {})
        return (-int(bool(row.get("complete"))), row.get("pareto_layer", math.inf),
                tuple(-float(x) for x in row.get("leximin_utility", ())), candidate.candidate_id)
    return spatial_priority(candidate)


def _site_epitope_key(positions, options_by_position):
    evidence = [options_by_position[position][0] for position in positions]
    priors = [option.epitope_prior for option in evidence]
    if any(option.epitope_evidence_status != "SCORED" or prior is None
           or not math.isfinite(float(prior)) for option, prior in zip(evidence, priors)):
        return None
    return (-sum(bool(option.epitope_supported) for option in evidence) / len(evidence),
            -sum(float(prior) for prior in priors) / len(priors))


def _geometry_key(geometry):
    fields = ("target_mean_coverage_fraction", "mean_marginal_area_a2",
              "mean_overlap_fraction", "direct_area_a2")
    values = [geometry.get(name) for name in fields]
    if any(value is None or not math.isfinite(float(value)) for value in values):
        return None
    return (-float(values[0]), -float(values[1]), float(values[2]), -float(values[3]))


def _triple_precedes_parents(triple, parents, principle, gain):
    """Use the strongest comparable exact double parent, never a weak surrogate."""
    if gain is None or gain < .05 or not parents:
        return False
    if principle == "mechanism_pareto":
        best_parent = rank_impacts(parents, assess_sensitivity=False)[0]
        return dominates(candidate_impact(triple), candidate_impact(best_parent))
    if principle == "sector_coverage":
        return float(triple.geometry_support.get("target_mean_coverage_fraction") or 0.0) >= max(
            float(parent.geometry_support.get("target_mean_coverage_fraction") or 0.0) for parent in parents) + .05
    if any(option.epitope_evidence_status != "SCORED" or option.epitope_prior is None
           for option in triple.mutations):
        return False
    if any(option.epitope_evidence_status != "SCORED" or option.epitope_prior is None
           for parent in parents for option in parent.mutations):
        return False
    return _principle_key(triple, principle) < min(_principle_key(parent, principle) for parent in parents)


def _exact_parent_mutations(triple):
    """Return the three exact AA-preserving double parents of a triple."""
    return [[option for option in triple.mutations if option is not omitted]
            for omitted in triple.mutations]


class SiteFirstQueue:
    def __init__(self):
        self._options = {}
        self._sectors = {}
        self._lanes = {}
        self._candidate = {}
        self._seen_candidate_ids = set()
        self._geometry = {}
        self._sequence = None
        self._closed = False
        self._ledger = {"schema_version": "ordinary_site_first_ledger_v3", "policy": POLICY,
                        "representative_policy": "role_ladder_3_6_7_v1",
                        "site_set_count": 0, "geometry_evaluations": 0,
                        "expanded_candidate_count": 0, "deep_candidate_count": 0,
                        "initial_candidate_count": 0,
                        "site_preparation_seconds": 0.0, "representative_expansion_seconds": 0.0,
                        "uncompleted_gaps": []}

    def close(self):
        self._closed = True

    def queue_states(self):
        return sorted(self._lanes)

    def _lane(self, sector, principle):
        try:
            return self._lanes[(str(sector), principle)]
        except KeyError:
            raise KeyError(f"Unknown site-first lane: {sector}/{principle}") from None

    def status(self, sector, principle):
        lane = self._lane(sector, principle)
        if lane["state"] != "ready":
            return {"status": lane["state"]}
        group = lane["groups"][lane["group_index"]] if lane["group_index"] < len(lane["groups"]) else None
        if group is not None and group["items"] is None:
            return {"status": "ready", "pending": None, "site_set_id": group["site_set_id"]}
        return {"status": "ready" if group and group["items"] else "exhausted",
                "pending": sum(x["state"] == "pending" for x in (group["items"] if group else [])),
                "site_set_id": group["site_set_id"] if group else None}

    def head(self, sector, principle, limit=4):
        if limit < 0:
            raise ValueError("limit must be nonnegative")
        if limit == 0:
            return []
        lane = self._lane(sector, principle)
        if lane["state"] != "ready" or lane["group_index"] >= len(lane["groups"]):
            return []
        self._skip_preblocked(sector, principle, lane)
        if lane["state"] != "ready" or lane["group_index"] >= len(lane["groups"]):
            return []
        group = lane["groups"][lane["group_index"]]
        if group["items"] is None:
            self._set_group_tier(group, principle, "initial")
            self._advance_empty_group_tier(group, principle)
        while not group["items"] and lane["group_index"] + 1 < len(lane["groups"]):
            lane["group_index"] += 1
            self._skip_preblocked(sector, principle, lane)
            if lane["state"] != "ready" or lane["group_index"] >= len(lane["groups"]):
                return []
            group = lane["groups"][lane["group_index"]]
            self._set_group_tier(group, principle, "initial")
            self._advance_empty_group_tier(group, principle)
        if not group["items"]:
            lane["state"] = "candidate_failures"
            self._ledger["uncompleted_gaps"].append({"sector_id": str(sector), "principle": principle,
                                                       "status": lane["state"]})
            return []
        return [self.materialize(row["id"]) for row in group["items"] if row["state"] == "pending"][:limit]

    def prefetch_groups(self, sector, principle):
        """Yield ranked pending candidates without changing the live queue."""
        lane = self._lane(sector, principle)
        if lane["state"] != "ready":
            return
        shadow = copy(self)
        shadow._candidate = {}
        shadow._seen_candidate_ids = set()
        shadow._ledger = deepcopy(self._ledger)
        for index in range(lane["group_index"], len(lane["groups"])):
            group = lane["groups"][index]
            if ((principle == "epitope_support" and group["site_priority"] is None)
                    or (group["order"] == 3 and not group["triple_gain_gate_passed"])):
                continue
            if index == lane["group_index"] and group["items"] is not None:
                rows = [deepcopy(self._candidate[row["id"]]) for row in group["items"]
                        if row["state"] == "pending" and row["id"] in self._candidate]
            else:
                temporary = deepcopy(group)
                temporary["items"] = None
                shadow._candidate.clear()
                items = shadow._materialize_group(temporary, principle, tier="initial")
                rows = [shadow._candidate[row["id"]] for row in items]
            if rows:
                yield rows

    def _skip_preblocked(self, sector, principle, lane):
        while lane["group_index"] < len(lane["groups"]):
            group = lane["groups"][lane["group_index"]]
            missing_epitope = principle == "epitope_support" and group["site_priority"] is None
            failed_gain = group["order"] == 3 and not group["triple_gain_gate_passed"]
            if not missing_epitope and not failed_gain:
                break
            if not group.get("preblocked_recorded"):
                self._ledger["uncompleted_gaps"].append({"sector_id": str(sector), "principle": principle,
                    "site_set_id": group["site_set_id"], "status": "epitope_prior_missing" if missing_epitope
                    else "regional_gain_below_or_unknown_threshold",
                    "measured_gain": group["gain"], "required_gain": .05})
                group["preblocked_recorded"] = True
            if missing_epitope:
                lane["missing_evidence"] = True
            lane["group_index"] += 1
        if lane["group_index"] >= len(lane["groups"]):
            lane["state"] = "missing_evidence" if lane.get("missing_evidence") else "candidate_failures"

    def queue_candidate_id(self, sector, principle, offset=0):
        if offset < 0:
            return None
        candidates = self.head(sector, principle, limit=offset + 1)
        return candidates[offset].candidate_id if len(candidates) > offset else None

    def prioritize_uncovered(self, covered_positions):
        """Reserve the next attempt for a new residue before retrying covered sites."""
        for (_, principle), lane in self._lanes.items():
            index = lane["group_index"]
            groups = lane["groups"]
            if lane["state"] != "ready" or index >= len(groups):
                continue
            current = groups[index]
            if not set(current["positions"]).issubset(covered_positions):
                continue
            for later in range(index + 1, len(groups)):
                group = groups[later]
                if (set(group["positions"]).issubset(covered_positions)
                        or (group["order"] == 3 and not group["triple_gain_gate_passed"])
                        or (principle == "epitope_support" and group["site_priority"] is None)):
                    continue
                # Suspended pending rows are regenerated from evidence when resumed.
                current["history"].extend(row for row in (current["items"] or [])
                                          if row["state"] != "pending")
                current["items"] = None
                current["tier"] = "initial"
                groups.insert(index, groups.pop(later))
                break

    def needs_exploration(self, covered_positions):
        return any(
            not set(group["positions"]).issubset(covered_positions)
            and (group["order"] != 3 or group["triple_gain_gate_passed"])
            and (principle != "epitope_support" or group["site_priority"] is not None)
            for (_, principle), lane in self._lanes.items() if lane["state"] == "ready"
            for group in lane["groups"][lane["group_index"]:]
        )

    def materialize(self, candidate_id):
        if candidate_id not in self._candidate:
            raise KeyError(candidate_id)
        candidate = deepcopy(self._candidate[candidate_id])
        return candidate

    def complete(self, sector, principle, candidate_id):
        lane, row = self._find(sector, principle, candidate_id)
        if row["state"] != "pending":
            raise KeyError(candidate_id)
        row["state"] = "complete"
        lane["state"] = "complete"
        self._prune_candidates()

    def record_outcome(self, sector, principle, candidate_id, outcome):
        lane, row = self._find(sector, principle, candidate_id)
        if row["state"] != "pending":
            raise KeyError(candidate_id)
        row["state"] = str(outcome)
        group = lane["groups"][lane["group_index"]]
        if all(item["state"] != "pending" for item in group["items"]):
            group["history"].extend(group["items"])
            # Failures open the next tier, then the next physical site set.
            if not self._advance_group_tier(group, principle):
                lane["group_index"] += 1
                self._prune_candidates()
                if lane["group_index"] >= len(lane["groups"]):
                    lane["state"] = "candidate_failures"
                    self._ledger["uncompleted_gaps"].append({"sector_id": str(sector), "principle": principle,
                                                               "status": lane["state"]})

    def _set_group_tier(self, group, principle, tier):
        group["tier"] = tier
        group["items"] = self._materialize_group(group, principle, tier=tier)
        group["expanded"] = {"initial": 0, "expanded": 1, "deep": 2}[tier]
        if tier == "expanded":
            self._ledger["expanded_candidate_count"] += len(group["items"])
        elif tier == "deep":
            self._ledger["deep_candidate_count"] += len(group["items"])

    def _advance_group_tier(self, group, principle):
        while True:
            current = group.get("tier", "initial")
            if current == "initial" and group["expansion"]:
                next_tier = "expanded"
            elif current == "expanded" and group.get("deep_options"):
                next_tier = "deep"
            else:
                return False
            self._set_group_tier(group, principle, next_tier)
            if group["items"]:
                self._prune_candidates()
                return True

    def _advance_empty_group_tier(self, group, principle):
        if group["items"]:
            return
        group["history"].extend(group["items"])
        self._advance_group_tier(group, principle)

    def _materialize_group(self, group, principle, *, tier):
        started = time.perf_counter()
        options = group[tier + "_options"]
        combos = itertools.product(*options)
        old = {row["id"] for row in group["history"] + (group["items"] or [])}
        rows = []
        for combo in combos:
            selected = [row[1] for row in combo]
            candidate = _candidate(self._sectors[group["sector_id"]], selected, len(selected))
            candidate.candidate_id = _identity(selected, group["sector_id"])
            if candidate.candidate_id in old:
                continue
            from .coverage_profile import sequence_reasons
            if sequence_reasons(self._sequence, [(o.residue.uniprot_pos, o.residue.one_letter, o.mutant_aa)
                                                  for o in selected]):
                continue
            candidate.geometry_support = deepcopy(group["geometry"])
            candidate.impact_evidence["regional_geometry_support"] = deepcopy(group["geometry"])
            allowed = set(self._sectors[group["sector_id"]].residue_positions)
            candidate.coverage_by_scenario, candidate.covered_positions_by_scenario = _coverage_from_distances(
                [_distances(self._graph, pos, allowed) for pos in group["positions"]], allowed)
            candidate.risk_components = {
                "regional_higher_order_gain": group["gain"] if candidate.order == 3 else None,
                "regional_lower_order_coverage": group.get("baseline_coverage"),
                "regional_candidate_coverage": group["geometry"].get("target_mean_coverage_fraction"),
                "regional_gain_baseline_candidate_id": group.get("baseline_site_set_id"),
                "regional_gain_region_ids": [group["sector_id"]],
                "regional_gain_basis": "mean_regional_area_coverage_9_12_15",
                "regional_gain_status": ("COMPUTED" if group.get("gain") is not None else "MISSING_GEOMETRY_GAIN")
                    if candidate.order == 3 else "NOT_APPLICABLE",
            }
            if candidate.order == 3:
                candidate.parent_candidate_ids = sorted(_identity([o for o in selected if o is not omitted],
                                                                    group["sector_id"])
                                                        for omitted in selected)
            candidate.impact_evidence.update({"protocol": POLICY, "site_set_id": group["site_set_id"],
                "source_site_set_id": group["site_set_id"], "representative_tier": tier,
                "approximate_pre_ranking": {"regional_gain": group["gain"], "geometry": group["geometry"]},
                "representative_roles": [row[0] for row in combo]})
            candidate.epitope_score = sum(float(getattr(o, "epitope_prior", 0.0) or 0.0) for o in selected)
            candidate.epitope_supported_sites = sum(bool(getattr(o, "epitope_supported", False)) for o in selected)
            candidate.impact_evidence["epitope_evidence_status"] = (
                "SCORED" if all(o.epitope_evidence_status == "SCORED" and o.epitope_prior is not None
                                for o in selected) else "MISSING")
            complete_candidate_impact([candidate], self._records, graph=self._graph)
            self._candidate[candidate.candidate_id] = candidate
            self._seen_candidate_ids.add(candidate.candidate_id)
            rows.append({"id": candidate.candidate_id, "state": "pending"})
        ranked = [self._candidate[row["id"]] for row in rows]
        if principle == "mechanism_pareto":
            ranked = rank_impacts(ranked, assess_sensitivity=False)
        else:
            ranked.sort(key=lambda c: _principle_key(c, principle))
        order = {c.candidate_id: i for i, c in enumerate(ranked)}
        rows.sort(key=lambda row: order[row["id"]])
        if tier == "initial" and not group.get("initial_materialized"):
            self._ledger["initial_candidate_count"] += len(rows)
            group["initial_materialized"] = True
        self._ledger["representative_expansion_seconds"] += time.perf_counter() - started
        return rows

    def _prune_candidates(self):
        active = {row["id"] for lane in self._lanes.values() if lane["state"] == "ready"
                  and lane["group_index"] < len(lane["groups"])
                  for row in (lane["groups"][lane["group_index"]]["items"] or [])
                  if row["state"] == "pending"}
        self._candidate = {key: value for key, value in self._candidate.items() if key in active}

    def mark_all_technical_failure(self):
        for lane in self._lanes.values():
            if lane["state"] == "ready":
                lane["state"] = "technical_failure"

    def mark_all_budget_limited(self, limit):
        for lane in self._lanes.values():
            if lane["state"] == "ready":
                lane["state"] = "budget_limited"
                lane["budget_limit"] = int(limit)

    def mark_queue_state(self, sector, principle, state, limit=None):
        lane = self._lane(sector, principle)
        if lane["state"] in {"complete", "not_applicable", "missing_evidence", "candidate_failures"}:
            return
        if state not in {"candidate_failures", "budget_limited", "technical_failure", "not_applicable"}:
            raise ValueError(state)
        lane["state"] = state
        if limit is not None:
            lane["budget_limit"] = int(limit)

    def ledger(self, validation=None):
        result = deepcopy(self._ledger)
        result["counts"] = {"lanes": len(self._lanes), "legal_candidates": len(self._seen_candidate_ids),
                             "resident_candidates": len(self._candidate),
                             "legal_candidate_count_scope": "materialized_initial_expanded_or_deep",
                             "ready_lanes": sum(x["state"] == "ready" for x in self._lanes.values()),
                             "completed_lanes": sum(x["state"] == "complete" for x in self._lanes.values())}
        rows = []
        for sector, principle in self.queue_states():
            lane = self._lanes[(sector, principle)]
            completed = [{"candidate_id": item["id"], "outcome": item["state"]}
                         for group in lane["groups"] for item in group["history"] + (group["items"] or [])
                         if item["state"] != "pending"]
            state = lane["state"]
            if state == "ready" and lane["group_index"] >= len(lane["groups"]):
                state = "exhausted"
            rows.append({"sector_id": sector, "principle": principle, "state": state,
                         "pending_count": self.status(sector, principle).get("pending"),
                         "site_set_id": self.status(sector, principle).get("site_set_id"),
                         "next_candidate_id": (lane["groups"][lane["group_index"]]["items"][0]["id"]
                                               if state == "ready" and lane["group_index"] < len(lane["groups"])
                                               and lane["groups"][lane["group_index"]]["items"] else None),
                         "completed": completed})
        result["queues"] = rows
        existing_gaps = {(gap.get("sector_id"), gap.get("principle"), gap.get("status"))
                         for gap in result["uncompleted_gaps"]}
        for row in rows:
            state = row["state"]
            if state not in {"candidate_failures", "budget_limited", "technical_failure", "missing_evidence", "exhausted"}:
                continue
            lane = self._lanes[(row["sector_id"], row["principle"])]
            status = "no_legal_candidates" if state == "candidate_failures" and not lane["groups"] else state
            key = (row["sector_id"], row["principle"], status)
            if key not in existing_gaps:
                gap = {"sector_id": row["sector_id"], "principle": row["principle"], "status": status}
                if lane.get("budget_limit") is not None:
                    gap["budget_limit"] = lane["budget_limit"]
                result["uncompleted_gaps"].append(gap)
                existing_gaps.add(key)
        result["principles"] = list(PRINCIPLES)
        result["spurs"] = dict(validation.spurs_ledger) if validation is not None else None
        result["foldx"] = dict(validation.foldx_ledger) if validation is not None else None
        result["status_counts"] = {
            "no_legal_candidates": sum(not lane["groups"] for lane in self._lanes.values()),
            "candidate_failures": sum(row["state"] == "candidate_failures" for row in rows),
            "not_yet_scored": sum(row["pending_count"] or 0 for row in rows),
            "pending_count_unknown_lanes": sum(row["pending_count"] is None and row["state"] == "ready" for row in rows),
            "technical_failures": sum(row["state"] == "technical_failure" for row in rows),
            "missing_evidence": sum(row["state"] == "missing_evidence" for row in rows),
        }
        result["timings"] = {"site_preparation_seconds": result["site_preparation_seconds"],
                              "representative_expansion_seconds": result["representative_expansion_seconds"]}
        if validation is not None:
            result["evaluated_candidates"] = [
                {"candidate_id": value.get("candidate_id"), "mutations": list(signature),
                 **{key: item for key, item in value.items() if key not in {"candidate_id", "mutations", "row", "audit"}}}
                for signature, value in sorted(validation.spurs_results.items())]
            result["foldx_results"] = [
                {"candidate_id": value.get("candidate_id"), "mutations": list(signature),
                 **{key: item for key, item in value.items() if key not in {"candidate_id", "mutations", "audit"}}}
                for signature, value in sorted(validation.foldx_results.items())]
        return result

    def _find(self, sector, principle, candidate_id):
        lane = self._lane(sector, principle)
        for group in lane["groups"]:
            for row in group["history"] + (group["items"] or []):
                if row["id"] == candidate_id:
                    return lane, row
        raise KeyError(candidate_id)


def build_site_first_queues(single, sectors, graph, areas, canonical_sequence, *, maximum_span=24.0,
                            records=None, epitope_applicable=True):
    """Build all legal local site sets and only their initial representative layer."""
    from .coverage_profile import sequence_reasons
    store = SiteFirstQueue()
    store._sequence = canonical_sequence
    store._graph = graph
    prepare_start = time.perf_counter()
    options = defaultdict(list)
    records_by_position = (dict(records) if isinstance(records, dict)
                           else {record.uniprot_pos: record for record in (records or [])})
    for item in single:
        pos = getattr(item.residue, "uniprot_pos", None)
        score = getattr(item, "spurs_ddg_kcal_mol", None)
        if (getattr(item, "eligible", False) and not getattr(item, "review_required", False)
                and not getattr(item, "hard_exclusions", ()) and pos in graph
                and score is not None and math.isfinite(float(score))):
            option = _option(item)
            option.epitope_prior = getattr(item, "epitope_prior", None)
            option.epitope_supported = bool(getattr(item, "epitope_supported", False))
            option.epitope_evidence_status = getattr(item, "epitope_evidence_status", "MISSING")
            record = records_by_position.get(pos)
            if record is not None:
                option.residue = replace(option.residue, one_letter=record.aa,
                    auth_seq_id=record.auth_seq_id, insertion_code=record.insertion_code or "",
                    label_seq_id=getattr(record, "label_seq_id", None))
            options[int(pos)].append(option)
    store._options = options
    store._records = records_by_position
    store._sectors = {str(s.sector_id): s for s in sectors}
    for sector in sorted(sectors, key=lambda s: str(s.sector_id)):
        allowed = set(sector.residue_positions) & set(graph)
        local = sorted(allowed & set(options))
        representatives = {p: _representatives(sector, options[p]) for p in local}
        site_sets = {}
        for order in (2, 3):
            for positions in itertools.combinations(local, order):
                ok, _span = _span_ok(positions, graph, allowed, maximum_span)
                if not ok:
                    continue
                # A site set is geometrically evaluated exactly once, independent of AA choice.
                pseudo = [_candidate(sector, [options[p][0] for p in positions], order)]
                if areas:
                    annotate_spatial_support(pseudo, graph, areas, target_positions=allowed, scope_id=sector.sector_id)
                geom = deepcopy(pseudo[0].geometry_support or {})
                store._geometry[(str(sector.sector_id), positions)] = geom
                store._ledger["geometry_evaluations"] += 1
                store._ledger["site_set_count"] += 1
                gain = None
                baseline_positions = None
                baseline_coverage = None
                if order == 3:
                    triple = geom.get("target_mean_coverage_fraction")
                    parents = [(positions, value.get("target_mean_coverage_fraction"))
                               for (sid, positions), value in store._geometry.items()
                               if sid == str(sector.sector_id) and len(positions) == 2]
                    parents = [(positions, float(value)) for positions, value in parents if value is not None]
                    if triple is not None and parents:
                        baseline_positions, baseline_coverage = max(parents, key=lambda row: row[1])
                        gain = float(triple) - baseline_coverage
                site_sets[positions] = {"site_set_id": _site_id(sector.sector_id, positions), "positions": positions,
                                        "sector_id": str(sector.sector_id), "order": order, "geometry": geom, "gain": gain,
                                        "triple_gain_gate_passed": order == 2 or (gain is not None and gain >= .05),
                                        "baseline_site_set_id": _site_id(sector.sector_id, baseline_positions) if baseline_positions else None,
                                        "baseline_coverage": baseline_coverage,
                                        "initial_options": [representatives[p][0] for p in positions],
                                        "expanded_options": [representatives[p][1] for p in positions],
                                        "deep_options": [representatives[p][2] for p in positions]}
        # Every principle gets an independent winner queue, while physical site sets and geometry are shared.
        for principle in PRINCIPLES:
            groups = []
            groups_by_positions = {}
            site_witnesses = {}
            all_mechanism_witnesses = []
            for positions, definition in site_sets.items():
                witnesses = []
                roles = ("charge_polarity", "size_aromatic", "combined")
                for index, role in enumerate(roles):
                    selected = [definition["initial_options"][p][index][1]
                                for p in range(len(positions)) if index < len(definition["initial_options"][p])]
                    if len(selected) != len(positions):
                        continue
                    if sequence_reasons(canonical_sequence, [(o.residue.uniprot_pos, o.residue.one_letter, o.mutant_aa) for o in selected]):
                        continue
                    candidate = _candidate(sector, selected, len(selected))
                    candidate.candidate_id = _identity(selected, sector.sector_id)
                    candidate.geometry_support = deepcopy(definition["geometry"])
                    candidate.impact_evidence.update({"protocol": POLICY, "site_set_id": definition["site_set_id"],
                        "source_site_set_id": definition["site_set_id"], "representative_tier": "site_witness",
                        "approximate_pre_ranking": {"regional_gain": definition["gain"], "geometry": definition["geometry"]},
                        "representative_roles": [role] * len(selected)})
                    candidate.epitope_score = sum(float(getattr(o, "epitope_prior", 0.0) or 0.0) for o in selected)
                    candidate.epitope_supported_sites = sum(bool(getattr(o, "epitope_supported", False)) for o in selected)
                    candidate.impact_evidence["epitope_evidence_status"] = (
                        "SCORED" if all(o.epitope_evidence_status == "SCORED" and o.epitope_prior is not None
                                        for o in selected) else "MISSING")
                    witnesses.append(candidate)
                if principle != "mechanism_pareto":
                    witnesses.sort(key=lambda c: _principle_key(c, principle))
                if witnesses:
                    site_witnesses[positions] = witnesses
                    all_mechanism_witnesses.extend(witnesses) if principle == "mechanism_pareto" else None
                if principle == "epitope_support":
                    definition["site_priority"] = _site_epitope_key(positions, options)
                    if definition["site_priority"] is None and epitope_applicable:
                        store._ledger["uncompleted_gaps"].append({"sector_id": str(sector.sector_id),
                            "principle": principle, "site_set_id": definition["site_set_id"],
                            "status": "epitope_prior_missing"})
                elif principle == "sector_coverage":
                    definition["site_priority"] = _geometry_key(definition["geometry"])
                else:
                    definition["site_priority"] = _principle_key(witnesses[0], principle) if witnesses else None
                definition["witness_count"] = len(witnesses)
                group = {**definition, "sector_id": str(sector.sector_id), "items": None,
                         "history": [], "expansion": True, "expanded": 0, "tier": "initial"}
                groups.append(group)
                groups_by_positions[positions] = group
            if principle == "mechanism_pareto" and all_mechanism_witnesses:
                global_order = rank_impacts(all_mechanism_witnesses, assess_sensitivity=False)
                ordinal = {id(candidate): index for index, candidate in enumerate(global_order)}
                for positions, witnesses in site_witnesses.items():
                    witnesses.sort(key=lambda candidate: ordinal[id(candidate)])
                    groups_by_positions[positions]["site_priority"] = _principle_key(witnesses[0], principle)
            # Pairs are visited before triples; among triples, geometry gain is descending.
            by_positions = {tuple(g["positions"]): g for g in groups}
            for group in groups:
                group["promoted"] = False
                if group["order"] != 3 or not group["triple_gain_gate_passed"] or not group["witness_count"]:
                    continue
                parents = [by_positions.get(pair) for pair in itertools.combinations(group["positions"], 2)]
                parents = [parent for parent in parents if parent]
                if parents and group["witness_count"]:
                    triple = site_witnesses[tuple(group["positions"])][0]
                    parent_candidates = []
                    for mutations in _exact_parent_mutations(triple):
                        parent_positions = tuple(sorted(_position(option) for option in mutations))
                        definition = site_sets.get(parent_positions)
                        if definition is None or sequence_reasons(canonical_sequence, [
                                (_position(option), option.residue.one_letter, option.mutant_aa)
                                for option in mutations]):
                            continue
                        parent = _candidate(sector, mutations, 2)
                        parent.candidate_id = _identity(mutations, sector.sector_id)
                        parent.geometry_support = deepcopy(definition["geometry"])
                        parent.epitope_score = sum(float(getattr(option, "epitope_prior", 0.0) or 0.0)
                                                   for option in mutations)
                        parent.epitope_supported_sites = sum(bool(option.epitope_supported) for option in mutations)
                        parent_candidates.append(parent)
                    if not parent_candidates:
                        continue
                    group["promoted"] = _triple_precedes_parents(
                        triple, parent_candidates, principle, group["gain"])
            groups.sort(key=lambda g: (not g["promoted"], g["order"], g["site_priority"] is None,
                                       g["site_priority"] or (),
                                       -(g["gain"] or 0), g["site_set_id"]))
            not_applicable = principle == "epitope_support" and not epitope_applicable
            store._lanes[(str(sector.sector_id), principle)] = {"groups": [] if not_applicable else groups, "group_index": 0,
                "state": "not_applicable" if not_applicable else "ready" if groups else "candidate_failures"}
    store._ledger["site_preparation_seconds"] = time.perf_counter() - prepare_start
    store._ledger["representative_expansion_seconds"] = 0.0
    return store
