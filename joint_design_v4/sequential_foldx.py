"""Bounded three-plus-two FoldX scheduling; original pair evidence is immutable."""
from __future__ import annotations
import math
import re
import statistics
from .models import PredictionResult
from .spatial_support import candidate_region_membership


def display_priority(rows):
    """Cover regional representatives before rotating through alternatives."""
    ids = [row.get("candidate_id") for row in rows]
    alternatives = [row.get("alternatives", []) for row in rows]
    for index in range(max((len(items) for items in alternatives), default=0)):
        ids.extend(items[index].get("candidate_id") for items in alternatives if index < len(items))
    return list(dict.fromkeys(candidate_id for candidate_id in ids if candidate_id))


def complete_display_pool(candidates):
    selected = []
    for candidate in candidates:
        repeat = (candidate.impact_evidence or {}).get("stability_repeat_audit", {})
        slots = repeat.get("foldx_repeat_slots") or []
        values = repeat.get("foldx_repeat_values_kcal_mol") or []
        try:
            complete = slots == list(range(5)) and len(values) == 5 and all(math.isfinite(float(v)) for v in values)
        except (TypeError, ValueError, OverflowError):
            complete = False
        if complete:
            selected.append(candidate)
    return selected


def _signature(row, position_map=None):
    mutations = row.get("mutations") or row.get("mutation_list") or []
    if isinstance(mutations, str):
        mutations = mutations.split(",")
    if position_map:
        mutations = [re.sub(r"\d+", lambda m: str(position_map.get(int(m[0]), int(m[0]))), str(x), count=1) for x in mutations]
    return tuple(sorted(str(x).strip() for x in mutations))


def merge_repeats(original, additional):
    """Merge distinct logical slots only; reject incomplete or conflicting evidence."""
    records = {_signature(row): dict(row) for row in original.records}
    for row in additional.records:
        key = _signature(row)
        if key not in records:
            raise ValueError("Supplemental FoldX returned an unrequested mutation identity")
        pairs = list(records[key].get("paired_models", [])) + list(row.get("paired_models", []))
        slots = [int(p["sampling_slot"]) for p in pairs]
        values = [float(p["ddg_kcal_mol"]) for p in pairs]
        if sorted(slots) != [0, 1, 2, 3, 4] or not all(math.isfinite(v) for v in values):
            raise ValueError("Supplemental FoldX requires five unique finite logical slots")
        pairs.sort(key=lambda p: int(p["sampling_slot"]))
        records[key].update(paired_models=pairs, ddg_runs_kcal_mol=[p["ddg_kcal_mol"] for p in pairs],
                            ddg=statistics.median(values), number_of_runs=5)
    return PredictionResult(status="success", records=list(records.values()),
        warnings=original.warnings + additional.warnings,
        provenance={**original.provenance, "supplemental_runs": [*original.provenance.get("supplemental_runs", []), additional.provenance]},
        runtime_seconds=original.runtime_seconds + additional.runtime_seconds)


def refine_foldx(candidates, sectors, initial, *, choose_panel, predict, recheck, position_map=None,
                 display_candidate_ids=()):
    available = {_signature(r, position_map) for r in initial.records}
    pool = [c for c in candidates if tuple(sorted(o.notation for o in c.mutations)) in available]
    if len(available) > 80:
        raise ValueError("Ordinary FoldX exceeds 80 unique combinations")
    audit = {"policy": "ordinary_v4_sequential_v1", "refinement_priority_policy": "display_representatives_first_v1",
        "unique_combinations": len(available),
        "initial_slots": len(available)*3, "supplemental_slots": 0, "attempted_refinement_candidate_ids": [], "batches": []}
    attempted = set()
    current_display_ids = []
    result = initial
    while len(attempted) < min(30, len(pool)):
        current_display_ids = list(dict.fromkeys(display_candidate_ids(candidates)
            if callable(display_candidate_ids) else display_candidate_ids))
        formal = [c for c in pool if c.formal_eligible]
        panel_ids = choose_panel(formal).selected_candidate_ids if formal else []
        regional_ids = []
        for sector in sectors:
            local = [c for c in formal if candidate_region_membership(c, sectors)["exclusive_region_id"] == sector.sector_id]
            if local:
                regional_ids.append(min(local, key=lambda c: ((c.impact_evidence or {}).get("selection_rank") or math.inf, c.candidate_id)).candidate_id)
        priority = list(dict.fromkeys(panel_ids + current_display_ids + regional_ids))
        remaining = [c for c in pool if c.candidate_id not in attempted]
        def key(c):
            values = (c.impact_evidence or {}).get("stability_repeat_audit", {}).get("foldx_repeat_values_kcal_mol", [])
            score = c.stability.foldx_ddg_kcal_mol if c.stability else None
            boundary = score is not None and (abs(score-2) <= .5 or (values and max(values)-min(values) > 1))
            return (0 if c.candidate_id in priority else 1 if boundary else 2,
                priority.index(c.candidate_id) if c.candidate_id in priority else 0,
                not c.formal_eligible, (c.impact_evidence or {}).get("selection_rank") or math.inf, c.candidate_id)
        batch = sorted(remaining, key=key)[:min(7, 30-len(attempted))]
        if not batch:
            break
        attempted.update(c.candidate_id for c in batch)
        audit["supplemental_slots"] += len(batch)*2
        audit["attempted_refinement_candidate_ids"].extend(c.candidate_id for c in batch)
        extra = predict(batch, len(audit["batches"])+1)
        audit["batches"].append({"candidate_ids": [c.candidate_id for c in batch], "status": extra.status, "errors": extra.errors})
        if extra.status == "success" and not extra.errors:
            expected = {tuple(sorted(o.notation for o in c.mutations)) for c in batch}
            if {_signature(r, position_map) for r in extra.records} != expected:
                raise ValueError("Supplemental FoldX response identity mismatch")
            result = merge_repeats(result, extra)
            recheck(result)
    audit["total_slots"] = audit["initial_slots"] + audit["supplemental_slots"]
    audit["display_candidate_ids"] = current_display_ids
    audit["display_candidates_not_refined"] = [candidate_id for candidate_id in current_display_ids
        if candidate_id not in attempted]
    pool_ids = {candidate.candidate_id for candidate in pool}
    audit["display_candidates_without_initial_foldx_evidence"] = [candidate_id for candidate_id in current_display_ids
        if candidate_id not in pool_ids]
    audit["display_candidates_blocked_by_refinement_budget"] = [candidate_id for candidate_id in current_display_ids
        if candidate_id in pool_ids and candidate_id not in attempted]
    audit["display_candidate_budget_status"] = (
        "INCOMPLETE_FIRST_PASS_SELECTION" if audit["display_candidates_without_initial_foldx_evidence"]
        else "INCOMPLETE_300_SLOT_CAP" if audit["display_candidates_blocked_by_refinement_budget"]
        else "COMPLETE")
    assert audit["total_slots"] <= 300
    result.provenance["sequential_budget"] = audit
    return result, audit
