"""Regional mutation combinations ranked by transparent chemical proxies.

This policy describes nominal side-chain changes and solvent exposure. It does
not estimate binding energy, affinity, or biological effect.
"""
from __future__ import annotations

from itertools import combinations, product
import hashlib
import math
from typing import Any

from .models import V4Candidate, V4MutationOption

POLICY = "regional_mechanism_pareto_v1"
RANKING_VERSION = "mechanism_channels_v2"
AXES = ("charge_delta", "volume_delta_a3", "hydropathy_delta",
        "donor_delta", "acceptor_delta", "aromatic_delta")
GROUPS = {
    "exposure": ("exposure_min",),
    "charge_polarity": ("charge_delta", "hydropathy_delta", "donor_delta", "acceptor_delta"),
    "size_aromatic": ("volume_delta_a3", "aromatic_delta"),
}
GROUP_NAMES = tuple(GROUPS)
DOUBLE_CAP = 20
TRIPLE_CAP = 10
DOUBLE_PARENT_CAP = 12


def _finite(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _mutations(candidate):
    return list(getattr(candidate, "mutations", ()) or ())


def _position(option):
    value = getattr(getattr(option, "residue", None), "uniprot_pos", None)
    return int(value) if value is not None else None


def _descriptors(option):
    from .coverage_profile import descriptors
    item = getattr(option, "chemical_descriptors", None)
    if isinstance(item, dict) and any(axis in item for axis in AXES):
        return item
    evidence = getattr(option, "impact_evidence", {}) or {}
    item = evidence.get("chemical_descriptors")
    if isinstance(item, dict) and any(axis in item for axis in AXES):
        return item
    residue = getattr(option, "residue", None)
    wt = getattr(residue, "one_letter", "X")
    mutant = getattr(option, "mutant_aa", "X")
    try:
        return descriptors(wt, mutant)
    except (KeyError, TypeError):
        return {axis: None for axis in AXES}


def candidate_impact(candidate) -> dict[str, Any]:
    """Return raw, auditable per-candidate impact axes without filling missing data."""
    mutations = _mutations(candidate)
    exposure_values, raw_mutations = [], []
    axis_values = {axis: [] for axis in AXES}
    missing = set()
    for option in mutations:
        evidence = getattr(option, "impact_evidence", {}) or {}
        exposure = evidence.get("sidechain_exposure", {})
        exposure_value = _finite(exposure.get("value")) if isinstance(exposure, dict) else None
        if exposure_value is None:
            missing.add("sidechain_exposure")
        else:
            exposure_values.append(exposure_value)
        descriptors = _descriptors(option)
        row = {"position": _position(option), "wt": getattr(getattr(option, "residue", None), "one_letter", None),
               "mutant": getattr(option, "mutant_aa", None), "sidechain_exposure": exposure,
               "chemical_descriptors": {axis: descriptors.get(axis) for axis in AXES}}
        raw_mutations.append(row)
        for axis in AXES:
            value = _finite(descriptors.get(axis))
            if value is None:
                missing.add(axis)
            else:
                axis_values[axis].append(abs(value))
    means = {axis: (sum(values) / len(values) if values else None)
             for axis, values in axis_values.items()}
    for axis, values in axis_values.items():
        if len(values) != len(mutations):
            means[axis] = None
    signed = {}
    for axis in ("charge_delta", "hydropathy_delta"):
        values = [_finite(row["chemical_descriptors"].get(axis)) for row in raw_mutations]
        known = bool(values) and all(value is not None for value in values)
        total = sum(values) if known else None
        absolute = sum(abs(value) for value in values) if known else None
        signed[axis] = {"status": "UNKNOWN" if not known else "NO_CHANGE" if absolute == 0 else "KNOWN",
                        "net_change": total, "absolute_sum": absolute,
                        "opposition_fraction": 1 - abs(total) / absolute if absolute else None}
    return {
        "policy": POLICY,
        "ranking_version": RANKING_VERSION,
        "order": len(mutations),
        "exposure_min": min(exposure_values) if len(exposure_values) == len(mutations) and exposure_values else None,
        "exposure_mean": sum(exposure_values) / len(exposure_values) if len(exposure_values) == len(mutations) and exposure_values else None,
        "mean_absolute_delta": means,
        "axes": dict(means),
        **means,
        "missing_mask": sorted(missing),
        "missing_axes": sorted(missing),
        "complete": not missing and bool(mutations),
        "signed_change_summary": signed,
        "aromatic_mutant_count": sum(row["mutant"] in "FWY" for row in raw_mutations),
        "introduced_aromatic_count": sum(row["mutant"] in "FWY" and row["wt"] not in "FWY" for row in raw_mutations),
        "mutations": raw_mutations,
    }


def _record(candidate):
    evidence = getattr(candidate, "impact_evidence", {}) or {}
    impact = evidence.get("regional_impact")
    return impact if isinstance(impact, dict) else candidate_impact(candidate)


def annotate_local_interactions(candidates, records):
    """WT C-alpha context only: neither mutant contacts nor a ranking penalty."""
    coordinates = {record.uniprot_pos: getattr(record, "ca", None) for record in records}
    for candidate in candidates:
        pairs = []
        for left, right in combinations(_mutations(candidate), 2):
            a, b = coordinates.get(_position(left)), coordinates.get(_position(right))
            valid = all(point is not None and len(point) == 3
                        and all(_finite(value) is not None for value in point) for point in (a, b))
            distance = math.dist(a, b) if valid else None
            changes = [_descriptors(option) for option in (left, right)]
            opposed = {}
            for axis in ("charge_delta", "hydropathy_delta"):
                values = [_finite(change.get(axis)) for change in changes]
                opposed[axis] = (values[0] * values[1] < 0
                                 if all(value is not None for value in values) else None)
            pairs.append({"positions": [_position(left), _position(right)], "wt_ca_distance_a": distance,
                          "opposed_changes": opposed,
                          "nearby_opposed_charge": (opposed["charge_delta"] and distance < 10
                              if opposed["charge_delta"] is not None and distance is not None else None),
                          "nearby_opposed_hydropathy": (opposed["hydropathy_delta"] and distance < 10
                              if opposed["hydropathy_delta"] is not None and distance is not None else None),
                          "ca_below_6_a": distance < 6 if distance is not None else None,
                          "ca_above_25_a": distance > 25 if distance is not None else None,
                          "aromatic_mutant_pair": all(option.mutant_aa in "FWY" for option in (left, right))})
        candidate.impact_evidence["local_interaction_diagnostics"] = {
            "basis": "WT_CA_COORDINATES", "review_only": True, "pairs": pairs,
            "nearby_ca_cutoff_a": 10.0,
            "interpretation": "WT CA distance is not CB distance, mutant side-chain contact, salt-bridge formation, steric clash or affinity evidence; signed net zero does not imply physical cancellation."}


def complete_candidate_impact(candidates, records, graph=None):
    """Attach descriptive chemical and structural context to every candidate.

    This is independent of regional mechanism ranking and never changes gates,
    eligibility, or recommendation roles. Unknown descriptors/exposure remain
    explicitly unknown rather than being imputed as zero.
    """
    from .coverage_profile import descriptors

    by_position = {record.uniprot_pos: record for record in (records.values() if isinstance(records, dict) else records or [])}
    for candidate in candidates:
        candidate.impact_evidence = dict(getattr(candidate, "impact_evidence", None) or {})
        for option in _mutations(candidate):
            position = _position(option)
            record = by_position.get(position)
            evidence = dict(getattr(option, "impact_evidence", None) or {})
            wt = getattr(getattr(option, "residue", None), "one_letter", None)
            mutant = getattr(option, "mutant_aa", None)
            try:
                chemical = descriptors(wt, mutant)
            except (KeyError, TypeError):
                chemical = {"version": "sidechain_descriptors_v1", "reference": "stability.AA_PROPERTIES; nominal sidechain functional-group counts",
                            **{axis: None for axis in AXES}, "protonation_status": "UNKNOWN_INVALID_RESIDUE"}
            exposure = getattr(record, "sidechain_exposure", None) if record is not None else None
            if not isinstance(exposure, dict) or not exposure:
                exposure = evidence.get("sidechain_exposure", {})
            evidence.update(sidechain_exposure=dict(exposure or {}), chemical_descriptors=chemical)
            evidence["site_change_summary"] = {
                "notation": getattr(option, "notation", None), "position": position,
                "wild_type": wt, "mutant": mutant,
                "basis": "stability.AA_PROPERTIES; nominal descriptor deltas",
                "protonation_status": chemical.get("protonation_status"),
            }
            option.impact_evidence = evidence

        raw = candidate_impact(candidate)
        # `POLICY` identifies a ranking/selection policy; descriptive evidence
        # must not make a generic candidate appear mechanism-ranked.
        raw.pop("policy", None)
        raw.pop("ranking_version", None)
        candidate.impact_evidence.update(raw)
        candidate.impact_evidence["descriptive_basis"] = "AA_PROPERTIES_and_observed_structure_context; not a qualification gate or selection policy"
        candidate.impact_evidence["mutation_site_details"] = [dict(option.impact_evidence.get("site_change_summary", {}),
            chemical_descriptors=dict(option.impact_evidence.get("chemical_descriptors", {})),
            sidechain_exposure=dict(option.impact_evidence.get("sidechain_exposure", {}))) for option in _mutations(candidate)]
        his_ranges = candidate_histidine_display_ranges(candidate, by_position)
        if his_ranges is not None:
            candidate.impact_evidence["display_chemistry_ranges"] = his_ranges
    annotate_local_interactions(candidates, list(by_position.values()))
    return candidates


def candidate_histidine_display_ranges(candidate, records):
    """Enumerate display-only His microstates without changing ranking inputs."""
    from .coverage_profile import AA_PROPERTIES, ACCEPTORS, DONORS, descriptors

    axes = AXES
    options = _mutations(candidate)
    his_sites = []
    for option in options:
        wt = getattr(getattr(option, "residue", None), "one_letter", None)
        mutant = getattr(option, "mutant_aa", None)
        side = "wild_type" if wt == "H" else "mutant" if mutant == "H" else None
        if side:
            his_sites.append({"notation": getattr(option, "notation", None), "position": _position(option),
                "histidine_side": side,
                "states": [{"state": "neutral", "charge": 0, "donor": 1, "acceptor": 1},
                           {"state": "protonated", "charge": 1, "donor": 2, "acceptor": 0}]})
    if not his_sites:
        return None

    state_choices = list(product(("neutral", "protonated"), repeat=len(his_sites)))
    scenarios = []
    axis_values = {axis: [] for axis in axes}
    site_axis_values = {}
    signed_values = {axis: {field: [] for field in ("net_change", "absolute_sum", "opposition_fraction")}
                     for axis in ("charge_delta", "hydropathy_delta")}
    interaction_values = {}
    for choice in state_choices:
        state_by_notation = {row["notation"]: state for row, state in zip(his_sites, choice)}
        per_site = []
        delta_by_notation = {}
        for option in options:
            wt = getattr(getattr(option, "residue", None), "one_letter", None)
            mutant = getattr(option, "mutant_aa", None)
            base = descriptors(wt, mutant)
            delta = dict(base)
            notation = getattr(option, "notation", None)
            state = state_by_notation.get(notation)
            if state:
                his_side = next(row["histidine_side"] for row in his_sites if row["notation"] == notation)
                his = {"neutral": (0, 1, 1), "protonated": (1, 2, 0)}[state]
                if his_side == "wild_type":
                    other = mutant
                    delta.update(charge_delta=AA_PROPERTIES[other][2] - his[0],
                                 donor_delta=DONORS.get(other, 0) - his[1],
                                 acceptor_delta=ACCEPTORS.get(other, 0) - his[2])
                else:
                    other = wt
                    delta.update(charge_delta=his[0] - AA_PROPERTIES[other][2],
                                 donor_delta=his[1] - DONORS.get(other, 0),
                                 acceptor_delta=his[2] - ACCEPTORS.get(other, 0))
            row = {"notation": notation, "position": _position(option),
                   "chemical_descriptors": {axis: delta.get(axis) for axis in axes}}
            per_site.append(row)
            delta_by_notation[notation] = row["chemical_descriptors"]
            for axis in axes:
                site_axis_values.setdefault(notation, {}).setdefault(axis, []).append(float(delta[axis]))
        mean_abs = {axis: sum(abs(float(row["chemical_descriptors"][axis])) for row in per_site) / len(per_site)
                    for axis in axes}
        for axis in axes:
            axis_values[axis].append(mean_abs[axis])
        signed_summary = {}
        for axis in ("charge_delta", "hydropathy_delta"):
            values = [float(row["chemical_descriptors"][axis]) for row in per_site]
            total, absolute = sum(values), sum(abs(value) for value in values)
            signed_summary[axis] = {"net_change": total, "absolute_sum": absolute,
                                    "opposition_fraction": 1 - abs(total) / absolute if absolute else None}
            for field, value in signed_summary[axis].items():
                if value is not None:
                    signed_values[axis][field].append(float(value))
        pairs = []
        for left, right in combinations(options, 2):
            left_pos, right_pos = _position(left), _position(right)
            left_record, right_record = records.get(left_pos), records.get(right_pos)
            a, b = getattr(left_record, "ca", None), getattr(right_record, "ca", None)
            valid = all(point is not None and len(point) == 3
                        and all(_finite(value) is not None for value in point) for point in (a, b))
            distance = math.dist(a, b) if valid else None
            dleft, dright = delta_by_notation[left.notation], delta_by_notation[right.notation]
            charge = dleft["charge_delta"] * dright["charge_delta"] < 0
            hydropathy = dleft["hydropathy_delta"] * dright["hydropathy_delta"] < 0
            pair = {"positions": [left_pos, right_pos], "wt_ca_distance_a": distance,
                    "nearby_opposed_charge": (bool(charge and distance < 10) if distance is not None else None),
                    "nearby_opposed_hydropathy": (bool(hydropathy and distance < 10) if distance is not None else None)}
            pairs.append(pair)
            key = tuple(pair["positions"])
            interaction_values.setdefault(key, {"wt_ca_distance_a": distance,
                "nearby_opposed_charge": [], "nearby_opposed_hydropathy": []})
            interaction_values[key]["nearby_opposed_charge"].append(pair["nearby_opposed_charge"])
            interaction_values[key]["nearby_opposed_hydropathy"].append(pair["nearby_opposed_hydropathy"])
        state_label = ";".join(f"{row['notation']}={state}" for row, state in zip(his_sites, choice))
        scenarios.append({"scenario_id": state_label,
            "histidine_states": [{"notation": row["notation"], "position": row["position"],
                                  "side": row["histidine_side"], "state": state}
                                 for row, state in zip(his_sites, choice)],
            "per_site_signed_deltas": per_site, "mean_absolute_delta": mean_abs,
            "signed_change_summary": signed_summary,
            "local_interaction_diagnostics": {"basis": "WT_CA_COORDINATES", "pairs": pairs}})

    def bounds(values):
        return {"minimum": min(values), "maximum": max(values)} if values else {"minimum": None, "maximum": None}
    ranges = {"axes": {axis: bounds(axis_values[axis]) for axis in axes},
        "per_site_signed_deltas": {notation: {axis: bounds(values) for axis, values in by_axis.items()}
                                   for notation, by_axis in site_axis_values.items()},
        "signed_change_summary": {axis: {field: bounds(values) for field, values in by_field.items()}
                                  for axis, by_field in signed_values.items()},
        "local_interactions": []}
    for positions, values in sorted(interaction_values.items()):
        row = {"positions": list(positions), "wt_ca_distance_a": values["wt_ca_distance_a"]}
        for field in ("nearby_opposed_charge", "nearby_opposed_hydropathy"):
            states = values[field]
            classification = ("DISTANCE_UNAVAILABLE" if not any(isinstance(value, bool) for value in states)
                else "PARTIALLY_UNAVAILABLE" if not all(isinstance(value, bool) for value in states)
                else "CONSISTENT_TRUE" if all(states) else "CONSISTENT_FALSE" if not any(states)
                else "STATE_DEPENDENT")
            row[field] = {"classification": classification, "scenario_values": states}
        ranges["local_interactions"].append(row)
    return {"schema_version": "display_chemistry_ranges_v1",
        "basis": "HIS_NEUTRAL_CHARGE0_DONOR1_ACCEPTOR1_VS_PROTONATED_CHARGE1_DONOR2_ACCEPTOR0_NO_PH_ASSUMPTION_v1",
        "protonation_assumption": "neutral and protonated states enumerated independently for each mutated histidine; no pH assignment",
        "microstates": his_sites, "scenario_count": len(scenarios), "scenarios": scenarios, "ranges": ranges}


def audit_single_thresholds(single, sectors, records, graph, parameters):
    """Describe REVIEW-band losses without admitting them to the strict pool."""
    by_position = {record.uniprot_pos: record for record in records}
    lower = float(parameters.get("acceptable_max_ddg", .5))
    upper = float(parameters.get("review_max_ddg", 1.0))
    audits = {}
    for sector in sectors:
        members = set(sector.residue_positions)
        rows = []
        for option in single:
            position = _position(option)
            value = _finite(option.spurs_ddg_kcal_mol)
            if position not in members or value is None or not lower < value <= upper:
                continue
            record = by_position.get(position)
            exposure = (getattr(option, "impact_evidence", {}) or {}).get("sidechain_exposure", {})
            rsasa = _finite(getattr(record, "rsasa", None))
            relative_sc = _finite(exposure.get("value"))
            otherwise = not option.hard_exclusions and not option.review_required and position in graph
            rows.append({"position": position, "native_aa": option.residue.one_letter,
                         "mutant_aa": option.mutant_aa, "spurs_ddg_kcal_mol": value,
                         "rsasa": rsasa, "relative_sidechain_exposure": relative_sc,
                         "otherwise_eligible": otherwise, "hard_exclusions": list(option.hard_exclusions),
                         "review_required": option.review_required})
        eligible = [row for row in rows if row["otherwise_eligible"]]
        audits[sector.sector_id] = {
            "model": "SPURS_SINGLE", "lower_exclusive_kcal_mol": lower, "upper_inclusive_kcal_mol": upper,
            "high_exposure_cutoff": .6, "review_only": True,
            "review_band_count": len(rows), "otherwise_eligible_count": len(eligible),
            "otherwise_eligible_position_count": len({row["position"] for row in eligible}),
            "high_rsasa_count": sum(row["rsasa"] is not None and row["rsasa"] > .6 for row in eligible),
            "high_sidechain_exposure_count": sum(row["relative_sidechain_exposure"] is not None
                and row["relative_sidechain_exposure"] > .6 for row in eligible),
            "rsasa_unknown_count": sum(row["rsasa"] is None for row in eligible),
            "sidechain_exposure_unknown_count": sum(row["relative_sidechain_exposure"] is None for row in eligible),
            "rows": sorted(rows, key=lambda row: (row["position"], row["mutant_aa"])),
            "interpretation": "Counts are substitutions, not confirmed safe mutations; no exposure-dependent threshold relaxation."}
    return audits


def dominates(left, right) -> bool:
    """Pareto dominance across exposure-min and six absolute chemical axes.

    Both inputs may be candidates or impact dictionaries. Cross-order comparisons
    are allowed for double/triple preference; missing masks must match and
    unknown axes are never imputed as zero.
    """
    a = _record(left) if not isinstance(left, dict) or "mean_absolute_delta" not in left else left
    b = _record(right) if not isinstance(right, dict) or "mean_absolute_delta" not in right else right
    if a.get("missing_mask") != b.get("missing_mask"):
        return False
    av = [a.get("exposure_min"), *(a.get("mean_absolute_delta", {}).get(axis) for axis in AXES)]
    bv = [b.get("exposure_min"), *(b.get("mean_absolute_delta", {}).get(axis) for axis in AXES)]
    comparable = [(x, y) for x, y in zip(av, bv) if x is not None and y is not None]
    return bool(comparable) and all(x >= y for x, y in comparable) and any(x > y for x, y in comparable)


def _skyline(rows, exposure_field="exposure_min", omit_group=None, omit_axis=None):
    """Exact first Pareto layer for one raw-axis sensitivity variant."""
    import numpy as np

    valid = []
    for candidate, impact in rows:
        values = [impact.get(exposure_field), *(impact.get("mean_absolute_delta", {}).get(a) for a in AXES)]
        keep = ([0] if omit_group != "exposure" else [])
        keep.extend(index + 1 for index, axis in enumerate(AXES)
                    if axis != omit_axis and (omit_group is None or axis not in GROUPS[omit_group]))
        valid.append((candidate, impact, tuple(None if values[i] is None else float(values[i]) for i in keep)))
    comparable_indices = [index for index in range(len(valid[0][2]))
                          if all(row[2][index] is not None for row in valid)]
    if not comparable_indices:
        return {id(candidate): 0 for candidate, _, _ in valid}
    vectors = sorted({tuple(values[i] for i in comparable_indices) for _, _, values in valid},
                     key=lambda vector: tuple(-value for value in vector))
    front = np.empty((len(vectors), len(comparable_indices)), dtype=np.float64)
    front_size = 0
    vector_layers = {}
    for vector in vectors:
        values = np.asarray(vector, dtype=np.float64)
        # Lexicographic descending order guarantees a later vector cannot
        # dominate an earlier one. Exact duplicate vectors were collapsed.
        dominated = front_size > 0 and bool(np.any(np.all(front[:front_size] >= values, axis=1)))
        vector_layers[vector] = 1 if dominated else 0
        if not dominated:
            front[front_size] = values
            front_size += 1
    layers = {id(candidate): vector_layers[tuple(values[i] for i in comparable_indices)]
              for candidate, _, values in valid}
    return layers


def _group_scores(impact, reference, exposure_field="exposure_min", normalization="ordinal", omit=None,
                  rank_cache=None, omit_axis=None):
    scores = {}
    for name, axes in GROUPS.items():
        if name == omit:
            continue
        vals = [impact.get(exposure_field) if name == "exposure" else impact.get("mean_absolute_delta", {}).get(axis)
                for axis in axes]
        axis_scores = []
        for axis, value in zip(axes, vals):
            if value is None or axis == omit_axis:
                continue
            ref_key = exposure_field if name == "exposure" else axis
            ref_vals = (reference.get("values", {}).get(ref_key, []) if isinstance(reference, dict) else
                        [item.get(exposure_field) if name == "exposure" else item.get("mean_absolute_delta", {}).get(axis)
                         for item in reference])
            cache_key = (id(reference), ref_key, normalization)
            if rank_cache is not None and cache_key in rank_cache:
                value_scores = rank_cache[cache_key]
            else:
                usable = [float(v) for v in ref_vals if v is not None and math.isfinite(float(v))]
                if normalization == "ordinal":
                    unique = sorted(set(usable))
                    value_scores = ({v: i / (len(unique) - 1) for i, v in enumerate(unique)}
                                    if len(unique) > 1 else {})
                else:
                    low, high = (min(usable), max(usable)) if usable else (0.0, 0.0)
                    value_scores = ({v: (v - low) / (high - low) for v in usable} if high > low else {})
                if rank_cache is not None:
                    rank_cache[cache_key] = value_scores
            if not value_scores:
                continue
            if float(value) in value_scores:
                axis_scores.append(value_scores[float(value)])
        if axis_scores:
            scores[name] = max(axis_scores)
    return scores


def _geometry_tiebreak(candidate):
    evidence = getattr(candidate, "impact_evidence", {}) or {}
    support = evidence.get("regional_geometry_support") or getattr(candidate, "geometry_support", {}) or {}
    return (float(support.get("target_mean_coverage_fraction") or 0.0),
            float(support.get("mean_marginal_area_a2") or 0.0),
            -float(support.get("mean_overlap_fraction") or 0.0))


def _risk(candidate):
    stability = getattr(candidate, "stability", None)
    combo = [_finite(getattr(stability, name, None)) for name in
             ("spurs_multi_ddg_kcal_mol", "foldx_ddg_kcal_mol")] if stability is not None else []
    known_combo = [value for value in combo if value is not None]
    if known_combo:
        return max(max(0.0, value) for value in known_combo)
    values = [_finite(getattr(option, "spurs_ddg_kcal_mol", None)) for option in _mutations(candidate)]
    return sum(max(0.0, value) for value in values if value is not None) if all(v is not None for v in values) else math.inf


def _risk_representative(candidate):
    """Lowest known combination risk, or lowest summed single SPURS score."""
    stability = getattr(candidate, "stability", None)
    combo = [_finite(getattr(stability, name, None)) for name in
             ("spurs_multi_ddg_kcal_mol", "foldx_ddg_kcal_mol")] if stability is not None else []
    known_combo = [value for value in combo if value is not None]
    if known_combo:
        return max(max(0.0, value) for value in known_combo)
    singles = [_finite(getattr(option, "spurs_ddg_kcal_mol", None)) for option in _mutations(candidate)]
    return sum(singles) if singles and all(value is not None for value in singles) else math.inf


def _reference_profile(items, order, missing_mask):
    values = {
        field: sorted({float(item[field]) for item in items
                       if item.get(field) is not None and math.isfinite(float(item[field]))})
        for field in ("exposure_min", "exposure_mean")
    }
    values.update({axis: sorted({float(item.get("mean_absolute_delta", {}).get(axis))
                                for item in items
                                if item.get("mean_absolute_delta", {}).get(axis) is not None
                                and math.isfinite(float(item["mean_absolute_delta"][axis]))})
                   for axis in AXES})
    return {"size": len(items), "order": order, "missing_mask": list(missing_mask), "values": values}


def rank_impacts(candidates, reference=None, assess_sensitivity=True):
    """Rank candidates using frozen-reference mechanism-group leximin scores."""
    candidates = list(candidates)
    if not candidates:
        return []
    frozen = None
    if reference is not None:
        frozen = [_record(item) if not isinstance(item, dict) else item for item in reference]
    if frozen is None:
        for candidate in candidates:
            evidence = getattr(candidate, "impact_evidence", {}) or {}
            saved = evidence.get("regional_reference_profiles")
            if isinstance(saved, dict) and saved:
                frozen = saved
                break
    if frozen is None:
        frozen = [candidate_impact(candidate) for candidate in candidates]
    elif not isinstance(frozen, dict):
        frozen = [dict(item) for item in frozen]
    by_order_mask = {}
    if isinstance(frozen, dict):
        for key, profile in frozen.items():
            order_text, _, mask_text = key.partition("|")
            mask = tuple(filter(None, mask_text.split(",")))
            by_order_mask[(int(order_text), mask)] = profile
        frozen_size = sum(int(profile.get("size", 0)) for profile in frozen.values())
    else:
        for item in frozen:
            by_order_mask.setdefault((item.get("order"), tuple(item.get("missing_mask", []))), []).append(item)
        frozen_size = len(frozen)
    impacts = {id(candidate): candidate_impact(candidate) for candidate in candidates}
    reference_profiles = {}
    if isinstance(frozen, dict):
        reference_profiles = dict(frozen)
    else:
        grouped_references = {}
        for item in frozen:
            grouped_references.setdefault((item.get("order"), tuple(item.get("missing_mask", []))), []).append(item)
        for (order, mask), members in grouped_references.items():
            reference_profiles[f"{order}|{','.join(mask)}"] = _reference_profile(members, order, mask)
    groups = {}
    rank_cache = {}
    variant_scores = {}
    for candidate in candidates:
        impact = impacts[id(candidate)]
        mask = tuple(impact.get("missing_mask", []))
        groups.setdefault((impact.get("order"), mask), []).append(candidate)
    layers = {}
    for key, members in groups.items():
        layers.update(_skyline([(candidate, impacts[id(candidate)]) for candidate in members]))
    for candidate in candidates:
        impact = impacts[id(candidate)]
        key = (impact.get("order"), tuple(impact.get("missing_mask", [])))
        profile_key = f"{impact.get('order')}|{','.join(impact.get('missing_mask', []))}"
        refs = reference_profiles.get(profile_key) or by_order_mask.get(key) or [impact]
        if key not in by_order_mask:
            by_order_mask[key] = refs
        scores = _group_scores(impact, refs, rank_cache=rank_cache)
        utilities = sorted(scores.values())
        if assess_sensitivity:
            minmax = _group_scores(impact, refs, normalization="minmax", rank_cache=rank_cache)
            mean_exposure = _group_scores(impact, refs, exposure_field="exposure_mean", rank_cache=rank_cache)
            without_aromaticity = _group_scores(impact, refs, omit_axis="aromatic_delta", rank_cache=rank_cache)
            leave_one = {name: _group_scores(impact, refs, omit=name, rank_cache=rank_cache) for name in GROUP_NAMES}
            sensitivity = {
                "status": "EVALUATED",
                "scope": "within_order_and_missing_mask",
                "evaluated_candidate_count": len(groups[key]),
                "ordinal_vs_minmax_leximin_agreement": utilities == sorted(minmax.values()),
                "min_vs_mean_exposure_leximin_agreement": utilities == sorted(mean_exposure.values()),
                "leave_one_group_out_leximin": {name: sorted(value.values()) for name, value in leave_one.items()},
                "robust": None,
            }
            variant_scores[id(candidate)] = {"ordinal": scores, "minmax": minmax,
                                             "mean_exposure": mean_exposure,
                                             "without_aromaticity": without_aromaticity,
                                             **{f"leave_one_out:{name}": score for name, score in leave_one.items()}}
        else:
            sensitivity = {"status": "NOT_EVALUATED_GENERATION_STAGE",
                           "scope": "within_order_and_missing_mask",
                           "evaluated_candidate_count": 0, "sensitive": None,
                           "robust": None, "variants": []}
        profile = refs if isinstance(refs, dict) else _reference_profile(refs, impact.get("order"), impact.get("missing_mask", []))
        reference_profiles.setdefault(profile_key, profile)
        row = {**impact, "pareto_layer": layers.get(id(candidate), 1),
               "group_aggregation": "maximum_axis_utility", "group_definitions": GROUPS,
               "mechanism_group_utilities": scores, "leximin_utility": utilities,
               "sensitivity": sensitivity, "reference_size": profile.get("size", len(refs)),
               "reference_values": profile.get("values", {})}
        evidence = dict(getattr(candidate, "impact_evidence", {}) or {})
        evidence.update(policy=POLICY, ranking_version=RANKING_VERSION, regional_impact=impact,
                        axes=impact["axes"], exposure_min=impact["exposure_min"],
                        exposure_mean=impact["exposure_mean"], missing_axes=impact["missing_axes"],
                        signed_change_summary=impact["signed_change_summary"],
                        aromatic_mutant_count=impact["aromatic_mutant_count"],
                        introduced_aromatic_count=impact["introduced_aromatic_count"],
                        sensitivity=sensitivity,
                        ranking=row, selection_rank=None)
        evidence["regional_reference_profiles"] = reference_profiles
        candidate.impact_evidence = evidence
    for key, members in groups.items() if assess_sensitivity else ():
        variants = ["ordinal", "minmax", "mean_exposure", "without_aromaticity", *(f"leave_one_out:{name}" for name in GROUP_NAMES)]
        variant_layers = {"ordinal": {id(item): layers[id(item)] for item in members},
                          "minmax": {id(item): layers[id(item)] for item in members}}
        variant_layers["mean_exposure"] = _skyline(
            [(item, impacts[id(item)]) for item in members], exposure_field="exposure_mean")
        variant_layers["without_aromaticity"] = _skyline(
            [(item, impacts[id(item)]) for item in members], omit_axis="aromatic_delta")
        for name in GROUP_NAMES:
            variant_layers[f"leave_one_out:{name}"] = _skyline(
                [(item, impacts[id(item)]) for item in members], omit_group=name)
        winner_by_variant = {}
        for variant in variants:
            ranked_variant = sorted(members, key=lambda candidate: (
                -int(candidate.impact_evidence["regional_impact"]["complete"]),
                variant_layers[variant][id(candidate)],
                tuple(-value for value in sorted(variant_scores[id(candidate)][variant].values())),
                tuple(-value for value in _geometry_tiebreak(candidate)),
                _risk(candidate), str(candidate.candidate_id)))
            winner_by_variant[variant] = ranked_variant[0].candidate_id if ranked_variant else None
        baseline_winner = winner_by_variant.get("ordinal")
        for candidate in members:
            sensitivity = candidate.impact_evidence["ranking"]["sensitivity"]
            sensitivity["selection_winners"] = winner_by_variant
            comparable = len(members) > 1
            sensitivity["status"] = "EVALUATED" if comparable else "INSUFFICIENT_COMPARATORS"
            sensitivity["sensitive"] = (any(value != baseline_winner for name, value in winner_by_variant.items()
                                             if name != "ordinal") if comparable else None)
            sensitivity["robust"] = not sensitivity["sensitive"] if comparable else None
            sensitivity["variants"] = [{"name": name, "winner_candidate_id": winner}
                                        for name, winner in winner_by_variant.items()]
            candidate.impact_evidence["sensitivity"] = sensitivity
    def sort_key(candidate):
        evidence = candidate.impact_evidence
        row = evidence["ranking"]
        return (-int(row["complete"]), len(tuple(row.get("missing_mask", []))),
                tuple(row.get("missing_mask", [])),
                row["pareto_layer"],
                tuple(-value for value in row["leximin_utility"]),
                tuple(-value for value in _geometry_tiebreak(candidate)),
                _risk(candidate), str(candidate.candidate_id))
    ordered = sorted(candidates, key=sort_key)
    for rank, candidate in enumerate(ordered):
        candidate.impact_evidence["selection_rank"] = rank
    return ordered


def _identity(mutations):
    signature = tuple(sorted((_position(option), option.mutant_aa) for option in mutations))
    return "V41_" + hashlib.sha256(repr(signature).encode()).hexdigest()[:16]


def _mechanism_groups(mutations):
    result = set()
    for option in mutations:
        values = _descriptors(option)
        result.update(axis for axis in AXES if values.get(axis) not in (None, 0))
    return result


def _option(single):
    residue = single.residue
    position = residue.uniprot_pos
    exposure = (getattr(single, "impact_evidence", {}) or {}).get("sidechain_exposure", {})
    exposure_value = _finite(exposure.get("value")) if isinstance(exposure, dict) else None
    accessibility = min(1.0, max(0.0, exposure_value)) if exposure_value is not None else 0.0
    return V4MutationOption(
        residue=residue, mutant_aa=single.mutant_aa,
        notation=f"{residue.one_letter}{position}{single.mutant_aa}",
        chemistry_class=single.chemistry_class, spurs_ddg_kcal_mol=single.spurs_ddg_kcal_mol,
        score=0.0, eligible=True, hard_exclusions=list(single.hard_exclusions),
        warnings=list(single.soft_warnings), accessibility_score=accessibility,
        chemical_disruption_score=0.0, stability_band=single.stability_band,
        structural_risk_flags=list(single.soft_warnings), review_required=False,
        impact_evidence={**dict(getattr(single, "impact_evidence", {}) or {}),
                         "chemical_descriptors": dict(getattr(single, "chemical_descriptors", {}) or {})},
    )


def _candidate(sector, mutations, order, parent_ids=()):
    risk_values = [_finite(getattr(option, "spurs_ddg_kcal_mol", None)) for option in mutations]
    candidate = V4Candidate(
        candidate_id=_identity(mutations), sector_id=sector.sector_id,
        mutations=list(mutations), order=order, connected=True,
        structural_risk=(sum(max(0.0, value) for value in risk_values)
                         if all(value is not None for value in risk_values) else math.inf),
        anchor_region_id=sector.sector_id, search_role="REGIONAL_PERTURBATION",
        qualification="FORMAL_POOL", parent_candidate_ids=list(parent_ids),
        recommendation_class="STABILITY_FALLBACK", epitope_score=0.0,
        impact_evidence={"policy": POLICY, "generation_origin": "regional_full_single_evidence",
                         "regional_provenance": [{"sector_id": sector.sector_id,
                                                  "parent_candidate_ids": list(parent_ids)}]},
    )
    candidate.impact_evidence["regional_impact"] = candidate_impact(candidate)
    return candidate


def _valid_sequence(sequence, options):
    from .coverage_profile import sequence_reasons
    if not sequence:
        return ["sequence-unavailable"]
    mutations = [(_position(option), option.residue.one_letter, option.mutant_aa) for option in options]
    return sequence_reasons(sequence, mutations)


def _span_ok(positions, graph, allowed, maximum_span):
    from .candidate_search import _distances
    if len(positions) < 2:
        return True, 0.0
    distances = [_distances(graph, position, allowed) for position in positions]
    spans = [distances[i].get(right, math.inf)
             for i, left in enumerate(positions) for right in positions[i + 1:]]
    span = max(spans, default=0.0)
    return math.isfinite(span) and span <= maximum_span, span


def _select(candidates, cap):
    ordered = list(candidates)
    selected, seen_sites, seen_mechanisms = [], set(), set()
    risk_representative = None
    if ordered and cap > 0:
        selected.append(ordered[0])
        selected[0].impact_evidence["regional_selection_role"] = "STRONG_IMPACT_SEED"
        seen_sites.add(tuple(sorted(_position(option) for option in ordered[0].mutations)))
        seen_mechanisms.update(_mechanism_groups(ordered[0].mutations))
        risk_representative = min(ordered, key=lambda item: (_risk_representative(item), str(item.candidate_id)))
        if len(selected) < cap and risk_representative is not ordered[0]:
            selected.append(risk_representative)
            risk_representative.impact_evidence["regional_selection_role"] = "LOWEST_RISK_REPRESENTATIVE"
            seen_sites.add(tuple(sorted(_position(option) for option in risk_representative.mutations)))
            seen_mechanisms.update(_mechanism_groups(risk_representative.mutations))
    # First retain one construct per physical site set, preferring the ranked
    # candidate that adds the most distinct mechanism groups to the sample.
    selected_ids = {id(item) for item in selected}
    remaining = [item for item in ordered if id(item) not in selected_ids]
    while remaining and len(selected) < cap:
        novel = [item for item in remaining
                 if tuple(sorted(_position(option) for option in item.mutations)) not in seen_sites]
        if not novel:
            break
        winner = min(novel, key=lambda item: (
            -len(_mechanism_groups(item.mutations) - seen_mechanisms),
            item.impact_evidence.get("selection_rank", math.inf),
            str(item.candidate_id)))
        selected.append(winner)
        seen_sites.add(tuple(sorted(_position(option) for option in winner.mutations)))
        seen_mechanisms.update(_mechanism_groups(winner.mutations))
        remaining.remove(winner)
    # Spare capacity retains different replacement chemistry within a site set.
    while remaining and len(selected) < cap:
        winner = min(remaining, key=lambda item: (
            -len(_mechanism_groups(item.mutations) - seen_mechanisms),
            item.impact_evidence.get("selection_rank", math.inf), str(item.candidate_id)))
        selected.append(winner)
        seen_mechanisms.update(_mechanism_groups(winner.mutations))
        remaining.remove(winner)
    return selected, len(seen_mechanisms), risk_representative.candidate_id if risk_representative else None


def generate_regional_candidates(single, sectors, graph, parameters, canonical_sequence, surface_areas=None):
    """Enumerate each region from all eligible, non-review single evidence.

    Pair generation uses all eligible replacements and all in-region site pairs.
    Triples extend a diverse retained pair set. Sequence and 24 A geodesic gates
    apply before retaining a construct.
    """
    from .candidate_search import _distances
    from .spatial_support import annotate_spatial_support
    from .coverage_profile import sequence_reasons
    from .models import V4Sector
    options_by_position = {}
    for item in single:
        position = getattr(item.residue, "uniprot_pos", None)
        if (not getattr(item, "eligible", False) or getattr(item, "review_required", False)
                or item.hard_exclusions or _finite(item.spurs_ddg_kcal_mol) is None
                or position is None or position not in graph):
            continue
        options_by_position.setdefault(position, []).append(_option(item))
    max_span = float(parameters.get("maximum_geodesic_span_a", 24.0))
    audit = {"policy": POLICY, "eligible_sites": len(options_by_position), "regions": []}
    result = []
    for sector in sorted(sectors, key=lambda item: item.sector_id):
        allowed = set(sector.residue_positions) & set(graph)
        positions = sorted(allowed & set(options_by_position))
        distances = {position: _distances(graph, position, allowed) for position in positions}
        local_options = {position: options_by_position[position] for position in positions}
        raw_pairs, rejected_span, rejected_sequence = [], 0, 0
        for left, right in combinations(positions, 2):
            span = distances[left].get(right, math.inf)
            if not math.isfinite(span) or span > max_span:
                rejected_span += len(local_options[left]) * len(local_options[right])
                continue
            for selected in product(local_options[left], local_options[right]):
                if _valid_sequence(canonical_sequence, selected):
                    rejected_sequence += 1
                    continue
                raw_pairs.append(_candidate(sector, selected, 2))
        if surface_areas:
            annotate_spatial_support(raw_pairs, graph, surface_areas, target_positions=allowed, scope_id=sector.sector_id)
            for candidate in raw_pairs:
                candidate.impact_evidence["regional_geometry_support"] = candidate.geometry_support
        # Freeze normalization on the full regional raw pair population before caps.
        pair_reference = [candidate_impact(item) for item in raw_pairs]
        pair_site_sets = {tuple(sorted(_position(option) for option in candidate.mutations))
                          for candidate in raw_pairs}
        ranked_pairs = rank_impacts(raw_pairs, reference=pair_reference, assess_sensitivity=False)
        retained_pairs, pair_mechanism_qualified, pair_risk_representative = _select(ranked_pairs, DOUBLE_CAP)
        # Parent choice prioritizes different site sets, then distinct mechanisms.
        parents = []
        seen_site_sets, seen_mechanisms = set(), set()
        for candidate in retained_pairs:
            mechanism = frozenset(_mechanism_groups(candidate.mutations))
            sites = tuple(sorted(_position(option) for option in candidate.mutations))
            if sites not in seen_site_sets or mechanism - seen_mechanisms:
                parents.append(candidate)
                seen_site_sets.add(sites)
                seen_mechanisms.update(mechanism)
            if len(parents) >= DOUBLE_PARENT_CAP:
                break
        raw_triples, triple_span_rejected, triple_sequence_rejected = [], 0, 0
        for parent in parents:
            parent_positions = tuple(sorted(_position(option) for option in parent.mutations))
            for added_position in positions:
                if added_position in parent_positions:
                    continue
                triple_positions = tuple(sorted((*parent_positions, added_position)))
                ok, _span = _span_ok(triple_positions, graph, allowed, max_span)
                if not ok:
                    triple_span_rejected += len(local_options[added_position])
                    continue
                for added in local_options[added_position]:
                    selected = [*parent.mutations, added]
                    if _valid_sequence(canonical_sequence, selected):
                        triple_sequence_rejected += 1
                        continue
                    raw_triples.append(_candidate(sector, selected, 3, [parent.candidate_id]))
        if surface_areas:
            annotate_spatial_support(raw_triples, graph, surface_areas, target_positions=allowed, scope_id=sector.sector_id)
            for candidate in raw_triples:
                candidate.impact_evidence["regional_geometry_support"] = candidate.geometry_support
        raw_triple_count_before_dedup = len(raw_triples)
        # A triple may extend multiple retained pairs. Collapse its physical
        # identity before freezing reference ranges and preserve every parent.
        triple_by_id = {}
        for candidate in raw_triples:
            previous = triple_by_id.get(candidate.candidate_id)
            if previous is None:
                triple_by_id[candidate.candidate_id] = candidate
            else:
                previous.parent_candidate_ids = sorted(set(previous.parent_candidate_ids + candidate.parent_candidate_ids))
                previous.impact_evidence["regional_provenance"][0]["parent_candidate_ids"] = previous.parent_candidate_ids
        raw_triples = list(triple_by_id.values())
        triple_reference = [candidate_impact(item) for item in raw_triples]
        triple_site_sets = {tuple(sorted(_position(option) for option in candidate.mutations))
                            for candidate in raw_triples}
        ranked_triples = rank_impacts(raw_triples, reference=triple_reference, assess_sensitivity=False)
        retained_triples, triple_mechanism_qualified, triple_risk_representative = _select(ranked_triples, TRIPLE_CAP)
        result.extend(retained_pairs)
        result.extend(retained_triples)
        retained_pair_site_sets = {tuple(sorted(_position(option) for option in candidate.mutations))
                                   for candidate in retained_pairs}
        retained_triple_site_sets = {tuple(sorted(_position(option) for option in candidate.mutations))
                                     for candidate in retained_triples}
        pair_quota = min(4, len(pair_site_sets))
        triple_quota = min(2, len(triple_site_sets))
        audit["regions"].append({
            "sector_id": sector.sector_id, "eligible_site_count": len(positions),
            "eligible_replacement_count": sum(len(local_options[p]) for p in positions),
            "raw_double_count": len(raw_pairs), "retained_double_count": len(retained_pairs),
            "double_available_site_set_count": len(pair_site_sets),
            "double_retained_site_set_count": len(retained_pair_site_sets),
            "double_site_set_quota_target": pair_quota,
            "double_requested_site_sets": 4,
            "double_site_set_deficit": max(0, 4 - len(retained_pair_site_sets)),
            "double_site_set_quota_met": len(retained_pair_site_sets) >= pair_quota,
            "double_mechanism_qualified_count": pair_mechanism_qualified,
            "double_lowest_risk_representative_id": pair_risk_representative,
            "double_span_rejected_count": rejected_span, "double_sequence_rejected_count": rejected_sequence,
            "double_reference_size": len(pair_reference), "triple_parent_count": len(parents),
            "raw_triple_count": raw_triple_count_before_dedup, "unique_triple_count": len(raw_triples),
            "retained_triple_count": len(retained_triples),
            "triple_available_site_set_count": len(triple_site_sets),
            "triple_retained_site_set_count": len(retained_triple_site_sets),
            "triple_site_set_quota_target": triple_quota,
            "triple_requested_site_sets": 2,
            "triple_site_set_deficit": max(0, 2 - len(retained_triple_site_sets)),
            "triple_site_set_quota_met": len(retained_triple_site_sets) >= triple_quota,
            "triple_mechanism_qualified_count": triple_mechanism_qualified,
            "triple_lowest_risk_representative_id": triple_risk_representative,
            "triple_span_rejected_count": triple_span_rejected,
            "triple_sequence_rejected_count": triple_sequence_rejected,
            "triple_reference_size": len(triple_reference),
            "surface_area_available": bool(surface_areas),
            "sequence_available": bool(canonical_sequence),
        })
    return result, audit


def merge_regional_candidates(existing, new):
    """Deduplicate by physical construct and preserve all regional provenance."""
    from copy import deepcopy
    combined = {}
    for candidate in list(existing) + list(new):
        signature = tuple(sorted((_position(option), option.mutant_aa) for option in candidate.mutations))
        current = combined.get(signature)
        if current is None:
            combined[signature] = candidate
            continue
        left = current.impact_evidence if isinstance(current.impact_evidence, dict) else {}
        right = candidate.impact_evidence if isinstance(candidate.impact_evidence, dict) else {}
        provenance = left.setdefault("regional_provenance", [])
        provenance.extend(row for row in right.get("regional_provenance", []) if row not in provenance)
        current.parent_candidate_ids = sorted(set(current.parent_candidate_ids + candidate.parent_candidate_ids))
        current.targeted_hypotheses = sorted(set(current.targeted_hypotheses + candidate.targeted_hypotheses))
        if right.get("policy") == POLICY:
            old_provenance = list(left.get("regional_provenance", []))
            left.update({key: value for key, value in right.items() if key != "regional_provenance"})
            left["regional_provenance"] = old_provenance + [row for row in right.get("regional_provenance", [])
                                                              if row not in old_provenance]
            incoming_by_key = {(_position(option), option.mutant_aa): option for option in candidate.mutations}
            for option in current.mutations:
                incoming = incoming_by_key.get((_position(option), option.mutant_aa))
                if incoming is not None and getattr(incoming, "impact_evidence", None):
                    option.impact_evidence = dict(incoming.impact_evidence)
        current.impact_evidence = left
    return sorted(combined.values(), key=lambda item: (item.order, item.candidate_id))
