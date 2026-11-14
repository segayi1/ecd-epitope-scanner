"""Auditable completeness checks for recommendation-card evidence."""
from __future__ import annotations

import math
import json
from itertools import combinations


AXES = ("charge_delta", "volume_delta_a3", "hydropathy_delta",
        "donor_delta", "acceptor_delta", "aromatic_delta")


def _repeat_protocol(analysis):
    protocol = analysis.get("execution_protocol") or {}
    if isinstance(protocol, str):
        protocol = {"name": protocol}
    name = protocol.get("name", protocol.get("policy", protocol.get("id", "")))
    if name in {"ordinary_adaptive_v1", "ordinary_three_principle_adaptive_v1", "ordinary_site_first_v2", "ordinary_site_first_perturbation_v1"}:
        policy = name if name in {"ordinary_site_first_v2", "ordinary_site_first_perturbation_v1"} else "ordinary_adaptive_v1"
        slots = protocol.get("foldx_repeat_slots")
        count = protocol.get("foldx_repeat_count", len(slots) if isinstance(slots, list) else None)
        if isinstance(slots, list) and slots:
            return policy, list(slots)
        if count in (3, 10):
            return policy, list(range(count))
        return policy, [list(range(3)), list(range(10))]
    return "five_slot_recommendation_cards_v1", list(range(5))


def write_card_acceptance(output_dir, analysis):
    from ..util import sha256_file, write_json
    from ..exporting import write_analysis_artifact

    if not analysis.get("design_summary", {}).get("sector_recommendations"):
        sidecar = output_dir / "sector_recommendations_v4.json"
        if sidecar.is_file():
            analysis.setdefault("design_summary", {})["sector_recommendations"] = json.loads(
                sidecar.read_text(encoding="utf-8")
            )
    recommendations = analysis.get("design_summary", {}).get("sector_recommendations", [])
    names = ("analysis.json", "panel_candidates.json", "sector_recommendations_v4.json")
    policy, expected_slots = _repeat_protocol(analysis)
    receipt = card_evidence_completeness(recommendations, {name: "0" * 64 for name in names}, expected_slots)
    analysis.setdefault("files", {})["card_evidence_completeness.json"] = "card_evidence_completeness.json"
    analysis.setdefault("design_summary", {})["display_evidence_acceptance"] = {
        "status": receipt["status"], "card_count": receipt["card_count"],
        "incomplete_card_count": sum(bool(row["missing_fields"]) for row in receipt["cards"]),
        "policy": policy,
        "diagnostics_only": receipt["card_count"] == 0,
    }
    write_analysis_artifact(output_dir, analysis)
    receipt["inputs_sha256"] = {name: sha256_file(output_dir / name) for name in names}
    write_json(output_dir / "card_evidence_completeness.json", receipt)
    return receipt


def _finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _finite_bounds(value):
    return isinstance(value, dict) and _finite(value.get("minimum")) and _finite(value.get("maximum"))


def _bounds_match(bound, values):
    return (_finite_bounds(bound) and bool(values)
            and math.isclose(bound["minimum"], min(values), rel_tol=1e-9, abs_tol=1e-9)
            and math.isclose(bound["maximum"], max(values), rel_tol=1e-9, abs_tol=1e-9))


def _complete_display_ranges(ranges, notations):
    if (ranges.get("schema_version") != "display_chemistry_ranges_v1"
            or ranges.get("basis") != "HIS_NEUTRAL_CHARGE0_DONOR1_ACCEPTOR1_VS_PROTONATED_CHARGE1_DONOR2_ACCEPTOR0_NO_PH_ASSUMPTION_v1"
            or "no ph" not in ranges.get("protonation_assumption", "").lower()):
        return False
    his_sites = ranges.get("microstates") or []
    scenarios = ranges.get("scenarios") or []
    if not his_sites or len(scenarios) != 2 ** len(his_sites) or ranges.get("scenario_count") != len(scenarios):
        return False
    if {row.get("notation") for row in his_sites} - set(notations):
        return False
    expected_pairs = len(notations) * (len(notations) - 1) // 2
    for state in his_sites:
        if {item.get("state") for item in state.get("states", [])} != {"neutral", "protonated"}:
            return False
    for scenario in scenarios:
        if len(scenario.get("histidine_states", [])) != len(his_sites):
            return False
        sites = scenario.get("per_site_signed_deltas") or []
        if {row.get("notation") for row in sites} != set(notations):
            return False
        if any(not all(_finite((row.get("chemical_descriptors") or {}).get(axis)) for axis in AXES) for row in sites):
            return False
        means = scenario.get("mean_absolute_delta") or {}
        if not all(_finite(means.get(axis)) for axis in AXES):
            return False
        signed = scenario.get("signed_change_summary") or {}
        if any(not _finite((signed.get(axis) or {}).get(field))
               for axis in ("charge_delta", "hydropathy_delta") for field in ("net_change", "absolute_sum")):
            return False
        pairs = ((scenario.get("local_interaction_diagnostics") or {}).get("pairs") or [])
        if len(pairs) != expected_pairs or any(
                not _finite(pair.get("wt_ca_distance_a")) or not isinstance(pair.get("nearby_opposed_charge"), bool)
                or not isinstance(pair.get("nearby_opposed_hydropathy"), bool) for pair in pairs):
            return False
    bound_axes = ((ranges.get("ranges") or {}).get("axes") or {})
    if any(not _bounds_match(bound_axes.get(axis), [scenario["mean_absolute_delta"][axis] for scenario in scenarios])
           for axis in AXES):
        return False
    signed_ranges = ((ranges.get("ranges") or {}).get("signed_change_summary") or {})
    for axis in ("charge_delta", "hydropathy_delta"):
        for field in ("net_change", "absolute_sum"):
            if not _bounds_match((signed_ranges.get(axis) or {}).get(field),
                    [scenario["signed_change_summary"][axis][field] for scenario in scenarios]):
                return False
        opposition = (signed_ranges.get(axis) or {}).get("opposition_fraction")
        opposition_values = [scenario["signed_change_summary"][axis]["opposition_fraction"] for scenario in scenarios]
        finite_opposition = [value for value in opposition_values if _finite(value)]
        if finite_opposition and not _bounds_match(opposition, finite_opposition):
            return False
        if not finite_opposition and not all(
                (scenario["signed_change_summary"][axis]["absolute_sum"] == 0
                 and scenario["signed_change_summary"][axis]["opposition_fraction"] is None)
                for scenario in scenarios):
            return False
    for notation in notations:
        site_ranges = ((ranges.get("ranges") or {}).get("per_site_signed_deltas") or {}).get(notation) or {}
        if any(not _bounds_match(site_ranges.get(axis), [
                next(row["chemical_descriptors"][axis] for row in scenario["per_site_signed_deltas"]
                     if row["notation"] == notation) for scenario in scenarios]) for axis in AXES):
            return False
    return True


def card_evidence_completeness(recommendations, inputs_sha256, expected_repeat_slots=None):
    """Audit primary and alternative cards without altering their eligibility."""
    required_inputs = {"analysis.json", "panel_candidates.json", "sector_recommendations_v4.json"}
    if not required_inputs <= set(inputs_sha256) or any(
            not isinstance(inputs_sha256[name], str) or len(inputs_sha256[name]) != 64
            or any(character not in "0123456789abcdefABCDEF" for character in inputs_sha256[name])
            for name in required_inputs):
        raise ValueError("Completeness receipt requires SHA-256 hashes for all three source evidence files")
    if expected_repeat_slots is None:
        expected_repeat_slot_sets = [list(range(5))]
    elif (
            isinstance(expected_repeat_slots, (list, tuple))
            and expected_repeat_slots
            and all(isinstance(value, (list, tuple)) for value in expected_repeat_slots)):
        expected_repeat_slot_sets = [list(value) for value in expected_repeat_slots]
    else:
        expected_repeat_slot_sets = [list(expected_repeat_slots)]
    cards = []
    for recommendation in recommendations:
        displayed = [(recommendation, "PRIMARY", "primary")]
        displayed.extend((row, "ALTERNATIVE", "alternatives") for row in recommendation.get("alternatives", []))
        displayed.extend((row, "ALL_FORMAL_LOCAL", "all_formal_local_candidates")
                         for row in recommendation.get("all_formal_local_candidates", []))
        for row, display_role, display_list in displayed:
            candidate_id = row.get("candidate_id")
            if not candidate_id:
                continue
            missing = []
            evidence = row.get("impact_evidence") or {}
            ranges = evidence.get("display_chemistry_ranges") or {}
            notations = row.get("mutation_notations") or []
            display_ranges_complete = _complete_display_ranges(ranges, notations) if ranges else False
            range_axes = (ranges.get("ranges") or {}).get("axes") or {}
            axes = evidence.get("axes") or {}
            for axis in AXES:
                if not _finite(axes.get(axis)) and not (display_ranges_complete and _finite_bounds(range_axes.get(axis))):
                    missing.append(f"impact_evidence.axes.{axis}")
            for axis in ("charge_delta", "hydropathy_delta"):
                signed = (evidence.get("signed_change_summary") or {}).get(axis) or {}
                signed_ranges = ((ranges.get("ranges") or {}).get("signed_change_summary") or {}).get(axis) or {}
                for field in ("net_change", "absolute_sum"):
                    if not _finite(signed.get(field)) and not (display_ranges_complete and _finite_bounds(signed_ranges.get(field))):
                        missing.append(f"impact_evidence.signed_change_summary.{axis}.{field}")
                opposition = signed.get("opposition_fraction")
                if not _finite(opposition) and not (display_ranges_complete and _finite_bounds(signed_ranges.get("opposition_fraction"))) \
                        and not (signed.get("absolute_sum") == 0 and opposition is None):
                    missing.append(f"impact_evidence.signed_change_summary.{axis}.opposition_fraction")
            sites = row.get("site_details") or []
            expected_sites = len(row.get("mutation_notations") or [])
            if len(sites) != expected_sites or expected_sites == 0:
                missing.append("site_details.count")
            for index, site in enumerate(sites):
                prefix = f"site_details[{index}]"
                if not site.get("notation") or site.get("native_aa") is None or site.get("mutant_aa") is None:
                    missing.append(f"{prefix}.identity")
                descriptors = site.get("chemical_descriptors") or {}
                site_ranges = (((ranges.get("ranges") or {}).get("per_site_signed_deltas") or {}).get(site.get("notation")) or {})
                for axis in AXES:
                    if not _finite(descriptors.get(axis)) and not (display_ranges_complete and _finite_bounds(site_ranges.get(axis))):
                        missing.append(f"{prefix}.chemical_descriptors.{axis}")
                exposure = site.get("sidechain_exposure") or {}
                for field in ("value", "numerator_a2", "denominator_a2"):
                    if not _finite(exposure.get(field)):
                        missing.append(f"{prefix}.sidechain_exposure.{field}")
                if exposure.get("status") != "KNOWN" or not exposure.get("method"):
                    missing.append(f"{prefix}.sidechain_exposure.provenance")
            scenarios = (row.get("regional_geometry_support") or {}).get("scenarios") or {}
            for radius in (9, 12, 15):
                scenario = scenarios.get(f"geodesic_{radius}") or {}
                if not _finite(scenario.get("coverage_fraction")):
                    missing.append(f"regional_geometry_support.geodesic_{radius}.coverage_fraction")
                if not isinstance(scenario.get("uncovered_positions"), list):
                    missing.append(f"regional_geometry_support.geodesic_{radius}.uncovered_positions")
            distances = row.get("pair_distances") or []
            expected_pairs = expected_sites * (expected_sites - 1) // 2
            if len(distances) != expected_pairs:
                missing.append("pair_distances.count")
            for index, pair in enumerate(distances):
                if not _finite(pair.get("geodesic_distance_a")):
                    missing.append(f"pair_distances[{index}].geodesic_distance_a")
            local = evidence.get("local_interaction_diagnostics") or {}
            local_pairs = local.get("pairs") or []
            local_ranges = {tuple(sorted(item.get("positions", []))): item
                            for item in ((ranges.get("ranges") or {}).get("local_interactions") or [])}
            if len(local_pairs) != expected_pairs:
                missing.append("impact_evidence.local_interaction_diagnostics.pairs.count")
            for index, pair in enumerate(local_pairs):
                if not _finite(pair.get("wt_ca_distance_a")):
                    missing.append(f"impact_evidence.local_interaction_diagnostics.pairs[{index}].wt_ca_distance_a")
                for field in ("nearby_opposed_charge", "nearby_opposed_hydropathy"):
                    range_pair = local_ranges.get(tuple(sorted(pair.get("positions", []))), {})
                    classification = (range_pair.get(field) or {}).get("classification")
                    if not isinstance(pair.get(field), bool) and (not display_ranges_complete or classification not in {
                            "CONSISTENT_TRUE", "CONSISTENT_FALSE", "STATE_DEPENDENT"}):
                        missing.append(f"impact_evidence.local_interaction_diagnostics.pairs[{index}].{field}")
            repeat = evidence.get("stability_repeat_audit") or {}
            slots = repeat.get("foldx_repeat_slots")
            if slots not in expected_repeat_slot_sets:
                missing.append("impact_evidence.stability_repeat_audit.foldx_repeat_slots")
            values = repeat.get("foldx_repeat_values_kcal_mol") or []
            expected_count = len(slots) if slots in expected_repeat_slot_sets else None
            if (expected_count is None or len(values) != expected_count
                    or not all(_finite(value) for value in values)):
                missing.append("impact_evidence.stability_repeat_audit.foldx_repeat_values_kcal_mol")
            stability = row.get("stability_evidence") or {}
            for field in ("spurs_multi_ddg_kcal_mol", "foldx_ddg_kcal_mol"):
                if not _finite(stability.get(field)):
                    missing.append(f"stability_evidence.{field}")
            cards.append({"parent_sector_id": recommendation.get("sector_id"),
                          "sector_id": row.get("sector_id", recommendation.get("sector_id")),
                          "candidate_id": candidate_id, "role": display_role, "display_list": display_list,
                          "missing_fields": sorted(set(missing))})
    all_missing = [{"candidate_id": card["candidate_id"], "sector_id": card["sector_id"],
                   "missing_fields": card["missing_fields"]} for card in cards if card["missing_fields"]]
    if not cards:
        all_missing.append({"candidate_id": None, "sector_id": None,
                            "missing_fields": ["recommendation_cards.count"]})
    return {"schema_version": "card_evidence_completeness_v1",
            "status": "PASS" if not all_missing else "INCOMPLETE", "card_count": len(cards),
            "inputs_sha256": dict(inputs_sha256), "cards": cards,
            "missing_fields": all_missing}
