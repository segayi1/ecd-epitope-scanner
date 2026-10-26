"""Independent engineering audit for ordinary sequential V4 FoldX task outputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def _candidate_signature(candidate):
    return tuple(sorted(str(m.get("notation", "")).upper() for m in candidate.get("mutations", [])))


def _mutations_from_worker(row, normalized_to_uniprot):
    import re
    output = []
    for item in row.get("mutations", []):
        match = re.fullmatch(r"([A-Z])(-?\d+)([A-Z])", str(item).upper())
        if not match:
            raise ValueError(f"Invalid paired mutation notation: {item}")
        wt, position, mutant = match.groups()
        position = int(position)
        output.append(f"{wt}{normalized_to_uniprot.get(position, position)}{mutant}")
    return tuple(sorted(output))


def _stability_state(spurs, foldx):
    if spurs is None or foldx is None:
        return "MISSING"
    if spurs <= 1 and foldx <= 2:
        return "PASS_PROXY_GATES"
    if spurs <= 2 and foldx <= 3:
        return "REVIEW_PROXY_BAND"
    return "OUTSIDE_PROXY_BANDS"


def _spurs_window(value):
    if value is None:
        return "MISSING"
    if value <= .5:
        return "LE_0_5"
    if value <= 1:
        return "GT_0_5_LE_1"
    if value <= 2:
        return "GT_1_LE_2"
    return "GT_2"


def _stage_diversity(rows, stages, interface_positions, accessible_positions):
    """Describe construct/site diversity and native-interface composition by funnel stage."""
    result = {}
    background_fraction = (len(interface_positions & accessible_positions) / len(accessible_positions)
                           if accessible_positions else None)
    for stage, members in stages.items():
        chemical = {tuple(sorted(row["mutations"])) for row in members}
        site_sets = {tuple(sorted(int(position) for position in row["site_set"])) for row in members}
        total_sites = sum(len(row["site_set"]) for row in members)
        interface_sites = sum(sum(int(position) in interface_positions for position in row["site_set"])
                              for row in members)
        interface_fraction = interface_sites / total_sites if total_sites else None
        result[stage] = {
            "candidate_count": len(members),
            "unique_chemical_identities": len(chemical),
            "unique_site_sets": len(site_sets),
            "chemical_identity_examples": [list(item) for item in sorted(chemical)[:100]],
            "site_set_examples": [list(item) for item in sorted(site_sets)[:100]],
            "mutation_site_count": total_sites,
            "native_interface_mutation_sites": interface_sites,
            "native_interface_mutation_fraction": interface_fraction,
            "accessible_surface_background_interface_fraction": background_fraction,
            "interface_enrichment_ratio_vs_accessible_surface": (
                interface_fraction / background_fraction
                if interface_fraction is not None and background_fraction not in (None, 0) else None),
        }
        for axis, key in (("hydropathy_delta", "hydropathy"),
                          ("charge_delta", "charge"), ("aromatic_delta", "aromatic")):
            sums, missing_candidates = [], 0
            for row in members:
                values = [descriptor.get(axis) for descriptor in row.get("mutation_descriptors", [])]
                if not values or any(_number(value) is None for value in values):
                    missing_candidates += 1
                    continue
                sums.append(sum(float(value) for value in values))
            result[stage][f"{key}_shift"] = {
                "candidate_denominator": len(members), "scored_candidate_count": len(sums),
                "missing_candidate_count": missing_candidates,
                "mean_sum_per_candidate": statistics.mean(sums) if sums else None,
                "median_sum_per_candidate": statistics.median(sums) if sums else None,
                "min_sum_per_candidate": min(sums) if sums else None,
                "max_sum_per_candidate": max(sums) if sums else None,
                "interpretation": "descriptive nominal mutation descriptor; not calibrated chemistry or binding evidence",
            }
    return result


def _verify_export_hashes(root: Path):
    manifest_path = root / "run_manifest.json"
    if not manifest_path.is_file():
        return {"status": "MISSING", "checked_files": 0, "errors": ["run_manifest.json missing"]}
    manifest = _json(manifest_path)
    errors, checked = [], 0
    files = manifest.get("files", {})
    for name, expected in files.items():
        relative = Path(name)
        path = (root / relative).resolve()
        if relative.is_absolute() or ".." in relative.parts or root.resolve() not in path.parents:
            errors.append(f"unsafe manifest path: {name}")
            continue
        if not path.is_file():
            errors.append(f"manifest artifact missing: {name}")
            continue
        checked += 1
        if _sha(path) != expected:
            errors.append(f"manifest hash mismatch: {name}")
    bundle = root / "ecd_epitope_design.zip"
    if not bundle.is_file():
        errors.append("ecd_epitope_design.zip missing")
    elif not errors:
        with zipfile.ZipFile(bundle) as archive:
            names = set(archive.namelist())
            for name, expected in files.items():
                if name not in names:
                    errors.append(f"bundle member missing: {name}")
                    continue
                if hashlib.sha256(archive.read(name)).hexdigest() != expected:
                    errors.append(f"bundle member hash mismatch: {name}")
            if "run_manifest.json" not in names:
                errors.append("bundle run_manifest.json missing")
            elif hashlib.sha256(archive.read("run_manifest.json")).hexdigest() != _sha(manifest_path):
                errors.append("bundle run_manifest.json differs from disk")
    return {"status": "PASS" if not errors else "FAIL", "checked_files": checked, "errors": errors}


def _verify_source_snapshot(root: Path):
    metadata_path = root / "run_source_snapshot_v1.json"
    if not metadata_path.is_file():
        return {"status": "MISSING", "checked_files": 0, "errors": ["run_source_snapshot_v1.json missing"]}
    metadata = _json(metadata_path)
    archive_spec = metadata.get("source_archive") or {}
    archive = root / str(archive_spec.get("path", ""))
    errors, checked = [], 0
    rows_digest = hashlib.sha256(json.dumps(metadata.get("source_files", []), sort_keys=True,
        separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()
    if metadata.get("source_files_sha256") != rows_digest:
        errors.append("source file-list fingerprint mismatch")
    snapshot_without_id = dict(metadata)
    recorded_snapshot_id = snapshot_without_id.pop("snapshot_id", None)
    calculated_snapshot_id = hashlib.sha256(json.dumps(snapshot_without_id, sort_keys=True,
        separators=(",", ":"), ensure_ascii=False, default=str).encode("utf-8")).hexdigest()
    if recorded_snapshot_id != calculated_snapshot_id:
        errors.append("source snapshot identity mismatch")
    if not archive.is_file():
        errors.append("executed source archive missing")
    else:
        if _sha(archive) != archive_spec.get("sha256"):
            errors.append("executed source archive hash mismatch")
        with zipfile.ZipFile(archive) as bundle:
            actual = set(bundle.namelist())
            expected = {row["path"]: row["sha256"] for row in metadata.get("source_files", [])}
            if actual != set(expected):
                errors.append("source archive member list differs from source snapshot")
            for name, digest in expected.items():
                if name in actual:
                    checked += 1
                    if hashlib.sha256(bundle.read(name)).hexdigest() != digest:
                        errors.append(f"source member hash mismatch: {name}")
    identity_path = root / "run_source_identity_v1.json"
    if identity_path.is_file() and archive_spec:
        identity = _json(identity_path).get("executed_source_snapshot", {})
        if identity.get("archive_sha256") != archive_spec.get("sha256"):
            errors.append("run source identity does not bind the archived source hash")
    return {"status": "PASS" if not errors else "FAIL", "checked_files": checked, "errors": errors,
            "snapshot_id": metadata.get("snapshot_id"), "archive_sha256": archive_spec.get("sha256")}


def _load_pairs(root: Path, normalized_to_uniprot):
    by_signature, errors, manifest_count = defaultdict(list), [], 0
    identities, worker_contracts = [], []
    for path in sorted(root.rglob("paired_model_manifest.json")):
        if {"paired_manifest_history", "prior_attempts"} & set(path.parts):
            continue
        manifest_count += 1
        manifest = _json(path)
        worker_contracts.append((manifest.get("binary_sha256"), manifest.get("input_structure_sha256"),
                                 manifest.get("repaired_structure_sha256")))
        contract = manifest.get("sampling_slot_contract") or {}
        if contract and (contract.get("protocol") != "time_rotabase_v3"
                         or contract.get("cohort") != "ordinary_v4_sequential_v1"
                         or contract.get("role") != "production"):
            errors.append(f"unexpected FoldX sampling contract: {path.relative_to(root)}")
        for record in manifest.get("records", []):
            signature = _mutations_from_worker(record, normalized_to_uniprot)
            requested = contract.get("requested_slots")
            record_slots = [int(pair["sampling_slot"]) for pair in record.get("paired_models", [])
                            if pair.get("sampling_slot") is not None]
            if requested is not None and sorted(record_slots) != sorted(int(value) for value in requested):
                errors.append(f"paired rows disagree with manifest slot contract: {path.relative_to(root)}")
            for pair in record.get("paired_models", []):
                slot = pair.get("sampling_slot")
                if slot is None:
                    errors.append(f"paired row lacks logical slot: {path.relative_to(root)}")
                    continue
                row = dict(pair)
                row["manifest"] = path.relative_to(root).as_posix()
                row["signature"] = signature
                row["sampling_slot"] = int(slot)
                value = _number(row.get("ddg_kcal_mol"))
                if value is None:
                    errors.append(f"nonfinite paired FoldX value: {path.relative_to(root)}")
                    continue
                row["ddg_kcal_mol"] = value
                identity = (row.get("batch_path"), row.get("group"), row.get("run"))
                row["identity"] = identity
                if identity in identities:
                    errors.append(f"paired model identity repeated across manifests: {identity}")
                identities.append(identity)
                batch_path = Path(str(row.get("batch_path", ""))).resolve()
                if root.resolve() not in batch_path.parents:
                    errors.append(f"paired batch path escapes task directory: {batch_path}")
                for field in ("wt_path", "mut_path", "pdb_list_path", "ddg_file_path"):
                    file_path, expected = row.get(field), row.get(field.replace("_path", "_sha256"))
                    if file_path and expected:
                        artifact = Path(file_path)
                        if root.resolve() not in artifact.resolve().parents:
                            errors.append(f"paired artifact path escapes task directory: {file_path}")
                        elif not artifact.is_file() or _sha(artifact) != expected:
                            errors.append(f"paired artifact hash mismatch: {field} {file_path}")
                for field, filename in (("input_sha256", "input.pdb"),
                                        ("mutant_file_sha256", "individual_list.txt")):
                    expected = row.get(field)
                    artifact = batch_path / filename
                    if expected and (not artifact.is_file() or _sha(artifact) != expected):
                        errors.append(f"paired batch artifact hash mismatch: {filename} {batch_path}")
                by_signature[signature].append(row)
    if worker_contracts and len(set(worker_contracts)) != 1:
        errors.append("FoldX binary/input/repaired-input identity changed across phase manifests")
    return by_signature, errors, manifest_count, len(set(identities)), len(identities), worker_contracts


def audit_task(root: Path):
    root = root.resolve()
    required = ("panel_candidates.json", "v4_manifest_v41.json", "selected_panel_v4.json",
                "sector_recommendations_v4.json", "residue_mapping.csv")
    missing = [name for name in required if not (root / name).is_file()]
    if missing:
        raise FileNotFoundError(f"{root}: missing required outputs {missing}")
    manifest = _json(root / "v4_manifest_v41.json")
    provenance = manifest.get("provenance", {})
    global_pool = provenance.get("global_candidate_pool", {})
    repeat_budget = global_pool.get("foldx_repeat_budget", {})
    allocation = global_pool.get("foldx_allocation", {})
    candidates = _json(root / "panel_candidates.json")
    selected_panel = _json(root / "selected_panel_v4.json")
    regional = _json(root / "sector_recommendations_v4.json")
    with (root / "residue_mapping.csv").open(newline="", encoding="utf-8-sig") as handle:
        normalized_to_uniprot = {int(row["normalized_position"]): int(row["uniprot_pos"])
                                 for row in csv.DictReader(handle)}
    pair_rows, pair_errors, pair_manifest_count, pair_identity_count, pair_row_count, worker_contracts = _load_pairs(root, normalized_to_uniprot)
    duplicate_pair_contents = 0
    pair_content_counts = Counter()
    for pairs in pair_rows.values():
        for pair in pairs:
            if pair.get("wt_sha256") and pair.get("mut_sha256"):
                pair_content_counts[(pair["wt_sha256"], pair["mut_sha256"])] += 1
    duplicate_pair_contents = sum(count - 1 for count in pair_content_counts.values() if count > 1)

    detail_rows, strata = [], {}
    for candidate in candidates:
        signature = _candidate_signature(candidate)
        order = int(candidate.get("order", len(signature)))
        stability = candidate.get("stability") or {}
        spurs = _number(stability.get("spurs_multi_ddg_kcal_mol"))
        audit = (candidate.get("impact_evidence") or {}).get("stability_repeat_audit", {})
        pairs = sorted(pair_rows.get(signature, []), key=lambda row: row["sampling_slot"])
        values_by_slot = {}
        for pair in pairs:
            slot = pair["sampling_slot"]
            if slot in values_by_slot:
                pair_errors.append(f"duplicate logical slot {slot} for {';'.join(signature)}")
            values_by_slot[slot] = pair["ddg_kcal_mol"]
        slots = sorted(values_by_slot)
        if slots and slots != list(range(min(slots), max(slots) + 1)):
            pair_errors.append(f"noncontiguous slots for {';'.join(signature)}: {slots}")
        foldx_3 = (statistics.median(values_by_slot[s] for s in (0, 1, 2))
                   if all(s in values_by_slot for s in (0, 1, 2)) else None)
        foldx_5 = (statistics.median(values_by_slot[s] for s in range(5))
                   if all(s in values_by_slot for s in range(5)) else None)
        value_list = [values_by_slot[s] for s in slots]
        reported_slots = audit.get("foldx_repeat_slots")
        if reported_slots is not None and sorted(int(value) for value in reported_slots) != slots:
            pair_errors.append(f"candidate repeat slots disagree with pair manifests for {candidate.get('candidate_id')}")
        reported_values = audit.get("foldx_repeat_values_kcal_mol")
        if reported_values is not None and len(reported_values) == len(value_list):
            if any(abs(float(left) - right) > 1e-8 for left, right in zip(reported_values, value_list)):
                pair_errors.append(f"candidate repeat values disagree with pair manifests for {candidate.get('candidate_id')}")
        elif reported_values not in (None, []) and value_list:
            pair_errors.append(f"candidate repeat value count disagrees with pair manifests for {candidate.get('candidate_id')}")
        formal = bool(candidate.get("formal_eligible"))
        reasons = candidate.get("formal_exclusion_reasons", []) or []
        structural_risks = [mutation.get("structural_risk_flags", []) for mutation in candidate.get("mutations", [])]
        consensus = stability.get("consensus_status")
        high_risk = bool(any(structural_risks) or consensus in {"COMBINATION_REJECT", "MODEL_DISAGREEMENT"})
        signed_gap = foldx_5 - spurs if foldx_5 is not None and spurs is not None else None
        signed_gap_3 = foldx_3 - spurs if foldx_3 is not None and spurs is not None else None
        state_3, state_5 = _stability_state(spurs, foldx_3), _stability_state(spurs, foldx_5)
        window = _spurs_window(spurs)
        row = {
            "candidate_id": candidate.get("candidate_id"), "mutations": list(signature), "order": order,
            "site_set": sorted({int(__import__("re").fullmatch(r"[A-Z](-?\d+)[A-Z]", item).group(1))
                                 for item in signature}),
            "mutation_descriptors": [
                (mutation.get("impact_evidence", {}) or {}).get("chemical_descriptors")
                or mutation.get("chemical_descriptors", {})
                for mutation in candidate.get("mutations", [])],
            "spurs_window": window, "spurs_multi_ddg": spurs,
            "foldx_observed_slots": slots, "foldx_repeat_count": len(values_by_slot),
            "foldx_repeat_values_kcal_mol": value_list,
            "foldx_3_slot_median_kcal_mol": foldx_3, "foldx_5_slot_median_kcal_mol": foldx_5,
            "foldx_min_kcal_mol": min(value_list) if value_list else None,
            "foldx_max_kcal_mol": max(value_list) if value_list else None,
            "foldx_range_kcal_mol": max(value_list)-min(value_list) if value_list else None,
            "foldx_sd_kcal_mol": statistics.stdev(value_list) if len(value_list) > 1 else None,
            "foldx_pair_identity_count": len({row["identity"] for row in pairs}),
            "foldx_pair_content_duplicate_count": sum(
                1 for pair in pairs if pair.get("wt_sha256") and pair.get("mut_sha256")
                and pair_content_counts[(pair["wt_sha256"], pair["mut_sha256"])] > 1),
            "foldx_minus_spurs_signed_gap_kcal_mol": signed_gap,
            "absolute_spurs_foldx_gap_kcal_mol": abs(signed_gap) if signed_gap is not None else None,
            "foldx_3_slot_minus_spurs_signed_gap_kcal_mol": signed_gap_3,
            "absolute_spurs_foldx_3_slot_gap_kcal_mol": abs(signed_gap_3) if signed_gap_3 is not None else None,
            "stability_state_at_3_slots": state_3, "stability_state_at_5_slots": state_5,
            "stability_state_changed_3_to_5": state_3 != state_5 and state_3 != "MISSING" and state_5 != "MISSING",
            "proxy_gate_recommendation_change": ("USE_FIVE_SLOT_REVIEW_STATE_FOR_FORMAL_GATING"
                if state_3 != state_5 and state_5 == "REVIEW_PROXY_BAND" else
                "REMOVE_FROM_PROXY_BANDS" if state_3 != state_5 and state_5 == "OUTSIDE_PROXY_BANDS" else
                "STABILITY_PROXY_STATE_UNCHANGED" if state_3 == state_5 and state_5 != "MISSING" else
                "INSUFFICIENT_REPEAT_EVIDENCE"),
            "formal_eligible_saved": formal, "formal_exclusion_reasons": reasons,
            "high_risk_evidence": high_risk, "formally_excluded_saved": not formal,
            "evidence_source": "candidate saved repeats plus immutable paired model manifests",
        }
        detail_rows.append(row)
        repeat_band = ("UNSCHEDULED" if not slots else "INITIAL_3" if slots == [0, 1, 2]
                       else "COMPLETE_5" if slots == [0, 1, 2, 3, 4] else "INCOMPLETE")
        key = f"order_{order}/{window}/{repeat_band}"
        bucket = strata.setdefault(key, {"order": order, "spurs_window": window, "repeat_band": repeat_band,
            "candidate_denominator": 0, "spurs_scored": 0, "foldx_scheduled": 0,
            "foldx_complete_5": 0, "formal_eligible": 0, "formally_excluded": 0,
            "high_risk_evidence": 0,
            "signed_gap_values": [], "absolute_gap_values": [],
            "signed_gap_3_slot_values": [], "absolute_gap_3_slot_values": [],
            "3_to_5_state_changes": 0,
            "denominator_note": "candidate rows in this saved task/order/SPURS-window/repeat stratum; cross-category flags can overlap"})
        bucket["candidate_denominator"] += 1
        bucket["spurs_scored"] += spurs is not None
        bucket["foldx_scheduled"] += bool(slots)
        bucket["foldx_complete_5"] += len(values_by_slot) == 5 and slots == list(range(5))
        bucket["formal_eligible"] += formal
        bucket["formally_excluded"] += not formal
        bucket["high_risk_evidence"] += high_risk
        if signed_gap is not None:
            bucket["signed_gap_values"].append(signed_gap)
            bucket["absolute_gap_values"].append(abs(signed_gap))
        if signed_gap_3 is not None:
            bucket["signed_gap_3_slot_values"].append(signed_gap_3)
            bucket["absolute_gap_3_slot_values"].append(abs(signed_gap_3))
        bucket["3_to_5_state_changes"] += row["stability_state_changed_3_to_5"]

    for bucket in strata.values():
        for field in ("signed_gap_values", "absolute_gap_values",
                      "signed_gap_3_slot_values", "absolute_gap_3_slot_values"):
            values = bucket[field]
            bucket[field] = {"n": len(values), "mean": statistics.mean(values) if values else None,
                "median": statistics.median(values) if values else None,
                "min": min(values) if values else None, "max": max(values) if values else None,
                "sample_sd": statistics.stdev(values) if len(values) > 1 else None}

    by_id = {row["candidate_id"]: row for row in detail_rows}
    panel_ids = selected_panel.get("selected_candidate_ids", [])
    panel_failures = [candidate_id for candidate_id in panel_ids
                      if candidate_id not in by_id or not by_id[candidate_id]["formal_eligible_saved"]
                      or by_id[candidate_id]["foldx_observed_slots"] != list(range(5))]
    regional_primary = [row for row in regional if row.get("recommendation_role") == "PRIMARY"
                        and row.get("candidate_id")]
    regional_failures = [row.get("candidate_id") for row in regional_primary
        if row.get("candidate_id") not in by_id or not row.get("formal_eligible")
        or by_id[row["candidate_id"]]["foldx_observed_slots"] != list(range(5))]
    orders = sorted({int(candidate.get("order", 0)) for candidate in candidates})
    analysis = _json(root / "analysis.json") if (root / "analysis.json").is_file() else {}
    analysis_residues = analysis.get("residues", []) or []
    interface_positions = {int(row["uniprot_pos"]) for row in analysis_residues
                           if row.get("uniprot_pos") is not None and row.get("interface") is True}
    accessible_positions = {int(row["uniprot_pos"]) for row in analysis_residues
                            if row.get("uniprot_pos") is not None and row.get("surface_accessible") is True}
    scheduled_ids = {candidate.get("candidate_id") for candidate in candidates
                     if (candidate.get("impact_evidence", {}).get("foldx_scheduling", {}) or {}).get("selected") is True}
    scheduled_signatures = {_candidate_signature(candidate) for candidate in candidates
                            if (candidate.get("impact_evidence", {}).get("foldx_scheduling", {}) or {}).get("selected") is True}
    formal_ids = {candidate.get("candidate_id") for candidate in candidates if candidate.get("formal_eligible")}
    panel_id_set = set(panel_ids)
    stages = {
        "candidate_pool": detail_rows,
        "foldx_scheduled": [row for row in detail_rows if row["candidate_id"] in scheduled_ids],
        "formal_eligible": [row for row in detail_rows if row["candidate_id"] in formal_ids],
        "main_panel": [row for row in detail_rows if row["candidate_id"] in panel_id_set],
    }
    diversity = _stage_diversity(detail_rows, stages, interface_positions, accessible_positions)
    policy_params = manifest.get("parameters", {})
    foldx_params = policy_params.get("foldx_confirmation", {})
    source_hashes = _verify_source_snapshot(root)
    export_hashes = _verify_export_hashes(root)
    checks = {
        "ordinary_sequential_policy_present": repeat_budget.get("policy") == "ordinary_v4_sequential_v1",
        "candidate_orders_only_double_triple": bool(orders) and set(orders) <= {2, 3},
        "unique_foldx_combinations_at_most_80": int(repeat_budget.get("unique_combinations", 10**9)) <= 80,
        "paired_manifest_unique_combos_match_budget": len(pair_rows) == int(repeat_budget.get("unique_combinations", -1)),
        "scheduled_candidate_signatures_match_pair_manifests": scheduled_signatures == set(pair_rows),
        "scheduled_signature_count_matches_budget": len(scheduled_signatures) == int(repeat_budget.get("unique_combinations", -1)),
        "foldx_requested_allocation_and_reallocation_within_total_and_exploration_caps": (
            allocation.get("requested", {}).get("local") == 60
            and allocation.get("requested", {}).get("global") == 20
            and allocation.get("requested", {}).get("total") == 80
            and allocation.get("requested", {}).get("local_strong_exploration", 10**9) <= 16
            and allocation.get("selected", {}).get("local_strong_exploration", 10**9) <= 16
            and allocation.get("selected", {}).get("total", 10**9) <= 80),
        "foldx_slot_budget_at_most_300": int(repeat_budget.get("total_slots", 10**9)) <= 300,
        "initial_three_slots_per_unique_combo": int(repeat_budget.get("initial_slots", -1)) == 3*int(repeat_budget.get("unique_combinations", 0)),
        "supplemental_at_most_30_combinations": int(repeat_budget.get("supplemental_slots", 10**9)) <= 60,
        "supplemental_slots_match_attempted_combos": (
            int(repeat_budget.get("supplemental_slots", -1)) ==
            2*len(repeat_budget.get("attempted_refinement_candidate_ids", []))
            and len(repeat_budget.get("attempted_refinement_candidate_ids", [])) <= 30),
        "total_slots_equal_phase_sum": int(repeat_budget.get("total_slots", -1)) == (
            int(repeat_budget.get("initial_slots", 0)) + int(repeat_budget.get("supplemental_slots", 0))),
        "configured_unique_cap_at_most_80": int(foldx_params.get("max_candidates", 10**9)) <= 80,
        "all_formal_primary_panel_candidates_have_five_slots": not panel_failures,
        "all_formal_regional_primaries_have_five_slots": not regional_failures,
        "source_snapshot_hashes_pass": source_hashes["status"] == "PASS",
        "export_manifest_and_zip_hashes_pass": export_hashes["status"] == "PASS",
        "paired_evidence_integrity_pass": not pair_errors,
        "paired_row_count_matches_slot_budget": pair_row_count == int(repeat_budget.get("total_slots", -1)),
    }
    request = analysis.get("request", {})
    protein = request.get("identifier") or request.get("protein_name") or root.name
    from audit_panel_repeat_sensitivity import audit_panel_repeat_sensitivity
    panel_repeat_sensitivity = (
        audit_panel_repeat_sensitivity(root)
        if checks["paired_evidence_integrity_pass"] and checks["source_snapshot_hashes_pass"]
        else {"status": "BLOCKED_UPSTREAM_AUDIT_FAILURE",
              "reasons": [name for name in ("paired_evidence_integrity_pass", "source_snapshot_hashes_pass")
                          if not checks[name]]}
    )
    counts = {"candidate_denominator": len(detail_rows),
        "missing_spurs_multi": sum(row["spurs_multi_ddg"] is None for row in detail_rows),
        "foldx_unscheduled": sum(not row["foldx_observed_slots"] for row in detail_rows),
        "foldx_incomplete_repeat_evidence": sum(bool(row["foldx_observed_slots"])
            and row["foldx_observed_slots"] != list(range(5)) for row in detail_rows),
        "foldx_complete_five": sum(row["foldx_observed_slots"] == list(range(5)) for row in detail_rows),
        "formally_excluded": sum(not row["formal_eligible_saved"] for row in detail_rows),
        "high_risk_evidence": sum(row["high_risk_evidence"] for row in detail_rows),
        "formal_eligible": sum(row["formal_eligible_saved"] for row in detail_rows),
        "selected_panel_count": len(panel_ids), "selected_panel_failures": panel_failures,
        "regional_primary_count": len(regional_primary), "regional_primary_failures": regional_failures,
        "pair_manifest_count": pair_manifest_count, "unique_pair_identity_count": pair_identity_count,
        "pair_row_count": pair_row_count,
        "duplicate_pair_content_count": duplicate_pair_contents}
    return {"schema_version": "p0_p3_batch_engineering_audit_v1", "protein": protein,
        "task_directory": str(root), "interpretation": "Engineering provenance and model-output audit only; repeated prediction variability is not an experimental confidence interval or stability calibration.",
        "budgets": {"foldx_allocation": allocation, "sequential_repeat_budget": repeat_budget,
                    "configured_foldx_confirmation": foldx_params,
                    "caps": {"unique_combinations": 80, "paired_run_slots": 300,
                             "supplemental_combinations": 30, "formal_repeat_slots": [0, 1, 2, 3, 4]}},
        "counts": counts, "checks": checks,
        "source_hash_audit": source_hashes, "export_hash_audit": export_hashes,
        "paired_integrity_errors": pair_errors,
        "foldx_worker_identity_contracts": worker_contracts,
        "panel_repeat_sensitivity": panel_repeat_sensitivity,
        "candidate_strata": list(strata.values()), "candidates": detail_rows,
        "selection_diversity_and_chemistry": diversity,
        "selection_diversity_caveat": "Stage comparisons describe the saved candidate funnel. Interface enrichment and nominal chemistry shifts are not causal, calibrated, or binding evidence.",
        "diagnostic_caveats": ["SPURS/FoldX thresholds are prespecified engineering gates, not experimentally calibrated limits.",
            "Absolute/signed model gaps are descriptive and do not gate eligibility.",
            "Missing, unscheduled, incomplete and high-risk flags use separate counts and may overlap.",
            "Repeat SD/range describe this run's sampling protocol only; they are not experimental uncertainty intervals."]}


def audit_batch(paths):
    tasks = [audit_task(Path(path)) for path in paths]
    return {"schema_version": "p0_p3_multi_task_engineering_audit_v1", "task_count": len(tasks),
            "all_checks_pass": all(all(task["checks"].values()) for task in tasks), "tasks": tasks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tasks", nargs="+", help="completed V4 task output directories")
    parser.add_argument("--json", dest="json_path", required=True, help="write audit JSON here")
    parser.add_argument("--csv", dest="csv_path", help="optional candidate-level audit CSV")
    args = parser.parse_args()
    result = audit_batch(args.tasks)
    output = Path(args.json_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    if args.csv_path:
        rows = [row | {"protein": task["protein"]}
                for task in result["tasks"] for row in task["candidates"]]
        fields = list(dict.fromkeys(key for row in rows for key in row))
        path = Path(args.csv_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps({"all_checks_pass": result["all_checks_pass"], "task_count": len(result["tasks"]),
                      "json": str(output.resolve())}, ensure_ascii=False))
    return 0 if result["all_checks_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
