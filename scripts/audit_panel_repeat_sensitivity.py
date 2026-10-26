"""Read-only panel sensitivity to saved FoldX 3-vs-5 repeats; never calls models."""
from __future__ import annotations

import copy
import hashlib
import json
import statistics
import sys
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT.parent) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT.parent))


STABILITY_REASON_TOKENS = ("SPURS", "FOLDX", "CONSENSUS", "MODEL_DISAGREEMENT", "COMBINATION_REJECT")
SELECTION_CODE = ("joint_design_v4/engine.py", "joint_design_v4/optimizer.py",
                 "joint_design_v4/identifiability.py", "joint_design_v4/models.py")


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def _source_identity(root: Path):
    meta_path, archive_path = root / "run_source_snapshot_v1.json", root / "executed_source_snapshot.zip"
    if not meta_path.is_file() or not archive_path.is_file():
        return {"status": "PENDING_SOURCE_ARCHIVE"}
    metadata = _read(meta_path)
    declared = {row["path"]: row["sha256"] for row in metadata.get("source_files", [])}
    current = {}
    for name in SELECTION_CODE:
        path = root.parents[1] / name
        if not path.is_file():
            path = Path(__file__).resolve().parents[1] / name
        current[name] = _sha_bytes(path.read_bytes()) if path.is_file() else None
    with zipfile.ZipFile(archive_path) as archive:
        archived = {name: _sha_bytes(archive.read(name)) for name in SELECTION_CODE if name in archive.namelist()}
    expected = {name: declared.get(name) for name in SELECTION_CODE}
    matches_snapshot = all(expected[name] and archived.get(name) == expected[name] == current[name]
                           for name in SELECTION_CODE)
    return {"status": "PASS" if matches_snapshot else "BLOCKED_SOURCE_MISMATCH",
            "snapshot_hashes": expected, "archive_hashes": archived, "current_hashes": current,
            "matches_snapshot": matches_snapshot}


def _status(spurs, foldx, parameters):
    if spurs is None or foldx is None:
        return "MISSING"
    spurs_parameters = parameters.get("spurs_multi", parameters)
    policy = parameters.get("stability_threshold_policy") or spurs_parameters.get("stability_threshold_policy")
    if policy != "ordinary_v4_directional_v1":
        return "UNSUPPORTED_STABILITY_POLICY"
    if spurs <= float(spurs_parameters.get("acceptable_max_ddg", 1.0)) and foldx <= float(
            parameters.get("foldx_confirmation", {}).get("combination_gate_kcal_mol", 2.0)):
        return "PASS_PROXY_GATES"
    if spurs <= float(spurs_parameters.get("review_max_ddg", 2.0)) and foldx <= float(
            parameters.get("foldx_confirmation", {}).get("review_max_kcal_mol", 3.0)):
        return "REVIEW_PROXY_BAND"
    return "OUTSIDE_PROXY_BANDS"


def _is_stability_reason(reason):
    upper = str(reason).upper()
    return any(token in upper for token in STABILITY_REASON_TOKENS) or "REPEAT_EVIDENCE_INCOMPLETE" in upper


def _deserialize_candidate(row):
    from ecd_epitope_scanner.joint_design_v4.models import (
        CombinationStabilityEvidence, ResidueKey, V4Candidate, V4MutationOption,
    )
    data = copy.deepcopy(row)
    data["mutations"] = [V4MutationOption(**{**mutation, "residue": ResidueKey(**mutation["residue"])})
                         for mutation in data.get("mutations", [])]
    if data.get("stability"):
        data["stability"] = CombinationStabilityEvidence(**data["stability"])
    return V4Candidate(**data)


def _choose_panel(rows, foldx_values, source_signatures, parameters, sectors, areas, hypotheses, budget):
    from ecd_epitope_scanner.joint_design_v4.optimizer import select_panel
    selected_candidates, statuses = [], {}
    for row in rows:
        signature = tuple(sorted(item.get("notation", "").upper() for item in row.get("mutations", [])))
        if signature not in foldx_values:
            continue
        spurs = (row.get("stability") or {}).get("spurs_multi_ddg_kcal_mol")
        foldx = foldx_values[signature]
        state = _status(spurs, foldx, parameters)
        statuses[row["candidate_id"]] = state
        reasons = [reason for reason in row.get("formal_exclusion_reasons", []) if not _is_stability_reason(reason)]
        eligible = state == "PASS_PROXY_GATES" and not reasons
        if not eligible:
            continue
        candidate = _deserialize_candidate(row)
        if candidate.stability is None:
            continue
        candidate.stability.foldx_ddg_kcal_mol = float(foldx)
        candidate.stability.consensus_status = state
        candidate.formal_exclusion_reasons = []
        candidate.formal_eligible = True
        candidate.risk_components = dict(candidate.risk_components or {})
        candidate.risk_components["spurs_multi_ddg_kcal_mol"] = spurs
        candidate.risk_components["foldx_ddg_kcal_mol"] = float(foldx)
        selected_candidates.append(candidate)
    panel_params = parameters.get("panel_optimizer", {})
    panel, selected = select_panel(
        selected_candidates,
        max(1, int(parameters.get("max_primary_constructs", 7))),
        int(panel_params.get("higher_order_max_primary", 2)),
        int(panel_params.get("recommendations_per_sector", 6)),
        int(panel_params.get("epitope_priority_per_sector", 3)),
        int(panel_params.get("stability_fallback_per_sector", 3)),
        sectors=sectors, surface_areas=areas, hypotheses=hypotheses,
        risk_penalty_weights=panel_params.get("risk_penalty_weights", {}),
        panel_budget=budget, stop_parameters=panel_params)
    return list(panel.selected_candidate_ids), statuses, len(selected_candidates)


def audit_panel_repeat_sensitivity(root: Path):
    """Compare panels on frozen model scores; all calculations are local selection only."""
    root = Path(root).resolve()
    required = ("panel_candidates.json", "run_manifest.json", "selected_panel_v4.json",
                "probe_sectors_v4.json", "surface_graph_v41.json", "inference_hypotheses_v411.json",
                "run_source_snapshot_v1.json", "executed_source_snapshot.zip")
    missing = [name for name in required if not (root / name).exists()]
    if missing:
        return {"status": "PENDING_REQUIRED_EVIDENCE", "missing": missing}
    identity = _source_identity(root)
    if identity["status"] != "PASS":
        return {"status": identity["status"], "selection_code_identity": identity,
                "interpretation": "Panel replay is withheld because the active selector code does not match the task source snapshot."}
    run_manifest = _read(root / "run_manifest.json")
    parameters = run_manifest["parameters"]["joint_v4_ml"]
    spurs_parameters = parameters.get("spurs_multi", parameters)
    stability_policy = parameters.get("stability_threshold_policy") or spurs_parameters.get(
        "stability_threshold_policy")
    if stability_policy != "ordinary_v4_directional_v1":
        return {"status": "UNSUPPORTED_STABILITY_POLICY", "selection_code_identity": identity,
                "policy": stability_policy}
    from ecd_epitope_scanner.joint_design_v4.models import InferenceHypothesis, V4Sector
    from ecd_epitope_scanner.joint_design_v4.engine import _map_mutation_notation
    import csv
    with (root / "residue_mapping.csv").open(newline="", encoding="utf-8-sig") as stream:
        mapping = {int(row["normalized_position"]): int(row["uniprot_pos"]) for row in csv.DictReader(stream)}
    pair_values = {}
    for path in sorted(root.rglob("paired_model_manifest.json")):
        if "paired_manifest_history" in path.parts:
            continue
        manifest = _read(path)
        for record in manifest.get("records", []):
            mutations = record.get("mutations") or record.get("mutation_list") or []
            if isinstance(mutations, str):
                mutations = [part.strip() for part in mutations.split(",")]
            signature = tuple(sorted(_map_mutation_notation(str(item), mapping) for item in mutations))
            for pair in record.get("paired_models", []):
                slot = pair.get("sampling_slot")
                value = pair.get("ddg_kcal_mol")
                if slot is not None and value is not None:
                    pair_values.setdefault(signature, {})[int(slot)] = float(value)
    candidates = _read(root / "panel_candidates.json")
    sectors = [V4Sector(**row) for row in _read(root / "probe_sectors_v4.json")]
    graph = _read(root / "surface_graph_v41.json")
    areas = {int(node["uniprot_pos"]): float(node["surface_area_a2"]) for node in graph["nodes"]}
    hypotheses = []
    for row in _read(root / "inference_hypotheses_v411.json"):
        row = dict(row)
        row["residue_weights"] = {int(position): value for position, value in row.get("residue_weights", {}).items()}
        hypotheses.append(InferenceHypothesis(**row))
    budget = _read(root / "panel_budget_decision_v411.json")
    budget = {key: budget[key] for key in ("minimum_constructs", "target_constructs", "maximum_constructs", "adaptive_budget")}
    signatures = {tuple(sorted(mutation.get("notation", "").upper() for mutation in row.get("mutations", [])))
                  for row in candidates}
    scheduled_three, fixed_five = {}, {}
    for signature in signatures:
        values = pair_values.get(signature, {})
        if all(slot in values for slot in (0, 1, 2)):
            scheduled_three[signature] = statistics.median(values[slot] for slot in (0, 1, 2))
        if all(slot in values for slot in range(5)):
            fixed_five[signature] = {
                "three": statistics.median(values[slot] for slot in (0, 1, 2)),
                "five": statistics.median(values[slot] for slot in range(5)),
            }
    if not fixed_five:
        return {"status": "PENDING_NO_COMPLETE_FIVE_SLOT_POOL", "selection_code_identity": identity,
                "phase1_scheduled_candidate_count": len(scheduled_three), "fixed_five_slot_candidate_count": 0}
    fixed3_values = {signature: item["three"] for signature, item in fixed_five.items()}
    fixed5_values = {signature: item["five"] for signature, item in fixed_five.items()}
    ids_fixed3, statuses_fixed3, eligible_fixed3 = _choose_panel(candidates, fixed3_values, signatures,
        parameters, sectors, areas, hypotheses, budget)
    ids_fixed5, statuses_fixed5, eligible_fixed5 = _choose_panel(candidates, fixed5_values, signatures,
        parameters, sectors, areas, hypotheses, budget)
    ids_full3, statuses_full3, eligible_full3 = _choose_panel(candidates, scheduled_three, signatures,
        parameters, sectors, areas, hypotheses, budget)
    saved_panel = _read(root / "selected_panel_v4.json").get("selected_candidate_ids", [])
    formal_set = {row["candidate_id"] for row in candidates if row.get("formal_eligible")}
    return {
        "status": "COMPLETE",
        "selection_code_identity": identity,
        "method": "Read-only selection replay with saved SPURS and saved FoldX slot values; no model or structure worker calls.",
        "fixed_completed_five_slot_pool": {
            "candidate_count": len(fixed_five),
            "eligible_count_at_three_slot_median": eligible_fixed3,
            "eligible_count_at_five_slot_median": eligible_fixed5,
            "three_slot_panel_ids": ids_fixed3, "five_slot_panel_ids": ids_fixed5,
            "intersection_ids": sorted(set(ids_fixed3) & set(ids_fixed5)),
            "jaccard": len(set(ids_fixed3) & set(ids_fixed5)) / max(1, len(set(ids_fixed3) | set(ids_fixed5))),
            "panel_changed_ids": {"three_only": sorted(set(ids_fixed3)-set(ids_fixed5)),
                                  "five_only": sorted(set(ids_fixed5)-set(ids_fixed3))},
            "stability_gate_counts_at_three_slots": {state: sum(value == state for value in statuses_fixed3.values())
                                                       for state in ("PASS_PROXY_GATES", "REVIEW_PROXY_BAND", "OUTSIDE_PROXY_BANDS", "MISSING")},
            "stability_gate_counts_at_five_slots": {state: sum(value == state for value in statuses_fixed5.values())
                                                      for state in ("PASS_PROXY_GATES", "REVIEW_PROXY_BAND", "OUTSIDE_PROXY_BANDS", "MISSING")},
            "saved_final_panel_ids": saved_panel,
            "five_slot_replay_exactly_matches_saved_final_panel": ids_fixed5 == saved_panel,
            "saved_formal_candidate_ids_subset_of_five_slot_available_pool": formal_set >= set(saved_panel),
        },
        "all_phase1_scheduled_vs_final_five_slot_pool": {
            "phase1_three_slot_scheduled_candidate_count": len(scheduled_three),
            "eligible_count_at_three_slot_median": eligible_full3,
            "eligible_count_at_final_five_slot_median": eligible_fixed5,
            "three_slot_scheduled_panel_ids": ids_full3,
            "final_five_slot_panel_ids": ids_fixed5,
            "intersection_ids": sorted(set(ids_full3) & set(ids_fixed5)),
            "jaccard": len(set(ids_full3) & set(ids_fixed5)) / max(1, len(set(ids_full3) | set(ids_fixed5))),
            "panel_changed_ids": {"three_slot_only": sorted(set(ids_full3)-set(ids_fixed5)),
                                  "five_slot_only": sorted(set(ids_fixed5)-set(ids_full3))},
            "interpretation": "This contrast mixes repeat updates with candidate-evidence availability; interpret separately from the fixed five-slot pool comparison.",
        },
        "interpretation": "These are saved-score panel-selection sensitivity diagnostics, not independent model validation or experimental confidence.",
    }
