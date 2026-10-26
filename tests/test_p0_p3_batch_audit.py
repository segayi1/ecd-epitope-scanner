import hashlib
import json
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_p0_p3_batch import _stage_diversity, audit_task
from verify_p0_p3_binder_task import _source_identity_assessment


def _write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")


def _pair(task, signature, slot, value, duplicate_content=False):
    folder = task / "pairs" / signature.replace(";", "_")
    folder.mkdir(parents=True, exist_ok=True)
    wt = folder / f"wt_{slot}.pdb"
    mut = folder / f"mut_{slot}.pdb"
    wt.write_bytes(b"same-wt")
    mut.write_bytes(b"same-mut" if duplicate_content or slot == 0 else f"mut-{slot}".encode())
    return {"sampling_protocol": "time_rotabase_v3", "sampling_role": "production",
        "sampling_cohort": "ordinary_v4_sequential_v1", "sampling_slot": slot,
        "batch_path": str(folder), "group": 1, "run": slot,
        "wt_path": str(wt), "wt_sha256": hashlib.sha256(wt.read_bytes()).hexdigest(),
        "mut_path": str(mut), "mut_sha256": hashlib.sha256(mut.read_bytes()).hexdigest(),
        "ddg_kcal_mol": value}


def _candidate(cid, position, spurs, values, formal, reasons=None, order=2):
    slots = list(range(len(values)))
    repeat_audit = {"foldx_repeat_slots": slots, "foldx_repeat_values_kcal_mol": values,
        "foldx_repeat_count": len(values), "foldx_evidence_complete": slots == list(range(5))}
    return {"candidate_id": cid, "order": order, "formal_eligible": formal,
        "formal_exclusion_reasons": reasons or [], "warnings": [], "mutations": [
            {"notation": f"A{position}V", "structural_risk_flags": []}],
        "stability": {"spurs_multi_ddg_kcal_mol": spurs,
                      "foldx_ddg_kcal_mol": sorted(values)[len(values)//2] if len(values) == 5 else None},
        "impact_evidence": {"stability_repeat_audit": repeat_audit}}


def _fixture(tmp_path):
    task = tmp_path / "TST"
    task.mkdir()
    candidates = [
        _candidate("c1", 1, .8, [1.0, 1.5, 2.1, 2.5, 3.0], True),
        _candidate("c2", 2, 2.4, [.1, .2, .3], False, ["FOLDX_REPEAT_EVIDENCE_INCOMPLETE"]),
        _candidate("c3", 3, None, [], False, ["SPURS_MISSING"]),
    ]
    _write_json(task / "panel_candidates.json", candidates)
    _write_json(task / "selected_panel_v4.json", {"selected_candidate_ids": ["c1"]})
    _write_json(task / "sector_recommendations_v4.json", [{"sector_id": "R1", "candidate_id": "c1",
        "recommendation_role": "PRIMARY", "formal_eligible": True}])
    _write_json(task / "analysis.json", {"request": {"identifier": "TST"}})
    (task / "residue_mapping.csv").write_text(
        "normalized_position,uniprot_pos,chain_id,auth_seq_id,insertion_code,label_seq_id,one_letter\n"
        "1,1,A,1,,,A\n2,2,A,2,,,A\n3,3,A,3,,,A\n", encoding="utf-8")
    initial_records = []
    for candidate, signature in zip(candidates[:2], ("A1V", "A2V")):
        phase1_values = candidate["impact_evidence"]["stability_repeat_audit"]["foldx_repeat_values_kcal_mol"][:3]
        pairs = [_pair(task, signature, slot, value, duplicate_content=(signature == "A1V" and slot == 1))
                 for slot, value in enumerate(phase1_values)]
        initial_records.append({"mutations": [signature], "paired_models": pairs})
    (task / "engine_runs" / "phase1").mkdir(parents=True)
    _write_json(task / "engine_runs" / "phase1" / "paired_model_manifest.json", {
        "sampling_slot_contract": {"protocol": "time_rotabase_v3", "role": "production",
            "cohort": "ordinary_v4_sequential_v1", "requested_slots": [0, 1, 2]},
        "records": [{"mutations": ["A1V"], "paired_models": initial_records[0]["paired_models"]},
                    {"mutations": ["A2V"], "paired_models": initial_records[1]["paired_models"]}]})
    _write_json(task / "engine_runs" / "phase2.json", {"records": []})
    phase2 = task / "engine_runs" / "phase2"
    phase2.mkdir()
    _write_json(phase2 / "paired_model_manifest.json", {
        "sampling_slot_contract": {"protocol": "time_rotabase_v3", "role": "production",
            "cohort": "ordinary_v4_sequential_v1", "requested_slots": [3, 4]},
        "records": [{"mutations": ["A1V"], "paired_models": [
            _pair(task, "A1V", 3, 2.5), _pair(task, "A1V", 4, 3.0)]}]})
    source_archive = task / "executed_source_snapshot.zip"
    with zipfile.ZipFile(source_archive, "w") as archive:
        archive.writestr("pipeline.py", b"source")
    source_digest = hashlib.sha256(source_archive.read_bytes()).hexdigest()
    source_metadata = {
        "schema_version": "run_source_snapshot_v1", "captured_at_utc": "fixture",
        "source_archive": {"path": source_archive.name, "sha256": source_digest},
        "source_files": [{"path": "pipeline.py", "sha256": hashlib.sha256(b"source").hexdigest()}]}
    source_metadata["source_files_sha256"] = hashlib.sha256(json.dumps(source_metadata["source_files"],
        sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()
    source_metadata["snapshot_id"] = hashlib.sha256(json.dumps(source_metadata, sort_keys=True,
        separators=(",", ":"), ensure_ascii=False, default=str).encode("utf-8")).hexdigest()
    _write_json(task / "run_source_snapshot_v1.json", source_metadata)
    _write_json(task / "run_source_identity_v1.json", {"executed_source_snapshot": {
        "archive_sha256": source_digest}})
    sequential = {"policy": "ordinary_v4_sequential_v1", "unique_combinations": 2,
        "initial_slots": 6, "supplemental_slots": 2, "total_slots": 8,
        "attempted_refinement_candidate_ids": ["c1"]}
    manifest = {"parameters": {"foldx_confirmation": {"max_candidates": 80}},
        "provenance": {"global_candidate_pool": {"foldx_repeat_budget": sequential,
            "foldx_allocation": {"requested": {"local": 60, "local_strong_exploration": 16,
                "global": 20, "total": 80},
                "selected": {"local": 61, "local_strong_exploration": 16, "global": 19, "total": 80},
                "shortage": {"local": 0, "local_strong_exploration": 0, "global": 1, "total": 0}}}}}
    _write_json(task / "v4_manifest_v41.json", manifest)
    names = ["panel_candidates.json", "selected_panel_v4.json", "sector_recommendations_v4.json",
        "analysis.json", "residue_mapping.csv", "run_source_snapshot_v1.json", "run_source_identity_v1.json",
        "executed_source_snapshot.zip", "v4_manifest_v41.json"]
    hashes = {name: hashlib.sha256((task / name).read_bytes()).hexdigest() for name in names}
    _write_json(task / "run_manifest.json", {"files": hashes})
    with zipfile.ZipFile(task / "ecd_epitope_design.zip", "w") as archive:
        for name in names:
            archive.write(task / name, name)
        archive.write(task / "run_manifest.json", "run_manifest.json")
    return task


def test_batch_audit_checks_slots_gates_hashes_and_separate_denominators(tmp_path):
    result = audit_task(_fixture(tmp_path))
    assert result["checks"]["all_formal_primary_panel_candidates_have_five_slots"]
    assert result["checks"]["all_formal_regional_primaries_have_five_slots"]
    assert result["checks"]["paired_evidence_integrity_pass"]
    assert result["checks"]["foldx_requested_allocation_and_reallocation_within_total_and_exploration_caps"]
    assert result["counts"]["foldx_unscheduled"] == 1
    assert result["counts"]["foldx_incomplete_repeat_evidence"] == 1
    assert result["counts"]["foldx_complete_five"] == 1
    assert result["counts"]["duplicate_pair_content_count"] >= 1
    complete = next(row for row in result["candidates"] if row["candidate_id"] == "c1")
    assert complete["stability_state_at_3_slots"] == "PASS_PROXY_GATES"
    assert complete["stability_state_at_5_slots"] == "REVIEW_PROXY_BAND"
    assert complete["stability_state_changed_3_to_5"]
    assert complete["foldx_5_slot_median_kcal_mol"] == 2.1
    assert abs(complete["foldx_3_slot_minus_spurs_signed_gap_kcal_mol"] - 0.7) < 1e-12


def test_batch_audit_rejects_four_point_candidates_and_changed_exports(tmp_path):
    task = _fixture(tmp_path)
    candidates = json.loads((task / "panel_candidates.json").read_text(encoding="utf-8"))
    candidates.append(_candidate("c4", 4, .2, [], False, order=4))
    _write_json(task / "panel_candidates.json", candidates)
    result = audit_task(task)
    assert not result["checks"]["candidate_orders_only_double_triple"]
    assert not result["checks"]["export_manifest_and_zip_hashes_pass"]


def test_local_global_reallocation_is_allowed_and_diversity_axes_are_separate():
    stages = {"candidate_pool": [
        {"candidate_id": "a", "mutations": ["A1V"], "site_set": [1],
         "mutation_descriptors": [{"hydropathy_delta": -1.0, "charge_delta": 0.0, "aromatic_delta": 0.0}]},
        {"candidate_id": "b", "mutations": ["A1L"], "site_set": [1],
         "mutation_descriptors": [{"hydropathy_delta": -2.0, "charge_delta": 1.0, "aromatic_delta": 1.0}]},
        {"candidate_id": "c", "mutations": ["A2V"], "site_set": [2],
         "mutation_descriptors": [{"hydropathy_delta": None, "charge_delta": 0.0, "aromatic_delta": 0.0}]},
    ]}
    # Two different mutation identities can share one site set; this is a separate diversity axis.
    summary = _stage_diversity(stages["candidate_pool"], stages, {1}, {1, 2})["candidate_pool"]
    assert summary["unique_chemical_identities"] == 3
    assert summary["unique_site_sets"] == 2
    assert summary["native_interface_mutation_fraction"] == 2 / 3
    assert summary["accessible_surface_background_interface_fraction"] == 0.5
    assert summary["interface_enrichment_ratio_vs_accessible_surface"] == (2 / 3) / 0.5
    assert summary["hydropathy_shift"]["scored_candidate_count"] == 2
    assert summary["charge_shift"]["scored_candidate_count"] == 3


def test_shared_process_source_identity_is_explicit_and_narrowly_inferred():
    deployed = {"pipeline.py": "p", "exporting.py": "old"}
    exact = _source_identity_assessment(deployed, deployed)
    assert exact["ok"] and exact["mode"] == "EXACT_TASK_SNAPSHOT"
    inferred = _source_identity_assessment(
        {"pipeline.py": "p", "exporting.py": "captured-from-disk"}, deployed,
        {"same_recovery_process": True, "anchor_snapshot_matches_deployment": True,
         "anchor_captured_at_process_start": True, "import_chain_loaded_before_task": True,
         "no_module_reload_path": True})
    assert inferred["ok"] and inferred["mode"] == "SAME_PROCESS_IMPORT_INFERRED"
    assert inferred["mismatches"] == {"exporting.py": {
        "task_snapshot": "captured-from-disk", "deployment": "old"}}
    insufficient = _source_identity_assessment(
        {"pipeline.py": "other", "exporting.py": "captured-from-disk"}, deployed,
        {"same_recovery_process": True, "anchor_snapshot_matches_deployment": True,
         "anchor_captured_at_process_start": True, "import_chain_loaded_before_task": True,
         "no_module_reload_path": True})
    assert not insufficient["ok"] and insufficient["mode"] == "UNRESOLVED"
