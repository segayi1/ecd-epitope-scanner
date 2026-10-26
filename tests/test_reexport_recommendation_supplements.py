import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest

from scripts.reexport_recommendation_supplements import (
    copy_refinement_manifests, digest, merge_pairs, update_target,
    validate_completed_supplement, write_card_evidence_receipt,
)


def pair(slot, value, identity):
    return {"sampling_slot": slot, "ddg_kcal_mol": value, "batch_path": f"batch-{identity}",
            "group": f"group-{identity}", "run": f"run-{identity}",
            "wt_sha256": f"wt-{identity}", "mut_sha256": f"mut-{identity}"}


def test_merge_pairs_preserves_original_evidence_and_adds_only_missing_slots():
    original = {"mutations": ["A1V"], "paired_models": [pair(i, 0.1 * i, f"old-{i}") for i in range(3)]}
    new_pairs = {"mutations": ["A1V"], "paired_models": [pair(i, 0.1 * i, f"new-{i}") for i in (3, 4)]}

    merged = merge_pairs([original], [new_pairs])[("A1V",)]

    assert merged["paired_models"][:3] == original["paired_models"]
    assert [item["sampling_slot"] for item in merged["paired_models"]] == [0, 1, 2, 3, 4]
    assert merged["number_of_runs"] == 5
    assert merged["ddg"] == 0.2


def test_update_target_clears_stale_missing_evidence_gates_but_keeps_real_gate_failure():
    candidate = {
        "candidate_id": "candidate-1",
        "mutations": [{"notation": "A1V"}],
        "stability": {"spurs_multi_ddg_kcal_mol": 0.2, "foldx_ddg_kcal_mol": None,
                      "provenance": {}},
        "formal_exclusion_reasons": ["FOLDX_REPEAT_EVIDENCE_INCOMPLETE", "COMBINATION_FOLDX_MISSING"],
        "formal_eligible": False,
    }
    merged = {("A1V",): {"paired_models": [pair(i, 4.0, f"id-{i}") for i in range(5)]}}

    update_target(candidate, merged, {"stability_threshold_policy": "ordinary_v4_directional_v1"},
                  {"combination_gate_kcal_mol": 2.0, "review_max_kcal_mol": 3.0})

    assert candidate["stability"]["consensus_status"] == "COMBINATION_REJECT"
    assert candidate["formal_exclusion_reasons"] == ["COMBINATION_REJECT"]
    assert candidate["formal_eligible"] is False
    assert candidate["impact_evidence"]["stability_repeat_audit"]["foldx_repeat_slots"] == [0, 1, 2, 3, 4]
    audit = candidate["impact_evidence"]["stability_repeat_audit"]
    assert audit["foldx_repeat_above_gate_count"] == 5
    assert audit["sampling_decision_sensitive"] is False
    assert audit["cross_execution_evidence_status"] == "NOT_ASSESSED"


def test_legacy_reexport_adds_sampling_risk_explanation_without_changing_median_qualification():
    candidate = {
        "candidate_id": "ctla4-e59y-k65f",
        "mutations": [{"notation": "E59Y"}, {"notation": "K65F"}],
        "stability": {"spurs_multi_ddg_kcal_mol": .209517, "foldx_ddg_kcal_mol": 1.21,
                      "provenance": {}},
        "formal_exclusion_reasons": [], "formal_eligible": True,
    }
    values = [2.99488, .276247, 3.4539, 1.21002, .710114]
    merged = {("E59Y", "K65F"): {"paired_models": [pair(i, value, f"legacy-{i}")
                                                   for i, value in enumerate(values)]}}

    update_target(candidate, merged, {}, {"combination_gate_kcal_mol": 2.0})

    assert candidate["stability"]["consensus_status"] == "CONSENSUS_PASS"
    assert candidate["formal_eligible"] is True
    assert candidate["impact_evidence"]["stability_repeat_audit"]["sampling_decision_sensitive"] is True
    assert candidate["recommendation_risk_label"] == "SAMPLING_SENSITIVE_REVIEW"
    assert "跨越配置门限" in candidate["recommendation_risk_explanation"]


def test_receipt_hashes_bind_to_final_output_files_not_historical_source(tmp_path):
    (tmp_path / "analysis.json").write_text(json.dumps({"files": {}}), encoding="utf-8")
    (tmp_path / "panel_candidates.json").write_text("[]", encoding="utf-8")
    (tmp_path / "sector_recommendations_v4.json").write_text("[]", encoding="utf-8")
    (tmp_path / "paired_model_manifest.json").write_text("{}", encoding="utf-8")
    (tmp_path / "run_manifest.json").write_text(json.dumps({"files": {}}), encoding="utf-8")
    historical_hash = "f" * 64

    receipt = write_card_evidence_receipt(tmp_path)

    for name in ("analysis.json", "panel_candidates.json", "sector_recommendations_v4.json"):
        assert receipt["inputs_sha256"][name] == digest(tmp_path / name)
        assert receipt["inputs_sha256"][name] != historical_hash
    assert receipt["card_count"] == 0
    assert receipt["status"] == "INCOMPLETE"
    assert receipt["missing_fields"][0]["missing_fields"] == ["recommendation_cards.count"]


def test_reexport_copies_active_refinement_manifests_into_bundle_inventory(tmp_path):
    source, output = tmp_path / "source", tmp_path / "output"
    active = source / "foldx_refinement" / "batch_01" / "paired_model_manifest.json"
    archived = source / "foldx_refinement" / "prior_attempts" / "old" / "paired_model_manifest.json"
    active.parent.mkdir(parents=True)
    archived.parent.mkdir(parents=True)
    output.mkdir()
    active.write_text('{"records": []}', encoding="utf-8")
    archived.write_text('{"records": []}', encoding="utf-8")
    (output / "analysis.json").write_text(json.dumps({"files": {}}), encoding="utf-8")
    (output / "run_manifest.json").write_text(json.dumps({"files": {}}), encoding="utf-8")

    copied = copy_refinement_manifests(source, output)

    name = "foldx_refinement/batch_01/paired_model_manifest.json"
    assert copied == [name]
    assert (output / name).read_bytes() == active.read_bytes()
    assert name in json.loads((output / "analysis.json").read_text())["files"]
    assert json.loads((output / "run_manifest.json").read_text())["files"][name] == digest(output / name)
    assert not (output / "foldx_refinement/prior_attempts/old/paired_model_manifest.json").exists()


def test_reexport_refuses_to_overwrite_different_refinement_manifest(tmp_path):
    source, output = tmp_path / "source", tmp_path / "output"
    source_manifest = source / "foldx_refinement" / "batch_01" / "paired_model_manifest.json"
    output_manifest = output / "foldx_refinement" / "batch_01" / "paired_model_manifest.json"
    source_manifest.parent.mkdir(parents=True)
    output_manifest.parent.mkdir(parents=True)
    source_manifest.write_text("source", encoding="utf-8")
    output_manifest.write_text("different", encoding="utf-8")

    with pytest.raises(ValueError, match="Existing refinement evidence differs"):
        copy_refinement_manifests(source, output)


def test_reexport_requires_finished_ledger_with_matching_group_completions(tmp_path):
    raw = tmp_path / "supplement"
    job = raw / "raw" / "job1"
    job.mkdir(parents=True)
    group_path = job / "group_slots_0_0.json"
    group_path.write_text('{"records": []}', encoding="utf-8")
    execution_id = "fixture-execution"
    plan = {"budget": {"new_unique_combinations": 1, "new_pair_slots": 1},
            "execution": {"status": "RAW_EVIDENCE_COMPLETE", "execution_id": execution_id,
                "reserved_unique_combinations": 1, "reserved_pair_slots": 1}}
    events = [
        {"event": "GROUP_RESERVED", "execution_id": execution_id, "job_id": "job1",
         "candidate_ids": ["c1"], "sampling_slots": [0], "unique_combinations": 1, "pair_slots": 1},
        {"event": "GROUP_COMPLETE", "execution_id": execution_id, "job_id": "job1",
         "candidate_ids": ["c1"], "sampling_slots": [0], "raw_records_path": str(group_path),
         "raw_records_sha256": digest(group_path)},
        {"event": "EXECUTION_COMPLETE", "execution_id": execution_id},
    ]
    (raw / "execution_ledger.jsonl").write_text("\n".join(json.dumps(row) for row in events), encoding="utf-8")

    assert validate_completed_supplement(raw, plan)["reserved_pair_slots"] == 1
    plan["execution"]["status"] = "RUNNING"
    with pytest.raises(ValueError, match="RAW_EVIDENCE_COMPLETE"):
        validate_completed_supplement(raw, plan)
