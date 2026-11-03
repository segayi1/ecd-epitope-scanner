from pathlib import Path
import json
import zipfile

import pytest

from ecd_epitope_scanner.joint_design_v4.impact_completeness import (
    card_evidence_completeness,
    write_card_acceptance,
)
from ecd_epitope_scanner.job_timing import timing_key
from ecd_epitope_scanner.exporting import write_scientific_outputs
from ecd_epitope_scanner.models import ProteinContext


HASHES = {
    "analysis.json": "a" * 64,
    "panel_candidates.json": "b" * 64,
    "sector_recommendations_v4.json": "c" * 64,
}


def _card(slots):
    return {
        "candidate_id": "candidate",
        "impact_evidence": {"stability_repeat_audit": {
            "foldx_repeat_slots": slots,
            "foldx_repeat_values_kcal_mol": [1.0] * len(slots),
        }},
    }


@pytest.mark.parametrize("count", [3, 5, 10])
def test_card_repeat_completeness_accepts_protocol_slots(count):
    slots = list(range(count))
    receipt = card_evidence_completeness([_card(slots)], HASHES, slots)
    missing = receipt["cards"][0]["missing_fields"]
    assert "impact_evidence.stability_repeat_audit.foldx_repeat_slots" not in missing
    assert "impact_evidence.stability_repeat_audit.foldx_repeat_values_kcal_mol" not in missing


def test_card_repeat_completeness_rejects_wrong_declared_slots():
    receipt = card_evidence_completeness([_card([0, 1, 2, 3, 4])], HASHES, list(range(10)))
    assert "impact_evidence.stability_repeat_audit.foldx_repeat_slots" in receipt["cards"][0]["missing_fields"]
    assert "impact_evidence.stability_repeat_audit.foldx_repeat_values_kcal_mol" in receipt["cards"][0]["missing_fields"]


def test_card_repeat_completeness_accepts_each_adaptive_card_length():
    receipt = card_evidence_completeness(
        [_card([0, 1, 2]), _card(list(range(10)))],
        HASHES,
        [list(range(3)), list(range(10))],
    )
    assert all(
        "impact_evidence.stability_repeat_audit.foldx_repeat_slots"
        not in card["missing_fields"]
        for card in receipt["cards"]
    )


def test_acceptance_uses_ordinary_adaptive_protocol_and_legacy_defaults(tmp_path):
    for name in ("panel_candidates.json", "sector_recommendations_v4.json"):
        (tmp_path / name).write_text("{}", encoding="utf-8")
    analysis = {
        "execution_protocol": {"name": "ordinary_adaptive_v1", "foldx_repeat_count": 3},
        "design_summary": {"sector_recommendations": [_card([0, 1, 2])]},
    }
    receipt = write_card_acceptance(tmp_path, analysis)
    assert receipt["status"] == "INCOMPLETE"  # Other card evidence is intentionally absent.
    assert analysis["design_summary"]["display_evidence_acceptance"]["policy"] == "ordinary_adaptive_v1"
    card_missing = receipt["cards"][0]["missing_fields"]
    assert "impact_evidence.stability_repeat_audit.foldx_repeat_slots" not in card_missing

    legacy = {"design_summary": {"sector_recommendations": [_card([0, 1, 2])]}}
    from ecd_epitope_scanner.joint_design_v4.impact_completeness import _repeat_protocol
    assert _repeat_protocol(legacy) == ("five_slot_recommendation_cards_v1", list(range(5)))


def test_timing_key_separates_adaptive_from_legacy():
    legacy = {"design_engine": "joint_v4_ml"}
    adaptive = {**legacy, "execution_protocol": "ordinary_adaptive_v1"}
    assert timing_key(legacy) != timing_key(adaptive)


def test_scientific_zip_includes_ledger_but_excludes_sqlite_candidate_index(tmp_path):
    output = tmp_path / "run"
    output.mkdir()
    (output / "structure_for_energy.pdb").write_text("END\\n", encoding="utf-8")
    (output / "structure.cif").write_text("data_test\\n", encoding="utf-8")
    ledger = {"candidate_count": 1, "candidates": [{"id": "ledger-only"}]}
    (output / "ordinary_candidate_ledger.json").write_text(
        json.dumps(ledger), encoding="utf-8"
    )
    (output / "ordinary_candidate_index.sqlite3").write_bytes(b"sqlite sentinel")
    protein = ProteinContext(
        accession="P00000", protein_name="test", organism="human", sequence="A",
        ecd_start=1, ecd_end=1, membrane_anchor_position=1, topology="single-pass",
    )

    write_scientific_outputs(output, {}, protein, [], [], [])

    manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    with zipfile.ZipFile(output / "ecd_epitope_design.zip") as archive:
        names = set(archive.namelist())
        assert "ordinary_candidate_ledger.json" in names
        assert "ordinary_candidate_index.sqlite3" not in names
        assert json.loads(archive.read("run_manifest.json"))["files"] == manifest["files"]
        assert "ordinary_candidate_ledger.json" in manifest["files"]
        assert "ordinary_candidate_index.sqlite3" not in manifest["files"]
