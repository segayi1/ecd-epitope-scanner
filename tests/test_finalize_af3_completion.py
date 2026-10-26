import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest

from workspace_data.af3_binder_batch_20261017.finalize_af3_completion import (
    active_pair_rows, validate_pairs, copy_refinement_evidence,
)


def test_validate_pairs_checks_pair_and_approved_paths(tmp_path):
    root = tmp_path / "approved"
    batch = root / "batch"
    batch.mkdir(parents=True)
    artifacts = {}
    for field in ("wt_path", "mut_path", "pdb_list_path"):
        path = batch / field
        path.write_text("evidence")
        artifacts[field] = str(path)
    pair = {"batch_path": str(batch), **artifacts}
    validated = []

    assert validate_pairs([pair], (root,), lambda row: validated.append(row)) == 1
    assert validated == [pair]


def test_validate_pairs_rejects_path_escape_before_pair_validation(tmp_path):
    approved = tmp_path / "approved"
    approved.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    validated = []
    pair = {"batch_path": str(outside)}

    with pytest.raises(ValueError, match="outside approved roots"):
        validate_pairs([pair], (approved,), lambda row: validated.append(row))
    assert not validated


def test_active_phase_slots_include_refinements_but_skip_archived_manifests(tmp_path):
    root = tmp_path / "job"
    refinement = root / "foldx_refinement" / "batch_01"
    history = root / "paired_manifest_history" / "old"
    prior = root / "prior_attempts" / "old"
    for path, count in ((root, 3), (refinement, 2), (history, 40), (prior, 60)):
        path.mkdir(parents=True, exist_ok=True)
        (path / "paired_model_manifest.json").write_text(json.dumps({
            "records": [{"paired_models": [{} for _ in range(count)]}]
        }), encoding="utf-8")

    phases = active_pair_rows(root)

    assert [(path.relative_to(root).as_posix(), len(rows))
            for path, _, rows in phases] == [
                ("paired_model_manifest.json", 3),
                ("foldx_refinement/batch_01/paired_model_manifest.json", 2),
            ]
    assert sum(len(rows) for _, _, rows in phases) == 5


def test_delivery_copies_refinement_phase_without_overwriting_augmented_root(tmp_path, monkeypatch):
    import ecd_epitope_scanner.scripts.reexport_recommendation_supplements as exporter
    monkeypatch.setattr(exporter, "write_card_evidence_receipt", lambda output: None)
    source, output = tmp_path / "source", tmp_path / "output"
    refinement = source / "foldx_refinement/batch_01"
    refinement.mkdir(parents=True)
    output.mkdir()
    (source / "paired_model_manifest.json").write_text('{"records": []}')
    (refinement / "paired_model_manifest.json").write_text('{"records": [{"paired_models": [{"sampling_slot": 3}]}]}')
    root_evidence = '{"records": [{"augmented": true}]}'
    (output / "paired_model_manifest.json").write_text(root_evidence)
    (output / "analysis.json").write_text('{"files": {}}')
    (output / "run_manifest.json").write_text('{"files": {}}')
    names = copy_refinement_evidence(source, output)
    assert names == ["foldx_refinement/batch_01/paired_model_manifest.json"]
    assert (output / names[0]).read_bytes() == (refinement / "paired_model_manifest.json").read_bytes()
    assert (output / "paired_model_manifest.json").read_text() == root_evidence
    assert names[0] in json.loads((output / "analysis.json").read_text())["files"]
    assert names[0] in json.loads((output / "run_manifest.json").read_text())["files"]
