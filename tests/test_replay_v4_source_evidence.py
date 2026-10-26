import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from ecd_epitope_scanner.foldx_source_evidence import (
    copy_foldx_source_evidence,
    verify_portable_source_evidence,
)
from ecd_epitope_scanner import exporting
from ecd_epitope_scanner.util import sha256_file
from ecd_epitope_scanner.surface_perturbation.electrostatics import write_atoms


def make_source(tmp_path, runs=5, layout="legacy", missing=None):
    source = tmp_path / "source_run"
    if layout == "legacy":
        batch = source / "engine_runs" / "foldx_batch_abc123"
    else:
        batch = source / "fresh_models_123" / "engine_runs" / "private_batch"
    batch.mkdir(parents=True)
    (batch / "individual_list.txt").write_text("AA1V;\n", encoding="ascii")
    (batch / "rotabase.txt").write_text("rotabase\n", encoding="ascii")
    names = []
    for slot in range(runs):
        suffix = f"_1_{slot}" if runs > 1 else "_1"
        names.extend([f"input{suffix}.pdb", f"WT_input{suffix}.pdb"])
    listed = []
    dif_rows = ["Pdb\ttotal energy"]
    for slot in range(runs):
        suffix = f"_1_{slot}" if runs > 1 else "_1"
        mut, wt = f"input{suffix}.pdb", f"WT_input{suffix}.pdb"
        listed.extend([mut, wt])
        dif_rows.append(f"{mut}\t{0.25 + slot / 10}")
        for name, residue in ((mut, "VAL"), (wt, "ALA")):
            if name == missing:
                continue
            write_atoms(batch / name, {("A", 1, "", "CA"): (np.array([0., 0., 0.]), residue, "C")})
    (batch / "PdbList_input.fxout").write_text("\n".join(listed) + "\n", encoding="utf-8")
    (batch / "Dif_input.fxout").write_text("\n".join(dif_rows) + "\n", encoding="utf-8")
    write_atoms(batch / "input.pdb", {("A", 1, "", "CA"): (np.array([0., 0., 0.]), "ALA", "C")})
    (batch / "Raw_input.fxout").write_text("raw scores\n", encoding="utf-8")
    if layout == "fresh":
        (batch / "buildmodel_request.json").write_text(json.dumps({
            "command": ["foldx", "--command=BuildModel", f"--numberOfRuns={runs}"],
            "number_of_runs": runs, "time_seed_rotabase": False, "worker_sha256": "worker-sha",
            "execution_layout": "legacy_batch_v1",
            "input_sha256": sha256_file(batch / "input.pdb"),
            "mutant_file_sha256": sha256_file(batch / "individual_list.txt"),
            "rotabase_sha256": sha256_file(batch / "rotabase.txt"),
            "output_sha256": {
                **{name: sha256_file(batch / name) for name in listed},
                "PdbList_input.fxout": sha256_file(batch / "PdbList_input.fxout"),
                "Dif_input.fxout": sha256_file(batch / "Dif_input.fxout"),
            },
        }), encoding="utf-8")
        invocation = source / "foldx_invocations" / "inv-001" / "source" / "integrations" / "workers"
        invocation.mkdir(parents=True)
        (invocation / "foldx_worker.py").write_text("# frozen worker\n", encoding="utf-8")
        invocation_root = source / "foldx_invocations" / "inv-001"
        (invocation_root / "request.json").write_text("{}", encoding="utf-8")
        (invocation_root / "execution_manifest.json").write_text("{}", encoding="utf-8")
    (source / "structure_for_energy.pdb").write_text("source structure\n", encoding="utf-8")
    manifest = {
        "worker_runs": {"foldx": {"status": "success", "provenance": {
            "worker": "foldx_worker", "number_of_runs": runs, "weights_sha256": "binary-hash",
            "execution_root": str(source)}}},
        "tool_commands": {"foldx": ["python", "integrations/workers/foldx_worker.py"]},
        "parameters": {"joint_v4_ml": {"foldx_confirmation": {"number_of_runs": runs}}},
    }
    (source / "v4_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return source, manifest


@pytest.mark.parametrize("runs", [1, 3, 5])
def test_portable_evidence_validates_explicit_pdb_list_slots_after_source_moves(tmp_path, runs):
    source, manifest = make_source(tmp_path, runs=runs)
    output = tmp_path / "replay"
    output.mkdir()

    evidence = copy_foldx_source_evidence(source, output, manifest)
    shutil.move(source, tmp_path / "archived-source")

    assert evidence["batches"][0]["pair_count"] == runs
    assert evidence["actual_foldx_cli"]["status"] == "UNKNOWN_NOT_RECORDED"
    assert evidence["engine_execution_status"] == "UNKNOWN_NOT_ARCHIVED"
    assert verify_portable_source_evidence(output, evidence)
    assert (output / evidence["batches"][0]["package_path"] / "PdbList_input.fxout").is_file()


def test_fresh_models_bundle_records_buildmodel_request_and_worker_snapshot(tmp_path):
    source, manifest = make_source(tmp_path, runs=3, layout="fresh")
    output = tmp_path / "replay"
    output.mkdir()

    evidence = copy_foldx_source_evidence(source, output, manifest)

    assert evidence["actual_foldx_cli"]["status"] == "RECORDED"
    request = evidence["actual_foldx_cli"]["buildmodel_requests"][0]
    assert request["time_seed_rotabase"] is False
    assert not Path(request["package_path"]).is_absolute()
    assert evidence["source_worker_script"]["status"] == "RECORDED"
    assert evidence["worker_execution_manifest_status"] == "RECORDED"
    assert evidence["historical_sampling"]["time_seed_rotabase"] == [False]
    assert verify_portable_source_evidence(output, evidence)


def test_portable_pair_paths_and_slot_identity_survive_source_move(tmp_path):
    from ecd_epitope_scanner.surface_perturbation.paired import scan_pairs, validate_pair

    source, manifest = make_source(tmp_path, runs=1, layout="fresh")
    scan_root = source / "fresh_models_123"
    rows, failures = scan_pairs(scan_root)
    assert not failures and len(rows) == 1
    row = rows[0]
    row.update(sampling_protocol="time_rotabase_v3", sampling_role="production",
               sampling_cohort="historical_cohort", sampling_slot=7,
               source_execution_identity={"batch_path": row["batch_path"],
                                          "group": row["group"], "run": row["run"]})
    (source / "paired_model_manifest.json").write_text(json.dumps({
        "schema_version": "1.0", "sampling_slot_contract": {"cohort": "historical_cohort"},
        "records": [{"mutations": ["A1V"], "ddg": row["ddg_kcal_mol"], "paired_models": [row]}],
    }), encoding="utf-8")
    output = tmp_path / "portable"
    output.mkdir()

    evidence = copy_foldx_source_evidence(source, output, manifest)
    shutil.move(source, tmp_path / "archived-source")
    portable = json.loads((output / evidence["portable_pair_manifest_path"]).read_text(encoding="utf-8"))
    pair = portable["records"][0]["paired_models"][0]

    assert pair["sampling_role"] == "production"
    assert pair["sampling_cohort"] == "historical_cohort"
    assert pair["sampling_slot"] == 7
    assert pair["source_execution_identity"]["batch_path"] == row["batch_path"]
    for key in ("wt_path", "mut_path", "pdb_list_path", "ddg_file_path", "batch_path",
                "input_path", "mutant_file_path", "buildmodel_request_path"):
        assert not Path(pair[key]).is_absolute()
        pair[key] = str(output / pair[key])
    validate_pair(pair)
    assert verify_portable_source_evidence(output, evidence)


def test_final_export_registers_and_archives_existing_foldx_evidence(tmp_path, monkeypatch):
    source, manifest = make_source(tmp_path, runs=1)
    (source / "run_manifest.json").write_text(json.dumps({"files": {"analysis.json": ""}}), encoding="utf-8")
    (source / "analysis.json").write_text(json.dumps({"files": {}}), encoding="utf-8")
    monkeypatch.setattr(exporting, "build_offline_html", lambda *_: "<html></html>")
    monkeypatch.setattr(exporting, "build_offline_site", lambda *_: [])

    exporting.finalize_scientific_bundle(source)

    final_analysis = json.loads((source / "analysis.json").read_text(encoding="utf-8"))
    run_manifest = json.loads((source / "run_manifest.json").read_text(encoding="utf-8"))
    evidence_files = [name for name in final_analysis["files"] if name.startswith("source_evidence/")]
    assert evidence_files
    assert all(run_manifest["files"][name] == sha256_file(source / name) for name in evidence_files)
    import zipfile
    with zipfile.ZipFile(source / "ecd_epitope_design.zip") as archive:
        assert set(evidence_files) <= set(archive.namelist())


@pytest.mark.parametrize("missing", ["input_1_2.pdb", "WT_input_1_2.pdb"])
def test_missing_mutant_or_wild_type_structure_is_rejected(tmp_path, missing):
    source, manifest = make_source(tmp_path, runs=3, missing=missing)
    output = tmp_path / "replay"
    output.mkdir()

    with pytest.raises(ValueError, match="Invalid FoldX PdbList pair evidence"):
        copy_foldx_source_evidence(source, output, manifest)


def test_missing_listed_slot_and_corrupted_evidence_are_rejected(tmp_path):
    source, manifest = make_source(tmp_path, runs=3)
    batch = next(source.rglob("PdbList_*.fxout"))
    batch.write_text("\n".join(batch.read_text().splitlines()[:-2]) + "\n", encoding="utf-8")
    output = tmp_path / "replay"
    output.mkdir()

    with pytest.raises(ValueError, match="Incomplete FoldX PdbList pairs"):
        copy_foldx_source_evidence(source, output, manifest)

    source, manifest = make_source(tmp_path / "second", runs=1)
    output = tmp_path / "second" / "replay"
    output.mkdir()
    evidence = copy_foldx_source_evidence(source, output, manifest)
    artifact = output / evidence["batches"][0]["package_path"] / "Dif_input.fxout"
    artifact.write_text("corrupt", encoding="utf-8")
    assert not verify_portable_source_evidence(output, evidence)


@pytest.mark.parametrize("bad_path", ["../outside.txt", "..\\outside.txt", "C:/outside.txt"])
def test_verifier_rejects_paths_outside_package_and_bad_json(tmp_path, bad_path):
    output = tmp_path / "package"
    output.mkdir()
    manifest = output / "evidence.json"
    manifest.write_text(json.dumps({"raw_artifacts": [{"path": bad_path, "sha256": "x"}]}), encoding="utf-8")
    evidence = {"manifest_path": "evidence.json", "manifest_sha256": __import__(
        "hashlib").sha256(manifest.read_bytes()).hexdigest()}
    assert not verify_portable_source_evidence(output, evidence)
    manifest.write_text("{bad", encoding="utf-8")
    assert not verify_portable_source_evidence(output, evidence)
