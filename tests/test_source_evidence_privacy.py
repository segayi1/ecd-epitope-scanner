import gzip
import json
import shutil
import zipfile
from pathlib import Path

from ecd_epitope_scanner import exporting
from ecd_epitope_scanner.foldx_source_evidence import (
    copy_foldx_source_evidence,
    verify_portable_source_evidence,
)
from ecd_epitope_scanner.integrations.paired import validate_pair
from ecd_epitope_scanner.util import sanitize_text_file, sha256_file


def _write_ca_pdb(path: Path, residue: str) -> None:
    path.write_text(
        f"ATOM      1  CA  {residue} A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n",
        encoding="ascii",
    )


def _source_with_private_receipts(tmp_path: Path):
    source = tmp_path / "source"
    batch = source / "engine_runs" / "foldx_batch"
    batch.mkdir(parents=True)
    (batch / "individual_list.txt").write_text("AA1V;\n", encoding="ascii")
    _write_ca_pdb(batch / "input_1.pdb", "VAL")
    _write_ca_pdb(batch / "WT_input_1.pdb", "ALA")
    _write_ca_pdb(batch / "input.pdb", "ALA")
    (batch / "PdbList_input.fxout").write_text("input_1.pdb\nWT_input_1.pdb\n", encoding="utf-8")
    (batch / "Dif_input.fxout").write_text("Pdb\ttotal energy\ninput_1.pdb\t0.25\n", encoding="utf-8")
    private_binary = ('C:' + '\\Users\\private-owner\\foldx\\foldx.exe')
    (batch / "buildmodel_request.json").write_text(json.dumps({
        "command": [private_binary, "--command=BuildModel"],
        "number_of_runs": 1,
        "time_seed_rotabase": False,
        "input_sha256": sha256_file(batch / "input.pdb"),
        "mutant_file_sha256": sha256_file(batch / "individual_list.txt"),
        "output_sha256": {
            name: sha256_file(batch / name)
            for name in ("input_1.pdb", "WT_input_1.pdb", "PdbList_input.fxout", "Dif_input.fxout")
        },
    }), encoding="utf-8")
    _write_ca_pdb(source / "structure_for_energy.pdb", "ALA")
    manifest = {
        "worker_runs": {"foldx": {"status": "success", "provenance": {
            "worker": "foldx_worker",
            "execution_root": str(source),
            "repair_cache_path": str(source / "foldx_repair_cache" / "private" / "structure_Repair.pdb"),
        }}},
        "tool_commands": {"foldx": [private_binary, "foldx_worker.py"]},
    }
    (source / "v4_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return source, manifest


def test_portable_foldx_evidence_redacts_copies_and_rebinds_hashes(tmp_path):
    source, manifest = _source_with_private_receipts(tmp_path)
    source_receipt_digest = sha256_file(next(source.rglob("buildmodel_request.json")))
    output = tmp_path / "portable"
    output.mkdir()

    evidence = copy_foldx_source_evidence(source, output, manifest)
    shutil.move(source, tmp_path / "moved-source")

    assert verify_portable_source_evidence(output, evidence)
    text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in output.rglob("*")
        if path.is_file() and path.suffix in {".json", ".jsonl", ".fxout", ".txt"}
    )
    assert "private-owner" not in text
    request = json.loads(next(output.rglob("buildmodel_request.json")).read_text(encoding="utf-8"))
    assert request["command"][0] == "<local>/foldx.exe"
    assert evidence["actual_foldx_cli"]["buildmodel_requests"][0]["command"][0] == "<local>/foldx.exe"
    sanitized = next(
        row for row in evidence["raw_artifacts"]
        if row["path"].endswith("buildmodel_request.json")
    )
    assert sanitized["sha256"] != sanitized["source_sha256"]

    portable = json.loads(next(output.rglob("portable_pairs.json")).read_text(encoding="utf-8"))
    pair = portable["records"][0]["paired_models"][0]
    assert pair["buildmodel_request_sha256"] == sha256_file(
        output / pair["buildmodel_request_path"])
    assert pair["buildmodel_request_sha256"] != source_receipt_digest
    for key in ("batch_path", "wt_path", "mut_path", "pdb_list_path", "ddg_file_path",
                "buildmodel_request_path"):
        if pair.get(key):
            pair[key] = str(output / pair[key])
    validate_pair(pair)


def test_finalizer_sanitizes_all_export_text_before_hashing(tmp_path, monkeypatch):
    output = tmp_path / "run"
    output.mkdir()
    private = ('C:' + '\\Users\\private-owner\\run')
    analysis = {"request": {"output_dir": private}, "files": {}, "parameters": {"binary": private}}
    (output / "analysis.json").write_text(json.dumps(analysis), encoding="utf-8")
    (output / "execution_trace.jsonl").write_text(
        json.dumps({"message": private}) + "\n", encoding="utf-8"
    )
    (output / "notes.csv").write_text(f"path\n{private}\n", encoding="utf-8")
    (output / "run_manifest.json").write_text(json.dumps({
        "files": {
            "analysis.json": "",
            "execution_trace.jsonl": "",
            "notes.csv": "",
        },
        "request": analysis["request"],
        "source": {},
        "parameters": analysis["parameters"],
    }), encoding="utf-8")
    monkeypatch.setattr(exporting, "build_offline_html", lambda *_: "<html></html>")
    monkeypatch.setattr(exporting, "build_offline_site", lambda *_: [])

    exporting.finalize_scientific_bundle(output)

    with zipfile.ZipFile(output / "ecd_epitope_design.zip") as archive:
        combined = b"\n".join(archive.read(name) for name in archive.namelist())
    assert b"private-owner" not in combined
    manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    for name, digest in manifest["files"].items():
        assert digest == sha256_file(output / name)


def test_gzip_sanitizer_streams_and_redacts_local_paths(tmp_path):
    private = ('C:' + '\\Users\\private-owner\\results.csv')
    artifact = tmp_path / "rows.csv.gz"
    with gzip.open(artifact, "wt", encoding="utf-8", newline="") as handle:
        handle.write(f"path,value\n{private},1\n")

    assert sanitize_text_file(artifact)
    with gzip.open(artifact, "rt", encoding="utf-8") as handle:
        assert handle.read() == "path,value\n<local-path>,1\n"
