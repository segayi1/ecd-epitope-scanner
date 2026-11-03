import io
import json
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import run_site_first_legacy_baseline as baseline


def make_archive(*members):
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w") as archive:
        for name, data in members:
            info = tarfile.TarInfo(name)
            payload = data.encode()
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    return stream.getvalue()


def test_git_archive_extraction_is_contained_and_rejects_traversal(tmp_path):
    destination = tmp_path / "snapshot"
    package = baseline.safe_extract_archive(
        make_archive(("ecd_epitope_scanner/pipeline.py", "# archived")), destination)
    assert package == destination / "ecd_epitope_scanner"
    assert (package / "pipeline.py").read_text() == "# archived"

    unsafe_root = tmp_path / "unsafe"
    with pytest.raises(ValueError, match="Unsafe archive member"):
        baseline.safe_extract_archive(make_archive(("../outside.txt", "no")), unsafe_root)
    assert not (tmp_path / "outside.txt").exists()


def test_request_preserves_legacy_defaults_and_uses_isolated_paths(tmp_path):
    structure = tmp_path / "structure.pdb"
    assets = tmp_path / "assets"
    binary = tmp_path / "foldx.exe"
    structure.write_text("PDB")
    binary.write_text("exe")
    source = {"metadata": {"execution_queue": {"old": True}, "joint_v4_ml": {
        "execution_queue": {"old": True}, "execution_protocol": {"name": "ordinary_site_first_v2"},
        "procaliper": {"device": "cuda"},
        "foldx_confirmation": {"number_of_runs": 5}}}}
    request = baseline.build_legacy_request(source, structure=structure, output=tmp_path / "output",
        baseline_root=tmp_path / "baseline", model_assets_root=assets, foldx_binary=binary)
    assert "execution_queue" not in request["metadata"]
    assert "execution_queue" not in request["metadata"]["joint_v4_ml"]
    joint = request["metadata"]["joint_v4_ml"]
    assert "execution_protocol" not in joint
    assert joint["procaliper"] == {"device": "cuda"}
    assert joint["foldx_confirmation"]["number_of_runs"] == 5
    assert joint["foldx_confirmation"]["binary"] == str(binary.resolve())
    assert joint["foldx_confirmation"]["short_binary_alias"] is True
    assert joint["runtime"]["cache_dir"] == str((tmp_path / "baseline" / "score_cache").resolve())
    assert request["model_assets_root"] == str(assets.resolve())


def test_preflight_resolves_commit_and_blocks_existing_baseline(tmp_path):
    source = tmp_path / "request.json"
    source.write_text(json.dumps({"input_type": "structure", "identifier": "IL7Ra_binder_AF2_35_AIB234.pdb",
        "chain": "A", "analysis_target": "binder", "design_engine": "joint_v4_ml", "ml_mode": "strict",
        "confidence_path": None}))
    assets = tmp_path / "assets"
    assets.mkdir()
    binary = tmp_path / "foldx.exe"
    names = ["foldx_binary", "spurs_single_checkpoint", "spurs_multi_checkpoint", "spurs_esm_weights",
             "discotope_model_archive"]
    identity_paths = {}
    for name in names:
        path = binary if name == "foldx_binary" else assets / name
        path.write_bytes(name.encode())
        identity_paths[name] = path
    acceptance = tmp_path / "acceptance"
    plan = baseline.preflight(root=assets, acceptance_dir=acceptance, source_request=source,
        structure=Path(baseline.SOURCE_STRUCTURE), foldx_binary=binary, model_paths_override=identity_paths,
        commit_resolver=lambda _root: "a" * 40, max_foldx_root_length=200)
    assert plan["baseline_commit"] == "a" * 40
    assert plan["short_foldx_work_root"] == str((acceptance / "baseline" / "fx").resolve())
    (acceptance / "baseline").mkdir(parents=True)
    with pytest.raises(FileExistsError, match="Refusing to reuse"):
        baseline.preflight(root=assets, acceptance_dir=acceptance, source_request=source,
            structure=Path(baseline.SOURCE_STRUCTURE), foldx_binary=binary, model_paths_override=identity_paths,
            commit_resolver=lambda _root: "a" * 40, max_foldx_root_length=200)


def test_invalid_archive_writes_needs_review_without_starting_pipeline(tmp_path):
    source = tmp_path / "request.json"
    source.write_text(json.dumps({"input_type": "structure", "identifier": "IL7Ra_binder_AF2_35_AIB234.pdb",
        "chain": "A", "analysis_target": "binder", "design_engine": "joint_v4_ml", "ml_mode": "strict",
        "confidence_path": None}))
    assets = tmp_path / "assets"
    assets.mkdir()
    identity_paths = {}
    for name in ("foldx_binary", "spurs_single_checkpoint", "spurs_multi_checkpoint", "spurs_esm_weights",
                 "discotope_model_archive"):
        path = tmp_path / name
        path.write_text(name)
        identity_paths[name] = path
    acceptance = tmp_path / "acceptance"
    report = baseline.execute(root=assets, acceptance_dir=acceptance, source_request=source,
        structure=Path(baseline.SOURCE_STRUCTURE), foldx_binary=identity_paths["foldx_binary"],
        model_paths_override=identity_paths, commit_resolver=lambda _root: "a" * 40,
        archive_provider=lambda _commit: make_archive(("../escape.py", "bad")), max_foldx_root_length=200)
    assert report["status"] == "NEEDS_REVIEW"
    assert "Unsafe archive member" in report["errors"][0]
    assert json.loads((acceptance / "baseline_report.json").read_text())["status"] == "NEEDS_REVIEW"
    assert not (tmp_path / "escape.py").exists()


def test_model_command_validation_requires_all_identity_files(tmp_path):
    binary = tmp_path / "foldx.exe"
    binary.write_text("exe")
    with pytest.raises(FileNotFoundError, match="Required legacy model identity files"):
        baseline.validate_model_paths(tmp_path, binary)


def test_actual_legacy_child_code_starts_and_emits_progress(tmp_path):
    import subprocess
    import sys
    package = tmp_path / 'ecd_epitope_scanner'
    package.mkdir()
    (package / '__init__.py').write_text('')
    (package / 'pipeline.py').write_text(
        'class PipelineRequest:\n    def __init__(self, **kwargs): pass\n'
        'def run_pipeline(request, callback):\n    callback("START", 0.1, "ready")\n    return {}\n')
    request = tmp_path / 'request.json'
    request.write_text('{}')
    completed = subprocess.run([sys.executable, '-c', baseline.CHILD_CODE, str(tmp_path), str(request)],
                               capture_output=True, text=True, check=True)
    assert [json.loads(line)['event'] for line in completed.stdout.splitlines()] == ['pipeline_start', 'progress', 'result']


def test_execute_tracks_pipeline_time_and_cold_cache_then_validates_deliverables(tmp_path):
    source = tmp_path / "request.json"
    source.write_text(json.dumps({"input_type": "structure", "identifier": "IL7Ra_binder_AF2_35_AIB234.pdb",
        "chain": "A", "analysis_target": "binder", "design_engine": "joint_v4_ml", "ml_mode": "strict",
        "confidence_path": None, "metadata": {}}))
    assets = tmp_path / "assets"
    assets.mkdir()
    identity_paths = {}
    for name in ("foldx_binary", "spurs_single_checkpoint", "spurs_multi_checkpoint", "spurs_esm_weights",
                 "discotope_model_archive"):
        path = tmp_path / name
        path.write_text(name)
        identity_paths[name] = path

    class FakeChild:
        stdout = io.StringIO('{"event":"pipeline_start"}\n{"event":"result","result_keys":[]}\n')

        def wait(self):
            return 0

    def fake_popen(args, **_kwargs):
        request_data = json.loads(Path(args[-1]).read_text())
        output = Path(request_data["output_dir"])
        output.mkdir(parents=True, exist_ok=True)
        analysis = {"execution_protocol": {"name": "legacy"}}
        (output / "analysis.json").write_text(json.dumps(analysis))
        (output / "structure_for_energy.pdb").write_text("normalized")
        (output / "card_evidence_completeness.json").write_text(json.dumps({"status": "PASS", "card_count": 1}))
        with zipfile.ZipFile(output / "ecd_epitope_design.zip", "w") as archive:
            archive.write(output / "analysis.json", arcname="analysis.json")
        return FakeChild()

    archive = make_archive(("ecd_epitope_scanner/pipeline.py", "# stub"),
                           ("ecd_epitope_scanner/config.py", '"procaliper": {"python": r"G:\\missing\\python.exe"}'))
    acceptance = tmp_path / "acceptance"
    report = baseline.execute(root=assets, acceptance_dir=acceptance, source_request=source,
        structure=Path(baseline.SOURCE_STRUCTURE), foldx_binary=identity_paths["foldx_binary"],
        model_paths_override=identity_paths, commit_resolver=lambda _root: "a" * 40,
        archive_provider=lambda _commit: archive, popen_factory=fake_popen, max_foldx_root_length=200)
    assert report["status"] == "COMPLETE"
    assert report["elapsed_pipeline_seconds"] is not None
    assert report["elapsed_total_seconds"] >= report["elapsed_pipeline_seconds"]
    assert report["input_sha256"] == baseline.EXPECTED_INPUT_SHA256
    assert report["normalized_structure"]["sha256"]
    assert report["cold_cache_evidence"]["all_possible_repair_cache_paths_absent_before_run"] is True
    assert report["score_cache"]["entries_at_start"] == 0
    assert report["validation"]["worker_model_identity"]["status"] == "NOT_DETERMINED"
    assert json.loads((acceptance / "baseline_report.json").read_text())["status"] == "COMPLETE"
