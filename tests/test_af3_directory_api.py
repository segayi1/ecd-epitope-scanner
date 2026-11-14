from __future__ import annotations

import io
import sys
import types
from pathlib import Path

from ecd_epitope_scanner.api import app
from fastapi.testclient import TestClient

client = TestClient(app)


def _setup(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.jobs import JobStore

    store = JobStore(tmp_path / "workspace")
    monkeypatch.setattr(api, "STORE", store)
    monkeypatch.setattr(api.RUNNER, "store", store)
    submitted = []
    monkeypatch.setattr(api.RUNNER, "submit", submitted.append)
    importer = types.ModuleType("ecd_epitope_scanner.af3_preprocessing")

    def import_af3_directory(directory, output_dir):
        manifest = output_dir / "af3_inputs" / "manifest.json"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text("{}", encoding="utf-8")
        return manifest

    importer.import_af3_directory = import_af3_directory
    monkeypatch.setitem(sys.modules, importer.__name__, importer)
    return api, store, submitted


def _files():
    return [
        ("files", ("seed-123_sample-0/job_seed-123_sample-0_model.cif", io.BytesIO(b"data_model"), "chemical/x-cif")),
        ("files", ("seed-123_sample-0/job_seed-123_sample-0_confidences.json", io.BytesIO(b'{"atom_plddts":[]}'), "application/json")),
    ]


def test_af3_directory_create_saves_manifest_and_binder_fields(monkeypatch, tmp_path):
    api, _, submitted = _setup(monkeypatch, tmp_path)
    response = client.post(
        "/api/jobs/af3-directory",
        data={"analysis_target": "binder", "chain": "A", "official_source_confirmed": "true",
              "execution_protocol": "ordinary_site_first_perturbation_v1"},
        files=_files(),
    )
    assert response.status_code == 202, response.text
    request = response.json()["request"]
    assert request["input_type"] == "af3_upload"
    assert request["analysis_target"] == "binder"
    assert request["chain"] == "A"
    assert request["metadata"]["joint_v4_ml"]["execution_protocol"]["name"] == "ordinary_site_first_perturbation_v1"
    assert request["af3_manifest_path"].replace("/", "\\").endswith("af3_inputs\\manifest.json")
    assert request["assay_profile"] is None
    assert request["metadata"]["official_af3_source_confirmed"] is True
    assert submitted == [response.json()["id"]]


def test_af3_directory_resume_requires_waiting_and_source_confirmation(monkeypatch, tmp_path):
    api, store, submitted = _setup(monkeypatch, tmp_path)
    from ecd_epitope_scanner.pipeline import PipelineRequest

    job = store.create(PipelineRequest(input_type="uniprot", identifier="Q9NZQ7"))
    store.update(job["id"], status="WAITING_FOR_AF3_UPLOAD", stage="WAITING_FOR_AF3_UPLOAD")
    rejected = client.post(f"/api/jobs/{job['id']}/af3-directory", files=_files())
    assert rejected.status_code == 422
    response = client.post(
        f"/api/jobs/{job['id']}/af3-directory",
        data={"official_source_confirmed": "true"},
        files=_files(),
    )
    assert response.status_code == 202, response.text
    request = store.request_for(job["id"])
    assert request.input_type == "af3_upload"
    assert request.accession == "Q9NZQ7"
    assert request.af3_manifest_path
    assert request.output_dir == job["workdir"]
    assert store.get(job["id"])["status"] == "QUEUED"
    assert request.metadata["official_af3_source_confirmed"] is True
    assert submitted == [job["id"]]


def test_af3_directory_rejects_traversal_before_creating_job(monkeypatch, tmp_path):
    _, store, _ = _setup(monkeypatch, tmp_path)
    files = [
        ("files", ("../escape.cif", io.BytesIO(b"data"), "chemical/x-cif")),
        ("files", ("confidence.json", io.BytesIO(b"{}"), "application/json")),
    ]
    response = client.post(
        "/api/jobs/af3-directory",
        data={"official_source_confirmed": "true"},
        files=files,
    )
    assert response.status_code == 422
    assert store.list() == []


def test_af3_directory_accepts_bare_server_filenames(monkeypatch, tmp_path):
    _, _, submitted = _setup(monkeypatch, tmp_path)
    response = client.post("/api/jobs/af3-directory",
        data={"official_source_confirmed": "true", "analysis_target": "binder"},
        files=[("files", ("model_0.cif", io.BytesIO(b"data"), "chemical/x-cif")),
               ("files", ("full_data_0.json", io.BytesIO(b"{}"), "application/json"))])
    assert response.status_code == 202, response.text
    assert submitted == [response.json()["id"]]


def test_af3_directory_rejects_normalized_duplicates_before_writes(monkeypatch, tmp_path):
    _, store, _ = _setup(monkeypatch, tmp_path)
    files = [
        ("files", ("seed-1_sample-0/model.cif", io.BytesIO(b"data"), "chemical/x-cif")),
        ("files", ("SEED-1_SAMPLE-0/MODEL.CIF", io.BytesIO(b"different"), "chemical/x-cif")),
        ("files", ("seed-1_sample-0/full_data.json", io.BytesIO(b"{}"), "application/json")),
    ]
    response = client.post(
        "/api/jobs/af3-directory",
        data={"official_source_confirmed": "true"},
        files=files,
    )
    assert response.status_code == 422
    assert store.list() == []


def test_cli_accepts_af3_directory_as_exclusive_input():
    from ecd_epitope_scanner.cli import build_parser

    args = build_parser().parse_args(["analyze", "--af3-directory", "raw_af3", "--output", "out"])
    assert args.af3_directory == "raw_af3"
    try:
        build_parser().parse_args(
            ["analyze", "--af3-directory", "raw_af3", "--structure", "model.cif", "--output", "out"]
        )
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("--af3-directory must be mutually exclusive with other inputs")


def test_cli_directory_dispatch_passes_manifest_and_selected_fields(monkeypatch, tmp_path):
    import ecd_epitope_scanner.af3_preprocessing as preprocessing
    import ecd_epitope_scanner.cli as cli

    captured = {}

    def import_directory(directory, output_dir):
        captured["source"] = directory
        manifest = output_dir / "af3_inputs" / "manifest.json"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text("{}", encoding="utf-8")
        return manifest

    def fake_run(request, progress_callback):
        captured["request"] = request
        return {"summary": {"ok": True}}

    monkeypatch.setattr(preprocessing, "import_af3_directory", import_directory)
    monkeypatch.setattr(cli, "run_pipeline", fake_run)
    output = tmp_path / "run"
    code = cli.main([
        "analyze", "--af3-directory", "source", "--output", str(output),
        "--analysis-target", "binder", "--chain", "A", "--ml-mode", "preview",
        "--surface-strategy", "surface_perturbation_v1_preview", "--max-sectors", "4",
    ])
    request = captured["request"]
    assert code == 0
    assert captured["source"] == Path("source")
    assert request.input_type == "af3_upload"
    assert request.af3_manifest_path == str(output / "af3_inputs" / "manifest.json")
    assert request.analysis_target == "binder" and request.chain == "A"
    assert request.ml_mode == "preview"
    assert request.surface_strategy == "surface_perturbation_v1_preview"
    assert request.max_sectors == 4


def test_display_projection_keeps_af3_model_diagnostics_for_rendering():
    from ecd_epitope_scanner.joint_design_v4.reporting import project_display_analysis

    report = {
        "method": "AF3_ENSEMBLE_CORE", "selected_model_id": "m001",
        "baselines": {"original_rank1": "m000", "ptm_best": 0.83},
        "models": [{"model_id": "m000", "eligible": False, "exclusion_reasons": ["outside_Q"], "f70": 0.8}],
        "model_diagnostics": [{"model_id": "m000", "has_clash": False, "weight": 0.0, "pae_status": "mapped", "sha256": {"structure": "hash"}}],
        "per_residue_diagnostics": {"10": {"mapped_plddt": [80.0], "low_confidence_models": []}},
    }
    display = project_display_analysis({"methods": {"af3_preprocessing": report}, "joint_v4_ml": {"selected_panel": {}}})
    compact = display["methods"]["af3_preprocessing"]
    for field in ("method", "selected_model_id", "baselines", "models", "model_diagnostics"):
        assert compact[field] == report[field]
    assert "per_residue_diagnostics" not in compact
    assert compact["per_residue_diagnostics_reference"] == {
        "file": "af3_selection_report.json", "json_path": ["per_residue_diagnostics"], "count": 1,
    }
    live = Path(__file__).resolve().parents[1] / "web" / "app.js"
    offline = live.with_name("offline_template.html")
    assert "state.result.methods?.af3_preprocessing" in live.read_text(encoding="utf-8")
    assert "Object.entries(af3.baselines || {})" in live.read_text(encoding="utf-8")
    assert "(D.methods||{}).af3_preprocessing" in offline.read_text(encoding="utf-8")
    for page in (live, offline):
        content = page.read_text(encoding="utf-8")
        assert "consensus_score" in content and "core_max_rmsd_flag" in content
