import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_site_first_single_acceptance.py"
spec = importlib.util.spec_from_file_location("site_first_single_acceptance", SCRIPT)
acceptance = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = acceptance
spec.loader.exec_module(acceptance)


def _fixtures(tmp_path, monkeypatch):
    structure = tmp_path / "input.pdb"
    structure.write_bytes(b"fixture structure")
    digest = hashlib.sha256(structure.read_bytes()).hexdigest()
    monkeypatch.setattr(acceptance, "EXPECTED_STRUCTURE_SHA256", digest)
    request = tmp_path / "request.json"
    request.write_text(json.dumps({
        "input_type": "structure", "identifier": "IL7Ra_binder_AF2_35_AIB234.pdb", "chain": "A",
        "analysis_target": "binder", "design_engine": "joint_v4_ml", "ml_mode": "strict",
        "structure_path": "old path", "confidence_path": None, "metadata": {"execution_queue": {"old": True}},
    }), encoding="utf-8")
    return structure, request


def test_preflight_is_read_only_and_validates_source(tmp_path, monkeypatch):
    structure, request = _fixtures(tmp_path, monkeypatch)
    output = tmp_path / "acceptance"
    plan = acceptance.preflight(request, structure, output)
    assert plan["status"] == "PREFLIGHT_PASS"
    assert plan["job_created"] is False
    assert not output.exists()


def test_execute_refuses_repeated_attempt_before_creating_job(tmp_path, monkeypatch):
    structure, request = _fixtures(tmp_path, monkeypatch)
    output = tmp_path / "acceptance"
    output.mkdir()
    (output / "state.json").write_text("{}", encoding="utf-8")
    called = []
    with pytest.raises(FileExistsError):
        acceptance.execute(acceptance_dir=output, source_request=request, structure=structure,
                           store_factory=lambda: called.append(True), data_dir=tmp_path / "data")
    assert not called


class FakeStore:
    def __init__(self):
        self.job = None

    def create(self, request, job_id):
        self.request = request
        self.job = {"id": job_id, "status": "QUEUED", "error": None, "workdir": request.output_dir}
        return self.job

    def get(self, job_id):
        assert self.job["id"] == job_id
        return self.job


class FakeRunner:
    def __init__(self, store):
        self.store = store
        self.executor = SimpleNamespace(shutdown=lambda **kwargs: None)

    def _run(self, job_id):
        workdir = Path(self.store.job["workdir"])
        workdir.mkdir(parents=True, exist_ok=True)
        (workdir / "analysis.json").write_text(json.dumps({"execution_protocol": {"name": "ordinary_site_first_v2"}, "regions": []}), encoding="utf-8")
        (workdir / "sector_recommendations_v4.json").write_text(json.dumps([{"candidate_id": "C1", "formal_eligible": True, "status": "PASSED"}]), encoding="utf-8")
        (workdir / "run_source_snapshot_v1.json").write_text(json.dumps({"source_files_sha256": "source-digest"}), encoding="utf-8")
        from ecd_epitope_scanner.joint_design_v4.run_audit import ExecutionTrace
        trace = ExecutionTrace(workdir / "execution_trace.jsonl")
        trace.record("STAGE_PROGRESS", "SPURS_SINGLE", .5, "fake model stage")
        self.store.job.update(status="COMPLETE", error=None)


def test_success_path_records_live_events_and_report(tmp_path, monkeypatch):
    structure, request = _fixtures(tmp_path, monkeypatch)
    acceptance_dir, data_dir = tmp_path / "acceptance", tmp_path / "data"
    store = FakeStore()
    def validate(workdir, analysis, expected_protocol):
        assert expected_protocol == "ordinary_site_first_v2"
        return {"status": "PASS", "checks": {"protocol": True}, "errors": []}
    report = acceptance.execute(acceptance_dir=acceptance_dir, source_request=request, structure=structure,
        store_factory=lambda: store, runner_factory=FakeRunner, validate_output=validate,
        source_digest=lambda: "source-digest", data_dir=data_dir,
        geometry_check=lambda workdir: {"status": "PASS", "region_count": 1})
    state = json.loads((acceptance_dir / "state.json").read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert state["latest_event"]["stage"] == "SPURS_SINGLE"
    assert (acceptance_dir / "report.md").is_file()
    assert json.loads((acceptance_dir / "report.json").read_text(encoding="utf-8"))["historical_reference"]["runtime_conclusion"].startswith("NOT_DETERMINED")
    assert store.job["id"] == report["job_id"]
    request_obj = store.request
    assert request_obj.metadata["joint_v4_ml"]["execution_protocol"]["name"] == "ordinary_site_first_v2"
    assert Path(request_obj.metadata["joint_v4_ml"]["runtime"]["cache_dir"]).parent == Path(store.job["workdir"])
    assert Path(request_obj.structure_path).read_bytes() == structure.read_bytes()
    assert report["geometry_validation"]["status"] == "PASS"
    cache = report["cache_evidence"]
    assert cache["before_run"]["cold_verified"] is True
    assert cache["before_run"]["foldx_repair_cache"]["path"] == str(Path(store.job["workdir"]) / "fx_shared" / "foldx_repair_cache")
    assert cache["before_run"]["foldx_repair_cache_resolution"].startswith("repair_work_root.parent")


def test_failed_job_writes_review_report_without_retry(tmp_path, monkeypatch):
    structure, request = _fixtures(tmp_path, monkeypatch)
    acceptance_dir, data_dir = tmp_path / "acceptance", tmp_path / "data"
    store = FakeStore()
    class FailingRunner(FakeRunner):
        def _run(self, job_id):
            self.store.job.update(status="FAILED", error="fake failure")
    report = acceptance.execute(acceptance_dir=acceptance_dir, source_request=request, structure=structure,
        store_factory=lambda: store, runner_factory=FailingRunner,
        source_digest=lambda: "source-digest", data_dir=data_dir)
    saved = json.loads((acceptance_dir / "report.json").read_text(encoding="utf-8"))
    assert report["status"] == "NEEDS_REVIEW"
    assert any("fake failure" in error for error in saved["errors"])
    assert store.job["id"] == report["job_id"]
    assert acceptance.preflight(request, structure, acceptance_dir)["execution_already_started"]


def test_validator_accepts_explicit_expected_protocol(tmp_path):
    sys.path.insert(0, str(ROOT.parent))
    from ecd_epitope_scanner.scripts.run_ordinary_adaptive_rerun import validate_case_output
    output = tmp_path / "run"
    output.mkdir()
    analysis = {"execution_protocol": {"name": "ordinary_site_first_v2"}}
    result = validate_case_output(output, analysis, expected_protocol="ordinary_site_first_v2")
    assert result["checks"]["ordinary_protocol"] is True
    legacy = validate_case_output(output, analysis)
    assert legacy["checks"]["ordinary_protocol"] is False


def test_region_principles_uses_primary_alternatives_and_gap_schema(tmp_path):
    rows = [{
        "sector_id": "VS001", "candidate_id": "C1", "status": "PASSED", "formal_eligible": True,
        "principles": ["epitope_support"],
        "alternatives": [{"candidate_id": "C2", "status": "PASSED", "formal_eligible": True,
                          "principles": ["mechanism_pareto", "sector_coverage"]}],
        "unavailable_principles": [],
    }, {
        "sector_id": "VS002", "candidate_id": None, "status": "NO_COMBINATION", "formal_eligible": False,
        "principles": [], "alternatives": [],
        "unavailable_principles": [{"principle": "epitope_support", "state": "NO_QUALIFIED_PRINCIPLE_CANDIDATE",
                                     "reason": "NO_SAFE_EXCLUSIVE_COMBINATION"}],
    }]
    (tmp_path / "sector_recommendations_v4.json").write_text(json.dumps(rows), encoding="utf-8")
    # The real output is keyed by sector on recommendation rows; analysis regions may not carry those rows.
    regions = acceptance._region_principles(tmp_path, {"regions": []})
    assert [row["region"] for row in regions] == ["VS001", "VS002"]
    assert regions[0]["principles"]["epitope_support"]["candidate_id"] == "C1"
    assert regions[0]["principles"]["mechanism_pareto"]["candidate_id"] == "C2"
    assert regions[0]["principles"]["sector_coverage"]["candidate_id"] == "C2"
    assert regions[1]["principles"]["epitope_support"]["gap"] == "NO_SAFE_EXCLUSIVE_COMBINATION"
    assert set(regions[1]["missing_principles"]) == {"mechanism_pareto", "sector_coverage"}


def test_region_principles_maps_real_display_ids_and_reports_technical_failures():
    job = ROOT / "workspace_data" / "d14576822c6f"
    analysis = json.loads((job / "analysis.json").read_text(encoding="utf-8"))
    report = acceptance._region_principles(job, analysis)
    assert [row["region"] for row in report] == [f"VS{i:03d}" for i in range(1, 8)]
    for row in report:
        assert row["principles"]["mechanism_pareto"]["status"] == "technical_failure"
        assert row["principles"]["sector_coverage"]["status"] == "technical_failure"
        assert set(row["missing_principles"]) == {"mechanism_pareto", "sector_coverage"}


def test_geometry_failure_downgrades_validation_and_report(tmp_path, monkeypatch):
    structure, request = _fixtures(tmp_path, monkeypatch)
    acceptance_dir, data_dir = tmp_path / "acceptance", tmp_path / "data"
    store = FakeStore()
    report = acceptance.execute(acceptance_dir=acceptance_dir, source_request=request, structure=structure,
        store_factory=lambda: store, runner_factory=FakeRunner,
        validate_output=lambda *args, **kwargs: {"status": "PASS", "checks": {}, "errors": []},
        source_digest=lambda: "source-digest", data_dir=data_dir,
        geometry_check=lambda workdir: (_ for _ in ()).throw(ValueError("disconnected region")))
    assert report["status"] == "NEEDS_REVIEW"
    assert report["geometry_validation"]["status"] == "NEEDS_REVIEW"
    assert report["validation"]["checks"]["independent_geometry"]["error"].endswith("disconnected region")


def test_independent_geometry_requires_the_verifier_diameter_contract(tmp_path):
    (tmp_path / "run_manifest.json").write_text(json.dumps({"parameters": {"joint_v4_ml": {
        "sectoring": {"maximum_geodesic_diameter_a": 26.0, "maximum_euclidean_diameter_a": 18.0}}}}), encoding="utf-8")
    with pytest.raises(ValueError, match="pinned to 30A/22A"):
        acceptance._independent_geometry_evidence(tmp_path, verifier=SimpleNamespace())


def test_independent_geometry_records_reused_verifier_evidence(tmp_path):
    (tmp_path / "run_manifest.json").write_text(json.dumps({"parameters": {"joint_v4_ml": {
        "sectoring": {"maximum_geodesic_diameter_a": 30.0, "maximum_euclidean_diameter_a": 22.0}}}}), encoding="utf-8")
    seen = []
    verifier = SimpleNamespace(_geometry_evidence=lambda root: seen.append(Path(root)) or {
        "status": "PASS", "region_count": 2, "region_sizes": [5, 6], "surface_residues": 11})
    result = acceptance._independent_geometry_evidence(tmp_path, verifier=verifier)
    assert seen == [tmp_path]
    assert result["status"] == "PASS"
    assert result["maximum_geodesic_diameter_a"] == 30
    assert result["maximum_ca_diameter_a"] == 22


def test_cache_evidence_uses_worker_repair_work_root_parent(tmp_path):
    workdir = tmp_path / "job"
    workdir.mkdir()
    evidence = acceptance._cache_evidence(workdir)
    assert evidence["score_cache"]["path"] == str(workdir / "site_first_score_cache")
    assert evidence["repair_work_root"]["path"] == str(workdir / "fx_shared" / "work")
    assert evidence["foldx_repair_cache"]["path"] == str(workdir / "fx_shared" / "foldx_repair_cache")
    assert evidence["cold_verified"] is True
    (workdir / "fx_shared" / "foldx_repair_cache").mkdir(parents=True)
    after = acceptance._cache_evidence(workdir)
    assert after["foldx_repair_cache"]["exists"] is True
    assert after["cold_verified"] is False
