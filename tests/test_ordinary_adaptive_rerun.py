import hashlib
import importlib.util
import json
from pathlib import Path
import sqlite3
import sys
import types

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/run_ordinary_adaptive_rerun.py"
spec = importlib.util.spec_from_file_location("ordinary_rerun", SCRIPT)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def frozen_case(tmp_path):
    structure = tmp_path / "input.cif"
    confidence = tmp_path / "confidence.json"
    structure.write_text("structure", encoding="utf-8")
    confidence.write_text("{}", encoding="utf-8")
    return {"case_id": "case1", "input_type": "af3_upload", "structure_path": str(structure),
            "structure_sha256": hashlib.sha256(structure.read_bytes()).hexdigest(),
            "confidence_path": str(confidence),
            "confidence_sha256": hashlib.sha256(confidence.read_bytes()).hexdigest(),
            "source_hashes": {"structure": "s", "confidence": "c"},
            "identifier": "sample", "chain": "A"}


def test_preflight_rejects_changed_frozen_input(tmp_path):
    case = frozen_case(tmp_path)
    case["confidence_path"] = None
    case["confidence_sha256"] = None
    case["input_type"] = "structure"
    Path(case["structure_path"]).write_text("changed", encoding="utf-8")
    with pytest.raises(ValueError, match="hash changed"):
        runner.preflight({"case_count": 13, "rerun_cases": [dict(case, case_id=f"case{i}") for i in range(13)]})


@pytest.mark.parametrize("failure", ["saved_failure", "new_failure", "validation", "source_changed", "success"])
def test_batch_stops_on_failed_case_and_records_interruption(tmp_path, monkeypatch, failure):
    case = frozen_case(tmp_path)
    cases = [dict(case, case_id=f"case{i}", name=f"case{i}", source_run_id=f"src{i}") for i in range(13)]
    cases[-1].pop("source_hashes")
    manifest = {"case_count": 13, "rerun_cases": cases,
                "protocol": {"engine": "ordinary_v4_sequential_v1"},
                "comparison_baselines": [], "comparison_dimensions": []}
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"status": "PASS", "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                                   "runner_sha256": runner.sha256(SCRIPT),
                                   "production_source_sha256": runner.production_source_sha256()}), encoding="utf-8")
    state_path = tmp_path / "state.json"
    state_path.write_text(json.dumps({"manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                                     "runner_sha256": runner.sha256(SCRIPT),
                                     "production_source_sha256": runner.production_source_sha256(),
                                     "jobs": {"case0": {"job_id": "old", "status": "FAILED"}} if failure == "saved_failure" else {}}), encoding="utf-8")
    monkeypatch.setattr(runner, "MANIFEST", manifest_path)
    monkeypatch.setattr(runner, "STATE", state_path)
    monkeypatch.setattr(runner, "OUT", tmp_path / "report")
    monkeypatch.setattr(runner, "ROOT", tmp_path)
    stores = []
    class Store:
        db_path = tmp_path / "jobs.sqlite3"
        def __init__(self, *_): self.created = []; stores.append(self)
        def create(self, req, job_id=None):
            job = {"id": job_id or f"j{len(self.created)}", "workdir": str(tmp_path / (job_id or f"j{len(self.created)}"))}
            self.created.append(job); return job
        def get(self, jid): return {"id": jid, "status": "COMPLETE" if failure in {"validation", "success"} else "FAILED",
                                    "error": None, "workdir": str(tmp_path / jid)}
    class Runner:
        def __init__(self, store): self.store = store; self.executor = types.SimpleNamespace(shutdown=lambda **_: None)
        def _run(self, jid): pass
    module = types.ModuleType("ecd_epitope_scanner")
    jobs = types.ModuleType("ecd_epitope_scanner.jobs"); jobs.JobStore = Store; jobs.JobRunner = Runner
    pipeline = types.ModuleType("ecd_epitope_scanner.pipeline")
    class Request:
        def __init__(self, **kwargs): self.kwargs = kwargs
    pipeline.PipelineRequest = Request
    monkeypatch.setitem(sys.modules, "ecd_epitope_scanner", module)
    monkeypatch.setitem(sys.modules, "ecd_epitope_scanner.jobs", jobs)
    monkeypatch.setitem(sys.modules, "ecd_epitope_scanner.pipeline", pipeline)
    monkeypatch.setattr(runner, "capture_post_validation", lambda *_: {
        "status": "PASS" if failure == "success" else "NEEDS_REVIEW", "errors": ["mock evidence missing"]})
    if failure == "source_changed":
        source_hashes = iter([json.loads(receipt.read_text())["production_source_sha256"], "changed"])
        monkeypatch.setattr(runner, "production_source_sha256", lambda: next(source_hashes))
        with pytest.raises(RuntimeError, match="source changed"):
            runner.execute(manifest, receipt, state_path, resume=True)
        saved = json.loads(state_path.read_text())
        assert saved["status"] == "INTERRUPTED" and "source changed" in saved["error"]
        assert stores[0].created == []
        report = json.loads((tmp_path / "report/report.json").read_text())
        assert report["status"] == "INTERRUPTED" and report["acceptance_gate"] == "NEEDS_REVIEW"
        assert "source changed" in report["error"]
        return
    result = runner.execute(manifest, receipt, state_path, resume=True)
    if failure == "success":
        assert result["status"] == "COMPLETE" and len(result["jobs"]) == 13
        assert len(stores[0].created) == 13
    else:
        assert result["status"] == "HALTED_FOR_REVIEW" and result["blocked_case_id"] == "case0"
        assert len(result["jobs"]) == 1
        assert len(stores[0].created) == (0 if failure == "saved_failure" else 1)
        report = json.loads((tmp_path / "report/report.json").read_text())
        assert report["status"] == "HALTED_FOR_REVIEW" and report["acceptance_gate"] == "NEEDS_REVIEW"


def test_preflight_rejects_duplicate_case_ids(tmp_path):
    with pytest.raises(ValueError, match="IDs must be unique"):
        runner.preflight({"case_count": 13, "rerun_cases": [frozen_case(tmp_path)] * 13})


def test_receipt_requires_production_source_hash(tmp_path, monkeypatch):
    manifest = {"case_count": 13, "rerun_cases": []}
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    monkeypatch.setattr(runner, "MANIFEST", path)
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"status": "PASS", "manifest_sha256": runner.sha256(path),
                                  "runner_sha256": runner.sha256(SCRIPT)}), encoding="utf-8")
    with pytest.raises(ValueError, match="production source SHA256"):
        runner.execute(manifest, receipt, tmp_path / "state.json")


def test_active_job_guard_runs_before_jobstore_initialization(tmp_path, monkeypatch):
    case = frozen_case(tmp_path)
    cases = [dict(case, case_id=f"case{i}", source_run_id=f"src{i}") for i in range(13)]
    manifest = {"case_count": 13, "rerun_cases": cases}
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"status": "PASS", "manifest_sha256": runner.sha256(manifest_path),
        "runner_sha256": runner.sha256(SCRIPT), "production_source_sha256": runner.production_source_sha256()}), encoding="utf-8")
    monkeypatch.setattr(runner, "ROOT", tmp_path)
    monkeypatch.setattr(runner, "MANIFEST", manifest_path)
    db = tmp_path / "workspace_data/jobs.sqlite3"
    db.parent.mkdir()
    with sqlite3.connect(db) as connection:
        connection.execute("CREATE TABLE jobs (id TEXT, status TEXT)")
        connection.execute("INSERT INTO jobs VALUES ('unrelated', 'RUNNING')")
    monkeypatch.setitem(sys.modules, "ecd_epitope_scanner", types.ModuleType("ecd_epitope_scanner"))
    jobs = types.ModuleType("ecd_epitope_scanner.jobs")
    jobs.JobStore = lambda *_: pytest.fail("JobStore initialized before active-job check")
    jobs.JobRunner = object
    monkeypatch.setitem(sys.modules, "ecd_epitope_scanner.jobs", jobs)
    with pytest.raises(RuntimeError, match="Active jobs"):
        runner.execute(manifest, receipt, tmp_path / "state.json")


def test_request_preserves_sequence_and_ecd_from_saved_request(tmp_path, monkeypatch):
    source = tmp_path / "request.json"
    source.write_text(json.dumps({"sequence": "ACDE", "ecd_start": 3, "ecd_end": 9,
                                 "metadata": {"disable_external_engines": True}}), encoding="utf-8")
    case = frozen_case(tmp_path)
    case["source_request"] = str(source)
    monkeypatch.setitem(sys.modules, "ecd_epitope_scanner.pipeline", types.SimpleNamespace(
        PipelineRequest=lambda **kwargs: kwargs))
    request = runner.request_for(case, Path(case["structure_path"]), Path(case["confidence_path"]))
    assert (request["sequence"], request["ecd_start"], request["ecd_end"]) == ("ACDE", 3, 9)
    assert request["input_type"] == "af3_upload"
    assert request["confidence_path"] == case["confidence_path"]
    assert request["metadata"] == {"ordinary_adaptive_rerun_case": "case1"}
    case.pop("identifier")
    case["name"] = "binder_name"
    assert runner.request_for(case, Path(case["structure_path"]), Path(case["confidence_path"]))["identifier"] == "binder_name"


def test_report_contains_case_level_comparison_and_marks_missing_baseline(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    case = {"case_id": "c1", "name": "sample", "source_run_id": "old1", "source_run_dir": str(source)}
    state = {"status": "COMPLETE", "manifest_sha256": "m", "jobs": {"c1": {
        "job_id": "new1", "status": "COMPLETE", "result_evidence": {"analysis": "missing.json",
        "region_count": 2, "design_count": 3}}}}
    monkeypatch.setattr(runner, "OUT", tmp_path / "out")
    runner.write_reports(state, {"rerun_cases": [case], "protocol": {"engine": "ordinary"},
                                 "comparison_baselines": [], "comparison_dimensions": []})
    report = json.loads((tmp_path / "out/report.json").read_text(encoding="utf-8"))
    assert report["cross_comparison"][0]["baseline_status"] == "MISSING"
    assert report["cross_comparison"][0]["rerun_status"] == "COMPLETE"


def test_cross_compare_matches_legacy_foldx_signatures_when_sequence_matches(tmp_path, monkeypatch):
    old = tmp_path / "paired_model_manifest.json"
    old_data = {"execution_identity": {"sequence": "ACDE"}, "records": [
        {"mutations": ["A1V", "C2D"], "ddg_runs_kcal_mol": [0.5, 0.7, 0.6]},
        {"mutations": ["D3E"], "ddg_runs_kcal_mol": [1.0, 1.2, 0.8]}]}
    old.write_text(json.dumps(old_data), encoding="utf-8")
    analysis_path = tmp_path / "analysis.json"
    analysis_path.write_text(json.dumps({"request": {"sequence": "ACDE"}, "designs": [
        {"mutations": [{"position": 1, "wt": "A", "mutant": "V"}, {"position": 2, "wt": "C", "mutant": "D"}],
         "stability": {"spurs_multi_ddg_kcal_mol": 0.4, "foldx_ddg_kcal_mol": 0.55}},
        {"mutations": [{"position": 4, "wt": "E", "mutant": "F"}], "stability": {}}],
        "regions": []}), encoding="utf-8")
    case = {"case_id": "c", "source_run_id": "s", "source_run_dir": str(tmp_path), "chain": "A",
        "structure_path": str(tmp_path / "structure.cif"), "request_contract": {},
        "baseline_evidence": {"legacy_summary": {"formal": 4, "foldx_unique": 80},
            "paired_manifests": [{"path": str(old), "sha256": runner.sha256(old)}]}}
    (tmp_path / "structure.cif").write_text("x", encoding="utf-8")
    monkeypatch.setattr(runner, "_chain_sequence", lambda *_: "ACDE")
    compared = runner.cross_compare_case(case, {"result_evidence": {"analysis": str(analysis_path)}})
    assert compared["shared_signatures"] == ["A1V|C2D"]
    assert compared["new_reported_signatures"] == ["E4F"]
    assert compared["legacy_scored_not_in_rerun"] == ["D3E"]
    assert compared["same_combo_scores"][0]["legacy_foldx_median"] == 0.6
    assert compared["same_combo_scores"][0]["rerun_spurs"] == 0.4
    case.update(name="sample", source_run_id="old", source_run_dir=str(tmp_path))
    monkeypatch.setattr(runner, "OUT", tmp_path / "report")
    state = {"status": "COMPLETE", "manifest_sha256": "m", "jobs": {"c": {
        "job_id": "new", "status": "COMPLETE", "result_evidence": {"analysis": str(analysis_path)}}}}
    runner.write_reports(state, {"rerun_cases": [case], "protocol": {"engine": "ordinary"},
        "comparison_baselines": [], "comparison_dimensions": []})
    csv_text = (tmp_path / "report/report.csv").read_text(encoding="utf-8-sig")
    markdown = (tmp_path / "report/report.md").read_text(encoding="utf-8")
    assert "A1V|C2D" in csv_text and "same_combo_scores" in csv_text
    assert "Legacy scored signatures absent" in markdown and "legacy_foldx_median" in markdown


def test_postvalidation_rejects_overbudget_ledger(tmp_path):
    (tmp_path / "ordinary_candidate_ledger.json").write_text(json.dumps({
        "spurs": {"unique_requests": 1001, "technical_failure": False},
        "foldx": {"attempted_slots": 301, "technical_failure": False}, "foldx_results": []}), encoding="utf-8")
    result = runner.validate_case_output(tmp_path, {"execution_protocol": {"name": "ordinary_three_principle_adaptive_v1"},
        "design_summary": {"sector_recommendations": []}})
    assert result["status"] == "NEEDS_REVIEW"
    assert not result["checks"]["spurs_budget"] and not result["checks"]["foldx_budget"]
    ranked = runner.validate_case_output(tmp_path, {"execution_protocol": {"name": "ordinary_three_principle_adaptive_v1"}},
                                         spurs_request_limit=10000)
    assert ranked["checks"]["spurs_budget"] and not ranked["checks"]["foldx_budget"]


def test_postvalidation_reads_real_recommendation_csv_large_provenance_field(tmp_path):
    import csv

    job = Path(__file__).resolve().parents[1] / "workspace_data" / "06710f1f861b"
    source_csv = job / "sector_recommendations_v4.csv"
    source_json = job / "sector_recommendations_v4.json"
    original_limit = csv.field_size_limit()
    csv.field_size_limit(16 * 1024 * 1024)
    try:
        with source_csv.open(encoding="utf-8-sig", newline="") as stream:
            largest_field = max(len(value or "") for row in csv.DictReader(stream) for value in row.values())
    finally:
        csv.field_size_limit(original_limit)
    assert largest_field > original_limit

    (tmp_path / source_csv.name).write_bytes(source_csv.read_bytes())
    (tmp_path / source_json.name).write_bytes(source_json.read_bytes())
    result = runner.validate_case_output(tmp_path, {
        "execution_protocol": {"name": "ordinary_three_principle_adaptive_v1"}})
    assert not any("sector recommendation CSV validation failed" in error
                   or "sector recommendation JSON/CSV membership mismatch" in error
                   for error in result["errors"])
    assert csv.field_size_limit() == original_limit


def test_postvalidation_rejects_real_site_first_zeroed_panel_signatures():
    job = Path(__file__).resolve().parents[1] / "workspace_data" / "fdbde4aa24d7"
    analysis = json.loads((job / "analysis.json").read_text(encoding="utf-8"))
    result = runner.validate_case_output(job, analysis, expected_protocol="ordinary_site_first_v2")
    assert result["status"] == "NEEDS_REVIEW"
    assert result["checks"]["panel_hypothesis_signatures"] is False
    assert any("panel signatures disagree with selected candidate responses" in error
               for error in result["errors"])


def test_site_first_panel_signature_validation_requires_artifacts(tmp_path):
    errors = runner._panel_hypothesis_signature_errors(tmp_path)
    assert any("selected_panel_v4.json" in error for error in errors)
    assert any("panel_identifiability_v41.json" in error for error in errors)
    assert any("panel_candidates.json" in error for error in errors)


def test_site_first_display_audit_rejects_sensitive_and_missing_evidence(tmp_path):
    audit = {"sampling_decision_diagnostic_status": "COMPLETE", "sampling_decision_sensitive": False,
             "recommendation_risk_label": "NO_WITHIN_RUN_GATE_CROSSING",
             "foldx_repeat_slot_values_kcal_mol": [.1, .2, .3]}
    candidate = {"candidate_id": "C1", "formal_eligible": True, "display_eligible": True,
                 "qualification": "FORMAL_POOL", "mutations": [{"review_required": False}],
                 "impact_evidence": {"policy": "regional_mechanism_pareto_v1", "missing_axes": [],
                                     "stability_repeat_audit": audit}}
    card = {"candidate_id": "C1", "status": "PASSED", "formal_eligible": True,
            "stability_repeat_diagnostic": dict(audit)}
    panel = {"selected_candidate_ids": ["C1"], "sector_recommendations": [card]}
    for name, value in (("panel_candidates.json", [candidate]), ("selected_panel_v4.json", panel),
                        ("sector_recommendations_v4.json", [card])):
        (tmp_path / name).write_text(json.dumps(value), encoding="utf-8")
    assert runner._site_first_display_errors(tmp_path) == []

    audit["foldx_repeat_slot_values_kcal_mol"] = [2.1, .2, .3] + [.2] * 7
    (tmp_path / "panel_candidates.json").write_text(json.dumps([candidate]), encoding="utf-8")
    assert any("cross" in error for error in runner._site_first_display_errors(tmp_path))
    audit["foldx_repeat_slot_values_kcal_mol"] = [.1, .2, .3]
    candidate["impact_evidence"]["missing_axes"] = ["charge_delta"]
    (tmp_path / "panel_candidates.json").write_text(json.dumps([candidate]), encoding="utf-8")
    assert any("evidence" in error for error in runner._site_first_display_errors(tmp_path))

    candidate["impact_evidence"]["missing_axes"] = []
    (tmp_path / "panel_candidates.json").write_text(json.dumps([candidate]), encoding="utf-8")
    card["alternatives"] = [{"candidate_id": "C2", "status": "EVIDENCE_LIMITED", "formal_eligible": True}]
    (tmp_path / "sector_recommendations_v4.json").write_text(json.dumps([card]), encoding="utf-8")
    assert any("C2" in error for error in runner._site_first_display_errors(tmp_path))


def test_site_first_display_accepts_only_explicit_valid_ten_slot_confirmation(tmp_path):
    def write_artifacts(audit, diagnostic=None):
        diagnostic = dict(audit) if diagnostic is None else diagnostic
        candidate = {"candidate_id": "C1", "formal_eligible": True, "display_eligible": True,
                     "qualification": "FORMAL_POOL", "mutations": [{"review_required": False}],
                     "impact_evidence": {"policy": "regional_mechanism_pareto_v1", "missing_axes": [],
                                         "stability_repeat_audit": audit}}
        card = {"candidate_id": "C1", "status": "PASSED", "formal_eligible": True,
                "stability_repeat_diagnostic": diagnostic}
        panel = {"selected_candidate_ids": ["C1"], "sector_recommendations": [card]}
        for name, value in (("panel_candidates.json", [candidate]), ("selected_panel_v4.json", panel),
                            ("sector_recommendations_v4.json", [card])):
            (tmp_path / name).write_text(json.dumps(value), encoding="utf-8")

    audit = {"sampling_decision_diagnostic_status": "COMPLETE", "sampling_decision_sensitive": True,
             "recommendation_risk_label": "TEN_SLOT_MEDIAN_CONFIRMED",
             "repeat_confirmation_policy": "ten_slot_median_v1",
             "repeat_confirmation_status": "CONFIRMED_PASS", "foldx_repeat_slots": list(range(10)),
             "foldx_repeat_slot_values_kcal_mol": [2.1, .2, .3] + [.2] * 7,
             "foldx_repeat_gate_kcal_mol": 2.0, "cross_execution_conflict": None}
    write_artifacts(audit)
    assert runner._site_first_display_errors(tmp_path) == []

    inconsistent_card = dict(audit, foldx_repeat_slot_values_kcal_mol=[.2] * 10)
    write_artifacts(audit, inconsistent_card)
    assert any("disagrees with raw evidence" in error
               for error in runner._site_first_display_errors(tmp_path))

    no_marker = dict(audit, repeat_confirmation_policy=None, repeat_confirmation_status=None)
    write_artifacts(no_marker)
    assert any("crosses the FoldX gate" in error for error in runner._site_first_display_errors(tmp_path))

    invalid = [
        dict(audit, foldx_repeat_slot_values_kcal_mol=[2.1, .2, .3] + [2.2] * 5 + [.2] * 2),
        dict(audit, foldx_repeat_slots=list(range(9)) + [8]),
        dict(audit, foldx_repeat_slots=list(range(9))),
    ]
    for bad_audit in invalid:
        write_artifacts(bad_audit)
        assert runner._site_first_display_errors(tmp_path)


@pytest.mark.parametrize("slots,scores,median", [
    ([0, 1], [.1, .2], .15),
    ([0, 1, 2], [2.1, .1, .1], .1),
    (list(range(10)), [2.1, 2.1, .1] + [.1] * 7, .1),
    ([0, 1, 3], [.1, .1, .1], .1),
])
def test_postvalidation_rejects_foldx_pass_with_invalid_protocol(tmp_path, slots, scores, median):
    (tmp_path / "ordinary_candidate_ledger.json").write_text(json.dumps({
        "spurs": {"unique_requests": 1, "technical_failure": False},
        "foldx": {"attempted_slots": 2, "technical_failure": False},
        "foldx_results": [{"status": "PASS", "candidate_id": "c", "slots": slots,
            "scores": scores, "median": median, "paired_models": []}]}), encoding="utf-8")
    result = runner.validate_case_output(tmp_path, {"execution_protocol": {"name": "ordinary_three_principle_adaptive_v1"},
        "design_summary": {"sector_recommendations": []}})
    assert result["status"] == "NEEDS_REVIEW"
    assert any("invalid PASS FoldX evidence" in error for error in result["errors"])


@pytest.mark.parametrize('string_protocol', [False, True])
def test_report_accepts_frozen_manifest_protocol_contract(tmp_path, monkeypatch, string_protocol):
    manifest = json.loads(runner.MANIFEST.read_text(encoding='utf-8'))
    if string_protocol:
        manifest['protocol'] = 'ordinary_site_first_v2'
    monkeypatch.setattr(runner, 'OUT', tmp_path)
    monkeypatch.setattr(runner, 'cross_compare_case', lambda case, *_: {
        'case_id': case['case_id'], 'source_run_id': case['source_run_id'],
        'baseline_status': 'UNKNOWN', 'mapping_status': 'UNKNOWN',
        'shared_signatures': [], 'new_reported_signatures': [],
        'legacy_scored_not_in_rerun': [], 'same_combo_scores': []})
    runner.write_reports({'status': 'RUNNING', 'jobs': {}, 'manifest_sha256': 'test'}, manifest)
    expected = 'ordinary_site_first_v2' if string_protocol else manifest['protocol']['engine']
    assert expected in (tmp_path / 'report.md').read_text(encoding='utf-8')


def test_report_failure_preserves_case_and_does_not_mask_execution_error(tmp_path, monkeypatch):
    def fail(*_):
        raise OSError('report disk unavailable')
    monkeypatch.setattr(runner, 'write_reports', fail)
    state = {'status': 'INTERRUPTED', 'error': 'source changed',
             'jobs': {'pd1': {'status': 'PASS', 'job_status': 'COMPLETE'}}}
    runner.publish_reports(state, {}, tmp_path / 'state.json')
    saved = json.loads((tmp_path / 'state.json').read_text(encoding='utf-8'))
    assert saved['jobs']['pd1']['status'] == 'PASS'
    assert saved['error'] == 'source changed'
    assert saved['report_error'] == 'OSError: report disk unavailable'


def test_regional_coverage_reports_gaps_separately_from_technical_status():
    ledger = {'queues': [
        {'sector_id': 'VS006', 'principle': 'mechanism_pareto', 'state': 'candidate_failures'},
        {'sector_id': 'VS004', 'principle': 'sector_coverage', 'state': 'budget_limited'},
        {'sector_id': 'VS004', 'principle': 'epitope_support', 'state': 'not_applicable'}]}
    result = runner.regional_principle_coverage(ledger)
    assert result['status'] == 'INCOMPLETE'
    assert [(x['sector_id'], x['state']) for x in result['gaps']] == [('VS006', 'candidate_failures'), ('VS004', 'budget_limited')]
    assert runner.regional_principle_coverage({})['status'] == 'UNKNOWN'
    assert runner.regional_principle_coverage({'queues': [None]})['status'] == 'UNKNOWN'
    assert runner.regional_principle_coverage({'queues':[{'sector_id':'S', 'principle':'p', 'state':'not_applicable'}, {'sector_id':'S', 'principle':'q', 'state':'complete'}]})['status'] == 'COMPLETE'
