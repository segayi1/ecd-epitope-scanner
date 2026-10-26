import json
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from write_p0_p3_upgrade_report import _offline_site_static_scan, build_report, render_markdown


def test_static_offline_scan_includes_standalone_report(tmp_path):
    for name in (
        "index.html", "offline_report.html", "assets/app.bundle.js", "assets/workbench.css", "assets/app.css",
        "assets/renderer.js", "assets/vendor/3Dmol-min.js", "assets/vendor/pdbe-molstar-plugin-3.10.1.js",
        "assets/vendor/pdbe-molstar-light-3.10.1.css", "assets/vendor/PDBe-Molstar-LICENSE",
    ):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("local content", encoding="utf-8")
    (tmp_path / "index.html").write_text("//(p-c)^2", encoding="utf-8")
    (tmp_path / "offline_report.html").write_text('<script src="https://remote.example/app.js"></script>', encoding="utf-8")

    present, remote_refs = _offline_site_static_scan(tmp_path)

    assert present
    assert remote_refs == ["https://remote.example/app.js"]

    (tmp_path / "offline_report.html").write_text("// vector comment", encoding="utf-8")
    assert _offline_site_static_scan(tmp_path) == (True, [])


def test_report_requires_each_independent_audit_before_claiming_completion(tmp_path, monkeypatch):
    import write_p0_p3_upgrade_report as report

    manifest = tmp_path / "rerun_manifest.json"
    manifest.write_text(json.dumps({"protocol": "ordinary_v4_sequential_v1", "jobs": [{
        "file": "binder_AIB1.pdb", "source_job_id": "old", "job_id": "new", "workdir": str(tmp_path / "new")
    }]}), encoding="utf-8")
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    (baseline / "case_summary.csv").write_text(
        "binder,source_job_id,chain_length_aa,foldx_candidates,panel_size,job_status\nAIB1,old,90,60,7,COMPLETE\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(report, "_summarize_job", lambda item, *_args: {
        "protein": "AIB1", "source_job_id": item["source_job_id"], "new_job_id": item["job_id"],
        "input_file": item["file"], "status": "COMPLETE_ARTIFACTS_PRESENT", "audit": {"status": "AUDIT_ERROR_REQUIRES_REVIEW"},
        "run_acceptance": {"status": "PASS"}, "bundle_acceptance": {"status": "PASS"},
        "offline_site_browser_audit": {"status": "POLICY_BLOCKED_FILE_URL"},
        "pipeline_elapsed_monotonic_seconds": None,
    })
    monkeypatch.setattr(report, "_external_benchmark", lambda _root: {
        "status": "RUNNING", "metrics_by_order": {"1": None, "2": {}, "3": None},
    })

    result = report.build_report(manifest, baseline, tmp_path / "report", root=tmp_path)

    assert not result["independent_audits_complete"]
    assert result["report_status"] == "SIX_TASK_ARTIFACTS_PRESENT_AUDIT_REVIEW_REQUIRED"


def test_report_waits_for_terminal_external_benchmark(tmp_path, monkeypatch):
    import write_p0_p3_upgrade_report as report

    manifest = tmp_path / "rerun_manifest.json"
    manifest.write_text(json.dumps({"protocol": "ordinary_v4_sequential_v1", "jobs": [{
        "file": "binder_AIB1.pdb", "source_job_id": "old", "job_id": "new", "workdir": str(tmp_path / "new")
    }]}), encoding="utf-8")
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    (baseline / "case_summary.csv").write_text(
        "binder,source_job_id,chain_length_aa,foldx_candidates,panel_size,job_status\nAIB1,old,90,60,7,COMPLETE\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(report, "_summarize_job", lambda item, *_args: {
        "protein": "AIB1", "source_job_id": item["source_job_id"], "new_job_id": item["job_id"],
        "input_file": item["file"], "status": "COMPLETE_ARTIFACTS_PRESENT",
        "audit": {"status": "AUDIT_COMPLETE"}, "run_acceptance": {"status": "PASS"},
        "bundle_acceptance": {"status": "PASS"},
        "offline_site_browser_audit": {"status": "POLICY_BLOCKED_FILE_URL"},
        "pipeline_elapsed_monotonic_seconds": None,
    })
    monkeypatch.setattr(report, "_external_benchmark", lambda _root: {"status": "PILOT_OR_RUNNING"})

    result = report.build_report(manifest, baseline, tmp_path / "report", root=tmp_path)

    assert result["external_benchmark_terminal"] is False
    assert result["report_status"] == "PENDING_EXTERNAL_BENCHMARK"


def test_markdown_reports_recovered_budget_guard_and_spurs_failures():
    markdown = render_markdown({
        "report_status": "COMPLETE_ENGINEERING_AUDIT_WITH_BENCHMARK_LIMITATION",
        "generated_at_local": "2026-10-17T00:00:00-07:00",
        "fresh_terminal_artifact_count": 0,
        "fresh_task_count": 0,
        "fresh_manifest": "rerun_manifest.json",
        "fresh_manifest_protocol": "ordinary_v4_sequential_v1",
        "old_baseline_report_dir": "baseline",
        "tasks": [],
        "offline_browser_policy_block": None,
        "external_benchmark": {
            "status": "COMPLETE_WITH_MODEL_FAILURES",
            "score_blind_candidate_count": 1,
            "prior_attempt_recovery": {
                "snapshot_file": "runner_failure.json",
                "root_cause": "candidate IDs mixed integer and string types",
            },
            "spurs": {
                "attempted_candidate_count": 1,
                "scored_candidate_count": 0,
                "failed_events": [{"candidate_ids": ["candidate-1"], "errors": [{
                    "message": "normalized chain mismatch",
                }]}],
            },
            "foldx": {
                "attempted_candidates": 1,
                "successful_candidates": 0,
                "pair_slots_reserved": 5,
                "approved_pair_slots": 240,
                "pair_proofs_validated": 0,
                "pair_validation_failures": 0,
                "missing_result_candidates": ["candidate-1"],
            },
            "metrics_by_order": {},
            "actual_spurs_training_overlap": "unknown",
            "PTMUL_reference_train_clusters": "unknown",
        },
    })

    recovery_at = markdown.index("runner_failure.json")
    failure_at = markdown.index("normalized chain mismatch")
    table_at = markdown.index("| 突变阶数 |")
    assert recovery_at < failure_at < table_at


def test_report_keeps_running_task_pending_and_never_fills_new_result_from_old_baseline(tmp_path):
    workspace = tmp_path / "workspace_data"
    task = workspace / "new-job"
    task.mkdir(parents=True)
    (task / "execution_trace.jsonl").write_text(json.dumps({
        "event": "V4_STAGE", "stage": "V4_SPURS_COMBINATION", "elapsed_monotonic_seconds": 12.5
    }) + "\n", encoding="utf-8")
    (task / "run_manifest.json").write_text("{}", encoding="utf-8")
    (task / "ecd_epitope_design.zip").write_bytes(b"archive already exists while the job is still running")
    with sqlite3.connect(workspace / "jobs.sqlite3") as conn:
        conn.execute("CREATE TABLE jobs (id TEXT PRIMARY KEY, status TEXT, stage TEXT)")
        conn.execute("INSERT INTO jobs VALUES ('new', 'RUNNING', 'V4_EXPORT')")
    browser_receipt = workspace / "p0_p3_upgrade_20261016" / "offline_browser_audit" / "acceptance.json"
    browser_receipt.parent.mkdir(parents=True)
    browser_receipt.write_text(json.dumps({"results": [{"job_id": "new",
        "status": "PASS_FILE_URL_NO_NETWORK", "opened_from": "unpacked/new/index.html",
        "network_requests": [], "browser_errors": [], "checks": {"regionCount": 7, "designCount": 12},
        "offline_draft": {"status": "LOCAL_DRAFT_NOT_ANALYZED", "downloads": [{}, {}, {}]}}]}), encoding="utf-8")
    block_record = workspace / "p0_p3_upgrade_20261016" / "offline_browser_audit" / "file_protocol_policy_block.json"
    block_record.parent.mkdir(parents=True, exist_ok=True)
    block_record.write_text(json.dumps({"status": "POLICY_BLOCKED", "reason": "file protocol blocked",
        "scope": "remaining task sites"}), encoding="utf-8")
    partial_audit_dir = workspace / "p0_p3_upgrade_20261016" / "offline_browser_audit"
    for suffix in (".png", "_local_draft.json", "_local_fasta.fasta", "_local_cds.fasta"):
        (partial_audit_dir / f"prior{suffix}").write_bytes(b"present")
    manifest = workspace / "rerun_manifest.json"
    manifest.write_text(json.dumps({"protocol": "ordinary_v4_sequential_v1", "status": "SUBMITTED",
        "source_manifest_sha256": "abc", "jobs": [{"file": "test_AIB1.pdb",
            "source_job_id": "old", "job_id": "new", "workdir": str(task)},
            {"file": "test_AIB2.pdb", "source_job_id": "old", "job_id": "prior",
             "workdir": str(workspace / "prior")},
            {"file": "test_AIB2.pdb", "source_job_id": "old", "job_id": "blocked",
             "workdir": str(workspace / "blocked")}]}), encoding="utf-8")
    baseline = tmp_path / "baseline"
    baseline.mkdir()
    (baseline / "case_summary.csv").write_text(
        "binder,source_job_id,chain_length_aa,foldx_candidates,job_status\nAIB1,old,90,60,COMPLETE\n",
        encoding="utf-8")
    result = build_report(manifest, baseline, workspace / "report", root=tmp_path)
    assert result["report_status"] == "PENDING_SIX_FRESH_TASKS"
    row = result["tasks"][0]
    assert row["status"] == "RUNNING:V4_EXPORT"
    assert row["job_store_status"] == {"id": "new", "status": "RUNNING", "stage": "V4_EXPORT"}
    assert row["chain_length_aa"] is None  # never mislabel the historical chain length as a new result
    assert row["historical_source_result"]["chain_length_aa"] == "90"
    assert row["audit"]["status"] == "PENDING_REQUIRED_OUTPUTS"
    assert row["offline_site_browser_audit"]["status"] == "PASS_FILE_URL_NO_NETWORK"
    assert row["offline_site_browser_audit"]["offline_draft_download_count"] == 3
    assert result["tasks"][1]["offline_site_browser_audit"]["status"] == "PRIOR_BROWSER_ARTIFACTS_ONLY"
    assert result["tasks"][1]["offline_site_browser_audit"]["draft_artifacts_present"]
    assert result["tasks"][2]["offline_site_browser_audit"]["status"] == "POLICY_BLOCKED_FILE_URL"
    assert result["offline_browser_policy_block"]["reason"] == "file protocol blocked"
    text = (workspace / "report" / "p0_p3_upgrade_report_zh.md").read_text(encoding="utf-8")
    assert "全部维持PENDING" in text
    assert "稳定性门槛的文献背景" in text
    assert "https://www.nature.com/articles/s41467-025-67609-4" in text
    assert "离线文件页的实时交互验收受浏览器安全策略阻止" in text
    assert "不能據" not in text
