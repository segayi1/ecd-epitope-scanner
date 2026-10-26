import importlib.util
from pathlib import Path

import pytest


def recovery_module():
    path = Path(__file__).resolve().parents[1] / "workspace_data/af3_binder_batch_20261017/recover_failed_af3.py"
    spec = importlib.util.spec_from_file_location("af3_recovery", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_recovery_only_selects_failed_terminal_jobs():
    module = recovery_module()
    ids = sorted(module.EXPECTED)
    rows = [{"id": job_id, "status": "FAILED" if job_id == ids[0] else "COMPLETE"} for job_id in ids]
    assert [row["id"] for row in module.select_failed(rows)] == [ids[0]]
    rows[0]["status"] = "RUNNING"
    with pytest.raises(ValueError, match="not terminal"):
        module.select_failed(rows)
    with pytest.raises(ValueError, match="identities"):
        module.select_failed(rows[1:])


def test_recovery_pair_cost_excludes_history_but_retains_identity(tmp_path):
    import json
    module = recovery_module()
    manifest = {"records": [{"paired_models": [{"batch_path": "kept", "group": 1, "run": 0}]}]}
    (tmp_path / "paired_model_manifest.json").write_text(json.dumps(manifest))
    history = tmp_path / "prior_attempts" / "old"
    history.mkdir(parents=True)
    manifest["records"][0]["paired_models"][0]["batch_path"] = "history_only"
    (history / "paired_model_manifest.json").write_text(json.dumps(manifest))
    assert module.pair_identities(tmp_path) == {("kept", 1, 0)}


def test_service_reload_refuses_active_jobs_before_stopping_service(tmp_path, monkeypatch):
    import json
    import sqlite3
    module = recovery_module()
    monkeypatch.setattr(module, "TASK", tmp_path)
    monkeypatch.setattr(module, "DB", tmp_path / "jobs.sqlite3")
    (tmp_path / "service_reload_target.json").write_text(json.dumps({"pid": 123, "start_ticks": 456}))
    with sqlite3.connect(module.DB) as conn:
        conn.execute("CREATE TABLE jobs (id TEXT, status TEXT)")
        conn.execute("INSERT INTO jobs VALUES ('active','RUNNING')")
    def unexpected_stop(*args, **kwargs):
        pytest.fail("An active task must prevent the service-stop command")
    monkeypatch.setattr(module.subprocess, "run", unexpected_stop)
    with pytest.raises(RuntimeError, match="another task is active"):
        module.reload_idle_service()
