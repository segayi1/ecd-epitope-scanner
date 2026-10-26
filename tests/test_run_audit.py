import hashlib
import json
import zipfile

from ecd_epitope_scanner.joint_design_v4.run_audit import ExecutionTrace, source_snapshot


def test_source_snapshot_archives_only_hashed_executable_sources(tmp_path):
    root = tmp_path / "repo"
    (root / "joint_design_v4").mkdir(parents=True)
    (root / "integrations" / "workers").mkdir(parents=True)
    (root / "vendor" / "spurs").mkdir(parents=True)
    (root / "web").mkdir()
    (root / "workspace_data").mkdir()
    files = {
        "pipeline.py": b"pipeline = True\n",
        "config.py": b"threshold = 1\n",
        "joint_design_v4/engine.py": b"run = True\n",
        "integrations/workers/spurs_gpu_worker.py": b"worker = True\n",
        "vendor/spurs/model.py": b"weights = False\n",
        "web/app.js": b"const app = true;\n",
        "workspace_data/old_result.json": b"historical = True\n",
    }
    for relative, content in files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)

    archive = tmp_path / "run" / "executed_source_snapshot.zip"
    snapshot = source_snapshot(root, {"mode": "ordinary"}, archive)
    rows = {row["path"]: row for row in snapshot["source_files"]}
    assert {"pipeline.py", "config.py", "joint_design_v4/engine.py",
            "integrations/workers/spurs_gpu_worker.py", "vendor/spurs/model.py", "web/app.js"} <= set(rows)
    assert "workspace_data/old_result.json" not in rows
    with zipfile.ZipFile(archive) as bundle:
        assert set(bundle.namelist()) == set(rows)
        for name, row in rows.items():
            assert hashlib.sha256(bundle.read(name)).hexdigest() == row["sha256"]
    assert snapshot["source_archive"]["sha256"] == hashlib.sha256(archive.read_bytes()).hexdigest()
    assert snapshot["effective_parameters"]["mode"] == "ordinary"


def test_execution_trace_is_append_only_and_reports_elapsed_wall_time(tmp_path):
    path = tmp_path / "execution_trace.jsonl"
    trace = ExecutionTrace(path)
    trace.progress("SASA", .25, "surface started")
    trace.progress("CANDIDATE_SEARCH", .5, "candidate search started")
    summary = trace.summary()
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    assert [row["sequence"] for row in rows] == [1, 2, 3]
    assert rows[-1]["stage"] == "CANDIDATE_SEARCH"
    assert summary["elapsed_wall_seconds"] >= 0
    assert summary["queue_delay_captured"] is False
    assert "CPU core-hours" in summary["interpretation"]
