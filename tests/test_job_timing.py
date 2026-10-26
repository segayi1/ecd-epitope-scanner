import asyncio
import datetime as dt
import json
from pathlib import Path
import shutil
import subprocess

import pytest
from fastapi.testclient import TestClient

from ecd_epitope_scanner import api, jobs
from ecd_epitope_scanner.job_timing import read_timing, remaining_time, timing_key
from ecd_epitope_scanner.pipeline import PipelineRequest

NOW = dt.datetime(2026, 10, 26, 12, tzinfo=dt.UTC)


def trace(workdir, execution, start, duration=None):
    rows = [dict(execution_id=execution, event="PIPELINE_STARTED", stage=None,
                 timestamp_utc=start.isoformat(), elapsed_monotonic_seconds=0),
            dict(execution_id=execution, event="V4_STAGE", stage="V4_FOLDX_FIRST_PASS",
                 timestamp_utc=(start + dt.timedelta(seconds=20)).isoformat(), elapsed_monotonic_seconds=20)]
    if duration:
        rows.append(dict(execution_id=execution, event="ANALYSIS_READY_FOR_EXPORT", stage="EXPORTING",
                         timestamp_utc=(start + dt.timedelta(seconds=duration)).isoformat(), elapsed_monotonic_seconds=duration))
    path = Path(workdir) / "execution_trace.jsonl"
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")
    return path


def test_stage_estimate_ignores_progress_and_excludes_duplicate_executions(tmp_path):
    current_dir = tmp_path / "current"
    current_dir.mkdir()
    current = read_timing(trace(current_dir, "current", NOW - dt.timedelta(seconds=80)).parent)
    references = []
    for index, duration in enumerate([300, 360, 420]):
        folder = tmp_path / str(index)
        folder.mkdir()
        reference = read_timing(trace(folder, str(index), NOW - dt.timedelta(days=1), duration).parent)
        reference["total_seconds"] = duration
        references.append(reference)
    result = remaining_time(current, references + references, NOW)
    assert result["status"] == "ESTIMATED" and result["reference_count"] == 3
    assert result["range_seconds"] == [210, 360]
    assert result["remaining_seconds"] == 300
    assert remaining_time(current, references[:2], NOW)["status"] == "INSUFFICIENT_HISTORY"
    overrun = remaining_time(current, references, NOW + dt.timedelta(hours=1))
    assert overrun["status"] == "OVERRUN" and overrun["remaining_seconds"] is None


def test_trace_uses_latest_attempt_and_tolerates_partial_corrupt_lines(tmp_path):
    path = trace(tmp_path, "old", NOW - dt.timedelta(days=1), 900)
    old = path.read_text()
    trace(tmp_path, "new", NOW - dt.timedelta(seconds=80))
    path.write_text(old + path.read_text() + '{"event":\n[]\n', encoding="utf-8")
    result = read_timing(tmp_path)
    assert result["execution_id"] == "new" and result["elapsed_seconds"] == 20
    assert "EXPORTING" not in result["stage_starts"]
    assert read_timing(tmp_path / "absent") is None


def test_timing_groups_separate_disabled_models_af3_and_binder():
    ordinary = {"input_type": "structure"}
    assert timing_key(ordinary) != timing_key({**ordinary, "analysis_target": "binder"})
    assert timing_key(ordinary) != timing_key({**ordinary, "metadata": {"disable_external_engines": True}})
    assert timing_key({"input_type": "af3_upload"}) != timing_key({"input_type": "af3_upload", "af3_manifest_path": "models.json"})


def test_api_and_sse_publish_estimate_without_changing_database(tmp_path, monkeypatch):
    store = jobs.JobStore(tmp_path)
    monkeypatch.setattr(api, "STORE", store)
    for index, duration in enumerate([300, 360, 420]):
        start = dt.datetime.now(dt.UTC) - dt.timedelta(days=index + 1)
        monkeypatch.setattr(jobs, "_now", lambda start=start: (start - dt.timedelta(seconds=10)).isoformat())
        job = store.create(PipelineRequest(input_type="structure"))
        trace(job["workdir"], f"history-{index}", start, duration)
        monkeypatch.setattr(jobs, "_now", lambda start=start, duration=duration: (start + dt.timedelta(seconds=duration)).isoformat())
        store.update(job["id"], status="COMPLETE")
    monkeypatch.setattr(jobs, "_now", lambda: dt.datetime.now(dt.UTC).isoformat())
    current = store.create(PipelineRequest(input_type="structure"))
    trace(current["workdir"], "current", dt.datetime.now(dt.UTC) - dt.timedelta(seconds=80))
    store.update(current["id"], status="RUNNING", stage="PROCALIPER_ANALYSIS", progress=.62)
    before = store.get(current["id"])
    shown = TestClient(api.app).get(f"/api/jobs/{current['id']}?view=display").json()
    assert shown["time_estimate"]["status"] == "ESTIMATED"
    assert shown["time_estimate"]["stage"] == "V4_FOLDX_FIRST_PASS"
    assert store.get(current["id"]) == before

    async def event():
        response = await api.job_events(current["id"])
        payload = await anext(response.body_iterator)
        await response.body_iterator.aclose()
        return json.loads(payload.removeprefix("data: "))

    assert asyncio.run(event())["time_estimate"]["reference_count"] == 3
    for status in ("QUEUED", "WAITING_FOR_AF3_UPLOAD", "FAILED", "COMPLETE"):
        assert store.time_estimate({**before, "status": status}) is None
    assert store.time_estimate({**before, "message": "Starting pipeline"})["status"] == "INSUFFICIENT_HISTORY"
    assert store.time_estimate({**before, "owner_id": "someone-else"})["reference_count"] == 0


def test_eta_renderer_handles_ranges_queue_overrun_and_terminal(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for frontend checks")
    app = (Path(__file__).resolve().parents[1] / "web/app.js").read_text(encoding="utf-8")
    function = "function setJob" + app.split("function setJob", 1)[1].split("function adoptResult", 1)[0]
    checks = r"""
const assert=require('node:assert/strict'),state={lang:'zh'},elements={};
const $=id=>elements[id]||=( {style:{}} ),t=(key,fallback)=>fallback;
setJob('RUNNING','FoldX',.62,null,{status:'ESTIMATED',stage:'V4_FOLDX_FIRST_PASS',range_seconds:[600,1200],reference_count:5});
assert(elements['job-eta'].textContent.includes('10–20 分钟'));assert(elements['job-eta'].textContent.includes('FoldX首轮'));
assert.equal(elements['job-eta'].hidden,false);
setJob('RUNNING','FoldX',.62,null,{status:'OVERRUN',reference_count:5});assert(elements['job-eta'].textContent.includes('超出历史'));
setJob('RUNNING','FoldX',.62);assert(elements['job-eta'].textContent.includes('暂无法估计'));
setJob('QUEUED','queue');assert(elements['job-eta'].textContent.includes('排队'));
setJob('COMPLETE','done',1);assert.equal(elements['job-eta'].hidden,true);
setJob('WAITING_FOR_AF3_UPLOAD','waiting');assert.equal(elements['job-eta'].hidden,true);
"""
    script = tmp_path / "eta.cjs"
    script.write_text(function + checks, encoding="utf-8")
    subprocess.run([node, str(script)], check=True, capture_output=True, text=True, timeout=30)
