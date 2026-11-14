import threading

from fastapi.testclient import TestClient

from ecd_epitope_scanner import api, jobs
from ecd_epitope_scanner.pipeline import PipelineRequest


def test_cancel_running_and_queued_jobs(tmp_path, monkeypatch):
    store = jobs.JobStore(tmp_path)
    runner = jobs.JobRunner(store)
    monkeypatch.setattr(api, "STORE", store)
    monkeypatch.setattr(api, "RUNNER", runner)
    started, release = threading.Event(), threading.Event()
    called = []

    def pipeline(request, progress):
        called.append(request.identifier)
        started.set()
        assert release.wait(10)
        progress("TEST", .5, "working")
        return {}

    monkeypatch.setattr(jobs, "run_pipeline", pipeline)
    first = store.create(PipelineRequest(input_type="uniprot", identifier="P16410"))
    second = store.create(PipelineRequest(input_type="uniprot", identifier="Q9NZQ7"))
    runner.submit(first["id"])
    assert started.wait(10)
    runner.submit(second["id"])
    client = TestClient(api.app)
    assert client.post(f"/api/jobs/{first['id']}/cancel").json()["status"] == "CANCELLING"
    assert client.post(f"/api/jobs/{second['id']}/cancel").json()["status"] == "CANCELLED"
    release.set()
    runner.executor.shutdown(wait=True)
    assert store.get(first["id"])["status"] == "CANCELLED"
    assert store.get(second["id"])["terminal"] is True
    assert called == ["P16410"]
    assert client.post(f"/api/jobs/{first['id']}/cancel").status_code == 200
    store.update(first["id"], status="COMPLETE")
    assert client.post(f"/api/jobs/{first['id']}/cancel").status_code == 409
    waiting = store.create(PipelineRequest(input_type="uniprot", identifier="Q15116"))
    store.update(waiting["id"], status="WAITING_FOR_AF3_UPLOAD")
    assert client.post(f"/api/jobs/{waiting['id']}/cancel").json()["status"] == "CANCELLED"


def test_completed_job_stores_projection_instead_of_large_candidate_records(tmp_path, monkeypatch):
    store = jobs.JobStore(tmp_path)
    runner = jobs.JobRunner(store)
    job = store.create(PipelineRequest(input_type="uniprot", identifier="Q9NZQ7"))
    result = {"joint_v4_ml": {"candidates": [{"candidate_id": "C1", "details": "x" * 100000}],
                               "selected_panel": {}}, "design_summary": {}}
    monkeypatch.setattr(jobs, "run_pipeline", lambda request, progress: result)
    runner._run(job["id"])
    saved = store.get(job["id"])
    assert saved["status"] == "COMPLETE"
    assert "candidates" not in saved["result"]["joint_v4_ml"]
    assert saved["result"]["joint_v4_ml"]["candidates_reference"]["ids"] == ["C1"]
    assert "candidates" in result["joint_v4_ml"]
    runner.executor.shutdown()
