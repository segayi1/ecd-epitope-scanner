from __future__ import annotations

import io
import json
import time

from fastapi.testclient import TestClient

from ecd_epitope_scanner.api import app


client = TestClient(app)


def test_health_and_workbench_assets_are_served():
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.json()["status"] == "ok"
    assert "ECD Epitope Workbench" in client.get("/").text
    assert client.get("/assets/demo/analysis.json").status_code == 200
    assert client.get("/assets/vendor/pdbe-molstar-plugin-3.10.1.js").status_code == 200


def test_job_list_is_a_small_summary_and_owner_scoped(tmp_path):
    from ecd_epitope_scanner.jobs import JobStore
    from ecd_epitope_scanner.pipeline import PipelineRequest

    store = JobStore(tmp_path / "workspace")
    job = store.create(PipelineRequest(input_type="uniprot", identifier="Q9NZQ7"))
    assert job["owner_id"] == "local"
    assert store.list(owner_id="another-user") == []
    # The HTTP list endpoint must not serialize full per-residue analyses.
    response = client.get("/api/jobs")
    assert response.status_code == 200
    assert all("result" not in item and "request" not in item for item in response.json())


def test_af3_resume_backfills_accession_from_identifier(monkeypatch, tmp_path):
    """P0 regression: a UniProt job must survive the AF3-upload resume.

    The API form stores the UniProt accession in `identifier` and leaves
    `accession` empty.  resume_with_af3 used to switch input_type to
    af3_upload without backfilling, so the pipeline fell into
    _synthetic_context and every resume failed with
    "requires sequence, ecd_start and ecd_end".
    """
    import ecd_epitope_scanner.api as api_module  # noqa: F401  (import sanity)
    from ecd_epitope_scanner.pipeline import PipelineRequest
    from ecd_epitope_scanner.util import WaitingForAF3

    seen_requests: list[PipelineRequest] = []

    def fake_run(request, progress_callback=None):
        if request.input_type == "uniprot":
            # First run: no eligible structure -> the WAITING_FOR_AF3_UPLOAD state.
            raise WaitingForAF3("No experimental structure passed the gate.", [])
        seen_requests.append(request)
        progress = progress_callback or (lambda *args: None)
        progress("COMPLETE", 1.0, "fake complete")
        return {"ok": True}

    monkeypatch.setattr("ecd_epitope_scanner.jobs.run_pipeline", fake_run)

    response = client.post(
        "/api/jobs",
        json={"input_type": "uniprot", "identifier": "Q9NZQ7"},
    )
    assert response.status_code == 202
    job_id = response.json()["id"]
    final = _wait_for_status(client, job_id, {"COMPLETE", "FAILED", "WAITING_FOR_AF3_UPLOAD"})
    assert final["status"] == "WAITING_FOR_AF3_UPLOAD"

    files = {
        "structure": ("model.cif", io.BytesIO(b"data_ cif placeholder\n"), "chemical/x-mmcif"),
        "confidence": (
            "confidences.json",
            io.BytesIO(json.dumps({"pae": [], "atom_plddts": []}).encode("utf-8")),
            "application/json",
        ),
    }
    data = {"official_source_confirmed": "true"}
    response = client.post(f"/api/jobs/{job_id}/af3", files=files, data=data)
    assert response.status_code == 202, response.text
    final = _wait_for_status(client, job_id, {"COMPLETE", "FAILED"})
    assert final["status"] == "COMPLETE", final.get("error")

    assert seen_requests, "resume must re-run the pipeline"
    assert seen_requests[-1].input_type == "af3_upload"
    assert seen_requests[-1].accession == "Q9NZQ7", (
        "resume must backfill accession from the UniProt identifier"
    )


def _wait_for_status(client: TestClient, job_id: str, statuses: set[str], timeout: float = 20.0) -> dict:
    deadline = time.time() + timeout
    last: dict = {}
    while time.time() < deadline:
        last = client.get(f"/api/jobs/{job_id}").json()
        if last["status"] in statuses:
            return last
        time.sleep(0.05)
    raise AssertionError(f"job {job_id} never reached {statuses}; last={last['status']}")
