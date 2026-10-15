import copy
import io

import pytest
from fastapi.testclient import TestClient

from ecd_epitope_scanner import api
from ecd_epitope_scanner.jobs import JobStore
from ecd_epitope_scanner.pipeline import PipelineRequest


def test_display_api_is_optional_and_preserves_full_result(monkeypatch, tmp_path):
    store = JobStore(tmp_path / "jobs")
    monkeypatch.setattr(api, "STORE", store)
    job = store.create(PipelineRequest(input_type="uniprot", identifier="Q9NZQ7"))
    analysis = {"protein": {"accession": "Q9NZQ7"}, "regions": [], "residues": [],
        "joint_v4_ml": {"candidates": [{"candidate_id": "C1"}],
            "selected_panel": {"selected_candidate_ids": [], "sector_recommendations": []}},
        "partitioning": {"selected_panel": {"selected_candidate_ids": []}}}
    store.update(job["id"], status="COMPLETE", stage="COMPLETE", result_json=analysis)
    original = copy.deepcopy(store.get(job["id"]))
    client = TestClient(api.app)
    full = client.get(f"/api/jobs/{job['id']}").json()["result"]
    shown = client.get(f"/api/jobs/{job['id']}?view=display").json()["result"]
    assert full["joint_v4_ml"]["candidates"] == analysis["joint_v4_ml"]["candidates"]
    assert "candidates" not in shown["joint_v4_ml"]
    assert shown["joint_v4_ml"]["candidates_reference"]["file"] == "panel_candidates.json"
    assert shown["structure_url"] == full["structure_url"]
    assert store.get(job["id"]) == original
    assert client.get(f"/api/jobs/{job['id']}?view=bad").status_code == 422


@pytest.mark.parametrize("field,value", [
    ("target_sector_residues", 0), ("target_sector_residues", True),
    ("target_sector_residues", 9), ("target_sector_residues", 16),
    ("maximum_sector_geodesic_diameter_a", float("nan")),
    ("maximum_sector_ca_diameter_a", float("inf")),
    ("maximum_sector_ca_diameter_a", -1),
])
def test_display_request_rejects_invalid_settings(field, value):
    with pytest.raises(ValueError):
        PipelineRequest(input_type="uniprot", identifier="Q9NZQ7", **{field: value})


def test_json_and_upload_propagate_display_settings(monkeypatch, tmp_path):
    monkeypatch.setattr(api, "STORE", JobStore(tmp_path / "jobs"))
    monkeypatch.setattr(api.RUNNER, "submit", lambda job_id: None)
    values = {"target_sector_residues": 12, "maximum_sector_geodesic_diameter_a": 30,
              "maximum_sector_ca_diameter_a": 22}
    client = TestClient(api.app)
    responses = [client.post("/api/jobs", json={"identifier": "Q9NZQ7", **values}),
        client.post("/api/jobs/upload", data=values,
            files={"structure": ("input.pdb", io.BytesIO(b"HEADER test\nEND\n"), "chemical/x-pdb")})]
    for response in responses:
        assert response.status_code == 202, response.text
        request = api.STORE.request_for(response.json()["id"])
        settings = request.metadata["joint_v4_ml"]["sectoring"]
        assert settings["target_sector_residues"] == 12
        assert settings["maximum_geodesic_diameter_a"] == 30
        assert settings["maximum_euclidean_diameter_a"] == 22
        assert api.STORE.request_for(response.json()["id"]).as_dict() == request.as_dict()
