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
    ("target_sector_residues", 6), ("target_sector_residues", 16),
    ("maximum_sector_geodesic_diameter_a", float("nan")),
    ("maximum_sector_ca_diameter_a", float("inf")),
    ("maximum_sector_ca_diameter_a", -1),
])
def test_display_request_rejects_invalid_settings(field, value):
    with pytest.raises(ValueError):
        PipelineRequest(input_type="uniprot", identifier="Q9NZQ7", **{field: value})


def test_display_request_defaults_to_latest_4_to_12_policy():
    from ecd_epitope_scanner.cli import build_parser

    request = PipelineRequest(input_type="uniprot", identifier="Q9NZQ7")
    job = api.JobCreate(identifier="Q9NZQ7")
    args = build_parser().parse_args(
        ["analyze", "--uniprot", "Q9NZQ7", "--output", "out"]
    )

    assert (request.minimum_regions, request.maximum_regions, request.max_sectors) == (4, 12, 12)
    assert job.max_sectors == 12
    assert (args.min_regions, args.max_regions, args.max_sectors) == (4, 12, 12)


def test_json_and_upload_accept_twelve_construct_budget(monkeypatch, tmp_path):
    monkeypatch.setattr(api, "STORE", JobStore(tmp_path / "jobs"))
    monkeypatch.setattr(api.RUNNER, "submit", lambda job_id: None)
    client = TestClient(api.app)
    values = {"minimum_constructs": 4, "target_constructs": 12, "maximum_constructs": 12}
    responses = [
        client.post("/api/jobs", json={"identifier": "P16410", **values}),
        client.post("/api/jobs/upload", data=values,
                    files={"structure": ("input.pdb", io.BytesIO(b"HEADER test\nEND\n"), "chemical/x-pdb")}),
    ]
    for response in responses:
        assert response.status_code == 202, response.text
        request = api.STORE.request_for(response.json()["id"])
        assert request.maximum_constructs == 12


def test_json_and_upload_propagate_explicit_protocol(monkeypatch, tmp_path):
    monkeypatch.setattr(api, "STORE", JobStore(tmp_path / "jobs"))
    monkeypatch.setattr(api.RUNNER, "submit", lambda job_id: None)
    client = TestClient(api.app)
    protocol = "ordinary_site_first_perturbation_v1"
    responses = [client.post("/api/jobs", json={"identifier": "P16410", "execution_protocol": protocol}),
        client.post("/api/jobs/upload", data={"execution_protocol": protocol},
            files={"structure": ("input.pdb", io.BytesIO(b"HEADER test\nEND\n"), "chemical/x-pdb")})]
    for response in responses:
        assert response.status_code == 202, response.text
        request = api.STORE.request_for(response.json()["id"])
        assert request.metadata["joint_v4_ml"]["execution_protocol"]["name"] == protocol
    assert client.post("/api/jobs", json={"identifier": "P16410", "ml_mode": "preview",
                                       "execution_protocol": protocol}).status_code == 422
    rejected = client.post("/api/jobs/upload", data={"execution_protocol": protocol, "ml_mode": "preview"},
        files={"structure": ("input.pdb", io.BytesIO(b"HEADER test\nEND\n"), "chemical/x-pdb")})
    assert rejected.status_code == 422


@pytest.mark.parametrize("target", [7, 8, 12])
def test_json_and_upload_propagate_display_settings(monkeypatch, tmp_path, target):
    monkeypatch.setattr(api, "STORE", JobStore(tmp_path / "jobs"))
    monkeypatch.setattr(api.RUNNER, "submit", lambda job_id: None)
    values = {"target_sector_residues": target, "maximum_sector_geodesic_diameter_a": 30,
              "maximum_sector_ca_diameter_a": 22}
    client = TestClient(api.app)
    responses = [client.post("/api/jobs", json={"identifier": "Q9NZQ7", **values}),
        client.post("/api/jobs/upload", data=values,
            files={"structure": ("input.pdb", io.BytesIO(b"HEADER test\nEND\n"), "chemical/x-pdb")})]
    for response in responses:
        assert response.status_code == 202, response.text
        request = api.STORE.request_for(response.json()["id"])
        settings = request.metadata["joint_v4_ml"]["sectoring"]
        assert settings["target_sector_residues"] == target
        assert settings["maximum_geodesic_diameter_a"] == 30
        assert settings["maximum_euclidean_diameter_a"] == 22
        assert api.STORE.request_for(response.json()["id"]).as_dict() == request.as_dict()
