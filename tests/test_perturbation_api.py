import pytest
from pydantic import ValidationError

from ecd_epitope_scanner.api import JobCreate
from ecd_epitope_scanner.pipeline import PipelineRequest


def test_explicit_protocol_and_construct_budget_are_preserved():
    payload = JobCreate(identifier="P16410", execution_protocol="ordinary_site_first_perturbation_v1",
                        maximum_constructs=12)
    request = PipelineRequest(**payload.model_dump())
    assert request.metadata["joint_v4_ml"]["execution_protocol"]["name"] == payload.execution_protocol
    assert request.maximum_constructs == 12
    assert PipelineRequest(input_type="uniprot", identifier="P16410").execution_protocol is None
    with pytest.raises(ValidationError):
        JobCreate(identifier="P16410", execution_protocol="unknown")


def test_uploaded_structure_preserves_selected_protocol(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from ecd_epitope_scanner import api
    from ecd_epitope_scanner.jobs import JobStore

    store = JobStore(tmp_path)
    monkeypatch.setattr(api, "STORE", store)
    monkeypatch.setattr(api.RUNNER, "submit", lambda job_id: None)
    response = TestClient(api.app).post(
        "/api/jobs/upload",
        data={"execution_protocol": "ordinary_site_first_perturbation_v1", "maximum_constructs": "12"},
        files={"structure": ("structure.pdb", b"END\n", "chemical/x-pdb")},
    )
    assert response.status_code == 202
    request = store.request_for(response.json()["id"])
    assert request.execution_protocol == "ordinary_site_first_perturbation_v1"
    assert request.maximum_constructs == 12
