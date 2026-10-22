from __future__ import annotations

import io

from fastapi.testclient import TestClient

from ecd_epitope_scanner.api import app


client = TestClient(app)


def _isolate_upload(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api_module
    from ecd_epitope_scanner.jobs import JobStore

    monkeypatch.setattr(api_module, "STORE", JobStore(tmp_path / "workspace"))
    monkeypatch.setattr(api_module.RUNNER, "submit", lambda job_id: None)


def _upload(data=None):
    return client.post(
        "/api/jobs/upload",
        data=data or {},
        files={"structure": ("binder.pdb", io.BytesIO(b"HEADER binder\nEND\n"), "chemical/x-pdb")},
    )


def test_binder_upload_persists_target_and_skips_antigen_assay_defaults(monkeypatch, tmp_path):
    _isolate_upload(monkeypatch, tmp_path)

    response = _upload({"analysis_target": "binder", "chain": "B"})

    assert response.status_code == 202, response.text
    request = response.json()["request"]
    assert request["analysis_target"] == "binder"
    assert request["assay_profile"] is None
    assert request["assay_format"] == "unspecified"
    assert request["chain"] == "B"


def test_structure_upload_defaults_to_antigen_mode(monkeypatch, tmp_path):
    _isolate_upload(monkeypatch, tmp_path)
    from ecd_epitope_scanner.assay_profiles import DEFAULT_ASSAY_PROFILE_ID

    response = _upload()

    assert response.status_code == 202, response.text
    request = response.json()["request"]
    assert request["analysis_target"] == "antigen"
    assert request["assay_profile"] == DEFAULT_ASSAY_PROFILE_ID
    assert request["assay_format"] == "cell_surface"


def test_structure_upload_rejects_unknown_analysis_target():
    response = _upload({"analysis_target": "receptor"})

    assert response.status_code == 422
