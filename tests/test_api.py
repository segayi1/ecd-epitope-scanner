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
