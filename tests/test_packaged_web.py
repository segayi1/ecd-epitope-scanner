import importlib

import pytest
from pydantic import ValidationError


def test_health_and_packaged_web_assets(monkeypatch, tmp_path):
    monkeypatch.setenv("ECD_DATA_DIR", str(tmp_path / "data"))
    api = importlib.import_module("ecd_epitope_scanner.api")
    api = importlib.reload(api)

    assert api.health()["status"] == "ok"
    from ecd_epitope_scanner import __version__
    from ecd_epitope_scanner.pipeline import _software_versions

    assert api.health()["version"] == api.app.version == __version__
    assert _software_versions()["ecd_epitope_workbench"].endswith(" " + __version__)
    index = (api.WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "PDBe Mol*" in index or "molstar" in index.lower()
    assert 'href="/experiments"' not in index
    assert 'surface_perturbation_v1_preview' not in index
    assert (api.WEB_ROOT / "vendor" / "PDBe-Molstar-LICENSE").is_file()
    assert any(getattr(route, "path", "") == "/assets" for route in api.app.routes)


def test_release_exposes_only_site_first_v2(monkeypatch, tmp_path):
    monkeypatch.setenv("ECD_DATA_DIR", str(tmp_path / "data"))
    api = importlib.import_module("ecd_epitope_scanner.api")
    api = importlib.reload(api)

    for asset in ("index.html", "app.js", "offline_template.html"):
        assert "surface_perturbation_v1_preview" not in (api.WEB_ROOT / asset).read_text(encoding="utf-8")
    with pytest.raises(ValidationError):
        api.JobCreate(identifier="Q9NZQ7", surface_strategy="surface_perturbation_v1_preview")
    with pytest.raises(ValidationError):
        api.JobCreate(identifier="Q9NZQ7", optimization_profile="coverage_scan_v2_preview")
