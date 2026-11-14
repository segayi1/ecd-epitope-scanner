import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_mutant_combination_uses_surface_showcase_contract():
    app = (ROOT / "web" / "app.js").read_text(encoding="utf-8")
    index = (ROOT / "web" / "index.html").read_text(encoding="utf-8")
    css = (ROOT / "web" / "app.css").read_text(encoding="utf-8")
    renderer = (ROOT / "web" / "renderer.js").read_text(encoding="utf-8")

    assert 'id="design-surface-callout"' in index
    assert ".design-surface-callout" in css
    assert 'await setRepresentation("molecular-surface")' in app
    assert 'representation: "molecular-surface"' in app
    assert 'data-representation="molecular-surface" aria-pressed="true"' in index
    assert 'new Set(["surface", "regions", "design"])' in app
    assert "function regionForDesign" in app
    assert "owner.parent_design_id" in app
    assert "function designSurfaceRecords" in app
    assert 'new URLSearchParams(window.location.search).get("job")' in app
    assert "add(records.envelope, DESIGN_SURFACE_COLORS.region)" in app
    assert "add(records.footprint, DESIGN_SURFACE_COLORS.footprint)" in app
    assert "add(records.mutations, DESIGN_SURFACE_COLORS.mutation);" in app
    assert "representationColor: DESIGN_SURFACE_COLORS.mutation" not in app
    assert "WT 位点" in app
    assert "instance.applyVisualStyles(instance.assemblyRef, visualStyle)" in renderer
    assert "visual.update(\n      { visualStyle" not in renderer
    assert '...panelBudgetInput()' in app
    assert "max_sectors: 12" in app
    assert 'form.append("max_sectors", "12")' in app
    assert 'const ceiling = 12;' in app
    assert 'class="design-series"' in app


def test_surface_showcase_translations_are_complete():
    keys = {
        "design.surfaceTitle",
        "design.surfaceExit",
        "design.surfaceCoordinateNote",
        "design.surfaceProtein",
        "design.surfaceRegion",
        "design.surfaceSites",
        "design.epitopePrioritized",
        "design.stabilityFallback",
        "design.combinations",
    }
    for language in ("zh", "en"):
        payload = json.loads(
            (ROOT / "web" / "i18n" / f"{language}.json").read_text(encoding="utf-8")
        )
        assert keys <= payload.keys()
        assert all(payload[key].strip() for key in keys)


@pytest.mark.parametrize(
    "page",
    (
        "index.html",
        "experiments.html",
        "surface_partitions.html",
        "ctla4_six_mutants.html",
        "pd1_antibodies.html",
    ),
)
def test_job_center_is_available_from_every_workbench_page(page):
    html = (ROOT / "web" / page).read_text(encoding="utf-8")
    assert '/assets/job_center.js?v=global-jobs-20261111-v1' in html
