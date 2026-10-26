import re
from pathlib import Path

from ecd_epitope_scanner.exporting import OFFLINE_SITE_ASSETS, build_offline_site


def test_task_site_is_local_complete_and_noncomputing(tmp_path):
    (tmp_path / "structure.cif").write_text("data_fixture\n#\n", encoding="utf-8")
    (tmp_path / "structure_for_energy.pdb").write_text("END\n", encoding="utf-8")
    names = build_offline_site(
        tmp_path,
        {"protein": {"protein_name": "fixture"}, "residues": [], "regions": [], "designs": []},
    )

    page = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert "window.APP_OFFLINE" in page
    assert "connect-src 'self' file: blob: data:" in page
    assert "type=\"module\"" not in page
    assert not re.search(r'(?:href|src)="(?:/|https?:|//)|/api/', page, re.IGNORECASE)
    for attr, url in re.findall(r'\b(href|src)="([^"]+)"', page):
        if url.startswith(("./", "#")):
            if url.startswith("./"):
                assert (tmp_path / url[2:]).is_file(), f"missing local {attr}: {url}"

    bundle = (tmp_path / "assets" / "app.bundle.js").read_text(encoding="utf-8")
    assert "if (state.offline) throw new Error" in bundle
    assert "LOCAL_DRAFT_NOT_ANALYZED" in bundle
    assert set(names) == {"index.html", "assets/app.bundle.js", *(f"assets/{name}" for name in OFFLINE_SITE_ASSETS)}


def test_standard_export_includes_active_refinement_evidence_only(tmp_path):
    import json
    import zipfile
    from ecd_epitope_scanner.exporting import write_scientific_outputs
    from ecd_epitope_scanner.models import ProteinContext
    (tmp_path / 'structure_for_energy.pdb').write_text('END\n', encoding='ascii')
    active = tmp_path / 'foldx_refinement/batch_01/paired_model_manifest.json'
    archived = tmp_path / 'foldx_refinement/prior_attempts/old/paired_model_manifest.json'
    (tmp_path / 'paired_model_manifest.json').write_text('{"records":[]}', encoding='utf-8')
    for path in (active, archived):
        path.parent.mkdir(parents=True)
        path.write_text('{"records":[]}', encoding='utf-8')
    protein = ProteinContext(accession='TEST', protein_name='test', organism='test',
                             sequence='AAAA', ecd_start=1, ecd_end=4,
                             membrane_anchor_position=None, topology='unspecified')
    files = write_scientific_outputs(tmp_path, {}, protein, [], [], [])
    name = active.relative_to(tmp_path).as_posix()
    assert name in files
    assert 'paired_model_manifest.json' in files
    assert archived.relative_to(tmp_path).as_posix() not in files
    manifest = json.loads((tmp_path / 'run_manifest.json').read_text())
    assert name in manifest['files']
    with zipfile.ZipFile(tmp_path / 'ecd_epitope_design.zip') as archive:
        assert archive.read(name) == active.read_bytes()
        assert archive.read('paired_model_manifest.json') == (tmp_path / 'paired_model_manifest.json').read_bytes()


def test_task_site_falls_back_to_saved_pdb_when_cif_is_absent(tmp_path):
    pdb = "HEADER    OFFLINE FIXTURE\nATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    (tmp_path / "structure_for_energy.pdb").write_text(pdb, encoding="ascii")
    build_offline_site(tmp_path, {"protein": {}, "residues": [], "regions": [], "designs": []})
    page = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert '"structure_data":"HEADER    OFFLINE FIXTURE' in page or '"structure_data": "HEADER    OFFLINE FIXTURE' in page
    assert "const format =" in (tmp_path / "assets" / "renderer.js").read_text(encoding="utf-8")


def test_site_and_standalone_report_embed_the_same_final_analysis_projection(tmp_path):
    import ecd_epitope_scanner.exporting as exporting
    from ecd_epitope_scanner.joint_design_v4.reporting import project_display_analysis

    (tmp_path / "structure_for_energy.pdb").write_text(
        "HEADER    OFFLINE FIXTURE\nATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n",
        encoding="ascii",
    )
    analysis = {"protein": {"protein_name": "offline projection fixture"}, "residues": [],
                "regions": [], "designs": []}
    exporting.build_offline_site(tmp_path, analysis)
    projection = exporting._json_inline(project_display_analysis(analysis))
    index = (tmp_path / "index.html").read_text(encoding="utf-8")
    report = exporting.build_offline_html(analysis, tmp_path)

    assert projection in index
    assert f"window.APP_DATA = {projection};" in report
    assert "offline projection fixture" in projection


def test_finalization_refreshes_offline_report_from_final_analysis(tmp_path, monkeypatch):
    import json
    import zipfile

    import ecd_epitope_scanner.exporting as exporting
    from ecd_epitope_scanner.util import sha256_file

    (tmp_path / "structure_for_energy.pdb").write_text("END\n", encoding="ascii")
    (tmp_path / "offline_report.html").write_text("stale", encoding="utf-8")
    analysis = {"files": {f"artifact_{i}.json": "hash" for i in range(81)},
                "run_diagnostics": {"elapsed_wall_seconds": 123.0}}
    (tmp_path / "analysis.json").write_text(json.dumps(analysis), encoding="utf-8")
    (tmp_path / "run_manifest.json").write_text(
        json.dumps({"files": {"analysis.json": "old", "offline_report.html": "old"}}),
        encoding="utf-8",
    )
    seen = {}

    def render(final_analysis, workdir):
        seen["analysis"] = final_analysis
        assert workdir == tmp_path
        return f"final artifacts={len(final_analysis['files'])};elapsed={final_analysis['run_diagnostics']['elapsed_wall_seconds']}"

    monkeypatch.setattr(exporting, "build_offline_html", render)
    exporting.finalize_scientific_bundle(tmp_path)

    report = (tmp_path / "offline_report.html").read_text(encoding="utf-8")
    assert seen["analysis"] == analysis
    assert report == "final artifacts=81;elapsed=123.0"
    manifest = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["files"]["offline_report.html"] == sha256_file(tmp_path / "offline_report.html")
    with zipfile.ZipFile(tmp_path / "ecd_epitope_design.zip") as archive:
        assert archive.read("offline_report.html").decode("utf-8") == report
