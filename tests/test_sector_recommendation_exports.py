import csv
import json
import shutil
import subprocess
import zipfile
from pathlib import Path

import pytest

from ecd_epitope_scanner import exporting
from ecd_epitope_scanner.joint_design_v4.models import JointV4Result, V4Panel
from ecd_epitope_scanner.joint_design_v4.reporting import write_v4_outputs
from ecd_epitope_scanner.models import ProteinContext
from ecd_epitope_scanner.util import sha256_file


ROOT = Path(__file__).resolve().parents[1]


def _recommendations():
    return [
        dict(sector_id="VS001", candidate_id="C1", status="PASSED", formal_eligible=True,
             failed_gates=[], gate_details=[], mutation_notations=["A1W"], design_id="D1",
             selected_in_panel=True, coverage_gain=0.5, coverage_gain_area_a2=50,
             regional_geometry_support=dict(basis="surface_geodesic_neighborhood_proxy",
                 area_basis="eligible_surface_area_a2", sector_area_a2=200, mean_coverage_fraction=0.75,
                 scenarios={f"geodesic_{radius}": dict(covered_area_a2=150, coverage_fraction=0.75,
                     uncovered_positions=list(range(10, 25))) for radius in (9, 12, 15)}),
             geometry_support=dict(basis="surface_geodesic_neighborhood_proxy", area_basis="eligible_surface_area_a2",
                 mean_marginal_area_a2=25, mean_overlap_fraction=0.4,
                 scenarios={f"geodesic_{radius}": dict(union_area_a2=100, best_lower_order_area_a2=75,
                    marginal_area_a2=25, mean_pair_overlap_fraction=0.4, max_pair_overlap_fraction=0.6)
                    for radius in (9, 12, 15)})),
        dict(sector_id="VS002", candidate_id="C2", status="NEAR_THRESHOLD", formal_eligible=False,
             failed_gates=["COMBINATION_REJECT"], mutation_notations=["C2W"], design_id="D2",
             selected_in_panel=False, coverage_gain=0.4, coverage_gain_area_a2=40,
             regional_geometry_support=dict(basis="surface_geodesic_neighborhood_proxy",
                 area_basis="eligible_surface_area_a2", sector_area_a2=0, mean_coverage_fraction=None,
                 scenarios={"geodesic_12": dict(covered_area_a2=0, coverage_fraction=None, uncovered_positions=[2])}),
             geometry_support=dict(basis="surface_geodesic_neighborhood_proxy", area_basis="residue_count_proxy",
                 scenarios={"geodesic_12": dict(marginal_area_a2=3, mean_pair_overlap_fraction=0.5)}),
             gate_details=[dict(gate="SPURS_MULTI", reason="Combination exceeds stability threshold",
                                observed=2.1, threshold=2.0, excess=0.1, evidence_missing=False),
                           dict(gate="FOLDX", reason="Confirmation unavailable", observed=None,
                                threshold=2.0, excess=None, evidence_missing=True)]),
        dict(sector_id="VS003", candidate_id=None, status="NO_COMBINATION", formal_eligible=False,
             failed_gates=["NO_CANDIDATES"], gate_details=[], mutation_notations=[], design_id=None,
             selected_in_panel=False, coverage_gain=0, coverage_gain_area_a2=0),
    ]


def test_sector_gate_evidence_survives_exports_manifest_and_zip(tmp_path):
    panel = V4Panel(["C1"], {}, {}, [], "COVERAGE_FIRST")
    panel.sector_recommendations = _recommendations()
    result = JointV4Result("4.0", "FORMAL_RECOMMENDATION", {}, [], [], [], [], [], panel)
    write_v4_outputs(tmp_path, result, parameters={}, contexts={})
    assert json.loads((tmp_path / "sector_recommendations_v4.json").read_text()) == _recommendations()
    with (tmp_path / "sector_recommendations_v4.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["status"] for row in rows] == ["PASSED", "NEAR_THRESHOLD", "NO_COMBINATION"]
    assert json.loads(rows[1]["gate_details"]) == _recommendations()[1]["gate_details"]
    assert json.loads(rows[0]["geometry_support"]) == _recommendations()[0]["geometry_support"]
    assert json.loads(rows[0]["regional_geometry_support"]) == _recommendations()[0]["regional_geometry_support"]
    assert json.loads(rows[1]["regional_geometry_support"])["mean_coverage_fraction"] is None
    report = (tmp_path / "v4_scientific_report.md").read_text(encoding="utf-8")
    assert "Confirmation unavailable" in report and "NEAR_THRESHOLD" in report
    assert "geodesic_15" in report and "not measured perturbation" in report
    protein = ProteinContext("TEST", "test", "test", "AC", 1, 2, 2, "single_pass")
    payload = {"joint_v4_ml": {"selected_panel": {"sector_recommendations": _recommendations()}}}
    (tmp_path / "structure_for_energy.pdb").write_text("END\n", encoding="ascii")
    files = exporting.write_scientific_outputs(tmp_path, payload, protein, [], [], [])
    manifest = json.loads((tmp_path / "run_manifest.json").read_text())
    with zipfile.ZipFile(tmp_path / "ecd_epitope_design.zip") as archive:
        for name in ("sector_recommendations_v4.json", "sector_recommendations_v4.csv"):
            assert name in files
            assert manifest["files"][name] == sha256_file(tmp_path / name)
            assert archive.read(name) == (tmp_path / name).read_bytes()


def test_sector_ui_uses_explicit_design_and_shows_gate_values(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required to execute UI rendering")
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    renderer = (ROOT / "web/recommendations.js").read_text(encoding="utf-8").replace("export function", "function", 1)
    code = renderer + app[app.index("const { sectorRecommendationText,"):app.index("function renderStep4()")]
    code += app[app.index("function renderStep5()"):app.index("function renderStep6()")]
    setup = """
const assert = require('assert');
const elements = {'region-list':{previousElementSibling:null}, 'design-list':{}};
const $ = id => elements[id];
const h = value => String(value ?? '');
const t = (key, fallback) => fallback;
const notation = item => item.wt + item.position + item.mutant;
const fmt = value => String(value);
const fmtPct = value => String(value * 100) + '%';
const document = {querySelectorAll:() => []};
const scoreInfo = () => ({real:false});
const renderExports = () => {};
const renderDesignSelect = () => {};
const regionsForDesign = design => [design.region_id];
"""
    setup += "const recommendations = " + json.dumps(_recommendations()) + ";\n"
    setup += """
const state = {lang:'zh', result:{residues:[], regions:[1,2,3].map(id=>({id,name:'Region '+id,residues:[]})),
  designs:[{id:'D1',region_id:1,design_role:'STAGE1_PROBE',mutations:[{wt:'A',position:1,mutant:'W'}]},
           {id:'CONTROL',region_id:1,design_role:'STAGE2_SINGLE',mutations:[{wt:'A',position:1,mutant:'E'}]},
           {id:'D2',region_id:1,design_role:'SECTOR_REVIEW_CANDIDATE',mutations:[{wt:'C',position:2,mutant:'W'}]}],
  joint_v4_ml:{selected_panel:{sector_recommendations:recommendations}, sectors:recommendations.map(r=>({sector_id:r.sector_id}))}}};
"""
    checks = """
renderStep3();
assert(elements['region-list'].innerHTML.includes('C2W'));
assert(!elements['region-list'].innerHTML.includes('A1E'));
renderStep5();
const output = elements['design-list'].innerHTML;
assert(output.includes('2.1'));
assert(!output.includes('geodesic_'));
assert(output.includes('9 Å') && output.includes('12 Å') && output.includes('15 Å'));
assert(output.includes('25 Å²'));
assert(!output.includes('3 Å²'));
assert(output.includes('40%'));
assert(output.includes('几何估计'));
assert(output.includes('75%'));
assert(output.includes('未计算／缺少证据'));
assert(output.includes('10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24'));
assert(output.includes('门槛'));
assert(output.includes('缺少证据'));
assert(output.includes('尚未解释的审查项'));
assert(output.includes('没有安全的专属组合'));
assert(output.includes('查看详细证据'));
assert(output.includes('候选追踪编号'));
assert(output.includes('data-review-design="D2"'));
assert(output.includes('data-design-residue="2"'));
delete state.result.joint_v4_ml.selected_panel.sector_recommendations;
renderStep3();
assert(elements['region-list'].innerHTML.includes('A1W'));
assert(!elements['region-list'].innerHTML.includes('A1E'));
"""
    script = tmp_path / "sector_ui.cjs"
    script.write_text(setup + code + checks, encoding="utf-8")
    subprocess.run([node, str(script)], check=True, capture_output=True, text=True)


def test_offline_sector_review_rendering_includes_failed_gate_details(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required to execute offline rendering")
    template = (ROOT / "web/offline_template.html").read_text(encoding="utf-8")
    renderer = (ROOT / "web/recommendations.js").read_text(encoding="utf-8").replace("export function", "function", 1)
    code = renderer + template[template.index("function stepDesign()"):template.index("function stepCustom()")]
    setup = "const D=" + json.dumps({"joint_v4_ml": {"selected_panel": {"sector_recommendations": _recommendations()}}}) + ";"
    setup += "D.designs=[]; D.regions=[]; D.residues=[]; D.joint_v4_ml.sectors=D.joint_v4_ml.selected_panel.sector_recommendations; const LANG='en';const designs=[]; const state={}; const esc=value=>String(value??''); const t=key=>key; const chip=()=>'';"
    checks = "const assert=require('assert');const output=stepDesign();assert(output.includes('2.100'));assert(output.includes('Missing evidence'));assert(output.includes('Unmapped review item'));assert(output.includes('data-review-design=\"D2\"'));assert(output.includes('No safe exclusive combination'));assert(!output.includes('geodesic_'));assert(output.includes('25.000 Å²'));assert(!output.includes('3.000 Å²'));assert(output.includes('geometric estimate'));"
    checks += "assert(output.includes('75%'));assert(output.includes('Not computed / missing evidence'));assert(output.includes('10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24'));assert(output.includes('View detailed evidence'));"
    script = tmp_path / "sector_offline.cjs"
    script.write_text(setup + code + checks, encoding="utf-8")
    subprocess.run([node, str(script)], check=True, capture_output=True, text=True)
