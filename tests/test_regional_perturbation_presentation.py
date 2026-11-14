import csv
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from ecd_epitope_scanner.joint_design_v4.models import (
    JointV4Result,
    ResidueKey,
    SingleMutationEvidence,
    SurfaceEvidence,
    V4Panel,
)
from ecd_epitope_scanner.joint_design_v4.reporting import write_v4_outputs


ROOT = Path(__file__).resolve().parents[1]
POLICY = "regional_mechanism_pareto_v1"


def _result():
    residue = ResidueKey("1", "A", 7, uniprot_pos=7, one_letter="K")
    exposure = {"value": 0.42, "numerator_a2": 12.6, "denominator_a2": 30.0,
                "status": "OK", "method": "sidechain_fraction_v1"}
    axes = {"charge_delta": 1, "volume_delta_a3": 14.2, "hydropathy_delta": 2.1,
            "donor_delta": 1, "acceptor_delta": 0, "aromatic_delta": 1}
    item = {
        "sector_id": "VS001", "candidate_id": "C1", "status": "PASSED",
        "formal_eligible": True, "recommendation_policy": POLICY,
        "principles": ["mechanism_pareto"], "mutation_notations": ["K7W"],
        "selected_in_panel": False, "design_id": None,
        "impact_evidence": {"policy": POLICY, "axes": axes, "exposure_min": 0.42,
                             "exposure_mean": 0.42, "epitope_prior": 0.8,
                             "ranking_version": "mechanism_channels_v2",
                             "signed_change_summary": {"charge_delta": {"net_change": 0, "absolute_sum": 4,
                                                                          "opposition_fraction": 1}}},
        "resource_audit": {"single_threshold_audit": {"review_only": True, "review_band_count": 4}},
        "sensitivity": {"sensitive": True, "selection_winners": {"ordinal": "C1", "minmax": "C2",
                       "mean_exposure": "C1", "leave_one_out:exposure": "C1"},
                        "variants": [{"name": "ordinal", "winner_candidate_id": "C1"},
                                     {"name": "minmax", "winner_candidate_id": "C2"},
                                     {"name": "mean_exposure", "winner_candidate_id": "C1"}]},
        "direct_area_a2": 120.0, "local_12_covered_area_a2": 88.0,
        "local_12_coverage_fraction": 0.75, "local_12_uncovered_area_a2": 30.0,
        "local_12_uncovered_positions": [10, 11],
        "site_details": [{"notation": "K7W", "position": 7, "sidechain_exposure": exposure,
                          "chemical_descriptors": {"charge_delta": -1, "volume_delta_a3": 14.2,
                                                     "hydropathy_delta": 2.1, "donor_delta": -1,
                                                     "acceptor_delta": 0, "aromatic_delta": 1},
                          "epitope_prior": 0.8, "epitope_supported": True}],
        "regional_geometry_support": {"scenarios": {
            f"geodesic_{radius}": {"coverage_fraction": 0.7, "covered_area_a2": 88.0,
                                    "uncovered_positions": [10, 11]}
            for radius in (9, 12, 15)}},
        "geometry_support": {"area_basis": "eligible_surface_area_a2", "scenarios": {
            f"geodesic_{radius}": {"marginal_area_a2": 20.0, "mean_pair_overlap_fraction": 0.2,
                                    "max_pair_overlap_fraction": 0.3}
            for radius in (9, 12, 15)}},
    }
    panel = V4Panel(["C1"], {}, {}, [], "TEST", sector_recommendations=[item])
    option = SingleMutationEvidence(
        residue=residue, mutant_aa="W", stability_band="LOW", chemistry_class="aromatic",
        eligible=True, chemical_descriptors={"charge_delta": -1}, impact_evidence={
            "policy": POLICY, "sidechain_exposure": exposure, "axes": axes},
    )
    return JointV4Result("4.1", "FORMAL_RECOMMENDATION", {},
                         [SurfaceEvidence(residue=residue)], [], [option], [], [], panel)


def test_regional_policy_evidence_exports_as_json_csv_and_report(tmp_path):
    result = _result()
    write_v4_outputs(tmp_path, result, parameters={}, contexts={})

    rows = json.loads((tmp_path / "sector_recommendations_v4.json").read_text(encoding="utf-8"))
    assert rows[0]["impact_evidence"]["axes"]["charge_delta"] == 1
    assert rows[0]["site_details"][0]["chemical_descriptors"]["charge_delta"] == -1
    assert rows[0]["local_12_uncovered_positions"] == [10, 11]
    with (tmp_path / "sector_recommendations_v4.csv").open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["recommendation_policy"] == POLICY
    assert json.loads(row["impact_evidence"])["axes"]["volume_delta_a3"] == 14.2
    assert json.loads(row["sensitivity"])["sensitive"] is True
    assert json.loads(row["resource_audit"])["single_threshold_audit"]["review_band_count"] == 4
    assert json.loads(row["impact_evidence"])["signed_change_summary"]["charge_delta"]["net_change"] == 0
    assert float(row["direct_area_a2"]) == 120.0
    assert json.loads(row["local_12_uncovered_positions"]) == [10, 11]

    with (tmp_path / "single_mutation_options.csv").open(encoding="utf-8", newline="") as handle:
        option = next(csv.DictReader(handle))
    assert json.loads(option["chemical_descriptors"])["charge_delta"] == -1
    assert json.loads(option["impact_evidence"])["sidechain_exposure"]["value"] == 0.42
    with (tmp_path / "surface_evidence.csv").open(encoding="utf-8", newline="") as handle:
        surface = next(csv.DictReader(handle))
    assert float(surface["sidechain_exposure_value"]) == 0.42
    assert float(surface["sidechain_exposure_numerator_a2"]) == 12.6
    report = (tmp_path / "v4_scientific_report.md").read_text(encoding="utf-8")
    assert "engineering representative for testing" in report
    assert "do not predict the maximum antibody-affinity loss" in report
    assert "not an independent mutation-effect axis" in report
    assert "mechanism_channels_v2" in report and "zero net change" in report


def test_online_and_offline_renderers_show_mechanism_without_changing_legacy(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required to execute recommendation rendering")
    renderer = (ROOT / "web/recommendations.js").read_text(encoding="utf-8").replace("export function", "function", 1)
    template = (ROOT / "web/offline_template.html").read_text(encoding="utf-8")
    online_app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    offline_step = template[template.index("function stepDesign(){"):template.index("function stepCustom(){")]
    setup = r"""
const assert=require('node:assert/strict');
const h=value=>String(value??'').replace(/[&<>\"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;',"'":'&#39;'}[c]));
const fmt=(value,d=3)=>Number(value).toFixed(d),fmtPct=value=>Math.round(value*100)+'%';
const notation=m=>m.wt+m.position+m.mutant;
const item={sector_id:'VS001',candidate_id:'C1',status:'EVIDENCE_LIMITED',formal_eligible:true,selected_in_panel:true,
 recommendation_policy:'regional_mechanism_pareto_v1',principles:['mechanism_pareto'],mutation_notations:['K7W'],
 impact_evidence:{policy:'regional_mechanism_pareto_v1',axes:{charge_delta:1,volume_delta_a3:14.2,hydropathy_delta:2.1,donor_delta:1,acceptor_delta:0,aromatic_delta:null},missing_axes:['aromatic_delta'],exposure_min:0.42},
sensitivity:{sensitive:true,selection_winners:{ordinal:'C1',minmax:'C2',mean_exposure:'C1','leave_one_out:exposure':'C1'},variants:[{name:'ordinal',winner_candidate_id:'C1'},{name:'minmax',winner_candidate_id:'C2'},{name:'mean_exposure',winner_candidate_id:'C1'}]},direct_area_a2:120,local_12_covered_area_a2:88,
 stability_evidence:{spurs_multi_ddg_kcal_mol:0.1,foldx_ddg_kcal_mol:0.2},
 gate_checks:[{gate:'combination_spurs',passed:true},{gate:'combination_foldx',passed:true}],
 local_12_coverage_fraction:0.75,local_12_uncovered_area_a2:30,local_12_uncovered_positions:[10,11],
site_details:[{notation:'K7W',position:7,risk_assessment:{conclusion:'存在需复核提示',notes:[{code:'salt-bridge',explanation:'尚未解释的保护或风险项，请查看原始证据。'}]},sidechain_exposure:{value:0.42,numerator_a2:12.6,denominator_a2:30,status:'OK',method:'sidechain_fraction_v1'},chemical_descriptors:{charge_delta:-1,volume_delta_a3:14.2,hydropathy_delta:2.1,donor_delta:-1,acceptor_delta:0,aromatic_delta:1},epitope_prior:0.8}],
 regional_geometry_support:{scenarios:{geodesic_9:{coverage_fraction:0.6,uncovered_positions:[10]},geodesic_12:{coverage_fraction:0.75,uncovered_positions:[10,11]},geodesic_15:{coverage_fraction:0.9,uncovered_positions:[]}}},
 geometry_support:{area_basis:'eligible_surface_area_a2',scenarios:{geodesic_9:{marginal_area_a2:20,mean_pair_overlap_fraction:0.2,max_pair_overlap_fraction:0.3},geodesic_12:{marginal_area_a2:22,mean_pair_overlap_fraction:0.2,max_pair_overlap_fraction:0.3},geodesic_15:{marginal_area_a2:24,mean_pair_overlap_fraction:0.2,max_pair_overlap_fraction:0.3}}}};
Object.assign(item.impact_evidence,{ranking_version:'mechanism_channels_v2',aromatic_mutant_count:1,introduced_aromatic_count:1,
 signed_change_summary:{charge_delta:{net_change:0,absolute_sum:4,opposition_fraction:1},hydropathy_delta:{net_change:null,absolute_sum:null,opposition_fraction:null}},
 local_interaction_diagnostics:{pairs:[{positions:[7,8],wt_ca_distance_a:8,nearby_opposed_charge:true,ca_below_6_a:false,ca_above_25_a:false}]}});
item.resource_audit={single_threshold_audit:{review_band_count:4,otherwise_eligible_count:3,lower_exclusive_kcal_mol:0.5,upper_inclusive_kcal_mol:1,high_exposure_cutoff:0.6,high_rsasa_count:1,high_sidechain_exposure_count:2}};
const state={lang:'zh',result:{regions:[{id:1,residues:[7]}],residues:[],designs:[],joint_v4_ml:{sectors:[{sector_id:'VS001'}]}}};
state.result.joint_v4_ml.epitope_evidence=[{residue:{uniprot_pos:7},raw_score:0.12345,calibrated_score:-1.23456,status:'PREDICTED'}];
const ui=createRecommendationRenderer({state,h,fmt,fmtPct,notation});
const online=ui.renderSectorRecommendations([item]);
for(const text of ['DiscoTope 原始分数','DiscoTope 校准分数','0.123','-1.235','并非突变后组合的预测分数']) assert(online.includes(text),text);
for(const text of ['本区强扰动工程代表','42%','120.000 Å²','88.000 Å²','75%','10, 11','表位倾向 P 按计算表面排序与局部峰值标注','平均绝对化学变化','变化幅度（非带符号）','14.200 Å³','选择敏感性','最小–最大归一化','组合轴按位点绝对变化汇总','逐位带符号化学变化','9 Å','12 Å','15 Å']) assert(online.includes(text),text);
assert(online.includes('已纳入正式主面板（稳定性模型通过；扰动证据不全）'));
assert(online.includes('稳定性模型门槛已通过；机制扰动证据不全，缺少：芳香性变化'));
assert(!online.includes('模型记录需核对'));
for(const text of ['各取最强轴效用','方向与近邻诊断（不参与排序）','F/W/Y替换数／新引入数','近邻异向','WT Cα','距离按 WT Cα 计算。','单点稳定性门槛边界审计','此审计不放宽门槛']) assert(online.includes(text),text);
state.lang='en';
const english=ui.renderSectorRecommendations([item]);
for(const text of ['maximum axis utility','F/W/Y mutants / introductions','Distances use WT Cα coordinates.','does not relax gates']) assert(english.includes(text),text);
state.lang='zh';
assert(online.includes('盐桥'));
assert(!online.includes('尚未解释的保护或风险项'));
assert(!online.includes('尚未解释的审查项'));
assert(!online.includes('没有安全的专属候选，不能建立计算资格'));
assert(!online.includes('亲和力下降最大化'));
const legacy={sector_id:'VS001',candidate_id:'OLD',status:'PASSED',mutation_notations:['K7W'],ownership:'EXCLUSIVE'};
const oldHtml=ui.renderSectorRecommendations([legacy]);
assert(!oldHtml.includes('本区强扰动工程代表'));
assert(!oldHtml.includes('regional_mechanism_pareto_v1'));
const noCandidate={sector_id:'VS001',candidate_id:null,status:'NO_COMBINATION',formal_eligible:false,
 recommendation_policy:'regional_mechanism_pareto_v1',principles:['mechanism_pareto'],mutation_notations:[],
 impact_evidence:{policy:'regional_mechanism_pareto_v1',missing_axes:['charge_delta']}};
const gapHtml=ui.renderSectorRecommendations([noCandidate]);
assert(gapHtml.includes('本区扰动建议缺口'));
assert(gapHtml.includes('本区没有可用候选组合，未选出工程测试代表。'));
assert(!gapHtml.includes('本区强扰动工程代表'));
assert(!gapHtml.includes('没有安全的专属候选，不能建立计算资格'));
const D={joint_v4_ml:{selected_panel:{sector_recommendations:[item]}}};
const LANG='zh', designs=[];
const esc=h, t=key=>key, chip=(label,value)=>`<span>${label}:${value}</span>`;
const offline=stepDesign();
assert(offline.includes('本区强扰动工程代表'));
"""
    # The offline entry point calls the same shared renderer; assert that contract and execute it.
    assert "createRecommendationRenderer" in offline_step
    script = tmp_path / "regional_renderer.cjs"
    assert 'import { createRecommendationRenderer } from "./recommendations.js?v=' in online_app
    script.write_text(setup + renderer + "\n" + offline_step + "\n"
                      + "assert.equal(typeof createRecommendationRenderer, 'function');\n"
                      + "console.log(online);\n", encoding="utf-8")
    completed = subprocess.run([node, str(script)], check=True, capture_output=True, text=True, encoding="utf-8")
    assert "本区强扰动工程代表" in completed.stdout
