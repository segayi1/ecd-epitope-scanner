"""Shared formal/offline renderer preserves unknowns and escapes stored content."""
from pathlib import Path
import shutil
import subprocess

import pytest


def test_renderer_unknowns_legacy_membership_and_html_escaping(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for the shared renderer")
    renderer = (Path(__file__).resolve().parents[1] / "web/recommendations.js").read_text(encoding="utf-8").replace("export function", "function", 1)
    checks = r"""
const assert = require('node:assert/strict');
const state={lang:'zh',result:{regions:[{id:1,residues:[1]},{id:2,residues:[2]}],residues:[],designs:[{id:'D',mutations:[{wt:'A',position:1,mutant:'W'},{wt:'S',position:2,mutant:'Y'}],mutant_ecd_sequence:'<img src=x onerror=alert(1)>'}],joint_v4_ml:{sectors:[{sector_id:'VS001'},{sector_id:'VS002'}]}}};
const h=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const ui=createRecommendationRenderer({state,h,fmt:(v,d=2)=>Number(v).toFixed(d),fmtPct:v=>Math.round(v*100)+'%',notation:m=>m.wt+m.position+m.mutant});
const item={sector_id:'VS001',candidate_id:'<script>bad()</script>',design_id:'D',status:'NEAR_THRESHOLD',mutation_notations:['A1W','S2Y'],gate_details:[{reason:'UNKNOWN_INTERNAL_REASON',evidence_missing:true}],regional_geometry_support:{scenarios:{geodesic_9:{coverage_fraction:null},geodesic_12:{coverage_fraction:0,uncovered_positions:[1]},geodesic_15:{coverage_fraction:1,uncovered_positions:[]}}},gate_checks:[{gate:'combination_spurs',passed:null,required:false,observed:0.2,threshold:1},{gate:'candidate_connectivity',passed:true,description:'组合连通'}]};
let html=ui.renderSectorRecommendations([item]);
assert(html.includes('历史建议'));
assert(html.includes('跨区组合'));
assert(html.includes('未计算／缺少证据'));
assert(html.includes('0%'));
assert(html.includes('仅记录，不作为准入门槛'));
assert(html.includes('组合连通性'));
assert(html.includes('&lt;script&gt;bad()&lt;/script&gt;'));
assert(html.includes('&lt;img src=x onerror=alert(1)&gt;'));
assert(!html.includes('<script>bad()'));
assert(!html.includes('<img src=x'));
assert(!html.includes('UNKNOWN_INTERNAL_REASON'));
assert(!html.includes('geodesic_9'));
item.primary_comparison={baseline_candidate_id:'BASE',baseline_metrics:{mean_coverage_fraction:0.25},candidate_metrics:{mean_coverage_fraction:0.5},coverage_difference:0.25,worst_threshold_utilization_difference:null};
html=ui.renderSectorRecommendations([item]);
assert(html.includes('与本区首选的定量比较'));
assert(html.includes('+25.00 个百分点'));
assert(html.includes('无法比较：基准或证据缺失'));
state.lang='en';html=ui.renderSectorRecommendations([item]);
assert(html.includes('Cross-region combination'));
assert(html.includes('Recorded, not a qualification gate'));
assert(html.includes('Not computed / missing evidence'));
state.lang='zh';item.provenance={structure_evidence:{entry_id:'1AH1',method:'X-RAY DIFFRACTION',resolution_a:2,experimental:true,plddt_applicable:false}};
item.site_details=[{notation:'A1W',position:1,epitope_prior:0.4,epitope_supported:false,plddt:null,warnings:['SASA_METHOD_SETTINGS_UNVERIFIED'],structural_risk_flags:['SASA_METHOD_SETTINGS_UNVERIFIED']}];
item.gate_checks=[{gate:'candidate_connectivity',passed:true,required:true,observed:null,threshold:null,description:'\ufffd\ufffd',evidence_missing:false}];
html=ui.renderSectorRecommendations([item]);
assert(html.includes('实验结构：pLDDT不适用'));
assert(html.includes('局部精度尚未评估'));
assert(html.includes('不适用：规则检查，无数值门槛'));
assert(!html.includes('\ufffd'));
assert(html.includes('表面积计算参数尚未完整核对'));
assert.equal(html.split('表面积计算参数尚未完整核对').length-1,1);
assert(html.includes('有先验评分，未达支持标记条件'));
assert(html.includes('追溯编号、文件指纹与运行路径'));
assert(html.includes('构建序列'));
item.gate_details=[{gate:'higher_order_gain',reason:'INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN',observed:0.011,threshold:0.05,excess:0.039,evidence_missing:false}];
html=ui.renderSectorRecommendations([item]);
assert(html.includes('距最低要求还差 0.039'));
item.gate_details=[{gate:'sasa_structure_context_classification_changed',reason:'SASA_STRUCTURE_CONTEXT_CLASSIFICATION_CHANGED',evidence_missing:false}];
item.gate_checks=[{gate:'sasa_structure_context_classification_changed',reason:'SASA_STRUCTURE_CONTEXT_CLASSIFICATION_CHANGED',passed:false,required:true,evidence_missing:false}];
html=ui.renderSectorRecommendations([item]);
assert(html.includes('完整输入结构与ECD计算范围的表面暴露分类不一致'));
assert(!html.includes('尚未解释的审查项'));
item.gate_details=[{reason:'INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN',measurement_basis:'mean_regional_area_coverage_9_12_15',observed:0.035,threshold:0.05}];
html=ui.renderSectorRecommendations([item]);
assert(html.includes('同区最佳低阶几何对照'));
assert(html.includes('3.5 个百分点'));
assert(html.includes('5.0 个百分点'));
item.gate_details=[{reason:'HIGHER_ORDER_COVERAGE_BASELINE_MISSING',evidence_missing:true}];
html=ui.renderSectorRecommendations([item]);
assert(!html.includes('HIGHER_ORDER_COVERAGE_BASELINE_MISSING'));
assert(!html.includes('尚未解释的审查项'));
item.provenance.structure_evidence.method='Solution NMR';item.provenance.structure_evidence.resolution_a=null;
html=ui.renderSectorRecommendations([item]);
assert(html.includes('NMR不采用晶体学分辨率'));
assert(html.includes('结构集合离散度'));
assert(!html.includes('局部电子密度'));
item.provenance.models={spurs_multi:{provenance:{version:'1.0.0',torch_version:'2.5.1+cu124',git_revision:'70cbb99af5cdf1f6da392389aaf93fd245610efb',version_identification:{basis:'current_environment_metadata_supplement'}}},discotope:{provenance:{version:'1.0',model_revision:'DiscoTope-3.0',torch_version:'2.5.1+cu124'}}};
html=ui.renderSectorRecommendations([item]);
assert(html.includes('软件／安装包版本 1.0.0'));
  item.stability_evidence={provenance:{models:{spurs_multi:{version:'1.0.1',torch_version:'2.5.1+cu124'}}}};
  html=ui.renderSectorRecommendations([item]);
  assert(html.includes('软件／安装包版本 1.0.1'));
  delete item.stability_evidence.provenance;
  html=ui.renderSectorRecommendations([item]);
assert(html.includes('模型系列 3.0 · 安装包版本 1.0'));
assert(html.includes('当前环境补充识别，非历史运行时版本记录'));
assert(!html.includes('软件／安装包版本 2.5.1'));
assert(!html.includes('模型系列 1.0'));
assert(html.includes('70cbb99af5cdf1f6da392389aaf93fd245610efb'));
"""
    script = tmp_path / "renderer.cjs"
    script.write_text(renderer + checks, encoding="utf-8")
    subprocess.run([node, str(script)], check=True, capture_output=True, text=True)


def test_cross_region_surface_union_and_readable_selector(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for design display helpers")
    app = (Path(__file__).resolve().parents[1] / "web/app.js").read_text(encoding="utf-8")
    code = app[app.index("function regionsForDesign("):app.index("function paintBeads()")]
    code += app[app.index("function renderDesignSelect()"):app.index("async function selectDesign(")]
    checks = """
const assert=require('node:assert/strict');
const state={selectedDesign:'INTERNAL_HASH',previewPlans:new Map(),result:{designs:[{id:'INTERNAL_HASH',source_engine:'joint_v4_ml',region_id:3,design_role:'SECTOR_REVIEW_CANDIDATE',mutations:[{wt:'A',position:1,mutant:'W'},{wt:'Q',position:3,mutant:'Y'}]}],regions:[{id:1,residues:[1],envelope_residues:[1,2],footprint_residues:[1]},{id:2,residues:[3],envelope_residues:[3,4],footprint_residues:[3]}],residues:[1,2,3,4].map(uniprot_pos=>({uniprot_pos}))}};
const selectedDesign=()=>state.result.designs[0];
const select={innerHTML:'',options:[{value:'INTERNAL_HASH'}]};const $=()=>select;
const h=value=>String(value);const notation=m=>m.wt+m.position+m.mutant;
const recommendationLabel=zh=>zh;
const records=designSurfaceRecords();
assert.deepEqual(records.envelope.map(r=>r.uniprot_pos),[1,2,3,4]);
assert.deepEqual(records.footprint.map(r=>r.uniprot_pos),[1,3]);
assert.deepEqual(records.mutations.map(r=>r.uniprot_pos),[1,3]);
renderDesignSelect();
assert(select.innerHTML.includes('value="INTERNAL_HASH"'));
assert(!select.innerHTML.includes('>INTERNAL_HASH'));
assert(select.innerHTML.includes('A1W + Q3Y'));
assert(select.innerHTML.includes('分区 1、分区 2'));
assert(select.innerHTML.includes('待审查建议'));
"""
    script = tmp_path / "cross_region_display.cjs"
    script.write_text(code + checks, encoding="utf-8")
    subprocess.run([node, str(script)], check=True, capture_output=True, text=True)
