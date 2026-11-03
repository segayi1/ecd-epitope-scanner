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
const item={sector_id:'VS001',candidate_id:'<script>bad()</script>',design_id:'D',status:'PASSED',formal_eligible:true,mutation_notations:['A1W','S2Y'],gate_details:[{reason:'UNKNOWN_INTERNAL_REASON',evidence_missing:true}],regional_geometry_support:{scenarios:{geodesic_9:{coverage_fraction:null},geodesic_12:{coverage_fraction:0,uncovered_positions:[1]},geodesic_15:{coverage_fraction:1,uncovered_positions:[]}}},gate_checks:[{gate:'combination_spurs',passed:null,required:false,observed:0.2,threshold:1},{gate:'candidate_connectivity',passed:true,description:'组合连通'}]};
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
assert(html.includes('data-select-recommendation="&lt;script&gt;bad()&lt;/script&gt;"'));
assert(html.includes('加入空间比较'));
state.recommendationSelections = new Map([
  ['C1', {id:'C1',label:'A1W + S2Y',status:'通过计算门槛',mutations:[{wt:'A',position:1,mutant:'W'},{wt:'S',position:2,mutant:'Y'}]}],
  ['C2', {id:'C2',label:'A1E + S2Y',status:'待补算',mutations:[{wt:'A',position:1,mutant:'E'},{wt:'S',position:2,mutant:'Y'}]}],
]);
const selectionHtml=ui.renderSectorRecommendations([item]);
assert(selectionHtml.includes('已选组合空间关系'));
assert(selectionHtml.includes('同位点替换冲突'));
assert(selectionHtml.includes('位点重合'));
assert(selectionHtml.includes('data-clear-recommendations'));
const repeatCandidate={...item,candidate_id:'REPEAT',status:'NEAR_THRESHOLD',mutation_notations:['A1W'],
  execution_protocol:'ordinary_three_principle_adaptive_v1',
  stability_evidence:{provenance:{repeat_audit:{foldx_repeat_count:3,foldx_repeat_slots:[0,1,2],foldx_repeat_values_kcal_mol:[0.1,0.2,0.3],foldx_evidence_complete:false}}},
  gate_details:[{gate:'foldx_repeat_evidence_incomplete',reason:'FOLDX_REPEAT_EVIDENCE_INCOMPLETE',evidence_missing:false}]};
const repeatHtml=ui.renderSectorRecommendations([repeatCandidate]);
assert(repeatHtml.includes('FoldX重复 3/3'));
assert(repeatHtml.includes('槽值 0: 0.100 / 1: 0.200 / 2: 0.300 kcal/mol'));
assert(repeatHtml.includes('待补算'));
repeatCandidate.execution_protocol={name:'ordinary_three_principle_adaptive_v1',foldx_repeat_count:3};
const adaptiveThreeHtml=ui.renderSectorRecommendations([repeatCandidate]);
assert(adaptiveThreeHtml.includes('FoldX重复 3/3'));
assert(adaptiveThreeHtml.includes('三次同侧通过') === false); // incomplete evidence remains visible
const adaptiveTen={...repeatCandidate,candidate_id:'ADAPTIVE-TEN',execution_protocol:{name:'ordinary_adaptive_v1',foldx_repeat_count:10},
  stability_repeat_diagnostic:{sampling_decision_diagnostic_status:'COMPLETE',sampling_decision_sensitive:false,foldx_repeat_count:10,
    foldx_repeat_slots:Array.from({length:10},(_,i)=>i),foldx_repeat_values_kcal_mol:Array(10).fill(0.2),foldx_repeat_gate_kcal_mol:2,
    foldx_repeat_min_kcal_mol:.2,foldx_repeat_max_kcal_mol:.2,foldx_repeat_above_gate_count:0,foldx_repeat_at_or_below_gate_count:10}};
const adaptiveTenHtml=ui.renderSectorRecommendations([adaptiveTen]);
assert(adaptiveTenHtml.includes('FoldX重复 10/10'));
assert(adaptiveTenHtml.includes('十槽中位数通过；采样不确定性'));
assert(adaptiveTenHtml.includes('十槽采样不确定性'));
const mixedLegacy={...repeatCandidate,candidate_id:'MIXED-LEGACY',execution_protocol:undefined};
const mixedHtml=ui.renderSectorRecommendations([mixedLegacy,repeatCandidate]);
assert(mixedHtml.includes('FoldX重复 3/5'));
assert(mixedHtml.includes('FoldX重复 3/3'));
state.result.joint_v4_ml.execution_protocol='ordinary_three_principle_adaptive_v1';
assert(ui.renderSectorRecommendations([mixedLegacy]).includes('FoldX重复 3/3'));
state.result.joint_v4_ml.execution_protocol={name:'ordinary_adaptive_v1',foldx_repeat_count:10};
assert(ui.renderSectorRecommendations([mixedLegacy]).includes('FoldX重复 3/10'));
delete state.result.joint_v4_ml.execution_protocol;
state.result.execution_protocol={name:'ordinary_three_principle_adaptive_v1'};
const runAdaptiveTen={...item,candidate_id:'RUN-ADAPTIVE-TEN',status:'PASSED',mutation_notations:['A1W'],
  stability_repeat_diagnostic:{sampling_decision_diagnostic_status:'COMPLETE',sampling_decision_sensitive:false,foldx_repeat_count:10,
    foldx_repeat_slots:Array.from({length:10},(_,i)=>i),foldx_repeat_values_kcal_mol:Array(10).fill(.2),foldx_repeat_gate_kcal_mol:2,
    foldx_repeat_min_kcal_mol:.2,foldx_repeat_max_kcal_mol:.2,foldx_repeat_above_gate_count:0,foldx_repeat_at_or_below_gate_count:10}};
const runAdaptiveTenHtml=ui.renderSectorRecommendations([runAdaptiveTen]);
assert(runAdaptiveTenHtml.includes('FoldX重复 10/10'));
assert(runAdaptiveTenHtml.includes('十槽采样不确定性'));
const auditLegacyThree={...runAdaptiveTen,candidate_id:'AUDIT-LEGACY-THREE',
  stability_repeat_diagnostic:{execution_protocol:'legacy',foldx_repeat_count:3,foldx_repeat_slots:[0,1,2],
    foldx_repeat_values_kcal_mol:[.1,.2,.3],foldx_evidence_complete:false}};
const auditLegacyHtml=ui.renderSectorRecommendations([auditLegacyThree]);
assert(auditLegacyHtml.includes('FoldX重复 3/5'));
delete state.result.execution_protocol;
assert(repeatHtml.includes('FoldX组合评分重复次数不足，待补算'));
assert(repeatHtml.includes('待补算'));
const oldCtla4Risk={...item,candidate_id:'CTLA4-LEGACY',status:'PASSED',design_id:'D',mutation_notations:['E59Y','K65F'],
  stability_repeat_diagnostic:{recommendation_risk_label:'SAMPLING_SENSITIVE_REVIEW',sampling_decision_diagnostic_status:'COMPLETE',
    sampling_decision_sensitive:true,foldx_repeat_count:5,foldx_repeat_slots:[0,1,2,3,4],
    foldx_repeat_values_kcal_mol:[2.99488,.276247,3.4539,1.21002,.710114],
    foldx_repeat_slot_values_kcal_mol:[2.99488,.276247,3.4539,1.21002,.710114],foldx_repeat_gate_kcal_mol:2,
    foldx_repeat_min_kcal_mol:.276247,foldx_repeat_max_kcal_mol:3.4539,foldx_repeat_above_gate_count:2,
    foldx_repeat_at_or_below_gate_count:3}};
const oldCtla4Html=ui.renderSectorRecommendations([oldCtla4Risk]);
assert(oldCtla4Html.includes('采样需复核'));
assert(oldCtla4Html.includes('门限 ≤ 2.000'));
assert(oldCtla4Html.includes('范围 0.276–3.454'));
assert(oldCtla4Html.includes('门上 2/5'));
assert(oldCtla4Html.includes('槽值 0: 2.995 / 1: 0.276 / 2: 3.454 / 3: 1.210 / 4: 0.710'));
assert(oldCtla4Html.includes('通过计算门槛'));
const incompleteRisk={...item,candidate_id:'UNKNOWN-REPEAT',status:'PASSED',mutation_notations:['A1W'],
  stability_repeat_diagnostic:{recommendation_risk_label:'REPEAT_EVIDENCE_INCOMPLETE',sampling_decision_diagnostic_status:'INCOMPLETE',
    sampling_decision_sensitive:null,foldx_repeat_count:4,foldx_repeat_slots:[0,1,2,3,4],foldx_repeat_values_kcal_mol:[1,1,1,1],
    foldx_repeat_slot_values_kcal_mol:[1,1,1,1,null],foldx_repeat_gate_kcal_mol:2,foldx_repeat_min_kcal_mol:1,
    foldx_repeat_max_kcal_mol:1,foldx_repeat_above_gate_count:0,foldx_repeat_at_or_below_gate_count:4}};
const incompleteHtml=ui.renderSectorRecommendations([incompleteRisk]);
assert(incompleteHtml.includes('重复证据不完整'));
assert(incompleteHtml.includes('4/5'));
assert(incompleteHtml.includes('槽值 0: 1.000 / 1: 1.000 / 2: 1.000 / 3: 1.000 / 4: 缺失'));
assert(!incompleteHtml.includes('低风险'));
const withinGate={...oldCtla4Risk,candidate_id:'WITHIN-GATE',recommendation_risk_label:'NO_WITHIN_RUN_GATE_CROSSING',
  stability_repeat_diagnostic:{...oldCtla4Risk.stability_repeat_diagnostic,recommendation_risk_label:'NO_WITHIN_RUN_GATE_CROSSING',sampling_decision_sensitive:false,
    foldx_repeat_values_kcal_mol:[.1,.2,.3,.4,.5],foldx_repeat_slot_values_kcal_mol:[.1,.2,.3,.4,.5],foldx_repeat_min_kcal_mol:.1,
    foldx_repeat_max_kcal_mol:.5,foldx_repeat_above_gate_count:0,foldx_repeat_at_or_below_gate_count:5}};
const withinGateHtml=ui.renderSectorRecommendations([withinGate]);
assert(withinGateHtml.includes('本批未见门限跨越'));
assert(!withinGateHtml.includes('采样需复核'));
assert(!withinGateHtml.includes('跨批稳定'));
const unknownLegacy={...item,candidate_id:'OLD-WITHOUT-DIAGNOSTIC',status:'PASSED',mutation_notations:['A1W'],
  stability_evidence:{foldx_ddg_kcal_mol:1.2,provenance:{repeat_audit:{foldx_repeat_count:5,foldx_repeat_slots:[0,1,2,3,4],
    foldx_repeat_values_kcal_mol:[.2,.4,.6,1.2,1.8],foldx_evidence_complete:true}}}};
const unknownLegacyHtml=ui.renderSectorRecommendations([unknownLegacy]);
assert(unknownLegacyHtml.includes('采样诊断未知'));
assert(unknownLegacyHtml.includes('门限未记录'));
assert(!unknownLegacyHtml.includes('低风险'));
assert(!unknownLegacyHtml.includes('跨批稳定'));
const crossRegion={...item,candidate_id:'CROSS-IMPACT',status:'PASSED',mutation_notations:['A1W','S2Y'],
  impact_evidence:{axes:{charge_delta:1,volume_delta_a3:2,hydropathy_delta:3,donor_delta:4,acceptor_delta:5,aromatic_delta:6},exposure_min:.31,
    signed_change_summary:{charge_delta:{status:'KNOWN',net_change:0,absolute_sum:2,opposition_fraction:1},hydropathy_delta:{status:'KNOWN',net_change:2,absolute_sum:2,opposition_fraction:0}},
    aromatic_mutant_count:2,introduced_aromatic_count:1,
    local_interaction_diagnostics:{basis:'WT_CA_COORDINATES',nearby_ca_cutoff_a:10,pairs:[{positions:[1,2],wt_ca_distance_a:8,
      opposed_changes:{charge_delta:true,hydropathy_delta:false},nearby_opposed_charge:true,nearby_opposed_hydropathy:false,ca_below_6_a:false,ca_above_25_a:false}] }},
  site_details:[{notation:'A1W',position:1,chemical_descriptors:{charge_delta:-1,volume_delta_a3:2,hydropathy_delta:3,donor_delta:4,acceptor_delta:5,aromatic_delta:6}},
    {notation:'S2Y',position:2,chemical_descriptors:{charge_delta:1,volume_delta_a3:2,hydropathy_delta:3,donor_delta:4,acceptor_delta:5,aromatic_delta:6}}]};
const crossHtml=ui.renderSectorRecommendations([crossRegion]);
for(const label of ['化学变化与暴露度','电荷变化','体积变化','疏水性变化','氢键供体变化','氢键受体变化','芳香性变化','电荷变化方向及距离分类','方向相反 · 10 Å内','方向相同','6–25 Å']) assert(crossHtml.includes(label),label);
for(const value of ['1.000','2.000','3.000','4.000','5.000','6.000','8.000 Å']) assert(crossHtml.includes(value),value);
assert(!crossHtml.includes('不据此推断碰撞') && !crossHtml.includes('未标记近邻异向电荷变化'));
const histidineRange={...crossRegion,candidate_id:'HIS-RANGE',mutation_notations:['H10A','A1W'],
  impact_evidence:{axes:{charge_delta:null,volume_delta_a3:2,hydropathy_delta:null,donor_delta:null,acceptor_delta:null,aromatic_delta:1},missing_axes:['charge_delta','hydropathy_delta','donor_delta','acceptor_delta'],
    display_chemistry_ranges:{schema_version:'display_chemistry_ranges_v1',basis:'HIS_NEUTRAL_CHARGE0_DONOR1_ACCEPTOR1_VS_PROTONATED_CHARGE1_DONOR2_ACCEPTOR0_NO_PH_ASSUMPTION_v1',
      protonation_assumption:'neutral/protonated states enumerated independently for each mutated His; no pH assignment',
      microstates:[{notation:'H10A',position:10,histidine_side:'wild_type',states:[{state:'neutral',charge:0,donor:1,acceptor:1},{state:'protonated',charge:1,donor:2,acceptor:0}]}],
      scenarios:[{scenario_id:'S0',histidine_states:[{notation:'H10A',side:'wild_type',state:'neutral'}],mean_absolute_delta:{charge_delta:0.5,volume_delta_a3:2,hydropathy_delta:1,donor_delta:1,acceptor_delta:1,aromatic_delta:1},
        per_site_signed_deltas:[{notation:'H10A',position:10,chemical_descriptors:{charge_delta:0,volume_delta_a3:0,hydropathy_delta:0,donor_delta:0,acceptor_delta:0,aromatic_delta:0}}],
        signed_change_summary:{charge_delta:{net_change:0,absolute_sum:1,opposition_fraction:1},hydropathy_delta:{net_change:0,absolute_sum:2,opposition_fraction:.5}},
        local_interaction_diagnostics:{pairs:[{positions:[10,1],wt_ca_distance_a:8,nearby_opposed_charge:true,nearby_opposed_hydropathy:false}]} }],
      ranges:{axes:{charge_delta:{minimum:.5,maximum:1},volume_delta_a3:{minimum:2,maximum:2},hydropathy_delta:{minimum:1,maximum:2},donor_delta:{minimum:1,maximum:2},acceptor_delta:{minimum:0,maximum:1},aromatic_delta:{minimum:1,maximum:1}},
        per_site_signed_deltas:{H10A:{charge_delta:{minimum:0,maximum:1},volume_delta_a3:{minimum:0,maximum:0},hydropathy_delta:{minimum:-1,maximum:0},donor_delta:{minimum:0,maximum:1},acceptor_delta:{minimum:-1,maximum:0},aromatic_delta:{minimum:0,maximum:0}}},
        signed_change_summary:{charge_delta:{net_change:{minimum:0,maximum:1},absolute_sum:{minimum:1,maximum:2},opposition_fraction:{minimum:.5,maximum:1}},hydropathy_delta:{net_change:{minimum:0,maximum:0},absolute_sum:{minimum:2,maximum:3},opposition_fraction:{minimum:0,maximum:.5}}},
        local_interactions:[{positions:[10,1],wt_ca_distance_a:8,nearby_opposed_charge:{classification:'STATE_DEPENDENT'},nearby_opposed_hydropathy:{classification:'CONSISTENT_FALSE'}}]}}},
  site_details:[{notation:'H10A',position:10,chemical_descriptors:{charge_delta:null,volume_delta_a3:0,hydropathy_delta:null,donor_delta:null,acceptor_delta:null,aromatic_delta:0}},
    {notation:'A1W',position:1,chemical_descriptors:{charge_delta:0,volume_delta_a3:2,hydropathy_delta:1,donor_delta:1,acceptor_delta:1,aromatic_delta:1}}]};
const hisHtml=ui.renderSectorRecommendations([histidineRange]);
for(const text of ['组氨酸微态情景范围','中性／质子化状态给出范围','质子化（电荷 1，供体 2，受体 0）','逐位带符号变化范围','signed变化汇总范围','是否近邻异向随微态变化']) assert(hisHtml.includes(text),text);
for(const value of ['0.500 Δe – 1.000 Δe','-1.000 – 0.000','50% – 100%','8.000 Å']) assert(hisHtml.includes(value),value);
assert.match(hisHtml, /电荷变化<\/td><td>0\.500 Δe – 1\.000 Δe/);
assert(hisHtml.includes('氢键供体变化</td><td>1.000 – 2.000'));
assert(hisHtml.includes('电荷变化: 0.000 Δe – 1.000 Δe'));
assert(hisHtml.includes('未计算／缺少证据'));
assert.equal(histidineRange.impact_evidence.axes.charge_delta,null);
assert.equal(histidineRange.site_details[0].chemical_descriptors.charge_delta,null);
const legacyCopy={...item,candidate_id:'OLD-COPY',mutation_notations:['A1W'],recommendation_reasons:['保留原覆盖优先方案作为面积对照；几何邻域不等于实际抗体受影响面积。']};
const legacyChem={...item,candidate_id:'OLD-CHEM',mutation_notations:['A1W'],recommendation_reasons:['保留charge_delta变化较大的机制备选；不能据此推断结合自由能变化。']};
const cleanCopyHtml=ui.renderSectorRecommendations([legacyCopy,legacyChem]);
assert(cleanCopyHtml.includes('保留原覆盖优先方案作为面积对照'));
assert(cleanCopyHtml.includes('按电荷变化较大的原则选择备选'));
assert(!cleanCopyHtml.includes('实际抗体受影响面积') && !cleanCopyHtml.includes('结合自由能变化') && !cleanCopyHtml.includes('变化变化'));
assert(ui.surfaceGateDescription({joint_v4_ml:{surface_gate:{policy:'native_total_rsasa_only',native_total_rsasa_threshold:0.25}}},'zh').includes('侧链 SASA 为描述证据，不是表面门'));
assert(ui.surfaceGateDescription({joint_v4_ml:{surface_gate:{policy:'native_total_rsasa_only',native_total_rsasa_threshold:0.25}}},'en').includes('native total rSASA ≥ 0.25'));
assert(ui.surfaceGateDescription({methods:{rsasa_threshold:0.25,sidechain_sasa_threshold_A2:20}},'zh').includes('侧链 SASA ≥ 20 Å²'));
const qualified={...item,candidate_id:'FORMAL1',status:'PASSED',ownership:'EXCLUSIVE',mutation_notations:['A1W','S2Y'],principles:[],gate_checks:[],gate_details:[],regional_geometry_support:{scenarios:{}},site_details:[]};
const primary={...item,candidate_id:null,design_id:null,status:'NEAR_THRESHOLD',recommendation_policy:'regional_mechanism_pareto_v1',mutation_notations:[],principles:[],formal_candidate_count:1,confirmed_alternative_count:0,all_formal_local_candidates:[qualified]};
html=ui.renderSectorRecommendations([primary]);
assert(html.includes('机制主代表状态: 待审查 · 未通过'));
assert(html.includes('本区正式合格候选: 1'));
assert(html.includes('全部正式展示合格本区候选（含主代表） · 1'));
assert(html.includes('FORMAL1'));
assert(html.includes('本区仍有正式合格候选'));
assert(!html.includes('当前未得到本区候选组合'));
const reviewableFormal={...qualified,candidate_id:'SENSITIVE-LOCAL',stability_repeat_diagnostic:{
  recommendation_risk_label:'SAMPLING_SENSITIVE_REVIEW',sampling_decision_diagnostic_status:'COMPLETE',sampling_decision_sensitive:true}};
const evidenceLimited={...qualified,candidate_id:'LIMITED',status:'EVIDENCE_LIMITED'};
const unknownRisk={...qualified,candidate_id:'UNKNOWN-RISK',stability_repeat_diagnostic:{recommendation_risk_label:'SAMPLING_DIAGNOSTIC_UNKNOWN'}};
const displayReview={...qualified,candidate_id:'DISPLAY-REVIEW',display_eligible:false,
  display_exclusion_reasons:['STRUCTURAL_CONFIDENCE_REVIEW_REQUIRED']};
const reviewCard={...qualified,candidate_id:'SENSITIVE-PRIMARY',stability_repeat_diagnostic:{
  recommendation_risk_label:'SAMPLING_SENSITIVE_REVIEW',sampling_decision_diagnostic_status:'COMPLETE',sampling_decision_sensitive:true},
  all_formal_local_candidates:[qualified,reviewableFormal,evidenceLimited,unknownRisk,displayReview]};
html=ui.renderSectorRecommendations([reviewCard]);
const formalHtml=html.split('class="regional-review-candidates"')[0];
assert(!formalHtml.includes('data-recommendation-role="primary"'));
assert(!formalHtml.includes('data-recommendation-role="alternative"'));
assert(html.includes('regional-review-candidates'));
assert(html.includes('SENSITIVE-PRIMARY') && html.includes('SENSITIVE-LOCAL') && html.includes('LIMITED') && html.includes('UNKNOWN-RISK') && html.includes('DISPLAY-REVIEW'));
assert(html.includes('全部正式展示合格本区候选（含主代表） · 1'));
const formalLocalSection=html.match(/<details class="formal-local-candidates">([\s\S]*?)<\/details>/)?.[1] || '';
assert(formalLocalSection.includes('FORMAL1'));
assert(!formalLocalSection.includes('SENSITIVE-LOCAL') && !formalLocalSection.includes('LIMITED') && !formalLocalSection.includes('UNKNOWN-RISK') && !formalLocalSection.includes('DISPLAY-REVIEW'));
assert(!html.includes('data-select-recommendation="SENSITIVE-PRIMARY"'));
assert(!html.includes('data-select-recommendation="SENSITIVE-LOCAL"'));
assert(!html.includes('data-select-recommendation="DISPLAY-REVIEW"'));
assert.equal(ui.recommendationSelectionRecord(reviewableFormal,null),null);
assert.equal(ui.recommendationSelectionRecord(evidenceLimited,null),null);
assert.equal(ui.recommendationSelectionRecord(unknownRisk,null),null);
assert.equal(ui.recommendationSelectionRecord(displayReview,null),null);
assert.equal(ui.recommendationSelectionRecord({...qualified,display_status:'DISPLAY_REVIEW'},null),null);
assert.equal(ui.recommendationSelectionRecord(qualified,null).id,'FORMAL1');
state.result.execution_protocol={name:'ordinary_site_first_v2'};
assert.equal(ui.recommendationSelectionRecord({...qualified,formal_eligible:undefined},null),null);
assert.equal(ui.recommendationSelectionRecord({...qualified,stability_repeat_diagnostic:undefined},null),null);
const siteFirstPass={...qualified,stability_repeat_diagnostic:{
  sampling_decision_diagnostic_status:'COMPLETE',sampling_decision_sensitive:false,
  cross_execution_conflict:null,cross_execution_evidence_status:'NOT_ASSESSED'}};
assert.equal(ui.recommendationSelectionRecord(siteFirstPass,null).id,'FORMAL1');
assert.equal(ui.recommendationSelectionRecord({...siteFirstPass,stability_repeat_diagnostic:{
  ...siteFirstPass.stability_repeat_diagnostic,sampling_decision_sensitive:true}},null),null);
assert.equal(ui.recommendationSelectionRecord({...siteFirstPass,stability_repeat_diagnostic:{
  ...siteFirstPass.stability_repeat_diagnostic,cross_execution_conflict:true}},null),null);
const policyEnabledThreeSlot={...qualified,candidate_id:'POLICY-THREE-SLOT',stability_repeat_diagnostic:{
  repeat_confirmation_policy:'ten_slot_median_v1',repeat_confirmation_status:'NOT_CONFIRMED',
  sampling_decision_diagnostic_status:'COMPLETE',sampling_decision_sensitive:false,foldx_repeat_count:3,
  foldx_repeat_slots:[0,1,2],foldx_repeat_slot_values_kcal_mol:[.2,.3,.4],foldx_repeat_gate_kcal_mol:2,
  cross_execution_conflict:false}};
assert.equal(ui.recommendationSelectionRecord(policyEnabledThreeSlot,null).id,'POLICY-THREE-SLOT');
const policyEnabledSensitiveTen={...qualified,candidate_id:'POLICY-SENSITIVE-TEN',stability_repeat_diagnostic:{
  repeat_confirmation_policy:'ten_slot_median_v1',repeat_confirmation_status:'NOT_CONFIRMED',
  recommendation_risk_label:'SAMPLING_SENSITIVE_REVIEW',sampling_decision_diagnostic_status:'COMPLETE',
  sampling_decision_sensitive:true,foldx_repeat_count:10,foldx_repeat_slots:Array.from({length:10},(_,i)=>i),
  foldx_repeat_slot_values_kcal_mol:[2.5,1,1,1,1,1,1,1,1,1],foldx_repeat_gate_kcal_mol:2,
  cross_execution_conflict:false}};
assert.equal(ui.recommendationSelectionRecord(policyEnabledSensitiveTen,null),null);
assert(ui.renderSectorRecommendations([policyEnabledSensitiveTen]).includes('采样需复核'));
const tenSlotAudit={repeat_confirmation_policy:'ten_slot_median_v1',repeat_confirmation_status:'CONFIRMED_PASS',
  recommendation_risk_label:'TEN_SLOT_MEDIAN_CONFIRMED',sampling_decision_diagnostic_status:'COMPLETE',
  sampling_decision_sensitive:true,foldx_repeat_count:10,foldx_repeat_slots:[0,1,2,3,4,5,6,7,8,9],
  foldx_repeat_slot_values_kcal_mol:[2.5,1,1,1,1,1,1,1,1,1],foldx_repeat_gate_kcal_mol:2,
  foldx_repeat_min_kcal_mol:1,foldx_repeat_max_kcal_mol:2.5,foldx_repeat_above_gate_count:1,
  foldx_repeat_at_or_below_gate_count:9,cross_execution_conflict:false};
const tenSlotPass={...qualified,candidate_id:'TEN-SLOT-PASS',stability_repeat_diagnostic:tenSlotAudit};
assert.equal(ui.recommendationSelectionRecord(tenSlotPass,null).id,'TEN-SLOT-PASS');
html=ui.renderSectorRecommendations([tenSlotPass]);
assert(html.includes('十槽复算确认通过；存在采样波动'));
assert(html.includes('范围 1.000–2.500'));
assert(html.includes('门上 1/10') && html.includes('门内 9/10'));
assert.equal(ui.recommendationSelectionRecord({...tenSlotPass,stability_repeat_diagnostic:{
  ...tenSlotAudit,foldx_repeat_slots:[1,2,3,4,5,6,7,8,9,10]}},null),null);
assert.equal(ui.recommendationSelectionRecord({...tenSlotPass,stability_repeat_diagnostic:{
  ...tenSlotAudit,foldx_repeat_slot_values_kcal_mol:[2.5,2.2,1,1,1,1,1,1,1,1]}},null),null);
assert.equal(ui.recommendationSelectionRecord({...tenSlotPass,stability_repeat_diagnostic:{
  ...tenSlotAudit,foldx_repeat_slot_values_kcal_mol:[2.5,1,1,2.1,2.2,2.3,2.4,2.6,1,1]}},null),null);
assert.equal(ui.recommendationSelectionRecord({...tenSlotPass,stability_repeat_diagnostic:{
  ...tenSlotAudit,cross_execution_conflict:true}},null),null);
delete state.result.execution_protocol;
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


def test_all_twelve_regions_render_qualified_options_and_separate_failed_reviews(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for the shared renderer")
    renderer = (Path(__file__).resolve().parents[1] / "web/recommendations.js").read_text(encoding="utf-8").replace("export function", "function", 1)
    checks = r"""
const assert=require('node:assert/strict');
const state={lang:'en',result:{regions:[],residues:[],designs:[],joint_v4_ml:{sectors:[]}}};
const rows=Array.from({length:12},(_,i)=>({sector_id:`VS${i+1}`,candidate_id:`PASS-${i+1}`,formal_eligible:true,status:'PASSED',
 mutation_notations:[`D${i*4+1}A`,`D${i*4+2}A`],principles:['mechanism_pareto','charge_extreme'],
 alternatives:[],all_formal_local_candidates:[],selected_in_panel:i<4}));
rows[0].alternatives=[{candidate_id:'FAILED-ALTERNATIVE',formal_eligible:false,status:'NEAR_THRESHOLD',mutation_notations:['D2A','D3A'],principles:['lower_structural_risk']}];
rows[1]={...rows[1],candidate_id:'FAILED-PRIMARY',status:'NEAR_THRESHOLD',formal_eligible:false,
 alternatives:[{candidate_id:'QUALIFIED-ALTERNATIVE',formal_eligible:true,status:'PASSED',mutation_notations:['D6A','D7A'],principles:['lower_structural_risk']}]};
const ui=createRecommendationRenderer({state,h:x=>String(x??''),fmt:x=>String(x),fmtPct:x=>String(x),notation:m=>m.wt+m.position+m.mutant});
const html=ui.renderSectorRecommendations(rows);
assert.equal((html.match(/class="recommendation-region"/g)||[]).length,12);
assert(html.includes('PASS-12'));assert(html.includes('QUALIFIED-ALTERNATIVE'));
assert(html.includes('without the global compact-panel limit'));
const first=html.split('class="recommendation-region"')[1];
assert(first.indexOf('class="regional-review-candidates"')<first.indexOf('FAILED-ALTERNATIVE'));
const second=html.split('class="recommendation-region"')[2];
assert(second.indexOf('QUALIFIED-ALTERNATIVE')<second.indexOf('class="regional-review-candidates"'));
assert(second.indexOf('class="regional-review-candidates"')<second.indexOf('FAILED-PRIMARY'));
"""
    script = tmp_path / "twelve_regions.cjs"
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
