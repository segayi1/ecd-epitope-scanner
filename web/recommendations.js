export function createRecommendationRenderer({ state, h, fmt, fmtPct, notation }) {
const recommendationColors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7A5195", "#B94E43"];
function recommendationRepeatAudit(item, design) {
  return item?.stability_repeat_diagnostic || item?.stability_evidence?.provenance?.repeat_audit
    || item?.impact_evidence?.stability_repeat_audit || design?.foldx_diagnostics?.stability_repeat_diagnostic
    || design?.foldx_diagnostics?.provenance?.repeat_audit;
}
function recommendationRepeatCount(item) {
  const audit = recommendationRepeatAudit(item);
  const protocol = item?.execution_protocol || audit?.execution_protocol || state.result?.execution_protocol
    || state.result?.joint_v4_ml?.execution_protocol || state.result?.request?.execution_protocol;
  const name = typeof protocol === "string" ? protocol : protocol?.name;
  if (["ordinary_three_principle_adaptive_v1", "ordinary_adaptive_v1", "ordinary_site_first_v2"].includes(name)) {
    const declared = Number(typeof protocol === "object" ? protocol.foldx_repeat_count : NaN);
    const actual = Number(audit?.foldx_repeat_count ?? audit?.foldx_repeat_slots?.length);
    return [3, 10].includes(declared) ? declared : [3, 10].includes(actual) ? actual : 3;
  }
  return 5;
}
function recommendationRiskCode(item, design) {
  const audit = recommendationRepeatAudit(item, design);
  if (recommendationTenSlotConfirmation(item, design)) return "TEN_SLOT_MEDIAN_CONFIRMED";
  if (recommendationHasTenSlotMarker(item, design)) return "SAMPLING_DIAGNOSTIC_UNKNOWN";
  const explicit = item?.recommendation_risk_label || audit?.recommendation_risk_label
    || design?.foldx_diagnostics?.recommendation_risk_label;
  if (explicit) return explicit;
  if (!audit) return item?.stability_evidence || item?.impact_evidence || design?.foldx_diagnostics
    ? "SAMPLING_DIAGNOSTIC_UNKNOWN" : null;
  if (audit.sampling_decision_diagnostic_status === "COMPLETE") {
    return audit.sampling_decision_sensitive === true ? "SAMPLING_SENSITIVE_REVIEW"
      : audit.sampling_decision_sensitive === false ? "NO_WITHIN_RUN_GATE_CROSSING" : "SAMPLING_DIAGNOSTIC_UNKNOWN";
  }
  if (audit.sampling_decision_diagnostic_status === "INCOMPLETE"
    || audit.foldx_evidence_complete === false || (Number.isInteger(audit.foldx_repeat_count) && audit.foldx_repeat_count < recommendationRepeatCount(item))
    || (Array.isArray(audit.foldx_repeat_slots) && audit.foldx_repeat_slots.length < recommendationRepeatCount(item))) return "REPEAT_EVIDENCE_INCOMPLETE";
  return "SAMPLING_DIAGNOSTIC_UNKNOWN";
}
function recommendationTenSlotConfirmation(item, design) {
  const audit = recommendationRepeatAudit(item, design);
  if (audit?.repeat_confirmation_policy !== "ten_slot_median_v1"
    || audit?.repeat_confirmation_status !== "CONFIRMED_PASS"
    || audit?.recommendation_risk_label !== "TEN_SLOT_MEDIAN_CONFIRMED"
    || audit?.sampling_decision_diagnostic_status !== "COMPLETE"
    || audit?.sampling_decision_sensitive !== true
    || audit?.foldx_repeat_gate_kcal_mol !== 2
    || audit?.cross_execution_conflict === true) return false;
  const slots = audit.foldx_repeat_slots;
  const values = audit.foldx_repeat_slot_values_kcal_mol;
  if (!Array.isArray(slots) || slots.length !== 10 || slots.some((slot, index) => slot !== index)
    || !Array.isArray(values) || values.length !== 10 || values.some((value) => typeof value !== "number" || !Number.isFinite(value))) return false;
  const sorted = [...values].sort((a, b) => a - b);
  return values.slice(0, 3).filter((value) => value > 2).length === 1
    && (sorted[4] + sorted[5]) / 2 <= 2;
}
function recommendationHasTenSlotMarker(item, design) {
  const audit = recommendationRepeatAudit(item, design);
  return (audit?.repeat_confirmation_status != null && audit.repeat_confirmation_status !== "NOT_CONFIRMED")
    || audit?.recommendation_risk_label === "TEN_SLOT_MEDIAN_CONFIRMED";
}
function recommendationRiskLabelText(item, design) {
  const label = recommendationRiskCode(item, design);
  const repeatCount = recommendationRepeatCount(item);
  return ({
    TEN_SLOT_MEDIAN_CONFIRMED: recommendationLabel("十槽复算确认通过；存在采样波动", "Ten-slot confirmation passed; sampling variation observed"),
    SAMPLING_SENSITIVE_REVIEW: recommendationLabel("采样需复核", "Sampling review needed"),
    REPEAT_EVIDENCE_INCOMPLETE: recommendationLabel("重复证据不完整", "Repeat evidence incomplete"),
    SAMPLING_DIAGNOSTIC_UNKNOWN: recommendationLabel("采样诊断未知", "Sampling diagnostic unknown"),
    NO_WITHIN_RUN_GATE_CROSSING: repeatCount === 3
      ? recommendationLabel("三次同侧通过", "Three repeats pass on the same side")
      : repeatCount === 10 ? recommendationLabel("十槽中位数通过；采样不确定性", "Ten-slot median pass / sampling uncertainty")
        : recommendationLabel("本批未见门限跨越", "No within-run gate crossing"),
  })[label] || "";
}
function recommendationNeedsReview(item, design) {
  if (recommendationTenSlotConfirmation(item, design)) return false;
  if (recommendationHasTenSlotMarker(item, design)) return true;
  return ["SAMPLING_SENSITIVE_REVIEW", "REPEAT_EVIDENCE_INCOMPLETE", "SAMPLING_DIAGNOSTIC_UNKNOWN"]
    .includes(recommendationRiskCode(item, design));
}
function recommendationDesign(item) {
  return (state.result?.designs || []).find((design) => design.id === item?.design_id);
}
function recommendationQualified(item) {
  const design = recommendationDesign(item);
  const protocol = item?.execution_protocol || state.result?.execution_protocol
    || state.result?.joint_v4_ml?.execution_protocol || state.result?.request?.execution_protocol;
  const protocolName = typeof protocol === "string" ? protocol : protocol?.name;
  const siteFirst = protocolName === "ordinary_site_first_v2";
  const audit = recommendationRepeatAudit(item, design);
  const confirmedTenSlot = recommendationTenSlotConfirmation(item, design);
  return (!siteFirst || (item?.formal_eligible === true
      && audit?.sampling_decision_diagnostic_status === "COMPLETE"
      && (audit?.sampling_decision_sensitive === false || confirmedTenSlot)
      && audit?.cross_execution_conflict !== true))
    && item?.formal_eligible !== false
    && item?.display_eligible !== false
    && item?.display_status !== "DISPLAY_REVIEW"
    && !(Array.isArray(item?.display_exclusion_reasons) && item.display_exclusion_reasons.length)
    && item?.status === "PASSED"
    && !recommendationNeedsReview(item, design);
}
function recommendationRepeatSummary(item, design) {
  const audit = recommendationRepeatAudit(item, design);
  if (!audit) return "";
  const expectedCount = recommendationRepeatCount(item);
  const values = Array.isArray(audit.foldx_repeat_values_kcal_mol)
    ? audit.foldx_repeat_values_kcal_mol.filter((value) => typeof value === "number" && Number.isFinite(value)) : [];
  const count = Number.isInteger(audit.foldx_repeat_count) ? audit.foldx_repeat_count : values.length;
  const min = Number.isFinite(audit.foldx_repeat_min_kcal_mol) ? audit.foldx_repeat_min_kcal_mol : values.length ? Math.min(...values) : null;
  const max = Number.isFinite(audit.foldx_repeat_max_kcal_mol) ? audit.foldx_repeat_max_kcal_mol : values.length ? Math.max(...values) : null;
  const gate = Number.isFinite(audit.foldx_repeat_gate_kcal_mol) ? audit.foldx_repeat_gate_kcal_mol : null;
  const above = Number.isInteger(audit.foldx_repeat_above_gate_count) ? audit.foldx_repeat_above_gate_count : null;
  const below = Number.isInteger(audit.foldx_repeat_at_or_below_gate_count) ? audit.foldx_repeat_at_or_below_gate_count : null;
  const range = min == null || max == null ? recommendationLabel("范围未计算", "range unavailable") : `${fmt(min, 3)}–${fmt(max, 3)}`;
  const gateText = gate == null ? recommendationLabel("门限未记录", "gate not recorded") : `${recommendationLabel("门限", "gate")} ≤ ${fmt(gate, 3)}`;
  const counts = above == null || below == null ? "" : ` · ${recommendationLabel("门上", "above gate")} ${above}/${expectedCount} · ${recommendationLabel("门内", "at/below gate")} ${below}/${expectedCount}`;
  const slots = Array.isArray(audit.foldx_repeat_slots) ? audit.foldx_repeat_slots : [];
  const alignedSlotValues = Array.isArray(audit.foldx_repeat_slot_values_kcal_mol)
    || slots.length === values.length;
  const perSlot = Array.isArray(audit.foldx_repeat_slot_values_kcal_mol)
    ? audit.foldx_repeat_slot_values_kcal_mol
    : values;
  const slotValues = perSlot.length ? ` · ${recommendationLabel("槽值", "slots")} ${perSlot.map((value, index) =>
    `${alignedSlotValues ? (slots[index] ?? index) : `${recommendationLabel("记录", "record")} ${index + 1}`}: ${typeof value === "number" && Number.isFinite(value) ? fmt(value, 3) : recommendationLabel("缺失", "missing")}`
  ).join(" / ")} kcal/mol` : "";
  const state = count < expectedCount ? ` · ${recommendationLabel("待补算", "computation pending")}` : "";
  const uncertainty = expectedCount === 10 ? ` · ${recommendationLabel("十槽采样不确定性", "10-slot sampling uncertainty")}` : "";
  return `${recommendationLabel("FoldX重复", "FoldX repeats")} ${count}/${expectedCount} · ${gateText} · ${recommendationLabel("范围", "range")} ${range}${counts}${slotValues}${state}${uncertainty}`;
}
function recommendationRepeatNote(item, design) {
  return recommendationRepeatSummary(item, design);
}
function recommendationDisplayStatus(item, design) {
  const risk = recommendationRiskLabelText(item, design);
  return recommendationNeedsReview(item, design) ? risk : recommendationStatus(item?.status);
}
function recommendationSelectionRecord(item, design) {
  if (!recommendationQualified(item)
    || recommendationNeedsReview(item, design)
    || item?.status !== "PASSED" || item?.formal_eligible === false) return null;
  const mutations = design?.mutations?.length ? design.mutations : (item.mutation_notations || []).map((value) => {
    const match = String(value).match(/^([A-Z])([0-9]+)([A-Z])$/);
    return match ? { wt: match[1], position: Number(match[2]), mutant: match[3] } : null;
  }).filter(Boolean);
  const id = item.candidate_id || item.design_id;
  const risk = recommendationRiskLabelText(item, design);
  const status = recommendationNeedsReview(item, design) ? risk : recommendationStatus(item.status || design?.status);
  return id && mutations.length ? { id: String(id), design_id: item.design_id || design?.id || null, label: mutations.map(notation).join(" + "), status, mutations } : null;
}
function selectionColor(index) { return recommendationColors[index % recommendationColors.length]; }
function cleanRecommendationCopy(value) {
  return String(value ?? "").replace(/，不表示已预测或最大化抗体亲和力下降。?/g, "。")
    .replace(/ · 不代表亲和力下降最大/g, "")
    .replace(/；不代表亲和力变化预测。?/g, "。")
    .replace(/；几何码不代表实测抗体响应码。?/g, "。")
    .replace(/，不能据此断言亲和力变化。?/g, "。")
    .replace(/，不能直接代表真实扰动范围或结合变化。?/g, "。")
    .replace(/；不能据此推断亲和力变化。?/g, "。")
    .replace(/；几何邻域不等于实际抗体受影响面积。?/g, "。")
    .replace(/；不能据此推断结合自由能变化。?/g, "。")
    .replace(/几何邻域不等于实际抗体受影响面积/g, "")
    .replace(/，不能推断亲和力变化。?/g, "。")
    .replace(/，不等于已证实结构失稳或结合变化/g, "")
    .replace(/，不是本次突变已造成损伤/g, "")
    .replace(/；这不等于已证实结构失稳或结合变化/g, "")
    .replace(/。化学轴是未校准的机制描述符。?/g, "。")
    .replace(/表位倾向\s*P?仅为[^。]*?概率。?/g, "表位倾向 P 按计算表面排序与局部峰值标注。")
    .replace(/epitope tendency\s*P? only [^.]*?probability\.?/gi, "Epitope tendency P uses computed surface ranking and local peaks.")
    .replace(/，不是亲和力目标或抗体结合概率/g, "")
    .replace(/；这些量不是结合自由能、亲和力或其变化的测量值。?/g, "。")
    .replace(/；这不是实验稳定性证明。?/g, "。")
    .replace(/。非加性与稳定性预测不是实验协同作用或结构破坏的证明。?/g, "。")
    .replace(/aromatic_delta/g, recommendationLabel("芳香性变化", "aromaticity change"))
    .replace(/charge_delta/g, recommendationLabel("电荷变化", "charge change"))
    .replace(/volume_delta_a3/g, recommendationLabel("体积变化", "volume change"))
    .replace(/hydropathy_delta/g, recommendationLabel("疏水性变化", "hydropathy change"))
    .replace(/donor_delta/g, recommendationLabel("氢键供体变化", "hydrogen-bond donor change"))
    .replace(/acceptor_delta/g, recommendationLabel("氢键受体变化", "hydrogen-bond acceptor change"))
    .replace(/变化变化/g, "变化")
    .replace(/保留(电荷变化|体积变化|疏水性变化|氢键供体变化|氢键受体变化|芳香性变化)较大的机制备选/g, "按$1较大的原则选择备选");
}
function selectionPanel() {
  const selections = [...(state.recommendationSelections?.values?.() || [])];
  const relations = new Map();
  selections.forEach((selection, index) => selection.mutations.forEach((mutation) => {
    const entries = relations.get(mutation.position) || [];
    entries.push({ ...mutation, label: selection.label, color: selectionColor(index) });
    relations.set(mutation.position, entries);
  }));
  const conflicts = [...relations.entries()].filter(([, entries]) => entries.length > 1).map(([position, entries]) => {
    const replacements = new Set(entries.map((entry) => `${entry.wt}>${entry.mutant}`));
    const type = replacements.size > 1 ? recommendationLabel("同位点替换冲突", "Conflicting substitutions") : recommendationLabel("位点重合", "Shared site");
    return `<li><strong>${position}</strong> · ${type}: ${entries.map((entry) => `<span class="recommendation-selection-site"><i style="background:${entry.color}"></i>${h(entry.label)} (${entry.wt}${position}→${entry.mutant})</span>`).join(" ")}</li>`;
  }).join("");
  return `<section class="recommendation-selection" aria-live="polite"><div class="recommendation-selection-head"><h3>${recommendationLabel("已选组合空间关系", "Selected combination relationships")} · ${selections.length}</h3>${selections.length ? `<button class="btn btn-mini" type="button" data-clear-recommendations>${recommendationLabel("清空", "Clear")}</button>` : ""}</div>
    ${selections.length ? `<ul class="recommendation-selection-list">${selections.map((selection, index) => `<li><span class="recommendation-selection-site"><i style="background:${selectionColor(index)}"></i><strong>${h(selection.label)}</strong><small>${h(selection.id)}</small></span><span class="tag">${h(selection.status)}</span><button class="btn btn-mini" type="button" data-remove-recommendation="${h(selection.id)}">${recommendationLabel("取消选择", "Remove")}</button></li>`).join("")}</ul>${conflicts ? `<h4>${recommendationLabel("位点关系", "Site relationships")}</h4><ul class="recommendation-selection-relations">${conflicts}</ul>` : `<p class="muted">${recommendationLabel("所选组合之间没有重合突变位点。", "Selected combinations have no shared mutation sites.")}</p>`}` : `<p class="muted">${recommendationLabel("选择推荐组合，查看各组合的空间位置、位点重合与替换冲突。", "Select recommendations to inspect spatial positions, shared sites and substitution conflicts.")}</p>`}</section>`;
}
function sectorRecommendationText(item) {
  return [recommendationStatus(item.status), ...(item.gate_details || []).map(gateExplanation)].join("；");
}

function recommendationLabel(zh, en) { return state.lang === "zh" ? zh : en; }
function surfaceGateDescription(analysis, lang = "zh") {
  if (!analysis) return lang === "zh" ? "载入任务后显示该次运行的实际表面门槛。" : "The run-specific surface gate appears after loading a task.";
  const gate = analysis?.joint_v4_ml?.surface_gate;
  const methods = analysis?.methods || {};
  const thresholds = analysis?.parameters?.scientific_thresholds || analysis?.parameters || {};
  if (gate?.policy === "native_total_rsasa_only") {
    const threshold = gate.native_total_rsasa_threshold;
    return lang === "zh"
      ? `本次表面门：原生 total rSASA ≥ ${threshold ?? "未记录"}；侧链 SASA 为描述证据，不是表面门。坐标与突变资格另行检查。`
      : `Surface gate for this run: native total rSASA ≥ ${threshold ?? "not recorded"}; side-chain SASA is descriptive, not a gate. Coordinate and mutation eligibility are checked separately.`;
  }
  if (methods.sidechain_sasa_role === "descriptive_only_not_surface_gate") {
    const threshold = methods.rsasa_threshold ?? thresholds.rsasa_threshold;
    return lang === "zh"
      ? `表面暴露按 total rSASA${threshold == null ? "" : ` ≥ ${threshold}`} 判断；侧链 SASA 仅作描述。坐标资格和突变审核独立。`
      : `Surface accessibility uses total rSASA${threshold == null ? "" : ` ≥ ${threshold}`}; side-chain SASA is descriptive. Coordinate and mutation eligibility are separate.`;
  }
  const rsasa = methods.rsasa_threshold ?? thresholds.rsasa_threshold;
  const sidechain = methods.sidechain_sasa_threshold_A2 ?? thresholds.sidechain_sasa_threshold_A2;
  return lang === "zh"
    ? `此历史运行使用 rSASA${rsasa == null ? "" : ` ≥ ${rsasa}`} 与侧链 SASA${sidechain == null ? "" : ` ≥ ${sidechain} Å²`} 联合门槛。`
    : `This historical run used a combined rSASA${rsasa == null ? "" : ` ≥ ${rsasa}`} and side-chain SASA${sidechain == null ? "" : ` ≥ ${sidechain} Å²`} gate.`;
}
function isRegionalMechanism(item) {
  return item.recommendation_policy === "regional_mechanism_pareto_v1"
    || item.impact_evidence?.policy === "regional_mechanism_pareto_v1";
}
function recommendationStatus(status) {
  return ({ PASSED: recommendationLabel("通过计算门槛", "Passed computational gates"),
    NEAR_THRESHOLD: recommendationLabel("待审查 · 未通过", "Review required · not qualified"),
    EVIDENCE_LIMITED: recommendationLabel("计算通过 · 扰动证据不全", "Computational gates passed · perturbation evidence incomplete"),
    NO_COMBINATION: recommendationLabel("没有安全的专属组合", "No safe exclusive combination") })[status]
    || recommendationLabel("资格尚未确认", "Qualification unknown");
}
function recommendationNumber(value, unit = "", percent = false) {
  return typeof value === "number" && Number.isFinite(value)
    ? `${percent ? fmtPct(value) : fmt(value, 3)}${unit ? ` ${unit}` : ""}`
    : recommendationLabel("未计算／缺少证据", "Not computed / missing evidence");
}
function recommendationRegion(id) {
  const index = state.result.joint_v4_ml?.sectors?.findIndex((sector) => sector.sector_id === id);
  return `${recommendationLabel("分区", "Region")} ${index >= 0 ? index + 1 : Number(String(id).match(/\d+/)?.[0]) || "—"}`;
}
function recommendationCode(code) {
  const labels = {
    mechanism_pareto: ["机制扰动代表", "Mechanism perturbation representative"],
    charge_extreme: ["电荷变化较大", "Larger charge change"],
    volume_extreme: ["体积变化较大", "Larger volume change"],
    hydropathy_extreme: ["疏水性变化较大", "Larger hydropathy change"],
    hbond_extreme: ["供受体变化较大", "Larger donor/acceptor change"],
    sector_coverage: ["本区专属首选", "Exclusive region primary"], lower_structural_risk: ["更低预测结构风险", "Lower predicted structural risk"],
    epitope_support: ["更强表位先验支持", "Stronger epitope prior"], cross_region_coverage: ["跨区覆盖补充", "Additional cross-region coverage"],
    combination_spurs: ["SPURS组合稳定性", "SPURS combination stability"], combination_foldx: ["FoldX组合稳定性", "FoldX combination stability"],
    model_agreement: ["两种模型的预测差异", "Difference between model predictions"], higher_order_gain: ["高阶组合新增收益", "Higher-order additional gain"],
    single_evidence: ["单点证据完整性", "Single-site evidence completeness"],
    formal_eligibility: ["完整计算资格", "Full computational qualification"],
    model_evidence_complete: ["组合双模型证据完整性", "Complete combination evidence from both models"],
    candidate_connectivity: ["组合连通性", "Combination connectivity"],
    protected_residue_policy: ["保护残基规则", "Protected-residue policy"],
    structural_confidence: ["结构位置可靠性", "Structural coordinate confidence"],
    SASA_METHOD_SETTINGS_UNVERIFIED: ["表面积计算参数尚未完整核对；暴露度存在方法不确定性，不等于已证实结构不稳定", "Surface-area settings have not been fully verified; accessibility has methodological uncertainty, not proven structural instability"],
    SASA_CONTEXT_SENSITIVE: ["表面积对结构计算范围敏感，应结合单链或复合物环境复核", "Surface area depends on the structural calculation context; review the single-chain or complex environment"],
    SASA_STRUCTURE_CONTEXT_CLASSIFICATION_CHANGED: ["完整输入结构与ECD计算范围的表面暴露分类不一致，需要复核结构环境；这不等于已证实结构不稳定", "Surface-accessibility classification differs between the full input structure and the ECD context; review the structural environment. This is not evidence of instability"],
    sasa_structure_context_classification_changed: ["结构计算范围的暴露分类一致性", "Accessibility consistency between structural contexts"],
    "beta-strand": ["记录将该位点标注为β链，需留意局部二级结构；此标注不是结构失稳证据", "This site is annotated as a beta strand; review local secondary structure. This annotation is not evidence of instability"],
    "salt-bridge": ["输入结构提示该位点参与盐桥，电荷替换可能改变局部联系；这不等于已证实结构失稳或结合变化", "The input structure annotates a salt bridge at this site; a charge substitution may alter local contacts. This does not establish instability or altered binding"],
    "engineered-neighborhood": ["该位点邻近输入结构中的工程化差异，需要与目标天然序列对照；不是本次突变已造成损伤", "This site is near engineered differences in the input structure; compare with the target native sequence. This is not evidence that the designed mutation caused damage"],
    COMBINATION_REJECT: ["组合未满足当前计算推荐门槛；具体超限或缺失证据见下方", "Combination does not meet the current computational recommendation gates; see specific excess or missing evidence below"],
    COMBINATION_SPURS_MISSING: ["缺少SPURS组合稳定性证据", "Missing SPURS combination evidence"],
    COMBINATION_FOLDX_MISSING: ["缺少FoldX组合评分，待补算", "Missing FoldX combination score; computation pending"],
    FOLDX_REPEAT_EVIDENCE_INCOMPLETE: ["FoldX组合评分重复次数不足，待补算", "Insufficient FoldX combination repeats; computation pending"],
    COMBINATION_SPURS_ABOVE_THRESHOLD: ["SPURS预测值超过当前门槛", "SPURS prediction exceeds the threshold"],
    COMBINATION_FOLDX_ABOVE_THRESHOLD: ["FoldX预测值超过当前门槛", "FoldX prediction exceeds the threshold"],
    MODEL_DISAGREEMENT: ["两种模型预测差异超过允许范围", "Model predictions disagree beyond the allowed range"],
    HIGHER_ORDER_REQUIRES_FOLDX: ["高阶组合仍需FoldX确认", "Higher-order combination needs FoldX confirmation"],
    INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN: ["高阶组合新增收益未满足最低要求", "Higher-order additional gain does not meet the minimum requirement"],
    HIGHER_ORDER_COVERAGE_BASELINE_MISSING: ["缺少同区低阶几何对照，无法确认新增覆盖收益；缺失不能当作零", "Missing a lower-order geometric comparison in the same region; additional coverage is unknown, not zero"],
    SINGLE_EVIDENCE_MISSING_OR_EXCLUDED: ["单点证据缺失或不符合准入条件", "Single-site evidence missing or excluded"],
    NO_CANDIDATE_WITH_MUTATION_IN_SECTOR: ["没有位点位于本区的安全候选", "No safe candidate mutates this region"],
    NO_QUALIFIED_PRINCIPLE_CANDIDATE: ["尚无通过本次计算门槛且符合此原则的已确认候选", "No confirmed candidate passes the gates for this principle"],
    NO_SAFE_EXCLUSIVE_CANDIDATE: ["本区没有安全的专属双点或三点候选", "No safe exclusive double or triple"],
    AUTO_EXCLUDED_MUTANT: ["替换氨基酸属于自动排除类型", "Substitution is automatically excluded"],
    AUTO_PROTECTED_NATIVE_RESIDUE: ["涉及需保护的原始残基", "Mutates a protected native residue"],
    DISCONNECTED_CANDIDATE: ["组合不满足连通要求", "Combination fails connectivity requirements"],
    FORMAL_ELIGIBILITY_NOT_ESTABLISHED: ["计算资格证据尚未建立", "Computational qualification not established"],
  };
  return labels[code] ? recommendationLabel(...labels[code])
    : recommendationLabel("尚未解释的审查项，请查看原始证据", "Unmapped review item; see original evidence");
}
function gateExplanation(gate) {
  const description = gate.description || gate.message;
  if (description && !description.includes("\ufffd")) return description;
  if (gate.reason) {
    if (gate.measurement_basis === "mean_regional_area_coverage_9_12_15" && Number.isFinite(gate.observed) && Number.isFinite(gate.threshold)) {
      return recommendationLabel("三点组合相对同区最佳低阶几何对照的平均估计覆盖增益不足", "Mean estimated coverage gain versus the best lower-order geometric comparison in the same region is insufficient")
        + `；${recommendationLabel("9/12/15 Å三尺度平均增加", "mean increase at 9/12/15 Å")} ${fmt(gate.observed * 100, 1)} ${recommendationLabel("个百分点", "percentage points")}，${recommendationLabel("最低要求", "required minimum")} ${fmt(gate.threshold * 100, 1)} ${recommendationLabel("个百分点", "percentage points")}`;
    }
    const reason = gate.gate === "higher_order_gain" && gate.measurement_basis === "information_gain" ? recommendationLabel("三点组合相对双点的信息增益不足", "Triple versus double information gain is insufficient") : gate.gate === "higher_order_gain" && gate.measurement_basis === "coverage_gain" ? recommendationLabel("三点组合相对双点的覆盖增益不足", "Triple versus double coverage gain is insufficient") : /^[A-Z][A-Z0-9_]+$/.test(gate.reason) ? recommendationCode(gate.reason) : gate.reason;
    if (typeof gate.observed === "number" && typeof gate.threshold === "number") {
      const minimum = gate.comparison === "minimum" || gate.gate === "higher_order_gain";
      const unit = gate.unit || (["combination_spurs", "combination_foldx", "model_agreement"].includes(gate.gate) ? "kcal/mol" : "");
      return `${reason}；${recommendationLabel("实际值", "observed")} ${recommendationNumber(gate.observed, unit)}，${recommendationLabel(minimum ? "最低要求" : "上限", minimum ? "required minimum" : "maximum")} ${recommendationNumber(gate.threshold, unit)}；${recommendationLabel(minimum ? "距最低要求还差" : "超出上限", minimum ? "below minimum by" : "exceeds maximum by")} ${recommendationNumber(gate.excess, unit)}`;
    }
    return reason;
  }
  const ruleConclusions = {
    candidate_connectivity: ["突变组合满足当前连通性规则", "Mutation combination satisfies the current connectivity rule"],
    protected_residue_policy: ["未触发保护残基或禁止替换的硬排除规则", "No protected-residue or prohibited-substitution hard exclusion was triggered"],
    single_evidence: ["所需单点证据齐全，且没有触发单点硬排除", "Required single-site evidence is complete without single-site hard exclusions"],
    structural_confidence: ["满足当前坐标准入规则；预测模型置信度分数未记录", "Current coordinate admission rules satisfied; predicted-model confidence score not recorded"],
    formal_eligibility: ["完整计算资格检查通过，仍需实验验证", "Full computational qualification passed; experimental validation remains necessary"],
    model_evidence_complete: ["SPURS与FoldX组合评分均已获得有限数值", "Finite SPURS and FoldX combination scores are available"],
  };
  if (gate.passed === true && ruleConclusions[gate.gate]) return recommendationLabel(...ruleConclusions[gate.gate]);
  return gate.required === false ? recommendationLabel("当前合同仅记录此项，不作为准入条件", "Recorded only; not required by this contract")
    : gate.passed === true ? recommendationLabel("满足当前准入条件", "Current qualification condition satisfied")
      : gate.evidence_missing ? recommendationLabel("缺少模型证据，不能确认此项", "Missing model evidence; cannot confirm this check")
        : recommendationLabel("不满足当前准入条件", "Current qualification condition not satisfied");
}
function recommendationTable(headers, rows) {
  return `<div class="recommendation-table"><table class="audit"><thead><tr>${headers.map((label) => `<th>${h(label)}</th>`).join("")}</tr></thead><tbody>${rows.map((values) => `<tr>${values.map((value) => `<td>${h(value)}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`;
}
function regionalMechanismSummary(item) {
  const exposure = item.impact_evidence?.exposure_min;
  const sensitivity = item.sensitivity || item.impact_evidence?.sensitivity || {};
  const sensitive = sensitivity.sensitive;
  const sensitivitySummary = regionalSensitivitySummary(sensitivity);
  const uncovered = item.local_12_uncovered_positions;
  const uncoveredText = Array.isArray(uncovered)
    ? uncovered.join(", ") || recommendationLabel("无", "None")
    : recommendationLabel("未知", "Unknown");
  return `<div class="regional-mechanism-summary">
    <p><strong>${recommendationLabel("本区强扰动工程代表", "Regional strong-perturbation engineering representative")}</strong></p>
    <p>${recommendationLabel("最小相对侧链暴露度", "Minimum relative side-chain exposure")}: ${recommendationNumber(exposure, "", true)} · ${recommendationLabel("选择敏感性", "Selection sensitivity")}: ${sensitive === true ? `${recommendationLabel("敏感", "Sensitive")}${sensitivitySummary ? `（${sensitivitySummary}）` : "：预设变体会改变代表项"}` : sensitive === false ? recommendationLabel("稳健：预设变体未改变代表项", "Stable: predefined variants did not change the representative") : recommendationLabel("未知", "Unknown")}</p>
    <p>${recommendationLabel("直接映射突变面积", "Direct mapped mutation area")}: ${recommendationNumber(item.direct_area_a2, "Å²")}</p>
    <p>${recommendationLabel("本区12 Å局部覆盖", "Regional 12 Å local coverage")}: ${recommendationNumber(item.local_12_covered_area_a2, "Å²")} · ${recommendationNumber(item.local_12_coverage_fraction, "", true)} · ${recommendationLabel("未覆盖面积", "uncovered area")}: ${recommendationNumber(item.local_12_uncovered_area_a2, "Å²")} · ${recommendationLabel("未覆盖位点", "uncovered positions")}: ${uncoveredText}</p>
    <p class="note">${recommendationLabel("表位倾向 P 按计算表面排序与局部峰值标注；化学轴记录电荷、体积、疏水性、氢键供受体及芳香性变化。", "Epitope tendency P uses computed surface ranking and local peaks; chemical axes record changes in charge, volume, hydropathy, hydrogen-bond donors and acceptors, and aromaticity.")}</p>
  </div>`;
}
const chemicalAxisDefinitions = [
  ["charge_delta", "电荷变化", "Charge change", "Δe"],
  ["volume_delta_a3", "体积变化", "Volume change", "Å³"],
  ["hydropathy_delta", "疏水性变化", "Hydropathy change", ""],
  ["donor_delta", "氢键供体变化", "H-bond donor change", ""],
  ["acceptor_delta", "氢键受体变化", "H-bond acceptor change", ""],
  ["aromatic_delta", "芳香性变化", "Aromaticity change", ""],
];
function candidateImpactSites(item) {
  return item.site_details?.length ? item.site_details
    : item.impact_evidence?.regional_impact?.mutations || item.impact_evidence?.mutations || [];
}
function candidateImpactAxes(item, sites = candidateImpactSites(item)) {
  const evidence = item.impact_evidence || {};
  const axes = { ...(evidence.regional_impact?.mean_absolute_delta || {}), ...(evidence.regional_impact?.axes || {}), ...(evidence.mean_absolute_delta || {}), ...(evidence.axes || {}) };
  const scenarioRanges = evidence.display_chemistry_ranges?.ranges?.axes || {};
  for (const [key] of chemicalAxisDefinitions) {
    if (typeof axes[key] === "number" && Number.isFinite(axes[key])) continue;
    if (scenarioRanges[key]) { axes[key] = scenarioRanges[key]; continue; }
    const values = sites.map((site) => site.chemical_descriptors?.[key]);
    if (values.length && values.every((value) => typeof value === "number" && Number.isFinite(value))) {
      axes[key] = values.reduce((sum, value) => sum + Math.abs(value), 0) / values.length;
    }
  }
  return axes;
}
function hasCandidateImpact(item, sites = candidateImpactSites(item)) {
  const evidence = item.impact_evidence || {};
  return Boolean(evidence.axes || evidence.mean_absolute_delta || evidence.regional_impact
    || evidence.signed_change_summary || evidence.local_interaction_diagnostics || evidence.exposure_min != null
    || evidence.display_chemistry_ranges
    || sites.some((site) => site.chemical_descriptors || site.sidechain_exposure));
}
function chemistryRangeValue(value, unit = "", percent = false) {
  if (!value || typeof value !== "object") return recommendationNumber(value, unit);
  const minimum = value.minimum ?? value.min;
  const maximum = value.maximum ?? value.max;
  if (minimum == null && maximum == null) return recommendationLabel("未计算", "Not computed");
  if (minimum != null && maximum != null && minimum === maximum) return recommendationNumber(minimum, unit, percent);
  return `${recommendationNumber(minimum, unit, percent)} – ${recommendationNumber(maximum, unit, percent)}`;
}
function displayChemistryRangeEvidence(item) {
  const ranges = item.impact_evidence?.display_chemistry_ranges;
  if (!ranges) return "";
  const axisRanges = ranges.ranges?.axes || {};
  const axes = chemicalAxisDefinitions.map(([key, zh, en, unit]) => [
    recommendationLabel(zh, en), chemistryRangeValue(axisRanges[key], unit),
  ]);
  const sites = ranges.ranges?.per_site_signed_deltas || {};
  const siteRows = Object.entries(sites).map(([notationText, siteAxes]) => [notationText,
    ...chemicalAxisDefinitions.map(([key]) => chemistryRangeValue(siteAxes?.[key]))]);
  const signedRanges = ranges.ranges?.signed_change_summary || {};
  const signedRows = ["charge_delta", "hydropathy_delta"].filter((axis) => signedRanges[axis]).map((axis) => {
    const value = signedRanges[axis];
    return [recommendationLabel(axis === "charge_delta" ? "电荷" : "疏水性", axis === "charge_delta" ? "Charge" : "Hydropathy"),
      chemistryRangeValue(value.net_change), chemistryRangeValue(value.absolute_sum), chemistryRangeValue(value.opposition_fraction, "", true)];
  });
  const stateLabels = { neutral: recommendationLabel("中性", "neutral"), protonated: recommendationLabel("质子化", "protonated") };
  const microstates = (ranges.microstates || []).map((site) => `${site.notation || site.position}: ${site.states.map((state) => `${stateLabels[state.state] || state.state}（电荷 ${state.charge}，供体 ${state.donor}，受体 ${state.acceptor}）`).join(" / ")}`);
  const scenarios = (ranges.scenarios || []).map((scenario) => {
    const stateSummary = (scenario.histidine_states || []).map((state) => `${state.notation || state.position} ${state.side === "wild_type" ? recommendationLabel("WT", "WT") : recommendationLabel("替换后", "mutant")} ${stateLabels[state.state] || state.state}`).join("；");
    return [scenario.scenario_id || "—", stateSummary || recommendationLabel("无组氨酸状态项", "No histidine state assignments"),
      ...chemicalAxisDefinitions.map(([key]) => recommendationNumber(scenario.mean_absolute_delta?.[key]))];
  });
  const interactionClass = (record) => {
    const classification = record?.classification;
    const labels = {
      CONSISTENT_TRUE: ["各情景均为近邻异向", "Opposed in every scenario"],
      CONSISTENT_FALSE: ["各情景均未见近邻异向", "No nearby opposition in any scenario"],
      STATE_DEPENDENT: ["是否近邻异向随微态变化", "Nearby opposition depends on microstate"],
    };
    return labels[classification] ? recommendationLabel(...labels[classification]) : recommendationLabel("未计算", "Not computed");
  };
  const interactions = ranges.ranges?.local_interactions || [];
  const interactionRows = interactions.map((pair) => [
    (pair.positions || []).join("–"), recommendationNumber(pair.wt_ca_distance_a, "Å"),
    interactionClass(pair.nearby_opposed_charge), interactionClass(pair.nearby_opposed_hydropathy),
  ]);
  return `<section class="chemistry-scenario-ranges"><h4>${recommendationLabel("组氨酸微态情景范围", "Histidine microstate scenario ranges")}</h4>
    <p class="note">${recommendationLabel("含组氨酸项按中性／质子化状态给出范围；各微态独立枚举，不指定pH或实际质子化状态。原始未知值不进入排序。", "Histidine-containing candidates show neutral/protonated scenarios; microstates are enumerated independently, without assigning pH or actual protonation. Original unknowns remain outside ranking.")}</p>
    ${microstates.length ? `<p>${microstates.join("；")}</p>` : ""}
    ${recommendationTable([recommendationLabel("化学轴", "Chemical axis"), recommendationLabel("情景范围（平均绝对变化）", "Scenario range (mean absolute change)")], axes)}
    ${siteRows.length ? `<h5>${recommendationLabel("逐位带符号变化范围", "Per-site signed-change ranges")}</h5>${recommendationTable([recommendationLabel("位点", "Site"), ...chemicalAxisDefinitions.map(([, zh, en]) => recommendationLabel(zh, en))], siteRows)}` : ""}
    ${signedRows.length ? `<h5>${recommendationLabel("signed变化汇总范围", "Signed-change summary ranges")}</h5>${recommendationTable([recommendationLabel("变化轴", "Axis"), recommendationLabel("净变化范围", "Net-change range"), recommendationLabel("绝对变化总和范围", "Absolute-sum range"), recommendationLabel("反向比例范围", "Opposition-fraction range")], signedRows)}` : ""}
    ${scenarios.length ? `<h5>${recommendationLabel("逐情景计算值", "Per-scenario computed values")}</h5>${recommendationTable([recommendationLabel("情景", "Scenario"), recommendationLabel("组氨酸状态", "Histidine states"), ...chemicalAxisDefinitions.map(([, zh, en]) => recommendationLabel(zh, en))], scenarios)}` : ""}
    ${interactionRows.length ? `<h5>${recommendationLabel("近邻方向分类范围", "Nearby direction-class ranges")}</h5>${recommendationTable([recommendationLabel("位点对", "Site pair"), "WT Cα", recommendationLabel("电荷", "Charge"), recommendationLabel("疏水性", "Hydropathy")], interactionRows)}` : ""}</section>`;
}
function regionalImpactEvidence(item, sites) {
  const evidence = item.impact_evidence || {};
  const axes = candidateImpactAxes(item, sites);
  const rows = [
    ...chemicalAxisDefinitions.map(([key, zh, en, unit]) => [recommendationLabel(zh, en), chemistryRangeValue(axes[key], unit)]),
    ...(evidence.exposure_min != null || evidence.regional_impact?.exposure_min != null
      ? [[recommendationLabel("最小相对侧链暴露度", "Minimum relative side-chain exposure"), recommendationNumber(evidence.exposure_min ?? evidence.regional_impact?.exposure_min, "", true)]] : []),
    ...(evidence.exposure_mean != null || evidence.regional_impact?.exposure_mean != null
      ? [[recommendationLabel("平均相对侧链暴露度", "Mean relative side-chain exposure"), recommendationNumber(evidence.exposure_mean ?? evidence.regional_impact?.exposure_mean, "", true)]] : []),
  ];
  const sensitivity = item.sensitivity || item.impact_evidence?.sensitivity || {};
  const sensitivitySummary = regionalSensitivitySummary(sensitivity, true);
  if (sensitivitySummary) rows.push([recommendationLabel("哪些预设变体改变了代表项", "Predefined variants that changed the representative"), sensitivitySummary]);
  if (Number.isInteger(sensitivity.evaluated_candidate_count)) rows.push([
    recommendationLabel("敏感性比较范围", "Sensitivity comparison scope"),
    recommendationLabel(`当前同阶、同缺失模式的 ${sensitivity.evaluated_candidate_count} 个候选${sensitivity.evaluated_candidate_count < 2 ? "；比较对象不足" : ""}。`,
      `${sensitivity.evaluated_candidate_count} current candidates of the same order and missing-evidence pattern${sensitivity.evaluated_candidate_count < 2 ? "; insufficient comparators" : ""}.`)]);
  const groupNote = item.impact_evidence?.ranking_version === "mechanism_channels_v2"
    ? `<p class="note">${recommendationLabel("当前机制折中：电荷／亲疏水／氢键组、体积／芳香组，各取最强轴效用；暴露度单独比较，组间补短板。此分组是工程偏好。", "Current mechanism compromise: maximum axis utility within charge/polarity/H-bond and size/aromaticity groups; exposure is separate, with leximin across groups. These groups express an engineering preference.")}</p>` : "";
  return groupNote + recommendationTable([recommendationLabel("平均绝对化学变化", "Mean absolute chemical change"), recommendationLabel("变化幅度（非带符号）", "Change magnitude (unsigned)")], rows) + displayChemistryRangeEvidence(item);
}
function regionalInteractionEvidence(item) {
  const evidence = item.impact_evidence || {};
  const regional = evidence.regional_impact || {};
  const signed = evidence.signed_change_summary || regional.signed_change_summary;
  const local = evidence.local_interaction_diagnostics || regional.local_interaction_diagnostics;
  const audit = item.resource_audit?.single_threshold_audit;
  if (!signed && !local && !audit) return "";
  const signedRanges = evidence.display_chemistry_ranges?.ranges?.signed_change_summary || {};
  const rows = signed ? [[recommendationLabel("电荷", "Charge"), signed.charge_delta, signedRanges.charge_delta],
    [recommendationLabel("疏水性", "Hydropathy"), signed.hydropathy_delta, signedRanges.hydropathy_delta]].map(([label, value, range]) => [label,
      chemistryRangeValue(value?.net_change ?? range?.net_change), chemistryRangeValue(value?.absolute_sum ?? range?.absolute_sum),
      chemistryRangeValue(value?.opposition_fraction ?? range?.opposition_fraction, "", true)]) : [];
  const rangePairs = evidence.display_chemistry_ranges?.ranges?.local_interactions || [];
  const direction = (pair, axis) => {
    const opposed = pair.opposed_changes?.[axis];
    const nearby = axis === "charge_delta" ? pair.nearby_opposed_charge : pair.nearby_opposed_hydropathy;
    if (typeof opposed === "boolean") return recommendationLabel(opposed ? (nearby === true ? "方向相反 · 10 Å内" : nearby === false ? "方向相反 · 超出10 Å" : "方向相反") : "方向相同", opposed ? (nearby === true ? "Opposed · within 10 Å" : nearby === false ? "Opposed · beyond 10 Å" : "Opposed") : "Same direction");
    if (nearby === true) return recommendationLabel("近邻异向", "Nearby and opposed");
    if (nearby === false) return recommendationLabel("未标记近邻异向", "No nearby opposition flagged");
    const positions = pair.positions || [pair.left, pair.right];
    const range = rangePairs.find((candidate) => (candidate.positions || []).join("|") === positions.join("|"));
    const classification = axis === "charge_delta" ? range?.nearby_opposed_charge : range?.nearby_opposed_hydropathy;
    if (classification) return interactionClass(classification);
    return recommendationLabel("未计算", "Not computed");
  };
  const interactionClass = (record) => ({
    CONSISTENT_TRUE: recommendationLabel("各情景均为近邻异向", "Opposed in every scenario"),
    CONSISTENT_FALSE: recommendationLabel("各情景均未见近邻异向", "No nearby opposition in any scenario"),
    STATE_DEPENDENT: recommendationLabel("是否近邻异向随微态变化", "Nearby opposition depends on microstate"),
  })[record?.classification];
  const distanceClass = (pair) => typeof pair.wt_ca_distance_a !== "number" ? recommendationLabel("未计算", "Not computed")
    : pair.ca_below_6_a ? recommendationLabel("小于6 Å", "Below 6 Å")
      : pair.ca_above_25_a ? recommendationLabel("大于25 Å", "Above 25 Å") : recommendationLabel("6–25 Å", "6–25 Å");
  const pairs = (local?.pairs || []).map(pair => [(pair.positions || [pair.left, pair.right]).join("–"),
    recommendationNumber(pair.wt_ca_distance_a, "Å"), direction(pair, "charge_delta"), direction(pair, "hydropathy_delta"), distanceClass(pair)]);
  return `<h4>${recommendationLabel("方向与近邻诊断（不参与排序）", "Signed changes and proximity diagnostics (not ranked)")}</h4>
    ${rows.length ? recommendationTable([recommendationLabel("变化轴", "Axis"), recommendationLabel("带符号净变化", "Signed net change"), recommendationLabel("绝对变化总和", "Absolute change sum"), recommendationLabel("代数反向比例", "Algebraic opposition fraction")], rows) : ""}
    ${signed ? `<p>${recommendationLabel("F/W/Y替换数／新引入数", "F/W/Y mutants / introductions")}: ${recommendationNumber(evidence.aromatic_mutant_count ?? regional.aromatic_mutant_count)} / ${recommendationNumber(evidence.introduced_aromatic_count ?? regional.introduced_aromatic_count)}</p>` : ""}
    ${pairs.length ? recommendationTable([recommendationLabel("位点对", "Site pair"), "WT Cα", recommendationLabel("电荷变化方向及距离分类", "Charge direction and distance class"), recommendationLabel("疏水变化方向及距离分类", "Hydropathy direction and distance class"), recommendationLabel("Cα距离分类", "Cα distance class")], pairs) : ""}
    ${local ? `<p class="note">${recommendationLabel("距离按 WT Cα 计算。", "Distances use WT Cα coordinates.")}</p>` : ""}
    ${audit?.review_band_count != null ? `<h4>${recommendationLabel("单点稳定性门槛边界审计", "Single-mutation stability threshold audit")}</h4><p>SPURS (${recommendationNumber(audit.lower_exclusive_kcal_mol)}, ${recommendationNumber(audit.upper_inclusive_kcal_mol)}] kcal/mol: ${recommendationNumber(audit.review_band_count)} · ${recommendationLabel("其余准入条件合格", "Otherwise eligible")}: ${recommendationNumber(audit.otherwise_eligible_count)} · RSASA &gt; ${recommendationNumber(audit.high_exposure_cutoff, "", true)}: ${recommendationNumber(audit.high_rsasa_count)} · ${recommendationLabel("相对侧链暴露度", "Relative side-chain exposure")} &gt; ${recommendationNumber(audit.high_exposure_cutoff, "", true)}: ${recommendationNumber(audit.high_sidechain_exposure_count)}</p><p class="note">${recommendationLabel("统计单位是替换方案；两种暴露度分别计数。此审计不放宽门槛，也不将REVIEW方案纳入正式推荐。", "Counts refer to substitutions, with the two exposure measures reported separately. This audit does not relax gates or admit REVIEW substitutions to formal recommendations.")}</p>` : ""}`;
}
function regionalSensitivityName(name) {
  const labels = {
    ordinal: ["基线：序数归一化", "Baseline: ordinal normalization"],
    minmax: ["最小–最大归一化", "Min–max normalization"],
    mean_exposure: ["平均侧链暴露度", "Mean side-chain exposure"],
    without_aromaticity: ["去除芳香性轴", "Omit aromaticity axis"],
    exposure: ["去除暴露度维度", "Omit exposure dimension"],
    charge: ["去除电荷维度", "Omit charge dimension"],
    volume: ["去除体积维度", "Omit volume dimension"],
    hydropathy_aromatic: ["去除疏水性／芳香性维度", "Omit hydropathy/aromaticity dimensions"],
    donor_acceptor: ["去除氢键供体／受体维度", "Omit donor/acceptor dimensions"],
    charge_polarity: ["去除电荷／亲疏水／氢键机制组", "Omit charge/polarity/H-bond group"],
    size_aromatic: ["去除体积／芳香机制组", "Omit size/aromaticity group"],
  };
  const key = String(name).replace(/^leave_one_out:/, "");
  return labels[key] ? recommendationLabel(...labels[key]) : recommendationLabel("预设敏感性变体", "Predefined sensitivity variant");
}
function regionalSensitivitySummary(sensitivity, detailed = false) {
  const winners = sensitivity.selection_winners || {};
  const variants = Array.isArray(sensitivity.variants) ? sensitivity.variants : Object.entries(winners).map(([name, winner_candidate_id]) => ({ name, winner_candidate_id }));
  const baseline = winners.ordinal ?? variants.find((variant) => variant.name === "ordinal")?.winner_candidate_id;
  const changed = variants.filter((variant) => variant.name !== "ordinal" && variant.winner_candidate_id != null
    && (baseline == null || variant.winner_candidate_id !== baseline));
  if (!changed.length) return "";
  return changed.map((variant) => `${regionalSensitivityName(variant.name)}${detailed && variant.winner_candidate_id ? ` → ${recommendationLabel("代表候选", "representative candidate")} ${variant.winner_candidate_id}` : ""}`).join("；");
}
function regionalMissingImpactNote(item) {
  const evidence = item.impact_evidence || {};
  const missing = evidence.missing_axes || evidence.missing_mask || evidence.regional_impact?.missing_axes || [];
  const labels = {
    charge_delta: ["电荷变化", "charge change"], volume_delta_a3: ["体积变化", "volume change"],
    hydropathy_delta: ["疏水性变化", "hydropathy change"], donor_delta: ["氢键供体变化", "H-bond donor change"],
    acceptor_delta: ["氢键受体变化", "H-bond acceptor change"], aromatic_delta: ["芳香性变化", "aromaticity change"],
  };
  const names = missing.map((key) => labels[key] ? recommendationLabel(...labels[key]) : null).filter(Boolean);
  return names.length
    ? recommendationLabel(`稳定性模型门槛已通过；机制扰动证据不全，缺少：${names.join("、")}。`, `Stability model gates passed; perturbation evidence is incomplete. Missing: ${names.join(", ")}.`)
    : recommendationLabel("稳定性模型门槛已通过；机制扰动证据未完整记录。", "Stability model gates passed; mechanism perturbation evidence is incomplete.");
}
function recommendationCoverage(item, detailed = false) {
  const support = item.regional_geometry_support?.scenarios || {};
  return recommendationTable([
    recommendationLabel("表面邻域距离", "Surface neighborhood"), recommendationLabel("本区估计覆盖", "Estimated region coverage"),
    recommendationLabel("未覆盖残基", "Uncovered residues"), ...(detailed ? [recommendationLabel("本区／区外面积", "Region / outside area")] : [])
  ], [9, 12, 15].map((radius) => {
    const values = support[`geodesic_${radius}`] || {};
    const missing = values.uncovered_positions ?? values.target_uncovered_positions;
    return [`${radius} Å`, recommendationNumber(values.coverage_fraction ?? values.target_coverage_fraction, "", true),
      Array.isArray(missing) ? (detailed ? (missing.join(", ") || recommendationLabel("无", "None")) : String(missing.length)) : recommendationLabel("未知", "Unknown"),
      ...(detailed ? [`${recommendationNumber(values.covered_area_a2 ?? values.target_covered_area_a2, "Å²")} / ${recommendationNumber(values.outside_covered_area_a2, "Å²")}`] : [])];
  }));
}
function recommendationDetails(item, design) {
  const mechanism = isRegionalMechanism(item);
  const detail = candidateImpactSites(item);
  const sites = detail.length ? detail : (design?.mutations || []).map((mutation) => ({ ...mutation, notation: notation(mutation) }));
  const showImpact = mechanism || hasCandidateImpact(item, sites);
  const showSiteImpact = mechanism || sites.some((site) => site.chemical_descriptors || site.sidechain_exposure);
  const rows = sites.map((site) => {
    const position = site.position ?? site.uniprot_pos;
    const residue = (state.result.residues || []).find((row) => row.uniprot_pos === position) || {};
    const assessment = site.risk_assessment;
    const riskText = assessment ? [assessment.conclusion, ...(assessment.notes || []).map((note) => note.explanation?.includes("尚未解释") ? recommendationCode(note.code) : note.explanation)].filter(Boolean).join("；") : [...new Set([...(site.risk_notes || site.hard_exclusions || []), ...(site.warnings || []), ...(site.structural_risk_flags || [])])].map((value) => /[\u4e00-\u9fff]/.test(value) && !value.includes("\ufffd") ? value : recommendationCode(value)).join("；") || recommendationLabel("未记录额外提示；不等于排除全部风险", "No additional note recorded; this does not rule out all risks");
    return [site.notation || `${site.native_aa || site.wildtype || site.wt || residue.one_letter || ""}${position ?? "?"}${site.mutant_aa || site.mutant || ""}`,
      recommendationRegion(site.sector_id || item.sector_id), recommendationNumber(site.rsasa ?? residue.rsasa),
      recommendationNumber(site.single_spurs_ddg_kcal_mol ?? site.spurs_ddg_kcal_mol ?? site.spurs_single_ddg_kcal_mol, "kcal/mol"),
      riskText];
  });
  const checks = item.gate_checks?.length ? item.gate_checks : item.gate_details || [];
  const gates = checks.map((gate) => {
    const rule = (gate.observed == null || typeof gate.observed === "boolean") && gate.threshold == null && !gate.evidence_missing;
    const notApplicable = recommendationLabel("不适用：规则检查，无数值门槛", "Not applicable: rule check without a numeric threshold");
    const proof = gate.evaluated_evidence || {};
    let ruleResult = gate.passed === true ? recommendationLabel("规则检查通过", "Rule check passed") : gate.passed === false ? recommendationLabel("规则检查未通过", "Rule check failed") : recommendationLabel("规则结论未知", "Rule outcome unknown");
    if (proof.mapped_sites != null && proof.mutation_count != null) ruleResult = `${proof.mapped_sites}/${proof.mutation_count} ${recommendationLabel("映射匹配", "sites mapped")}；${proof.coordinate_available_sites ?? "—"}/${proof.mutation_count} ${recommendationLabel("坐标可用", "coordinates available")}`;
    else if (proof.single_scores_present != null && proof.mutation_count != null) ruleResult = `${proof.single_scores_present}/${proof.mutation_count} ${recommendationLabel("已有单点评分", "single scores available")}；${recommendationLabel("硬排除项", "hard exclusions")}: ${(proof.hard_exclusions || []).length}`;
    else if (typeof proof.connected === "boolean") ruleResult = proof.connected ? recommendationLabel("组合连通", "Combination connected") : recommendationLabel("组合不连通", "Combination disconnected");
    else if (proof.native_residues && proof.mutant_residues) ruleResult = `${recommendationLabel("原始残基", "Native residues")}: ${proof.native_residues.join("/")}；${recommendationLabel("替换残基", "Substitutions")}: ${proof.mutant_residues.join("/")}`;
    return [gate.gate === "higher_order_gain" && gate.measurement_basis === "information_gain" ? recommendationLabel("三点组合相对双点的信息增益", "Triple versus double information gain") : gate.gate === "higher_order_gain" && gate.measurement_basis === "coverage_gain" ? recommendationLabel("三点组合相对双点的覆盖增益", "Triple versus double coverage gain") : recommendationCode(gate.gate),
    gate.required === false ? recommendationLabel("仅记录，不作为准入门槛", "Recorded, not a qualification gate") : gate.evidence_missing ? recommendationLabel("缺少证据", "Missing evidence") : gate.passed === true ? recommendationLabel("满足", "Satisfied") : gate.passed === false || item.gate_details?.includes(gate) ? recommendationLabel("未满足", "Not satisfied") : recommendationLabel("未知", "Unknown"),
    rule ? ruleResult : typeof gate.observed === "boolean" ? (gate.observed ? recommendationLabel("是", "Yes") : recommendationLabel("否", "No")) : recommendationNumber(gate.observed, gate.unit || ""),
    rule ? notApplicable : `${gate.comparison === "minimum" ? "≥ " : "≤ "}${recommendationNumber(gate.threshold, gate.unit || "")}`,
    rule ? notApplicable : recommendationNumber(gate.excess, gate.unit || ""), gateExplanation(gate)];
  });
  const spatial = item.regional_geometry_support?.scenarios || {};
  const geometry = item.geometry_support?.area_basis === "eligible_surface_area_a2" ? item.geometry_support.scenarios || {} : {};
  const spatialRows = [9, 12, 15].map((radius) => {
    const value = spatial[`geodesic_${radius}`] || {};
    const full = geometry[`geodesic_${radius}`] || {};
    const marginal = value.marginal_area_a2 ?? full.marginal_area_a2;
    const leave = value.leave_one_out_gain_a2 || full.leave_one_out_gain_a2 || {};
    return [`${radius} Å`, recommendationNumber(marginal, "Å²"), recommendationNumber(value.mean_pair_overlap_fraction ?? full.mean_pair_overlap_fraction, "", true),
      recommendationNumber(value.max_pair_overlap_fraction ?? full.max_pair_overlap_fraction, "", true),
      Object.entries(leave).map(([position, gain]) => `${position}: ${recommendationNumber(gain, "Å²")}`).join("；") || recommendationLabel("未计算", "Not computed")];
  });
  const sequence = design?.mutant_ecd_sequence || item.mutant_ecd_sequence;
  const coding = design?.mutant_coding_sequence || item.mutant_coding_sequence;
  const provenance = item.provenance || {};
  const structure = provenance.structure_evidence || item.structure_evidence || state.result.structure_evidence || {};
  const nmr = /NMR/i.test(structure.method || "");
  const evidence = item.stability_evidence || design?.foldx_diagnostics || {};
  const source = provenance.source || {};
  const identity = provenance.structure_identity || state.result.source || {};
  const sourceRows = [
    [recommendationLabel("结构与链", "Structure and chain"), `${structure.entry_id || identity.identifier || state.result.structure_selection?.selected_identifier || recommendationLabel("来源未记录", "Source not recorded")} · ${recommendationLabel("链", "chain")} ${identity.chain_id || state.result.structure_selection?.selected_chain || state.result.source?.chain_selection?.chain_id || "—"}`],
    [recommendationLabel("实验方法与质量指标", "Experimental method and quality metrics"), `${structure.method || recommendationLabel("方法未记录", "Method not recorded")} · ${nmr ? recommendationLabel("NMR不采用晶体学分辨率", "NMR does not use crystallographic resolution") : typeof structure.resolution_a === "number" ? recommendationNumber(structure.resolution_a, "Å") : recommendationLabel("分辨率未记录／不适用", "Resolution not recorded / not applicable")}`],
    [recommendationLabel("残基映射", "Residue mapping"), recommendationLabel("突变采用UniProt位置；结构定位使用作者链号、残基号和插入码。", "Mutations use UniProt positions; structure uses author chain, residue number and insertion code.")],
    [recommendationLabel("非加性诊断", "Nonadditivity diagnostic"), recommendationNumber(evidence.epistasis_estimate, "kcal/mol")],
  ];
  const technicalRows = [[recommendationLabel("候选追踪编号", "Candidate trace identifier"), item.candidate_id || "—"]];
  if (state.result.source?.structure_sha256) technicalRows.push([recommendationLabel("输入结构指纹", "Input structure fingerprint"), state.result.source.structure_sha256]);
    const models = { ...(source.models || {}), ...(provenance.models || {}), ...(evidence.provenance?.models || {}), ...(evidence.provenance || {}) };
  const toolNames = { spurs_single: "SPURS · " + recommendationLabel("单点", "single"), spurs_multi: "SPURS · " + recommendationLabel("组合", "combination"), foldx: "FoldX", discotope: "DiscoTope", procaliper: "ProCaliper" };
  for (const [tool, record] of Object.entries(models)) {
    if (!record || typeof record !== "object" || !toolNames[tool]) continue;
    const version = record.tool_inventory?.version || record.provenance?.version || record.version;
    const metadata = record.provenance || record;
    const revision = metadata.model_revision || record.tool_inventory?.model_revision || record.model_revision;
    const series = typeof revision === "string" ? revision.match(/^DiscoTope-(\d+(?:\.\d+)*)$/i)?.[1] : null;
    let versionText = version ? `${recommendationLabel("软件／安装包版本", "Software / installed package version")} ${version}` : recommendationLabel("软件版本未记录；计算来源已保存", "Software version not recorded; computation provenance retained");
    if (tool === "discotope") versionText = `${series ? `${recommendationLabel("模型系列", "Model series")} ${series}` : recommendationLabel("模型系列未记录", "Model series not recorded")} · ${version ? `${recommendationLabel("安装包版本", "Installed package version")} ${version}` : recommendationLabel("安装包版本未记录", "Installed package version not recorded")}`;
    if (metadata.version_identification || record.version_identification) versionText += ` · ${recommendationLabel("当前环境补充识别，非历史运行时版本记录", "Supplement identified from current environment, not a historical runtime version record")}`;
    sourceRows.push([toolNames[tool], versionText]);
    const records = [record, ...(Object.values(record).filter((value) => value && typeof value === "object" && !Array.isArray(value)))];
    for (const model of records) {
      for (const [key, value] of Object.entries(model)) {
        if (/sha256|fingerprint|revision|weights|path|binary|executable/i.test(key) && typeof value === "string") technicalRows.push([`${toolNames[tool]} · ${recommendationLabel("追溯标识", "trace identifier")}`, value]);
      }
    }
  }
  const siteEvidence = sites.map((site) => {
    const confidence = site.coordinate_confidence || structure;
    const confidenceText = confidence.coordinate_available === false ? recommendationLabel("坐标不可用，需要审查", "Coordinates unavailable; review required")
      : confidence.mapping_matches === false ? recommendationLabel("残基映射未匹配，需要审查", "Residue mapping mismatch; review required")
        : confidence.experimental === true && confidence.plddt_applicable === false ? recommendationLabel("实验结构：pLDDT不适用；局部精度尚未评估", "Experimental structure: pLDDT is not applicable; local precision not evaluated")
          : typeof site.plddt === "number" ? `pLDDT ${recommendationNumber(site.plddt)}`
            : confidence.plddt_applicable === true ? recommendationLabel("预测结构：缺少pLDDT证据", "Predicted structure: missing pLDDT evidence")
              : recommendationLabel("来源类型待核对；没有pLDDT记录，不能推断局部可靠性", "Source type unconfirmed; no pLDDT record, local reliability cannot be inferred");
    const exposure = site.sidechain_exposure || {};
    const signedDescriptors = site.chemical_descriptors || {};
    const siteRanges = item.impact_evidence?.display_chemistry_ranges?.ranges?.per_site_signed_deltas?.[site.notation] || {};
    const signedText = chemicalAxisDefinitions.map(([key, zh, en, unit]) => `${recommendationLabel(zh, en)}: ${chemistryRangeValue(signedDescriptors[key] ?? siteRanges[key], unit)}`).join(" · ");
    return [site.notation || String(site.position ?? "—"),
    `${site.chain_id || "—"}: ${site.auth_seq_id ?? "—"}${site.insertion_code || ""}`,
    recommendationNumber(site.sasa_a2, "Å²"), recommendationNumber(site.sidechain_sasa_a2, "Å²"),
    confidenceText,
    site.epitope_assessment?.conclusion || (site.epitope_supported == null ? recommendationLabel("支持标记未知", "Support flag unknown") : site.epitope_supported ? recommendationLabel("达到表位先验支持标记条件", "Epitope prior support flag satisfied") : typeof site.epitope_prior === "number" ? recommendationLabel("有先验评分，未达支持标记条件", "Prior score present; support flag condition not met") : recommendationLabel("未标记支持，先验评分缺失", "Not flagged as supported; prior score missing")),
    recommendationNumber(site.epitope_prior),
    ...(showSiteImpact ? [`${recommendationNumber(exposure.value, "", true)} · ${recommendationNumber(exposure.numerator_a2, "Å²")} / ${recommendationNumber(exposure.denominator_a2, "Å²")} · ${exposure.status || recommendationLabel("状态未知", "status unknown")} · ${exposure.method || recommendationLabel("方法未知", "method unknown")}`, signedText] : [])];
  });
  const comparison = item.primary_comparison;
  const comparisons = comparison ? [
    [recommendationLabel("三尺度平均本区覆盖", "Mean region coverage"), "mean_coverage_fraction", "coverage_difference", true],
    [recommendationLabel("最不利稳定性阈值占用", "Worst stability threshold utilization"), "worst_threshold_utilization", "worst_threshold_utilization_difference", true],
    [recommendationLabel("有表位先验支持的位点比例", "Epitope-supported site fraction"), "epitope_supported_fraction", "epitope_supported_fraction_difference", true],
    [recommendationLabel("平均表位先验评分", "Mean epitope prior score"), "mean_epitope_prior", "mean_epitope_prior_difference", false],
  ].map(([label, key, deltaKey, percent]) => {
    const difference = comparison[deltaKey];
    const delta = typeof difference === "number" && Number.isFinite(difference)
      ? `${difference > 0 ? "+" : ""}${fmt(difference * (percent ? 100 : 1), 2)}${percent ? recommendationLabel(" 个百分点", " percentage points") : ""}`
      : recommendationLabel("无法比较：基准或证据缺失", "Cannot compare: missing baseline or evidence");
    return [label, recommendationNumber(comparison.baseline_metrics?.[key], "", percent), recommendationNumber(comparison.candidate_metrics?.[key], "", percent), delta];
  }) : [];
  return `<details class="recommendation-details"><summary>${recommendationLabel("查看详细证据", "View detailed evidence")}</summary>
    ${showImpact ? `<h4>${recommendationLabel("化学变化与暴露度", "Chemical changes and exposure")}</h4>${regionalImpactEvidence(item, sites)}${showSiteImpact ? `<p class="note">${recommendationLabel("组合轴按位点绝对变化汇总；位点表列出带符号变化。", "Combination axes summarize per-site absolute changes; the site table lists signed changes.")}</p>` : ""}` : ""}
    ${showImpact ? regionalInteractionEvidence(item) : ""}
    <h4>${recommendationLabel("突变位点与单点预测", "Mutation sites and single-site predictions")}</h4>
    ${recommendationTable([recommendationLabel("突变", "Mutation"), recommendationLabel("归属", "Membership"), recommendationLabel("相对溶剂暴露度（RSASA）", "Relative solvent exposure (RSASA)"), "SPURS ΔΔG", recommendationLabel("保护／风险结论", "Protection / risk conclusion")], rows)}
    ${recommendationTable([recommendationLabel("突变", "Mutation"), recommendationLabel("作者链与残基编号", "Author chain / residue identifier"), recommendationLabel("总表面积", "Total surface area"), recommendationLabel("侧链表面积", "Side-chain surface area"), recommendationLabel("坐标可靠性证据", "Coordinate reliability evidence"), recommendationLabel("表位先验标记", "Epitope prior flag"), recommendationLabel("先验评分", "Prior score"), ...(showSiteImpact ? [recommendationLabel("相对侧链暴露度与计算来源", "Relative side-chain exposure and calculation provenance"), recommendationLabel("逐位带符号化学变化", "Per-site signed chemical changes")] : [])], siteEvidence)}
    <p class="note">${recommendationLabel("表位先验是计算表面排序与局部峰值证据，不是抗体结合概率；先验评分存在不等于达到支持标记条件。", "Epitope priors describe computed surface ranking and local peaks, not antibody binding probability; a prior score can exist without satisfying the support flag.")}</p>
    <h4>${recommendationLabel("覆盖范围与完整未覆盖列表", "Coverage and complete uncovered list")}</h4>${recommendationCoverage(item, true)}
    <p>${recommendationLabel("直接突变位点表面积", "Direct mutated-site surface area")}: ${recommendationNumber(item.coverage_gain_area_a2, "Å²")} (${recommendationNumber(item.coverage_gain, "", true)})</p>
    <h4>${recommendationLabel("空间互补性", "Spatial complementarity")}</h4>
    ${recommendationTable([recommendationLabel("邻域", "Neighborhood"), recommendationLabel("新增独立面积", "Added independent area"), recommendationLabel("平均重叠", "Mean overlap"), recommendationLabel("最大重叠", "Max overlap"), recommendationLabel("移除各位点的面积损失", "Area lost when each site is removed")], spatialRows)}
    <p>${recommendationLabel("位点间表面距离", "Pairwise surface distance")}: ${(item.pair_distances || item.geometry_support?.pair_distances || []).map((pair) => `${(pair.positions || [pair.left, pair.right]).join("–")}: ${recommendationNumber(pair.geodesic_distance_a ?? pair.distance_a, "Å")}`).join("；") || recommendationLabel("未计算", "Not computed")}</p>
    <h4>${recommendationLabel("组合准入检查", "Combination qualification checks")}</h4>
    ${gates.length ? recommendationTable([recommendationLabel("审查项", "Check"), recommendationLabel("结论", "Result"), recommendationLabel("观察值／规则结果", "Observed / rule outcome"), recommendationLabel("门槛", "Threshold"), recommendationLabel("超限量", "Excess"), recommendationLabel("说明", "Explanation")], gates) : `<p>${recommendationLabel("此历史结果没有保存完整通过项，请查看下载证据", "This historical result has no complete passed checks; see downloaded evidence")}</p>`}
    <h4>${recommendationLabel("推荐依据", "Recommendation rationale")}</h4><p>${h(cleanRecommendationCopy((item.recommendation_reasons || []).join("；")) || recommendationLabel("本区覆盖优先；计算资格独立审核", "Region coverage first; qualification checked independently"))}</p>
    ${comparison ? `<h4>${recommendationLabel("与本区首选的定量比较", "Quantitative comparison with region primary")}</h4>${recommendationTable([recommendationLabel("指标", "Metric"), recommendationLabel("本区首选", "Region primary"), recommendationLabel("此备选", "This alternative"), recommendationLabel("备选减首选", "Alternative minus primary")], comparisons)}<p class="note">${h(comparison.explanation || recommendationLabel("覆盖和先验越高越符合对应原则；稳定性阈值占用越低表示预测风险越低。各自资格仍独立审查。", "Higher coverage and prior match their respective principles; lower threshold utilization indicates lower predicted risk. Qualification remains independently checked."))}</p>` : ""}
    <h4>${recommendationLabel("证据来源与构建序列", "Evidence provenance and construct sequences")}</h4>${recommendationTable([recommendationLabel("项目", "Item"), recommendationLabel("详情", "Details")], sourceRows)}
    <p class="note">${h(structure.experimental === true ? (nmr ? recommendationLabel("NMR实验结构的局部精度需结合结构集合离散度、实验约束与残基覆盖评估。本次没有计算这些局部指标；pLDDT与晶体学分辨率不适用。", "Local precision of NMR structures uses ensemble dispersion, experimental restraints and residue coverage. These local metrics were not evaluated; pLDDT and crystallographic resolution do not apply.") : recommendationLabel("实验来源已记录；分辨率是整体指标。本次没有评价局部电子密度、占有率或B因子。", "Experimental source is recorded; resolution is a global metric. Local density, occupancy and B factors were not evaluated.")) : recommendationLabel("预测结构的pLDDT描述局部模型置信度。", "pLDDT describes local model confidence."))}</p>
    <details class="recommendation-technical"><summary>${recommendationLabel("追溯编号、文件指纹与运行路径", "Trace identifiers, fingerprints and run paths")}</summary>${recommendationTable([recommendationLabel("记录", "Record"), recommendationLabel("原始标识", "Original identifier")], technicalRows)}${item.candidate_id ? `<button class="btn btn-mini" type="button" data-copy-candidate="${h(item.candidate_id)}">${recommendationLabel("复制候选编号", "Copy candidate identifier")}</button>` : ""}</details>
    <details class="recommendation-sequences"><summary>${recommendationLabel("构建序列", "Construct sequences")} · ${sequence ? sequence.length : "—"} aa · ${coding ? coding.length : "—"} nt</summary><p>${recommendationLabel("突变ECD氨基酸序列", "Mutant ECD amino acid sequence")}</p><pre class="recommendation-sequence">${h(sequence || recommendationLabel("未生成", "Not generated"))}</pre>
    <p>${recommendationLabel("已有编码序列", "Existing coding sequence")}</p><pre class="recommendation-sequence">${h(coding || recommendationLabel("未生成；不推测缺失编码序列", "Not generated; missing coding sequence is not inferred"))}</pre></details>
    
  </details>`;
}

function renderRecommendationCard(item, isPrimary = false, formalLocalCandidate = false) {
  const mechanism = isRegionalMechanism(item);
  const design = (state.result.designs || []).find((value) => value.id === item.design_id);
  const principles = item.principles || [isPrimary ? "sector_coverage" : "cross_region_coverage"];
  const gateFailures = (item.gate_details || []).map((gate) => ({ ...(item.gate_checks || []).find((check) => check.gate === gate.gate), ...gate }));
  const hasCandidateMutations = Boolean(item.candidate_id)
    && ((item.mutation_notations || []).length > 0 || (design?.mutations || []).length > 0);
  const title = formalLocalCandidate ? recommendationLabel("正式合格本区组合", "Formally qualified local combination") : isPrimary ? (mechanism ? hasCandidateMutations
    ? recommendationLabel("本区强扰动工程代表", "Regional strong-perturbation engineering representative")
    : recommendationLabel("本区扰动建议缺口", "Regional perturbation recommendation gap")
    : item.ownership || !item.candidate_id ? recommendationLabel("本区专属首选", "Exclusive region primary") : recommendationLabel("本区历史建议", "Historical region suggestion")) : principles.map(recommendationCode).join(" · ");
  const mutations = item.mutation_notations || design?.mutations?.map(notation) || [];
  const hasMechanismRepresentative = mechanism && isPrimary && Boolean(item.candidate_id) && mutations.length > 0;
  const mechanismGap = mechanism && isPrimary && !hasMechanismRepresentative;
  const exclusive = item.ownership === "EXCLUSIVE";
  const stability = item.stability_evidence || design?.foldx_diagnostics || {};
  const repeatNote = recommendationRepeatNote(item, design);
  const riskLabel = recommendationRiskLabelText(item, design);
  const repeatNeedsReview = recommendationNeedsReview(item, design);
  const statusLabel = repeatNeedsReview ? riskLabel : recommendationStatus(item.status);
  const modelGateNames = ["combination_spurs", "combination_foldx"];
  const modelGateChecks = modelGateNames.map((name) => (item.gate_checks || []).find((gate) => gate.gate === name));
  const modelGatesPassed = modelGateChecks.every((gate) => gate?.passed === true)
    || (!(item.gate_checks || []).length && stability.consensus_status === "CONSENSUS_PASS");
  const hasFiniteModelScores = [stability.spurs_multi_ddg_kcal_mol, stability.foldx_ddg_kcal_mol]
    .every((value) => typeof value === "number" && Number.isFinite(value));
  const acceptedModelStatus = item.status === "PASSED"
    || (item.status === "EVIDENCE_LIMITED" && item.formal_eligible === true);
  const bothModels = hasFiniteModelScores && acceptedModelStatus && modelGatesPassed;
  const panelLabel = item.selected_in_panel ? (bothModels
    ? item.status === "EVIDENCE_LIMITED"
      ? recommendationLabel("已纳入正式主面板（稳定性模型通过；扰动证据不全）", "In formal primary panel (stability models passed; perturbation evidence incomplete)")
      : recommendationLabel("已纳入正式主面板（双模型通过）", "In formal primary panel (both models passed)")
    : recommendationLabel("已纳入正式主面板；模型记录需核对", "In formal primary panel; model records need checking"))
    : recommendationLabel("分区独立推荐／审查 · 资格见状态", "Independent regional recommendation/review · see qualification status");
  const positions = design?.mutations?.map((mutation) => mutation.position) || mutations.map((mutation) => Number(mutation.match(/\d+/)?.[0]));
  const owners = new Set((state.result.regions || []).filter((region) => positions.some((position) => (region.residues || []).includes(position))).map((region) => region.id));
  const allMapped = positions.length > 0 && positions.every((position) => (state.result.regions || []).some((region) => (region.residues || []).includes(position)));
  const ownership = item.ownership ? (exclusive ? recommendationLabel("专属组合", "Exclusive combination") : recommendationLabel("跨区组合", "Cross-region combination")) : !allMapped ? recommendationLabel("历史组合 · 归属证据不完整", "Historical combination · incomplete membership evidence") : owners.size === 1 ? recommendationLabel("位点仅在一区 · 历史建议", "Sites in one region · historical suggestion") : recommendationLabel("跨区组合 · 历史建议", "Cross-region combination · historical suggestion");
  const rationale = formalLocalCandidate
    ? recommendationLabel("该组合通过本次正式计算资格门槛，并作为完整本区候选保留。", "This combination passed the formal computational gates for this run and is retained as a complete local candidate.")
    : hasMechanismRepresentative
    ? recommendationLabel("此项按本区机制扰动规则选出，用作工程测试代表。", "Selected under the regional mechanism perturbation policy as an engineering test representative.")
    : mechanismGap ? Number(item.formal_candidate_count) > 0
      ? recommendationLabel("机制排序未选出主代表；本区仍有正式合格组合，请展开完整候选列表。", "The mechanism ranking selected no primary representative; formally qualified local combinations remain in the full candidate list.")
      : recommendationLabel("本区没有可用候选组合，未选出工程测试代表。", "No usable candidate combination is available for this region; no engineering test representative was selected.")
      : item.recommendation_reasons?.[0] || (!item.ownership ? recommendationLabel("沿用历史已评分建议；本次未重新实施专属资源保障。", "Uses a historical scored suggestion; exclusive resource guarantees were not rerun.") : isPrimary ? recommendationLabel("在安全专属候选中，先检查计算资格，再优先本区邻域覆盖。", "Among safe exclusive candidates, qualification precedes region neighborhood coverage.") : recommendationLabel("按此设计原则比较合格候选；资格与覆盖独立报告。", "Compare qualified candidates for this design principle; qualification and coverage are reported independently."));
  const warningText = gateFailures.length ? gateFailures.map(gateExplanation).join("；")
    : (item.failed_gates || []).map(recommendationCode).join("；");
  const selection = recommendationSelectionRecord(item, design);
  const selected = selection && state.recommendationSelections?.has?.(selection.id);
  return `<article class="design-card recommendation-card" data-recommendation-role="${isPrimary ? "primary" : "alternative"}" data-recommendation-sector="${h(item.sector_id)}">
    <h4>${h(title)}</h4><div>${mutations.map((mutation) => `<span class="mut-chip" data-design-residue="${Number(mutation.match(/\d+/)?.[0])}">${h(mutation)}</span>`).join(" ")}</div>
    <div class="recommendation-tags"><span class="tag ${repeatNeedsReview || item.status !== "PASSED" ? "warn" : ""}">${h(statusLabel)}</span>${repeatNeedsReview ? `<span class="tag">${h(recommendationStatus(item.status))}</span>` : !riskLabel ? "" : `<span class="tag">${h(riskLabel)}</span>`}<span class="tag">${h(ownership)}</span><span class="tag">${h(panelLabel)}</span><span class="tag">${mutations.length} ${recommendationLabel("个突变", "mutations")}</span></div>
    <p>${h(cleanRecommendationCopy(rationale))}</p>${mechanismGap ? `<p class="recommendation-note">${h(Number(item.formal_candidate_count) > 0 ? recommendationLabel("机制主代表未选出；本区仍有正式合格候选，见完整列表。", "No mechanism primary was selected; formally qualified local candidates remain in the full list.") : recommendationLabel("推荐缺口：当前未得到本区候选组合。", "Recommendation gap: no regional candidate combination is available."))}</p>` : ""}${item.status === "EVIDENCE_LIMITED" && bothModels ? `<p class="recommendation-note">${h(cleanRecommendationCopy(regionalMissingImpactNote(item)))}</p>` : ""}${warningText ? `<p class="recommendation-warning">${h(warningText)}</p>` : ""}
    ${mutations.length ? `<p class="muted">SPURS ΔΔG: ${recommendationNumber(stability.spurs_multi_ddg_kcal_mol, "kcal/mol")} · FoldX ΔΔG: ${recommendationNumber(stability.foldx_ddg_kcal_mol, "kcal/mol")}${repeatNote ? ` · ${h(repeatNote)}` : ""}</p>` : ""}
    ${hasMechanismRepresentative ? regionalMechanismSummary(item) : mechanismGap ? "" : mutations.length ? recommendationCoverage(item) : ""}
    ${(item.principles || []).filter((value) => value !== "sector_coverage").length && isPrimary ? `<p class="muted">${recommendationLabel("同时符合", "Also selected for")}: ${h(item.principles.filter((value) => value !== "sector_coverage").map(recommendationCode).join(" · "))}</p>` : ""}
    ${item.design_id ? `<div class="actions-row"><button class="btn btn-mini" type="button" data-review-design="${h(item.design_id)}">${recommendationLabel("定位结构与序列", "Locate structure and sequence")}</button>${selection ? `<button class="btn btn-mini" type="button" data-select-recommendation="${h(selection.id)}" aria-pressed="${Boolean(selected)}">${selected ? recommendationLabel("取消选择", "Remove selection") : recommendationLabel("加入空间比较", "Add to spatial comparison")}</button>` : ""}<button class="btn btn-mini" type="button" data-copy-design="${h(item.design_id)}">${recommendationLabel("复制到自定义", "Copy to custom")}</button></div>` : ""}
    ${!item.design_id && selection ? `<div class="actions-row"><button class="btn btn-mini" type="button" data-select-recommendation="${h(selection.id)}" aria-pressed="${Boolean(selected)}">${selected ? recommendationLabel("取消选择", "Remove selection") : recommendationLabel("加入空间比较", "Add to spatial comparison")}</button></div>` : ""}
    ${item.candidate_id ? recommendationDetails(item, design) : ""}</article>`;
}

function renderSectorRecommendations(recommendations) {
  return `<section class="design-series region-recommendations"><h3>${recommendationLabel("按分区推荐组合", "Recommendations by region")}</h3>
    <p class="note">${recommendationLabel("各区独立保留不同设计原则下通过门槛的推荐，不受全局精简面板数量限制。同一组合满足多个原则时合并标签；缺少合格方案时列明缺口。表面邻域覆盖按表面测地距离估计。", "Each region retains qualified recommendations for its design principles without the global compact-panel limit. Matching combinations share principle labels; missing qualified options are listed as gaps. Surface neighborhood coverage uses surface geodesic distance.")}</p>
    ${recommendations.map((item) => {
      const sectorIndex = (state.result.joint_v4_ml.sectors || []).findIndex((sector) => sector.sector_id === item.sector_id);
      const region = (state.result.regions || [])[sectorIndex] || {};
      const summary = item.sector_summary || {};
      const localCandidates = Array.isArray(item.all_formal_local_candidates) ? item.all_formal_local_candidates : null;
      const qualifiedLocalCandidates = localCandidates?.filter(recommendationQualified) || [];
      const displayedReviewIds = new Set([item, ...(item.alternatives || [])]
        .filter((candidate) => candidate.candidate_id)
        .map((candidate) => candidate.candidate_id));
      const reviewLocalCandidates = (localCandidates || []).filter((candidate) =>
        !recommendationQualified(candidate) && !displayedReviewIds.has(candidate.candidate_id));
      const alternatives = (item.alternatives || []).filter(recommendationQualified);
      const reviewAlternatives = (item.alternatives || []).filter((candidate) => !recommendationQualified(candidate));
      const primaryReview = Boolean(item.candidate_id) && !recommendationQualified(item);
      const area = summary.surface_area_a2 ?? item.regional_geometry_support?.sector_area_a2 ?? region.total_sasa ?? region.surface_area;
      const audit = item.resource_audit;
      return `<section class="recommendation-region"><h3><span class="swatch" style="background:${h(region.color || "#668899")}"></span> ${h(recommendationRegion(item.sector_id))}</h3>
        <p class="muted">${summary.residue_count ?? region.residues?.length ?? "—"} ${recommendationLabel("个残基", "residues")} · ${recommendationNumber(area, "Å²")} · ${summary.mutable_site_count ?? region.mutable_residues?.length ?? "—"} ${recommendationLabel("个可变异位点", "mutable sites")}</p>
        ${item.formal_candidate_count != null || item.confirmed_alternative_count != null ? `<p class="note">${recommendationLabel("机制主代表状态", "Mechanism primary status")}: ${h(recommendationDisplayStatus(item))}${recommendationRepeatNote(item) ? ` · ${h(recommendationRepeatNote(item))}` : ""} · ${recommendationLabel("本区正式合格候选", "Formally qualified local candidates")}: ${item.formal_candidate_count ?? recommendationLabel("未记录", "not recorded")} · ${recommendationLabel("已确认替代", "Confirmed alternatives")}: ${item.confirmed_alternative_count ?? recommendationLabel("未记录", "not recorded")}</p>` : ""}
        ${primaryReview ? `<p class="note">${recommendationLabel("本区主推荐尚无通过门槛的方案；其他原则的合格推荐见下方。", "No qualified primary is available for this region; qualified options for other principles appear below.")}</p>` : renderRecommendationCard(item, true)}<h4>${recommendationLabel("其他设计原则下的合格推荐", "Qualified recommendations for other principles")}</h4>
        ${alternatives.map((alternative) => renderRecommendationCard({ ...alternative, sector_id: item.sector_id })).join("") || `<p class="muted">${recommendationLabel("没有不同的合格备选；同一合格组合的原则合并显示。", "No distinct qualified alternative; matching principles share one qualified combination.")}</p>`}
        ${primaryReview || reviewAlternatives.length || reviewLocalCandidates.length ? `<details class="regional-review-candidates"><summary>${recommendationLabel("未通过正式展示条件／待复核的审查详情", "Review details: not qualified for formal display or review required")}</summary>${primaryReview ? renderRecommendationCard(item, true) : ""}${reviewAlternatives.map((candidate) => renderRecommendationCard({ ...candidate, sector_id: item.sector_id })).join("")}${reviewLocalCandidates.map((candidate) => renderRecommendationCard({ ...candidate, sector_id: item.sector_id }, false, true)).join("")}</details>` : ""}
        ${Array.isArray(localCandidates) ? `<details class="formal-local-candidates"><summary>${recommendationLabel("全部正式展示合格本区候选（含主代表）", "All qualified local candidates for formal display (including the primary)")} · ${qualifiedLocalCandidates.length}</summary>${qualifiedLocalCandidates.map((candidate) => renderRecommendationCard({ ...candidate, sector_id: item.sector_id }, false, true)).join("") || `<p class="muted">${recommendationLabel("当前没有可正式展示的本区组合。", "No local combinations qualify for formal display.")}</p>`}</details>` : ""}
        ${(item.unavailable_principles || []).map((value) => `<p class="note">${h(recommendationCode(value.principle))}: ${h(cleanRecommendationCopy(/[\u4e00-\u9fff]/.test(value.reason || "") ? value.reason : recommendationCode(value.reason)))}</p>${value.gate_details?.length ? `<details class="principle-review-details"><summary>${recommendationLabel("未通过门槛／待补证据", "Failed gates / missing evidence")}</summary><code>${h(value.candidate_id)}</code><ul>${value.gate_details.map((gate) => `<li>${h(gateExplanation(gate))}</li>`).join("")}</ul></details>` : ""}`).join("")}
        ${audit ? `<details><summary>${recommendationLabel("搜索与验证情况", "Search and validation resources")}</summary>${recommendationResourceAudit(audit)}</details>` : `<p class="note">${recommendationLabel("此历史任务未记录新的分区资源审计。", "This historical job has no new region resource audit.")}</p>`}
      </section>`;
    }).join("")}${selectionPanel()}</section>`;
}

function recommendationResourceAudit(audit) {
  const labels = { pool_candidate_count: ["进入候选池", "In candidate pool"], pool_exclusive_count: ["入池专属候选", "Exclusive pool candidates"], pool_exclusive_groups: ["入池不同专属位点组", "Distinct exclusive site sets in pool"], spurs_scored_count: ["完成SPURS评分", "SPURS scored"], foldx_scheduled_count: ["安排FoldX", "Scheduled FoldX"], foldx_success_count: ["FoldX评分成功", "FoldX successfully scored"], formal_eligible_count: ["最终合格组合", "Qualified combinations"], formal_exclusive_count: ["最终合格专属组合", "Qualified exclusive combinations"], generated_exclusive_candidates: ["生成专属候选", "Generated exclusive candidates"], confirmable_exclusive_groups: ["可验证专属位点组", "Confirmable exclusive site sets"], scheduled_exclusive_groups: ["已安排专属验证位点组", "Scheduled exclusive site sets"], scheduled_associated_groups: ["已安排关联验证位点组", "Scheduled associated site sets"] };
  const rows = Object.entries(labels).filter(([key]) => audit[key] != null).map(([key, label]) => [recommendationLabel(...label), String(audit[key])]);
  for (const [key, label] of [["generated_exclusive_site_sets", "生成的专属位点组"], ["pool_exclusive_site_sets", "入池专属位点组"], ["quota_deficits", "基础配额缺口"]]) {
    if (audit[key]) rows.push([recommendationLabel(label, key === "quota_deficits" ? "Base quota deficit" : "Exclusive site sets"), [2,3].map((order) => `${order} ${recommendationLabel("点", "sites")}: ${audit[key][String(order)] ?? "—"}`).join(" · ")]);
  }
  return recommendationTable([recommendationLabel("阶段", "Stage"), recommendationLabel("数量", "Count")], rows) + `<p>${h((audit.reasons || []).map((reason) => /[\u4e00-\u9fff]/.test(reason) ? reason : recommendationCode(reason)).join("；"))}</p>`;
}
  return { sectorRecommendationText, renderSectorRecommendations, recommendationLabel, surfaceGateDescription,
    recommendationSelectionRecord, recommendationRiskLabelText, recommendationNeedsReview,
    recommendationRepeatSummary, selectionColor };
}
