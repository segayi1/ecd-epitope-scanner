export function createRecommendationRenderer({ state, h, fmt, fmtPct, notation }) {
function sectorRecommendationText(item) {
  return [recommendationStatus(item.status), ...(item.gate_details || []).map(gateExplanation)].join("；");
}

function recommendationLabel(zh, en) { return state.lang === "zh" ? zh : en; }
function recommendationStatus(status) {
  return ({ PASSED: recommendationLabel("通过计算门槛", "Passed computational gates"),
    NEAR_THRESHOLD: recommendationLabel("待审查 · 未通过", "Review required · not qualified"),
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
    "engineered-neighborhood": ["该位点邻近输入结构中的工程化差异，需要与目标天然序列对照；不是本次突变已造成损伤", "This site is near engineered differences in the input structure; compare with the target native sequence. This is not evidence that the designed mutation caused damage"],
    COMBINATION_REJECT: ["组合未满足当前计算推荐门槛；具体超限或缺失证据见下方", "Combination does not meet the current computational recommendation gates; see specific excess or missing evidence below"],
    COMBINATION_SPURS_MISSING: ["缺少SPURS组合稳定性证据", "Missing SPURS combination evidence"],
    COMBINATION_FOLDX_MISSING: ["缺少FoldX组合稳定性证据", "Missing FoldX combination evidence"],
    COMBINATION_SPURS_ABOVE_THRESHOLD: ["SPURS预测值超过当前门槛", "SPURS prediction exceeds the threshold"],
    COMBINATION_FOLDX_ABOVE_THRESHOLD: ["FoldX预测值超过当前门槛", "FoldX prediction exceeds the threshold"],
    MODEL_DISAGREEMENT: ["两种模型预测差异超过允许范围", "Model predictions disagree beyond the allowed range"],
    HIGHER_ORDER_REQUIRES_FOLDX: ["高阶组合仍需FoldX确认", "Higher-order combination needs FoldX confirmation"],
    INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN: ["高阶组合新增收益未满足最低要求", "Higher-order additional gain does not meet the minimum requirement"],
    HIGHER_ORDER_COVERAGE_BASELINE_MISSING: ["缺少同区低阶几何对照，无法确认新增覆盖收益；缺失不能当作零", "Missing a lower-order geometric comparison in the same region; additional coverage is unknown, not zero"],
    SINGLE_EVIDENCE_MISSING_OR_EXCLUDED: ["单点证据缺失或不符合准入条件", "Single-site evidence missing or excluded"],
    NO_CANDIDATE_WITH_MUTATION_IN_SECTOR: ["没有位点位于本区的安全候选", "No safe candidate mutates this region"],
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
    structural_confidence: ["满足当前坐标准入规则；不意味着存在预测模型置信度分数", "Current coordinate admission rules satisfied; this does not imply a predicted-model confidence score"],
    formal_eligibility: ["完整计算资格检查通过，仍需实验验证", "Full computational qualification passed; experimental validation remains necessary"],
    model_evidence_complete: ["SPURS与FoldX组合评分均已获得有限数值；这不是实验稳定性证明", "Finite SPURS and FoldX combination scores are available; this does not prove experimental stability"],
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
  const detail = item.site_details || [];
  const sites = detail.length ? detail : (design?.mutations || []).map((mutation) => ({ ...mutation, notation: notation(mutation) }));
  const rows = sites.map((site) => {
    const position = site.position ?? site.uniprot_pos;
    const residue = (state.result.residues || []).find((row) => row.uniprot_pos === position) || {};
    const assessment = site.risk_assessment;
    const riskText = assessment ? [assessment.conclusion, ...(assessment.notes || []).map((note) => note.explanation?.includes("尚未解释") ? recommendationCode(note.code) : note.explanation)].filter(Boolean).join("；") : [...new Set([...(site.risk_notes || site.hard_exclusions || []), ...(site.warnings || []), ...(site.structural_risk_flags || [])])].map((value) => /[\u4e00-\u9fff]/.test(value) && !value.includes("\ufffd") ? value : recommendationCode(value)).join("；") || recommendationLabel("未记录额外提示；不等于排除全部风险", "No additional note recorded; this does not rule out all risks");
    return [site.notation || `${site.native_aa || site.wildtype || residue.one_letter || ""}${position ?? "?"}${site.mutant_aa || site.mutant || ""}`,
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
    return [site.notation || String(site.position ?? "—"),
    `${site.chain_id || "—"}: ${site.auth_seq_id ?? "—"}${site.insertion_code || ""}`,
    recommendationNumber(site.sasa_a2, "Å²"), recommendationNumber(site.sidechain_sasa_a2, "Å²"),
    confidenceText,
    site.epitope_assessment?.conclusion || (site.epitope_supported == null ? recommendationLabel("支持标记未知", "Support flag unknown") : site.epitope_supported ? recommendationLabel("达到表位先验支持标记条件", "Epitope prior support flag satisfied") : typeof site.epitope_prior === "number" ? recommendationLabel("有先验评分，未达支持标记条件", "Prior score present; support flag condition not met") : recommendationLabel("未标记支持，先验评分缺失", "Not flagged as supported; prior score missing")),
    recommendationNumber(site.epitope_prior)];
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
    <h4>${recommendationLabel("突变位点与单点预测", "Mutation sites and single-site predictions")}</h4>
    ${recommendationTable([recommendationLabel("突变", "Mutation"), recommendationLabel("归属", "Membership"), recommendationLabel("相对溶剂暴露度（RSASA）", "Relative solvent exposure (RSASA)"), "SPURS ΔΔG", recommendationLabel("保护／风险结论", "Protection / risk conclusion")], rows)}
    ${recommendationTable([recommendationLabel("突变", "Mutation"), recommendationLabel("作者链与残基编号", "Author chain / residue identifier"), recommendationLabel("总表面积", "Total surface area"), recommendationLabel("侧链表面积", "Side-chain surface area"), recommendationLabel("坐标可靠性证据", "Coordinate reliability evidence"), recommendationLabel("表位先验标记", "Epitope prior flag"), recommendationLabel("先验评分", "Prior score")], siteEvidence)}
    <p class="note">${recommendationLabel("表位先验是计算表面排序与局部峰值证据，不是抗体结合概率；先验评分存在不等于达到支持标记条件。", "Epitope priors describe computed surface ranking and local peaks, not antibody binding probability; a prior score can exist without satisfying the support flag.")}</p>
    <h4>${recommendationLabel("覆盖范围与完整未覆盖列表", "Coverage and complete uncovered list")}</h4>${recommendationCoverage(item, true)}
    <p>${recommendationLabel("直接突变位点表面积", "Direct mutated-site surface area")}: ${recommendationNumber(item.coverage_gain_area_a2, "Å²")} (${recommendationNumber(item.coverage_gain, "", true)})</p>
    <h4>${recommendationLabel("空间互补性", "Spatial complementarity")}</h4>
    ${recommendationTable([recommendationLabel("邻域", "Neighborhood"), recommendationLabel("新增独立面积", "Added independent area"), recommendationLabel("平均重叠", "Mean overlap"), recommendationLabel("最大重叠", "Max overlap"), recommendationLabel("移除各位点的面积损失", "Area lost when each site is removed")], spatialRows)}
    <p>${recommendationLabel("位点间表面距离", "Pairwise surface distance")}: ${(item.pair_distances || item.geometry_support?.pair_distances || []).map((pair) => `${(pair.positions || [pair.left, pair.right]).join("–")}: ${recommendationNumber(pair.geodesic_distance_a ?? pair.distance_a, "Å")}`).join("；") || recommendationLabel("未计算", "Not computed")}</p>
    <h4>${recommendationLabel("组合准入检查", "Combination qualification checks")}</h4>
    ${gates.length ? recommendationTable([recommendationLabel("审查项", "Check"), recommendationLabel("结论", "Result"), recommendationLabel("观察值／规则结果", "Observed / rule outcome"), recommendationLabel("门槛", "Threshold"), recommendationLabel("超限量", "Excess"), recommendationLabel("说明", "Explanation")], gates) : `<p>${recommendationLabel("此历史结果没有保存完整通过项，请查看下载证据", "This historical result has no complete passed checks; see downloaded evidence")}</p>`}
    <h4>${recommendationLabel("推荐依据", "Recommendation rationale")}</h4><p>${h((item.recommendation_reasons || []).join("；") || recommendationLabel("本区覆盖优先；计算资格独立审核", "Region coverage first; qualification checked independently"))}</p>
    ${comparison ? `<h4>${recommendationLabel("与本区首选的定量比较", "Quantitative comparison with region primary")}</h4>${recommendationTable([recommendationLabel("指标", "Metric"), recommendationLabel("本区首选", "Region primary"), recommendationLabel("此备选", "This alternative"), recommendationLabel("备选减首选", "Alternative minus primary")], comparisons)}<p class="note">${h(comparison.explanation || recommendationLabel("覆盖和先验越高越符合对应原则；稳定性阈值占用越低表示预测风险越低。各自资格仍独立审查。", "Higher coverage and prior match their respective principles; lower threshold utilization indicates lower predicted risk. Qualification remains independently checked."))}</p>` : ""}
    <h4>${recommendationLabel("证据来源与构建序列", "Evidence provenance and construct sequences")}</h4>${recommendationTable([recommendationLabel("项目", "Item"), recommendationLabel("详情", "Details")], sourceRows)}
    <p class="note">${h(structure.experimental === true ? (nmr ? recommendationLabel("NMR实验结构的局部精度需结合结构集合离散度、实验约束与残基覆盖评估。本次没有计算这些局部指标；pLDDT与晶体学分辨率不适用，不能将它们的缺失解释为低置信度。", "Local precision of NMR structures requires ensemble dispersion, experimental restraints and residue coverage. These local metrics were not evaluated; pLDDT and crystallographic resolution do not apply, so their absence is not low confidence.") : recommendationLabel("实验来源不能直接证明每个侧链都可靠。分辨率是整体指标；本次没有评价局部电子密度、占有率或B因子，不能将缺少pLDDT解释为低置信度。", "Experimental origin does not prove every side chain is reliable. Resolution is a global metric; local density, occupancy and B factors were not evaluated. Absence of pLDDT does not mean low confidence.")) : recommendationLabel("预测结构的pLDDT只描述局部模型置信度；来源不明或分数缺失时不能当作高置信度。", "Predicted-structure pLDDT describes local model confidence; unknown sources or missing scores do not imply high confidence."))}</p>
    <details class="recommendation-technical"><summary>${recommendationLabel("追溯编号、文件指纹与运行路径", "Trace identifiers, fingerprints and run paths")}</summary>${recommendationTable([recommendationLabel("记录", "Record"), recommendationLabel("原始标识", "Original identifier")], technicalRows)}${item.candidate_id ? `<button class="btn btn-mini" type="button" data-copy-candidate="${h(item.candidate_id)}">${recommendationLabel("复制候选编号", "Copy candidate identifier")}</button>` : ""}</details>
    <details class="recommendation-sequences"><summary>${recommendationLabel("构建序列", "Construct sequences")} · ${sequence ? sequence.length : "—"} aa · ${coding ? coding.length : "—"} nt</summary><p>${recommendationLabel("突变ECD氨基酸序列", "Mutant ECD amino acid sequence")}</p><pre class="recommendation-sequence">${h(sequence || recommendationLabel("未生成", "Not generated"))}</pre>
    <p>${recommendationLabel("已有编码序列", "Existing coding sequence")}</p><pre class="recommendation-sequence">${h(coding || recommendationLabel("未生成；不推测缺失编码序列", "Not generated; missing coding sequence is not inferred"))}</pre></details>
    <p class="note">${recommendationLabel("非加性与稳定性预测不是实验协同作用或结构破坏的证明。", "Nonadditivity and stability predictions do not demonstrate experimental synergy or structural damage.")}</p>
  </details>`;
}

function renderRecommendationCard(item, isPrimary = false) {
  const design = (state.result.designs || []).find((value) => value.id === item.design_id);
  const principles = item.principles || [isPrimary ? "sector_coverage" : "cross_region_coverage"];
  const gateFailures = (item.gate_details || []).map((gate) => ({ ...(item.gate_checks || []).find((check) => check.gate === gate.gate), ...gate }));
  const title = isPrimary ? (item.ownership || !item.candidate_id ? recommendationLabel("本区专属首选", "Exclusive region primary") : recommendationLabel("本区历史建议", "Historical region suggestion")) : principles.map(recommendationCode).join(" · ");
  const mutations = item.mutation_notations || design?.mutations?.map(notation) || [];
  const exclusive = item.ownership === "EXCLUSIVE";
  const stability = item.stability_evidence || design?.foldx_diagnostics || {};
  const bothModels = [stability.spurs_multi_ddg_kcal_mol, stability.foldx_ddg_kcal_mol].every((value) => typeof value === "number" && Number.isFinite(value))
    && item.status === "PASSED" && ((item.gate_checks || []).length ? ["combination_spurs", "combination_foldx"].every((name) => item.gate_checks.some((gate) => gate.gate === name && gate.passed === true)) : stability.consensus_status === "CONSENSUS_PASS");
  const panelLabel = item.selected_in_panel ? (bothModels ? recommendationLabel("已纳入正式主面板（双模型通过）", "In formal primary panel (both models passed)") : recommendationLabel("已纳入正式主面板；模型记录需核对", "In formal primary panel; model records need checking")) : recommendationLabel("分区建议 · 未纳入正式主面板", "Region suggestion · outside formal primary panel");
  const positions = design?.mutations?.map((mutation) => mutation.position) || mutations.map((mutation) => Number(mutation.match(/\d+/)?.[0]));
  const owners = new Set((state.result.regions || []).filter((region) => positions.some((position) => (region.residues || []).includes(position))).map((region) => region.id));
  const allMapped = positions.length > 0 && positions.every((position) => (state.result.regions || []).some((region) => (region.residues || []).includes(position)));
  const ownership = item.ownership ? (exclusive ? recommendationLabel("专属组合", "Exclusive combination") : recommendationLabel("跨区组合", "Cross-region combination")) : !allMapped ? recommendationLabel("历史组合 · 归属证据不完整", "Historical combination · incomplete membership evidence") : owners.size === 1 ? recommendationLabel("位点仅在一区 · 历史建议", "Sites in one region · historical suggestion") : recommendationLabel("跨区组合 · 历史建议", "Cross-region combination · historical suggestion");
  const rationale = item.recommendation_reasons?.[0] || (!item.ownership ? recommendationLabel("沿用历史已评分建议；本次未重新实施专属资源保障。", "Uses a historical scored suggestion; exclusive resource guarantees were not rerun.") : isPrimary ? recommendationLabel("在安全专属候选中，先检查计算资格，再优先本区邻域覆盖。", "Among safe exclusive candidates, qualification precedes region neighborhood coverage.") : recommendationLabel("按此设计原则比较合格候选；资格与覆盖独立报告。", "Compare qualified candidates for this design principle; qualification and coverage are reported independently."));
  return `<article class="design-card recommendation-card" data-recommendation-role="${isPrimary ? "primary" : "alternative"}" data-recommendation-sector="${h(item.sector_id)}">
    <h4>${h(title)}</h4><div>${mutations.map((mutation) => `<span class="mut-chip" data-design-residue="${Number(mutation.match(/\d+/)?.[0])}">${h(mutation)}</span>`).join(" ")}</div>
    <div class="recommendation-tags"><span class="tag ${item.status === "PASSED" ? "" : "warn"}">${h(recommendationStatus(item.status))}</span><span class="tag">${h(ownership)}</span><span class="tag">${h(panelLabel)}</span><span class="tag">${mutations.length} ${recommendationLabel("个突变", "mutations")}</span></div>
    <p>${h(rationale)}</p>${item.status !== "PASSED" ? `<p class="recommendation-warning">${h(gateFailures.length ? gateFailures.map(gateExplanation).join("；") : (item.failed_gates || []).map(recommendationCode).join("；") || recommendationLabel("没有安全的专属候选，不能建立计算资格。", "No safe exclusive candidate establishes qualification."))}</p>` : ""}
    ${mutations.length ? `<p class="muted">SPURS ΔΔG: ${recommendationNumber(stability.spurs_multi_ddg_kcal_mol, "kcal/mol")} · FoldX ΔΔG: ${recommendationNumber(stability.foldx_ddg_kcal_mol, "kcal/mol")}</p>` : ""}
    ${mutations.length ? recommendationCoverage(item) : ""}
    ${(item.principles || []).filter((value) => value !== "sector_coverage").length && isPrimary ? `<p class="muted">${recommendationLabel("同时符合", "Also selected for")}: ${h(item.principles.filter((value) => value !== "sector_coverage").map(recommendationCode).join(" · "))}</p>` : ""}
    ${item.design_id ? `<div class="actions-row"><button class="btn btn-mini" type="button" data-review-design="${h(item.design_id)}">${recommendationLabel("定位结构与序列", "Locate structure and sequence")}</button><button class="btn btn-mini" type="button" data-copy-design="${h(item.design_id)}">${recommendationLabel("复制到自定义", "Copy to custom")}</button></div>` : ""}
    ${item.candidate_id ? recommendationDetails(item, design) : ""}</article>`;
}

function renderSectorRecommendations(recommendations) {
  return `<section class="design-series region-recommendations"><h3>${recommendationLabel("按分区推荐组合", "Recommendations by region")}</h3>
    <p class="note">${recommendationLabel("表面邻域覆盖是几何估计，不能直接代表真实扰动范围或结合变化。待审查建议不计入合格主面板；主面板数量仍由独立预算控制。", "Surface neighborhood coverage is a geometric estimate, not measured perturbation or binding change. Review suggestions are outside the qualified primary panel; its separate budget still applies.")}</p>
    ${recommendations.map((item) => {
      const sectorIndex = (state.result.joint_v4_ml.sectors || []).findIndex((sector) => sector.sector_id === item.sector_id);
      const region = (state.result.regions || [])[sectorIndex] || {};
      const summary = item.sector_summary || {};
      const area = summary.surface_area_a2 ?? item.regional_geometry_support?.sector_area_a2 ?? region.total_sasa ?? region.surface_area;
      const audit = item.resource_audit;
      return `<section class="recommendation-region"><h3><span class="swatch" style="background:${h(region.color || "#668899")}"></span> ${h(recommendationRegion(item.sector_id))}</h3>
        <p class="muted">${summary.residue_count ?? region.residues?.length ?? "—"} ${recommendationLabel("个残基", "residues")} · ${recommendationNumber(area, "Å²")} · ${summary.mutable_site_count ?? region.mutable_residues?.length ?? "—"} ${recommendationLabel("个可变异位点", "mutable sites")}</p>
        ${renderRecommendationCard(item, true)}<h4>${recommendationLabel("其他设计原则", "Alternative design principles")}</h4>
        ${(item.alternatives || []).map((alternative) => renderRecommendationCard({ ...alternative, sector_id: item.sector_id })).join("") || `<p class="muted">${recommendationLabel("没有不同的安全备选；相同组合的原则已合并到首选。", "No different safe alternative; matching principles are merged into the primary.")}</p>`}
        ${(item.unavailable_principles || []).map((value) => `<p class="note">${h(recommendationCode(value.principle))}: ${h(/[\u4e00-\u9fff]/.test(value.reason || "") ? value.reason : recommendationCode(value.reason))}</p>`).join("")}
        ${audit ? `<details><summary>${recommendationLabel("搜索与验证情况", "Search and validation resources")}</summary>${recommendationResourceAudit(audit)}</details>` : `<p class="note">${recommendationLabel("此历史任务未记录新的分区资源审计。", "This historical job has no new region resource audit.")}</p>`}
      </section>`;
    }).join("")}</section>`;
}

function recommendationResourceAudit(audit) {
  const labels = { pool_candidate_count: ["进入候选池", "In candidate pool"], pool_exclusive_count: ["入池专属候选", "Exclusive pool candidates"], pool_exclusive_groups: ["入池不同专属位点组", "Distinct exclusive site sets in pool"], spurs_scored_count: ["完成SPURS评分", "SPURS scored"], foldx_scheduled_count: ["安排FoldX", "Scheduled FoldX"], foldx_success_count: ["FoldX评分成功", "FoldX successfully scored"], formal_eligible_count: ["最终合格组合", "Qualified combinations"], formal_exclusive_count: ["最终合格专属组合", "Qualified exclusive combinations"], generated_exclusive_candidates: ["生成专属候选", "Generated exclusive candidates"], confirmable_exclusive_groups: ["可验证专属位点组", "Confirmable exclusive site sets"], scheduled_exclusive_groups: ["已安排专属验证位点组", "Scheduled exclusive site sets"], scheduled_associated_groups: ["已安排关联验证位点组", "Scheduled associated site sets"] };
  const rows = Object.entries(labels).filter(([key]) => audit[key] != null).map(([key, label]) => [recommendationLabel(...label), String(audit[key])]);
  for (const [key, label] of [["generated_exclusive_site_sets", "生成的专属位点组"], ["pool_exclusive_site_sets", "入池专属位点组"], ["quota_deficits", "基础配额缺口"]]) {
    if (audit[key]) rows.push([recommendationLabel(label, key === "quota_deficits" ? "Base quota deficit" : "Exclusive site sets"), [2,3].map((order) => `${order} ${recommendationLabel("点", "sites")}: ${audit[key][String(order)] ?? "—"}`).join(" · ")]);
  }
  return recommendationTable([recommendationLabel("阶段", "Stage"), recommendationLabel("数量", "Count")], rows) + `<p>${h((audit.reasons || []).map((reason) => /[\u4e00-\u9fff]/.test(reason) ? reason : recommendationCode(reason)).join("；"))}</p>`;
}
  return { sectorRecommendationText, renderSectorRecommendations, recommendationLabel };
}
