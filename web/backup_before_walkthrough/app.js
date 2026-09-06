const $ = (id) => document.getElementById(id);
const AA = "ACDEFGHIKLMNPQRSTVWY".split("");
const CODONS = {
  A: "GCC", R: "CGC", N: "AAC", D: "GAC", C: "TGC", Q: "CAG", E: "GAG",
  G: "GGC", H: "CAC", I: "ATC", L: "CTG", K: "AAG", M: "ATG", F: "TTC",
  P: "CCC", S: "AGC", T: "ACC", W: "TGG", Y: "TAC", V: "GTG"
};

const state = {
  result: null,
  jobId: null,
  eventSource: null,
  viewer: null,
  viewerReady: false,
  selectedResidue: null,
  selectedRegion: null,
  selectedDesign: "wt",
  sequenceMode: "aa",
  representation: "cartoon",
  layers: new Set(["folding", "regions", "design"]),
  layerFocus: "folding",
  previewPlans: new Map(),
  spinning: false,
  customMutations: [],
  demo: false,
  engine: "molstar",       // molstar (cartoon) | beads (3Dmol)
  beadsReady: false,
  lang: "zh",
  i18n: {},
};

const STEP_META = {
  ecd: { label: "已映射 ECD", color: "#9AA8B2" },
  surface: { label: "充分表面可及", color: "#1565C0" },
  folding: { label: "折叠风险筛选通过", color: "#2E7D32" },
  interface: { label: "天然功能/结构界面", color: "#E4572E" },
  regions: { label: "结构语义分区", color: "#6A4C93" },
  design: { label: "当前组合突变", color: "#E4572E" },
};

function h(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

/* ---- i18n (static chrome; scientific wording from the backend stays as-is) */
function t(key, fallback) {
  const entry = state.i18n?.[state.lang]?.[key]
    ?? state.i18n?.zh?.[key]
    ?? state.i18n?.en?.[key];
  if (entry == null) {
    if (fallback != null) return fallback;
    console.warn("i18n missing key:", key);
    return key;
  }
  return entry;
}

async function loadI18n() {
  for (const lang of ["zh", "en"]) {
    try {
      state.i18n[lang] = await api(`/assets/i18n/${lang}.json`);
    } catch (error) {
      console.error(`i18n load failed for ${lang}`, error);
      state.i18n[lang] = {};
    }
  }
  const saved = localStorage.getItem("ecd.lang");
  state.lang = saved === "en" || saved === "zh" ? saved : "zh";
}

function applyStaticI18n() {
  document.querySelectorAll("[data-i18n]").forEach((node) => {
    node.textContent = t(node.dataset.i18n, node.textContent);
  });
  $("btn-lang").textContent = state.lang === "zh" ? "EN" : "中文";
}

function toggleLang() {
  state.lang = state.lang === "zh" ? "en" : "zh";
  localStorage.setItem("ecd.lang", state.lang);
  applyStaticI18n();
  if (state.result) { renderStageRail(); renderLegend(); }
}

/* ---- data-driven interface label (P0 of the "no CD80 hard-coding" rule) */
function interfacePartnerChains() {
  if (!state.result) return [];
  const chains = new Set();
  for (const residue of state.result.residues) {
    if (!residue.interface) continue;
    for (const partner of residue.interface_partners || []) {
      chains.add(String(partner).split(":")[0]);
    }
  }
  return [...chains].sort();
}

function interfaceLabel() {
  const chains = interfacePartnerChains();
  if (!chains.length) return t("interface.none", "无蛋白伙伴（apo/单体）");
  return t("interface.with", "界面") + " " + chains.join("/");
}

function notation(mutation) {
  return `${mutation.wt}${mutation.position}${mutation.mutant}`;
}

function showToast(message) {
  const toast = $("toast");
  toast.textContent = message;
  toast.hidden = false;
  clearTimeout(showToast.timer);
  showToast.timer = setTimeout(() => { toast.hidden = true; }, 4000);
}

async function api(url, options = {}) {
  const response = await fetch(url, options);
  const contentType = response.headers.get("content-type") || "";
  const data = contentType.includes("json") ? await response.json() : await response.text();
  if (!response.ok) {
    const detail = typeof data === "object" ? data.detail ?? data.error ?? data : data;
    const message = typeof detail === "object" ? detail.message ?? JSON.stringify(detail) : detail;
    const error = new Error(message || `HTTP ${response.status}`);
    error.detail = detail;
    throw error;
  }
  return data;
}

function setJobState(status, message, progress = 0, error = null) {
  const badge = $("job-state");
  const names = {
    QUEUED: "排队中", RUNNING: "运行中", COMPLETE: "完成",
    WAITING_FOR_AF3_UPLOAD: "等待 AF3", FAILED: "失败", DEMO: "冻结示例"
  };
  const classes = {
    QUEUED: "running", RUNNING: "running", COMPLETE: "complete",
    WAITING_FOR_AF3_UPLOAD: "waiting", FAILED: "failed", DEMO: "neutral"
  };
  badge.textContent = names[status] || status || "未运行";
  badge.className = `status-badge ${classes[status] || "neutral"}`;
  $("job-message").textContent = message || "";
  const percent = Math.max(0, Math.min(100, Math.round((progress || 0) * 100)));
  $("job-percent").textContent = `${percent}%`;
  $("job-progress").style.width = `${percent}%`;
  $("job-error").hidden = !error;
  $("job-error").textContent = error || "";
  $("run-button").disabled = status === "RUNNING" || status === "QUEUED";
  const exportButton = $("export-offline");
  if (exportButton) exportButton.disabled = !(status === "COMPLETE");
}

function queryForResidue(residue) {
  if (residue.label_seq_id != null && residue.label_chain_id) {
    return { label_asym_id: residue.label_chain_id, label_seq_id: residue.label_seq_id };
  }
  const query = { auth_asym_id: residue.chain_id, auth_seq_id: residue.auth_seq_id };
  if (residue.insertion_code) query.pdbx_PDB_ins_code = residue.insertion_code;
  return query;
}

async function ensureViewer(structureUrl) {
  const target = $("viewer");
  try {
    state.viewerReady = false;
    if (!state.viewer) {
      state.viewer = new window.PDBeMolstarPlugin();
      await state.viewer.render(target, {
        customData: { url: structureUrl, format: "cif", binary: false },
        defaultPreset: "default",
        visualStyle: { polymer: state.representation, het: "ball-and-stick", carbs: "carbohydrate" },
        hideStructure: ["water"],
        bgColor: "#f4f0e7",
        highlightColor: "#FFD166",
        selectColor: "#FFD166",
        granularity: "residue",
        hideControls: true,
        hideCanvasControls: ["all"],
        sequencePanel: false,
        leftPanel: false,
        rightPanel: false,
        logPanel: false,
        pdbeLink: false,
        loadingOverlay: true,
        reactive: true,
        subscribeEvents: false,
      });
    } else {
      await state.viewer.visual.update({
        customData: { url: structureUrl, format: "cif", binary: false },
        alphafoldView: false,
      }, true);
    }
    state.viewerReady = true;
    $("viewer-placeholder").classList.add("hidden");
    await installViewerTooltips();
    await applyViewerLayer();
    // Keep the primary work surface on the mapped ECD rather than the full-length
    // chain, which may include flexible tails that make the antigen appear tiny.
    await state.viewer.visual.focus(state.result.residues.map(queryForResidue));
  } catch (error) {
    console.error(error);
    $("viewer-placeholder").classList.remove("hidden");
    $("viewer-placeholder").querySelector("strong").textContent = "结构画布载入失败";
    $("viewer-placeholder").querySelector("span").textContent = error.message;
  }
}

async function installViewerTooltips() {
  if (!state.viewerReady || !state.result) return;
  const data = state.result.residues.map((residue) => ({
    ...queryForResidue(residue),
    tooltip: `${residue.aa}${residue.uniprot_pos} | rSASA ${Number(residue.rsasa).toFixed(3)} | side-chain SASA ${Number(residue.sidechain_sasa).toFixed(1)} A2 | ${residue.interface ? "functional interface" : "non-interface"}`,
  }));
  await state.viewer.visual.tooltips({ data });
}

function selectedDesign() {
  if (!state.result || state.selectedDesign === "wt") return null;
  return state.previewPlans.get(state.selectedDesign)
    || state.result.designs.find((design) => design.id === state.selectedDesign)
    || null;
}

async function applyViewerLayer() {
  if (!state.viewerReady || !state.result) return;
  const data = [];
  const append = (records, color, extra = {}) => records.forEach((residue) => {
    data.push({ ...queryForResidue(residue), color, ...extra });
  });
  if (state.layers.has("surface")) append(state.result.residues.filter((item) => item.surface_accessible), STEP_META.surface.color);
  if (state.layers.has("folding")) append(state.result.residues.filter((item) => item.folding_screen_pass), STEP_META.folding.color);
  if (state.layers.has("interface")) append(state.result.residues.filter((item) => item.interface), STEP_META.interface.color, { sideChain: true, representation: "ball-and-stick", representationColor: STEP_META.interface.color });
  if (state.layers.has("regions")) {
    for (const region of state.result.regions) {
      const positions = new Set(region.residues);
      append(state.result.residues.filter((item) => positions.has(item.uniprot_pos)), region.color);
    }
  }
  if (state.layers.has("design")) {
    const positions = new Set((selectedDesign()?.mutations || []).map((item) => item.position));
    append(state.result.residues.filter((item) => positions.has(item.uniprot_pos)), STEP_META.design.color, { sideChain: true, representation: "ball-and-stick", representationColor: STEP_META.design.color });
  }
  await state.viewer.visual.select({ data, nonSelectedColor: "#7c909d", keepRepresentations: false });
  if (state.selectedResidue) {
    await state.viewer.visual.highlight({
      data: [queryForResidue(state.selectedResidue)], color: "#FFD166", focus: false
    });
  }
  renderLegend();
  /* keep the beads engine in sync with layer/selection changes */
  if (state.engine === "beads") paintBeads();
}

function renderLegend() {
  if (!state.result) {
    $("legend").innerHTML = '<span class="legend-item">等待结构结果</span>';
    return;
  }
  const items = [];
  for (const key of ["surface", "folding", "interface", "design"]) {
    if (!state.layers.has(key)) continue;
    const label = key === "interface" ? interfaceLabel() : STEP_META[key].label;
    items.push(`<span class="legend-item"><i class="legend-swatch" style="background:${STEP_META[key].color}"></i>${h(label)}</span>`);
  }
  if (state.layers.has("regions")) {
    items.push(...state.result.regions.map((region) => `<span class="legend-item"><i class="legend-swatch" style="background:${h(region.color)}"></i>分区 ${region.id}</span>`));
  }
  items.push('<span class="legend-item"><i class="legend-swatch" style="background:#7c909d"></i>ECD 骨架</span>');
  $("legend").innerHTML = items.join("");
}

function criteriaText() {
  if (!state.result) return "载入结果后显示本步骤的定量标准。";
  const m = state.result.methods;
  const texts = {
    ecd: `结构中可映射到 UniProt ${state.result.protein.ecd_start}-${state.result.protein.ecd_end} 的全部 ECD 残基；缺失坐标不参与设计。`,
    surface: `正向筛选：rSASA ≥ ${m.rsasa_threshold} 且侧链 SASA ≥ ${m.sidechain_sasa_threshold_A2} Å²；SASA 引擎：${m.sasa_method}。`,
    folding: "在表面标准基础上，排除二硫键 Cys、已注释糖基化/PTM、无法映射残基和低置信度 AF3 位点。这里表示低预测风险，不是实验确认。",
    interface: `任一目标重原子距非目标链或异质残基 ≤5 Å。功能界面正常参与设计，但永久保留功能风险标记。`,
    regions: `${state.result.partitioning.method}；目标 5–8 区，实际 ${state.result.regions.length} 区，合格表面覆盖 ${formatPercent(state.result.partitioning.coverage)}。`,
    design: "每区先从最多 4 个关键位点的候选替换中，用 FoldX 单点结果排除不稳定替换；最终推荐组合再进行 FoldX 组合复核。预览子组合必须复核后才可称为低风险。",
  };
  return texts[state.layerFocus] || "当前图层可与其他图层叠加显示。";
}

function setInspectorPanel(panel) {
  document.querySelectorAll(".secondary-nav button").forEach((button) => {
    const active = button.dataset.panel === panel;
    button.classList.toggle("active", active);
  });
  document.querySelectorAll(".inspector-panel").forEach((content) => {
    content.classList.toggle("active", content.dataset.panelContent === panel);
  });
}

function setPrimaryPage(panel) {
  document.querySelectorAll(".primary-pages button").forEach((button) => {
    const active = button.dataset.primary === panel;
    button.classList.toggle("active", active);
    button.setAttribute("aria-selected", String(active));
  });
  setInspectorPanel(panel);
}

async function setRepresentation(representation) {
  state.representation = representation;
  document.querySelectorAll("[data-representation]").forEach((button) => {
    button.setAttribute("aria-pressed", String(button.dataset.representation === representation));
  });
  if (!state.viewerReady) return;
  await state.viewer.visual.update({ visualStyle: { polymer: representation }, hideStructure: ["water"] }, false);
  await state.viewer.canvas.applySettings({ color: "#f4f0e7", lighting: representation === "molecular-surface" ? "matte" : "flat" });
  await applyViewerLayer();
}

function formatPercent(value) {
  return `${Math.round(Number(value || 0) * 100)}%`;
}

function formatNumber(value, digits = 2) {
  return value == null ? "—" : Number(value).toFixed(digits);
}

function renderSummary() {
  const s = state.result.summary;
  const values = [
    [s.mapped_ecd_residues, "ECD 已映射"],
    [s.surface_accessible, "充分表面可及"],
    [s.folding_screen_pass, "折叠风险通过"],
    [s.functional_interface, "功能/结构界面"],
    [s.regions, "结构分区"],
    [s.recommended_constructs, "推荐组合"],
  ];
  $("summary-stats").innerHTML = values.map(([value, label]) =>
    `<div class="stat"><strong>${h(value)}</strong><span>${h(label)}</span></div>`
  ).join("");
}

function renderStageRail() {
  if (!state.result) return;
  const stages = [
    ["ecd", "S1", "结构 / ECD", "overview"],
    ["surface", "S2", "表面可及", "layers"],
    ["folding", "S3", "FoldX 单点筛选", "layers"],
    ["interface", "S4", "界面标识", "layers"],
    ["regions", "S5", "结构分区", "regions"],
    ["design", "S6", "组合设计", "design"],
  ];
  $("stage-rail").innerHTML = stages.map(([key, step, label, panel], index) => `
    <button type="button" class="stage-card ${state.layerFocus === key ? "active" : ""}" data-stage="${key}" data-stage-panel="${panel}" aria-label="${step} ${label}">
      <span class="stage-number">${step}</span><span class="stage-label">${label}</span><span class="stage-state">${index < 4 ? "已计算" : key === "regions" ? `${state.result.regions.length} 区` : `${state.result.designs.length} 条`}</span>
    </button>`).join("");
  document.querySelectorAll("[data-stage]").forEach((button) => button.addEventListener("click", () => {
    state.layerFocus = button.dataset.stage;
    if (button.dataset.stage === "ecd") state.layers.add("folding");
    else state.layers.add(button.dataset.stage);
    if (button.dataset.stagePanel === "overview") setPrimaryPage("overview");
    else if (button.dataset.stagePanel === "regions") setPrimaryPage("regions");
    else if (button.dataset.stagePanel === "design") setPrimaryPage("design");
    else setInspectorPanel("layers");
    renderLayerControls();
    applyViewerLayer();
    renderStageRail();
  }));
}

function renderSource() {
  const source = state.result.source;
  const selection = state.result.structure_selection;
  const rows = [
    ["蛋白", `${state.result.protein.protein_name} (${state.result.protein.accession})`],
    ["ECD", `${state.result.protein.ecd_start}-${state.result.protein.ecd_end} · ${state.result.protein.ecd_length} aa`],
    ["拓扑", state.result.protein.topology],
    ["结构", `${source.identifier || "user upload"} · ${source.provider || source.kind}`],
    ["链与覆盖", `${selection.selected_chain} · ${formatPercent(selection.observed_ecd_coverage)}`],
    ["内部缺口", `${selection.max_internal_gap} aa`],
    ["序列一致性", formatPercent(selection.aligned_sequence_identity)],
    ["模型类型", source.explicit_af3 ? "明确 AF3 来源" : source.category || source.kind],
  ];
  const list = $("source-detail");
  list.classList.remove("empty-state");
  list.innerHTML = rows.map(([key, value]) => `<dt>${h(key)}</dt><dd>${h(value)}</dd>`).join("") +
    (source.page_url ? `<dt>来源页</dt><dd><a href="${h(source.page_url)}" target="_blank" rel="noreferrer">打开权威来源</a></dd>` : "");
}

function renderMethods() {
  const stability = state.result.stability;
  $("method-detail").innerHTML = [
    `<strong>SASA</strong>　${h(state.result.methods.sasa_method)}`,
    `<strong>分区</strong>　${h(state.result.partitioning.method)}`,
    `<strong>稳定性</strong>　${h(stability.primary)}；${h(stability.proxy_units)}`,
    `<strong>结论边界</strong>　计算筛选不能替代表达、构象抗体、功能和稳定性实测。`,
  ].join("<br>");
  const warnings = state.result.warnings?.length ? state.result.warnings : ["没有运行级警告。"];
  $("warning-list").innerHTML = warnings.map((warning) => `<li>${h(warning)}</li>`).join("");
}

function renderLayerControls() {
  if (!state.result) return;
  const entries = [
    ["surface", "充分表面可及", "rSASA 与侧链 SASA 双阈值", STEP_META.surface.color],
    ["folding", "折叠风险通过", "保护位点排除后的候选", STEP_META.folding.color],
    ["interface", "功能/结构界面", "保留功能风险标记", STEP_META.interface.color],
    ["regions", "结构语义分区", "按区域颜色显示", STEP_META.regions.color],
    ["design", "当前突变组合", "橙色球棍叠加", STEP_META.design.color],
  ];
  $("layer-controls").innerHTML = entries.map(([key, label, detail, color]) => `
    <div class="layer-row">
      <label><input type="checkbox" data-layer="${key}" ${state.layers.has(key) ? "checked" : ""}>
      <i style="background:${color}"></i>
      <span><strong>${label}</strong><small>${detail}</small></span></label>
      <button type="button" data-layer-focus="${key}">定位</button>
    </div>`).join("");
  document.querySelectorAll("[data-layer]").forEach((input) => input.addEventListener("change", () => {
    if (input.checked) state.layers.add(input.dataset.layer);
    else state.layers.delete(input.dataset.layer);
    state.layerFocus = input.dataset.layer;
    $("criterion").textContent = criteriaText();
    applyViewerLayer();
  }));
  document.querySelectorAll("[data-layer-focus]").forEach((button) => button.addEventListener("click", () => {
    const key = button.dataset.layerFocus;
    state.layers.add(key);
    state.layerFocus = key;
    renderLayerControls();
    $("criterion").textContent = criteriaText();
    applyViewerLayer();
  }));
  $("criterion").textContent = criteriaText();
}

function aminoAcidName(aa) {
  return ({ A: "Ala", R: "Arg", N: "Asn", D: "Asp", C: "Cys", Q: "Gln", E: "Glu", G: "Gly", H: "His", I: "Ile", L: "Leu", K: "Lys", M: "Met", F: "Phe", P: "Pro", S: "Ser", T: "Thr", W: "Trp", Y: "Tyr", V: "Val" })[aa] || aa;
}

function secondaryName(code) {
  return ({ H: "α-螺旋", E: "β-链", C: "环 / 无规则" })[code] || "未注释";
}

function designForResidue(residue) {
  const regions = state.result.regions || [];
  const region = regions.find((item) => item.id === residue.region_id);
  const design = state.result.designs.find((item) => item.region_id === residue.region_id);
  if (!region) return { text: "未分配至自动分区", design: null, region: null };
  const kind = residue.interface ? "functional-interface patch" : region.structural_character;
  return { text: `分区 R${region.id}（${kind}）`, design, region };
}

function scanChoicesForResidue(residue) {
  for (const design of state.result.designs || []) {
    const scan = design.saturation_scan?.[String(residue.uniprot_pos)];
    if (scan) return scan;
  }
  return [];
}

function hasFoldXScore(design) {
  return Boolean(design?.stability_engine && /foldx/i.test(design.stability_engine) && design.combined_stability_score != null);
}

function designScoreText(design) {
  if (!design) return "未生成自动组合";
  if (hasFoldXScore(design)) return `FoldX ΔΔG ${formatNumber(design.combined_stability_score, 3)} kcal/mol`;
  return "FoldX 组合复核待执行";
}

function mutationDiffersInCharge(mutation) {
  const positive = new Set(["R", "K", "H"]);
  const negative = new Set(["D", "E"]);
  const charge = (aa) => positive.has(aa) ? 1 : negative.has(aa) ? -1 : 0;
  return charge(mutation.wt) !== charge(mutation.mutant);
}

function makePreviewPlan(design, suffix, label, mutations) {
  const wtAa = state.result.protein.ecd_sequence;
  const wtDna = state.result.protein.ecd_coding_sequence || [...wtAa].map((aa) => CODONS[aa]).join("");
  const start = state.result.protein.ecd_start;
  const aa = [...wtAa];
  const dna = wtDna.match(/.{1,3}/g) || [];
  mutations.forEach((mutation) => {
    const index = mutation.position - start;
    if (index < 0 || index >= aa.length || aa[index] !== mutation.wt) return;
    aa[index] = mutation.mutant;
    if (dna[index]) dna[index] = CODONS[mutation.mutant];
  });
  return {
    ...design,
    id: `${design.id}_${suffix}`,
    name: `${design.name} · ${label}`,
    mutations: mutations.map((mutation) => ({ ...mutation })),
    mutant_ecd_sequence: aa.join(""),
    mutant_coding_sequence: dna.join(""),
    stability_engine: "FoldX final-combination review pending",
    combined_stability_score: null,
    status: "FOLDX_REVIEW_REQUIRED",
    preview: true,
  };
}

function regionPlanOptions(design) {
  if (!design?.mutations?.length) return [];
  const all = design.mutations;
  const minimum = all.slice(0, Math.min(2, all.length));
  const charged = all.filter(mutationDiffersInCharge);
  const chargePlan = charged.length >= 2 ? charged : all.slice(0, Math.min(2, all.length));
  return [
    { id: design.id, label: `推荐组合 · ${all.length} 位点 · ${designScoreText(design)}`, design },
    { id: `${design.id}_minimal`, label: `最小扰动预览 · ${minimum.length} 位点 · 待 FoldX 组合复核`, design: makePreviewPlan(design, "minimal", "最小扰动预览", minimum) },
    { id: `${design.id}_charge`, label: `电荷扰动预览 · ${chargePlan.length} 位点 · 待 FoldX 组合复核`, design: makePreviewPlan(design, "charge", "电荷扰动预览", chargePlan) },
  ];
}

function filteredResidues() {
  if (!state.result) return [];
  const filter = $("residue-filter").value;
  const query = $("residue-search").value.trim().toLowerCase();
  return state.result.residues.filter((residue) => {
    let matches = true;
    if (filter === "folding") matches = residue.folding_screen_pass;
    if (filter === "surface") matches = residue.surface_accessible;
    if (filter === "interface") matches = residue.interface;
    if (filter === "protected") matches = residue.hard_exclusions.length > 0;
    if (!matches || !query) return matches;
    const text = [
      residue.uniprot_pos, `${residue.aa}${residue.uniprot_pos}`, residue.domain,
      residue.interface ? "interface" : "", ...(residue.risk_flags || []), ...(residue.hard_exclusions || [])
    ].join(" ").toLowerCase();
    return text.includes(query);
  });
}

function renderResidueList() {
  if (!state.result) return;
  const residues = filteredResidues();
  $("residue-list").innerHTML = residues.map((residue) => {
    const active = state.selectedResidue?.uniprot_pos === residue.uniprot_pos ? " active" : "";
    const status = residue.folding_screen_pass ? "通过" : residue.surface_accessible ? "风险" : "排除";
    return `<button type="button" class="residue-row${active}" data-residue="${residue.uniprot_pos}" role="option" aria-selected="${Boolean(active.trim())}">
      <span class="residue-id">${h(residue.aa)}${residue.uniprot_pos}</span>
      <span class="residue-main"><strong>${h(residue.domain)} · 分区 ${residue.region_id ?? "—"}</strong><span>${h(status)}${residue.interface ? " · 功能界面" : ""}${residue.risk_flags.length ? ` · ${h(residue.risk_flags.join(", "))}` : ""}</span></span>
      <span class="residue-score">rSASA ${formatNumber(residue.rsasa, 2)}</span>
    </button>`;
  }).join("") || '<div class="empty-state" style="padding:14px">当前筛选没有残基。</div>';
  document.querySelectorAll("[data-residue]").forEach((button) => {
    button.addEventListener("click", () => selectResidue(Number(button.dataset.residue), true));
  });
}

function renderResidueDetail() {
  const residue = state.selectedResidue;
  if (!residue) {
    $("residue-detail").className = "empty-state";
    $("residue-detail").textContent = "在结构、序列或列表中选择残基。";
    return;
  }
  const tags = [];
  const mapping = designForResidue(residue);
  const singleChoices = scanChoicesForResidue(residue)
    .filter((choice) => choice.allowed)
    .sort((left, right) => (left.foldx_ddg ?? left.stability_score) - (right.foldx_ddg ?? right.stability_score))
    .slice(0, 5);
  tags.push(`<span class="tag ${residue.folding_screen_pass ? "pass" : "block"}">${residue.folding_screen_pass ? "折叠风险通过" : "未通过自动筛选"}</span>`);
  if (residue.interface) tags.push('<span class="tag interface">功能/结构界面</span>');
  (residue.risk_flags || []).forEach((tag) => tags.push(`<span class="tag">${h(tag)}</span>`));
  (residue.hard_exclusions || []).forEach((tag) => tags.push(`<span class="tag block">${h(tag)}</span>`));
  const options = AA.filter((aa) => aa !== residue.aa).map((aa) => `<option value="${aa}">${aa}</option>`).join("");
  $("residue-detail").className = "";
  $("residue-detail").innerHTML = `
    <div class="detail-title"><div><h3 class="residue-audit-title">${h(residue.aa)}${residue.uniprot_pos} <span>· ${aminoAcidName(residue.aa)} · 位置 ${residue.uniprot_pos}</span></h3><div class="empty-state">PDB ${h(residue.chain_id)}:${residue.auth_seq_id}${h(residue.insertion_code)} · ${h(residue.residue_name)}</div></div><button id="focus-residue" type="button" class="text-button">结构定位</button></div>
    <div class="tag-row">${tags.join("")}</div>
    <table class="audit-table"><tbody>
      <tr><th>状态</th><td>${h(mapping.text)}${mapping.design ? `；${h(mapping.design.status)}` : ""}</td></tr>
      <tr><th>可及性</th><td>rSASA=${formatNumber(residue.rsasa, 2)}（总 SASA ${formatNumber(residue.sasa, 0)} Å²；侧链 ${formatNumber(residue.sidechain_sasa, 0)} Å²）</td></tr>
      <tr><th>局部堆积</th><td>Cα 接触数 ${residue.contacts}；二级结构 ${secondaryName(residue.secondary_structure)}${residue.salt_bridge ? "；盐桥风险" : ""}</td></tr>
      <tr><th>界面</th><td>${residue.interface ? `与 ${h((residue.interface_partners || []).join(", ") || "非目标链")} 形成接触` : "未检测到结构界面接触"}</td></tr>
      <tr><th>结构保护</th><td>${h([...(residue.hard_exclusions || []), ...(residue.risk_flags || [])].join("；") || "无特殊基序")}</td></tr>
      <tr><th>置信度</th><td>${residue.plddt == null ? "实验结构 / 未提供 pLDDT" : `pLDDT=${formatNumber(residue.plddt, 1)}`}</td></tr>
    </tbody></table>
    <strong>判定理由</strong><ul class="audit-reasons">${(residue.reasons || []).map((reason) => `<li>${h(reason)}</li>`).join("") || "<li>没有额外判定信息。</li>"}</ul>
    <div class="single-mutation-options"><h3>候选单点替换</h3>${singleChoices.length ? singleChoices.map((choice) => `<div class="candidate-mutation-row"><code>${h(choice.notation || `${choice.wt}${choice.position}${choice.mutant}`)}</code><small>${choice.foldx_ddg != null ? `FoldX ΔΔG ${formatNumber(choice.foldx_ddg, 3)} kcal/mol` : choice.stability_class || "预筛"}</small><button type="button" data-add-choice="${choice.mutant}">加入</button></div>`).join("") : '<span class="empty-state">当前结构未保存该位点的单点替换扫描。</span>'}</div>
    <div class="mutation-adder"><select id="custom-aa">${options}</select><button id="add-custom-site" type="button" class="secondary-button">加入自定义构建</button></div>`;
  $("focus-residue").addEventListener("click", () => focusSelectedResidue());
  $("add-custom-site").addEventListener("click", () => addCustomMutation(residue, $("custom-aa").value));
  document.querySelectorAll("[data-add-choice]").forEach((button) => button.addEventListener("click", () => addCustomMutation(residue, button.dataset.addChoice)));
}

async function selectResidue(position, focus = false) {
  if (!state.result) return;
  const residue = state.result.residues.find((item) => item.uniprot_pos === position);
  if (!residue) return;
  state.selectedResidue = residue;
  renderResidueList();
  renderResidueDetail();
  renderSequences();
  if (state.viewerReady) {
    await state.viewer.visual.highlight({ data: [queryForResidue(residue)], color: "#FFD166", focus: false });
    if (focus) await focusSelectedResidue();
  }
}

async function focusSelectedResidue() {
  if (!state.viewerReady || !state.selectedResidue) return;
  const query = queryForResidue(state.selectedResidue);
  await state.viewer.visual.focus([query]);
  await state.viewer.visual.interactivityFocus({ data: [query] });
}

function renderRegions() {
  if (!state.result) return;
  const designByRegion = new Map(state.result.designs.map((design) => [design.region_id, design]));
  $("region-list").innerHTML = state.result.regions.map((region) => {
    const design = designByRegion.get(region.id);
    const mutations = design?.mutations?.map(notation).join(" · ") || "无自动组合";
    const options = regionPlanOptions(design);
    const current = options.some((option) => option.id === state.selectedDesign) ? state.selectedDesign : "wt";
    const planPicker = options.length ? `<label class="region-plan-picker"><span>此区方案</span><select data-region-plan="${region.id}"><option value="wt">不选择</option>${options.map((option) => `<option value="${h(option.id)}" ${option.id === current ? "selected" : ""}>${h(option.label)}</option>`).join("")}</select></label>` : "";
    return `<article class="region-item" data-region-item="${region.id}">
      <button class="region-header" type="button" data-region="${region.id}" aria-expanded="false">
        <i class="region-color" style="background:${h(region.color)}"></i>
        <span><strong>分区 ${region.id} · ${h(region.name)}</strong><span>${h(region.structural_character)} · ${h(mutations)}</span></span>
        <span class="region-count">${region.residues.length} 位点</span>
      </button>
      <div class="region-body">
        <div class="region-metrics">
          <div>侧链表面积<strong>${formatNumber(region.surface_area, 0)} Å²</strong></div>
          <div>测地直径<strong>${formatNumber(region.geodesic_diameter, 1)} Å</strong></div>
          <div>界面比例<strong>${formatPercent(region.interface_fraction)}</strong></div>
        </div>
        <p>${h(region.rationale)}</p>
        ${planPicker}
        <div class="design-actions"><button type="button" data-show-region="${region.id}">聚焦此区</button>${design?.mutations?.length ? `<button type="button" data-copy-design="${h(design.id)}">复制推荐组合</button>` : ""}</div>
      </div>
    </article>`;
  }).join("") || '<div class="empty-state" style="padding:14px">没有可辩护的结构分区。</div>';
  document.querySelectorAll("[data-region]").forEach((button) => {
    button.addEventListener("click", () => {
      const id = Number(button.dataset.region);
      const item = document.querySelector(`[data-region-item="${id}"]`);
      const willOpen = !item.classList.contains("open");
      document.querySelectorAll("[data-region-item]").forEach((node) => node.classList.remove("open"));
      if (willOpen) item.classList.add("open");
      button.setAttribute("aria-expanded", String(willOpen));
      state.selectedRegion = willOpen ? state.result.regions.find((region) => region.id === id) : null;
      state.layers.add("regions");
      state.layerFocus = "regions";
      renderLayerControls();
      applyViewerLayer();
    });
  });
  document.querySelectorAll("[data-show-region]").forEach((button) => {
    button.addEventListener("click", async () => {
      state.selectedRegion = state.result.regions.find((region) => region.id === Number(button.dataset.showRegion));
      state.layers.add("regions");
      state.layerFocus = "regions";
      renderLayerControls();
      await applyViewerLayer();
      if (state.viewerReady && state.selectedRegion) {
        const positions = new Set(state.selectedRegion.residues);
        await state.viewer.visual.focus(state.result.residues.filter((residue) => positions.has(residue.uniprot_pos)).map(queryForResidue));
      }
    });
  });
  document.querySelectorAll("[data-copy-design]").forEach((button) => {
    button.addEventListener("click", () => copyDesignToCustom(button.dataset.copyDesign));
  });
  document.querySelectorAll("[data-region-plan]").forEach((select) => {
    select.addEventListener("change", () => {
      const regionId = Number(select.dataset.regionPlan);
      const option = regionPlanOptions(designByRegion.get(regionId)).find((item) => item.id === select.value);
      if (option?.design?.preview) state.previewPlans.set(option.id, option.design);
      selectDesign(select.value);
    });
  });
}

function renderDesigns() {
  if (!state.result) return;
  $("design-list").innerHTML = state.result.designs.map((design) => {
    const recommended = design.status.includes("RECOMMENDED") || design.status.includes("CONFIRMED");
    const chips = (design.mutations || []).map((mutation) =>
      `<button type="button" class="mutation-chip" data-design-residue="${mutation.position}" title="定位 ${h(notation(mutation))}">${h(notation(mutation))}</button>`
    ).join("") || '<span class="empty-state">未生成自动组合</span>';
    return `<article class="design-item">
      <div class="design-heading"><div><strong>${h(design.id)} · ${h(design.name)}</strong><span>${h(design.stability_engine)} · ${h(designScoreText(design))}</span></div><span class="design-status ${recommended ? "" : "review"}">${h(design.status)}</span></div>
      <div class="mutation-chips">${chips}</div>
      <div class="design-actions"><button type="button" data-view-design="${h(design.id)}">结构与序列审查</button>${design.mutations?.length ? `<button type="button" data-copy-design="${h(design.id)}">复制到自定义</button>` : ""}</div>
    </article>`;
  }).join("");
  document.querySelectorAll("[data-design-residue]").forEach((button) => {
    button.addEventListener("click", () => selectResidue(Number(button.dataset.designResidue), true));
  });
  document.querySelectorAll("[data-view-design]").forEach((button) => {
    button.addEventListener("click", () => selectDesign(button.dataset.viewDesign));
  });
  document.querySelectorAll("#panel-design [data-copy-design]").forEach((button) => {
    button.addEventListener("click", () => copyDesignToCustom(button.dataset.copyDesign));
  });
}

function renderDesignSelector() {
  const select = $("sequence-design");
  const previews = [...state.previewPlans.values()];
  select.innerHTML = '<option value="wt">WT ECD</option>' + [...state.result.designs, ...previews].map((design) =>
    `<option value="${h(design.id)}">${h(design.id)} · ${h(design.mutations.map(notation).join("; ") || design.status)}</option>`
  ).join("");
  select.value = [...select.options].some((option) => option.value === state.selectedDesign) ? state.selectedDesign : "wt";
}

function selectDesign(id) {
  state.selectedDesign = id;
  renderDesignSelector();
  $("sequence-design").value = id;
  if (id !== "wt") {
    state.layers.add("design");
    state.layerFocus = "design";
  }
  renderLayerControls();
  applyViewerLayer();
  renderSequences();
  if (window.innerWidth <= 860) setMobileView("sequence");
}

function renderSequenceTrack(target, sequence, mutations, mode, isMutant) {
  const mutationPositions = new Set((mutations || []).map((item) => item.position));
  const selectedPosition = state.selectedResidue?.uniprot_pos;
  const start = state.result.protein.ecd_start;
  const protectedPositions = new Set(
    state.result.residues.filter((item) => item.hard_exclusions.length).map((item) => item.uniprot_pos)
  );
  if (mode === "aa") {
    target.innerHTML = [...sequence].map((aa, index) => {
      const position = start + index;
      const classes = ["residue-token"];
      if (isMutant && mutationPositions.has(position)) classes.push("mutated");
      if (protectedPositions.has(position)) classes.push("protected");
      if (selectedPosition === position) classes.push("selected");
      return `<span class="${classes.join(" ")}" data-seq-residue="${position}" title="${aa}${position}">${aa}</span>`;
    }).join("");
  } else {
    const codons = sequence.match(/.{1,3}/g) || [];
    target.innerHTML = codons.map((codon, index) => {
      const position = start + index;
      const classes = ["residue-token"];
      if (isMutant && mutationPositions.has(position)) classes.push("mutated");
      if (protectedPositions.has(position)) classes.push("protected");
      if (selectedPosition === position) classes.push("selected");
      return `<span class="${classes.join(" ")}" style="min-width:31px" data-seq-residue="${position}" title="${position}: ${codon}">${codon}</span>`;
    }).join(" ");
  }
}

function renderSequences() {
  if (!state.result) return;
  const design = selectedDesign();
  const wtAa = state.result.protein.ecd_sequence;
  const wtDna = state.result.protein.ecd_coding_sequence || [...wtAa].map((aa) => CODONS[aa]).join("");
  const mutantAa = design?.mutant_ecd_sequence || wtAa;
  const mutantDna = design?.mutant_coding_sequence || wtDna;
  const mutations = design?.mutations || [];
  state.lastSequences = {
    wtAa, wtDna, mutAa: mutantAa, mutDna: mutantDna,
    hasCds: Boolean(design?.mutant_coding_sequence || state.result.protein.ecd_coding_sequence),
  };
  renderSequenceTrack($("wt-sequence"), state.sequenceMode === "aa" ? wtAa : wtDna, [], state.sequenceMode, false);
  renderSequenceTrack($("mutant-sequence"), state.sequenceMode === "aa" ? mutantAa : mutantDna, mutations, state.sequenceMode, true);
  $("wt-sequence-meta").textContent = `${state.result.protein.ecd_start}-${state.result.protein.ecd_end} · ${state.sequenceMode === "aa" ? `${wtAa.length} aa` : `${wtDna.length} nt`}`;
  $("mutant-sequence-meta").textContent = design ? `${design.id} · ${mutations.map(notation).join("; ") || design.status}` : "WT 对照";
  document.querySelectorAll("[data-seq-residue]").forEach((token) => {
    token.addEventListener("mouseenter", () => {
      const residue = state.result.residues.find((item) => item.uniprot_pos === Number(token.dataset.seqResidue));
      if (!residue || state.selectedResidue?.uniprot_pos === residue.uniprot_pos) return;
      setInspectorPanel("residues");
      selectResidue(residue.uniprot_pos, false);
    });
    token.addEventListener("click", () => selectResidue(Number(token.dataset.seqResidue), true));
  });
}

function addCustomMutation(residue, mutant) {
  const existing = state.customMutations.find((item) => item.position === residue.uniprot_pos);
  if (existing) existing.mutant = mutant;
  else state.customMutations.push({ position: residue.uniprot_pos, wt: residue.aa, mutant });
  state.customMutations.sort((a, b) => a.position - b.position);
  renderCustomBuilder();
  setInspectorPanel("design");
}

function copyDesignToCustom(id) {
  const design = state.previewPlans.get(id) || state.result.designs.find((item) => item.id === id);
  if (!design) return;
  state.customMutations = design.mutations.map((item) => ({ position: item.position, wt: item.wt, mutant: item.mutant }));
  $("custom-name").value = `${design.id} edited`;
  renderCustomBuilder();
  setInspectorPanel("design");
}

function renderCustomBuilder() {
  const container = $("custom-mutations");
  if (!state.customMutations.length) {
    container.className = "custom-mutations empty-state";
    container.textContent = "从位点详情或区域设计中添加突变。";
    $("save-custom").disabled = true;
    return;
  }
  container.className = "custom-mutations";
  container.innerHTML = state.customMutations.map((mutation) => {
    const options = AA.filter((aa) => aa !== mutation.wt).map((aa) =>
      `<option value="${aa}" ${aa === mutation.mutant ? "selected" : ""}>${aa}</option>`
    ).join("");
    return `<div class="custom-row"><code>${h(mutation.wt)}${mutation.position}</code><select data-custom-select="${mutation.position}" aria-label="${mutation.wt}${mutation.position} 替换氨基酸">${options}</select><button type="button" data-custom-remove="${mutation.position}">移除</button></div>`;
  }).join("");
  $("save-custom").disabled = state.demo || !state.jobId;
  $("custom-feedback").textContent = state.demo ? "冻结示例可编辑预览；提交真实作业后可保存导出。" : "";
  document.querySelectorAll("[data-custom-select]").forEach((select) => {
    select.addEventListener("change", () => {
      const mutation = state.customMutations.find((item) => item.position === Number(select.dataset.customSelect));
      if (mutation) mutation.mutant = select.value;
    });
  });
  document.querySelectorAll("[data-custom-remove]").forEach((button) => {
    button.addEventListener("click", () => {
      state.customMutations = state.customMutations.filter((item) => item.position !== Number(button.dataset.customRemove));
      renderCustomBuilder();
    });
  });
}

async function saveCustomDesign() {
  if (!state.jobId || !state.customMutations.length) return;
  $("save-custom").disabled = true;
  $("custom-feedback").textContent = "正在验证并生成 ECD 序列…";
  try {
    const response = await api(`/api/jobs/${state.jobId}/designs`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        name: $("custom-name").value,
        mutations: state.customMutations.map(({ position, mutant }) => ({ position, mutant })),
        force_high_risk: $("force-risk").checked,
      }),
    });
    state.result = response.job.result;
    state.selectedDesign = response.design.id;
    renderAll();
    $("custom-feedback").textContent = `${response.design.id} 已保存并进入导出文件。`;
    showToast("自定义构建已保存");
  } catch (error) {
    const issues = error.detail?.issues;
    $("custom-feedback").textContent = issues ? issues.join("；") : error.message;
  } finally {
    $("save-custom").disabled = false;
  }
}

function renderExports() {
  const container = $("export-links");
  const files = Object.keys(state.result.files || {});
  if (!files.length) {
    container.className = "export-links empty-state";
    container.textContent = "当前结果没有可下载文件。";
    return;
  }
  const labels = {
    "ecd_epitope_design.zip": "完整结果 ZIP", "analysis.json": "分析 JSON",
    "residues.csv": "逐残基 CSV", "designs.csv": "构建 CSV",
    "sequences.fasta": "氨基酸 FASTA", "coding_sequences.fasta": "核苷酸 FASTA",
    "sequences.gb": "GenBank", "methods_report.md": "方法报告",
    "run_manifest.json": "复现 manifest", "mapping.csv": "编号映射 CSV",
  };
  const order = ["ecd_epitope_design.zip", "methods_report.md", "designs.csv", "sequences.fasta", "coding_sequences.fasta", "sequences.gb", "residues.csv", "mapping.csv", "analysis.json", "run_manifest.json"];
  const sorted = [...files].sort((a, b) => order.indexOf(a) - order.indexOf(b));
  container.className = "export-links";
  container.innerHTML = sorted.map((file) => {
    const url = state.demo ? `/assets/demo/${file}` : `/api/jobs/${state.jobId}/files/${file}`;
    return `<a href="${h(url)}" download>${h(labels[file] || file)}</a>`;
  }).join("");
}

function renderAll() {
  if (!state.result) return;
  if (!state.result.protein.ecd_coding_sequence) {
    state.result.protein.ecd_coding_sequence = [...state.result.protein.ecd_sequence].map((aa) => CODONS[aa]).join("");
  }
  state.selectedRegion = null;
  if (!state.result.designs.some((item) => item.id === state.selectedDesign)) state.selectedDesign = "wt";
  renderSummary();
  renderStageRail();
  renderSource();
  renderMethods();
  renderLayerControls();
  renderResidueList();
  renderResidueDetail();
  renderRegions();
  renderDesigns();
  renderDesignSelector();
  renderSequences();
  renderCustomBuilder();
  renderExports();
  const structureUrl = state.result.structure_url || "/assets/demo/structure.cif";
  ensureViewer(structureUrl);
}

async function loadDemo() {
  try {
    setJobState("RUNNING", "载入冻结 PD-L1 审查结果", 0.4);
    const result = await api("/assets/demo/analysis.json");
    result.structure_url = "/assets/demo/ecd_structure.cif";
    state.result = result;
    state.previewPlans.clear();
    state.jobId = null;
    state.demo = true;
    state.selectedResidue = result.residues.find((item) => item.folding_screen_pass) || result.residues[0];
    state.selectedDesign = result.designs.find((item) => item.mutations.length)?.id || "wt";
    setJobState("DEMO", "PD-L1 冻结示例；提交真实作业后可保存自定义组合", 1);
    renderAll();
  } catch (error) {
    setJobState("FAILED", "示例载入失败", 0, error.message);
  }
}

function watchJob(jobId) {
  state.eventSource?.close();
  const source = new EventSource(`/api/jobs/${jobId}/events`);
  state.eventSource = source;
  source.onmessage = async (event) => {
    const update = JSON.parse(event.data);
    setJobState(update.status, update.message, update.progress, update.error);
    if (["COMPLETE", "FAILED", "WAITING_FOR_AF3_UPLOAD"].includes(update.status)) {
      source.close();
      const job = await api(`/api/jobs/${jobId}`);
      if (job.status === "COMPLETE") {
        state.result = job.result;
        state.previewPlans.clear();
        state.demo = false;
        state.selectedResidue = job.result.residues.find((item) => item.folding_screen_pass) || job.result.residues[0];
        state.selectedDesign = job.result.designs.find((item) => item.mutations.length)?.id || "wt";
        renderAll();
      } else if (job.status === "WAITING_FOR_AF3_UPLOAD") {
        $("af3-dialog").showModal();
      }
    }
  };
  source.onerror = () => {
    if (source.readyState === EventSource.CLOSED) return;
    setJobState("FAILED", "作业事件连接中断", 0, "请刷新任务状态或重新提交。" );
  };
}

async function submitRun(event) {
  event.preventDefault();
  const inputType = document.querySelector('input[name="input-type"]:checked').value;
  const payload = {
    input_type: inputType,
    identifier: $("identifier").value.trim(),
    chain: $("chain").value.trim() || null,
    strict_af3: true,
    allow_partial: false,
  };
  try {
    setJobState("QUEUED", "正在提交作业", 0.01);
    const job = await api("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    state.jobId = job.id;
    state.demo = false;
    watchJob(job.id);
  } catch (error) {
    setJobState("FAILED", "提交失败", 0, error.message);
  }
}

async function submitStructureUpload(event) {
  event.preventDefault();
  const file = $("upload-structure").files[0];
  if (!file) return;
  const form = new FormData();
  form.append("structure", file);
  for (const [name, id] of [
    ["accession", "upload-accession"], ["chain", "upload-chain"], ["sequence", "upload-sequence"],
    ["ecd_start", "upload-ecd-start"], ["ecd_end", "upload-ecd-end"],
  ]) {
    const value = $(id).value.trim();
    if (value) form.append(name, value);
  }
  try {
    setJobState("QUEUED", "正在上传结构", 0.01);
    const job = await api("/api/jobs/upload", { method: "POST", body: form });
    $("upload-dialog").close();
    state.jobId = job.id;
    state.demo = false;
    watchJob(job.id);
  } catch (error) {
    setJobState("FAILED", "结构上传失败", 0, error.message);
  }
}

async function submitAf3(event) {
  event.preventDefault();
  if (!state.jobId) return;
  const form = new FormData();
  form.append("structure", $("af3-structure").files[0]);
  form.append("confidence", $("af3-confidence").files[0]);
  form.append("official_source_confirmed", String($("af3-confirm").checked));
  try {
    const job = await api(`/api/jobs/${state.jobId}/af3`, { method: "POST", body: form });
    $("af3-dialog").close();
    setJobState(job.status, job.message, job.progress);
    watchJob(state.jobId);
  } catch (error) {
    setJobState("WAITING_FOR_AF3_UPLOAD", "AF3 文件验证失败", 0.34, error.message);
  }
}

function handleMolstarEvent(event) {
  const detail = event.eventData || event.detail || {};
  if (!state.result || detail.auth_seq_id == null) return;
  const residue = state.result.residues.find((item) =>
    item.auth_seq_id === Number(detail.auth_seq_id) &&
    (!detail.auth_asym_id || item.chain_id === detail.auth_asym_id) &&
    (!detail.ins_code || item.insertion_code === detail.ins_code)
  );
  if (residue) selectResidue(residue.uniprot_pos, false);
}

function setMobileView(view) {
  document.body.classList.remove("mobile-structure", "mobile-residues", "mobile-regions", "mobile-sequence");
  document.body.classList.add(`mobile-${view}`);
  document.querySelectorAll("[data-mobile-view]").forEach((button) => button.classList.toggle("active", button.dataset.mobileView === view));
  if (view === "residues") setInspectorPanel("residues");
  if (view === "regions") setInspectorPanel("regions");
}

function bindEvents() {
  $("run-form").addEventListener("submit", submitRun);
  $("load-demo").addEventListener("click", loadDemo);
  $("open-upload").addEventListener("click", () => $("upload-dialog").showModal());
  $("close-upload").addEventListener("click", () => $("upload-dialog").close());
  $("upload-form").addEventListener("submit", submitStructureUpload);
  $("close-af3").addEventListener("click", () => $("af3-dialog").close());
  $("af3-form").addEventListener("submit", submitAf3);
  $("viewer").addEventListener("PDB.molstar.click", handleMolstarEvent);
  $("viewer").addEventListener("PDB.molstar.mouseover", handleMolstarEvent);
  document.querySelectorAll(".primary-pages button").forEach((button) => button.addEventListener("click", () => setPrimaryPage(button.dataset.primary)));
  document.querySelectorAll(".secondary-nav button").forEach((button) => button.addEventListener("click", () => setInspectorPanel(button.dataset.panel)));
  document.querySelectorAll("[data-representation]").forEach((button) => button.addEventListener("click", () => setRepresentation(button.dataset.representation)));
  document.querySelectorAll("[data-engine]").forEach((button) => button.addEventListener("click", () => setRenderEngine(button.dataset.engine)));
  $("btn-lang").addEventListener("click", toggleLang);
  $("export-offline").addEventListener("click", () => {
    if (!state.jobId || state.demo) { showToast(t("msg.exportNeedJob", "提交真实作业后可导出")); return; }
    window.open(`/api/jobs/${state.jobId}/export/offline?lang=${state.lang}`, "_blank");
  });
  const copyText = (text) => {
    if (!text) return;
    if (navigator.clipboard?.writeText) navigator.clipboard.writeText(text).then(() => showToast(t("msg.copied", "已复制")), () => showToast(t("msg.copyFail", "复制失败")));
    else showToast(t("msg.copyFail", "复制失败"));
  };
  $("copy-wt-aa").addEventListener("click", () => copyText(state.lastSequences?.wtAa));
  $("copy-wt-nt").addEventListener("click", () => copyText(state.lastSequences?.wtDna));
  $("copy-mut-aa").addEventListener("click", () => copyText(state.lastSequences?.mutAa));
  $("copy-mut-nt").addEventListener("click", () => copyText(state.lastSequences?.mutDna));
  $("residue-filter").addEventListener("change", renderResidueList);
  $("residue-search").addEventListener("input", renderResidueList);
  $("sequence-design").addEventListener("change", () => selectDesign($("sequence-design").value));
  document.querySelectorAll('input[name="sequence-mode"]').forEach((input) => input.addEventListener("change", () => {
    state.sequenceMode = input.value;
    renderSequences();
  }));
  document.querySelectorAll('input[name="input-type"]').forEach((input) => input.addEventListener("change", () => {
    const uniprot = input.value === "uniprot" && input.checked;
    $("identifier-label").textContent = uniprot ? "蛋白 accession" : "PDB ID";
    $("identifier").value = uniprot ? "Q9NZQ7" : "4ZQK";
  }));
  $("reset-view").addEventListener("click", () => state.viewer?.visual.reset({ camera: true }));
  $("spin-view").addEventListener("click", async () => {
    state.spinning = !state.spinning;
    $("spin-view").setAttribute("aria-pressed", String(state.spinning));
    await state.viewer?.visual.toggleSpin(state.spinning, false);
  });
  $("clear-focus").addEventListener("click", async () => {
    state.selectedResidue = null;
    await state.viewer?.visual.clearHighlight();
    await state.viewer?.visual.reset({ camera: true });
    renderResidueList(); renderResidueDetail(); renderSequences();
  });
  $("clear-custom").addEventListener("click", () => { state.customMutations = []; renderCustomBuilder(); });
  $("save-custom").addEventListener("click", saveCustomDesign);
  document.querySelectorAll("[data-mobile-view]").forEach((button) => button.addEventListener("click", () => setMobileView(button.dataset.mobileView)));
}

/* ---- beads engine (3Dmol) wiring; BeadsRenderer itself lives in renderer.js.
 * These must stay inside the module: they close over module-scope state. */
function setRenderEngine(engine) {
  if (!state.result) return;
  if (engine === "beads" && !state.jobId) {
    showToast(t("msg.beadsNeedJob", "冻结示例无作业结构文件；提交真实作业后可用珠子渲染"));
    document.querySelectorAll("[data-engine]").forEach((button) => {
      button.setAttribute("aria-pressed", String(button.dataset.engine === "molstar"));
    });
    return;
  }
  state.engine = engine;
  document.querySelectorAll("[data-engine]").forEach((button) => {
    button.setAttribute("aria-pressed", String(button.dataset.engine === engine));
  });
  const molstarHost = $("viewer");
  const beadsHost = $("beads-viewer");
  if (!beadsHost) return;
  if (engine === "beads") {
    molstarHost.style.visibility = "hidden";
    beadsHost.hidden = false;
    if (!state.beadsReady) {
      BeadsRenderer.init(
        beadsHost,
        `/api/jobs/${state.jobId}/files/structure_for_energy.pdb`,
        state.result.residues
      ).then(() => {
        state.beadsReady = true;
        BeadsRenderer.setPickCallbacks(
          (residue) => { if (residue) selectResidue(residue.uniprot_pos, true); },
          (residue) => { if (residue) setInspectorPanel("residues"); }
        );
        paintBeads();
      }).catch((error) => {
        console.error(error);
        showToast((t("msg.beadsFail", "珠子渲染载入失败：")) + error.message);
        setRenderEngine("molstar");
      });
    } else {
      paintBeads();
    }
  } else {
    beadsHost.hidden = true;
    molstarHost.style.visibility = "visible";
    BeadsRenderer.destroy();
    state.beadsReady = false;
  }
}

/* Colour semantics for beads.  Mirrors applyViewerLayer():
 * design/custom > interface > region colour > folding pass > surface > grey. */
function beadsColorFor(residue) {
  const design = selectedDesign();
  const designPositions = new Set((design?.mutations || []).map((m) => m.position));
  if (state.customMutations.some((m) => m.position === residue.uniprot_pos)) return "#8a2f1d";
  if (designPositions.has(residue.uniprot_pos)) return STEP_META.design.color;
  if (state.layers.has("interface") && residue.interface) return STEP_META.interface.color;
  if (state.layers.has("regions") && residue.region_id != null) {
    const region = state.result.regions.find((item) => item.id === residue.region_id);
    if (region) return region.color;
  }
  if (state.layers.has("folding") && residue.folding_screen_pass) return STEP_META.folding.color;
  if (state.layers.has("surface") && residue.surface_accessible) return STEP_META.surface.color;
  return "#7c909d";
}

function paintBeads() {
  if (!state.beadsReady) return;
  const selected = new Set(
    state.customMutations.map((m) => m.position)
      .concat((selectedDesign()?.mutations || []).map((m) => m.position))
  );
  const sidechain = new Set(
    [...selected].concat(
      state.layers.has("interface")
        ? state.result.residues.filter((r) => r.interface).map((r) => r.uniprot_pos)
        : []
    )
  );
  BeadsRenderer.setColors(beadsColorFor, { selected, sidechain });
}

(async function boot() {
  await loadI18n();
  applyStaticI18n();
  bindEvents();
  renderLegend();
  setMobileView("structure");
  loadDemo();
})();
