import { createRecommendationRenderer } from "./recommendations.js?v=regional-qualified-20261025-v1";

/* ECD Epitope Workbench — frontend (rewrite v4: full-page i18n + interaction fixes)
 *
 * Contract with the backend (ecd_epitope_scanner api.py — unchanged):
 *   POST /api/jobs            {input_type, identifier, chain, strict_af3, allow_partial}
 *   POST /api/jobs/upload     multipart(structure, analysis_target, accession?, chain?, sequence?, ecd_start?, ecd_end?)
 *   GET  /api/jobs/{id}       -> {status, message, progress, result}
 *   GET  /api/jobs/{id}/events        SSE progress
 *   POST /api/jobs/{id}/af3           multipart(structure, confidence, official_source_confirmed)
 *   POST /api/jobs/{id}/designs       {name, mutations[{position,mutant}], force_high_risk}
 *   GET  /api/jobs/{id}/files/{name}  result files
 *   GET  /api/jobs/{id}/export/offline?lang=zh|en
 *   GET  /assets/*                    static frontend + frozen demo
 *
 * Honest-scoring rule (requirement #3): the folding-tolerance step renders
 * ONLY real engine numbers.  `score_origin` (written by the backend) decides:
 *   foldx / evoef2 -> kcal/mol values are shown;
 *   proxy / absent -> a status card explains that no energy score exists.
 * Proxy ranking numbers are never rendered as scores anywhere.
 *
 * Interaction rules (v4):
 *   - Mol* pick is CLICK-ONLY.  Binding mouseover to the pick handler used to
 *     rebuild the residue-detail panel on every hover, which tore the custom-
 *     design selects out of the DOM mid-interaction ("can't pick a residue").
 *   - renderResidueDetail() skips rebuilds when nothing it displays changed,
 *     so stray re-renders can never fight an open <select>.
 *   - Adding a custom mutation stays on the current step (toast instead of a
 *     forced jump to step 6).
 *   - The whole page is i18n: static shell via data-i18n*, dynamic panels via
 *     t(); the panel scaffold rebinds its own controls so a language switch
 *     can rebuild it safely.
 */

"use strict";

const $ = (id) => document.getElementById(id);
const AA = "ACDEFGHIKLMNPQRSTVWY".split("");
const CODONS = {
  A: "GCC", R: "CGC", N: "AAC", D: "GAC", C: "TGC", Q: "CAG", E: "GAG",
  G: "GGC", H: "CAC", I: "ATC", L: "CTG", K: "AAG", M: "ATG", F: "TTC",
  P: "CCC", S: "AGC", T: "ACC", W: "TGG", Y: "TAC", V: "GTG",
};

const state = {
  result: null,
  jobId: null,
  demo: false,
  eventSource: null,
  lang: "zh",
  i18n: {},
  step: 1,
  engine: "molstar",
  beadsReady: false,
  // Surface-first display: the primary scientific object is the solvent-
  // exposed ECD surface, while cartoon/stick remain explicit alternatives.
  representation: "molecular-surface",
  spinning: false,
  layers: new Set(["surface", "regions", "design"]),
  selectedFootprint: null,
  selectedRegion: null,
  layerFocus: "folding",
  selectedResidue: null,
  selectedDesign: "wt",
  recommendationSelections: new Map(),
  resultJobId: null,
  showPartner: false,
  designSurfaceMode: false,
  designFocusRequest: 0,
  previewMode: "construct", // construct | custom; custom never shadows a clicked design
  previewPlans: new Map(),
  customMutations: [],
  lastSequences: null,
  offline: Boolean(window.APP_OFFLINE),
};

const STEP_LIST = ["step.surface", "step.tolerance", "step.regions", "step.interface", "step.design", "step.custom", "step.overview"];
const LAYER_META = {
  surface:   { labelKey: "layer.surface", color: "#1565C0" },
  folding:   { labelKey: "layer.folding", color: "#2E7D32" },
  interface: { labelKey: "layer.interface", color: "#E4572E" },
  regions:   { labelKey: "layer.regions", color: "#6A4C93" },
  design:    { labelKey: "layer.design", color: "#E4572E" },
};

/* ---------------- tiny utils ---------------- */

const h = (v) => String(v ?? "")
  .replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;")
  .replaceAll('"', "&quot;").replaceAll("'", "&#039;");

const fmt = (v, d = 2) => (v == null ? "—" : Number(v).toFixed(d));
const fmtPct = (v) => `${Math.round(Number(v || 0) * 100)}%`;

function t(key, fallback) {
  const entry = state.i18n?.[state.lang]?.[key] ?? state.i18n?.zh?.[key] ?? state.i18n?.en?.[key];
  if (entry == null) {
    if (fallback != null) return fallback;
    console.warn("i18n missing key:", key);
    return key;
  }
  return entry;
}

async function api(url, options = {}) {
  if (state.offline) throw new Error("离线结果包不连接服务器，也不执行分析计算。" );
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

async function loadI18n() {
  if (state.offline) {
    state.i18n = window.APP_OFFLINE.i18n || {};
    const saved = (() => { try { return localStorage.getItem("ecd.lang"); } catch { return null; } })();
    state.lang = saved === "en" || saved === "zh" ? saved : (window.APP_OFFLINE.language || "zh");
    return;
  }
  for (const lang of ["zh", "en"]) {
    try { state.i18n[lang] = await api(`/assets/i18n/${lang}.json?v=job-eta-20261026-v1`); }
    catch (e) { console.error(`i18n load failed for ${lang}`, e); state.i18n[lang] = {}; }
  }
  const saved = localStorage.getItem("ecd.lang");
  state.lang = saved === "en" || saved === "zh" ? saved : "zh";
}

function applyStaticI18n() {
  document.documentElement.lang = state.lang === "en" ? "en" : "zh-CN";
  document.title = t("app.title", "ECD Epitope Workbench");
  document.querySelectorAll("[data-i18n]").forEach((n) => { n.textContent = t(n.dataset.i18n, n.textContent); });
  document.querySelectorAll("[data-i18n-placeholder]").forEach((n) => {
    n.placeholder = t(n.dataset.i18nPlaceholder, n.placeholder);
  });
  document.querySelectorAll("[data-i18n-title]").forEach((n) => {
    n.title = t(n.dataset.i18nTitle, n.title);
  });
  document.querySelectorAll("[data-i18n-aria]").forEach((n) => {
    n.setAttribute("aria-label", t(n.dataset.i18nAria, n.getAttribute("aria-label") || ""));
  });
  $("btn-lang").textContent = state.lang === "zh" ? "EN" : "中文";
}

function toggleLang() {
  state.lang = state.lang === "zh" ? "en" : "zh";
  try { localStorage.setItem("ecd.lang", state.lang); } catch { /* file:// storage may be disabled */ }
  applyStaticI18n();
  renderPanelScaffold();   // rebuilds the six step sections in the new language
  if (state.result) renderAll(false);
  else renderLegend();
  renderSeqDockHead();
  updateUploadTargetUI();
}

function showToast(message) {
  const el = $("toast");
  el.textContent = message;
  el.hidden = false;
  clearTimeout(showToast.timer);
  showToast.timer = setTimeout(() => { el.hidden = true; }, 4000);
}

/* ---------------- honest scoring (requirement #3) ---------------- */

function scoreOrigin(design) {
  if (!design) return "none";
  if (design.score_origin) return design.score_origin;
  const real = design.stability_engine && /foldx|evoef2/i.test(design.stability_engine)
    && design.combined_stability_score != null;
  return real ? "foldx" : "proxy";
}

/* Returns a display contract for one design's combination score.
   `real=false` means: render a STATUS, never a number. */
function scoreInfo(design) {
  const origin = scoreOrigin(design);
  if (origin === "foldx" || origin === "evoef2") {
    return {
      real: true, origin,
      value: design.combined_stability_score,
      engine: design.stability_engine,
      label: origin === "foldx" ? t("score.foldx", "FoldX 5.1 ΔΔG（真实计算值）") : t("score.evoef2", "EvoEF2 能量（真实计算值）"),
    };
  }
  return {
    real: false, origin: "proxy", value: null,
    engine: design?.stability_engine || "builtin structural-risk proxy",
    label: t("score.noReal", "未获得真实能量评分"),
  };
}

function singleChoiceScore(choice) {
  /* per-substitution: only FoldX numbers count as scores */
  return choice.foldx_ddg != null ? { real: true, value: choice.foldx_ddg } : { real: false, value: null };
}

/* ---------------- job lifecycle ---------------- */

function setJob(status, message, progress = 0, error = null, estimate = null) {
  const names = {
    QUEUED: t("status.QUEUED", "排队中"), RUNNING: t("status.RUNNING", "运行中"),
    COMPLETE: t("status.COMPLETE", "完成"), WAITING_FOR_AF3_UPLOAD: t("status.WAITING", "等待 AF3"),
    FAILED: t("status.FAILED", "失败"), DEMO: t("status.DEMO", "冻结示例"),
  };
  const cls = {
    QUEUED: "running", RUNNING: "running", COMPLETE: "complete",
    WAITING_FOR_AF3_UPLOAD: "waiting", FAILED: "failed", DEMO: "",
  };
  const badge = $("job-state");
  badge.textContent = names[status] || status || t("status.idle", "未运行");
  badge.className = `badge ${cls[status] || ""}`;
  $("job-message").textContent = message || "";
  const percent = Math.max(0, Math.min(100, Math.round((progress || 0) * 100)));
  $("job-percent").textContent = `${percent}%`;
  $("job-progress").style.width = `${percent}%`;
  $("job-error").hidden = !error;
  $("job-error").textContent = error || "";
  $("run-button").disabled = status === "RUNNING" || status === "QUEUED";
  $("btn-export-offline").disabled = status !== "COMPLETE";
  const eta = $("job-eta");
  eta.hidden = !["RUNNING", "QUEUED"].includes(status);
  const stages = {
    V4_FOLDX_FIRST_PASS: ["FoldX首轮", "FoldX first pass"], V4_FOLDX_REFINEMENT: ["FoldX追加确认", "FoldX refinement"],
    V4_SPURS_SINGLE: ["单点评分", "Single-mutation scoring"], V4_SPURS_COMBINATION: ["组合评分", "Combination scoring"],
    V4_COMBINATION_SEARCH: ["组合搜索", "Combination search"], EXPORTING: ["导出结果", "Exporting"],
  };
  const stage = stages[estimate?.stage]?.[state.lang === "en" ? 1 : 0];
  eta.textContent = status === "QUEUED" ? t("eta.queued", "排队中，开始后估时")
    : estimate?.status === "ESTIMATED" ? (() => {
      const [low, high] = estimate.range_seconds.map((seconds) => Math.ceil(seconds / 60));
      const range = high <= 1 ? t("eta.underMinute", "不足1分钟") : low === high ? `${high}` : `${low}–${high}`;
      return `${t("eta.remaining", "预计剩余")} ${range}${high <= 1 ? "" : ` ${t("eta.minutes", "分钟")}`}${stage ? ` · ${stage}` : ""}`;
    })() : estimate?.status === "OVERRUN" ? t("eta.overrun", "已超出历史耗时范围，剩余时间待更新")
    : t("eta.pending", "正在收集耗时信息，暂无法估计剩余时间");
  eta.title = `${t("eta.basis", "参考同类已完成任务的阶段耗时，非完成承诺")} · ${estimate?.reference_count || 0} ${t("eta.references", "个参考任务")}`;
}

function adoptResult(result, { demo, jobId }) {
  if (state.resultJobId !== jobId) state.recommendationSelections?.clear();
  state.resultJobId = jobId;
  state.result = result;
  state.showPartner = false;
  state.selectedRegion = null;
  state.selectedFootprint = null;
  state.reliabilityView = false;
  state.designSurfaceMode = false;
  state.designFocusRequest += 1;
  state.previewPlans.clear();
  state.demo = demo;
  state.jobId = jobId;
  state.selectedResidue = result.residues.find((r) => r.folding_screen_pass) || result.residues[0] || null;
  state.selectedDesign = "wt";
  state.previewMode = "construct";
  if (!result.protein.ecd_coding_sequence) {
    result.protein.ecd_coding_sequence = [...result.protein.ecd_sequence].map((aa) => CODONS[aa]).join("");
  }
  renderAll(false);
  ensureViewers(result.structure_url || (state.offline ? "" : "/assets/demo/ecd_structure.cif"));
}

function taskName() {
  const protein = state.result?.protein || {};
  return protein.protein_name
    ? `${protein.protein_name}${protein.accession ? ` (${protein.accession})` : ""}`
    : (protein.accession || "ECD Epitope Workbench");
}

function applyTaskIdentity() {
  const title = taskName();
  $("app-brand-title").textContent = "ECD Epitope Workbench";
  $("workspace-title").textContent = state.result ? title : t("hud.canvasTitle", "结构审查画布");
  $("workspace-context").textContent = state.result
    ? `${state.result.protein.ecd_start}–${state.result.protein.ecd_end} · ${state.result.regions.length} ${t("workspace.regions", "分区")}`
    : t("workspace.loading", "正在载入结构");
  document.title = state.result ? `${title} | ECD Epitope Workbench` : t("app.title", "ECD Epitope Workbench");
}

let resultLoadId = 0;

function loadDemo() {
  setJob("READY", "", 0);
}

async function loadSavedJob(jobId) {
  const loadId = ++resultLoadId;
  try {
    const job = await api(`/api/jobs/${encodeURIComponent(jobId)}?view=display`);
    if (loadId !== resultLoadId) return;
    setJob(job.status, job.message || "", job.progress || 0, job.error || "", job.time_estimate);
    if (job.status === "COMPLETE" && job.result) {
      adoptResult(job.result, { demo: false, jobId });
    } else if (job.status === "RUNNING" || job.status === "QUEUED") {
      state.jobId = jobId;
      watchJob(jobId);
    }
  } catch (error) {
    if (loadId !== resultLoadId) return;
    console.error(error);
    showToast(t("msg.savedJobFail", "无法载入指定任务：") + error.message);
    loadDemo();
  }
}

function watchJob(jobId) {
  resultLoadId += 1;
  state.eventSource?.close();
  const source = new EventSource(`/api/jobs/${jobId}/events`);
  state.eventSource = source;
  source.onmessage = async (event) => {
    if (state.eventSource !== source) return;
    try {
      const u = JSON.parse(event.data);
      setJob(u.status, u.message, u.progress, u.error, u.time_estimate);
      if (["COMPLETE", "FAILED", "WAITING_FOR_AF3_UPLOAD"].includes(u.status)) {
        source.close();
        const job = await api(`/api/jobs/${jobId}?view=display`);
        if (state.eventSource !== source) return;
        if (job.status === "COMPLETE") {
          adoptResult(job.result, { demo: false, jobId });
        } else if (job.status === "WAITING_FOR_AF3_UPLOAD") {
          $("af3-dialog").showModal();
        }
      }
    } catch (e) {
      console.error("bad SSE event", e);
    }
  };
  source.onerror = () => {
    if (state.eventSource !== source) return;
    if (source.readyState === EventSource.CLOSED) return;
    setJob("FAILED", t("msg.sseLost", "作业事件连接中断"), 0, t("msg.sseLostHint", "请刷新任务状态或重新提交。"));
  };
}

function panelBudgetInput() {
  const minimum = Number($("budget-minimum")?.value || 4);
  const target = Number($("budget-target")?.value || 7);
  const maximum = Number($("budget-maximum")?.value || 8);
  const ceiling = 8;
  if (![minimum, target, maximum].every(Number.isInteger) || minimum < 1 || minimum > target || target > maximum || maximum > ceiling) {
    throw new Error(`构建体预算需要满足 1 ≤ 最少目标 ≤ 常见目标 ≤ 硬上限 ≤ ${ceiling}`);
  }
  return {minimum_constructs: minimum, target_constructs: target,
          maximum_constructs: maximum, adaptive_budget: $("adaptive-budget")?.checked ?? true};
}

function displayPartitionInput() {
  const rawTarget = $("sector-target").value.trim();
  const target = rawTarget ? Number(rawTarget) : null;
  const geodesic = Number($("sector-geodesic").value);
  const ca = Number($("sector-ca").value);
  if ((target !== null && (!Number.isInteger(target) || target < 7 || target > 15)) || ![geodesic, ca].every((value) => Number.isFinite(value) && value > 0)) {
    throw new Error(t("form.partitionInvalid", "分区目标留空自动选择，或填写7–15之间的整数；两个直径上限必须为正数。"));
  }
  return { ...(target === null ? {} : {target_sector_residues: target}), maximum_sector_geodesic_diameter_a: geodesic, maximum_sector_ca_diameter_a: ca };
}

async function submitRun(event) {
  event.preventDefault();
  const inputType = document.querySelector('input[name="input-type"]:checked').value;
  const designEngine = "joint_v4_ml";
  try {
    $("run-feedback").hidden = true;
    setJob("QUEUED", t("msg.submitting", "提交作业…"), 0.02);
    const job = await api("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        input_type: inputType,
        identifier: $("identifier").value.trim(),
        chain: $("chain").value.trim() || null,
        strict_af3: true,
        allow_partial: false,
        design_engine: designEngine,
        ml_mode: $("ml-mode").value,
        max_sectors: 12,
        ...displayPartitionInput(),
        ...panelBudgetInput(),
      }),
    });
    state.jobId = job.id;
    state.demo = false;
    $("run-dialog").close();
    setJob(job.status, job.message, job.progress, job.error, job.time_estimate);
    watchJob(job.id);
  } catch (error) {
    setJob("FAILED", t("msg.submitFail", "作业提交失败"), 0, error.message);
    $("run-feedback").hidden = false;
    $("run-feedback").textContent = error.message;
  }
}

async function submitUpload(event) {
  event.preventDefault();
  const folderMode = $("upload-source-mode").value === "af3-directory";
  const structure = $("upload-structure").files[0];
  const folderFiles = [...$("upload-af3-directory").files];
  if (folderMode ? (!folderFiles.length || !$("upload-af3-confirm").checked) : !structure) return;
  $("upload-dialog").close();
  const form = new FormData();
  const analysisTarget = $("upload-analysis-target").value;
  if (folderMode) {
    folderFiles.forEach((file) => form.append("files", file, file.webkitRelativePath || file.name));
    form.append("official_source_confirmed", "true");
  } else form.append("structure", structure);
  form.append("analysis_target", analysisTarget);
  if ($("upload-accession").value.trim()) form.append("accession", $("upload-accession").value.trim());
  if ($("upload-chain2").value.trim()) form.append("chain", $("upload-chain2").value.trim());
  if (analysisTarget === "antigen") {
    if ($("upload-sequence").value.trim()) form.append("sequence", $("upload-sequence").value.trim());
    if ($("upload-ecd-start").value) form.append("ecd_start", $("upload-ecd-start").value);
    if ($("upload-ecd-end").value) form.append("ecd_end", $("upload-ecd-end").value);
  }
  form.append("design_engine", "joint_v4_ml");
  form.append("ml_mode", $("upload-ml-mode").value);
  form.append("max_sectors", "12");
  try {
    Object.entries(displayPartitionInput()).forEach(([key, value]) => form.append(key, String(value)));
    Object.entries(panelBudgetInput()).forEach(([key, value]) => form.append(key, String(value)));
    setJob("QUEUED", t("msg.uploading", "上传结构并创建作业…"), 0.02);
    const job = await api(folderMode ? "/api/jobs/af3-directory" : "/api/jobs/upload", { method: "POST", body: form });
    state.jobId = job.id;
    state.demo = false;
    setJob(job.status, job.message, job.progress, job.error, job.time_estimate);
    watchJob(job.id);
  } catch (error) {
    setJob("FAILED", t("msg.uploadFail", "上传失败"), 0, error.message);
  }
}

function updateUploadTargetUI() {
  const folderMode = $("upload-source-mode").value === "af3-directory";
  $("upload-structure-field").hidden = folderMode;
  $("upload-structure").required = !folderMode;
  $("upload-af3-directory-field").hidden = !folderMode;
  $("upload-af3-directory").required = folderMode;
  $("upload-af3-confirm-field").hidden = !folderMode;
  $("upload-af3-confirm").required = folderMode;
  const binder = $("upload-analysis-target").value === "binder";
  $("upload-target-hint").textContent = t(
    binder ? "upload.binderHint" : "upload.hint",
    binder ? "Analyze the binder surface. Monomers use their only chain automatically; for complexes, enter the binder chain ID. Antigen membrane orientation is not applied, and B-factors are not assumed to be pLDDT." : "Only an antigen ECD PDB/mmCIF is required."
  );
  $("upload-accession-label").textContent = t(binder ? "upload.accessionBinder" : "upload.accession", binder ? "UniProt accession (optional)" : "UniProt accession");
  $("upload-chain-label").textContent = t(binder ? "upload.chainBinder" : "upload.chain", binder ? "Binder chain (optional)" : "Target chain");
  $("upload-chain2").placeholder = t(binder ? "upload.chainBinderHint" : "form.chainPlaceholder", binder ? "e.g. A; blank for a monomer" : "auto");
  $("upload-accession").placeholder = t(binder ? "upload.accessionBinderHint" : "upload.accessionHint", binder ? "optional; blank is OK" : "recommended");
  $("upload-antigen-advanced").hidden = binder;
}

async function submitAf3(event) {
  event.preventDefault();
  if ($("af3-source-mode").value === "directory") {
    const files = [...$("af3-directory").files];
    if (!files.length || !$("af3-confirm").checked) {
      toast(t("msg.af3FilesRequired", "请选择官方 AF3 文件夹并确认来源。"));
      return;
    }
    const form = new FormData();
    files.forEach((file) => form.append("files", file, file.webkitRelativePath || file.name));
    form.append("official_source_confirmed", "true");
    await submitAf3Form(form, "/af3-directory");
    return;
  }
  const structure = $("af3-structure").files[0];
  const confidence = $("af3-confidence").files[0];
  if (!structure || !confidence || !$("af3-confirm").checked) {
    toast(t("msg.af3FilesRequired", "请选择结构、置信度文件并确认来源。"));
    return;
  }
  const form = new FormData();
  form.append("structure", structure);
  form.append("confidence", confidence);
  // FastAPI declares this as a required Form field.  The checkbox alone is
  // browser state; it must be explicitly included in multipart data.
  form.append("official_source_confirmed", "true");
  await submitAf3Form(form, "/af3");
}

async function submitAf3Form(form, endpoint) {
  const submit = $("af3-submit");
  const originalLabel = submit.textContent;
  try {
    submit.disabled = true;
    submit.textContent = t("msg.af3Uploading", "正在上传…");
    const job = await api(`/api/jobs/${state.jobId}${endpoint}`, { method: "POST", body: form });
    $("af3-dialog").close();
    setJob(job.status, job.message, job.progress, job.error, job.time_estimate);
    watchJob(state.jobId);
  } catch (error) {
    setJob("WAITING_FOR_AF3_UPLOAD", t("msg.af3Fail", "AF3 文件验证失败"), 0.34, error.message);
    toast(`${t("msg.af3Fail", "AF3 文件验证失败")}: ${error.message}`);
  } finally {
    submit.disabled = false;
    submit.textContent = originalLabel;
  }
}

function updateAf3SourceUI() {
  const directory = $("af3-source-mode").value === "directory";
  $("af3-pair-fields").hidden = directory;
  $("af3-structure").required = !directory;
  $("af3-confidence").required = !directory;
  $("af3-directory-field").hidden = !directory;
  $("af3-directory").required = directory;
}

/* ---------------- viewport (3D) ---------------- */

function residueQuery(residue) {
  // ECD export rebuilds label_seq_id from 1; records retain the input AF3
  // labels. Author IDs and insertion codes survive the export unchanged.
  if (residue.auth_seq_id == null || !residue.chain_id) {
    return { label_asym_id: residue.label_chain_id, label_seq_id: residue.label_seq_id };
  }
  const q = { auth_asym_id: residue.chain_id, auth_seq_id: residue.auth_seq_id };
  if (residue.insertion_code) q.pdbx_PDB_ins_code = residue.insertion_code;
  return q;
}

let molstarReady = false;
let viewerLoadId = 0;

function viewerError(error, stage, blocking = false) {
  const labels = {
    fetch: t("hud.stageFetch", "获取结构文件"), parse: t("hud.stageParse", "检查结构格式"),
    render: t("hud.stageRender", "解析并建立三维结构"), structure: t("hud.stageStructure", "等待结构就绪"),
    webgl: t("hud.stageWebgl", "初始化图形画布"), paint: t("hud.stagePaint", "分区着色"),
    focus: t("hud.stageFocus", "定位结构"), tooltip: t("hud.stageTooltip", "安装位点提示"),
  };
  const detail = `${labels[error.stage || stage] || labels.render}：${error.message}`;
  console.error("ECD structure viewer", { jobId: state.jobId, structureUrl: state.result?.structure_url, stage: error.stage || stage, error });
  if (!blocking) { showToast(detail); return; }
  $("viewport-empty").classList.remove("hidden");
  $("viewport-empty").querySelector("strong").textContent = t("hud.canvasFail", "结构画布载入失败");
  $("viewport-empty").querySelector("span").textContent = detail;
  $("load-demo").hidden = true;
  $("retry-viewer").hidden = false;
}

async function ensureViewers(structureUrl) {
  const loadId = ++viewerLoadId;
  molstarReady = false;
  state.beadsReady = false;
  BeadsView.destroy();
  if (state.engine === "beads") setEngine("molstar");
  $("viewport-empty").classList.add("hidden");
  $("retry-viewer").hidden = true;
  let stage = "render";
  try {
    // Keep partner chains available on initial load and retry; scope the target below.
    if (interfacePartnerChains().length) structureUrl = structureUrl.replace(/\/ecd_structure\.cif$/, "/structure.cif");
    const instance = await MolstarView.init($("viewer-molstar"), structureUrl, state.representation,
      state.offline ? window.APP_OFFLINE.structure_data : null);
    if (loadId !== viewerLoadId || !instance) return;
    if (state.result.residues.length) {
      const chains = [...new Set(state.result.residues.map((residue) => residue.chain_id).filter(Boolean))];
      await MolstarView.configureChains(chains, interfacePartnerChains(), state.representation,
        state.result.residues.map(residueQuery));
    }
    molstarReady = true;
    stage = "paint";
    await paintMolstar();
    if (loadId !== viewerLoadId) return;
    stage = "focus";
    const design = state.designSurfaceMode ? selectedDesign() : null;
    const focusRecords = designSurfaceRecords(design).envelope;
    await MolstarView.focus(
      (focusRecords.length ? focusRecords : state.result.residues).map(residueQuery)
    );
    if (loadId !== viewerLoadId) return;
    try { await installMolstarTooltips(); }
    catch (error) { if (loadId === viewerLoadId) viewerError(error, "tooltip"); }
  } catch (error) {
    if (loadId !== viewerLoadId) return;
    viewerError(error, stage, !molstarReady);
  }
}

async function installMolstarTooltips() {
  if (!molstarReady || !state.result) return;
  const records = state.result.residues.map((residue) => ({
    ...residueQuery(residue),
    tooltip: `${residue.aa}${residue.uniprot_pos} | rSASA ${Number(residue.rsasa).toFixed(3)} | ${t("hud.sidechain", "侧链")} ${Number(residue.sidechain_sasa).toFixed(1)} Å² | ${residue.interface ? t("tag.interface", "蛋白界面") : t("hud.nonInterface", "非界面")}`,
  }));
  await MolstarView.tooltips(records);
}

const RELIABILITY_COLORS = {STRUCTURED: "#0072B2", FLEXIBLE_STRUCTURED: "#56B4E9", UNCERTAIN: "#E69F00", CONDITIONALLY_FOLDED: "#CC79A7", DISORDERED: "#D55E00", UNRESOLVED: "#999999"};
function reliabilityColor(residue) {
  const row = state.result?.structure_reliability?.residues?.find((r) => r.residue_id === residue.uniprot_pos);
  return RELIABILITY_COLORS[row?.state] || "#999999";
}
function beadsColorFor(residue) {
  if (state.reliabilityView) return reliabilityColor(residue);
  const selectedFootprint = state.result?.surface_v2?.footprints?.find(
    (item) => item.footprint_id === state.selectedFootprint
  );
  if (selectedFootprint?.residues?.includes(residue.uniprot_pos)) return "#D81B60";
  const design = selectedDesign();
  const designPos = new Set((design?.mutations || []).map((m) => m.position));
  if (state.previewMode === "custom" && state.customMutations.some((m) => m.position === residue.uniprot_pos)) return "#8a2f1d";
  if (designPos.has(residue.uniprot_pos)) return LAYER_META.design.color;
  if (state.layers.has("interface") && residue.interface) return LAYER_META.interface.color;
  if (state.layers.has("regions")) {
    const region = state.result.regions.find((r) => r.id === residue.region_id);
    if (region) return state.selectedRegion && state.selectedRegion !== region.id ? "#C8CDD1" : region.color;
    return "#7c909d";
  }
  if (state.layers.has("folding") && residue.folding_screen_pass) return LAYER_META.folding.color;
  if (state.layers.has("surface") && residue.surface_accessible) return LAYER_META.surface.color;
  return "#7c909d";
}

const DESIGN_SURFACE_COLORS = {
  background: "#EEECE6", region: "#D69E2E", footprint: "#F0E76A", mutation: "#D55E00",
};

function regionsForDesign(design) {
  if (!design || !state.result) return [];
  let owner = design;
  const visited = new Set();
  while (owner && owner.region_id == null && owner.parent_design_id && !visited.has(owner.id)) {
    visited.add(owner.id);
    owner = state.result.designs.find((item) => item.id === owner.parent_design_id);
  }
  if (owner?.source_engine === "joint_v4_ml" && owner.mutations?.length) {
    const positions = new Set(owner.mutations.map((mutation) => mutation.position));
    return state.result.regions.filter((region) => region.residues.some((position) => positions.has(position)));
  }
  return state.result.regions.filter((region) => region.id === owner?.region_id);
}

function regionForDesign(design) {
  return regionsForDesign(design)[0] || null;
}

function designSurfaceRecords(design = selectedDesign()) {
  const regions = regionsForDesign(design);
  const region = regions[0];
  if (!design || !region) return { region: null, envelope: [], footprint: [], mutations: [] };
  const envelopePositions = new Set(regions.flatMap((item) => item.envelope_residues?.length ? item.envelope_residues : item.residues));
  const footprintPositions = new Set(regions.flatMap((item) => item.footprint_residues?.length ? item.footprint_residues : item.residues));
  const mutationPositions = new Set((design.mutations || []).map((item) => item.position));
  return {
    region,
    envelope: state.result.residues.filter((item) => envelopePositions.has(item.uniprot_pos)),
    footprint: state.result.residues.filter((item) => footprintPositions.has(item.uniprot_pos)),
    mutations: state.result.residues.filter((item) => mutationPositions.has(item.uniprot_pos)),
  };
}

function paintBeads() {
  if (!state.beadsReady) return;
  const selected = new Set(state.previewMode === "custom"
    ? state.customMutations.map((m) => m.position)
    : (selectedDesign()?.mutations || []).map((m) => m.position));
  const sidechain = new Set(
    [...selected].concat(
      state.layers.has("interface")
        ? state.result.residues.filter((r) => r.interface).map((r) => r.uniprot_pos)
        : []
    )
  );
  BeadsView.setColors(beadsColorFor, { selected, sidechain });
}

async function paintMolstar() {
  if (!molstarReady || !state.result) return;
  const loadId = viewerLoadId;
  const data = [];
  const add = (records, color, extra = {}) => records.forEach((residue) => {
    data.push({ ...residueQuery(residue), color, ...extra });
  });
  if (state.recommendationSelections.size) {
    const byPosition = new Map();
    [...state.recommendationSelections.values()].forEach((selection, index) => selection.mutations.forEach((mutation) => {
      const variants = byPosition.get(mutation.position) || new Map();
      variants.set(`${mutation.wt}>${mutation.mutant}`, selectionColor(index));
      byPosition.set(mutation.position, variants);
    }));
    for (const [position, variants] of byPosition) {
      const residue = state.result.residues.find((item) => item.uniprot_pos === position);
      if (residue) add([residue], variants.size > 1 ? "#D81B60" : variants.values().next().value);
    }
    await MolstarView.color(data, "#c8d2d5", { keepRepresentations: true });
    if (loadId !== viewerLoadId) return;
    renderLegend();
    return;
  }
  const showcase = state.designSurfaceMode && state.previewMode === "construct" ? selectedDesign() : null;
  if (showcase?.mutations?.length) {
    const records = designSurfaceRecords(showcase);
    add(records.envelope, DESIGN_SURFACE_COLORS.region);
    add(records.footprint, DESIGN_SURFACE_COLORS.footprint);
    add(records.mutations, DESIGN_SURFACE_COLORS.mutation);
    await MolstarView.color(data, DESIGN_SURFACE_COLORS.background);
    if (loadId !== viewerLoadId) return;
    renderLegend();
    renderDesignSurfaceCallout();
    return;
  }
  if (state.layers.has("surface") && !state.layers.has("regions")) add(state.result.residues.filter((r) => r.surface_accessible), LAYER_META.surface.color);
  if (state.layers.has("folding") && !state.layers.has("regions")) add(state.result.residues.filter((r) => r.folding_screen_pass), LAYER_META.folding.color);
  if (state.layers.has("interface")) {
    add(state.result.residues.filter((r) => r.interface), LAYER_META.interface.color,
      { sideChain: true, representation: "ball-and-stick", representationColor: LAYER_META.interface.color });
  }
  if (state.layers.has("regions")) {
    for (const region of state.result.regions) {
      const positions = new Set(region.envelope_residues?.length ? region.envelope_residues : region.residues);
      const members = state.result.residues.filter((r) => region.partition_mode === "surface_v2"
        ? r.region_id === region.id : positions.has(r.uniprot_pos));
      add(members, state.selectedRegion && state.selectedRegion !== region.id ? "#C8CDD1" : region.color);
    }
  }
  const selectedFootprint = state.result.surface_v2?.footprints?.find(
    (item) => item.footprint_id === state.selectedFootprint
  );
  if (selectedFootprint) {
    const positions = new Set(selectedFootprint.residues || []);
    add(state.result.residues.filter((r) => positions.has(r.uniprot_pos)), "#D81B60");
  }
  if (state.layers.has("design")) {
    const positions = new Set(state.previewMode === "custom"
      ? state.customMutations.map((m) => m.position)
      : (selectedDesign()?.mutations || []).map((m) => m.position));
    add(state.result.residues.filter((r) => positions.has(r.uniprot_pos)), LAYER_META.design.color,
      { sideChain: true, representation: "ball-and-stick", representationColor: LAYER_META.design.color });
  }
  if (state.reliabilityView) {
    data.length = 0;
    state.result.residues.forEach((r) => add([r], reliabilityColor(r)));
  }
  await MolstarView.color(data, "#7c909d");
  if (loadId !== viewerLoadId) return;
  if (state.selectedResidue) await MolstarView.highlight(residueQuery(state.selectedResidue));
  renderLegend();
  if (state.engine === "beads") paintBeads();
}

function renderDesignSurfaceCallout() {
  const container = $("design-surface-callout");
  if (!container) return;
  const design = state.designSurfaceMode && state.previewMode === "construct" ? selectedDesign() : null;
  const region = regionForDesign(design);
  const regions = regionsForDesign(design);
  if (!design?.mutations?.length || !region) {
    container.hidden = true;
    container.innerHTML = "";
    return;
  }
  container.hidden = false;
  container.innerHTML = `<div class="design-surface-head">
    <div><strong>${h(t("design.surfaceTitle", "突变组合表面视图"))}</strong>
      <span>${regions.length > 1 ? h(recommendationLabel("跨区组合", "Cross-region combination")) + " · " : ""}${h(regions.map((item) => `${recommendationLabel("分区", "Region")} ${item.id}`).join("、"))}</span></div>
    <span class="grow"></span><button type="button" class="btn btn-mini" id="exit-design-surface">${h(t("design.surfaceExit", "退出聚焦"))}</button>
  </div>
  <div class="design-surface-labels">${design.mutations.map((mutation) =>
    `<button type="button" class="design-surface-label" data-surface-mutation="${mutation.position}">${h(notation(mutation))}</button>`
  ).join("")}</div>
  <div class="note" style="margin:0">${h(t("design.surfaceCoordinateNote", "表面来自输入结构；红色表面标记 WT 位点，标签表示设计替换，并非突变后原子构象。"))}</div>
  <div class="design-surface-key">
    <span><i style="background:${DESIGN_SURFACE_COLORS.background}"></i>${h(t("design.surfaceProtein", "其余表面"))}</span>
    <span><i style="background:${DESIGN_SURFACE_COLORS.region}"></i>${h(t("design.surfaceRegion", "区域"))}</span>
    <span><i style="background:${DESIGN_SURFACE_COLORS.footprint}"></i>${recommendationLabel("局部足迹", "Local footprint")}</span>
    <span><i style="background:${DESIGN_SURFACE_COLORS.mutation}"></i>${h(t("design.surfaceSites", "突变位点"))}</span>
  </div>`;
  $("exit-design-surface").addEventListener("click", exitDesignSurfaceMode);
  container.querySelectorAll("[data-surface-mutation]").forEach((button) => {
    button.addEventListener("click", () => selectResidue(Number(button.dataset.surfaceMutation), true));
  });
}

async function exitDesignSurfaceMode() {
  state.designSurfaceMode = false;
  state.designFocusRequest += 1;
  renderDesignSurfaceCallout();
  await paintMolstar();
  if (molstarReady) await MolstarView.focus(state.result.residues.map(residueQuery));
}

async function setRepresentation(representation) {
  state.representation = representation;
  document.querySelectorAll("[data-representation]").forEach((b) => {
    b.setAttribute("aria-pressed", String(b.dataset.representation === representation));
  });
  if (state.engine === "molstar" && molstarReady) {
    const loadId = viewerLoadId;
    try {
      await MolstarView.setRepresentation(representation);
      if (loadId === viewerLoadId) await paintMolstar();
    } catch (error) { if (loadId === viewerLoadId) viewerError(error, "render"); }
  }
}

function setEngine(engine) {
  if (!state.result) return;
  if (engine === "beads" && state.showPartner) {
    state.showPartner = false;
    MolstarView.setPartnerVisible(false).catch((error) => viewerError(error, "render"));
    renderStep4();
  }
  if (engine === "beads" && !state.jobId) {
    showToast(t("msg.beadsNeedJob", "冻结示例无作业结构文件；提交真实作业后可用珠子渲染"));
    document.querySelectorAll("[data-engine]").forEach((b) => {
      b.setAttribute("aria-pressed", String(b.dataset.engine === "molstar"));
    });
    return;
  }
  state.engine = engine;
  document.querySelectorAll("[data-engine]").forEach((b) => {
    b.setAttribute("aria-pressed", String(b.dataset.engine === engine));
  });
  const mol = $("viewer-molstar");
  const beads = $("viewer-beads");
  if (engine === "beads") {
    mol.style.visibility = "hidden";
    beads.hidden = false;
    if (!state.beadsReady) {
      const loadId = viewerLoadId;
      BeadsView.init(beads, state.offline ? "" : `/api/jobs/${state.jobId}/files/structure_for_energy.pdb`,
        state.result.residues, state.offline ? window.APP_OFFLINE.structure_pdb : null)
        .then((instance) => {
          if (loadId !== viewerLoadId || !instance) return;
          state.beadsReady = true;
          BeadsView.setCallbacks(
            (residue) => { if (residue) selectResidue(residue.uniprot_pos, false); },
            (residue) => showBeadsTip(residue)
          );
          paintBeads();
        })
        .catch((error) => {
          if (loadId !== viewerLoadId) return;
          console.error(error);
          showToast(t("msg.beadsFail", "珠子渲染载入失败：") + error.message);
          setEngine("molstar");
        });
    } else {
      paintBeads();
    }
  } else {
    beads.hidden = true;
    mol.style.visibility = "visible";
    BeadsView.destroy();
    state.beadsReady = false;
    hideBeadsTip();
  }
}

let beadsTipXY = { x: 0, y: 0 };
function showBeadsTip(residue) {
  const tip = $("residue-tip");
  if (!residue) { tip.hidden = true; return; }
  const ddg = foldxDdgFor(residue.uniprot_pos);
  tip.textContent = `${residue.aa}${residue.uniprot_pos} · rSASA ${fmt(residue.rsasa, 2)}`
    + (ddg != null ? ` · FoldX ΔΔG ${fmt(ddg, 3)}` : "");
  tip.style.left = `${beadsTipXY.x}px`;
  tip.style.top = `${beadsTipXY.y}px`;
  tip.hidden = false;
}
function hideBeadsTip() { $("residue-tip").hidden = true; }

function renderLegend() {
  if (!state.result) { $("legend").innerHTML = ""; return; }
  const items = [];
  if (state.designSurfaceMode && selectedDesign()?.mutations?.length) {
    items.push(legendItem(DESIGN_SURFACE_COLORS.background, t("design.surfaceProtein", "其余表面")));
    items.push(legendItem(DESIGN_SURFACE_COLORS.region, t("design.surfaceRegion", "所选区域")));
    items.push(legendItem(DESIGN_SURFACE_COLORS.footprint, "Footprint"));
    items.push(legendItem(DESIGN_SURFACE_COLORS.mutation, t("design.surfaceSites", "突变位点")));
  } else {
  const design = selectedDesign();
  const designCount = (design?.mutations || []).length + state.customMutations.length;
  if (designCount) items.push(legendItem(LAYER_META.design.color, `${t("legend.selected", "已选突变")} (${designCount})`));
  if (state.layers.has("interface")) items.push(legendItem(LAYER_META.interface.color, interfaceLabel()));
  if (state.layers.has("regions")) {
    state.result.regions.forEach((region) => items.push(legendItem(region.color, `${recommendationLabel("分区", "Region")} ${region.id}`, region.name)));
  }
  items.push(legendItem("#7c909d", t("hud.backbone", "主链")));
  }
  if (state.reliabilityView) {
    items.length = 0;
    Object.entries(RELIABILITY_COLORS).forEach(([name, color]) => items.push(legendItem(color, name)));
  }
  const collapsed = (() => { try { return localStorage.getItem("ecd.legend") === "collapsed"; } catch { return false; } })();
  $("legend").innerHTML =
    `<button type="button" class="legend-toggle">${collapsed ? t("hud.legendCollapsed", "图例 ▸") : t("hud.legendOpen", "图例 ▾")}</button>`
    + (collapsed ? "" : items.join(""));
  const toggle = $("legend").querySelector(".legend-toggle");
  toggle.addEventListener("click", () => {
    try { localStorage.setItem("ecd.legend", collapsed ? "open" : "collapsed"); } catch { /* optional preference */ }
    renderLegend();
  });
}

const legendItem = (color, text, title = text) =>
  `<span class="legend-item" title="${h(title)}"><i class="legend-swatch" style="background:${h(color)}"></i>${h(text)}</span>`;

function interfacePartnerChains() {
  if (!state.result) return [];
  const chains = new Set();
  for (const residue of state.result.residues) {
    if (!residue.interface) continue;
    for (const partner of residue.interface_partners || []) chains.add(String(partner).split(":")[0]);
  }
  return [...chains].sort();
}

function interfaceLabel() {
  const chains = interfacePartnerChains();
  if (!chains.length) return t("interface.none", "无蛋白伙伴（apo/单体）");
  return `${t("interface.with", "界面")} ${chains.join("/")}`;
}

/* ---------------- steps ---------------- */

function renderSteps() {
  $("steps").innerHTML = STEP_LIST.map((key, i) => {
    const n = i + 1;
    const cls = `step-pill${state.step === n ? " active" : ""}${n < state.step ? " done" : ""}`;
    return `<button type="button" class="${cls}" data-step="${n}"${state.step === n ? ' aria-current="step"' : ""}><span class="n">${n}</span>${h(t(key, key.split(".")[1]))}</button>`;
  }).join("") + `<div class="step-controls"><button type="button" class="btn btn-quiet" id="previous-step" ${state.step === 1 ? "disabled" : ""}>${h(t("workspace.previous", "上一步"))}</button><button type="button" class="btn" id="next-step" ${state.step === STEP_LIST.length ? "disabled" : ""}>${h(t("workspace.next", "下一步"))} →</button></div>`;
  document.querySelectorAll("#steps [data-step]").forEach((b) => {
    b.addEventListener("click", () => setStep(Number(b.dataset.step)));
  });
  $("previous-step").addEventListener("click", () => setStep(state.step - 1));
  $("next-step").addEventListener("click", () => setStep(state.step + 1));
}

function setStep(n) {
  if (!Number.isInteger(n) || n < 1 || n > STEP_LIST.length) return;
  const changed = state.step !== n;
  state.step = n;
  if (n !== 4 && state.showPartner) {
    state.showPartner = false;
    MolstarView.setPartnerVisible(false).catch((error) => viewerError(error, "render"));
    renderStep4();
  }
  document.body.dataset.workspaceStep = String(n);
  renderSteps();
  document.querySelectorAll("#panel [data-step-content]").forEach((s) => {
    s.hidden = Number(s.dataset.stepContent) !== n;
  });
  if (changed) $("panel").scrollTop = 0;
}

function renderPanelScaffold() {
  $("panel").innerHTML = `
  <section data-step-content="1">
    <h2>${h(t("step.surface", "表面可及性"))}</h2>
    <div class="note" id="surface-gate-note">${h(surfaceGateDescription(state.result, state.lang))}</div>
    <div id="surface-stats" class="statgrid" style="margin-bottom:10px"></div>
    <div id="source-card"></div>
    <div class="filter-row">
      <label class="fld">${h(t("s1.search", "搜索"))} <input id="residue-search" type="search" placeholder="${h(t("s1.searchPlaceholder", "如 123、K123 或 interface"))}"></label>
      <label class="fld">${h(t("s1.filter", "筛选"))}
        <select id="residue-filter">
          <option value="folding">${h(t("filter.folding", "折叠风险通过"))}</option>
          <option value="surface">${h(t("filter.surface", "表面可及"))}</option>
          <option value="interface">${h(t("filter.interface", "功能界面"))}</option>
          <option value="protected">${h(t("filter.protected", "受保护/排除"))}</option>
          <option value="all">${h(t("filter.all", "全部 ECD"))}</option>
        </select></label>
    </div>
    <div id="residue-list"></div>
    <div style="margin-top:12px"><h2>${h(t("s1.detailTitle", "位点详情"))}</h2></div>
    <div id="residue-detail" class="empty">${h(t("s1.detailEmpty", "在结构、序列或列表中选择残基。"))}</div>
  </section>

  <section data-step-content="2" hidden>
    <h2>${h(t("step.tolerance", "折叠耐受"))} · ${h(t("fold.realScores", "真实能量评分"))}</h2>
    <div class="note">${h(t("fold.note", "本模块只显示 FoldX / EvoEF2 实际计算得到的能量分数。引擎未产出评分时不显示任何数字，只显示状态。"))}</div>
    <div id="foldx-status"></div>
    <div id="foldx-diagnostics"></div>
    <div id="foldx-table"></div>
    <div class="card" style="margin-top:12px">
      <h3>${h(t("s2.methods", "方法"))}</h3>
      <div id="methods-card" class="note" style="margin:0"></div>
    </div>
    <div class="card">
      <h3>${h(t("s2.warnings", "运行警告"))}</h3>
      <ul id="warning-list" class="note" style="margin:0;padding-left:18px"><li>${h(t("s2.awaiting", "等待分析结果。"))}</li></ul>
    </div>
  </section>

  <section data-step-content="3" hidden>
    <h2>${h(t("step.regions", "结构分区"))}</h2>
    <div class="note">${h(t("region.note", "Fab 足迹尺度补丁：测地半径 13 Å 且欧氏距离 ≤18 Å，允许跨结构域；点击分区头可聚焦该区。"))}</div>
    <div id="region-list"></div>
  </section>

  <section data-step-content="4" hidden>
    <h2>${h(t("step.interface", "蛋白-蛋白界面"))}（${h(t("s4.measured", "本结构实测"))}）</h2>
    <div class="note">${h(t("iface.note", "界面残基来自当前结构中实测的 5 Å 链伙伴接触——它是这条结构里实际存在的结合伙伴，不是预测。若选中结构为单体/apo，则没有可显示的界面。"))}</div>
    <div id="iface-partners"></div>
    <div class="card" style="margin-top:12px">
      <h3>${h(t("s4.layers", "结构叠加图层"))}</h3>
      <div id="layer-controls"></div>
      <div id="criterion" class="note" style="margin:8px 0 0"></div>
    </div>
  </section>

  <section data-step-content="5" hidden>
    <h2>${h(t("step.design", "推荐组合"))}</h2>
    <div class="note">${h(t("design.note", "每个补丁一条 Stage-1 探针（3-6 点，带锚点）；列表同时含备选探针与 Stage-2 单点/双点预设（命中后才需订购）。HIGH 补丁有两条探针覆盖。"))}</div>
    <div id="design-list"></div>
    <div class="card" style="margin-top:12px">
      <h3>${h(t("s5.exports", "导出文件"))}</h3>
      <div id="export-links" class="export-links empty">${h(t("s5.exportsEmpty", "完成分析后提供 JSON、CSV、FASTA、GenBank、方法报告和 ZIP。"))}</div>
    </div>
  </section>

  <section data-step-content="6" hidden>
    <h2>${h(t("step.custom", "自主设计"))}</h2>
    <div class="note">${h(t("custom.note", "在位点详情（第 1 步）选择残基与目标氨基酸后加入组合；点击 3D 或底部序列字母可定位残基。"))}</div>
    <div class="card">
      <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px">
        <h3 style="margin:0">${h(t("custom.current", "我的组合"))}</h3>
        <span><button class="btn btn-mini" id="preview-custom" type="button">预览自定义草稿</button>
        <button class="btn btn-mini" id="clear-custom" type="button">${h(t("custom.clear", "清空"))}</button></span>
      </div>
      <label class="fld" style="margin-bottom:8px"><span>${h(t("custom.name", "构建名称"))}</span>
        <input id="custom-name" value="Custom construct"></label>
      <div id="custom-mutations" class="empty">${h(t("custom.empty", "还没有选残基。"))}</div>
      <label class="chk"><input id="force-risk" type="checkbox">
        <span>${h(t("custom.forceRisk", "显式保留高风险替换（受保护位点仍不可强制）"))}</span></label>
      <button id="save-custom" class="btn primary" type="button" style="width:100%" disabled>${h(t("custom.save", "确认并保存此组合"))}</button>
      <p id="custom-feedback" class="note" style="margin:8px 0 0"></p>
    </div>
  </section>

  <section data-step-content="7" hidden>
    <h2>${h(t("step.overview", "任务详情"))}</h2>
    <div id="task-overview"></div>
  </section>`;
  bindPanelControls();
  setStep(state.step);
}

function bindPanelControls() {
  $("residue-filter").addEventListener("change", renderResidueList);
  $("residue-search").addEventListener("input", renderResidueList);
  $("design-select").addEventListener("change", () => selectDesign($("design-select").value));
  $("clear-custom").addEventListener("click", () => {
    state.customMutations = [];
    state.previewMode = "construct";
    renderCustomBuilder();
    renderSeqDock();
    repaintViews();
  });
  $("preview-custom").addEventListener("click", () => {
    if (!state.customMutations.length) return;
    state.previewMode = "custom";
    state.designSurfaceMode = false;
    state.designFocusRequest += 1;
    renderDesignSurfaceCallout();
    renderSeqDock();
    repaintViews();
  });
  $("save-custom").addEventListener("click", saveCustomDesign);
}

/* ---------------- step renderers ---------------- */

function renderAll(repaint = true) {
  if (!state.result) return;
  renderSteps();
  renderStep1();
  renderStep2();
  renderStep3();
  renderStep4();
  renderStep5();
  renderStep6();
  renderStep7();
  applyTaskIdentity();
  renderLegend();
  renderSeqDock();
  if (repaint) {
    repaintViews();
  }
}

function renderStep1() {
  $("surface-gate-note").textContent = surfaceGateDescription(state.result, state.lang);
  const s = state.result.summary;
  const values = [
    [s.mapped_ecd_residues, t("s1.stat.mapped", "ECD 已映射")],
    [s.surface_accessible, t("s1.stat.surface", "表面可及")],
    [s.folding_screen_pass, t("s1.stat.folding", "折叠通过")],
    [s.functional_interface, t("s1.stat.interface", "功能界面")],
    [s.regions, t("s1.stat.regions", "补丁")],
    [s.recommended_constructs, t("s1.stat.designs", "推荐组合")],
  ];
  $("surface-stats").innerHTML = values.map(([v, l]) =>
    `<div class="stat"><strong>${h(v)}</strong><span>${h(l)}</span></div>`).join("");

  const source = state.result.source;
  const sel = state.result.structure_selection;
  const rows = [
    [t("s1.src.protein", "蛋白"), `${state.result.protein.protein_name} (${state.result.protein.accession})`],
    [t("s1.src.ecd", "ECD"), `${state.result.protein.ecd_start}-${state.result.protein.ecd_end} · ${state.result.protein.ecd_length} aa`],
    [t("s1.src.structure", "结构"), `${source.identifier || "user upload"} · ${source.provider || source.kind}`],
    [t("s1.src.chain", "链与覆盖"), `${sel.selected_chain} · ${fmtPct(sel.observed_ecd_coverage)}`],
    [t("s1.src.identity", "序列一致性"), fmtPct(sel.aligned_sequence_identity)],
  ];
  $("source-card").innerHTML = `<div class="card"><h3>${h(t("s1.src.title", "结构来源"))}</h3><dl class="kv">${
    rows.map(([k, v]) => `<dt>${h(k)}</dt><dd>${h(v)}</dd>`).join("")
  }${source.page_url ? `<dt>${h(t("s1.src.page", "来源页"))}</dt><dd><a href="${h(source.page_url)}" target="_blank" rel="noreferrer">${h(t("s1.src.open", "打开权威来源"))}</a></dd>` : ""}</dl></div>`;

  renderResidueList();
  renderResidueDetail(true);
}

function renderStep2() {
  const design = selectedDesign();
  const statusEl = $("foldx-status");
  const diagEl = $("foldx-diagnostics");
  const tableEl = $("foldx-table");

  const stability = state.result.stability || {};
  const diag = state.result.design_summary?.foldx_diagnostics || null;

  /* --- engine availability line (facts, not scores) --- */
  diagEl.innerHTML = `<div class="card"><h3>${h(t("s2.engineStatus", "本次运行引擎状态"))}</h3><dl class="kv">
    <dt>FoldX</dt><dd>${stability.foldx_available ? t("common.yes", "是") : t("common.no", "否")}</dd>
    <dt>EvoEF2</dt><dd>${stability.evoef2_available ? t("common.yes", "是") : t("common.no", "否")}</dd>
    ${diag ? `<dt>${h(t("s2.batch", "组合批量"))}</dt><dd>${diag.batch_scores_returned}/${diag.batch_constructs} ${h(t("s2.success", "成功"))}` +
      `（${h(t("s2.gates", "门限"))} ${fmt(diag.single_gate_kcal_mol?.default ?? diag.single_gate_kcal_mol)} / ${fmt(diag.combination_gate_kcal_mol)} kcal/mol）</dd>` : ""}
  </dl></div>`;

  /* --- combination score: real number or explicit status --- */
  if (!design) {
    statusEl.innerHTML = `<div class="status-card"><div class="big">${h(t("fold.noDesign", "未选择构建体"))}</div>
      <div class="note" style="margin:0">${h(t("fold.pickDesign", "到第 5 步选择一个推荐组合后，这里显示它的真实 FoldX 评分。"))}</div></div>`;
    tableEl.innerHTML = "";
    renderMethods();
    return;
  }

  const info = scoreInfo(design);
  if (info.real) {
    statusEl.innerHTML = `<div class="card">
      <div class="muted" style="font-size:11.5px">${h(design.id)} · ${h(design.name)}</div>
      <div><span class="score-big score-real">${fmt(info.value, 3)}</span>
        <span class="score-unit">kcal/mol · ${h(t("s2.comboMedian", "组合 ΔΔG（中位数）"))}</span></div>
      <div class="note" style="margin:6px 0 0">${h(info.engine)} · score_origin=${h(info.origin)}</div>
    </div>`;
  } else {
    const reasons = (design.warnings || []).filter((w) => /foldx|proxy|score/i.test(w));
    statusEl.innerHTML = `<div class="status-card">
      <div class="big">${h(t("fold.pendingTitle", "尚未获得真实能量评分"))}</div>
      <div class="note" style="margin:0 0 6px">${h(t("fold.pendingBody",
        "该构建体的组合与单点能量尚未由 FoldX/EvoEF2 计算产出，因此这里不显示任何分数——显示数字反而是误导。"))}</div>
      <dl class="kv">
        <dt>score_origin</dt><dd>${h(info.origin)}</dd>
        <dt>${h(t("s2.engineRecord", "稳定性引擎记录"))}</dt><dd>${h(info.engine)}</dd>
        ${reasons.length ? `<dt>${h(t("s2.reasons", "原因"))}</dt><dd>${reasons.map((w) => h(w)).join("；")}</dd>` : ""}
      </dl>
      <div class="note" style="margin:6px 0 0">${h(t("fold.howto",
        "提交真实作业并确保 FoldX 引擎可用（vendor/foldx 存在或设置 FOLDX_PATH）后，此处会显示 kcal/mol 真值。"))}</div>
    </div>`;
  }

  /* --- single-point table: real numbers only, otherwise "未评分" --- */
  const rows = [];
  Object.keys(design.saturation_scan || {}).forEach((pos) => {
    (design.saturation_scan[pos] || []).forEach((choice) => {
      if (!choice.allowed) return;
      const sc = singleChoiceScore(choice);
      rows.push(`<tr>
        <td><b>${h(choice.wt)}${choice.position}</b></td>
        <td>${h(choice.mutant)}</td>
        <td>${sc.real ? fmt(sc.value, 3) : `<span class="muted">${h(t("s2.unscored", "未评分"))}</span>`}</td>
        <td>${h(choice.stability_class && sc.real ? choice.stability_class : "")}</td>
      </tr>`);
    });
  });
  tableEl.innerHTML = rows.length
    ? `<div class="card"><h3>${h(t("s2.passedTable", "通过单点门的替换（ΔΔG 列仅显示 FoldX 真值）"))}</h3>
       <table class="audit"><tbody><tr><th>${h(t("tol.site", "位点"))}</th><th>${h(t("s2.substitution", "替换"))}</th><th>ΔΔG</th><th>${h(t("tol.class", "分类"))}</th></tr>${rows.join("")}</tbody></table></div>`
    : `<div class="card"><div class="empty">${h(t("s2.noPassed", "该构建体没有通过单点门的替换记录。"))}</div></div>`;

  renderMethods();
}

function renderMethods() {
  const stability = state.result.stability || {};
  const rows = [
    `<strong>SASA</strong>　${h(state.result.methods.sasa_method)}`,
    `<strong>${h(t("s2.partition", "分区"))}</strong>　${h(state.result.partitioning.method)}`,
    `<strong>${h(t("s2.engineLabel", "稳定性引擎"))}</strong>　${h(stability.primary || t("s2.notRecorded", "未记录"))}；${h(stability.proxy_units || "")}`,
    `<strong>${h(t("s2.boundary", "结论边界"))}</strong>　${h(t("s2.boundaryBody", "计算预筛不等于实验确认；表位归因仍需表达、构象抗体、稳定性与功能读出。"))}`,
  ];
  const af3 = state.result.methods?.af3_preprocessing;
  if (af3) {
    const baselines = Object.entries(af3.baselines || {}).map(([key, value]) =>
      `${h(key)}: ${h(Array.isArray(value) ? value.join(", ") : typeof value === "object" ? value?.model_id ?? "—" : value ?? "—")}`
    ).join(" · ");
    const models = (af3.models || []).map((model) =>
      `${h(model.model_id)}${model.model_id === af3.selected_model_id ? ` (${h(t("af3.selected", "selected"))})` : ""}: ${h(model.eligible ? t("af3.eligible", "eligible") : (model.exclusion_reasons || []).join(", ") || t("af3.ineligible", "ineligible"))}`
    );
    rows.push(`<strong>${h(t("af3.preprocessing", "AF3 model preprocessing"))}</strong>　${h(af3.method || af3.status || "—")} · ${h(t("af3.selectedModel", "Selected model"))}: ${h(af3.selected_model_id || "—")}${baselines ? `<br><span>${h(t("af3.baselines", "Baselines"))}: ${baselines}</span>` : ""}${models.length ? `<br><span>${h(t("af3.models", "Models"))}: ${models.join(" · ")}</span>` : ""}`);
    const columns = ["model_id", "coverage", "f70", "mean_plddt", "ptm", "ranking_score", "consensus_score"];
    rows.push(`<div class="table-wrap"><table class="audit"><thead><tr>${columns.map(key => `<th>${h(key)}</th>`).join("")}</tr></thead><tbody>${(af3.models || []).map(model => `<tr>${columns.map(key => `<td>${h(typeof model[key] === "number" ? model[key].toFixed(4) : model[key] ?? "—")}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`);
    rows.push(h([...(af3.warnings || []), af3.core_max_rmsd_flag ? "Reliable core pairwise RMSD > 2 Å" : ""].filter(Boolean).join(" · ")));
  }
  $("methods-card").innerHTML = rows.join("<br>");
  const warnings = state.result.warnings?.length ? state.result.warnings : [t("s2.noWarnings", "没有运行级警告。")];
  $("warning-list").innerHTML = warnings.map((w) => `<li>${h(w)}</li>`).join("");
  $("criterion").textContent = criteriaText();
}

function criteriaText() {
  if (!state.result) return t("s2.criteriaEmpty", "载入结果后显示本步骤的定量标准。");
  const m = state.result.methods;
  return {
    surface: `${surfaceGateDescription(state.result, state.lang)}${m.sasa_method ? `; SASA: ${m.sasa_method}` : ""}`,
    folding: t("crit.folding", "在表面标准上排除二硫键 Cys、注释糖基化/PTM 与低置信度位点；表示低预测风险，不是实验确认。"),
    interface: t("crit.interface", "任一目标重原子距非目标链或异质残基 ≤5 Å；功能界面保留风险标记。"),
    regions: t("crit.regions", `${state.result.partitioning.method}；实际 ${state.result.regions.length} 补丁，覆盖 ${fmtPct(state.result.partitioning.coverage)}。`),
    design: t("crit.design", `FoldX 单点门与暴露度联动（+0.46/+0.60/+0.80）；组合门 +2.0 kcal/mol；探针内两两间距 ≥5.5 Å 且最大展开 ≤${state.result.partitioning?.partition_mode === "compact_surface_face" ? 15 : 22} Å。`),
  }[state.layerFocus] || t("crit.default", "当前图层可与其他图层叠加显示。");
}

const { sectorRecommendationText, renderSectorRecommendations, recommendationLabel, surfaceGateDescription, recommendationSelectionRecord, recommendationRiskLabelText, recommendationNeedsReview, recommendationRepeatSummary, selectionColor } = createRecommendationRenderer({ state, h, fmt, fmtPct, notation: (mutation) => notation(mutation) });


function renderStep3() {
  const policy = state.result.partitioning?.display_partition_policy;
  if (state.result.joint_v4_ml) {
    const note = $("region-list").previousElementSibling;
    if (note?.classList.contains("note")) note.textContent = state.lang === "zh"
      ? "颜色表示空间连续、直径受限的显示区域；区域数可因几何约束增加。重叠足迹是假设层，不是已验证表位。AF3 位置不可靠的末端保留在结构中，但不参与三维分区；序列相距远不等于空间相距远。"
      : "Colors represent connected display regions with bounded diameters; geometry may require extra regions. Overlapping footprints are hypotheses, not validated epitopes. Unreliable AF3 termini remain in the structure and are omitted from 3D partitioning; sequence distance is not spatial distance.";
  }
  if (state.result.partitioning?.partition_mode === "surface_v2") {
    const note = $("region-list").previousElementSibling;
    if (note?.classList.contains("note")) note.textContent = state.lang === "zh"
      ? "表面区域由主尺度足迹聚合而成，基础颜色按残基唯一归属显示。点击区域可突出显示，再次点击恢复全部区域；下方可查看重叠足迹。"
      : "Regions group primary footprints. Each residue has one base region color. Click a region to highlight it; click again to restore all regions. Overlapping footprints remain available below.";
  }
  if (policy) {
    const note = $("region-list").previousElementSibling;
    const target = policy.target_sector_residues;
    const surface = policy.maximum_geodesic_diameter_a;
    const ca = policy.maximum_euclidean_diameter_a;
    const surfaceLabel = recommendationLabel("全表面图测地直径", "Full-surface graph geodesic diameter");
    if (note?.classList.contains("note")) note.textContent += state.lang === "zh"
      ? ` 当前任务：${state.result.regions.length}区；${target != null ? `目标每区约${target}个合格表面残基（软目标）；` : ""}每区至少${policy.minimum_sector_residues}残基、${policy.minimum_sector_area_a2} Å²；${surfaceLabel}≤${surface} Å，Cα直径≤${ca} Å。显示分区数量不是构建体预算。`
      : ` This run: ${state.result.regions.length} regions; ${target != null ? `about ${target} eligible surface residues per region (soft target); ` : ""}at least ${policy.minimum_sector_residues} residues and ${policy.minimum_sector_area_a2} Å² per region; ${surfaceLabel} ≤${surface} Å, Cα diameter ≤${ca} Å. Region count is independent of the construct budget.`;
    if (note?.classList.contains("note") && policy.preferred_minimum_sector_count != null) note.textContent += state.lang === "zh"
      ? ` 总分区数以${policy.preferred_minimum_sector_count}–${policy.preferred_maximum_sector_count}个为软范围。`
      : ` Preferred region count: ${policy.preferred_minimum_sector_count}–${policy.preferred_maximum_sector_count} (soft range).`;
  }
  const recommendations = state.result.joint_v4_ml?.selected_panel?.sector_recommendations || [];
  const designByRegion = new Map(state.result.designs.filter((d) => !d.design_role || d.design_role === "STAGE1_PROBE").map((d) => [d.region_id, d]));
  const legacyHtml = state.result.regions.map((region) => {
    const sector = state.result.joint_v4_ml?.sectors?.[state.result.regions.indexOf(region)];
    const recommendation = recommendations.find((item) => item.sector_id === sector?.sector_id);
    const design = recommendation ? state.result.designs.find((item) => item.id === recommendation.design_id) : designByRegion.get(region.id);
    const mutations = design?.mutations?.map(notation).join(" · ") || t("region.noAuto", "无自动组合");
    const tier = region.tier ? ` · ${region.tier}` : "";
    const accessibility = region.predicted_accessibility === "reduced"
      ? ` · ${t("region.reducedAccessibility", "预测细胞表面可及性降低")} (${fmt(region.membrane_accessibility_score)})`
      : "";
    const metric = policy?.region_metrics?.find((row) => row.sector_id === sector?.sector_id);
    const exception = policy?.size_exceptions?.find((row) => row.sector_id === sector?.sector_id);
    const geometry = metric ? (state.lang === "zh"
      ? `全表面图测地直径 ${fmt(metric.full_graph_geodesic_diameter_a)} Å · 区内路径直径 ${fmt(metric.induced_graph_geodesic_diameter_a)} Å（诊断） · Cα直径 ${fmt(metric.ca_diameter_a)} Å`
      : `Full-surface graph geodesic diameter ${fmt(metric.full_graph_geodesic_diameter_a)} Å · Within-region path diameter ${fmt(metric.induced_graph_geodesic_diameter_a)} Å (diagnostic) · Cα diameter ${fmt(metric.ca_diameter_a)} Å`) : "";
    const exceptions = { COMPONENT_SIZE_BELOW_PREFERRED_RANGE: recommendationLabel("该连通表面块本身较小", "The connected surface component is small"), LOCAL_BALANCE_LIMIT: recommendationLabel("在当前连通性与直径约束下，局部调整仍未达到目标大小", "Local adjustment could not reach the target size under the connectivity and diameter limits"), SEARCH_BUDGET_EXHAUSTED: recommendationLabel("调整搜索预算已用尽，未证明无法得到更均衡的分区", "Adjustment search budget exhausted; a more balanced partition has not been ruled out") };
    const sizeNote = exception ? `${recommendationLabel("大小例外", "Size exception")}：${exceptions[exception.reason] || recommendationLabel("原因待核对", "Reason pending review")}` : "";
    return `<div role="button" tabindex="0" class="rowline click${state.selectedRegion === region.id ? " on" : ""}" data-focus-region="${region.id}" aria-pressed="${state.selectedRegion === region.id}">
      <span class="swatch" style="background:${h(region.color)}"></span>
      <span><b>${h(region.partition_mode === "compact_surface_face" ? t("region.surfaceFace", "表面区域") : t("region.patch", "补丁"))} ${region.id} · ${h(region.name)}</b><br>
      <span class="muted">${h(region.structural_character)} · ${h(region.membrane_position || "")}${h(tier)}${h(accessibility)} · ${h(mutations)}${recommendation ? `<br>${h(sectorRecommendationText(recommendation))}` : ""}</span></span>
      <span class="grow"></span><span class="muted">${(region.envelope_residues?.length ?? region.residues.length)} ${h(t("region.surfaceSites", "表面残基"))}<br>${(region.mutable_residues?.length ?? region.residues.length)} ${h(t("region.mutableSites", "可突变"))}</span>
    </div><div class="note" style="margin:0 0 10px 20px">${h(region.rationale)}${geometry ? ` · ${h(geometry)}` : ""}${sizeNote ? `<br>${h(sizeNote)}` : ""}${region.footprint_residues?.length ? ` · ${h(t("region.footprint", "局部足迹"))} ${region.footprint_residues.length} ${h(t("region.sites", "位点"))}` : ""}${region.candidate_footprints?.length ? ` · ${region.candidate_footprints.length} ${h(t("region.candidateFootprints", "个候选足迹"))}` : ""}${region.partition_basis ? ` · ${h(t("region.basis", "依据"))}: ${h(region.partition_basis)}` : ""}</div>`;
  }).join("") || `<div class="empty">${h(t("region.none", "没有可辩护的结构分区。"))}</div>`;
  const v2 = state.result.surface_v2;
  const footprintHtml = v2?.footprints?.length ? `<div class="box" style="margin-top:12px">
    <h3>Surface V2 footprints</h3>
    <div class="note">${v2.qc?.deduplicated_footprint_count ?? v2.footprints.length} overlapping footprints · ${v2.qc?.surface_region_count ?? 0} surface regions · coverage ${fmtPct(v2.qc?.coverage ?? 0)}</div>
    <label class="field"><span>Footprint</span><select id="footprint-select">
      <option value="">—</option>${v2.footprints.map((item) => `<option value="${h(item.footprint_id)}">${h(item.footprint_id)} · ${h(item.generation_mode)} · ${item.residues?.length || 0} residues · ${fmt(item.total_surface_area)} Å²</option>`).join("")}
    </select></label>
  </div>` : "";
  const reliability = state.result.structure_reliability;
  const reliabilityHtml = reliability ? `<details class="box"><summary>结构可靠性 / Reliability · ${reliability.idr_candidates.length} sequence candidates</summary>
    <label><input id="reliability-view" type="checkbox" ${state.reliabilityView ? "checked" : ""}> 可靠性着色（覆盖分区颜色）</label>
    <div class="note">多证据规则分流，不是无序概率；IDR 窗口仅供审核，不自动运行 FoldX。${reliability.qc.enabled ? "" : "分流已禁用。"}</div>
    <div class="note">PAE: ${h(reliability.qc.pae_status)} · IUPred: ${h(reliability.qc.sequence_predictor?.status || "unavailable")} · boundary guard: ${reliability.qc.boundary_positions.length} residues</div>
    ${Object.entries(RELIABILITY_COLORS).map(([name, color]) => `<span style="color:${color}">■</span> ${name} (${reliability.qc.state_counts[name] || 0})`).join(" · ")}
    ${reliability.segments.map((s) => `<div>${s.start}–${s.end}: ${h(s.state)} · ${s.use_for_surface_mesh ? "保留几何" : "不进入三维分区"}</div>`).join("")}
    <details><summary>IDR / conditional sequence windows</summary>${reliability.idr_candidates.map((c) => `<div>${c.start}–${c.end} ${h(c.sequence)} · ${h(c.candidate_type)} · REVIEW_REQUIRED</div>`).join("") || "无候选"}</details>
  </details>` : "";
  $("region-list").innerHTML = legacyHtml + footprintHtml + reliabilityHtml;
  $("reliability-view")?.addEventListener("change", (event) => { state.reliabilityView = event.target.checked; repaintViews(); });
  document.querySelectorAll("[data-focus-region]").forEach((el) => {
    el.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") { event.preventDefault(); el.click(); }
    });
    el.addEventListener("click", async () => {
      const id = Number(el.dataset.focusRegion);
      state.selectedRegion = state.selectedRegion === id ? null : id;
      state.selectedFootprint = null;
      state.designSurfaceMode = false;
      state.reliabilityView = false;
      state.layers.add("regions");
      document.querySelectorAll("[data-focus-region]").forEach((button) => {
        const selected = Number(button.dataset.focusRegion) === state.selectedRegion;
        button.setAttribute("aria-pressed", String(selected));
        button.classList.toggle("on", selected);
      });
      renderDesignSurfaceCallout();
      repaintViews();
      if (molstarReady) {
        const positions = state.selectedRegion ? new Set(state.result.regions.find((region) => region.id === id).residues) : null;
        try { await MolstarView.focus(state.result.residues.filter((residue) => !positions || positions.has(residue.uniprot_pos)).map(residueQuery)); }
        catch (error) { viewerError(error, "focus"); }
      }
    });
  });
  const footprintSelect = $("footprint-select");
  if (footprintSelect) {
    footprintSelect.value = state.selectedFootprint || "";
    footprintSelect.addEventListener("change", async () => {
      state.selectedFootprint = footprintSelect.value || null;
      state.selectedRegion = null;
      repaintViews();
      const footprint = v2.footprints.find((item) => item.footprint_id === state.selectedFootprint);
      if (footprint && molstarReady) {
        const positions = new Set(footprint.residues || []);
        await MolstarView.focus(state.result.residues.filter((r) => positions.has(r.uniprot_pos)).map(residueQuery));
      }
    });
  }
}

function renderStep4() {
  const chains = interfacePartnerChains();
  const contacts = state.result.residues.filter((r) => r.interface);
  const partnerToggle = chains.length
    ? `<button type="button" class="btn" id="toggle-protein-partner" aria-pressed="${state.showPartner}">${h(state.showPartner ? t("iface.hidePartner", "隐藏伙伴链") : t("iface.showPartner", "显示伙伴链"))}</button>` : "";
  $("iface-partners").innerHTML = chains.length
    ? chains.map((chain) => {
      const count = new Set(contacts.flatMap((r) => (r.interface_partners || [])
        .filter((p) => p.split(":")[0] === chain))).size;
      return `<div class="rowline"><span class="swatch" style="background:${LAYER_META.interface.color}"></span>
        <span><b>${h(t("iface.with", "与"))} ${h(chain)}</b><br><span class="muted">${count} ${h(t("iface.contacts", "个接触残基"))}</span></span></div>`;
    }).join("") + partnerToggle
    : `<div class="empty">${h(t("iface.none", "当前结构无蛋白伙伴，界面叠加不可用。"))}</div>`;
  $("toggle-protein-partner")?.addEventListener("click", async () => {
    if (!molstarReady || state.step !== 4) return;
    if (state.engine === "beads") setEngine("molstar");
    state.showPartner = !state.showPartner;
    try { await MolstarView.setPartnerVisible(state.showPartner); renderStep4(); }
    catch (error) { state.showPartner = state.step === 4 && !state.showPartner; renderStep4(); viewerError(error, "render"); }
  });

  const entries = [
    ["surface", t("layer.surface", "充分表面可及"), t("layer.surfaceDetail", "rSASA 与侧链 SASA 双阈值"), LAYER_META.surface.color],
    ["folding", t("layer.folding", "折叠风险通过"), t("layer.foldingDetail", "保护位点排除后的候选"), LAYER_META.folding.color],
    ["interface", t("layer.interface", "功能/结构界面"), t("layer.interfaceDetail", "保留功能风险标记"), LAYER_META.interface.color],
    ["regions", t("layer.regions", "Fab 足迹补丁"), t("layer.regionsDetail", "按补丁颜色显示"), LAYER_META.regions.color],
    ["design", t("layer.design", "当前突变组合"), t("layer.designDetail", "叠加显示并标注侧链"), LAYER_META.design.color],
  ];
  $("layer-controls").innerHTML = entries.map(([key, label, detail, color]) => `
    <div class="rowline">
      <label style="display:flex;gap:8px;align-items:center;flex:1;cursor:pointer">
        <input type="checkbox" data-layer="${key}" ${state.layers.has(key) ? "checked" : ""}>
        <span class="swatch" style="background:${color}"></span>
        <span><b>${h(label)}</b><br><span class="muted" style="font-size:11.5px">${h(detail)}</span></span>
      </label>
      <button type="button" class="btn btn-mini" data-layer-focus="${key}">${h(t("layer.locate", "定位"))}</button>
    </div>`).join("");
  document.querySelectorAll("[data-layer]").forEach((input) => input.addEventListener("change", () => {
    if (input.checked) state.layers.add(input.dataset.layer);
    else state.layers.delete(input.dataset.layer);
    state.layerFocus = input.dataset.layer;
    $("criterion").textContent = criteriaText();
    repaintViews();
  }));
  document.querySelectorAll("[data-layer-focus]").forEach((button) => button.addEventListener("click", () => {
    state.layers.add(button.dataset.layerFocus);
    state.layerFocus = button.dataset.layerFocus;
    renderStep4();
    repaintViews();
  }));
  $("criterion").textContent = criteriaText();
}

function renderRunDiagnostics() {
  const joint = state.result?.joint_v4_ml;
  if (!joint?.compute_budget) return "";
  const allocation = joint.compute_budget.foldx_allocation;
  const repeats = joint.compute_budget.foldx_repeat_budget;
  const compare = joint.selected_panel?.interface_comparison;
  const identification = joint.selected_panel?.identifiability_diagnostics;
  const row = (label, value) => `<tr><th>${h(label)}</th><td>${h(String(value ?? "—"))}</td></tr>`;
  const panelText = compare ? `<details><summary>${h(recommendationLabel("全表面 / 界面聚焦同池对照", "Full-surface / interface-focused comparison"))}</summary>
    <p>${h(recommendationLabel("只改变选择目标，不新增模型计算；面积和邻域覆盖为几何代理。", "Changes only the selection objective; no new inference. Area and neighborhood coverage are geometric proxies."))}</p>
    <table><tbody>${row("全表面候选 / Full-surface", (compare.full_surface_candidate_ids || []).join(", "))}${row("界面聚焦候选 / Interface-focused", (compare.candidate_ids || []).join(", "))}${row("交集 / Intersection", (compare.intersection_candidate_ids || []).join(", "))}${row("全表面：界面直接面积 Å²", compare.full_surface_metrics?.direct_interface_area_a2)}${row("界面聚焦：界面直接面积 Å²", compare.direct_interface_area_a2)}${row("合并去重突变构建 / Unique mutant constructs", compare.procurement?.deduplicated_mutant_constructs)}${row("WT对照（实验对照另计）", compare.procurement?.wt_controls)}</tbody></table>
    <details><summary>9/12/15 Å、分区覆盖与证据 / Coverage and evidence</summary><pre>${h(JSON.stringify(compare, null, 2))}</pre></details></details>` : "";
  return `<section class="result-section"><details><summary>${h(recommendationLabel("计算预算、重复诊断与目标对照", "Compute budget, repeat diagnostics and objective comparison"))}</summary>
    <table><tbody>${row("FoldX唯一组合 / Unique combinations", allocation?.selected?.total)}${row("局部 / Local", allocation?.selected?.local)}${row("全局 / Global", allocation?.selected?.global)}${row("局部强扰动探索 / Local exploration", allocation?.selected?.local_strong_exploration)}${row("前三次配对槽 / Initial pair slots", repeats?.initial_slots)}${row("追加两次配对槽 / Additional pair slots", repeats?.supplemental_slots)}${row("总配对槽（上限300）/ Total pair slots", repeats?.total_slots)}</tbody></table>
    ${panelText}<details><summary>${h(recommendationLabel("近邻与精确重复响应诊断", "Nearby and exactly identical response diagnostics"))}</summary><pre>${h(JSON.stringify(identification || {}, null, 2))}</pre></details>
    <details><summary>${h(recommendationLabel("逐候选模型差、重复值及资格", "Per-candidate differences, repetitions and eligibility"))}</summary><pre>${h(JSON.stringify(joint.model_evidence_diagnostics || {}, null, 2))}</pre></details>
    <p>${h(recommendationLabel("模型差和重复离散度仅作诊断，不是实验置信区间或折叠/结合保证。", "Model differences and repeat spread are diagnostics, not experimental confidence intervals or folding/binding guarantees."))}</p></details></section>`;
}

function renderStep5() {
  const budget = state.result.design_summary?.panel_budget_decision;
  const budgetReasons = {
    NO_FEASIBLE_NONREDUNDANT_CANDIDATES: "没有更多合格且提供新信息的候选",
    NO_POSITIVE_NET_INFORMATION_GAIN: "继续增加构建体没有正的净信息增益",
    NO_POSITIVE_COVERAGE_GAIN: "继续增加构建体没有新增覆盖",
    COVERAGE_GAIN_BELOW_STOPPING_POLICY: "新增覆盖不足，按自适应规则停止",
    INFORMATION_GAIN_BELOW_STOPPING_POLICY: "新增信息不足，按自适应规则停止",
    MAXIMUM_CONSTRUCT_BUDGET_REACHED: "已达到构建体硬上限",
  };
  const coveragePolicy = (state.result.joint_v4_ml?.selected_panel?.sector_recommendations || []).length > 0;
  $("design-list").previousElementSibling.textContent = coveragePolicy
    ? recommendationLabel("各分区按设计原则独立给出通过门槛的推荐，数量不受全局精简面板预算限制。选择组合可联动结构与序列；未通过项单列审查详情。", "Each region independently recommends combinations that pass the gates under each design principle, without the global compact-panel limit. Select a combination to inspect its structure and sequence; failed candidates appear in review details.")
    : t("design.note", "选择推荐组合以查看结构位置、序列及评分证据。");
  const minimumNote = coveragePolicy ? "；可用合格候选不足，未达到最少目标" : "；合格候选或信息增益不足，未达到最少目标";
  const budgetText = budget ? `入选 ${budget.selected_constructs} 个；目标 ${budget.minimum_constructs}–${budget.maximum_constructs}，常见目标 ${budget.target_constructs}。${budgetReasons[budget.stop_reason] || budget.stop_reason}${budget.minimum_target_met ? "" : minimumNote}` : "";
  const designs = [...state.result.designs].sort((left, right) =>
    (left.region_id ?? 999) - (right.region_id ?? 999) || left.id.localeCompare(right.id));
  const renderCard = (design, patch) => {
    const info = scoreInfo(design);
    const chips = (design.mutations || []).map((mutation) =>
      `<span class="mut-chip" data-design-residue="${mutation.position}" title="${h(t("design.locate", "定位"))} ${h(notation(mutation))}">${h(notation(mutation))}</span>`
    ).join("") || `<span class="empty">${h(t("design.noAuto", "未生成自动组合"))}</span>`;
    const scoreLine = info.real
      ? `<span class="score-real">${fmt(info.value, 3)} kcal/mol</span> · ${h(info.engine)}`
      : `<span class="tag warn">${h(t("design.proxyHidden", "无真实能量评分 · 不显示代理分"))}</span>`;
    const role = design.design_role
      ? `<span class="tag">${h(t("role." + design.design_role, design.design_role))}</span>`
      : "";
    const recommendationClass = design.foldx_diagnostics?.recommendation_class;
    const classLabel = recommendationClass === "EPITOPE_PRIORITIZED"
      ? t("design.epitopePrioritized", "DiscoTope 优先")
      : recommendationClass === "STABILITY_FALLBACK"
        ? t("design.stabilityFallback", "结构稳定回退") : "";
    const classTag = classLabel ? `<span class="tag">${h(classLabel)}</span>` : "";
    const riskLabel = recommendationRiskLabelText({}, design);
    const repeatNeedsReview = recommendationNeedsReview({}, design);
    const repeatSummary = recommendationRepeatSummary({}, design);
    const statusTag = repeatNeedsReview
      ? `<span class="tag warn">${h(riskLabel)}</span>`
      : `<span class="st">${h(design.status)}</span>${riskLabel ? `<span class="tag">${h(riskLabel)}</span>` : ""}`;
    const actionable = Boolean(design.mutations?.length);
    const cardAction = actionable ? ` data-view-design="${h(design.id)}" tabindex="0" role="button"` : "";
    return `<div class="design-card ${actionable ? "clickable" : ""} ${state.selectedDesign === design.id ? "is-selected" : ""}"${cardAction}>
      <div class="design-head">
        <div><strong>${h(t("region.patch", "补丁"))} ${h(patch)} · ${h(design.id)}</strong><br><span class="st">${h(design.source_engine === "joint_v4_ml" ? (design.mutations || []).map(notation).join(" + ") : design.name)} · ${scoreLine}</span></div>
        <span class="st">${classTag} ${role} ${statusTag}</span>
      </div>
      <div>${chips}</div>
      ${repeatSummary ? `<div class="muted">${h(repeatSummary)}</div>` : ""}
      <div class="actions-row">
        ${design.mutations?.length ? `<button type="button" class="btn btn-mini" data-review-design="${h(design.id)}">${h(t("design.review", "审查此构建体"))}</button>
        <button type="button" class="btn btn-mini" data-copy-design="${h(design.id)}">${h(t("design.copyCustom", "复制到自定义"))}</button>` : ""}
      </div>
    </div>`;
  };
  const byPatch = new Map();
  designs.forEach((design) => {
    const regions = regionsForDesign(design);
    for (const patch of regions.length ? regions.map((region) => region.id) : ["—"]) {
      if (!byPatch.has(patch)) byPatch.set(patch, []);
      byPatch.get(patch).push(design);
    }
  });
  const recommendations = state.result.joint_v4_ml?.selected_panel?.sector_recommendations || [];
  const sectorHtml = recommendations.length ? renderSectorRecommendations(recommendations) : "";
  const legacySeries = recommendations.length ? "" : ([...byPatch.entries()].map(([patch, series]) => `
    <section class="design-series">
      <h3>${h(t("region.patch", "补丁"))} ${h(patch)} · ${series.length} ${h(t("design.combinations", "个推荐组合"))}</h3>
      ${series.map((design) => renderCard(design, patch)).join("")}
    </section>`).join("") || `<div class="empty">${h(t("design.none", "没有推荐组合。"))}</div>`);
  const budgetHtml = budgetText ? coveragePolicy
    ? `<details class="compact-panel-budget"><summary>${h(recommendationLabel("附加全局精简面板（不限制逐区推荐）", "Additional compact panel (does not limit regional recommendations)"))}</summary><p>${h(budgetText)}</p></details>`
    : `<div class="note">${h(budgetText)}</div>` : "";
  $("design-list").innerHTML = sectorHtml + legacySeries + budgetHtml + renderRunDiagnostics();
  document.querySelectorAll("[data-copy-candidate]").forEach((el) => {
    el.addEventListener("click", async () => {
      try { await navigator.clipboard.writeText(el.dataset.copyCandidate); showToast(recommendationLabel("追踪编号已复制", "Trace identifier copied")); }
      catch { showToast(recommendationLabel("复制失败，请手动选择编号", "Copy failed; select the identifier manually")); }
    });
  });
  document.querySelectorAll("[data-design-residue]").forEach((el) => {
    el.addEventListener("click", (event) => { event.stopPropagation(); selectResidue(Number(el.dataset.designResidue), true); });
  });
  document.querySelectorAll("[data-view-design]").forEach((el) => {
    el.addEventListener("click", () => selectDesign(el.dataset.viewDesign));
    el.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") { event.preventDefault(); selectDesign(el.dataset.viewDesign); }
    });
  });
  document.querySelectorAll("[data-review-design]").forEach((el) => {
    el.addEventListener("click", (event) => { event.stopPropagation(); selectDesign(el.dataset.reviewDesign); });
  });
  document.querySelectorAll("[data-select-recommendation]").forEach((el) => {
    el.addEventListener("click", (event) => { event.stopPropagation(); toggleRecommendationSelection(el.dataset.selectRecommendation); });
  });
  document.querySelectorAll("[data-remove-recommendation]").forEach((el) => {
    el.addEventListener("click", () => toggleRecommendationSelection(el.dataset.removeRecommendation));
  });
  document.querySelectorAll("[data-clear-recommendations]").forEach((el) => el.addEventListener("click", clearRecommendationSelections));
  document.querySelectorAll("[data-copy-design]").forEach((el) => {
    el.addEventListener("click", (event) => { event.stopPropagation(); copyDesignToCustom(el.dataset.copyDesign); });
  });
  renderExports();
  renderDesignSelect();
}

function recommendationItems() {
  return (state.result?.joint_v4_ml?.selected_panel?.sector_recommendations || []).flatMap((item) => [
    item,
    ...(item.alternatives || []).map((candidate) => ({ ...candidate, sector_id: item.sector_id })),
    ...(item.all_formal_local_candidates || []).map((candidate) => ({ ...candidate, sector_id: item.sector_id })),
  ]);
}

async function toggleRecommendationSelection(id) {
  if (state.recommendationSelections.has(id)) state.recommendationSelections.delete(id);
  else {
    const item = recommendationItems().find((candidate) => String(candidate.candidate_id || candidate.design_id) === id);
    if (!item) return;
    const design = state.result.designs.find((candidate) => candidate.id === item.design_id);
    const selection = recommendationSelectionRecord(item, design);
    if (!selection) return;
    state.recommendationSelections.set(selection.id, selection);
  }
  renderStep5();
  if (!molstarReady) return;
  await paintMolstar();
  if (!state.recommendationSelections.size) return;
  const positions = new Set();
  for (const selection of state.recommendationSelections.values()) {
    const design = state.result.designs.find((candidate) => candidate.id === selection.design_id);
    const records = design ? designSurfaceRecords(design).envelope : [];
    (records.length ? records : state.result.residues.filter((residue) => selection.mutations.some((mutation) => mutation.position === residue.uniprot_pos)))
      .forEach((residue) => positions.add(residue.uniprot_pos));
  }
  await MolstarView.focus(state.result.residues.filter((residue) => positions.has(residue.uniprot_pos)).map(residueQuery));
}

function clearRecommendationSelections() {
  state.recommendationSelections.clear();
  renderStep5();
  if (molstarReady) paintMolstar().catch((error) => viewerError(error, "paint"));
}

function renderStep6() {
  renderCustomBuilder();
}

function renderStep7() {
  const result = state.result;
  const protein = result.protein || {};
  const source = result.source || {};
  const selection = result.structure_selection || {};
  const summary = result.summary || {};
  const partition = result.partitioning || {};
  const request = result.request || {};
  const v4 = result.joint_v4_ml || {};
  const rangeLabel = request.analysis_target === "binder"
    ? t("overview.analysisChainRange", "分析链范围") : t("overview.ecd", "ECD 范围");
  const files = Object.keys(result.files || {});
  const row = (label, value) => `<dt>${h(label)}</dt><dd>${h(value ?? "—")}</dd>`;
  $("task-overview").innerHTML = `<div class="card task-overview-card">
    <h3>${h(taskName())}</h3><dl class="kv task-kv">
    ${row(t("overview.accession", "UniProt / 输入标识"), protein.accession || request.identifier)}
    ${row(rangeLabel, protein.ecd_start && protein.ecd_end ? `${protein.ecd_start}–${protein.ecd_end} (${protein.ecd_sequence?.length || summary.mapped_ecd_residues || "—"} aa)` : "—")}
    ${row(t("overview.topology", "拓扑注释"), protein.topology)}
    ${row(t("overview.assayProfile", "固定实验方位"), result.assay_profile?.display_name || "—")}
    ${row(t("overview.orientation", "PDB 方位证据"), result.assay_geometry?.status === "relative_axis_only" ? t("overview.relativeAxis", "C 端锚点相对轴（不是物理膜平面）") : "—")}
    ${row(t("overview.structure", "使用的结构"), selection.selected_identifier || source.identifier)}
    ${row(t("overview.provider", "结构来源 / 类型"), [source.provider, source.kind].filter(Boolean).join(" / "))}
    ${row(t("overview.chain", "目标链"), selection.selected_chain || selection.chain || "—")}
    ${row(t("overview.confidence", "结构置信度"), result.methods?.confidence_method || (request.input_type === "af3_upload" ? "AF3 confidence JSON" : "—"))}
    ${row(t("overview.run", "运行信息"), result.created_at || "—")}
    ${row(t("overview.structureFile", "结构文件名"), source.identifier || selection.selected_identifier || "—")}
    ${row(t("overview.summary", "结果摘要"), `${summary.mapped_ecd_residues ?? "—"} ECD residues · ${partition.actual_regions ?? result.regions.length} patches · ${result.designs.length} ${state.lang === "zh" ? "设计记录（含控制与逐区建议）" : "design records (including controls and sector suggestions)"}`)}
    ${v4.recommendation_status ? row("Joint V4-ML 状态", v4.recommendation_status) : ""}
    ${v4.sectors ? row("V4 补丁 / 主面板组合", `${v4.sectors.length} / ${(v4.selected_panel?.selected_candidate_ids || []).length}`) : ""}
    ${v4.selected_panel?.sector_recommendations?.length ? row(state.lang === "zh" ? "逐区建议（通过 / 待审查 / 无组合）" : "Sector suggestions (passed / review / no combination)", ["PASSED", "NEAR_THRESHOLD", "NO_COMBINATION"].map((status) => v4.selected_panel.sector_recommendations.filter((item) => item.status === status).length).join(" / ")) : ""}
    ${row(t("overview.coverage", "表面覆盖率"), partition.coverage != null ? fmtPct(partition.coverage) : "—")}
    ${row(t("overview.files", "导出文件"), files.length ? `${files.length} files: ${files.join(", ")}` : "—")}
    </dl></div>`;
}

/* ---------------- residues ---------------- */

function filteredResidues() {
  if (!state.result) return [];
  const filter = $("residue-filter").value;
  const query = $("residue-search").value.trim().toLowerCase();
  return state.result.residues.filter((residue) => {
    let matches = true;
    if (filter === "folding") matches = residue.folding_screen_pass;
    if (filter === "surface") matches = residue.surface_accessible;
    if (filter === "interface") matches = residue.interface;
    if (filter === "protected") matches = residue.hard_exclusions?.length > 0;
    if (!matches || !query) return matches;
    const text = [
      residue.uniprot_pos, `${residue.aa}${residue.uniprot_pos}`, residue.domain,
      residue.interface ? "interface" : "", ...(residue.risk_flags || []), ...(residue.hard_exclusions || []),
    ].join(" ").toLowerCase();
    return text.includes(query);
  });
}

function renderResidueList() {
  if (!state.result || !$("residue-list")) return;
  const list = filteredResidues();
  const statusWord = (residue) => residue.folding_screen_pass
    ? t("rd.pass", "通过")
    : residue.surface_accessible ? t("rd.risk", "风险") : t("rd.excluded", "排除");
  $("residue-list").innerHTML = list.map((residue) => {
    const active = state.selectedResidue?.uniprot_pos === residue.uniprot_pos ? " active" : "";
    return `<button type="button" class="resrow${active}" data-residue="${residue.uniprot_pos}">
      <span class="id">${h(residue.aa)}${residue.uniprot_pos}</span>
      <span class="main">${h(residue.domain)} · ${h(t("region.patch", "补丁"))} ${residue.region_id ?? "—"} · ${h(statusWord(residue))}${residue.interface ? ` · ${h(t("rd.interfaceShort", "界面"))}` : ""}${residue.engineered_neighborhood ? ` · ${h(t("rd.engineered", "结构差异邻域"))}` : ""}</span>
      <span class="muted">rSASA ${fmt(residue.rsasa)}</span>
    </button>`;
  }).join("") || `<div class="empty">${h(t("rd.listEmpty", "当前筛选没有残基。"))}</div>`;
  document.querySelectorAll("[data-residue]").forEach((b) => {
    b.addEventListener("click", () => selectResidue(Number(b.dataset.residue), true));
  });
}

function foldxDdgFor(position) {
  for (const design of state.result?.designs || []) {
    const scan = design.saturation_scan?.[String(position)];
    if (scan) {
      const hit = scan.find((c) => c.foldx_ddg != null);
      return hit ? hit.foldx_ddg : null;
    }
  }
  return null;
}

/* Rebuild guard: the detail panel owns <select> elements the user interacts
   with; rebuilding it on every stray re-render used to rip the open dropdown
   out of the DOM.  Only rebuild when something displayed actually changed. */
let detailRenderKey = null;

function renderResidueDetail(force = false) {
  const residue = state.selectedResidue;
  const el = $("residue-detail");
  if (!residue) {
    detailRenderKey = "empty";
    el.className = "empty";
    el.textContent = t("s1.detailEmpty", "在结构、序列或列表中选择残基。");
    return;
  }
  const key = [
    residue.uniprot_pos, state.lang, state.result.designs.length, state.customMutations.length,
  ].join("|");
  if (!force && key === detailRenderKey) return;
  detailRenderKey = key;

  el.className = "";
  const tags = [];
  tags.push(`<span class="tag ${residue.folding_screen_pass ? "ok" : "bad"}">${residue.folding_screen_pass ? h(t("rd.foldingPass", "折叠风险通过")) : h(t("rd.foldingFail", "未通过自动筛选"))}</span>`);
  if (residue.interface) tags.push(`<span class="tag warn">${h(t("rd.interface", "功能/结构界面"))}</span>`);
  if (residue.glycan_shielded) tags.push(`<span class="tag warn">${h(t("rd.glycan", "糖链遮蔽"))}</span>`);
  if (residue.disulfide) tags.push(`<span class="tag bad">${h(t("rd.disulfide", "二硫键"))}</span>`);
  (residue.risk_flags || []).forEach((f) => tags.push(`<span class="tag">${h(f)}</span>`));
  (residue.hard_exclusions || []).forEach((f) => tags.push(`<span class="tag bad">${h(f)}</span>`));

  const choices = scanChoicesForResidue(residue)
    .filter((c) => c.allowed)
    .sort((a, b) => (a.foldx_ddg ?? 99) - (b.foldx_ddg ?? 99))
    .slice(0, 5);
  const options = AA.filter((aa) => aa !== residue.aa).map((aa) => `<option value="${aa}">${aa}</option>`).join("");

  el.innerHTML = `
    <div style="display:flex;justify-content:space-between;gap:8px;align-items:flex-start">
      <h3 style="margin:0">${h(residue.aa)}${residue.uniprot_pos}
        <span class="muted" style="font-size:11.5px">PDB ${h(residue.chain_id)}:${residue.auth_seq_id}${h(residue.insertion_code || "")} · ${h(residue.residue_name)}</span></h3>
      <button class="btn btn-mini" id="focus-residue" type="button">${h(t("rd.locate", "结构定位"))}</button>
    </div>
    <div style="margin:6px 0">${tags.join("")}</div>
    <dl class="kv">
      <dt>rSASA</dt><dd>${fmt(residue.rsasa, 3)}（${h(t("rd.total", "总"))} ${fmt(residue.sasa, 0)} Å² / ${h(t("rd.sidechain", "侧链"))} ${fmt(residue.sidechain_sasa, 0)} Å²）</dd>
      <dt>表面暴露</dt><dd>${residue.surface_accessible ? "通过" : "未通过"}</dd>
      <dt>坐标资格</dt><dd>${residue.coordinate_eligible === true ? "可用于几何" : residue.coordinate_eligible === false ? "不用于三维足迹" : "历史运行未记录"} · ${h(residue.coordinate_status || "历史运行未记录")}</dd>
      <dt>位点初筛</dt><dd>${h(residue.mutation_qualification || "历史运行未记录")} ${h((residue.review_reasons || []).join("；"))}；具体替换仍需独立审核</dd>
      <dt>${h(t("rd.local", "局部"))}</dt><dd>Cα ${h(t("rd.contacts", "接触"))} ${residue.contacts} · ${h(secondaryName(residue.secondary_structure))}${residue.salt_bridge ? ` · ${h(t("rd.saltBridge", "盐桥"))}` : ""}</dd>
      <dt>${h(t("rd.interfaceLabel", "界面"))}</dt><dd>${residue.interface ? `${h(t("iface.with", "与"))} ${h((residue.interface_partners || []).join(", ") || t("rd.nonTarget", "非目标链"))} ${h(t("rd.contact", "接触"))}` : h(t("rd.none", "无"))}</dd>
      <dt>${h(t("insp.region", "分区"))}</dt><dd>${h(designForResidue(residue).text)}</dd>
    </dl>
    <div style="margin-top:8px"><h3>${h(t("rd.choices", "候选替换（ΔΔG 仅显示 FoldX 真值）"))}</h3>
      ${choices.length ? choices.map((c) => {
        const sc = singleChoiceScore(c);
        return `<div class="rowline"><code>${h(c.notation || `${c.wt}${c.position}${c.mutant}`)}</code>
          <span class="muted">${sc.real ? `FoldX ${fmt(sc.value, 3)} kcal/mol` : h(t("s2.unscored", "未评分"))}</span>
          <span class="grow"></span>
          <button type="button" class="btn btn-mini" data-add-choice="${h(c.mutant)}">${h(t("insp.add", "加入组合"))}</button></div>`;
      }).join("") : `<div class="empty">${h(t("rd.noScan", "当前结构未保存该位点的替换扫描。"))}</div>`}
    </div>
    <div class="rowline" style="margin-top:6px">
      <select id="custom-aa" class="fld-inline">${options}</select>
      <button id="add-custom-site" class="btn btn-mini" type="button">${h(t("rd.addCustom", "加入自定义构建"))}</button>
    </div>`;
  $("focus-residue").addEventListener("click", () => focusSelectedResidue());
  $("add-custom-site").addEventListener("click", () => addCustomMutation(residue, $("custom-aa").value));
  document.querySelectorAll("[data-add-choice]").forEach((b) => {
    b.addEventListener("click", () => addCustomMutation(residue, b.dataset.addChoice));
  });
}

async function selectResidue(position, focus = false) {
  if (!state.result) return;
  const residue = state.result.residues.find((item) => item.uniprot_pos === position);
  if (!residue) return;
  state.selectedResidue = residue;
  renderResidueList();
  renderResidueDetail();
  renderSeqDock();
  if (state.engine === "molstar" && molstarReady) {
    /* A surface showcase must remain surface-only. Outside that mode, retain
       the side-chain marker used by the residue inspector. */
    await MolstarView.highlight(state.designSurfaceMode
      ? residueQuery(residue)
      : {
        ...residueQuery(residue),
        sideChain: true,
        representation: "ball-and-stick",
        representationColor: "#FFD166",
      });
    if (focus) await MolstarView.focus([residueQuery(residue)]);
  }
}

async function focusSelectedResidue() {
  if (state.engine !== "molstar" || !molstarReady || !state.selectedResidue) return;
  const q = residueQuery(state.selectedResidue);
  await MolstarView.focus([q]);
}

/* ---------------- designs ---------------- */

function selectedDesign() {
  if (!state.result || state.selectedDesign === "wt") return null;
  return state.previewPlans.get(state.selectedDesign)
    || state.result.designs.find((d) => d.id === state.selectedDesign)
    || null;
}

function designForResidue(residue) {
  const region = state.result.regions.find((item) => item.id === residue.region_id);
  if (!region) return { text: t("rd.unassigned", "未分配至自动分区") };
  return { text: `${t("region.patch", "补丁")} R${region.id}（${region.tier || ""}）` };
}

function scanChoicesForResidue(residue) {
  for (const design of state.result.designs || []) {
    const scan = design.saturation_scan?.[String(residue.uniprot_pos)];
    if (scan) return scan;
  }
  return [];
}

function renderDesignSelect() {
  const select = $("design-select");
  if (!select || !state.result) return;
  const previews = [...state.previewPlans.values()];
  select.innerHTML = `<option value="wt">WT ECD</option>` +
    [...state.result.designs, ...previews].map((design) => {
      const regions = regionsForDesign(design);
      const regionLabel = regions.length ? regions.map((region) => `${recommendationLabel("分区", "Region")} ${region.id}`).join("、") : recommendationLabel("全域／自定义", "Global / custom");
      const role = design.design_role === "SECTOR_REVIEW_CANDIDATE" ? recommendationLabel("待审查建议", "Review suggestion")
        : design.design_role === "STAGE2_SINGLE" ? recommendationLabel("单点控制", "Single-site control")
          : design.design_role === "STAGE2_PAIR" ? recommendationLabel("双点控制", "Pair control")
            : design.design_role === "SECTOR_QUALIFIED_CANDIDATE" ? recommendationLabel("逐区建议", "Region suggestion")
              : recommendationLabel("组合构建体", "Combination construct");
      return `<option value="${h(design.id)}">${h(design.mutations.map(notation).join(" + ") || recommendationLabel("无突变", "No mutation"))} · ${h(regionLabel)} · ${h(role)}</option>`;
    }).join("");
  select.value = [...select.options].some((o) => o.value === state.selectedDesign) ? state.selectedDesign : "wt";
}

async function selectDesign(id) {
  const requestId = ++state.designFocusRequest;
  state.selectedDesign = id;
  state.previewMode = "construct";
  state.reliabilityView = false;
  const design = selectedDesign();
  state.designSurfaceMode = Boolean(design?.mutations?.length && regionForDesign(design));
  renderDesignSelect();
  $("design-select").value = id;
  if (id !== "wt") state.layers.add("design");
  renderStep2();
  renderStep5();
  renderSeqDock();
  renderDesignSurfaceCallout();
  if (!state.designSurfaceMode) {
    repaintViews();
    if (id === "wt" && molstarReady) await MolstarView.reset();
    return;
  }
  if (state.engine !== "molstar") setEngine("molstar");
  await setRepresentation("molecular-surface");
  if (requestId !== state.designFocusRequest || !molstarReady) return;
  const records = designSurfaceRecords(design);
  await MolstarView.focus(records.envelope.length ? records.envelope.map(residueQuery) : records.mutations.map(residueQuery));
}

/* ---------------- persistent sequence dock ---------------- */

function currentDesignSequences() {
  const design = selectedDesign();
  const wtAa = state.result.protein.ecd_sequence;
  const wtDna = state.result.protein.ecd_coding_sequence || [...wtAa].map((aa) => CODONS[aa]).join("");
  const map = {};
  if (state.previewMode === "custom" && state.customMutations.length) {
    state.customMutations.forEach((m) => { map[m.position] = m.mutant; });
  } else if (design) {
    design.mutations.forEach((m) => { map[m.position] = m.mutant; });
  }
  const out = wtAa.split("");
  const dnaOut = wtDna.match(/.{1,3}/g) || [];
  Object.keys(map).forEach((p) => {
    const index = Number(p) - state.result.protein.ecd_start;
    if (index >= 0 && index < out.length) {
      out[index] = map[p];
      if (dnaOut[index]) dnaOut[index] = CODONS[map[p]];
    }
  });
  return {
    wtAa, wtDna,
    designAa: out.join(""),
    designDna: dnaOut.join(""),
    map,
    designLabel: state.previewMode === "custom" && state.customMutations.length
      ? `${t("seqdock.custom", "自定义")} · ${state.customMutations.length} ${t("region.sites", "位点")}`
      : design ? `${design.mutations.map(notation).join(" + ") || recommendationLabel("无突变", "No mutation")} · ${regionsForDesign(design).map((region) => `${recommendationLabel("分区", "Region")} ${region.id}`).join("、") || recommendationLabel("全域／自定义", "Global / custom")}` : t("seqdock.wtControl", "WT 对照"),
    hasRealCds: Boolean(design?.mutant_coding_sequence || state.result.protein.ecd_coding_sequence),
  };
}

function renderSeqDock() {
  renderSeqDockHead();
  if (!state.result) {
    $("wt-sequence").innerHTML = `<span class="muted">${h(t("seqdock.loadFirst", "载入结果后显示序列比对。"))}</span>`;
    $("design-sequence").innerHTML = "";
    return;
  }
  const seqs = currentDesignSequences();
  state.lastSequences = seqs;
  const start = state.result.protein.ecd_start;
  const mode = document.querySelector('input[name="seq-mode"]:checked').value;
  const design = selectedDesign();
  $("sequence-summary").textContent = `${seqs.wtAa.length} aa · ${Object.keys(seqs.map).length} ${t("workspace.mutations", "处替换")}`;
  $("sequence-summary").title = seqs.designLabel;

  function track(sequence, isDesign) {
    const mutationPositions = isDesign ? new Set(Object.keys(seqs.map).map(Number)) : new Set();
    const protectedPositions = new Set(
      state.result.residues.filter((r) => r.hard_exclusions?.length).map((r) => r.uniprot_pos)
    );
    const selectedPosition = state.selectedResidue?.uniprot_pos;
    if (mode === "aa") {
      return [...sequence].map((aa, i) => {
        const pos = start + i;
        const cls = ["tok"];
        if (isDesign && mutationPositions.has(pos)) cls.push("mut");
        if (protectedPositions.has(pos)) cls.push("prot");
        if (selectedPosition === pos) cls.push("sel");
        const ddg = isDesign && mutationPositions.has(pos) ? foldxDdgFor(pos) : null;
        const tip = ddg != null ? ` · FoldX ΔΔG ${fmt(ddg, 3)}` : "";
        return `<span class="${cls.join(" ")}" role="button" tabindex="0" aria-label="${aa}${pos}${h(tip)}" data-seq-residue="${pos}" title="${aa}${pos}${tip}">${aa}</span>`;
      }).join("");
    }
    const codons = sequence.match(/.{1,3}/g) || [];
    return codons.map((codon, i) => {
      const pos = start + i;
      const cls = ["tok"];
      if (isDesign && mutationPositions.has(pos)) cls.push("mut");
      if (selectedPosition === pos) cls.push("sel");
      return `<span class="${cls.join(" ")}" role="button" tabindex="0" aria-label="${pos}: ${codon}" data-seq-residue="${pos}" title="${pos}: ${codon}" style="min-width:31px">${codon}</span>`;
    }).join(" ");
  }

  $("wt-sequence").innerHTML = track(mode === "aa" ? seqs.wtAa : seqs.wtDna, false);
  $("design-sequence").innerHTML = track(mode === "aa" ? seqs.designAa : seqs.designDna, true);
  $("wt-meta").textContent = `${start}-${start + seqs.wtAa.length - 1} · ${mode === "aa" ? `${seqs.wtAa.length} aa` : `${seqs.wtDna.length} nt`}`;
  $("design-meta").textContent = seqs.designLabel;
  document.querySelectorAll("#design-sequence [data-seq-residue], #wt-sequence [data-seq-residue]").forEach((token) => {
    /* focus=true: the camera centres on the residue so it can never hide
       behind the rest of the molecule (the "clicked but can't see it" bug). */
    token.addEventListener("click", () => selectResidue(Number(token.dataset.seqResidue), true));
    token.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") { event.preventDefault(); token.click(); }
    });
  });
}

function renderSeqDockHead() {
  const expanded = !$("sequence-content").hidden;
  $("toggle-sequence").setAttribute("aria-expanded", String(expanded));
  $("toggle-sequence").textContent = expanded
    ? t("workspace.collapseSequence", "收起比对") : t("workspace.expandSequence", "展开比对");
}

/* ---------------- custom builder / exports ---------------- */

function addCustomMutation(residue, mutant) {
  const existing = state.customMutations.find((item) => item.position === residue.uniprot_pos);
  if (existing) existing.mutant = mutant;
  else state.customMutations.push({ position: residue.uniprot_pos, wt: residue.aa, mutant });
  state.customMutations.sort((a, b) => a.position - b.position);
  renderCustomBuilder();
  renderSeqDock();
  repaintViews();
  /* Stay on the current step: a forced jump to step 6 made it impossible to
     add a second mutation from the detail panel.  The toast tells the user
     where to find the construct. */
  showToast(`${residue.aa}${residue.uniprot_pos}→${mutant} ${t("msg.customAdded", "已加入组合（第 6 步查看/保存）")}`);
}

function copyDesignToCustom(id) {
  const design = state.previewPlans.get(id) || state.result.designs.find((d) => d.id === id);
  if (!design) return;
  state.customMutations = design.mutations.map((m) => ({ position: m.position, wt: m.wt, mutant: m.mutant }));
  $("custom-name").value = `${design.id} edited`;
  renderCustomBuilder();
  renderSeqDock();
  repaintViews();
  setStep(6);
}

function renderCustomBuilder() {
  const container = $("custom-mutations");
  if (!state.customMutations.length) {
    container.className = "empty";
    container.textContent = t("custom.empty", "还没有选残基。");
    $("save-custom").disabled = true;
    return;
  }
  container.className = "";
  container.innerHTML = state.customMutations.map((mutation) => {
    const options = AA.filter((aa) => aa !== mutation.wt).map((aa) =>
      `<option value="${aa}" ${aa === mutation.mutant ? "selected" : ""}>${aa}</option>`).join("");
    return `<div class="custom-row"><code>${h(mutation.wt)}${mutation.position}</code>
      <select data-custom-select="${mutation.position}">${options}</select>
      <button type="button" class="btn btn-mini" data-custom-remove="${mutation.position}">${h(t("insp.remove", "移除"))}</button></div>`;
  }).join("");
  $("save-custom").disabled = state.offline ? !state.customMutations.length : state.demo || !state.jobId;
  $("custom-feedback").textContent = state.offline
    ? t("offline.draftNote", "离线草稿只导出本地序列与突变清单；不会计算、写回原任务或生成稳定性结论。")
    : state.demo ? t("custom.demoNote", "冻结示例可编辑预览；提交真实作业后可保存导出。") : "";
  document.querySelectorAll("[data-custom-select]").forEach((select) => {
    select.addEventListener("change", () => {
      const mutation = state.customMutations.find((item) => item.position === Number(select.dataset.customSelect));
      if (mutation) mutation.mutant = select.value;
      renderSeqDock();
      repaintViews();
    });
  });
  document.querySelectorAll("[data-custom-remove]").forEach((button) => {
    button.addEventListener("click", () => {
      state.customMutations = state.customMutations.filter((item) => item.position !== Number(button.dataset.customRemove));
      renderCustomBuilder();
      renderSeqDock();
      repaintViews();
    });
  });
}

async function saveCustomDesign() {
  if (state.offline) {
    if (!state.customMutations.length) return;
    renderSeqDock();
    const seqs = state.lastSequences || currentDesignSequences();
    const name = ($("custom-name").value || "offline_draft").trim().replace(/[^\w.-]+/g, "_");
    const draft = {
      schema: "ecd-epitope-offline-draft-v1",
      status: "LOCAL_DRAFT_NOT_ANALYZED",
      note: "No model or server computation was run. This draft is not added to the original task results.",
      protein: state.result.protein,
      mutations: state.customMutations.map(({ position, wt, mutant }) => ({ position, wt, mutant, notation: `${wt}${position}${mutant}` })),
      sequences: { amino_acid: seqs.designAa, coding_dna: seqs.designDna },
    };
    const download = (filename, content, type) => {
      const link = document.createElement("a");
      link.href = URL.createObjectURL(new Blob([content], { type }));
      link.download = filename;
      link.click();
      setTimeout(() => URL.revokeObjectURL(link.href), 1000);
    };
    download(`${name}.draft.json`, JSON.stringify(draft, null, 2), "application/json");
    download(`${name}.fasta`, `>${name}_local_draft_not_analyzed\n${seqs.designAa.match(/.{1,80}/g)?.join("\n") || ""}\n`, "text/plain");
    download(`${name}.cds.fasta`, `>${name}_local_draft_not_analyzed\n${seqs.designDna.match(/.{1,90}/g)?.join("\n") || ""}\n`, "text/plain");
    $("custom-feedback").textContent = t("offline.draftSaved", "本地草稿已下载；未运行模型，也未修改原结果包。请按实验设计复核序列。" );
    return;
  }
  if (!state.jobId || !state.customMutations.length) return;
  $("save-custom").disabled = true;
  $("custom-feedback").textContent = t("custom.saving", "正在验证并生成 ECD 序列…");
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
    renderAll(false);
    $("custom-feedback").textContent = t("custom.saved", `${response.design.id} 已保存并进入导出文件。`);
    showToast(t("custom.savedToast", "自定义构建已保存"));
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
    container.className = "export-links empty";
    container.textContent = t("s5.noFiles", "当前结果没有可下载文件。");
    return;
  }
  const labels = {
    "ecd_epitope_design.zip": t("export.zip", "完整结果 ZIP"), "analysis.json": t("export.analysis", "分析 JSON"),
    "residues.csv": t("export.residues", "逐残基 CSV"), "designs.csv": t("export.designs", "构建 CSV"),
    "panel_coverage_matrix.csv": t("export.coverage", "覆盖矩阵 CSV"),
    "sequences.fasta": t("export.aaFasta", "氨基酸 FASTA"), "coding_sequences.fasta": t("export.cdsFasta", "核苷酸 FASTA"),
    "full_length_protein_sequences.fasta": t("export.fullProtein", "全长蛋白 FASTA"),
    "full_length_coding_sequences.gb": t("export.fullCds", "全长编码 GenBank"),
    "sequences.gb": t("export.genbank", "GenBank"), "methods_report.md": t("export.methods", "方法报告"),
    "run_manifest.json": t("export.manifest", "复现 manifest"), "mapping.csv": t("export.mapping", "编号映射 CSV"),
  };
  const order = ["ecd_epitope_design.zip", "methods_report.md", "designs.csv", "panel_coverage_matrix.csv", "sequences.fasta",
    "coding_sequences.fasta", "full_length_protein_sequences.fasta", "full_length_coding_sequences.gb", "sequences.gb",
    "residues.csv", "mapping.csv", "analysis.json", "run_manifest.json"];
  const sorted = [...files].sort((a, b) => order.indexOf(a) - order.indexOf(b));
  container.className = "export-links";
  container.innerHTML = sorted.map((file) => {
    const url = state.offline ? `./${file}` : state.demo ? `/assets/demo/${file}` : `/api/jobs/${state.jobId}/files/${file}`;
    return `<a href="${h(url)}" download>${h(labels[file] || file)}</a>`;
  }).join("");
}

/* ---------------- repaint + events ---------------- */

function repaintViews() {
  if (state.engine === "molstar") {
    const loadId = viewerLoadId;
    paintMolstar().catch((error) => { if (loadId === viewerLoadId) viewerError(error, "paint"); });
  }
  else paintBeads();
  renderLegend();
}

function secondaryName(code) {
  return ({
    H: t("sse.helix", "α-螺旋"), E: t("sse.strand", "β-链"), C: t("sse.loop", "环 / 无规则"),
  })[code] || t("sse.none", "未注释");
}

const notation = (mutation) => `${mutation.wt}${mutation.position}${mutation.mutant}`;

function bindEvents() {
  $("open-run").addEventListener("click", () => $("run-dialog").showModal());
  $("close-run").addEventListener("click", () => $("run-dialog").close());
  $("toggle-sequence").addEventListener("click", () => {
    $("sequence-content").hidden = !$("sequence-content").hidden;
    renderSeqDockHead();
  });
  $("steps").addEventListener("keydown", (event) => {
    if (!event.target.matches("[data-step]") || !["ArrowLeft", "ArrowRight"].includes(event.key)) return;
    event.preventDefault();
    setStep(Math.min(STEP_LIST.length, Math.max(1, state.step + (event.key === "ArrowRight" ? 1 : -1))));
    document.querySelector(`#steps [data-step="${state.step}"]`).focus();
  });
  $("run-form").addEventListener("submit", submitRun);
  $("load-demo").addEventListener("click", loadDemo);
  $("retry-viewer").addEventListener("click", () => {
    if (state.result) ensureViewers(state.result.structure_url || "/assets/demo/ecd_structure.cif");
  });
  $("open-upload").addEventListener("click", () => { $("run-dialog").close(); updateUploadTargetUI(); $("upload-dialog").showModal(); });
  $("close-upload").addEventListener("click", () => $("upload-dialog").close());
  $("upload-form").addEventListener("submit", submitUpload);
  $("upload-source-mode").addEventListener("change", updateUploadTargetUI);
  $("upload-analysis-target").addEventListener("change", updateUploadTargetUI);
  $("close-af3").addEventListener("click", () => $("af3-dialog").close());
  $("af3-form").addEventListener("submit", submitAf3);
  $("af3-source-mode").addEventListener("change", updateAf3SourceUI);
  updateAf3SourceUI();
  /* CLICK ONLY.  Binding mouseover here used to re-render the residue panel on
     every hover, destroying <select> elements mid-interaction (the "custom
     design can't pick anything" bug).  Tooltips cover hover already. */
  $("viewer-molstar").addEventListener("PDB.molstar.click", handleMolstarPick);
  document.querySelectorAll("[data-engine]").forEach((b) => b.addEventListener("click", () => setEngine(b.dataset.engine)));
  document.querySelectorAll("[data-representation]").forEach((b) => b.addEventListener("click", () => setRepresentation(b.dataset.representation)));
  $("btn-reset-view").addEventListener("click", () => {
    if (state.engine === "beads") BeadsView.reset();
    else MolstarView.reset();
  });
  $("btn-spin").addEventListener("click", async () => {
    state.spinning = !state.spinning;
    $("btn-spin").setAttribute("aria-pressed", String(state.spinning));
    if (state.engine === "molstar") await MolstarView.spin(state.spinning);
  });
  $("btn-lang").addEventListener("click", toggleLang);
  $("btn-export-offline").addEventListener("click", () => {
    if (!state.jobId || state.demo) { showToast(t("msg.exportNeedJob", "提交真实作业并完成后才能导出离线报告")); return; }
    window.open(`/api/jobs/${state.jobId}/export/offline?lang=${state.lang}`, "_blank");
  });
  const copyText = (text) => {
    if (!text) return;
    if (navigator.clipboard?.writeText) {
      navigator.clipboard.writeText(text).then(() => showToast(t("msg.copied", "已复制到剪贴板")), () => showToast(t("msg.copyFail", "复制失败，请手动复制")));
    } else showToast(t("msg.copyFail", "复制失败，请手动复制"));
  };
  $("copy-wt-aa").addEventListener("click", () => copyText(state.lastSequences?.wtAa));
  $("copy-wt-nt").addEventListener("click", () => copyText(state.lastSequences?.wtDna));
  $("copy-design-aa").addEventListener("click", () => copyText(state.lastSequences?.designAa));
  $("copy-design-nt").addEventListener("click", () => copyText(state.lastSequences?.designDna));
  document.querySelectorAll('input[name="seq-mode"]').forEach((input) => input.addEventListener("change", renderSeqDock));
  document.querySelectorAll('input[name="input-type"]').forEach((input) => input.addEventListener("change", () => {
    const uniprot = input.value === "uniprot" && input.checked;
    $("identifier-label").textContent = uniprot ? t("form.identifier", "蛋白 accession") : t("form.pdbId", "PDB ID");
    $("identifier").value = uniprot ? "Q9NZQ7" : "4ZQK";
  }));
  $("viewer-beads").addEventListener("mousemove", (e) => {
    const rect = $("viewer-beads").getBoundingClientRect();
    beadsTipXY = { x: e.clientX - rect.left, y: e.clientY - rect.top };
  });
  $("viewer-beads").addEventListener("mouseleave", hideBeadsTip);
}

function handleMolstarPick(event) {
  const detail = event.eventData || event.detail || {};
  if (!state.result || detail.auth_seq_id == null) return;
  const residue = state.result.residues.find((item) =>
    item.auth_seq_id === Number(detail.auth_seq_id) &&
    (!detail.auth_asym_id || item.chain_id === detail.auth_asym_id) &&
    (!detail.ins_code || item.insertion_code === detail.ins_code)
  );
  if (residue) selectResidue(residue.uniprot_pos, false);
}

(async function boot() {
  await loadI18n();
  applyStaticI18n();
  renderPanelScaffold();
  bindEvents();
  renderLegend();
  if (state.offline) {
    document.body.classList.add("offline-mode");
    const snapshot = window.APP_OFFLINE.analysis;
    snapshot.structure_url = "";
    adoptResult(snapshot, { demo: false, jobId: null });
    setJob("COMPLETE", t("offline.loaded", "离线任务结果已载入；仅查看与本地序列草稿，不运行计算。"), 1);
    return;
  }
  const savedJobId = new URLSearchParams(window.location.search).get("job");
  if (savedJobId) loadSavedJob(savedJobId);
  else setJob("READY", "", 0);

  /* keep the 3D canvas in sync when the viewport box changes without a
     window resize (sequence dock expanding, panel resize) — otherwise the
     structure gets clipped at the stale canvas edge and the page background
     shows through. */
  let resizeTimer = null;
  new ResizeObserver(() => {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(() => {
      if (state.engine === "beads") BeadsView.resize();
      else MolstarView.resize();
    }, 120);
  }).observe($("viewport"));
})();
