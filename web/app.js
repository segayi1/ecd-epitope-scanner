/* ECD Epitope Workbench — frontend (rewrite v4: full-page i18n + interaction fixes)
 *
 * Contract with the backend (ecd_epitope_scanner api.py — unchanged):
 *   POST /api/jobs            {input_type, identifier, chain, strict_af3, allow_partial}
 *   POST /api/jobs/upload     multipart(structure, accession?, chain?, sequence?, ecd_start?, ecd_end?)
 *   GET  /api/jobs/{id}       -> {status, message, progress, result}
 *   GET  /api/jobs/{id}/events        SSE progress
 *   POST /api/jobs/{id}/af3           multipart(structure, confidence)
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
  representation: "cartoon",
  spinning: false,
  layers: new Set(["folding", "regions", "design"]),
  layerFocus: "folding",
  selectedResidue: null,
  selectedDesign: "wt",
  previewPlans: new Map(),
  customMutations: [],
  lastSequences: null,
};

const STEP_LIST = ["step.surface", "step.tolerance", "step.regions", "step.interface", "step.design", "step.custom"];
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
  for (const lang of ["zh", "en"]) {
    try { state.i18n[lang] = await api(`/assets/i18n/${lang}.json`); }
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
  localStorage.setItem("ecd.lang", state.lang);
  applyStaticI18n();
  renderPanelScaffold();   // rebuilds the six step sections in the new language
  if (state.result) renderAll(false);
  else renderLegend();
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

function setJob(status, message, progress = 0, error = null) {
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
}

function adoptResult(result, { demo, jobId }) {
  state.result = result;
  state.previewPlans.clear();
  state.demo = demo;
  state.jobId = jobId;
  state.selectedResidue = result.residues.find((r) => r.folding_screen_pass) || result.residues[0] || null;
  state.selectedDesign = result.designs.find((d) => d.mutations?.length)?.id || "wt";
  if (!result.protein.ecd_coding_sequence) {
    result.protein.ecd_coding_sequence = [...result.protein.ecd_sequence].map((aa) => CODONS[aa]).join("");
  }
  renderAll();
  const url = result.structure_url || "/assets/demo/ecd_structure.cif";
  ensureViewers(url);
}

function loadDemo() {
  setJob("RUNNING", t("msg.loadingDemo", "载入冻结 PD-L1 审查结果"), 0.4);
  return api("/assets/demo/analysis.json").then((result) => {
    result.structure_url = "/assets/demo/ecd_structure.cif";
    adoptResult(result, { demo: true, jobId: null });
    setJob("DEMO", t("msg.demoReady", "PD-L1 冻结示例；提交真实作业后可保存自定义组合与导出"), 1);
  }).catch((e) => setJob("FAILED", t("msg.demoFail", "示例载入失败"), 0, e.message));
}

function watchJob(jobId) {
  state.eventSource?.close();
  const source = new EventSource(`/api/jobs/${jobId}/events`);
  state.eventSource = source;
  source.onmessage = async (event) => {
    try {
      const u = JSON.parse(event.data);
      setJob(u.status, u.message, u.progress, u.error);
      if (["COMPLETE", "FAILED", "WAITING_FOR_AF3_UPLOAD"].includes(u.status)) {
        source.close();
        const job = await api(`/api/jobs/${jobId}`);
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
    if (source.readyState === EventSource.CLOSED) return;
    setJob("FAILED", t("msg.sseLost", "作业事件连接中断"), 0, t("msg.sseLostHint", "请刷新任务状态或重新提交。"));
  };
}

async function submitRun(event) {
  event.preventDefault();
  const inputType = document.querySelector('input[name="input-type"]:checked').value;
  try {
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
      }),
    });
    state.jobId = job.id;
    state.demo = false;
    setJob(job.status, job.message, job.progress);
    watchJob(job.id);
  } catch (error) {
    setJob("FAILED", t("msg.submitFail", "作业提交失败"), 0, error.message);
  }
}

async function submitUpload(event) {
  event.preventDefault();
  const structure = $("upload-structure").files[0];
  if (!structure) return;
  $("upload-dialog").close();
  const form = new FormData();
  form.append("structure", structure);
  if ($("upload-accession").value.trim()) form.append("accession", $("upload-accession").value.trim());
  if ($("upload-chain2").value.trim()) form.append("chain", $("upload-chain2").value.trim());
  if ($("upload-sequence").value.trim()) form.append("sequence", $("upload-sequence").value.trim());
  if ($("upload-ecd-start").value) form.append("ecd_start", $("upload-ecd-start").value);
  if ($("upload-ecd-end").value) form.append("ecd_end", $("upload-ecd-end").value);
  try {
    setJob("QUEUED", t("msg.uploading", "上传结构并创建作业…"), 0.02);
    const job = await api("/api/jobs/upload", { method: "POST", body: form });
    state.jobId = job.id;
    state.demo = false;
    setJob(job.status, job.message, job.progress);
    watchJob(job.id);
  } catch (error) {
    setJob("FAILED", t("msg.uploadFail", "上传失败"), 0, error.message);
  }
}

async function submitAf3(event) {
  event.preventDefault();
  if (!$("af3-confirm").checked) return;
  const form = new FormData();
  form.append("structure", $("af3-structure").files[0]);
  form.append("confidence", $("af3-confidence").files[0]);
  try {
    const job = await api(`/api/jobs/${state.jobId}/af3`, { method: "POST", body: form });
    $("af3-dialog").close();
    setJob(job.status, job.message, job.progress);
    watchJob(state.jobId);
  } catch (error) {
    setJob("WAITING_FOR_AF3_UPLOAD", t("msg.af3Fail", "AF3 文件验证失败"), 0.34, error.message);
  }
}

/* ---------------- viewport (3D) ---------------- */

function residueQuery(residue) {
  if (residue.label_seq_id != null && residue.label_chain_id) {
    return { label_asym_id: residue.label_chain_id, label_seq_id: residue.label_seq_id };
  }
  const q = { auth_asym_id: residue.chain_id, auth_seq_id: residue.auth_seq_id };
  if (residue.insertion_code) q.pdbx_PDB_ins_code = residue.insertion_code;
  return q;
}

let molstarReady = false;

async function ensureViewers(structureUrl) {
  $("viewport-empty").classList.add("hidden");
  try {
    await MolstarView.init($("viewer-molstar"), structureUrl, state.representation);
    molstarReady = true;
    await paintMolstar();
    await MolstarView.focus(state.result.residues.map(residueQuery));
    await installMolstarTooltips();
  } catch (error) {
    console.error(error);
    $("viewport-empty").classList.remove("hidden");
    $("viewport-empty").querySelector("strong").textContent = t("hud.canvasFail", "结构画布载入失败");
    $("viewport-empty").querySelector("span").textContent = error.message;
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

function beadsColorFor(residue) {
  const design = selectedDesign();
  const designPos = new Set((design?.mutations || []).map((m) => m.position));
  if (state.customMutations.some((m) => m.position === residue.uniprot_pos)) return "#8a2f1d";
  if (designPos.has(residue.uniprot_pos)) return LAYER_META.design.color;
  if (state.layers.has("interface") && residue.interface) return LAYER_META.interface.color;
  if (state.layers.has("regions") && residue.region_id != null) {
    const region = state.result.regions.find((r) => r.id === residue.region_id);
    if (region) return region.color;
  }
  if (state.layers.has("folding") && residue.folding_screen_pass) return LAYER_META.folding.color;
  if (state.layers.has("surface") && residue.surface_accessible) return LAYER_META.surface.color;
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
  BeadsView.setColors(beadsColorFor, { selected, sidechain });
}

async function paintMolstar() {
  if (!molstarReady || !state.result) return;
  const data = [];
  const add = (records, color, extra = {}) => records.forEach((residue) => {
    data.push({ ...residueQuery(residue), color, ...extra });
  });
  if (state.layers.has("surface")) add(state.result.residues.filter((r) => r.surface_accessible), LAYER_META.surface.color);
  if (state.layers.has("folding")) add(state.result.residues.filter((r) => r.folding_screen_pass), LAYER_META.folding.color);
  if (state.layers.has("interface")) {
    add(state.result.residues.filter((r) => r.interface), LAYER_META.interface.color,
      { sideChain: true, representation: "ball-and-stick", representationColor: LAYER_META.interface.color });
  }
  if (state.layers.has("regions")) {
    for (const region of state.result.regions) {
      const positions = new Set(region.residues);
      add(state.result.residues.filter((r) => positions.has(r.uniprot_pos)), region.color);
    }
  }
  if (state.layers.has("design")) {
    const positions = new Set(
      (selectedDesign()?.mutations || []).map((m) => m.position)
        .concat(state.customMutations.map((m) => m.position))
    );
    add(state.result.residues.filter((r) => positions.has(r.uniprot_pos)), LAYER_META.design.color,
      { sideChain: true, representation: "ball-and-stick", representationColor: LAYER_META.design.color });
  }
  await MolstarView.color(data, "#7c909d");
  if (state.selectedResidue) await MolstarView.highlight(residueQuery(state.selectedResidue));
  renderLegend();
  if (state.engine === "beads") paintBeads();
}

async function setRepresentation(representation) {
  state.representation = representation;
  document.querySelectorAll("[data-representation]").forEach((b) => {
    b.setAttribute("aria-pressed", String(b.dataset.representation === representation));
  });
  if (state.engine === "molstar" && molstarReady) {
    await MolstarView.setRepresentation(representation);
    await paintMolstar();
  }
}

function setEngine(engine) {
  if (!state.result) return;
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
      BeadsView.init(beads, `/api/jobs/${state.jobId}/files/structure_for_energy.pdb`, state.result.residues)
        .then(() => {
          state.beadsReady = true;
          BeadsView.setCallbacks(
            (residue) => { if (residue) selectResidue(residue.uniprot_pos, false); },
            (residue) => showBeadsTip(residue)
          );
          paintBeads();
        })
        .catch((error) => {
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
  const design = selectedDesign();
  const designCount = (design?.mutations || []).length + state.customMutations.length;
  if (designCount) items.push(legendItem(LAYER_META.design.color, `${t("legend.selected", "已选突变")} (${designCount})`));
  if (state.layers.has("interface")) items.push(legendItem(LAYER_META.interface.color, interfaceLabel()));
  if (state.layers.has("regions")) {
    state.result.regions.forEach((region) => items.push(legendItem(region.color, region.name)));
  }
  items.push(legendItem("#7c909d", t("hud.backbone", "主链")));
  const collapsed = localStorage.getItem("ecd.legend") === "collapsed";
  $("legend").innerHTML =
    `<button type="button" class="legend-toggle">${collapsed ? t("hud.legendCollapsed", "图例 ▸") : t("hud.legendOpen", "图例 ▾")}</button>`
    + (collapsed ? "" : items.join(""));
  const toggle = $("legend").querySelector(".legend-toggle");
  toggle.addEventListener("click", () => {
    localStorage.setItem("ecd.legend", collapsed ? "open" : "collapsed");
    renderLegend();
  });
}

const legendItem = (color, text) =>
  `<span class="legend-item"><i class="legend-swatch" style="background:${h(color)}"></i>${h(text)}</span>`;

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
    return `<button type="button" class="${cls}" data-step="${n}"><span class="n">${n}</span>${h(t(key, key.split(".")[1]))}</button>`;
  }).join("");
  document.querySelectorAll("#steps [data-step]").forEach((b) => {
    b.addEventListener("click", () => setStep(Number(b.dataset.step)));
  });
}

function setStep(n) {
  state.step = n;
  renderSteps();
  document.querySelectorAll("#panel [data-step-content]").forEach((s) => {
    s.hidden = Number(s.dataset.stepContent) !== n;
  });
}

function renderPanelScaffold() {
  $("panel").innerHTML = `
  <section data-step-content="1">
    <h2>${h(t("step.surface", "表面可及性"))}</h2>
    <div class="note">${h(t("surface.note", "表面闸门：rSASA ≥ 0.25 且侧链 SASA ≥ 20 Å²。只有过闸残基进入设计候选。"))}</div>
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
        <button class="btn btn-mini" id="clear-custom" type="button">${h(t("custom.clear", "清空"))}</button>
      </div>
      <label class="fld" style="margin-bottom:8px"><span>${h(t("custom.name", "构建名称"))}</span>
        <input id="custom-name" value="Custom construct"></label>
      <div id="custom-mutations" class="empty">${h(t("custom.empty", "还没有选残基。"))}</div>
      <label class="chk"><input id="force-risk" type="checkbox">
        <span>${h(t("custom.forceRisk", "显式保留高风险替换（受保护位点仍不可强制）"))}</span></label>
      <button id="save-custom" class="btn primary" type="button" style="width:100%" disabled>${h(t("custom.save", "确认并保存此组合"))}</button>
      <p id="custom-feedback" class="note" style="margin:8px 0 0"></p>
    </div>
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
    renderCustomBuilder();
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
  renderLegend();
  renderSeqDock();
  if (repaint) {
    const url = state.result.structure_url || "/assets/demo/ecd_structure.cif";
    ensureViewers(url);
  }
}

function renderStep1() {
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
  $("methods-card").innerHTML = [
    `<strong>SASA</strong>　${h(state.result.methods.sasa_method)}`,
    `<strong>${h(t("s2.partition", "分区"))}</strong>　${h(state.result.partitioning.method)}`,
    `<strong>${h(t("s2.engineLabel", "稳定性引擎"))}</strong>　${h(stability.primary || t("s2.notRecorded", "未记录"))}；${h(stability.proxy_units || "")}`,
    `<strong>${h(t("s2.boundary", "结论边界"))}</strong>　${h(t("s2.boundaryBody", "计算预筛不等于实验确认；表位归因仍需表达、构象抗体、稳定性与功能读出。"))}`,
  ].join("<br>");
  const warnings = state.result.warnings?.length ? state.result.warnings : [t("s2.noWarnings", "没有运行级警告。")];
  $("warning-list").innerHTML = warnings.map((w) => `<li>${h(w)}</li>`).join("");
  $("criterion").textContent = criteriaText();
}

function criteriaText() {
  if (!state.result) return t("s2.criteriaEmpty", "载入结果后显示本步骤的定量标准。");
  const m = state.result.methods;
  return {
    surface: t("crit.surface", `rSASA ≥ ${m.rsasa_threshold} 且侧链 SASA ≥ ${m.sidechain_sasa_threshold_A2} Å²；SASA 引擎：${m.sasa_method}。`),
    folding: t("crit.folding", "在表面标准上排除二硫键 Cys、注释糖基化/PTM 与低置信度位点；表示低预测风险，不是实验确认。"),
    interface: t("crit.interface", "任一目标重原子距非目标链或异质残基 ≤5 Å；功能界面保留风险标记。"),
    regions: t("crit.regions", `${state.result.partitioning.method}；实际 ${state.result.regions.length} 补丁，覆盖 ${fmtPct(state.result.partitioning.coverage)}。`),
    design: t("crit.design", "FoldX 单点门与暴露度联动（+0.46/+0.60/+0.80）；组合门 +2.0 kcal/mol；探针内两两间距 ≥5.5 Å 且最大展开 ≤22 Å。"),
  }[state.layerFocus] || t("crit.default", "当前图层可与其他图层叠加显示。");
}

function renderStep3() {
  const designByRegion = new Map(state.result.designs.map((d) => [d.region_id, d]));
  $("region-list").innerHTML = state.result.regions.map((region) => {
    const design = designByRegion.get(region.id);
    const mutations = design?.mutations?.map(notation).join(" · ") || t("region.noAuto", "无自动组合");
    const tier = region.tier ? ` · ${region.tier}` : "";
    return `<div class="rowline click" data-focus-region="${region.id}">
      <span class="swatch" style="background:${h(region.color)}"></span>
      <span><b>${h(t("region.patch", "补丁"))} ${region.id} · ${h(region.name)}</b><br>
      <span class="muted">${h(region.structural_character)} · ${h(region.membrane_position || "")}${h(tier)} · ${h(mutations)}</span></span>
      <span class="grow"></span><span class="muted">${region.residues.length} ${h(t("region.sites", "位点"))}</span>
    </div><div class="note" style="margin:0 0 10px 20px">${h(region.rationale)}</div>`;
  }).join("") || `<div class="empty">${h(t("region.none", "没有可辩护的结构分区。"))}</div>`;
  document.querySelectorAll("[data-focus-region]").forEach((el) => {
    el.addEventListener("click", () => {
      state.selectedRegion = Number(el.dataset.focusRegion);
      state.layers.add("regions");
      repaintViews();
    });
  });
}

function renderStep4() {
  const chains = interfacePartnerChains();
  const contacts = state.result.residues.filter((r) => r.interface);
  $("iface-partners").innerHTML = chains.length
    ? chains.map((chain) => {
      const count = new Set(contacts.flatMap((r) => (r.interface_partners || [])
        .filter((p) => p.split(":")[0] === chain))).size;
      return `<div class="rowline"><span class="swatch" style="background:${LAYER_META.interface.color}"></span>
        <span><b>${h(t("iface.with", "与"))} ${h(chain)}</b><br><span class="muted">${count} ${h(t("iface.contacts", "个接触残基"))}</span></span></div>`;
    }).join("")
    : `<div class="empty">${h(t("iface.none", "当前结构无蛋白伙伴，界面叠加不可用。"))}</div>`;

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

function renderStep5() {
  $("design-list").innerHTML = state.result.designs.map((design) => {
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
    return `<div class="design-card">
      <div class="design-head">
        <div><strong>${h(design.id)} · ${h(design.name)}</strong><br><span class="st">${scoreLine}</span></div>
        <span class="st">${role} ${h(design.status)}</span>
      </div>
      <div>${chips}</div>
      <div class="actions-row">
        ${design.mutations?.length ? `<button type="button" class="btn btn-mini" data-view-design="${h(design.id)}">${h(t("design.review", "审查此构建体"))}</button>
        <button type="button" class="btn btn-mini" data-copy-design="${h(design.id)}">${h(t("design.copyCustom", "复制到自定义"))}</button>` : ""}
      </div>
    </div>`;
  }).join("") || `<div class="empty">${h(t("design.none", "没有推荐组合。"))}</div>`;
  document.querySelectorAll("[data-design-residue]").forEach((el) => {
    el.addEventListener("click", () => selectResidue(Number(el.dataset.designResidue), true));
  });
  document.querySelectorAll("[data-view-design]").forEach((el) => {
    el.addEventListener("click", () => selectDesign(el.dataset.viewDesign));
  });
  document.querySelectorAll("[data-copy-design]").forEach((el) => {
    el.addEventListener("click", () => copyDesignToCustom(el.dataset.copyDesign));
  });
  renderExports();
  renderDesignSelect();
}

function renderStep6() {
  renderCustomBuilder();
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
      <span class="main">${h(residue.domain)} · ${h(t("region.patch", "补丁"))} ${residue.region_id ?? "—"} · ${h(statusWord(residue))}${residue.interface ? ` · ${h(t("rd.interfaceShort", "界面"))}` : ""}</span>
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
    /* sideChain ball-and-stick so the picked residue stays visible even when
       no design layer covers it; yellow halo marks the exact site */
    await MolstarView.highlight({
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
    [...state.result.designs, ...previews].map((d) =>
      `<option value="${h(d.id)}">${h(d.id)} · ${h(d.mutations.map(notation).join("; ") || d.status)}</option>`
    ).join("");
  select.value = [...select.options].some((o) => o.value === state.selectedDesign) ? state.selectedDesign : "wt";
}

function selectDesign(id) {
  state.selectedDesign = id;
  renderDesignSelect();
  $("design-select").value = id;
  if (id !== "wt") state.layers.add("design");
  renderStep2();
  renderStep5();
  repaintViews();
  renderSeqDock();
}

/* ---------------- persistent sequence dock ---------------- */

function currentDesignSequences() {
  const design = selectedDesign();
  const wtAa = state.result.protein.ecd_sequence;
  const wtDna = state.result.protein.ecd_coding_sequence || [...wtAa].map((aa) => CODONS[aa]).join("");
  const map = {};
  if (state.customMutations.length) {
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
    designLabel: state.customMutations.length
      ? `${t("seqdock.custom", "自定义")} · ${state.customMutations.length} ${t("region.sites", "位点")}`
      : design ? `${design.id} · ${design.mutations.map(notation).join("; ") || design.status}` : t("seqdock.wtControl", "WT 对照"),
    hasRealCds: Boolean(design?.mutant_coding_sequence || state.result.protein.ecd_coding_sequence),
  };
}

function renderSeqDock() {
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
        return `<span class="${cls.join(" ")}" data-seq-residue="${pos}" title="${aa}${pos}${tip}">${aa}</span>`;
      }).join("");
    }
    const codons = sequence.match(/.{1,3}/g) || [];
    return codons.map((codon, i) => {
      const pos = start + i;
      const cls = ["tok"];
      if (isDesign && mutationPositions.has(pos)) cls.push("mut");
      if (selectedPosition === pos) cls.push("sel");
      return `<span class="${cls.join(" ")}" data-seq-residue="${pos}" title="${pos}: ${codon}" style="min-width:31px">${codon}</span>`;
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
  });
}

function renderSeqDockHead() { /* headers are static data-i18n elements now */ }

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
  $("save-custom").disabled = state.demo || !state.jobId;
  $("custom-feedback").textContent = state.demo ? t("custom.demoNote", "冻结示例可编辑预览；提交真实作业后可保存导出。") : "";
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
    const url = state.demo ? `/assets/demo/${file}` : `/api/jobs/${state.jobId}/files/${file}`;
    return `<a href="${h(url)}" download>${h(labels[file] || file)}</a>`;
  }).join("");
}

/* ---------------- repaint + events ---------------- */

function repaintViews() {
  if (state.engine === "molstar") paintMolstar();
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
  $("run-form").addEventListener("submit", submitRun);
  $("load-demo").addEventListener("click", loadDemo);
  $("open-upload").addEventListener("click", () => $("upload-dialog").showModal());
  $("close-upload").addEventListener("click", () => $("upload-dialog").close());
  $("upload-form").addEventListener("submit", submitUpload);
  $("close-af3").addEventListener("click", () => $("af3-dialog").close());
  $("af3-form").addEventListener("submit", submitAf3);
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
  loadDemo();

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
