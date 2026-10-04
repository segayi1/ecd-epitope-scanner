"use strict";

const element = (id) => document.getElementById(id);
let structureEffectViewer = null;

async function request(url, options) {
  const response = await fetch(url, options);
  const body = await response.json();
  if (!response.ok) {
    const detail = body.detail;
    throw new Error(typeof detail === "string" ? detail : detail?.message || detail?.code || `HTTP ${response.status}`);
  }
  return body.data;
}

function describeQcFindings(findings) {
  if (!findings?.length) return "无自动标记";
  return findings.map(item => {
    const details = [];
    if (item.channel) details.push(item.channel);
    if (item.gate_id) details.push(item.gate_id);
    const fraction = item.fraction_at_boundary ?? item.fraction_at_rail;
    if (Number.isFinite(fraction)) details.push(`${(100 * fraction).toFixed(2)}%`);
    if (Number.isFinite(item.count)) details.push(`n=${item.count}`);
    if (item.windows?.length) details.push(`时间窗 ${item.windows.map(window => `${window.event_start}-${window.event_end}`).join(",")}`);
    return `${item.code}${details.length ? ` (${details.join("; ")})` : ""}`;
  }).join("; ");
}

function appendWellEffectReadout(effectDocument, antibodies, root) {
  const wells = (effectDocument?.well_results || []).filter(item => Number.isFinite(Number(item.fraction)));
  if (!wells.length) return;
  const section = document.createElement("section");
  const title = document.createElement("h4"); title.textContent = "逐孔阳性比例";
  const note = document.createElement("p"); note.className = "muted";
  note.textContent = "每个点代表一个已审核样本孔，按批次和抗体分面；事件数不是独立重复数。技术重复汇总只纳入 PASS/PASS_WITH_OVERRIDE 孔，悬停可查看样本、计数和 QC。";
  section.append(title, note);
  const colors = {MUTANT: "#bd4d43", REFERENCE: "#347ca0", MUTANT_NC: "#d19332", REFERENCE_NC: "#657985"};
  const labels = {MUTANT: "突变抗原", REFERENCE: "WT 抗原", MUTANT_NC: "突变抗原 NC", REFERENCE_NC: "WT NC"};
  const groups = new Map();
  for (const item of wells) {
    const key = [item.probe_id, item.batch_id, item.antibody_construct_id].join("\u001f");
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(item);
  }
  for (const [key, rows] of [...groups.entries()].sort(([a], [b]) => a.localeCompare(b))) {
    const [probeId, batchId, antibodyId] = key.split("\u001f");
    const figure = document.createElement("figure");
    const caption = document.createElement("figcaption");
    caption.textContent = `${probeId} · ${batchId} · ${antibodies.find(item => item.id === antibodyId)?.name || antibodyId}`;
    const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    svg.setAttribute("viewBox", "0 0 720 320"); svg.setAttribute("role", "img");
    svg.setAttribute("aria-label", "逐孔阳性比例图");
    const add = (tag, attrs, text = "") => {
      const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
      for (const [name, value] of Object.entries(attrs)) node.setAttribute(name, String(value));
      if (text) node.textContent = text;
      svg.append(node); return node;
    };
    const conditions = ["MUTANT", "REFERENCE", "MUTANT_NC", "REFERENCE_NC"];
    const centers = [120, 285, 450, 615];
    for (const fraction of [0, .25, .5, .75, 1]) {
      const y = 245 - 200 * fraction;
      add("line", {x1: 60, x2: 680, y1: y, y2: y, stroke: "#d8e2e8"});
      add("text", {x: 52, y: y + 4, "text-anchor": "end", fill: "#526874", "font-size": 12}, `${Math.round(fraction * 100)}%`);
    }
    conditions.forEach((condition, index) => {
      add("text", {x: centers[index], y: 280, "text-anchor": "middle", fill: "#243b49", "font-size": 12}, labels[condition]);
      const points = rows.filter(item => item.condition === condition)
        .sort((a, b) => a.replicate_index - b.replicate_index || a.sample_id.localeCompare(b.sample_id));
      points.forEach((item, pointIndex) => {
        const x = centers[index] + (pointIndex - (points.length - 1) / 2) * Math.min(12, 72 / points.length);
        const point = add("circle", {
          cx: x, cy: 245 - 200 * Math.max(0, Math.min(1, Number(item.fraction))), r: 6,
          fill: colors[condition] || "#657985", opacity: ["HOLD", "EXCLUDE"].includes(item.qc_status) ? .35 : .9,
          stroke: item.qc_status === "PASS_WITH_OVERRIDE" ? "#172b39" : "white", "stroke-width": 1.5,
        });
        const tooltip = document.createElementNS("http://www.w3.org/2000/svg", "title");
        tooltip.textContent = `${item.sample_id} · ${item.positive_count}/${item.parent_count} · ${item.qc_status}${item.qc_reason ? ` · ${item.qc_reason}` : ""}`;
        point.append(tooltip);
      });
    });
    figure.append(caption, svg); section.append(figure);
  }
  const technical = new Map();
  for (const item of wells.filter(row => row.replicate_type === "TECHNICAL")) {
    const key = [item.probe_id, item.batch_id, item.antibody_construct_id, item.condition].join("\u001f");
    if (!technical.has(key)) technical.set(key, {accepted: [], excluded: 0});
    const summary = technical.get(key);
    if (["PASS", "PASS_WITH_OVERRIDE"].includes(item.qc_status)) summary.accepted.push(Number(item.fraction));
    else summary.excluded += 1;
  }
  const technicalTitle = document.createElement("h4"); technicalTitle.textContent = "技术重复孔离散度";
  section.append(technicalTitle);
  if (!technical.size) {
    const empty = document.createElement("p"); empty.className = "muted";
    empty.textContent = "本次没有技术重复孔；不计算孔内离散度。"; section.append(empty);
  } else {
    const table = document.createElement("table"), head = document.createElement("tr");
    for (const value of ["探针", "批次", "条件", "纳入技术孔", "未纳入技术孔", "比例中位数", "范围"]) {
      const cell = document.createElement("th"); cell.textContent = value; head.append(cell);
    }
    table.append(head);
    for (const [key, summary] of [...technical.entries()].sort(([a], [b]) => a.localeCompare(b))) {
      const values = summary.accepted.sort((a, b) => a - b);
      const [probeId, batchId, , condition] = key.split("\u001f");
      const middle = values.length >> 1;
      const median = values.length % 2 ? values[middle] : (values[middle - 1] + values[middle]) / 2;
      const row = document.createElement("tr");
      for (const value of [probeId, batchId, labels[condition] || condition, String(values.length), String(summary.excluded),
        values.length ? `${(100 * median).toFixed(2)}%` : "无通过孔",
        values.length ? `${(100 * values[0]).toFixed(2)}–${(100 * values.at(-1)).toFixed(2)}%` : "—"]) {
        const cell = document.createElement("td"); cell.textContent = value; row.append(cell);
      }
      table.append(row);
    }
    section.append(table);
  }
  root.append(section);
}

function appendFlowJoPercentageReplicateReadout(effect) {
  const observations = Array.isArray(effect?.replicate_observations)
    ? effect.replicate_observations : [];
  if (!observations.length) return null;

  const details = document.createElement("details");
  details.className = "flowjo-percentage-replicate-details";
  const summary = document.createElement("summary");
  summary.textContent = `逐孔分布与 QC（${observations.length} 个孔）`;
  const note = document.createElement("p");
  note.className = "muted";
  note.textContent = "每点对应一个结果表观测孔；点色表示该孔的人工 QC 决定。技术重复展示孔间差异，不代表独立生物学样本量。坐标只呈现输入百分比，不计算误差区间或显著性。";
  details.append(summary, note);

  const groupLabels = {MUTANT: "变体", REFERENCE: "WT", MUTANT_NC: "NC"};
  const qcColors = {PASS: "#0072B2", PASS_WITH_OVERRIDE: "#E69F00", HOLD: "#777777", EXCLUDE: "#B8B8B8"};
  const normalized = observations.map(point => {
    let group = point.comparison_group;
    if (!group) {
      const condition = String(point.condition_label || "").trim();
      if (condition === String(effect.probe_label || "")) group = "MUTANT";
      else if (condition === String(effect.reference_label || "")) group = "REFERENCE";
      else if (/\bNC\b/i.test(condition)) group = "MUTANT_NC";
    }
    const fraction = Number(point.fraction);
    return {...point, comparison_group: group, fraction};
  });
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("viewBox", "0 0 390 255");
  svg.setAttribute("role", "img");
  svg.setAttribute("aria-label", `${effect.panel_label} ${effect.probe_label} 逐孔百分比散点图`);
  svg.style.maxWidth = "440px";
  svg.style.height = "auto";
  const add = (tag, attrs, text = "") => {
    const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
    for (const [name, value] of Object.entries(attrs)) node.setAttribute(name, String(value));
    if (text) node.textContent = text;
    svg.append(node);
    return node;
  };
  const yTop = 32, yBottom = 207, centers = {MUTANT: 105, REFERENCE: 205, MUTANT_NC: 305};
  for (const tick of [0, 25, 50, 75, 100]) {
    const y = yBottom - (yBottom - yTop) * tick / 100;
    add("line", {x1: 52, x2: 360, y1: y, y2: y, stroke: "#d8e2e8"});
    add("text", {x: 45, y: y + 4, "text-anchor": "end", fill: "#526874", "font-size": 11}, `${tick}%`);
  }
  for (const [group, label] of Object.entries(groupLabels)) {
    add("text", {x: centers[group], y: 232, "text-anchor": "middle", fill: "#243b49", "font-size": 12}, label);
    const points = normalized.filter(point => point.comparison_group === group)
      .sort((a, b) => String(a.well_id).localeCompare(String(b.well_id))
        || Number(a.replicate_index || 0) - Number(b.replicate_index || 0));
    const jitter = Math.min(13, 100 / Math.max(1, points.length - 1));
    points.forEach((point, index) => {
      if (!Number.isFinite(point.fraction) || point.fraction < 0 || point.fraction > 1) return;
      const status = String(point.qc_status || "");
      const x = centers[group] + (index - (points.length - 1) / 2) * jitter;
      const y = yBottom - (yBottom - yTop) * point.fraction;
      const circle = add("circle", {
        cx: x, cy: y, r: 5, fill: qcColors[status] || "#657985",
        opacity: status === "HOLD" || status === "EXCLUDE" ? 0.55 : 0.95,
        stroke: "white", "stroke-width": 1.5,
      });
      const tooltip = document.createElementNS("http://www.w3.org/2000/svg", "title");
      tooltip.textContent = `${groupLabels[group]} · ${point.well_id || point.sample_id} · ${point.displayed_value || `${(100 * point.fraction).toFixed(2)}%`} · QC ${status || "未记录"}`;
      circle.append(tooltip);
    });
  }
  const legend = document.createElement("p");
  legend.className = "muted";
  legend.textContent = "颜色：蓝色 PASS、橙色 PASS_WITH_OVERRIDE、灰色 HOLD/EXCLUDE。每个分组按 0–100% 同一坐标范围显示。";
  details.append(svg, legend);

  const table = document.createElement("table");
  table.setAttribute("aria-label", "结果表逐孔观测与 QC 状态");
  const head = document.createElement("tr");
  for (const label of ["比较组", "孔", "结果表原值", "人工 QC", "纳入效应", "样本 ID"]) {
    const cell = document.createElement("th"); cell.textContent = label; head.append(cell);
  }
  table.append(head);
  for (const point of normalized) {
    const row = document.createElement("tr");
    const accepted = point.accepted_for_summary ?? ["PASS", "PASS_WITH_OVERRIDE"].includes(point.qc_status);
    const displayed = point.displayed_value ?? (Number.isFinite(point.percentage)
      ? `${point.percentage}%` : Number.isFinite(point.fraction) ? `${(100 * point.fraction).toFixed(2)}%` : "—");
    const values = [groupLabels[point.comparison_group] || "未分组", point.well_id || "—",
      displayed, point.qc_status || "未记录", accepted ? "是" : "否", point.sample_id || "—"];
    for (const value of values) { const cell = document.createElement("td"); cell.textContent = String(value); row.append(cell); }
    table.append(row);
  }
  details.append(table);
  return details;
}

async function showProjects() {
  const list = element("project-list");
  list.replaceChildren();
  try {
    const projects = await request("/api/v4/experiment-projects");
    if (!projects.length) { list.textContent = "尚无实验项目。"; return; }
    for (const project of projects) {
      const card = document.createElement("article");
      const title = document.createElement("button");
      title.type = "button";
      title.textContent = `${project.name} · ${project.project_id}`;
      title.addEventListener("click", () => showProject(project.project_id));
      const status = document.createElement("p");
      status.className = project.status === "COMPLETE" ? "muted" : "error";
      status.textContent = `状态：${project.status} · 来源：${project.source_job_id} · 修订：${project.revision}`;
      card.append(title, status);
      list.append(card);
    }
  } catch (error) { list.textContent = `读取失败：${error.message}`; }
}

async function showProject(projectId) {
  const detail = element("project-detail");
  const content = element("project-content");
  detail.hidden = false;
  content.textContent = "核验文件中…";
  try {
    const project = await request(`/api/v4/experiment-projects/${encodeURIComponent(projectId)}`);
    content.replaceChildren();
    const summary = document.createElement("p");
    summary.textContent = `${project.name} · ${project.status} · 创建时间：${project.created_at}`;
    content.append(summary);
    if (project.integrity_errors.length) {
      const errors = document.createElement("p");
      errors.className = "error";
      errors.textContent = project.integrity_errors.join("；");
      content.append(errors);
    }
    const hashes = document.createElement("dl");
    for (const [kind, value] of Object.entries(project.frozen_hashes || {})) {
      const key = document.createElement("dt"); key.textContent = kind;
      const hash = document.createElement("dd"); hash.className = "hash"; hash.textContent = value;
      hashes.append(key, hash);
    }
    content.append(hashes);
    const files = document.createElement("ul");
    for (const artifact of project.artifacts) {
      const item = document.createElement("li");
      const link = document.createElement("a");
      link.href = `/api/v4/experiment-projects/${encodeURIComponent(projectId)}/files/${encodeURIComponent(artifact.id)}`;
      link.textContent = `${artifact.kind} · ${artifact.size_bytes} bytes`;
      item.append(link); files.append(item);
    }
    content.append(files);
    const support = document.createElement("section");
    const supportTitle = document.createElement("h3"); supportTitle.textContent = "支持与审计";
    const includeOriginals = document.createElement("input");
    includeOriginals.type = "checkbox";
    includeOriginals.id = "support-include-originals";
    const includeLabel = document.createElement("label");
    includeLabel.htmlFor = includeOriginals.id;
    includeLabel.append(includeOriginals, document.createTextNode(" 将原始 FCS/WSP/FlowJo 结果表加入支持包"));
    const supportLink = document.createElement("a");
    supportLink.textContent = "下载项目支持包（含审计记录与校验文件）";
    supportLink.download = `experiment_${projectId}_support.zip`;
    const updateSupportUrl = () => {
      const query = includeOriginals.checked ? "?include_originals=true" : "";
      supportLink.href = `/api/v4/experiment-projects/${encodeURIComponent(projectId)}/support-bundle${query}`;
    };
    includeOriginals.addEventListener("change", updateSupportUrl);
    updateSupportUrl();
    support.append(supportTitle, includeLabel, supportLink);
    content.append(support);
    await renderProjectTools(projectId, project, content);
  } catch (error) { content.textContent = `读取失败：${error.message}`; }
}

async function renderProjectTools(projectId, project, root) {
  const base = `/api/v4/experiment-projects/${encodeURIComponent(projectId)}`;
  const [batches, samples, recommendations, constructs, antibodies, summaryMeasurements, gateRevisions, gateAnalyses, qcReviews, effectAnalyses, paraplumeRuns, flowjoPercentageImports] = await Promise.all([
    request(`${base}/batches`), request(`${base}/samples`),
    request(`${base}/recommended-antigens`), request(`${base}/antigen-constructs`),
    request(`${base}/antibody-constructs`), request(`${base}/summary-measurements`),
    request(`${base}/gates`), request(`${base}/gate-analyses`), request(`${base}/qc-reviews`),
    request(`${base}/effect-analyses`), request(`${base}/paraplume-runs`),
    request(`${base}/flowjo-percentage-imports`)
  ]);
  const data = document.createElement("section");
  const constructTitle = document.createElement("h3"); constructTitle.textContent = "抗原构建体"; data.append(constructTitle);
  const constructList = document.createElement("ul");
  for (const construct of constructs) {
    const li = document.createElement("li");
    const changes = (construct.mutations || []).map(item => item.position
      ? `${item.wt}${item.position}${item.mut}` : `${item.kind} ${item.position ?? "?"}`).join(", ") || "序列已登记";
    const riskFlags = (construct.risk_flags || []).join(", ") || "无风险标记";
    li.textContent = `${construct.name} · ${construct.source} · ${construct.mapping_status} · ${changes} · ${riskFlags} · ${construct.sequence_sha256.slice(0, 12)}`;
    constructList.append(li);
  }
  if (!constructs.length) {
    const empty = document.createElement("li"); empty.textContent = "尚无抗原构建体。"; constructList.append(empty);
  }
  data.append(constructList);
  const customForm = document.createElement("form");
  customForm.innerHTML = '<label>自定义名称 <input name="name" required maxlength="120"></label><label>UniProt 位置 <input name="position" type="number" min="1" required></label><label>WT <input name="wt" required maxlength="1" pattern="[ACDEFGHIKLMNPQRSTVWY]"></label><label>突变 <input name="mut" required maxlength="1" pattern="[ACDEFGHIKLMNPQRSTVWY]"></label><button>登记单点构建体</button>';
  customForm.addEventListener("submit", async (event) => {
    event.preventDefault(); const fields = new FormData(customForm);
    const button = customForm.querySelector("button"); button.disabled = true;
    try {
      await request(`${base}/antigen-constructs`, {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({name: fields.get("name"), source: "CUSTOM",
          mutations: [{position: Number(fields.get("position")), wt: fields.get("wt"), mut: fields.get("mut")}],
          expected_revision: project.revision})});
      await showProject(projectId);
    } catch (error) { alert(`登记构建体失败：${error.message}`); button.disabled = false; }
  });
  data.append(customForm);
  const fullSequenceForm = document.createElement("form");
  fullSequenceForm.innerHTML = '<label>完整序列构建体名称 <input name="name" required maxlength="120"></label><label>完整 ECD 序列 <textarea name="sequence" required></textarea></label><button name="preview" type="button">预览与 WT 差异</button><button name="save" disabled>确认并登记完整序列</button>';
  const fullSequencePreview = document.createElement("pre"); fullSequencePreview.textContent = "请先预览差异；预览不会创建 revision。";
  let previewedFullSequence = null;
  const normalizeFullSequence = () => fullSequenceForm.elements.sequence.value.replace(/\s+/g, "").toUpperCase();
  fullSequenceForm.elements.sequence.addEventListener("input", () => {
    previewedFullSequence = null;
    fullSequenceForm.elements.save.disabled = true;
    fullSequencePreview.textContent = "序列已改变，请重新预览。";
  });
  fullSequenceForm.elements.preview.addEventListener("click", async () => {
    const button = fullSequenceForm.elements.preview; button.disabled = true;
    try {
      const sequence = normalizeFullSequence();
      const preview = await request(`${base}/antigen-constructs/preview`, {method: "POST",
        headers: {"Content-Type": "application/json"}, body: JSON.stringify({full_sequence: sequence})});
      const changes = (preview.mutations || []).map(item => item.kind === "SUBSTITUTION"
        ? `${item.wt}${item.position}${item.mut}${item.mapped ? "" : " (未映射)"}`
        : `${item.kind} ${item.position ?? item.after_position ?? "?"}${item.mapped ? "" : " (未映射)"}`);
      fullSequencePreview.textContent = `差异 ${changes.length} 项：${changes.join("；") || "无"}\n映射 ${preview.mapping_status}；风险 ${(preview.risk_flags || []).join("、") || "无"}\n序列 SHA-256 ${preview.sequence_sha256}`;
      previewedFullSequence = sequence;
      fullSequenceForm.elements.save.disabled = false;
    } catch (error) {
      fullSequencePreview.textContent = `预览失败：${error.message}`;
    } finally { button.disabled = false; }
  });
  fullSequenceForm.addEventListener("submit", async (event) => {
    event.preventDefault();
    const sequence = normalizeFullSequence();
    if (!previewedFullSequence || sequence !== previewedFullSequence) {
      fullSequencePreview.textContent = "登记前必须对当前序列重新预览。";
      return;
    }
    const button = fullSequenceForm.elements.save; button.disabled = true;
    try {
      await request(`${base}/antigen-constructs`, {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({name: fullSequenceForm.elements.name.value, source: "CUSTOM",
          full_sequence: sequence, expected_revision: project.revision})});
      await showProject(projectId);
    } catch (error) { alert(`登记完整序列失败：${error.message}`); button.disabled = false; }
  });
  data.append(fullSequenceForm, fullSequencePreview);
  const wtButton = document.createElement("button"); wtButton.type = "button";
  const hasWtConstruct = constructs.some(item => item.sequence === project.manifest.reference_sequence);
  wtButton.textContent = hasWtConstruct ? "已登记 WT 参考构建体" : "登记冻结参考序列为 WT 构建体";
  wtButton.disabled = hasWtConstruct;
  wtButton.addEventListener("click", async () => {
    wtButton.disabled = true;
    try {
      await request(`${base}/antigen-constructs`, {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({name: "WT reference", source: "CUSTOM",
          full_sequence: project.manifest.reference_sequence, expected_revision: project.revision})});
      await showProject(projectId);
    } catch (error) { alert(`登记 WT 参考构建体失败：${error.message}`); wtButton.disabled = false; }
  });
  data.append(wtButton);
  const antibodyTitle = document.createElement("h3"); antibodyTitle.textContent = "抗体构建体"; data.append(antibodyTitle);
  const paraplumeLicenseNote = document.createElement("p");
  paraplumeLicenseNote.className = "muted";
  paraplumeLicenseNote.textContent = "当前用途：仅本机非商业研究，不分发（由用户确认）；尚无与 MIT 项目组合的书面许可。Paraplume 自定义许可含不与不兼容许可软件组合的限制，且未说明仅限分发场景，因此本机集成调用的许可适用范围仍待澄清。商业产品/服务集成另需事先书面授权。";
  data.append(paraplumeLicenseNote);
  const antibodyList = document.createElement("ul");
  const predictionsByConstruct = new Map(paraplumeRuns.map(item => [item.antibody_construct_id, item]));
  for (const antibody of antibodies) {
    const li = document.createElement("li");
    li.textContent = `${antibody.name} · ${antibody.role} · ${antibody.mutation_class} · VH ${antibody.heavy_mutations.length} changes · VL ${antibody.light_mutations.length} changes`;
    const prediction = predictionsByConstruct.get(antibody.id);
    if (prediction) {
      const device = prediction.provenance?.actual_device || "device unknown";
      const track = (chain, scores) => [...scores]
        .map((probability, index) => ({position: index + 1, probability}))
        .sort((a, b) => b.probability - a.probability)
        .slice(0, 5)
        .map(item => `${chain}${item.position} ${item.probability.toFixed(3)}`).join(" · ");
      const result = document.createElement("p");
      result.className = "muted";
      result.textContent = `PREDICTED_UNREVIEWED · Paraplume 2.0.0 · ${device} · 抗体侧预测，不是抗原表位证据。VH: ${track("H", prediction.vh_probabilities)}；VL: ${track("L", prediction.vl_probabilities)}`;
      li.append(result);
    }
    const predictButton = document.createElement("button");
    predictButton.type = "button";
    predictButton.textContent = prediction ? "重新运行完整 Paraplume" : "运行完整 Paraplume";
    predictButton.addEventListener("click", async () => {
      predictButton.disabled = true;
      predictButton.textContent = "正在排队…";
      try {
        const queued = await request(`${base}/paraplume`, {method: "POST", headers: {"Content-Type": "application/json"},
          body: JSON.stringify({antibody_construct_ids: [antibody.id], model_version: "2.0.0", model_size: "large", expected_revision: project.revision})});
        let task;
        const terminal = new Set(["COMPLETE", "FAILED", "INTERRUPTED"]);
        const progressMessages = new Map([
          ["Checking full-model cache and preparing inference", "正在核验完整模型缓存并准备推理"],
          ["Preparing the full large model", "正在准备完整 large 模型"],
          ["Loading the full Paraplume model and preparing five embeddings", "正在加载完整 Paraplume 模型并准备 5 种嵌入"],
          ["Loading the full Paraplume model", "正在加载完整 Paraplume 模型"],
          ["Computing IgBert residue embeddings", "正在计算 IgBert 残基嵌入"],
          ["Completed IgBert residue embeddings", "IgBert 残基嵌入完成"],
          ["Computing IgT5 residue embeddings", "正在计算 IgT5 残基嵌入"],
          ["Completed IgT5 residue embeddings", "IgT5 残基嵌入完成"],
          ["Computing ESM2 residue embeddings", "正在计算 ESM2 残基嵌入"],
          ["Completed ESM2 residue embeddings", "ESM2 残基嵌入完成"],
          ["Computing AntiBERTy residue embeddings", "正在计算 AntiBERTy 残基嵌入"],
          ["Completed AntiBERTy residue embeddings", "AntiBERTy 残基嵌入完成"],
          ["Computing ProtT5 residue embeddings", "正在计算 ProtT5 残基嵌入"],
          ["Completed ProtT5 residue embeddings", "ProtT5 残基嵌入完成"],
          ["GPU memory was insufficient; clearing allocations and retrying the full model on CPU", "GPU 显存不足，正在清理并用 CPU 重跑完整模型"],
          ["Validating residue probabilities and model provenance", "正在核对残基概率和模型来源"],
          ["Saving validated residue results and provenance", "正在保存已核验的残基结果和来源信息"],
        ]);
        const showTaskProgress = update => {
          const percent = Math.round(100 * Math.max(0, Math.min(1, Number(update.progress) || 0)));
          const message = update.message || update.output?.message || "正在执行完整模型";
          const localizedMessage = progressMessages.get(message) || message;
          predictButton.textContent = `${percent}% · ${localizedMessage}`;
        };
        if ("EventSource" in window && queued.events_url) {
          task = await new Promise(resolve => {
            const stream = new EventSource(queued.events_url);
            let settled = false;
            const finish = value => {
              if (settled) return;
              settled = true; stream.close(); resolve(value);
            };
            stream.onmessage = event => {
              try {
                const update = JSON.parse(event.data);
                showTaskProgress(update);
                if (terminal.has(update.status)) finish(update);
              } catch { finish(null); }
            };
            stream.onerror = () => finish(null);
          });
        }
        if (!task) {
          do {
            await new Promise(resolve => setTimeout(resolve, 2000));
            task = await request(`/api/v4/experiment-tasks/${encodeURIComponent(queued.task_id)}`);
            showTaskProgress(task);
          } while (!terminal.has(task.status));
        }
        if (task.status !== "COMPLETE") throw new Error(task.output?.message || task.error_code || task.status);
        await showProject(projectId);
      } catch (error) {
        alert(`Paraplume 任务失败：${error.message}`);
        predictButton.disabled = false;
        predictButton.textContent = prediction ? "重新运行完整 Paraplume" : "运行完整 Paraplume";
      }
    });
    li.append(predictButton);
    antibodyList.append(li);
  }
  if (!antibodies.length) { const li = document.createElement("li"); li.textContent = "尚无抗体构建体"; antibodyList.append(li); }
  data.append(antibodyList);
  const antibodyForm = document.createElement("form");
  antibodyForm.innerHTML = '<label>构建体名称 <input name="name" required maxlength="120"></label><label>类型 <select name="role"><option value="PARENT">亲本</option><option value="VARIANT">变体</option></select></label><label>亲本构建体 <select name="parent"><option value="">选择亲本（变体必填）</option></select></label><label>VH 序列 <textarea name="vh" required></textarea></label><label>VL 序列 <textarea name="vl" required></textarea></label><button>登记抗体构建体</button>';
  const parentSelect = antibodyForm.elements.parent;
  for (const antibody of antibodies.filter(item => item.role === "PARENT")) {
    const option = document.createElement("option"); option.value = antibody.id; option.textContent = antibody.name;
    parentSelect.append(option);
  }
  const updateAntibodyRole = () => {
    const variant = antibodyForm.elements.role.value === "VARIANT";
    parentSelect.required = variant;
    parentSelect.closest("label").hidden = !variant;
  };
  antibodyForm.elements.role.addEventListener("change", updateAntibodyRole); updateAntibodyRole();
  antibodyForm.addEventListener("submit", async (event) => {
    event.preventDefault(); const fields = new FormData(antibodyForm);
    const button = antibodyForm.querySelector("button"); button.disabled = true;
    try {
      await request(`${base}/antibody-constructs`, {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({name: fields.get("name"), role: fields.get("role"),
          parent_construct_id: fields.get("parent") || null,
          vh_sequence: fields.get("vh"), vl_sequence: fields.get("vl"), expected_revision: project.revision})});
      await showProject(projectId);
    } catch (error) { alert(`登记抗体构建体失败：${error.message}`); button.disabled = false; }
  });
  data.append(antibodyForm);
  const recommendedTitle = document.createElement("h4"); recommendedTitle.textContent = "V4 推荐设计（来源于冻结设计结果）"; data.append(recommendedTitle);
  const recommendationList = document.createElement("ul");
  for (const recommendation of recommendations) {
    const li = document.createElement("li");
    const label = document.createElement("span");
    label.textContent = `${recommendation.name || recommendation.id} · ${recommendation.status} · ${(recommendation.mutations || []).map(item => `${item.wt}${item.position}${item.mutant || item.mut}`).join(", ")}`;
    li.append(label);
    if (String(recommendation.status || "").startsWith("RECOMMENDED_")) {
      const add = document.createElement("button"); add.type = "button"; add.textContent = "加入构建列表";
      add.addEventListener("click", async () => {
        add.disabled = true;
        try {
          await request(`${base}/antigen-constructs`, {method: "POST", headers: {"Content-Type": "application/json"},
            body: JSON.stringify({name: recommendation.name || recommendation.id, source: "RECOMMENDED",
              recommended_construct_id: recommendation.id, expected_revision: project.revision})});
          await showProject(projectId);
        } catch (error) { alert(`加入推荐构建体失败：${error.message}`); add.disabled = false; }
      });
      li.append(" ", add);
    }
    recommendationList.append(li);
  }
  if (!recommendations.length) {
    const empty = document.createElement("li");
    empty.textContent = "此项目的旧版冻结快照未包含推荐设计；自定义构建体仍可登记。";
    recommendationList.append(empty);
  }
  data.append(recommendationList);
  const batchTitle = document.createElement("h3"); batchTitle.textContent = "实验批次"; data.append(batchTitle);
  const batchWarningLabels = {
    ACTIVATION_SAMPLE_MISSING: "缺少 activation 样本",
    NC_CONTROL_MISSING: "缺少同批 NC",
    WT_CONSTRUCT_MISSING: "尚未登记 WT 参考构建体",
    WT_REFERENCE_SAMPLE_MISSING: "缺少同批 WT 参考样本",
    CONSTRUCT_LINK_MISSING: "有样本未关联抗原或抗体构建体",
  };
  for (const batch of batches) {
    const line = document.createElement("p");
    const dose = batch.dose == null ? "未记录剂量" : `${batch.dose} ${batch.dose_unit || "(单位未填)"}`;
    const independence = batch.independent_batch ? "独立批次" : "非独立/技术批次";
    line.textContent = `${batch.name} · ${batch.date} · ${batch.assay} · ${batch.cell_line} · ${dose} · ${independence}`;
    data.append(line);
    if ((batch.design_warnings || []).length) {
      const warning = document.createElement("p"); warning.className = "error";
      warning.textContent = `批次设计未完整：${batch.design_warnings.map(code => batchWarningLabels[code] || code).join("；")}`;
      data.append(warning);
    }
  }
  const batchForm = document.createElement("form");
  batchForm.innerHTML = '<label>批次名称 <input name="name" required></label><label>日期 <input name="date" type="date" required></label><label>实验类型 <select name="assay"><option value="ANTIGEN_MYC">抗原 Myc</option><option value="FAB_MYC">Fab Myc</option><option value="ACTIVATION">激活</option><option value="OTHER">其他</option></select></label><label>细胞系 <input name="cell_line" value="unspecified" required></label><label>剂量 <input name="dose" type="number" min="0" step="any"></label><label>剂量单位 <input name="dose_unit" placeholder="例如 ug/mL"></label><label>阳性门 ID <input name="positive_gate_id"></label><label><input name="independent_batch" type="checkbox" checked> 独立实验批次</label><label>批次备注 <input name="notes" maxlength="1000"></label><button>新增批次</button>';
  batchForm.addEventListener("submit", async (event) => {
    event.preventDefault(); const fields = new FormData(batchForm);
    try {
      await request(`${base}/batches`, {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({name: fields.get("name"), date: fields.get("date"), assay: fields.get("assay"),
          cell_line: fields.get("cell_line"), dose: fields.get("dose") === "" ? null : Number(fields.get("dose")),
          dose_unit: fields.get("dose_unit") || null, positive_gate_id: fields.get("positive_gate_id") || null,
          independent_batch: fields.get("independent_batch") === "on",
          settings: fields.get("notes") ? {notes: fields.get("notes")} : {}, expected_revision: project.revision})});
      await showProject(projectId);
    } catch (error) { alert(`新增批次失败：${error.message}`); }
  });
  data.append(batchForm);

  const summaryTitle = document.createElement("h4"); summaryTitle.textContent = "历史 Myc 汇总（summary-only）"; data.append(summaryTitle);
  const summaryNotice = document.createElement("p"); summaryNotice.className = "muted";
  summaryNotice.textContent = "仅登记历史汇总值及来源，不创建 FCS 样本、不标记为同期表达对照，也不会进入 activation 效应计算。";
  data.append(summaryNotice);
  const summaryList = document.createElement("ul");
  const antigenNames = new Map(constructs.map(item => [item.id, item.name]));
  const antibodyNames = new Map(antibodies.map(item => [item.id, item.name]));
  for (const measurement of summaryMeasurements) {
    const li = document.createElement("li");
    const name = measurement.assay === "ANTIGEN_MYC"
      ? antigenNames.get(measurement.construct_id) : antibodyNames.get(measurement.construct_id);
    li.textContent = `HISTORICAL_SUMMARY_ONLY · ${measurement.measured_on} · ${measurement.assay} · ${name || measurement.construct_id} · ${measurement.metric} ${measurement.value} ${measurement.unit || ""} · 来源 ${measurement.source_label} · 审核 ${measurement.reviewer}`;
    summaryList.append(li);
  }
  if (!summaryMeasurements.length) {
    const empty = document.createElement("li"); empty.textContent = "尚无历史 Myc 汇总记录。"; summaryList.append(empty);
  }
  data.append(summaryList);
  const summaryForm = document.createElement("form");
  summaryForm.innerHTML = '<label>类型 <select name="assay"><option value="ANTIGEN_MYC">抗原 Myc</option><option value="FAB_MYC">Fab Myc</option></select></label><label>构建体 <select name="construct" required></select></label><label>测量日期 <input name="measured_on" type="date" required></label><label>指标 <select name="metric"><option value="POSITIVE_FRACTION">阳性比例 (0–1)</option><option value="MEDIAN_FLUORESCENCE">中位荧光</option><option value="GEOMETRIC_MEAN">几何均值</option><option value="OTHER">其他</option></select></label><label>数值 <input name="value" type="number" min="0" step="any" required></label><label>单位 <input name="unit" maxlength="80"></label><label>来源标签 <input name="source_label" maxlength="240" required placeholder="记录本、文件或报告编号"></label><label>审核人 <input name="reviewer" maxlength="120" required></label><label>备注 <input name="notes" maxlength="1000"></label><button>登记历史汇总</button>';
  const updateSummaryConstructs = () => {
    const choices = summaryForm.elements.assay.value === "ANTIGEN_MYC" ? constructs : antibodies;
    summaryForm.elements.construct.replaceChildren();
    for (const construct of choices) {
      const option = document.createElement("option"); option.value = construct.id; option.textContent = construct.name;
      summaryForm.elements.construct.append(option);
    }
    summaryForm.elements.construct.disabled = !choices.length;
    summaryForm.querySelector("button").disabled = !choices.length;
  };
  summaryForm.elements.assay.addEventListener("change", updateSummaryConstructs); updateSummaryConstructs();
  summaryForm.addEventListener("submit", async (event) => {
    event.preventDefault(); const fields = new FormData(summaryForm);
    const button = summaryForm.querySelector("button"); button.disabled = true;
    try {
      await request(`${base}/summary-measurements`, {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({assay: fields.get("assay"), construct_id: fields.get("construct"),
          measured_on: fields.get("measured_on"), metric: fields.get("metric"), value: Number(fields.get("value")),
          unit: fields.get("unit"), source_label: fields.get("source_label"), reviewer: fields.get("reviewer"),
          notes: fields.get("notes"), expected_revision: project.revision})});
      await showProject(projectId);
    } catch (error) { alert(`历史 Myc 汇总未保存：${error.message}`); button.disabled = false; }
  });
  data.append(summaryForm);

  const samplesTitle = document.createElement("h3"); samplesTitle.textContent = `已确认样本（${samples.length}）`; data.append(samplesTitle);
  const sampleList = document.createElement("ul");
  for (const sample of samples) {
    const li = document.createElement("li");
    li.textContent = `${sample.condition?.sample_label || sample.id} · ${sample.role} · ${sample.replicate_type} ${sample.replicate_index}`;
    sampleList.append(li);
  }
  data.append(sampleList);

  const flowJoPreviews = project.artifacts.filter(item => item.kind === "flowjo_import_preview")
    .sort((left, right) => String(right.created_at).localeCompare(String(left.created_at)));
  let unresolvedFlowJoRebuild = null;
  if (flowJoPreviews.length) {
    const previewTitle = document.createElement("h3"); previewTitle.textContent = "FlowJo 完整事件门计数"; data.append(previewTitle);
    const disclaimer = document.createElement("p");
    disclaimer.className = "muted";
    const mycOnly = samples.length > 0 && samples.every(sample =>
      batches.find(batch => batch.batch_id === sample.batch_id)?.assay === "ANTIGEN_MYC");
    disclaimer.textContent = mycOnly
      ? "这是导入 workspace 时按完整事件重算的描述性门计数。本项目为 Myc 染色，不能据此推导激活效应或 WT 比较。"
      : "这是导入 workspace 时按完整事件重算的描述性门计数；计数本身不等于功能效应，需结合批次设计、对照和 QC 审核解释。";
    data.append(disclaimer);
    try {
      const response = await fetch(`${base}/files/${encodeURIComponent(flowJoPreviews[0].id)}`);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const preview = await response.json();
      const records = preview.records || {};
      const flowjoWaiting = records.status === "WAITING_FOR_USER" || records.unresolved_nodes?.length > 0;
      const cachedCountMismatch = (records.gate_count_verification?.mismatched_gate_count || 0) > 0;
      const approximateEllipsoid = (records.samples || []).some(sample =>
        (sample.gates || []).some(gate =>
          gate.mapping_note === "FLOWKIT_APPROXIMATES_FLOWJO_ELLIPSOID_AS_POLYGON"));
      if (approximateEllipsoid) {
        const ellipseNote = document.createElement("p"); ellipseNote.className = "muted";
        ellipseNote.textContent = "FlowJo 椭圆门由 FlowKit 转换为近似多边形后分类；计数采用该近似边界。";
        data.append(ellipseNote);
      }
      if (!flowjoWaiting && cachedCountMismatch) {
        const countNote = document.createElement("p"); countNote.className = "muted";
        countNote.textContent = `FlowKit 在同一批次按一致的门层级、门型、通道、门内/门外与 Boolean 逻辑、变换和补偿设置重算；与 FlowJo 保存的百分比不同属于方法差异提示，不会阻止计数进入 QC。${records.gate_count_verification.mismatched_gate_count}/${records.gate_count_verification.compared_gate_count} 条门计数与缓存值不同。门的具体边界不要求跨样本完全相同。`;
        data.append(countNote);
      }
      if (flowjoWaiting) {
        const unresolvedNodes = records.unresolved_nodes || [];
        const inconsistentMethod = unresolvedNodes.some(node =>
          node.reason_code === "FLOWJO_GATE_METHOD_INCONSISTENT");
        const unmatchedSample = unresolvedNodes.some(node =>
          node.reason_code === "WSP_SAMPLE_UNMATCHED");
        const transformedEllipsoid = unresolvedNodes.some(node =>
          node.reason_code === "FLOWJO_ELLIPSOID_TRANSFORM_UNSUPPORTED");
        const unsupportedTransform = unresolvedNodes.some(node =>
          node.reason_code?.startsWith("FLOWJO_TRANSFORM_"));
        const unsupportedGate = unresolvedNodes.some(node =>
          node.reason_code?.startsWith("FLOWJO_GATE_"));
        const unresolved = document.createElement("p"); unresolved.className = "error";
        unresolved.textContent = inconsistentMethod
          ? "同一批次的门层级、门型、通道、门内/门外选择、Boolean 组合、变换或补偿设置不一致；本次门计数不会进入 QC 或效应分析。请检查分批方式，或使用下方内部门控完整重建并留下审核记录。"
          : unmatchedSample
            ? "部分 FCS 无法唯一对应到 FlowJo workspace 样本；本次门计数不会进入 QC 或效应分析。请检查文件分组或使用下方内部门控完整重建。"
          : transformedEllipsoid
            ? "FlowKit 1.3.2 不支持带变换引用的 FlowJo 椭圆门，已停止重算并标记待处理；可导入已圈门结果表，或人工确认完整门树后再重建。"
          : unsupportedTransform
            ? "FlowJo workspace 含 FlowKit 1.3.2 不支持或参数异常的通道变换，已停止重算并列出待处理项；可导入已圈门结果表，或在 FlowJo 调整变换后重新导出。"
          : unsupportedGate
            ? "FlowKit 1.3.2 无法忠实重放这些 FlowJo 门型或门结构，已停止重算；可直接导入已圈门结果表，或在 FlowJo 改用受支持的门后重新导出。"
          : "该 FlowJo workspace 含有未解析门；计数不会用于 QC 或效应分析。可用下方内部门控完整重建，并留下审核记录。";
        data.append(unresolved);
        const nodeList = document.createElement("ul");
        for (const node of unresolvedNodes) {
          const item = document.createElement("li");
          const parent = (node.parent_path || []).join(" / ") || "root";
          item.textContent = node.reason_code === "FLOWJO_GATE_METHOD_INCONSISTENT"
            ? `${node.name || "门控方法"} · ${node.method_variants} 种方法变体 · ${node.samples_checked} 个样本`
            : node.reason_code === "WSP_SAMPLE_UNMATCHED"
              ? `${node.sample_id || "样本未知"} · 无法唯一对应 workspace 样本`
            : node.reason_code?.startsWith("FLOWJO_TRANSFORM_")
              ? `${node.sample_id || "样本未知"} · 变换 ${node.transform_type || "未知"} · 通道 ${(node.channels || []).join(", ") || "未知"} · ${node.reason_code}`
            : `${node.sample_id || "样本未知"} · ${node.name || "未命名门"} · ${node.source_gate_type || "未知类型"} · 通道 ${(node.channels || []).join(", ") || "未知"} · ${node.reason_code || "未解析"} · 父路径 ${parent}`;
          nodeList.append(item);
        }
        data.append(nodeList);
        const sourceFlowjoImportId = flowJoPreviews[0].relative_path.split("/")[2];
        const linkedRebuild = gateRevisions.filter(revision =>
          revision.flowjo_rebuild?.flowjo_import_id === sourceFlowjoImportId).at(-1);
        const completedRebuild = linkedRebuild && gateAnalyses.find(analysis =>
          analysis.gate_revision_id === linkedRebuild.id && analysis.status === "COMPLETE");
        if (completedRebuild) {
          const rebuilt = document.createElement("p"); rebuilt.className = "muted";
          rebuilt.textContent = `已由内部门控修订 v${linkedRebuild.version} 完整重建并重算；审核人 ${linkedRebuild.flowjo_rebuild.reviewer}。原始 FlowJo 未解析记录仍保留。`;
          data.append(rebuilt);
        } else if (Array.isArray(preview.flowjo_input_artifact_ids)
          && preview.flowjo_input_artifact_ids.length) {
          unresolvedFlowJoRebuild = {
            flowjo_import_id: sourceFlowjoImportId,
            input_artifact_ids: preview.flowjo_input_artifact_ids,
            nodes: unresolvedNodes,
          };
          const instructions = document.createElement("p"); instructions.className = "muted";
          instructions.textContent = "重建时须在下方填写完整门树（包含未解析节点的等价门），确认覆盖此 WSP 的层级逻辑，并填写审核人和理由。分析将限定使用该 WSP 绑定的 FCS 文件。";
          data.append(instructions);
        } else {
          const reimport = document.createElement("p"); reimport.className = "error";
          reimport.textContent = "此旧版预览未记录其 FCS artifact 清单，无法安全绑定重建结果。请重新导入 WSP 后再处理。";
          data.append(reimport);
        }
      }
      const sampleByLabel = new Map();
      for (const sample of samples) {
        const well = sample.condition?.well_id;
        if (well) sampleByLabel.set(String(well).toUpperCase(), sample);
        if (sample.condition?.sample_label) sampleByLabel.set(sample.condition.sample_label, sample);
      }
      const table = document.createElement("table");
      const head = document.createElement("tr");
      for (const value of ["样本", "角色", "门", "父门事件", "门内事件", "占父门比例"]) {
        const cell = document.createElement("th"); cell.textContent = value; head.append(cell);
      }
      table.append(head); data.append(table);
      const bySample = new Map();
      for (const gate of records.gate_counts || []) {
        if (!bySample.has(gate.sample_id)) bySample.set(gate.sample_id, []);
        bySample.get(gate.sample_id).push(gate);
      }
      for (const [sampleId, gates] of bySample) {
        const sampleLabel = sampleId.replace(/\.fcs$/i, "");
        const sample = sampleByLabel.get(sampleLabel.toUpperCase()) || sampleByLabel.get(sampleLabel);
        const parentCounts = new Map(gates.map(gate => [gate.gate_name, gate.count]));
        for (const gate of gates) {
          const row = document.createElement("tr");
          const parentCount = gate.parent_gate === "root" ? null : parentCounts.get(gate.parent_gate);
          const cells = [sampleLabel, sample?.role || "未匹配", gate.gate_name,
            parentCount == null ? "全部事件" : String(parentCount), String(gate.count),
            parentCount == null ? "—" : `${(100 * gate.count / parentCount).toFixed(2)}%`];
          for (const value of cells) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
          table.append(row);
        }
      }
      data.append(table);
      const qcTitle = document.createElement("h4"); qcTitle.textContent = "导入样本自动 QC 提示"; data.append(qcTitle);
      const qcTable = document.createElement("table");
      const qcHead = document.createElement("tr");
      for (const value of ["样本", "角色", "QC 状态", "标记", "时间窗异常", "检查事件"]) {
        const cell = document.createElement("th"); cell.textContent = value; qcHead.append(cell);
      }
      qcTable.append(qcHead); data.append(qcTable);
      const flowjoImportId = flowJoPreviews[0].relative_path.split("/")[2];
      const flowjoQcReviews = qcReviews.filter(item =>
        item.source_type === "FLOWJO_IMPORT" && item.source_id === flowjoImportId);
      for (const sampleRecord of records.samples || []) {
        const qc = sampleRecord.qc || {
          status: flowjoWaiting ? "WAITING_FOR_USER" : "QC_NOT_AVAILABLE", findings: [], time_windows: [],
        };
        const sampleLabel = sampleRecord.sample_id.replace(/\.fcs$/i, "");
        const sample = sampleByLabel.get(sampleLabel.toUpperCase()) || sampleByLabel.get(sampleLabel);
        const row = document.createElement("tr");
        const values = [sampleLabel, sample?.role || "未匹配", qc.status,
          describeQcFindings(qc.findings),
          String((qc.time_windows || []).filter(window => window.rate_outlier).length),
          String(qc.events_examined ?? sampleRecord.event_count ?? "—")];
        for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
        qcTable.append(row);
        if (sample && !flowjoWaiting) {
          const previous = flowjoQcReviews.filter(item => item.sample_id === sample.id).at(-1);
          if (previous) {
            const reviewed = document.createElement("p"); reviewed.className = "muted";
            reviewed.textContent = `${sampleLabel} 上次人工 QC：${previous.status} · ${previous.reviewer} · ${previous.reason || "未填写理由"} · v${previous.version}`;
            data.append(reviewed);
          }
          const reviewForm = document.createElement("form");
          reviewForm.innerHTML = '<label>FlowJo 样本 QC 决定 <select name="status" required><option value="">请选择</option><option value="PASS">PASS</option><option value="PASS_WITH_OVERRIDE">PASS_WITH_OVERRIDE</option><option value="HOLD">HOLD</option><option value="EXCLUDE">EXCLUDE</option></select></label><label>审核人 <input name="reviewer" required maxlength="120"></label><label>理由 <input name="reason" maxlength="1000"></label><button>保存 FlowJo 样本 QC</button>';
          if (previous) {
            reviewForm.elements.status.value = previous.status;
            reviewForm.elements.reviewer.value = previous.reviewer;
            reviewForm.elements.reason.value = previous.reason;
          }
          const setReasonRequirement = () => {
            const status = reviewForm.elements.status.value;
            const flagsRequireReason = (qc.findings || []).some(item => item.severity === "REVIEW");
            reviewForm.elements.reason.required = ["PASS_WITH_OVERRIDE", "HOLD", "EXCLUDE"].includes(status)
              || (status === "PASS" && flagsRequireReason);
          };
          reviewForm.elements.status.addEventListener("change", setReasonRequirement); setReasonRequirement();
          reviewForm.addEventListener("submit", async (event) => {
            event.preventDefault(); const fields = new FormData(reviewForm);
            const button = reviewForm.querySelector("button"); button.disabled = true;
            try {
              await request(`${base}/qc-reviews`, {method: "POST", headers: {"Content-Type": "application/json"},
                body: JSON.stringify({flowjo_import_id: flowjoImportId, sample_id: sample.id,
                  status: fields.get("status"), reviewer: fields.get("reviewer"), reason: fields.get("reason"),
                  expected_revision: project.revision})});
              await showProject(projectId);
            } catch (error) { alert(`FlowJo QC 审核未保存：${error.message}`); button.disabled = false; }
          });
          data.append(reviewForm);
        }
      }
      const qcPolicy = document.createElement("p"); qcPolicy.className = "muted";
      qcPolicy.textContent = `自动 QC 仅提出审核提示，不会删除事件。策略：${JSON.stringify(records.qc_policy || {})}`;
      data.append(qcPolicy);
    } catch (error) {
      const failed = document.createElement("p"); failed.className = "error";
      failed.textContent = `FlowJo 预览读取失败：${error.message}`; data.append(failed);
    }
  }

  const percentageTitle = document.createElement("h3");
  percentageTitle.textContent = "已圈门 FlowJo 百分比结果";
  data.append(percentageTitle);
  const percentageNote = document.createElement("p"); percentageNote.className = "muted";
  percentageNote.textContent = "可直接导入已圈门的逐孔结果表，不读取或重算 FCS，也不要求导入 WSP。每份表独立承载其人工 QC；若关联同批 WSP，只作额外方法来源。系统保留表头门路径和指标，不推断或改写其分母语义；表格未提供事件分子/分母时不会补造计数。每次导入都需声明该表所依据的补偿状态；每孔仍需人工 QC 后才能计算描述性效应。";
  data.append(percentageNote);
  const percentageForm = document.createElement("form");
  percentageForm.innerHTML = '<label>批次 <select name="batch_id" required></select></label><label>用于生成此表的补偿状态 <select name="compensation_status" required><option value="UNKNOWN">未知</option><option value="NOT_APPLIED">未应用补偿</option><option value="APPLIED">已应用补偿</option></select></label><label>补偿设置说明（可选） <input name="compensation_note" maxlength="500"></label><label>FlowJo 导出结果表 <input name="table" type="file" accept=".xls,.xlsx,.csv,.tsv,.txt" required></label><label>无百分号数字尺度 <select name="numeric_scale"><option value="AUTO">自动（带 % 或 Excel 百分比格式）</option><option value="PERCENT">百分数 0–100</option><option value="FRACTION">比例 0–1</option></select></label><label>文本小数分隔符 <select name="decimal_separator"><option value=".">点（如 2.1%）</option><option value=",">逗号（如 2,1%）</option></select></label><button>导入已圈门结果</button>';
  for (const batch of batches.filter(item => item.assay === "ACTIVATION")) {
    const option = document.createElement("option"); option.value = batch.batch_id;
    option.textContent = `${batch.name} · ${batch.date}`; percentageForm.elements.batch_id.append(option);
  }
  percentageForm.querySelector("button").disabled = !batches.some(item => item.assay === "ACTIVATION");
  percentageForm.addEventListener("submit", async (event) => {
    event.preventDefault(); const button = percentageForm.querySelector("button"); button.disabled = true;
    try {
      const upload = new FormData();
      upload.append("batch_id", percentageForm.elements.batch_id.value);
      upload.append("numeric_scale", percentageForm.elements.numeric_scale.value);
      upload.append("decimal_separator", percentageForm.elements.decimal_separator.value);
      upload.append("compensation_status", percentageForm.elements.compensation_status.value);
      upload.append("compensation_note", percentageForm.elements.compensation_note.value);
      upload.append("expected_revision", String(project.revision));
      upload.append("table", percentageForm.elements.table.files[0], percentageForm.elements.table.files[0].name);
      const imported = await request(`${base}/flowjo-percentage-imports`, {method: "POST", body: upload});
      const compensationLabel = {NOT_APPLIED: "未应用补偿", APPLIED: "已应用补偿", UNKNOWN: "补偿状态未知"}[imported.result.compensation?.status] || "补偿状态未知";
      alert(`已保存 ${imported.result.imported_sample_count}/${imported.result.expected_sample_count} 个逐孔百分比；${compensationLabel}。事件分母未知；此结果表本身可作为人工 QC 来源${imported.flowjo_import_id ? "，并保留同批 FlowJo workspace 的可选关联。" : "，无需导入或处理 FCS/FlowJo workspace。"}`);
      await showProject(projectId);
    } catch (error) { alert(`FlowJo 百分比表导入失败：${error.message}`); button.disabled = false; }
  });
  data.append(percentageForm);

  let flowjoPercentageQcSaveInProgress = false;
  const flowjoPercentageQcFormButtons = [];
  for (const imported of flowjoPercentageImports) {
    const result = imported.result || {};
    const card = document.createElement("article");
    const title = document.createElement("h4");
    title.textContent = `${result.batch_name || imported.batch_id} · ${result.source_filename || "FlowJo 结果表"}`;
    card.append(title);
    const provenance = document.createElement("p"); provenance.className = "muted";
    const compensationLabels = {NOT_APPLIED: "未应用（用户报告）", APPLIED: "已应用（用户报告）", UNKNOWN: "未知"};
    const compensation = result.compensation || {status: "UNKNOWN"};
    const compensationLabel = compensationLabels[compensation.status] || "未知";
    const compensationDetail = compensation.note ? ` · 补偿说明：${compensation.note}` : "";
    provenance.textContent = `门路径：${result.gate_label} · ${result.metric} · ${result.imported_sample_count}/${result.expected_sample_count} 孔 · 仅百分比，分母未知 · 补偿：${compensationLabel}${compensationDetail} · SHA-256 ${result.source_sha256}`;
    card.append(provenance);
    const missingSampleIds = result.missing_sample_ids || [];
    if (missingSampleIds.length) {
      const incomplete = document.createElement("p"); incomplete.className = "error";
      const examples = missingSampleIds.slice(0, 6).join(", ");
      incomplete.textContent = `结果表缺少 ${missingSampleIds.length} 个已登记孔位${examples ? `（例如 ${examples}）` : ""}。此表不可生成效应；请补全并重新导入，缺失值不会按阴性或零值处理。`;
      card.append(incomplete);
    }
    const percentageImportId = imported.id || result.import_id;
    const sourceReviews = qcReviews.filter(item =>
      item.source_type === "FLOWJO_PERCENTAGE_IMPORT" && item.source_id === percentageImportId);
    const latestBySample = new Map();
    for (const review of sourceReviews) {
      const old = latestBySample.get(review.sample_id);
      if (!old || review.version > old.version) latestBySample.set(review.sample_id, review);
    }
    const observations = result.observations || [];
    let missingReviews = observations.filter(item => !latestBySample.has(item.sample_id)).length;
    const percentageReviewForms = [];
    let effectButton = null;
    let effectReviewStatus = null;
    let refreshEffectButton = null;
    const effects = await request(`${base}/flowjo-percentage-imports/${encodeURIComponent(percentageImportId)}/effect-analyses`);
    const latestEffectAnalysis = effects.at(-1);
    const latestPercentageEffect = latestEffectAnalysis?.result;
    const effectReviewIds = new Map((latestPercentageEffect?.effects || []).flatMap(item =>
      (item.replicate_observations || []).map(point => [point.sample_id, point.qc_review_id])));
    const isPercentageEffectStale = () => !latestPercentageEffect || !effectReviewIds.size
      || [...effectReviewIds].some(([sampleId, reviewId]) => latestBySample.get(sampleId)?.id !== reviewId);
    const reviewProgress = document.createElement("p");
    reviewProgress.className = "muted";
    reviewProgress.setAttribute("role", "status");
    const updateReviewProgress = () => {
      reviewProgress.textContent = `人工 QC：${observations.length - missingReviews}/${observations.length} 孔已审核 · ${missingReviews} 孔待审核`;
      if (effectButton) {
        effectButton.disabled = missingReviews > 0 || missingSampleIds.length > 0;
        effectButton.textContent = missingReviews
          ? `完成人工 QC 后计算标签级效应（尚有 ${missingReviews} 孔未审核）`
          : missingSampleIds.length
            ? `先补全结果表（缺少 ${missingSampleIds.length} 孔）`
            : "计算面板内描述性效应";
      }
      if (effectReviewStatus && refreshEffectButton) {
        const stale = isPercentageEffectStale();
        effectReviewStatus.textContent = stale
          ? "最新人工 QC 与此只读效应版本不一致；请基于当前审核重新计算。"
          : "此只读效应版本与当前纳入的人工 QC 决定一致。";
        refreshEffectButton.disabled = !stale || missingReviews > 0 || missingSampleIds.length > 0;
      }
    };
    updateReviewProgress();
    let reviewSaveInProgress = false;
    const table = document.createElement("table");
    const head = document.createElement("tr");
    for (const label of ["孔位", "临时样本标签", "面板/条件", "技术重复", "表格导出百分比", "事件分母", "人工 QC"]) {
      const cell = document.createElement("th"); cell.textContent = label; head.append(cell);
    }
    table.append(head);
    for (const observation of observations) {
      const row = document.createElement("tr");
      const sample = samples.find(item => item.id === observation.sample_id);
      const values = [observation.well_id, sample?.condition?.sample_label || observation.sample_id,
        `${observation.panel_label || "—"} / ${observation.condition_label || "—"}`,
        `${observation.replicate_type || "—"} ${observation.replicate_index ?? "—"}`,
        observation.displayed_value || `${Number(observation.percentage).toFixed(2)} %`, "未知"];
      for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
      const qcCell = document.createElement("td");
      const previous = latestBySample.get(observation.sample_id);
      const reviewForm = document.createElement("form"); reviewForm.className = "flowjo-percentage-qc-form";
      reviewForm.innerHTML = '<select name="status" aria-label="人工 QC 状态" required><option value="">待审核</option><option value="PASS">PASS</option><option value="PASS_WITH_OVERRIDE">PASS_WITH_OVERRIDE</option><option value="HOLD">HOLD</option><option value="EXCLUDE">EXCLUDE</option></select><input name="reviewer" aria-label="审核人" placeholder="审核人" required maxlength="120"><input name="reason" aria-label="理由" placeholder="理由（覆盖/保持/排除必填）" maxlength="1000"><button>保存</button><span class="flowjo-percentage-qc-save-status" role="status"></span>';
      const saveButton = reviewForm.querySelector("button");
      const saveStatus = reviewForm.querySelector(".flowjo-percentage-qc-save-status");
      flowjoPercentageQcFormButtons.push(saveButton);
      percentageReviewForms.push(reviewForm);
      if (previous) {
        reviewForm.elements.status.value = previous.status;
        reviewForm.elements.reviewer.value = previous.reviewer;
        reviewForm.elements.reason.value = previous.reason;
        saveStatus.textContent = `${previous.status} · v${previous.version} · ${previous.reviewer}`;
      }
      for (const field of [reviewForm.elements.status, reviewForm.elements.reviewer, reviewForm.elements.reason]) {
        field.addEventListener("input", () => { if (saveStatus.textContent) saveButton.textContent = "保存修改"; });
        field.addEventListener("change", () => { if (saveStatus.textContent) saveButton.textContent = "保存修改"; });
      }
      const updateReviewRequirement = () => {
        reviewForm.elements.reason.required = ["PASS_WITH_OVERRIDE", "HOLD", "EXCLUDE"].includes(reviewForm.elements.status.value);
      };
      reviewForm.elements.status.addEventListener("change", updateReviewRequirement); updateReviewRequirement();
      reviewForm.addEventListener("submit", async (event) => {
        event.preventDefault();
        if (reviewSaveInProgress || flowjoPercentageQcSaveInProgress) return;
        reviewSaveInProgress = true;
        flowjoPercentageQcSaveInProgress = true;
        flowjoPercentageQcFormButtons.forEach(button => { button.disabled = true; });
        saveButton.textContent = "保存中…";
        const hadPrevious = latestBySample.has(observation.sample_id);
        try {
          const saved = await request(`${base}/qc-reviews`, {method: "POST", headers: {"Content-Type": "application/json"},
            body: JSON.stringify({flowjo_percentage_import_id: percentageImportId, sample_id: observation.sample_id,
              status: reviewForm.elements.status.value, reviewer: reviewForm.elements.reviewer.value,
              reason: reviewForm.elements.reason.value, expected_revision: project.revision})});
          project.revision = saved.project_revision;
          latestBySample.set(observation.sample_id, {
            id: saved.qc_review_id, sample_id: saved.sample_id, status: saved.status,
            reviewer: saved.reviewer, reason: saved.reason, version: saved.version,
            project_revision: saved.project_revision,
          });
          if (!hadPrevious) missingReviews -= 1;
          saveStatus.textContent = `${saved.status} · v${saved.version} · ${saved.reviewer}`;
          saveButton.textContent = "已保存";
          updateReviewProgress();
        } catch (error) {
          alert(`结果表样本 QC 未保存：${error.message}`);
          saveButton.textContent = hadPrevious ? "保存修改" : "保存审核";
        } finally {
          reviewSaveInProgress = false;
          flowjoPercentageQcSaveInProgress = false;
          flowjoPercentageQcFormButtons.forEach(button => { button.disabled = false; });
        }
      });
      qcCell.append(reviewForm);
      row.append(qcCell); table.append(row);
    }
    const reviewerTools = document.createElement("div");
    reviewerTools.className = "flowjo-percentage-qc-reviewer-tools";
    const reviewerLabel = document.createElement("label");
    reviewerLabel.textContent = "本表默认审核人（只填空白项）";
    const reviewerDefault = document.createElement("input");
    reviewerDefault.type = "text";
    reviewerDefault.maxLength = 120;
    reviewerDefault.required = true;
    reviewerDefault.setAttribute("aria-label", "本表默认审核人");
    reviewerDefault.placeholder = "填写姓名后用于快速填充，不会自动提交 QC";
    reviewerLabel.append(reviewerDefault);
    const fillReviewers = document.createElement("button");
    fillReviewers.type = "button";
    fillReviewers.textContent = "填入本表空白审核人";
    fillReviewers.addEventListener("click", () => {
      const reviewer = reviewerDefault.value.trim();
      if (!reviewer) {
        reviewerDefault.reportValidity();
        reviewerDefault.focus();
        return;
      }
      const blankRows = percentageReviewForms.filter(form => !form.elements.reviewer.value.trim());
      for (const form of blankRows) form.elements.reviewer.value = reviewer;
      reviewerToolsStatus.textContent = `已填入 ${blankRows.length} 个空白审核人；未更改 QC 状态，也未保存。`;
    });
    flowjoPercentageQcFormButtons.push(fillReviewers);
    const reviewerToolsStatus = document.createElement("span");
    reviewerToolsStatus.className = "muted";
    reviewerToolsStatus.setAttribute("role", "status");
    reviewerTools.append(reviewerLabel, fillReviewers, reviewerToolsStatus);
    const worksheetTools = document.createElement("div");
    worksheetTools.className = "flowjo-qc-worksheet-tools";
    const downloadWorksheet = document.createElement("button");
    downloadWorksheet.type = "button";
    downloadWorksheet.textContent = "下载离线 QC 工作表";
    downloadWorksheet.addEventListener("click", async () => {
      downloadWorksheet.disabled = true;
      try {
        const response = await fetch(`${base}/flowjo-percentage-imports/${encodeURIComponent(percentageImportId)}/qc-worksheet`);
        if (!response.ok) {
          const body = await response.json().catch(() => ({}));
          throw new Error(body.detail?.message || body.detail || `HTTP ${response.status}`);
        }
        const blob = await response.blob();
        const objectUrl = URL.createObjectURL(blob);
        const anchor = document.createElement("a");
        anchor.href = objectUrl;
        anchor.download = `flowjo-qc-review-${percentageImportId.slice(0, 12)}.xlsx`;
        anchor.click();
        URL.revokeObjectURL(objectUrl);
        worksheetStatus.textContent = "已下载当前结果表的人工 QC 工作表。";
      } catch (error) {
        worksheetStatus.textContent = `下载失败：${error.message}`;
      } finally {
        downloadWorksheet.disabled = false;
      }
    });
    const worksheetInput = document.createElement("input");
    worksheetInput.type = "file";
    worksheetInput.accept = ".xlsx,application/vnd.openxmlformats-officedocument.spreadsheetml.sheet";
    worksheetInput.setAttribute("aria-label", "填写后的人工 QC 工作表");
    const importWorksheet = document.createElement("button");
    importWorksheet.type = "button";
    importWorksheet.textContent = "验证并导入人工 QC";
    importWorksheet.disabled = true;
    worksheetInput.addEventListener("change", () => {
      importWorksheet.disabled = !worksheetInput.files?.length;
    });
    const worksheetStatus = document.createElement("span");
    worksheetStatus.className = "muted";
    worksheetStatus.setAttribute("role", "status");
    importWorksheet.addEventListener("click", async () => {
      const file = worksheetInput.files?.[0];
      if (!file) return;
      const confirmed = confirm(
        "导入后会把工作表中填写的人工 QC 决定逐孔写入审计记录。空白决定保持未审核；任何来源信息不匹配都会整批拒绝。继续？",
      );
      if (!confirmed) return;
      importWorksheet.disabled = true;
      downloadWorksheet.disabled = true;
      worksheetInput.disabled = true;
      const form = new FormData();
      form.set("worksheet", file, file.name);
      try {
        const result = await request(
          `${base}/flowjo-percentage-imports/${encodeURIComponent(percentageImportId)}/qc-worksheet`,
          {method: "POST", body: form},
        );
        project.revision = result.project_revision;
        worksheetStatus.textContent = result.already_applied
          ? "这份工作表此前已成功导入；没有重复创建审核记录。"
          : `已导入 ${result.applied_count} 个 QC 决定；${result.unchanged_count} 行未新增决定。`;
        await showProject(projectId);
      } catch (error) {
        worksheetStatus.textContent = `未写入任何 QC：${error.message}`;
        importWorksheet.disabled = false;
        worksheetInput.disabled = false;
        downloadWorksheet.disabled = false;
      }
    });
    const worksheetNote = document.createElement("span");
    worksheetNote.className = "muted";
    worksheetNote.textContent = "XLSX 工作表只允许填写黄色 decision、reviewer、reason 栏；导入逐行校验来源，全部通过后才保存。";
    worksheetTools.append(downloadWorksheet, worksheetInput, importWorksheet, worksheetNote, worksheetStatus);
    card.append(reviewProgress, reviewerTools, worksheetTools, table);
    if (latestPercentageEffect?.effects?.length) {
      const effectTable = document.createElement("table");
      const effectHead = document.createElement("tr");
      for (const label of ["固定面板", "抗原临时标签", "重复类型", "WT 中位数", "突变体中位数", "NC", "突变−WT", "NC 校正下降", "合格孔 n（变体/WT/NC）", "逐孔分布与 QC"]) {
        const cell = document.createElement("th"); cell.textContent = label; effectHead.append(cell);
      }
      effectTable.append(effectHead);
      for (const item of latestPercentageEffect.effects) {
        const effect = item.batch_effect;
        const cellValues = [item.panel_label, item.probe_label, item.replicate_type,
          effect ? `${(100 * effect.reference_fraction).toFixed(2)}%` : "—",
          effect ? `${(100 * effect.mutant_fraction).toFixed(2)}%` : "—",
          effect ? `${(100 * effect.mutant_nc_fraction).toFixed(2)}%` : "—",
          effect ? `${effect.difference_pp.toFixed(2)} pp` : item.effect_status,
          effect?.decline_pct == null ? "无效窗口" : `${effect.decline_pct.toFixed(2)}%`,
          `${item.replicate_count_by_group.mutant}/${item.replicate_count_by_group.reference}/${item.replicate_count_by_group.nc}`];
        const row = document.createElement("tr");
        for (const value of cellValues) { const cell = document.createElement("td"); cell.textContent = String(value); row.append(cell); }
        const replicateCell = document.createElement("td");
        const replicateView = appendFlowJoPercentageReplicateReadout(item);
        if (replicateView) replicateCell.append(replicateView);
        else replicateCell.textContent = "此旧版结果未保存逐孔明细";
        row.append(replicateCell);
        effectTable.append(row);
      }
      card.append(effectTable);
      const effectNote = document.createElement("p"); effectNote.className = "muted";
      const effectCompensation = latestPercentageEffect.provenance?.compensation || compensation;
      const effectCompensationNote = effectCompensation.status === "NOT_APPLIED"
        ? "输入按用户报告的未补偿状态导出；串色可能影响门内分类。"
        : effectCompensation.status === "APPLIED"
          ? "输入按用户报告的已补偿状态导出；补偿矩阵未由百分比表独立验证。"
          : "补偿状态未知，不能假定不同补偿条件可直接比较。";
      effectNote.textContent = `仅为同批、同面板内同类型重复孔的中位数描述；输入为表格导出百分比，事件分母未知。${effectCompensationNote}未计算事件级误差、跨批显著性或跨抗体汇总，也不作结构定位。CD80 配体对照不纳入抗体面板效应。`;
      card.append(effectNote);
      effectReviewStatus = document.createElement("p");
      effectReviewStatus.className = "muted";
      effectReviewStatus.setAttribute("role", "status");
      refreshEffectButton = document.createElement("button");
      refreshEffectButton.type = "button";
      refreshEffectButton.textContent = "基于当前 QC 重新计算描述性效应";
      refreshEffectButton.addEventListener("click", async () => {
        refreshEffectButton.disabled = true;
        try {
          await request(`${base}/flowjo-percentage-imports/${encodeURIComponent(percentageImportId)}/effect-analyses`, {
            method: "POST", headers: {"Content-Type": "application/json"},
            body: JSON.stringify({expected_revision: project.revision}),
          });
          await showProject(projectId);
        } catch (error) {
          alert(`描述性效应未生成：${error.message}`);
          refreshEffectButton.disabled = false;
        }
      });
      card.append(effectReviewStatus, refreshEffectButton);
      updateReviewProgress();
    } else {
      effectButton = document.createElement("button");
      effectButton.type = "button";
      effectButton.addEventListener("click", async () => {
        effectButton.disabled = true;
        try {
          await request(`${base}/flowjo-percentage-imports/${encodeURIComponent(imported.id)}/effect-analyses`, {
            method: "POST", headers: {"Content-Type": "application/json"},
            body: JSON.stringify({expected_revision: project.revision}),
          });
          await showProject(projectId);
        } catch (error) { alert(`描述性效应未生成：${error.message}`); effectButton.disabled = false; }
      });
      card.append(effectButton);
      updateReviewProgress();
    }
    data.append(card);
  }

  const fcsTitle = document.createElement("h3"); fcsTitle.textContent = "FCS 导入与样本匹配"; data.append(fcsTitle);
  const fcsForm = document.createElement("form");
  fcsForm.innerHTML = '<label>批次 <select name="batch_id" required></select></label><label>FCS 文件 <input name="files" type="file" accept=".fcs" multiple required></label><button>解析并暂存</button>';
  const batchSelect = fcsForm.elements.batch_id;
  for (const batch of batches) {
    const option = document.createElement("option");
    option.value = batch.batch_id; option.textContent = `${batch.name} · ${batch.date}`;
    batchSelect.append(option);
  }
  fcsForm.querySelector("button").disabled = !batches.length;
  fcsForm.addEventListener("submit", async (event) => {
    event.preventDefault(); const button = fcsForm.querySelector("button"); button.disabled = true;
    try {
      const form = new FormData(); form.append("batch_id", fcsForm.elements.batch_id.value);
      form.append("expected_revision", String(project.revision));
      for (const file of fcsForm.elements.files.files) form.append("files", file, file.name);
      const imported = await request(`${base}/samples/import`, {method: "POST", body: form});
      renderPendingMatches(projectId, imported, data, {constructs, antibodies});
    } catch (error) { alert(`FCS 导入失败：${error.message}`); }
    finally { button.disabled = false; }
  });
  data.append(fcsForm);

  const plotTitle = document.createElement("h3"); plotTitle.textContent = "事件分布图（抽样显示）";
  data.append(plotTitle);
  const plotNote = document.createElement("p"); plotNote.className = "muted";
  plotNote.textContent = "图中使用原始、未补偿、未变换事件，只用于查看分布。门控计数和 QC 仍使用完整事件。横纵轴按全样本第 1–99 百分位显示，超出范围的点会裁切。";
  data.append(plotNote);
  const plotSamples = samples.filter(sample => sample.artifact_id && sample.fcs_metadata?.channels?.length);
  const plotForm = document.createElement("form");
  plotForm.innerHTML = '<label>样本 <select name="sample" required></select></label><label>X 通道 <select name="x" required></select></label><label>Y 通道 <select name="y" required></select></label><label>显示方式 <select name="mode"><option value="scatter">散点</option><option value="density">密度</option></select></label><label>抽样点数（100–10000） <input name="maximum_points" type="number" min="100" max="10000" step="100" value="3000" required></label><button type="submit">读取并绘图</button>';
  const plotSampleSelect = plotForm.elements.sample;
  for (const sample of plotSamples) {
    const option = document.createElement("option"); option.value = sample.id;
    option.textContent = `${sample.condition?.sample_label || sample.id} · ${sample.role}`;
    plotSampleSelect.append(option);
  }
  const updatePlotChannels = () => {
    const sample = plotSamples.find(item => item.id === plotSampleSelect.value);
    const channels = (sample?.fcs_metadata?.channels || []).map(channel => ({
      value: channel.pnn || channel.pns,
      label: channel.pns && channel.pns !== channel.pnn
        ? `${channel.pnn || channel.pns} · ${channel.pns}` : channel.pnn || channel.pns,
    })).filter(channel => channel.value);
    for (const axis of [plotForm.elements.x, plotForm.elements.y]) {
      axis.replaceChildren();
      for (const channel of channels) {
        const option = document.createElement("option"); option.value = channel.value;
        option.textContent = channel.label; axis.append(option);
      }
    }
    const xFsc = channels.find(channel => channel.value === "FSC-A");
    const ySsc = channels.find(channel => channel.value === "SSC-A");
    if (xFsc) plotForm.elements.x.value = xFsc.value;
    if (ySsc) plotForm.elements.y.value = ySsc.value;
    else if (channels.length > 1 && plotForm.elements.x.value === plotForm.elements.y.value) {
      plotForm.elements.y.selectedIndex = 1;
    }
    plotForm.querySelector("button").disabled = channels.length < 2;
  };
  plotSampleSelect.addEventListener("change", updatePlotChannels);
  updatePlotChannels();
  if (!plotSamples.length) {
    const option = document.createElement("option"); option.textContent = "尚无已匹配 FCS 样本";
    option.value = ""; plotSampleSelect.append(option); plotForm.querySelector("button").disabled = true;
  }
  data.append(plotForm);
  const plotStatus = document.createElement("p"); plotStatus.className = "muted";
  plotStatus.id = "sample-plot-status"; plotStatus.setAttribute("role", "status");
  plotStatus.textContent = "选择样本和通道后读取事件抽样。";
  data.append(plotStatus);
  const plotCanvas = document.createElement("canvas");
  plotCanvas.width = 900; plotCanvas.height = 560; plotCanvas.setAttribute("role", "img");
  plotCanvas.setAttribute("aria-label", "FCS 事件散点图或密度图");
  plotCanvas.style.width = "100%"; plotCanvas.style.maxWidth = "900px";
  plotCanvas.style.height = "auto"; plotCanvas.style.display = "block";
  data.append(plotCanvas);
  let latestPlotData = null;
  const drawEventPlot = (result, mode) => {
    const context = plotCanvas.getContext("2d");
    const width = plotCanvas.width, height = plotCanvas.height;
    const margin = {left: 88, right: 24, top: 28, bottom: 82};
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;
    const xDomain = result.domain.x, yDomain = result.domain.y;
    const xSpan = Math.max(Number.EPSILON, xDomain[1] - xDomain[0]);
    const ySpan = Math.max(Number.EPSILON, yDomain[1] - yDomain[0]);
    const toCanvas = point => [
      margin.left + ((point[0] - xDomain[0]) / xSpan) * plotWidth,
      margin.top + plotHeight - ((point[1] - yDomain[0]) / ySpan) * plotHeight,
    ];
    context.clearRect(0, 0, width, height);
    context.fillStyle = "#fff"; context.fillRect(0, 0, width, height);
    context.font = "13px system-ui, sans-serif"; context.textAlign = "right";
    for (let tick = 0; tick <= 4; tick++) {
      const x = margin.left + plotWidth * tick / 4;
      const y = margin.top + plotHeight - plotHeight * tick / 4;
      context.strokeStyle = "#e2e8ec"; context.lineWidth = 1;
      context.beginPath(); context.moveTo(x, margin.top); context.lineTo(x, margin.top + plotHeight); context.stroke();
      context.beginPath(); context.moveTo(margin.left, y); context.lineTo(margin.left + plotWidth, y); context.stroke();
      context.fillStyle = "#536772";
      context.textAlign = "center"; context.fillText((xDomain[0] + xSpan * tick / 4).toPrecision(4), x, margin.top + plotHeight + 20);
      context.textAlign = "right"; context.fillText((yDomain[0] + ySpan * tick / 4).toPrecision(4), margin.left - 9, y + 4);
    }
    const bins = 64, counts = new Uint32Array(bins * bins);
    const visiblePoints = [];
    for (const point of result.points || []) {
      if (!Array.isArray(point) || point.length < 2 || !point.every(Number.isFinite)) continue;
      if (point[0] < xDomain[0] || point[0] > xDomain[1] || point[1] < yDomain[0] || point[1] > yDomain[1]) continue;
      visiblePoints.push(point);
      const bx = Math.min(bins - 1, Math.floor((point[0] - xDomain[0]) / xSpan * bins));
      const by = Math.min(bins - 1, Math.floor((point[1] - yDomain[0]) / ySpan * bins));
      counts[by * bins + bx]++;
    }
    context.save(); context.beginPath(); context.rect(margin.left, margin.top, plotWidth, plotHeight); context.clip();
    if (mode === "density") {
      let maximum = 0; for (const count of counts) maximum = Math.max(maximum, count);
      const cellWidth = plotWidth / bins, cellHeight = plotHeight / bins;
      for (let by = 0; by < bins; by++) for (let bx = 0; bx < bins; bx++) {
        const count = counts[by * bins + bx]; if (!count) continue;
        const intensity = Math.sqrt(count / maximum);
        const red = Math.round(235 - 205 * intensity);
        const green = Math.round(247 - 112 * intensity);
        const blue = Math.round(252 - 48 * intensity);
        context.fillStyle = `rgb(${red},${green},${blue})`;
        context.fillRect(margin.left + bx * cellWidth, margin.top + plotHeight - (by + 1) * cellHeight,
          cellWidth + 0.5, cellHeight + 0.5);
      }
    } else {
      context.fillStyle = "rgba(7,95,141,0.22)";
      for (const point of visiblePoints) {
        const [x, y] = toCanvas(point); context.beginPath(); context.arc(x, y, 1.8, 0, Math.PI * 2); context.fill();
      }
    }
    context.restore();
    context.strokeStyle = "#536772"; context.lineWidth = 1.2;
    context.strokeRect(margin.left, margin.top, plotWidth, plotHeight);
    context.fillStyle = "#172b39"; context.textAlign = "center"; context.font = "14px system-ui, sans-serif";
    context.fillText(`${result.x_channel}（原始值）`, margin.left + plotWidth / 2, height - 32);
    context.save(); context.translate(23, margin.top + plotHeight / 2); context.rotate(-Math.PI / 2);
    context.fillText(`${result.y_channel}（原始值）`, 0, 0); context.restore();
    plotStatus.textContent = `${result.sample?.label || result.file_name} · ${result.sample?.role || "样本"}：完整事件 ${result.event_count.toLocaleString()}，有效配对事件 ${result.finite_event_count.toLocaleString()}，抽样 ${result.sampled_event_count.toLocaleString()}，当前范围内显示 ${visiblePoints.length.toLocaleString()}。抽样方法：${result.sampling_method}。`;
    plotStatus.className = "muted";
  };
  plotForm.addEventListener("submit", async event => {
    event.preventDefault(); const button = plotForm.querySelector("button"); button.disabled = true;
    plotStatus.textContent = "正在读取原始事件并生成确定性抽样…"; plotStatus.className = "muted";
    try {
      const result = await request(`${base}/samples/${encodeURIComponent(plotSampleSelect.value)}/plot-data`, {
        method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({schema_version: "1.0", x_channel: plotForm.elements.x.value,
          y_channel: plotForm.elements.y.value, maximum_points: Number(plotForm.elements.maximum_points.value)})
      });
      latestPlotData = result;
      drawEventPlot(latestPlotData, plotForm.elements.mode.value);
    } catch (error) { plotStatus.textContent = `绘图数据读取失败：${error.message}`; plotStatus.className = "error"; }
    finally { button.disabled = false; }
  });
  plotForm.elements.mode.addEventListener("change", () => {
    if (latestPlotData) drawEventPlot(latestPlotData, plotForm.elements.mode.value);
  });

  const fcsArtifacts = project.artifacts.filter(item => item.kind === "raw_fcs");
  const wspForm = document.createElement("form");
  wspForm.innerHTML = `<label>FlowJo workspace <input name="workspace" type="file" accept=".wsp" required></label><button ${fcsArtifacts.length ? "" : "disabled"}>导入 WSP 并预览门控</button>`;
  wspForm.addEventListener("submit", async (event) => {
    event.preventDefault(); const form = new FormData();
    form.append("workspace", wspForm.elements.workspace.files[0]);
    form.append("fcs_artifact_ids", JSON.stringify(fcsArtifacts.map(item => item.id)));
    form.append("expected_revision", String(project.revision));
    try {
      const imported = await request(`${base}/samples/import-flowjo`, {method: "POST", body: form});
      const counts = imported.preview.records.gate_counts || [];
      alert(`WSP ${imported.status}：${counts.length} 个完整事件门控计数；未解析门 ${imported.preview.records.unresolved_nodes.length}`);
      await showProject(projectId);
    } catch (error) { alert(`WSP 导入失败：${error.message}`); }
  });
  data.append(wspForm);

  const matrixArtifacts = project.artifacts.filter(item => item.kind === "compensation_matrix");
  const matrixTitle = document.createElement("h3"); matrixTitle.textContent = "外部补偿矩阵"; data.append(matrixTitle);
  const matrixNote = document.createElement("p"); matrixNote.className = "muted";
  matrixNote.textContent = "接受带行列通道标签的 .csv，或包含 detectors、matrix（可选 fluorochromes）的 .json。矩阵会按只读版本登记并绑定到门控修订。";
  data.append(matrixNote);
  const matrixList = document.createElement("ul");
  for (const artifact of matrixArtifacts) {
    const item = document.createElement("li");
    item.textContent = `${artifact.id.slice(0, 12)} · ${(artifact.detectors || []).join(", ")} · SHA-256 ${artifact.sha256.slice(0, 12)}`;
    matrixList.append(item);
  }
  if (!matrixArtifacts.length) {
    const item = document.createElement("li"); item.textContent = "尚未登记外部补偿矩阵。"; matrixList.append(item);
  }
  data.append(matrixList);
  const matrixForm = document.createElement("form");
  matrixForm.innerHTML = '<label>矩阵文件 (.csv 或 .json) <input name="matrix" type="file" accept=".csv,.json,text/csv,application/json" required></label><button>登记只读矩阵</button>';
  matrixForm.addEventListener("submit", async (event) => {
    event.preventDefault(); const button = matrixForm.querySelector("button"); button.disabled = true;
    try {
      const form = new FormData();
      form.append("matrix", matrixForm.elements.matrix.files[0]);
      form.append("expected_revision", String(project.revision));
      const registered = await request(`${base}/compensation-matrices`, {method: "POST", body: form});
      alert(`矩阵已登记：${registered.source_filename} · ${registered.detectors.join(", ")}`);
      await showProject(projectId);
    } catch (error) { alert(`矩阵登记失败：${error.message}`); button.disabled = false; }
  });
  data.append(matrixForm);

  const gateEditorTitle = document.createElement("h3"); gateEditorTitle.textContent = "自定义门控与修订"; data.append(gateEditorTitle);
  const gateNotice = document.createElement("p"); gateNotice.className = "muted";
  gateNotice.textContent = unresolvedFlowJoRebuild
    ? "此修订将作为 FlowJo 未解析门的人工重建版本，原始 WSP 与待处理节点记录会保留。内部编辑器支持无补偿、FCS 内补偿和已登记的外部矩阵；自动 QC 只标记异常，不会自动删除事件。"
    : "内部编辑器支持无补偿、FCS 内补偿和已登记的外部矩阵，门界限作用于线性通道；FlowJo workspace 门请用上方导入流程。自动 QC 只标记异常，不会自动删除时间窗或事件。";
  data.append(gateNotice);
  const gateHistory = document.createElement("ul");
  for (const revision of gateRevisions) {
    const item = document.createElement("li");
    const rebuild = revision.flowjo_rebuild
      ? ` · 重建 FlowJo ${revision.flowjo_rebuild.flowjo_import_id.slice(0, 8)} · ${revision.flowjo_rebuild.reviewer}` : "";
    item.textContent = `门控修订 v${revision.version} · ${revision.name} · ${revision.compensation_decision}${rebuild}`;
    gateHistory.append(item);
  }
  if (!gateRevisions.length) {
    const item = document.createElement("li"); item.textContent = "尚无内部门控修订。"; gateHistory.append(item);
  }
  data.append(gateHistory);
  if (gateAnalyses.length) {
    const latestAnalysis = gateAnalyses[gateAnalyses.length - 1];
    try {
      const response = await fetch(`${base}/files/${encodeURIComponent(latestAnalysis.result_artifact_id)}`);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const report = await response.json();
      const table = document.createElement("table");
      const head = document.createElement("tr");
      for (const value of ["样本", "角色", "门", "父门事件", "门内事件", "占父门比例"]) {
        const cell = document.createElement("th"); cell.textContent = value; head.append(cell);
      }
      table.append(head);
      const exclusionForm = document.createElement("form");
      const exclusionNote = document.createElement("p"); exclusionNote.className = "muted";
      exclusionNote.textContent = "自动 QC 只标记采集速率异常。操作人可勾选完整事件时间窗并填写审核人和理由；系统会保持原 FCS 不变，按排除后的完整事件重新门控并保存新分析版本。";
      exclusionForm.append(exclusionNote);
      for (const sampleResult of report.samples || []) {
        const qc = sampleResult.qc || {status: "QC_NOT_AVAILABLE", findings: []};
        const qcSummary = document.createElement("p");
        qcSummary.className = qc.status === "REVIEW_REQUIRED" ? "error" : "muted";
        qcSummary.textContent = `自动 QC：${qc.status} · 本版检查 ${qc.events_examined ?? sampleResult.event_count} 个事件 · 人工排除 ${qc.events_removed ?? 0} 个事件 · ${describeQcFindings(qc.findings)}`;
        data.append(qcSummary);
        const windows = qc.time_windows || [];
        if (windows.length) {
          const fieldset = document.createElement("fieldset");
          const legend = document.createElement("legend");
          legend.textContent = sampleResult.sample.condition?.sample_label || sampleResult.sample.sample_id;
          fieldset.append(legend);
          for (const window of windows) {
            const label = document.createElement("label");
          label.style.cssText = "display:flex; flex:1 1 100%; align-items:center; gap:6px;";
            const checkbox = document.createElement("input"); checkbox.type = "checkbox";
            checkbox.value = JSON.stringify({event_start: window.event_start, event_end: window.event_end});
            checkbox.dataset.sampleId = sampleResult.sample.sample_id;
            const span = window.time_span == null ? "—" : Number(window.time_span).toPrecision(4);
            const rate = window.event_rate == null ? "—" : Number(window.event_rate).toPrecision(4);
            label.append(checkbox, document.createTextNode(
              ` 事件 [${window.event_start}, ${window.event_end}) · 时间跨度 ${span} · 采集率 ${rate}${window.rate_outlier ? " · 自动标记异常" : ""}`
            ));
            fieldset.append(label);
          }
          exclusionForm.append(fieldset);
        }
        for (const gate of sampleResult.gate_counts || []) {
          const row = document.createElement("tr");
          const values = [sampleResult.sample.condition?.sample_label || sampleResult.sample.sample_id,
            sampleResult.sample.role, gate.gate_name, String(gate.parent_count), String(gate.count),
            gate.parent_fraction == null ? "—" : `${(100 * gate.parent_fraction).toFixed(2)}%`];
          for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
          table.append(row);
        }
        const review = qcReviews.filter(item => item.source_type === "GATE_ANALYSIS"
          && item.source_id === latestAnalysis.id
          && item.sample_id === sampleResult.sample.sample_id).at(-1);
        const reviewForm = document.createElement("form");
        reviewForm.innerHTML = '<label>QC 决定 <select name="status" required><option value="">请选择</option><option value="PASS">PASS</option><option value="PASS_WITH_OVERRIDE">PASS_WITH_OVERRIDE</option><option value="HOLD">HOLD</option><option value="EXCLUDE">EXCLUDE</option></select></label><label>审核人 <input name="reviewer" required maxlength="120"></label><label>理由 <input name="reason" maxlength="1000"></label><button>保存样本 QC 审核</button>';
        if (review) {
          reviewForm.elements.status.value = review.status;
          reviewForm.elements.reviewer.value = review.reviewer;
          reviewForm.elements.reason.value = review.reason;
          const previous = document.createElement("p"); previous.className = "muted";
          previous.textContent = `上次决定：${review.status} · ${review.reviewer} · ${review.reason || "未填写理由"} · v${review.version}`;
          data.append(previous);
        }
        const updateReason = () => {
          const status = reviewForm.elements.status.value;
          const flagsRequireReason = (sampleResult.qc?.findings || []).some(item => item.severity === "REVIEW");
          reviewForm.elements.reason.required = ["PASS_WITH_OVERRIDE", "HOLD", "EXCLUDE"].includes(status)
            || (status === "PASS" && flagsRequireReason);
        };
        reviewForm.elements.status.addEventListener("change", updateReason); updateReason();
        reviewForm.addEventListener("submit", async (event) => {
          event.preventDefault(); const fields = new FormData(reviewForm);
          const button = reviewForm.querySelector("button"); button.disabled = true;
          try {
            await request(`${base}/qc-reviews`, {method: "POST", headers: {"Content-Type": "application/json"},
              body: JSON.stringify({gate_analysis_id: latestAnalysis.id, sample_id: sampleResult.sample.sample_id,
                status: fields.get("status"), reviewer: fields.get("reviewer"), reason: fields.get("reason"),
                expected_revision: project.revision})});
            await showProject(projectId);
          } catch (error) { alert(`QC 审核未保存：${error.message}`); button.disabled = false; }
        });
        data.append(reviewForm);
      }
      const previousExclusions = report.event_exclusions;
      if (previousExclusions) {
        const history = document.createElement("p"); history.className = "muted";
        const excluded = (previousExclusions.windows || []).reduce((sum, item) => sum + (item.events_removed || 0), 0);
        history.textContent = `本分析继承已审核时间窗排除：${excluded} 个事件 · 审核人 ${previousExclusions.reviewer} · ${previousExclusions.reason}。再次选择会累积排除，并映射回原始 FCS 事件索引。`;
        data.append(history);
      }
      if (exclusionForm.querySelector('input[type="checkbox"]')) {
        const reviewerLabel = document.createElement("label"); reviewerLabel.textContent = "审核人";
        const reviewerInput = document.createElement("input"); reviewerInput.name = "reviewer"; reviewerInput.required = true; reviewerInput.maxLength = 120;
        reviewerLabel.append(reviewerInput);
        const reasonLabel = document.createElement("label"); reasonLabel.textContent = "时间窗排除理由";
        const reasonInput = document.createElement("input"); reasonInput.name = "reason"; reasonInput.required = true; reasonInput.maxLength = 1000;
        reasonLabel.append(reasonInput);
        const excludeButton = document.createElement("button"); excludeButton.textContent = "排除所选窗口并重新分析";
        exclusionForm.append(reviewerLabel, reasonLabel, excludeButton);
        exclusionForm.addEventListener("submit", async (event) => {
          event.preventDefault();
          const selectedBySample = new Map();
          for (const checkbox of exclusionForm.querySelectorAll('input[type="checkbox"]:checked')) {
            if (!selectedBySample.has(checkbox.dataset.sampleId)) selectedBySample.set(checkbox.dataset.sampleId, []);
            selectedBySample.get(checkbox.dataset.sampleId).push(JSON.parse(checkbox.value));
          }
          if (!selectedBySample.size) { alert("请至少选择一个完整时间窗。"); return; }
          excludeButton.disabled = true;
          try {
            const next = await request(`${base}/gate-analyses`, {method: "POST", headers: {"Content-Type": "application/json"},
              body: JSON.stringify({gate_revision_id: report.gate_revision_id,
                fcs_artifact_ids: report.input_artifact_ids, qc_policy: report.qc_policy,
                source_gate_analysis_id: latestAnalysis.id,
                window_exclusions: [...selectedBySample.entries()].map(([sample_id, windows]) => ({sample_id, windows})),
                exclusion_reviewer: reviewerInput.value, exclusion_reason: reasonInput.value,
                expected_revision: project.revision})});
            const removed = next.samples.reduce((sum, sample) => sum + (sample.events_removed || 0), 0);
            alert(`已人工排除 ${removed} 个事件并完成新门控分析；原始 FCS 与来源分析均保留。`);
            await showProject(projectId);
          } catch (error) { alert(`时间窗排除分析未完成：${error.message}`); excludeButton.disabled = false; }
        });
      } else {
        const unavailable = document.createElement("p"); unavailable.className = "muted";
        unavailable.textContent = "该分析没有可供人工选择的采集时间窗（例如缺少 Time 通道或事件数不足）；未排除任何事件。";
        exclusionForm.append(unavailable);
      }
      data.append(exclusionForm);
      if (constructs.length) {
        const effectTitle = document.createElement("h4");
        effectTitle.textContent = "Activation 效应分析"; data.append(effectTitle);
        const effectForm = document.createElement("form");
        const comparisonModeLabel = document.createElement("label");
        comparisonModeLabel.textContent = "比较类型";
        const comparisonMode = document.createElement("select"); comparisonMode.name = "comparisonMode";
        comparisonMode.innerHTML = '<option value="ANTIGEN_EFFECT">固定抗体：抗原突变 vs WT</option><option value="ANTIBODY_EFFECT">固定 WT 抗原：亲本 Fab vs 变体 Fab</option><option value="DOUBLE_MUTANT_INTERACTION">亲本/变体抗体 × WT/突变抗原四格交互</option>';
        comparisonModeLabel.append(comparisonMode); effectForm.append(comparisonModeLabel);
        const gateLabel = document.createElement("label"); gateLabel.textContent = "阳性门";
        const gateSelect = document.createElement("select"); gateSelect.name = "positiveGate";
        const gateChoices = [...new Map((report.samples || []).flatMap(item => item.gate_counts || [])
          .map(item => [item.gate_id, item.gate_name])).entries()];
        for (const [id, name] of gateChoices) {
          const option = document.createElement("option"); option.value = id; option.textContent = `${name} (${id})`;
          gateSelect.append(option);
        }
        gateLabel.append(gateSelect); effectForm.append(gateLabel);
        const assignmentTitle = document.createElement("p");
        assignmentTitle.className = "muted";
        assignmentTitle.textContent = "逐样本标记 MUTANT、REFERENCE 及对应 NC；未选择的样本不参加本次分析。只接受 ACTIVATION 批次和人工 PASS/PASS_WITH_OVERRIDE QC。";
        effectForm.append(assignmentTitle);
        const effectButton = document.createElement("button"); effectButton.textContent = "保存并计算批次效应";
        const assignments = [];
        for (const sampleResult of report.samples || []) {
          const row = document.createElement("div"); row.style.flex = "1 1 100%";
          const label = document.createElement("span");
          const sample = sampleResult.sample;
          const sampleAntigen = constructs.find(item => item.id === sample.antigen_construct_id);
          const sampleAntibody = antibodies.find(item => item.id === sample.antibody_construct_id);
          const temporaryIdentity = sample.condition?.temporary_identity;
          label.textContent = `${sample.condition?.sample_label || sample.sample_id} · ${sample.role} · 抗原 ${sampleAntigen?.name || temporaryIdentity?.antigen_label || "未关联"} · 面板 ${sampleAntibody?.name || temporaryIdentity?.panel_label || "未关联"}${temporaryIdentity ? "（临时标签）" : ""} `;
          const condition = document.createElement("select"); condition.setAttribute("aria-label", "效应条件");
          condition.innerHTML = '<option value="">不参加</option><option value="MUTANT">突变体</option><option value="REFERENCE">WT 参考</option><option value="MUTANT_NC">突变体 NC</option><option value="REFERENCE_NC">WT NC</option>';
          const probe = document.createElement("select"); probe.setAttribute("aria-label", "突变探针");
          const blankProbe = document.createElement("option"); blankProbe.value = ""; blankProbe.textContent = "选择突变构建体"; probe.append(blankProbe);
          for (const construct of constructs) {
            const option = document.createElement("option"); option.value = construct.id; option.textContent = construct.name;
            probe.append(option);
          }
          row.append(label, condition, probe); effectForm.append(row);
          assignments.push({sampleResult, condition, probe});
        }
        const updateComparisonMode = () => {
          const factorial = comparisonMode.value === "DOUBLE_MUTANT_INTERACTION";
          const antibodyMode = comparisonMode.value === "ANTIBODY_EFFECT";
          assignmentTitle.textContent = factorial
            ? "四格交互需在每个独立批次中，分别为亲本和变体抗体指定 MUTANT、REFERENCE 及各自匹配 NC；探针均选择同一 WT 背景下的抗原变体。将只报告描述性交互线索，不拆到位点。"
            : antibodyMode
              ? "固定一个 WT 抗原；亲本和变体 Fab 分别选择 activation 与各自 NC，同批比较 NC 校正阳性率窗口差。该读数不等同于亲和力或表面表达。"
              : "逐样本标记 MUTANT、REFERENCE 及对应 NC；未选择的样本不参加本次分析。只接受 ACTIVATION 批次和人工 PASS/PASS_WITH_OVERRIDE QC。";
          effectButton.textContent = factorial ? "保存并计算四格交互"
            : antibodyMode ? "保存并计算抗体变体效应" : "保存并计算批次效应";
          for (const assignment of assignments) {
            const selectedCondition = assignment.condition.value;
            const conditionOptions = antibodyMode
              ? [["", "不参加"], ["REFERENCE", "WT activation"], ["REFERENCE_NC", "WT NC"]]
              : [["", "不参加"], ["MUTANT", "突变体"], ["REFERENCE", "WT 参考"],
                ["MUTANT_NC", "突变体 NC"], ["REFERENCE_NC", "WT NC"]];
            assignment.condition.replaceChildren(...conditionOptions.map(([value, text]) => {
              const option = document.createElement("option"); option.value = value; option.textContent = text;
              return option;
            }));
            assignment.condition.value = conditionOptions.some(([value]) => value === selectedCondition)
              ? selectedCondition : "";
            const selectedProbe = assignment.probe.value;
            const validConstructs = antibodyMode
              ? constructs.filter(item => item.sequence === project.manifest.reference_sequence) : constructs;
            assignment.probe.replaceChildren();
            const blankProbe = document.createElement("option"); blankProbe.value = "";
            blankProbe.textContent = antibodyMode ? "选择固定 WT 抗原" : "选择突变构建体";
            assignment.probe.append(blankProbe);
            for (const construct of validConstructs) {
              const option = document.createElement("option"); option.value = construct.id; option.textContent = construct.name;
              assignment.probe.append(option);
            }
            assignment.probe.value = validConstructs.some(item => item.id === selectedProbe)
              ? selectedProbe : "";
          }
        };
        comparisonMode.addEventListener("change", updateComparisonMode);
        updateComparisonMode();
        const windowLabel = document.createElement("label"); windowLabel.textContent = "最小 WT 对照窗口（比例）";
        const windowInput = document.createElement("input"); windowInput.name = "minWindow";
        windowInput.type = "number"; windowInput.min = "0"; windowInput.max = "1";
        windowInput.step = "0.01"; windowInput.value = "0.10"; windowLabel.append(windowInput);
        effectForm.append(windowLabel, effectButton);
        effectForm.addEventListener("submit", async (event) => {
          event.preventDefault(); effectButton.disabled = true;
          const selected = assignments.filter(item => item.condition.value && item.probe.value).map(item => ({
            sample_id: item.sampleResult.sample.sample_id, probe_id: item.probe.value,
            condition: item.condition.value,
          }));
          try {
            const factorial = comparisonMode.value === "DOUBLE_MUTANT_INTERACTION";
            const antibodyMode = comparisonMode.value === "ANTIBODY_EFFECT";
            if (selected.length < (factorial ? 8 : 4)) {
              throw new Error(factorial
                ? "四格交互至少需要亲本/变体抗体各自的 MUTANT、REFERENCE 和两类 NC 样本"
                : antibodyMode
                  ? "抗体效应至少需要同一批亲本和变体 Fab 的 WT activation 与各自 NC（共 4 个样本）"
                  : "至少选择同批 MUTANT、REFERENCE 和两类 NC 样本");
            }
            if (selected.some(item => {
              const source = assignments.find(assignment => assignment.sampleResult.sample.sample_id === item.sample_id)?.sampleResult.sample;
              return !source?.antigen_construct_id || !source?.antibody_construct_id;
            })) {
              throw new Error("临时标签尚未对应到真实抗原和抗体构建体，当前不能计算效应或结构归属。");
            }
            await request(`${base}/effect-analyses`, {method: "POST", headers: {"Content-Type": "application/json"},
              body: JSON.stringify({gate_analysis_id: latestAnalysis.id, positive_gate_id: gateSelect.value,
                samples: selected, comparison_mode: comparisonMode.value,
                min_reference_window: Number(windowInput.value), expected_revision: project.revision})});
            await showProject(projectId);
          } catch (error) { alert(`效应分析未完成：${error.message}`); effectButton.disabled = false; }
        });
        data.append(effectForm);
      }
      const qcNotice = document.createElement("p"); qcNotice.className = "muted";
      qcNotice.textContent = `自动 QC 策略：非根门至少 500 事件、20 个采集时间窗、事件率偏离中位数 5 MAD、上限饱和事件超过 1% 时提出复核；下限堆积作提示。不会自动排除事件。人工 QC 修订留存决定和理由；此版本不改变门控计数。${report.qc_policy ? ` 参数：${JSON.stringify(report.qc_policy)}` : ""}`;
      data.append(qcNotice);
    } catch (error) {
      const failed = document.createElement("p"); failed.className = "error";
      failed.textContent = `门控分析读取失败：${error.message}`; data.append(failed);
    }
  }
  let latestEffectDocument = null;
  if (effectAnalyses.length) {
    const latestEffect = effectAnalyses[effectAnalyses.length - 1];
    const title = document.createElement("h3"); title.textContent = "已保存效应分析"; data.append(title);
    try {
      const response = await fetch(`${base}/files/${encodeURIComponent(latestEffect.result_artifact_id)}`);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const result = await response.json();
      latestEffectDocument = result;
      const factorialResult = result.comparison_mode === "DOUBLE_MUTANT_INTERACTION";
      const antibodyResult = result.comparison_mode === "ANTIBODY_EFFECT";
      if (!antibodyResult) {
        const table = document.createElement("table");
        const head = document.createElement("tr");
        const batchHeaders = ["探针", ...(factorialResult ? ["抗体构建体"] : []), "批次", "突变率", "WT 率", "百分点差", "保留率", "下降比例", "QC flags"];
        for (const value of batchHeaders) {
          const cell = document.createElement("th"); cell.textContent = value; head.append(cell);
        }
        table.append(head);
        for (const item of result.batch_effects || []) {
          const effect = item.batch_effect;
          const row = document.createElement("tr");
          const antibodyName = antibodies.find(antibody => antibody.id === item.antibody_construct_id)?.name
            || item.antibody_construct_id || "—";
          const values = [item.probe_id, ...(factorialResult ? [antibodyName] : []), effect.batch_id, `${(100 * effect.mutant_fraction).toFixed(2)}%`,
            `${(100 * effect.reference_fraction).toFixed(2)}%`, `${effect.difference_pp.toFixed(2)} pp`,
            effect.retention == null ? "—" : effect.retention.toFixed(3),
            effect.decline_pct == null ? "—" : `${effect.decline_pct.toFixed(2)}%`,
            (effect.qc_flags || []).join(", ") || "—"];
          for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
          table.append(row);
        }
        data.append(table);
      }
      if (result.batch_effects?.length) {
        const controlTitle = document.createElement("h4"); controlTitle.textContent = "每批四条件阳性率与 NC 校正窗口";
        const controlNote = document.createElement("p"); controlNote.className = "muted";
        controlNote.textContent = "条件比例由该批 QC 通过样本计算；校正窗口以百分点表示。逐批查看 WT 与突变抗原各自的 NC 基线。";
        data.append(controlTitle, controlNote);
        const controlTable = document.createElement("table");
        controlTable.setAttribute("aria-label", "每批四条件对照窗口");
        const controlHead = document.createElement("tr");
        const controlHeaders = ["探针", ...(factorialResult ? ["抗体构建体"] : []), "批次",
          "突变抗原", "WT 抗原", "突变抗原 NC", "WT NC", "突变抗原−NC 窗口", "WT−NC 窗口", "QC flags"];
        for (const value of controlHeaders) { const cell = document.createElement("th"); cell.textContent = value; controlHead.append(cell); }
        controlTable.append(controlHead);
        for (const item of result.batch_effects) {
          const effect = item.batch_effect;
          const row = document.createElement("tr");
          const antibodyName = antibodies.find(antibody => antibody.id === item.antibody_construct_id)?.name
            || item.antibody_construct_id || "—";
          const values = [item.probe_id, ...(factorialResult ? [antibodyName] : []), effect.batch_id,
            `${(100 * effect.mutant_fraction).toFixed(2)}%`, `${(100 * effect.reference_fraction).toFixed(2)}%`,
            `${(100 * effect.mutant_nc_fraction).toFixed(2)}%`, `${(100 * effect.reference_nc_fraction).toFixed(2)}%`,
            Number.isFinite(effect.mutant_control_window_pp) ? `${effect.mutant_control_window_pp.toFixed(2)} pp` : "—",
            Number.isFinite(effect.reference_control_window_pp) ? `${effect.reference_control_window_pp.toFixed(2)} pp` : "—",
            (effect.qc_flags || []).join(", ") || "—"];
          for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
          controlTable.append(row);
        }
        data.append(controlTable);
      }
      if (antibodyResult) {
        const antibodyTitle = document.createElement("h4");
        antibodyTitle.textContent = "亲本与变体 Fab：各自 NC 校正的 WT 阳性率窗口";
        const antibodyNote = document.createElement("p"); antibodyNote.className = "muted";
        antibodyNote.textContent = "每批比较两个抗体各自的 WT activation−NC 百分点窗口；正负方向只描述功能读数差异，不代表亲和力、折叠、表面表达或抗原空间位点。汇总以有效独立批次为 n，少于 3 批不显示区间。";
        data.append(antibodyTitle, antibodyNote);
        const antibodyTable = document.createElement("table");
        antibodyTable.setAttribute("aria-label", "亲本与变体抗体各自 NC 校正的 WT 窗口");
        const antibodyHead = document.createElement("tr");
        for (const value of ["固定 WT 抗原", "批次", "亲本构建体", "亲本 WT / NC", "亲本窗口",
          "变体构建体", "变体 WT / NC", "变体窗口", "变体−亲本", "QC flags"]) {
          const cell = document.createElement("th"); cell.textContent = value; antibodyHead.append(cell);
        }
        antibodyTable.append(antibodyHead);
        for (const item of result.antibody_effects || []) {
          const row = document.createElement("tr");
          const parentName = antibodies.find(value => value.id === item.parent_antibody_construct_id)?.name
            || item.parent_antibody_construct_id || "亲本";
          const variantName = antibodies.find(value => value.id === item.variant_antibody_construct_id)?.name
            || item.variant_antibody_construct_id || "变体";
          const probeName = constructs.find(value => value.id === item.probe_id)?.name || item.probe_id;
          const values = [probeName, item.batch_id, parentName,
            `${(100 * item.parent_wt_fraction).toFixed(2)}% / ${(100 * item.parent_wt_nc_fraction).toFixed(2)}%`,
            `${item.parent_control_window_pp.toFixed(2)} pp`, variantName,
            `${(100 * item.variant_wt_fraction).toFixed(2)}% / ${(100 * item.variant_wt_nc_fraction).toFixed(2)}%`,
            `${item.variant_control_window_pp.toFixed(2)} pp`, `${item.variant_minus_parent_pp.toFixed(2)} pp`,
            (item.qc_flags || []).join(", ") || "—"];
          for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
          antibodyTable.append(row);
        }
        data.append(antibodyTable);
        const summaryTable = document.createElement("table");
        summaryTable.setAttribute("aria-label", "跨批次抗体变体效应汇总");
        const summaryHead = document.createElement("tr");
        for (const value of ["固定 WT 抗原", "亲本构建体", "变体构建体", "有效 / 总批次",
          "变体−亲本窗口中位数", "有效批次 bootstrap 95% 区间", "正 / 负 / 零批次", "QC flags"]) {
          const cell = document.createElement("th"); cell.textContent = value; summaryHead.append(cell);
        }
        summaryTable.append(summaryHead);
        for (const summary of result.antibody_summaries || []) {
          const row = document.createElement("tr");
          const interval = summary.interval_pp;
          const directions = [summary.positive_batch_count, summary.negative_batch_count, summary.zero_batch_count];
          const probeName = constructs.find(value => value.id === summary.probe_id)?.name || summary.probe_id;
          const parentName = antibodies.find(value => value.id === summary.parent_antibody_construct_id)?.name
            || summary.parent_antibody_construct_id || "亲本";
          const variantName = antibodies.find(value => value.id === summary.variant_antibody_construct_id)?.name
            || summary.variant_antibody_construct_id || "变体";
          const values = [probeName, parentName, variantName, `${summary.valid_batch_count} / ${summary.batch_count}`,
            summary.median_variant_minus_parent_pp == null
              ? "无有效批次" : `${summary.median_variant_minus_parent_pp.toFixed(2)} pp`,
            interval ? `${interval[0].toFixed(2)} 至 ${interval[1].toFixed(2)} pp` : "不足 3 个有效批次",
            directions.every(Number.isInteger) ? directions.join(" / ") : "旧版未记录",
            (summary.qc_flags || []).join(", ") || "—"];
          for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
          summaryTable.append(row);
        }
        data.append(summaryTable);
      }
      if (!factorialResult && result.surface_projection?.probes?.length) {
        const summaryTitle = document.createElement("h4"); summaryTitle.textContent = "跨批次效应汇总";
        const summaryNote = document.createElement("p"); summaryNote.className = "muted";
        summaryNote.textContent = "每批是一个独立重复；原始百分点差中位数保留所有批次，下降中位数、方向计数和 percentile bootstrap 区间只根据有效批次计算。少于 3 个有效批次时不显示区间。";
        data.append(summaryTitle, summaryNote);
        const summaryTable = document.createElement("table");
        summaryTable.setAttribute("aria-label", "跨批次效应汇总");
        const summaryHead = document.createElement("tr");
        for (const value of ["探针", "独立批次 / 有效下降批次", "原始百分点差中位数", "下降比例中位数", "有效批次 bootstrap 95% 区间", "正 / 负 / 零批次", "结构归属"]) {
          const cell = document.createElement("th"); cell.textContent = value; summaryHead.append(cell);
        }
        summaryTable.append(summaryHead);
        for (const summary of result.surface_projection.probes) {
          const row = document.createElement("tr");
          const interval = summary.decline_interval_pct;
          const directionCounts = [summary.positive_decline_batch_count,
            summary.negative_decline_batch_count, summary.zero_decline_batch_count];
          const values = [summary.probe_id,
            `${summary.valid_decline_batch_count} / ${summary.batch_count}`,
            `${summary.difference_pp.toFixed(2)} pp`,
            summary.decline_pct == null ? "无有效批次" : `${summary.decline_pct.toFixed(2)}%`,
            interval ? `${interval[0].toFixed(2)}% 至 ${interval[1].toFixed(2)}%` : "不足 3 个有效批次",
            directionCounts.every(Number.isInteger) ? directionCounts.join(" / ") : "旧版未记录",
            summary.attribution];
          for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
          summaryTable.append(row);
        }
        data.append(summaryTable);
      }
      if (factorialResult) {
        const interactionTitle = document.createElement("h4");
        interactionTitle.textContent = "双侧交互：百分点差中之差"; data.append(interactionTitle);
        const interactionNote = document.createElement("p"); interactionNote.className = "muted";
        interactionNote.textContent = "每批展示亲本和变体抗体各自的突变抗原/WT 四格原始阳性率与差中之差。汇总只把独立批次作为 n；交互量是耦合线索，不拆分到抗原位点，也不作为亲和力。";
        data.append(interactionNote);
        const batchInteractionTable = document.createElement("table");
        const interactionHead = document.createElement("tr");
        for (const value of ["探针", "批次", "亲本 Ab Mut/WT", "变体 Ab Mut/WT", "交互 DID", "QC flags"]) {
          const cell = document.createElement("th"); cell.textContent = value; interactionHead.append(cell);
        }
        batchInteractionTable.append(interactionHead);
        for (const item of result.factorial_interactions || []) {
          const parentCells = item.four_cells.parent_antibody;
          const variantCells = item.four_cells.variant_antibody;
          const row = document.createElement("tr");
          const values = [item.probe_id, item.batch_id,
            `${(100 * parentCells.mutant_antigen).toFixed(2)}% / ${(100 * parentCells.wt_antigen).toFixed(2)}%`,
            `${(100 * variantCells.mutant_antigen).toFixed(2)}% / ${(100 * variantCells.wt_antigen).toFixed(2)}%`,
            `${item.interaction_pp.toFixed(2)} pp`, (item.qc_flags || []).join(", ") || "—"];
          for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
          batchInteractionTable.append(row);
        }
        data.append(batchInteractionTable);
        const summaryTable = document.createElement("table");
        const summaryHead = document.createElement("tr");
        for (const value of ["探针", "独立批次数", "DID 中位数", "批次 bootstrap 95% 区间", "正/负/零批次"]) {
          const cell = document.createElement("th"); cell.textContent = value; summaryHead.append(cell);
        }
        summaryTable.append(summaryHead);
        for (const summary of result.factorial_summaries || []) {
          const interval = summary.interaction_interval_pp;
          const row = document.createElement("tr");
          const values = [summary.probe_id, String(summary.batch_count),
            `${summary.median_interaction_pp.toFixed(2)} pp`,
            interval ? `${interval[0].toFixed(2)} 至 ${interval[1].toFixed(2)} pp` : "不足 3 个独立批次",
            `${summary.positive_batch_count} / ${summary.negative_batch_count} / ${summary.zero_batch_count}`];
          for (const value of values) { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); }
          summaryTable.append(row);
        }
        data.append(summaryTable);
      }
      const heatmapRows = result.batch_effects || [];
      const heatmapBatchIds = [...new Set(heatmapRows.map(item => item.batch_effect?.batch_id).filter(Boolean))].sort();
      const heatmapProbeIds = [...new Set(heatmapRows.map(item => item.probe_id))].sort();
      if (!factorialResult && heatmapBatchIds.length && heatmapProbeIds.length) {
        const heatmapTitle = document.createElement("h4");
        heatmapTitle.textContent = "探针 × 独立批次下降热图"; data.append(heatmapTitle);
        const heatmapNote = document.createElement("p"); heatmapNote.className = "muted";
        heatmapNote.textContent = "每格展示该探针在该独立批次中的下降比例；无效参考窗口显示为无效，不填补、不合并重复孔。悬停可查看完整指标和 QC flags。";
        data.append(heatmapNote);
        const heatmap = document.createElement("table");
        heatmap.setAttribute("aria-label", "逐探针独立批次效应热图");
        const heatmapHead = document.createElement("tr");
        const probeHeading = document.createElement("th"); probeHeading.scope = "col";
        probeHeading.textContent = "探针"; heatmapHead.append(probeHeading);
        for (const batchId of heatmapBatchIds) {
          const batch = batches.find(item => item.id === batchId);
          const heading = document.createElement("th"); heading.scope = "col";
          heading.textContent = batch?.name || batchId.slice(0, 10);
          heading.title = batch ? `${batch.name} · ${batch.date}` : batchId;
          heatmapHead.append(heading);
        }
        heatmap.append(heatmapHead);
        for (const probeId of heatmapProbeIds) {
          const row = document.createElement("tr");
          const label = document.createElement("th"); label.scope = "row"; label.textContent = probeId;
          row.append(label);
          for (const batchId of heatmapBatchIds) {
            const record = heatmapRows.find(item => item.probe_id === probeId
              && item.batch_effect?.batch_id === batchId);
            const effect = record?.batch_effect;
            const cell = document.createElement("td");
            if (!effect || effect.decline_pct == null) {
              cell.textContent = "无效"; cell.style.backgroundColor = "#e9eef1";
            } else {
              const decline = Number(effect.decline_pct);
              cell.textContent = `${decline.toFixed(1)}%`;
              if (decline < 0) cell.style.backgroundColor = "#d9e8f4";
              else {
                const intensity = Math.max(0, Math.min(100, decline)) / 100;
                const base = [250, 241, 235], target = [177, 44, 38];
                cell.style.backgroundColor = `rgb(${base.map((value, index) => Math.round(value * (1 - intensity) + target[index] * intensity)).join(",")})`;
              }
            }
            cell.title = effect
              ? `突变率 ${(100 * effect.mutant_fraction).toFixed(2)}% · WT 率 ${(100 * effect.reference_fraction).toFixed(2)}% · 差 ${effect.difference_pp.toFixed(2)} pp · 保留率 ${effect.retention == null ? "无效" : effect.retention.toFixed(3)} · ${effect.qc_flags?.join(", ") || "无 QC flag"}`
              : "该探针在此批次没有效应记录";
            row.append(cell);
          }
          heatmap.append(row);
        }
        data.append(heatmap);
      }
      const note = document.createElement("p"); note.className = "muted";
      note.textContent = "效应是指定门和条件下的功能阳性率变化；不是亲和力、接触概率或组合突变单个位点的效应。结果绑定至冻结结构哈希和 QC 审核记录。";
      data.append(note);
    } catch (error) {
      const failed = document.createElement("p"); failed.className = "error";
      failed.textContent = `效应分析读取失败：${error.message}`; data.append(failed);
    }
  }
  if (latestEffectDocument) appendWellEffectReadout(latestEffectDocument, antibodies, data);
  await renderStructureEffects(
    projectId, project, constructs, antibodies, latestEffectDocument, summaryMeasurements, qcReviews, data,
  );

  const matchedArtifacts = unresolvedFlowJoRebuild
    ? unresolvedFlowJoRebuild.input_artifact_ids
    : [...new Set(samples.filter(sample => sample.artifact_id).map(sample => sample.artifact_id))];
  const channelsBySample = samples.filter(sample => matchedArtifacts.includes(sample.artifact_id)).map(sample =>
    (sample.fcs_metadata?.channels || []).map(channel => channel.pnn || channel.pns).filter(Boolean));
  const channelChoices = channelsBySample.length
    ? channelsBySample.slice(1).reduce((common, channels) => common.filter(channel => channels.includes(channel)), channelsBySample[0])
    : [];
  const gateForm = document.createElement("form");
  gateForm.innerHTML = '<label>门名称 <input name="name" value="Gate1" required maxlength="120"></label><label>门类型 <select name="kind"><option value="RECTANGLE_2D">二维矩形</option><option value="RANGE_1D">一维范围</option><option value="POLYGON_2D">二维多边形</option></select></label><label>X 通道 <select name="x" required></select></label><label data-gate-input="xmin">X 最小值 <input name="xmin" type="number" step="any" required></label><label data-gate-input="xmax">X 最大值 <input name="xmax" type="number" step="any" required></label><label data-gate-input="y">Y 通道 <select name="y" required></select></label><label data-gate-input="ymin">Y 最小值 <input name="ymin" type="number" step="any" required></label><label data-gate-input="ymax">Y 最大值 <input name="ymax" type="number" step="any" required></label><label data-gate-input="vertices" hidden>顶点坐标 JSON <textarea name="vertices" placeholder="[[0,0],[100,0],[100,100],[0,100]]"></textarea></label><label>该门变换 ID（可空）<input name="transformId" maxlength="120"></label><label>本修订的变换定义 JSON <textarea name="transforms" placeholder="[{&quot;transform_id&quot;:&quot;logicle-v1&quot;,&quot;kind&quot;:&quot;LOGICLE&quot;,&quot;parameters&quot;:{&quot;param_t&quot;:262144,&quot;param_w&quot;:0.5,&quot;param_m&quot;:4.5,&quot;param_a&quot;:0}}]"></textarea><span class="muted">门边界使用变换后的通道坐标。支持 LOGICLE、ASINH、HYPERLOG、LINEAR、LOG、WSP_BIEXP、WSP_LOG。</span></label><label>补偿来源 <select name="compensation" required><option value="">请明确选择</option><option value="NONE">无补偿</option><option value="FCS_EMBEDDED">FCS 内补偿矩阵</option></select></label><label>完整层级门树 JSON（可选，填写后使用此门树） <textarea name="gateTree" placeholder="[{&quot;gate_id&quot;:&quot;live&quot;,&quot;parent_gate_id&quot;:&quot;root&quot;,&quot;name&quot;:&quot;Live&quot;,&quot;kind&quot;:&quot;POLYGON_2D&quot;,&quot;channels&quot;:[&quot;FSC-A&quot;,&quot;SSC-A&quot;],&quot;parameters&quot;:{&quot;vertices&quot;:[[0,0],[1,0],[1,1],[0,1]]}}]"></textarea></label><button>保存修订并对已匹配 FCS 执行</button>';
  const externalOption = document.createElement("option");
  externalOption.value = "EXTERNAL_MATRIX"; externalOption.textContent = "外部补偿矩阵";
  gateForm.elements.compensation.append(externalOption);
  const matrixLabel = document.createElement("label");
  matrixLabel.dataset.gateInput = "matrix"; matrixLabel.hidden = true;
  matrixLabel.append(document.createTextNode("已登记矩阵 "));
  const matrixSelect = document.createElement("select"); matrixSelect.name = "compensationArtifact";
  for (const artifact of matrixArtifacts) {
    const option = document.createElement("option"); option.value = artifact.id;
    option.textContent = `${artifact.id.slice(0, 12)} · ${(artifact.detectors || []).join(", ")}`;
    matrixSelect.append(option);
  }
  matrixLabel.append(matrixSelect);
  gateForm.insertBefore(matrixLabel, gateForm.querySelector("button"));
  for (const axis of ["x", "y"]) {
    const select = gateForm.elements[axis];
    for (const channel of channelChoices) {
      const option = document.createElement("option"); option.value = channel; option.textContent = channel; select.append(option);
    }
  }
  if (channelChoices.includes("FSC-A")) gateForm.elements.x.value = "FSC-A";
  if (channelChoices.includes("SSC-A")) gateForm.elements.y.value = "SSC-A";
  if (unresolvedFlowJoRebuild) {
    const rebuildFieldset = document.createElement("fieldset");
    rebuildFieldset.innerHTML = '<legend>FlowJo 未解析门重建审核</legend><label>审核人 <input name="flowjoReviewer" required maxlength="120"></label><label>重建依据 <textarea name="flowjoReason" required maxlength="1000" placeholder="说明如何将所列节点完整重建到下方门树"></textarea></label><label><input name="flowjoAcknowledged" type="checkbox" required>我已核对未解析节点，并确认下方完整门树覆盖此 workspace 的层级逻辑</label>';
    gateForm.insertBefore(rebuildFieldset, gateForm.querySelector("button"));
    gateForm.elements.gateTree.required = true;
  }
  const updateGateFields = () => {
    const kind = gateForm.elements.kind.value;
    gateForm.querySelector('[data-gate-input="y"]').hidden = kind === "RANGE_1D";
    gateForm.querySelector('[data-gate-input="ymin"]').hidden = kind !== "RECTANGLE_2D";
    gateForm.querySelector('[data-gate-input="ymax"]').hidden = kind !== "RECTANGLE_2D";
    gateForm.querySelector('[data-gate-input="vertices"]').hidden = kind !== "POLYGON_2D";
    gateForm.elements.vertices.required = kind === "POLYGON_2D";
    gateForm.elements.y.required = kind !== "RANGE_1D";
    gateForm.elements.ymin.required = kind === "RECTANGLE_2D";
    gateForm.elements.ymax.required = kind === "RECTANGLE_2D";
    gateForm.elements.xmin.closest("label").hidden = kind === "POLYGON_2D";
    gateForm.elements.xmax.closest("label").hidden = kind === "POLYGON_2D";
    const externalCompensation = gateForm.elements.compensation.value === "EXTERNAL_MATRIX";
    matrixLabel.hidden = !externalCompensation;
    matrixSelect.required = externalCompensation;
    matrixSelect.disabled = !externalCompensation || !matrixArtifacts.length;
  };
  gateForm.elements.kind.addEventListener("change", updateGateFields);
  gateForm.elements.compensation.addEventListener("change", updateGateFields);
  updateGateFields();
  const matchedSampleArtifactIds = new Set(samples.filter(sample => sample.artifact_id).map(sample => sample.artifact_id));
  const rebuildFilesMatched = !unresolvedFlowJoRebuild
    || matchedArtifacts.every(artifactId => matchedSampleArtifactIds.has(artifactId));
  const gateButton = gateForm.querySelector("button");
  gateButton.disabled = !matchedArtifacts.length || channelChoices.length < 2 || !rebuildFilesMatched
    || (gateForm.elements.compensation.value === "EXTERNAL_MATRIX" && !matrixArtifacts.length);
  if (!rebuildFilesMatched) {
    gateNotice.textContent += " WSP 绑定的 FCS 中仍有文件未完成样本匹配，请先匹配这些文件。";
  }
  gateForm.addEventListener("submit", async (event) => {
    event.preventDefault(); gateButton.disabled = true;
    const fields = new FormData(gateForm); let savedRevision = false;
    try {
      const kind = fields.get("kind");
      let channels = [fields.get("x")];
      let parameters;
      if (kind === "RECTANGLE_2D") {
        channels.push(fields.get("y"));
        parameters = {bounds: [
          {min: Number(fields.get("xmin")), max: Number(fields.get("xmax"))},
          {min: Number(fields.get("ymin")), max: Number(fields.get("ymax"))}
        ]};
      } else if (kind === "RANGE_1D") {
        parameters = {min: Number(fields.get("xmin")), max: Number(fields.get("xmax"))};
      } else {
        channels.push(fields.get("y"));
        parameters = {vertices: JSON.parse(fields.get("vertices"))};
      }
      const transformId = fields.get("transformId").trim() || null;
      const transforms = fields.get("transforms").trim() ? JSON.parse(fields.get("transforms")) : [];
      if (fields.get("compensation") === "EXTERNAL_MATRIX" && !fields.get("compensationArtifact")) {
        throw new Error("请先登记并选择外部补偿矩阵");
      }
      if (unresolvedFlowJoRebuild && !fields.get("gateTree").trim()) {
        throw new Error("重建 FlowJo 未解析门时必须提供完整层级门树 JSON");
      }
      const gateTree = fields.get("gateTree").trim()
        ? JSON.parse(fields.get("gateTree"))
        : [{gate_id: "Gate1", parent_gate_id: "root", name: fields.get("name"), kind,
            channels, parameters, ...(transformId ? {transform_id: transformId} : {})}];
      const saved = await request(`${base}/gates`, {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({name: fields.get("name"), compensation_decision: fields.get("compensation"),
          ...(fields.get("compensation") === "EXTERNAL_MATRIX"
            ? {compensation_artifact_id: fields.get("compensationArtifact")} : {}),
          transforms, expected_revision: project.revision, gate_tree: gateTree,
          ...(unresolvedFlowJoRebuild ? {flowjo_rebuild: {
            flowjo_import_id: unresolvedFlowJoRebuild.flowjo_import_id,
            reviewer: fields.get("flowjoReviewer"), reason: fields.get("flowjoReason")
          }} : {})})});
      savedRevision = true;
      const analysis = await request(`${base}/gate-analyses`, {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({gate_revision_id: saved.gate_revision_id, fcs_artifact_ids: matchedArtifacts,
          expected_revision: saved.project_revision})});
      alert(unresolvedFlowJoRebuild
        ? `已对绑定的 ${analysis.samples.length} 个样本按完整内部门树重建并重算；审核记录已关联原始 FlowJo 导入。`
        : `已对 ${analysis.samples.length} 个样本执行完整事件门控；结果已保存为只读分析 artifact。`);
      await showProject(projectId);
    } catch (error) {
      alert(`门控修订或执行失败：${error.message}`);
      if (savedRevision) await showProject(projectId);
      else gateButton.disabled = false;
    }
  });
  data.append(gateForm);
  const exportSection = document.createElement("section");
  const exportTitle = document.createElement("h3"); exportTitle.textContent = "GRASP IR 导出";
  const exportNotice = document.createElement("p"); exportNotice.className = "muted";
  exportNotice.textContent = "实验 IR 只开放给具有有效下降值的真实单点探针；Paraplume IR 仍标记为预测辅助证据。当前未登记复合物结构，因此不提供 RPR。";
  const mappingByUniprot = new Map();
  const mappingArtifact = project.artifacts.find(item => item.kind === "mapping");
  if (mappingArtifact) {
    try {
      const response = await fetch(`${base}/files/${encodeURIComponent(mappingArtifact.id)}`);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);
      const lines = (await response.text()).trim().split(/\r?\n/);
      const columns = lines[0].split(",");
      for (const line of lines.slice(1)) {
        const fields = line.split(",");
        const row = Object.fromEntries(columns.map((key, index) => [key, fields[index] || ""]));
        if (row.uniprot_pos && row.normalized_position) mappingByUniprot.set(Number(row.uniprot_pos), row);
      }
    } catch (error) { exportNotice.textContent += ` 映射表读取失败：${error.message}`; }
  }
  const effectProbeOptions = [];
  const latestEffectId = effectAnalyses.at(-1)?.id;
  for (const probe of latestEffectDocument?.surface_projection?.probes || []) {
    const position = probe.mutated_positions?.length === 1 ? Number(probe.mutated_positions[0]) : null;
    const mapped = position ? mappingByUniprot.get(position) : null;
    if (probe.attribution === "SINGLE_SITE" && probe.decline_pct != null
        && probe.valid_decline_batch_count > 0 && mapped?.one_letter) {
      effectProbeOptions.push({probe, position, mapped});
    }
  }
  const exportForm = document.createElement("form");
  exportForm.innerHTML = '<label>证据类型 <select name="source"><option value="predicted_ir">Paraplume 预测 IR</option><option value="experimental_ir">实验单点 IR</option></select></label><label>抗体构建体 <select name="antibody" required></select></label><label data-prediction-field>链 <select name="chain"><option value="2">VH</option><option value="3">VL</option></select></label><label data-prediction-field>位置（1 起始）<input name="position" type="number" min="1"></label><label data-effect-field hidden>实验探针 <select name="probe"></select></label><label>审阅人 <input name="reviewer" required maxlength="120"></label><label>选择理由 <input name="rationale" required maxlength="1000"></label><button>校验并生成 ZIP</button>';
  const exportAntibody = exportForm.elements.antibody;
  for (const antibody of antibodies) {
    const run = predictionsByConstruct.get(antibody.id);
    const option = document.createElement("option"); option.value = antibody.id;
    option.textContent = `${antibody.name} · ${run?.prediction_status || "无 Paraplume 结果"}`; exportAntibody.append(option);
  }
  const probeSelect = exportForm.elements.probe;
  for (const item of effectProbeOptions) {
    const option = document.createElement("option"); option.value = item.probe.probe_id;
    option.textContent = `${item.probe.probe_id} · UniProt ${item.position} → ECD ${item.mapped.normalized_position} · 下降 ${Number(item.probe.decline_pct).toFixed(1)}%`;
    probeSelect.append(option);
  }
  const sourceSelect = exportForm.elements.source;
  if (!effectProbeOptions.length) sourceSelect.querySelector('[value="experimental_ir"]').disabled = true;
  if (!paraplumeRuns.length) sourceSelect.value = "experimental_ir";
  const exportPreview = document.createElement("pre"); exportPreview.className = "muted";
  const exportResult = document.createElement("p"); exportResult.setAttribute("role", "status");
  const updateExportPreview = () => {
    const predicted = sourceSelect.value === "predicted_ir";
    exportForm.querySelectorAll("[data-prediction-field]").forEach(label => { label.hidden = !predicted; });
    exportForm.querySelectorAll("[data-effect-field]").forEach(label => { label.hidden = predicted; });
    exportForm.elements.position.required = predicted;
    if (exportButton) exportButton.disabled = predicted
      ? !predictionsByConstruct.has(exportAntibody.value)
      : !effectProbeOptions.length;
    if (!predicted) {
      const item = effectProbeOptions.find(candidate => candidate.probe.probe_id === probeSelect.value);
      exportPreview.textContent = item
        ? `FASTA chain 1 · ${item.mapped.normalized_position} · ${item.mapped.one_letter}\n1-${item.mapped.normalized_position}-${item.mapped.one_letter}\nSINGLE_SITE · decline ${Number(item.probe.decline_pct).toFixed(2)}% · ${item.probe.valid_decline_batch_count} valid batches`
        : "没有可用的已验证实验单点探针。";
      return;
    }
    const antibody = antibodies.find(item => item.id === exportAntibody.value);
    const chain = Number(exportForm.elements.chain.value);
    const sequence = antibody ? (chain === 2 ? antibody.vh_sequence : antibody.vl_sequence) : "";
    const position = Number(exportForm.elements.position.value);
    exportPreview.textContent = Number.isInteger(position) && position >= 1 && position <= sequence.length
      ? `FASTA chain ${chain} · ${position} · ${sequence[position - 1]}\n${chain}-${position}-${sequence[position - 1]}\nPREDICTED_UNREVIEWED`
      : "位置须落在所选链序列范围内。";
  };
  exportForm.addEventListener("input", updateExportPreview);
  exportForm.addEventListener("change", updateExportPreview);
  const exportButton = exportForm.querySelector("button");
  exportButton.disabled = !antibodies.length || (!effectProbeOptions.length && !paraplumeRuns.length);
  updateExportPreview();
  if (!effectProbeOptions.length && !paraplumeRuns.length) exportNotice.textContent += " 尚无可导出的已验证实验单点或 Paraplume 结果。";
  sourceSelect.addEventListener("change", updateExportPreview);
  exportForm.addEventListener("submit", async event => {
    event.preventDefault();
    const predicted = sourceSelect.value === "predicted_ir";
    const antibody = antibodies.find(item => item.id === exportAntibody.value);
    const run = predictionsByConstruct.get(antibody?.id);
    const chain = Number(exportForm.elements.chain.value);
    const sequence = antibody ? (chain === 2 ? antibody.vh_sequence : antibody.vl_sequence) : "";
    const position = Number(exportForm.elements.position.value);
    const selectedProbe = effectProbeOptions.find(item => item.probe.probe_id === probeSelect.value);
    if (predicted && (!run || !Number.isInteger(position) || position < 1 || position > sequence.length)) return;
    if (!predicted && !selectedProbe) return;
    exportButton.disabled = true; exportResult.textContent = "正在验证来源和 FASTA 编号…";
    try {
      const result = await request(`${base}/constraint-exports`, {method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({antibody_construct_id: antibody?.id || antibodies[0]?.id,
          experimental_ir: predicted ? [] : [{constraint: `1-${selectedProbe.mapped.normalized_position}-${selectedProbe.mapped.one_letter}`,
            evidence_id: latestEffectId, rationale: exportForm.elements.rationale.value,
            details: {probe_id: selectedProbe.probe.probe_id}}],
          predicted_ir: predicted ? [{constraint: `${chain}-${position}-${sequence[position - 1]}`, evidence_id: run.run_id,
            rationale: exportForm.elements.rationale.value, details: {chain_name: chain === 2 ? "VH" : "VL"}}] : [],
          manual_rpr: [], merge_sources: [predicted ? "predicted_ir" : "experimental_ir"], reviewer: exportForm.elements.reviewer.value,
          expected_revision: project.revision})});
      project.revision = result.project_revision;
      const link = document.createElement("a");
      link.href = `/api/v4/experiment-projects/${encodeURIComponent(projectId)}/files/${encodeURIComponent(result.artifact.id)}`;
      link.textContent = `下载 GRASP ZIP · ${result.manifest.constraint_count} 条 IR · ${result.manifest.status}`;
      exportResult.replaceChildren(link);
    } catch (error) {
      exportResult.className = "error"; exportResult.textContent = `导出失败：${error.message}`;
      exportButton.disabled = false;
    }
  });
  exportSection.append(exportTitle, exportNotice, exportForm, exportPreview, exportResult);
  data.append(exportSection);
  root.append(data);
}

async function renderStructureEffects(
  projectId, project, constructs, antibodies, effectDocument, summaryMeasurements, qcReviews, root,
) {
  const section = document.createElement("section");
  section.id = "structure-effects";
  const title = document.createElement("h3"); title.textContent = "冻结结构上的实验构建体与效应";
  const description = document.createElement("p"); description.className = "muted";
  description.textContent = effectDocument
    ? "结构颜色只呈现已审核的构建体/分区/单点效应。多位点探针保留整体效应，不拆成单残基因果；渐淡层仅供可视化。"
    : "此视图只标出已登记构建体的突变位置；没有 activation 效应分析时不显示功能效应。";
  section.append(title, description);
  const controls = document.createElement("form");
  const layerLabel = document.createElement("label"); layerLabel.textContent = "图层";
  const layer = document.createElement("select");
  layer.innerHTML = '<option value="constructs">构建体突变位点</option>';
  if (effectDocument?.surface_projection) {
    layer.insertAdjacentHTML("beforeend", '<option value="direct">直接观测（整探针）</option><option value="sectors">分区效应</option><option value="sites">单点敏感性</option><option value="fade">显示渐淡（VISUALIZATION_ONLY）</option>');
  }
  layerLabel.append(layer); controls.append(layerLabel);
  const screenshotButton = document.createElement("button");
  screenshotButton.type = "button"; screenshotButton.textContent = "导出当前结构图 PNG";
  screenshotButton.disabled = true; controls.append(screenshotButton);
  const evidenceExportButton = document.createElement("button");
  evidenceExportButton.type = "button"; evidenceExportButton.textContent = "导出结构效应证据 JSON";
  evidenceExportButton.disabled = true; controls.append(evidenceExportButton);
  const freezeEvidenceButton = document.createElement("button");
  freezeEvidenceButton.type = "button";
  freezeEvidenceButton.textContent = "冻结来源的结构可靠性/SASA证据";
  freezeEvidenceButton.hidden = true;
  const host = document.createElement("div"); host.className = "structure-effect-view";
  const detail = document.createElement("p"); detail.className = "structure-effect-detail muted";
  detail.textContent = "读取冻结结构与残基映射…";
  const evidenceSummary = document.createElement("div");
  evidenceSummary.className = "structure-evidence-summary";
  const antibodyTrackSection = document.createElement("section");
  antibodyTrackSection.className = "antibody-mutation-tracks";
  const antibodyTrackTitle = document.createElement("h4"); antibodyTrackTitle.textContent = "VH/VL 抗体突变轨道";
  const antibodyTrackNote = document.createElement("p"); antibodyTrackNote.className = "muted";
  antibodyTrackNote.textContent = "位置按各自 VH/VL 序列编号；抗体轨道独立于抗原三维结构，不投影到抗原残基。";
  antibodyTrackSection.append(antibodyTrackTitle, antibodyTrackNote);
  const variantsWithChanges = (antibodies || []).filter(item =>
    item.role === "VARIANT" && ((item.heavy_mutations || []).length || (item.light_mutations || []).length));
  if (!variantsWithChanges.length) {
    const emptyTracks = document.createElement("p"); emptyTracks.className = "muted";
    emptyTracks.textContent = "尚无带 VH/VL 突变的抗体变体。"; antibodyTrackSection.append(emptyTracks);
  }
  const antibodyById = new Map((antibodies || []).map(item => [item.id, item]));
  for (const antibody of variantsWithChanges) {
    const parent = antibodyById.get(antibody.parent_id);
    const variantTitle = document.createElement("h5");
    variantTitle.textContent = `${antibody.name} · 相对亲本 ${parent?.name || "未知"} · ${antibody.mutation_class}`;
    antibodyTrackSection.append(variantTitle);
    for (const [chain, sequenceKey, mutationsKey] of [
      ["VH", "vh_sequence", "heavy_mutations"], ["VL", "vl_sequence", "light_mutations"],
    ]) {
      const sequence = String(antibody[sequenceKey] || "");
      const mutations = antibody[mutationsKey] || [];
      const row = document.createElement("div"); row.className = "antibody-mutation-track-row";
      const label = document.createElement("strong");
      label.textContent = `${chain} · ${sequence.length} aa · ${mutations.length} substitutions`;
      const scroll = document.createElement("div"); scroll.className = "sequence-track-scroll";
      const lane = document.createElement("div"); lane.className = "sequence-track-lane";
      lane.style.minWidth = `${Math.max(600, sequence.length * 5)}px`;
      for (let position = 10; position <= sequence.length; position += 10) {
        const tick = document.createElement("span"); tick.className = "sequence-track-tick";
        tick.textContent = String(position);
        tick.style.left = `${(position - 1) / Math.max(1, sequence.length - 1) * 100}%`;
        lane.append(tick);
      }
      for (const mutation of mutations) {
        const position = Number(mutation.position);
        if (!Number.isInteger(position) || position < 1 || position > sequence.length) continue;
        const marker = document.createElement("button"); marker.type = "button";
        marker.className = "antibody-track-marker";
        marker.style.left = `${(position - 1) / Math.max(1, sequence.length - 1) * 100}%`;
        marker.title = `${chain}${mutation.position} ${mutation.wt}→${mutation.mut}`;
        marker.setAttribute("aria-label", `${antibody.name} ${chain}${mutation.position} ${mutation.wt} to ${mutation.mut}`);
        marker.addEventListener("click", () => {
          detail.textContent = `${antibody.name} · ${chain} ${mutation.wt}${mutation.position}${mutation.mut} · 相对 ${parent?.name || "未知亲本"}；此位点属于 VH/VL 序列轨道，不是抗原结构位置。`;
        });
        lane.append(marker);
      }
      scroll.append(lane); row.append(label, scroll); antibodyTrackSection.append(row);
    }
  }
  section.append(evidenceSummary, controls, host, antibodyTrackSection, detail);
  root.append(section);
  if (structureEffectViewer?.clear) structureEffectViewer.clear();
  structureEffectViewer = null;
  try {
    if (!window.$3Dmol) throw new Error("3D 结构查看器未载入");
    const context = await request(`/api/v4/experiment-projects/${encodeURIComponent(projectId)}/structure-context`);
    const structureEvidence = context.structure_evidence;
    const appendStatus = (label, value, note) => {
      const card = document.createElement("article");
      const heading = document.createElement("strong"); heading.textContent = label;
      const status = document.createElement("span"); status.textContent = value;
      const explanation = document.createElement("small"); explanation.textContent = note;
      card.append(heading, status, explanation); evidenceSummary.append(card);
    };
    if (structureEvidence?.available) {
      const counts = structureEvidence.reliability.state_counts || {};
      appendStatus(
        "残基可靠性",
        `${structureEvidence.reliability.predicted ? "预测模型" : "非预测结构"} · ${structureEvidence.reliability.records.length} 位点`,
        `${structureEvidence.reliability.classification_version} · ${Object.entries(counts).map(([state, count]) => `${state} ${count}`).join(" · ")} · ${structureEvidence.provenance.residue_reliability.sha256.slice(0, 12)}`,
      );
      const sasaCounts = structureEvidence.sasa.status_counts || {};
      appendStatus(
        "SASA 方法比较",
        `一致 ${sasaCounts.AGREE || 0} · 分歧 ${sasaCounts.SASA_DISAGREEMENT || 0} · 设置未核 ${sasaCounts.METHOD_SETTINGS_UNVERIFIED || 0} · 范围未核 ${sasaCounts.COMPARISON_CONTEXT_UNVERIFIED || 0} · 未比较 ${sasaCounts.NOT_COMPARED || 0}`,
        `ProCaliper 对照同一份规范化 ECD 输入计算；不是功能或表位结论 · ${structureEvidence.provenance.sasa_comparison.sha256.slice(0, 12)}`,
      );
      const sasaMethodCounts = structureEvidence.sasa.method_status_counts || {};
      const sasaContextCounts = structureEvidence.sasa.context_status_counts || {};
      appendStatus(
        "SASA 结构范围影响",
        `稳定 ${sasaContextCounts.CONTEXT_STABLE || 0} · 敏感 ${sasaContextCounts.CONTEXT_SENSITIVE || 0} · 未比较 ${sasaContextCounts.NOT_COMPARED || 0} · 旧记录范围未核 ${sasaContextCounts.LEGACY_SCOPE_UNVERIFIED || 0}`,
        `Biopython 版本相同 ${sasaMethodCounts.MATCHED || 0} · 版本不同 ${sasaMethodCounts.VERSION_DIFFERS || 0} · 设置未核 ${sasaMethodCounts.SETTINGS_DIFFER || 0} · 来源不全 ${sasaMethodCounts.PROVENANCE_INCOMPLETE || 0}。比较原始选定链与规范化 ECD；敏感表示两种输入范围的 rSASA 分类或数值明显变化。`,
      );
    } else {
      appendStatus("残基可靠性 / SASA", "尚未冻结", "仅接受来源作业清单校验并冻结到本项目的证据。");
      freezeEvidenceButton.hidden = false;
      controls.append(freezeEvidenceButton);
      freezeEvidenceButton.addEventListener("click", async () => {
        freezeEvidenceButton.disabled = true;
        try {
          const result = await request(`/api/v4/experiment-projects/${encodeURIComponent(projectId)}/structure-evidence`, {
            method: "POST", headers: {"Content-Type": "application/json"},
            body: JSON.stringify({schema_version: "1.0", expected_revision: project.revision}),
          });
          project.revision = result.project_revision;
          section.remove();
          await renderStructureEffects(
            projectId, project, constructs, antibodies, effectDocument, summaryMeasurements, qcReviews, root,
          );
        } catch (error) {
          freezeEvidenceButton.disabled = false;
          detail.className = "structure-effect-detail error";
          detail.textContent = `冻结结构证据失败：${error.message}`;
        }
      });
    }
    const mycMeasurements = (summaryMeasurements || []).filter(item =>
      item.assay === "ANTIGEN_MYC" || item.assay === "FAB_MYC");
    appendStatus(
      "历史 Myc 汇总",
      mycMeasurements.length ? `${mycMeasurements.length} 条 summary-only 记录` : "未登记",
      mycMeasurements.length
        ? mycMeasurements.map(item => `${item.assay} ${item.metric} ${item.value}${item.unit ? ` ${item.unit}` : ""}`).join(" · ")
        : "独立历史汇总，不与结构可靠性或 SASA 合并。",
    );
    const manualReviewCounts = (qcReviews || []).reduce((counts, item) => {
      counts[item.status] = (counts[item.status] || 0) + 1;
      return counts;
    }, {});
    const reviewLabel = Object.entries(manualReviewCounts).length
      ? Object.entries(manualReviewCounts).map(([status, count]) => `${status} ${count}`).join(" · ")
      : "未记录人工 QC 审核";
    appendStatus(
      "人工 QC / 豁免",
      reviewLabel,
      "审核状态单独呈现；PASS_WITH_OVERRIDE 保留为独立类别。",
    );
    const unmappedConstructs = (constructs || []).filter(item => item.mapping_status !== "MAPPED").length;
    const unmappedProbes = effectDocument?.unmapped_probe_ids?.length || 0;
    appendStatus(
      "UNMAPPED 记录",
      `构建体 ${unmappedConstructs} · 探针 ${unmappedProbes}`,
      "未映射条目保留在数据层，不补到抗原三维表面。",
    );
    if (structureEvidence?.available) {
      layer.insertAdjacentHTML("beforeend", '<option value="reliability">残基可靠性类别</option><option value="sasa">SASA 一致/分歧</option>');
    }
    const structureResponse = await fetch(context.structure.url);
    if (!structureResponse.ok) throw new Error(`结构文件读取失败：HTTP ${structureResponse.status}`);
    const structureText = await structureResponse.text();
    const viewer = window.$3Dmol.createViewer(host, {backgroundColor: "#f4f0e7", antialias: true});
    structureEffectViewer = viewer;
    screenshotButton.addEventListener("click", () => {
      const image = viewer.pngURI();
      if (typeof image !== "string" || !image.startsWith("data:image/png")) {
        detail.textContent = "结构图导出失败：查看器没有返回 PNG 图像。";
        return;
      }
      const link = document.createElement("a");
      link.href = image;
      link.download = `structure_${projectId}_${layer.value}_${context.frozen_hashes.structure.slice(0, 12)}.png`;
      link.click();
    });
    evidenceExportButton.disabled = false;
    evidenceExportButton.addEventListener("click", () => {
      const evidenceDocument = {
        schema_version: "1.0",
        exported_at: new Date().toISOString(),
        project: {
          project_id: projectId,
          revision: project.revision,
          frozen_hashes: context.frozen_hashes,
          structure_layer: layer.value,
        },
        structure_context: {
          mapping: context.mapping,
          surface_graph: context.surface_graph,
          structure_evidence: structureEvidence,
        },
        antigen_constructs: constructs,
        antibody_constructs: antibodies || [],
        effect_analysis: effectDocument,
        historical_myc_summary_only: mycMeasurements,
        manual_qc_reviews: qcReviews || [],
        limitations: [
          "Historical Myc summaries are separate from concurrent activation controls.",
          "Reliability and SASA categories are source evidence, not functional effect estimates.",
          "VH/VL mutation tracks use antibody sequence coordinates and are not projected onto antigen structure.",
        ],
      };
      const blob = new Blob([JSON.stringify(evidenceDocument, null, 2)], {type: "application/json"});
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = `structure_effect_evidence_${projectId}_${layer.value}_${context.frozen_hashes.structure.slice(0, 12)}.json`;
      link.click();
      URL.revokeObjectURL(url);
    });
    viewer.addModel(structureText, context.structure.format === "pdb" ? "pdb" : "cif");
    const byPosition = new Map(context.mapping.map(item => [item.uniprot_pos, item]));
    const graph = context.surface_graph;
    const nodes = graph.nodes || [];
    const sectorByPosition = new Map(nodes.map(item => [item.uniprot_pos, item.sector_id]));
    const reliabilityByPosition = new Map((structureEvidence?.reliability?.records || [])
      .map(item => [item.uniprot_pos, item]));
    const sasaByPosition = new Map((structureEvidence?.sasa?.records || [])
      .map(item => [item.uniprot_pos, item]));
    const effectProbes = effectDocument?.surface_projection?.probes || [];
    const sectorEffects = effectDocument?.surface_projection?.sector_decline_pct || {};
    const singleProbeIds = effectDocument?.surface_projection?.single_site_probe_ids || {};
    const probeById = new Map(effectProbes.map(item => [item.probe_id, item]));
    const mutationsByPosition = new Map();
    for (const construct of constructs) {
      for (const mutation of construct.mutations || []) {
        if (mutation.kind !== "SUBSTITUTION") continue;
        const position = Number(mutation.position);
        if (!mutationsByPosition.has(position)) mutationsByPosition.set(position, []);
        mutationsByPosition.get(position).push({construct: construct.name, ...mutation});
      }
    }
    const adjacency = new Map(nodes.map(item => [item.uniprot_pos, []]));
    for (const edge of graph.edges || []) {
      adjacency.get(edge.source)?.push([edge.target, edge.distance_a]);
      adjacency.get(edge.target)?.push([edge.source, edge.distance_a]);
    }
    const colorForDecline = (value) => {
      const fraction = Math.max(0, Math.min(100, Number(value) || 0)) / 100;
      const base = [190, 202, 210], target = [190, 47, 47];
      return `rgb(${base.map((item, index) => Math.round(item * (1 - fraction) + target[index] * fraction)).join(",")})`;
    };
    const selector = (mapping) => ({chain: mapping.chain_id, resi: mapping.auth_seq_id, icode: mapping.insertion_code || ""});
    const shortestDistances = (start) => {
      const distances = new Map([[start, 0]]), queue = [[0, start]];
      while (queue.length) {
        queue.sort((a, b) => b[0] - a[0]);
        const [distance, current] = queue.pop();
        if (distance !== distances.get(current)) continue;
        for (const [next, step] of adjacency.get(current) || []) {
          const candidate = distance + step;
          if (candidate > 12 || candidate >= (distances.get(next) ?? Infinity)) continue;
          distances.set(next, candidate); queue.push([candidate, next]);
        }
      }
      return distances;
    };
    const paint = () => {
      viewer.setStyle({}, {cartoon: {color: "#7c909d"}});
      const mode = layer.value;
      const colorByPosition = new Map();
      if (mode === "constructs" || mode === "direct") {
        for (const [position, records] of mutationsByPosition) colorByPosition.set(position, "#d78b00");
        if (mode === "direct") {
          for (const probe of effectProbes) {
            const color = probe.attribution === "CROSS_SECTOR_UNASSIGNED" ? "#7d4caf" : "#cf443b";
            for (const position of probe.mutated_positions || []) colorByPosition.set(position, color);
          }
        }
      } else if (mode === "sectors") {
        for (const node of nodes) {
          const decline = sectorEffects[node.sector_id];
          if (decline != null) colorByPosition.set(node.uniprot_pos, colorForDecline(decline));
        }
      } else if (mode === "sites") {
        for (const [positionText, probeIds] of Object.entries(singleProbeIds)) {
          const rows = probeIds.map(id => probeById.get(id)).filter(Boolean);
          const declines = rows.map(item => item.decline_pct).filter(value => value != null);
          if (declines.length) colorByPosition.set(Number(positionText), colorForDecline(Math.max(...declines)));
        }
      } else if (mode === "reliability") {
        const colors = {
          STRUCTURED: "#3f7f68", FLEXIBLE_STRUCTURED: "#d79b2b",
          UNCERTAIN: "#925da8", DISORDERED: "#68747c", UNRESOLVED: "#b3b9bd",
          CONDITIONALLY_FOLDED: "#3f78a3",
        };
        for (const [position, record] of reliabilityByPosition) {
          if (colors[record.state]) colorByPosition.set(position, colors[record.state]);
        }
      } else if (mode === "sasa") {
        for (const [position, record] of sasaByPosition) {
          const color = record.status === "SASA_DISAGREEMENT" ? "#bd4d43"
            : record.native_context_status === "CONTEXT_SENSITIVE" ? "#d79b2b"
            : record.status === "AGREE" ? "#43876b" : "#aeb7bd";
          colorByPosition.set(position, color);
        }
      } else if (mode === "fade") {
        const intensity = new Map();
        for (const probe of effectProbes) {
          if (probe.attribution !== "SINGLE_SITE" || probe.decline_pct == null || probe.decline_pct <= 0) continue;
          const source = probe.mutated_positions?.[0];
          if (!adjacency.has(source)) continue;
          for (const [position, distance] of shortestDistances(source)) {
            const weight = Math.min(1, probe.decline_pct / 100) * Math.exp(-(distance ** 2) / (2 * 5 ** 2));
            intensity.set(position, Math.max(intensity.get(position) || 0, weight));
          }
        }
        for (const [position, value] of intensity) colorByPosition.set(position, colorForDecline(value * 100));
      }
      for (const [position, color] of colorByPosition) {
        const mapping = byPosition.get(Number(position));
        if (mapping) viewer.setStyle(selector(mapping), {cartoon: {color}});
      }
      viewer.setClickable({}, true, (atom) => {
        const mapping = context.mapping.find(item => item.chain_id === atom.chain
          && item.auth_seq_id === Number(atom.resi)
          && (item.insertion_code || "") === (atom.icode || ""));
        if (!mapping) return;
        const mutations = mutationsByPosition.get(mapping.uniprot_pos) || [];
        const probes = effectProbes.filter(item => (item.mutated_positions || []).includes(mapping.uniprot_pos));
        const reliability = reliabilityByPosition.get(mapping.uniprot_pos);
        const sasa = sasaByPosition.get(mapping.uniprot_pos);
        const facts = [`UniProt ${mapping.uniprot_pos} ${mapping.aa || ""}`, `${mapping.chain_id}:${mapping.auth_seq_id}${mapping.insertion_code || ""}`];
        if (reliability) facts.push(`结构可靠性 ${reliability.state}${reliability.plddt == null ? "" : ` · pLDDT ${reliability.plddt.toFixed(1)}`}${reliability.mean_core_pae == null ? "" : ` · core PAE ${reliability.mean_core_pae.toFixed(1)}`}`);
        if (sasa) facts.push(`SASA 方法 ${sasa.status}${sasa.absolute_delta == null ? "" : ` · Δ ${sasa.absolute_delta.toFixed(3)}`} · 结构范围 ${sasa.native_context_status}${sasa.native_context_absolute_delta == null ? "" : ` · Δ ${sasa.native_context_absolute_delta.toFixed(3)}`}`);
        if (mutations.length) facts.push(`构建体 ${mutations.map(item => `${item.construct} ${item.wt}${item.position}${item.mut}`).join("; ")}`);
        if (probes.length) facts.push(`探针 ${probes.map(item => {
          const replacements = (effectDocument?.probe_mutations?.[item.probe_id] || [])
            .filter(mutation => mutation.position === mapping.uniprot_pos)
            .map(mutation => `${mutation.wt}${mutation.position}${mutation.mut}`).join(",");
          return `${replacements || item.attribution}${item.decline_pct == null ? "" : ` · 整探针下降 ${item.decline_pct.toFixed(1)}%`}`;
        }).join("; ")}`);
        detail.textContent = facts.join(" · ");
      });
      viewer.render();
      const modeLabels = {
        constructs: "橙色表示登记构建体的替换位置；不表示功能效应。",
        direct: "红色表示有效应结果的整探针突变位点；紫色表示跨分区组合。组合效应不分摊到单点。",
        sectors: "颜色表示同一分区内不同探针的效应中位数；灰色为无已测分区效应。",
        sites: "仅真正单点探针着色；颜色取同一位置已测替换中的最大下降值以便定位，点击残基查看各替换数据。",
        reliability: "按来源分类状态着色；这是模型可靠性路由信息，不是校准概率。点击残基查看 pLDDT/PAE。",
        sasa: "绿色表示规范化 ECD 输入上的 SASA 结果一致；红色表示方法结果分歧；黄色表示原始链与规范化 ECD 的结构范围敏感；灰色表示未比较或旧记录范围未知。",
        fade: "单点沿冻结表面图最短路渐淡，5 Å 带宽、12 Å 截断；仅供可视化。",
      };
      detail.textContent = modeLabels[mode];
    };
    layer.addEventListener("change", paint);
    viewer.zoomTo(); viewer.render(); paint();
    screenshotButton.disabled = false;
    if (!effectDocument) detail.textContent = "当前没有效应分析；只显示构建体位点，不推断功能或表位效应。";
  } catch (error) {
    detail.className = "structure-effect-detail error";
    detail.textContent = `结构图读取失败：${error.message}`;
  }
  if (window.location.hash === "#structure-effects") section.scrollIntoView({block: "start"});
}

function renderPendingMatches(projectId, imported, root, context = {}) {
  const constructs = context.constructs || [];
  const antibodies = context.antibodies || [];
  const section = document.createElement("section");
  const title = document.createElement("h3"); title.textContent = "确认文件与样本对应关系"; section.append(title);
  const guidance = document.createElement("p"); guidance.className = "muted";
  guidance.textContent = "Activation 样本可关联已登记构建体，或先用明确的临时标签登记；临时标签不会伪造序列，也不能用于结构归属或效应分析。重复类型和编号决定孔/批次汇总方式。未选择角色的文件保持未匹配。";
  section.append(guidance);
  const form = document.createElement("form");
  const rows = [];
  for (const file of imported.files) {
    const row = document.createElement("label");
    row.append(document.createTextNode(`${file.original_filename} · ${file.event_count} events · ${file.sha256.slice(0, 12)} `));
    const role = document.createElement("select");
    role.innerHTML = '<option value="">保留为未匹配</option><option value="ACTIVATION">Activation</option><option value="NC">NC</option><option value="ANTIGEN_MYC">Antigen Myc</option><option value="FAB_MYC">Fab Myc</option><option value="OTHER_CONTROL">其他对照</option>';
    const sample = document.createElement("input"); sample.value = file.original_filename.replace(/\.fcs$/i, ""); sample.setAttribute("aria-label", "样本或构建体临时标签");
    const antigen = document.createElement("select"); antigen.setAttribute("aria-label", "抗原构建体");
    antigen.innerHTML = '<option value="">未关联抗原</option>';
    for (const construct of constructs) {
      const option = document.createElement("option"); option.value = construct.id; option.textContent = construct.name;
      antigen.append(option);
    }
    const antibody = document.createElement("select"); antibody.setAttribute("aria-label", "抗体构建体");
    antibody.innerHTML = '<option value="">未关联抗体</option>';
    for (const construct of antibodies) {
      const option = document.createElement("option"); option.value = construct.id; option.textContent = construct.name;
      antibody.append(option);
    }
    const replicateType = document.createElement("select"); replicateType.setAttribute("aria-label", "重复类型");
    replicateType.innerHTML = '<option value="SINGLE">单孔</option><option value="TECHNICAL">技术重复</option><option value="BIOLOGICAL">生物重复</option>';
    const replicateIndex = document.createElement("input"); replicateIndex.type = "number"; replicateIndex.min = "1";
    replicateIndex.value = "1"; replicateIndex.setAttribute("aria-label", "重复编号");
    const temporaryToggle = document.createElement("input"); temporaryToggle.type = "checkbox";
    temporaryToggle.setAttribute("aria-label", "使用临时身份标签");
    const temporaryLabel = document.createElement("span");
    temporaryLabel.append(temporaryToggle, document.createTextNode(" 暂用标签登记（序列未核实）"));
    const temporaryFields = document.createElement("span"); temporaryFields.style.display = "none";
    const antigenLabel = document.createElement("input"); antigenLabel.placeholder = "抗原标签，如 CTLA4-E1";
    antigenLabel.setAttribute("aria-label", "临时抗原标签");
    const antigenRole = document.createElement("select"); antigenRole.setAttribute("aria-label", "临时抗原类别");
    antigenRole.innerHTML = '<option value="REFERENCE">WT 参考标签</option><option value="VARIANT">变体标签</option><option value="NONE">无抗原（NC）</option>';
    const panelLabel = document.createElement("input"); panelLabel.placeholder = "面板标签，如 L32 或 CD80";
    panelLabel.setAttribute("aria-label", "临时面板标签");
    const panelRole = document.createElement("select"); panelRole.setAttribute("aria-label", "面板试剂类别");
    panelRole.innerHTML = '<option value="ANTIBODY">抗体</option><option value="POSITIVE_LIGAND_CONTROL">配体阳性对照</option>';
    temporaryFields.append(antigenLabel, antigenRole, panelLabel, panelRole);
    temporaryToggle.addEventListener("change", () => {
      temporaryFields.style.display = temporaryToggle.checked ? "inline-flex" : "none";
      antigenLabel.required = temporaryToggle.checked && antigenRole.value !== "NONE";
      panelLabel.required = temporaryToggle.checked;
    });
    antigenRole.addEventListener("change", () => {
      antigenLabel.required = temporaryToggle.checked && antigenRole.value !== "NONE";
      antigenLabel.disabled = temporaryToggle.checked && antigenRole.value === "NONE";
      if (antigenRole.value === "NONE") antigenLabel.value = "";
    });
    row.append(role, sample, antigen, antibody, replicateType, replicateIndex, temporaryLabel, temporaryFields);
    form.append(row); rows.push({file, role, sample, antigen, antibody, replicateType, replicateIndex,
      temporaryToggle, antigenLabel, antigenRole, panelLabel, panelRole});
  }
  const confirm = document.createElement("button"); confirm.textContent = "保存已选择匹配"; form.append(confirm);
  form.addEventListener("submit", async (event) => {
    event.preventDefault();
    const project = await request(`/api/v4/experiment-projects/${encodeURIComponent(projectId)}`);
    const matches = rows.filter(row => row.role.value).map(row => ({
      upload_id: row.file.upload_id, role: row.role.value,
      antigen_construct_id: row.antigen.value || null, antibody_construct_id: row.antibody.value || null,
      replicate_type: row.replicateType.value, replicate_index: Number(row.replicateIndex.value),
      ...(row.temporaryToggle.checked ? {temporary_identity: {
        status: "TEMPORARY_LABEL", source: "USER_CONFIRMED_LABELS",
        antigen_label: row.antigenRole.value === "NONE" ? null : row.antigenLabel.value.trim(),
        antigen_role: row.antigenRole.value, panel_label: row.panelLabel.value.trim(),
        panel_role: row.panelRole.value,
      }} : {}),
      condition: {sample_label: row.sample.value, source_filename: row.file.original_filename}
    }));
    try {
      if (!matches.length) throw new Error("至少选择一个文件角色；未选择的文件会继续保持未匹配");
      const incomplete = matches.find(match =>
        (match.role === "ACTIVATION" && ((!match.antigen_construct_id || !match.antibody_construct_id) && !match.temporary_identity))
        || (match.role === "ANTIGEN_MYC" && !match.antigen_construct_id)
        || (match.role === "FAB_MYC" && !match.antibody_construct_id));
      if (incomplete) throw new Error("Activation 需要构建体关联或完整临时标签；Antigen/Fab Myc 必须关联对应构建体");
      const result = await request(`/api/v4/experiment-projects/${encodeURIComponent(projectId)}/samples/confirm`, {
        method: "POST", headers: {"Content-Type": "application/json"},
        body: JSON.stringify({upload_session_id: imported.upload_session_id, matches, expected_revision: project.revision})
      });
      alert(`已保存 ${result.sample_ids.length} 个样本，${result.unmatched_upload_ids.length} 个文件保留未匹配。`);
      await showProject(projectId);
    } catch (error) { alert(`样本匹配失败：${error.message}`); }
  });
  section.append(form); root.append(section);
}

element("create-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  const button = element("create-button");
  const message = element("create-message");
  button.disabled = true;
  message.className = "";
  message.textContent = "创建并校验中…";
  try {
    const project = await request("/api/v4/experiment-projects", {
      method: "POST", headers: {"Content-Type": "application/json"},
      body: JSON.stringify({schema_version: "1.0", source_job_id: element("source-job").value.trim(), name: element("project-name").value.trim()})
    });
    message.textContent = `已创建：${project.project_id}`;
    await showProjects();
    await showProject(project.project_id);
  } catch (error) { message.className = "error"; message.textContent = `创建失败：${error.message}`; }
  finally { button.disabled = false; }
});

const initialProjectId = new URLSearchParams(window.location.search).get("project_id");
if (initialProjectId) { void showProjects(); void showProject(initialProjectId); }
else showProjects();
