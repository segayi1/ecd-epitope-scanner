/* Browser smoke test for the experiment workspace. Expects a local server. */
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { pathToFileURL } from "node:url";

const server = process.env.ECD_WORKBENCH_URL || "http://127.0.0.1:8765";
const playwrightPath = process.env.PLAYWRIGHT_CORE ||
  "C:/Users/陈淏文/.workbuddy/binaries/node/workspace/node_modules/playwright-core/index.js";
const chromePath = process.env.CHROME_PATH ||
  "C:/Program Files/Google/Chrome/Application/chrome.exe";
if (!fs.existsSync(playwrightPath)) throw new Error(`Playwright module not found: ${playwrightPath}`);
if (!fs.existsSync(chromePath)) throw new Error(`Chrome executable not found: ${chromePath}`);
const playwright = await import(pathToFileURL(playwrightPath).href);
const browserType = playwright.chromium || playwright.default?.chromium;
const browser = await browserType.launch({executablePath: chromePath, headless: true});
const page = await browser.newPage({viewport: {width: 1440, height: 1000}});
const failures = [];
const check = (condition, label) => {
  console.log(`${condition ? "PASS" : "FAIL"} ${label}`);
  if (!condition) failures.push(label);
};
const pageErrors = [];
page.on("pageerror", error => pageErrors.push(error.message));
page.on("console", message => {
  if (message.type() === "error") pageErrors.push(message.text());
});

try {
  const projectsResponse = await page.request.get(`${server}/api/v4/experiment-projects`);
  check(projectsResponse.ok(), "project list API responds successfully");
  const projects = await projectsResponse.json();
  const project = projects.data.find(item => item.project_id === process.env.ECD_PROJECT_ID) || projects.data[0];
  check(Boolean(project), "at least one experiment project is available");
  if (!project) throw new Error("No experiment project is available for browser verification");

  await page.route(`**/api/v4/experiment-projects/${project.project_id}/structure-context`, async route => {
    const response = await route.fetch();
    const body = await response.json();
    const positions = body.data.mapping.slice(0, 3).map((item, index) => ({
      uniprot_pos: item.uniprot_pos,
      state: ["STRUCTURED", "UNCERTAIN", "DISORDERED"][index],
      plddt: [82, 57, 41][index], mean_local_pae: 2.1, mean_core_pae: 8.5,
      quality_flags: [], use_for_surface_mesh: index !== 2,
      use_for_3d_footprint: index === 0, use_for_epitope_screening: true,
    }));
    body.data.structure_evidence = {
      available: true,
      provenance: {
        residue_reliability: {id: "browser-reliability", sha256: "a".repeat(64)},
        sasa_comparison: {id: "browser-sasa", sha256: "b".repeat(64)},
        reliability_qc: {id: "browser-qc", sha256: "c".repeat(64)},
      },
      reliability: {
        predicted: true, classification_version: "reliability_rules_v1",
        state_counts: {STRUCTURED: 1, UNCERTAIN: 1, DISORDERED: 1}, records: positions,
      },
      sasa: {
        records: positions.map((item, index) => ({
          uniprot_pos: item.uniprot_pos, native_rsasa: 0.5,
          normalized_context_rsasa: 0.2, procaliper_rsasa: 0.2,
          absolute_delta: index === 1 ? 0.12 : 0.01,
          native_context_absolute_delta: index === 0 ? 0.3 : 0.01,
          native_context_status: index === 0 ? "CONTEXT_SENSITIVE" : "CONTEXT_STABLE",
          native_context_classification_changed: index === 0,
          procaliper_classification_changed: index === 1,
          status: index === 1 ? "SASA_DISAGREEMENT" : "AGREE",
        })),
        status_counts: {AGREE: 2, SASA_DISAGREEMENT: 1, NOT_COMPARED: 0},
        context_status_counts: {CONTEXT_STABLE: 2, CONTEXT_SENSITIVE: 1, NOT_COMPARED: 0},
      },
    };
    await route.fulfill({response, json: body});
  });
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/antibody-constructs`, route =>
    route.fulfill({json: {data: [
      {id: "browser-parent-antibody", name: "Track parent", role: "PARENT", parent_id: null,
        mutation_class: "PARENT", heavy_mutations: [], light_mutations: [],
        vh_sequence: "ACDEFG", vl_sequence: "HIKLMN"},
      {id: "browser-variant-antibody", name: "Track variant", role: "VARIANT",
        parent_id: "browser-parent-antibody", mutation_class: "SINGLE",
        heavy_mutations: [{position: 3, wt: "D", mut: "V"}], light_mutations: [],
        vh_sequence: "ACVEFG", vl_sequence: "HIKLMN"},
    ]}}));

  await page.goto(`${server}/experiments`, {waitUntil: "domcontentloaded"});
  await page.waitForSelector("#project-list article button", {timeout: 15000});
  const projectButton = page.locator("#project-list article button").filter({hasText: project.project_id}).first();
  await projectButton.click();
  await page.waitForFunction(() => {
    const detail = document.querySelector("#project-content");
    return detail?.textContent.includes("GRASP IR 导出");
  }, {timeout: 15000});

  check(await page.locator("#project-detail").isVisible(), "project detail panel is visible");
  check(await page.locator("#project-content").getByText("GRASP IR 导出").isVisible(), "GRASP export panel is rendered");
  check(await page.locator("#structure-effects").getByText("残基可靠性", {exact: true}).count() === 1,
    "residue reliability summary is rendered separately");
  check(await page.locator("#structure-effects").getByText("SASA 方法比较", {exact: true}).count() === 1,
    "SASA comparison summary is rendered separately");
  check(await page.locator("#structure-effects").getByText("SASA 结构范围影响", {exact: true}).count() === 1,
    "source-context sensitivity is summarized separately from method agreement");
  check(await page.locator("#structure-effects").getByText("历史 Myc 汇总", {exact: true}).count() === 1,
    "historical Myc summary remains a separate status");
  check(await page.locator("#structure-effects").getByText("人工 QC / 豁免", {exact: true}).count() === 1,
    "manual QC review status remains separate");
  check(await page.locator("#structure-effects").getByRole("heading", {name: "VH/VL 抗体突变轨道"}).count() === 1,
    "antibody mutations are displayed on independent VH/VL sequence tracks");
  const antibodyMarker = page.locator("#structure-effects .antibody-track-marker").first();
  check(await antibodyMarker.count() === 1, "variant substitution is marked on its sequence track");
  await antibodyMarker.click();
  check((await page.locator(".structure-effect-detail").innerText()).includes("此位点属于 VH/VL 序列轨道"),
    "antibody track click does not map sequence positions onto the antigen structure");
  const structureLayer = page.locator("#structure-effects select").first();
  const structureLayerOptions = await structureLayer.locator("option").evaluateAll(options => options.map(item => item.value));
  check(structureLayerOptions.includes("reliability") && structureLayerOptions.includes("sasa"),
    "categorical reliability and SASA structure layers are available when evidence is frozen");
  await structureLayer.selectOption("reliability");
  check((await page.locator(".structure-effect-detail").innerText()).includes("按来源分类状态着色"),
    "reliability layer explains its categorical interpretation");
  await structureLayer.selectOption("sasa");
  check((await page.locator(".structure-effect-detail").innerText()).includes("结构范围敏感"),
    "SASA layer distinguishes context-sensitive residues from algorithm disagreement");

  const freezePage = await browser.newPage({viewport: {width: 1280, height: 900}});
  let freezeRequestRevision = null;
  let mockedEvidenceFrozen = false;
  let freezeContextRequests = 0;
  await freezePage.route(`**/api/v4/experiment-projects/${project.project_id}/structure-context`, async route => {
    freezeContextRequests += 1;
    const response = await route.fetch();
    const body = await response.json();
    body.data.structure_evidence = mockedEvidenceFrozen ? {
      available: true,
      provenance: {
        residue_reliability: {id: "freeze-reliability", sha256: "d".repeat(64)},
        sasa_comparison: {id: "freeze-sasa", sha256: "e".repeat(64)},
        reliability_qc: {id: "freeze-qc", sha256: "f".repeat(64)},
      },
      reliability: {predicted: true, classification_version: "reliability_rules_v1",
        state_counts: {STRUCTURED: body.data.mapping.length}, records: []},
      sasa: {records: [], status_counts: {AGREE: body.data.mapping.length},
        context_status_counts: {CONTEXT_STABLE: body.data.mapping.length}},
    } : null;
    await route.fulfill({response, json: body});
  });
  await freezePage.route(`**/api/v4/experiment-projects/${project.project_id}/structure-evidence`, async route => {
    freezeRequestRevision = route.request().postDataJSON().expected_revision;
    mockedEvidenceFrozen = true;
    await route.fulfill({json: {data: {project_revision: freezeRequestRevision + 1, already_frozen: false}}});
  });
  await freezePage.goto(`${server}/experiments`, {waitUntil: "domcontentloaded"});
  await freezePage.locator("#project-list article button").filter({hasText: project.project_id}).first().click();
  const freezeEvidenceButton = freezePage.getByRole("button", {name: "冻结来源的结构可靠性/SASA证据"});
  try {
    await freezeEvidenceButton.waitFor({state: "visible", timeout: 15000});
  } catch (error) {
    console.error(`Structure-context mock requests: ${freezeContextRequests}`);
    console.error((await freezePage.locator("#project-content").innerText()).slice(-2500));
    throw error;
  }
  await freezeEvidenceButton.click();
  await freezePage.waitForFunction(() => [...document.querySelectorAll("#structure-effects select option")]
    .some(option => option.value === "reliability"), {timeout: 15000});
  check(freezeRequestRevision === project.revision,
    "evidence freeze UI submits the current project revision and refreshes after success");
  await freezePage.close();

  const structurePngButton = page.getByRole("button", {name: "导出当前结构图 PNG"});
  await page.waitForFunction(() => {
    const button = [...document.querySelectorAll("#structure-effects button")]
      .find(item => item.textContent.includes("导出当前结构图 PNG"));
    return button && !button.disabled;
  }, {timeout: 20000});
  const structurePngDownloadPromise = page.waitForEvent("download");
  await structurePngButton.click();
  const structurePngDownload = await structurePngDownloadPromise;
  const structurePngPath = path.join(os.tmpdir(), `ecd-structure-${project.project_id}.png`);
  await structurePngDownload.saveAs(structurePngPath);
  check(structurePngDownload.suggestedFilename().endsWith(".png")
    && fs.statSync(structurePngPath).size > 1000,
  "current frozen structure layer can be exported as a PNG");
  fs.unlinkSync(structurePngPath);
  const evidenceDownloadPromise = page.waitForEvent("download");
  await page.getByRole("button", {name: "导出结构效应证据 JSON"}).click();
  const evidenceDownload = await evidenceDownloadPromise;
  const evidencePath = path.join(os.tmpdir(), `ecd-evidence-${project.project_id}.json`);
  await evidenceDownload.saveAs(evidencePath);
  const exportedEvidence = JSON.parse(fs.readFileSync(evidencePath, "utf8"));
  check(evidenceDownload.suggestedFilename().endsWith(".json")
    && exportedEvidence.project.project_id === project.project_id
    && exportedEvidence.structure_context.structure_evidence.available
    && exportedEvidence.limitations.some(item => item.includes("VH/VL")),
  "structure effect evidence JSON includes frozen provenance and interpretation limits");
  fs.unlinkSync(evidencePath);
  const supportLink = page.getByRole("link", {name: /下载项目支持包/});
  const supportCheckbox = page.locator("#support-include-originals");
  check(await supportLink.count() === 1, "project support bundle link is rendered");
  check(!(await supportLink.getAttribute("href")).includes("include_originals"), "support bundle excludes raw FCS/WSP by default");
  await supportCheckbox.check();
  check((await supportLink.getAttribute("href")).includes("include_originals=true"), "support bundle can explicitly include raw FCS/WSP");
  await supportCheckbox.uncheck();
  const projectStatus = await page.locator("#project-content").innerText();
  check(projectStatus.includes("COMPLETE"), "project integrity status is visible");
  check(await page.getByText("历史 Myc 汇总（summary-only）").count() === 1,
    "historical Myc summary-only panel is rendered");
  check(await page.getByRole("heading", {name: "事件分布图（抽样显示）"}).count() === 1,
    "display-only event plot panel is rendered");
  const sampleResponse = await page.request.get(`${server}/api/v4/experiment-projects/${project.project_id}/samples`);
  const sampleRows = (await sampleResponse.json()).data;
  const plotSample = sampleRows.find(item => item.artifact_id && item.fcs_metadata?.channels?.length);
  if (plotSample) {
    const plotForm = page.locator("#project-content form").filter({has: page.getByRole("button", {name: "读取并绘图"})});
    await plotForm.locator('select[name="sample"]').selectOption(plotSample.id);
    const availableChannels = await plotForm.locator('select[name="x"] option').evaluateAll(options => options.map(option => option.value));
    const xChannel = availableChannels.includes("FSC-A") ? "FSC-A" : availableChannels[0];
    const yChannel = availableChannels.includes("SSC-A") ? "SSC-A" : availableChannels.find(value => value !== xChannel);
    await plotForm.locator('select[name="x"]').selectOption(xChannel);
    await plotForm.locator('select[name="y"]').selectOption(yChannel);
    await plotForm.getByRole("button", {name: "读取并绘图"}).click();
    await page.waitForFunction(() => document.querySelector("#sample-plot-status")?.textContent.includes("完整事件"), {timeout: 20000});
    check((await page.locator("#sample-plot-status").innerText()).includes("抽样"),
      "event plot reports full and sampled event counts separately");
    const plotCanvas = page.locator('canvas[aria-label="FCS 事件散点图或密度图"]');
    check(await plotCanvas.evaluate(canvas => canvas.width === 900 && canvas.height === 560),
      "sampled scatter plot renders to its canvas");
    await plotForm.locator('select[name="mode"]').selectOption("density");
    check(await plotCanvas.evaluate(canvas => canvas.getContext("2d") !== null),
      "density view redraws from the same display sample");
  } else {
    check(false, "project has a matched FCS sample for event plot verification");
  }
  check(await page.getByRole("button", {name: "预览与 WT 差异"}).count() === 1,
    "full-sequence antigen preview action is rendered");
  const exportButton = page.getByRole("button", {name: "校验并生成 ZIP"});
  check(await exportButton.count() === 1, "single export action is rendered");

  const antibodies = await page.request.get(`${server}/api/v4/experiment-projects/${project.project_id}/antibody-constructs`);
  const runs = await page.request.get(`${server}/api/v4/experiment-projects/${project.project_id}/paraplume-runs`);
  const antibodyRows = (await antibodies.json()).data;
  const runRows = (await runs.json()).data;
  if (!antibodyRows.length && !runRows.length) {
    check(await exportButton.isDisabled(), "export is disabled when the project has no reviewed IR candidates");
    check(projectStatus.includes("尚无可导出的已验证实验单点或 Paraplume 结果"), "empty export state explains why export is unavailable");
  }

  const testAntibodyId = "browser-fixture-antibody";
  const testEffectId = "browser-fixture-effect";
  const testEffectArtifactId = "browser-fixture-effect-document";
  const testGateAnalysisId = "browser-fixture-gate-analysis";
  const testGateAnalysisArtifactId = "browser-fixture-gate-document";
  let showExcludedGateAnalysis = false;
  let submittedWindowExclusion = null;
  let showFactorialResult = false;
  let showAntibodyResult = false;
  const projectResponse = await page.request.get(`${server}/api/v4/experiment-projects/${project.project_id}`);
  const projectDetails = (await projectResponse.json()).data;
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/antigen-constructs`, async route => {
    const response = await route.fetch();
    const body = await response.json();
    const rows = body.data || [];
    const reference = projectDetails.manifest.reference_sequence;
    if (!rows.some(item => item.sequence === reference)) {
      rows.push({id: "browser-wt-antigen", name: "Browser WT antigen", source: "FIXTURE",
        sequence: reference, sequence_sha256: "w".repeat(64), mapping_status: "MAPPED",
        mutations: [], risk_flags: []});
    }
    if (!rows.some(item => item.id === "browser-mutant-antigen")) {
      const alternate = [..."ACDEFGHIKLMNPQRSTVWY"].find(letter => letter !== reference[0]);
      rows.push({id: "browser-mutant-antigen", name: "Browser mutant antigen", source: "FIXTURE",
        sequence: `${alternate}${reference.slice(1)}`,
        sequence_sha256: "m".repeat(64), mapping_status: "MAPPED",
        mutations: [{kind: "SUBSTITUTION", position: 1, wt: reference[0], mut: alternate, mapped: true}],
        risk_flags: []});
    }
    await route.fulfill({response, json: {...body, data: rows}});
  });
  const testGateSample = sampleRows.find(item => item.artifact_id);
  const testGateFcs = projectDetails.artifacts.find(item => item.id === testGateSample?.artifact_id);
  const fullSequenceForm = page.locator("#project-content form").filter({
    has: page.getByRole("button", {name: "预览与 WT 差异"}),
  });
  await fullSequenceForm.locator('input[name="name"]').fill("Browser WT preview");
  await fullSequenceForm.locator('textarea[name="sequence"]').fill(projectDetails.manifest.reference_sequence);
  await fullSequenceForm.getByRole("button", {name: "预览与 WT 差异"}).click();
  await page.waitForFunction(() => {
    const form = document.querySelector('textarea[name="sequence"]')?.closest("form");
    return form?.nextElementSibling?.textContent.includes("差异 0 项");
  }, {timeout: 10000});
  check((await fullSequenceForm.locator("xpath=following-sibling::pre[1]").innerText()).includes("差异 0 项"),
    "full-sequence preview compares against the frozen WT without writing");
  check(!(await fullSequenceForm.getByRole("button", {name: "确认并登记完整序列"}).isDisabled()),
    "full-sequence save is enabled after a successful preview");
  const mappingArtifact = projectDetails.artifacts.find(item => item.kind === "mapping");
  const mappingTextResponse = await page.request.get(`${server}/api/v4/experiment-projects/${project.project_id}/files/${mappingArtifact.id}`);
  const mappingLines = (await mappingTextResponse.text()).trim().split(/\r?\n/);
  const mappingColumns = mappingLines[0].split(",");
  const firstMappedRow = mappingLines.slice(1).map(line => {
    const values = line.split(",");
    return Object.fromEntries(mappingColumns.map((column, index) => [column, values[index] || ""]));
  }).find(row => row.uniprot_pos && row.normalized_position && row.one_letter);
  if (!firstMappedRow) throw new Error("No mapped reference residue is available for the experimental IR fixture");
  const testProbeId = "browser-fixture-single-site";
  let submittedExport = null;
  const gateAnalysisSummary = {
    id: testGateAnalysisId, result_artifact_id: testGateAnalysisArtifactId, status: "COMPLETE",
    gate_revision_id: "browser-fixture-gate-revision", input_artifact_ids: [testGateFcs?.id],
  };
  const gateAnalysisDocument = {
    analysis_id: testGateAnalysisId, gate_revision_id: gateAnalysisSummary.gate_revision_id,
    input_artifact_ids: [testGateFcs?.id], qc_policy: {time_window_count: 2},
    samples: [{artifact_id: testGateFcs?.id, input_sha256: testGateFcs?.sha256,
      sample: {sample_id: testGateSample?.id, role: testGateSample?.role || "NC",
        condition: {sample_label: "browser-H1"}},
      event_count: 100, source_event_count: 100, events_removed: 0,
      gate_counts: [{gate_id: "browser-positive", gate_name: "Browser Positive",
        parent_count: 100, count: 60, parent_fraction: 0.6}],
      qc: {status: "REVIEW_REQUIRED", events_examined: 100, events_removed: 0,
        findings: [{code: "TIME_EVENT_RATE_OUTLIER", severity: "REVIEW", windows: [{event_start: 0, event_end: 50}]}],
        time_windows: [
          {event_start: 0, event_end: 50, time_span: 10, event_rate: 4.9, rate_outlier: true},
          {event_start: 50, event_end: 100, time_span: 50, event_rate: 1, rate_outlier: false},
        ]},
    }],
  };
  const excludedGateAnalysisDocument = {
    ...gateAnalysisDocument,
    analysis_id: "browser-fixture-excluded-gate-analysis",
    event_exclusions: {source_gate_analysis_id: testGateAnalysisId, reviewer: "browser reviewer",
      reason: "browser time-window review", windows: [{sample_id: testGateSample?.id,
        ranges: [{event_start: 0, event_end: 50}], events_removed: 50}]},
    samples: gateAnalysisDocument.samples.map(item => ({...item, event_count: 50,
      source_event_count: 100, events_removed: 50, qc: {...item.qc, events_examined: 50, events_removed: 50,
        time_windows: [{event_start: 0, event_end: 50, time_span: 50, event_rate: 1, rate_outlier: false}]}})),
  };
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/gate-analyses`, async route => {
    if (route.request().method() === "POST") {
      submittedWindowExclusion = route.request().postDataJSON();
      showExcludedGateAnalysis = true;
      await route.fulfill({status: 201, json: {data: {samples: [{events_removed: 50}]}}});
      return;
    }
    await route.fulfill({json: {data: [showExcludedGateAnalysis
      ? {...gateAnalysisSummary, id: excludedGateAnalysisDocument.analysis_id}
      : gateAnalysisSummary]}});
  });
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/files/${testGateAnalysisArtifactId}`, route =>
    route.fulfill({json: showExcludedGateAnalysis ? excludedGateAnalysisDocument : gateAnalysisDocument}));
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/antibody-constructs`, route =>
    route.fulfill({json: {data: [{
      id: testAntibodyId, name: "Browser fixture", role: "PARENT", mutation_class: "PARENT",
      heavy_mutations: [], light_mutations: [], vh_sequence: "ACDE", vl_sequence: "FGHI",
    }]}}));
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/paraplume-runs`, route =>
    route.fulfill({json: {data: [{
      run_id: "browser-fixture-run", antibody_construct_id: testAntibodyId,
      prediction_status: "PREDICTED_UNREVIEWED", vh_probabilities: [0.1, 0.2, 0.3, 0.4],
      vl_probabilities: [0.1, 0.7, 0.2, 0.3], provenance: {actual_device: "fixture"},
    }]}}));
  let submittedParaplume = null;
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/paraplume`, async route => {
    submittedParaplume = route.request().postDataJSON();
    await route.fulfill({status: 202, json: {data: {
      task_id: "browser-progress-task",
      events_url: "/api/v4/experiment-tasks/browser-progress-task/events",
    }}});
  });
  await page.route("**/api/v4/experiment-tasks/browser-progress-task/events", route => {
    const events = [
      {task_id: "browser-progress-task", status: "RUNNING", progress: 0.14,
        message: "Loading the full Paraplume model", output: null},
      {task_id: "browser-progress-task", status: "RUNNING", progress: 0.30,
        message: "Computing IgBert residue embeddings", output: null},
      {task_id: "browser-progress-task", status: "RUNNING", progress: 0.86,
        message: "Completed ProtT5 residue embeddings", output: null},
      {task_id: "browser-progress-task", status: "COMPLETE", progress: 1,
        message: null, output: {record_count: 1}},
    ];
    return route.fulfill({status: 200,
      headers: {"Content-Type": "text/event-stream; charset=utf-8", "Cache-Control": "no-cache"},
      body: events.map(item => `data: ${JSON.stringify(item)}\n\n`).join("")});
  });
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/effect-analyses`, route =>
    route.fulfill({json: {data: [{id: testEffectId, result_artifact_id: testEffectArtifactId, status: "COMPLETE"}]}}));
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/files/${testEffectArtifactId}`, route =>
    route.fulfill({json: showAntibodyResult ? {
      comparison_mode: "ANTIBODY_EFFECT", batch_effects: [], surface_projection: null,
      antibody_effects: [{probe_id: "browser-wt-antigen", batch_id: "browser-batch-1",
        parent_antibody_construct_id: "browser-parent-antibody",
        variant_antibody_construct_id: "browser-variant-antibody",
        parent_wt_fraction: 0.8, parent_wt_nc_fraction: 0.1, parent_control_window_pp: 70,
        variant_wt_fraction: 0.6, variant_wt_nc_fraction: 0.1, variant_control_window_pp: 50,
        variant_minus_parent_pp: -20, valid: true, qc_flags: []}],
      antibody_summaries: [{probe_id: "browser-wt-antigen", batch_count: 1, valid_batch_count: 1,
        median_variant_minus_parent_pp: -20, interval_pp: null,
        positive_batch_count: 0, negative_batch_count: 1, zero_batch_count: 0, qc_flags: []}],
      well_results: [], probe_mutations: {}, unmapped_probe_ids: [],
    } : showFactorialResult ? {
      comparison_mode: "DOUBLE_MUTANT_INTERACTION",
      batch_effects: [
        {probe_id: testProbeId, antibody_construct_id: "browser-parent-antibody",
          batch_effect: {batch_id: "browser-batch-1", mutant_fraction: 0.4,
            reference_fraction: 0.7, mutant_nc_fraction: 0.1, reference_nc_fraction: 0.1,
            mutant_control_window_pp: 30, reference_control_window_pp: 60,
            difference_pp: -30, retention: 0.5, decline_pct: 50, qc_flags: []}},
        {probe_id: testProbeId, antibody_construct_id: "browser-variant-antibody",
          batch_effect: {batch_id: "browser-batch-1", mutant_fraction: 0.6,
            reference_fraction: 0.8, mutant_nc_fraction: 0.1, reference_nc_fraction: 0.1,
            mutant_control_window_pp: 50, reference_control_window_pp: 70,
            difference_pp: -20, retention: 0.5, decline_pct: 50, qc_flags: []}},
      ],
      factorial_interactions: [{probe_id: testProbeId, batch_id: "browser-batch-1", interaction_pp: 10,
        four_cells: {
          parent_antibody: {mutant_antigen: 0.4, wt_antigen: 0.7},
          variant_antibody: {mutant_antigen: 0.6, wt_antigen: 0.8},
        }, qc_flags: []}],
      factorial_summaries: [{probe_id: testProbeId, batch_count: 1, median_interaction_pp: 10,
        interaction_interval_pp: null, positive_batch_count: 1, negative_batch_count: 0, zero_batch_count: 0}],
      surface_projection: null, unmapped_probe_ids: [], well_results: [], probe_mutations: {},
    } : {
      batch_effects: [
        {probe_id: testProbeId, batch_effect: {batch_id: "browser-batch-1", mutant_fraction: 0.4,
          reference_fraction: 0.8, mutant_nc_fraction: 0.1, reference_nc_fraction: 0.1,
          mutant_control_window_pp: 30, reference_control_window_pp: 70,
          difference_pp: -40, retention: 0.5, decline_pct: 50, qc_flags: []}},
        {probe_id: testProbeId, batch_effect: {batch_id: "browser-batch-2", mutant_fraction: 0.3,
          reference_fraction: 0.8, mutant_nc_fraction: 0.1, reference_nc_fraction: 0.1,
          mutant_control_window_pp: 20, reference_control_window_pp: 70,
          difference_pp: -50, retention: null, decline_pct: null,
          qc_flags: ["INSUFFICIENT_REFERENCE_WINDOW"]}},
      ], well_results: ["MUTANT", "REFERENCE", "MUTANT_NC", "REFERENCE_NC"].flatMap((condition, index) => {
        const fraction = [0.4, 0.8, 0.1, 0.1][index];
        return [
          {sample_id: `tech-${condition}-1`, probe_id: testProbeId, batch_id: "browser-batch-1",
            antibody_construct_id: testAntibodyId, condition, fraction, positive_count: fraction * 100,
            parent_count: 100, qc_status: "PASS", qc_reason: "", replicate_type: "TECHNICAL", replicate_index: 1},
          {sample_id: `tech-${condition}-2`, probe_id: testProbeId, batch_id: "browser-batch-1",
            antibody_construct_id: testAntibodyId, condition, fraction: fraction + 0.02, positive_count: fraction * 100 + 2,
            parent_count: 100, qc_status: "PASS_WITH_OVERRIDE", qc_reason: "fixture override",
            replicate_type: "TECHNICAL", replicate_index: 2},
          {sample_id: `held-${condition}-3`, probe_id: testProbeId, batch_id: "browser-batch-1",
            antibody_construct_id: testAntibodyId, condition, fraction: fraction + 0.05, positive_count: fraction * 100 + 5,
            parent_count: 100, qc_status: "HOLD", qc_reason: "held fixture well",
            replicate_type: "TECHNICAL", replicate_index: 3},
          {sample_id: `bio-${condition}-1`, probe_id: testProbeId, batch_id: "browser-batch-2",
            antibody_construct_id: testAntibodyId, condition, fraction, positive_count: fraction * 100,
            parent_count: 100, qc_status: "PASS", qc_reason: "", replicate_type: "BIOLOGICAL", replicate_index: 1},
        ];
      }), probe_mutations: {},
      surface_projection: {
        probes: [{probe_id: testProbeId, attribution: "SINGLE_SITE",
          mutated_positions: [Number(firstMappedRow.uniprot_pos)], decline_pct: 31.4,
          decline_interval_pct: [25.2, 38.1], difference_pp: -42.5,
          positive_decline_batch_count: 2, negative_decline_batch_count: 0, zero_decline_batch_count: 0,
          valid_decline_batch_count: 2, batch_count: 2}],
        sector_decline_pct: {}, single_site_probe_ids: {
          [firstMappedRow.uniprot_pos]: [testProbeId],
        },
      },
    }}));
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/constraint-exports`, async route => {
    submittedExport = route.request().postDataJSON();
    await route.fulfill({status: 201, json: {data: {
      project_revision: project.revision + 1,
      artifact: {id: "browser-fixture-export"},
      manifest: {constraint_count: 1, status: "IR_ONLY"},
    }}});
  });
  await page.reload({waitUntil: "domcontentloaded"});
  await page.waitForSelector("#project-list article button", {timeout: 15000});
  await page.locator("#project-list article button").filter({hasText: project.project_id}).first().click();
  await page.waitForFunction(() => document.querySelector("#project-content")?.textContent.includes("GRASP IR 导出"), {timeout: 15000});
  const windowForm = page.locator("#project-content form").filter({
    has: page.getByRole("button", {name: "排除所选窗口并重新分析"}),
  });
  check(await windowForm.count() === 1
    && (await windowForm.innerText()).includes("自动标记异常"),
  "reviewer can choose a complete acquisition time window flagged by automatic QC");
  await windowForm.locator('input[type="checkbox"]').first().check();
  await windowForm.locator('input[name="reviewer"]').fill("browser reviewer");
  await windowForm.locator('input[name="reason"]').fill("browser time-window review");
  await windowForm.getByRole("button", {name: "排除所选窗口并重新分析"}).click();
  await page.waitForFunction(() => document.querySelector("#project-content")?.textContent.includes("本分析继承已审核时间窗排除"), {timeout: 15000});
  check(submittedWindowExclusion?.source_gate_analysis_id === testGateAnalysisId
    && submittedWindowExclusion?.exclusion_reviewer === "browser reviewer"
    && submittedWindowExclusion?.exclusion_reason === "browser time-window review"
    && submittedWindowExclusion?.window_exclusions?.[0]?.windows?.[0]?.event_start === 0,
  "window exclusion creates a reviewer-attributed reanalysis and shows its inherited audit record");
  const effectHeatmap = page.locator('table[aria-label="逐探针独立批次效应热图"]');
  check(await effectHeatmap.count() === 1
    && await effectHeatmap.locator("tr").count() === 2
    && (await effectHeatmap.innerText()).includes("无效"),
  "effect heatmap separates independent batches and keeps invalid windows explicit");
  const controlDashboard = page.locator('table[aria-label="每批四条件对照窗口"]');
  check(await controlDashboard.count() === 1 && await controlDashboard.locator("tr").count() === 3
    && (await controlDashboard.innerText()).includes("突变抗原 NC")
    && (await controlDashboard.innerText()).includes("70.00 pp"),
  "batch control dashboard shows all four rates and NC-adjusted windows, including invalid windows");
  const batchSummary = page.locator('table[aria-label="跨批次效应汇总"]');
  check(await batchSummary.count() === 1 && await batchSummary.locator("tr").count() === 2
    && (await batchSummary.innerText()).includes("不足") === false
    && (await batchSummary.innerText()).includes("25.20% 至 38.10%")
    && (await batchSummary.innerText()).includes("2 / 0 / 0"),
  "probe summaries show independent-batch median, bootstrap range, and direction counts");
  check(await page.locator('svg[aria-label="逐孔阳性比例图"]').count() === 2
    && await page.locator('svg[aria-label="逐孔阳性比例图"] circle').count() === 16,
  "per-well charts keep samples in separate batch and antibody panels");
  const wellReadoutText = await page.locator("#project-content").innerText();
  check(await page.locator("#project-content").getByText("技术重复孔离散度", {exact: true}).count() === 1
    && wellReadoutText.includes("纳入技术孔") && wellReadoutText.includes("41.00%")
    && wellReadoutText.includes("40.00–42.00%") && wellReadoutText.includes("PASS/PASS_WITH_OVERRIDE"),
  "technical summaries use accepted wells only and show median/range separately from biological batches");
  const exportForm = page.locator("#project-content form").filter({has: page.locator('input[name="rationale"]')});
  await exportForm.locator('select[name="antibody"]').selectOption(testAntibodyId);
  await exportForm.locator('select[name="chain"]').selectOption("3");
  await exportForm.locator('input[name="position"]').fill("2");
  const preview = await page.locator("#project-content pre").filter({hasText: "PREDICTED_UNREVIEWED"}).innerText();
  check(preview.includes("3-2-G"), "predicted IR preview uses selected chain, local residue, and FASTA amino acid");
  await exportForm.locator('input[name="reviewer"]').fill("browser reviewer");
  await exportForm.locator('input[name="rationale"]').fill("browser fixture review");
  await exportForm.getByRole("button", {name: "校验并生成 ZIP"}).click();
  await page.getByRole("link", {name: /下载 GRASP ZIP/}).waitFor({timeout: 10000});
  check(submittedExport?.predicted_ir?.[0]?.constraint === "3-2-G", "UI submits the previewed IR constraint");
  check(submittedExport?.predicted_ir?.[0]?.evidence_id === "browser-fixture-run", "UI submits the Paraplume run provenance");

  await exportForm.locator('select[name="source"]').selectOption("experimental_ir");
  const experimentalPreview = await page.locator("#project-content pre").filter({hasText: "SINGLE_SITE"}).innerText();
  const expectedExperimentalIr = `1-${firstMappedRow.normalized_position}-${firstMappedRow.one_letter}`;
  check(experimentalPreview.includes(expectedExperimentalIr), "experimental IR preview uses the frozen UniProt-to-ECD mapping");
  await exportForm.getByRole("button", {name: "校验并生成 ZIP"}).click();
  await page.getByRole("link", {name: /下载 GRASP ZIP/}).waitFor({timeout: 10000});
  check(submittedExport?.experimental_ir?.[0]?.constraint === expectedExperimentalIr, "UI submits the mapped experimental IR");
  check(submittedExport?.experimental_ir?.[0]?.evidence_id === testEffectId
    && submittedExport.experimental_ir[0].details.probe_id === testProbeId,
  "UI retains the effect analysis and probe provenance");

  showFactorialResult = true;
  await page.reload({waitUntil: "domcontentloaded"});
  await page.waitForSelector("#project-list article button", {timeout: 15000});
  await page.locator("#project-list article button").filter({hasText: project.project_id}).first().click();
  await page.waitForFunction(() => document.querySelector("#project-content")?.textContent.includes("双侧交互：百分点差中之差"), {timeout: 15000});
  check((await page.locator("#project-content").innerText()).includes("10.00 pp")
    && (await page.locator("#project-content").innerText()).includes("不足 3 个独立批次")
    && await page.locator('table[aria-label="逐探针独立批次效应热图"]').count() === 0,
  "factorial interaction display preserves raw four-cell comparison without site attribution or antigen-effect heatmap");

  const effectMode = page.locator('#project-content select[name="comparisonMode"]');
  const modeValues = await effectMode.locator("option").evaluateAll(options => options.map(option => option.value));
  check(modeValues.includes("ANTIBODY_EFFECT"), "standalone antibody-effect mode is offered in the workbench");
  await effectMode.selectOption("ANTIBODY_EFFECT");
  const antibodyConditionValues = await page.locator('#project-content select[aria-label="效应条件"]').first()
    .locator("option").evaluateAll(options => options.map(option => option.value));
  const antibodyProbeNames = await page.locator('#project-content select[aria-label="突变探针"]').first()
    .locator("option").allTextContents();
  check(antibodyConditionValues.join(",") === ",REFERENCE,REFERENCE_NC"
    && antibodyProbeNames[0] === "选择固定 WT 抗原"
    && await page.getByRole("button", {name: "保存并计算抗体变体效应"}).count() === 1,
  "antibody-effect mode restricts conditions to WT activation/NC and probes to the fixed WT construct");

  showFactorialResult = false;
  showAntibodyResult = true;
  await page.reload({waitUntil: "domcontentloaded"});
  await page.waitForSelector("#project-list article button", {timeout: 15000});
  await page.locator("#project-list article button").filter({hasText: project.project_id}).first().click();
  await page.waitForFunction(() => document.querySelector("#project-content")?.textContent.includes("亲本与变体 Fab：各自 NC 校正的 WT 阳性率窗口"), {timeout: 15000});
  const antibodyEffectTable = page.locator('table[aria-label="亲本与变体抗体各自 NC 校正的 WT 窗口"]');
  const antibodySummaryTable = page.locator('table[aria-label="跨批次抗体变体效应汇总"]');
  check(await antibodyEffectTable.count() === 1 && await antibodyEffectTable.locator("tr").count() === 2
    && (await antibodyEffectTable.innerText()).includes("70.00 pp")
    && (await antibodyEffectTable.innerText()).includes("-20.00 pp")
    && await antibodySummaryTable.count() === 1
    && (await antibodySummaryTable.innerText()).includes("不足 3 个有效批次")
    && await page.locator('table[aria-label="逐探针独立批次效应热图"]').count() === 0,
  "antibody-effect results show own-NC windows and valid-batch summary without antigen mapping");

  await page.evaluate(() => {
    window.__paraplumeProgressSeen = [];
    const root = document.querySelector("#project-content");
    const observer = new MutationObserver(() => {
      for (const button of root.querySelectorAll("button")) {
        window.__paraplumeProgressSeen.push(button.textContent);
      }
    });
    observer.observe(root, {subtree: true, childList: true, characterData: true});
  });
  await page.getByRole("button", {name: "重新运行完整 Paraplume"}).click();
  await page.waitForFunction(() => window.__paraplumeProgressSeen
    ?.some(text => text.includes("正在计算 IgBert 残基嵌入")), {timeout: 10000});
  await page.waitForFunction(() => document.querySelector("#project-content")?.textContent.includes("GRASP IR 导出"), {timeout: 15000});
  const paraplumeProgressSeen = await page.evaluate(() => window.__paraplumeProgressSeen || []);
  check(submittedParaplume?.model_size === "large"
    && paraplumeProgressSeen.some(text => text.includes("14% · 正在加载完整 Paraplume 模型"))
    && paraplumeProgressSeen.some(text => text.includes("86% · ProtT5 残基嵌入完成")),
  "Paraplume SSE displays embedding-stage progress and submits the full large model");

  const unresolvedPreviewId = "browser-fixture-unresolved-preview";
  const unresolvedImportId = "browser-fixture-unresolved-import";
  const firstFcsArtifactId = projectDetails.artifacts.find(item => item.kind === "raw_fcs")?.id;
  const mismatchPreviewId = "browser-fixture-mismatch-preview";
  const mismatchImportId = "browser-fixture-mismatch-import";
  const mismatchPage = await browser.newPage({viewport: {width: 1280, height: 900}});
  await mismatchPage.route(`**/api/v4/experiment-projects/${project.project_id}`, async route => {
    const response = await route.fetch();
    const body = await response.json();
    body.data.artifacts = [...body.data.artifacts, {
      id: mismatchPreviewId, kind: "flowjo_import_preview",
      relative_path: `analyses/flowjo_imports/${mismatchImportId}/preview.json`,
      sha256: "mismatch-fixture", size_bytes: 1, media_type: "application/json", created_at: "9998-01-01T00:00:00Z",
    }];
    await route.fulfill({response, json: body});
  });
  await mismatchPage.route(`**/api/v4/experiment-projects/${project.project_id}/files/${mismatchPreviewId}`, route => route.fulfill({json: {
    schema_version: "1.0", status: "COMPLETE", flowjo_input_artifact_ids: [firstFcsArtifactId],
    records: {
      status: "COMPLETE", unresolved_nodes: [],
      gate_count_verification: {status: "MISMATCH", compared_gate_count: 3,
        matched_gate_count: 0, mismatched_gate_count: 3, unmatched_gate_count: 0},
      gate_method_consistency: {status: "CONSISTENT", samples_checked: 1, method_variants: 1,
        gate_coordinates_compared: false},
      gate_counts: [{sample_id: "H1.fcs", gate_name: "Positive", gate_path: ["root"],
        parent_gate: "root", count: 75, absolute_fraction: 0.75, parent_fraction: 0.75}],
      samples: [{sample_id: "H1.fcs", event_count: 100, gates: [],
        qc: {status: "REVIEW_REQUIRED", findings: [], events_examined: 100, time_windows: []}}],
    },
  }}));
  await mismatchPage.goto(`${server}/experiments`, {waitUntil: "domcontentloaded"});
  await mismatchPage.locator("#project-list article button").filter({hasText: project.project_id}).first().click();
  await mismatchPage.waitForFunction(() => document.querySelector("#project-content")?.textContent
    .includes("与 FlowJo 保存的百分比不同属于方法差异提示"), {timeout: 15000});
  check((await mismatchPage.locator("#project-content").innerText())
    .includes("不会阻止计数进入 QC"),
  "FlowJo cache count differences stay visible but do not block consistent-batch counts");
  check(await mismatchPage.locator("#project-content").getByRole("button", {
    name: "保存修订并对已匹配 FCS 执行",
  }).count() === 1,
  "manual gate revision remains available after an informational count difference");
  await mismatchPage.close();

  await page.route(`**/api/v4/experiment-projects/${project.project_id}`, route => route.fulfill({json: {data: {
    ...projectDetails,
    artifacts: [...projectDetails.artifacts, {
      id: unresolvedPreviewId, kind: "flowjo_import_preview",
      relative_path: `analyses/flowjo_imports/${unresolvedImportId}/preview.json`,
      sha256: "fixture", size_bytes: 1, media_type: "application/json", created_at: "9999-01-01T00:00:00Z",
    }],
  }}}));
  await page.route(`**/api/v4/experiment-projects/${project.project_id}/files/${unresolvedPreviewId}`, route => route.fulfill({json: {
    schema_version: "1.0", status: "WAITING_FOR_USER",
    flowjo_input_artifact_ids: [firstFcsArtifactId],
    records: {status: "WAITING_FOR_USER", gate_counts: [], samples: [{
      sample_id: "H1.fcs", event_count: 13863, gates: [],
    }], unresolved_nodes: [{
      sample_id: "H1.fcs", name: "Lymphocytes", parent_path: ["root"],
      source_gate_type: "EllipsoidGate", channels: ["FSC-A", "SSC-A"],
    }]},
  }}));
  await page.reload({waitUntil: "domcontentloaded"});
  await page.waitForSelector("#project-list article button", {timeout: 15000});
  await page.locator("#project-list article button").filter({hasText: project.project_id}).first().click();
  await page.waitForFunction(() => document.querySelector("#project-content")?.textContent.includes("EllipsoidGate"), {timeout: 15000});
  const rebuildGateForm = page.locator("#project-content form").filter({
    has: page.getByRole("button", {name: "保存修订并对已匹配 FCS 执行"}),
  });
  check((await page.locator("#project-content").innerText()).includes("计数不会用于 QC 或效应分析"),
    "unresolved FlowJo counts are visibly blocked from interpretation");
  check(await page.getByRole("button", {name: "保存 FlowJo 样本 QC"}).count() === 0,
    "unresolved FlowJo preview does not offer an unusable QC review form");
  check(await rebuildGateForm.locator('input[name="flowjoAcknowledged"][required]').count() === 1,
    "manual WSP rebuild requires an explicit completeness acknowledgement");
  check(await rebuildGateForm.locator('input[name="flowjoReviewer"][required]').count() === 1
    && await rebuildGateForm.locator('textarea[name="flowjoReason"][required]').count() === 1,
  "manual WSP rebuild records reviewer and reason");
  check(await rebuildGateForm.locator('textarea[name="gateTree"][required]').count() === 1,
    "manual WSP rebuild requires a complete internal gate tree");

  const mappingPage = await browser.newPage({viewport: {width: 1440, height: 1000}});
  let temporaryMatchPayload = null;
  let resolveTemporaryMatch;
  const temporaryMatchRequest = new Promise(resolve => { resolveTemporaryMatch = resolve; });
  mappingPage.on("dialog", dialog => dialog.accept());
  await mappingPage.route(`**/api/v4/experiment-projects/${project.project_id}/samples/import`, route =>
    route.fulfill({json: {data: {upload_session_id: "browser-temporary-label-session", files: [{
      upload_id: "browser-temporary-label-upload", original_filename: "E1.fcs", event_count: 100,
      sha256: "d".repeat(64),
    }]}}}));
  await mappingPage.route(`**/api/v4/experiment-projects/${project.project_id}/samples/confirm`, async route => {
    temporaryMatchPayload = route.request().postDataJSON();
    resolveTemporaryMatch();
    await route.fulfill({json: {data: {upload_session_id: temporaryMatchPayload.upload_session_id,
      sample_ids: ["browser-temporary-label-sample"], unmatched_upload_ids: [],
      project_revision: (project.revision || 13) + 1, status: "CONFIRMED"}}});
  });
  await mappingPage.goto(`${server}/experiments`, {waitUntil: "domcontentloaded"});
  await mappingPage.locator("#project-list article button").filter({hasText: project.project_id}).first().click();
  await mappingPage.waitForFunction(() => document.querySelector("#project-content")?.textContent
    .includes("FCS 导入与样本匹配"), {timeout: 15000});
  await mappingPage.locator('#project-content input[name="files"]').setInputFiles({
    name: "E1.fcs", mimeType: "application/octet-stream", buffer: Buffer.from("browser fixture"),
  });
  await mappingPage.getByRole("button", {name: "解析并暂存"}).click();
  const temporarySection = mappingPage.locator("#project-content section").filter({hasText: "确认文件与样本对应关系"}).last();
  await temporarySection.waitFor({state: "visible", timeout: 10000});
  await temporarySection.locator("select").first().selectOption("ACTIVATION");
  await temporarySection.getByLabel("使用临时身份标签").check();
  await temporarySection.getByLabel("临时抗原标签").fill("CTLA4-E1");
  await temporarySection.getByLabel("临时抗原类别").selectOption("VARIANT");
  await temporarySection.getByLabel("临时面板标签").fill("L32");
  await temporarySection.getByRole("button", {name: "保存已选择匹配"}).click();
  await temporaryMatchRequest;
  check(temporaryMatchPayload?.matches?.length === 1
    && temporaryMatchPayload.matches[0].antigen_construct_id === null
    && temporaryMatchPayload.matches[0].antibody_construct_id === null
    && temporaryMatchPayload.matches[0].temporary_identity?.status === "TEMPORARY_LABEL"
    && temporaryMatchPayload.matches[0].temporary_identity?.antigen_label === "CTLA4-E1"
    && temporaryMatchPayload.matches[0].temporary_identity?.panel_label === "L32",
  "sample matching preserves explicit temporary labels without inventing construct sequences");
  await temporarySection.waitFor({state: "detached", timeout: 10000});
  await mappingPage.waitForFunction(() => {
    const root = document.querySelector("#project-content");
    const fileInput = root?.querySelector('input[name="files"]');
    return root?.textContent.includes("FCS 导入与样本匹配")
      && fileInput?.getClientRects().length > 0;
  }, {timeout: 30000});
  await mappingPage.close();

  const screenshot = path.join(os.tmpdir(), `ecd-experiments-${project.project_id}.png`);
  await page.screenshot({path: screenshot, fullPage: true});
  check(fs.statSync(screenshot).size > 1000, "experiment workspace screenshot was captured");
  check(pageErrors.length === 0, `browser has no uncaught JavaScript errors (${pageErrors.join("; ")})`);
} finally {
  await browser.close();
}

if (failures.length) process.exitCode = 1;
else console.log("EXPERIMENT WORKSPACE BROWSER CHECKS PASSED");
