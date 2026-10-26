import assert from "node:assert/strict";
import crypto from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const modulePath = process.argv[2];
const jobIds = process.argv.slice(3);
if (!modulePath || !jobIds.length) {
  throw new Error("Usage: node audit_offline_task_sites.mjs <playwright-index.mjs> <job-id> [...]");
}
const pw = await import(pathToFileURL(path.resolve(modulePath)).href);
const chromium = pw.chromium || pw.default?.chromium;
const repo = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const root = path.join(repo, "workspace_data");
const offlineSiteRoot = process.env.OFFLINE_SITE_ROOT ? path.resolve(process.env.OFFLINE_SITE_ROOT) : root;
const auditDir = process.env.OFFLINE_AUDIT_DIR ? path.resolve(process.env.OFFLINE_AUDIT_DIR) :
  path.join(root, "p0_p3_upgrade_20261016", "offline_browser_audit");
fs.mkdirSync(auditDir, { recursive: true });
const browser = await chromium.launch({
  executablePath: "C:/Program Files/Google/Chrome/Application/chrome.exe",
  headless: true,
});

try {
  const results = [];
  for (const jobId of jobIds) {
    const context = await browser.newContext({ viewport: { width: 1440, height: 900 }, acceptDownloads: true });
    const page = await context.newPage();
    const errors = [];
    const networkRequests = [];
    page.on("pageerror", error => errors.push(error.message));
    page.on("request", request => {
      if (/^https?:\/\//i.test(request.url())) networkRequests.push(request.url());
    });
    await page.route(/^https?:\/\//i, route => route.abort());
    const site = path.join(offlineSiteRoot, jobId, "index.html");
    assert.ok(fs.existsSync(site), `${jobId}: offline package root index is missing`);
    await page.goto(pathToFileURL(site).href, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => window.APP_OFFLINE?.analysis &&
      document.querySelectorAll("#steps [data-step]").length === 7, { timeout: 30000 });
    await page.waitForFunction(() => {
      const plugin = window.__ecdMolstar?.plugin;
      const structures = plugin?.managers?.structure?.hierarchy?.current?.structures;
      const canvas = document.querySelector("#viewer-molstar canvas");
      return Boolean(structures?.length && typeof MolstarView !== "undefined" && MolstarView.isReady() &&
        canvas && canvas.clientWidth > 0 && canvas.clientHeight > 0);
    }, { timeout: 60000 });

    const checks = await page.evaluate(() => ({
      offlineSnapshotLoaded: Boolean(window.APP_OFFLINE?.analysis),
      offlineNoticeVisible: getComputedStyle(document.querySelector(".offline-notice")).display !== "none",
      newAnalysisHidden: getComputedStyle(document.querySelector("#open-run")).display === "none",
      offlineExportHidden: getComputedStyle(document.querySelector("#btn-export-offline")).display === "none",
      steps: document.querySelectorAll("#steps [data-step]").length,
      title: document.title,
      regionCount: window.APP_OFFLINE.analysis.regions?.length ??
        window.APP_OFFLINE.analysis.partitioning?.sector_count ??
        window.APP_OFFLINE.analysis.surface_display_regions_v41?.regions?.length ??
        window.APP_OFFLINE.analysis.surface_regions?.length ?? null,
      designCount: window.APP_OFFLINE.analysis.designs?.length ?? null,
      structureEmbedded: typeof window.APP_OFFLINE.structure_pdb === "string" &&
        window.APP_OFFLINE.structure_pdb.includes("ATOM"),
      structureViewerReady: Boolean(window.__ecdMolstar?.plugin?.managers?.structure?.hierarchy?.current?.structures?.length &&
        typeof MolstarView !== "undefined" && MolstarView.isReady() &&
        document.querySelector("#viewer-molstar canvas")?.clientWidth > 0),
      externalResources: Array.from(document.querySelectorAll("script[src],link[rel=stylesheet][href],img[src],iframe[src],source[src]"))
        .map(node => node.src || node.href).filter(url => /^https?:\/\//i.test(url)),
    }));
    assert.equal(checks.offlineSnapshotLoaded, true, `${jobId}: snapshot not loaded`);
    assert.equal(checks.offlineNoticeVisible, true, `${jobId}: offline notice hidden`);
    assert.equal(checks.newAnalysisHidden, true, `${jobId}: new-analysis control visible`);
    assert.equal(checks.offlineExportHidden, true, `${jobId}: server-export control visible`);
    assert.equal(checks.steps, 7, `${jobId}: expected seven result steps`);
    assert.ok(checks.regionCount > 0, `${jobId}: offline region data missing`);
    assert.equal(checks.structureEmbedded, true, `${jobId}: structure not embedded`);
    assert.deepEqual(checks.externalResources, [], `${jobId}: absolute network resource referenced`);
    const expectedAnalysis = await page.evaluate(() => JSON.parse(JSON.stringify(window.APP_OFFLINE.analysis)));

    for (let step = 1; step <= 7; step += 1) {
      await page.locator(`#steps [data-step="${step}"]`).click();
      assert.equal(await page.locator(`#steps [data-step="${step}"]`).getAttribute("aria-current"), "step",
        `${jobId}: step ${step} did not activate`);
      assert.equal(await page.locator(`#panel [data-step-content="${step}"]`).isVisible(), true,
        `${jobId}: step ${step} panel is not visible`);
    }
    await page.locator('#steps [data-step="1"]').click();
    await page.locator("[data-residue]").first().click();
    await page.locator("#custom-aa").waitFor({ state: "attached" });
    await page.locator("#custom-aa").selectOption(await page.locator("#custom-aa").evaluate(select => select.options[0].value));
    await page.locator("#add-custom-site").click();
    await page.locator('#steps [data-step="6"]').click();
    const downloads = [];
    page.on("download", download => downloads.push(download));
    await page.locator("#save-custom").click();
    for (let attempt = 0; attempt < 40 && downloads.length < 3; attempt += 1) await page.waitForTimeout(100);
    assert.equal(downloads.length, 3, `${jobId}: expected draft JSON and both FASTA downloads`);
    const savedDrafts = [];
    let draft;
    for (const download of downloads) {
      const originalName = download.suggestedFilename();
      const label = originalName.endsWith(".draft.json") ? "draft.json" : originalName.endsWith(".cds.fasta") ? "cds.fasta" : "fasta";
      const destination = path.join(auditDir, `${jobId}_local_${label}${label === "fasta" ? ".fasta" : ""}`);
      await download.saveAs(destination);
      const content = fs.readFileSync(destination);
      savedDrafts.push({ file: path.basename(destination), sha256: crypto.createHash("sha256").update(content).digest("hex") });
      if (label === "draft.json") draft = JSON.parse(content.toString("utf8"));
    }
    assert.deepEqual(savedDrafts.map(row => row.file).sort(), [
      `${jobId}_local_cds.fasta`, `${jobId}_local_draft.json`, `${jobId}_local_fasta.fasta`,
    ].sort(), `${jobId}: unexpected draft download set`);
    assert.equal(draft?.schema, "ecd-epitope-offline-draft-v1", `${jobId}: draft schema missing`);
    assert.equal(draft.status, "LOCAL_DRAFT_NOT_ANALYZED", `${jobId}: draft was not clearly marked unscored`);
    assert.equal(draft.mutations.length, 1, `${jobId}: local mutation missing`);
    assert.ok(draft.sequences.amino_acid && draft.sequences.coding_dna, `${jobId}: AA/CDS sequence missing`);
    assert.deepEqual(Object.keys(draft).sort(), ["mutations", "note", "protein", "schema", "sequences", "status"],
      `${jobId}: draft unexpectedly contains analysis fields`);
    assert.equal(await page.locator("#custom-feedback").textContent().then(text => text.includes("未运行模型")), true,
      `${jobId}: offline draft feedback does not disclose that no model was run`);
    await page.screenshot({ path: path.join(auditDir, `${jobId}.png`), fullPage: true });
    assert.deepEqual(networkRequests, [], `${jobId}: attempted an HTTP(S) request`);
    assert.deepEqual(errors, [], `${jobId}: browser runtime errors`);
    const reportPath = path.join(offlineSiteRoot, jobId, "offline_report.html");
    assert.ok(fs.existsSync(reportPath), `${jobId}: standalone offline report is missing`);
    const reportPage = await context.newPage();
    const reportErrors = [];
    const reportRequests = [];
    reportPage.on("pageerror", error => reportErrors.push(error.message));
    reportPage.on("request", request => {
      if (/^https?:\/\//i.test(request.url())) reportRequests.push(request.url());
    });
    await reportPage.route(/^https?:\/\//i, route => route.abort());
    await reportPage.goto(pathToFileURL(reportPath).href, { waitUntil: "domcontentloaded" });
    await reportPage.waitForFunction(() => {
      const canvas = document.querySelector("#viewport canvas");
      return Boolean(window.APP_DATA && typeof window.APP_PDB === "string" &&
        canvas && canvas.clientWidth > 0 && canvas.clientHeight > 0);
    }, { timeout: 30000 });
    const reportChecks = await reportPage.evaluate(expected => {
      const stable = value => Array.isArray(value) ? value.map(stable) :
        value && typeof value === "object" ? Object.fromEntries(Object.keys(value).sort().map(key => [key, stable(value[key])])) : value;
      const comparable = value => {
        const result = JSON.parse(JSON.stringify(value));
        delete result.files;
        delete result.structure_url;
        const trace = result.run_diagnostics?.execution_trace;
        if (trace) {
          delete trace.elapsed_wall_seconds;
          delete trace.stage_elapsed_wall_seconds;
        }
        return stable(result);
      };
      const differences = (left, right, path = "$") => {
        if (JSON.stringify(comparable(left)) === JSON.stringify(comparable(right))) return [];
        if (Array.isArray(left) && Array.isArray(right)) {
          return left.length === right.length ? left.flatMap((item, index) => differences(item, right[index], `${path}[${index}]`)) : [`${path}.length`];
        }
        if (left && right && typeof left === "object" && typeof right === "object") {
          return [...new Set([...Object.keys(left), ...Object.keys(right)])].flatMap(key =>
            Object.hasOwn(left, key) && Object.hasOwn(right, key) ? differences(left[key], right[key], `${path}.${key}`) : [`${path}.${key}`]);
        }
        return [path];
      };
      const analysisMatchesSite = JSON.stringify(comparable(window.APP_DATA)) === JSON.stringify(comparable(expected));
      return {
        analysisMatchesSite,
        analysisDifferences: analysisMatchesSite ? [] : differences(window.APP_DATA, expected).slice(0, 20),
        comparisonScope: "display/scientific payload; excludes transport URL, artifact hashes and elapsed timings",
        reportKeys: Object.keys(window.APP_DATA || {}).sort(),
        siteKeys: Object.keys(expected || {}).sort(),
        embeddedStructure: window.APP_PDB.includes("ATOM"),
        viewerCanvasReady: Boolean(document.querySelector("#viewport canvas")?.clientWidth > 0),
        title: document.title,
        externalResources: Array.from(document.querySelectorAll("script[src],link[rel=stylesheet][href],img[src],iframe[src],source[src]"))
          .map(node => node.src || node.href).filter(url => /^https?:\/\//i.test(url)),
      };
    }, expectedAnalysis);
    console.log(`${jobId} offline_report_checks=${JSON.stringify(reportChecks)}`);
    assert.equal(reportChecks.analysisMatchesSite, true, `${jobId}: standalone report analysis differs from the site`);
    assert.equal(reportChecks.embeddedStructure, true, `${jobId}: standalone report structure is not embedded`);
    assert.equal(reportChecks.viewerCanvasReady, true, `${jobId}: standalone report viewer did not render`);
    assert.deepEqual(reportChecks.externalResources, [], `${jobId}: standalone report references a network resource`);
    await reportPage.locator("#btn-reset").click();
    assert.deepEqual(reportRequests, [], `${jobId}: standalone report attempted an HTTP(S) request`);
    assert.deepEqual(reportErrors, [], `${jobId}: standalone report browser runtime errors`);
    await reportPage.screenshot({ path: path.join(auditDir, `${jobId}_report.png`), fullPage: true });
    results.push({ job_id: jobId, status: "PASS_FILE_URL_NO_NETWORK", opened_from: path.relative(repo, site), checks, network_requests: networkRequests,
      browser_errors: errors, offline_report: { status: "PASS_FILE_URL_NO_NETWORK", checks: reportChecks,
        network_requests: reportRequests, browser_errors: reportErrors },
      offline_draft: { status: draft.status, mutations: draft.mutations, downloads: savedDrafts },
      screenshot: `${jobId}.png` });
    await context.close();
  }
  const receipt = { schema: "offline_task_site_browser_audit_v1", results };
  fs.writeFileSync(path.join(auditDir, "acceptance.json"), JSON.stringify(receipt, null, 2));
  console.log(JSON.stringify(receipt, null, 2));
} finally {
  await browser.close();
}
