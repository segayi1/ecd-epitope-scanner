import assert from "node:assert/strict";
import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL, fileURLToPath } from "node:url";
import { createHash } from "node:crypto";
import { execFileSync } from "node:child_process";

// Read-only replay of a frozen result through the current online and offline UI.
// This intercepts job APIs and is evidence of renderer compatibility only; it
// must never be reported as a live-job API acceptance.
const [outputDir, evidenceDir, modulePath = process.env.PLAYWRIGHT_CORE || "playwright-core",
  baseUrl = "http://127.0.0.1:8765"] = process.argv.slice(2);
assert(outputDir && evidenceDir, "Frozen output and a new evidence directory are required");
assert.notEqual(path.resolve(outputDir), path.resolve(evidenceDir), "Keep browser evidence outside frozen job output");
const playwright = await import(path.isAbsolute(modulePath) ? pathToFileURL(modulePath).href : modulePath);
const { chromium } = playwright.default || playwright;
const read = (file) => fs.readFile(path.join(outputDir, file));
const analysisBytes = await read("analysis.json");
const analysis = JSON.parse(analysisBytes.toString("utf8"));
const recommendations = analysis.joint_v4_ml?.selected_panel?.sector_recommendations || [];
assert(recommendations.length > 0, "Frozen result has no regional recommendations");
assert(recommendations.every((item) => item.recommendation_policy === "regional_mechanism_pareto_v1"),
  "Frozen result does not use the regional mechanism policy in every region");
const offlinePath = path.join(outputDir, "offline_report.html");
await fs.access(offlinePath);
await fs.mkdir(evidenceDir, { recursive: false });
// Match the real ?view=display API rather than injecting the full scientific audit into the browser.
const packageRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const displayPath = path.join(evidenceDir, "display_analysis.json");
execFileSync(path.join(packageRoot, ".venv", "Scripts", "python.exe"), ["-c",
  "import json,sys; from pathlib import Path; sys.path.insert(0,str(Path(sys.argv[1]).parent)); from ecd_epitope_scanner.joint_design_v4.reporting import project_display_analysis; a=json.loads(Path(sys.argv[2]).read_text(encoding='utf8')); Path(sys.argv[3]).write_text(json.dumps(project_display_analysis(a),ensure_ascii=False),encoding='utf8')",
  packageRoot, path.resolve(outputDir, "analysis.json"), path.resolve(displayPath)]);
const displayAnalysis = JSON.parse(await fs.readFile(displayPath, "utf8"));
const offlineBytes = await fs.readFile(offlinePath);
const browser = await chromium.launch({ channel: "chrome", headless: true });
const pageErrors = [];
const checks = [];
const frozenJobId = `frozen-regional-${process.pid}`;
displayAnalysis.structure_url = `/api/jobs/${frozenJobId}/files/ecd_structure.cif`;
const frozenJob = { id: frozenJobId, status: "COMPLETE", progress: 1,
  message: "Read-only frozen result replay", result: displayAnalysis };

async function inspect(page, label, viewportWidth) {
  await page.waitForSelector(label.startsWith("offline") ? '#steps [data-step="5"]' : "#design-list .design-series", { state: "attached", timeout: 90000 });
  await page.locator('#steps [data-step="5"]').click();
  const regions = page.locator(".recommendation-region");
  assert.equal(await regions.count(), recommendations.length, `${label}: regional card count`);
  const panel = page.locator(".design-series");
  const summary = await panel.innerText();
  assert(summary.includes("本区强扰动工程代表"), `${label}: missing engineering representative label`);
  assert(summary.includes("不代表亲和力下降最大"), `${label}: missing affinity claim boundary`);
  assert(summary.includes("最小相对侧链暴露度"), `${label}: missing relative side-chain exposure`);
  assert(summary.includes("12 Å局部覆盖"), `${label}: missing 12 Å coverage summary`);
  assert(summary.includes("表位倾向 P 仅为区域／位点先验注释"), `${label}: missing P annotation boundary`);
  for (const item of recommendations) {
    const card = panel.locator(`[data-recommendation-sector="${item.sector_id}"][data-recommendation-role="primary"]`);
    assert.equal(await card.count(), 1, `${label}: one primary card for ${item.sector_id}`);
    const text = await card.innerText();
    assert(text.includes("本区强扰动工程代表"), `${label}: representative label for ${item.sector_id}`);
    assert(text.includes("直接映射突变面积"), `${label}: direct area for ${item.sector_id}`);
    assert(text.includes("未覆盖面积"), `${label}: uncovered area for ${item.sector_id}`);
    if (item.candidate_id) {
      const details = card.locator("details.recommendation-details");
      await details.locator(":scope > summary").click();
      const detailsText = await details.innerText();
      for (const labelText of ["区内平均绝对化学变化", "变化幅度（非带符号）", "逐位带符号化学变化",
        "氢键供体变化", "氢键受体变化"]) {
        assert(detailsText.includes(labelText), `${label}: ${item.sector_id} missing ${labelText}`);
      }
      for (const radius of [9, 12, 15]) assert(detailsText.includes(`${radius} Å`), `${label}: missing ${radius} Å detail`);
      if (item.impact_evidence?.ranking_version === "mechanism_channels_v2") {
        for (const title of ["各取最强轴效用", "带符号净变化", "净变化为零不代表局部作用抵消", "WT Cα", "不能替代Cβ"]) {
          assert(detailsText.includes(title), `${label}: missing v2 diagnostic ${title}`);
        }
        if (item.resource_audit?.single_threshold_audit) assert(detailsText.includes("单点稳定性门槛边界审计"));
      }
    }
  }
  assert(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 2),
    `${label}: horizontal viewport overflow at ${viewportWidth}px`);
  await page.screenshot({ path: path.join(evidenceDir, `${label}_${viewportWidth}.png`), fullPage: true });
  checks.push({ view: label, viewport_width: viewportWidth, region_count: recommendations.length });
}

try {
  for (const viewport of [{ width: 1500, height: 1100 }, { width: 390, height: 844 }]) {
    const online = await browser.newPage({ viewport });
    online.on("pageerror", (error) => pageErrors.push(`online: ${error.message}`));
    await online.route(`**/api/jobs/${frozenJobId}*`, (route) => route.fulfill({ json: frozenJob }));
    await online.route(`**/api/jobs/${frozenJobId}/files/*`, async (route) => {
      const name = path.basename(new URL(route.request().url()).pathname);
      const bytes = await read(name);
      await route.fulfill({ body: bytes, contentType: name.endsWith(".cif") ? "chemical/x-cif" : "text/plain" });
    });
    await online.goto(`${baseUrl}/?job=${frozenJobId}`, { waitUntil: "domcontentloaded" });
    await inspect(online, `online_frozen`, viewport.width);
    await online.close();

    const offline = await browser.newPage({ viewport });
    offline.on("pageerror", (error) => pageErrors.push(`offline: ${error.message}`));
    await offline.goto(pathToFileURL(offlinePath).href, { waitUntil: "domcontentloaded" });
    await inspect(offline, `offline_frozen`, viewport.width);
    await offline.close();
  }
  assert.deepEqual(pageErrors, [], "No uncaught browser page errors");
  const hash = (bytes) => createHash("sha256").update(bytes).digest("hex");
  await fs.writeFile(path.join(evidenceDir, "browser_verification.json"), JSON.stringify({
    status: "PASS", policy: "regional_mechanism_pareto_v1", verification_kind: "frozen_result_renderer_replay",
    mocked_routes: true, live_job_api_claim: false, base_url: baseUrl,
    online_payload: "project_display_analysis", display_payload_bytes: (await fs.stat(displayPath)).size,
    output_dir: path.resolve(outputDir), analysis_sha256: hash(analysisBytes), offline_report_sha256: hash(offlineBytes),
    region_count: recommendations.length, checks, page_errors: pageErrors,
  }, null, 2), "utf8");
} finally {
  await browser.close();
}
