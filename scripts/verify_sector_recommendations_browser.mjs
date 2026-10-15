import assert from "node:assert/strict";
import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";

// Usage: node script.mjs <replay directory> <new evidence directory> [playwright module path] [server URL]
const [replayDir, evidenceDir, modulePath = "playwright", baseUrl = "http://127.0.0.1:8765", liveJobId] = process.argv.slice(2);
assert(replayDir && evidenceDir, "Replay and new evidence directories are required");
assert.notEqual(path.resolve(replayDir), path.resolve(evidenceDir), "Preserve replay artifacts");
const playwright = await import(path.isAbsolute(modulePath) ? pathToFileURL(modulePath).href : modulePath);
const { chromium } = playwright.default || playwright;
const result = JSON.parse(await fs.readFile(path.join(replayDir, "analysis.json"), "utf8"));
const recommendations = result.joint_v4_ml.selected_panel.sector_recommendations;
assert.equal(recommendations.length, result.regions.length);
await fs.mkdir(evidenceDir, { recursive: false });
const browser = await chromium.launch({ channel: "chrome", headless: true });
const page = await browser.newPage({ viewport: { width: 1500, height: 1100 } });
const errors = [];
page.on("pageerror", (error) => errors.push(error.message));
const verificationJob = liveJobId || "sector-review-verification";
const structureName = path.basename((result.structure_url || "structure.cif").split("?")[0]);
result.structure_url = `/api/jobs/${verificationJob}/files/${structureName}`;
if (liveJobId) {
  const response = await page.request.get(`${baseUrl}/api/jobs/${liveJobId}`);
  assert.equal(response.status(), 200);
  const job = await response.json();
  assert.equal(job.status, "COMPLETE");
  assert.deepEqual(job.result.joint_v4_ml.selected_panel, result.joint_v4_ml.selected_panel);
  const jobs = await (await page.request.get(`${baseUrl}/api/jobs`)).json();
  assert(jobs.some((item) => item.id === liveJobId), "Live task absent from task list");
  for (const name of ["ecd_structure.cif", "ecd_epitope_design.zip"]) {
    const file = await page.request.get(`${baseUrl}/api/jobs/${liveJobId}/files/${name}`);
    assert.equal(file.status(), 200);
    assert.deepEqual(await file.body(), await fs.readFile(path.join(replayDir, name)));
  }
} else {
await page.route(`**/api/jobs/${verificationJob}`, (route) => route.fulfill({
  json: { id: verificationJob, status: "COMPLETE", message: "Read-only replay verification", progress: 1, result },
}));
await page.route(`**/api/jobs/${verificationJob}/files/*`, async (route) => {
  const name = path.basename(new URL(route.request().url()).pathname);
  const body = await fs.readFile(path.join(replayDir, name));
  await route.fulfill({ body, contentType: name.endsWith(".json") ? "application/json" : "text/plain" });
});
}
try {
  await page.goto(`${baseUrl}/?job=${verificationJob}`, { waitUntil: "domcontentloaded" });
  await page.waitForSelector("#design-list .design-series", { state: "attached", timeout: 90000 });
  await page.locator('#steps [data-step="5"]').click();
  const text = await page.locator("#design-list").innerText();
  for (const item of recommendations) {
    assert(text.includes(item.sector_id), `Missing sector ${item.sector_id}`);
    for (const gate of item.failed_gates || []) assert(text.includes(gate), `Missing failed gate ${gate}`);
    for (const detail of item.gate_details || []) {
      if (detail.reason) assert(text.includes(detail.reason), `Missing reason ${detail.reason}`);
      if (detail.observed != null) assert(text.includes("模型值"), "Prediction mislabeled");
      if (detail.evidence_missing) assert(text.includes("证据缺失"), "Missing evidence not disclosed");
    }
  }
  assert(text.includes("geodesic_0"));
  if (recommendations.some((item) => item.regional_geometry_support?.scenarios)) {
    assert(text.includes("本区代理覆盖"));
    assert(text.includes("未覆盖残基"));
    for (const radius of [9, 12, 15]) assert(text.includes(`geodesic_${radius}`));
    assert(text.includes("不保证100%真实扰动"));
  }
  assert(!text.includes("正式主面板已选"));
  const region6 = recommendations.find((item) => item.sector_id === "VS006");
  if (region6?.design_id && region6.status === "PASSED") {
    const series = page.locator("#design-list .design-series").filter({ has: page.locator("h3", { hasText: /^补丁 6 ·/ }) });
    assert.equal(await series.count(), 1, "Actual region 6 has no construct section");
    assert((await series.innerText()).includes(region6.design_id), "Region 6 construct placed under generation window");
    await series.locator(`[data-review-design="${region6.design_id}"]`).first().click();
    assert.equal(await page.locator("#design-select").inputValue(), region6.design_id);
  }
  const review = recommendations.find((item) => item.design_id && item.status === "NEAR_THRESHOLD");
  if (review) {
    await page.locator(`[data-review-design="${review.design_id}"]`).first().click();
    assert.equal(await page.locator("#design-select").inputValue(), review.design_id);
  }
  await page.screenshot({ path: path.join(evidenceDir, "sector_recommendations.png"), fullPage: true });
  await page.locator('#steps [data-step="7"]').click();
  assert((await page.locator('[data-step-content="7"]').innerText()).includes("逐区建议"));
  assert.deepEqual(errors, []);
  await fs.writeFile(path.join(evidenceDir, "browser_verification.json"), JSON.stringify({
    status: "PASS", live_job_id: liveJobId || null, mocked_routes: !liveJobId,
    url: `${baseUrl}/?job=${verificationJob}`, replayDir: path.resolve(replayDir), sector_count: recommendations.length,
    selected_primary_count: result.joint_v4_ml.selected_panel.selected_candidate_ids.length,
    review_count: recommendations.filter((item) => item.status === "NEAR_THRESHOLD").length,
    exercised_review_design: review?.design_id || null, page_errors: errors,
  }, null, 2));
} finally {
  await browser.close();
}
