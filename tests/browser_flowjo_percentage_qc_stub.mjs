/* Synthetic UI submission smoke; every QC write is intercepted and never reaches the API. */
import assert from "node:assert/strict";
import fs from "node:fs";
import {pathToFileURL} from "node:url";

const server = process.env.ECD_WORKBENCH_URL || "http://127.0.0.1:8766";
const projectId = process.env.ECD_PROJECT_ID;
if (!projectId) throw new Error("ECD_PROJECT_ID is required");
const playwrightPath = process.env.PLAYWRIGHT_CORE ||
  "C:/Users/陈淏文/.workbuddy/binaries/node/workspace/node_modules/playwright-core/index.js";
const chromePath = process.env.CHROME_PATH ||
  "C:/Program Files/Google/Chrome/Application/chrome.exe";
if (!fs.existsSync(playwrightPath)) throw new Error(`Playwright module not found: ${playwrightPath}`);
if (!fs.existsSync(chromePath)) throw new Error(`Chrome executable not found: ${chromePath}`);

const apiBase = `${server}/api/v4/experiment-projects/${encodeURIComponent(projectId)}`;
const [projectResponse, importsResponse, reviewsResponse] = await Promise.all([
  fetch(apiBase), fetch(`${apiBase}/flowjo-percentage-imports`), fetch(`${apiBase}/qc-reviews`),
]);
assert.equal(projectResponse.status, 200, "project should be readable");
assert.equal(importsResponse.status, 200, "percentage imports should be readable");
assert.equal(reviewsResponse.status, 200, "QC review ledger should be readable");
const project = (await projectResponse.json()).data;
const imports = (await importsResponse.json()).data;
const initialReviews = (await reviewsResponse.json()).data;
assert.ok(imports.length > 0, "the project should contain an imported result table");
const imported = imports[0];
const observation = imported.result.observations[0];
assert.ok(observation?.sample_id, "a result-table observation should map to a sample");
const staleEffect = {
  analysis_id: "synthetic-old-effect",
  result: {
    effects: [{
      panel_label: "synthetic-panel", probe_label: "synthetic-probe", replicate_type: "TECHNICAL",
      batch_effect: {reference_fraction: 0.5, mutant_fraction: 0.4, mutant_nc_fraction: 0.05,
        difference_pp: -10, decline_pct: 20},
      replicate_count_by_group: {mutant: 1, reference: 1, nc: 1},
      replicate_observations: [{sample_id: observation.sample_id, well_id: observation.well_id,
        fraction: observation.fraction, percentage: observation.percentage, comparison_group: "MUTANT",
        qc_status: "PASS", accepted_for_summary: true, qc_review_id: "older-review"}],
    }],
    provenance: {compensation: {status: "NOT_APPLIED"}},
  },
};

const playwright = await import(pathToFileURL(playwrightPath).href);
const chromium = playwright.chromium || playwright.default?.chromium;
const browser = await chromium.launch({executablePath: chromePath, headless: true});
const page = await browser.newPage({viewport: {width: 1440, height: 1000}});
const pageErrors = [];
const pageWrites = [];
let capturedReview = null;
const capturedReviews = [];
const syntheticReviews = [];
page.on("pageerror", error => pageErrors.push(error.message));
page.on("request", request => {
  if (!["GET", "HEAD"].includes(request.method())) pageWrites.push(request.method());
});

await page.route(`**/api/v4/experiment-projects/${projectId}/qc-reviews`, async route => {
  if (route.request().method() === "POST") {
    capturedReview = JSON.parse(route.request().postData() || "{}");
    assert.equal(capturedReview.expected_revision, project.revision + capturedReviews.length,
      "every QC decision should use the latest revision from the preceding decision");
    capturedReviews.push(capturedReview);
    const syntheticReview = {
      id: `synthetic-browser-review-${capturedReviews.length}`,
      qc_review_id: `synthetic-browser-review-${capturedReviews.length}`,
      source_type: "FLOWJO_PERCENTAGE_IMPORT", source_id: imported.id,
      sample_id: capturedReview.sample_id, status: capturedReview.status,
      reviewer: capturedReview.reviewer, reason: capturedReview.reason,
      version: 1, project_revision: capturedReview.expected_revision + 1,
    };
    syntheticReviews.push(syntheticReview);
    await route.fulfill({
      status: 201, contentType: "application/json",
      body: JSON.stringify({schema_version: "1.0", data: syntheticReview}),
    });
    return;
  }
  if (!syntheticReviews.length) {
    await route.continue();
    return;
  }
  const response = await route.fetch();
  const payload = await response.json();
  payload.data = [...(payload.data || []), ...syntheticReviews];
  await route.fulfill({response, json: payload});
});
await page.route(`**/api/v4/experiment-projects/${projectId}/flowjo-percentage-imports/${imported.id}/effect-analyses`, route =>
  route.fulfill({json: {data: [staleEffect]}}));

try {
  const response = await page.goto(
    `${server}/experiments?project_id=${encodeURIComponent(projectId)}`,
    {waitUntil: "domcontentloaded"},
  );
  assert.equal(response?.status(), 200, "experiment workspace should load");
  const form = page.locator(".flowjo-percentage-qc-form").first();
  await form.waitFor({state: "visible", timeout: 20000});
  assert.equal(await page.getByText("最新人工 QC 与此只读效应版本不一致；请基于当前审核重新计算。", {exact: true}).count(), 1,
    "the workspace should identify an effect that cites an older QC review");
  const reviewerTools = page.locator(".flowjo-percentage-qc-reviewer-tools").first();
  const firstTableCard = reviewerTools.locator("xpath=..");
  const firstTableForms = firstTableCard.locator(".flowjo-percentage-qc-form");
  await reviewerTools.getByLabel("本表默认审核人").fill("Synthetic batch reviewer");
  await reviewerTools.getByRole("button", {name: "填入本表空白审核人"}).click();
  await page.waitForFunction((expected) => document.querySelector(".flowjo-percentage-qc-reviewer-tools [role='status']")
    ?.textContent.includes(`已填入 ${expected} 个空白审核人`), imported.result.observations.length, {timeout: 20000});
  assert.equal(await firstTableForms.count(), imported.result.observations.length,
    "reviewer autofill should be scoped to the selected source table");
  assert.equal(await firstTableForms.evaluateAll(forms => forms.every(form =>
    form.querySelector('input[name="reviewer"]').value === "Synthetic batch reviewer"
      && form.querySelector('select[name="status"]').value === "")), true,
  "reviewer autofill should change only reviewer fields and leave QC decisions untouched");
  assert.equal(capturedReviews.length, 0, "reviewer autofill must not submit a QC decision");
  assert.deepEqual(pageWrites, [], "filling reviewer names alone must issue no write request");
  await form.locator('select[name="status"]').selectOption("HOLD");
  assert.equal(await form.locator('input[name="reason"]').evaluate(input => input.required), true,
    "HOLD requires a reviewer reason");
  await form.locator("button").click();
  assert.equal(capturedReview, null, "native validation must prevent HOLD without a reason");

  await form.locator('input[name="reason"]').fill("Synthetic UI route stub; no real review decision");
  await page.evaluate(() => { window.__qcReviewFirstForm = document.querySelector(".flowjo-percentage-qc-form"); });
  await form.locator("button").click();
  await page.waitForFunction(() => {
    const first = document.querySelector(".flowjo-percentage-qc-form");
    return first?.querySelector(".flowjo-percentage-qc-save-status")?.textContent.includes("HOLD · v1")
      && first?.querySelector('select[name="status"]')?.value === "HOLD"
      && first?.querySelector('input[name="reviewer"]')?.value === "Synthetic batch reviewer";
  }, {timeout: 20000});
  assert.equal(await page.evaluate(() => window.__qcReviewFirstForm === document.querySelector(".flowjo-percentage-qc-form")),
    true, "saving QC should update the row in place without rebuilding the project view");

  assert.deepEqual(capturedReview, {
    flowjo_percentage_import_id: imported.id,
    sample_id: observation.sample_id,
    status: "HOLD",
    reviewer: "Synthetic batch reviewer",
    reason: "Synthetic UI route stub; no real review decision",
    expected_revision: project.revision,
  }, "the UI should bind the decision to the exact table, sample, and current revision");
  const secondForm = page.locator(".flowjo-percentage-qc-form").nth(1);
  await secondForm.locator('select[name="status"]').selectOption("PASS");
  await secondForm.locator("button").click();
  await page.waitForFunction(() => document.querySelectorAll(".flowjo-percentage-qc-save-status")[1]
    ?.textContent.includes("PASS · v1"), {timeout: 20000});
  assert.equal(capturedReviews[1]?.expected_revision, project.revision + 1,
    "the next row should use the revision returned by the preceding QC write");
  assert.equal(capturedReviews[1]?.sample_id, imported.result.observations[1].sample_id,
    "the second QC write should stay bound to its own sample");
  const observationCount = imported.result.observations.length;
  for (let index = 2; index < observationCount; index += 1) {
    const rowForm = page.locator(".flowjo-percentage-qc-form").nth(index);
    await rowForm.locator('select[name="status"]').selectOption("PASS");
    await rowForm.locator("button").click();
    await page.waitForFunction((rowIndex) => document.querySelectorAll(".flowjo-percentage-qc-save-status")[rowIndex]
      ?.textContent.includes("PASS · v1"), index, {timeout: 20000});
  }
  assert.equal(capturedReviews.length, observationCount, "each result-table well should be reviewable in sequence");
  assert.equal(await page.locator("[role='status']").filter({hasText: `人工 QC：${observationCount}/${observationCount}`}).count(), 1,
    "the result table should show that all observations have a human QC decision");
  assert.equal(await page.locator("button").filter({hasText: "基于当前 QC 重新计算描述性效应"}).isEnabled(), true,
    "recalculation should become available after every result-table observation has a review");
  assert.equal(await page.evaluate(() => window.__qcReviewFirstForm === document.querySelector(".flowjo-percentage-qc-form")),
    true, "later QC saves should preserve the in-place review table");
  assert.deepEqual(pageWrites, Array(observationCount).fill("POST"),
    "QC saves should use only the intercepted writes, without page reloads");
  assert.deepEqual(pageErrors, [], `browser page errors: ${pageErrors.join("; ")}`);

  const persistedReviews = (await (await fetch(`${apiBase}/qc-reviews`)).json()).data;
  assert.equal(persistedReviews.length, initialReviews.length,
    "the route stub must leave the real QC ledger unchanged");
  assert.equal(persistedReviews.some(item => item.reviewer === "Synthetic batch reviewer"), false,
    "the synthetic review must not be persisted");
  console.log("PASS synthetic result-table QC UI: reason validation, exact source/sample/revision binding, route-intercepted write, no database change");
} finally {
  await browser.close();
}
