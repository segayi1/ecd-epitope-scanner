/* Read-only browser smoke for table-only QC, without a FlowJo workspace dependency. */
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

const playwright = await import(pathToFileURL(playwrightPath).href);
const chromium = playwright.chromium || playwright.default?.chromium;
const browser = await chromium.launch({executablePath: chromePath, headless: true});
const page = await browser.newPage({viewport: {width: 1440, height: 1000}});
const pageErrors = [];
const writes = [];
page.on("pageerror", error => pageErrors.push(error.message));
page.on("request", request => {
  if (!["GET", "HEAD"].includes(request.method())) writes.push(`${request.method()} ${request.url()}`);
});
await page.route(`**/api/v4/experiment-projects/${projectId}/flowjo-percentage-imports`, async route => {
  const response = await route.fetch();
  const payload = await response.json();
  for (const item of payload.data || []) {
    item.flowjo_import_id = null;
    if (item.result) item.result.flowjo_import_id = null;
  }
  await route.fulfill({response, json: payload});
});
await page.route(
  `**/api/v4/experiment-projects/${projectId}/flowjo-percentage-imports/*/effect-analyses`,
  async route => route.fulfill({
    status: 200,
    contentType: "application/json",
    body: JSON.stringify({data: [{result: {effects: [{
      panel_label: "L32",
      probe_label: "CTLA4-E1",
      reference_label: "CTLA4-WT",
      replicate_type: "TECHNICAL",
      batch_effect: {
        reference_fraction: 0.70,
        mutant_fraction: 0.35,
        mutant_nc_fraction: 0.05,
        difference_pp: -35,
        decline_pct: 53.85,
      },
      replicate_count_by_group: {mutant: 1, reference: 1, nc: 1},
      replicate_observations: [
        {sample_id: "synthetic-a1", well_id: "A1", comparison_group: "REFERENCE",
          condition_label: "CTLA4-WT", fraction: 0.70, percentage: 70,
          displayed_value: "70.0 %", qc_status: "PASS_WITH_OVERRIDE", accepted_for_summary: true},
        {sample_id: "synthetic-a2", well_id: "A2", comparison_group: "MUTANT",
          condition_label: "CTLA4-E1", fraction: 0.35, percentage: 35,
          displayed_value: "35.0 %", qc_status: "PASS", accepted_for_summary: true},
        {sample_id: "synthetic-a3", well_id: "A3", comparison_group: "MUTANT_NC",
          condition_label: "NC", fraction: 0.05, percentage: 5,
          displayed_value: "5.0 %", qc_status: "HOLD", accepted_for_summary: false},
      ],
    }]}}]})
  }),
);

try {
  const response = await page.goto(
    `${server}/experiments?project_id=${encodeURIComponent(projectId)}`,
    {waitUntil: "domcontentloaded"},
  );
  assert.equal(response?.status(), 200, "experiment workspace should load");
  await page.waitForFunction(
    () => document.querySelector("#project-content")?.innerText.includes("已圈门 FlowJo 百分比结果"),
    {timeout: 20000},
  );

  const rows = await page.locator(".flowjo-percentage-qc-form").count();
  const statuses = await page.locator(
    ".flowjo-percentage-qc-form select[name=status]",
  ).evaluateAll(elements => elements.map(element => element.value));
  const newImportCompensation = await page.locator(
    'form select[name="compensation_status"]',
  ).inputValue();
  const newImportDecimalSeparator = await page.locator(
    'form select[name="decimal_separator"]',
  ).inputValue();
  const text = await page.locator("#project-content").innerText();
  if (!rows) console.error(`Rendered content: ${text.slice(-3000)}\nPage errors: ${pageErrors.join("; ")}`);
  assert.ok(rows > 0, "per-well manual QC forms should render");
  assert.ok(statuses.every(value => value === ""), "the read-only smoke must not prefill QC decisions");
  assert.equal(newImportCompensation, "UNKNOWN", "new table imports must not inherit compensation state");
  assert.equal(newImportDecimalSeparator, ".", "new table imports must state their text decimal separator");
  const decimalSeparatorSelector = page.locator('form select[name="decimal_separator"]');
  await decimalSeparatorSelector.selectOption(",");
  assert.equal(await decimalSeparatorSelector.inputValue(), ",", "the comma-decimal input mode should be selectable");
  await decimalSeparatorSelector.selectOption(".");
  assert.match(text, /21\.0\s*%/, "the imported W212 B5 value should be visible");
  assert.match(text, /事件分母/, "the unknown event denominator should be stated");
  assert.match(text, /未补偿/, "the uncompensated acquisition state should be stated");
  const replicateSummary = page.locator(".flowjo-percentage-replicate-details summary").first();
  await replicateSummary.click();
  const replicateDetails = page.locator(".flowjo-percentage-replicate-details").first();
  assert.equal(await replicateDetails.locator("svg circle").count(), 3,
    "the plot should show every synthetic well observation");
  const replicateText = await replicateDetails.innerText();
  assert.match(replicateText, /35\.0\s*%/, "the original per-well table value should be retained");
  assert.match(replicateText, /PASS_WITH_OVERRIDE/, "override QC should remain visible");
  assert.match(replicateText, /HOLD/, "held observations should remain visible");
  assert.match(replicateText, /技术重复展示孔间差异/, "the plot should distinguish technical wells from biological n");
  assert.equal(writes.length, 0, `read-only browser navigation made writes: ${writes.join(", ")}`);
  assert.deepEqual(pageErrors, [], `browser console/page errors: ${pageErrors.join("; ")}`);
  console.log(`PASS read-only result-table workbench without WSP: ${rows} QC forms, no review decisions, no writes`);
} finally {
  await browser.close();
}
