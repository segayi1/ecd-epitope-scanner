import assert from "node:assert/strict";
import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";

const [evidenceDir, modulePath = "playwright", jobId = "be18bf6290fb", baseUrl = "http://127.0.0.1:8765"] = process.argv.slice(2);
assert(evidenceDir, "A new evidence directory is required");
await fs.mkdir(evidenceDir, { recursive: false });
const { chromium } = await import(path.isAbsolute(modulePath) ? pathToFileURL(modulePath).href : modulePath);
const browser = await chromium.launch({ channel: "chrome", headless: true });
const errors = [];
const checks = [];
try {
  for (const viewport of [{ width: 1366, height: 768 }, { width: 1440, height: 900 }, { width: 1920, height: 1080 }, { width: 390, height: 844 }]) {
    const context = await browser.newContext({ viewport, permissions: ["clipboard-read", "clipboard-write"] });
    const page = await context.newPage();
    page.on("pageerror", (error) => errors.push({ width: viewport.width, message: error.message }));
    const response = await page.request.get(`${baseUrl}/api/jobs/${jobId}?view=display`);
    assert.equal(response.status(), 200);
    const job = await response.json();
    assert.equal(job.status, "COMPLETE");
    await page.goto(`${baseUrl}/?job=${jobId}`, { waitUntil: "domcontentloaded" });
    await page.waitForFunction(() => window.__ecdMolstar?.assemblyRef &&
      window.__ecdMolstar.plugin.managers.structure.hierarchy.current.structures.length > 0 &&
      document.querySelector("#viewport-empty").classList.contains("hidden"), null, { timeout: 90000 });
    await page.waitForTimeout(600);
    await page.evaluate(() => { window.__uiTestPlugin = window.__ecdMolstar; });
    assert.equal(await page.locator("#viewer-molstar canvas").count(), 1);
    assert.equal(await page.locator("#design-select").inputValue(), "wt");
    assert(await page.locator("#sequence-content").isHidden());
    const initialCamera = await page.evaluate(() => Array.from(window.__ecdMolstar.plugin.canvas3d.camera.getSnapshot().target));
    const initial = await page.locator("#viewport").boundingBox();
    const stage = await page.locator(".page").boundingBox();
    if (viewport.width > 1000) {
      assert(initial.width / stage.width >= 0.63, "Structure steps must reserve most of the work area for the canvas");
      assert(initial.height >= viewport.height * 0.53, "Laptop canvas is too short");
    }
    await page.screenshot({ path: path.join(evidenceDir, `structure_${viewport.width}.png`) });
    for (let step = 1; step <= 7; step += 1) {
      await page.locator(`#steps [data-step="${step}"]`).click();
      assert(await page.locator(`[data-step-content="${step}"]`).isVisible());
      assert.equal(await page.locator("#steps [aria-current='step']").count(), 1);
      assert.equal(await page.locator("#design-select").inputValue(), "wt");
      assert(await page.evaluate(() => window.__uiTestPlugin === window.__ecdMolstar));
      assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), `Overflow on step ${step}`);
    }
    const finalCamera = await page.evaluate(() => Array.from(window.__ecdMolstar.plugin.canvas3d.camera.getSnapshot().target));
    assert(initialCamera.every((value, index) => Math.abs(value - finalCamera[index]) < 0.001), "Step navigation changed the camera target");
    await page.locator('#steps [data-step="5"]').click();
    const resultCanvas = await page.locator("#viewport").boundingBox();
    if (viewport.width > 1000) assert(resultCanvas.width < initial.width - 40, "Result panel did not expand");
    const candidate = page.locator('[data-review-design]').first();
    const candidateId = await candidate.getAttribute("data-review-design");
    await candidate.click();
    assert.equal(await page.locator("#design-select").inputValue(), candidateId);
    const mutationCount = await page.locator("#design-sequence .mut").count();
    assert(mutationCount >= 2, "Construct selection did not update the sequence");
    await page.waitForTimeout(200);
    await page.waitForFunction(() => window.__ecdMolstar.plugin.isBusy === false, null, { timeout: 90000 });
    await page.waitForTimeout(600);
    await page.locator("#panel").evaluate((element) => { element.scrollTop = 0; });
    await page.screenshot({ path: path.join(evidenceDir, `recommendations_${viewport.width}.png`) });
    await page.locator('#steps [data-step="3"]').click();
    const region = page.locator('[data-focus-region]').first();
    await region.focus();
    await page.keyboard.press("Enter");
    assert.equal(await region.getAttribute("aria-pressed"), "true");
    assert(await page.locator("#design-surface-callout").isHidden(), "Region selection remained masked by construct focus");
    await page.locator("#toggle-sequence").click();
    assert.equal(await page.locator("#toggle-sequence").getAttribute("aria-expanded"), "true");
    assert.equal(await page.locator("#wt-sequence .tok").count(), job.result.protein.ecd_sequence.length);
    const residue = page.locator("#wt-sequence .tok").nth(10);
    const position = await residue.getAttribute("data-seq-residue");
    await residue.focus();
    await page.keyboard.press("Enter");
    await page.locator('#steps [data-step="1"]').click();
    assert((await page.locator("#residue-detail").innerText()).includes(position));
    await page.locator("#copy-wt-aa").click();
    assert.equal(await page.evaluate(() => navigator.clipboard.readText()), job.result.protein.ecd_sequence);
    await page.locator("#toggle-sequence").click();
    await page.locator("#open-run").click();
    assert(await page.locator("#run-dialog").isVisible());
    await page.locator(".run-advanced > summary").click();
    assert(await page.locator("#surface-strategy").isVisible());
    await page.keyboard.press("Escape");
    assert(await page.locator("#run-dialog").isHidden());
    await page.locator("#open-run").click();
    await page.locator("#open-upload").click();
    assert(await page.locator("#upload-dialog").isVisible());
    assert(await page.locator("#run-dialog").isHidden());
    await page.keyboard.press("Escape");
    await page.locator("#btn-lang").click();
    assert.equal(await page.locator("html").getAttribute("lang"), "en");
    assert.equal(await page.locator("#toggle-sequence").innerText(), "Expand alignment");
    assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1));
    assert(await page.evaluate(() => window.__uiTestPlugin === window.__ecdMolstar));
    let resultSwitch = "covered_at_1366";
    if (viewport.width === 1366) {
      await page.evaluate(() => {
        window.__uiDisposals = 0;
        const plugin = window.__ecdMolstar.plugin;
        const dispose = plugin.dispose.bind(plugin);
        plugin.dispose = () => { window.__uiDisposals += 1; return dispose(); };
      });
      await page.locator("#load-demo").evaluate((button) => button.click());
      await page.waitForFunction(() => window.__ecdMolstar && window.__ecdMolstar !== window.__uiTestPlugin && document.querySelector("#viewport-empty").classList.contains("hidden"), null, { timeout: 90000 });
      assert.equal(await page.evaluate(() => window.__uiDisposals), 1);
      assert.equal(await page.locator("#viewer-molstar canvas").count(), 1);
      assert.equal(await page.locator("#design-select").inputValue(), "wt");
      resultSwitch = "PASS_live_saved_result_to_frozen_demo";
    }
    checks.push({ viewport, canvas: initial, result_canvas: resultCanvas, mutation_count: mutationCount, steps: 7, copy: "PASS", keyboard: "PASS", single_instance: true, result_switch: resultSwitch });
    await context.close();
  }
  assert.deepEqual(errors, []);
  await fs.writeFile(path.join(evidenceDir, "verification.json"), JSON.stringify({ status: "PASS", job_id: jobId, real_api: true, mocked_routes: false, checks, errors }, null, 2));
  console.log(JSON.stringify({ status: "PASS", checks: checks.length, errors }));
} finally { await browser.close(); }
