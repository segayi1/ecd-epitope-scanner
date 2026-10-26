import assert from "node:assert/strict";
import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";

// Uses real job and file endpoints; never intercepts requests.
const [jobId, evidenceDir, modulePath = "playwright", baseUrl = "http://127.0.0.1:8765"] = process.argv.slice(2);
assert(jobId && evidenceDir, "Live job ID and a new evidence directory are required");
const playwright = await import(path.isAbsolute(modulePath) ? pathToFileURL(modulePath).href : modulePath);
const { chromium } = playwright.default || playwright;
const browser = await chromium.launch({ channel: "chrome", headless: true });
await fs.mkdir(evidenceDir, { recursive: false });
const errors = [];
const checks = [];
try {
  for (const viewport of [{ width: 1500, height: 1100 }, { width: 390, height: 844 }]) {
    const page = await browser.newPage({ viewport });
    page.on("pageerror", (error) => errors.push(error.message));
    const response = await page.request.get(`${baseUrl}/api/jobs/${jobId}`);
    assert.equal(response.status(), 200);
    const job = await response.json();
    assert.equal(job.status, "COMPLETE");
    const recommendations = job.result.joint_v4_ml.selected_panel.sector_recommendations;
    for (const candidateId of job.result.joint_v4_ml.selected_panel.selected_candidate_ids) {
      const construct = job.result.designs.find((design) => design.foldx_diagnostics?.candidate_id === candidateId && design.design_role === "STAGE1_PROBE");
      assert(construct, "Formal primary construct is absent");
      for (const key of ["spurs_multi_ddg_kcal_mol", "foldx_ddg_kcal_mol"]) assert(Number.isFinite(construct.foldx_diagnostics[key]), "Formal primary lacks finite combination model evidence");
    }
    assert.equal(recommendations.length, job.result.regions.length);
    await page.goto(`${baseUrl}/?job=${jobId}`, { waitUntil: "domcontentloaded" });
    await page.waitForSelector(".recommendation-region", { state: "attached", timeout: 90000 });
    await page.locator('#steps [data-step="5"]').click();
    assert.equal(await page.locator(".recommendation-region").count(), recommendations.length);
    const summary = await page.locator("#design-list").innerText();
    assert(!/VS\d{3}|geodesic_\d+|V41_[a-f\d]{12,}|COMBINATION_[A-Z_]+/.test(summary), "Internal code in default panel");
    assert(summary.includes("本区估计覆盖") || summary.includes("12 Å局部覆盖") || !recommendations.some((item) => item.candidate_id));
    assert.equal(await page.locator("#design-list .design-series").count(), 1, "Duplicate legacy patch panels");
    for (const item of recommendations) {
      const card = page.locator(`[data-recommendation-sector="${item.sector_id}"][data-recommendation-role="primary"]`);
      assert.equal(await card.count(), 1);
      const heading = await card.locator("h4").first().innerText();
      assert(heading.includes("专属首选") || heading.includes("本区强扰动工程代表"));
      for (const mutation of item.mutation_notations || []) assert((await card.innerText()).includes(mutation));
      if (item.candidate_id) {
        const details = card.locator("details.recommendation-details");
        await details.locator(":scope > summary").click();
        const text = await details.innerText();
        for (const title of ["突变位点", "完整未覆盖列表", "空间互补性", "组合准入检查", "推荐依据", "证据来源"]) assert(text.includes(title));
        for (const gate of item.gate_checks || []) {
          if (typeof gate.observed === "number") assert(text.includes(gate.observed.toFixed(3)), "Missing actual model value");
        }
        assert(!text.includes("尚未解释的审查项"), "Uninterpreted warning in detailed panel");
        assert(!text.includes("\ufffd"), "Damaged explanation in detailed panel");
        if (item.provenance?.structure_evidence?.experimental) {
          assert(text.includes("实验结构：pLDDT不适用"), "Experimental structure mislabeled as missing predicted confidence");
          assert(text.includes("局部精度尚未评估"), "Experimental origin incorrectly promoted to local certainty");
          if (/NMR/i.test(item.provenance.structure_evidence.method || "")) {
            assert(text.includes("结构集合离散度"), "NMR precision described using the wrong experimental method");
            assert(text.includes("NMR不采用晶体学分辨率"), "Missing NMR resolution treated as absent quality evidence");
          }
        }
        if ((item.gate_checks || []).some((gate) => gate.check_type === "rule")) {
          assert(text.includes("不适用：规则检查，无数值门槛"), "Rule check shown as missing numeric evidence");
        }
        assert(!/[a-f\d]{64}|[A-Z]:\\/.test(text), "Technical hashes or paths clutter readable details");
        assert(await details.locator(".recommendation-technical:not([open])").count(), "Technical provenance is not collapsed");
        assert(await details.locator(".recommendation-sequences:not([open])").count(), "Construct sequences are not collapsed");
        assert((await card.innerText()).includes(item.selected_in_panel ? "已纳入正式主面板" : "分区建议 · 未纳入正式主面板"), "Region qualification confused with formal panel selection");
        await details.locator(":scope > summary").click();
      }
      if (item.design_id) {
        await card.locator(`[data-review-design="${item.design_id}"]`).click();
        assert.equal(await page.locator("#design-select").inputValue(), item.design_id);
      }
    }
    const altButton = page.locator('[data-recommendation-role="alternative"] [data-review-design]').first();
    if (await altButton.count()) {
      const designId = await altButton.getAttribute("data-review-design");
      await altButton.click();
      assert.equal(await page.locator("#design-select").inputValue(), designId);
    }
    const firstDetails = page.locator(".recommendation-details > summary").first();
    if (await firstDetails.count()) {
      await firstDetails.focus();
      await page.keyboard.press("Enter");
      assert(await firstDetails.evaluate((element) => element.parentElement.open), "Keyboard expansion failed");
    }
    assert(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 2), "Viewport overflow");
    await page.screenshot({ path: path.join(evidenceDir, `panel_${viewport.width}.png`), fullPage: true });
    checks.push({ width: viewport.width, region_count: recommendations.length, selected_primary_count: job.result.joint_v4_ml.selected_panel.selected_candidate_ids.length });
    await page.close();
  }
  assert.deepEqual(errors, []);
  await fs.writeFile(path.join(evidenceDir, "browser_verification.json"), JSON.stringify({ status: "PASS", job_id: jobId, mocked_routes: false, url: `${baseUrl}/?job=${jobId}`, checks, page_errors: errors }, null, 2));
} finally { await browser.close(); }
