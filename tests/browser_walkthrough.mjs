/* Real-browser walkthrough verification for the integrated online UI.
   Boots nothing itself: expects the workbench server on 127.0.0.1:8765.
   Checks (frozen PD-L1 demo):
     1. the 6 step pills sit on ONE horizontal row
     2. every step switches the right panel
     3. the bottom sequence bar is visible with 4 copy buttons
     4. clicking a sequence token fills the residue detail
     5. full-page screenshots for visual review
   Run: node tests/browser_walkthrough.mjs
   Screenshots land in test_output/ui_check/. */

import fs from "fs";
import path from "path";
import { pathToFileURL } from "url";

const pw = await import(
  pathToFileURL("C:/Users/陈淏文/.workbuddy/binaries/node/workspace/node_modules/playwright-core/index.js").href
);
const chromium = pw.chromium || (pw.default && pw.default.chromium);

const outDir = "E:/user/New project/ecd_epitope_scanner/test_output/ui_check";
fs.mkdirSync(outDir, { recursive: true });

const browser = await chromium.launch({
  executablePath: "C:/Program Files/Google/Chrome/Application/chrome.exe",
  headless: true,
});
const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
const fails = [];
const ok = (cond, label) => {
  console.log((cond ? "  PASS  " : "  FAIL  ") + label);
  if (!cond) fails.push(label);
};

await page.goto("http://127.0.0.1:8765/", { waitUntil: "domcontentloaded" });
await page.waitForSelector("#surface-stats .stat", { timeout: 20000 });
await page.waitForTimeout(2500); // let Mol* settle

/* 1. step pills on one row */
const boxes = await page.$$eval("#steps .step-pill", (els) =>
  els.map((el) => { const r = el.getBoundingClientRect(); return { x: r.x, y: r.y, w: r.width }; })
);
ok(boxes.length === 6, `6 step pills rendered (got ${boxes.length})`);
const rowTolerance = 4;
const sameRow = boxes.every((b) => Math.abs(b.y - boxes[0].y) <= rowTolerance);
ok(sameRow, "step pills share one row (horizontal nav)");
const insideViewport = boxes.every((b) => b.x >= 0 && b.x + b.w <= 1441);
ok(insideViewport, "step pills inside viewport width");

await page.screenshot({ path: path.join(outDir, "01_full.png"), fullPage: false });

/* 2. each step switches the panel */
const expectations = [
  [1, "#residue-list"],
  [2, "#foldx-table"],
  [3, "#region-list"],
  [4, "#iface-partners"],
  [5, "#design-list"],
  [6, "#custom-mutations"],
];
for (const [n, selector] of expectations) {
  await page.click(`#steps [data-step="${n}"]`);
  await page.waitForTimeout(150);
  const section = await page.$(`[data-step-content="${n}"]`);
  const visible = await section.isVisible();
  const target = await page.$(selector);
  const targetVisible = target ? await target.isVisible() : false;
  ok(visible && targetVisible, `step ${n} visible with ${selector}`);
  if (n === 3) await page.screenshot({ path: path.join(outDir, "02_step3_regions.png") });
  if (n === 5) await page.screenshot({ path: path.join(outDir, "03_step5_designs.png") });
}

/* 2.5 persistent seqdock: visible on EVERY step */
const dockPerStep = [];
for (const n of [1, 2, 3, 4, 5, 6]) {
  await page.click(`#steps [data-step="${n}"]`);
  await page.waitForTimeout(120);
  const dock = await page.$(".seqdock");
  const vis = dock ? await dock.isVisible() : false;
  const box = dock ? await dock.boundingBox() : null;
  dockPerStep.push(vis && box && box.height > 60);
}
ok(dockPerStep.every(Boolean), `persistent seqdock visible on all steps (${dockPerStep})`);
await page.click('#steps [data-step="1"]');

/* 3. sequence bar visible with copy buttons */
await page.click('#steps [data-step="1"]');
const seqBar = await page.$(".seqdock");
ok(seqBar && await seqBar.isVisible(), "persistent sequence dock visible");
for (const id of ["copy-wt-aa", "copy-wt-nt", "copy-design-aa", "copy-design-nt"]) {
  const btn = await page.$(`#${id}`);
  ok(btn && await btn.isVisible(), `copy button #${id} visible`);
}
const seqTokens = await page.$$eval("#wt-sequence .tok", (els) => els.length);
ok(seqTokens > 100, `WT sequence rendered (${seqTokens} tokens)`);

/* 4. clicking a sequence token fills residue detail */
await page.click('#wt-sequence [data-seq-residue="100"]');
await page.waitForTimeout(200);
const detail = await page.$eval("#residue-detail", (el) => el.textContent);
ok(detail.length > 50 && !detail.includes("在结构、序列或列表中选择残基"),
  "residue detail filled after sequence click");

/* 5. back to top, final full screenshot */
await page.click('#steps [data-step="1"]');
await page.evaluate(() => window.scrollTo(0, 0));
await page.waitForTimeout(200);
await page.screenshot({ path: path.join(outDir, "04_final.png") });

await browser.close();
console.log(fails.length ? `\nFAILURES: ${fails.length}` : "\nALL CHECKS PASSED");
process.exit(fails.length ? 1 : 0);
