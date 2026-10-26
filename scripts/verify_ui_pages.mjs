import assert from "node:assert/strict";
import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";
import { execFileSync } from "node:child_process";
import { createHash } from "node:crypto";

const [evidenceDir, modulePath = "playwright", baseUrl = "http://127.0.0.1:8765"] = process.argv.slice(2);
assert(evidenceDir);
await fs.mkdir(evidenceDir, { recursive: false });
const { chromium } = await import(path.isAbsolute(modulePath) ? pathToFileURL(modulePath).href : modulePath);
const browser = await chromium.launch({ channel: "chrome", headless: true });
const checks = [], errors = [], writes = [], assets = [];
try {
  for (const file of ["web/surface_partitions.html", "web/ctla4_six_mutants.html"]) {
    const original = execFileSync("git", ["show", `HEAD:${file}`], { encoding: "utf8", maxBuffer: 32 * 1024 * 1024 });
    const current = await fs.readFile(file, "utf8");
    const scripts = (html) => [...html.matchAll(/<script\b[^>]*>([\s\S]*?)<\/script>/gi)].map((match) => match[1].replaceAll("\r\n", "\n")).join("\n");
    const hash = (text) => createHash("sha256").update(text).digest("hex");
    assert.equal(hash(scripts(current)), hash(scripts(original)), `${file}: embedded science or viewer scripts changed`);
    assets.push({ file, scientific_scripts_sha256: hash(scripts(current)), unchanged_from_head: true });
  }
  for (const viewport of [{ width: 1366, height: 768 }, { width: 1440, height: 900 }, { width: 390, height: 844 }]) {
    for (const [name, url] of [
      ["experiments", "/experiments"],
      ["pd1", "/pd1-antibodies"],
      ["surface", "/assets/surface_partitions.html"],
      ["ctla4", "/assets/ctla4_six_mutants.html"],
    ]) {
      const page = await browser.newPage({ viewport });
      page.on("pageerror", (error) => errors.push({ name, width: viewport.width, message: error.message }));
      await page.route("**/api/**", (route) => {
        if (["GET", "HEAD"].includes(route.request().method())) return route.continue();
        writes.push(route.request().url());
        return route.abort();
      });
      const response = await page.goto(`${baseUrl}${url}`, { waitUntil: "domcontentloaded" });
      assert.equal(response.status(), 200);
      if (name === "experiments") {
        await page.locator('#project-list article[data-project-id="e0f2c26650f04681"] button').click();
        await page.waitForSelector(".operation-details");
        await page.waitForSelector("#structure-effects canvas");
        assert.equal(await page.locator("#project-list [aria-current='true']").count(), 1);
        const disclosure = page.locator(".operation-details > summary").first();
        await disclosure.focus();
        await page.keyboard.press("Enter");
        assert(await disclosure.evaluate((element) => element.parentElement.open));
        await page.keyboard.press("Enter");
        assert.equal(await page.locator(".operation-details .flowjo-percentage-qc-form").count(), 0);
        const dimensions = await page.evaluate(() => [...document.querySelectorAll("#project-content table")].map((element) => ({ width: element.clientWidth, scroll: element.scrollWidth, overflow: getComputedStyle(element).overflowX })));
        assert(dimensions.filter((item) => item.scroll > item.width + 1).every((item) => item.overflow === "auto"), "Wide tables must remain scrollable");
        assert.equal(await page.locator("#project-content").evaluate((element) => getComputedStyle(element).overflowX), "visible");
        await page.locator(".structure-effect-view").scrollIntoViewIfNeeded();
      } else if (name === "pd1") {
        await page.waitForFunction(() => document.querySelector("#status").textContent.includes("结构已对齐"));
        await page.locator("#focus").selectOption("pembro");
        await page.waitForFunction(() => document.querySelector("#status").textContent.includes("结构已对齐"));
      } else if (name === "ctla4") {
        await page.waitForFunction(() => document.querySelector("#status").textContent.includes("已加载"));
      }
      await page.waitForSelector("canvas");
      await page.waitForTimeout(500);
      assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), `${name}: horizontal overflow`);
      const canvas = page.locator(name === "experiments" ? "#structure-effects canvas" : "canvas");
      assert.equal(await canvas.count(), 1, `${name}: duplicate structure canvas`);
      await page.screenshot({ path: path.join(evidenceDir, `${name}_${viewport.width}.png`) });
      checks.push({ name, viewport, document_width: await page.evaluate(() => document.documentElement.scrollWidth), canvas: await canvas.boundingBox() });
      console.log(`PASS ${name} ${viewport.width}`);
      await page.close();
    }
  }
  const offlineBytes = Buffer.from(await (await fetch(`${baseUrl}/api/jobs/be18bf6290fb/export/offline?lang=zh`)).arrayBuffer());
  const offlinePath = path.resolve(evidenceDir, "offline_report.html");
  await fs.writeFile(offlinePath, offlineBytes);
  const offline = await browser.newPage({ viewport: { width: 1366, height: 768 } });
  const external = [];
  offline.on("pageerror", (error) => errors.push({ name: "offline", message: error.message }));
  await offline.route(/^https?:/, (route) => { external.push(route.request().url()); return route.abort(); });
  await offline.goto(pathToFileURL(offlinePath).href, { waitUntil: "domcontentloaded" });
  await offline.waitForSelector("canvas");
  await offline.locator('#steps [data-step="5"]').click();
  assert((await offline.locator("#panel").innerText()).includes("待审查项不计入合格主面板"));
  assert.equal(await offline.locator("#steps [data-step]").count(), 7);
  assert.deepEqual(external, []);
  assert(await offline.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1));
  await offline.screenshot({ path: path.join(evidenceDir, "offline_recommendations_1366.png") });
  assets.push({ file: "offline_report.html", sha256: createHash("sha256").update(offlineBytes).digest("hex"), bytes: offlineBytes.length, external_requests: external });
  await offline.close();
  assert.deepEqual(errors, []);
  assert.deepEqual(writes, []);
  await fs.writeFile(path.join(evidenceDir, "verification.json"), JSON.stringify({ status: "PASS", checks, assets, page_errors: errors, writes }, null, 2));
} finally { await browser.close(); }
