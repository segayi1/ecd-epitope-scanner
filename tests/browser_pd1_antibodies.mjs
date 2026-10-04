import { pathToFileURL } from "url";

const pw = await import(pathToFileURL("C:/Users/陈淏文/.workbuddy/binaries/node/workspace/node_modules/playwright-core/index.js").href);
const chromium = pw.chromium || pw.default?.chromium;
const browser = await chromium.launch({executablePath:"C:/Program Files/Google/Chrome/Application/chrome.exe",headless:true});
const page = await browser.newPage({viewport:{width:1280,height:800}});
const errors = [];
page.on("pageerror", error => errors.push(error.message));
try {
  const response = await page.goto("http://127.0.0.1:8765/pd1-antibodies", {waitUntil:"domcontentloaded"});
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("结构已对齐"), {timeout:30000});
  await page.waitForTimeout(12000);
  const result = await page.evaluate(() => ({
    status:document.querySelector("#status").textContent,
    canvas:document.querySelector("#mol canvas")?.width ?? 0,
    controls:document.querySelectorAll("[data-key]").length
  }));
  console.log(JSON.stringify({http:response.status(),...result,errors}));
  if (response.status() !== 200 || !result.canvas || result.controls !== 4 || errors.length) process.exitCode = 1;
  await page.screenshot({path:"test_output/pd1_antibodies.png"});
  await page.locator('#focus').selectOption('pembro');
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("结构已对齐"), {timeout:30000});
  await page.screenshot({path:"test_output/pd1_pembro_epitope.png"});
  await page.locator('#focus').selectOption('all');
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("结构已对齐"), {timeout:30000});
  await page.locator('[data-key="pembro"]').uncheck();
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("结构已对齐"), {timeout:30000});
  await page.locator('[data-key="pembro"]').check();
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("结构已对齐"), {timeout:30000});
  await page.locator('#focus').selectOption('clone19');
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("结构已对齐"), {timeout:30000});
  await page.screenshot({path:"test_output/pd1_clone19_epitope.png"});
  await page.locator('[data-key="clone19"]').uncheck();
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("结构已对齐"), {timeout:30000});
  await page.screenshot({path:"test_output/pd1_clone19_contact_only.png"});
} finally { await browser.close(); }
