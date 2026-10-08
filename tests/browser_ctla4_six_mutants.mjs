import { pathToFileURL } from "url";

const pw = await import(pathToFileURL("C:/Users/陈淏文/.workbuddy/binaries/node/workspace/node_modules/playwright-core/index.js").href);
const chromium = pw.chromium || pw.default?.chromium;
const browser = await chromium.launch({executablePath:"C:/Program Files/Google/Chrome/Application/chrome.exe",headless:true});
const page = await browser.newPage({viewport:{width:1500,height:950}});
const errors = [];
const structureResponses = [];
page.on("pageerror", error => errors.push(error.message));
page.on("response", response => {
  if (response.url().includes("/structures/")) structureResponses.push({url:response.url(),status:response.status()});
});
try {
  const home = await page.goto("http://127.0.0.1:8765/", {waitUntil:"domcontentloaded"});
  const homeLink = page.locator('a[href="/assets/ctla4_six_mutants.html"]');
  if (home.status() !== 200 || await homeLink.count() !== 1) throw new Error("CTLA-4 structure link missing from the workbench");
  const response = await page.goto("http://127.0.0.1:8765/assets/ctla4_six_mutants.html", {waitUntil:"domcontentloaded"});
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("已加载"), {timeout:30000});
  await page.screenshot({path:"test_output/ctla4_six_mutants.png",fullPage:true});
  const labels = ["A61W/M134R","T67R/M134R","Q117R/E155W","S55F/S101R","D78A/M122W","S55F/I102R"];
  const wtContentCounts = [];
  for (const label of labels) {
    const button = page.locator(".candidate").filter({hasText:label}).first();
    await button.click();
    await page.waitForFunction(expected => document.querySelector("#detail-name")?.textContent.includes(expected), label);
    await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("Mut · 槽位 1 已加载"), {timeout:30000});
    if (!await page.locator("#viewer canvas").count()) throw new Error("3D canvas missing for " + label);
    const counts = await page.locator("#file-counts").textContent();
    const expectedWt = label === "Q117R/E155W" ? "3/5" : label === "D78A/M122W" ? "1/5" : "5/5";
    if (!counts.includes("WT 文件 " + expectedWt + " 唯一")) throw new Error("WT content count mismatch for " + label + ": " + counts);
    wtContentCounts.push(counts.match(/WT 文件 ([0-9]+\/5) 唯一/)?.[1]);
  }
  await page.locator("#slot-select").selectOption("4");
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("槽位 5 已加载"), {timeout:30000});
  await page.locator('[data-kind="wt"]').click();
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("WT · 槽位 5 已加载"), {timeout:30000});
  await page.locator('[data-representation="surface"]').click();
  await page.waitForFunction(() => document.querySelector("#status")?.textContent.includes("WT · 槽位 5 已加载"), {timeout:60000});
  const result = await page.evaluate(() => ({
    title:document.title,
    badge:document.querySelector("#identity-badge")?.textContent,
    candidates:document.querySelectorAll(".candidate").length,
    canvas:document.querySelector("#viewer canvas")?.width ?? 0,
    wtDiversityWarning:document.querySelector("#file-counts")?.textContent,
    status:document.querySelector("#status")?.textContent
  }));
  const failures = structureResponses.filter(item => item.status !== 200);
  console.log(JSON.stringify({homeHttp:home.status(),http:response.status(),...result,wtContentCounts,structureLoads:structureResponses.length,failedStructureLoads:failures.length,errors}));
  if (response.status() !== 200 || result.badge !== "六项身份核验通过" || result.candidates !== 6 || !result.canvas || failures.length || errors.length) process.exitCode = 1;
} finally {
  await browser.close();
}
