import assert from "node:assert/strict";
import fs from "node:fs/promises";
import {pathToFileURL} from "node:url";
const pw = await import(pathToFileURL("C:/Users/陈淏文/.workbuddy/binaries/node/workspace/node_modules/playwright-core/index.js").href);
const chromium = pw.chromium || pw.default?.chromium;
const runDir = process.argv[2];
const base = process.argv[3] || "http://127.0.0.1:8765";
if (!runDir) throw new Error("Provide a validated V3 run directory");
const analysis = JSON.parse(await fs.readFile(`${runDir}/analysis.json`,"utf8"));
const browser = await chromium.launch({executablePath:"C:/Program Files/Google/Chrome/Application/chrome.exe",headless:true});
const page = await browser.newPage({viewport:{width:1600,height:1000}});
const errors=[];
page.on("pageerror",e=>errors.push(e.message));
try {
  // Read-only browser fixture: expose validated files through intercepted
  // responses. No historical job or production recommendation is modified.
  await page.route("**/api/jobs/v411-browser-validation",route=>route.fulfill({json:{
    id:"v411-browser-validation",status:"COMPLETE",stage:"COMPLETE",terminal:true,
    progress:1,message:"V3 validation",result:analysis,
  }}));
  await page.route("**/api/jobs/v411-browser-validation/result",route=>route.fulfill({json:analysis}));
  await page.route("**/api/jobs/v411-browser-validation/files/**",async route=>{
    const name=decodeURIComponent(new URL(route.request().url()).pathname.split("/").at(-1));
    if (!name || name.includes("..")) return route.abort();
    try {await route.fulfill({body:await fs.readFile(`${runDir}/${name}`)});} catch {await route.abort();}
  });
  await page.goto(`${base}/?job=v411-browser-validation`,{waitUntil:"domcontentloaded"});
  await page.waitForSelector("#residue-filter",{timeout:30000});
  const budgetControls=await page.locator("#adaptive-budget").isChecked();
  assert.equal(budgetControls,true);
  // Switch to combination review and assert its independently audited count.
  const buttons=page.locator("[data-step='5']");
  if (await buttons.count()) await buttons.click();
  const budgetText=await page.locator("#design-list").innerText();
  assert(budgetText.includes(`入选 ${analysis.design_summary.panel_budget_decision.selected_constructs} 个`));
  await page.evaluate(()=>document.querySelector("[data-step='1']")?.click());
  await page.locator("#residue-filter").selectOption("all");
  await page.waitForTimeout(1000);
  assert.equal(errors.length,0,errors.join("\n"));
  await page.screenshot({path:`${runDir}/browser_v411.png`,fullPage:true});
  await fs.writeFile(`${runDir}/browser_v411_validation.json`,JSON.stringify({budgetControls,budgetText,errors,read_only_fixture:true},null,2));
  console.log(JSON.stringify({budgetControls,budgetText,errors}));
} finally {await browser.close();}
