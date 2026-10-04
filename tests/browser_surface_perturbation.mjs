/* Browser submission contract; API writes are intercepted, no real job created. */
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
const playwright = await import(pathToFileURL('C:/Users/陈淏文/.workbuddy/binaries/node/workspace/node_modules/playwright-core/index.js').href);
const browser=await (playwright.chromium||playwright.default.chromium).launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true});
const page=await browser.newPage();const errors=[],requests=[];
page.on('pageerror',error=>errors.push(error.message));
const web=path.resolve('web');
await page.route('http://surface-preview.test/**',async route=>{
  const url=new URL(route.request().url());
  if(url.pathname==='/assets/demo/analysis.json'){
    const fixture=JSON.parse(fs.readFileSync(path.join(web,'demo','analysis.json'),'utf8'));
    fixture.joint_v4_ml={...(fixture.joint_v4_ml||{}),surface_perturbation:{
      status:'BROWSER_FIXTURE_NO_PHYSICAL_EVIDENCE',actual_constructs:0,
      panel_contract:'combo_2_3_hard80_v1',goal_status:'COVERAGE_BELOW_TARGET',search_status:'CALIBRATION_EXPANSION_REQUIRED',
      verified_local_constructs:[{candidate_id:'fixture_combo',mutation_notations:['T82W','Q117R'],model_repeats:5,
        site_evidence:{complete:true,sites:[{notation:'T82W',verified_direct_area_a2:4},{notation:'Q117R',verified_direct_area_a2:6}]},
        spurs_ddg_kcal_mol:0,foldx_ddg_kcal_mol:0,actual_changed_area_a2:10,
        hard_scan_center_area_a2:20,interpretation:'BROWSER_FIXTURE_ONLY'}],
      coverage_curve:Array.from({length:12},(_,i)=>({K:i+1,actual_constructs:0,
        C_hard_all:0,C_hard_designable:0,C_extended_all:null,C_extended_designable:null,
        largest_blind_patch_a2:null}))}};
    // Match the native pipeline/API location, not an invented viewer schema.
    fixture.design_summary={...(fixture.design_summary||{}),surface_perturbation:fixture.joint_v4_ml.surface_perturbation};
    delete fixture.joint_v4_ml.surface_perturbation;
    await route.fulfill({status:200,json:fixture});return;
  }
  if(url.pathname==='/api/jobs' && route.request().method()==='POST'){
    requests.push(route.request().postDataJSON());
    await route.fulfill({status:202,json:{id:'intercepted',status:'FAILED',message:'UI-only',progress:0}});return;
  }
  const relative=url.pathname==='/'?'index.html':url.pathname.replace(/^\/assets\//,'');
  const file=path.resolve(web,relative);
  if(file.startsWith(web+path.sep)&&fs.existsSync(file)&&fs.statSync(file).isFile()){
    const mime={'.html':'text/html','.js':'application/javascript','.json':'application/json','.css':'text/css'}[path.extname(file)]||'application/octet-stream';
    await route.fulfill({status:200,contentType:mime,body:fs.readFileSync(file)});
  }else await route.fulfill({status:404,body:'missing'});
});
try{
  await page.goto('http://surface-preview.test/',{waitUntil:'domcontentloaded'});
  await page.waitForTimeout(1500);
  await page.waitForFunction(()=>document.querySelector('#task-overview')?.textContent.includes('BROWSER_FIXTURE_NO_PHYSICAL_EVIDENCE'));
  assert.equal(await page.locator('#task-overview tbody tr').count(),13);
  assert.match(await page.locator('#task-overview').innerText(),/局部构建体/);
  assert.match(await page.locator('#task-overview').innerText(),/COVERAGE_BELOW_TARGET/);
  assert.match(await page.locator('#task-overview').innerText(),/T82W/);
  assert.match(await page.locator('#task-overview').innerText(),/Q117R/);
  assert.match(await page.locator('#task-overview').innerText(),/逐位直接变化/);
  assert.match(await page.locator('#task-overview').innerText(),/未知/);
  assert.doesNotMatch(await page.locator('#task-overview').innerText(),/NaN/);
  assert.equal(await page.locator('#surface-strategy').inputValue(),'legacy');
  await page.locator('#surface-strategy').selectOption('surface_perturbation_v1_preview');
  await page.locator('summary').first().click();
  await page.locator('#budget-maximum').fill('12');
  await page.locator('#run-button').click();
  await page.waitForTimeout(300);
  assert.equal(requests.length,1);assert.equal(requests[0].surface_strategy,'surface_perturbation_v1_preview');
  assert.equal(requests[0].maximum_constructs,12);
  await page.locator('#surface-strategy').selectOption('legacy');
  await page.locator('#upload-surface-strategy').selectOption('surface_perturbation_v1_preview',{force:true});
  await page.locator('#run-button').click();await page.waitForTimeout(200);
  assert.equal(requests.length,1,'legacy must reject 12 locally');
  await page.locator('#budget-maximum').fill('8');
  await page.locator('#run-button').click();await page.waitForTimeout(200);
  assert.equal(requests.length,2,'legacy accepts its independent ceiling');
  assert.equal(requests[1].surface_strategy,'legacy');
  assert.equal(requests[1].maximum_constructs,8);
  await page.locator('#optimization-profile').selectOption('coverage_scan_v2_preview');
  await page.locator('#budget-maximum').fill('12');
  await page.locator('#run-button').click();await page.waitForTimeout(200);
  assert.equal(requests.length,3);
  assert.equal(requests[2].optimization_profile,'coverage_scan_v2_preview');
  assert.equal(requests[2].maximum_constructs,12);
  assert.deepEqual(errors,[]);
  console.log(JSON.stringify({status:'PASS',intercepted:requests.length,errors}));
}finally{await browser.close();}
