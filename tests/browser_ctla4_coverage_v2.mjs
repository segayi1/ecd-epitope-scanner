/* Read-only real-result check. Every API write is rejected and counted. */
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import http from 'node:http';
import {pathToFileURL} from 'node:url';
const result=path.resolve(process.argv[2]);
const bytes=fs.readFileSync(path.join(result,'analysis.json'));
const analysis=JSON.parse(bytes);
const physical=analysis.design_summary?.surface_perturbation;
assert.ok(physical,'Actual API schema must contain physical summary');
const playwright=await import(pathToFileURL('C:/Users/陈淏文/.workbuddy/binaries/node/workspace/node_modules/playwright-core/index.js').href);
const browser=await (playwright.chromium||playwright.default.chromium).launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true});
const page=await browser.newPage();const writes=[],errors=[];const web=path.resolve('web');
page.on('pageerror',e=>errors.push(e.message));
const server=http.createServer((request,response)=>{
  const url=new URL(request.url,'http://127.0.0.1');
  if(request.method!=='GET'){writes.push(url.pathname);response.writeHead(405);response.end('Read only');return;}
  if(url.pathname==='/assets/demo/analysis.json'){response.writeHead(200,{'Content-Type':'application/json'});response.end(bytes);return;}
  if(url.pathname==='/assets/demo/ecd_structure.cif'){
    response.writeHead(200,{'Content-Type':'text/plain'});fs.createReadStream(path.join(result,'ecd_structure.cif')).pipe(response);return;
  }
  const relative=url.pathname==='/'?'index.html':url.pathname.replace(/^\/assets\//,'');
  const file=path.resolve(web,relative);
  if(file.startsWith(web+path.sep)&&fs.existsSync(file)&&fs.statSync(file).isFile()){
    const mime={'.html':'text/html','.js':'application/javascript','.json':'application/json','.css':'text/css'}[path.extname(file)]||'application/octet-stream';
    response.writeHead(200,{'Content-Type':mime});fs.createReadStream(file).pipe(response);
  }else {response.writeHead(404);response.end('missing');}
});
await new Promise(resolve=>server.listen(0,'127.0.0.1',resolve));
try{
  await page.goto(`http://127.0.0.1:${server.address().port}/`,{waitUntil:'domcontentloaded'});
  await page.waitForFunction(status=>document.querySelector('#task-overview')?.textContent.includes(status),physical.status,{timeout:60000});
  const text=await page.locator('#task-overview').innerText();
  assert.ok(text.includes(physical.goal_status));
  assert.ok(text.includes(physical.search_status));
  assert.doesNotMatch(text,/NaN/);assert.deepEqual(writes,[]);assert.deepEqual(errors,[]);
  const report={status:'PASS',analysis_sha256:crypto.createHash('sha256').update(bytes).digest('hex'),
    physical_status:physical.status,goal_status:physical.goal_status,search_status:physical.search_status,
    api_writes:writes,errors,interpretation:'Browser displays actual scientific state; rendering PASS is not scientific acceptance'};
  if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(report,null,2));
  console.log(JSON.stringify(report));
}catch(error){
  console.error(JSON.stringify({status:'FAIL',error:String(error),errors,writes,
    rendered_text:(await page.locator('body').innerText().catch(()=>'' )).slice(-12000)}));
  throw error;
}finally{await browser.close();await new Promise(resolve=>server.close(resolve));}
