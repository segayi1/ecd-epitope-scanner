"""Viewer lifecycle checks use small JS doubles; final browser checks use real APIs."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run_js(tmp_path, source):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for viewer lifecycle checks")
    script = tmp_path / "viewer.cjs"
    script.write_text(source, encoding="utf-8")
    subprocess.run([node, str(script)], check=True, capture_output=True, text=True, timeout=30)


def test_renderer_disposes_once_and_ignores_old_load(tmp_path):
    source = (ROOT / "web/renderer.js").read_text(encoding="utf-8")
    checks = r"""
const assert = require('node:assert/strict');
const instances=[];let releaseOld;
global.window={PDBeMolstarPlugin:class {
  constructor(){instances.push(this);this.disposals=0;}
  async render(){this.plugin={dispose:()=>this.disposals++,canvas3d:{setProps(){}},managers:{structure:{hierarchy:{current:{structures:[{}]}}}}};this.assemblyRef='loaded';
    if(instances.indexOf(this)===0)await new Promise(resolve=>releaseOld=resolve);
  }
}};
global.fetch=async url=>({ok:url!=='bad-http',status:404,text:async()=>url==='bad-format'?'not CIF':'# generated structure\ndata_test\n#'});
(async()=>{
  const host={innerHTML:''};const old=MolstarView.init(host,'old','cartoon');
  await new Promise(resolve=>setTimeout(resolve,0));
  const current=await MolstarView.init(host,'new','cartoon');
  assert.equal(instances.length,2);assert.equal(instances[0].disposals,1);
  releaseOld();assert.equal(await old,null);
  assert.equal(instances[0].disposals,1);assert.equal(window.__ecdMolstar,current);assert(MolstarView.isReady());
  MolstarView.destroy();assert.equal(current.disposals,1);assert(!MolstarView.isReady());assert.equal(window.__ecdMolstar,null);
  await assert.rejects(MolstarView.init(host,'bad-http'),e=>e.stage==='fetch'&&e.message==='HTTP 404');
  await assert.rejects(MolstarView.init(host,'bad-format'),e=>e.stage==='parse');
})().catch(e=>{console.error(e);process.exitCode=1;});
"""
    run_js(tmp_path, source + checks)


def test_adoption_initializes_once_and_refresh_avoids_init(tmp_path):
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    adoption = app.split("function adoptResult", 1)[1].split("function taskName", 1)[0]
    checks = r"""
const assert=require('node:assert/strict');
const state={previewPlans:new Map(),designFocusRequest:0};let renders=[],loads=[];
const renderAll=value=>renders.push(value),ensureViewers=url=>loads.push(url),regionForDesign=()=>null;
adoptResult({residues:[],designs:[],protein:{ecd_coding_sequence:'AAA'},structure_url:'real.cif'},{demo:false,jobId:'NEW'});
assert.deepEqual(renders,[false]);assert.deepEqual(loads,['real.cif']);assert.equal(state.jobId,'NEW');
"""
    run_js(tmp_path, "function adoptResult" + adoption + checks)
    render_all = app.split("function renderAll", 1)[1].split("function renderStep1", 1)[0]
    assert "ensureViewers" not in render_all
    assert "repaintViews();" in render_all
    assert '?view=display`' in app


def test_stale_app_load_and_optional_tooltip_failure(tmp_path):
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    lifecycle = app.split("let molstarReady = false;", 1)[1].split("async function installMolstarTooltips", 1)[0]
    checks = r"""
const assert=require('node:assert/strict');let paints=0,focuses=0,releaseOld,initCount=0;
const empty={classList:{add(){this.hidden=true;},remove(){this.hidden=false;}},querySelector(){return {textContent:''};}};
const buttons={'viewport-empty':empty,'retry-viewer':{},'load-demo':{}};const $=id=>buttons[id];
const state={engine:'molstar',jobId:'new',result:{residues:[],structure_url:'new.cif'}};
const t=(key,fallback)=>fallback,toasts=[],showToast=s=>toasts.push(s),setEngine=()=>{};
const BeadsView={destroy(){}};
const MolstarView={init:async()=>{if(++initCount===1)return new Promise(r=>releaseOld=r);return {};},focus:async()=>focuses++};
const paintMolstar=async()=>paints++,selectedDesign=()=>null,designSurfaceRecords=()=>({envelope:[]}),residueQuery=x=>x;
const installMolstarTooltips=async()=>{throw new Error('optional tooltip failed');};
(async()=>{
 const old=ensureViewers('old.cif');await ensureViewers('new.cif');releaseOld({});await old;
 assert.equal(paints,1);assert.equal(focuses,1);assert(molstarReady);assert(empty.classList.hidden);assert.equal(toasts.length,1);assert(toasts[0].includes('安装位点提示'));
 MolstarView.init=async()=>{throw Object.assign(new Error('HTTP 404'),{stage:'fetch'});};
 await ensureViewers('missing.cif');assert(!molstarReady);assert(!empty.classList.hidden);assert.equal(buttons['retry-viewer'].hidden,false);
 MolstarView.init=async()=>({});await ensureViewers('missing.cif');assert(molstarReady);assert(empty.classList.hidden);assert.equal(buttons['retry-viewer'].hidden,true);
})().catch(e=>{console.error(e);process.exitCode=1;});
"""
    run_js(tmp_path, "let molstarReady = false;" + lifecycle + checks)


def test_viewer_and_partition_translations():
    keys = {"hud.retryViewer", "hud.stageFetch", "hud.stageParse", "hud.stageRender", "hud.stageStructure", "hud.stageWebgl", "hud.stagePaint", "hud.stageFocus", "hud.stageTooltip", "form.partitionSettings", "form.partitionTarget", "form.partitionSurface", "form.partitionCa", "form.partitionInvalid"}
    for language in ("zh", "en"):
        payload = json.loads((ROOT / "web/i18n" / f"{language}.json").read_text(encoding="utf-8"))
        assert keys <= payload.keys()


def test_display_partition_inputs_are_validated(tmp_path):
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    function = app.split("function displayPartitionInput", 1)[1].split("async function submitRun", 1)[0]
    checks = r"""
const assert=require('node:assert/strict');
const fields={'sector-target':{value:'12'},'sector-geodesic':{value:'30'},'sector-ca':{value:'22'}};
const $=id=>fields[id],t=(key,fallback)=>fallback;
assert.deepEqual(displayPartitionInput(),{target_sector_residues:12,maximum_sector_geodesic_diameter_a:30,maximum_sector_ca_diameter_a:22});
for(const [id,value] of [['sector-target','3'],['sector-target','257'],['sector-target','4.5'],['sector-geodesic','0'],['sector-ca','bad']]){
 const previous=fields[id].value;fields[id].value=value;assert.throws(displayPartitionInput);fields[id].value=previous;
}
"""
    run_js(tmp_path, "function displayPartitionInput" + function + checks)
