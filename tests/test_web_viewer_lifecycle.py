"""Viewer lifecycle checks use small JS doubles; final browser checks use real APIs."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_binder_partner_components_hide_only_context_and_keep_viewer(tmp_path):
    source = (ROOT / "web/renderer.js").read_text(encoding="utf-8")
    checks = r"""
const assert=require('node:assert/strict');
const representation=()=>({hidden:false,visible:true});
const original={cell:{state:{},obj:{label:'Polymer'}},representations:[representation()]},components=[original],focusCalls=[];
global.window={PDBeMolstarPlugin:class {
 constructor(){this.assemblyRef='assembly';this.plugin={dispose(){},canvas3d:{setProps(){}},managers:{
 structure:{hierarchy:{current:{structures:[{components}]},toggleVisibility(items,action){items.forEach(c=>{c.cell.state.isHidden=action==='hide';c.representations.forEach(r=>{r.hidden=action==='hide';r.visible=action!=='hide';});});}}},
 camera:{focusLoci:(loci)=>focusCalls.push(loci)}},builders:{structure:{
 tryCreateComponentFromExpression:async(ref,expression,key,options)=>{assert.equal(ref,'assembly');assert(expression.args['chain-test']);const c={cell:{state:{},obj:{label:options.label}},representations:[representation()]};components.push(c);return c;},
 representation:{addRepresentation:async()=>{}}}}};this.visual={select:async()=>{},reset:async()=>{throw Error('full structure reset');}};this.canvas={applySettings:async()=>{}};}
 async applyVisualStyles(){components.forEach(c=>{c.representations=[representation()];});}
 async render(){} getLociForParams(records){return records;}
}};
(async()=>{
 const instance=await MolstarView.init({innerHTML:''},'unused','cartoon','data_test\n#');
 await MolstarView.configureChains(['A'],['B'],'cartoon');
 assert.equal(original.cell.state.isHidden,true);assert.equal(components[2].cell.state.isHidden,true);
 await MolstarView.setRepresentation('molecular-surface');
 assert.equal(original.representations[0].visible,false);assert.equal(components[2].representations[0].visible,false);
 assert.equal(components[1].representations[0].visible,true);
 await MolstarView.setPartnerVisible(true);assert.equal(components[2].cell.state.isHidden,false);
 assert.deepEqual(focusCalls.at(-1),[{auth_asym_id:'A'},{auth_asym_id:'B'}]);
 await MolstarView.color([{auth_seq_id:7}],'#fff');
 await MolstarView.setRepresentation('cartoon');
 assert.equal(original.representations[0].visible,false);assert.equal(components[2].representations[0].visible,true);
 await MolstarView.setPartnerVisible(false);assert.equal(components[2].cell.state.isHidden,true);
 await MolstarView.setRepresentation('molecular-surface');
 assert.equal(components[2].representations[0].visible,false);
 assert.deepEqual(focusCalls.at(-1),[{auth_asym_id:'A'}]);assert.equal(window.__ecdMolstar,instance);
 MolstarView.destroy();
})().catch(e=>{console.error(e);process.exitCode=1;});
"""
    run_js(tmp_path, source + checks)


def run_js(tmp_path, source):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for viewer lifecycle checks")
    script = tmp_path / "viewer.cjs"
    script.write_text(source, encoding="utf-8")
    subprocess.run([node, str(script)], check=True, capture_output=True, text=True, timeout=30)


def test_leaving_interface_step_hides_partner_before_mutant_selection(tmp_path):
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    function = app.split("function setStep(n)", 1)[1].split("function renderPanelScaffold", 1)[0]
    checks = r"""
const assert=require('node:assert/strict');
const state={step:4,showPartner:true},STEP_LIST=Array(7),visibility=[],errors=[];
const document={body:{dataset:{}},querySelectorAll:()=>[]},$=()=>({scrollTop:10});
let interfaceRenders=0;
const renderSteps=()=>{},renderStep4=()=>interfaceRenders++,viewerError=e=>errors.push(e);
const MolstarView={setPartnerVisible:async visible=>visibility.push(visible)};
setStep(5);assert.equal(state.showPartner,false);assert.deepEqual(visibility,[false]);
for(const step of [1,2,3,5,6,7]){
 state.step=4;state.showPartner=true;setStep(step);
 assert.equal(state.showPartner,false);assert.equal(visibility.at(-1),false);
}
assert.equal(interfaceRenders,7);
setStep(4);assert.equal(state.showPartner,false);
state.showPartner=true;setStep(4);assert.equal(state.showPartner,true);
assert.equal(visibility.length,7);
for(const step of [0,8,2.5]){setStep(step);assert.equal(state.step,4);assert.equal(state.showPartner,true);}
assert.deepEqual(errors,[]);
"""
    run_js(tmp_path, "function setStep(n)" + function + checks)


def test_mutant_focus_rebuilds_representation_before_painting(tmp_path):
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    function = app.split("async function selectDesign(id)", 1)[1].split("/* ---------------- persistent sequence dock", 1)[0]
    checks = r"""
const assert=require('node:assert/strict');
const state={designFocusRequest:0,layers:new Set(),engine:'molstar'},calls=[];
let molstarReady=true;
const selectedDesign=()=>({mutations:[{position:1}]}),regionForDesign=()=>({id:1});
const renderDesignSelect=()=>{},$=()=>({}),renderStep2=()=>{},renderStep5=()=>{},renderSeqDock=()=>{},renderDesignSurfaceCallout=()=>{};
const repaintViews=()=>{throw Error('paint started before representations finished rebuilding');};
const setRepresentation=async()=>{calls.push('representation and paint');};
const designSurfaceRecords=()=>({envelope:[],mutations:[{auth_seq_id:1}]}),residueQuery=x=>x;
const MolstarView={focus:async()=>calls.push('focus')},setEngine=()=>{};
(async()=>{await selectDesign('mutant');assert.deepEqual(calls,['representation and paint','focus']);})().catch(e=>{console.error(e);process.exitCode=1;});
"""
    run_js(tmp_path, "async function selectDesign(id)" + function + checks)


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


def test_recommendation_focus_reuses_viewer_and_transitions_with_reduced_motion(tmp_path):
    source = (ROOT / "web/renderer.js").read_text(encoding="utf-8")
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    checks = r"""
const assert=require('node:assert/strict');
const instances=[],focusCalls=[],selectCalls=[];let reducedMotion=false;
global.window={matchMedia:()=>({matches:reducedMotion}),PDBeMolstarPlugin:class {
  constructor(){instances.push(this);this.assemblyRef='same-assembly';this.disposals=0;
    this.plugin={dispose:()=>this.disposals++,canvas3d:{setProps(){}},managers:{
      structure:{hierarchy:{current:{structures:[{}]}}},camera:{focusLoci:(loci,options)=>focusCalls.push({loci,options})}}};
    this.visual={select:async options=>selectCalls.push(options)};this.canvas={applySettings:async()=>{}};
  }
  async render(){}
  getLociForParams(records){return {records,structure:'same-assembly'};}
}};
global.fetch=async()=>({ok:true,text:async()=>'data_test\n#'});
(async()=>{
 const host={innerHTML:''};
 const instance=await MolstarView.init(host,'structure.cif','cartoon');
 const canvas=instance.plugin.canvas3d;
 await MolstarView.color([{auth_seq_id:7}], '#c8d2d5', {keepRepresentations:true});
 await MolstarView.focus([{auth_seq_id:7}]);
 assert.equal(focusCalls[0].options.durationMs,600);
 assert.deepEqual(focusCalls[0].loci.records,[{auth_seq_id:7}]);
 reducedMotion=true;await MolstarView.focus([{auth_seq_id:8}]);
 assert.equal(focusCalls[1].options.durationMs,0);
 assert.equal(selectCalls[0].keepRepresentations,true);
 assert.equal(instances.length,1);assert.equal(instance.plugin.canvas3d,canvas);assert.equal(instance.disposals,0);
 assert.equal(typeof instance.visual.focus,'undefined');
})().catch(e=>{console.error(e);process.exitCode=1;});
"""
    run_js(tmp_path, source + checks)
    selection_paint = app.split("if (state.recommendationSelections.size)", 1)[1].split("const showcase", 1)[0]
    assert "keepRepresentations: true" in selection_paint
    assert "representation:" not in selection_paint
    assert "ensureViewers(" not in selection_paint


def test_adoption_initializes_once_and_refresh_avoids_init(tmp_path):
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    adoption = app.split("function adoptResult", 1)[1].split("function taskName", 1)[0]
    checks = r"""
const assert=require('node:assert/strict');
const state={previewPlans:new Map(),designFocusRequest:0,resultJobId:'OLD',recommendationSelections:new Map([['OLD-CANDIDATE',{}]])};let renders=[],loads=[];
const renderAll=value=>renders.push(value),ensureViewers=url=>loads.push(url),regionForDesign=()=>null;
adoptResult({residues:[],designs:[],protein:{ecd_coding_sequence:'AAA'},structure_url:'real.cif'},{demo:false,jobId:'NEW'});
assert.deepEqual(renders,[false]);assert.deepEqual(loads,['real.cif']);assert.equal(state.jobId,'NEW');assert.equal(state.resultJobId,'NEW');assert.equal(state.recommendationSelections.size,0);
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
const t=(key,fallback)=>fallback,toasts=[],showToast=s=>toasts.push(s),setEngine=()=>{},interfacePartnerChains=()=>[];
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
fields['sector-target'].value='';
assert.deepEqual(displayPartitionInput(),{maximum_sector_geodesic_diameter_a:30,maximum_sector_ca_diameter_a:22});
fields['sector-target'].value='8';
assert.equal(displayPartitionInput().target_sector_residues,8);
fields['sector-target'].value='7';
assert.equal(displayPartitionInput().target_sector_residues,7);
for(const [id,value] of [['sector-target','6'],['sector-target','257'],['sector-target','4.5'],['sector-geodesic','0'],['sector-ca','bad']]){
 const previous=fields[id].value;fields[id].value=value;assert.throws(displayPartitionInput);fields[id].value=previous;
}
"""
    run_js(tmp_path, "function displayPartitionInput" + function + checks)
