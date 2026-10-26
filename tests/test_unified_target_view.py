"""Both analysis targets use the same residue-scoped viewer and partner controls."""
from pathlib import Path
import json
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run_js(tmp_path, source):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for viewer checks")
    script = tmp_path / "unified.cjs"
    script.write_text(source, encoding="utf-8")
    subprocess.run([node, str(script)], check=True, capture_output=True, text=True, timeout=30)


@pytest.mark.parametrize("target", ["antigen", "binder"])
def test_target_adoption_resets_partner_and_reuses_shared_loader(tmp_path, target):
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    adoption = "function adoptResult" + app.split("function adoptResult", 1)[1].split("function taskName", 1)[0]
    run_js(tmp_path, f"const target='{target}';\n" + adoption + r"""
const assert=require('node:assert/strict');
const state={previewPlans:new Map(),designFocusRequest:0},loads=[];
const renderAll=()=>{},ensureViewers=url=>loads.push(url);
let partners=['B'];const interfacePartnerChains=()=>partners;
const result={analysis_target:target,residues:[],protein:{ecd_coding_sequence:'AAA'},structure_url:'/output/job/ecd_structure.cif'};
adoptResult(result,{demo:false,jobId:'job'});
assert.equal(loads.pop(),'/output/job/ecd_structure.cif');assert.equal(state.showPartner,false);
partners=[];adoptResult(result,{demo:false,jobId:'job'});
assert.equal(loads.pop(),'/output/job/ecd_structure.cif');
""")


@pytest.mark.parametrize("target", ["antigen", "binder"])
def test_viewer_initialization_scopes_both_targets_to_exact_residues(tmp_path, target):
    app = (ROOT / "web/app.js").read_text(encoding="utf-8")
    function = "async function ensureViewers" + app.split("async function ensureViewers", 1)[1].split("async function installMolstarTooltips", 1)[0]
    run_js(tmp_path, f"const target='{target}';\n" + function + r"""
const assert=require('node:assert/strict');
const residues=[{chain_id:'A',auth_seq_id:51,insertion_code:'B'}],calls=[];
const state={result:{analysis_target:target,residues},engine:'molstar',representation:'cartoon',offline:false};
let viewerLoadId=0,molstarReady=false;
const $=()=>({classList:{add(){}}}),BeadsView={destroy(){}};
const interfacePartnerChains=()=>['B'],residueQuery=r=>({auth_asym_id:r.chain_id,auth_seq_id:r.auth_seq_id,pdbx_PDB_ins_code:r.insertion_code});
const urls=[],MolstarView={init:async(host,url)=>{urls.push(url);return {};},configureChains:async(...args)=>calls.push(args),focus:async()=>{}};
const paintMolstar=async()=>{},designSurfaceRecords=()=>({envelope:[]}),installMolstarTooltips=async()=>{};
const viewerError=e=>{throw e;};
(async()=>{await ensureViewers('/output/job/ecd_structure.cif');assert.equal(molstarReady,true);
await ensureViewers('/output/job/ecd_structure.cif');
assert.deepEqual(urls,['/output/job/structure.cif','/output/job/structure.cif']);
assert.deepEqual(calls[0],[['A'],['B'],'cartoon',[{auth_asym_id:'A',auth_seq_id:51,pdbx_PDB_ins_code:'B'}]]);
})().catch(e=>{console.error(e);process.exitCode=1;});
""")


@pytest.mark.parametrize("partners", [[], ["B"]])
def test_exact_target_bundle_survives_representation_reload_and_reset(tmp_path, partners):
    source = (ROOT / "web/renderer.js").read_text(encoding="utf-8")
    run_js(tmp_path, "const partners=" + json.dumps(partners) + ";\n" + source + r"""
const assert=require('node:assert/strict');
const records=[{auth_asym_id:'A',auth_seq_id:51,pdbx_PDB_ins_code:'B'}],components=[],focusCalls=[],bundles=[];
const make=label=>({label,hidden:false});let full;
global.window={PDBeMolstarPlugin:class {
 constructor(){this.assemblyRef='assembly';this.plugin={dispose(){},canvas3d:{setProps(){}},managers:{
 structure:{hierarchy:{current:{structures:[{components}]},toggleVisibility(items,action){items.forEach(c=>c.hidden=action==='hide');}}},
 camera:{focusLoci:(records)=>focusCalls.push(records)}},builders:{structure:{
 tryCreateComponent:async(ref,params)=>{assert.equal(params.type.name,'bundle');assert.deepEqual(params.type.params,{records});
   const c=make(params.label);components.push(c);return c;},
 tryCreateComponentFromExpression:async()=>{const c=make('Protein partner');components.push(c);return c;},
 representation:{addRepresentation:async()=>{}}}}};this.visual={update:async()=>{components.length=0;full=make('Polymer');components.push(full);}};this.canvas={applySettings:async()=>{}};}
 async render(){full=make('Polymer');components.push(full);}
 getBundle(records){bundles.push(records);return {records};} getLociForParams(records){return records;}
}};
(async()=>{
 await MolstarView.init({innerHTML:''},'unused','cartoon','data_test\n#');
 await MolstarView.configureChains(['A'],partners,'cartoon',records);
 assert(full.hidden);assert.equal(components[1].label,'Analysis target');assert.equal(components[1].hidden,false);
 await MolstarView.setPartnerVisible(true);
 await MolstarView.setRepresentation('molecular-surface');
 assert(full.hidden);assert.equal(bundles.length,2);assert.equal(components[1].hidden,false);
 if(partners.length)assert.equal(components[2].hidden,false);
 await MolstarView.reset();assert.deepEqual(focusCalls.at(-1),[...records,...partners.map(auth_asym_id=>({auth_asym_id}))]);
 await MolstarView.setPartnerVisible(false);await MolstarView.reset();assert.deepEqual(focusCalls.at(-1),records);
 if(partners.length)assert.equal(components[2].hidden,true);
})().catch(e=>{console.error(e);process.exitCode=1;});
""")
