"""Frozen CTLA4 staged validation with one charged case ledger and owned storage.

Development/manual gaps between invocations are excluded by closed intervals.
Every launched scientific interval, including failures, is charged.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ecd_epitope_scanner.pipeline import PipelineRequest, run_pipeline
from ecd_epitope_scanner.retrieval import parse_uniprot_entry
from ecd_epitope_scanner.scripts.replay_joint_v41_run import replay_request_data
from ecd_epitope_scanner.scripts.surface_phase_guard import stop_owned_tree
from ecd_epitope_scanner.surface_perturbation.case_resources import (
    GIB, EvidenceStore, StorageGuard, affordable_batch, owned_bytes,
)
from ecd_epitope_scanner.surface_perturbation.runtime import RunBudget, read_budget
from ecd_epitope_scanner.util import sha256_file, write_json


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def physical_stage_engineering_status(physical):
    if not isinstance(physical,dict) or not physical.get('status') or physical['status']=='INCOMPLETE':
        return 'ENGINEERING_FAILURE'
    return 'COMPLETED'


def stage_attempt_number(case,prefix,phases=('calibration','production')):
    index=1
    while True:
        stems=[f'{prefix}_{index}'+('_'+phase if phase else '') for phase in phases]
        occupied=any((case/'stages'/(stem+suffix)).exists()
            for stem in stems for suffix in ('_request.json','_source_manifest.json','_report.json','_result.json','.log'))
        if not occupied:return index
        index+=1


def case_limit(case):
    amendment=case/'resource_contract_amendment.json'
    return None if amendment.is_file() and read(amendment).get('time_budget')=='UNLIMITED' else 14400


def case_remaining(case):
    limit=case_limit(case)
    return math.inf if limit is None else limit-read_budget(case/'case_budget.json')['elapsed_seconds']


def case_storage_limit(case):
    amendment=case/'resource_contract_amendment.json'
    return int(read(amendment)['storage_byte_budget']) if amendment.is_file() else 20*GIB


def activate_resource_amendment(case):
    """Apply the explicit user amendment only at a closed, locked stage boundary."""
    amendment_path=case/'resource_contract_amendment.json'
    amendment=read(amendment_path)
    if amendment.get('time_budget')!='UNLIMITED' or not amendment.get('authorization'):
        raise ValueError('Explicit unlimited-time authorization is missing')
    old=read_budget(case/'case_budget.json')
    if old.get('active'):raise ValueError('Amendment requires a closed stage boundary')
    if amendment.get('storage_policy')=='AVAILABLE_DISK_WITH_HEADROOM':
        import shutil
        reserve=10*GIB
        amendment['minimum_free_disk_bytes']=reserve
        amendment['storage_byte_budget']=owned_bytes(case)+max(0,shutil.disk_usage(case).free-reserve)
    budget=RunBudget(case/'case_budget.json',None)
    try:
        if old.get('limit_seconds') is not None:
            snapshot=case/'case_budget_before_unlimited_amendment.json'
            if snapshot.exists():raise ValueError('Historical amendment ledger already exists')
            write_json(snapshot,old)
        budget.limit=math.inf
    finally:budget.save(False)
    amendment.update(activation='ACTIVE',activated_wall_time=time.time(),
        elapsed_seconds_at_activation=read_budget(case/'case_budget.json')['elapsed_seconds'])
    write_json(amendment_path,amendment)


def worker(case, assets, config):
    source=assets/'workspace_data/8c068ad701a6'
    analysis, entry, old = (read(source/name) for name in ('analysis.json','uniprot.json','request.json'))
    output=Path(config.get('results_root',case/'results'));output.mkdir(parents=True,exist_ok=True)
    if case not in output.resolve().parents:raise ValueError('Results must remain inside the owned case')
    profile=config.get('optimization_profile','coverage_scan_v2_preview')
    write_json(output/'frozen_input_manifest.json',read(case/'T0_preflight.json'))
    snapshot=Path(__file__).resolve().parents[1]/'source_overlay_manifest.json'
    write_json(output/'launch_source_manifest.json',{'source_snapshot':read(snapshot) if snapshot.is_file() else None,
        'python':sys.executable,'configuration':config})
    data=replay_request_data(source,output,analysis,old)
    binaries=sorted((assets/'vendor/foldx').glob('foldx_*.exe'))
    # The installed Windows executable crashes before main with its long path;
    # a byte-identical short-named task-owned copy passed the actual startup probe.
    execution_binary=case/'fx/foldx.exe'
    tool_pin=case/'T0_frozen_tools.json'
    expected_binary=read(tool_pin if tool_pin.is_file() else case/'T0_preflight.json')['foldx_binary_sha256']
    if (not execution_binary.is_file() or sha256_file(execution_binary)!=expected_binary
            or sha256_file(binaries[-1])!=expected_binary):
        raise ValueError('Task-owned FoldX execution alias missing or changed')
    upstream={name:{'model_assets_root':str(assets)} for name in ('spurs_single','spurs_multi','discotope')}
    frozen_assets=case/'T0_frozen_model_assets.json'
    if frozen_assets.is_file():
        for name in ('spurs_single','spurs_multi'):
            upstream[name]['frozen_model_assets_path']=str(frozen_assets)
    upstream['runtime']={'cache_dir':str(case/'upstream_cache')}
    upstream['foldx_confirmation']={'binary':str(execution_binary),'engine_work_root':str(case/'fx'),
        'max_candidates':config['confirmation_count'],'paired_cache_roots':[str(case)]}
    upstream['candidate_search']={'global_shortlist':config.get('pool_cap',32),
        'triple_parent_ids':config.get('triple_parent_ids',[]),'triple_role':config.get('triple_role','TRIPLE_COMPRESSION'),
        'allow_remote_triples':config.get('allow_remote_triples',False),
        'mandatory_candidate_signatures':config.get('mandatory_candidate_signatures',[]),
        'conditional_pair_backups':config.get('conditional_pair_backups',[])}
    for key in ('triple_roles_by_candidate_id','triple_parents_by_candidate_id'):
        upstream['candidate_search'][key]=config.get(key,{})
    upstream['menu_supplements']=config.get('menu_supplements',[])
    upstream['full_replacement_search']=config.get('full_replacement_search',False)
    upstream['candidate_stage_limit']=config.get('stage_limit')
    upstream['surface_perturbation']={
        'optimization_profile':profile,'phase':config['phase'],
        'case_time_budget_seconds':config['seconds'], 'case_storage_root':str(case),
        'case_byte_budget':case_storage_limit(case),'cumulative_budget_path':str(case/'ledgers'/f"{config['stage']}.json"),
        'calibration_dir':str(output/'surface_perturbation/calibration'),
        'paired_cache_roots':[str(case)],'tool_cache':str(case/'physical_cache'),
        'paired_manifest_cache_roots':[str(case/'evidence_index')],
        'calibration_candidates':config['confirmation_count'],
        'confirmation_batch':config['confirmation_count'],
        'compression_baseline_ids':config.get('compression_baseline_ids'),
        'calibration_repeat_schedule':config.get('schedule',[5]),
        'calibration_extension_candidate_ids':config.get('calibration_extension_candidate_ids'),
        'maximum_iterations':2,
        'final_panel_audit':config.get('final_panel_audit',False),
    }
    data.update(optimization_profile=profile, assay_profile=None,
        max_primary_constructs=None, minimum_constructs=4,target_constructs=7,maximum_constructs=12,
        adaptive_budget=True, export_mode='checkpoint', candidate_ids=config.get('candidate_ids'),
        calibration_repeat_schedule=config.get('schedule',[5]),
        metadata={**data.get('metadata',{}),'joint_v4_ml':upstream})
    conservation={r['uniprot_pos']:r['conservation_score'] for r in analysis['residues'] if r.get('conservation_score') is not None}
    def protein(accession,client,output_dir,ecd_override=None):
        if accession!='P16410':raise ValueError('Frozen accession mismatch')
        write_json(output_dir/'uniprot.json',entry)
        return parse_uniprot_entry(entry,ecd_override)
    with patch('ecd_epitope_scanner.pipeline.fetch_protein_context',protein), patch(
        'ecd_epitope_scanner.conservation.fetch_conservation',return_value=conservation or None):
        result=run_pipeline(PipelineRequest(**data),lambda stage,value,message:print(stage,message,flush=True))
    physical=result['design_summary'].get('surface_perturbation')
    write_json(case/'stages'/f"{config['stage']}_result.json", {
        'engineering_status':physical_stage_engineering_status(physical),'physical':physical,
        'configuration':config})


def launch(case,assets,stage,config,seconds):
    (case/'stages').mkdir(exist_ok=True)
    budget=RunBudget(case/'case_budget.json',case_limit(case))
    # Write admissions retain fresh inventories; only the outer supervisor
    # amortizes whole-case traversal. Disk headroom is checked on every poll.
    storage=StorageGuard(case,inventory_interval=30 if case_limit(case) is None else 2)
    if stage.startswith('T3'):
        seconds=min(seconds,1200-budget.elapsed)
    elif stage.startswith(('T4','T5')):
        seconds=min(seconds,3600-budget.elapsed)
    config={**config,'stage':stage,'seconds':min(seconds,budget.remaining(),7200)}
    path=case/'stages'/f'{stage}_request.json';write_json(path,config)
    source_manifest=Path(__file__).resolve().parents[1]/'source_overlay_manifest.json'
    if source_manifest.is_file():
        write_json(case/'stages'/f'{stage}_source_manifest.json',read(source_manifest))
    process=None;start=time.monotonic();started_wall=time.time();before=owned_bytes(case);status='COMPLETED'
    try:
        storage.check(reserved=256*1024**2)
        with (case/'stages'/f'{stage}.log').open('w',encoding='utf-8') as log:
            process=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),
                '--case-root',str(case),'--assets-root',str(assets),'--worker',str(path)],
                cwd=Path(__file__).resolve().parents[2],stdout=log,stderr=subprocess.STDOUT)
            while process.poll() is None:
                storage.check()
                if time.monotonic()-start >= config['seconds'] or budget.remaining()<=0:
                    status='INCOMPLETE_RUNTIME_EXHAUSTED';stop_owned_tree(process);break
                time.sleep(.25)
            if process.returncode not in (0,None):status='ENGINEERING_FAILURE'
    except Exception as exc:
        status='ENGINEERING_FAILURE';config['failure']=str(exc)
        if process is not None and process.poll() is None:stop_owned_tree(process)
    finally:
        # Storage/checkpoint work belongs to this scientific stage as well.
        try:
            store=EvidenceStore(case)
            if config.get('optimization_profile')=='coverage_scan_v3_preview':
                store.checkpoint(stage,source_roots=[Path(config['results_root']),case/'fx'])
            else:store.checkpoint(stage)
            storage.check(solving=False,force_inventory=True)
        finally:
            budget.save(False)
    closed_elapsed=read_budget(case/'case_budget.json')['elapsed_seconds']
    report={'stage':stage,'engineering_status':status,'charged_seconds':closed_elapsed-budget.carried,
        'case_elapsed_seconds':closed_elapsed,'case_elapsed_before_stage':budget.carried,
        'worker_and_checkpoint_wall_seconds':time.monotonic()-start,
        'charge_basis':'persisted case ledger delta including admission and checkpoint',
        'new_bytes':max(0,owned_bytes(case)-before),
        'peak_owned_bytes':storage.peak,'configuration':config}
    summary=Path(config.get('results_root',case/'results'))/'surface_perturbation_preview.json'
    if summary.is_file() and summary.stat().st_mtime>=started_wall:
        report['physical']=read(summary)
        if report['engineering_status']=='COMPLETED':
            report['engineering_status']=physical_stage_engineering_status(report['physical'])
    else:
        report['physical_status']='NOT_PRODUCED_FOR_THIS_ATTEMPT'
        if report['engineering_status']=='COMPLETED':report['engineering_status']='ENGINEERING_FAILURE'
    if report.get('physical',{}).get('status')=='INCOMPLETE':
        report['engineering_status']='ENGINEERING_FAILURE'
    if report.get('physical',{}).get('errors'):
        failures=report['physical']['errors']
        if any('FoldX supplemental repeats failed' in row.get('reason','') for row in failures):
            report['engineering_status']='ENGINEERING_FAILURE'
    write_json(case/'stages'/f'{stage}_report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='physical'},ensure_ascii=False),flush=True)
    return report


def preflight(case,assets):
    import shutil
    source=assets/'workspace_data/8c068ad701a6'
    expected={'structure.cif':'c3550fdbd6a113891080ad7ead61076c9da8620af66fb417a715654ea5c589bc',
        'uniprot.json':'b9e0717afdf6c9c807b055a44368fb7a34bfc89c5f06f01292d1890db3db4ff9'}
    hashes={name:sha256_file(source/name) for name in expected}
    if hashes!=expected:raise ValueError('CTLA4 frozen input identity changed')
    if read(source/'analysis.json')['protein']['accession']!='P16410':raise ValueError('Wrong input protein')
    StorageGuard(case).check()
    binaries=sorted((assets/'vendor/foldx').glob('foldx_*.exe'))
    if not binaries:raise ValueError('FoldX binary missing')
    alias=case/'fx/foldx.exe';alias.parent.mkdir(exist_ok=True)
    digest=sha256_file(binaries[-1])
    if alias.exists() and sha256_file(alias)!=digest:raise ValueError('FoldX execution alias changed')
    if not alias.exists():shutil.copy2(binaries[-1],alias)
    probe=subprocess.run([str(alias),'--help'],cwd=alias.parent,capture_output=True,timeout=30)
    if probe.returncode:raise ValueError('FoldX startup preflight failed: '+str(probe.returncode))
    write_json(case/'T0_preflight.json',{'engineering_status':'PASS','scientific_computation':0,
        'input_type':'experimental','chain':'A','accession':'P16410','hashes':hashes,
        'assets_root':str(assets),'source_snapshot':str(Path(__file__).resolve().parents[1]),
        'foldx_source_binary':str(binaries[-1]),'foldx_execution_binary':str(alias),
        'foldx_binary_sha256':digest,'foldx_startup_exit_code':probe.returncode,
        'case_time_budget_seconds':14400,'case_byte_budget':20*GIB,
        'contract':{'target_hard_coverage':.8,'scan_radius_a':12,'compression_loss':.01,
                    'pool_rule':'clamp(4*affordable_groups,32,256)','batch_max':8}})


def replay_saved_models(case,assets):
    from Bio.PDB import MMCIFParser, PDBIO
    from dataclasses import fields
    from ecd_epitope_scanner.models import ResidueRecord
    from ecd_epitope_scanner.surface_perturbation.paired import scan_pairs, validate_pair, atoms
    from ecd_epitope_scanner.surface_perturbation.fields import aligned_atoms
    fixtures=[('positive_double','ctla4_v41_repair_20260930/verified_pipeline','foldx_batch_00f0dfc2',2),
        ('positive_triple','v411_upgrade_20260930/ctla4_v3_strict','foldx_batch_4735fa0b',5),
        ('disagreement_triple','ctla4_v41_repair_20260930/verified_pipeline','foldx_batch_15714a7d',6)]
    reports=[];raw=[]
    names={f.name for f in fields(ResidueRecord)}
    native_path=case/'T2_native_author.pdb'
    writer=PDBIO();writer.set_structure(MMCIFParser(QUIET=True).get_structure('native',
        str(assets/'workspace_data/8c068ad701a6/structure.cif')))
    writer.save(str(native_path))
    for label,relative,batch,group in fixtures:
        root=assets/'workspace_data'/relative
        all_pairs,failures=scan_pairs(root)
        if failures:raise ValueError('Saved pair identity failures: '+str(failures[:3]))
        pairs=[p for p in all_pairs if Path(p['batch_path']).name==batch and p['group']==group]
        if len(pairs)!=5:raise ValueError('Saved fixture requires five real pairs')
        records=[ResidueRecord(**{k:v for k,v in row.items() if k in names}) for row in read(root/'analysis.json')['residues']]
        native=atoms(native_path)
        audits=[]
        for pair in pairs:
            validate_pair(pair)
            _,_,audit=aligned_atoms(pair,native,records)
            if audit['core_rmsd_a']>.5:raise ValueError('Saved core registration failed')
            audits.append(audit)
        manifest=read(root/'v4_manifest.json')
        binary_hash=manifest['tools']['foldx']['weights_sha256']
        for pair in pairs:
            pair.update(binary_sha256=binary_hash,normalized_input_sha256=sha256_file(root/'prediction_context.pdb'))
        raw.append({'mutations':pairs[0]['mutations'],'paired_models':pairs})
        reports.append({'fixture':label,'engineering_status':'PASS','pairs':5,
            'core_rmsd_a':[a['core_rmsd_a'] for a in audits],
            'ddg_runs':[p['ddg_kcal_mol'] for p in pairs],
            'scientific_status':'RAW_MODELS_ONLY_CURRENT_QUALIFICATION_NOT_INHERITED'})
    write_json(case/'T2_saved_model_replay.json',{'engineering_status':'PASS','new_foldx_calls':0,'fixtures':reports})
    # Exact numeric-input and binary identities are checked again by the worker.
    folder=case/'legacy_raw_models';folder.mkdir(exist_ok=True)
    groups={}
    for row in raw:
        pair=row['paired_models'][0]
        groups.setdefault((pair['binary_sha256'],pair['normalized_input_sha256']),[]).append(row)
    for i,((binary,normalized),rows) in enumerate(groups.items()):
        write_json(folder/str(i)/'paired_model_manifest.json',{'schema_version':'1.0',
            'binary_sha256':binary,'input_structure_sha256':normalized,'records':rows})
    return reports


def run(case,assets):
    preflight(case,assets)
    replay_saved_models(case,assets)
    base={'phase':'calibration','confirmation_count':1,'stage_limit':1,'pool_cap':32,'schedule':[5]}
    reports=[]
    first=launch(case,assets,'T3_calibration',base,1200);reports.append(first)
    if first['engineering_status']!='COMPLETED':return reports
    candidates=read(case/'results/panel_candidates.json')
    confirmed=[c['candidate_id'] for c in candidates if c.get('stability',{}).get('consensus_status')=='CONSENSUS_PASS']
    if not confirmed:
        write_json(case/'blocked.json',{'status':'T3_NO_STABLE_DOUBLE','reason':'No passing double in one-combination pilot; no scientific pass inferred'})
        return reports
    remaining=1200-first['charged_seconds']
    if remaining<=0:return reports
    production=launch(case,assets,'T3_production',{**base,'phase':'production','candidate_ids':confirmed},remaining);reports.append(production)
    if production['engineering_status']!='COMPLETED':return reports
    cost=max(1,sum(r['charged_seconds'] for r in reports));byte_cost=max(1,sum(r['new_bytes'] for r in reports))
    catalog=read(case/'results/candidate_pool_v2.json')['candidates']
    ids=[c['candidate_id'] for c in catalog[:4]]
    pilot_remaining=3600-read_budget(case/'case_budget.json')['elapsed_seconds']
    if pilot_remaining<=0:return reports
    for phase in ('calibration','production'):
        report=launch(case,assets,'T4_'+phase,{**base,'phase':phase,'candidate_ids':ids,
            'stage_limit':None,'confirmation_count':4},pilot_remaining)
        reports.append(report)
        pilot_remaining=3600-read_budget(case/'case_budget.json')['elapsed_seconds']
        if report['engineering_status']!='COMPLETED' or pilot_remaining<=0:return reports
    reports.extend(triple_test(case,assets,base,ids,pilot_remaining))
    tested=set(ids)
    while True:
        remaining=case_remaining(case)
        count,cap=affordable_batch(remaining,case_storage_limit(case)-owned_bytes(case),cost,byte_cost)
        pending=[c['candidate_id'] for c in catalog if c['candidate_id'] not in tested]
        if count==0 or not pending:break
        next_ids=rank_pending(case,catalog,pending)[:count];all_ids=sorted(tested | set(next_ids))
        started=time.monotonic();before=owned_bytes(case)
        for phase in ('calibration','production'):
            report=launch(case,assets,f'T6_{len(tested)}_'+phase,{**base,'phase':phase,
                'candidate_ids':all_ids,'stage_limit':None,'confirmation_count':len(all_ids),
                'pool_cap':cap},remaining-1800)
            reports.append(report)
            remaining=case_remaining(case)
            if report['engineering_status']!='COMPLETED':return reports
        cost=max(cost,(time.monotonic()-started)/len(next_ids))
        byte_cost=max(byte_cost,max(1,owned_bytes(case)-before)/len(next_ids))
        tested.update(next_ids)
        catalog=read(case/'results/candidate_pool_v2.json')['candidates']
        if reports[-1].get('physical',{}).get('goal_status')=='ACHIEVED':break
    return reports


def rank_pending(case,catalog,ids):
    import numpy as np
    from scipy.sparse import load_npz
    from ecd_epitope_scanner.surface_perturbation.coverage import largest_blind
    folder=case/'results/surface_perturbation'
    with np.load(folder/'canonical_reference.npz',allow_pickle=False) as data:
        area,owner=data['area'],data['owner']
    covered=np.zeros(len(area),bool)
    matrix_path=folder/'coverage_matrix.npz'
    if matrix_path.is_file():
        with np.load(matrix_path,allow_pickle=False) as data:
            if len(data['hard']):covered=np.any(data['hard'],axis=0)
    graph=load_npz(folder/'reference_graph.npz')
    selected={c['candidate_id']:c for c in catalog if c['candidate_id'] in ids}
    scores={}
    for cid,c in selected.items():
        hit=np.isin(owner,c['covered_positions_by_scenario'].get('geodesic_12',[]))
        scores[cid]=(-float(area[hit&~covered].sum()),largest_blind(graph,covered|hit,area),c['structural_risk'],cid)
    ordered=sorted(ids,key=lambda cid:scores[cid])
    write_json(case/'next_batch_priority.json',{'candidate_ids':ordered,'scores':scores,
        'interpretation':'Geometric proxy scheduling against real official hard blind mask; never a qualification or coverage measurement'})
    return ordered


def choose_menu_supplement(case,existing):
    """Expose one eligible backup to the measured blind surface, with provenance."""
    import numpy as np
    from ecd_epitope_scanner.joint_design_v4.coverage_profile import vector
    from ecd_epitope_scanner.joint_design_v4.candidate_search import _distances
    rows=read(case/'results/mutation_menu_v2.json')['rows']
    with np.load(case/'results/surface_perturbation/canonical_reference.npz',allow_pickle=False) as data:
        area,owner=data['area'],data['owner']
    covered=np.zeros(len(area),bool)
    matrix=case/'results/surface_perturbation/coverage_matrix.npz'
    if matrix.is_file():
        with np.load(matrix,allow_pickle=False) as data:
            if len(data['hard']):covered=np.any(data['hard'],axis=0)
    if covered.all():return None
    graph={}
    for e in read(case/'results/surface_graph_v4.json')['edges']:
        graph.setdefault(e['source'],{})[e['target']]=e['distance_a']
        graph.setdefault(e['target'],{})[e['source']]=e['distance_a']
    selected={}
    for row in rows:
        if row['menu_selected']:selected.setdefault(row['residue']['uniprot_pos'],[]).append(row)
    previous={(r['position'],r['mutant']) for r in existing}
    expanded_sites={r['position'] for r in existing if r['reason']=='MEASURED_HARD_BLIND_ZONE_BACKUP_MECHANISM'}
    supports={}
    choices=[]
    for row in rows:
        p=row['residue']['uniprot_pos'];aa=row['mutant_aa']
        if (not row['eligible'] or row.get('review_required') or row['menu_selected']
                or p not in graph or (p,aa) in previous):continue
        if p not in supports:
            positions=[q for q,d in _distances(graph,p,set(graph)).items() if d<=12]
            supports[p]=float(area[np.isin(owner,positions)&~covered].sum())
        if supports[p]<=0:continue
        v=vector(row['chemical_descriptors'])
        diversity=min((math.dist(v,vector(s['chemical_descriptors'])) for s in selected.get(p,[])),default=0.)
        choices.append(((p in expanded_sites,-supports[p],-diversity,row['spurs_ddg_kcal_mol'],p,aa),row))
    if not choices:return None
    score,row=min(choices,key=lambda item:item[0]);p=row['residue']['uniprot_pos']
    return {'position':p,'mutant':row['mutant_aa'],'reason':'MEASURED_HARD_BLIND_ZONE_BACKUP_MECHANISM',
        'scheduling_proxy_blind_area_a2':-score[1],'chemical_distance_to_current_menu':-score[2],
        'official_hard_coverage_at_request':float(area[covered].sum()/area.sum())}


def completed_triple_roles(case):
    """Exclude only actually completed production attempts, with their role."""
    import hashlib
    completed=set()
    history=read(case/'screened_candidate_history.json') if (case/'screened_candidate_history.json').is_file() else {}
    for path in (case/'stages').glob('*_report.json'):
        report=read(path);config=report.get('configuration',{})
        if report.get('engineering_status')!='COMPLETED' or config.get('phase')!='production':continue
        role=config.get('triple_role')
        if not role:continue
        cid=config.get('triple_candidate_id')
        if not cid:
            parents=[history.get(p,{}).get('candidate',{}) for p in config.get('triple_parent_ids',[])]
            notations={m['notation'] for parent in parents for m in parent.get('mutations',[])}
            if len(notations)==3:
                cid='CV2_'+hashlib.sha256('|'.join(sorted(notations)).encode()).hexdigest()[:16]
        if cid:completed.add((cid,role))
    return completed


def triple_test(case,assets,base,parent_ids,seconds,*,remote=False,proposal_only=False):
    from itertools import combinations
    from types import SimpleNamespace as NS
    import hashlib
    from ecd_epitope_scanner.surface_perturbation.compression import merge_parents
    from ecd_epitope_scanner.joint_design_v4.candidate_search import _distances
    snapshot=read(case/'results/surface_graph_v4.json');graph={}
    for edge in snapshot['edges']:
        graph.setdefault(edge['source'],{})[edge['target']]=edge['distance_a']
        graph.setdefault(edge['target'],{})[edge['source']]=edge['distance_a']
    candidates=read(case/'results/panel_candidates.json')
    pairs=read(case/'results/surface_perturbation/paired_manifest.json')
    eligible=[]
    for c in candidates:
        if c['candidate_id'] not in parent_ids or c['order']!=2 or len(pairs.get(c['candidate_id'],[]))<5:
            continue
        if not c.get('stability') or c['stability']['consensus_status']!='CONSENSUS_PASS':continue
        eligible.append(NS(**{**c,'mutations':[NS(**{**m,'residue':NS(**m['residue'])}) for m in c['mutations']]}))
    from ecd_epitope_scanner.joint_design_v4.coverage_profile import sequence_reasons
    sequence=read(case/'results/uniprot.json')['sequence']['value']
    summary=read(case/'results/surface_perturbation_preview.json')
    site_evidence=summary.get('candidate_site_evidence',{})
    role='TRIPLE_COMPRESSION' if summary.get('goal_status')=='ACHIEVED' else 'TRIPLE_EXPANSION'
    if role=='TRIPLE_COMPRESSION':
        eligible=[c for c in eligible if c.candidate_id in summary.get('selected_candidate_ids',[])]
    completed_roles=completed_triple_roles(case)
    if proposal_only and (case/'v3_triple_attempts.json').exists():
        completed_roles.update((row['candidate_id'],row['role']) for row in read(case/'v3_triple_attempts.json'))
    proposals=[];geometry_rejected=[];sequence_rejected=[]
    for left,right in combinations(eligible,2):
        merged=merge_parents(left,right)
        if merged is not None:
            reasons=sequence_reasons(sequence,[(m.residue.uniprot_pos,m.residue.one_letter,m.mutant_aa) for m in merged])
            if reasons:
                sequence_rejected.append({'parent_ids':[left.candidate_id,right.candidate_id],'reasons':reasons})
                continue
            positions=[m.residue.uniprot_pos for m in merged]
            distances={p:_distances(graph,p,set(graph)) for p in positions}
            span=max(distances[p].get(q,float('inf')) for p,q in combinations(positions,2))
            signature='|'.join(sorted(m.notation for m in merged))
            cid='CV2_'+hashlib.sha256(signature.encode()).hexdigest()[:16]
            if (cid,role) in completed_roles:continue
            if remote and span<=24:continue
            if not remote and span>24:
                geometry_rejected.append({'parent_ids':[left.candidate_id,right.candidate_id],'geodesic_span_a':span})
                continue
            supported=sum(bool(site_evidence.get(c.candidate_id,{}).get('complete')) for c in (left,right))
            proposals.append(((-supported,span,left.candidate_id,right.candidate_id),(left,right,merged)))
    if not proposals:
        write_json(case/'T5_decision.json',{'status':'SKIPPED_NO_LEGAL_SHARED_REPLACEMENT_PARENTS',
            'measured_parent_ids':[c.candidate_id for c in eligible],
            'geometry_rejected':geometry_rejected,
            'sequence_rejected':sequence_rejected,
            'interpretation':'No synthetic real triple or control evidence is inferred from absent parents'})
        return []
    left,right,mutations=min(proposals,key=lambda item:item[0])[1]
    signature='|'.join(sorted(m.notation for m in mutations))
    triple_id='CV2_'+hashlib.sha256(signature.encode()).hexdigest()[:16]
    catalog=read(case/'results/candidate_pool_v2.json')['candidates']
    controls=[];control_signatures=[]
    for options in combinations(mutations,2):
        key=tuple(sorted(o.notation for o in options))
        control_signatures.append(list(key))
        # Force the exact genotype into the coverage double pool before its
        # bounded truncation, so a cap change cannot switch CV2 to V41 identity.
        controls.append('CV2_'+hashlib.sha256('|'.join(key).encode()).hexdigest()[:16])
    retained=bounded_retained_ids(pairs,summary,[],[left.candidate_id,right.candidate_id,*controls])
    active_ids=sorted(retained|set(controls)|{triple_id})
    if len(active_ids)>256:raise ValueError('Triple confirmation exceeds active pool bound')
    if remote:
        remote_path=case/'remote_triple_confirmation_ids.json'
        admitted=read(remote_path) if remote_path.is_file() else []
        if triple_id not in admitted:
            if len(admitted)>=2:
                write_json(case/'remote_triple_stop.json',{'status':'ROUTE_REMOTE_CONFIRMATION_LIMIT_REACHED','candidate_id':triple_id})
                return []
            write_json(remote_path,[*admitted,triple_id])
    config={**base,'stage_limit':None,'candidate_ids':active_ids,
        'confirmation_count':len(active_ids),'pool_cap':max(base.get('pool_cap',32),len(active_ids)),
        'allow_remote_triples':remote,
        'triple_parent_ids':[left.candidate_id,right.candidate_id],'triple_role':role,
        'triple_candidate_id':triple_id}
    config['mandatory_candidate_signatures']=[*base.get('mandatory_candidate_signatures',[]),*control_signatures]
    if role=='TRIPLE_COMPRESSION':config['compression_baseline_ids']=summary['selected_candidate_ids']
    write_json(case/'T5_decision.json',{'status':'LEGAL_PARENT_PROPOSAL','role':role,
        'geometry_rejected':geometry_rejected,'sequence_rejected':sequence_rejected,
        'exact_double_control_ids':controls,'configuration':config})
    if proposal_only:return config
    reports=[]
    start=time.monotonic()
    stage_prefix=('T5' if read_budget(case/'case_budget.json')['elapsed_seconds']<3600 else 'T6_triple')+('_remote' if remote else '')
    attempt=stage_attempt_number(case,stage_prefix)
    for phase in ('calibration','production'):
        remaining=seconds-(time.monotonic()-start)
        if remaining<=0:break
        stage=f'{stage_prefix}_{attempt}_{phase}'
        # Preserve attempt provenance if a later proposal revisits this stage.
        if (case/'stages'/f'{stage}_report.json').exists():stage+='_'+triple_id
        reports.append(launch(case,assets,stage,{**config,'phase':phase},remaining))
        if reports[-1]['engineering_status']!='COMPLETED':break
        preserve_pair_index(case)
        if phase=='production':record_completed_snapshot(case,reports[-1])
    return reports


def review_final_triples(case,assets,base,parent_ids,seconds):
    """Review novel local proposals, then at most two remote proposals.

    Stop after one accepted compression so separate phases cannot reset the
    cumulative one-percentage-point loss allowance against a new baseline.
    """
    reports=[];start=time.monotonic()
    for remote in (False,True):
        while True:
            remaining=seconds-(time.monotonic()-start)
            if remaining<=0:return reports
            summary=read(case/'results/surface_perturbation_preview.json')
            if summary.get('goal_status')=='ACHIEVED' and len(summary.get('selected_candidate_ids',[]))<=4:
                write_json(case/'triple_route_review.json',{'status':'NO_CONSTRUCT_REDUCTION_ABOVE_MINIMUM',
                    'formal_constructs':len(summary.get('selected_candidate_ids',[]))})
                return reports
            if any(a.get('accepted') for a in summary.get('compression_audits',[])):return reports
            batch=triple_test(case,assets,base,parent_ids,remaining,remote=remote)
            if not batch:break
            reports.extend(batch)
            if batch[-1]['engineering_status']!='COMPLETED':return reports
    write_json(case/'triple_route_review.json',{'status':'RULE_GENERATED_PROPOSALS_REVIEWED',
        'remote_confirmation_limit':2,'interpretation':'Not exhaustive enumeration of all triples'})
    return reports


def continue_after_pilot(case,assets):
    """Continue a verified five-pair pilot; never infer completion from exit alone."""
    base=read(case/'stages/T3_saved_control_production_request.json')
    base={k:v for k,v in base.items() if k not in {'stage','seconds','phase'}}
    output=case/'results'
    pairs=read(output/'surface_perturbation/paired_manifest.json')
    if not pairs or any(len(rows)<5 for rows in pairs.values()):
        raise ValueError('T3 five-pair pilot has not completed')
    from ecd_epitope_scanner.surface_perturbation.paired import validate_pair
    for rows in pairs.values():
        for pair in rows[:5]:validate_pair(pair)
    reports=[read(case/'stages'/name) for name in ('T3_saved_control_calibration_report.json','T3_saved_control_production_report.json')]
    if any(r['engineering_status']!='COMPLETED' for r in reports):raise ValueError('Pilot engineering gate failed')
    catalog=read(output/'candidate_pool_v2.json')['candidates']
    ids=list(pairs)
    boundary=[c['candidate_id'] for c in catalog if c['candidate_id'] not in ids
        and any('RSASA_BOUNDARY_020_025' in m.get('warnings',[]) for m in c['mutations'])]
    if boundary:ids.append(boundary[0])
    ids += [c['candidate_id'] for c in catalog if c['candidate_id'] not in ids][:4-len(ids)]
    remaining=3600-read_budget(case/'case_budget.json')['elapsed_seconds']
    for phase in ('calibration','production'):
        report=launch(case,assets,'T4_'+phase,{**base,'phase':phase,'candidate_ids':ids,
            'stage_limit':None,'confirmation_count':4},remaining)
        reports.append(report)
        remaining=3600-read_budget(case/'case_budget.json')['elapsed_seconds']
        if report['engineering_status']!='COMPLETED' or remaining<=0:return reports
    reports.extend(triple_test(case,assets,base,ids,remaining))
    cost=max(1,sum(r['charged_seconds'] for r in reports[:2]))
    byte_cost=max(1,sum(r['new_bytes'] for r in reports[:2]))
    tested=set(ids)
    while True:
        remaining=case_remaining(case)
        count,cap=affordable_batch(remaining,case_storage_limit(case)-owned_bytes(case),cost,byte_cost)
        pending=[c['candidate_id'] for c in catalog if c['candidate_id'] not in tested and c['order']==2]
        if count==0 or not pending:break
        additions=rank_pending(case,catalog,pending)[:count]
        all_ids=sorted(tested|set(additions))
        start=time.monotonic();before=owned_bytes(case)
        for phase in ('calibration','production'):
            report=launch(case,assets,f'T6_{len(tested)}_'+phase,{**base,'phase':phase,
                'candidate_ids':all_ids,'stage_limit':None,'confirmation_count':len(all_ids),'pool_cap':max(cap,len(all_ids))},remaining-1800)
            reports.append(report)
            remaining=case_remaining(case)
            if report['engineering_status']!='COMPLETED':return reports
        cost=max(cost,(time.monotonic()-start)/len(additions))
        byte_cost=max(byte_cost,max(1,owned_bytes(case)-before)/len(additions))
        tested.update(additions);catalog=read(output/'candidate_pool_v2.json')['candidates']
        if reports[-1].get('physical',{}).get('goal_status')=='ACHIEVED':
            reports.extend(triple_test(case,assets,base,list(tested),max(0,remaining-1800)))
            break
    return reports


def audit_attempt_resolution(reports,retired):
    """Keep failed attempts charged; certify only explicitly recovered nodes."""
    resolved=[];outstanding=[]
    for index,report in enumerate(reports):
        if report['engineering_status']=='COMPLETED':continue
        config=report.get('configuration',{})
        required=set(config.get('candidate_ids',[]))-set(retired)
        recovery=next((later for later in reports[index+1:]
            if later['engineering_status']=='COMPLETED'
            and later.get('configuration',{}).get('phase')==config.get('phase')
            and (bool(config.get('candidate_ids')) or all(later.get('configuration',{}).get(key)==config.get(key)
                for key in ('stage_limit','confirmation_count','pool_cap','schedule')))
            and required<=set(later.get('configuration',{}).get('candidate_ids',[]))),None)
        row={'failed_stage':report['stage'],'required_candidate_ids':sorted(required),
            'retired_candidate_ids':sorted(set(config.get('candidate_ids',[]))&set(retired))}
        if recovery is None:outstanding.append(row)
        else:resolved.append({**row,'recovered_by':recovery['stage']})
    return {'status':'PASS' if not outstanding else 'INCOMPLETE',
        'resolved_failures':resolved,'outstanding_failures':outstanding,
        'interpretation':'Failure charges and original reports remain; a later successful phase must cover the same non-retired nodes'}


def write_case_provenance(case,output,guard):
    """Portable current source and stage evidence; assets remain pinned references."""
    import zipfile
    source=Path(__file__).resolve().parents[1]
    overlay=read(source/'source_overlay_manifest.json')
    files={}
    for name,digest in overlay['files'].items():
        path=source/name
        if sha256_file(path)!=digest:raise ValueError('Final source snapshot changed: '+name)
        files['source/ecd_epitope_scanner/'+name]=path
    files['source/ecd_epitope_scanner/source_overlay_manifest.json']=source/'source_overlay_manifest.json'
    replay={}
    from ecd_epitope_scanner.surface_perturbation.paired import validate_pair
    groups=[]
    for manifest_path in (case/'legacy_raw_models').glob('*/paired_model_manifest.json'):
        for record_index,record in enumerate(read(manifest_path)['records']):
            groups.append((manifest_path.parent.name+'_'+str(record_index),record['paired_models'],'RAW_SAVED_MODEL_REPLAY_ONLY'))
    physical=read(output/'surface_perturbation/paired_manifest.json')
    seen={(r['batch_path'],r['group'],r['run']) for rows in physical.values() for r in rows}
    screened=[]
    index_path=case/'evidence_index/surface_perturbation/paired_manifest.json'
    if index_path.is_file():
        files['case/evidence_index/surface_perturbation/paired_manifest.json']=index_path
        historical=[]
        for rows in read(index_path).values():
            for row in rows:
                identity=(row['batch_path'],row['group'],row['run'])
                if identity not in seen:historical.append(row);seen.add(identity)
        if historical:groups.append(('historical_physical_pairs',historical,'RAW_HISTORICAL_PAIR_CURRENT_QUALIFICATION_NOT_INHERITED'))
    for manifest_path in output.rglob('paired_model_manifest.json'):
        for record in read(manifest_path).get('records',[]):
            for row in record.get('paired_models',[]):
                key=(row['batch_path'],row['group'],row['run'])
                if key not in seen:screened.append(row);seen.add(key)
    if screened:groups.append(('screened_unqualified',screened,'RAW_SCREENED_MODEL_NOT_QUALIFIED'))
    for identity,rows,qualification in groups:
        replay[identity]=[]
        for index,row in enumerate(rows):
            validate_pair(row);item=dict(row);prefix=f'saved_model_replay/{identity}/{index:03d}'
            for key in ('wt_path','mut_path','pdb_list_path','ddg_file_path'):
                if row.get(key):
                    path=Path(row[key]);name=prefix+'/'+path.name;files[name]=path;item[key]=name
            for filename,key in (('input.pdb','input_path'),('individual_list.txt','mutant_file_path')):
                path=Path(row['batch_path'])/filename
                if not path.is_file():raise ValueError('Saved replay dependency missing: '+str(path))
                files[prefix+'/'+filename]=path;item[key]=prefix+'/'+filename
            item['batch_path']=prefix;item['qualification']=qualification
            replay[identity].append(item)
    if (case/'T2_native_author.pdb').is_file():files['case/T2_native_author.pdb']=case/'T2_native_author.pdb'
    for path in case.glob('*.json'):
        prefix='case/historical_snapshots/' if path.name in {'validation_report.json','resume_checkpoint.json'} else 'case/'
        files[prefix+path.name]=path
    for path in case.glob('T1*.log'):files['case/'+path.name]=path
    # Preserve the exact previous frozen source before a boundary upgrade,
    # alongside its independently checked per-file manifest.
    for path in case.glob('source_snapshot_*.zip'):
        files['case/source_snapshots/'+path.name]=path
    for directory in ('stages','ledgers'):
        for path in (case/directory).glob('*'):
            if path.is_file() and path.suffix in {'.json','.log'}:
                files['case/'+directory+'/'+path.name]=path
    guard.check(reserved=sum(p.stat().st_size for p in files.values()),solving=False)
    destination=output/'validation_provenance.zip';pending=destination.with_suffix('.zip.pending')
    try:
        with zipfile.ZipFile(pending,'w',zipfile.ZIP_DEFLATED) as archive:
            identities={}
            for name,path in sorted(files.items()):
                content=path.read_bytes();archive.writestr(name,content)
                import hashlib
                identities[name]=hashlib.sha256(content).hexdigest()
            content=json.dumps(replay,sort_keys=True).encode()
            archive.writestr('saved_model_replay_pairs.json',content)
            identities['saved_model_replay_pairs.json']=hashlib.sha256(content).hexdigest()
            archive.writestr('provenance_manifest.json',json.dumps({'files':identities,
                'scope':'Actual current frozen source, tests and stage history; original tool/weight identities are referenced, not bundled'},sort_keys=True))
        pending.replace(destination)
    finally:pending.unlink(missing_ok=True)
    return destination


def finalize(case,reports):
    """One final portable export; streaming independent checks avoid extraction copies."""
    import hashlib
    import zipfile
    from ecd_epitope_scanner.surface_perturbation.export import evidence_archive
    from ecd_epitope_scanner.exporting import finalize_scientific_bundle
    from ecd_epitope_scanner.surface_perturbation.archive_audit import validate_archive
    from ecd_epitope_scanner.scripts.audit_run_bundle import audit_run_bundle
    from ecd_epitope_scanner.scripts.audit_ctla4_case import audit as audit_case
    output=case/'results';physical=output/'surface_perturbation'
    budget=RunBudget(case/'case_budget.json',case_limit(case))
    try:
        if budget.remaining()<1:
            return
        guard=StorageGuard(case)
        retired=read(case/'retired_candidate_ids.json') if (case/'retired_candidate_ids.json').is_file() else {}
        attempts=audit_attempt_resolution(reports,retired)
        report={'stages':reports,'attempt_resolution':attempts,
            'engineering_status':attempts['status'],
            'scientific_status':read(output/'surface_perturbation_preview.json').get('goal_status','NOT_EVALUATED')
                if (output/'surface_perturbation_preview.json').is_file() else 'NOT_EVALUATED',
            'interpretation':'Only official goal_status can establish the 80 percent target; shadow stays separate'}
        pairs_path=physical/'paired_manifest.json'
        if not pairs_path.is_file():return
        semantic=audit_case(output)
        write_json(output/'independent_case_semantic_audit.json',semantic)
        checked_objects=0
        for path in (case/'objects').glob('*/*'):
            if path.is_file():
                if sha256_file(path)!=path.name:raise ValueError('Case evidence object changed: '+str(path))
                checked_objects+=1
        report['verified_content_objects']=checked_objects
        write_case_provenance(case,output,guard)
        # Reserve the physical archive and its enclosing final bundle peak.
        size=sum(p.stat().st_size for p in physical.rglob('*') if p.is_file())
        for rows in read(pairs_path).values():
            for pair in rows:
                size+=sum(Path(pair[k]).stat().st_size for k in ('wt_path','mut_path','pdb_list_path'))
        guard.check(reserved=2*size,solving=False)
        # Preserve obsolete archives from the failed early finalization rather
        # than silently overwriting their historical evidence.
        obsolete=case/'historical_exports';obsolete.mkdir(exist_ok=True)
        for name in ('surface_perturbation_evidence.zip','ecd_epitope_design.zip'):
            old=output/name
            if old.is_file():old.replace(obsolete/(sha256_file(old)+'_'+name))
        archive=evidence_archive(output,read(pairs_path),physical/'calibration')
        with zipfile.ZipFile(archive) as bundle:
            unpacked=sum(item.file_size for item in bundle.infolist())
        guard.check(reserved=unpacked,solving=False)
        audit_root=case/'archive_audit';audit_root.mkdir(exist_ok=True)
        portable=validate_archive(archive,temporary_root=audit_root)
        provenance=output/'validation_provenance.zip'
        with zipfile.ZipFile(provenance) as bundle:
            unpacked=sum(item.file_size for item in bundle.infolist())
        guard.check(reserved=unpacked,solving=False)
        provenance_audit=validate_archive(provenance,temporary_root=audit_root,
            manifest_name='provenance_manifest.json',pairs_name='saved_model_replay_pairs.json')
        with zipfile.ZipFile(archive) as bundle:
            manifest=json.loads(bundle.read('evidence_manifest.json'))['files']
            checked=0
            for name,digest in manifest.items():
                calculated=hashlib.sha256()
                with bundle.open(name) as stream:
                    while data:=stream.read(1024**2):calculated.update(data)
                if calculated.hexdigest()!=digest:raise ValueError('Independent archive hash mismatch: '+name)
                checked+=1
        write_json(output/'independent_bundle_audit.json',{'engineering_status':'PASS','files_checked':checked,
            'portable_model_ddg_null_audit':portable,
            'portable_source_stage_saved_model_audit':provenance_audit,
            'evidence_archive_sha256':sha256_file(archive),'test':'Every archived object streamed and independently hashed'})
        if (output/'run_manifest.json').is_file():
            manifest=read(output/'run_manifest.json')
            manifest['files']['surface_perturbation_evidence.zip']=sha256_file(archive)
            manifest['files']['independent_bundle_audit.json']=sha256_file(output/'independent_bundle_audit.json')
            for name in ('candidate_sequence_replay.csv','independent_case_semantic_audit.json','validation_provenance.zip'):
                manifest['files'][name]=sha256_file(output/name)
            write_json(output/'run_manifest.json',manifest)
            finalize_scientific_bundle(output)
            report['final_bundle_audit']=audit_run_bundle(output)
            report['final_bundle_bytes']=(output/'ecd_epitope_design.zip').stat().st_size
            report['final_bundle_sha256']=sha256_file(output/'ecd_epitope_design.zip')
        guard.check(solving=False)
    except Exception as exc:
        if 'report' in locals():
            report['engineering_status']='INCOMPLETE'
            report['delivery_failure']={'type':type(exc).__name__,'message':str(exc)}
        raise
    finally:
        budget.save(False)
        if 'report' in locals():
            report['case_elapsed_seconds']=read_budget(case/'case_budget.json')['elapsed_seconds']
            report['owned_bytes']=owned_bytes(case)
            write_json(case/'validation_report.json',report)


def completed_candidate_ids(reports):
    identities=set()
    for report in reports:
        config=report.get('configuration',{})
        if report.get('engineering_status')=='COMPLETED' and config.get('phase')=='production':
            identities.update(config.get('candidate_ids') or report.get('physical',{}).get('candidate_site_evidence',{}))
    return identities


def preserve_pair_index(case):
    """Cumulative raw pair identities survive a bounded active scoring pool."""
    destination=case/'evidence_index/surface_perturbation/paired_manifest.json'
    index=read(destination) if destination.is_file() else {}
    current=read(case/'results/surface_perturbation/paired_manifest.json')
    for cid,rows in current.items():
        merged={(row['batch_path'],row['group'],row['run']):row for row in index.get(cid,[])}
        for row in rows:
            identity=(row['batch_path'],row['group'],row['run'])
            if identity in merged and merged[identity]!=row:
                raise ValueError('Historical paired identity changed')
            merged[identity]=row
        index[cid]=list(merged.values())
    write_json(destination,index)
    return index


def bounded_retained_ids(pairs,summary,improvable,parent_ids=(),maximum=248):
    ids=set(pairs)
    if len(ids)<=maximum:return ids
    required=set(summary.get('independent_K_verified_candidate_ids',[]))
    required.update(summary.get('shadow_positive_lower_bound',{}).get('candidate_ids',[]))
    required.update(cid for cid,row in summary.get('candidate_site_evidence',{}).items() if row.get('complete'))
    required.update(improvable);required.update(parent_ids)
    required &= ids
    if len(required)>maximum:raise ValueError('Necessary active evidence exceeds bounded pool capacity')
    return required|set(sorted(ids-required)[:maximum-len(required)])


def record_completed_snapshot(case,report):
    if report.get('engineering_status')!='COMPLETED' or report.get('configuration',{}).get('phase')!='production':
        raise ValueError('Only completed production can record processed search identities')
    path=case/'screened_candidate_history.json'
    history=read(path) if path.is_file() else {}
    origins={row['candidate_id']:row for row in report.get('configuration',{}).get('conditional_pair_backups',[])}
    for candidate in read(case/'results/panel_candidates.json'):
        cid=candidate['candidate_id'];old=history.get(cid,{})
        history[cid]={**old,'last_completed_stage':report['stage'],'candidate':candidate}
        if cid in origins:history[cid]['conditional_backup_origin']=origins[cid]
        for key in ('candidate_site_evidence','candidate_local_calibration'):
            if cid in report.get('physical',{}).get(key,{}):history[cid][key]=report['physical'][key][cid]
    write_json(path,history)
    preserve_pair_index(case)
    return history


def choose_failure_backup(case,completed,*,scientific_v3=False):
    """One stability-first menu alternative after an actual prescreen failure."""
    from ecd_epitope_scanner.joint_design_v4.models import ResidueKey,V4MutationOption
    from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
    from ecd_epitope_scanner.models import to_dict
    output=case/'results'
    history=read(case/'screened_candidate_history.json')
    pool=read(output/'panel_candidates.json')
    if scientific_v3:pool=list({c['candidate_id']:c for c in [*pool,*[row['candidate'] for row in history.values()]]}.values())
    failed=[];failure_reasons={}
    for candidate in pool:
        cid=candidate['candidate_id'];row=history.get(cid,{})
        if candidate['order']!=2:continue
        status=(candidate.get('stability') or {}).get('consensus_status')
        if status in {'COMBINATION_REJECT','MODEL_DISAGREEMENT','STRUCTURAL_ENERGY_CONFLICT'}:
            failed.append(candidate);failure_reasons[cid]='OBSERVED_COMBINATION_STABILITY_FAILURE'
        elif scientific_v3 and row.get('candidate_site_evidence',{}).get('complete') is False:
            failed.append(candidate);failure_reasons[cid]='MEASURED_SITE_EVIDENCE_MISSING_OR_UNRESOLVED_ALTERNATIVE'
    if not failed:return None
    visits={}
    for row in history.values():
        candidate=row['candidate']
        if candidate['order']==2:
            pair=tuple(sorted(m['residue']['uniprot_pos'] for m in candidate['mutations']))
            visits[pair]=visits.get(pair,0)+1
    saved=read(output/'surface_graph_v4.json');graph={node['uniprot_pos']:{} for node in saved['nodes']}
    for edge in saved['edges']:
        a,b,d=edge['source'],edge['target'],edge['distance_a'];graph[a][b]=d;graph[b][a]=d
    rows=read(output/'mutation_menu_v2.json')['rows'];choices=[]
    for parent in failed:
        pair=tuple(sorted(m['residue']['uniprot_pos'] for m in parent['mutations']))
        # Already tried alternatives remain evidence, not an obligation to
        # enumerate all nine genotypes of this position pair.
        if visits.get(pair,0)>1:continue
        options=[V4MutationOption(ResidueKey(**row['residue']),row['mutant_aa'],
            row['residue']['one_letter']+str(row['residue']['uniprot_pos'])+row['mutant_aa'],
            row['chemistry_class'],row['spurs_ddg_kcal_mol'],0,True)
            for row in rows if row['residue']['uniprot_pos'] in pair and row['eligible']
                and not row.get('review_required') and row['menu_selected']]
        pool=generate_candidates([],options,graph,{'coverage_scan_v2':True,
            'canonical_sequence':read(output/'uniprot.json')['sequence']['value'],
            'global_shortlist':16,'exclude_candidate_ids':sorted(completed)})
        if not pool:continue
        winner=min(pool,key=lambda c:(c.structural_risk,sum(o.spurs_ddg_kcal_mol for o in c.mutations),c.candidate_id))
        candidate=to_dict(winner)
        candidate['conditional_backup_parent_id']=parent['candidate_id']
        candidate['conditional_backup_reason']=failure_reasons[parent['candidate_id']] if scientific_v3 else 'OBSERVED_COMBINATION_STABILITY_PRESCREEN_FAILURE'
        choices.append(candidate)
    if not choices:return None
    ranked=rank_pending(case,choices,[c['candidate_id'] for c in choices])
    return next(c for c in choices if c['candidate_id']==ranked[0])


def replenish_saved_frontier(case,excluded,cap):
    """Rebuild lightweight unmeasured doubles from the exact saved menu/graph."""
    from ecd_epitope_scanner.joint_design_v4.models import ResidueKey,V4MutationOption
    from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
    from ecd_epitope_scanner.models import to_dict
    output=case/'results'
    rows=read(output/'mutation_menu_v2.json')['rows']
    amendment=case/'resource_contract_amendment.json'
    scope=read(amendment).get('exhaustion_scope') if amendment.is_file() else None
    full=scope=='ALL_ELIGIBLE_LOCAL_DOUBLES'
    seen_pairs=[]
    if scope=='ELIGIBLE_LOCAL_SITE_PAIRS':
        history=read(case/'screened_candidate_history.json') if (case/'screened_candidate_history.json').is_file() else {}
        seen_pairs=[sorted(m['residue']['uniprot_pos'] for m in row['candidate']['mutations'])
            for cid,row in history.items() if cid in excluded and row['candidate']['order']==2]
    options=[V4MutationOption(ResidueKey(**r['residue']),r['mutant_aa'],
        r['residue']['one_letter']+str(r['residue']['uniprot_pos'])+r['mutant_aa'],
        r['chemistry_class'],r['spurs_ddg_kcal_mol'],0,True) for r in rows
        if r['eligible'] and not r.get('review_required') and (full or r['menu_selected'])]
    saved_graph=read(output/'surface_graph_v4.json')
    graph={node['uniprot_pos']:{} for node in saved_graph.get('nodes',[])}
    for edge in saved_graph['edges']:
        a,b,d=edge['source'],edge['target'],edge['distance_a']
        graph.setdefault(a,{})[b]=d;graph.setdefault(b,{})[a]=d
    diagnostic={}
    pool=generate_candidates([],options,graph,{'coverage_scan_v2':True,
        'canonical_sequence':read(output/'uniprot.json')['sequence']['value'],
        'global_shortlist':cap,'exclude_candidate_ids':sorted(excluded),
        'site_pair_priority':scope=='ELIGIBLE_LOCAL_SITE_PAIRS','exclude_site_pairs':seen_pairs},diagnostics=diagnostic)
    if scope=='ELIGIBLE_LOCAL_SITE_PAIRS':
        contract_path=case/'site_pair_search_contract.json'
        contract={'scope':scope,'graph_sha256':sha256_file(output/'surface_graph_v4.json'),
            'initial_menu_sha256':sha256_file(output/'mutation_menu_v2.json'),
            'maximum_geodesic_span_a':24,'site_pairs':diagnostic['all_eligible_site_pairs'],
            'interpretation':'One representative per safe eligible local position pair; alternatives are not exhaustively solved'}
        if contract_path.is_file():
            previous=read(contract_path)
            if previous['graph_sha256']!=contract['graph_sha256'] or previous['site_pairs']!=contract['site_pairs']:
                raise ValueError('Frozen site-pair search scope changed')
        else:write_json(contract_path,contract)
    write_json(case/'frontier_replenishment.json',{'menu_sha256':sha256_file(output/'mutation_menu_v2.json'),
        'graph_sha256':sha256_file(output/'surface_graph_v4.json'),'diagnostics':diagnostic,
        'excluded_candidate_ids':sorted(excluded),'candidates':to_dict(pool),
        'interpretation':'Unmeasured lightweight search only; saved evidence is retained independently'})
    return to_dict(pool)


def resume_batches(case,assets,max_batches=None):
    """Resume the last incomplete batch with its exact identities and ledger.

    Completed pilots are not scheduled again. Historical failed attempts remain
    visible; each retry has a distinct request/log/report identity.
    """
    if case_limit(case) is None and read(case/'resource_contract_amendment.json').get('activation')!='ACTIVE':
        raise ValueError('Resource amendment must be activated at a closed stage boundary')
    stages=case/'stages'
    reports=[read(p) for p in sorted(stages.glob('*_report.json'),key=lambda p:p.stat().st_mtime)]
    production=[r for r in reports if r['engineering_status']=='COMPLETED'
        and r['configuration'].get('phase')=='production']
    if not production:raise ValueError('No completed production checkpoint to resume')
    previous=production[-1]
    completed=completed_candidate_ids(production)
    history=record_completed_snapshot(case,previous)
    retired_path=case/'retired_candidate_ids.json'
    retired=read(retired_path) if retired_path.is_file() else {}
    pending=next((r for r in reversed(reports) if r['stage'].startswith('T6_')
        and r['engineering_status']!='COMPLETED'
        and not (set(r['configuration'].get('candidate_ids',[]))-set(retired))<=completed),None)
    base=dict(previous['configuration'])
    if pending:
        base=dict(pending['configuration'])
    for key in ('stage','seconds','phase','candidate_ids','confirmation_count','failure'):
        base.pop(key,None)
    base['stage_limit']=None
    if case_limit(case) is None:
        base['full_replacement_search']=read(case/'resource_contract_amendment.json').get('exhaustion_scope')=='ALL_ELIGIBLE_LOCAL_DOUBLES'
    physical=read(case/'results/surface_perturbation_preview.json')
    unknown_audit={}
    for old in reversed(reports):
        unknown_audit=old.get('physical',{}).get('calibration_unknowns',{})
        if unknown_audit:break
    pair_counts={cid:len(rows) for cid,rows in read(case/'results/surface_perturbation/paired_manifest.json').items()}
    calibration_counts={r['candidate']:r.get('calibration_runs',0)
        for r in old.get('physical',{}).get('trajectories',[]) if 'candidate' in r}
    def extendable(cid,row):
        modes=row.get('observed_failure_modes',{})
        first_failures=modes.get('first_five_topology_rejected_runs',
            modes.get('topology_rejected_runs',0) if calibration_counts.get(cid)==5 else 0)
        return row.get('appendable') and pair_counts.get(cid,0)<20 and first_failures<2
    improvable=[cid for cid,row in unknown_audit.items() if extendable(cid,row)]
    if physical.get('goal_status')!='ACHIEVED' and improvable:
        base['schedule']=[5,10,20]
        # New five-model cohorts are reviewed in the following batch; they do
        # not automatically consume twenty repeats during initial calibration.
        base['calibration_extension_candidate_ids']=improvable[:2]
        write_json(case/'null_extension_decision.json',{
            'candidate_ids_with_improvable_statistics':improvable,
            'schedule':[5,10,20],'same_protocol_cohort_cap':20,
            'interpretation':'Engine appends only after rechecking finite statistical uncertainty; deterministic unknowns remain excluded'})
    cost=max(1.,sum(r['charged_seconds'] for r in reports
        if r['stage'].startswith('T3_saved_control')))
    byte_cost=max(1.,sum(r['new_bytes'] for r in reports
        if r['stage'].startswith('T3_saved_control')))
    batch_number=0
    if pending is None and physical.get('goal_status')=='ACHIEVED':
        return reports+review_final_triples(case,assets,base,list(completed),max(0,case_remaining(case)-1800))
    while max_batches is None or batch_number<max_batches:
        remaining=case_remaining(case)
        count,cap=affordable_batch(remaining,case_storage_limit(case)-owned_bytes(case),cost,byte_cost)
        if not count:break
        catalog=read(case/'results/candidate_pool_v2.json')['candidates']
        if pending:
            ids=[cid for cid in pending['configuration']['candidate_ids'] if cid not in retired]
            additions=set(ids)-completed
            cap=max(cap,pending['configuration']['pool_cap'],len(ids))
        else:
            site_priority=case_limit(case) is None and read(case/'resource_contract_amendment.json').get('exhaustion_scope')=='ELIGIBLE_LOCAL_SITE_PAIRS'
            pool=[] if site_priority else [c['candidate_id'] for c in catalog if c['order']==2 and c['candidate_id'] not in completed]
            if len(pool)<count:
                # All prior screened identities remain excluded, but only
                # measured physical evidence needs the next stage's whitelist.
                search_budget=RunBudget(case/'case_budget.json',case_limit(case))
                try:catalog=replenish_saved_frontier(case,completed|set(retired),cap)
                finally:search_budget.save(False)
                pool=[c['candidate_id'] for c in catalog]
            if not pool:
                if site_priority:
                    contract=read(case/'site_pair_search_contract.json')
                    seen={tuple(sorted(m['residue']['uniprot_pos'] for m in row['candidate']['mutations']))
                        for row in history.values() if row['candidate']['order']==2}
                    remaining_pairs=set(map(tuple,contract['site_pairs']))-seen
                    write_json(case/'site_pair_exhaustion_audit.json',{
                        'status':'EXHAUSTED' if not remaining_pairs else 'INCOMPLETE',
                        'total_site_pairs':len(contract['site_pairs']),
                        'processed_site_pairs':len(set(map(tuple,contract['site_pairs']))&seen),
                        'remaining_site_pairs':[list(pair) for pair in sorted(remaining_pairs)],
                        'contract_sha256':sha256_file(case/'site_pair_search_contract.json')})
                    if remaining_pairs:raise ValueError('Empty pool does not establish site-pair exhaustion')
                    reports.extend(review_final_triples(case,assets,base,list(completed),max(0,remaining-1800)))
                break
            backup=choose_failure_backup(case,completed|set(retired)) if site_priority and count>=2 else None
            additions=set(rank_pending(case,catalog,pool)[:count-(1 if backup else 0)])
            base['conditional_pair_backups']=[]
            if backup:
                catalog.append(backup);additions.add(backup['candidate_id'])
                base['conditional_pair_backups']=[{'candidate_id':backup['candidate_id'],
                    'parent_id':backup['conditional_backup_parent_id'],'reason':backup['conditional_backup_reason']}]
            mandatory={tuple(sorted(s)) for s in base.get('mandatory_candidate_signatures',[])}
            controls=[cid for cid,row in history.items()
                if tuple(sorted(m['notation'] for m in row['candidate']['mutations'])) in mandatory]
            retained=bounded_retained_ids(read(case/'results/surface_perturbation/paired_manifest.json'),
                physical,improvable,[*base.get('triple_parent_ids',[]),*controls])
            ids=sorted(retained|additions)
            cap=max(cap,len(ids))
            if cap>256:raise ValueError('Retained physical evidence exceeds hard candidate pool capacity')
        # Re-audit known context conflicts before expensive adapters. Retired
        # nodes remain explicit evidence, and are never counted as measured.
        import csv
        context_path=case/'results/sasa_comparison.csv'
        with context_path.open(encoding='utf-8-sig',newline='') as stream:
            conflicting={int(row['position']) for row in csv.DictReader(stream)
                if row.get('native_context_classification_changed')=='True'}
        by_id={c['candidate_id']:c for c in catalog}
        for cid in ids:
            if cid in by_id and any(m['residue']['uniprot_pos'] in conflicting for m in by_id[cid]['mutations']):
                retired[cid]={'reason':'NATIVE_NORMALIZED_SASA_THRESHOLD_CONFLICT',
                    'mutations':[m['notation'] for m in by_id[cid]['mutations']],
                    'sasa_comparison_sha256':sha256_file(context_path)}
        if retired:write_json(retired_path,retired)
        ids=[cid for cid in ids if cid not in retired]
        completed-=set(retired);additions=set(ids)-completed
        origins={cid:row['conditional_backup_origin'] for cid,row in history.items()
            if cid in ids and row.get('conditional_backup_origin')}
        origins.update({row['candidate_id']:row for row in base.get('conditional_pair_backups',[]) if row['candidate_id'] in ids})
        base['conditional_pair_backups']=list(origins.values())
        # Backup variants are lightweight pool additions, never positive evidence.
        supplements=base.get('menu_supplements',[])
        if len(supplements)<5:
            supplement=choose_menu_supplement(case,supplements)
            if supplement:
                base['menu_supplements']=[*supplements,supplement]
                write_json(case/f'menu_expansion_{len(supplements)}.json',supplement)
        # Never overwrite an earlier attempt, including failed logs/configs.
        prefix=f'T6_resume_{len(completed)}'
        attempt=stage_attempt_number(case,prefix)
        if (case/'frontier_replenishment.json').is_file():
            write_json(stages/f'{prefix}_{attempt}_frontier_contract.json',{
                'frontier':read(case/'frontier_replenishment.json'),
                'resource_contract':read(case/'resource_contract_amendment.json') if (case/'resource_contract_amendment.json').is_file() else None,
                'estimated_seconds_per_group':cost,'estimated_bytes_per_group':byte_cost})
        started=time.monotonic();before=owned_bytes(case)
        for phase in ('calibration','production'):
            remaining=case_remaining(case)
            report=launch(case,assets,f'{prefix}_{attempt}_{phase}',{
                **base,'phase':phase,'candidate_ids':ids,'confirmation_count':len(ids),'pool_cap':cap},
                max(0,remaining-1800))
            reports.append(report)
            if report['engineering_status']!='COMPLETED':return reports
            preserve_pair_index(case)
        completed.update(ids);pending=None;batch_number+=1
        history=record_completed_snapshot(case,reports[-1])
        physical=reports[-1].get('physical',{})
        pair_counts={cid:len(rows) for cid,rows in read(case/'results/surface_perturbation/paired_manifest.json').items()}
        latest_calibration=reports[-2].get('physical',{})
        calibration_counts={row['candidate']:row.get('calibration_runs',0)
            for row in latest_calibration.get('trajectories',[]) if 'candidate' in row}
        improvable=[cid for cid,row in latest_calibration.get('calibration_unknowns',{}).items() if extendable(cid,row)]
        base['calibration_extension_candidate_ids']=improvable[:2]
        cost=max(cost,(time.monotonic()-started)/max(1,len(additions)))
        byte_cost=max(byte_cost,max(1,owned_bytes(case)-before)/max(1,len(additions)))
        if reports[-1].get('physical',{}).get('goal_status')=='ACHIEVED':
            reports.extend(review_final_triples(case,assets,base,list(completed),max(0,remaining-1800)))
            break
    return reports


def paired_execution_identities(pairs):
    """Keep independent executions distinct even when their model bytes match."""
    identities={}
    for cid,rows in pairs.items():
        executions=set()
        evidence=set()
        for pair in rows:
            execution=(pair.get('batch_path'),pair['group'],pair['run'])
            if not execution[0]:
                raise ValueError('Missing paired model batch identity')
            if execution in executions:
                raise ValueError('Duplicate paired model execution identity')
            executions.add(execution)
            evidence.add((*execution,pair['wt_sha256'],pair['mut_sha256']))
        identities[cid]=evidence
    return identities


def final_panel_audit(case,assets,reports):
    """Audit the saved formal panel without adding models or extending null cohorts."""
    summary=read(case/'results/surface_perturbation_preview.json')
    if not summary.get('selected_candidate_ids'):return reports
    completed=[r for r in reports if r.get('engineering_status')=='COMPLETED'
        and r.get('configuration',{}).get('phase')=='production']
    if not completed:raise ValueError('No completed production for final APBS audit')
    config=dict(completed[-1]['configuration'])
    for key in ('stage','seconds','failure'):config.pop(key,None)
    ids=sorted(read(case/'results/surface_perturbation/paired_manifest.json'))
    config.update(phase='production',candidate_ids=ids,confirmation_count=len(ids),
        pool_cap=max(32,len(ids)),schedule=[5],calibration_extension_candidate_ids=[],
        final_panel_audit=True,stage_limit=None)
    index=stage_attempt_number(case,'T6_final_apbs',('',))
    remaining=case_remaining(case)
    # Keep five minutes for portable export and independent verification.
    if remaining<=300:return reports
    before=paired_execution_identities(read(case/'results/surface_perturbation/paired_manifest.json'))
    report=launch(case,assets,f'T6_final_apbs_{index}',config,remaining-300)
    reports.append(report)
    after=read(case/'results/surface_perturbation/paired_manifest.json')
    if before!=paired_execution_identities(after):raise ValueError('Final APBS stage changed paired model identities')
    write_json(case/'final_apbs_model_reuse_audit.json',{'status':'PASS',
        'new_paired_models':0,'candidate_count':len(after),
        'pair_count':sum(map(len,after.values())),'stage':report['stage']})
    return reports


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-root',type=Path,required=True)
    parser.add_argument('--assets-root',type=Path,required=True)
    parser.add_argument('--worker',type=Path)
    parser.add_argument('--preflight-only',action='store_true')
    parser.add_argument('--continue-after-pilot',action='store_true')
    parser.add_argument('--resume-batches',action='store_true')
    parser.add_argument('--max-new-batches',type=int)
    parser.add_argument('--finalize-current',action='store_true')
    parser.add_argument('--activate-amendment',action='store_true')
    args=parser.parse_args();case=args.case_root.resolve();assets=args.assets_root.resolve()
    if case==assets or assets in case.parents:raise ValueError('Task output must be separate from source assets')
    case.mkdir(parents=True,exist_ok=True)
    if args.worker:worker(case,assets,read(args.worker))
    elif args.preflight_only:preflight(case,assets)
    elif args.activate_amendment:activate_resource_amendment(case)
    elif args.finalize_current:
        reports=[read(p) for p in sorted((case/'stages').glob('*_report.json'),key=lambda p:p.stat().st_mtime)]
        reports=final_panel_audit(case,assets,reports)
        finalize(case,reports)
    else:
        reports=(resume_batches(case,assets,args.max_new_batches) if args.resume_batches else
            continue_after_pilot(case,assets) if args.continue_after_pilot else run(case,assets))
        if args.max_new_batches is not None or reports[-1]['engineering_status']!='COMPLETED':
            write_json(case/'resume_checkpoint.json',{'stage':reports[-1]['stage'],
                'engineering_status':reports[-1]['engineering_status'],
                'case_budget':read_budget(case/'case_budget.json'),
                'interpretation':'Intermediate checkpoint; no new final archive created'})
        else:
            finalize(case,final_panel_audit(case,assets,reports))


if __name__=='__main__':main()
