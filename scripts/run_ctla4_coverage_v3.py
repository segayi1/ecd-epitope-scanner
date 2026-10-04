"""Coverage-v3 on the existing charged case, with an isolated result history."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import time
from zipfile import ZipFile,ZIP_DEFLATED

import numpy as np

from ecd_epitope_scanner.scripts import run_ctla4_coverage_v2 as old
from ecd_epitope_scanner.util import write_json,sha256_file
from ecd_epitope_scanner.surface_perturbation.runtime import read_budget,RunBudget
from ecd_epitope_scanner.surface_perturbation.case_resources import StorageGuard


def initialize(case,assets):
    if read_budget(case/'case_budget.json').get('active'):
        raise ValueError('v3 requires a closed production boundary')
    route=case/'routes/coverage_v3';output=route/'results'
    if not output.exists():
        output.mkdir(parents=True)
        for name in ('mutation_menu_v2.json','surface_graph_v4.json','uniprot.json','sasa_comparison.csv','candidate_pool_v2.json'):
            shutil.copy2(case/'results'/name,output/name)
        for name in ('resource_contract_amendment.json','site_pair_search_contract.json'):
            shutil.copy2(case/name,route/name)
        write_json(route/'source_contract.json',{'profile':'coverage_scan_v3_preview',
            'shared_case_budget':str(case/'case_budget.json'),'old_ledger_sha256':sha256_file(case/'case_budget.json'),
            'frozen_scope_sha256':sha256_file(case/'site_pair_search_contract.json'),
            'historical_results':str(case/'results'),'assets_root':str(assets),
            'interpretation':'old classifications not inherited; no refund or ledger reset'})
    return route


def replay(case,route):
    from ecd_epitope_scanner.joint_design_v4.stability_v3 import decision
    source=case/'audits/scientific_methods_20261004/research_data.json'
    data=old.read(source)
    restored=[r for r in data['stability_17'] if decision(r['spurs'],r['foldx'])=='CONSENSUS_PASS']
    assert len(restored)==17
    with (case/'audits/screening_96/screening_audit_96.csv').open(encoding='utf-8-sig',newline='') as stream:
        rows=list(csv.DictReader(stream))
    mild=[r for r in rows if r['stability_status']=='COMBINATION_REJECT' and .5<float(r['spurs_multi_ddg'])<=1]
    assert len(mild)==39 and all(decision(float(r['spurs_multi_ddg']),None)=='SPURS_ONLY_REVIEW_REQUIRED' for r in mild)
    conflict=[r for r in rows if r['stability_status']=='MODEL_DISAGREEMENT' and float(r['foldx_ddg'])>2]
    assert len(conflict)==4 and all(decision(float(r['spurs_multi_ddg']),float(r['foldx_ddg']))=='STRUCTURAL_ENERGY_CONFLICT' for r in conflict)
    write_json(route/'T2_stability_replay.json',{'engineering_status':'PASS','restored_confirmation_opportunities':17,
        'mild_reinstated':39,'retained_conflicts':4,'new_scientific_calls':0,
        'source_sha256':sha256_file(source),'classification_inheritance':False})


def audit_current(route):
    from ecd_epitope_scanner.surface_perturbation.sampling_v3 import validate_roles
    from ecd_epitope_scanner.surface_perturbation.paired import validate_pair
    output=route/'results';folder=output/'surface_perturbation'
    summary=old.read(output/'surface_perturbation_preview.json')
    pairs=old.read(folder/'paired_manifest.json')
    for cid,rows in pairs.items():
        for p in rows:validate_pair(p)
        production=validate_roles(rows,'production',5,'production_v3')
        # Scientific exclusions may have only the first medium-cost model.
        if summary.get('cascade',{}).get(cid,{}).get('null_models')==20:
            assert len(production)==5
            groups={p['sampling_cohort'] for p in rows if p['sampling_role']=='null_train'}
            assert len(groups)==1
            cohort=next(iter(groups))
            assert len(validate_roles(rows,'null_train',10,cohort))==10
            assert len(validate_roles(rows,'null_holdout',10,cohort))==10
    if (folder/'coverage_matrix.npz').exists():
        with np.load(folder/'coverage_matrix.npz',allow_pickle=False) as arrays:
            assert not np.any(arrays['hard']&~arrays['upper'])
            actual=float(arrays['area'][np.any(arrays['hard'],axis=0)].sum()/arrays['area'].sum())
            assert np.isclose(actual,summary['panel_information']['C_screen_lower_bound'],atol=1e-12)
            assert np.isclose(arrays['area'].sum(),6874.363277945836,atol=1e-8,rtol=0)
    write_json(route/'current_evidence_audit.json',{'engineering_status':'PASS',
        'scientific_status':summary.get('hard_target_status'),'candidate_count':len(pairs),
        'actual_pairs':sum(map(len,pairs.values())),'summary_sha256':sha256_file(output/'surface_perturbation_preview.json')})
    return summary


def run_stage(case,assets,route,label,ids,**extra):
    from itertools import combinations
    import hashlib
    path=route/'results/panel_candidates.json'
    catalog={c['candidate_id']:c for c in old.read(path)} if path.exists() else {}
    history_path=route/'screened_candidate_history.json'
    if history_path.exists():
        for cid,row in old.read(history_path).items():catalog.setdefault(cid,row['candidate'])
    parents=set(extra.get('triple_parent_ids',[]));mandatory={tuple(sorted(s)) for s in extra.get('mandatory_candidate_signatures',[])}
    roles=dict(extra.get('triple_roles_by_candidate_id',{}));parent_map=dict(extra.get('triple_parents_by_candidate_id',{}))
    expanded=set(ids)
    for cid in ids:
        c=catalog.get(cid,{})
        if c.get('order')!=3:continue
        current_parents=extra.get('triple_parent_ids') if cid==extra.get('triple_candidate_id') else c.get('parent_candidate_ids',[])
        parent_map[cid]=current_parents;parents.update(current_parents);expanded.update(current_parents)
        roles[cid]=extra.get('triple_role') if cid==extra.get('triple_candidate_id') else c['triple_role']
        for options in combinations(c['mutations'],2):
            signature=tuple(sorted(m['notation'] for m in options));mandatory.add(signature)
            expanded.add('CV2_'+hashlib.sha256('|'.join(signature).encode()).hexdigest()[:16])
    if roles:
        extra.update(triple_parent_ids=sorted(parents),mandatory_candidate_signatures=[list(s) for s in sorted(mandatory)],
            triple_roles_by_candidate_id=roles,triple_parents_by_candidate_id=parent_map,
            allow_remote_triples=bool(extra.get('allow_remote_triples') or any(r['remote'] for r in
                (old.read(route/'v3_triple_attempts.json') if (route/'v3_triple_attempts.json').exists() else []) if r['candidate_id'] in roles)))
    ids=sorted(expanded)
    if len(ids)>256:raise ValueError('Retained triples/controls exceed active evidence bound')
    stage=f'V3_{label}_{old.stage_attempt_number(case,"V3_"+label,phases=("",))}'
    config={'phase':'production','optimization_profile':'coverage_scan_v3_preview','results_root':str(route/'results'),
        'confirmation_count':len(ids),'candidate_ids':ids,'stage_limit':None,'pool_cap':max(32,len(ids)),
        'schedule':[5],'full_replacement_search':False,**extra}
    result=old.launch(case,assets,stage,config,7200)
    if result['engineering_status']=='COMPLETED':
        audit_current(route)
        old.record_completed_snapshot(route,result)
        # record_completed_snapshot references only this route's result data.
        write_json(route/'last_completed_report.json',result)
    else:
        write_json(route/'pending_stage.json',{'ids':ids,'label':label,'extra':extra,'failed_stage':stage})
    return result


def finalize(case,assets,route):
    output=route/'results';summary=old.read(output/'surface_perturbation_preview.json')
    if summary.get('selected_candidate_ids'):
        result=run_stage(case,assets,route,'FINAL_APBS',summary['selected_candidate_ids'],final_panel_audit=True)
        if result['engineering_status']!='COMPLETED':return
    from ecd_epitope_scanner.surface_perturbation.export import evidence_archive
    from ecd_epitope_scanner.exporting import finalize_scientific_bundle
    budget=RunBudget(case/'case_budget.json',None)
    try:
        StorageGuard(case).check()
        pairs=old.read(output/'surface_perturbation/paired_manifest.json')
        evidence_archive(output,pairs,output/'surface_perturbation/calibration')
        from .audit_ctla4_coverage_v3 import audit
        write_json(output/'independent_v3_semantic_audit.json',audit(route))
        source=Path(__file__).resolve().parents[1]
        pin=old.read(source/'source_overlay_manifest.json')
        with ZipFile(output/'implementation_source_v3.zip','w',ZIP_DEFLATED) as archive:
            for name,digest in pin['files'].items():
                path=source/name
                if sha256_file(path)!=digest:raise ValueError('Frozen source changed before delivery')
                archive.write(path,name)
            archive.write(source/'source_overlay_manifest.json','source_overlay_manifest.json')
        receipts=[old.read(p) for p in (case/'stages').glob('V3_*_report.json')]
        write_json(output/'v3_cost_and_execution_report.json',{
            'charged_seconds':sum(r['charged_seconds'] for r in receipts),
            'case_cumulative_seconds':budget.elapsed,'role_counts':audit(route)['role_counts'],
            'stages':[{k:r.get(k) for k in ('stage','engineering_status','charged_seconds','new_bytes','peak_owned_bytes')} for r in receipts],
            'actual_evidence_pairs':sum(map(len,pairs.values())),
            'interpretation':'role counts are retained logical pairs, not a sum of reused worker calls; elapsed includes failures'})
        # Checkpoint exports could not list an archive that did not yet exist.
        # Explicitly bind the newly generated evidence, source and audits.
        manifest=old.read(output/'run_manifest.json')
        for name in ('surface_perturbation_evidence.zip','implementation_source_v3.zip',
                     'independent_v3_semantic_audit.json','v3_cost_and_execution_report.json'):
            manifest['files'][name]=sha256_file(output/name)
        write_json(output/'run_manifest.json',manifest)
        StorageGuard(case).check()
        finalize_scientific_bundle(output)
        audit_current(route)
        write_json(route/'independent_final_package_audit.json',audit(route,packages=True))
    finally:budget.save(False)
    write_json(route/'delivery_status.json',{'scientific_status':old.read(output/'surface_perturbation_preview.json').get('hard_target_status'),
        'case_budget':read_budget(case/'case_budget.json'),'search':old.read(route/'search_state.json') if (route/'search_state.json').exists() else {},
        'interpretation':'archive integrity and scientific target are separate'})


def triple_review(case,assets,route,label):
    """Generate from measured concrete doubles; charge only original case."""
    summary=old.read(route/'results/surface_perturbation_preview.json')
    if summary.get('hard_target_status')=='ACHIEVED' and summary.get('actual_constructs',0)<=4:
        write_json(route/'T5_decision.json',{'status':'SKIPPED_NO_CONSTRUCT_COUNT_BENEFIT','actual_constructs':summary['actual_constructs']})
        return True
    retained=[cid for cid,row in summary.get('candidate_local_calibration',{}).items()
              if row.get('positive_lower_bound_qualified')]
    attempts_path=route/'v3_triple_attempts.json'
    attempts=old.read(attempts_path) if attempts_path.exists() else []
    for remote in (False,True):
        config=old.triple_test(route,assets,{},retained,7200,remote=remote,proposal_only=True)
        if not config:continue
        ids=config.pop('candidate_ids')
        for key in ('confirmation_count','pool_cap'):config.pop(key,None)
        result=run_stage(case,assets,route,label,ids,**config)
        if result['engineering_status']!='COMPLETED':return False
        attempts.append({'candidate_id':config['triple_candidate_id'],'role':config['triple_role'],
                         'remote':remote,'stage':result['stage']})
        write_json(attempts_path,attempts)
        write_json(route/'T5_gate.json',{'engineering_status':'PASS','scientific_status':result['physical'].get('hard_target_status'),
                                      'actual_attempts':attempts})
        return True
    write_json(route/'T5_gate.json',{'engineering_status':'PASS','scientific_status':'NO_LEGAL_TRIPLE_PROPOSAL',
        'interpretation':'No missing parents or controls are fabricated','actual_attempts':attempts})
    return True


def rank_pending(route,pool):
    """Deterministic lexicographic proxies, never exclusion or formal codes."""
    from scipy.sparse import load_npz
    from ecd_epitope_scanner.surface_perturbation.coverage import largest_blind
    from ecd_epitope_scanner.joint_design_v4.coverage_profile import vector
    folder=route/'results/surface_perturbation'
    with np.load(folder/'canonical_reference.npz',allow_pickle=False) as data:area,owner=data['area'],data['owner']
    with np.load(folder/'coverage_matrix.npz',allow_pickle=False) as data:lower=data['hard'];upper=data['upper']
    covered=np.any(lower,axis=0) if len(lower) else np.zeros(len(area),bool)
    graph=load_npz(folder/'reference_graph.npz')
    summary=old.read(route/'results/surface_perturbation_preview.json')
    costs=[r['physical_seconds'] for r in summary.get('cascade',{}).values() if r.get('physical_seconds')]
    cost=max(costs,default=0)
    menus=old.read(route/'results/mutation_menu_v2.json')['rows']
    descriptors={(r['residue']['uniprot_pos'],r['mutant_aa']):vector(r['chemical_descriptors']) for r in menus}
    codes=np.where(lower,1,np.where(~upper,0,-1)).T
    groups,inverse=np.unique(codes,axis=0,return_inverse=True)
    weights=np.bincount(inverse,weights=area,minlength=len(groups))
    scores={}
    for candidate in pool:
        hit=np.isin(owner,candidate['covered_positions_by_scenario'].get('geodesic_12',[]))
        # Proxy splits within existing ternary groups: this is not a verified
        # response code, and cannot justify a pruning decision.
        yes=np.bincount(inverse,weights=area*hit,minlength=len(groups))
        split=float(np.sum(yes*(weights-yes)))
        values=[descriptors[(m['residue']['uniprot_pos'],m['mutant_aa'])] for m in candidate['mutations']]
        novelty=float(np.linalg.norm(np.array(values[0])-np.array(values[1])))
        scores[candidate['candidate_id']]=(-float(area[hit&~covered].sum()),largest_blind(graph,covered|hit,area),
            -split,-novelty,candidate['structural_risk'],cost,candidate['candidate_id'])
    ordered=sorted(scores,key=scores.get)
    write_json(route/'next_batch_priority.json',{'candidate_ids':ordered,'scores':scores,
        'interpretation':'priority_proxy only; low-ranked pairs remain pending; no upper-bound pruning'})
    return ordered


def next_batch_size(case):
    costs=[state['physical_seconds'] for path in (case/'stages').glob('V3_*_report.json')
           for state in old.read(path).get('physical',{}).get('cascade',{}).values() if state.get('physical_seconds')]
    maximum=max(costs,default=1200.)
    return max(1,min(8,int(6600/maximum)))


def run(case,assets,pilot_only=False):
    route=initialize(case,assets);replay(case,route)
    history_path=route/'screened_candidate_history.json'
    # A boundary restart can occur after the charged worker receipt but before
    # the coordinator records its gate. Recover that exact completed stage.
    receipts=[old.read(p) for p in (case/'stages').glob('V3_T4*_report.json')]
    if not (route/'T4_gate.json').exists() and receipts:
        latest=max(receipts,key=lambda r:r['case_elapsed_seconds'])
        if latest['engineering_status']=='COMPLETED':
            audit_current(route);old.record_completed_snapshot(route,latest)
            write_json(route/'T4_gate.json',{'engineering_status':'PASS','stage':latest['stage'],
                'recovered_closed_boundary':True,'candidate_ids':latest['configuration']['candidate_ids']})
    receipts=[old.read(p) for p in (case/'stages').glob('V3_T5*_report.json')]
    if not (route/'T5_gate.json').exists() and receipts:
        latest=max(receipts,key=lambda r:r['case_elapsed_seconds'])
        if latest['engineering_status']=='COMPLETED':
            audit_current(route);old.record_completed_snapshot(route,latest)
            config=latest['configuration']
            attempt={'candidate_id':config['triple_candidate_id'],'role':config['triple_role'],
                     'remote':bool(config.get('allow_remote_triples')),'stage':latest['stage']}
            path=route/'v3_triple_attempts.json';attempts=old.read(path) if path.exists() else []
            if not any((r['candidate_id'],r['role'])==(attempt['candidate_id'],attempt['role']) for r in attempts):attempts.append(attempt)
            write_json(path,attempts)
            write_json(route/'T5_gate.json',{'engineering_status':'PASS','stage':latest['stage'],
                'recovered_closed_boundary':True,'actual_attempts':attempts})
    history=old.read(history_path) if history_path.is_file() else {}
    # Preserve immutable stage exports; the live registry must reflect an
    # already-recorded five-model conflict instead of a stale prescreen flag.
    reaudit=[]
    for cid,row in history.items():
        candidate=row['candidate'];stability=candidate.get('stability') or {}
        if stability.get('consensus_status')=='STRUCTURAL_ENERGY_CONFLICT' and candidate.get('formal_eligible'):
            reaudit.append({'candidate_id':cid,'previous_stage':row['last_completed_stage'],
                            'reason':'FIVE_MODEL_CONFLICT_FLAG_SYNCHRONIZED','previous_candidate':candidate.copy()})
            candidate['formal_eligible']=False
            candidate['formal_exclusion_reasons']=sorted(set(candidate.get('formal_exclusion_reasons',[]))|{'STRUCTURAL_ENERGY_CONFLICT'})
    if reaudit:
        write_json(route/'stability_flag_reaudit.json',{'new_scientific_calls':0,'changes':reaudit,
                   'interpretation':'working qualification corrected; original stage exports preserved'})
        write_json(history_path,history)
    # A successful recorded stage is never recomputed after restart. Failed
    # stages retain logs and model slots; run_stage resumes their missing nodes.
    pending=route/'pending_stage.json'
    if pending.exists():
        job=old.read(pending)
        result=run_stage(case,assets,route,job['label'],job['ids'],**job['extra'])
        if result['engineering_status']!='COMPLETED':return
        pending.unlink()
        history=old.read(history_path)
    if not (route/'T3_gate.json').exists():
        source_history=old.read(case/'screened_candidate_history.json')
        cid=next(cid for cid,row in source_history.items() if
            {m['notation'] for m in row['candidate']['mutations']}=={'D78A','M122W'})
        result=run_stage(case,assets,route,'T3',[cid])
        if result['engineering_status']!='COMPLETED':return
        write_json(route/'T3_gate.json',{'engineering_status':'PASS','scientific_status':result['physical'].get('hard_target_status'),
            'stage':result['stage'],'interpretation':'unknown is permitted as a scientific result, not mislabeled PASS'})
        history=old.read(history_path)
    if pilot_only:return
    if not (route/'T4_gate.json').exists():
        pool=old.replenish_saved_frontier(route,set(history),256)
        with (route/'results/sasa_comparison.csv').open(encoding='utf-8-sig',newline='') as stream:
            boundary={int(row['position']) for row in csv.DictReader(stream) if .2<=float(row['native_rsasa'])<.25
                      and row['native_context_classification_changed']!='True'}
        ranked=rank_pending(route,pool);by_id={c['candidate_id']:c for c in pool}
        boundary_ids=[cid for cid in ranked if any(m['residue']['uniprot_pos'] in boundary for m in by_id[cid]['mutations'])]
        additions=(boundary_ids[:1]+[cid for cid in ranked if cid not in boundary_ids[:1]])[:3]
        result=run_stage(case,assets,route,'T4',sorted(set(history)|set(additions)))
        if result['engineering_status']!='COMPLETED':return
        write_json(route/'T4_gate.json',{'engineering_status':'PASS','stage':result['stage'],'new_position_pairs':3,
            'boundary_exposure_candidate_ids':boundary_ids[:1],'candidate_ids':additions})
        history=old.read(history_path)
    if not (route/'T5_gate.json').exists():
        if not triple_review(case,assets,route,'T5'):return
    batch=0
    while True:
        StorageGuard(case).check()
        summary=old.read(route/'results/surface_perturbation_preview.json')
        if summary.get('hard_target_status')=='ACHIEVED':
            if not triple_review(case,assets,route,'T6_COMPRESSION'):return
            write_json(route/'search_state.json',{'status':'TARGET_REACHED','scope_exhausted':False});break
        history=old.read(history_path)
        pool=old.replenish_saved_frontier(route,set(history),256)
        if not pool:
            seen={tuple(sorted(m['residue']['uniprot_pos'] for m in row['candidate']['mutations'])) for row in history.values() if row['candidate']['order']==2}
            frozen=set(map(tuple,old.read(case/'site_pair_search_contract.json')['site_pairs']))
            missing=frozen-seen
            write_json(route/'search_state.json',{'status':'EXHAUSTED' if not missing else 'INCOMPLETE',
                'total':len(frozen),'processed':len(seen&frozen),'remaining':[list(p) for p in sorted(missing)],
                'scope':'position pairs, not all replacements or triples'})
            if missing:raise ValueError('Empty candidate pool is not scope exhaustion')
            if not triple_review(case,assets,route,'T6_FINAL_TRIPLE'):return
            break
        # Keep formal masks active; unresolved historical cohorts remain in
        # the cumulative evidence index, never assigned fresh null retries.
        accepted_triples=set(summary.get('selected_candidate_ids',[]))|{r['candidate_id'] for r in summary.get('triple_audits',[]) if r.get('accepted')}
        orders={cid:row['candidate']['order'] for cid,row in history.items()}
        retained={cid for cid,row in summary.get('candidate_local_calibration',{}).items()
                  if row.get('positive_lower_bound_qualified') and (orders.get(cid)==2 or cid in accepted_triples)}
        backup=old.choose_failure_backup(route,set(history),scientific_v3=True)
        batch_size=next_batch_size(case)
        additions=rank_pending(route,pool)[:max(1,batch_size-int(backup is not None))]
        if batch_size==1:backup=None
        extra={}
        if backup:
            additions.append(backup['candidate_id'])
            extra={'mandatory_candidate_signatures':[[m['notation'] for m in backup['mutations']]],
                   'conditional_pair_backups':[{'candidate_id':backup['candidate_id'],
                        'parent_id':backup['conditional_backup_parent_id'],'reason':backup['conditional_backup_reason']}]}
        ids=sorted(retained|set(additions))
        if len(ids)>256:raise ValueError('Necessary evidence exceeds bounded active pool')
        batch+=1
        result=run_stage(case,assets,route,f'T6_{len(history)}',ids,**extra)
        if result['engineering_status']!='COMPLETED':return
        if batch%4==0 and not triple_review(case,assets,route,f'T6_TRIPLE_{len(history)}'):return
    finalize(case,assets,route)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-root',type=Path,required=True)
    parser.add_argument('--assets-root',type=Path,required=True)
    parser.add_argument('--pilot-only',action='store_true')
    args=parser.parse_args()
    run(args.case_root.resolve(),args.assets_root.resolve(),args.pilot_only)


if __name__=='__main__':main()
