"""Isolated coverage-v3 contract; legacy physical qualification is unchanged."""
from dataclasses import asdict
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse import save_npz

from ..util import write_json, sha256_file
from .paired import atoms, validate_pair
from .tools import PhysicalTools, fingerprint, make_grid
from .reference import canonical_reference
from .fields import compute_pair, reference_id, aligned_atoms
from .scoped_fields import complete_scope
from .coverage import scan_support, budget_curve
from .coverage_upper import hard_scan_upper_bound
from .sampling_v3 import validate_roles
from .null_v3 import obtain
from .classification import classify_fields
from .combo_contract import site_evidence, actual_pair_controls
from .calibration_scope import local_scope
from .panel_v3 import select, metrics, compression_information, feasible_target
from .runtime import ACTIVE_BUDGET, RunBudget
from .local_constructs import local_construct, export_local_constructs

CONTRACT='positive_lower_bound_information_v3'


def coupling(components):
    contact=components.get('minimum_native_heavy_atom_distance_a')
    overlap=components.get('proxy_support_overlap_fraction')
    if contact is None or overlap is None:return 'UNKNOWN'
    if contact<=4.5:return 'CONTACT_COUPLED_PROXY'
    if overlap==0:return 'SPATIALLY_SEPARATED_PROXY'
    return 'INTERMEDIATE_PROXY'


def _restore(path):
    return json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {}


def run(**kwargs):
    from .engine import direct_support, direct_site_masks, cached_numerical_bounds
    output=Path(kwargs['output']);settings=kwargs['settings'];folder=output/'surface_perturbation'
    folder.mkdir(parents=True,exist_ok=True)
    budget=ACTIVE_BUDGET.get() or RunBudget(folder/'production_budget.json',7200)
    tools=PhysicalTools(settings.get('tool_cache',folder/'tool_cache'),settings,budget)
    from .deadlines import guard_mesh
    from .query_cache import memoize_mesh
    ses=tools.ses
    tools.ses=lambda *a,**k:memoize_mesh(guard_mesh(ses(*a,**k),budget),tools)
    records=kwargs['records'];native=atoms(Path(kwargs['structure']))
    all_pairs=_restore(folder/'paired_manifest.json')
    # Import only matching v3 roles; old physical/manifests stay historical.
    rows={tuple(sorted(r['mutations'])):r for r in kwargs['paired_records']}
    summaries={};classified={};lower={};upper={};direct={};scopes={};errors=[];qualified=[];final=[]
    summary={'strategy':'surface_perturbation_v1_preview','phase':settings.get('phase','production'),
        'status':'NO_VERIFIED_PANEL','panel_contract':CONTRACT,'goal_status':'NOT_EVALUATED',
        'hard_target_status':'NOT_EVALUATED','search_status':'OPEN','candidate_site_evidence':{},
        'method_contract':{'version':CONTRACT,'scan_center_geodesic_radius_a':12.,
                           'production_repeats':5,'null_train_WT':10,'null_holdout_WT':10},
        'interpretation':'computational scan opportunity and geometric codes; not antibody response'}
    try:
        reliable={r.uniprot_pos for r in records if r.coordinate_eligible and not r.review_reasons}
        designable={o.residue.uniprot_pos for o in kwargs['options'] if o.eligible and not o.review_required
                    and o.mutant_aa not in 'CGP' and o.residue.one_letter not in 'CGP'}
        mesh=tools.ses(native)
        reference=canonical_reference(mesh,native,records,reliable,designable,spacing=2.)
        reference['context_identity']=fingerprint({'native_input':sha256_file(Path(kwargs['structure'])),
            'source':settings.get('reference_input',{}),
            'mapping':[(r.chain_id,r.auth_seq_id,r.insertion_code,r.uniprot_pos,r.coordinate_eligible,r.review_reasons,r.domain) for r in records]})
        area=reference['area'];graph=reference['graph'];n=len(area)
        summary.update(reference_area_all_a2=float(area.sum()),reference_area_designable_a2=float(area[reference['designable']].sum()))
        save_npz(folder/'reference_graph.npz',graph)
        np.savez_compressed(folder/'canonical_reference.npz',**{k:v for k,v in reference.items() if isinstance(v,np.ndarray)})
        np.savez_compressed(folder/'native_WT_SES_mesh.npz',vertices=mesh.vertices,faces=mesh.faces,triangle_area=mesh.area,triangle_normals=mesh.normals)
        write_json(folder/'reference_contract.json',{'reference_id':reference_id(reference),
            'area_all_a2':float(area.sum()),'input_sha256':sha256_file(Path(kwargs['structure'])),
            'all_denominator':'fixed reliable external WT surface, including protected positions',
            'reference_spacing_a':2.,'scan_radius_a':12.,'basis':'SES_NATIVE'})
        # Candidate IDs are upstream-whitelisted. Cheap scores only order this
        # queue; every deferred candidate remains explicitly pending.
        candidates=sorted(kwargs['candidates'],key=lambda c:(
            c.order,not {0,1,2,3,4}<={p.get('sampling_slot') for p in all_pairs.get(c.candidate_id,[])
                if p.get('sampling_role')=='production' and p.get('sampling_cohort')=='production_v3'},
            -len(c.covered_positions_by_scenario.get('geodesic_12',[])),c.structural_risk,c.candidate_id))
        capacity=int(settings.get('confirmation_batch',8));spent=0
        for candidate_index,candidate in enumerate(candidates):
            cid=candidate.candidate_id
            from itertools import combinations
            native_sites=[];ca=[]
            for option in candidate.mutations:
                r=option.residue;key=(r.chain_id,int(r.auth_seq_id),(r.insertion_code or '').strip())
                native_sites.append(np.array([v[0] for k,v in native.items() if k[:3]==key and v[2]!='H']))
                ca.append(native.get((*key,'CA')))
            distances=[float(np.linalg.norm(a[:,None]-b[None,:],axis=2).min()) for a,b in combinations(native_sites,2) if len(a) and len(b)]
            if distances:candidate.risk_components['minimum_native_heavy_atom_distance_a']=min(distances)
            ca_distances=[float(np.linalg.norm(a[0]-b[0])) for a,b in combinations(ca,2) if a is not None and b is not None]
            candidate.risk_components['CA_distances_a']=ca_distances
            candidate.coupling_class=coupling(candidate.risk_components)
            state={'candidate_id':cid,'engineering_status':'PENDING','scientific_status':'UNKNOWN',
                   'coupling_class':candidate.coupling_class,'priority_kind':'proxy_not_upper_bound'}
            state['coupling_diagnostics']=candidate.risk_components
            summaries[cid]=state
            source=rows.get(tuple(sorted(kwargs['notation'](candidate))),{}).get('paired_models',[])
            previous=all_pairs.setdefault(cid,[])
            known={(p['batch_path'],p['group'],p['run']) for p in previous}
            previous.extend(p for p in source if (p['batch_path'],p['group'],p['run']) not in known
                            and p.get('sampling_protocol')=='time_rotabase_v3')
            production=validate_roles(previous,'production',5,'production_v3')
            if candidate.order==3 and actual_pair_controls(candidate,[c for c in candidates
                    if c.order==2 and c.candidate_id in lower and c.formal_eligible]) is None:
                state.update(engineering_status='COMPLETED',scientific_status='EXACT_DOUBLE_CONTROLS_PENDING',
                             expensive_confirmation='NOT_PAID_MISSING_AUDITED_CONTROLS')
                continue
            if len(production)<5 and spent>=capacity:
                state['scientific_status']='DEFERRED_NOT_EXCLUDED';continue
            started=time.monotonic()
            try:
                if len(production)<5:spent+=1
                for count in (1,3,5):
                    if len(production)<count:
                        production+=kwargs['supplement'](candidate,count-len(production),slot_start=len(production))
                        ids={(p['batch_path'],p['group'],p['run']) for p in previous}
                        previous.extend(p for p in production if (p['batch_path'],p['group'],p['run']) not in ids)
                        write_json(folder/'paired_manifest.json',all_pairs)
                if len(validate_roles(production,'production',5,'production_v3'))!=5:
                    raise ValueError('five production sampling slots missing')
                for p in production:validate_pair(p)
                if not kwargs['recheck'](candidate,production):
                    state.update(engineering_status='COMPLETED',scientific_status='STABILITY_EXCLUDED');continue
                registration=[aligned_atoms(p,native,records)[2]['core_rmsd_a'] for p in production]
                state['production_core_rmsd_a']=registration
                state['production_content_diversity']={
                    'WT':len({p['wt_sha256'] for p in production}),
                    'mutant':len({p['mut_sha256'] for p in production}),
                    'paired_content':len({(p['wt_sha256'],p['mut_sha256']) for p in production})}
                if any(not np.isfinite(v) or v>.5 for v in registration):
                    candidate.formal_eligible=False
                    candidate.formal_exclusion_reasons.append('PAIRED_PHYSICAL_MODEL_QC_FAILED')
                    state.update(engineering_status='COMPLETED',scientific_status='PAIRED_PHYSICAL_MODEL_QC_FAILED')
                    continue
                support=np.array([direct_support(p,native,records,reference) for p in production])
                required=support.any(axis=0)
                upper[cid]=hard_scan_upper_bound(support,graph)
                site_support=np.array([direct_site_masks(p,native,records,reference) for p in production])
                if site_support.shape[:2]!=(5,candidate.order) or not np.all(np.any(site_support.sum(axis=0)>=4,axis=1)):
                    state.update(engineering_status='COMPLETED',scientific_status='NO_POSSIBLE_SITE_QUORUM');continue
                runs=[complete_scope(p,compute_pair(p,native,records,reference,tools,folder/'fields'),native,records,reference,tools,required,folder/'fields') for p in production]
                numerical=cached_numerical_bounds(production[0],native,records,reference,tools,runs[0],folder/'fields',required_mask=required)
                protocol,calibration=obtain(candidate,production,previous,kwargs['supplement'],native,records,
                    reference,tools,required,numerical,folder/'calibration',folder/'fields')
                state['null_content_diversity']={role:{
                    'logical_WT':10,'distinct_WT':len({p['wt_sha256'] for p in previous
                        if p.get('sampling_role')==role and p.get('sampling_cohort')==protocol})}
                        for role in ('null_train','null_holdout')}
                write_json(folder/'paired_manifest.json',all_pairs)
                result=classify_fields(runs,production,native,records,reference,calibration,protocol,run_support=support)
                classified[cid]=result;lower[cid]=result['hard'];direct[cid]=result['direct']
                if np.any(lower[cid]&~upper[cid]):raise ValueError('LOWER_OUTSIDE_PROVEN_UPPER')
                evidence=site_evidence(candidate,result['site_changes'],area)
                from .witness_v3 import witnesses
                witness=witnesses(runs,site_support,calibration)
                for site,changes in enumerate(result['site_changes']):
                    supported=np.any([value!=0 for name,value in witness.items() if name.startswith(f'site{site}__')],axis=0)
                    if np.any(changes & ~supported):raise ValueError('POSITIVE_WITHOUT_SAME_FOUR_WITNESS')
                np.savez_compressed(folder/f'{cid}_quorum_witnesses.npz',**witness)
                summary['candidate_site_evidence'][cid]=evidence
                scope=local_scope(calibration,required,area)
                scope.update(qualification_contract=CONTRACT,positive_lower_bound_qualified=bool(evidence['complete'] and lower[cid].any()),
                    holdout_validated=True,hard_scan_bound_resolved=not bool(np.any(upper[cid]&~lower[cid])),
                    measurement_local_complete=not bool(np.any(~result['known']&required)),
                    footprint_characterization_status='COMPLETE' if scope['complete'] and not np.any(~result['known']&required) else 'PARTIAL')
                scopes[cid]=scope
                candidate.positive_lower_bound_qualified=scope['positive_lower_bound_qualified']
                candidate.footprint_characterization_status=scope['footprint_characterization_status']
                np.savez_compressed(folder/f'{cid}_coverage.npz',direct=direct[cid],hard=lower[cid],upper=upper[cid],
                    unknown=~result['known'],direct_support_runs=support,direct_support=result['repeat_support'])
                np.savez_compressed(folder/f'{cid}_direct_site_evidence.npz',site_changes=result['site_changes'])
                np.savez_compressed(folder/f'{cid}_channel_resolution.npz',**result['channel_known'])
                np.savez_compressed(folder/f'{cid}_direction_evidence.npz',**result['vector_evidence'])
                write_json(folder/f'{cid}_field_audit.json',{'candidate_id':cid,'protocol':protocol,
                    'production_execution_ids':[(p['batch_path'],p['group'],p['run']) for p in production],
                    'production_runs':5,'null_runs':20,'scope':scope,
                    'physical_area_a2':float(area[direct[cid]].sum()),'scan_area_a2':float(area[lower[cid]].sum()),
                    'unknown_area_a2':float(area[~result['known']].sum())})
                if scope['positive_lower_bound_qualified']:qualified.append(candidate)
                state.update(engineering_status='COMPLETED',scientific_status='POSITIVE_LOWER_BOUND_QUALIFIED'
                    if scope['positive_lower_bound_qualified'] else 'UNRESOLVED_OR_NO_SITE_POSITIVE',
                    physical_seconds=time.monotonic()-started,production_models=5,null_models=20)
            except (RuntimeError,ValueError) as exc:
                state.update(engineering_status='FAILED',scientific_status='UNKNOWN',reason=str(exc))
                errors.append({'candidate':cid,'reason':str(exc),'classification':'ENGINEERING_FAILURE_NOT_NEGATIVE'})
                write_json(folder/'paired_manifest.json',all_pairs)
            write_json(folder/'cascade_checkpoint.json',summaries)
            cold_or_compression_triple=any(c.order==3 and (c.triple_role=='TRIPLE_COMPRESSION' or
                not {0,1,2,3,4}<={p.get('sampling_slot') for p in all_pairs.get(c.candidate_id,[])
                    if p.get('sampling_role')=='production' and p.get('sampling_cohort')=='production_v3'}) for c in candidates)
            witness_candidates=[c for c in qualified if c.order==2]
            if not errors and not cold_or_compression_triple and len(witness_candidates)>=4:
                ids=[c.candidate_id for c in witness_candidates]
                witness=feasible_target(np.array([lower[cid] for cid in ids]),np.array([upper[cid] for cid in ids]),
                    area,graph,[c.structural_risk for c in witness_candidates],ids)
                if witness is not None:
                    summary['early_target_witness_ids']=[ids[i] for i in witness]
                    for remaining in candidates[candidate_index+1:]:
                        rcid=remaining.candidate_id
                        pending=all_pairs.setdefault(rcid,[]);known={(p['batch_path'],p['group'],p['run']) for p in pending}
                        source=rows.get(tuple(sorted(kwargs['notation'](remaining))),{}).get('paired_models',[])
                        pending.extend(p for p in source if (p['batch_path'],p['group'],p['run']) not in known
                                       and p.get('sampling_protocol')=='time_rotabase_v3')
                        summaries[rcid]={'engineering_status':'COMPLETED','scientific_status':'TARGET_REACHED_NOT_CONFIRMED',
                            'expensive_confirmation':'NOT_PAID_VERIFIED_FEASIBLE_PANEL_ALREADY_AVAILABLE',
                            'interpretation':'not a scientific rejection or scope exhaustion'}
                    break
        # Controls need stable, measured intervals; a positive-only parent mask
        # is never substituted for its unresolved upper footprint.
        doubles=[c for c in qualified if c.order==2]
        eligible=list(doubles)
        summary['triple_audits']=[]
        for c in (c for c in qualified if c.order==3):
            controls=actual_pair_controls(c,[x for x in candidates if x.order==2 and x.candidate_id in lower and x.formal_eligible])
            audit={'candidate_id':c.candidate_id,'role':c.triple_role,'accepted':False,'reasons':[]}
            if controls is None:audit['reasons'].append('EXACT_DOUBLE_CONTROLS_PENDING')
            elif c.triple_role=='TRIPLE_COMPRESSION':eligible.append(c)
            else:
                gains={p.candidate_id:float(area[lower[c.candidate_id]&~upper[p.candidate_id]].sum()) for p in controls}
                gain=min(gains.values())
                baseline=np.any([lower[p.candidate_id] for p in doubles],axis=0) if doubles else np.zeros(n,bool)
                audit.update(exact_double_control_ids=[p.candidate_id for p in controls],
                    gain_vs_pair_upper_a2=gains,minimum_required_gain_a2=.05*area.sum(),
                    actual_added_vs_double_pool_a2=float(area[lower[c.candidate_id]&~baseline].sum()))
                if gain>=.05*area.sum() and np.any(lower[c.candidate_id]&~baseline):
                    c.warnings.append('REAL_TRIPLE_EXTRA_GAIN_VERIFIED');eligible.append(c);audit['accepted']=True
                else:audit['reasons'].append('TRIPLE_GAIN_LOWER_BOUND_NOT_DEMONSTRATED')
            summary['triple_audits'].append(audit)
        regular=[c for c in eligible if c.triple_role!='TRIPLE_COMPRESSION']
        lo=np.array([lower[c.candidate_id] for c in regular],bool).reshape(-1,n)
        up=np.array([upper[c.candidate_id] for c in regular],bool).reshape(-1,n)
        indices,summary['panel_optimizer']=select(lo,up,area,graph,[c.structural_risk for c in regular],[c.candidate_id for c in regular],solver_seconds=min(5.,budget.remaining()))
        witness_ids=summary.get('early_target_witness_ids')
        if witness_ids and metrics(lo[indices],up[indices],area)['C_screen_lower_bound']<.8:
            lookup={c.candidate_id:i for i,c in enumerate(regular)}
            proposal=[lookup[cid] for cid in witness_ids]
            if 4<=len(proposal)<=12 and metrics(lo[proposal],up[proposal],area)['C_screen_lower_bound']>=.8:
                indices=proposal;summary['panel_optimizer']['verified_witness_fallback']=True
        final=[regular[i] for i in indices]
        baseline=np.any([lower[c.candidate_id] for c in final],axis=0) if final else np.zeros(n,bool)
        from .compression import compression_audit
        summary['compression_audits']=[]
        for c in (c for c in eligible if c.triple_role=='TRIPLE_COMPRESSION'):
            ids=[x.candidate_id for x in final]
            if not set(c.parent_candidate_ids)<=set(ids):continue
            audit=compression_audit(ids,c.parent_candidate_ids,c.candidate_id,lower,area,graph,baseline)
            information=compression_information(ids,audit['after_ids'],lower,upper,area)
            audit['information']=information; audit['accepted'] &= information['accepted']
            audit['reasons']+=information['reasons'];summary['compression_audits'].append(audit)
            if audit['accepted']:
                c.warnings.append('REAL_TRIPLE_COMPRESSION_VERIFIED')
                final=[x for x in final if x.candidate_id not in c.parent_candidate_ids]+[c]
        lo=np.array([lower[c.candidate_id] for c in final],bool).reshape(-1,n)
        up=np.array([upper[c.candidate_id] for c in final],bool).reshape(-1,n)
        di=np.array([direct[c.candidate_id] for c in final],bool).reshape(-1,n)
        summary['panel_information']=metrics(lo,up,area,graph,di)
        target=4<=len(final)<=12 and summary['panel_information']['C_screen_lower_bound']>=.8
        summary['hard_target_status']='ACHIEVED' if target else 'BELOW_TARGET'
        summary['goal_status']='ACHIEVED' if target else 'COVERAGE_BELOW_TARGET'
        summary['status']='TARGET_HARD_VERIFIED_SEARCH_OPEN' if target else 'CALIBRATION_EXPANSION_REQUIRED'
        extended=np.full(lo.shape,np.nan);apbs_reports={};completed=[]
        # Persist the accepted lower bound before the independent final audit.
        summary['candidate_local_calibration']=scopes
        summary['verified_local_constructs']=[local_construct(asdict(c),validate_roles(all_pairs[c.candidate_id],'production',5,'production_v3'),
            scopes[c.candidate_id],area,lower[c.candidate_id],direct[c.candidate_id]) for c in final]
        np.savez_compressed(folder/'coverage_matrix.npz',hard=lo,upper=up,direct=di,extended=extended,
            candidate_ids=np.array([c.candidate_id for c in final]),area=area,designable=reference['designable'])
        requested=bool(final) and (target or settings.get('final_panel_audit',False))
        if requested:
            summary['final_apbs_evidence']={'status':'INCOMPLETE','execution_status':'STARTED',
                'preparation_status':'NOT_COMPLETED','numerical_resolution_status':'NOT_EVALUATED'}
            from .electrostatics import prepare_candidate,audit_candidate
            prepared={};xyz=[]
            for c in final:
                production=validate_roles(all_pairs[c.candidate_id],'production',5,'production_v3')
                prepared[c.candidate_id]=prepare_candidate(c.candidate_id,production,native,records,tools,folder)
                xyz.extend(prepared[c.candidate_id][2])
            grid=make_grid(xyz,.5)
            for i,c in enumerate(final):
                production=validate_roles(all_pairs[c.candidate_id],'production',5,'production_v3')
                changes,missing,audit=audit_candidate(c.candidate_id,production,native,records,reference,tools,folder,
                    convergence=i==0,prepared_data=prepared[c.candidate_id],shared_grid=grid)
                completed.append(c.candidate_id);apbs_reports[c.candidate_id]=audit
                extended[i]=changes;extended[i,missing]=np.nan
            from .final_evidence import require_apbs_audits
            summary['final_apbs_evidence']=require_apbs_audits([c.candidate_id for c in final],completed,reports=apbs_reports)
        else:
            summary['final_apbs_evidence']={'status':'NOT_REQUESTED_TARGET_NOT_MET','execution_status':'NOT_REQUESTED',
                'preparation_status':'NOT_REQUESTED','numerical_resolution_status':'NOT_EVALUATED'}
        summary['coverage_curve']=budget_curve(lo,extended,list(range(len(final))),area,reference['designable'],graph)
        np.savez_compressed(folder/'coverage_matrix.npz',hard=lo,upper=up,direct=di,extended=extended,
            candidate_ids=np.array([c.candidate_id for c in final]),area=area,designable=reference['designable'])
        np.savez_compressed(folder/'verified_pool_v3.npz',hard=np.array([lower[c.candidate_id] for c in regular],bool).reshape(-1,n),
            upper=np.array([upper[c.candidate_id] for c in regular],bool).reshape(-1,n),candidate_ids=np.array([c.candidate_id for c in regular]),area=area)
        summary['verified_local_constructs']=[local_construct(asdict(c),validate_roles(all_pairs[c.candidate_id],'production',5,'production_v3'),
            scopes[c.candidate_id],area,lower[c.candidate_id],direct[c.candidate_id]) for c in final]
        summary['candidate_local_calibration']=scopes
    except (TimeoutError,RuntimeError,ValueError) as exc:
        summary['status']='INCOMPLETE';errors.append({'reason':str(exc)})
        if summary.get('final_apbs_evidence',{}).get('execution_status')=='STARTED':
            summary['final_apbs_evidence'].update(execution_status='FAILED',reason=str(exc))
    finally:
        summary.update(selected_candidate_ids=[c.candidate_id for c in final],actual_constructs=len(final),
            errors=errors,tool_calls=tools.calls,tool_activity=tools.activity,elapsed_seconds=budget.elapsed,
            cascade=summaries,trajectories=[])
        if any(v.get('engineering_status')=='FAILED' for v in summaries.values()):summary['status']='INCOMPLETE'
        write_json(folder/'paired_manifest.json',all_pairs)
        write_json(output/'surface_perturbation_preview.json',summary)
        export_local_constructs(output,summary.get('verified_local_constructs',[]),summary['status'])
        if settings.get('export_mode')=='deliverable':
            from .export import evidence_archive
            evidence_archive(output,all_pairs,folder/'calibration')
        budget.save(ACTIVE_BUDGET.get() is not None)
    return final,summary
