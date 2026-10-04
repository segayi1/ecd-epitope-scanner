"""Opt-in physical preview. Missing calibration never becomes verified coverage.

The upstream engine owns substitution/stability gates; this engine owns SES
coordinates, physical repeats, matched nulls and measured coverage. FoldX
supplements are requested by a callback, never inferred from single fields.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import save_npz
from scipy.sparse.csgraph import connected_components
from ..util import write_json, sha256_file
from . import STRATEGY_ID
from .paired import atoms, validate_pair
from .reference import canonical_reference
from .fields import compute_pair, aligned_atoms, reference_id, IMPLEMENTATION_ID
from .null_store import FIELD_CHANNELS, numerical_bounds
from .scoped_fields import protocol_identity,save_calibration,complete_scope
from .calibration_jobs import usable_calibration,isolated_job,clear_failure
from .calibration import resolved_changes,calibrate_channel
from .coverage import scan_support, nested_panel, budget_curve, independent_bounds
from .tools import PhysicalTools, fingerprint
from .runtime import RunBudget, ACTIVE_BUDGET
from .buffers import repeated_field
from .chemistry import CHANNELS,DONORS,ACCEPTORS,RINGS
from .directions import VECTOR_GROUPS,resolve_vector


def load_calibration(folder,protocol):
    """Recompute qualification from matched raw nulls, including imported ones."""
    entries=usable_calibration(folder,protocol)
    for name,value in entries.items():
        if 'null_blocks' in value and 'numerical_error' in value:
            fresh=calibrate_channel(value['null_blocks'],value['numerical_error'],protocol,name,
                required_mask=value.get('required_mask'))
            # Retain source area-weighted diagnostics, not unit-weighted area
            # summaries. Candidate-local current areas are recomputed separately.
            keys=('status','valid','threshold','q95','q99','block_count','requested_scope_converged')
            entries[name]={**value,**{k:fresh[k] for k in keys if k in fresh},'raw_null_revalidated':True}
    return entries


def atom_property(value,key):
    _,resname,element=value;name=key[3]
    return (element,name in DONORS.get(resname,set()),name in ACCEPTORS.get(resname,set()),
            name in RINGS.get(resname,()))


def direct_support(pair, native, records, reference, distance=4.0):
    """Local query support from changed atoms at the mutated residues only.

    This conservative geometric support is an engineering contact test. It
    does not assert a universal radius of antibody binding disruption.
    """
    wt, mut, audit = aligned_atoms(pair, native, records)
    if audit['core_rmsd_a']>.5:
        raise ValueError('REFERENCE_UNRELIABLE: core RMSD exceeds 0.5 A engineering guard')
    residues = {tuple(x) for x in audit['mutated_author_residues']}
    changed = []
    for key in wt.keys() | mut.keys():
        if key[:3] not in residues: continue
        if key[3] in {'N','CA','C','O','OXT'}: continue
        left, right = wt.get(key), mut.get(key)
        if left is None or right is None or atom_property(left,key)!=atom_property(right,key) or np.linalg.norm(left[0]-right[0]) > .001:
            if left is not None: changed.append(left[0])
            if right is not None: changed.append(right[0])
    if not changed: return np.zeros(len(reference['area']), bool)
    separation, _ = cKDTree(np.array(changed)).query(reference['points'])
    return separation <= distance


def direct_site_masks(pair,native,records,reference):
    wt,mut,audit=aligned_atoms(pair,native,records)
    result=[]
    for residue in audit['mutated_author_residues']:
        coords=[]
        for key in wt.keys()|mut.keys():
            if key[:3]!=tuple(residue) or key[3] in {'N','CA','C','O','OXT'}: continue
            left,right=wt.get(key),mut.get(key)
            if left is None or right is None or atom_property(left,key)!=atom_property(right,key) or np.linalg.norm(left[0]-right[0])>.001:
                if left is not None: coords.append(left[0])
                if right is not None: coords.append(right[0])
        if coords: result.append(cKDTree(np.asarray(coords)).query(reference['points'])[0]<=4.)
    return result


def supported_run_values(values,support):
    """Require perturbation and direct support in the same repeated model.

    Zero is a classification sentinel below every nonnegative threshold; raw
    fields remain untouched and their missing-evidence masks are exported.
    """
    values=np.asarray(values,float);support=np.asarray(support,bool)
    if values.shape[:2]!=support.shape: raise ValueError('direct support/run shape mismatch')
    return np.where(support[...,None] if values.ndim==3 else support,values,0.)


def cached_numerical_bounds(pair,native,records,reference,tools,run,folder,*,required_mask=None):
    """Reuse only the exact numerical sensitivity calculation for this pair."""
    contract={'protocol':protocol_identity([pair],reference,tools),'wt':pair['wt_sha256'],
              'mut':pair['mut_sha256'],'kind':'matched_pair_numerical_bounds_1'}
    if 'query_scope_sha256' in run:contract['query_scope_sha256']=str(run['query_scope_sha256'])
    numerical_run=run
    if required_mask is not None:
        required_mask=np.asarray(required_mask,bool)
        if required_mask.shape!=(len(reference['area']),):raise ValueError('numerical scope/reference mismatch')
        contract['numerical_query_policy']='required_direct_support_only_1'
        contract['numerical_mask_sha256']=fingerprint(required_mask.tolist())
        numerical_run={**run,'wt_sdf':np.where(required_mask,run['wt_sdf'],np.nan)}
    path=Path(folder)/('numerical_'+fingerprint(contract)+'.npz');manifest=path.with_suffix('.json')
    from .legacy_fields import verified_legacy_id
    legacy_id=verified_legacy_id(getattr(tools,'config',{}))
    if legacy_id and not path.exists():
        legacy_contract={**contract,'protocol':protocol_identity([pair],reference,tools,implementation_id=legacy_id)}
        legacy_path=Path(folder)/('numerical_'+fingerprint(legacy_contract)+'.npz')
        if legacy_path.exists() and legacy_path.with_suffix('.json').exists():
            metadata=json.loads(legacy_path.with_suffix('.json').read_text(encoding='utf-8'))
            if metadata.get('contract')!=legacy_contract or sha256_file(legacy_path)!=metadata.get('sha256'):
                raise ValueError('legacy numerical sensitivity fingerprint mismatch')
            import shutil
            shutil.copy2(legacy_path,path)
            write_json(manifest,{'contract':contract,'sha256':sha256_file(path),
                'raw_value_migration':{'source':str(legacy_path),'legacy_implementation_id':legacy_id,
                    'scope':'unchanged numerical_bounds AST and physical source; thresholds recomputed'}})
    if path.is_file() and manifest.is_file():
        metadata=json.loads(manifest.read_text(encoding='utf-8'))
        if metadata.get('contract')!=contract or metadata.get('sha256')!=sha256_file(path):
            raise ValueError('numerical sensitivity cache fingerprint mismatch')
        with np.load(path,allow_pickle=False) as data:
            values={name:data[name].copy() for name in data.files}
        if set(values)!=set(FIELD_CHANNELS) or any(v.shape!=(len(reference['area']),) for v in values.values()):
            raise ValueError('numerical sensitivity cache reference/channel mismatch')
        return values
    values=numerical_bounds(pair,native,records,reference,tools,numerical_run)
    path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,**values)
    write_json(manifest,{'contract':contract,'sha256':sha256_file(path)})
    return values


def seed_raw_field_cache(source,destination,reference,tools):
    """Copy content-verified raw fields across calibration/production phases.

    Classification and coverage are always recomputed. The field loader still
    checks its complete expected contract; no approximate matching is allowed.
    """
    import shutil
    source,destination=Path(source),Path(destination)
    if not source.is_dir():return []
    from .legacy_fields import verified_legacy_id
    legacy_id=verified_legacy_id(getattr(tools,'config',{}))
    destination.mkdir(parents=True,exist_ok=True);reference_hash=reference_id(reference)
    nano_hash=tools.identity(tools.nano);copied=[]
    for manifest in source.glob('field_*.json'):
        metadata=json.loads(manifest.read_text(encoding='utf-8'));contract=metadata.get('contract',{})
        accepted_ids={IMPLEMENTATION_ID}|({legacy_id} if legacy_id else set())
        if contract.get('reference')!=reference_hash or contract.get('implementation_sha256') not in accepted_ids or contract.get('nano')!=nano_hash:continue
        if manifest.name!='field_'+fingerprint(contract)+'.json':raise ValueError('raw field cache filename identity mismatch')
        array=manifest.with_suffix('.npz')
        if not array.is_file() or sha256_file(array)!=metadata.get('sha256'):
            raise ValueError('raw field cache source fingerprint mismatch')
        if contract.get('implementation_sha256')!=IMPLEMENTATION_ID:
            new_contract={**contract,'implementation_sha256':IMPLEMENTATION_ID}
            name='field_'+fingerprint(new_contract)
            target=destination/(name+'.npz');target_manifest=destination/(name+'.json')
            if target.exists() or target_manifest.exists():continue
            shutil.copy2(array,target)
            write_json(target_manifest,{**metadata,'contract':new_contract,
                'raw_value_migration':{'source':str(array),'legacy_implementation_id':legacy_id,
                    'scope':'physical numerical source byte-identical; old calibration acceptance not reused'}})
            copied.append({'field':target.name,'sha256':metadata['sha256'],'source':str(array),'raw_value_migration':True})
        else:
            if (destination/array.name).exists() or (destination/manifest.name).exists():continue
            shutil.copy2(array,destination/array.name);shutil.copy2(manifest,destination/manifest.name)
            copied.append({'field':array.name,'sha256':metadata['sha256'],'source':str(array)})
    return copied


def run_physical_preview(*, structure, records, options, candidates, paired_records,
                         notation, supplement, recheck, output, settings=None, register_pairs=None):
    if (settings or {}).get('optimization_profile') == 'coverage_scan_v3_preview':
        from .engine_v3 import run
        return run(structure=structure,records=records,options=options,candidates=candidates,
                   paired_records=paired_records,notation=notation,supplement=supplement,recheck=recheck,
                   output=output,settings=settings,register_pairs=register_pairs)
    settings = settings or {}; output = Path(output)
    from .combo_contract import enabled,waiting_confirmation,site_evidence,goal_state
    combo=enabled(settings)
    folder = output/'surface_perturbation'; folder.mkdir(exist_ok=True)
    phase = settings.get('phase', 'production')
    limit=settings.get('calibration_limit_seconds',14400) if phase=='calibration' else min(7200,settings.get('limit_seconds',7200))
    budget = ACTIVE_BUDGET.get() or RunBudget(folder/f'{phase}_budget.json', limit)
    tools = PhysicalTools(settings.get('tool_cache', folder/'tool_cache'), settings, budget)
    from .deadlines import guard_mesh
    from .query_cache import memoize_mesh
    original_ses=tools.ses
    tools.ses=lambda *args,**kwargs:memoize_mesh(guard_mesh(original_ses(*args,**kwargs),budget),tools) if hasattr(tools,'cache') else guard_mesh(original_ses(*args,**kwargs),budget)
    native = atoms(Path(structure)); calibration_dir = Path(settings.get('calibration_dir', folder/'calibration'))
    trajectories, errors, final, all_pairs = [], [], [], {}
    summary = {'strategy': STRATEGY_ID, 'phase': phase, 'status': 'INCOMPLETE',
               'panel_contract':settings.get('panel_contract'),
               'goal_status':'NOT_EVALUATED','search_status':'NOT_EVALUATED',
               'interpretation': 'computational preview; not functional antibody validation',
               'initial_budget': int(np.clip(math.ceil(int(settings.get('ecd_length',len(records)))/30),4,12)),
               'method_contract':{'version':'paired_WT_native_SES_preview_1.2',
                 'query_scope_contract':'complete_actual_candidate_support_1',
                 'classification_sha256':fingerprint({name:sha256_file(Path(__file__).parent/name)
                    for name in ('engine.py','classification.py','final_evidence.py','directions.py','coverage.py','electrostatics.py',
                                 'calibration_scope.py','convergence.py','coverage_upper.py','scoped_fields.py')}),
                 'difference':'SDF(Mut)-SDF(matchedWT) in native WT coordinates',
                 'ecd_length':int(settings.get('ecd_length',len(records))),'mapped_residues':len(records),
                 'direct_atom_support_a':4.,'local_query_a':7.,'shell_a':[-2.,4.],
                 'scan_center_geodesic_radius_a':12.,'core_alignment_guard_a':.5,
                 'directions':'joint vector null magnitude; same four repeats above threshold with pairwise positive cosine',
                 'direct_repeat_evidence':'direct atom support and threshold exceedance in the same four of five models',
                 'direction_minimum_cosine':float(settings.get('direction_minimum_cosine',0.)),
                 'chemical_bandwidth_a':1.,'triple_extra_area_fraction':float(settings.get('triple_extra_area_fraction',.05)),
                 'apbs_grid_resolution':{'version':'grid_box_resolution_preview_1',
                     'absolute_tolerance_kt_e':float(settings.get('apbs_grid_absolute_tolerance_kt_e',.05)),
                     'relative_tolerance':float(settings.get('apbs_grid_relative_tolerance',.1))},
                 'coverage':'verified area-weighted WT footprint-center opportunity; excludes unresolved evidence',
                 'runtime':'cumulative charged wall time; interrupted unclosed intervals conservatively bounded'}}
    try:
        budget.remaining()
        reliable = {r.uniprot_pos for r in records if r.coordinate_eligible and not r.review_reasons}
        # The denominator freezes upstream strict replacements, not selected candidates.
        from ..stability import AUTO_EXCLUDED_MUTANTS
        def policy_allowed(option):
            return option.mutant_aa not in AUTO_EXCLUDED_MUTANTS and option.residue.one_letter not in 'CGP'
        designable = {o.residue.uniprot_pos for o in options
                      if o.eligible and not o.review_required and policy_allowed(o)}
        summary['method_contract']['mutation_policy'] = {
            'version': 'shared_auto_exclusions_1',
            'excluded_mutants': sorted(AUTO_EXCLUDED_MUTANTS),
            'protected_native_residues': list('CGP'),
            'qualification': 'also requires existing residue protection and formal stability gates'}
        native_mesh=tools.ses(native)
        reference = canonical_reference(native_mesh,native,records,reliable,designable,
                                        spacing=float(settings.get('reference_spacing_a',2.)))
        reference['context_identity']=fingerprint({'native_input':sha256_file(Path(structure)),
            'source':settings.get('reference_input',{}),
            'mapping':[(r.chain_id,r.auth_seq_id,r.insertion_code,r.uniprot_pos,
                        r.coordinate_eligible,r.review_reasons,r.domain) for r in records]})
        summary['reference_area_all_a2']=float(reference['area'].sum())
        summary['reference_area_designable_a2']=float(reference['area'][reference['designable']].sum())
        save_npz(folder/'reference_graph.npz',reference['graph'])
        np.savez_compressed(folder/'native_WT_SES_mesh.npz',vertices=native_mesh.vertices,faces=native_mesh.faces,
                            triangle_area=native_mesh.area,triangle_normals=native_mesh.normals)
        _,reference['component']=connected_components(reference['graph'],directed=False)
        np.savez_compressed(folder/'canonical_reference.npz',**{k:v for k,v in reference.items() if isinstance(v,np.ndarray)})
        write_json(folder/'reference_contract.json',{'reference_id':reference_id(reference),
                   'input_sha256':sha256_file(Path(structure)),
                   'context_identity':reference['context_identity'],
                   'source_confidence_identity':settings.get('reference_input',{}),
                   'samples':len(reference['area']),'area_all_a2':float(reference['area'].sum()),
                   'area_designable_a2':float(reference['area'][reference['designable']].sum()),
                   'reference_spacing_a':settings.get('reference_spacing_a',2.),
                   'basis':'native SES, independent of SAS',
                   'all_denominator':'reliable external ECD, including protected and unmutable areas',
                   'designable_positions':sorted(designable),
                   'excluded_area_a2':reference['excluded_area_a2'],
                   'internal_cavity_area_a2':reference['internal_cavity_area_a2']})
        summary['raw_field_cache_seeded']=seed_raw_field_cache(calibration_dir.parent/'fields',folder/'fields',reference,tools)
        # A new, explicitly requested validation batch may use historical raw
        # physical values, but never historical calibration or acceptance.
        for source in settings.get('raw_field_cache_roots',[]):
            summary['raw_field_cache_seeded']+=seed_raw_field_cache(source,folder/'fields',reference,tools)
        summary['raw_field_cache_seeded']+=seed_raw_field_cache(folder/'fields',folder/'fields',reference,tools)
        rows = {tuple(sorted(row['mutations'])):row for row in paired_records}
        from .resume_pairs import restore_pairs,capture_manifests
        resume_manifests=capture_manifests(folder,calibration_dir.parent)
        resume_manifests += [Path(p)/'surface_perturbation/paired_manifest.json' for p in settings.get('paired_manifest_cache_roots',[])]
        opposite=output.parent/('production' if phase=='calibration' else 'calibration')/'surface_perturbation'
        resume_manifests+=capture_manifests(opposite,opposite)
        pool = [c for c in candidates if c.order < 4 and (c.formal_eligible or (combo and waiting_confirmation(c)))
                and all(policy_allowed(o) for o in c.mutations)]
        if combo:
            pool=[c for c in pool if c.order in (2,3)]
            pending=[c for c in pool if not c.formal_eligible]
            # A finite new-confirmation batch; later attempts restore exact
            # supplemented models before advancing the pending queue.
            cached=set()
            for manifest_path in resume_manifests:
                if manifest_path.is_file():cached.update(json.loads(manifest_path.read_text(encoding='utf-8')))
            cached_pending=[c for c in pending if c.candidate_id in cached]
            unseen=[c for c in pending if c.candidate_id not in cached]
            batch=unseen[:int(settings.get('confirmation_batch',8))]
            summary['model_confirmation_pending_candidate_ids']=[c.candidate_id for c in unseen[len(batch):]]
            pool=[c for c in pool if c.formal_eligible]+cached_pending+batch
        fields = {}; proxy = {}; hard = {}; unknown = {}; actual = {}; ablations={}; q99_masks={}; two_masks={}; buffers={}; local_support={}
        # Retain earlier repeat order and null-model sources before overwrite.
        def ensure(candidate, count):
            budget.remaining(); cid=candidate.candidate_id; signature=tuple(sorted(notation(candidate)))
            if cid not in all_pairs:
                initial=list(rows.get(signature,{}).get('paired_models',[]))
                all_pairs[cid]=restore_pairs(cid,signature,initial,resume_manifests)
                summary['restored_model_pairs']=summary.get('restored_model_pairs',0)+max(0,len(all_pairs[cid])-len(initial))
            pairs=all_pairs[cid]
            # Cold production starts have not populated all_pairs in the outer
            # loop yet. Decide restoration only after restoring legal identities.
            if (count==1 and phase!='calibration' and len(pairs)>=5
                    and settings.get('optimization_profile')=='coverage_scan_v2_preview'):
                count=5
            if register_pairs is not None:register_pairs(pairs)
            # Deduplicate run identity, never deduplicate coordinate-identical repeats.
            identities={(p['batch_path'],p['group'],p['run']) for p in pairs}
            if len(pairs)<count:
                for p in supplement(candidate,count-len(pairs)):
                    key=(p['batch_path'],p['group'],p['run'])
                    if key not in identities: pairs.append(p);identities.add(key)
            if len(pairs)<count: raise ValueError('paired repeats incomplete')
            for pair in pairs: validate_pair(pair)
            runs=fields.setdefault(cid,[])
            # Reuse all legal model evidence, but calculate only the repeats
            # requested at this stage. A cache with 20 calibration models must
            # not force 20 expensive fields during the one-run shortlist.
            for pair in pairs[len(runs):count]:
                budget.remaining();started=time.perf_counter()
                value,reused=repeated_field(pair,buffers.setdefault(cid,{}),
                    lambda:compute_pair(pair,native,records,reference,tools,folder/'fields'))
                runs.append(value)
                if reused:summary['exact_buffer_reuses']=summary.get('exact_buffer_reuses',0)+1
                if hasattr(tools,'activity'):tools.activity.append({'stage':'paired_surface_field',
                    'candidate_id':cid,'run':pair['run'],'runtime_seconds':time.perf_counter()-started})
            if not recheck(candidate,pairs[:min(5,len(runs))]):
                candidate.formal_eligible=False
                hard[cid]=np.zeros(len(reference['area']),bool)
                return
            support=direct_support(pairs[0],native,records,reference)
            run_support=np.array([direct_support(p,native,records,reference) for p in pairs[:len(runs)]])
            local_support[cid]=np.any(run_support,axis=0)
            # Evaluate every actual requested support cell in every repeat.
            # Source caches remain immutable; a separate scope-bound cache
            # stores real supplements, including their physical unknowns.
            originals=[repeated_field(p,buffers.setdefault(cid,{}),
                lambda p=p:compute_pair(p,native,records,reference,tools,folder/'fields'))[0]
                for p in pairs[:len(runs)]]
            runs=[complete_scope(p,raw,native,records,reference,tools,local_support[cid],folder/'fields')
                  for p,raw in zip(pairs,originals)]
            fields[cid]=runs
            # Uncalibrated fields arrange computation ONLY, never verified coverage.
            provisional=support & np.isfinite(runs[0]['shape_sdf_a']) & (np.abs(runs[0]['shape_sdf_a'])>0)
            proxy[cid]=scan_support(reference['graph'],provisional,reference['owner'])
            protocol=protocol_identity(pairs,reference,tools)
            calibration=load_calibration(calibration_dir,protocol)
            from .classification import classify_fields
            classified=classify_fields(runs,pairs,native,records,reference,calibration,protocol,
                minimum_cosine=float(settings.get('direction_minimum_cosine',0.)),support=support,run_support=run_support)
            known=classified['known'];repeat_support=classified['repeat_support'];changed=classified['direct']
            run_support=classified['run_support']
            np.savez_compressed(folder/f'{cid}_direction_evidence.npz',**classified['vector_evidence'])
            np.savez_compressed(folder/f'{cid}_channel_resolution.npz',**classified['channel_known'])
            np.savez_compressed(folder/f'{cid}_direct_site_evidence.npz',site_changes=classified['site_changes'])
            actual[cid]=changed;unknown[cid]=~known;hard[cid]=classified['hard']
            if combo:
                summary.setdefault('candidate_site_evidence',{})[cid]=site_evidence(candidate,classified['site_changes'],reference['area'])
            q99_masks[cid]=classified['q99_hard'];two_masks[cid]=classified['two_direct_sites']
            ablations[cid]=classified['groups']
            np.savez_compressed(folder/f'{cid}_coverage.npz',direct=changed,hard=hard[cid],unknown=~known,
                                direct_support=repeat_support,direct_support_runs=run_support,priority_proxy=proxy[cid],
                                two_direct_sites=two_masks[cid])
            write_json(folder/f'{cid}_field_audit.json',{'candidate_id':cid,'runs':len(runs),
                       'calibration_protocol':protocol,'calibrated_channels':list(calibration),
                       'query_scope_sha256':str(runs[0]['query_scope_sha256']),
                       'supplemented_query_points_by_run':[int(r['supplemented_query_points']) for r in runs],
                       'known_mask_semantics':'at least one resolvable hard channel; does not imply every channel is resolved',
                       'contact_site_rule':'direct support and threshold exceedance must share the same four repeats for each residue',
                       'channel_resolved_local_points':{name:int(np.sum(mask&support)) for name,mask in classified['channel_known'].items()},
                       'valid_local_points':int(np.sum(known&support)),
                       'actual_changed_area_a2':float(reference['area'][changed].sum()),
                       'unknown_area_a2':float(reference['area'][~known].sum()),
                       'status':('VERIFIED_LOCAL_CHANNELS' if support.any() and known[support].all() else
                                 'PARTIALLY_VERIFIED_LOCAL_CHANNELS' if np.any(known&support) else 'UNVERIFIED') if len(runs)>=5 else 'UNVERIFIED',
                       'units':{name:'angstrom_squared' if name.endswith('_a2') else 'angstrom' if name.endswith('_a') else 'dimensionless_descriptor' for name in FIELD_CHANNELS}})

        for candidate in pool:
            try:
                # Restore already measured legal repeats before classifying.
                # A one-repeat provisional zero/unknown must not hide existing
                # five-repeat evidence or trigger redundant model production.
                ensure(candidate,1) # ensure restores the legal cache before choosing repeat count.
            except (RuntimeError,ValueError) as exc:
                candidate.formal_eligible=False
                candidate.formal_exclusion_reasons=list(getattr(candidate,'formal_exclusion_reasons',[]))+['PAIRED_PHYSICAL_MODEL_QC_FAILED']
                errors.append({'candidate':candidate.candidate_id,'reason':str(exc),'qualification':'physical model evidence rejected, not zero perturbation'})
        measured=[c for c in pool if c.candidate_id in proxy and c.formal_eligible]
        from .benchmark import controls
        summary['baseline_controls']=controls(measured,records,reference)
        if phase=='calibration':
            # Separate job/cost; matched per-candidate environment avoids pooling
            # incompatible reference neighborhoods into a permissive threshold.
            calibration_pool=measured
            expansion_file=output.parent/'production/surface_perturbation/blindspot_calibration_requests.json'
            expanding=bool(settings.get('calibrate_expansion_requests') and expansion_file.is_file())
            if expanding:
                requested_ids=set(json.loads(expansion_file.read_text(encoding='utf-8')).get('candidate_ids',[]))
                calibration_pool=[c for c in measured if c.candidate_id in requested_ids]
                summary['calibration_requested_candidate_ids']=sorted(requested_ids)
            buckets=[calibration_pool];appendable_ids={c.candidate_id for c in calibration_pool};deferred_complete=True
            if expanding or combo:
                from .calibration_schedule import calibration_buckets
                evidence={c.candidate_id:load_calibration(calibration_dir,protocol_identity(all_pairs[c.candidate_id],reference,tools)) for c in calibration_pool}
                failed={c.candidate_id for c in calibration_pool if
                    (calibration_dir/'task_failures'/(protocol_identity(all_pairs[c.candidate_id],reference,tools)+'.json')).is_file()}
                fresh,retry,deferred=calibration_buckets(calibration_pool,evidence,local_support,reference['area'],failed)
                buckets=[fresh,retry];appendable_ids={c.candidate_id for c in fresh+retry}
                from .calibration_scope import local_scope
                deferred_complete=all(local_scope(evidence[cid],local_support[cid],reference['area'])['complete'] for cid in deferred)
                summary['calibration_deferred_candidate_ids']=deferred
                summary['calibration_deferral_contract']='matched complete local null or same-scope 20-model cap; unresolved channels stay unknown and continue blocking global convergence'
            selected=[];capacity=int(settings.get('calibration_candidates',12))
            for bucket in buckets:
                count=min(len(bucket),max(0,capacity-len(selected)))
                priority=nested_panel(np.array([proxy[c.candidate_id] for c in bucket],bool).reshape(len(bucket),len(reference['area'])),
                                      reference['area'],[c.structural_risk for c in bucket],reference['graph'],
                                      [c.candidate_id for c in bucket],count)
                priority += [i for i in range(len(bucket)) if i not in priority][:count-len(priority)]
                selected.extend(bucket[i] for i in priority)
            count=len(selected)
            summary['calibration_selected_candidate_ids']=[c.candidate_id for c in selected]
            summary['calibration_remaining_appendable_candidate_ids']=sorted(appendable_ids-{c.candidate_id for c in selected})
            for candidate in selected:
                cid=candidate.candidate_id;protocol=protocol_identity(all_pairs[cid],reference,tools)
                repeat_count=0;numerical=None;numerical_scope=None;runs=None;calibration=None
                with isolated_job(calibration_dir,protocol,cid,trajectories,errors,lambda:repeat_count):
                    try:
                        for repeat_count in settings.get('calibration_repeat_schedule') or (5,10,20):
                            extension_ids=settings.get('calibration_extension_candidate_ids')
                            if (repeat_count>5 and settings.get('optimization_profile')=='coverage_scan_v2_preview'
                                    and extension_ids is not None and cid not in extension_ids):
                                summary.setdefault('calibration_extension_deferred_candidate_ids',[]).append(cid)
                                break
                            ensure(candidate,repeat_count)
                            pairs=all_pairs[cid];runs=fields[cid]
                            if numerical is None or numerical_scope!=str(runs[0]['query_scope_sha256']):
                                numerical=cached_numerical_bounds(pairs[0],native,records,reference,tools,runs[0],folder/'fields',
                                                                 required_mask=local_support[cid])
                                numerical_scope=str(runs[0]['query_scope_sha256'])
                            calibration=save_calibration(calibration_dir,pairs,runs,reference,tools,numerical,
                                                         required_mask=local_support[cid])
                            from .calibration_scope import local_scope
                            if local_scope(calibration,local_support[cid],reference['area'])['complete']:break
                            if settings.get('optimization_profile') == 'coverage_scan_v2_preview':
                                from .unknowns import calibration_unknowns
                                unknown_audit=calibration_unknowns(calibration,local_support[cid],runs=runs)
                                summary.setdefault('calibration_unknowns',{})[cid]=unknown_audit
                                if not unknown_audit['appendable']:break
                        clear_failure(calibration_dir,protocol)
                        ensure(candidate,len(runs))
                        trajectories.append({'candidate':cid,'calibration_runs':len(runs),
                            'local_scope_converged':local_scope(calibration,local_support[cid],reference['area'])['complete'],
                            'local_scope':local_scope(calibration,local_support[cid],reference['area']),
                            'statuses':{k:v.get('status') for k,v in calibration.items()}})
                    finally:
                        fields[cid]=fields[cid][:1];buffers.pop(cid,None)
                        runs=None;calibration=None;numerical=None
            selected_complete=bool(count) and all(r['local_scope_converged'] for r in trajectories)
            scope_complete=selected_complete
            if expanding or combo:
                scope_complete=bool(calibration_pool) and deferred_complete and not summary['calibration_remaining_appendable_candidate_ids'] and (not count or selected_complete)
                summary['calibration_expansion_scope_converged']=scope_complete
            summary['status']='LOCAL_CALIBRATION_COMPLETE' if scope_complete else 'CALIBRATION_INSUFFICIENT'
        else:
            previous=None
            upper_bounds={}
            def actual_upper(candidate):
                cid=candidate.candidate_id
                if cid not in upper_bounds and len(all_pairs[cid])>=5:
                    from .coverage_upper import hard_scan_upper_bound
                    masks=np.array([direct_support(p,native,records,reference) for p in all_pairs[cid][:5]])
                    upper_bounds[cid]=hard_scan_upper_bound(masks,reference['graph'])
                return upper_bounds.get(cid)
            for iteration in range(int(settings.get('maximum_iterations',8))):
                budget.remaining()
                measured=[c for c in measured if c.formal_eligible]
                if not measured: break
                priorities=np.array([hard[c.candidate_id] if len(fields[c.candidate_id])>=5 else proxy[c.candidate_id] for c in measured])
                order=nested_panel(priorities,reference['area'],[c.structural_risk for c in measured],
                                   reference['graph'],[c.candidate_id for c in measured],12)
                proposed=[measured[i] for i in order]
                identities=[c.candidate_id for c in proposed]
                for candidate in proposed:
                    try: ensure(candidate,3);ensure(candidate,5)
                    except (RuntimeError,ValueError) as exc:
                        candidate.formal_eligible=False;errors.append({'candidate':candidate.candidate_id,'reason':str(exc)})
                from .calibration_scope import local_scope
                scope_reports={c.candidate_id:local_scope(
                    load_calibration(calibration_dir,protocol_identity(all_pairs[c.candidate_id],reference,tools)),
                    local_support[c.candidate_id],reference['area']) for c in measured}
                for candidate in measured:
                    cid=candidate.candidate_id;bound=actual_upper(candidate)
                    measured_five=len(fields[cid])>=5
                    closed=bool(measured_five and bound is not None and not np.any(bound&~hard[cid]))
                    scope_reports[cid].update(hard_scan_bound_resolved=closed,
                        measurement_local_complete=bool(measured_five and not np.any(unknown[cid]&local_support[cid])),
                        hard_scan_unknown_upper_area_a2=float(reference['area'][bound&~hard[cid]].sum()) if measured_five and bound is not None else None)
                summary['candidate_local_calibration']=scope_reports
                verified=[c for c in measured if c.formal_eligible and len(fields[c.candidate_id])>=5
                          and hard[c.candidate_id].any() and (
                              scope_reports[c.candidate_id]['hard_scan_bound_resolved'] or
                              (scope_reports[c.candidate_id]['complete'] and scope_reports[c.candidate_id]['measurement_local_complete']))]
                if combo:
                    verified=[c for c in verified if summary.get('candidate_site_evidence',{}).get(c.candidate_id,{}).get('complete')]
                if settings.get('optimization_profile')=='coverage_scan_v2_preview':
                    remote_controls=[c for c in verified if c.order==2
                        and getattr(c,'search_role','')=='EXACT_TRIPLE_PAIR_CONTROL'
                        and c.risk_components.get('geodesic_span_a',0)>24]
                    summary['remote_double_control_only_candidate_ids']=[c.candidate_id for c in remote_controls]
                    verified=[c for c in verified if c not in remote_controls]
                if settings.get('optimization_profile') == 'coverage_scan_v2_preview':
                    shadow=[c for c in measured if c.formal_eligible and len(fields[c.candidate_id])>=5
                        and summary.get('candidate_site_evidence',{}).get(c.candidate_id,{}).get('complete')]
                    shadow_matrix=np.array([hard[c.candidate_id] for c in shadow],bool).reshape(len(shadow),len(reference['area']))
                    shadow_order=nested_panel(shadow_matrix,reference['area'],[c.structural_risk for c in shadow],reference['graph'],[c.candidate_id for c in shadow],12)
                    summary['shadow_positive_lower_bound']={
                        'status':'SHADOW_ONLY_NOT_FORMAL_QUALIFICATION',
                        'candidate_ids':[shadow[i].candidate_id for i in shadow_order],
                        'coverage_curve':budget_curve(shadow_matrix,np.full(shadow_matrix.shape,np.nan),shadow_order,
                            reference['area'],reference['designable'],reference['graph']),
                        'interpretation':'Same measured fields and original channel rules; full-scope qualification remains mandatory for official results'}
                    np.savez_compressed(folder/'shadow_positive_matrix.npz',hard=shadow_matrix,
                        candidate_ids=np.array([c.candidate_id for c in shadow]),area=reference['area'])
                summary['triple_counterpart_pending_candidate_ids']=[]
                if settings.get('optimization_profile') == 'coverage_scan_v2_preview':
                    doubles=[c for c in verified if c.order==2]
                    double_matrix=np.array([hard[c.candidate_id] for c in doubles],bool).reshape(len(doubles),len(reference['area']))
                    double_order=nested_panel(double_matrix,reference['area'],[c.structural_risk for c in doubles],reference['graph'],[c.candidate_id for c in doubles],12)
                    double_baseline_ids=[doubles[i].candidate_id for i in double_order]
                    summary['triple_expansion_audits']=[]
                for candidate in list(verified):
                    if candidate.order!=3: continue
                    mutations={o.notation for o in candidate.mutations}
                    lower=[c for c in verified if c.order<3 and {o.notation for o in c.mutations}<mutations]
                    if combo:
                        from .combo_contract import actual_pair_controls
                        resolved_controls=[c for c in measured if c.formal_eligible and len(fields[c.candidate_id])>=5
                            and (scope_reports[c.candidate_id]['hard_scan_bound_resolved'] or
                                (scope_reports[c.candidate_id]['complete'] and scope_reports[c.candidate_id]['measurement_local_complete']))]
                        lower=actual_pair_controls(candidate,resolved_controls)
                        if lower is None:
                            verified.remove(candidate)
                            summary.setdefault('triple_counterpart_pending_candidate_ids',[]).append(candidate.candidate_id)
                            continue
                    if getattr(candidate,'triple_role',None) == 'TRIPLE_COMPRESSION':
                        # All own/site/stability/scope and three exact controls passed above.
                        # Panel replacement is audited against a frozen double baseline below.
                        candidate.warnings=[w for w in candidate.warnings if w!='REAL_HIGHER_ORDER_GAIN_PENDING']
                        continue
                    advantage=min((float(reference['area'][hard[candidate.candidate_id]&~hard[c.candidate_id]].sum()) for c in lower),default=0.)
                    if advantage<=0 or advantage<float(settings.get('triple_extra_area_fraction',.05))*reference['area'].sum():
                        verified.remove(candidate)
                        candidate.formal_eligible=False
                        candidate.formal_exclusion_reasons.append('REAL_TRIPLE_GAIN_NOT_DEMONSTRATED')
                    else:
                        if settings.get('optimization_profile') == 'coverage_scan_v2_preview':
                            from .compression import expansion_audit
                            audit=expansion_audit(candidate.candidate_id,double_baseline_ids,hard,reference['area'])
                            summary['triple_expansion_audits'].append(audit)
                            if not audit['accepted']:
                                verified.remove(candidate)
                                candidate.formal_eligible=False
                                candidate.formal_exclusion_reasons.extend(audit['reasons'])
                                continue
                        candidate.warnings=[w for w in candidate.warnings if w!='REAL_HIGHER_ORDER_GAIN_PENDING']
                        if 'REAL_TRIPLE_EXTRA_GAIN_VERIFIED' not in candidate.warnings:
                            candidate.warnings.append('REAL_TRIPLE_EXTRA_GAIN_VERIFIED')
                matrix=np.array([hard[c.candidate_id] for c in verified],bool).reshape(len(verified),-1) if verified else np.zeros((0,len(reference['area'])),bool)
                regular=[i for i,c in enumerate(verified) if getattr(c,'triple_role',None)!='TRIPLE_COMPRESSION']
                verified_order=[regular[i] for i in nested_panel(matrix[regular],reference['area'],[verified[i].structural_risk for i in regular],reference['graph'],[verified[i].candidate_id for i in regular],12)]
                final=[verified[i] for i in verified_order]
                if settings.get('compression_baseline_ids'):
                    from .compression import compression_audit
                    baseline_ids=settings['compression_baseline_ids']
                    baseline_mask=np.any([hard[cid] for cid in baseline_ids],axis=0)
                    panel_ids=list(baseline_ids)
                    by_id={c.candidate_id:c for c in verified}
                    summary['compression_audits']=[]
                    if not set(baseline_ids)<=by_id.keys():
                        raise ValueError('COMPRESSION_BASELINE_NOT_CURRENTLY_VERIFIED')
                    for triple in sorted((c for c in verified if getattr(c,'triple_role',None)=='TRIPLE_COMPRESSION'),key=lambda c:c.candidate_id):
                        audit=compression_audit(panel_ids,triple.parent_candidate_ids,triple.candidate_id,
                            hard,reference['area'],reference['graph'],baseline_mask)
                        summary['compression_audits'].append(audit)
                        if audit['accepted']:
                            panel_ids=audit['after_ids']
                            if 'REAL_TRIPLE_COMPRESSION_VERIFIED' not in triple.warnings:
                                triple.warnings.append('REAL_TRIPLE_COMPRESSION_VERIFIED')
                    final=[by_id[cid] for cid in panel_ids]
                    verified_order=[verified.index(c) for c in final]
                else:
                    final=[c for c in final if getattr(c,'triple_role',None)!='TRIPLE_COMPRESSION']
                current=[c.candidate_id for c in final]
                from ..models import to_dict
                from .local_constructs import local_construct,export_local_constructs
                summary['verified_local_constructs']=[local_construct(
                    to_dict(c),all_pairs[c.candidate_id],scope_reports[c.candidate_id],
                    reference['area'],hard[c.candidate_id],actual[c.candidate_id]) for c in final]
                if combo:
                    for item in summary['verified_local_constructs']:
                        item['site_evidence']=summary['candidate_site_evidence'][item['candidate_id']]
                trajectories.append({'iteration':iteration,'proposed':identities,'verified':current,
                                     'verification_removed':sorted(set(identities)-set(current))})
                # Expand the real verification pool at blind patches before
                # declaring panel identity convergence. A low proxy score never
                # permanently excludes a candidate with matched calibration.
                from .expansion import blindspot_priority
                covered=np.any(matrix[verified_order],axis=0) if verified_order else np.zeros(len(reference['area']),bool)
                pending=[];uncalibrated=[];excluded_by_upper=[]
                for candidate in measured:
                    cid=candidate.candidate_id
                    if not candidate.formal_eligible:continue
                    bound=actual_upper(candidate)
                    if bound is not None and not np.any(bound&~covered):
                        excluded_by_upper.append(cid)
                        continue
                    if not scope_reports[cid]['complete'] or (len(fields[cid])>=5 and not scope_reports[cid]['measurement_local_complete'] and not scope_reports[cid]['hard_scan_bound_resolved']):uncalibrated.append(candidate)
                    elif len(fields[cid])<5:pending.append(candidate)
                summary['hard_blindspot_upper_bound_excluded_candidate_ids']=excluded_by_upper
                summary['upper_bound_scope']='five actual paired models; four-of-five direct atomic support and the identical 12 A canonical graph scan rule; excludes additional hard blindspot area only, not risk/identifiability or K-optimality'
                np.savez_compressed(folder/'hard_scan_upper_bounds.npz',
                    **{cid:bound for cid,bound in upper_bounds.items()})
                ordered,priority_diagnostics=blindspot_priority(reference,covered,pending,proxy)
                requested,uncalibrated_diagnostics=blindspot_priority(reference,covered,uncalibrated,proxy)
                write_json(folder/'blindspot_calibration_requests.json',{
                    'scope':'existing globally deduplicated qualified candidate pool; no unbounded new chemistry enumeration',
                    'candidate_ids':[c.candidate_id for c in requested],
                    'priorities':uncalibrated_diagnostics,'reason':'matched null missing; cannot enter verified panel'})
                batch=ordered[:max(1,int(settings.get('blindspot_verification_batch',2)))]
                if batch:
                    expansion=[]
                    for candidate in batch:
                        cid=candidate.candidate_id
                        try:
                            ensure(candidate,3);ensure(candidate,5)
                            expansion.append({'candidate_id':cid,**priority_diagnostics[cid],
                                'real_field_gain_a2':float(reference['area'][hard[cid]&~covered].sum()),
                                'qualification':'pending next iteration stability/triple gain/panel recheck'})
                        except (RuntimeError,ValueError) as exc:
                            candidate.formal_eligible=False;errors.append({'candidate':cid,'blindspot_verification':str(exc)})
                    trajectories[-1]['blindspot_verification_expansion']=expansion
                    previous=current
                    continue
                from .convergence import convergence_state
                blockers=[c.candidate_id for c in requested]
                if combo:
                    blockers+=summary.get('model_confirmation_pending_candidate_ids',[])+summary.get('triple_counterpart_pending_candidate_ids',[])
                state=convergence_state(current,previous,blockers)
                if state=='CALIBRATION_EXPANSION_REQUIRED':
                    summary['status']='CALIBRATION_EXPANSION_REQUIRED'
                    summary['calibration_expansion_candidate_ids']=[c.candidate_id for c in requested]
                    summary['verified_within_calibrated_pool_candidate_ids']=current
                    summary['convergence_scope']='candidate pool incomplete; proxy is not an exclusion bound'
                    export_local_constructs(output,summary['verified_local_constructs'],summary['status'])
                    break
                if state=='VERIFIED_PANEL':
                    export_local_constructs(output,summary['verified_local_constructs'],state)
                    summary['convergence_scope']='hard blindspot search closed by validated fields or five-model upper bounds; no global K, risk or identifiability optimum claim'
                    summary['status']='VERIFIED_PANEL'; break
                previous=current
                if not proposed: break
            summary['search_status']=summary['status']
            if summary['status']!='VERIFIED_PANEL' and not combo: final=[]
            if final:
                summary['apbs_audit_target_candidate_ids']=[c.candidate_id for c in final]
                summary['apbs_base_solve_ceiling']=2*5*len(final)
                draft_matrix=np.array([hard[c.candidate_id] for c in final],bool)
                np.savez_compressed(folder/'verified_hard_draft_matrix.npz',hard=draft_matrix,
                    candidate_ids=np.array([c.candidate_id for c in final]),area=reference['area'],designable=reference['designable'])
                write_json(folder/'verified_hard_draft_panel.json',{
                    'status':'HARD_EVIDENCE_COMPLETE_APBS_AUDIT_PENDING_NOT_FINAL',
                    'candidate_ids':[c.candidate_id for c in final],
                    'coverage_curve':budget_curve(draft_matrix,np.full(draft_matrix.shape,np.nan),
                        list(range(len(final))),reference['area'],reference['designable'],reference['graph']),
                    'interpretation':'real verified hard evidence checkpoint; a timed-out case remains incomplete'})
                summary['baseline_controls']=controls(final,records,reference)
                matrix=np.array([hard[c.candidate_id] for c in final]); extended=np.full(matrix.shape,np.nan)
                completed_apbs=[];apbs_reports={}
                hard_curve=budget_curve(matrix,extended,list(range(len(final))),reference['area'],reference['designable'],reference['graph'])
                from .final_evidence import apbs_panel_ready
                audit_ready=apbs_panel_ready(combo=combo,construct_count=len(final),
                    hard_coverage=hard_curve[-1]['C_hard_all'],
                    final_panel_audit=settings.get('final_panel_audit',False))
                summary['final_panel_audit_requested']=bool(settings.get('final_panel_audit',False))
                if settings.get('audit_apbs',True) and audit_ready:
                    from .electrostatics import audit_candidate,prepare_candidate
                    from .tools import make_grid
                    prepared_panel={};panel_xyz=[]
                    for candidate in final:
                        try:
                            data=prepare_candidate(candidate.candidate_id,all_pairs[candidate.candidate_id][:5],native,records,tools,folder)
                            prepared_panel[candidate.candidate_id]=data;panel_xyz.extend(data[2])
                        except (RuntimeError,ValueError) as exc: errors.append({'candidate':candidate.candidate_id,'apbs_preparation':str(exc)})
                    shared_grid=make_grid(panel_xyz,.5) if panel_xyz else None
                    summary['shared_panel_APBS_grid']=shared_grid
                    charge={'K':1,'R':1,'D':-1,'E':-1}
                    sensitive=max(range(len(final)),key=lambda i:sum(abs(charge.get(o.mutant_aa,0)-charge.get(o.residue.one_letter,0)) for o in final[i].mutations))
                    summary['grid_check_candidate_id']=final[sensitive].candidate_id
                    summary['grid_check_selection']='largest formal substitution charge-change proxy; not hidden contact labels'
                    for i,candidate in enumerate(final):
                        if candidate.candidate_id not in prepared_panel: continue
                        try:
                            changes,missing,audit=audit_candidate(candidate.candidate_id,all_pairs[candidate.candidate_id][:5],native,records,reference,tools,folder,convergence=i==sensitive,
                                prepared_data=prepared_panel[candidate.candidate_id],shared_grid=shared_grid)
                            completed_apbs.append(candidate.candidate_id)
                            apbs_reports[candidate.candidate_id]=audit
                            extended[i]=changes.astype(float);extended[i,missing]=np.nan
                        except (RuntimeError,ValueError) as exc: errors.append({'candidate':candidate.candidate_id,'apbs':str(exc)})
                from .final_evidence import require_apbs_audits
                summary['apbs_completed_candidate_ids']=completed_apbs
                summary['apbs_incomplete_candidate_ids']=sorted(set(c.candidate_id for c in final)-set(completed_apbs))
                if settings.get('final_panel_audit',False) and audit_ready and set(completed_apbs)!=set(c.candidate_id for c in final):
                    summary['final_apbs_evidence']={'status':'APBS_EXECUTION_INCOMPLETE',
                        'execution_status':'INCOMPLETE','preparation_status':'UNRESOLVED',
                        'numerical_resolution_status':'NOT_EVALUATED',
                        'completed_candidate_ids':completed_apbs,
                        'missing_candidate_ids':summary['apbs_incomplete_candidate_ids'],
                        'interpretation':'Final audit failed; retained hard evidence does not imply APBS acceptance'}
                else:
                    summary['final_apbs_evidence']=(require_apbs_audits([c.candidate_id for c in final],completed_apbs,
                    requested=settings.get('audit_apbs',True),reports=apbs_reports) if audit_ready else
                    {'status':'NOT_REQUESTED_TARGET_NOT_MET','execution_status':'NOT_REQUESTED',
                     'preparation_status':'NOT_REQUESTED','numerical_resolution_status':'NOT_EVALUATED',
                     'interpretation':'hard draft only; no APBS on intermediate candidate panels'})
                summary['coverage_curve']=budget_curve(matrix,extended,list(range(len(final))),reference['area'],reference['designable'],reference['graph'])
                verified_matrix=np.array([hard[c.candidate_id] for c in verified],bool)
                summary['independent_K_verified_candidate_ids']=[c.candidate_id for c in verified]
                summary['independent_K_hard_bounds']=[independent_bounds(verified_matrix,reference['area'],k,seconds=min(5.,budget.remaining())) for k in range(1,13)]
                np.savez_compressed(folder/'independent_K_verified_pool.npz',hard=verified_matrix,
                                    candidate_ids=np.array([c.candidate_id for c in verified]),area=reference['area'])
                summary['physical_ablations']={label:budget_curve(np.array([ablations[c.candidate_id][label] for c in final]),
                    np.full(matrix.shape,np.nan),list(range(len(final))),reference['area'],reference['designable'],reference['graph'])
                    for label in ('atom_occupancy_shape','shape_chemistry','shape_chemistry_ses')}
                summary['sensitivity_curves']={label:budget_curve(np.array([masks[c.candidate_id] for c in final]),np.full(matrix.shape,np.nan),
                    list(range(len(final))),reference['area'],reference['designable'],reference['graph'])
                    for label,masks in (('empirical_Q99',q99_masks),('two_direct_change_residues',two_masks))}
                np.savez_compressed(folder/'coverage_matrix.npz',hard=matrix,extended=extended,
                                    candidate_ids=np.array([c.candidate_id for c in final]),area=reference['area'],designable=reference['designable'])
                summary['nested_verified_candidate_ids']=[c.candidate_id for c in final]
                record_by_position={r.uniprot_pos:r for r in records}
                within_partition=[i for i,c in enumerate(verified) if len({record_by_position[o.residue.uniprot_pos].region_id for o in c.mutations})==1]
                partition_matrix=verified_matrix[within_partition]
                partition_order=nested_panel(partition_matrix,reference['area'],[verified[i].structural_risk for i in within_partition],
                    reference['graph'],[verified[i].candidate_id for i in within_partition],12)
                summary['route_controls']={'hard_partition_same_measured_pool':budget_curve(partition_matrix,
                    np.full(partition_matrix.shape,np.nan),partition_order,reference['area'],reference['designable'],reference['graph']),
                    'overlapping_windows_global_optimization':summary['coverage_curve'],
                    'scope':'measured candidate filter comparison; not exhaustive independent searches',
                    'unpartitioned_new_candidate_search':'NOT_MEASURED_IN_THIS_POOL'}
            else:
                if summary['status']!='CALIBRATION_EXPANSION_REQUIRED':
                    summary['status']='NO_VERIFIED_PANEL'
                empty=np.zeros((0,len(reference['area'])),bool)
                summary['coverage_curve']=budget_curve(empty,empty.astype(float),[],reference['area'],reference['designable'],reference['graph'])
    except (TimeoutError,RuntimeError,ValueError) as exc:
        summary['status']='INCOMPLETE';errors.append({'reason':str(exc)})
        final=[]
    finally:
        budget.save(ACTIVE_BUDGET.get() is not None)
        nested_count=len(final)
        requested=settings.get('primary_budget',{})
        primary_count=min(len(final),int(requested.get('maximum_constructs',12)))
        if combo and phase=='production':
            from .combo_contract import primary_count as combo_count
            primary_count=combo_count(summary.get('coverage_curve',[]),len(final),summary['initial_budget'],requested.get('maximum_constructs',12)) if final else 0
        elif requested.get('adaptive_budget',True): primary_count=min(primary_count,summary['initial_budget'])
        primary=final[:primary_count]
        if settings.get('optimization_profile')=='coverage_scan_v2_preview' and phase=='production':
            summary.setdefault('final_apbs_evidence',{
                'status':'NOT_REQUESTED_EMPTY_OFFICIAL_PANEL' if not primary else 'NOT_REQUESTED_TARGET_NOT_MET',
                'execution_status':'NOT_REQUESTED','preparation_status':'NOT_REQUESTED',
                'numerical_resolution_status':'NOT_EVALUATED',
                'interpretation':'No APBS acceptance is inferred before a qualifying target panel exists'})
        if combo and phase=='production':
            coverage=summary.get('coverage_curve',[])
            hard_all=coverage[primary_count-1]['C_hard_all'] if primary_count and coverage else None
            summary['goal_status']=goal_state(primary,hard_all,summary.get('candidate_site_evidence',{}),
                audit_complete=summary.get('final_apbs_evidence',{}).get('status')=='APBS_EXECUTION_COMPLETE',fatal=summary['status']=='INCOMPLETE')
            if summary['goal_status']=='ACHIEVED' and summary['status']!='VERIFIED_PANEL':
                summary['status']='TARGET_PANEL_VERIFIED_SEARCH_OPEN'
        summary['selected_candidate_ids']=[c.candidate_id for c in primary]
        summary.update(errors=errors,trajectories=trajectories,tool_calls=tools.calls,
                       tool_activity=getattr(tools,'activity',[]),
                       elapsed_seconds=budget.elapsed,actual_constructs=primary_count,
                       nested_verified_constructs=nested_count,
                       missing_constructs_to_twelve=12-nested_count,
                       apbs_base_solve_ceiling=summary.get('apbs_base_solve_ceiling',2*5*nested_count))
        # Local evidence survives incomplete global closure, but not a fatal
        # validation error. Always replace the current attempt's exports.
        if summary['status']=='INCOMPLETE':
            summary['verified_local_constructs']=[]
        from .local_constructs import export_local_constructs
        export_local_constructs(output,summary.get('verified_local_constructs',[]),summary['status'])
        write_json(folder/'paired_manifest.json',all_pairs)
        write_json(output/'surface_perturbation_preview.json',summary)
        from .export import evidence_archive
        if settings.get('export_mode', 'deliverable') == 'deliverable':
            evidence_archive(output,all_pairs,calibration_dir)
    return primary,summary
