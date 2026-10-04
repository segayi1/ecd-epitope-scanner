"""Bounded engineering replay of saved pairs; never a production budget reset."""
import argparse,json,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.sparse import load_npz
from ..util import write_json,sha256_file
from ..surface_perturbation.paired import atoms,validate_pair
from ..surface_perturbation.fields import reference_id,IMPLEMENTATION_ID
from ..surface_perturbation.scoped_fields import complete_scope,query_present
from ..surface_perturbation.tools import PhysicalTools
from ..surface_perturbation.runtime import RunBudget
from ..surface_perturbation.deadlines import guard_mesh
from ..surface_perturbation.query_cache import memoize_mesh
from ..surface_perturbation.power_guard import prevent_sleep
from ..surface_perturbation.legacy_fields import verified_legacy_id


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('case',type=Path)
    parser.add_argument('--diagnosis',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    cohorts=json.loads(args.diagnosis.read_text(encoding='utf-8'))[args.case.name]['cohorts']
    cohort=min(cohorts,key=lambda r:r['all_WT_queries_present_local_area_fraction'])
    phase=args.case/'calibration';folder=phase/'surface_perturbation'
    records=[SimpleNamespace(**r) for r in json.loads((phase/'analysis.json').read_text(encoding='utf-8'))['residues']]
    native=atoms(phase/'structure_for_energy.pdb')
    with np.load(folder/'canonical_reference.npz',allow_pickle=False) as data:reference={k:data[k].copy() for k in data.files}
    reference['graph']=load_npz(folder/'reference_graph.npz')
    reference['context_identity']=json.loads((folder/'reference_contract.json').read_text(encoding='utf-8'))['context_identity']
    manifest_path=Path(cohort['null_manifest']).resolve()
    meta=json.loads(manifest_path.read_text(encoding='utf-8'))
    if reference_id(reference)!=meta['reference_id']:raise ValueError('saved reference mismatch')
    with np.load(manifest_path.with_suffix('.npz'),allow_pickle=False) as data:required=data['shape_sdf_a__required_mask'].copy()
    candidates=json.loads((folder/'paired_manifest.json').read_text(encoding='utf-8'))
    def signature(pair):return (pair['wt_sha256'],pair['mut_sha256'],pair['pdb_list_sha256'])
    matched=[(cid,rows[:20]) for cid,rows in candidates.items()
        if [signature(p) for p in rows[:20]]==[signature(p) for p in meta['pairs']]]
    if len(matched)!=1:raise ValueError('saved null cohort identity is not unique')
    cid,pairs=matched[0]
    budget=RunBudget(args.output/'engineering_budget.json',300)
    # Same existing physical cache; no FoldX, SPURS or APBS is requested.
    original_settings=json.loads((args.case.parents[1]/'frozen_settings_20261002_v2.json').read_text(encoding='utf-8'))
    legacy=verified_legacy_id(original_settings)
    accepted_sources={IMPLEMENTATION_ID}|({legacy} if legacy else set())
    cache=Path(original_settings['reuse_tool_cache_root'])/args.case.name/'physical_cache'
    tools=PhysicalTools(cache,{},budget)
    ses=tools.ses;tools.ses=lambda *a,**kw:memoize_mesh(guard_mesh(ses(*a,**kw),budget),tools)
    report={'scope':'saved-pair query allocation engineering replay only; not calibrated thresholds, a final panel or five-case acceptance',
        'candidate_id':cid,'before_WT_query_fraction':cohort['all_WT_queries_present_local_area_fraction'],
        'original_protocol':meta['protocol_id'],'pairs':[]}
    masks=[];started=time.perf_counter()
    try:
        with prevent_sleep():
            for pair,path_string in zip(pairs,cohort['raw_manifests']):
                budget.remaining();validate_pair(pair)
                path=Path(path_string).resolve();raw_meta=json.loads(path.read_text(encoding='utf-8'));array=path.with_suffix('.npz')
                before=sha256_file(array)
                if before!=raw_meta['sha256'] or raw_meta['contract']['implementation_sha256'] not in accepted_sources:
                    raise ValueError('source raw field identity mismatch')
                with np.load(array,allow_pickle=False) as data:raw={k:data[k].copy() for k in data.files}
                result=complete_scope(pair,raw,native,records,reference,tools,required,args.output/'scoped_fields')
                masks.append(query_present(result))
                if before!=sha256_file(array):raise ValueError('original raw field was modified')
                report['pairs'].append({'wt_sha256':pair['wt_sha256'],'mut_sha256':pair['mut_sha256'],
                    'source_sha256_unchanged':before,'supplemented_points':int(result['supplemented_query_points']),
                    'source_implementation':raw_meta['contract']['implementation_sha256'],
                    'topology_rejection_preserved':bool(result['topology_changed'])==bool(raw['topology_changed'])})
        known=np.all(masks,axis=0);area=reference['area']
        report['after_WT_query_fraction']=float(area[required&known].sum()/area[required].sum())
        report['status']='PASS' if known[required].all() else 'UNRESOLVED_QUERY_ALLOCATION'
    finally:
        report['runtime_seconds']=time.perf_counter()-started;report['tool_calls']=tools.calls
        budget.save(False);write_json(args.output/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='pairs'}))


if __name__=='__main__':main()
