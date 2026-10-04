"""Small, read-only diagnostics of saved evidence and pure calibration logic."""
from __future__ import annotations
import collections
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
CAMPAIGN = ROOT/'workspace_data/surface_combo_four_case_20261003/batch_20261003_155853'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        while chunk:=handle.read(1024*1024):h.update(chunk)
    return h.hexdigest()


def main():
    snapshots=json.loads((OUT/'campaign_snapshot.json').read_text(encoding='utf-8'))
    report={'scope':'saved artifacts and synthetic pure-function diagnostics; no production/solver runs', 'phases':{}}
    for key,value in snapshots.items():
        if '/' not in key:continue
        s=value.get('surface_perturbation_preview.json',{})
        if not s:continue
        timing=collections.defaultdict(lambda:{'calls':0,'seconds':0.0})
        for row in s.get('tool_activity',[]):
            name=row.get('stage') or row.get('tool')
            timing[name]['calls']+=1;timing[name]['seconds']+=row.get('runtime_seconds',0)
        scopes=s.get('candidate_local_calibration',{})
        report['phases'][key]={
            'status':s.get('status'),'actual_constructs':s.get('actual_constructs'),
            'goal_status':s.get('goal_status'),'elapsed_seconds':s.get('elapsed_seconds'),
            'tool_calls':s.get('tool_calls'),'nested_timing_do_not_sum':dict(timing),
            'scopes':len(scopes),'complete_scopes':sum(v['complete'] for v in scopes.values()),
            'union_full_but_incomplete_scopes':sum(v['calibrated_local_area_fraction']==1 and not v['complete'] for v in scopes.values()),
            'site_evidence_complete':sum(v.get('complete',False) for v in s.get('candidate_site_evidence',{}).values()),
            'calibration_trajectories':[{'candidate':v.get('candidate'),'runs':v.get('calibration_runs'),
                'complete':v.get('local_scope_converged')} for v in s.get('trajectories',[]) if 'calibration_runs' in v],
            'budget':value.get('surface_'+key.split('/')[1]+'_cumulative_budget.json')}
    folder=CAMPAIGN/'CTLA4/production'
    pool=json.loads((folder/'panel_candidates.json').read_text(encoding='utf-8'))
    report['ctla4_pool']={'count':len(pool),'orders':dict(collections.Counter(str(x['order']) for x in pool)),
        'formal_eligible':sum(x['formal_eligible'] for x in pool),
        'exclusions_nonexclusive':dict(collections.Counter(y for x in pool for y in x['formal_exclusion_reasons']))}
    hierarchy=json.loads((folder/'selected_panel_v4.json').read_text(encoding='utf-8'))['hierarchy_resolution']
    report['hierarchy']={k:{'hypotheses':len(v['signatures']),'pairs':len(v['pairwise_distances']),
        'near_pairs':len(v['near_pairs'])} for k,v in hierarchy.items() if isinstance(v,dict)}
    report['hierarchy']['identical_compact_payload_bytes']=len(json.dumps(hierarchy,separators=(',',':')))
    other=json.loads((folder/'hypothesis_hierarchy_v411.json').read_text(encoding='utf-8'))
    report['hierarchy']['selected_panel_and_hierarchy_payload_equal']=hierarchy==other['resolution_by_level']
    physical=collections.defaultdict(lambda:[0,0])
    txt_names=collections.defaultdict(lambda:[0,0])
    for base,dirs,files in os.walk(ROOT/'workspace_data'):
        if 'physical_cache' not in Path(base).parts:continue
        for name in files:
            path=Path(base)/name;n=path.stat().st_size
            physical[path.suffix.lower()][0]+=n;physical[path.suffix.lower()][1]+=1
            if path.suffix.lower()=='.txt':txt_names[name][0]+=n;txt_names[name][1]+=1
    report['physical_cache_extensions']=dict(physical)
    report['physical_cache_txt']=sorted(txt_names.items(),key=lambda x:x[1][0],reverse=True)[:12]
    left=CAMPAIGN/'CTLA4/calibration/surface_perturbation/fields'
    right=folder/'surface_perturbation/fields'
    shared=0; shared_bytes=0; changed=[];same_file_count=0
    for a in left.glob('*.npz'):
        b=right/a.name
        if b.is_file() and a.stat().st_size==b.stat().st_size:
            if digest(a)==digest(b):
                shared+=1;shared_bytes+=b.stat().st_size
                same_file_count+=int(a.samefile(b))
            else:changed.append(a.name)
    report['ctla4_cross_phase_arrays']={'sha256_verified_identical_files':shared,'duplicate_bytes':shared_bytes,
        'same_underlying_file_count':same_file_count,'same_size_changed':changed}
    spec=importlib.util.spec_from_file_location('audit_calibration',ROOT/'surface_perturbation/calibration.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    repro=[]
    for n in (5,10,20):
        blocks=np.zeros((n-1,2));blocks[:,1]=np.nan
        result=module.calibrate_channel(blocks,np.zeros(2),'synthetic','channel',required_mask=np.ones(2,bool))
        assert result['valid'].tolist()==[True,False]
        assert not result['requested_scope_converged']
        repro.append({'models':n,'valid_mask':result['valid'].tolist(),'complete':result['requested_scope_converged']})
    report['deterministic_unknown_repeat_repro']=repro
    report['iid_illustration_not_actual_null_statistics']={str(n):1-.95**(n-1) for n in (5,10,20)}
    reference=folder/'surface_perturbation/canonical_reference.npz'
    with np.load(reference,allow_pickle=False) as data:area=data['area'].copy()
    with np.load(folder/'surface_perturbation/verified_hard_draft_matrix.npz',allow_pickle=False) as data:
        report['draft_matrix_fields']=list(data.files)
        hard=data['hard'];names=data['candidate_ids'].tolist()
        report['independent_saved_draft_coverage']={'constructs':len(names),'candidate_ids':names,
            'C_hard_all':float(area[np.any(hard,axis=0)].sum()/area.sum()),'area_a2':float(area.sum())}
    sources=['pipeline.py','models.py','exporting.py','scripts/benchmark_surface_perturbation.py',
        'surface_perturbation/engine.py','surface_perturbation/tools.py','surface_perturbation/geometry.py',
        'surface_perturbation/query_cache.py','surface_perturbation/fields.py','surface_perturbation/calibration.py',
        'surface_perturbation/calibration_scope.py','surface_perturbation/null_store.py','surface_perturbation/export.py',
        'surface_perturbation/combo_contract.py','joint_design_v4/reporting.py','joint_design_v4/identifiability.py',
        'integrations/workers/foldx_worker.py','integrations/workers/spurs_gpu_worker.py']
    report['source_sha256']={s:digest(ROOT/s) for s in sources}
    launch=json.loads((folder/'launch_source_manifest.json').read_text(encoding='utf-8'))
    frozen=launch.get('overlay',{}).get('files',{})
    report['current_vs_executed_source']={s:{'current':h,'executed':frozen.get(s),'equal':h==frozen.get(s)}
                                           for s,h in report['source_sha256'].items()}
    (OUT/'analysis.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['source_sha256','physical_cache_txt','phases']},ensure_ascii=False))


if __name__=='__main__':main()
