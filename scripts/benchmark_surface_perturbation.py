"""Frozen five-case replay. Calibration is a separately costed phase.

Run from repository parent with the project venv. Outputs never overwrite
historical runs; production/retry time is cumulative within each output case.
"""
from __future__ import annotations
import argparse
import json
import sys
import math
import time
import shutil
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.pipeline import PipelineRequest, run_pipeline
from ecd_epitope_scanner.retrieval import parse_uniprot_entry
from ecd_epitope_scanner.util import write_json, sha256_file
from ecd_epitope_scanner.surface_perturbation.case_status import case_status
from ecd_epitope_scanner.scripts.replay_joint_v41_run import replay_request_data

CASES={'CTLA4':'8c068ad701a6','PD1':'fa85119485cd','PDL1':'95bfd808b533','B7-H6':'2763d972f2b7','EGFR':'26655601df08'}
_overlay_path=Path(__file__).resolve().parents[1]/'source_overlay_manifest.json'
LAUNCH_SOURCE_MANIFEST={'python':sys.executable,'argv':sys.argv,
    'overlay':json.loads(_overlay_path.read_text(encoding='utf-8')) if _overlay_path.is_file() else None,
    'interpretation':'source snapshot captured at benchmark process startup; shared workers separately audited'}


def preview_request(data,ecd_length):
    target=max(4,min(12,math.ceil(ecd_length/30)))
    return {**data,'design_engine':'joint_v4_ml','surface_strategy':'surface_perturbation_v1_preview',
            'max_primary_constructs':None,'minimum_constructs':4,'target_constructs':target,
            'maximum_constructs':12,'adaptive_budget':True}


def frozen_worker_parameters(root,settings):
    """Resolve installed dependencies against the saved-input workspace.

    A detached source snapshot intentionally excludes vendor executables.
    Fail before upstream computation if the chosen dependency is unavailable.
    """
    upstream=dict(settings.get('upstream_parameters',{}))
    foldx=dict(upstream.get('foldx_confirmation',{}))
    if foldx.get('binary'):
        binary=Path(foldx['binary'])
        if not binary.is_absolute():binary=root/binary
    else:
        installed=sorted((root/'vendor/foldx').glob('foldx_*.exe'))
        if not installed:raise ValueError('FROZEN_DEPENDENCY_MISSING: FoldX in source workspace')
        binary=installed[-1]
    if not binary.is_file():raise ValueError('FROZEN_DEPENDENCY_MISSING: '+str(binary))
    digest=sha256_file(binary)
    expected=settings.get('foldx_binary_sha256')
    if expected and expected!=digest:raise ValueError('FROZEN_DEPENDENCY_CHANGED: FoldX binary')
    foldx['binary']=str(binary.resolve())
    foldx.setdefault('engine_work_root',str((root/'workspace_data/foldx_short').resolve()))
    upstream['foldx_confirmation']=foldx
    return upstream,{'foldx_binary':str(binary.resolve()),'foldx_binary_sha256':digest}


def run_case(name,root,destination,phase,settings):
    upstream,dependencies=frozen_worker_parameters(root,settings)
    source=root/'workspace_data'/CASES[name]
    analysis=json.loads((source/'analysis.json').read_text(encoding='utf-8'))
    entry=json.loads((source/'uniprot.json').read_text(encoding='utf-8'))
    request=json.loads((source/'request.json').read_text(encoding='utf-8'))
    output=destination/name/phase;output.mkdir(parents=True,exist_ok=True)
    snapshot=output/'attempt_snapshots'/str(time.time_ns())
    for filename in ('launch_source_manifest.json','frozen_input_manifest.json','surface_perturbation_preview.json','phase_result.json','case_acceptance.json',
                     'surface_perturbation_local_constructs.json','surface_perturbation_local_constructs.csv',
                     'surface_perturbation_evidence.zip','ecd_epitope_design.zip','run_manifest.json',
                     'surface_perturbation/reference_contract.json','surface_perturbation/canonical_reference.npz',
                     'surface_perturbation/reference_graph.npz'):
        if (output/filename).is_file():
            (snapshot/filename).parent.mkdir(parents=True,exist_ok=True);shutil.copy2(output/filename,snapshot/filename)
    write_json(output/'launch_source_manifest.json',LAUNCH_SOURCE_MANIFEST)
    data=replay_request_data(source,output,analysis,request)
    # Historical B7-H6 requests selected the retired V3 engine. Keep their
    # inputs, while explicitly testing the newly requested V4 preview method.
    data=preview_request(data,len(analysis['protein']['ecd_sequence']))
    budget_root=Path(settings.get('cumulative_budget_root',destination))
    calibration_root=budget_root/name/'calibration'
    field_roots=settings.get('raw_field_cache_roots',[])
    if isinstance(field_roots,dict):field_roots=field_roots.get(name,[])
    manifest_roots=settings.get('paired_manifest_cache_roots',[])
    if isinstance(manifest_roots,dict):manifest_roots=manifest_roots.get(name,[])
    null_roots=settings.get('null_cache_roots',[])
    if isinstance(null_roots,dict):null_roots=null_roots.get(name,[])
    imported=[]
    for null_root in null_roots:
        cache=Path(null_root)
        for path in cache.rglob('*'):
            if not path.is_file() or path.suffix not in {'.json','.npz'}:continue
            target=calibration_root/'surface_perturbation/calibration'/path.relative_to(cache)
            if target.exists():continue
            target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
            imported.append({'source':str(path),'target':str(target),'sha256':sha256_file(target)})
    if imported:write_json(output/'raw_null_import_manifest.json',{'files':imported,'scope':'raw matched evidence only; current engine recomputes calibration qualification'})
    data['metadata']={**data.get('metadata',{}),'joint_v4_ml':{
        **upstream,
        'runtime':{'cache_dir':str(destination/name/'upstream_cache')},
        'surface_perturbation':{**settings,'phase':phase,'calibration_dir':str(calibration_root/'surface_perturbation'/'calibration'),
            'raw_field_cache_roots':field_roots,
            'paired_manifest_cache_roots':manifest_roots,
            'cumulative_budget_path':str(budget_root/name/phase/f'surface_{phase}_cumulative_budget.json'),
            'paired_cache_roots':settings.get('paired_cache_roots',[])+[str(output)]+([str(calibration_root)] if phase=='production' else []),
            'tool_cache':str(Path(settings.get('reuse_tool_cache_root',destination))/name/'physical_cache')}}}
    frozen={r['uniprot_pos']:r['conservation_score'] for r in analysis['residues'] if r.get('conservation_score') is not None}
    write_json(output/'frozen_input_manifest.json',{'case':name,'source':str(source),
        'structure_type':request['input_type'],'replay_input_type':data['input_type'],
        'source_request_sha256':sha256_file(source/'request.json'),
        'source_design_engine':request.get('design_engine'),'replay_design_engine':data['design_engine'],
        'model_sha256':sha256_file(Path(data['structure_path'])),
        'confidence_sha256':sha256_file(Path(data['confidence_path'])) if data.get('confidence_path') else None,
        'uniprot_sha256':sha256_file(source/'uniprot.json'),'phase':phase,
        'worker_dependencies':dependencies})
    def protein(accession,client,output_dir,ecd_override=None):
        if accession!=analysis['protein']['accession']: raise ValueError('frozen accession mismatch')
        write_json(output_dir/'uniprot.json',entry);return parse_uniprot_entry(entry,ecd_override)
    with patch('ecd_epitope_scanner.pipeline.fetch_protein_context',protein),patch('ecd_epitope_scanner.conservation.fetch_conservation',return_value=frozen or None):
        return run_pipeline(PipelineRequest(**data),lambda stage,value,message:print(f'{name} {phase} {stage}: {message}',flush=True))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=list(CASES),action='append')
    parser.add_argument('--phase',choices=['calibration','production','both'],required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--settings',type=Path)
    parser.add_argument('--candidate-route',choices=['hard_partition','unpartitioned_local','overlapping_windows'])
    parser.add_argument('--budget-root',type=Path,help='shared case ledger; required for alternate physical route comparisons')
    parser.add_argument('--source-root',type=Path)
    args=parser.parse_args();root=args.source_root.resolve() if args.source_root else Path(__file__).resolve().parents[1]
    settings=json.loads(args.settings.read_text(encoding='utf-8')) if args.settings else {}
    if args.candidate_route in {'hard_partition','unpartitioned_local'} and not args.budget_root:
        parser.error('alternate routes require --budget-root; changing output folders must not reset the case budget')
    if args.budget_root:settings['cumulative_budget_root']=str(args.budget_root.resolve())
    if args.candidate_route:settings['candidate_route']=args.candidate_route
    results={}
    for name in args.case or CASES:
        results[name]={}
        for phase in ('calibration','production') if args.phase=='both' else (args.phase,):
            try:
                started=time.time()
                run_case(name,root,args.output.resolve(),phase,settings)
                summary_path=args.output/name/phase/'surface_perturbation_preview.json'
                if not summary_path.is_file() or summary_path.stat().st_mtime < started:
                    raise RuntimeError('NO_CURRENT_ATTEMPT_SUMMARY: preserved prior outputs are not current results')
                summary=json.loads(summary_path.read_text(encoding='utf-8'))
                status=case_status(args.output/name/phase,summary)
                results[name][phase]={'status':status,'physical_status':summary['status'],
                                      'actual_constructs':summary.get('actual_constructs'),
                                      'goal_status':summary.get('goal_status'),'search_status':summary.get('search_status'),
                                      'elapsed_seconds':summary.get('elapsed_seconds'),'summary':str(summary_path)}
            except Exception as exc:
                status='INCOMPLETE_RUNTIME_EXHAUSTED' if isinstance(exc,TimeoutError) and str(exc)=='CUMULATIVE_RUNTIME_EXHAUSTED' else 'INCOMPLETE'
                results[name][phase]={'status':status,'reason':str(exc)}
                print(name,repr(exc),flush=True)
            write_json(args.output/name/phase/'phase_result.json',results[name][phase])
        # Disjoint case lanes must not race on the same temporary report file.
        report_name=f'{args.phase}_{name}_cases.json' if args.case and len(args.case)==1 else f'{args.phase}_cases.json'
        write_json(args.output/report_name,results)
if __name__=='__main__': main()
