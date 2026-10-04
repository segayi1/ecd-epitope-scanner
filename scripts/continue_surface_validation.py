"""Serial frozen phases with fresh process boundaries and cumulative case ledgers."""
from __future__ import annotations
import argparse,json,subprocess,time
from pathlib import Path
from ..util import write_json,sha256_file
from .benchmark_surface_perturbation import CASES
from ..surface_perturbation.case_status import case_status
from .surface_phase_guard import run_guarded
from ..surface_perturbation.runtime import read_budget


def verify_frozen_source(source,checkout):
    """Concurrent lanes share a read-only snapshot, never a live overlay write."""
    manifest=json.loads((checkout/'source_overlay_manifest.json').read_text(encoding='utf-8'))
    if Path(manifest['source']).resolve()!=source.resolve():
        raise ValueError('frozen source belongs to a different workspace')
    for name,digest in manifest['files'].items():
        if sha256_file(source/name)!=digest or sha256_file(checkout/name)!=digest:
            raise ValueError('FROZEN_SOURCE_CHANGED: '+name)


def current_bundle(folder,started):
    return all(path.is_file() and path.stat().st_mtime>=started for path in
               (folder/'run_manifest.json',folder/'ecd_epitope_design.zip'))


def resumed_phases(case_folder,start_phase):
    if start_phase=='calibration':return ['calibration','production']
    folder=Path(case_folder)/'calibration'
    ledger=read_budget(folder/'surface_calibration_cumulative_budget.json')
    launch=folder/'launch_source_manifest.json';summary=folder/'surface_perturbation_preview.json'
    if (not ledger or ledger.get('active') or not launch.is_file() or not summary.is_file()
        or summary.stat().st_mtime<launch.stat().st_mtime
        or not current_bundle(folder,launch.stat().st_mtime)):
        raise ValueError('PRODUCTION_RESUME_REQUIRES_COMPLETED_CURRENT_CALIBRATION_DELIVERY')
    # Partial null evidence is allowed here. The physical engine still rechecks
    # protocol, local scope, real repeats, stability and unresolved candidates.
    return ['production']


def expansion_has_work(identities,previous_sets,calibration_summary):
    # Unchanged requests do not mean all requests were attempted when a batch
    # cap left additional candidates unprocessed. Scientific status stays unresolved.
    requested=set(identities)
    remaining=set(calibration_summary.get('calibration_remaining_appendable_candidate_ids',[]))
    if requested&remaining:return True
    known=set(calibration_summary.get('calibration_requested_candidate_ids',[]))
    known.update(calibration_summary.get('calibration_selected_candidate_ids',[]))
    known.update(calibration_summary.get('calibration_deferred_candidate_ids',[]))
    if 'calibration_remaining_appendable_candidate_ids' in calibration_summary and known:
        return bool(requested-known)
    return identities not in previous_sets


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,required=True);parser.add_argument('--checkout',type=Path,required=True)
    parser.add_argument('--settings',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--case',action='append',choices=list(CASES))
    parser.add_argument('--results-file',type=Path,help='separate coordinator file for disjoint case lanes')
    parser.add_argument('--adopt-calibration-pid',type=int,help='wait for a verified running first-case worker without restarting it')
    parser.add_argument('--retry-incomplete-calibration',action='store_true',help='one repaired-source retry within the same calibration ledger')
    parser.add_argument('--frozen-source',action='store_true',help='verify and reuse an immutable source overlay; safe for disjoint concurrent lanes')
    parser.add_argument('--start-phase',choices=['calibration','production'],default='calibration',
                        help='first-case continuation after a completed calibration delivery; physical gates remain unchanged')
    args=parser.parse_args()
    if args.start_phase=='production' and args.adopt_calibration_pid:
        parser.error('production continuation cannot adopt a live calibration worker')
    source=args.source_root.resolve();checkout=args.checkout.resolve();destination=args.output.resolve()
    python=source/'.venv/Scripts/python.exe';results_path=args.results_file or destination/'continuation_results.json'
    results=json.loads(results_path.read_text(encoding='utf-8')) if results_path.exists() else {}
    labels=source/'workspace_data/surface_perturbation_20261001/hidden_contact_labels'
    label_manifest=json.loads((labels/'manifest.json').read_text(encoding='utf-8'))
    adoption_used=False
    cases=args.case or ['PD1','PDL1','B7-H6','EGFR']
    for case in cases:
        results.setdefault(case,{})
        phases=resumed_phases(destination/case,args.start_phase if case==cases[0] else 'calibration')
        retries=0;expansion_sets=[]
        for phase_index,phase in enumerate(phases):
            adopting=bool(args.adopt_calibration_pid and not adoption_used and case==(args.case or ['PD1'])[0] and phase=='calibration')
            if adopting:
                import psutil
                worker=psutil.Process(args.adopt_calibration_pid);command=worker.cmdline()
                if not any('benchmark_surface_perturbation' in word for word in command) or '--case' not in command or command[command.index('--case')+1]!=case or '--phase' not in command or command[command.index('--phase')+1]!=phase:
                    raise ValueError('adoption requires the matching live benchmark worker')
                started=worker.create_time();print(f'ADOPT {case} {phase}: worker {worker.pid}, preserved live process and ledger',flush=True)
                adoption_used=True
                while True:
                    try:worker.wait(timeout=30);break
                    except psutil.TimeoutExpired:pass
                exit_code=None # The OS does not expose an unrelated process exit code.
            else:
                if args.frozen_source:verify_frozen_source(source,checkout)
                else:
                    subprocess.run([str(python),str(source/'scripts/freeze_surface_validation.py'),str(checkout)],
                        cwd=source,check=True)
                started=time.time();print(f'BEGIN {case} {phase}: preserved cumulative ledger; new frozen process',flush=True)
                configuration=json.loads(args.settings.read_text(encoding='utf-8'))
                exit_code=run_guarded([str(python),'-m','ecd_epitope_scanner.scripts.benchmark_surface_perturbation',
                    '--case',case,'--phase',phase,'--source-root',str(source),'--settings',str(args.settings.resolve()),
                    '--output',str(destination)],cwd=checkout.parent,folder=destination/case/phase,phase=phase,
                    limit=configuration.get('calibration_limit_seconds',14400) if phase=='calibration' else min(7200,configuration.get('limit_seconds',7200)))
            folder=destination/case/phase;summary_path=folder/'surface_perturbation_preview.json'
            current=summary_path.exists() and summary_path.stat().st_mtime>=started
            summary=json.loads(summary_path.read_text(encoding='utf-8')) if current else {}
            status=case_status(folder,summary) if current else 'INCOMPLETE_NO_CURRENT_SUMMARY'
            driver_path=folder/'phase_result.json'
            if driver_path.exists() and driver_path.stat().st_mtime>=started:
                driver=json.loads(driver_path.read_text(encoding='utf-8'))
                if status!='INCOMPLETE_RUNTIME_EXHAUSTED' and driver.get('status','').startswith('INCOMPLETE'):status=driver['status']
            results[case][phase]={'status':status,'physical_status':summary.get('status'),
                'goal_status':summary.get('goal_status'),'search_status':summary.get('search_status'),
                'driver_exit_code':exit_code,'adopted_worker_pid':args.adopt_calibration_pid if adopting else None,
                'phase_wall_seconds':time.time()-started,
                'summary':str(summary_path),'settings_sha256':sha256_file(args.settings.resolve())}
            audit_root=destination/'posthoc_audits'/case/phase;audit_root.mkdir(parents=True,exist_ok=True)
            if current and (folder/'surface_perturbation_evidence.zip').exists():
                subprocess.run([str(python),'-m','ecd_epitope_scanner.scripts.audit_surface_archive',
                    str(folder/'surface_perturbation_evidence.zip'),'--output',str(audit_root/'archive.json')],cwd=source.parent)
            if current and current_bundle(folder,started):
                subprocess.run([str(python),'-m','ecd_epitope_scanner.scripts.audit_run_bundle',
                    str(folder),'--output',str(audit_root/'run_bundle.json')],cwd=source.parent)
            elif current:
                write_json(audit_root/'run_bundle.json',{'status':'UNAVAILABLE_STALE_OR_MISSING_BUNDLE',
                    'reason':'Failed attempt retained historical artifacts; not a current delivery'})
            if current and phase=='production':
                label=labels/f'{case}.json'
                subprocess.run([str(python),'-m','ecd_epitope_scanner.scripts.audit_surface_contacts',str(folder),
                    str(label),'--label-sha256',label_manifest['label_files_sha256'][label.name],
                    '--output',str(audit_root/'contact_sanity.json')],cwd=source.parent)
            write_json(results_path,results)
            print(f'END {case} {phase}: {results[case][phase]["status"]}',flush=True)
            if phase=='production' and status=='CALIBRATION_EXPANSION_REQUIRED':
                identities=tuple(sorted(summary.get('calibration_expansion_candidate_ids',[])))
                budgets=[destination/case/p/f'surface_{p}_cumulative_budget.json' for p in ('calibration','production')]
                remaining=all(path.exists() and (lambda b: not b.get('active') and b['elapsed_seconds']<b['limit_seconds'])(json.loads(path.read_text(encoding='utf-8'))) for path in budgets)
                cal_summary_path=destination/case/'calibration/surface_perturbation_preview.json'
                cal_summary=json.loads(cal_summary_path.read_text(encoding='utf-8')) if cal_summary_path.is_file() else {}
                confirmation_pending=bool(summary.get('model_confirmation_pending_candidate_ids'))
                if remaining and (confirmation_pending or expansion_has_work(identities,expansion_sets,cal_summary)) and len(expansion_sets)<12:
                    expansion_sets.append(identities)
                    phases[phase_index+1:phase_index+1]=['calibration','production']
                    print(f'EXPAND {case}: matched local calibration required; same cumulative budgets',flush=True)
                else:
                    print(f'UNRESOLVED {case}: exhausted budget or unchanged requests after maximum-repeat calibration; no global convergence claim',flush=True)
            if phase=='calibration' and args.retry_incomplete_calibration and retries==0 and status.startswith('INCOMPLETE'):
                ledger=folder/'surface_calibration_cumulative_budget.json'
                budget=json.loads(ledger.read_text(encoding='utf-8')) if ledger.exists() else {}
                if not budget.get('active') and budget.get('elapsed_seconds',0)<budget.get('limit_seconds',0):
                    retries+=1;phases.insert(phase_index+1,'calibration')
                    print(f'RETRY {case} calibration once with repaired source; cumulative time preserved',flush=True)

if __name__=='__main__':main()
