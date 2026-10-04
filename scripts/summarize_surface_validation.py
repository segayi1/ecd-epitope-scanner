"""Read-only current-attempt acceptance, never reuse historical coverage as current."""
import argparse
import json
import time
from pathlib import Path
from ..surface_perturbation.case_status import case_status
from ..util import sha256_file, write_json


def read(path):
    return json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {}


def phase_report(folder):
    folder=Path(folder)
    summary_path=folder/'surface_perturbation_preview.json'
    launch=folder/'launch_source_manifest.json'
    summary=read(summary_path)
    current=summary_path.is_file() and (not launch.is_file() or summary_path.stat().st_mtime>=launch.stat().st_mtime)
    phase=folder.name
    ledger=read(folder/f'surface_{phase}_cumulative_budget.json')
    elapsed=ledger.get('elapsed_seconds')
    if ledger.get('active'):
        elapsed+=max(0.,time.time()-ledger['checkpoint_wall_time'])
        status='RUNNING_NOT_ACCEPTED'
    elif elapsed is not None and elapsed>=ledger['limit_seconds']:
        status='INCOMPLETE_RUNTIME_EXHAUSTED'
    elif not current:
        status='INCOMPLETE_NO_CURRENT_SUMMARY' if launch.is_file() else 'NOT_STARTED'
    else:
        status=case_status(folder,summary)
    driver_path=folder/'phase_result.json'
    driver=read(driver_path)
    if not ledger.get('active') and driver_path.is_file() and (not launch.is_file() or driver_path.stat().st_mtime>=launch.stat().st_mtime):
        # A generic exception must not erase the more specific budget failure.
        if status!='INCOMPLETE_RUNTIME_EXHAUSTED' and driver.get('status','').startswith('INCOMPLETE'):
            status=driver['status']
    accepted=phase=='production' and current and status=='VERIFIED_PANEL'
    result={'status':status,'physical_status':summary.get('status') if current else None,
        'summary_is_current':current,'accepted_panel':accepted,'elapsed_seconds':elapsed,
        'limit_seconds':ledger.get('limit_seconds'),'driver_reason':driver.get('reason'),
        'summary_sha256':sha256_file(summary_path) if summary_path.is_file() else None,
        'historical_summary_retained':bool(summary and not current),
        'actual_constructs':summary.get('actual_constructs') if accepted else None,
        'calibration_expansion_candidate_count':len(summary.get('calibration_expansion_candidate_ids',[])) if current else None,
        'apbs_evidence':summary.get('final_apbs_evidence') if current else None,
        'coverage':None}
    if accepted:
        curve=summary.get('coverage_curve',[])
        primary=next((row for row in curve if row.get('actual_constructs')==summary.get('actual_constructs')),None)
        result['coverage']={key:primary.get(key) for key in
            ('C_hard_all','C_hard_designable','C_extended_all','C_extended_designable')} if primary else None
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result={'interpretation':'Current source/attempt, cumulative budget and policy determine acceptance; no coverage for unfinished or historical panels.',
        'cases':{case:{phase:phase_report(args.root/case/phase) for phase in ('calibration','production')}
                 for case in ('CTLA4','PD1','PDL1','B7-H6','EGFR')}}
    write_json(args.output,result)
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':main()
