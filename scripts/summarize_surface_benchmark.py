"""Read current five-case state without treating preserved outputs as fresh."""
from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from ecd_epitope_scanner.scripts.benchmark_surface_perturbation import CASES
from ecd_epitope_scanner.util import write_json
from ecd_epitope_scanner.surface_perturbation.case_status import case_status


def read(path):
    try:return json.loads(path.read_text(encoding='utf-8'))
    except (OSError,ValueError):return None


def summarize(root):
    report={}
    for name in CASES:
        report[name]={}
        for phase in ('calibration','production'):
            folder=root/name/phase
            budget=read(folder/f'surface_{phase}_cumulative_budget.json') or {}
            summary_path=folder/'surface_perturbation_preview.json'
            summary=read(summary_path) or {}
            manifest=folder/'frozen_input_manifest.json'
            current=bool(summary and manifest.is_file() and summary_path.stat().st_mtime>=manifest.stat().st_mtime)
            active=budget.get('active',False)
            acceptance=read(folder/'case_acceptance.json') or {}
            elapsed=budget.get('elapsed_seconds')
            if active and elapsed is not None: elapsed+=max(0.,time.time()-budget['checkpoint_wall_time'])
            status='ACTIVE_OR_INTERRUPTED' if active else case_status(folder,summary) if current else 'NO_CURRENT_SUMMARY'
            accepted=current and status=='VERIFIED_PANEL'
            report[name][phase]={'status':status,'physical_status':summary.get('status') if current else None,
                'summary_is_current':current,'elapsed_seconds':elapsed,
                'field_files':len(list((folder/'surface_perturbation'/'fields').glob('*.json'))),
                'null_contexts':len(list((folder/'surface_perturbation'/'calibration').glob('*.json'))),
                'actual_constructs':summary.get('actual_constructs') if accepted else None,
                'nested_verified_constructs':summary.get('nested_verified_constructs') if accepted else None,
                'checkpoint_nested_constructs':summary.get('nested_verified_constructs') if current and not active else None,
                'tool_calls':summary.get('tool_calls') if current and not active else None}
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    p.add_argument('--output',type=Path);args=p.parse_args();result=summarize(args.root)
    if args.output:write_json(args.output,result)
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
