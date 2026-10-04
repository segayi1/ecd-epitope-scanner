"""Stop an owned coordinator only at a verified closed production boundary.

Scientific workers are never terminated. Suspending the coordinator briefly
closes the spawn race while checking children and the cumulative ledger lock.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import psutil

from ecd_epitope_scanner.util import sha256_file, write_json


def closed_production(case, earliest):
    reports=list((case/'stages').glob('*_report.json'))
    if not reports or not earliest.is_file():return None
    latest=max(reports,key=lambda path:path.stat().st_mtime_ns)
    if latest.stat().st_mtime_ns<earliest.stat().st_mtime_ns:return None
    report=json.loads(latest.read_text(encoding='utf-8'))
    if report.get('engineering_status')!='COMPLETED':return None
    if report.get('configuration',{}).get('phase')!='production':return None
    return latest


def stop_at_boundary(case, pid, earliest, timeout):
    process=psutil.Process(pid)
    identity=process.create_time()
    command=process.cmdline()
    owned=('--resume-batches' in command or 'ecd_epitope_scanner.scripts.run_ctla4_coverage_v3' in command)
    case_arg=Path(command[command.index('--case-root')+1]).resolve() if '--case-root' in command else None
    if not owned or case_arg!=case:
        raise ValueError('PID is not the explicitly owned case coordinator')
    deadline=time.monotonic()+timeout
    while time.monotonic()<deadline:
        if not process.is_running():
            return {'status':'COORDINATOR_ALREADY_TERMINAL','pid':pid}
        if closed_production(case,earliest) is None:
            time.sleep(.2);continue
        suspended=False
        try:
            if process.create_time()!=identity:raise ValueError('Coordinator PID identity changed')
            process.suspend();suspended=True
            if any(child.is_running() for child in process.children(recursive=True)):
                continue
            report_path=closed_production(case,earliest)
            if report_path is None:continue
            budget_path=case/'case_budget.json'
            with budget_path.with_suffix('.lock').open('r+b') as lock:
                try:
                    if os.name=='nt':
                        import msvcrt
                        msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
                    else:
                        import fcntl
                        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                except OSError:
                    continue
                budget=json.loads(budget_path.read_text(encoding='utf-8'))
                if budget.get('active'):continue
                digest=sha256_file(budget_path)
                process.kill();suspended=False
                process.wait(timeout=10)
                return {'status':'PASS','pid':pid,'coordinator_create_time':identity,
                    'reason':'ENGINEERING_RESTART_AT_VERIFIED_CLOSED_PRODUCTION_BOUNDARY',
                    'production_report':str(report_path),'production_report_sha256':sha256_file(report_path),
                    'ledger_sha256':digest,'ledger':budget,'scientific_workers_terminated':0}
        finally:
            if suspended:
                process.resume()
                # Also yield on the early-continue paths above, so a missed
                # boundary cannot repeatedly suspend an active supervisor.
                time.sleep(.2)
        time.sleep(.2)
    return {'status':'NO_BOUNDARY_OBSERVED','pid':pid,'scientific_workers_terminated':0}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-root',type=Path,required=True)
    parser.add_argument('--pid',type=int,required=True)
    parser.add_argument('--earliest-production-report',type=Path,required=True)
    parser.add_argument('--timeout-seconds',type=float,default=3600)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    case=args.case_root.resolve()
    earliest=args.earliest_production_report.resolve()
    if earliest.parent!=(case/'stages').resolve():raise ValueError('Report must belong to this case')
    report=stop_at_boundary(case,args.pid,earliest,args.timeout_seconds)
    write_json(args.report or case/'coordinator_boundary_restart.json',report)
    print(json.dumps(report,ensure_ascii=False))


if __name__=='__main__':main()
