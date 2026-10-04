import json
from types import SimpleNamespace

import pytest

from ecd_epitope_scanner.scripts import restart_coordinator_at_boundary as boundary


def prepared(case):
    stages=case/'stages';stages.mkdir()
    report=stages/'production_report.json'
    report.write_text(json.dumps({'engineering_status':'COMPLETED','configuration':{'phase':'production'}}))
    (case/'case_budget.json').write_text(json.dumps({'active':False,'elapsed_seconds':123.,'limit_seconds':None}))
    (case/'case_budget.lock').write_bytes(b'1')
    return report


class Process:
    def __init__(self,case,children=()):
        self.case=case;self.child_rows=children;self.suspended=0;self.resumed=0;self.killed=0
    def create_time(self):return 17.
    def cmdline(self):return ['python','--resume-batches',str(self.case)]
    def is_running(self):return True
    def children(self,recursive):return self.child_rows
    def suspend(self):self.suspended+=1
    def resume(self):self.resumed+=1
    def kill(self):self.killed+=1
    def wait(self,timeout):return 0


def mock_process(monkeypatch,process):
    clock=[0.]
    def sleep(seconds):clock[0]+=seconds
    monkeypatch.setattr(boundary,'psutil',SimpleNamespace(Process=lambda pid:process))
    monkeypatch.setattr(boundary,'time',SimpleNamespace(monotonic=lambda:clock[0],sleep=sleep))


def test_latest_calibration_report_cannot_establish_production_boundary(tmp_path):
    report=prepared(tmp_path)
    assert boundary.closed_production(tmp_path,report)==report
    newer=tmp_path/'stages/new_calibration_report.json'
    newer.write_text(json.dumps({'engineering_status':'COMPLETED','configuration':{'phase':'calibration'}}))
    import os
    stamp=report.stat().st_mtime_ns+10_000_000
    os.utime(newer,ns=(stamp,stamp))
    assert boundary.closed_production(tmp_path,report) is None


def test_live_scientific_child_prevents_coordinator_stop_and_is_never_killed(tmp_path,monkeypatch):
    report=prepared(tmp_path)
    child=SimpleNamespace(is_running=lambda:True)
    process=Process(tmp_path,[child]);mock_process(monkeypatch,process)
    result=boundary.stop_at_boundary(tmp_path,42,report,.1)
    assert result['status']=='NO_BOUNDARY_OBSERVED'
    assert process.killed==0 and process.suspended==process.resumed==1


def test_active_ledger_prevents_stop_even_after_worker_exit(tmp_path,monkeypatch):
    report=prepared(tmp_path)
    (tmp_path/'case_budget.json').write_text(json.dumps({'active':True}))
    process=Process(tmp_path);mock_process(monkeypatch,process)
    result=boundary.stop_at_boundary(tmp_path,42,report,.1)
    assert result['status']=='NO_BOUNDARY_OBSERVED'
    assert process.killed==0 and process.suspended==process.resumed==1


def test_locked_ledger_prevents_stop_despite_closed_json_state(tmp_path,monkeypatch):
    import os
    report=prepared(tmp_path)
    process=Process(tmp_path);mock_process(monkeypatch,process)
    with (tmp_path/'case_budget.lock').open('r+b') as lock:
        if os.name=='nt':
            import msvcrt
            msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
        else:
            import fcntl
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        result=boundary.stop_at_boundary(tmp_path,42,report,.1)
    assert result['status']=='NO_BOUNDARY_OBSERVED'
    assert process.killed==0 and process.suspended==process.resumed==1


def test_verified_closed_boundary_stops_only_owned_coordinator_and_preserves_ledger(tmp_path,monkeypatch):
    report=prepared(tmp_path)
    ledger=(tmp_path/'case_budget.json').read_bytes()
    process=Process(tmp_path);mock_process(monkeypatch,process)
    result=boundary.stop_at_boundary(tmp_path,42,report,.1)
    assert result['status']=='PASS' and result['scientific_workers_terminated']==0
    assert process.killed==1 and process.resumed==0
    assert (tmp_path/'case_budget.json').read_bytes()==ledger
    assert result['ledger']['elapsed_seconds']==123.


def test_unowned_pid_is_rejected_before_any_process_control(tmp_path,monkeypatch):
    report=prepared(tmp_path)
    process=Process(tmp_path)
    process.cmdline=lambda:['python','unrelated.py']
    mock_process(monkeypatch,process)
    with pytest.raises(ValueError,match='explicitly owned'):
        boundary.stop_at_boundary(tmp_path,42,report,.1)
    assert process.suspended==process.resumed==process.killed==0
