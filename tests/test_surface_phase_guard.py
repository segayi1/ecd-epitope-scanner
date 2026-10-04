import json
from ecd_epitope_scanner.scripts.surface_phase_guard import charged_elapsed,guard_elapsed,run_guarded


def test_unclosed_pause_is_charged_not_refunded():
    assert charged_elapsed({'active':True,'elapsed_seconds':40,'checkpoint_wall_time':100},200)==140
    assert charged_elapsed({'active':False,'elapsed_seconds':40},200)==40
    assert guard_elapsed({'active':False,'elapsed_seconds':140},
                         {'active':False,'elapsed_seconds':40},100,200)==140


def test_exhausted_phase_never_launches_worker(tmp_path,monkeypatch):
    from ecd_epitope_scanner.scripts import surface_phase_guard as guard
    path=tmp_path/'surface_production_cumulative_budget.json'
    path.write_text(json.dumps({'active':False,'elapsed_seconds':7201,'limit_seconds':7200}))
    def forbidden(*args,**kwargs):raise AssertionError('worker must not start')
    monkeypatch.setattr(guard.subprocess,'Popen',forbidden)
    assert run_guarded(['unused'],cwd=tmp_path,folder=tmp_path,phase='production',limit=7200)==124
    assert not json.loads(path.read_text())['active']
    assert json.loads((tmp_path/'case_acceptance.json').read_text())['status']=='INCOMPLETE_RUNTIME_EXHAUSTED'


def test_blocked_child_without_ledger_is_stopped_and_launch_time_charged(tmp_path):
    import sys
    result=run_guarded([sys.executable,'-c','import time; time.sleep(30)'],
        cwd=tmp_path,folder=tmp_path,phase='production',limit=.1)
    assert result==124
    ledger=json.loads((tmp_path/'surface_production_cumulative_budget.json').read_text())
    assert ledger['elapsed_seconds']>=.1 and not ledger['active']


def test_absent_null_and_exhausted_calibration_prevents_unverifiable_production(tmp_path,monkeypatch):
    from ecd_epitope_scanner.scripts import surface_phase_guard as guard
    calibration=tmp_path/'calibration';calibration.mkdir()
    (calibration/'surface_calibration_cumulative_budget.json').write_text(json.dumps(
        {'active':False,'elapsed_seconds':14401,'limit_seconds':14400}))
    def forbidden(*args,**kwargs):raise AssertionError('cannot complete formal evidence')
    monkeypatch.setattr(guard.subprocess,'Popen',forbidden)
    production=tmp_path/'production'
    assert run_guarded(['unused'],cwd=tmp_path,folder=production,phase='production',limit=7200)==125
    assert json.loads((production/'case_acceptance.json').read_text())['status']=='INCOMPLETE_CALIBRATION_BUDGET_EXHAUSTED'
    # A saved null array means absence was not proven; do not infer invalidity
    # or waive its scientific requalification from this cheap preflight.
    null=calibration/'surface_perturbation/calibration';null.mkdir(parents=True)
    (null/'saved.npz').write_bytes(b'identity must be checked by the engine')
    assert not guard.exhausted_missing_calibration(production)


def test_budget_reader_retries_sharing_denial_without_losing_time(tmp_path,monkeypatch):
    from pathlib import Path
    from ecd_epitope_scanner.surface_perturbation.runtime import read_budget
    p=tmp_path/'budget.json';p.write_text(json.dumps({'elapsed_seconds':123,'active':True}))
    original=Path.read_text;attempts=[]
    def transient(path,*args,**kwargs):
        if path==p:
            attempts.append(1)
            if len(attempts)<3:raise PermissionError('native sharing denial')
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'read_text',transient)
    assert read_budget(p)=={'elapsed_seconds':123,'active':True}
    assert len(attempts)==3


def test_guard_read_failure_stops_owned_worker_and_charges_launch(tmp_path,monkeypatch):
    import sys,psutil
    from ecd_epitope_scanner.scripts import surface_phase_guard as guard
    original=guard.read_budget;calls=[];started=[];popen=guard.subprocess.Popen
    def start(*args,**kwargs):
        p=popen(*args,**kwargs);started.append(p.pid);return p
    def fail(path,default=None):
        calls.append(1)
        if len(calls)>1:raise PermissionError('persistent budget read failure')
        return original(path,default)
    monkeypatch.setattr(guard,'read_budget',fail)
    monkeypatch.setattr(guard.subprocess,'Popen',start)
    assert run_guarded([sys.executable,'-c','import time; time.sleep(30)'],
        cwd=tmp_path,folder=tmp_path,phase='production',limit=60)==126
    assert not psutil.pid_exists(started[0])
    ledger=json.loads((tmp_path/'surface_production_cumulative_budget.json').read_text())
    assert not ledger['active'] and ledger['elapsed_seconds']>0
    report=json.loads((tmp_path/'phase_result.json').read_text())
    assert report['status']=='INCOMPLETE_GUARD_FAILURE' and report['owned_worker_stopped']


def test_normal_exit_persists_launch_overhead_and_keeps_repeat_budget(tmp_path):
    import sys
    path=tmp_path/'surface_production_cumulative_budget.json'
    path.write_text(json.dumps({'active':False,'elapsed_seconds':7.,'limit_seconds':60}))
    code=("import json,time;from pathlib import Path;time.sleep(.3);"
          f"Path({str(path)!r}).write_text(json.dumps({{'active':False,'elapsed_seconds':7.,"
          "'limit_seconds':60,'checkpoint_wall_time':time.time()}))")
    assert run_guarded([sys.executable,'-c',code],cwd=tmp_path,folder=tmp_path,
                       phase='production',limit=60)==0
    first=json.loads(path.read_text())
    assert not first['active'] and first['elapsed_seconds']>=7.3
    assert run_guarded([sys.executable,'-c','import time;time.sleep(.3)'],
                      cwd=tmp_path,folder=tmp_path,phase='production',limit=60)==0
    second=json.loads(path.read_text())
    assert not second['active'] and second['elapsed_seconds']>=first['elapsed_seconds']+.3


def test_completed_worker_cannot_accept_understated_internal_clock(tmp_path,monkeypatch):
    from ecd_epitope_scanner.scripts import surface_phase_guard as guard
    class Finished:
        returncode=0
        def poll(self):return 0
    ticks=iter([100.,101.,101.])
    monkeypatch.setattr(guard.time,'time',lambda:next(ticks,101.))
    monkeypatch.setattr(guard.subprocess,'Popen',lambda *a,**k:Finished())
    assert guard._run_guarded(['already exited'],cwd=tmp_path,folder=tmp_path,
                              phase='production',limit=.5)==124
    budget=json.loads((tmp_path/'surface_production_cumulative_budget.json').read_text())
    assert budget['elapsed_seconds']>=1. and not budget['active']
    assert json.loads((tmp_path/'phase_result.json').read_text())['status']=='INCOMPLETE_RUNTIME_EXHAUSTED'
