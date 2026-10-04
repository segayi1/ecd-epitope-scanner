"""External wall-time guard for frozen phase processes, including blocked workers."""
import subprocess
import time
from pathlib import Path
from ..surface_perturbation.runtime import RunBudget,read_budget
from ..util import write_json


def charged_elapsed(ledger,now):
    return float(ledger.get('elapsed_seconds',0.))+(
        max(0.,now-ledger['checkpoint_wall_time']) if ledger.get('active') else 0.)


def guard_elapsed(ledger,initial,start,now):
    # The closed child ledger already includes its execution: use the launch
    # envelope as a lower bound, never add the same interval a second time.
    return max(charged_elapsed(ledger,now),charged_elapsed(initial,start)+max(0.,now-start))


def exhausted_missing_calibration(folder):
    calibration=Path(folder).parent/'calibration'
    path=calibration/'surface_calibration_cumulative_budget.json'
    if not path.is_file():return False
    ledger=read_budget(path)
    return (not ledger.get('active') and ledger['elapsed_seconds']>=ledger['limit_seconds']
            and not any((calibration/'surface_perturbation/calibration').glob('*.npz')))


def stop_owned_tree(process):
    import psutil
    try:parent=psutil.Process(process.pid);children=parent.children(recursive=True)
    except psutil.NoSuchProcess:return
    # This exact child was started here; never terminate an unrelated WSL service.
    for child in [parent,*children]:
        try:child.terminate()
        except psutil.NoSuchProcess:pass
    _,alive=psutil.wait_procs([parent,*children],timeout=2)
    for child in alive:
        try:child.kill()
        except psutil.NoSuchProcess:pass


def settle_completed_phase(code, *, ledger_path, initial, start, folder, phase, limit):
    """Persist the launch envelope even when a worker exits normally."""
    try:
        observed=guard_elapsed(read_budget(ledger_path,initial),initial,start,time.time())
        budget=RunBudget(ledger_path,limit)
        budget.carried=max(budget.carried,observed)
        budget.save(False)
        if budget.elapsed<limit:return code
        report={'status':'INCOMPLETE_RUNTIME_EXHAUSTED','phase':phase,
                'reason':'Completed worker exceeded the cumulative launch envelope.',
                'elapsed_seconds':budget.elapsed,'limit_seconds':limit}
        result=124
    except (OSError,RuntimeError) as exc:
        report={'status':'INCOMPLETE_GUARD_FAILURE','phase':phase,'reason':str(exc),
                'owned_worker_stopped':True}
        result=126
    write_json(Path(folder)/'case_acceptance.json',report)
    write_json(Path(folder)/'phase_result.json',report)
    return result


def run_guarded(command,*,cwd,folder,phase,limit):
    from ..surface_perturbation.power_guard import prevent_sleep
    with prevent_sleep() as report:
        write_json(Path(folder)/'power_guard.json',report)
        return _run_guarded(command,cwd=cwd,folder=folder,phase=phase,limit=limit)


def _run_guarded(command,*,cwd,folder,phase,limit):
    folder=Path(folder);ledger_path=folder/f'surface_{phase}_cumulative_budget.json'
    initial=read_budget(ledger_path)
    if phase=='production' and exhausted_missing_calibration(folder):
        report={'status':'INCOMPLETE_CALIBRATION_BUDGET_EXHAUSTED','phase':phase,
            'reason':'No matched null arrays and no remaining calibration budget; formal physical validation cannot complete.'}
        write_json(folder/'case_acceptance.json',report)
        write_json(folder/'phase_result.json',report)
        return 125
    limit=min(float(limit),float(initial.get('limit_seconds',limit)))
    start=time.time()
    process=None
    if charged_elapsed(initial,start)<limit:
        process=subprocess.Popen(command,cwd=cwd)
        try:
            while process.poll() is None:
                ledger=read_budget(ledger_path,initial)
                # Before the child opens its budget, launch overhead still counts.
                elapsed=guard_elapsed(ledger,initial,start,time.time())
                if elapsed>=limit:break
                time.sleep(.25)
            else:return settle_completed_phase(process.returncode,ledger_path=ledger_path,
                initial=initial,start=start,folder=folder,phase=phase,limit=limit)
        except Exception as exc:
            stop_owned_tree(process)
            # A failing guard must never leave an unmanaged scientific worker.
            # Leave an unreadable active ledger intact: restart still charges it.
            report={'status':'INCOMPLETE_GUARD_FAILURE','phase':phase,'reason':str(exc),
                    'launch_elapsed_lower_bound':guard_elapsed(initial,initial,start,time.time()),
                    'owned_worker_stopped':True}
            try:
                budget=RunBudget(ledger_path,limit)
                budget.carried=max(budget.carried,report['launch_elapsed_lower_bound'])
                budget.save(False)
            except (OSError,RuntimeError) as close_error:
                report['ledger_closure_error']=str(close_error)
            write_json(folder/'phase_result.json',report)
            return 126
        if process.poll() is not None:return settle_completed_phase(process.returncode,
            ledger_path=ledger_path,initial=initial,start=start,folder=folder,phase=phase,limit=limit)
        stop_owned_tree(process)
    # Explicitly close an interrupted ledger; unclosed wall time is never refunded.
    ledger=read_budget(ledger_path,initial)
    observed=guard_elapsed(ledger,initial,start,time.time())
    budget=RunBudget(ledger_path,limit)
    # Include launch/wait time even if a blocked child never opened its ledger.
    budget.carried=max(budget.carried,observed)
    budget.save(False)
    report={'status':'INCOMPLETE_RUNTIME_EXHAUSTED','phase':phase,
        'reason':'External cumulative wall-time guard; physical checkpoints retained, no acceptance.',
        'elapsed_seconds':budget.elapsed,'limit_seconds':limit}
    write_json(folder/'case_acceptance.json',report)
    write_json(folder/'phase_result.json',report)
    return 124
