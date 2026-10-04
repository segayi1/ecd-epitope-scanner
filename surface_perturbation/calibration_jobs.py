"""Isolate unresolved null tasks without accepting partial failed evidence."""
import time
from contextlib import contextmanager
from pathlib import Path
from ..util import write_json
from .null_store import load_calibration


def usable_calibration(folder,protocol):
    if (Path(folder)/'task_failures'/f'{protocol}.json').exists():return {}
    return load_calibration(folder,protocol)


def clear_failure(folder,protocol):
    path=Path(folder)/'task_failures'/f'{protocol}.json'
    if path.exists():
        history=Path(folder)/'task_failure_history';history.mkdir(exist_ok=True)
        path.rename(history/f'{protocol}_{time.time_ns()}.json')


@contextmanager
def isolated_job(folder,protocol,candidate_id,trajectories,errors,requested_repeats):
    try:yield
    except (RuntimeError,ValueError) as exc:
        row={'candidate':candidate_id,'protocol_id':protocol,'status':'CALIBRATION_TASK_FAILED',
             'reason':str(exc),'requested_repeats':requested_repeats(),
             'local_scope_converged':False,'hard_evidence':'unknown until a successful matched recalibration'}
        path=Path(folder)/'task_failures'/f'{protocol}.json';path.parent.mkdir(parents=True,exist_ok=True)
        write_json(path,row);trajectories.append(row);errors.append({'candidate':candidate_id,'reason':str(exc),'stage':'null_calibration'})
