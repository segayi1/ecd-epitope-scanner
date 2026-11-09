"""Persist cumulative runtime across normal restarts; bounded subprocesses."""
from __future__ import annotations

import json
import time
import os
import threading
import math
from contextvars import ContextVar
from pathlib import Path
from ..util import write_json

ACTIVE_BUDGET = ContextVar('surface_perturbation_budget', default=None)


def read_budget(path, default=None):
    """Retry transient Windows sharing denials without substituting empty evidence."""
    path=Path(path)
    for attempt in range(12):
        try:
            return json.loads(path.read_text(encoding='utf-8')) if path.exists() else (default or {})
        except PermissionError:
            if attempt==11:raise
            time.sleep(min(.005*2**attempt,.1))


class RunBudget:
    def __init__(self, path: Path, limit=7200.):
        self.write_lock=threading.RLock()
        self.path = path
        path.parent.mkdir(parents=True,exist_ok=True)
        self.lock=open(path.with_suffix('.lock'),'a+b')
        try:
            if os.fstat(self.lock.fileno()).st_size==0:
                self.lock.write(b'1');self.lock.flush()
            self.lock.seek(0)
            if os.name=='nt':
                import msvcrt
                msvcrt.locking(self.lock.fileno(),msvcrt.LK_NBLCK,1)
            else:
                import fcntl
                fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except OSError:
            self.lock.close()
            raise RuntimeError('CASE_ALREADY_RUNNING: cumulative budget is locked')
        old = read_budget(path)
        self.carried = float(old.get("elapsed_seconds", 0))
        # An unclosed interrupted session conservatively charges elapsed wall time.
        if old.get("active"):
            self.carried += max(0., time.time()-old["checkpoint_wall_time"])
        def seconds(value): return math.inf if value is None else float(value)
        self.limit = min(seconds(old.get("limit_seconds", limit)), seconds(limit))
        self.started = time.monotonic()
        self.save()

    @property
    def elapsed(self): return self.carried + max(0., time.monotonic()-self.started)

    def save(self, active=True):
        with self.write_lock:
            checkpoint=self.path.with_name(self.path.name+f'.checkpoint.{os.getpid()}.tmp')
            write_json(checkpoint, {"limit_seconds": None if math.isinf(self.limit) else self.limit, "elapsed_seconds": self.elapsed,
                                   "active": active, "checkpoint_wall_time": time.time()})
            for attempt in range(12):
                try:
                    checkpoint.replace(self.path)
                    break
                except PermissionError:
                    # Windows readers/indexers can briefly deny rename. Keep
                    # the previous complete ledger; never truncate it in place.
                    if attempt==11:raise
                    time.sleep(min(.005*2**attempt,.1))
            if not active and not self.lock.closed:
                self.lock.close()

    def remaining(self):
        from .case_resources import ACTIVE_STORAGE
        if ACTIVE_STORAGE.get() is not None:
            ACTIVE_STORAGE.get().check()
        self.save()
        remaining = self.limit-self.elapsed
        if remaining <= 0: raise TimeoutError("CUMULATIVE_RUNTIME_EXHAUSTED")
        return remaining
