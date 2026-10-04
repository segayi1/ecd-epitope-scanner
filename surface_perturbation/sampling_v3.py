"""Audited native time seeding; no invented numeric seed or IID claim."""
import json
import os
import subprocess
import time
from pathlib import Path

from ..util import write_json

PROTOCOL = 'time_rotabase_v3'


def seeded_run(command, cwd, timeout):
    # One cross-process start lock per volume/tool asset tree. Hold it only
    # until Popen, not for the FoldX solve. The OS releases locks after a crash.
    lock_path = Path(command[0]).resolve().parent / 'time_seed_start.lock'
    state_path = lock_path.with_suffix('.json')
    started = time.monotonic()
    with lock_path.open('a+b') as handle:
        if os.fstat(handle.fileno()).st_size == 0:
            handle.write(b'1'); handle.flush()
        while True:
            try:
                handle.seek(0)
                if os.name == 'nt':
                    import msvcrt
                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except OSError:
                if time.monotonic() - started > timeout:
                    raise subprocess.TimeoutExpired(command, timeout)
                time.sleep(.05)
        previous = json.loads(state_path.read_text()) if state_path.exists() else {}
        delay = max(0., 2. - (time.time() - previous.get('last_start', 0)))
        if delay: time.sleep(delay)
        stamp = time.time()
        process = subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True)
        write_json(state_path, {'last_start': stamp, 'pid': process.pid})
        write_json(Path(cwd) / 'sampling_execution.json', {
            'protocol': PROTOCOL, 'command': command, 'start_wall_time': stamp,
            'pid': process.pid, 'internal_seed': None,
            'seed_source': 'FoldX timeSeedRotabase; internal seed not reported',
            'minimum_start_spacing_seconds': 2.,
            'interpretation': 'sidechain sampling, not biological or backbone independence'})
    try:
        stdout, stderr = process.communicate(timeout=max(.001, timeout-(time.monotonic()-started)))
    except subprocess.TimeoutExpired:
        process.kill(); process.communicate()
        raise
    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


def validate_roles(pairs, role, count, cohort=None):
    selected = [p for p in pairs if p.get('sampling_protocol') == PROTOCOL
                and p.get('sampling_role') == role
                and (cohort is None or p.get('sampling_cohort') == cohort)]
    slots = [p.get('sampling_slot') for p in selected]
    ids = [(p['batch_path'], p['group'], p['run']) for p in selected]
    if len(set(ids)) != len(ids) or len(set(slots)) != len(slots):
        raise ValueError('duplicate sampling execution or logical slot')
    ordered = sorted(selected, key=lambda p: p['sampling_slot'])
    if [p['sampling_slot'] for p in ordered] != list(range(len(ordered))):
        raise ValueError('sampling slots must be a fixed prefix')
    return ordered[:count]
