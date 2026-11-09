"""Post-export runtime evidence; kept outside the frozen scientific archive."""
import ctypes
import math
import os
import platform
from pathlib import Path
from time import perf_counter

from ..util import read_json, sha256_file, write_json


def peak_process_bytes():
    if os.name == 'nt':
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_ = [('cb', wintypes.DWORD), ('faults', wintypes.DWORD)] + [
                (name, ctypes.c_size_t) for name in ('peak', 'working', 'peak_pool', 'pool',
                                                     'peak_nonpaged', 'nonpaged', 'pagefile', 'peak_pagefile')]
        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        query = ctypes.WinDLL('psapi').GetProcessMemoryInfo
        query.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        if not query(wintypes.HANDLE(-1), ctypes.byref(counters), counters.cb):
            raise ctypes.WinError()
        return counters.peak
    import resource
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value * (1 if platform.system() == 'Darwin' else 1024))


def write_performance(output, request, analysis, started, export_started, cache_was_empty):
    output = Path(output)
    ledger_path = output / 'ordinary_candidate_ledger.json'
    ledger = read_json(ledger_path) if ledger_path.is_file() else {}
    fingerprint = {}
    for field in ('structure_path', 'confidence_path'):
        source = getattr(request, field, None)
        fingerprint[field] = sha256_file(Path(source)) if source and Path(source).is_file() else None
    fingerprint.update({key: getattr(request, key) for key in ('chain', 'input_type', 'analysis_target')})
    receipt = {'schema': 'site_first_performance_v1', 'protocol': analysis.get('execution_protocol'),
        'input': fingerprint, 'hardware': {'platform': platform.platform(), 'machine': platform.machine(),
                                         'processor': platform.processor(), 'logical_cpus': os.cpu_count()},
        'cache_condition': request.metadata.get('benchmark_cache_condition', 'UNSPECIFIED'),
        'models_disabled': bool(request.metadata.get('disable_external_engines')),
        'score_cache_initially_empty': cache_was_empty,
        'total_seconds': perf_counter() - started, 'export_seconds': perf_counter() - export_started,
        'peak_process_memory_bytes': peak_process_bytes(),
        'peak_memory_scope': 'process lifetime; use a fresh process for each benchmark case',
        'disk_bytes_before_receipt': sum(p.stat().st_size for p in output.rglob('*') if p.is_file()),
        'candidate_timings': ledger.get('timings', {}),
        'spurs': ledger.get('spurs', {}), 'foldx': ledger.get('foldx', {}),
        'worker_metrics': ledger.get('worker_metrics', {}), 'counts': ledger.get('counts', {}),
        'scientific_acceptance': 'NOT_ASSESSED_BY_RUNTIME_RECEIPT',
        'scope': 'post-export sidecar; does not modify scientific manifest or ZIP'}
    write_json(output / 'site_first_performance.json', receipt)
    return receipt


def compare_performance(baseline, optimized):
    def protocol_name(row):
        protocol = row.get('protocol') or {}
        return protocol.get('name') if isinstance(protocol, dict) else protocol
    if protocol_name(baseline) != 'legacy' or protocol_name(optimized) != 'ordinary_site_first_v2':
        raise ValueError('Benchmark requires legacy baseline and site-first optimized protocol')
    for field in ('input', 'hardware', 'cache_condition'):
        if baseline[field] != optimized[field]:
            raise ValueError(f'Unmatched benchmark {field}')
    if baseline['cache_condition'] not in {'cold', 'warm'}:
        raise ValueError('Benchmark cache condition must be explicitly cold or warm')
    if not all(row.get('models_disabled') is False for row in (baseline, optimized)):
        raise ValueError('Both benchmark runs must enable real model engines')
    if not baseline['input'].get('structure_path'):
        raise ValueError('Benchmark requires a verified saved structure hash')
    if baseline['input'].get('input_type') == 'af3_upload' and not baseline['input'].get('confidence_path'):
        raise ValueError('AF3 benchmark requires the confidence file hash')
    if not all(row.get('score_cache_initially_empty') is True for row in (baseline, optimized)):
        raise ValueError('Both runs must start without target score caches')
    if not all(isinstance(row['total_seconds'], (int, float)) and math.isfinite(row['total_seconds'])
               and row['total_seconds'] > 0 for row in (baseline, optimized)):
        raise ValueError('Benchmark times must be positive')
    ratio = optimized['total_seconds'] / baseline['total_seconds']
    return {'runtime_gate': 'PROVISIONAL_PASS' if ratio <= .5 else 'NOT_MET', 'time_ratio': ratio,
            'reduction_fraction': 1 - ratio, 'scientific_acceptance': 'SEPARATE_REQUIRED',
            'remaining_verification': 'Verify actual GPU/model identities, repair cache condition, real model execution and scientific evidence independently.'}
