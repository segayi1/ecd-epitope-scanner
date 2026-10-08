"""Bounded independent APBS solves and exact repeated solvent-query reuse."""
from __future__ import annotations
import ctypes
import hashlib
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from ..util import sha256_file,write_json
from .tools import fingerprint
from .pqr_surface import pqr_surface


def available_memory():
    try:
        if os.name=='nt':
            class Memory(ctypes.Structure):
                _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in
                    ('total','available','page_total','page_available','virtual_total','virtual_available','extended')]
            value=Memory();value.length=ctypes.sizeof(value)
            if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(value)):return 0
            return value.available
        return os.sysconf('SC_AVPHYS_PAGES')*os.sysconf('SC_PAGE_SIZE')
    except (ValueError,OSError,AttributeError):return 0


def potential_pair(tools,wt,mut,grid):
    # Identical content is one job, preventing simultaneous writers to one cache.
    paths=[Path(wt),Path(mut)];digests=[sha256_file(p) for p in paths]
    unique=dict(zip(digests,paths))
    configured=max(1,min(2,int(getattr(tools,'config',{}).get('apbs_workers',2))))
    estimate=int(np.prod(grid['counts']))*256+256*1024**2
    memory=available_memory()
    if getattr(tools,'use_wsl',False):
        try:
            output=subprocess.run(['wsl.exe','-d',tools.distro,'--exec','/usr/bin/cat','/proc/meminfo'],
                capture_output=True,timeout=10,check=True).stdout.decode()
            linux_available=int(next(row.split()[1] for row in output.splitlines() if row.startswith('MemAvailable:')))*1024
            memory=min(memory,linux_available)
        except (OSError,subprocess.SubprocessError,ValueError,StopIteration):memory=0
    workers=max(1,min(configured,len(unique),int(memory*.5//estimate)))
    # Only the APBS subprocess receives thread controls; do not alter global env.
    if getattr(tools,'use_wsl',False) and not getattr(tools,'_apbs_thread_control',False):
        original_run=tools.run
        def run(argv,*args,**kwargs):
            if argv[0]==tools.apbs:
                argv=['/usr/bin/env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',*argv]
            return original_run(argv,*args,**kwargs)
        tools.run=run;tools._apbs_thread_control=True
    if hasattr(tools,'identity'):tools.identity(tools.apbs) # identity cache populated before threads
    started=time.perf_counter()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures={digest:pool.submit(tools.potential,path,grid) for digest,path in unique.items()}
        values={digest:future.result() for digest,future in futures.items()}
    if hasattr(tools,'activity'):tools.activity.append({'stage':'apbs_pair','workers':workers,
        'unique_jobs':len(unique),'estimated_bytes_per_job':estimate,'available_memory_bytes':memory,
        'subprocess_threads':1 if getattr(tools,'use_wsl',False) else 'unchanged',
        'runtime_seconds':time.perf_counter()-started})
    return tuple(values[d] for d in digests)


def solvent_distances(tools,pqr_atoms,queries,*,return_audit=False):
    queries=np.asarray(queries,float)
    module=Path(__file__).parent
    contract={'atoms':fingerprint([(list(k),v) for k,v in sorted(pqr_atoms.items())]),
        'queries':hashlib.sha256(queries.tobytes()).hexdigest(),'query_shape':list(queries.shape),
        'geometry':sha256_file(module/'geometry.py'),'surface':sha256_file(module/'pqr_surface.py'),
        'nano':tools.identity(tools.nano),'version':'exact_solvent_sdf_cache_1'}
    folder=tools.cache/('solvent_queries_'+fingerprint(contract));folder.mkdir(exist_ok=True)
    data=folder/'distances.npz';manifest=folder/'complete.json'
    if data.exists() and manifest.exists():
        audit=json.loads(manifest.read_text())
        if audit['contract']!=contract or audit['sha256']!=sha256_file(data):raise ValueError('solvent query cache changed')
        with np.load(data) as saved:distance=saved['distance'];valid=saved['valid']
        if distance.shape!=(len(queries),) or valid.shape!=distance.shape:raise ValueError('solvent cache shape changed')
        tools.calls['cache_hits']+=1
        result_audit={key:audit[key] for key in ('status','reason') if key in audit}
        result_audit.setdefault('status','COMPLETE')
        return (distance,valid,result_audit) if return_audit else (distance,valid)
    started=time.perf_counter()
    try:
        distance,valid=pqr_surface(tools,pqr_atoms).signed_distance(queries)
        result_audit={'status':'COMPLETE'}
    except RuntimeError as exc:
        if str(exc)!='PQR SES is not closed':raise
        distance=np.full(len(queries),np.nan);valid=np.zeros(len(queries),bool)
        result_audit={'status':'UNRESOLVED','reason':str(exc)}
    np.savez_compressed(data,distance=distance,valid=valid)
    write_json(manifest,{'contract':contract,'sha256':sha256_file(data),**result_audit})
    if hasattr(tools,'activity'):tools.activity.append({'stage':'solvent_sdf','runtime_seconds':time.perf_counter()-started,
        'queries':len(queries),**result_audit})
    return (distance,valid,result_audit) if return_audit else (distance,valid)
