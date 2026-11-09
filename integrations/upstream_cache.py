"""SPURS reuse only with measured checkpoint, code and input identities."""
from __future__ import annotations
import hashlib
import json
import subprocess
import time
from uuid import uuid4
from pathlib import Path
from ..util import sha256_file,write_json
from ..joint_design_v4.models import PredictionResult
def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
from ..integrations.wsl_bridge import _decode_output
from .runtime import ACTIVE_BUDGET

_RUN_FILE_DIGESTS = {}


def _measured_digest(path, scope):
    """Reuse a measured digest only inside one task and unchanged file identity."""
    stat = path.stat()
    identity = (str(path.resolve()), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    key = (scope, identity)
    if scope and key in _RUN_FILE_DIGESTS:
        return _RUN_FILE_DIGESTS[key], 0
    value = sha256_file(path)
    after = path.stat()
    if (after.st_size, after.st_mtime_ns, after.st_ctime_ns) != identity[1:]:
        raise ValueError('SPURS model changed while fingerprinting')
    if scope:
        _RUN_FILE_DIGESTS[key] = value
    return value, stat.st_size


def predict_cached(adapter,request):
    started = time.perf_counter()
    source_root=Path(__file__).resolve().parents[1]
    root=Path(request.parameters.get('model_assets_root') or source_root)
    model=root/'models'/'spurs'/('spurs_multi' if request.mutations else 'spurs')
    required=[model/'checkpoints'/'best.ckpt',root/'models'/'spurs'/'esm2_t33_650M_UR50D.pt']
    frozen=request.parameters.get('frozen_model_assets_path')
    if frozen and not all(p.is_file() for p in required):
        raise ValueError('Frozen SPURS model asset missing')
    if (adapter.cache_dir is None and not frozen) or not all(p.is_file() for p in required):
        return adapter.predict(request)
    scope = str(getattr(request, 'output_dir', '') or '')
    measurements = {p: _measured_digest(p, scope) for p in required}
    weights={p.relative_to(root).as_posix():measurements[p][0] for p in required}
    hashed_bytes = sum(v[1] for v in measurements.values())
    if frozen:
        expected=json.loads(Path(frozen).read_text(encoding='utf-8'))['files']
        if any(expected.get(name)!=digest for name,digest in weights.items()):
            raise ValueError('Frozen SPURS model asset identity changed')
    if adapter.cache_dir is None:return adapter.predict(request)
    code=hashlib.sha256()
    for p in sorted((root/'vendor'/'spurs').rglob('*.py')):
        code.update(p.relative_to(root).as_posix().encode());code.update(p.read_bytes())
    for name in ('spurs_gpu_worker.py','gpu_runtime.py'):
        code.update((source_root/'integrations'/'workers'/name).read_bytes())
    for name in ('_adapter.py', 'base.py', 'spurs_adapter.py', 'wsl_bridge.py'):
        path = source_root/'integrations'/name
        if path.is_file():
            code.update(name.encode());code.update(path.read_bytes())
    configs={p.relative_to(root).as_posix():sha256_file(p) for p in model.rglob('*.yaml')}
    for p in (root/'vendor'/'spurs').rglob('*'):
        if p.is_file() and p.suffix in {'.yaml','.yml','.json','.toml'}:
            configs[p.relative_to(root).as_posix()]=sha256_file(p)
    budget=ACTIVE_BUDGET.get()
    try:
        distro=request.parameters.get('wsl_distro','Ubuntu-22.04')
        python=request.parameters.get('wsl_python','python')
        # SPURS uses fair-esm; an optional, absent transformers distribution
        # must not disable all reuse. Capture the installed environment instead
        # of requiring packages that this worker does not import.
        script=("import importlib.metadata as m,json,sys,torch;print(json.dumps({"
                "'packages':{d.metadata['Name']:d.version for d in m.distributions()},"
                "'python':sys.version,'torch':torch.__version__,'cuda_build':torch.version.cuda,"
                "'gpu_capability':torch.cuda.get_device_capability(0) if torch.cuda.is_available() else None,"
                "'deterministic_algorithms':torch.are_deterministic_algorithms_enabled()},sort_keys=True))")
        from ..integrations.process_runner import managed_spurs_environment
        versions = managed_spurs_environment(adapter, request)
        if versions is None:
            runtime=subprocess.run(['wsl.exe','-d',distro,'--exec',python,'-c',script],capture_output=True,
                                   timeout=min(30,budget.remaining()) if budget else 30,check=True)
            versions=json.loads(_decode_output(runtime.stdout).strip())
            if not isinstance(versions, dict) or not versions.get('python') or not versions.get('torch'):
                raise ValueError('SPURS runtime identity incomplete')
            gpu=subprocess.run(['wsl.exe','-d',distro,'--exec','/usr/lib/wsl/lib/nvidia-smi',
                                '--query-gpu=driver_version,name','--format=csv,noheader'],capture_output=True,
                               timeout=min(30,budget.remaining()) if budget else 30,check=True)
            versions['gpu_driver']=_decode_output(gpu.stdout).strip()
    except (OSError,ValueError,subprocess.SubprocessError):
        return adapter.predict(request)
    contract={'structure':sha256_file(Path(request.structure_path)),'chain':request.chain_id,
              'sequence':request.sequence,'structure_type':request.structure_type,'mutations':request.mutations,
              'parameters':{k:v for k,v in request.parameters.items() if k not in {'frozen_model_assets_path','stability_contract'}},
              'device':request.device,'weights':weights,
              'source_sha256':code.hexdigest(),'model_configuration':configs}
    contract['runtime']=versions
    contract['worker_protocol']=('json-spurs-v3-persistent-no-oom-retry' if
        versions.get('_persistent_worker') else 'json-spurs-v2-wt-encoding')
    folder=Path(adapter.cache_dir)/'verified_spurs';folder.mkdir(parents=True,exist_ok=True)
    path=folder/(fingerprint(contract)+'.json')
    cache_problem = None
    preflight = {'fingerprint_seconds':time.perf_counter()-started,
                 'weights_bytes_hashed':hashed_bytes,'fingerprint_scope':'task_file_identity',
                 'cache_contract_sha256':fingerprint(contract)}
    if path.exists():
        try:
            value=json.loads(path.read_text(encoding='utf-8'))
            if value['contract']!=contract: raise ValueError('SPURS cache identity mismatch')
            if value.get('result_sha256')!=fingerprint(value['result']): raise ValueError('SPURS cached result changed')
            data=value['result'];result=adapter.checked_result(request,PredictionResult(**data))
            if not result.ok: raise ValueError('SPURS cache payload failed validation')
            result.provenance.update(cache_hit=True, **preflight)
            return result
        except (KeyError, TypeError, ValueError, OSError):
            cache_problem = 'CACHE_PAYLOAD_INVALID'
            path.rename(path.with_name(path.stem+'.invalid_'+uuid4().hex+'.json'))
    result=adapter.predict(request)
    result.provenance.update(preflight)
    if cache_problem:
        result.provenance['cache_miss_reason']=cache_problem
    if result.status=='success' and not result.errors and result.provenance.get('model_files_sha256')==weights:
        from dataclasses import asdict
        # CPU fallback is preserved as evidence but not reused as a CUDA result.
        actual_device=result.provenance.get('actual_device')
        # Successful nvidia-smi above identifies the CUDA runtime selected by
        # auto. CPU fallback cannot be reused as that automatic CUDA result.
        if actual_device==request.device or (request.device=='auto' and actual_device=='cuda'):
            write_json(path,{'contract':contract,'result':asdict(result),'result_sha256':fingerprint(asdict(result))})
    return result
