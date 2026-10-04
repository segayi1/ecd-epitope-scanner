"""SPURS reuse only with measured checkpoint, code and input identities."""
from __future__ import annotations
import hashlib
import json
import subprocess
from pathlib import Path
from ..util import sha256_file,write_json
from ..joint_design_v4.models import PredictionResult
from .tools import fingerprint
from ..integrations.wsl_bridge import _decode_output
from .runtime import ACTIVE_BUDGET


def predict_cached(adapter,request):
    source_root=Path(__file__).resolve().parents[1]
    root=Path(request.parameters.get('model_assets_root') or source_root)
    model=root/'models'/'spurs'/('spurs_multi' if request.mutations else 'spurs')
    required=[model/'checkpoints'/'best.ckpt',root/'models'/'spurs'/'esm2_t33_650M_UR50D.pt']
    frozen=request.parameters.get('frozen_model_assets_path')
    if frozen and not all(p.is_file() for p in required):
        raise ValueError('Frozen SPURS model asset missing')
    if (adapter.cache_dir is None and not frozen) or not all(p.is_file() for p in required):
        return adapter.predict(request)
    weights={p.relative_to(root).as_posix():sha256_file(p) for p in required}
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
    configs={p.relative_to(root).as_posix():sha256_file(p) for p in model.rglob('*.yaml')}
    for p in (root/'vendor'/'spurs').rglob('*'):
        if p.is_file() and p.suffix in {'.yaml','.yml','.json','.toml'}:
            configs[p.relative_to(root).as_posix()]=sha256_file(p)
    budget=ACTIVE_BUDGET.get()
    try:
        distro=request.parameters.get('wsl_distro','Ubuntu-22.04')
        python=request.parameters.get('wsl_python','/home/segayi1/miniconda3/envs/spurs_v4/bin/python')
        # SPURS uses fair-esm; an optional, absent transformers distribution
        # must not disable all reuse. Capture the installed environment instead
        # of requiring packages that this worker does not import.
        script="import importlib.metadata as m,json;print(json.dumps({d.metadata['Name']:d.version for d in m.distributions()},sort_keys=True))"
        runtime=subprocess.run(['wsl.exe','-d',distro,'--exec',python,'-c',script],capture_output=True,
                               timeout=min(30,budget.remaining()) if budget else 30,check=True)
        versions=json.loads(_decode_output(runtime.stdout).strip())
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
    folder=Path(adapter.cache_dir)/'verified_spurs';folder.mkdir(parents=True,exist_ok=True)
    path=folder/(fingerprint(contract)+'.json')
    if path.exists():
        value=json.loads(path.read_text(encoding='utf-8'))
        if value['contract']!=contract: raise ValueError('SPURS cache identity mismatch')
        if value.get('result_sha256')!=fingerprint(value['result']): raise ValueError('SPURS cached result changed')
        data=value['result'];result=PredictionResult(**data)
        result.provenance.update(cache_hit=True,cache_contract_sha256=fingerprint(contract))
        return adapter.checked_result(request,result)
    result=adapter.predict(request)
    if result.status=='success' and not result.errors and result.provenance.get('model_files_sha256')==weights:
        from dataclasses import asdict
        # CPU fallback is preserved as evidence but not reused as a CUDA result.
        actual_device=result.provenance.get('actual_device')
        # Successful nvidia-smi above identifies the CUDA runtime selected by
        # auto. CPU fallback cannot be reused as that automatic CUDA result.
        if actual_device==request.device or (request.device=='auto' and actual_device=='cuda'):
            write_json(path,{'contract':contract,'result':asdict(result),'result_sha256':fingerprint(asdict(result))})
    return result
