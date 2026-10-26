import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from ecd_epitope_scanner.joint_design_v4.models import PredictionResult
from ecd_epitope_scanner.surface_perturbation import upstream_cache as cache


@pytest.fixture
def cached_run(tmp_path, monkeypatch):
    root = tmp_path/'root'
    module = root/'surface_perturbation/upstream_cache.py'
    module.parent.mkdir(parents=True); module.write_text('')
    monkeypatch.setattr(cache, '__file__', str(module))
    cache._RUN_FILE_DIGESTS.clear()
    for name in ('spurs_gpu_worker.py', 'gpu_runtime.py'):
        p = root/'integrations/workers'/name
        p.parent.mkdir(parents=True, exist_ok=True);p.write_text(name)
    weights = {}
    for name in ('models/spurs/spurs/checkpoints/best.ckpt', 'models/spurs/esm2_t33_650M_UR50D.pt'):
        p = root/name;p.parent.mkdir(parents=True, exist_ok=True);p.write_bytes(name.encode())
        weights[name] = hashlib.sha256(p.read_bytes()).hexdigest()
    monkeypatch.setattr(cache.subprocess, 'run', lambda args, **kwargs: SimpleNamespace(
        stdout=b'{"python":"3.11","torch":"2.5","cuda_build":"12.1"}' if '-c' in args else b'driver,gpu'))
    structure = tmp_path/'wt.pdb';structure.write_text('WT')
    request = SimpleNamespace(structure_path=str(structure),chain_id='A',sequence='AC',
                              structure_type='pdb',mutations=[],parameters={},device='cuda',output_dir=str(tmp_path/'job'))
    class Adapter:
        cache_dir=tmp_path/'cache'
        calls=0
        actual='cuda'
        def predict(self, request):
            self.calls+=1
            return PredictionResult(status='success',records=[{'position':1,'ddg':.2}],provenance={
                'actual_device':self.actual,'model_files_sha256':weights})
        def checked_result(self, request, result):
            return result
    return root, Adapter(), request


def test_cache_corruption_is_preserved_and_recomputed(cached_run):
    _, adapter, request = cached_run
    cache.predict_cached(adapter,request)
    path=next((adapter.cache_dir/'verified_spurs').glob('*.json'))
    payload=json.loads(path.read_text());payload['result']['records'][0]['ddg']=100
    path.write_text(json.dumps(payload))
    result=cache.predict_cached(adapter,request)
    assert adapter.calls==2 and result.records[0]['ddg']==.2
    assert result.provenance['cache_miss_reason']=='CACHE_PAYLOAD_INVALID'
    assert len(list(path.parent.glob('*.invalid_*.json')))==1
    assert cache.predict_cached(adapter,request).provenance['cache_hit']


def test_changed_parameters_and_adapter_source_do_not_hit(cached_run):
    root, adapter, request = cached_run
    cache.predict_cached(adapter,request)
    hit=cache.predict_cached(adapter,request)
    assert hit.provenance['cache_hit'] and hit.provenance['weights_bytes_hashed']==0
    request.parameters['reuse_wt_encoding']=False
    cache.predict_cached(adapter,request)
    assert adapter.calls==2
    p=root/'integrations/_adapter.py';p.write_text('protocol change')
    cache.predict_cached(adapter,request)
    assert adapter.calls==3


def test_missing_runtime_identity_and_cpu_fallback_are_not_cached(cached_run,monkeypatch):
    _, adapter, request = cached_run
    adapter.actual='cpu'
    cache.predict_cached(adapter,request);cache.predict_cached(adapter,request)
    assert adapter.calls==2
    adapter.actual='cuda'
    monkeypatch.setattr(cache.subprocess,'run',lambda *a,**kw:SimpleNamespace(stdout=b'{}'))
    cache.predict_cached(adapter,request);cache.predict_cached(adapter,request)
    assert adapter.calls==4
