import importlib.util
import io
import json
import sys
import pytest
from pathlib import Path

WORKERS = Path(__file__).resolve().parents[1] / 'integrations' / 'workers'
sys.path.insert(0, str(WORKERS))
spec = importlib.util.spec_from_file_location('spurs_persistent_worker', WORKERS / 'spurs_persistent_worker.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_model_cache_checks_file_contents(tmp_path):
    model = tmp_path / 'spurs_multi'
    for path in (model / '.hydra/config.yaml', model / 'checkpoints/best.ckpt', tmp_path / 'esm2_t33_650M_UR50D.pt'):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'one')
    calls = []
    cache = module.ModelCache(lambda path, device: calls.append((path, device)) or object())
    first = cache(model, 'cuda')
    assert cache(model, 'cuda') is first and len(calls) == 1
    (model / 'checkpoints/best.ckpt').write_bytes(b'two')
    with pytest.raises(ValueError, match='model files changed'):
        cache(model, 'cuda')
    with pytest.raises(ValueError, match='model files changed'):
        cache.verify()
    assert len(calls) == 1


def test_stream_preserves_request_identity_and_clears_on_failure(monkeypatch):
    def fake_main():
        request = json.load(sys.stdin)
        print(json.dumps({'status': request['status'], 'records': request.get('records', [])}))
    monkeypatch.setattr(module.worker, 'main', fake_main)
    payload = [{'request_id': 'a', 'request': {'status': 'success', 'records': [1]}},
               {'request_id': 'b', 'request': {'status': 'failed'}},
               {'request_id': 'c', 'request': {'status': 'success', 'records': [3]}}]
    output = io.StringIO()
    module.serve(io.StringIO(''.join(json.dumps(x) + '\n' for x in payload)), output)
    rows = [json.loads(line) for line in output.getvalue().splitlines()]
    assert [x['request_id'] for x in rows] == ['a', 'b']
    assert rows[0]['result']['records'] == [1]
