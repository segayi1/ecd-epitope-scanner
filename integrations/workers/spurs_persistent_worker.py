"""Experimental line-oriented SPURS worker with content-verified model reuse."""
import contextlib
import hashlib
import io
import json
import os
import subprocess
import uuid
import sys
from pathlib import Path

import spurs_gpu_worker as worker


class ModelCache:
    def __init__(self, load):
        self.load = load
        self.key = None
        self.value = None
        self.files = {}

    def clear(self):
        self.key = self.value = None

    def verify(self):
        for path, (_, expected) in self.files.items():
            with path.open('rb') as stream:
                if hashlib.file_digest(stream, 'sha256').hexdigest() != expected:
                    raise ValueError('model files changed during persistent process')

    def __call__(self, model_dir, device):
        paths = [Path(model_dir) / '.hydra/config.yaml', Path(model_dir) / 'checkpoints/best.ckpt',
                 Path(model_dir).parent / 'esm2_t33_650M_UR50D.pt']
        digests = []
        for path in paths:
            stat = path.stat()
            identity = (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
            previous = self.files.get(path)
            if previous is not None and previous[0] != identity:
                raise ValueError('model files changed during persistent process')
            if previous is None:
                with path.open('rb') as stream:
                    digest = hashlib.file_digest(stream, 'sha256').hexdigest()
                self.files[path] = (identity, digest)
            digests.append(self.files[path][1])
        key = (str(Path(model_dir).resolve()), device, *digests)
        if key != self.key:
            self.clear()
            value = self.load(model_dir, device)
            self.key, self.value = key, value
        return self.value


def serve(stdin, stdout):
    cache = ModelCache(worker._load_multi)
    original_load, original_cleanup = worker._load_multi, worker.cleanup_cuda
    original_hash = worker.sha256_model_files
    file_hashes = {}
    def measured_files(paths, root):
        identities = tuple((str(p), p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_ctime_ns)
                           for p in paths)
        key = tuple(map(str, paths)), str(root)
        if key in file_hashes:
            previous, value = file_hashes[key]
            if identities != previous:
                raise ValueError('model files changed during persistent process')
            return value
        value = original_hash(paths, root)
        file_hashes[key] = (identities, value)
        return value
    def verify_files():
        cache.verify()
        for (paths, root), (_, expected) in file_hashes.items():
            if original_hash([Path(p) for p in paths], Path(root)) != expected:
                raise ValueError('model files changed during persistent process')
    worker.sha256_model_files = measured_files
    code_identity = None
    environment = None
    instance = uuid.uuid4().hex
    request_count = 0
    def runtime_identity():
        import importlib.metadata as metadata
        import torch
        driver = subprocess.run(['/usr/lib/wsl/lib/nvidia-smi',
            '--query-gpu=driver_version,name', '--format=csv,noheader'],
            capture_output=True, text=True, check=True, timeout=15)
        return {'packages': {d.metadata['Name']: d.version for d in metadata.distributions()},
                'python': sys.version, 'torch': torch.__version__, 'cuda_build': torch.version.cuda,
                'gpu_capability': torch.cuda.get_device_capability(0) if torch.cuda.is_available() else None,
                'deterministic_algorithms': torch.are_deterministic_algorithms_enabled(),
                'gpu_driver': driver.stdout.strip()}
    def source_identity(root):
        digest = hashlib.sha256()
        paths = sorted((root / 'vendor/spurs').rglob('*.py')) + [Path(__file__), Path(worker.__file__), Path(worker.__file__).with_name('gpu_runtime.py')]
        for path in paths:
            digest.update(str(path).encode())
            digest.update(path.read_bytes())
        return digest.hexdigest()
    def cleanup(torch):
        cache.clear()
        original_cleanup(torch)
    worker._load_multi, worker.cleanup_cuda = cache, cleanup
    try:
        for line in stdin:
            try:
                envelope = json.loads(line)
                request_id = envelope['request_id']
                request = envelope['request']
                if not isinstance(request, dict):
                    raise ValueError('request must be an object')
                if request.get('_finalize'):
                    verify_files()
                    if environment is not None and runtime_identity() != environment:
                        raise ValueError('runtime environment changed during persistent process')
                    if code_identity and source_identity(Path(code_identity[0])) != code_identity[1]:
                        raise ValueError('worker source changed during persistent process')
                    stdout.write(json.dumps({'request_id': request_id, 'result': {'status': 'success', 'verified': True}}) + '\n')
                    stdout.flush()
                    break
                request.setdefault('parameters', {}).update(cuda_oom_retries=0, allow_cpu_fallback_on_oom=False)
                root = Path(request.get('parameters', {}).get('model_assets_root') or worker._project_root()).resolve()
                current_identity = (str(root), source_identity(root))
                if code_identity is None:
                    code_identity = current_identity
                if current_identity != code_identity:
                    raise ValueError('worker source changed during persistent process')
                if request.get('_environment'):
                    if environment is None:
                        with contextlib.redirect_stdout(sys.stderr):
                            environment = runtime_identity()
                    stdout.write(json.dumps({'request_id': request_id, 'result': {
                        'status': 'success', 'environment': environment}}) + '\n')
                    stdout.flush()
                    continue
                buffer = io.StringIO()
                old_stdin = sys.stdin
                try:
                    sys.stdin = io.StringIO(json.dumps(request))
                    with contextlib.redirect_stdout(buffer):
                        worker.main()
                finally:
                    sys.stdin = old_stdin
                result = json.loads(buffer.getvalue())
                request_count += 1
                result.setdefault('provenance', {})['worker_lifecycle'] = {
                    'pid': os.getpid(), 'instance_id': instance, 'request_count': request_count,
                    'scheduler_version': 'strict_cohort_persistent_v1'}
                if source_identity(root) != code_identity[1]:
                    raise ValueError('worker source changed during prediction')
                if request.get('mutations') and result.get('status') == 'success' and cache.key:
                    weights = result.get('provenance', {}).get('model_files_sha256', {})
                    checkpoint = 'models/spurs/spurs_multi/checkpoints/best.ckpt'
                    esm = 'models/spurs/esm2_t33_650M_UR50D.pt'
                    if weights.get(checkpoint) != cache.key[-2] or weights.get(esm) != cache.key[-1]:
                        raise ValueError('model files changed during prediction')
                if result.get('status') != 'success':
                    cache.clear()
                response = {'request_id': request_id, 'result': result}
            except Exception as exc:
                cache.clear()
                response = {'request_id': None, 'error': str(exc)}
            stdout.write(json.dumps(response, sort_keys=True) + '\n')
            stdout.flush()
            if response.get('error') or response.get('result', {}).get('status') != 'success':
                break
    finally:
        cache.clear()
        worker._load_multi, worker.cleanup_cuda = original_load, original_cleanup
        worker.sha256_model_files = original_hash


if __name__ == '__main__':
    if '--control-file' in sys.argv:
        if os.getsid(0) != os.getpid():
            os.setsid()
        if os.getpgrp() != os.getpid():
            raise RuntimeError('persistent worker must own its process group')
        control = Path(sys.argv[sys.argv.index('--control-file') + 1])
        temporary = control.with_suffix('.pending')
        temporary.write_text(json.dumps({'pid': os.getpid(), 'start_ticks': Path('/proc/self/stat').read_text().split()[21]}))
        temporary.replace(control)
    serve(sys.stdin, sys.stdout)
