import subprocess
from types import SimpleNamespace
from ecd_epitope_scanner.integrations import wsl_bridge


def test_path_conversion_cannot_consume_worker_requests(monkeypatch):
    def run(*args, **kwargs):
        assert kwargs['stdin'] is subprocess.DEVNULL
        return SimpleNamespace(returncode=0, stdout=b'/mnt/e/test\n', stderr=b'')
    monkeypatch.setattr(wsl_bridge.subprocess, 'run', run)
    assert wsl_bridge._wsl_path('E:/test', 'Ubuntu') == '/mnt/e/test'
