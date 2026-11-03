from ecd_epitope_scanner.integrations import persistent_runner, process_runner
from ecd_epitope_scanner.joint_design_v4.models import PredictionRequest


def test_persistent_scope_reuses_and_closes_worker(monkeypatch, tmp_path):
    made = []

    class FakeRunner:
        def __init__(self, command, timeout, log_path):
            self.command = command
            self.closed = False
            self.requests = []
            made.append(self)

        def run(self, request):
            if request.get('_finalize'):
                return {'verified': True}
            self.requests.append(request['request_id'])
            return {'status': 'success', 'records': [], 'provenance': {}}

        def close(self):
            self.closed = True

    monkeypatch.setattr(persistent_runner, 'PersistentRunner', FakeRunner)
    command = ['python', 'wsl_bridge.py', '--worker', 'spurs_gpu_worker.py']
    worker = process_runner.JsonWorkerRunner(command)
    with process_runner.persistent_spurs_scope():
        for request_id in ('one', 'two'):
            request = PredictionRequest('v1', request_id, 'x.pdb', 'A', mutations=[['A1W']])
            assert worker.run(request, log_path=tmp_path / 'log').provenance['persistent_worker']
        assert len(made) == 1 and made[0].requests == ['one', 'two']
        assert not made[0].closed
    assert made[0].closed
