import sys
import pytest

from ecd_epitope_scanner.integrations.base import ToolError
from ecd_epitope_scanner.integrations.persistent_runner import PersistentRunner


def test_multiple_requests_use_one_process(tmp_path):
    worker = tmp_path / 'worker.py'
    worker.write_text('import sys,json,os\nfor line in sys.stdin:\n x=json.loads(line); print(json.dumps({"request_id":x["request_id"],"result":{"pid":os.getpid(),"value":x["request"]["value"]}}),flush=True)\n')
    with PersistentRunner([sys.executable, str(worker)], 5) as runner:
        first = runner.run({'value': 1})
        second = runner.run({'value': 2})
        assert first['pid'] == second['pid'] and (first['value'], second['value']) == (1, 2)
    assert runner.process.poll() is not None


def test_identity_mismatch_stops_worker(tmp_path):
    worker = tmp_path / 'worker.py'
    worker.write_text('import sys,json\nfor line in sys.stdin: print(json.dumps({"request_id":"wrong","result":{}}),flush=True)\n')
    runner = PersistentRunner([sys.executable, str(worker)], 5)
    with pytest.raises(ToolError):
        runner.run({})
    assert runner.process.poll() is not None


@pytest.mark.parametrize('source', ['import sys\nsys.exit(0)\n',
    'import sys,time\nfor line in sys.stdin: time.sleep(30)\n'])
def test_eof_and_timeout_stop_process(tmp_path, source):
    worker = tmp_path / 'worker.py'
    worker.write_text(source)
    runner = PersistentRunner([sys.executable, str(worker)], .1)
    with pytest.raises(ToolError):
        runner.run({})
    assert runner.process.poll() is not None
