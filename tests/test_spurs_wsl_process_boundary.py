import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from ecd_epitope_scanner.integrations.persistent_runner import PersistentRunner


@pytest.mark.skipif(sys.platform != "win32" or not shutil.which("wsl.exe"), reason="requires WSL")
def test_bridge_reuses_worker_and_isolates_cases(tmp_path):
    python = "/home/segayi1/miniconda3/envs/spurs_v4/bin/python"
    probe = subprocess.run(["wsl.exe", "-d", "Ubuntu-22.04", "--exec", "test", "-x", python],
                           capture_output=True, timeout=15)
    if probe.returncode:
        pytest.skip("SPURS WSL environment unavailable")
    root = Path(__file__).resolve().parents[1]
    bridge = root / "integrations" / "wsl_bridge.py"
    worker = tmp_path / "spurs_gpu_worker.py"
    worker.write_text("", encoding="utf-8")
    (tmp_path / "spurs_persistent_worker.py").write_text('''
import json, os, sys
if os.getsid(0) != os.getpid():
    os.setsid()
control = sys.argv[sys.argv.index("--control-file") + 1]
with open(control, "w") as f:
    json.dump({"pid": os.getpid(), "start_ticks": open("/proc/self/stat").read().split()[21]}, f)
for line in sys.stdin:
    envelope = json.loads(line)
    request = envelope["request"]
    result = {"verified": True} if request.get("_finalize") else {"pid": os.getpid(), "value": request["value"]}
    print(json.dumps({"request_id": envelope["request_id"], "result": result}), flush=True)
''', encoding="utf-8")
    command = [sys.executable, str(bridge), "--distro", "Ubuntu-22.04", "--python",
               python, "--worker", str(worker), "--persistent"]
    pids = []
    for _ in range(2):
        with PersistentRunner(command, 30, tmp_path / "stderr.log") as runner:
            first = runner.run({"value": 1})
            second = runner.run({"value": 2})
            assert first["pid"] == second["pid"]
            assert (first["value"], second["value"]) == (1, 2)
            assert runner.run({"_finalize": True})["verified"]
            pids.append(first["pid"])
        assert runner.process.poll() is not None
        gone = subprocess.run(["wsl.exe", "-d", "Ubuntu-22.04", "--exec", "kill", "-0", str(pids[-1])],
                              capture_output=True, timeout=15)
        assert gone.returncode != 0
    assert pids[0] != pids[1]
