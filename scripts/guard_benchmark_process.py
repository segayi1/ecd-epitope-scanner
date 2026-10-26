from __future__ import annotations

"""Hold a scoped Windows execution request only while one benchmark PID runs."""
import argparse
import ctypes
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
from ecd_epitope_scanner.surface_perturbation.power_guard import prevent_sleep  # noqa: E402


def alive(pid: int) -> bool:
    if os.name != "nt":
        try:
            os.kill(pid, 0)
            return True
        except OSError:
            return False
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.OpenProcess.argtypes = [ctypes.c_uint32, ctypes.c_int, ctypes.c_uint32]
    kernel.OpenProcess.restype = ctypes.c_void_p
    kernel.GetExitCodeProcess.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint32)]
    kernel.GetExitCodeProcess.restype = ctypes.c_int
    kernel.CloseHandle.argtypes = [ctypes.c_void_p]
    kernel.CloseHandle.restype = ctypes.c_int
    handle = kernel.OpenProcess(0x1000, 0, pid)
    if not handle:
        return False
    try:
        code = ctypes.c_uint32()
        return bool(kernel.GetExitCodeProcess(handle, ctypes.byref(code))) and code.value == 259
    finally:
        kernel.CloseHandle(handle)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pid", type=int, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--receipt", type=Path, required=True)
    args = p.parse_args()
    manifest = args.manifest.resolve()
    receipt = args.receipt.resolve()
    request_started = time.time()
    report = {"scope": "benchmark runner PID lifetime only", "target_pid": args.pid,
              "manifest": str(manifest), "request_started_epoch": request_started}
    terminal = {"COMPLETE", "FAILED", "PAUSED", "CANCELLED"}
    with prevent_sleep() as request:
        report["power_request"] = request
        report["status"] = "ACTIVE"
        receipt.parent.mkdir(parents=True, exist_ok=True)
        receipt.write_text(json.dumps(report, indent=2), encoding="utf-8")
        while alive(args.pid):
            try:
                state = json.loads(manifest.read_text(encoding="utf-8")).get("status")
            except (OSError, ValueError):
                state = None
            if state in terminal:
                break
            time.sleep(5)
        report["status"] = "RELEASED_RUNNER_EXITED_OR_TERMINAL"
        report["request_released_epoch"] = time.time()
        report["active_seconds"] = report["request_released_epoch"] - request_started
        receipt.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
