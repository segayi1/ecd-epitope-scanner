from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


def _decode_output(value: bytes) -> str:
    # WSL service errors are localized UTF-16LE on Windows and often appear
    # on stdout. Linux programs normally return UTF-8.
    if b"\x00" in value[:256]:
        # A WSL UTF-16 warning may precede the worker's UTF-8 traceback in
        # the same pipe. Decode each prefix line without corrupting the tail.
        end = value.find(b"\n\x00")
        if end >= 0:
            prefix = value[:end + 2].decode("utf-16-le", errors="replace").strip()
            tail = _decode_output(value[end + 2:])
            return prefix + ("\n" + tail if tail else "")
        return value.decode("utf-16-le", errors="replace").strip()
    return value.decode("utf-8", errors="replace").strip()


class WSLBridgeError(RuntimeError):
    def __init__(self, operation: str, distro: str, completed):
        stdout = _decode_output(completed.stdout or b"")
        stderr = _decode_output(completed.stderr or b"")
        self.code = "WSL_ACCESS_DENIED" if "E_ACCESSDENIED" in stdout + stderr else "WSL_EXECUTION_FAILED"
        self.details = {
            "operation": operation, "distro": distro,
            "returncode": completed.returncode,
            "stdout": stdout[-4000:], "stderr": stderr[-4000:],
        }
        message = stdout or stderr or f"exit code {completed.returncode}"
        if self.code == "WSL_ACCESS_DENIED":
            message += (
                " WSL rejected the Windows process identity. Check the launching "
                "Python executable's Windows integrity label and run the application "
                "with a normal user token; low-integrity Python cannot access WSL."
            )
        super().__init__(f"WSL {operation} ({distro}) failed: {message}")


def _wsl_path(value: str, distro: str) -> str:
    if not value:
        return value
    completed = subprocess.run(
        ["wsl.exe", "-d", distro, "--exec", "wslpath", "-a", value],
        capture_output=True,
        check=False,
        timeout=30,
    )
    if completed.returncode:
        raise WSLBridgeError("path conversion", distro, completed)
    path = completed.stdout.decode("utf-8", errors="strict").strip()
    if not path.startswith("/"):
        raise RuntimeError(f"WSL returned an invalid Linux path: {path!r}")
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--distro", default="Ubuntu-22.04")
    parser.add_argument("--python", required=True)
    parser.add_argument("--worker", required=True)
    args = parser.parse_args()

    started = time.perf_counter()
    try:
        request = json.load(sys.stdin)
        for key in ("structure_path", "output_dir"):
            if request.get(key):
                request[key] = _wsl_path(str(Path(request[key]).resolve()), args.distro)
        if request.get("parameters", {}).get("model_assets_root"):
            request["parameters"]["model_assets_root"] = _wsl_path(
                str(Path(request["parameters"]["model_assets_root"]).resolve()), args.distro)
        worker = _wsl_path(str(Path(args.worker).resolve()), args.distro)
        completed = subprocess.run(
            ["wsl.exe", "-d", args.distro, "--exec", args.python, worker],
            input=json.dumps(request, sort_keys=True).encode("utf-8"),
            capture_output=True,
            check=False,
        )
        if completed.returncode:
            raise WSLBridgeError("worker execution", args.distro, completed)
        if completed.stderr:
            sys.stderr.buffer.write((_decode_output(completed.stderr) + "\n").encode("utf-8"))
        sys.stdout.buffer.write(completed.stdout)
    except (OSError, RuntimeError, subprocess.TimeoutExpired, ValueError) as exc:
        print(json.dumps({
            "status": "failed", "records": [], "warnings": [],
            "errors": [{
                "code": getattr(exc, "code", "WSL_ENVIRONMENT_INVALID"),
                "message": str(exc), "details": getattr(exc, "details", {}),
            }],
            "provenance": {"bridge": "wsl_bridge", "distro": args.distro},
            "runtime_seconds": time.perf_counter() - started,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
