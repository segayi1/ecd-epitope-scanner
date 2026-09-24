from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def _wsl_path(value: str, distro: str) -> str:
    if not value:
        return value
    completed = subprocess.run(
        ["wsl.exe", "-d", distro, "--", "wslpath", "-a", value],
        capture_output=True,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr.decode("utf-8", errors="replace").strip())
    return completed.stdout.decode("utf-8", errors="strict").strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--distro", default="Ubuntu-22.04")
    parser.add_argument("--python", required=True)
    parser.add_argument("--worker", required=True)
    args = parser.parse_args()

    request = json.load(sys.stdin)
    for key in ("structure_path", "output_dir"):
        if request.get(key):
            request[key] = _wsl_path(str(Path(request[key]).resolve()), args.distro)
    worker = _wsl_path(str(Path(args.worker).resolve()), args.distro)
    completed = subprocess.run(
        ["wsl.exe", "-d", args.distro, "--", args.python, worker],
        input=json.dumps(request, sort_keys=True).encode("utf-8"),
        capture_output=True,
        check=False,
    )
    if completed.stderr:
        sys.stderr.buffer.write(completed.stderr)
    if completed.stdout:
        sys.stdout.buffer.write(completed.stdout)
    return int(completed.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
