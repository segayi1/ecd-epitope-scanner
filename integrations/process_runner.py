from __future__ import annotations

import json
import os
import subprocess
import time
from contextlib import contextmanager
from contextvars import ContextVar
from collections.abc import Sequence
from pathlib import Path

from ..joint_design_v4.models import PredictionRequest, PredictionResult
from .base import ToolError

_PERSISTENT = ContextVar('spurs_persistent_runners', default=None)


def _managed_runner(command, timeout_seconds, log_path=None):
    runners = _PERSISTENT.get()
    if (runners is None or len(command) < 2 or Path(command[1]).name != 'wsl_bridge.py'
            or not any(Path(arg).name == 'spurs_gpu_worker.py' for arg in command)):
        return None
    from .persistent_runner import PersistentRunner
    key = tuple(command)
    if key not in runners:
        runners[key] = PersistentRunner([*command, '--persistent'], timeout_seconds, log_path)
    return runners[key]


def managed_spurs_environment(adapter, request):
    if _PERSISTENT.get() is None:
        return None
    runner = _managed_runner(adapter.command, adapter.timeout_seconds)
    if runner is None:
        return None
    return runner.run({'_environment': True, 'parameters': request.parameters})['environment'] | {'_persistent_worker': True}


@contextmanager
def persistent_spurs_scope():
    runners = {}
    token = _PERSISTENT.set(runners)
    try:
        yield
        for runner in runners.values():
            response = runner.run({'_finalize': True})
            if not response.get('verified'):
                raise ToolError('TOOL_CRASHED', 'Persistent worker final identity verification failed')
    finally:
        try:
            for runner in runners.values():
                runner.close()
        finally:
            _PERSISTENT.reset(token)


class JsonWorkerRunner:
    """Run an isolated worker using argv and a strict stdout JSON contract."""

    def __init__(self, command: Sequence[str] | None, timeout_seconds: float = 300):
        self.command = list(command or [])
        self.timeout_seconds = float(timeout_seconds)

    def run(self, request: PredictionRequest, *, log_path: Path | None = None) -> PredictionResult:
        if not self.command:
            raise ToolError("TOOL_NOT_CONFIGURED", "No external worker command was configured")
        runners = _PERSISTENT.get()
        if (runners is not None and len(self.command) >= 2
                and Path(self.command[1]).name == 'wsl_bridge.py'
                and any(Path(arg).name == 'spurs_gpu_worker.py' for arg in self.command)):
            runner = _managed_runner(self.command, self.timeout_seconds, log_path)
            payload = runner.run(request.as_dict())
            if payload.get('status') not in {'success', 'partial', 'failed'}:
                raise ToolError('OUTPUT_SCHEMA_INVALID', 'Worker JSON lacks a valid status')
            return PredictionResult(status=payload['status'], records=list(payload.get('records', [])),
                warnings=list(payload.get('warnings', [])), errors=list(payload.get('errors', [])),
                provenance=dict(payload.get('provenance', {})) | {'persistent_worker': True},
                runtime_seconds=float(payload.get('runtime_seconds', 0)))
        started = time.perf_counter()
        try:
            completed = subprocess.run(
                self.command,
                input=json.dumps(request.as_dict(), sort_keys=True),
                text=True,
                encoding="utf-8",
                errors="replace",
                capture_output=True,
                timeout=self.timeout_seconds,
                check=False,
                env=dict(os.environ),
            )
        except subprocess.TimeoutExpired as exc:
            raise ToolError("TOOL_TIMEOUT", f"Worker exceeded {self.timeout_seconds:g}s") from exc
        except OSError as exc:
            raise ToolError("TOOL_ENVIRONMENT_INVALID", str(exc)) from exc
        runtime = time.perf_counter() - started
        if log_path:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(completed.stderr or "", encoding="utf-8")
        if completed.returncode:
            raise ToolError(
                "TOOL_CRASHED",
                f"Worker exited with code {completed.returncode}",
                details={"stderr": (completed.stderr or "")[-4000:]},
            )
        try:
            payload = json.loads(completed.stdout)
        except json.JSONDecodeError as exc:
            raise ToolError(
                "OUTPUT_SCHEMA_INVALID",
                "Worker stdout was not a single JSON document",
                details={"stdout": completed.stdout[-4000:]},
            ) from exc
        if not isinstance(payload, dict) or payload.get("status") not in {
            "success",
            "partial",
            "failed",
        }:
            raise ToolError("OUTPUT_SCHEMA_INVALID", "Worker JSON lacks a valid status")
        return PredictionResult(
            status=payload["status"],
            records=list(payload.get("records", [])),
            warnings=list(payload.get("warnings", [])),
            errors=list(payload.get("errors", [])),
            provenance=dict(payload.get("provenance", {})),
            runtime_seconds=float(payload.get("runtime_seconds", runtime)),
        )
