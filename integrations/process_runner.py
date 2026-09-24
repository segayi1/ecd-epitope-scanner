from __future__ import annotations

import json
import os
import subprocess
import time
from collections.abc import Sequence
from pathlib import Path

from ..joint_design_v4.models import PredictionRequest, PredictionResult
from .base import ToolError


class JsonWorkerRunner:
    """Run an isolated worker using argv and a strict stdout JSON contract."""

    def __init__(self, command: Sequence[str] | None, timeout_seconds: float = 300):
        self.command = list(command or [])
        self.timeout_seconds = float(timeout_seconds)

    def run(self, request: PredictionRequest, *, log_path: Path | None = None) -> PredictionResult:
        if not self.command:
            raise ToolError("TOOL_NOT_CONFIGURED", "No external worker command was configured")
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
