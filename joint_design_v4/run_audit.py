"""Small run-local provenance and elapsed-time trace helpers."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

def source_snapshot(root: Path, effective_parameters: dict[str, Any] | None = None,
                    archive_path: Path | None = None) -> dict[str, Any]:
    root = root.resolve()
    source_files = list(root.glob("*.py"))
    for directory in ("joint_design_v4", "joint_design", "integrations", "surface_patch",
                      "surface_perturbation", "vendor/spurs"):
        source_files.extend((root / directory).rglob("*.py"))
    for suffix in ("*.html", "*.js", "*.css"):
        source_files.extend((root / "web").rglob(suffix))
    lock_names = {"pyproject.toml", "requirements.txt", "requirements.lock", "uv.lock", "poetry.lock",
                  "Pipfile.lock", "environment.yml", "environment.yaml"}
    source_files.extend(path for path in root.iterdir() if path.is_file() and path.name in lock_names)
    rows, contents = [], {}
    for path in sorted(set(source_files)):
        if not path.is_file() or path.is_symlink():
            continue
        relative = path.relative_to(root).as_posix()
        content = path.read_bytes()
        contents[relative] = content
        rows.append({"path": relative, "sha256": hashlib.sha256(content).hexdigest(),
                     "size_bytes": len(content)})
    git = {"revision": None, "dirty": None, "dirty_paths": [], "status": "UNAVAILABLE"}
    try:
        revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True,
                                  text=True, timeout=5, check=False)
        status = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"], cwd=root,
                                capture_output=True, text=True, timeout=5, check=False)
        if revision.returncode == 0 and status.returncode == 0:
            git = {"revision": revision.stdout.strip(), "dirty": bool(status.stdout.strip()),
                   "dirty_paths": [line[3:] for line in status.stdout.splitlines() if len(line) >= 4],
                   "status": "CAPTURED"}
    except (OSError, subprocess.TimeoutExpired):
        pass
    parameters = effective_parameters or {}
    parameters_digest = hashlib.sha256(json.dumps(parameters, sort_keys=True, separators=(",", ":"),
                                         ensure_ascii=False, default=str).encode("utf-8")).hexdigest()
    digest_input = json.dumps(rows, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    archive = None
    if archive_path is not None:
        archive_path = archive_path.resolve()
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = archive_path.with_name(archive_path.name + ".tmp")
        try:
            with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
                for relative, content in sorted(contents.items()):
                    bundle.writestr(relative, content)
            temporary.replace(archive_path)
        finally:
            temporary.unlink(missing_ok=True)
        archive = {"path": archive_path.name, "sha256": hashlib.sha256(archive_path.read_bytes()).hexdigest(),
                   "size_bytes": archive_path.stat().st_size, "member_count": len(contents)}
    snapshot = {
        "schema_version": "run_source_snapshot_v1",
        "captured_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_files": rows,
        "source_files_sha256": hashlib.sha256(digest_input).hexdigest(),
        "source_archive": archive,
        "effective_parameters_sha256": parameters_digest,
        "effective_parameters": parameters,
        "git": git,
        "runtime": {"python": sys.version, "platform": platform.platform()},
        "interpretation": "Source, configuration and dependency identity at run start; missing historical identities are not inferred.",
    }
    snapshot["snapshot_id"] = hashlib.sha256(
        json.dumps(snapshot, sort_keys=True, separators=(",", ":"), ensure_ascii=False, default=str).encode("utf-8")
    ).hexdigest()
    return snapshot


class ExecutionTrace:
    """Append-only UTC plus monotonic elapsed events; elapsed time is not CPU time."""

    def __init__(self, path: Path):
        self.execution_id = __import__("uuid").uuid4().hex
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.start_ns = time.monotonic_ns()
        self.last_ns = self.start_ns
        self.last_stage: str | None = None
        self.stage_started_ns = self.start_ns
        self.stage_elapsed: dict[str, float] = {}
        self.sequence = 0
        self.record("PIPELINE_STARTED", None, None, "Pipeline entered; queue delay is not captured here.")

    def record(self, event: str, stage: str | None, fraction: float | None, message: str,
               **details: Any) -> dict[str, Any]:
        now = time.monotonic_ns()
        if stage and self.last_stage and stage != self.last_stage:
            self.stage_elapsed[self.last_stage] = self.stage_elapsed.get(self.last_stage, 0.0) + (
                now - self.stage_started_ns) / 1e9
        if stage and stage != self.last_stage:
            self.stage_started_ns = now
            self.last_stage = stage
        self.sequence += 1
        row = {"execution_id": self.execution_id, "sequence": self.sequence, "event": event, "stage": stage,
               "progress_fraction": fraction, "message": message,
               "timestamp_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_monotonic_seconds": round((now - self.start_ns) / 1e9, 6),
               "interval_monotonic_seconds": round((now - self.last_ns) / 1e9, 6),
               **details}
        encoded = (json.dumps(row, ensure_ascii=False, sort_keys=True, default=str) + "\n").encode("utf-8")
        descriptor = os.open(self.path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        try:
            os.write(descriptor, encoded)
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        self.last_ns = now
        return row

    def progress(self, stage: str, fraction: float, message: str) -> None:
        self.record("STAGE_PROGRESS", stage, fraction, message)

    def summary(self) -> dict[str, Any]:
        elapsed = (time.monotonic_ns() - self.start_ns) / 1e9
        stage_elapsed = dict(self.stage_elapsed)
        if self.last_stage:
            stage_elapsed[self.last_stage] = stage_elapsed.get(self.last_stage, 0.0) + (
                time.monotonic_ns() - self.stage_started_ns) / 1e9
        return {"trace_file": self.path.name, "elapsed_wall_seconds": round(elapsed, 6),
                "stage_elapsed_wall_seconds": {key: round(value, 6) for key, value in stage_elapsed.items()},
                "queue_delay_captured": bool(getattr(self, "queue_delay_captured", False)),
                "interpretation": "Monotonic elapsed wall time; does not represent CPU core-hours or GPU utilization."}
