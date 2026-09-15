from __future__ import annotations

import datetime as dt
import json
import sqlite3
import shutil
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from .models import TERMINAL_STATUSES
from .pipeline import PipelineRequest, run_pipeline
from .util import WaitingForAF3


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


class JobStore:
    def __init__(self, root: Path):
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        self.db_path = self.root / "jobs.sqlite3"
        self._lock = threading.RLock()
        self._initialize()

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.db_path, timeout=30)
        connection.row_factory = sqlite3.Row
        return connection

    def _initialize(self) -> None:
        with self._connect() as connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS jobs (
                    id TEXT PRIMARY KEY,
                    owner_id TEXT NOT NULL DEFAULT 'local',
                    status TEXT NOT NULL,
                    stage TEXT NOT NULL,
                    progress REAL NOT NULL,
                    message TEXT NOT NULL,
                    request_json TEXT NOT NULL,
                    result_json TEXT,
                    error TEXT,
                    workdir TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                )
                """
            )
            # SQLite has no portable ADD COLUMN IF NOT EXISTS.  Migrate
            # workspaces created before job ownership was introduced.
            columns = {row[1] for row in connection.execute("PRAGMA table_info(jobs)")}
            if "owner_id" not in columns:
                connection.execute("ALTER TABLE jobs ADD COLUMN owner_id TEXT NOT NULL DEFAULT 'local'")
            # Jobs left QUEUED/RUNNING by a previous process can never finish:
            # the single worker thread died with the old interpreter.  Mark them
            # failed so SSE streams terminate and the UI does not spin forever.
            connection.execute(
                """
                UPDATE jobs
                SET status='FAILED', stage='FAILED', error='interrupted by service restart',
                    message='Pipeline was interrupted by a service restart', updated_at=?
                WHERE status IN ('QUEUED', 'RUNNING')
                """,
                (_now(),),
            )
            connection.commit()

    def create(self, request: PipelineRequest, job_id: str | None = None, owner_id: str = "local") -> dict[str, Any]:
        job_id = job_id or uuid.uuid4().hex[:12]
        workdir = self.root / job_id
        workdir.mkdir(parents=True, exist_ok=True)
        request.output_dir = str(workdir)
        now = _now()
        request_data = request.as_dict() | {"original_cds": request.original_cds}
        with self._lock, self._connect() as connection:
            connection.execute(
                """
                INSERT INTO jobs
                (id,owner_id,status,stage,progress,message,request_json,result_json,error,workdir,created_at,updated_at)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
                """,
                (
                    job_id,
                    owner_id,
                    "QUEUED",
                    "QUEUED",
                    0.0,
                    "Job queued",
                    json.dumps(request_data, ensure_ascii=False),
                    None,
                    None,
                    str(workdir),
                    now,
                    now,
                ),
            )
            connection.commit()
        return self.get(job_id)

    def update(self, job_id: str, **values: Any) -> dict[str, Any]:
        allowed = {
            "status",
            "stage",
            "progress",
            "message",
            "request_json",
            "result_json",
            "error",
        }
        payload = {key: value for key, value in values.items() if key in allowed}
        if not payload:
            return self.get(job_id)
        payload["updated_at"] = _now()
        assignments = ", ".join(f"{key}=?" for key in payload)
        parameters = [
            json.dumps(value, ensure_ascii=False)
            if key in {"request_json", "result_json"} and not isinstance(value, str)
            else value
            for key, value in payload.items()
        ]
        with self._lock, self._connect() as connection:
            cursor = connection.execute(
                f"UPDATE jobs SET {assignments} WHERE id=?", parameters + [job_id]
            )
            if not cursor.rowcount:
                raise KeyError(job_id)
            connection.commit()
        return self.get(job_id)

    def get(self, job_id: str, owner_id: str = "local") -> dict[str, Any]:
        with self._connect() as connection:
            row = connection.execute("SELECT * FROM jobs WHERE id=? AND owner_id=?", (job_id, owner_id)).fetchone()
        if row is None:
            raise KeyError(job_id)
        data = dict(row)
        for key in ("request_json", "result_json"):
            raw = data.pop(key)
            data[key.removesuffix("_json")] = json.loads(raw) if raw else None
        data["terminal"] = data["status"] in TERMINAL_STATUSES
        return data

    def list(self, limit: int = 30, owner_id: str = "local") -> list[dict[str, Any]]:
        with self._connect() as connection:
            rows = connection.execute(
                "SELECT id FROM jobs WHERE owner_id=? ORDER BY created_at DESC LIMIT ?", (owner_id, limit)
            ).fetchall()
        return [self.get(row["id"], owner_id) for row in rows]

    def request_for(self, job_id: str, owner_id: str = "local") -> PipelineRequest:
        return PipelineRequest(**self.get(job_id, owner_id)["request"])

    def cleanup_preview(self, include_test_output: bool = False) -> list[dict[str, Any]]:
        """Return removable successful-run diagnostics without changing data."""
        candidates: list[dict[str, Any]] = []
        for job in self.list(limit=10_000):
            if job["status"] != "COMPLETE":
                continue
            path = Path(job["workdir"]) / "engine_runs"
            if path.exists():
                size = sum(item.stat().st_size for item in path.rglob("*") if item.is_file())
                candidates.append({"job_id": job["id"], "path": str(path), "bytes": size})
        if include_test_output:
            test_output = self.root.parent / "test_output"
            if test_output.exists():
                size = sum(item.stat().st_size for item in test_output.rglob("*") if item.is_file())
                candidates.append({"job_id": "test_output", "path": str(test_output), "bytes": size})
        return candidates

    def cleanup_successful_engine_runs(self) -> int:
        removed = 0
        for item in self.cleanup_preview():
            path = Path(item["path"])
            # This path comes from a known job workdir, not user input.
            shutil.rmtree(path)
            removed += int(item["bytes"])
        return removed


class JobRunner:
    def __init__(self, store: JobStore):
        self.store = store
        self.executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="ecd-pipeline")
        self._active: set[str] = set()
        self._lock = threading.Lock()

    def submit(self, job_id: str) -> None:
        with self._lock:
            if job_id in self._active:
                return
            self._active.add(job_id)
        self.executor.submit(self._run, job_id)

    def _run(self, job_id: str) -> None:
        try:
            request = self.store.request_for(job_id)
            self.store.update(
                job_id,
                status="RUNNING",
                stage="FETCHING",
                progress=0.01,
                message="Starting pipeline",
                error=None,
            )
            # Final artefacts and source data remain; FoldX diagnostics are
            # reproducible and are discarded only after a successful export.
            engine_runs = Path(self.store.get(job_id)["workdir"]) / "engine_runs"
            if engine_runs.exists():
                shutil.rmtree(engine_runs)

            def progress(stage: str, value: float, message: str) -> None:
                status = "COMPLETE" if stage == "COMPLETE" else "RUNNING"
                self.store.update(
                    job_id,
                    status=status,
                    stage=stage,
                    progress=max(0.0, min(1.0, float(value))),
                    message=message,
                )

            result = run_pipeline(request, progress)
            self.store.update(
                job_id,
                status="COMPLETE",
                stage="COMPLETE",
                progress=1.0,
                message="Analysis complete",
                result_json=result,
                error=None,
            )
        except WaitingForAF3 as exc:
            self.store.update(
                job_id,
                status="WAITING_FOR_AF3_UPLOAD",
                stage="WAITING_FOR_AF3_UPLOAD",
                progress=0.34,
                message=str(exc),
                result_json={"candidates": exc.candidates},
                error=None,
            )
        except Exception as exc:
            self.store.update(
                job_id,
                status="FAILED",
                stage="FAILED",
                message="Pipeline failed",
                error=f"{type(exc).__name__}: {exc}",
            )
        finally:
            with self._lock:
                self._active.discard(job_id)

    def resume_with_af3(self, job_id: str, structure_path: Path, confidence_path: Path) -> None:
        request = self.store.request_for(job_id)
        # P0 fix: a UniProt job stores its accession in `identifier` (the API
        # form does not populate `accession`).  Without this backfill the
        # af3_upload branch falls into _synthetic_context and fails with
        # "requires sequence, ecd_start and ecd_end" -- every AF3 resume died.
        if not request.accession and request.input_type == "uniprot":
            request.accession = request.identifier
        request.input_type = "af3_upload"
        request.structure_path = str(structure_path)
        request.confidence_path = str(confidence_path)
        request.allow_partial = False
        request.output_dir = self.store.get(job_id)["workdir"]
        self.store.update(
            job_id,
            status="QUEUED",
            stage="QUEUED",
            progress=0.0,
            message="Official AF3 files uploaded; job queued for resume",
            request_json=request.as_dict() | {"original_cds": request.original_cds},
            result_json=None,
            error=None,
        )
        self.submit(job_id)
