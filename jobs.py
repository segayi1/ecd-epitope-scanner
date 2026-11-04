from __future__ import annotations

import datetime as dt
import json
import os
import shutil
import sqlite3
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from .models import TERMINAL_STATUSES
from .pipeline import PipelineRequest, run_pipeline
from .util import ClosingSQLiteConnection, WaitingForAF3


def _now() -> str:
    return dt.datetime.now(dt.UTC).isoformat()


def _process_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        open_process = kernel32.OpenProcess
        open_process.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        open_process.restype = wintypes.HANDLE
        wait_for_single_object = kernel32.WaitForSingleObject
        wait_for_single_object.argtypes = [wintypes.HANDLE, wintypes.DWORD]
        wait_for_single_object.restype = wintypes.DWORD
        close_handle = kernel32.CloseHandle
        close_handle.argtypes = [wintypes.HANDLE]
        close_handle.restype = wintypes.BOOL

        # os.kill(pid, 0) is not a reliable existence check on Windows and can
        # surface ERROR_INVALID_PARAMETER as a SystemError during restart scans.
        handle = open_process(0x00100000, False, pid)  # SYNCHRONIZE
        if not handle:
            error = ctypes.get_last_error()
            if error == 87:  # ERROR_INVALID_PARAMETER: no process has this PID
                return False
            # Access denied or an unfamiliar error is treated as possibly alive
            # so recovery never interrupts a task based on uncertain evidence.
            return True
        try:
            return wait_for_single_object(handle, 0) == 0x102  # WAIT_TIMEOUT
        finally:
            close_handle(handle)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except (OSError, SystemError):
        return False
    return True


class JobStore:
    def __init__(self, root: Path):
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        self.db_path = self.root / "jobs.sqlite3"
        self._lock = threading.RLock()
        self._timing_cache = {}
        self._initialize()
        self._summary_migration = threading.Thread(target=self._backfill_summaries, daemon=True)
        self._summary_migration.start()

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(
            self.db_path, timeout=30, factory=ClosingSQLiteConnection
        )
        connection.row_factory = sqlite3.Row
        return connection

    def _initialize(self) -> None:
        with self._connect() as connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS jobs (
                    id TEXT PRIMARY KEY,
                    owner_id TEXT NOT NULL DEFAULT 'local',
                    owner_pid INTEGER,
                    status TEXT NOT NULL,
                    stage TEXT NOT NULL,
                    progress REAL NOT NULL,
                    message TEXT NOT NULL,
                    request_json TEXT NOT NULL,
                    result_json TEXT,
                    summary_json TEXT,
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
            if "owner_pid" not in columns:
                connection.execute("ALTER TABLE jobs ADD COLUMN owner_pid INTEGER")
            if "summary_json" not in columns:
                connection.execute("ALTER TABLE jobs ADD COLUMN summary_json TEXT")
            # Jobs left QUEUED/RUNNING by a previous process can never finish:
            # the single worker thread died with the old interpreter.  Mark them
            # failed so SSE streams terminate and the UI does not spin forever.
            for row in connection.execute("SELECT id,owner_pid FROM jobs WHERE status IN ('QUEUED', 'RUNNING')"):
                if row["owner_pid"] is not None and _process_alive(row["owner_pid"]):
                    continue
                connection.execute(
                    """UPDATE jobs SET status='FAILED', stage='FAILED',
                       error='interrupted by service restart',
                       message='Pipeline was interrupted by a service restart', updated_at=? WHERE id=?""",
                    (_now(), row["id"]),
                )
            connection.commit()

    def _backfill_summaries(self):
        """Migrate old rows off the request path, committing one compact row at a time."""
        with self._connect() as connection:
            ids = [row[0] for row in connection.execute(
                "SELECT id FROM jobs WHERE summary_json IS NULL AND result_json IS NOT NULL")]
        for job_id in ids:
            with self._connect() as connection:
                # Multiple JSON paths share a single parse inside SQLite; no blob enters Python.
                try:
                    row = connection.execute("""SELECT json_extract(result_json,
                        '$.protein.accession', '$.summary.regions', '$.summary.recommended_constructs')
                        FROM jobs WHERE id=? AND summary_json IS NULL""", (job_id,)).fetchone()
                except sqlite3.OperationalError as exc:
                    if "malformed JSON" in str(exc):
                        continue
                    raise
                if row is None or row[0] is None:
                    continue
                accession, regions, constructs = json.loads(row[0])
                summary = {"protein": {"accession": accession}, "summary": {
                    "regions": regions, "recommended_constructs": constructs}}
                with self._lock:
                    connection.execute("UPDATE jobs SET summary_json=? WHERE id=? AND summary_json IS NULL",
                                       (json.dumps(summary, ensure_ascii=False), job_id))
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
                (id,owner_id,owner_pid,status,stage,progress,message,request_json,result_json,error,workdir,created_at,updated_at)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
                """,
                (
                    job_id,
                    owner_id,
                    os.getpid(),
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
        if "result_json" in payload:
            result = payload["result_json"]
            if isinstance(result, str):
                try:
                    result = json.loads(result)
                except (TypeError, ValueError):
                    result = None
            result = result if isinstance(result, dict) else {}
            payload["summary_json"] = json.dumps({
                "protein": {"accession": (result.get("protein") or {}).get("accession")},
                "summary": {
                    "regions": (result.get("summary") or {}).get("regions"),
                    "recommended_constructs": (result.get("summary") or {}).get(
                        "recommended_constructs"
                    ),
                },
            }, ensure_ascii=False)
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
        data.pop("summary_json", None)
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

    def summaries(self, limit: int = 30, owner_id: str = "local") -> list[dict[str, Any]]:
        """Return list-page fields without parsing request or result JSON."""
        with self._connect() as connection:
            rows = connection.execute(
                """SELECT id,status,stage,progress,message,summary_json,error,created_at,updated_at
                   FROM jobs WHERE owner_id=? ORDER BY created_at DESC LIMIT ?""",
                (owner_id, limit),
            ).fetchall()
        result = []
        for row in rows:
            data = dict(row)
            raw = data.pop("summary_json")
            data["summary"] = json.loads(raw) if raw else {}
            data["terminal"] = data["status"] in TERMINAL_STATUSES
            result.append(data)
        return result

    def request_for(self, job_id: str, owner_id: str = "local") -> PipelineRequest:
        return PipelineRequest(**self.get(job_id, owner_id)["request"])

    def time_estimate(self, job: dict[str, Any]) -> dict[str, Any] | None:
        if job["status"] != "RUNNING":
            return None
        from .job_timing import read_timing, remaining_time, timing_key
        key = (job.get("owner_id", "local"), timing_key(job["request"]))
        cached_at, references = self._timing_cache.get(key, (0., []))
        if time.monotonic() - cached_at >= 60 or key not in self._timing_cache:
            with self._connect() as connection:
                rows = connection.execute(
                    "SELECT request_json,workdir,created_at,updated_at FROM jobs WHERE owner_id=? AND status='COMPLETE' ORDER BY created_at DESC LIMIT 80",
                    (key[0],),
                ).fetchall()
            references = []
            for row in rows:
                try:
                    if timing_key(json.loads(row["request_json"])) != key[1]:
                        continue
                    reference = read_timing(row["workdir"])
                    if reference:
                        if dt.datetime.fromisoformat(row["created_at"]) > reference["started_at"]:
                            continue  # A published replay may contain its source run's older trace.
                        finished = dt.datetime.fromisoformat(row["updated_at"])
                        reference["total_seconds"] = max(reference["elapsed_seconds"],
                            (finished - reference["started_at"]).total_seconds())
                        references.append(reference)
                except (ValueError, TypeError, KeyError):
                    continue
            self._timing_cache[key] = (time.monotonic(), references)
        current = read_timing(job["workdir"])
        if job.get("message") == "Starting pipeline" or (
                job["stage"] == "FETCHING" and current and current["stage"] not in {None, "STARTING", "FETCHING"}):
            current = None  # A resumed worker has not replaced the old execution trace yet.
        return remaining_time(current, references, dt.datetime.now(dt.UTC))

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
            job = self.store.get(job_id)
            request.metadata = {**request.metadata, "execution_queue": {
                "job_id": job_id, "created_at": job["created_at"],
                "worker_started_at": dt.datetime.now(dt.UTC).isoformat()}}
            # Keep raw model paths stable; snapshot metadata before a resumed attempt.
            workdir = Path(job["workdir"]).resolve()
            engine_runs = workdir / "engine_runs"
            if engine_runs.exists():
                history = workdir / "prior_attempts" / uuid.uuid4().hex
                if not engine_runs.resolve().is_relative_to(workdir) or not history.resolve().is_relative_to(workdir):
                    raise ValueError("Engine-run history must remain within the task directory")
                history.mkdir(parents=True)
                sources = [p for p in workdir.iterdir() if p.is_file()]
                sources += [p for p in engine_runs.rglob("*") if p.is_file() and p.suffix in {".json", ".jsonl", ".log", ".txt"}]
                for source in sources:
                    target = history / source.relative_to(workdir)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, target)

            def progress(stage: str, value: float, message: str) -> None:
                # Publish terminal state only with the result below.
                status = "RUNNING"
                self.store.update(
                    job_id,
                    status=status,
                    stage=stage,
                    progress=max(0.0, min(1.0, float(value))),
                    message=message,
                )

            result = run_pipeline(request, progress)
            design_summary = result.get("design_summary") or {}
            recommendation_status = design_summary.get("recommendation_status")
            completion_message = "Analysis complete"
            if recommendation_status == "BLOCKED_EXTERNAL_EVIDENCE":
                outcome = (design_summary.get("recommendation_diagnostics") or {}).get(
                    "outcome", "FORMAL_RECOMMENDATIONS_BLOCKED"
                )
                completion_message = f"Analysis complete; no formal panel: {outcome}. See recommendation diagnostics."
            self.store.update(
                job_id,
                status="COMPLETE",
                stage="COMPLETE",
                progress=1.0,
                message=completion_message,
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

    def resume_with_af3_manifest(self, job_id: str, manifest_path: Path) -> None:
        request = self.store.request_for(job_id)
        if not request.accession and request.input_type == "uniprot":
            request.accession = request.identifier
        request.input_type = "af3_upload"
        request.af3_manifest_path = str(manifest_path)
        request.metadata["official_af3_source_confirmed"] = True
        request.structure_path = None
        request.confidence_path = None
        request.allow_partial = False
        request.output_dir = self.store.get(job_id)["workdir"]
        self.store.update(
            job_id,
            status="QUEUED",
            stage="QUEUED",
            progress=0.0,
            message="AF3 folder uploaded; job queued for resume",
            request_json=request.as_dict() | {"original_cds": request.original_cds},
            result_json=None,
            error=None,
        )
        self.submit(job_id)
