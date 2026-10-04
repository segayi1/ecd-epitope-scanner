"""Task orchestration for the isolated full-model Paraplume worker."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import threading
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from ..experiment_projects import ExperimentStore
from ..util import read_json

ROOT = Path(__file__).resolve().parents[1]
MODEL_ENV = ROOT / ".tool_envs" / "paraplume"
WORKER = ROOT / "integrations" / "workers" / "paraplume_worker.py"
MODEL_PYTHON = MODEL_ENV / "Scripts" / "python.exe"
MODEL_CACHE = ROOT / ".tool_cache" / "hf" / "hub"


class ParaplumeConfigurationError(RuntimeError):
    code = "PARAPLUME_MODEL_UNAVAILABLE"


class ParaplumeTaskError(RuntimeError):
    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def model_fingerprint() -> str:
    if not MODEL_PYTHON.is_file() or not WORKER.is_file():
        raise ParaplumeConfigurationError("The isolated Paraplume environment or worker is missing")
    release = read_json(ROOT / "tooling" / "manifest.json")
    tool = next((item for item in release["tools"] if item.get("name") == "Paraplume"), None)
    if not tool or tool.get("version") != "2.0.0":
        raise ParaplumeConfigurationError("Paraplume 2.0.0 is not registered in the tool manifest")
    weights = read_json(ROOT / ".tool_cache" / "paraplume_weight_manifest.json")
    repositories = weights.get("repositories") or {}
    if not repositories or any(not item.get("files") for item in repositories.values()):
        raise ParaplumeConfigurationError("One or more full-model weights are not verified in cache")
    hub_root = ROOT / ".tool_cache" / "hf" / "hub"
    for repository, entry in repositories.items():
        local_name = "models--" + repository.replace("/", "--")
        repository_root = hub_root / local_name
        ref = repository_root / "refs" / "main"
        if not ref.is_file():
            raise ParaplumeConfigurationError(f"The cached revision for {repository} is missing")
        snapshot = repository_root / "snapshots" / ref.read_text(encoding="utf-8").strip()
        if not snapshot.is_dir():
            raise ParaplumeConfigurationError(f"The cached snapshot for {repository} is missing")
        for asset in entry["files"]:
            path = snapshot / asset["file"]
            if not path.is_file() or path.stat().st_size != int(asset["size"]):
                raise ParaplumeConfigurationError(
                    f"Full-model weight is missing or has the wrong size: {repository}/{asset['file']}"
                )
    checkpoint = (
        MODEL_ENV
        / "Lib"
        / "site-packages"
        / "paraplume"
        / "trained_models"
        / "large_1_0"
        / "checkpoint.pt"
    )
    summary = checkpoint.with_name("summary_dict.json")
    if not checkpoint.is_file() or not summary.is_file():
        raise ParaplumeConfigurationError("The installed full-model checkpoint is missing")
    digest = hashlib.sha256()
    digest.update(checkpoint.read_bytes())
    identity = {
        "version": tool["version"],
        "wheel_sha256": tool["wheel"]["sha256"],
        "model_checkpoint_sha256": digest.hexdigest(),
        "model_summary_sha256": hashlib.sha256(summary.read_bytes()).hexdigest(),
        "repositories": {
            name: {
                "commit": value["commit"],
                "files": [
                    {key: file[key] for key in ("file", "size", "sha256_or_blob")}
                    for file in value["files"]
                ],
            }
            for name, value in sorted(repositories.items())
        },
    }
    return hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _run_worker(
    constructs: list[dict[str, str]], fingerprint: str, progress_callback=None
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    request = {
        "schema_version": "1.0",
        "model_size": "large",
        "model_fingerprint": fingerprint,
        "constructs": constructs,
    }
    env = dict(os.environ)
    env.update(
        {
            "HF_HUB_CACHE": str(MODEL_CACHE),
            "HF_HUB_OFFLINE": "1",
            "TRANSFORMERS_OFFLINE": "1",
            "PYTHONIOENCODING": "utf-8",
        }
    )
    try:
        completed = subprocess.Popen(
            [str(MODEL_PYTHON), str(WORKER)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            env=env,
            cwd=ROOT,
        )
    except OSError as exc:
        raise ParaplumeTaskError("PARAPLUME_WORKER_UNAVAILABLE", str(exc)) from exc

    stderr_tail: deque[str] = deque(maxlen=250)
    result_payloads: list[dict[str, Any]] = []
    progress_errors: list[Exception] = []

    def read_stdout() -> None:
        assert completed.stdout is not None
        for line in completed.stdout:
            try:
                message = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not isinstance(message, dict):
                continue
            if message.get("event") == "progress":
                if progress_callback is not None:
                    try:
                        progress = float(message.get("progress"))
                        detail = str(message.get("message") or "Model inference is running")[:500]
                        if not 0 <= progress <= 1:
                            raise ValueError("Worker progress must be between 0 and 1")
                        progress_callback(progress, detail)
                    except Exception as exc:
                        progress_errors.append(exc)
                continue
            result_payloads.append(message)

    def read_stderr() -> None:
        assert completed.stderr is not None
        for line in completed.stderr:
            stderr_tail.append(line.rstrip())

    stdout_reader = threading.Thread(target=read_stdout, name="paraplume-stdout", daemon=True)
    stderr_reader = threading.Thread(target=read_stderr, name="paraplume-stderr", daemon=True)
    stdout_reader.start()
    stderr_reader.start()
    try:
        assert completed.stdin is not None
        completed.stdin.write(json.dumps(request, separators=(",", ":")))
        completed.stdin.close()
        completed.wait(timeout=1800)
    except subprocess.TimeoutExpired as exc:
        completed.kill()
        completed.wait()
        stdout_reader.join()
        stderr_reader.join()
        raise ParaplumeTaskError("PARAPLUME_TIMEOUT", "Full-model inference exceeded 30 minutes") from exc
    except OSError:
        # A worker that exits during input delivery reports its structured failure on stdout.
        completed.wait()
    stdout_reader.join()
    stderr_reader.join()
    if progress_errors:
        raise ParaplumeTaskError(
            "PARAPLUME_PROGRESS_FAILED", f"Could not save model progress: {progress_errors[0]}"
        )
    stderr = "\n".join(stderr_tail)[-3000:]
    if not result_payloads:
        raise ParaplumeTaskError(
            "PARAPLUME_OUTPUT_INVALID",
            f"Worker returned no result (exit {completed.returncode}): {stderr}",
        )
    payload = result_payloads[-1]
    if completed.returncode or payload.get("status") != "success":
        errors = payload.get("errors") or []
        message = errors[0].get("message", "Full-model inference failed") if errors else "Full-model inference failed"
        raise ParaplumeTaskError(
            errors[0].get("code", "PARAPLUME_INFERENCE_FAILED") if errors else "PARAPLUME_INFERENCE_FAILED",
            f"{message}\n{stderr}".strip(),
        )
    provenance = payload.get("provenance")
    records = payload.get("records")
    if not isinstance(provenance, dict) or provenance.get("model_version") != "2.0.0":
        raise ParaplumeTaskError("PARAPLUME_OUTPUT_INVALID", "Worker provenance is incomplete")
    if provenance.get("model_fingerprint") != fingerprint:
        raise ParaplumeTaskError("PARAPLUME_OUTPUT_INVALID", "Worker model fingerprint differs from task")
    if not isinstance(records, list) or len(records) != len(constructs):
        raise ParaplumeTaskError("PARAPLUME_OUTPUT_INVALID", "Worker returned an unexpected residue record count")
    return records, provenance


class ParaplumeRunner:
    def __init__(self, store: ExperimentStore):
        self.store = store
        self.executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="ecd-paraplume")
        self._active: set[str] = set()
        self._lock = threading.Lock()

    def submit(self, task_id: str) -> None:
        with self._lock:
            if task_id in self._active:
                return
            self._active.add(task_id)
        self.executor.submit(self._run, task_id)

    def _run(self, task_id: str) -> None:
        try:
            task = self.store.update_experiment_task(
                task_id,
                status="RUNNING",
                progress=0.05,
                output={"message": "Checking full-model cache and preparing inference"},
            )
            input_records = task["input"]["constructs"]
            cached_values: dict[str, dict[str, Any]] = {}
            provenance_by_key: dict[str, dict[str, Any]] = {}
            cache_hits = set()
            for item in input_records:
                cache_key = item["cache_key"]
                cached = self.store.get_paraplume_cache(cache_key)
                if cached is not None:
                    cached_values[cache_key] = cached["result"]
                    provenance_by_key[cache_key] = cached["provenance"]
                    cache_hits.add(cache_key)
            input_by_construct = {
                item["antibody_construct_id"]: item for item in input_records
            }
            missing_inputs = self.store.paraplume_construct_inputs(task_id)
            for item in missing_inputs:
                saved = input_by_construct[item["construct_id"]]
                if hashlib.sha256(item["vh_sequence"].encode("ascii")).hexdigest() != saved["vh_sha256"]:
                    raise ParaplumeTaskError("PARAPLUME_INPUT_CHANGED", "VH sequence changed after task creation")
                if hashlib.sha256(item["vl_sequence"].encode("ascii")).hexdigest() != saved["vl_sha256"]:
                    raise ParaplumeTaskError("PARAPLUME_INPUT_CHANGED", "VL sequence changed after task creation")
            uncached = [
                item
                for item in missing_inputs
                if input_by_construct[item["construct_id"]]["cache_key"] not in cache_hits
            ]
            if uncached:
                progress_state = {"value": 0.12}

                def save_progress(progress: float, message: str) -> None:
                    progress_state["value"] = max(
                        progress_state["value"], min(0.92, max(0.12, progress))
                    )
                    self.store.update_experiment_task(
                        task_id,
                        status="RUNNING",
                        progress=progress_state["value"],
                        output={"stage": "INFERENCE", "message": message},
                    )

                self.store.update_experiment_task(
                    task_id,
                    status="RUNNING",
                    progress=0.12,
                    output={"stage": "INFERENCE", "message": "Preparing the full large model"},
                )
                worker_records, provenance = _run_worker(
                    uncached, task["input"]["model_fingerprint"],
                    progress_callback=save_progress,
                )
                by_construct = {record["construct_id"]: record for record in worker_records}
                for item in input_records:
                    construct_id = item["antibody_construct_id"]
                    if construct_id in by_construct:
                        cache_key = item["cache_key"]
                        cached_values[cache_key] = by_construct[construct_id]
                        provenance_by_key[cache_key] = provenance
            self.store.update_experiment_task(
                task_id,
                status="RUNNING",
                progress=0.96,
                output={"stage": "PERSISTING", "message": "Saving validated residue results and provenance"},
            )
            self.store.complete_paraplume_task(
                task_id,
                results_by_cache_key=cached_values,
                provenance_by_cache_key=provenance_by_key,
                cache_hits=cache_hits,
            )
        except ParaplumeTaskError as exc:
            self.store.update_experiment_task(
                task_id,
                status="FAILED",
                progress=1.0,
                output={"stage": "FAILED", "message": str(exc)},
                error_code=exc.code,
            )
        except Exception as exc:
            self.store.update_experiment_task(
                task_id,
                status="FAILED",
                progress=1.0,
                output={"stage": "FAILED", "message": f"{type(exc).__name__}: {exc}"},
                error_code="PARAPLUME_TASK_FAILED",
            )
        finally:
            with self._lock:
                self._active.discard(task_id)
