from __future__ import annotations

import json
import os
import time
import zipfile
from io import BytesIO, StringIO
from pathlib import Path

import pytest
from ecd_epitope_scanner.experiment_projects import ExperimentStore
from ecd_epitope_scanner.integrations import paraplume_service
from ecd_epitope_scanner.integrations.workers.paraplume_worker import validate_request
from ecd_epitope_scanner.jobs import JobStore
from ecd_epitope_scanner.util import sha256_file
from test_experiment_projects import _completed_v4_job


def _fixture(tmp_path: Path):
    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    project = store.create(_completed_v4_job(jobs)["id"], "Paraplume task fixture")
    antibody = store.add_antibody_construct(
        project["project_id"],
        name="Test antibody",
        role="PARENT",
        vh_sequence="ACDE",
        vl_sequence="FGHI",
        parent_construct_id=None,
        expected_revision=project["revision"],
    )
    return store, store.get(project["project_id"]), antibody


def _crash_paraplume_completion(
    root: str,
    project_id: str,
    antibody_construct_id: str,
    model_fingerprint: str,
    task_record_path: str,
    crash_point: str,
) -> None:
    import ecd_epitope_scanner.experiment_projects as project_module

    workspace = Path(root)
    store = ExperimentStore(workspace, JobStore(workspace))
    project = store.get(project_id)
    task = store.create_paraplume_task(
        project_id,
        antibody_construct_ids=[antibody_construct_id],
        model_fingerprint=model_fingerprint,
        expected_revision=project["revision"],
    )
    store.update_experiment_task(task["id"], status="RUNNING", progress=0.9)
    cache_key = task["input"]["constructs"][0]["cache_key"]
    project = store.get(project_id)
    Path(task_record_path).write_text(
        json.dumps({
            "task_id": task["id"],
            "cache_key": cache_key,
            "project_revision": project["revision"],
            "manifest_artifacts": project["manifest"]["artifacts"],
        }),
        encoding="utf-8",
    )

    if crash_point == "artifact_published":
        original_chmod = Path.chmod

        def chmod_then_exit(path: Path, mode: int, *args, **kwargs):
            result = original_chmod(path, mode, *args, **kwargs)
            if path.name == "paratope_predictions.json":
                os._exit(91)
            return result

        Path.chmod = chmod_then_exit
    elif crash_point == "manifest_published":
        original_write_json = project_module.write_json

        def write_then_exit(path: Path, document: dict):
            original_write_json(path, document)
            if Path(path).name == "project_manifest.json":
                os._exit(92)

        project_module.write_json = write_then_exit
    else:
        raise AssertionError(f"Unknown Paraplume crash point: {crash_point}")

    store.complete_paraplume_task(
        task["id"],
        results_by_cache_key={cache_key: {
            "vh_length": 4,
            "vl_length": 4,
            "vh_probabilities": [0.1, 0.2, 0.3, 0.4],
            "vl_probabilities": [0.4, 0.3, 0.2, 0.1],
        }},
        provenance_by_cache_key={cache_key: {
            "model_version": "2.0.0",
            "actual_device": "cpu",
            "model_fingerprint": model_fingerprint,
        }},
        cache_hits=set(),
    )


def _crash_constraint_export(
    root: str, project_id: str, export_payload: dict, crash_point: str
) -> None:
    import ecd_epitope_scanner.experiment_projects as project_module

    workspace = Path(root)
    store = ExperimentStore(workspace, JobStore(workspace))
    if crash_point == "artifact_published":
        original_chmod = Path.chmod

        def chmod_then_exit(path: Path, mode: int, *args, **kwargs):
            result = original_chmod(path, mode, *args, **kwargs)
            if path.suffix == ".zip":
                os._exit(93)
            return result

        Path.chmod = chmod_then_exit
    elif crash_point == "manifest_published":
        original_write_json = project_module.write_json

        def write_then_exit(path: Path, document: dict):
            original_write_json(path, document)
            if Path(path).name == "project_manifest.json":
                os._exit(94)

        project_module.write_json = write_then_exit
    else:
        raise AssertionError(f"Unknown GRASP export crash point: {crash_point}")
    store.create_constraint_export(project_id, **export_payload)


def test_worker_streams_structured_embedding_progress(monkeypatch):
    fingerprint = "f" * 64
    worker_output = [
        json.dumps({"event": "progress", "progress": 0.34, "message": "Computing IgBert"}),
        json.dumps({
            "status": "success", "records": [{"construct_id": "c1"}],
            "provenance": {"model_version": "2.0.0", "model_fingerprint": fingerprint},
        }),
    ]

    class FakeProcess:
        def __init__(self, *args, **kwargs):
            self.stdin = StringIO()
            self.stdout = StringIO("\n".join(worker_output) + "\n")
            self.stderr = StringIO("worker log\n")
            self.returncode = 0

        def wait(self, timeout=None):
            return self.returncode

    monkeypatch.setattr(paraplume_service.subprocess, "Popen", FakeProcess)
    progress = []
    records, provenance = paraplume_service._run_worker(
        [{"construct_id": "c1", "vh_sequence": "ACD", "vl_sequence": "EFG"}],
        fingerprint,
        progress_callback=lambda value, message: progress.append((value, message)),
    )

    assert records == [{"construct_id": "c1"}]
    assert provenance["model_version"] == "2.0.0"
    assert progress == [(0.34, "Computing IgBert")]


def test_worker_request_requires_full_model_and_valid_sequences():
    request = {
        "schema_version": "1.0",
        "model_size": "large",
        "model_fingerprint": "a" * 64,
        "constructs": [{"construct_id": "c1", "vh_sequence": "ACD", "vl_sequence": "EFG"}],
    }
    assert validate_request(request)[0]["construct_id"] == "c1"
    with pytest.raises(ValueError, match="full large model"):
        validate_request({**request, "model_size": "small"})
    with pytest.raises(ValueError, match="standard amino acids"):
        validate_request({**request, "constructs": [{**request["constructs"][0], "vh_sequence": "ACX"}]})


def test_paraplume_task_is_immutable_cached_and_manifested(tmp_path, monkeypatch):
    store, project, antibody = _fixture(tmp_path)
    fingerprint = "b" * 64
    first = store.create_paraplume_task(
        project["project_id"],
        antibody_construct_ids=[antibody["construct_id"]],
        model_fingerprint=fingerprint,
        expected_revision=project["revision"],
    )
    store.update_experiment_task(first["id"], status="RUNNING", progress=0.1)
    cache_key = first["input"]["constructs"][0]["cache_key"]
    scores = {
        "construct_id": antibody["construct_id"],
        "vh_length": 4,
        "vl_length": 4,
        "vh_probabilities": [0.1, 0.2, 0.3, 0.4],
        "vl_probabilities": [0.4, 0.3, 0.2, 0.1],
    }
    provenance = {"model_version": "2.0.0", "actual_device": "cpu", "model_fingerprint": fingerprint}
    store.complete_paraplume_task(
        first["id"],
        results_by_cache_key={cache_key: scores},
        provenance_by_cache_key={cache_key: provenance},
        cache_hits=set(),
    )
    task = store.get_experiment_task(first["id"])
    assert task["status"] == "COMPLETE"
    assert store.get(project["project_id"])["status"] == "COMPLETE"
    records = store.list_paraplume_runs(project["project_id"])
    assert records[0]["prediction_status"] == "PREDICTED_UNREVIEWED"
    assert records[0]["vh_probabilities"] == scores["vh_probabilities"]

    second = store.create_paraplume_task(
        project["project_id"],
        antibody_construct_ids=[antibody["construct_id"]],
        model_fingerprint=fingerprint,
        expected_revision=store.get(project["project_id"])["revision"],
    )
    monkeypatch.setattr(
        paraplume_service,
        "_run_worker",
        lambda *args, **kwargs: pytest.fail("cache hit should skip full-model execution"),
    )
    paraplume_service.ParaplumeRunner(store)._run(second["id"])
    assert store.get_experiment_task(second["id"])["status"] == "COMPLETE"
    latest_run = store.list_paraplume_runs(project["project_id"])[-1]
    assert latest_run["cache_hit"] is True
    result_artifact_id = latest_run["artifact_id"]
    path, artifact = store.artifact_path(project["project_id"], result_artifact_id)
    assert sha256_file(path) == artifact["sha256"]


@pytest.mark.parametrize(
    ("crash_point", "exit_code"),
    [("artifact_published", 91), ("manifest_published", 92)],
)
def test_paraplume_result_recovers_after_hard_crash(
    tmp_path, crash_point, exit_code
):
    import multiprocessing

    store, project, antibody = _fixture(tmp_path)
    model_fingerprint = "9" * 64
    task_record = tmp_path / "crashed_task.json"
    process = multiprocessing.get_context("spawn").Process(
        target=_crash_paraplume_completion,
        args=(
            str(tmp_path), project["project_id"], antibody["construct_id"],
            model_fingerprint, str(task_record), crash_point,
        ),
    )
    process.start()
    process.join(timeout=30)
    assert not process.is_alive()
    assert process.exitcode == exit_code
    crashed = json.loads(task_record.read_text(encoding="utf-8"))
    orphan = (
        tmp_path / "experiments" / project["project_id"] / "analyses"
        / "paraplume_runs" / crashed["task_id"] / "paratope_predictions.json"
    )
    assert orphan.is_file()

    recovered = ExperimentStore(tmp_path, JobStore(tmp_path))
    recovered_project = recovered.get(project["project_id"])
    assert recovered_project["status"] == "COMPLETE"
    assert recovered_project["revision"] == crashed["project_revision"]
    assert recovered_project["manifest"]["artifacts"] == crashed["manifest_artifacts"]
    interrupted = recovered.get_experiment_task(crashed["task_id"])
    assert interrupted["status"] == "INTERRUPTED"
    assert interrupted["error_code"] == "TASK_INTERRUPTED_BY_RESTART"
    assert recovered.list_paraplume_runs(project["project_id"]) == []

    retry = recovered.create_paraplume_task(
        project["project_id"],
        antibody_construct_ids=[antibody["construct_id"]],
        model_fingerprint=model_fingerprint,
        expected_revision=recovered_project["revision"],
    )
    recovered.update_experiment_task(retry["id"], status="RUNNING", progress=0.9)
    cache_key = retry["input"]["constructs"][0]["cache_key"]
    recovered.complete_paraplume_task(
        retry["id"],
        results_by_cache_key={cache_key: {
            "vh_length": 4,
            "vl_length": 4,
            "vh_probabilities": [0.1, 0.2, 0.3, 0.4],
            "vl_probabilities": [0.4, 0.3, 0.2, 0.1],
        }},
        provenance_by_cache_key={cache_key: {
            "model_version": "2.0.0",
            "actual_device": "cpu",
            "model_fingerprint": model_fingerprint,
        }},
        cache_hits=set(),
    )
    runs = recovered.list_paraplume_runs(project["project_id"])
    assert len(runs) == 1
    result_path, artifact = recovered.artifact_path(project["project_id"], runs[0]["artifact_id"])
    assert result_path != orphan
    assert sha256_file(result_path) == artifact["sha256"]
    assert orphan.is_file()


def test_paraplume_cache_rejects_non_probability_values(tmp_path):
    store, project, antibody = _fixture(tmp_path)
    task = store.create_paraplume_task(
        project["project_id"],
        antibody_construct_ids=[antibody["construct_id"]],
        model_fingerprint="c" * 64,
        expected_revision=project["revision"],
    )
    store.update_experiment_task(task["id"], status="RUNNING", progress=0.1)
    cache_key = task["input"]["constructs"][0]["cache_key"]
    with pytest.raises(ValueError, match="finite values"):
        store.complete_paraplume_task(
            task["id"],
            results_by_cache_key={
                cache_key: {
                    "vh_length": 4,
                    "vl_length": 4,
                    "vh_probabilities": [0.1, 0.2, 0.3, 1.1],
                    "vl_probabilities": [0.1, 0.2, 0.3, 0.4],
                }
            },
            provenance_by_cache_key={
                cache_key: {"model_version": "2.0.0", "model_fingerprint": "c" * 64}
            },
            cache_hits=set(),
        )


def test_restart_marks_only_orphaned_paraplume_tasks_interrupted(tmp_path, monkeypatch):
    import ecd_epitope_scanner.experiment_projects as experiment_projects

    live_pid = os.getpid()
    monkeypatch.setattr(experiment_projects, "_process_alive", lambda pid: pid == live_pid)
    store, project, antibody = _fixture(tmp_path)
    task = store.create_paraplume_task(
        project["project_id"],
        antibody_construct_ids=[antibody["construct_id"]],
        model_fingerprint="e" * 64,
        expected_revision=project["revision"],
    )
    store.update_experiment_task(task["id"], status="RUNNING", progress=0.2)
    same_process = ExperimentStore(tmp_path, JobStore(tmp_path))
    assert same_process.get_experiment_task(task["id"])["status"] == "RUNNING"
    monkeypatch.setattr(experiment_projects, "_process_alive", lambda pid: False)
    restarted = ExperimentStore(tmp_path, JobStore(tmp_path))
    recovered = restarted.get_experiment_task(task["id"])
    assert recovered["status"] == "INTERRUPTED"
    assert recovered["error_code"] == "TASK_INTERRUPTED_BY_RESTART"


def test_constraint_export_links_prediction_and_keeps_ir_only(tmp_path):
    store, project, antibody = _fixture(tmp_path)
    task = store.create_paraplume_task(
        project["project_id"], antibody_construct_ids=[antibody["construct_id"]],
        model_fingerprint="f" * 64, expected_revision=project["revision"],
    )
    store.update_experiment_task(task["id"], status="RUNNING", progress=0.1)
    cache_key = task["input"]["constructs"][0]["cache_key"]
    store.complete_paraplume_task(
        task["id"],
        results_by_cache_key={cache_key: {
            "vh_length": 4, "vl_length": 4,
            "vh_probabilities": [0.1, 0.2, 0.3, 0.4],
            "vl_probabilities": [0.4, 0.3, 0.2, 0.1],
        }},
        provenance_by_cache_key={cache_key: {"model_version": "2.0.0", "model_fingerprint": "f" * 64}},
        cache_hits=set(),
    )
    run = store.list_paraplume_runs(project["project_id"])[0]
    result = store.create_constraint_export(
        project["project_id"],
        antibody_construct_id=antibody["construct_id"],
        experimental_ir=[],
        predicted_ir=[{
            "constraint": "2-2-C", "evidence_id": run["run_id"],
            "rationale": "Selected for review", "details": {},
        }],
        manual_rpr=[], merge_sources=["predicted_ir"], reviewer="reviewer",
        expected_revision=store.get(project["project_id"])["revision"],
    )
    assert result["manifest"]["status"] == "IR_ONLY"
    assert result["manifest"]["rpr_count"] == 0
    assert store.get(project["project_id"])["status"] == "COMPLETE"
    exports = store.list_constraint_exports(project["project_id"])
    assert exports[0]["status"] == "IR_ONLY"
    path, artifact = store.artifact_path(project["project_id"], result["artifact"]["id"])
    assert sha256_file(path) == artifact["sha256"]
    with zipfile.ZipFile(BytesIO(path.read_bytes())) as bundle:
        assert "IR-only" in bundle.read("export_manifest.json").decode()
        assert "PREDICTED_UNREVIEWED" in bundle.read("evidence_records.json").decode()
    restarted = ExperimentStore(tmp_path, JobStore(tmp_path))
    assert restarted.get(project["project_id"])["status"] == "COMPLETE"
    assert restarted.list_constraint_exports(project["project_id"])[0]["sha256"] == artifact["sha256"]
    with pytest.raises(ValueError, match="RPR_COMPLEX_STRUCTURE_REQUIRED"):
        store.create_constraint_export(
            project["project_id"], antibody_construct_id=antibody["construct_id"],
            experimental_ir=[], predicted_ir=[],
            manual_rpr=[{"constraint": "1-1-A, 2-1-A, 8.0, 0.05", "evidence_id": "x",
                         "rationale": "manual complex evidence", "details": {}}],
            merge_sources=["manual_rpr"], reviewer="reviewer",
            expected_revision=store.get(project["project_id"])["revision"],
        )


@pytest.mark.parametrize(
    ("crash_point", "exit_code"),
    [("artifact_published", 93), ("manifest_published", 94)],
)
def test_grasp_export_recovers_after_hard_crash(tmp_path, crash_point, exit_code):
    import multiprocessing

    store, project, antibody = _fixture(tmp_path)
    task = store.create_paraplume_task(
        project["project_id"],
        antibody_construct_ids=[antibody["construct_id"]],
        model_fingerprint="a" * 64,
        expected_revision=project["revision"],
    )
    store.update_experiment_task(task["id"], status="RUNNING", progress=0.9)
    cache_key = task["input"]["constructs"][0]["cache_key"]
    store.complete_paraplume_task(
        task["id"],
        results_by_cache_key={cache_key: {
            "vh_length": 4,
            "vl_length": 4,
            "vh_probabilities": [0.1, 0.2, 0.3, 0.4],
            "vl_probabilities": [0.4, 0.3, 0.2, 0.1],
        }},
        provenance_by_cache_key={cache_key: {
            "model_version": "2.0.0",
            "actual_device": "cpu",
            "model_fingerprint": "a" * 64,
        }},
        cache_hits=set(),
    )
    run = store.list_paraplume_runs(project["project_id"])[0]
    export_payload = {
        "antibody_construct_id": antibody["construct_id"],
        "experimental_ir": [],
        "predicted_ir": [{
            "constraint": "2-2-C",
            "evidence_id": run["run_id"],
            "rationale": "Synthetic recovery fixture",
            "details": {},
        }],
        "manual_rpr": [],
        "merge_sources": ["predicted_ir"],
        "reviewer": "synthetic test reviewer",
        "expected_revision": store.get(project["project_id"])["revision"],
    }
    registered_artifacts_before_crash = store.get(project["project_id"])["manifest"]["artifacts"]
    process = multiprocessing.get_context("spawn").Process(
        target=_crash_constraint_export,
        args=(str(tmp_path), project["project_id"], export_payload, crash_point),
    )
    process.start()
    process.join(timeout=30)
    assert not process.is_alive()
    assert process.exitcode == exit_code

    export_dir = (
        tmp_path / "experiments" / project["project_id"] / "analyses" / "constraint_exports"
    )
    orphan_files = list(export_dir.glob("*.zip"))
    assert len(orphan_files) == 1
    orphan = orphan_files[0]
    assert orphan.is_file()

    recovered = ExperimentStore(tmp_path, JobStore(tmp_path))
    restored = recovered.get(project["project_id"])
    assert restored["status"] == "COMPLETE"
    assert restored["revision"] == export_payload["expected_revision"]
    assert restored["manifest"]["artifacts"] == registered_artifacts_before_crash
    assert recovered.list_constraint_exports(project["project_id"]) == []

    result = recovered.create_constraint_export(
        project["project_id"], **export_payload
    )
    exports = recovered.list_constraint_exports(project["project_id"])
    assert len(exports) == 1
    result_path, artifact = recovered.artifact_path(project["project_id"], result["artifact"]["id"])
    assert result_path != orphan
    assert sha256_file(result_path) == artifact["sha256"]
    assert orphan.is_file()


def test_paraplume_api_queues_and_returns_prediction_artifact(tmp_path, monkeypatch):
    import ecd_epitope_scanner.api as api
    from fastapi.testclient import TestClient

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    project = store.create(_completed_v4_job(jobs)["id"], "Paraplume API fixture")
    antibody = store.add_antibody_construct(
        project["project_id"],
        name="API antibody",
        role="PARENT",
        vh_sequence="ACDE",
        vl_sequence="FGHI",
        parent_construct_id=None,
        expected_revision=project["revision"],
    )
    fingerprint = "d" * 64
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    monkeypatch.setattr(api, "model_fingerprint", lambda: fingerprint)

    worker_progress = []

    def fake_worker(constructs, model_fingerprint, progress_callback=None):
        if progress_callback is not None:
            progress_callback(0.34, "Computing IgBert residue embeddings")
            progress_callback(0.48, "Completed IgBert residue embeddings")
            worker_progress.extend([0.34, 0.48])
        return (
            [
                {
                    "construct_id": item["construct_id"],
                    "vh_length": len(item["vh_sequence"]),
                    "vl_length": len(item["vl_sequence"]),
                    "vh_probabilities": [0.5] * len(item["vh_sequence"]),
                    "vl_probabilities": [0.25] * len(item["vl_sequence"]),
                }
                for item in constructs
            ],
            {
                "model_version": "2.0.0",
                "model_fingerprint": model_fingerprint,
                "actual_device": "cpu",
            },
        )

    monkeypatch.setattr(paraplume_service, "_run_worker", fake_worker)
    monkeypatch.setattr(api, "PARAPLUME_RUNNER", paraplume_service.ParaplumeRunner(store))
    client = TestClient(api.app)
    response = client.post(
        f"/api/v4/experiment-projects/{project['project_id']}/paraplume",
        json={
            "antibody_construct_ids": [antibody["construct_id"]],
            "expected_revision": store.get(project["project_id"])["revision"],
        },
    )
    assert response.status_code == 202, response.text
    task_id = response.json()["data"]["task_id"]
    deadline = time.monotonic() + 5
    task = None
    while time.monotonic() < deadline:
        task = client.get(f"/api/v4/experiment-tasks/{task_id}").json()["data"]
        if task["status"] in {"COMPLETE", "FAILED"}:
            break
        time.sleep(0.02)
    assert task["status"] == "COMPLETE", task
    assert worker_progress == [0.34, 0.48]
    assert client.get(
        f"/api/v4/experiment-projects/{project['project_id']}/paraplume-runs"
    ).json()["data"][0]["prediction_status"] == "PREDICTED_UNREVIEWED"
    run = store.list_paraplume_runs(project["project_id"])[0]
    current = store.get(project["project_id"])
    exported = client.post(
        f"/api/v4/experiment-projects/{project['project_id']}/constraint-exports",
        json={
            "antibody_construct_id": antibody["construct_id"],
            "experimental_ir": [],
            "predicted_ir": [{
                "constraint": "2-2-C", "evidence_id": run["run_id"],
                "rationale": "Reviewed candidate", "details": {},
            }],
            "manual_rpr": [], "merge_sources": ["predicted_ir"],
            "reviewer": "reviewer", "expected_revision": current["revision"],
        },
    )
    assert exported.status_code == 201, exported.text
    export_data = exported.json()["data"]
    downloaded = client.get(
        f"/api/v4/experiment-projects/{project['project_id']}/files/"
        f"{export_data['artifact']['id']}"
    )
    assert downloaded.status_code == 200
    assert downloaded.content[:2] == b"PK"
    no_complex = client.post(
        f"/api/v4/experiment-projects/{project['project_id']}/constraint-exports",
        json={
            "antibody_construct_id": antibody["construct_id"], "experimental_ir": [],
            "predicted_ir": [], "manual_rpr": [{
                "constraint": "1-1-A, 2-1-A, 8.0, 0.05", "evidence_id": "manual",
                "rationale": "No complex structure is available", "details": {},
            }], "merge_sources": ["manual_rpr"], "reviewer": "reviewer",
            "expected_revision": store.get(project["project_id"])["revision"],
        },
    )
    assert no_complex.status_code == 422
    assert no_complex.json()["detail"]["code"] == "RPR_COMPLEX_STRUCTURE_REQUIRED"
