from __future__ import annotations

import json
import multiprocessing
import os
from pathlib import Path

import ecd_epitope_scanner.jobs as jobs_module
import pytest
from ecd_epitope_scanner.jobs import JobStore
from ecd_epitope_scanner.pipeline import PipelineRequest
from ecd_epitope_scanner.util import sha256_file
from fastapi.testclient import TestClient


def _crash_after_project_snapshot_publish(root: str, source_job_id: str) -> None:
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    workspace = Path(root)
    store = ExperimentStore(workspace, JobStore(workspace))

    def terminate_before_database_registration():
        os._exit(73)

    store._connect = terminate_before_database_registration
    store.create(source_job_id, "Crash recovery snapshot")


def _crash_during_structure_evidence_publish(
    root: str, project_id: str, crash_point: str
) -> None:
    import ecd_epitope_scanner.experiment_projects as project_module

    workspace = Path(root)
    store = project_module.ExperimentStore(workspace, JobStore(workspace))
    if crash_point == "first_artifact":
        replace = project_module.os.replace

        def replace_then_exit(source, target):
            replace(source, target)
            if Path(target).name == "residue_reliability.csv":
                os._exit(75)

        project_module.os.replace = replace_then_exit
    elif crash_point == "updated_manifest":
        write_json = project_module.write_json

        def write_then_exit(path, document):
            write_json(path, document)
            if Path(path).name == "project_manifest.json":
                os._exit(76)

        project_module.write_json = write_then_exit
    else:
        raise AssertionError(f"Unknown crash point: {crash_point}")
    store.freeze_source_structure_evidence(project_id, expected_revision=1)


def _crash_flowjo_import_after_workspace_is_read_only(
    root: str, project_id: str, workspace_path: str, preview: dict
) -> None:
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    workspace = Path(root)
    store = ExperimentStore(workspace, JobStore(workspace))
    original_chmod = Path.chmod

    def chmod_then_exit(path: Path, mode: int, *args, **kwargs):
        result = original_chmod(path, mode, *args, **kwargs)
        if path.parent.name == "flowjo" and path.suffix == ".wsp":
            os._exit(77)
        return result

    Path.chmod = chmod_then_exit
    store.register_flowjo_import(
        project_id,
        Path(workspace_path),
        "crash-recovery.wsp",
        preview,
        expected_revision=1,
    )


def _crash_fcs_import_after_artifact_is_read_only(
    root: str, project_id: str, source_path: str, metadata: dict
) -> None:
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    workspace = Path(root)
    store = ExperimentStore(workspace, JobStore(workspace))
    original_chmod = Path.chmod

    def chmod_then_exit(path: Path, mode: int, *args, **kwargs):
        result = original_chmod(path, mode, *args, **kwargs)
        if path.parent.name == "fcs" and path.suffix == ".fcs":
            os._exit(80)
        return result

    Path.chmod = chmod_then_exit
    store.import_fcs_files(
        project_id,
        [{"path": source_path, "original_filename": "crash-recovery.fcs", "metadata": metadata}],
        batch_id=None,
        expected_revision=1,
    )


def _crash_compensation_matrix_registration(
    root: str, project_id: str, matrix_document: dict, crash_point: str
) -> None:
    import ecd_epitope_scanner.experiment_projects as project_module

    workspace = Path(root)
    store = project_module.ExperimentStore(workspace, JobStore(workspace))
    if crash_point == "artifact_published":
        original_chmod = Path.chmod

        def chmod_then_exit(path: Path, mode: int, *args, **kwargs):
            result = original_chmod(path, mode, *args, **kwargs)
            if path.parent.name == "matrices" and path.suffix == ".json":
                os._exit(81)
            return result

        Path.chmod = chmod_then_exit
    elif crash_point == "manifest_published":
        original_write_json = project_module.write_json

        def write_then_exit(path: Path, document: dict):
            original_write_json(path, document)
            if Path(path).name == "project_manifest.json":
                os._exit(82)

        project_module.write_json = write_then_exit
    else:
        raise AssertionError(f"Unknown compensation matrix crash point: {crash_point}")
    store.register_compensation_matrix(project_id, matrix_document, expected_revision=1)


def _gate_artifact_recovery_fixture(root: Path) -> dict:
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Gate artifact recovery")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id, name="Synthetic gate batch", date="2026-09-29", assay="ANTIGEN_MYC",
        cell_line="unspecified", dose=None, dose_unit=None, positive_gate_id="positive",
        independent_batch=True, settings={}, expected_revision=1,
    )
    source = root / "gate-recovery.fcs"
    source.write_bytes(b"synthetic gate recovery FCS")
    imported = store.import_fcs_files(
        project_id,
        [{"path": str(source), "original_filename": "H1.fcs", "metadata": {
            "event_count": 10, "channels": [{"pnn": "FSC-A"}], "metadata": {},
            "compensation_source": "USER_CONFIRMED_NOT_APPLIED",
        }}],
        batch_id=batch["batch_id"], expected_revision=2,
    )
    upload = imported["files"][0]
    store.confirm_fcs_matches(
        project_id, imported["upload_session_id"], [{
            "upload_id": upload["upload_id"], "batch_id": batch["batch_id"], "role": "NC",
            "replicate_type": "SINGLE", "replicate_index": 1,
            "condition": {"sample_label": "H1"},
        }], expected_revision=3,
    )
    gate = store.save_gate_revision(
        project_id, name="Synthetic recovery gate", gates=[{
            "gate_id": "positive", "name": "Positive", "parent_gate_id": "root",
            "kind": "RANGE_1D", "channels": ["FSC-A"],
            "parameters": {"min": 0, "max": 100},
        }], compensation_decision="NONE", compensation_artifact_id=None,
        parent_revision_id=None, expected_revision=4,
    )
    return {
        "store": store,
        "project_id": project_id,
        "gate_revision_id": gate["gate_revision_id"],
        "expected_revision": gate["project_revision"],
        "input_artifact_ids": [upload["artifact_id"]],
        "worker_result": {
            "samples": [{
                "event_count": 10,
                "gate_counts": [{
                    "gate_id": "positive", "gate_name": "Positive", "gate_path": ["root"],
                    "parent_gate": "root", "parent_count": 10, "count": 5,
                    "absolute_fraction": 0.5, "parent_fraction": 0.5,
                }],
            }],
            "provenance": {"worker": "synthetic_recovery_fixture"},
        },
    }


def _crash_gate_analysis_save(
    root: str, project_id: str, gate_revision_id: str, input_artifact_ids: list[str],
    worker_result: dict, expected_revision: int, crash_point: str,
) -> None:
    import ecd_epitope_scanner.experiment_projects as project_module

    workspace = Path(root)
    store = project_module.ExperimentStore(workspace, JobStore(workspace))
    if crash_point == "artifact_published":
        original_chmod = Path.chmod

        def chmod_then_exit(path: Path, mode: int, *args, **kwargs):
            result = original_chmod(path, mode, *args, **kwargs)
            if path.name == "gate_counts.json":
                os._exit(83)
            return result

        Path.chmod = chmod_then_exit
    elif crash_point == "manifest_published":
        original_write_json = project_module.write_json

        def write_then_exit(path: Path, document: dict):
            original_write_json(path, document)
            if Path(path).name == "project_manifest.json":
                os._exit(84)

        project_module.write_json = write_then_exit
    else:
        raise AssertionError(f"Unknown gate analysis crash point: {crash_point}")
    store.save_gate_analysis(
        project_id, gate_revision_id=gate_revision_id,
        input_artifact_ids=input_artifact_ids, worker_result=worker_result,
        expected_revision=expected_revision,
    )


def _crash_effect_analysis_save(
    root: str, project_id: str, gate_analysis_id: str, expected_revision: int,
    crash_point: str,
) -> None:
    import ecd_epitope_scanner.experiment_projects as project_module

    workspace = Path(root)
    store = project_module.ExperimentStore(workspace, JobStore(workspace))
    if crash_point == "artifact_published":
        original_chmod = Path.chmod

        def chmod_then_exit(path: Path, mode: int, *args, **kwargs):
            result = original_chmod(path, mode, *args, **kwargs)
            if path.name == "effect_results.json":
                os._exit(85)
            return result

        Path.chmod = chmod_then_exit
    elif crash_point == "manifest_published":
        original_write_json = project_module.write_json

        def write_then_exit(path: Path, document: dict):
            original_write_json(path, document)
            if Path(path).name == "project_manifest.json":
                os._exit(86)

        project_module.write_json = write_then_exit
    else:
        raise AssertionError(f"Unknown effect analysis crash point: {crash_point}")
    store.save_effect_analysis(
        project_id, gate_analysis_id=gate_analysis_id, positive_gate_id="positive",
        document={"synthetic_storage_recovery_probe": True},
        expected_revision=expected_revision,
    )


def _crash_flowjo_percentage_effect_save(
    root: str,
    project_id: str,
    import_id: str,
    document: dict,
    expected_revision: int,
    crash_point: str,
) -> None:
    import ecd_epitope_scanner.experiment_projects as project_module

    workspace = Path(root)
    store = project_module.ExperimentStore(workspace, JobStore(workspace))
    if crash_point == "artifact_published":
        original_chmod = Path.chmod

        def chmod_then_exit(path: Path, mode: int, *args, **kwargs):
            result = original_chmod(path, mode, *args, **kwargs)
            if path.name == "effect_results.json":
                os._exit(95)
            return result

        Path.chmod = chmod_then_exit
    elif crash_point == "manifest_published":
        original_write_json = project_module.write_json

        def write_then_exit(path: Path, value: dict):
            original_write_json(path, value)
            if Path(path).name == "project_manifest.json":
                os._exit(96)

        project_module.write_json = write_then_exit
    else:
        raise AssertionError(f"Unknown FlowJo percentage effect crash point: {crash_point}")
    store.save_flowjo_percentage_effect_analysis(
        project_id,
        import_id=import_id,
        document=document,
        expected_revision=expected_revision,
    )


def _crash_flowjo_percentage_import(
    root: str, project_id: str, source_path: str, parsed_table: dict, crash_point: str
) -> None:
    import ecd_epitope_scanner.experiment_projects as project_module

    workspace = Path(root)
    store = project_module.ExperimentStore(workspace, JobStore(workspace))
    if crash_point == "source_published":
        original_chmod = Path.chmod

        def chmod_then_exit(path: Path, mode: int, *args, **kwargs):
            result = original_chmod(path, mode, *args, **kwargs)
            if path.parent.name == "flowjo_result_tables":
                os._exit(78)
            return result

        Path.chmod = chmod_then_exit
    elif crash_point == "manifest_published":
        original_write_json = project_module.write_json

        def write_then_exit(path: Path, document: dict):
            original_write_json(path, document)
            if Path(path).name == "project_manifest.json":
                os._exit(79)

        project_module.write_json = write_then_exit
    else:
        raise AssertionError(f"Unknown FlowJo percentage crash point: {crash_point}")
    store.register_flowjo_percentage_import(
        project_id,
        batch_id=parsed_table["batch_id"],
        source_path=Path(source_path),
        parsed_table=parsed_table["table"],
        expected_revision=parsed_table["expected_revision"],
    )


def _concurrent_add_assay_batch(
    root: str, project_id: str, name: str, start_event, results
) -> None:
    from ecd_epitope_scanner.experiment_projects import ExperimentStore, RevisionConflict

    workspace = Path(root)
    store = ExperimentStore(workspace, JobStore(workspace))
    if not start_event.wait(timeout=20):
        results.put(("ERROR", "start barrier timed out"))
        return
    try:
        created = store.add_batch(
            project_id,
            name=name,
            date="2026-09-28",
            assay="ACTIVATION",
            cell_line="Jurkat",
            dose=1.0,
            dose_unit="ug/mL",
            positive_gate_id="positive",
            independent_batch=True,
            settings={},
            expected_revision=1,
        )
    except RevisionConflict:
        results.put(("REVISION_CONFLICT", name))
    except Exception as exc:
        results.put(("ERROR", f"{type(exc).__name__}: {exc}"))
    else:
        results.put(("CREATED", created["batch_id"]))


def _completed_v4_job(jobs: JobStore, *, graph: bool = True) -> dict:
    job = jobs.create(PipelineRequest(input_type="structure", design_engine="joint_v4_ml"))
    source = Path(job["workdir"])
    (source / "sequences.fasta").write_text(">WT_ECD\nACD\n", encoding="utf-8")
    (source / "structure_for_energy.pdb").write_text("MODEL        1\nENDMDL\n", encoding="utf-8")
    (source / "residue_mapping.csv").write_text(
        "normalized_position,uniprot_pos,chain_id,auth_seq_id,insertion_code,label_seq_id,one_letter\n"
        "1,10,A,10,,1,A\n2,11,A,11,,2,C\n3,12,A,12,,3,D\n",
        encoding="utf-8",
    )
    (source / "probe_sectors_v4.json").write_text(
        json.dumps([{"sector_id": "VS001", "residue_positions": [10, 11]}]), encoding="utf-8"
    )
    if graph:
        (source / "surface_graph_v4.json").write_text(
            json.dumps(
                {
                    "schema_version": "1.0",
                    "position_system": "uniprot_pos",
                    "mapping_sha256": sha256_file(source / "residue_mapping.csv"),
                    "structure_sha256": sha256_file(source / "structure_for_energy.pdb"),
                    "nodes": [
                        {"uniprot_pos": 10, "sector_id": "VS001"},
                        {"uniprot_pos": 11, "sector_id": "VS001"},
                    ],
                    "edges": [{"source": 10, "target": 11, "distance_a": 3.8}],
                }
            ),
            encoding="utf-8",
        )
    (source / "v4_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "4.0",
                "engine": "joint_v4_ml",
                "software": {"revision": "fixture"},
                "provenance": {
                    "surface_graph": {
                        "snapshot_sha256": (
                            sha256_file(source / "surface_graph_v4.json") if graph else "missing"
                        )
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    (source / "analysis.json").write_text(
        json.dumps({"protein": {"ecd_sequence": "ACD"}}), encoding="utf-8"
    )
    (source / "designs.json").write_text(json.dumps([{
        "id": "V4_01", "name": "Recommended A10C", "status": "RECOMMENDED_LOW_PREDICTED_RISK",
        "mutations": [{"position": 10, "wt": "A", "mutant": "C"}],
        "mutant_ecd_sequence": "CCD", "warnings": [],
    }]), encoding="utf-8")
    return jobs.update(job["id"], status="COMPLETE", stage="COMPLETE")


def _add_structure_evidence(job: dict) -> dict[str, str]:
    source = Path(job["workdir"])
    reliability = (
        "residue_id,state,plddt,mean_local_pae,mean_core_pae,quality_flags,"
        "use_for_surface_mesh,use_for_3d_footprint,use_for_epitope_screening\n"
        '10,STRUCTURED,81.0,2.1,8.5,[],True,True,True\n'
        '11,UNCERTAIN,55.0,6.2,16.0,"[""TRANSITION""]",True,False,True\n'
        '12,DISORDERED,42.0,7.0,18.0,[],False,False,True\n'
    )
    sasa = (
        "position,native_rsasa,procaliper_rsasa,absolute_delta,status\n"
        "10,0.2,0.21,0.01,AGREE\n"
        "11,0.3,0.38,0.08,SASA_DISAGREEMENT\n"
        "12,0.4,0.42,0.02,AGREE\n"
    )
    (source / "residue_reliability.csv").write_text(reliability, encoding="utf-8")
    (source / "sasa_comparison.csv").write_text(sasa, encoding="utf-8")
    (source / "reliability_qc.json").write_text(json.dumps({
        "classification_version": "reliability_rules_v1",
        "predicted": True,
        "state_counts": {"STRUCTURED": 1, "UNCERTAIN": 1, "DISORDERED": 1},
    }), encoding="utf-8")
    declared = {
        filename: sha256_file(source / filename)
        for filename in (
            "residue_reliability.csv", "sasa_comparison.csv", "reliability_qc.json"
        )
    }
    (source / "run_manifest.json").write_text(json.dumps({"files": declared}), encoding="utf-8")
    return declared


def _add_unmapped_residue_evidence(job: dict, *, state: str = "UNRESOLVED") -> None:
    source = Path(job["workdir"])
    (source / "sequences.fasta").write_text(">WT_ECD|P16410:10-13\nACDE\n", encoding="utf-8")
    (source / "analysis.json").write_text(
        json.dumps({"protein": {"ecd_sequence": "ACDE"}}), encoding="utf-8"
    )
    (source / "designs.json").write_text(json.dumps([{
        "id": "V4_01", "name": "Recommended A10C", "status": "RECOMMENDED_LOW_PREDICTED_RISK",
        "mutations": [{"position": 10, "wt": "A", "mutant": "C"}],
        "mutant_ecd_sequence": "CCDE", "warnings": [],
    }]), encoding="utf-8")
    with (source / "residue_reliability.csv").open("a", encoding="utf-8") as handle:
        handle.write(f"13,{state},,,,[],False,False,False\n")
    (source / "reliability_qc.json").write_text(json.dumps({
        "classification_version": "reliability_rules_v1",
        "predicted": True,
        "state_counts": {
            "STRUCTURED": 1,
            "UNCERTAIN": 1,
            "DISORDERED": 1,
            state: 1,
        },
    }), encoding="utf-8")
    run_manifest = json.loads((source / "run_manifest.json").read_text(encoding="utf-8"))
    run_manifest["files"] = {
        filename: sha256_file(source / filename)
        for filename in (
            "residue_reliability.csv", "sasa_comparison.csv", "reliability_qc.json"
        )
    }
    (source / "run_manifest.json").write_text(json.dumps(run_manifest), encoding="utf-8")


def test_unresolved_unmapped_terminal_residue_is_preserved_in_project_snapshot(tmp_path):
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    job = _completed_v4_job(jobs)
    _add_structure_evidence(job)
    _add_unmapped_residue_evidence(job)
    store = ExperimentStore(tmp_path, jobs)

    project = store.create(job["id"], "Unresolved terminal residue")

    assert project["status"] == "COMPLETE"
    assert project["manifest"]["reference_sequence"] == "ACDE"
    reliability_artifact = next(
        item for item in project["artifacts"] if item["kind"] == "residue_reliability"
    )
    reliability_path, _ = store.artifact_path(project["project_id"], reliability_artifact["id"])
    assert "13,UNRESOLVED" in reliability_path.read_text(encoding="utf-8")


def test_unmapped_reliability_record_must_be_unresolved_and_excluded_from_geometry(tmp_path):
    from ecd_epitope_scanner.experiment_projects import ExperimentStore, SnapshotError

    jobs = JobStore(tmp_path)
    job = _completed_v4_job(jobs)
    _add_structure_evidence(job)
    _add_unmapped_residue_evidence(job, state="STRUCTURED")
    store = ExperimentStore(tmp_path, jobs)

    with pytest.raises(SnapshotError, match="must be UNRESOLVED"):
        store.create(job["id"], "Invalid unmapped reliability")
    assert store.list() == []


def test_snapshot_api_restart_and_integrity(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    monkeypatch.setattr(api, "STORE", jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    client = TestClient(api.app)
    created = client.post(
        "/api/v4/experiment-projects", json={"source_job_id": job["id"], "name": "PD-L1 activation"}
    )
    assert created.status_code == 201, created.text
    project = created.json()["data"]
    assert project["status"] == "COMPLETE"
    assert len(project["artifacts"]) == 7
    assert project["frozen_hashes"]["surface_graph"] == sha256_file(
        Path(job["workdir"]) / "surface_graph_v4.json"
    )
    project_id = project["project_id"]
    assert client.get(f"/api/v4/experiment-projects/{project_id}").status_code == 200
    for artifact in project["artifacts"]:
        result = client.get(f"/api/v4/experiment-projects/{project_id}/files/{artifact['id']}")
        assert result.status_code == 200
        assert len(result.content) == artifact["size_bytes"]
    assert (
        client.get(f"/api/v4/experiment-projects/{project_id}/files/not-registered").status_code
        == 404
    )

    restarted = ExperimentStore(root, JobStore(root))
    assert restarted.get(project_id)["status"] == "COMPLETE"
    assert len(restarted.list()) == 1
    assert JobStore(root).get(job["id"])["status"] == "COMPLETE"
    monkeypatch.setattr(api, "EXPERIMENTS", restarted)
    assert client.get("/api/v4/experiment-projects").json()["data"][0]["status"] == "COMPLETE"

    mapping = root / "experiments" / project_id / "reference" / "residue_mapping.csv"
    mapping.write_text("tampered", encoding="utf-8")
    assert (
        client.get(f"/api/v4/experiment-projects/{project_id}").json()["data"]["status"]
        == "INTEGRITY_ERROR"
    )
    artifact_id = next(a["id"] for a in project["artifacts"] if a["kind"] == "mapping")
    blocked = client.get(f"/api/v4/experiment-projects/{project_id}/files/{artifact_id}")
    assert blocked.status_code == 409
    assert blocked.json()["detail"]["code"] == "INTEGRITY_ERROR"


def test_snapshot_publish_is_removed_when_database_registration_fails(monkeypatch, tmp_path):
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    job = _completed_v4_job(jobs)
    connect = store._connect
    failed = False

    def fail_first_database_connection():
        nonlocal failed
        if not failed:
            failed = True
            raise RuntimeError("injected failure after snapshot publish")
        return connect()

    monkeypatch.setattr(store, "_connect", fail_first_database_connection)
    with pytest.raises(RuntimeError, match="injected failure"):
        store.create(job["id"], "Interrupted snapshot")

    assert store.list() == []
    assert list((root / "experiments").iterdir()) == []
    restarted = ExperimentStore(root, jobs)
    assert restarted.list() == []


def test_restart_recovers_project_after_process_crashes_between_publish_and_commit(tmp_path):
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    original_store = ExperimentStore(root, jobs)
    job = _completed_v4_job(jobs)
    process = multiprocessing.get_context("spawn").Process(
        target=_crash_after_project_snapshot_publish,
        args=(str(root), job["id"]),
    )
    process.start()
    process.join(timeout=30)
    if process.is_alive():
        process.terminate()
        process.join(timeout=5)
        pytest.fail("snapshot crash-injection child did not terminate")
    assert process.exitcode == 73
    assert original_store.list() == []
    orphaned_directory = next(
        path for path in (root / "experiments").iterdir()
        if path.is_dir() and len(path.name) == 16
    )

    restarted = ExperimentStore(root, JobStore(root))
    projects = restarted.list()
    assert len(projects) == 1
    recovered = restarted.get(orphaned_directory.name)
    assert recovered["status"] == "COMPLETE"
    assert recovered["revision"] == 1
    assert {item["kind"] for item in recovered["artifacts"]} >= {
        "sequence", "structure", "mapping", "sectors", "surface_graph",
        "v4_manifest", "recommendations",
    }


def test_frozen_structure_context_serves_mapped_positions_without_local_paths(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Structure context")
    client = TestClient(api.app)
    response = client.get(
        f"/api/v4/experiment-projects/{project['project_id']}/structure-context"
    )
    assert response.status_code == 200, response.text
    context = response.json()["data"]
    assert context["structure"]["format"] == "pdb"
    assert context["structure"]["artifact_id"] in context["structure"]["url"]
    assert [item["uniprot_pos"] for item in context["mapping"]] == [10, 11, 12]
    assert context["surface_graph"]["nodes"][0]["sector_id"] == "VS001"
    assert str(tmp_path) not in response.text


def test_structure_evidence_is_explicitly_frozen_verified_and_exposed_by_position(
    monkeypatch, tmp_path
):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Structure evidence")
    project_id = project["project_id"]
    assert store.read_structure_evidence(project_id) is None

    hashes = _add_structure_evidence(job)
    client = TestClient(api.app)
    stale = client.post(
        f"/api/v4/experiment-projects/{project_id}/structure-evidence",
        json={"expected_revision": 2},
    )
    assert stale.status_code == 409
    assert not any(item["kind"] == "residue_reliability" for item in store.get(project_id)["artifacts"])

    frozen = client.post(
        f"/api/v4/experiment-projects/{project_id}/structure-evidence",
        json={"expected_revision": 1},
    )
    assert frozen.status_code == 200, frozen.text
    assert frozen.json()["data"]["project_revision"] == 2
    repeated = client.post(
        f"/api/v4/experiment-projects/{project_id}/structure-evidence",
        json={"expected_revision": 1},
    )
    assert repeated.status_code == 200
    assert repeated.json()["data"]["already_frozen"] is True
    assert repeated.json()["data"]["project_revision"] == 2
    current = store.get(project_id)
    evidence_artifacts = {
        item["kind"]: item for item in current["artifacts"]
        if item["kind"] in {"residue_reliability", "sasa_comparison", "reliability_qc"}
    }
    assert set(evidence_artifacts) == {Path(name).stem for name in hashes}
    assert {
        f"{kind}.json" if kind == "reliability_qc" else f"{kind}.csv": item["sha256"]
        for kind, item in evidence_artifacts.items()
    } == hashes
    assert all(item["read_only"] for item in evidence_artifacts.values())

    context = client.get(
        f"/api/v4/experiment-projects/{project_id}/structure-context"
    )
    assert context.status_code == 200, context.text
    evidence = context.json()["data"]["structure_evidence"]
    assert evidence["reliability"]["predicted"] is True
    assert evidence["reliability"]["records"][1]["state"] == "UNCERTAIN"
    assert evidence["reliability"]["records"][1]["quality_flags"] == ["TRANSITION"]
    assert evidence["sasa"]["records"][1]["status"] == "SASA_DISAGREEMENT"
    assert evidence["sasa"]["status_counts"] == {"AGREE": 2, "SASA_DISAGREEMENT": 1}
    assert str(tmp_path) not in context.text
    restarted = ExperimentStore(tmp_path, JobStore(tmp_path))
    assert restarted.read_structure_evidence(project_id) == evidence

    fresh_job = _completed_v4_job(jobs)
    _add_structure_evidence(fresh_job)
    fresh_project = restarted.create(fresh_job["id"], "Evidence present at creation")
    assert restarted.read_structure_evidence(fresh_project["project_id"])["available"] is True


@pytest.mark.parametrize(
    ("crash_point", "exit_code"),
    [("first_artifact", 75), ("updated_manifest", 76)],
)
def test_structure_evidence_freeze_recovers_after_process_crash(
    tmp_path, crash_point, exit_code
):
    import sqlite3

    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Interrupted structure evidence")
    _add_structure_evidence(job)
    process = multiprocessing.get_context("spawn").Process(
        target=_crash_during_structure_evidence_publish,
        args=(str(root), project["project_id"], crash_point),
    )
    process.start()
    process.join(timeout=30)
    if process.is_alive():
        process.terminate()
        process.join(timeout=5)
        pytest.fail("structure evidence crash-injection child did not terminate")
    assert process.exitcode == exit_code

    restarted = ExperimentStore(root, JobStore(root))
    recovered_project = restarted.get(project["project_id"])
    assert recovered_project["status"] == "COMPLETE"
    assert recovered_project["revision"] == 1
    completed = restarted.freeze_source_structure_evidence(
        project["project_id"], expected_revision=1
    )
    assert completed["project_revision"] == 2
    assert restarted.read_structure_evidence(project["project_id"])["available"] is True
    assert {
        item["kind"] for item in restarted.get(project["project_id"])["artifacts"]
    } >= {"residue_reliability", "sasa_comparison", "reliability_qc"}
    with sqlite3.connect(root / "experiments.sqlite3") as db:
        repair_count = db.execute(
            "SELECT count(*) FROM audit_events WHERE project_id=? AND action='MANIFEST_RECOVERED'",
            (project["project_id"],),
        ).fetchone()[0]
    assert repair_count == (1 if crash_point == "updated_manifest" else 0)
    ExperimentStore(root, JobStore(root))
    with sqlite3.connect(root / "experiments.sqlite3") as db:
        assert db.execute(
            "SELECT count(*) FROM audit_events WHERE project_id=? AND action='MANIFEST_RECOVERED'",
            (project["project_id"],),
        ).fetchone()[0] == repair_count


def test_restart_does_not_repair_changed_snapshot_identity_fields(tmp_path):
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Identity tamper")
    manifest_path = (
        root / "experiments" / project["project_id"] / "project_manifest.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["name"] = "Unverified replacement"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    restarted = ExperimentStore(root, jobs)
    recovered = restarted.get(project["project_id"])
    assert recovered["status"] == "INTEGRITY_ERROR"
    assert "manifest identity differs from database" in recovered["integrity_errors"]


def test_concurrent_processes_cannot_commit_the_same_project_revision(tmp_path):
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Cross-process revision conflict")
    context = multiprocessing.get_context("spawn")
    start_event = context.Event()
    results = context.Queue()
    workers = [
        context.Process(
            target=_concurrent_add_assay_batch,
            args=(str(root), project["project_id"], name, start_event, results),
        )
        for name in ("Concurrent A", "Concurrent B")
    ]
    for worker in workers:
        worker.start()
    start_event.set()
    for worker in workers:
        worker.join(timeout=30)
        if worker.is_alive():
            worker.terminate()
            worker.join(timeout=5)
            pytest.fail("concurrent revision writer did not terminate")
        assert worker.exitcode == 0

    outcome_rows = [results.get(timeout=5) for _ in workers]
    outcomes = [row[0] for row in outcome_rows]
    assert sorted(outcomes) == ["CREATED", "REVISION_CONFLICT"], outcome_rows
    restarted = ExperimentStore(root, JobStore(root))
    batches = restarted.list_entities(project["project_id"], "assay_batches")
    assert len(batches) == 1
    assert restarted.get(project["project_id"])["revision"] == 2
    assert restarted.get(project["project_id"])["status"] == "COMPLETE"
    results.close()
    results.join_thread()


def test_api_concurrent_revision_writes_return_one_success_and_one_conflict(monkeypatch, tmp_path):
    import threading
    from concurrent.futures import ThreadPoolExecutor

    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Concurrent API revision")
    barrier = threading.Barrier(2)

    def create_batch(name: str):
        payload = {
            "name": name,
            "date": "2026-09-28",
            "assay": "ACTIVATION",
            "cell_line": "Jurkat",
            "dose": 1.0,
            "dose_unit": "ug/mL",
            "positive_gate_id": "positive",
            "independent_batch": True,
            "settings": {},
            "expected_revision": 1,
        }
        with TestClient(api.app) as client:
            barrier.wait(timeout=10)
            return client.post(
                f"/api/v4/experiment-projects/{project['project_id']}/batches",
                json=payload,
            )

    with ThreadPoolExecutor(max_workers=2) as pool:
        responses = list(pool.map(create_batch, ("API concurrent A", "API concurrent B")))

    assert sorted(response.status_code for response in responses) == [201, 409]
    conflict = next(response for response in responses if response.status_code == 409)
    assert conflict.json()["detail"]["code"] == "REVISION_CONFLICT"
    assert len(store.list_entities(project["project_id"], "assay_batches")) == 1
    assert store.get(project["project_id"])["revision"] == 2


def test_flowjo_import_reuses_readonly_workspace_after_hard_crash(tmp_path):
    import stat

    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "FlowJo file publish recovery")
    project_id = project["project_id"]
    workspace = tmp_path / "source-workspace.wsp"
    workspace.write_text("<Workspace fixture='read-only recovery' />", encoding="utf-8")
    preview = {
        "schema_version": "1.0",
        "status": "COMPLETE",
        "records": {
            "status": "COMPLETE", "samples": [], "gate_counts": [],
            "unresolved_nodes": [],
        },
        "flowjo_input_artifact_ids": [],
    }

    context = multiprocessing.get_context("spawn")
    worker = context.Process(
        target=_crash_flowjo_import_after_workspace_is_read_only,
        args=(str(tmp_path), project_id, str(workspace), preview),
    )
    worker.start()
    worker.join(timeout=30)
    if worker.is_alive():
        worker.terminate()
        worker.join(timeout=5)
        pytest.fail("FlowJo workspace publish crash worker did not terminate")
    assert worker.exitcode == 77

    folder = tmp_path / "experiments" / project_id
    target = folder / "raw" / "flowjo" / f"{sha256_file(workspace)}.wsp"
    assert target.is_file()
    assert not (target.stat().st_mode & stat.S_IWRITE)
    manifest = json.loads((folder / "project_manifest.json").read_text(encoding="utf-8"))
    assert all(item["kind"] != "raw_flowjo" for item in manifest["artifacts"])

    restarted = ExperimentStore(tmp_path, JobStore(tmp_path))
    assert restarted.get(project_id)["status"] == "COMPLETE"
    assert restarted.get(project_id)["revision"] == 1
    imported = restarted.register_flowjo_import(
        project_id, workspace, "crash-recovery.wsp", preview, expected_revision=1
    )
    assert imported["status"] == "COMPLETE"
    assert imported["project_revision"] == 2
    assert sha256_file(target) == sha256_file(workspace)
    final_store = ExperimentStore(tmp_path, JobStore(tmp_path))
    assert final_store.get(project_id)["status"] == "COMPLETE"
    assert {
        item["kind"] for item in final_store.get(project_id)["artifacts"]
    } >= {"raw_flowjo", "flowjo_import_preview"}


@pytest.mark.parametrize(
    ("crash_point", "exit_code", "expected_manifest_repairs"),
    [("source_published", 78, 0), ("manifest_published", 79, 1)],
)
def test_flowjo_percentage_import_reuses_verified_source_after_hard_crash(
    tmp_path, crash_point, exit_code, expected_manifest_repairs
):
    import sqlite3
    import stat

    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "FlowJo percentage import recovery")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id, name="Activation", date="2026-09-28", assay="ACTIVATION",
        cell_line="Jurkat", dose=None, dose_unit=None, positive_gate_id="positive",
        independent_batch=True, settings={}, expected_revision=1,
    )
    fcs = tmp_path / "A1.fcs"
    fcs.write_bytes(b"activation sample A1")
    imported_fcs = store.import_fcs_files(
        project_id,
        [{"path": str(fcs), "original_filename": "A1.fcs", "metadata": {
            "event_count": 10, "channels": [], "metadata": {},
            "compensation_source": "USER_CONFIRMED_NOT_APPLIED",
        }}],
        batch_id=batch["batch_id"], expected_revision=2,
    )
    upload = imported_fcs["files"][0]
    store.confirm_fcs_matches(
        project_id, imported_fcs["upload_session_id"], [{
            "upload_id": upload["upload_id"], "batch_id": batch["batch_id"],
            "role": "ACTIVATION", "replicate_type": "TECHNICAL", "replicate_index": 1,
            "temporary_identity": {
                "status": "TEMPORARY_LABEL", "source": "USER_CONFIRMED_LABELS",
                "antigen_label": "CTLA4-WT", "antigen_role": "REFERENCE",
                "panel_label": "L32", "panel_role": "ANTIBODY",
            },
            "condition": {"well_id": "A1", "condition_label": "WT", "panel_label": "L32"},
        }],
        expected_revision=3,
    )
    source = tmp_path / "flowjo-results.csv"
    source.write_text("well,percent\nA1,72.0%\n", encoding="utf-8")
    parsed_table = {
        "batch_id": batch["batch_id"],
        "expected_revision": 4,
        "table": {
            "source_filename": source.name,
            "sheet_name": "",
            "gate_label": "Target-positive | Freq. of Parent",
            "metric": "FLOWJO_FREQ_OF_PARENT",
            "observations": [{
                "well_id": "A1", "fraction": 0.72, "raw_value": "72.0 %",
                "display_precision_decimal_places": 1, "source_row": 2,
            }],
        },
    }

    process = multiprocessing.get_context("spawn").Process(
        target=_crash_flowjo_percentage_import,
        args=(str(root), project_id, str(source), parsed_table, crash_point),
    )
    process.start()
    process.join(timeout=30)
    if process.is_alive():
        process.terminate()
        process.join(timeout=5)
        pytest.fail("FlowJo percentage import crash worker did not terminate")
    assert process.exitcode == exit_code

    source_hash = sha256_file(source)
    source_target = (
        root / "experiments" / project_id / "raw" / "flowjo_result_tables"
        / f"{source_hash}.csv"
    )
    assert source_target.is_file()
    assert sha256_file(source_target) == source_hash
    assert not source_target.stat().st_mode & stat.S_IWRITE

    restarted = ExperimentStore(root, JobStore(root))
    recovered = restarted.get(project_id)
    assert recovered["status"] == "COMPLETE"
    assert recovered["revision"] == 4
    with sqlite3.connect(root / "experiments.sqlite3") as db:
        repair_count = db.execute(
            "SELECT count(*) FROM audit_events WHERE project_id=? AND action='MANIFEST_RECOVERED'",
            (project_id,),
        ).fetchone()[0]
    assert repair_count == expected_manifest_repairs

    completed = restarted.register_flowjo_percentage_import(
        project_id,
        batch_id=batch["batch_id"],
        source_path=source,
        parsed_table=parsed_table["table"],
        expected_revision=4,
    )
    assert completed["project_revision"] == 5
    assert completed["source_artifact"]["sha256"] == source_hash
    assert completed["status"] == "IMPORTED_PENDING_MANUAL_QC"
    assert len(restarted.list_flowjo_percentage_imports(project_id)) == 1
    final_project = ExperimentStore(root, JobStore(root)).get(project_id)
    assert final_project["status"] == "COMPLETE"
    assert final_project["revision"] == 5


def test_structure_evidence_rejects_partial_manifest_declaration_without_snapshot_changes(
    tmp_path,
):
    from ecd_epitope_scanner.experiment_projects import ExperimentStore, SnapshotError

    jobs = JobStore(tmp_path)
    job = _completed_v4_job(jobs)
    declared = _add_structure_evidence(job)
    source = Path(job["workdir"])
    (source / "run_manifest.json").write_text(json.dumps({"files": {
        "residue_reliability.csv": declared["residue_reliability.csv"],
    }}), encoding="utf-8")
    store = ExperimentStore(tmp_path, jobs)
    with pytest.raises(SnapshotError, match="only part"):
        store.create(job["id"], "Invalid evidence source")
    assert store.list() == []


def test_structure_evidence_rejects_source_hash_mismatch_without_snapshot_changes(tmp_path):
    from ecd_epitope_scanner.experiment_projects import ExperimentStore, SnapshotError

    jobs = JobStore(tmp_path)
    job = _completed_v4_job(jobs)
    _add_structure_evidence(job)
    source = Path(job["workdir"])
    with (source / "sasa_comparison.csv").open("a", encoding="utf-8") as handle:
        handle.write("tampered\n")
    store = ExperimentStore(tmp_path, jobs)
    with pytest.raises(SnapshotError, match="hash differs"):
        store.create(job["id"], "Hash mismatch evidence source")
    assert store.list() == []


def test_factorial_effect_api_reports_four_cells_and_batch_level_interaction(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api

    project_id = "factorial-fixture"
    result_path = tmp_path / "gate-analysis.json"
    samples = []
    sample_results = []
    reviews = []
    conditions = [
        ("MUTANT", "ACTIVATION", "antigen-mut", {"ab-parent": 40, "ab-variant": 60}),
        ("REFERENCE", "ACTIVATION", "antigen-wt", {"ab-parent": 70, "ab-variant": 80}),
        ("MUTANT_NC", "NC", "antigen-mut", {"ab-parent": 5, "ab-variant": 5}),
        ("REFERENCE_NC", "NC", "antigen-wt", {"ab-parent": 5, "ab-variant": 5}),
    ]
    for antibody_id in ("ab-parent", "ab-variant"):
        for index, (condition, role, antigen_id, positives) in enumerate(conditions):
            sample_id = f"{antibody_id}-{condition}"
            samples.append({
                "id": sample_id, "batch_id": "batch-1", "role": role,
                "antigen_construct_id": antigen_id, "antibody_construct_id": antibody_id,
                "replicate_type": "BIOLOGICAL", "replicate_index": index + 1,
            })
            sample_results.append({
                "sample": {"sample_id": sample_id},
                "gate_counts": [{"gate_id": "positive", "count": positives[antibody_id], "parent_count": 100}],
            })
            reviews.append({
                "id": f"review-{sample_id}", "sample_id": sample_id,
                "status": "PASS", "reason": "", "version": 1,
            })
    result_path.write_text(json.dumps({
        "gate_revision_id": "gate-revision-1", "gate_revision_version": 1,
        "samples": sample_results,
    }), encoding="utf-8")
    graph_dir = tmp_path / project_id / "reference"
    graph_dir.mkdir(parents=True)
    (graph_dir / "surface_graph_v4.json").write_text(json.dumps({
        "nodes": [{"uniprot_pos": 1, "sector_id": "VS001"}], "edges": [],
    }), encoding="utf-8")

    class FakeExperimentStore:
        projects_root = tmp_path

        def get(self, requested_project_id):
            assert requested_project_id == project_id
            return {
                "project_id": project_id, "status": "COMPLETE", "revision": 4,
                "manifest": {"reference_sequence": "ACD"},
                "frozen_hashes": {"structure": "s" * 64},
                "artifacts": [{"id": "gate-artifact", "sha256": "g" * 64}],
            }

        def list_gate_analyses(self, requested_project_id):
            assert requested_project_id == project_id
            return [{"id": "gate-analysis-1", "status": "COMPLETE", "result_artifact_id": "gate-artifact"}]

        def artifact_path(self, requested_project_id, artifact_id):
            assert requested_project_id == project_id and artifact_id == "gate-artifact"
            return result_path, {"id": artifact_id}

        def list_entities(self, requested_project_id, kind):
            assert requested_project_id == project_id
            if kind == "samples":
                return samples
            if kind == "assay_batches":
                return [{"id": "batch-1", "assay": "ACTIVATION", "independent_batch": True,
                         "cell_line": "Jurkat", "dose": 1.0, "dose_unit": "ug/mL"}]
            if kind == "antigen_constructs":
                return [
                    {"id": "antigen-wt", "sequence": "ACD", "mutations": []},
                    {"id": "antigen-mut", "sequence": "CCD", "mutations": [
                        {"kind": "SUBSTITUTION", "position": 1, "wt": "A", "mut": "C"},
                    ]},
                ]
            if kind == "antibody_constructs":
                return [
                    {"id": "ab-parent", "name": "Parent", "role": "PARENT", "parent_id": None},
                    {"id": "ab-variant", "name": "Variant", "role": "VARIANT", "parent_id": "ab-parent"},
                ]
            raise AssertionError(kind)

        def list_qc_reviews(self, requested_project_id, source_type, source_id):
            assert (requested_project_id, source_type, source_id) == (project_id, "GATE_ANALYSIS", "gate-analysis-1")
            return reviews

        def save_effect_analysis(self, *args, **kwargs):
            return {"analysis_id": "effect-analysis-1", "project_revision": 5,
                    "result_artifact": {"id": "effect-artifact"}}

    monkeypatch.setattr(api, "EXPERIMENTS", FakeExperimentStore())
    client = TestClient(api.app)
    payload = {
        "gate_analysis_id": "gate-analysis-1", "positive_gate_id": "positive",
        "comparison_mode": "DOUBLE_MUTANT_INTERACTION",
        "samples": [
            {"sample_id": sample["id"], "probe_id": "antigen-mut", "condition": condition[0]}
            for antibody_id in ("ab-parent", "ab-variant")
            for condition in conditions
            for sample in samples if sample["id"] == f"{antibody_id}-{condition[0]}"
        ],
        "expected_revision": 4,
    }
    response = client.post(
        f"/api/v4/experiment-projects/{project_id}/effect-analyses", json=payload
    )
    assert response.status_code == 201, response.text
    data = response.json()["data"]
    assert data["comparison_mode"] == "DOUBLE_MUTANT_INTERACTION"
    assert data["surface_projection"] is None
    assert data["factorial_interactions"][0]["parent_antigen_effect_pp"] == pytest.approx(-30)
    assert data["factorial_interactions"][0]["variant_antigen_effect_pp"] == pytest.approx(-20)
    assert data["factorial_interactions"][0]["interaction_pp"] == pytest.approx(10)
    assert data["factorial_interactions"][0]["four_cells"]["variant_antibody"]["mutant_antigen"] == pytest.approx(.6)
    assert data["factorial_summaries"][0]["median_interaction_pp"] == pytest.approx(10)
    assert "not assigned to individual antigen residues" in data["limitations"][1]

    fixed_antibody_payload = dict(payload)
    fixed_antibody_payload.pop("comparison_mode")
    rejected = client.post(
        f"/api/v4/experiment-projects/{project_id}/effect-analyses", json=fixed_antibody_payload
    )
    assert rejected.status_code == 422
    assert "one fixed antibody" in rejected.json()["detail"]["message"]


def test_antibody_effect_api_uses_fixed_wt_and_antibody_specific_nc(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api

    project_id = "antibody-effect-fixture"
    result_path = tmp_path / "gate-analysis.json"
    samples = []
    sample_results = []
    reviews = []
    conditions = [
        ("REFERENCE", "ACTIVATION", 80, "ab-parent"),
        ("REFERENCE_NC", "NC", 10, "ab-parent"),
        ("REFERENCE", "ACTIVATION", 60, "ab-variant"),
        ("REFERENCE_NC", "NC", 10, "ab-variant"),
    ]
    for index, (condition, role, positive, antibody_id) in enumerate(conditions):
        sample_id = f"{antibody_id}-{condition}"
        samples.append({
            "id": sample_id, "batch_id": "batch-1", "role": role,
            "antigen_construct_id": "antigen-wt", "antibody_construct_id": antibody_id,
            "replicate_type": "BIOLOGICAL", "replicate_index": index + 1,
        })
        sample_results.append({
            "sample": {"sample_id": sample_id},
            "gate_counts": [{"gate_id": "positive", "count": positive, "parent_count": 100}],
        })
        reviews.append({
            "id": f"review-{sample_id}", "sample_id": sample_id,
            "status": "PASS", "reason": "", "version": 1,
        })
    result_path.write_text(json.dumps({
        "gate_revision_id": "gate-revision-1", "gate_revision_version": 1,
        "samples": sample_results,
    }), encoding="utf-8")
    graph_dir = tmp_path / project_id / "reference"
    graph_dir.mkdir(parents=True)
    (graph_dir / "surface_graph_v4.json").write_text(json.dumps({
        "nodes": [{"uniprot_pos": 1, "sector_id": "VS001"}], "edges": [],
    }), encoding="utf-8")

    class FakeExperimentStore:
        projects_root = tmp_path

        def get(self, requested_project_id):
            assert requested_project_id == project_id
            return {
                "project_id": project_id, "status": "COMPLETE", "revision": 4,
                "manifest": {"reference_sequence": "ACD"},
                "frozen_hashes": {"structure": "s" * 64},
                "artifacts": [{"id": "gate-artifact", "sha256": "g" * 64}],
            }

        def list_gate_analyses(self, requested_project_id):
            assert requested_project_id == project_id
            return [{"id": "gate-analysis-1", "status": "COMPLETE", "result_artifact_id": "gate-artifact"}]

        def artifact_path(self, requested_project_id, artifact_id):
            assert requested_project_id == project_id and artifact_id == "gate-artifact"
            return result_path, {"id": artifact_id}

        def list_entities(self, requested_project_id, kind):
            assert requested_project_id == project_id
            if kind == "samples":
                return samples
            if kind == "assay_batches":
                return [{"id": "batch-1", "assay": "ACTIVATION", "independent_batch": True,
                         "cell_line": "Jurkat", "dose": 1.0, "dose_unit": "ug/mL"}]
            if kind == "antigen_constructs":
                return [{"id": "antigen-wt", "sequence": "ACD", "mutations": []}]
            if kind == "antibody_constructs":
                return [
                    {"id": "ab-parent", "name": "Parent", "role": "PARENT", "parent_id": None},
                    {"id": "ab-variant", "name": "Variant", "role": "VARIANT", "parent_id": "ab-parent"},
                ]
            raise AssertionError(kind)

        def list_qc_reviews(self, requested_project_id, source_type, source_id):
            assert (requested_project_id, source_type, source_id) == (
                project_id, "GATE_ANALYSIS", "gate-analysis-1"
            )
            return reviews

        def save_effect_analysis(self, requested_project_id, **kwargs):
            assert requested_project_id == project_id
            self.document = kwargs["document"]
            return {"analysis_id": "antibody-effect-1", "project_revision": 5,
                    "result_artifact": {"id": "effect-artifact"}}

    store = FakeExperimentStore()
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    client = TestClient(api.app)
    payload = {
        "gate_analysis_id": "gate-analysis-1", "positive_gate_id": "positive",
        "comparison_mode": "ANTIBODY_EFFECT",
        "samples": [
            {"sample_id": sample["id"], "probe_id": "antigen-wt", "condition": condition}
            for condition, _, _, antibody_id in conditions
            for sample in samples if sample["id"] == f"{antibody_id}-{condition}"
        ],
        "expected_revision": 4,
    }
    response = client.post(
        f"/api/v4/experiment-projects/{project_id}/effect-analyses", json=payload
    )
    assert response.status_code == 201, response.text
    data = response.json()["data"]
    assert data["comparison_mode"] == "ANTIBODY_EFFECT"
    assert data["batch_effects"] == []
    assert data["surface_projection"] is None
    assert data["antibody_effects"][0]["parent_control_window_pp"] == pytest.approx(70)
    assert data["antibody_effects"][0]["variant_control_window_pp"] == pytest.approx(50)
    assert data["antibody_effects"][0]["variant_minus_parent_pp"] == pytest.approx(-20)
    assert data["antibody_summaries"][0]["valid_batch_count"] == 1
    assert data["antibody_summaries"][0]["interval_pp"] is None
    assert store.document["antibody_effects"][0]["variant_minus_parent_pp"] == pytest.approx(
        data["antibody_effects"][0]["variant_minus_parent_pp"]
    )
    assert "surface expression" in data["limitations"][1]

    duplicate_activation = {
        **samples[2], "id": "ab-variant-REFERENCE-duplicate", "replicate_index": 5,
    }
    samples.append(duplicate_activation)
    sample_results.append({
        "sample": {"sample_id": duplicate_activation["id"]},
        "gate_counts": [{"gate_id": "positive", "count": 61, "parent_count": 100}],
    })
    reviews.append({
        "id": f"review-{duplicate_activation['id']}", "sample_id": duplicate_activation["id"],
        "status": "PASS", "reason": "", "version": 1,
    })
    result_path.write_text(json.dumps({
        "gate_revision_id": "gate-revision-1", "gate_revision_version": 1,
        "samples": sample_results,
    }), encoding="utf-8")
    incomplete = {
        **payload,
        "samples": [sample for sample in payload["samples"] if sample["sample_id"] != "ab-variant-REFERENCE_NC"]
        + [{"sample_id": duplicate_activation["id"], "probe_id": "antigen-wt", "condition": "REFERENCE"}],
    }
    rejected = client.post(
        f"/api/v4/experiment-projects/{project_id}/effect-analyses", json=incomplete
    )
    assert rejected.status_code == 422
    assert "Missing REFERENCE_NC wells" in rejected.json()["detail"]["message"]


def test_incomplete_source_does_not_create_project(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs, graph=False)
    result = TestClient(api.app).post(
        "/api/v4/experiment-projects", json={"source_job_id": job["id"], "name": "Incomplete"}
    )
    assert result.status_code == 422
    assert result.json()["detail"]["code"] == "PROJECT_SNAPSHOT_INCOMPLETE"
    assert store.list() == []


def test_job_store_preserves_live_worker_and_recovers_orphaned_owner(tmp_path, monkeypatch):
    import sqlite3

    live_pid = os.getpid()
    monkeypatch.setattr(jobs_module, "_process_alive", lambda pid: pid == live_pid)
    jobs = JobStore(tmp_path)
    job = jobs.create(PipelineRequest(input_type="uniprot", identifier="Q9NZQ7"))
    jobs.update(job["id"], status="RUNNING", stage="FETCHING")
    assert JobStore(tmp_path).get(job["id"])["status"] == "RUNNING"
    with sqlite3.connect(tmp_path / "jobs.sqlite3") as db:
        db.execute("UPDATE jobs SET owner_pid=? WHERE id=?", (99999999, job["id"]))
    assert JobStore(tmp_path).get(job["id"])["status"] == "FAILED"


def test_support_bundle_audit_restarts_and_includes_originals_only_on_request(
    monkeypatch, tmp_path
):
    import io
    import zipfile

    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    monkeypatch.setattr(api, "STORE", jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Support bundle")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id, name="Support batch", date="2026-09-21", assay="ACTIVATION",
        cell_line="Jurkat", dose=1.0, dose_unit="ug/mL", positive_gate_id="positive",
        independent_batch=True, settings={}, expected_revision=1,
    )
    source = tmp_path / "original.fcs"
    source.write_bytes(b"synthetic original flow data")
    store.import_fcs_files(
        project_id,
        [{
            "path": str(source), "original_filename": source.name,
            "metadata": {
                "event_count": 3, "channels": [{"pnn": "FSC-A"}],
                "metadata": {}, "compensation_source": "UNKNOWN_COMPENSATION_SOURCE",
            },
        }],
        batch_id=batch["batch_id"], expected_revision=2,
    )
    revision = store.get(project_id)["revision"]

    without_originals = store.build_support_bundle(project_id)
    with zipfile.ZipFile(io.BytesIO(without_originals)) as bundle:
        names = bundle.namelist()
        manifest = json.loads(bundle.read("support_bundle_manifest.json"))
        state = json.loads(bundle.read("project_state.json"))
        records = json.loads(bundle.read("database_records.json"))
        assert "project_manifest.json" in names
        assert "audit_events" in records and "project_revisions" in records
        assert state["revision"] == revision
        assert manifest["include_originals"] is False
        assert len(manifest["omitted_originals"]) == 1
        assert not any(name.endswith("original.fcs") for name in names)

    with_originals = store.build_support_bundle(project_id, include_originals=True)
    with zipfile.ZipFile(io.BytesIO(with_originals)) as bundle:
        manifest = json.loads(bundle.read("support_bundle_manifest.json"))
        raw_path = next(
            item["relative_path"]
            for item in manifest["included_artifacts"]
            if item["kind"] == "raw_fcs"
        )
        assert f"artifacts/{raw_path}" in bundle.namelist()
        assert bundle.read(f"artifacts/{raw_path}") == source.read_bytes()
        assert manifest["omitted_originals"] == []

    restarted = ExperimentStore(root, JobStore(root))
    monkeypatch.setattr(api, "EXPERIMENTS", restarted)
    response = TestClient(api.app).get(
        f"/api/v4/experiment-projects/{project_id}/support-bundle"
    )
    assert response.status_code == 200
    assert response.headers["content-type"].startswith("application/zip")
    assert response.content == without_originals
    included_response = TestClient(api.app).get(
        f"/api/v4/experiment-projects/{project_id}/support-bundle?include_originals=true"
    )
    assert included_response.status_code == 200
    with zipfile.ZipFile(io.BytesIO(included_response.content)) as bundle:
        included_manifest = json.loads(bundle.read("support_bundle_manifest.json"))
        included_raw_path = next(
            item["relative_path"]
            for item in included_manifest["included_artifacts"]
            if item["kind"] == "raw_fcs"
        )
        assert bundle.read(f"artifacts/{included_raw_path}") == source.read_bytes()
    assert restarted.get(project_id)["revision"] == revision


def test_constructs_and_batches_are_versioned(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Fixture")
    project_id = project["project_id"]
    client = TestClient(api.app)
    base = f"/api/v4/experiment-projects/{project_id}"

    antigen = client.post(
        f"{base}/antigen-constructs",
        json={
            "name": "A10C",
            "mutations": [{"position": 10, "wt": "A", "mut": "C"}],
            "expected_revision": 1,
        },
    )
    assert antigen.status_code == 201, antigen.text
    antigen_data = antigen.json()["data"]
    assert antigen_data["mapping_status"] == "MAPPED"
    assert antigen_data["project_revision"] == 2
    duplicate = client.post(
        f"{base}/antigen-constructs",
        json={
            "name": "stale",
            "mutations": [{"position": 10, "wt": "A", "mut": "C"}],
            "expected_revision": 1,
        },
    )
    assert duplicate.status_code == 409
    wrong_wt = client.post(
        f"{base}/antigen-constructs",
        json={
            "name": "bad",
            "mutations": [{"position": 10, "wt": "C", "mut": "A"}],
            "expected_revision": 2,
        },
    )
    assert wrong_wt.status_code == 422

    parent = client.post(
        f"{base}/antibody-constructs",
        json={
            "name": "parent",
            "role": "PARENT",
            "vh_sequence": "ACD",
            "vl_sequence": "EFG",
            "expected_revision": 2,
        },
    )
    assert parent.status_code == 201, parent.text
    parent_id = parent.json()["data"]["construct_id"]
    variant = client.post(
        f"{base}/antibody-constructs",
        json={
            "name": "variant",
            "role": "VARIANT",
            "parent_construct_id": parent_id,
            "vh_sequence": "AAD",
            "vl_sequence": "EFG",
            "expected_revision": 3,
        },
    )
    assert variant.status_code == 201, variant.text
    assert variant.json()["data"]["mutation_class"] == "SINGLE"

    batch = client.post(
        f"{base}/batches",
        json={
            "name": "Myc 260921",
            "date": "2026-09-21",
            "assay": "ANTIGEN_MYC",
            "cell_line": "unknown",
            "expected_revision": 4,
        },
    )
    assert batch.status_code == 201, batch.text
    assert client.get(f"{base}").json()["data"]["revision"] == 5
    assert len(client.get(f"{base}/antigen-constructs").json()["data"]) == 1
    assert len(client.get(f"{base}/antibody-constructs").json()["data"]) == 2
    assert len(client.get(f"{base}/batches").json()["data"]) == 1
    invalid_summary = client.post(f"{base}/summary-measurements", json={
        "assay": "ANTIGEN_MYC", "construct_id": antigen_data["construct_id"],
        "measured_on": "2026-09-21", "metric": "POSITIVE_FRACTION", "value": 1.2,
        "unit": "fraction", "source_label": "historical notebook", "reviewer": "operator",
        "expected_revision": 5,
    })
    assert invalid_summary.status_code == 422
    summary = client.post(f"{base}/summary-measurements", json={
        "assay": "ANTIGEN_MYC", "construct_id": antigen_data["construct_id"],
        "measured_on": "2026-09-21", "metric": "POSITIVE_FRACTION", "value": 0.72,
        "unit": "fraction", "source_label": "historical notebook", "reviewer": "operator",
        "notes": "No paired activation batch", "expected_revision": 5,
    })
    assert summary.status_code == 201, summary.text
    assert summary.json()["data"]["interpretation"] == (
        "HISTORICAL_SUMMARY_ONLY_NOT_CONCURRENT_CONTROL"
    )
    summaries = client.get(f"{base}/summary-measurements").json()["data"]
    assert summaries[0]["historical"] is True
    assert summaries[0]["construct_id"] == antigen_data["construct_id"]
    restarted = ExperimentStore(tmp_path, JobStore(tmp_path))
    assert restarted.get(project_id)["revision"] == 6
    assert restarted.list_summary_measurements(project_id)[0]["value"] == pytest.approx(0.72)


def test_frozen_recommendation_can_be_added_as_antigen_construct(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Fixture")
    base = f"/api/v4/experiment-projects/{project['project_id']}"
    client = TestClient(api.app)

    recommendations = client.get(f"{base}/recommended-antigens")
    assert recommendations.status_code == 200, recommendations.text
    assert recommendations.json()["data"][0]["id"] == "V4_01"
    added = client.post(
        f"{base}/antigen-constructs",
        json={
            "name": "Frozen V4 recommendation",
            "source": "RECOMMENDED",
            "recommended_construct_id": "V4_01",
            "expected_revision": 1,
        },
    )
    assert added.status_code == 201, added.text
    record = added.json()["data"]
    assert record["source"] == "RECOMMENDED"
    assert record["recommended_construct_id"] == "V4_01"
    assert record["sequence"] == "CCD"
    assert record["sequence_sha256"]
    assert store.get(project["project_id"])["revision"] == 2

    invalid = client.post(
        f"{base}/antigen-constructs",
        json={"name": "Mismatch", "source": "RECOMMENDED", "expected_revision": 2},
    )
    assert invalid.status_code == 422


def test_full_sequence_preview_is_non_mutating_and_matches_saved_construct(
    monkeypatch, tmp_path
):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    project = store.create(_completed_v4_job(jobs)["id"], "Full sequence preview")
    project_id = project["project_id"]
    base = f"/api/v4/experiment-projects/{project_id}"
    client = TestClient(api.app)

    preview = client.post(
        f"{base}/antigen-constructs/preview", json={"full_sequence": "AAD"}
    )
    assert preview.status_code == 200, preview.text
    preview_data = preview.json()["data"]
    assert preview_data["mutations"] == [{
        "kind": "SUBSTITUTION", "position": 11, "wt": "C", "mut": "A", "mapped": True,
    }]
    assert store.get(project_id)["revision"] == 1

    saved = client.post(f"{base}/antigen-constructs", json={
        "name": "full sequence C11A", "source": "CUSTOM", "full_sequence": "AAD",
        "expected_revision": 1,
    })
    assert saved.status_code == 201, saved.text
    saved_data = saved.json()["data"]
    assert saved_data["mutations"] == preview_data["mutations"]
    assert saved_data["sequence_sha256"] == preview_data["sequence_sha256"]
    assert store.get(project_id)["revision"] == 2

    identical = client.post(
        f"{base}/antigen-constructs/preview", json={"full_sequence": "ACD"}
    )
    assert identical.status_code == 200
    assert identical.json()["data"]["mutations"] == []
    assert store.get(project_id)["revision"] == 2


def test_fcs_import_is_content_addressed_and_match_is_a_revision(tmp_path):
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "FCS import")
    batch = store.add_batch(
        project["project_id"],
        name="Myc",
        date="2026-09-21",
        assay="ANTIGEN_MYC",
        cell_line="unspecified",
        dose=None,
        dose_unit=None,
        positive_gate_id=None,
        independent_batch=True,
        settings={},
        expected_revision=1,
    )
    first = tmp_path / "first.fcs"
    duplicate = tmp_path / "copy.fcs"
    first.write_bytes(b"same fcs bytes")
    duplicate.write_bytes(b"same fcs bytes")
    metadata = {
        "event_count": 42,
        "channels": [{"pnn": "RL1-A"}],
        "metadata": {},
        "compensation_source": "UNKNOWN_COMPENSATION_SOURCE",
    }
    imported = store.import_fcs_files(
        project["project_id"],
        [
            {"path": str(first), "original_filename": "H1.fcs", "metadata": metadata},
            {"path": str(duplicate), "original_filename": "H1-copy.fcs", "metadata": metadata},
        ],
        batch_id=batch["batch_id"],
        expected_revision=2,
    )
    assert imported["project_revision"] == 3
    assert imported["files"][0]["artifact_id"] == imported["files"][1]["artifact_id"]
    assert (
        len(
            [
                item
                for item in store.get(project["project_id"])["artifacts"]
                if item["kind"] == "raw_fcs"
            ]
        )
        == 1
    )
    confirmed = store.confirm_fcs_matches(
        project["project_id"],
        imported["upload_session_id"],
        [
            {
                "upload_id": imported["files"][0]["upload_id"],
                "role": "ANTIGEN_MYC",
                "batch_id": batch["batch_id"],
                "replicate_type": "SINGLE",
                "replicate_index": 1,
                "condition": {"sample_label": "H1"},
            }
        ],
        expected_revision=3,
    )
    assert confirmed["status"] == "PARTIAL"
    assert len(confirmed["sample_ids"]) == 1
    assert len(confirmed["unmatched_upload_ids"]) == 1
    samples = store.list_entities(project["project_id"], "samples")
    assert len(samples) == 1
    assert samples[0]["fcs_metadata"]["event_count"] == 42
    assert samples[0]["artifact_id"] == imported["files"][0]["artifact_id"]
    assert store.get(project["project_id"])["status"] == "COMPLETE"


def test_fcs_import_reuses_verified_readonly_content_after_hard_crash(tmp_path):
    import sqlite3
    import stat

    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "FCS artifact publish recovery")
    project_id = project["project_id"]
    source = tmp_path / "crash-recovery.fcs"
    source.write_bytes(b"crash-recovery raw FCS content")
    metadata = {
        "event_count": 12,
        "channels": [{"pnn": "GFP-A"}],
        "metadata": {},
        "compensation_source": "USER_CONFIRMED_NOT_APPLIED",
    }

    worker = multiprocessing.get_context("spawn").Process(
        target=_crash_fcs_import_after_artifact_is_read_only,
        args=(str(root), project_id, str(source), metadata),
    )
    worker.start()
    worker.join(timeout=30)
    if worker.is_alive():
        worker.terminate()
        worker.join(timeout=5)
        pytest.fail("FCS import crash worker did not terminate")
    assert worker.exitcode == 80

    digest = sha256_file(source)
    target = root / "experiments" / project_id / "raw" / "fcs" / f"{digest}.fcs"
    assert target.is_file()
    assert sha256_file(target) == digest
    assert not target.stat().st_mode & stat.S_IWRITE
    assert all(item["kind"] != "raw_fcs" for item in store.get(project_id)["artifacts"])

    restarted = ExperimentStore(root, JobStore(root))
    assert restarted.get(project_id)["status"] == "COMPLETE"
    assert restarted.get(project_id)["revision"] == 1
    imported = restarted.import_fcs_files(
        project_id,
        [{"path": str(source), "original_filename": "crash-recovery.fcs", "metadata": metadata}],
        batch_id=None,
        expected_revision=1,
    )
    assert imported["project_revision"] == 2
    assert imported["files"][0]["sha256"] == digest
    assert imported["files"][0]["event_count"] == 12
    assert sha256_file(target) == digest
    with sqlite3.connect(root / "experiments.sqlite3") as db:
        assert db.execute(
            "SELECT count(*) FROM artifacts WHERE project_id=? AND kind='raw_fcs'",
            (project_id,),
        ).fetchone()[0] == 1
        assert db.execute(
            "SELECT count(*) FROM import_sessions WHERE project_id=?",
            (project_id,),
        ).fetchone()[0] == 1
    assert restarted.get(project_id)["status"] == "COMPLETE"


@pytest.mark.parametrize(
    ("crash_point", "exit_code", "expected_manifest_repairs"),
    [("artifact_published", 81, 0), ("manifest_published", 82, 1)],
)
def test_compensation_matrix_registration_recovers_after_hard_crash(
    tmp_path, crash_point, exit_code, expected_manifest_repairs
):
    import hashlib
    import json
    import sqlite3
    import stat

    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    jobs = JobStore(root)
    store = ExperimentStore(root, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Compensation matrix publish recovery")
    project_id = project["project_id"]
    matrix_document = {
        "schema_version": "1.0",
        "source_filename": "synthetic-spill.csv",
        "import_format": "csv",
        "detectors": ["FL1-A", "FL2-A"],
        "fluorochromes": ["FITC-A", "PE-A"],
        "matrix": [[1.0, 0.08], [0.11, 1.0]],
    }
    payload = json.dumps(matrix_document, ensure_ascii=False, indent=2).encode("utf-8")
    digest = hashlib.sha256(payload).hexdigest()

    worker = multiprocessing.get_context("spawn").Process(
        target=_crash_compensation_matrix_registration,
        args=(str(root), project_id, matrix_document, crash_point),
    )
    worker.start()
    worker.join(timeout=30)
    if worker.is_alive():
        worker.terminate()
        worker.join(timeout=5)
        pytest.fail("Compensation matrix crash worker did not terminate")
    assert worker.exitcode == exit_code

    target = (
        root / "experiments" / project_id / "compensation" / "matrices" / f"{digest}.json"
    )
    assert target.is_file()
    assert sha256_file(target) == digest
    assert not target.stat().st_mode & stat.S_IWRITE

    restarted = ExperimentStore(root, JobStore(root))
    recovered = restarted.get(project_id)
    assert recovered["status"] == "COMPLETE"
    assert recovered["revision"] == 1
    assert all(item["kind"] != "compensation_matrix" for item in recovered["artifacts"])
    with sqlite3.connect(root / "experiments.sqlite3") as db:
        repair_count = db.execute(
            "SELECT count(*) FROM audit_events WHERE project_id=? AND action='MANIFEST_RECOVERED'",
            (project_id,),
        ).fetchone()[0]
    assert repair_count == expected_manifest_repairs

    registered = restarted.register_compensation_matrix(
        project_id, matrix_document, expected_revision=1
    )
    assert registered["project_revision"] == 2
    assert registered["artifact"]["sha256"] == digest
    assert restarted.get(project_id)["status"] == "COMPLETE"
    assert sum(
        item["kind"] == "compensation_matrix" for item in restarted.get(project_id)["artifacts"]
    ) == 1


@pytest.mark.parametrize(
    ("crash_point", "exit_code", "expected_manifest_repairs"),
    [("artifact_published", 83, 0), ("manifest_published", 84, 1)],
)
def test_gate_analysis_artifact_recovers_after_hard_crash(
    tmp_path, crash_point, exit_code, expected_manifest_repairs
):
    import sqlite3
    import stat

    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    fixture = _gate_artifact_recovery_fixture(root)
    project_id = fixture["project_id"]
    process = multiprocessing.get_context("spawn").Process(
        target=_crash_gate_analysis_save,
        args=(
            str(root), project_id, fixture["gate_revision_id"],
            fixture["input_artifact_ids"], fixture["worker_result"],
            fixture["expected_revision"], crash_point,
        ),
    )
    process.start()
    process.join(timeout=30)
    if process.is_alive():
        process.terminate()
        process.join(timeout=5)
        pytest.fail("Gate analysis crash worker did not terminate")
    assert process.exitcode == exit_code

    output_root = root / "experiments" / project_id / "analyses" / "gate_runs"
    orphaned_files = list(output_root.rglob("gate_counts.json"))
    assert len(orphaned_files) == 1
    assert not orphaned_files[0].stat().st_mode & stat.S_IWRITE

    restarted = ExperimentStore(root, JobStore(root))
    assert restarted.get(project_id)["status"] == "COMPLETE"
    assert restarted.get(project_id)["revision"] == fixture["expected_revision"]
    assert restarted.list_gate_analyses(project_id) == []
    assert all(item["kind"] != "gate_analysis" for item in restarted.get(project_id)["artifacts"])
    with sqlite3.connect(root / "experiments.sqlite3") as db:
        repair_count = db.execute(
            "SELECT count(*) FROM audit_events WHERE project_id=? AND action='MANIFEST_RECOVERED'",
            (project_id,),
        ).fetchone()[0]
    assert repair_count == expected_manifest_repairs

    retry = restarted.save_gate_analysis(
        project_id, gate_revision_id=fixture["gate_revision_id"],
        input_artifact_ids=fixture["input_artifact_ids"],
        worker_result=fixture["worker_result"],
        expected_revision=fixture["expected_revision"],
    )
    assert retry["project_revision"] == fixture["expected_revision"] + 1
    assert len(restarted.list_gate_analyses(project_id)) == 1
    result_path, _ = restarted.artifact_path(
        project_id, retry["result_artifact"]["id"]
    )
    assert json.loads(result_path.read_text(encoding="utf-8"))["status"] == "COMPLETE"
    assert restarted.get(project_id)["status"] == "COMPLETE"


@pytest.mark.parametrize(
    ("crash_point", "exit_code", "expected_manifest_repairs"),
    [("artifact_published", 85, 0), ("manifest_published", 86, 1)],
)
def test_effect_analysis_artifact_recovers_after_hard_crash(
    tmp_path, crash_point, exit_code, expected_manifest_repairs
):
    import sqlite3
    import stat

    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    root = tmp_path / "workspace"
    fixture = _gate_artifact_recovery_fixture(root)
    store = fixture["store"]
    project_id = fixture["project_id"]
    gate_analysis = store.save_gate_analysis(
        project_id, gate_revision_id=fixture["gate_revision_id"],
        input_artifact_ids=fixture["input_artifact_ids"],
        worker_result=fixture["worker_result"],
        expected_revision=fixture["expected_revision"],
    )
    expected_revision = gate_analysis["project_revision"]

    process = multiprocessing.get_context("spawn").Process(
        target=_crash_effect_analysis_save,
        args=(str(root), project_id, gate_analysis["analysis_id"], expected_revision, crash_point),
    )
    process.start()
    process.join(timeout=30)
    if process.is_alive():
        process.terminate()
        process.join(timeout=5)
        pytest.fail("Effect analysis crash worker did not terminate")
    assert process.exitcode == exit_code

    output_root = root / "experiments" / project_id / "analyses" / "effect_runs"
    orphaned_files = list(output_root.rglob("effect_results.json"))
    assert len(orphaned_files) == 1
    assert not orphaned_files[0].stat().st_mode & stat.S_IWRITE

    restarted = ExperimentStore(root, JobStore(root))
    assert restarted.get(project_id)["status"] == "COMPLETE"
    assert restarted.get(project_id)["revision"] == expected_revision
    assert restarted.list_effect_analyses(project_id) == []
    assert all(item["kind"] != "effect_analysis" for item in restarted.get(project_id)["artifacts"])
    with sqlite3.connect(root / "experiments.sqlite3") as db:
        repair_count = db.execute(
            "SELECT count(*) FROM audit_events WHERE project_id=? AND action='MANIFEST_RECOVERED'",
            (project_id,),
        ).fetchone()[0]
    assert repair_count == expected_manifest_repairs

    retry = restarted.save_effect_analysis(
        project_id, gate_analysis_id=gate_analysis["analysis_id"],
        positive_gate_id="positive", document={"synthetic_storage_recovery_probe": True},
        expected_revision=expected_revision,
    )
    assert retry["project_revision"] == expected_revision + 1
    assert len(restarted.list_effect_analyses(project_id)) == 1
    result_path, _ = restarted.artifact_path(project_id, retry["result_artifact"]["id"])
    assert json.loads(result_path.read_text(encoding="utf-8"))[
        "synthetic_storage_recovery_probe"
    ] is True
    assert restarted.get(project_id)["status"] == "COMPLETE"


def test_activation_samples_can_be_registered_with_explicit_temporary_labels(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Temporary activation labels")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id,
        name="Activation label-only batch",
        date="2026-09-28",
        assay="ACTIVATION",
        cell_line="Jurkat",
        dose=None,
        dose_unit=None,
        positive_gate_id="positive",
        independent_batch=True,
        settings={},
        expected_revision=1,
    )
    files = []
    for name in ("WT.fcs", "NC.fcs"):
        path = tmp_path / name
        path.write_bytes(name.encode("ascii"))
        files.append({
            "path": str(path),
            "original_filename": name,
            "metadata": {
                "event_count": 100,
                "channels": [{"pnn": "FSC-A"}],
                "metadata": {"fil": name},
                "compensation_source": "UNKNOWN_COMPENSATION_SOURCE",
            },
        })
    imported = store.import_fcs_files(
        project_id, files, batch_id=batch["batch_id"], expected_revision=2
    )
    matches = [
        {
            "upload_id": imported["files"][0]["upload_id"],
            "role": "ACTIVATION",
            "replicate_type": "TECHNICAL",
            "replicate_index": 1,
            "temporary_identity": {
                "status": "TEMPORARY_LABEL",
                "source": "USER_CONFIRMED_LABELS",
                "antigen_label": "CTLA4-WT",
                "antigen_role": "REFERENCE",
                "panel_label": "L32",
                "panel_role": "ANTIBODY",
            },
            "condition": {"sample_label": "WT · L32"},
        },
        {
            "upload_id": imported["files"][1]["upload_id"],
            "role": "NC",
            "replicate_type": "TECHNICAL",
            "replicate_index": 1,
            "temporary_identity": {
                "status": "TEMPORARY_LABEL",
                "source": "USER_CONFIRMED_LABELS",
                "antigen_label": None,
                "antigen_role": "NONE",
                "panel_label": "L32",
                "panel_role": "ANTIBODY",
            },
            "condition": {"sample_label": "NC · L32"},
        },
    ]
    with pytest.raises(ValueError, match="require antigen and antibody constructs"):
        store.confirm_fcs_matches(
            project_id,
            imported["upload_session_id"],
            [{key: value for key, value in matches[0].items() if key != "temporary_identity"}],
            expected_revision=3,
        )
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    response = TestClient(api.app).post(
        f"/api/v4/experiment-projects/{project_id}/samples/confirm",
        json={
            "upload_session_id": imported["upload_session_id"],
            "matches": matches,
            "expected_revision": 3,
        },
    )
    assert response.status_code == 201, response.text
    confirmed = response.json()["data"]
    samples = store.list_entities(project_id, "samples")
    assert confirmed["status"] == "CONFIRMED"
    assert len(samples) == 2
    assert samples[0]["antigen_construct_id"] is None
    assert samples[0]["antibody_construct_id"] is None
    assert samples[0]["condition"]["temporary_identity"]["antigen_label"] == "CTLA4-WT"
    assert store.list_entities(project_id, "antigen_constructs") == []
    assert store.list_entities(project_id, "antibody_constructs") == []
    batch_record = store.list_entities(project_id, "assay_batches")[0]
    assert "TEMPORARY_IDENTITIES_UNVERIFIED" in batch_record["design_warnings"]
    assert "WT_IDENTITY_UNVERIFIED" in batch_record["design_warnings"]
    assert "CONSTRUCT_LINK_MISSING" not in batch_record["design_warnings"]


def test_gate_revisions_are_validated_immutable_and_audited(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Gate fixture")
    client = TestClient(api.app)
    url = f"/api/v4/experiment-projects/{project['project_id']}/gates"
    gate_tree = [
        {
            "gate_id": "Live",
            "kind": "POLYGON_2D",
            "channels": ["FSC-A", "SSC-A"],
            "transform_id": "scatter-logicle-v1",
            "parameters": {"vertices": [[0, 0], [2, 0], [2, 2], [0, 2]]},
        }
    ]
    transforms = [{
        "transform_id": "scatter-logicle-v1", "kind": "LOGICLE",
        "parameters": {"param_t": 262144, "param_w": 0.5, "param_m": 4.5, "param_a": 0},
    }]
    missing_decision = client.post(
        url,
        json={
            "name": "initial",
            "gate_tree": gate_tree,
            "transforms": transforms,
            "compensation_decision": "UNKNOWN",
            "expected_revision": 1,
        },
    )
    assert missing_decision.status_code == 422
    assert missing_decision.json()["detail"]["code"] == "UNKNOWN_COMPENSATION_SOURCE"
    saved = client.post(
        url,
        json={
            "name": "initial",
            "gate_tree": gate_tree,
            "transforms": transforms,
            "compensation_decision": "NONE",
            "expected_revision": 1,
        },
    )
    assert saved.status_code == 201, saved.text
    assert saved.json()["data"]["version"] == 1
    listed = client.get(url).json()["data"]
    assert len(listed) == 1
    assert listed[0]["transforms"] == transforms
    stale = client.post(
        url,
        json={
            "name": "stale",
            "gate_tree": gate_tree,
            "transforms": transforms,
            "compensation_decision": "NONE",
            "expected_revision": 1,
        },
    )
    assert stale.status_code == 409


def test_gate_analysis_api_records_complete_event_counts_and_inputs(monkeypatch, tmp_path):
    from types import SimpleNamespace

    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Gate analysis fixture")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id, name="Fixture", date="2026-09-21", assay="ANTIGEN_MYC",
        cell_line="unspecified", dose=None, dose_unit=None, positive_gate_id=None,
        independent_batch=True, settings={}, expected_revision=1,
    )
    fcs = tmp_path / "H1.fcs"
    fcs.write_bytes(b"registered test FCS")
    imported = store.import_fcs_files(
        project_id,
        [{"path": str(fcs), "original_filename": "H1.fcs", "metadata": {
            "event_count": 42, "channels": [{"pnn": "FSC-A"}, {"pnn": "SSC-A"}],
            "metadata": {}, "compensation_source": "UNKNOWN_COMPENSATION_SOURCE",
        }}],
        batch_id=batch["batch_id"], expected_revision=2,
    )
    upload = imported["files"][0]
    client = TestClient(api.app)
    base = f"/api/v4/experiment-projects/{project_id}"
    gate = client.post(f"{base}/gates", json={
        "name": "Uncompensated preview", "compensation_decision": "NONE",
        "expected_revision": 3,
        "gate_tree": [{
            "gate_id": "cells", "name": "Cells", "parent_gate_id": "root",
            "kind": "RECTANGLE_2D", "channels": ["FSC-A", "SSC-A"],
            "parameters": {"bounds": [{"min": 0, "max": 100}, {"min": 0, "max": 100}]},
        }],
    })
    assert gate.status_code == 201, gate.text
    gate_id = gate.json()["data"]["gate_revision_id"]
    unmatched = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": gate_id, "fcs_artifact_ids": [upload["artifact_id"]],
        "expected_revision": 4,
    })
    assert unmatched.status_code == 422
    assert unmatched.json()["detail"]["code"] == "SAMPLE_MAPPING_REQUIRED"

    confirmed = store.confirm_fcs_matches(
        project_id, imported["upload_session_id"], [{
            "upload_id": upload["upload_id"], "batch_id": batch["batch_id"], "role": "NC",
            "replicate_type": "SINGLE", "replicate_index": 1,
            "condition": {"sample_label": "H1"},
        }], expected_revision=4,
    )
    assert confirmed["project_revision"] == 5
    plot_worker_result = {"schema_version": "1.0", "status": "COMPLETE", "records": {
        "file_name": "H1.fcs", "sample_id": "H1.fcs", "event_count": 42,
        "finite_event_count": 41, "omitted_nonfinite_count": 1, "sampled_event_count": 2,
        "maximum_points": 100, "sampling_method": "DETERMINISTIC_EVEN_INDEX",
        "display_only": True, "source": "RAW_UNTRANSFORMED",
        "x_channel": "FSC-A", "y_channel": "SSC-A",
        "domain": {"x": [1.0, 40.0], "y": [2.0, 39.0], "quantiles": [0.01, 0.99]},
        "points": [[4.0, 6.0], [32.0, 31.0]],
    }, "provenance": {"flowkit_version": "1.3.2", "worker": "flowkit_worker"}}
    plot_calls = []

    def sample_plot_worker(command, **kwargs):
        plot_calls.append((command, json.loads(kwargs["input"])))
        return SimpleNamespace(returncode=0, stdout=json.dumps(plot_worker_result), stderr="")

    monkeypatch.setattr(api.subprocess, "run", sample_plot_worker)
    sample_id = confirmed["sample_ids"][0]
    plot = client.post(f"{base}/samples/{sample_id}/plot-data", json={
        "x_channel": "FSC-A", "y_channel": "SSC-A", "maximum_points": 100,
    })
    assert plot.status_code == 200, plot.text
    plot_data = plot.json()["data"]
    assert plot_data["display_only"] is True
    assert plot_data["source"] == "RAW_UNTRANSFORMED"
    assert plot_data["event_count"] == 42
    assert plot_data["points"] == [[4.0, 6.0], [32.0, 31.0]]
    registered_path, _ = store.artifact_path(project_id, upload["artifact_id"])
    assert plot_calls[0][1]["path"] == str(registered_path)
    assert plot_calls[0][1]["operation"] == "sample_events"
    unknown_plot_channel = client.post(f"{base}/samples/{sample_id}/plot-data", json={
        "x_channel": "arbitrary-path-or-channel", "y_channel": "SSC-A",
    })
    assert unknown_plot_channel.status_code == 422
    assert len(plot_calls) == 1

    gate = client.post(f"{base}/gates", json={
        "name": "Uncompensated preview rerun", "compensation_decision": "NONE",
        "parent_gate_revision_id": gate_id, "expected_revision": 5,
        "transforms": [{
            "transform_id": "fsc-logicle-v1", "kind": "LOGICLE",
            "parameters": {"param_t": 262144, "param_w": 0.5, "param_m": 4.5, "param_a": 0},
        }],
        "gate_tree": [{
            "gate_id": "cells", "name": "Cells", "parent_gate_id": "root",
            "kind": "RECTANGLE_2D", "channels": ["FSC-A", "SSC-A"],
            "transform_id": "fsc-logicle-v1",
            "parameters": {"bounds": [{"min": 0, "max": 100}, {"min": 0, "max": 100}]},
        }],
    })
    assert gate.status_code == 201, gate.text
    gate_id = gate.json()["data"]["gate_revision_id"]
    worker_result = {"schema_version": "1.0", "status": "COMPLETE", "records": {
        "samples": [{"file_name": "H1.fcs", "sample_id": "H1.fcs", "event_count": 42,
                     "gate_counts": [{"gate_id": "cells", "gate_name": "Cells",
                                      "gate_path": ["root"], "parent_gate": "root",
                                      "parent_count": 42, "count": 40,
                                      "absolute_fraction": 40 / 42, "parent_fraction": 40 / 42}]}]
    }, "provenance": {"flowkit_version": "1.3.2", "worker": "flowkit_worker"}}
    monkeypatch.setattr(api.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(
        returncode=0, stdout=json.dumps(worker_result), stderr=""
    ))
    result = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": gate_id, "fcs_artifact_ids": [upload["artifact_id"]],
        "expected_revision": 6,
    })
    assert result.status_code == 201, result.text
    record = result.json()["data"]
    assert record["project_revision"] == 7
    assert record["samples"][0]["sample"]["role"] == "NC"
    assert record["samples"][0]["gate_counts"][0]["count"] == 40

    saved = client.get(f"{base}/files/{record['result_artifact']['id']}")
    assert saved.status_code == 200
    document = saved.json()
    assert document["gate_revision_id"] == gate_id
    assert document["transforms"][0]["transform_id"] == "fsc-logicle-v1"
    assert document["samples"][0]["input_sha256"] == upload["sha256"]
    assert document["samples"][0]["sample"]["condition"]["sample_label"] == "H1"
    assert client.get(f"{base}/gate-analyses").json()["data"][0]["status"] == "COMPLETE"
    sample_id = record["samples"][0]["sample"]["sample_id"]
    override_without_reason = client.post(f"{base}/qc-reviews", json={
        "gate_analysis_id": record["analysis_id"], "sample_id": sample_id,
        "status": "PASS_WITH_OVERRIDE", "reviewer": "operator", "expected_revision": 7,
    })
    assert override_without_reason.status_code == 422
    override = client.post(f"{base}/qc-reviews", json={
        "gate_analysis_id": record["analysis_id"], "sample_id": sample_id,
        "status": "PASS_WITH_OVERRIDE", "reason": "No same-day conformation control",
        "reviewer": "operator", "expected_revision": 7,
    })
    assert override.status_code == 201, override.text
    assert override.json()["data"]["project_revision"] == 8
    passed = client.post(f"{base}/qc-reviews", json={
        "gate_analysis_id": record["analysis_id"], "sample_id": sample_id,
        "status": "PASS", "reviewer": "operator", "expected_revision": 8,
    })
    assert passed.status_code == 201, passed.text
    assert passed.json()["data"]["version"] == 2
    reviews = client.get(
        f"{base}/qc-reviews?source_type=GATE_ANALYSIS&source_id={record['analysis_id']}"
    )
    assert reviews.status_code == 200
    assert [item["version"] for item in reviews.json()["data"]] == [1, 2]
    assert reviews.json()["data"][0]["reason"] == "No same-day conformation control"

    workspace = tmp_path / "unresolved.wsp"
    workspace.write_text("<Workspace />", encoding="utf-8")
    unresolved_node = {
        "sample_id": "H1.fcs", "name": "Lymphocytes", "parent_path": ["root"],
        "source_gate_type": "EllipsoidGate", "kind": None,
        "channels": ["FSC-A", "SSC-A"], "status": "UNRESOLVED_FLOWJO_NODE",
    }
    flowjo_import = store.register_flowjo_import(
        project_id, workspace, "unresolved.wsp", {
            "schema_version": "1.0", "status": "COMPLETE",
            "records": {
                "status": "WAITING_FOR_USER", "samples": [], "gate_counts": [],
                "unresolved_nodes": [unresolved_node],
            },
            "flowjo_input_artifact_ids": [upload["artifact_id"]],
        }, expected_revision=store.get(project_id)["revision"],
    )
    assert flowjo_import["status"] == "WAITING_FOR_USER"
    rebuild_tree = [{
        "gate_id": "lymphocytes", "name": "Lymphocytes", "parent_gate_id": "root",
        "kind": "RECTANGLE_2D", "channels": ["FSC-A", "SSC-A"],
        "parameters": {"bounds": [{"min": 0, "max": 100}, {"min": 0, "max": 100}]},
    }]
    rebuilt = client.post(f"{base}/gates", json={
        "name": "Manual FlowJo reconstruction", "compensation_decision": "NONE",
        "expected_revision": flowjo_import["project_revision"], "gate_tree": rebuild_tree,
        "flowjo_rebuild": {
            "flowjo_import_id": flowjo_import["flowjo_import_id"],
            "reviewer": "operator", "reason": "Recreated the reviewed node as a rectangular gate",
        },
    })
    assert rebuilt.status_code == 201, rebuilt.text
    rebuild_revision = rebuilt.json()["data"]
    assert rebuild_revision["flowjo_rebuild"]["unresolved_node_count"] == 1
    assert rebuild_revision["flowjo_rebuild"]["input_artifact_ids"] == [upload["artifact_id"]]
    assert store.list_gate_revisions(project_id)[-1]["flowjo_rebuild"]["reviewer"] == "operator"
    with pytest.raises(ValueError, match="same registered FCS artifacts"):
        store.validate_gate_analysis_inputs(
            project_id, rebuild_revision["gate_revision_id"], ["different-artifact-id"]
        )
    rebuilt_analysis = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": rebuild_revision["gate_revision_id"],
        "fcs_artifact_ids": [upload["artifact_id"]],
        "expected_revision": rebuild_revision["project_revision"],
    })
    assert rebuilt_analysis.status_code == 201, rebuilt_analysis.text
    rebuilt_artifact = client.get(
        f"{base}/files/{rebuilt_analysis.json()['data']['result_artifact']['id']}"
    ).json()
    assert rebuilt_artifact["flowjo_rebuild"]["flowjo_import_id"] == flowjo_import["flowjo_import_id"]
    assert rebuilt_artifact["flowjo_rebuild"]["reason"] == "Recreated the reviewed node as a rectangular gate"
    assert ExperimentStore(tmp_path, JobStore(tmp_path)).get(project_id)["status"] == "COMPLETE"


def test_external_compensation_matrix_upload_gate_binding_and_result_provenance(monkeypatch, tmp_path):
    from types import SimpleNamespace

    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore, RevisionConflict

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    project = store.create(_completed_v4_job(jobs)["id"], "External matrix integration")
    project_id = project["project_id"]
    base = f"/api/v4/experiment-projects/{project_id}"
    batch = store.add_batch(
        project_id, name="Matrix fixture", date="2026-09-21", assay="ANTIGEN_MYC",
        cell_line="unspecified", dose=None, dose_unit=None, positive_gate_id=None,
        independent_batch=True, settings={}, expected_revision=1,
    )
    fcs = tmp_path / "H1.fcs"
    fcs.write_bytes(b"registered test FCS")
    imported = store.import_fcs_files(
        project_id,
        [{"path": str(fcs), "original_filename": "H1.fcs", "metadata": {
            "event_count": 42,
            "channels": [{"pnn": label} for label in ("FSC-A", "SSC-A", "RL1-A", "YL2-A")],
            "metadata": {}, "compensation_source": "UNKNOWN_COMPENSATION_SOURCE",
        }}],
        batch_id=batch["batch_id"], expected_revision=2,
    )
    upload = imported["files"][0]
    confirmed = store.confirm_fcs_matches(
        project_id, imported["upload_session_id"], [{
            "upload_id": upload["upload_id"], "batch_id": batch["batch_id"], "role": "NC",
            "replicate_type": "SINGLE", "replicate_index": 1,
            "condition": {"sample_label": "H1"},
        }], expected_revision=3,
    )
    assert confirmed["project_revision"] == 4

    gate_tree = [{
        "gate_id": "positive", "name": "Positive", "parent_gate_id": "root",
        "kind": "RANGE_1D", "channels": ["RL1-A"],
        "parameters": {"min": 0, "max": 100},
    }]
    client = TestClient(api.app)
    wrong_kind = client.post(f"{base}/gates", json={
        "name": "Invalid matrix link", "compensation_decision": "EXTERNAL_MATRIX",
        "compensation_artifact_id": upload["artifact_id"], "expected_revision": 4,
        "gate_tree": gate_tree,
    })
    assert wrong_kind.status_code == 422

    source_matrix = {
        "schema_version": "1.0", "detectors": ["RL1-A", "YL2-A"],
        "fluorochromes": ["Alexa Fluor 647-A", "mCherry-A"],
        "matrix": [[1.0, 0.12], [0.08, 1.0]],
    }
    matrix_upload = client.post(
        f"{base}/compensation-matrices",
        files={"matrix": ("spill.json", json.dumps(source_matrix).encode(), "application/json")},
        data={"expected_revision": "4"},
    )
    assert matrix_upload.status_code == 201, matrix_upload.text
    matrix_record = matrix_upload.json()["data"]
    matrix_artifact_id = matrix_record["artifact_id"]
    matrix_artifact = next(
        item for item in store.get(project_id)["artifacts"] if item["id"] == matrix_artifact_id
    )
    assert matrix_artifact["read_only"] == 1
    assert matrix_artifact["detectors"] == ["RL1-A", "YL2-A"]

    duplicate_matrix_upload = client.post(
        f"{base}/compensation-matrices",
        files={"matrix": ("spill.json", json.dumps(source_matrix).encode(), "application/json")},
        data={"expected_revision": str(matrix_record["project_revision"])},
    )
    assert duplicate_matrix_upload.status_code == 201, duplicate_matrix_upload.text
    duplicate_record = duplicate_matrix_upload.json()["data"]
    assert duplicate_record["artifact_id"] == matrix_artifact_id
    assert duplicate_record["already_registered"] is True
    assert duplicate_record["project_revision"] == matrix_record["project_revision"]
    assert sum(
        item["kind"] == "compensation_matrix" for item in store.get(project_id)["artifacts"]
    ) == 1

    gate = client.post(f"{base}/gates", json={
        "name": "External matrix revision", "compensation_decision": "EXTERNAL_MATRIX",
        "compensation_artifact_id": matrix_artifact_id,
        "expected_revision": matrix_record["project_revision"], "gate_tree": gate_tree,
    })
    assert gate.status_code == 201, gate.text
    gate_revision = gate.json()["data"]

    worker_result = {"schema_version": "1.0", "status": "COMPLETE", "records": {
        "samples": [{"file_name": "H1.fcs", "sample_id": "H1.fcs", "event_count": 42,
                     "gate_counts": [{"gate_id": "positive", "gate_name": "Positive",
                                      "gate_path": ["root"], "parent_gate": "root",
                                      "parent_count": 42, "count": 11,
                                      "absolute_fraction": 11 / 42,
                                      "parent_fraction": 11 / 42}]}]
    }, "provenance": {"flowkit_version": "1.3.2", "worker": "flowkit_worker"}}
    worker_requests = []

    def fake_worker(command, **kwargs):
        worker_requests.append(json.loads(kwargs["input"]))
        return SimpleNamespace(returncode=0, stdout=json.dumps(worker_result), stderr="")

    monkeypatch.setattr(api.subprocess, "run", fake_worker)
    analysis = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": gate_revision["gate_revision_id"],
        "fcs_artifact_ids": [upload["artifact_id"]],
        "expected_revision": gate_revision["project_revision"],
    })
    assert analysis.status_code == 201, analysis.text
    assert worker_requests[0]["compensation_matrix"]["matrix"] == source_matrix["matrix"]
    result_artifact = analysis.json()["data"]["result_artifact"]["id"]
    saved_result = client.get(f"{base}/files/{result_artifact}").json()
    assert saved_result["compensation"] == {
        "decision": "EXTERNAL_MATRIX", "artifact_id": matrix_artifact_id,
        "artifact_sha256": matrix_artifact["sha256"], "detectors": ["RL1-A", "YL2-A"],
    }
    assert saved_result["samples"][0]["gate_counts"][0]["count"] == 11
    restarted = ExperimentStore(tmp_path, JobStore(tmp_path))
    assert restarted.get(project_id)["status"] == "COMPLETE"

    gate_result_files = list(
        (tmp_path / "experiments" / project_id / "analyses" / "gate_runs").rglob("gate_counts.json")
    )
    with pytest.raises(RevisionConflict):
        store.save_gate_analysis(
            project_id,
            gate_revision_id=gate_revision["gate_revision_id"],
            input_artifact_ids=[upload["artifact_id"]],
            worker_result={
                "samples": worker_result["records"]["samples"],
                "provenance": worker_result["provenance"],
            },
            expected_revision=gate_revision["project_revision"],
        )
    assert list(
        (tmp_path / "experiments" / project_id / "analyses" / "gate_runs").rglob("gate_counts.json")
    ) == gate_result_files

    effect_result_root = tmp_path / "experiments" / project_id / "analyses" / "effect_runs"
    with pytest.raises(RevisionConflict):
        store.save_effect_analysis(
            project_id,
            gate_analysis_id=analysis.json()["data"]["analysis_id"],
            positive_gate_id="positive",
            document={"synthetic_cleanup_probe": True},
            expected_revision=gate_revision["project_revision"],
        )
    assert list(effect_result_root.rglob("effect_results.json")) == []


def test_flowjo_qc_review_is_bound_to_import_preview_and_sample(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "FlowJo QC review")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id, name="Myc", date="2026-09-21", assay="ANTIGEN_MYC",
        cell_line="unspecified", dose=None, dose_unit=None, positive_gate_id=None,
        independent_batch=True, settings={}, expected_revision=1,
    )
    fcs = tmp_path / "H1.fcs"
    fcs.write_bytes(b"FCS fixture")
    imported = store.import_fcs_files(
        project_id,
        [{"path": str(fcs), "original_filename": "H1.fcs", "metadata": {
            "event_count": 42, "channels": [{"pnn": "FSC-A"}], "metadata": {},
            "compensation_source": "UNKNOWN_COMPENSATION_SOURCE",
        }}],
        batch_id=batch["batch_id"], expected_revision=2,
    )
    upload = imported["files"][0]
    confirmed = store.confirm_fcs_matches(
        project_id, imported["upload_session_id"], [{
            "upload_id": upload["upload_id"], "batch_id": batch["batch_id"], "role": "NC",
            "replicate_type": "SINGLE", "replicate_index": 1,
            "condition": {"sample_label": "H1"},
        }], expected_revision=3,
    )
    preview = {"schema_version": "1.0", "status": "COMPLETE", "records": {
        "samples": [{"sample_id": "H1.fcs", "qc": {"findings": [{
            "code": "DETECTOR_UPPER_SATURATION", "severity": "REVIEW",
        }]}}], "unresolved_nodes": [], "gate_counts": [{
            "sample_id": "H1.fcs", "gate_name": "Positive", "count": 11,
        }],
        "gate_count_verification": {
            "status": "MISMATCH", "compared_gate_count": 1, "matched_gate_count": 0,
            "mismatched_gate_count": 1, "unmatched_gate_count": 0,
        },
        "gate_method_consistency": {
            "status": "CONSISTENT", "samples_checked": 1, "method_variants": 1,
        },
    }}
    client = TestClient(api.app)
    base = f"/api/v4/experiment-projects/{project_id}"
    monkeypatch.setattr(api, "_import_uploaded_wsp", lambda path, paths: preview)
    flowjo_response = client.post(
        f"{base}/samples/import-flowjo",
        data={"fcs_artifact_ids": json.dumps([upload["artifact_id"]]), "expected_revision": "4"},
        files={"workspace": ("workspace.wsp", b"<Workspace />", "application/xml")},
    )
    assert flowjo_response.status_code == 201, flowjo_response.text
    flowjo = flowjo_response.json()["data"]
    assert flowjo["status"] == "COMPLETE"
    preview_response = client.get(f"{base}/files/{flowjo['preview_artifact_id']}")
    assert preview_response.json()["flowjo_input_artifact_ids"] == [upload["artifact_id"]]
    assert preview_response.json()["records"]["gate_count_verification"]["status"] == "MISMATCH"
    assert preview_response.json()["records"]["gate_counts"][0]["count"] == 11
    pass_without_reason = client.post(f"{base}/qc-reviews", json={
        "flowjo_import_id": flowjo["flowjo_import_id"],
        "sample_id": confirmed["sample_ids"][0], "status": "PASS",
        "reviewer": "operator", "expected_revision": 5,
    })
    assert pass_without_reason.status_code == 422
    review = client.post(f"{base}/qc-reviews", json={
        "flowjo_import_id": flowjo["flowjo_import_id"],
        "sample_id": confirmed["sample_ids"][0], "status": "HOLD",
        "reason": "Inspect scatter boundary accumulation", "reviewer": "operator",
        "expected_revision": 5,
    })
    assert review.status_code == 201, review.text
    record = review.json()["data"]
    assert record["source_type"] == "FLOWJO_IMPORT"
    assert record["source_id"] == flowjo["flowjo_import_id"]
    result = client.get(f"{base}/qc-reviews?source_type=FLOWJO_IMPORT&source_id={flowjo['flowjo_import_id']}")
    assert [item["status"] for item in result.json()["data"]] == ["HOLD"]


def test_flowjo_percentage_api_requires_qc_and_records_descriptive_effect(monkeypatch, tmp_path):
    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "FlowJo percentage effect")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id, name="Synthetic activation", date="2026-09-28", assay="ACTIVATION",
        cell_line="Jurkat", dose=None, dose_unit=None, positive_gate_id="positive",
        independent_batch=True, settings={}, expected_revision=1,
    )

    conditions = [
        ("A1", "ACTIVATION", "CTLA4-WT", "REFERENCE", "L32", "ANTIBODY"),
        ("A2", "ACTIVATION", "CTLA4-E1", "VARIANT", "L32", "ANTIBODY"),
        ("A3", "NC", None, "NONE", "L32", "ANTIBODY"),
        ("B1", "OTHER_CONTROL", None, "NONE", "CD80", "POSITIVE_LIGAND_CONTROL"),
    ]
    uploads = []
    matches = []
    for well_id, role, antigen_label, antigen_role, panel_label, panel_role in conditions:
        fcs = tmp_path / f"{well_id}.fcs"
        fcs.write_bytes(f"synthetic-{well_id}".encode("ascii"))
        uploads.append({
            "path": str(fcs),
            "original_filename": fcs.name,
            "metadata": {
                "event_count": 100,
                "channels": [],
                "metadata": {},
                "compensation_source": "USER_CONFIRMED_NOT_APPLIED",
            },
        })
        matches.append({
            "role": role,
            "replicate_type": "TECHNICAL",
            "replicate_index": 1,
            "temporary_identity": {
                "status": "TEMPORARY_LABEL",
                "source": "USER_CONFIRMED_LABELS",
                "antigen_label": antigen_label,
                "antigen_role": antigen_role,
                "panel_label": panel_label,
                "panel_role": panel_role,
            },
            "condition": {
                "well_id": well_id,
                "sample_label": f"{panel_label} {antigen_label or role} {well_id}",
                "condition_label": antigen_label or role,
                "panel_label": panel_label,
            },
        })

    imported_fcs = store.import_fcs_files(
        project_id, uploads, batch_id=batch["batch_id"], expected_revision=2
    )
    for match, imported in zip(matches, imported_fcs["files"], strict=True):
        match["upload_id"] = imported["upload_id"]
        match["batch_id"] = batch["batch_id"]
    confirmed = store.confirm_fcs_matches(
        project_id, imported_fcs["upload_session_id"], matches, expected_revision=3
    )

    # A user-supplied, already-gated percentage table must work without a
    # FlowJo workspace. This API path does not call the FlowKit worker.
    client = TestClient(api.app)
    base = f"/api/v4/experiment-projects/{project_id}"
    partial_table = (
        b"Well,Target-positive | Freq. of Parent\n"
        b"A1,70.0%\nA2,35.0%\nA3,5.0%\n"
    )
    partial_import = client.post(
        f"{base}/flowjo-percentage-imports",
        data={
            "batch_id": batch["batch_id"],
            "numeric_scale": "AUTO",
            "expected_revision": "4",
        },
        files={"table": ("partial-flowjo.csv", partial_table, "text/csv")},
    )
    assert partial_import.status_code == 201, partial_import.text
    partial_data = partial_import.json()["data"]
    assert len(partial_data["result"]["missing_sample_ids"]) == 1
    assert partial_data["result"]["compensation"]["status"] == "UNKNOWN"
    partial_effect = client.post(
        f"{base}/flowjo-percentage-imports/{partial_data['import_id']}/effect-analyses",
        json={"expected_revision": store.get(project_id)["revision"]},
    )
    assert partial_effect.status_code == 422
    assert partial_effect.json()["detail"]["code"] == "QC_HOLD"
    assert "omits 1 registered samples" in partial_effect.json()["detail"]["message"]

    invalid_separator = client.post(
        f"{base}/flowjo-percentage-imports",
        data={
            "batch_id": batch["batch_id"],
            "decimal_separator": "AUTO",
            "expected_revision": str(store.get(project_id)["revision"]),
        },
        files={"table": ("invalid-separator.csv", partial_table, "text/csv")},
    )
    assert invalid_separator.status_code == 422
    assert invalid_separator.json()["detail"]["code"] == "FLOWJO_DECIMAL_SEPARATOR_INVALID"
    assert store.get(project_id)["revision"] == 5

    table = (
        b"Well;Target-positive | Freq. of Parent\n"
        b"A1;70,0%\nA2;35,0%\nA3;5,0%\nB1;90,0%\n"
    )
    imported_table = client.post(
        f"{base}/flowjo-percentage-imports",
        data={
            "batch_id": batch["batch_id"],
            "numeric_scale": "AUTO",
            "decimal_separator": ",",
            "compensation_status": "APPLIED",
            "compensation_note": "Synthetic declared setup",
            "expected_revision": str(store.get(project_id)["revision"]),
        },
        files={"table": ("flowjo.csv", table, "text/csv")},
    )
    assert imported_table.status_code == 201, imported_table.text
    table_data = imported_table.json()["data"]
    percentage_import_id = table_data["import_id"]
    assert table_data["flowjo_import_id"] is None
    assert table_data["result"]["missing_sample_ids"] == []
    assert table_data["result"]["parser_provenance"]["numeric_scale"] == "AUTO"
    assert table_data["result"]["parser_provenance"]["decimal_separator"] == ","
    assert table_data["result"]["observations"][0]["fraction"] == pytest.approx(0.70)
    assert table_data["result"]["observations"][0]["displayed_value"] == "70,0%"
    assert table_data["result"]["compensation"] == {
        "status": "APPLIED",
        "source": "USER_REPORTED_RESULT_TABLE_CONTEXT",
        "note": "Synthetic declared setup",
    }
    assert all(item["positive_count"] is None and item["parent_count"] is None
               for item in table_data["result"]["observations"])

    blocked = client.post(
        f"{base}/flowjo-percentage-imports/{percentage_import_id}/effect-analyses",
        json={"expected_revision": store.get(project_id)["revision"]},
    )
    assert blocked.status_code == 422
    assert blocked.json()["detail"]["code"] == "QC_HOLD"
    assert store.get(project_id)["revision"] == 6

    outside_import = client.post(f"{base}/qc-reviews", json={
        "flowjo_percentage_import_id": percentage_import_id,
        "sample_id": "not-in-the-table", "status": "PASS",
        "reviewer": "synthetic-test-reviewer", "expected_revision": store.get(project_id)["revision"],
    })
    assert outside_import.status_code == 404

    import io
    from openpyxl import load_workbook

    worksheet_response = client.get(
        f"{base}/flowjo-percentage-imports/{percentage_import_id}/qc-worksheet"
    )
    assert worksheet_response.status_code == 200
    assert worksheet_response.headers["content-type"].startswith(
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )
    workbook = load_workbook(io.BytesIO(worksheet_response.content), data_only=False)
    worksheet = workbook["QC Review"]
    headers = {cell.value: cell.column for cell in worksheet[1]}
    assert worksheet.max_row == 5
    assert worksheet.cell(2, headers["displayed_value"]).value == "70,0%"
    assert worksheet.cell(2, headers["current_qc_review_id"]).value == ""

    invalid_worksheet = io.BytesIO()
    worksheet.cell(2, headers["decision"], "HOLD")
    worksheet.cell(2, headers["reviewer"], "synthetic-test-reviewer")
    # Missing HOLD rationale invalidates the whole worksheet before any write.
    for row_number in range(3, worksheet.max_row + 1):
        worksheet.cell(row_number, headers["decision"], "PASS")
        worksheet.cell(row_number, headers["reviewer"], "synthetic-test-reviewer")
    workbook.save(invalid_worksheet)
    invalid_result = client.post(
        f"{base}/flowjo-percentage-imports/{percentage_import_id}/qc-worksheet",
        files={"worksheet": ("review.xlsx", invalid_worksheet.getvalue(),
                             "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")},
    )
    assert invalid_result.status_code == 422
    assert invalid_result.json()["detail"]["code"] == "QC_WORKSHEET_INVALID"
    assert store.get(project_id)["revision"] == 6
    assert store.list_qc_reviews(project_id, "FLOWJO_PERCENTAGE_IMPORT", percentage_import_id) == []

    worksheet.cell(2, headers["decision"], "PASS")
    worksheet.cell(5, headers["decision"], "PASS_WITH_OVERRIDE")
    worksheet.cell(5, headers["reason"], "Synthetic worksheet override")
    valid_worksheet = io.BytesIO()
    workbook.save(valid_worksheet)
    valid_bytes = valid_worksheet.getvalue()
    imported_reviews = client.post(
        f"{base}/flowjo-percentage-imports/{percentage_import_id}/qc-worksheet",
        files={"worksheet": ("review.xlsx", valid_bytes,
                             "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")},
    )
    assert imported_reviews.status_code == 201, imported_reviews.text
    imported_summary = imported_reviews.json()["data"]
    assert imported_summary["applied_count"] == 4
    assert imported_summary["unchanged_count"] == 0
    assert imported_summary["project_revision"] == 10
    assert imported_summary["already_applied"] is False

    retried_reviews = client.post(
        f"{base}/flowjo-percentage-imports/{percentage_import_id}/qc-worksheet",
        files={"worksheet": ("review.xlsx", valid_bytes,
                             "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")},
    )
    assert retried_reviews.status_code == 201
    assert retried_reviews.json()["data"]["already_applied"] is True
    assert store.get(project_id)["revision"] == 10
    review_ledger = store.list_qc_reviews(project_id, "FLOWJO_PERCENTAGE_IMPORT", percentage_import_id)
    assert len(review_ledger) == 4
    assert {item["status"] for item in review_ledger} == {"PASS", "PASS_WITH_OVERRIDE"}

    refreshed_worksheet = client.get(
        f"{base}/flowjo-percentage-imports/{percentage_import_id}/qc-worksheet"
    )
    unchanged_upload = client.post(
        f"{base}/flowjo-percentage-imports/{percentage_import_id}/qc-worksheet",
        files={"worksheet": ("review.xlsx", refreshed_worksheet.content,
                             "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")},
    )
    assert unchanged_upload.status_code == 201
    assert unchanged_upload.json()["data"]["applied_count"] == 0
    assert unchanged_upload.json()["data"]["project_revision"] == 10

    effect_response = client.post(
        f"{base}/flowjo-percentage-imports/{percentage_import_id}/effect-analyses",
        json={"expected_revision": store.get(project_id)["revision"]},
    )
    assert effect_response.status_code == 201, effect_response.text
    result = effect_response.json()["data"]
    assert len(result["effects"]) == 1
    effect = result["effects"][0]
    assert effect["panel_label"] == "L32"
    assert effect["replicate_type"] == "TECHNICAL"
    assert effect["batch_effect"]["difference_pp"] == pytest.approx(-35)
    assert effect["batch_effect"]["decline_pct"] == pytest.approx(100 * (1 - 0.30 / 0.65))
    assert effect["replicate_count_by_group"] == {"mutant": 1, "reference": 1, "nc": 1}
    assert {row["well_id"]: row["comparison_group"] for row in effect["replicate_observations"]} == {
        "A1": "REFERENCE", "A2": "MUTANT", "A3": "MUTANT_NC",
    }
    assert all(row["accepted_for_summary"] for row in effect["replicate_observations"])
    stored = client.get(
        f"{base}/flowjo-percentage-imports/{percentage_import_id}/effect-analyses"
    )
    assert stored.status_code == 200
    provenance = stored.json()["data"][-1]["result"]["provenance"]
    assert provenance["measurement_basis"] == "PERCENTAGE_ONLY"
    assert provenance["denominator_available"] is False
    assert provenance["compensation"]["status"] == "APPLIED"
    assert any("无法独立验证所用补偿矩阵" in item
               for item in stored.json()["data"][-1]["result"]["limitations"])
    assert provenance["qc_source_type"] == "FLOWJO_PERCENTAGE_IMPORT"
    assert provenance["qc_source_id"] == percentage_import_id
    assert provenance["excluded_control_panels"][0]["panel_label"] == "CD80"
    assert provenance["excluded_control_panels"][0]["panel_role"] == "POSITIVE_LIGAND_CONTROL"

    import multiprocessing

    effect_document = dict(stored.json()["data"][-1]["result"])
    for key in (
        "schema_version", "analysis_id", "project_id",
        "flowjo_percentage_import_id", "status", "created_at",
    ):
        effect_document.pop(key, None)
    effect_root = (
        tmp_path / "experiments" / project_id / "analyses" / "flowjo_percentage_effects"
    )
    for crash_point, exit_code in (("artifact_published", 95), ("manifest_published", 96)):
        current = store.get(project_id)
        registered_artifacts = current["manifest"]["artifacts"]
        existing_effect_count = len(
            store.list_flowjo_percentage_effect_analyses(project_id, percentage_import_id)
        )
        existing_effect_paths = list(effect_root.rglob("effect_results.json"))
        process = multiprocessing.get_context("spawn").Process(
            target=_crash_flowjo_percentage_effect_save,
            args=(
                str(tmp_path), project_id, percentage_import_id, effect_document,
                current["revision"], crash_point,
            ),
        )
        process.start()
        process.join(timeout=30)
        assert not process.is_alive()
        assert process.exitcode == exit_code

        recovered = ExperimentStore(tmp_path, JobStore(tmp_path))
        recovered_project = recovered.get(project_id)
        assert recovered_project["status"] == "COMPLETE"
        assert recovered_project["revision"] == current["revision"]
        assert recovered_project["manifest"]["artifacts"] == registered_artifacts
        assert len(
            recovered.list_flowjo_percentage_effect_analyses(project_id, percentage_import_id)
        ) == existing_effect_count
        assert len(list(effect_root.rglob("effect_results.json"))) == len(existing_effect_paths) + 1

        recovered.save_flowjo_percentage_effect_analysis(
            project_id,
            import_id=percentage_import_id,
            document=effect_document,
            expected_revision=recovered_project["revision"],
        )
        assert len(
            recovered.list_flowjo_percentage_effect_analyses(project_id, percentage_import_id)
        ) == existing_effect_count + 1
        assert len(list(effect_root.rglob("effect_results.json"))) == len(existing_effect_paths) + 2
        store = recovered


def test_qc_review_schema_migration_preserves_gate_analysis_reviews(tmp_path):
    import sqlite3

    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "QC migration")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id, name="Fixture", date="2026-09-21", assay="OTHER",
        cell_line="unspecified", dose=None, dose_unit=None, positive_gate_id=None,
        independent_batch=True, settings={}, expected_revision=1,
    )
    fcs = tmp_path / "sample.fcs"
    fcs.write_bytes(b"migration sample")
    imported = store.import_fcs_files(
        project_id,
        [{"path": str(fcs), "original_filename": "sample.fcs", "metadata": {
            "event_count": 4, "channels": [], "metadata": {},
            "compensation_source": "UNKNOWN_COMPENSATION_SOURCE",
        }}],
        batch_id=batch["batch_id"], expected_revision=2,
    )
    upload = imported["files"][0]
    confirmed = store.confirm_fcs_matches(
        project_id, imported["upload_session_id"], [{
            "upload_id": upload["upload_id"], "batch_id": batch["batch_id"], "role": "OTHER_CONTROL",
            "replicate_type": "SINGLE", "replicate_index": 1,
            "condition": {"sample_label": "sample"},
        }], expected_revision=3,
    )
    sample_id = confirmed["sample_ids"][0]
    with sqlite3.connect(store.db_path) as db:
        db.execute("DROP TABLE qc_reviews")
        db.executescript("""
            CREATE TABLE qc_reviews (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                gate_analysis_id TEXT NOT NULL REFERENCES gate_analyses(id),
                sample_id TEXT NOT NULL REFERENCES samples(id), version INTEGER NOT NULL,
                status TEXT NOT NULL, reason TEXT NOT NULL, reviewer TEXT NOT NULL,
                project_revision INTEGER NOT NULL, created_at TEXT NOT NULL,
                UNIQUE(gate_analysis_id,sample_id,version)
            );
        """)
        db.execute(
            "INSERT INTO qc_reviews VALUES (?,?,?,?,?,?,?,?,?,?)",
            ("legacy-review", project_id, "legacy-analysis", sample_id, 1, "PASS_WITH_OVERRIDE",
             "Legacy reason", "operator", 4, "2026-09-21T00:00:00+00:00"),
        )
        db.execute("DELETE FROM schema_migrations WHERE version=8")
        db.execute("DELETE FROM schema_migrations WHERE version=14")
        db.commit()

    migrated = ExperimentStore(tmp_path, JobStore(tmp_path))
    reviews = migrated.list_qc_reviews(project_id)
    assert len(reviews) == 1
    assert reviews[0]["source_type"] == "GATE_ANALYSIS"
    assert reviews[0]["source_id"] == "legacy-analysis"
    assert reviews[0]["reason"] == "Legacy reason"


def test_activation_effect_to_grasp_export_restarts_with_frozen_evidence(
    monkeypatch, tmp_path
):
    import hashlib
    import io
    import zipfile
    from types import SimpleNamespace

    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Activation effects")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id, name="Activation 1", date="2026-09-21", assay="ACTIVATION",
        cell_line="Jurkat", dose=1.0, dose_unit="ug/mL", positive_gate_id="positive",
        independent_batch=True, settings={}, expected_revision=1,
    )
    incomplete_batch = store.list_entities(project_id, "assay_batches")[0]
    assert incomplete_batch["design_status"] == "INCOMPLETE"
    assert set(incomplete_batch["design_warnings"]) == {
        "ACTIVATION_SAMPLE_MISSING", "NC_CONTROL_MISSING", "WT_CONSTRUCT_MISSING"
    }
    antibody = store.add_antibody_construct(
        project_id, name="Parent Fab", role="PARENT", vh_sequence="ACD", vl_sequence="EFG",
        parent_construct_id=None, expected_revision=2,
    )
    reference = store.add_antigen_construct(
        project_id, name="WT reference", mutations=None, full_sequence="ACD", expected_revision=3,
    )
    mutant = store.add_antigen_construct(
        project_id, name="H1 temporary", mutations=[{"position": 11, "wt": "C", "mut": "A"}],
        full_sequence=None, expected_revision=4,
    )
    positive_counts = [50, 80, 10, 10]
    files = []
    for label in ("mutant", "reference", "mutant-nc", "reference-nc"):
        path = tmp_path / f"{label}.fcs"
        path.write_bytes(f"synthetic FCS {label}".encode())
        files.append({"path": str(path), "original_filename": f"{label}.fcs", "metadata": {
            "event_count": 100, "channels": [{"pnn": "FSC-A"}, {"pnn": "SSC-A"}],
            "metadata": {}, "compensation_source": "UNKNOWN_COMPENSATION_SOURCE",
        }})
    imported = store.import_fcs_files(project_id, files, batch_id=batch["batch_id"], expected_revision=5)
    with pytest.raises(ValueError, match="require antigen and antibody"):
        store.confirm_fcs_matches(
            project_id,
            imported["upload_session_id"],
            [{
                "upload_id": imported["files"][0]["upload_id"],
                "batch_id": batch["batch_id"],
                "role": "ACTIVATION",
                "replicate_type": "BIOLOGICAL",
                "replicate_index": 1,
                "condition": {"sample_label": "incomplete"},
            }],
            expected_revision=6,
        )
    roles = ["ACTIVATION", "ACTIVATION", "NC", "NC"]
    antigens = [mutant["construct_id"], reference["construct_id"], mutant["construct_id"], reference["construct_id"]]
    matches = [{
        "upload_id": upload["upload_id"], "batch_id": batch["batch_id"], "role": role,
        "antigen_construct_id": antigen, "antibody_construct_id": antibody["construct_id"],
        "replicate_type": "BIOLOGICAL", "replicate_index": index + 1,
        "condition": {"sample_label": upload["original_filename"]},
    } for index, (upload, role, antigen) in enumerate(zip(imported["files"], roles, antigens, strict=True))]
    confirmed = store.confirm_fcs_matches(
        project_id, imported["upload_session_id"], matches, expected_revision=6,
    )
    ready_batch = store.list_entities(project_id, "assay_batches")[0]
    assert ready_batch["design_status"] == "READY"
    assert ready_batch["design_warnings"] == []
    client = TestClient(api.app)
    base = f"/api/v4/experiment-projects/{project_id}"
    gate = client.post(f"{base}/gates", json={
        "name": "Activation positive", "compensation_decision": "NONE", "expected_revision": 7,
        "gate_tree": [{
            "gate_id": "positive", "name": "Positive", "parent_gate_id": "root",
            "kind": "RANGE_1D", "channels": ["FSC-A"], "parameters": {"min": 0, "max": 100},
        }],
    })
    assert gate.status_code == 201, gate.text
    gate_revision_id = gate.json()["data"]["gate_revision_id"]
    worker_result = {"status": "COMPLETE", "records": {"samples": [
        {"file_name": upload["original_filename"], "sample_id": upload["original_filename"],
         "event_count": 100, "gate_counts": [{
             "gate_id": "positive", "gate_name": "Positive", "gate_path": ["root"],
             "parent_gate": "root", "parent_count": 100, "count": positive,
             "absolute_fraction": positive / 100, "parent_fraction": positive / 100,
         }], "qc": {"events_examined": 100, "events_removed": 0, "findings": []}}
        for upload, positive in zip(imported["files"], positive_counts, strict=True)
    ]}, "provenance": {"flowkit_version": "1.3.2"}}
    monkeypatch.setattr(api.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(
        returncode=0, stdout=json.dumps(worker_result), stderr=""
    ))
    gated = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": gate_revision_id,
        "fcs_artifact_ids": [upload["artifact_id"] for upload in imported["files"]],
        "expected_revision": 8,
    })
    assert gated.status_code == 201, gated.text
    gate_analysis_id = gated.json()["data"]["analysis_id"]
    revision = gated.json()["data"]["project_revision"]
    for sample_id in confirmed["sample_ids"]:
        reviewed = client.post(f"{base}/qc-reviews", json={
            "gate_analysis_id": gate_analysis_id, "sample_id": sample_id,
            "status": "PASS", "reviewer": "reviewer", "expected_revision": revision,
        })
        assert reviewed.status_code == 201, reviewed.text
        revision = reviewed.json()["data"]["project_revision"]

    effect = client.post(f"{base}/effect-analyses", json={
        "gate_analysis_id": gate_analysis_id, "positive_gate_id": "positive",
        "samples": [
            {"sample_id": confirmed["sample_ids"][0], "probe_id": mutant["construct_id"], "condition": "MUTANT"},
            {"sample_id": confirmed["sample_ids"][1], "probe_id": mutant["construct_id"], "condition": "REFERENCE"},
            {"sample_id": confirmed["sample_ids"][2], "probe_id": mutant["construct_id"], "condition": "MUTANT_NC"},
            {"sample_id": confirmed["sample_ids"][3], "probe_id": mutant["construct_id"], "condition": "REFERENCE_NC"},
        ],
        "expected_revision": revision,
    })
    assert effect.status_code == 201, effect.text
    data = effect.json()["data"]
    assert data["batch_effects"][0]["batch_effect"]["difference_pp"] == pytest.approx(-30)
    assert data["batch_effects"][0]["batch_effect"]["retention"] == pytest.approx(4 / 7)
    assert data["batch_effects"][0]["batch_effect"]["decline_pct"] == pytest.approx(300 / 7)
    assert data["surface_projection"]["probes"][0]["attribution"] == "SINGLE_SITE"
    assert data["surface_projection"]["probes"][0]["valid_decline_batch_count"] == 1
    saved = client.get(f"{base}/files/{data['result_artifact']['id']}")
    assert saved.status_code == 200
    assert saved.json()["well_results"][0]["parent_count"] == 100
    assert saved.json()["frozen_hashes"] == store.get(project_id)["frozen_hashes"]

    restarted_after_effect = ExperimentStore(tmp_path, JobStore(tmp_path))
    monkeypatch.setattr(api, "EXPERIMENTS", restarted_after_effect)
    assert restarted_after_effect.get(project_id)["status"] == "COMPLETE"
    export_payload = {
        "antibody_construct_id": antibody["construct_id"],
        "experimental_ir": [{
            "constraint": "1-2-C",
            "evidence_id": data["analysis_id"],
            "rationale": "Synthetic single-site activation decline",
            "details": {"probe_id": mutant["construct_id"]},
        }],
        "predicted_ir": [],
        "manual_rpr": [],
        "merge_sources": ["experimental_ir"],
        "reviewer": "integration reviewer",
        "expected_revision": data["project_revision"],
    }
    exported = client.post(f"{base}/constraint-exports", json=export_payload)
    assert exported.status_code == 201, exported.text
    export = exported.json()["data"]
    assert export["manifest"]["status"] == "IR_ONLY"
    assert export["manifest"]["constraint_count"] == 1
    export_dir = tmp_path / "experiments" / project_id / "analyses" / "constraint_exports"
    existing_exports = list(export_dir.glob("*.zip"))
    stale_export = client.post(f"{base}/constraint-exports", json=export_payload)
    assert stale_export.status_code == 409, stale_export.text
    assert list(export_dir.glob("*.zip")) == existing_exports
    first_download = client.get(export["download_url"])
    assert first_download.status_code == 200
    assert hashlib.sha256(first_download.content).hexdigest() == export["artifact"]["sha256"]
    with zipfile.ZipFile(io.BytesIO(first_download.content)) as bundle:
        assert bundle.read("grasp_constraints.txt") == b"1-2-C\n"
        evidence = json.loads(bundle.read("evidence_records.json"))
        record = evidence["experimental_ir"][0]
        assert record["evidence_id"] == data["analysis_id"]
        assert record["details"]["probe_id"] == mutant["construct_id"]
        assert record["details"]["uniprot_position"] == 11
        assert record["details"]["ecd_position"] == 2
        assert record["details"]["effect_artifact_sha256"] == data["result_artifact"]["sha256"]

    restarted_after_export = ExperimentStore(tmp_path, JobStore(tmp_path))
    monkeypatch.setattr(api, "EXPERIMENTS", restarted_after_export)
    assert restarted_after_export.get(project_id)["status"] == "COMPLETE"
    listed = client.get(f"{base}/constraint-exports")
    assert listed.status_code == 200
    assert listed.json()["data"][0]["sha256"] == export["artifact"]["sha256"]
    second_download = client.get(export["download_url"])
    assert second_download.status_code == 200
    assert second_download.content == first_download.content


def test_manual_time_window_exclusion_creates_audited_filtered_gate_analysis(
    monkeypatch, tmp_path
):
    from types import SimpleNamespace

    import ecd_epitope_scanner.api as api
    from ecd_epitope_scanner.experiment_projects import ExperimentStore

    jobs = JobStore(tmp_path)
    store = ExperimentStore(tmp_path, jobs)
    monkeypatch.setattr(api, "EXPERIMENTS", store)
    job = _completed_v4_job(jobs)
    project = store.create(job["id"], "Manual event-window exclusion")
    project_id = project["project_id"]
    batch = store.add_batch(
        project_id, name="Window fixture", date="2026-09-21", assay="ANTIGEN_MYC",
        cell_line="unspecified", dose=None, dose_unit=None, positive_gate_id=None,
        independent_batch=True, settings={}, expected_revision=project["revision"],
    )
    fcs = tmp_path / "H1.fcs"
    original_bytes = b"registered test FCS; exclusion must not mutate it"
    fcs.write_bytes(original_bytes)
    fcs2 = tmp_path / "H2.fcs"
    original_bytes_2 = b"second registered test FCS; exclusion must not mutate it"
    fcs2.write_bytes(original_bytes_2)
    imported = store.import_fcs_files(
        project_id,
        [{"path": str(path), "original_filename": path.name, "metadata": {
            "event_count": 42, "channels": [{"pnn": "FSC-A"}, {"pnn": "SSC-A"}],
            "metadata": {}, "compensation_source": "UNKNOWN_COMPENSATION_SOURCE",
        }} for path in (fcs, fcs2)],
        batch_id=batch["batch_id"], expected_revision=store.get(project_id)["revision"],
    )
    confirmed = store.confirm_fcs_matches(
        project_id, imported["upload_session_id"], [{
            "upload_id": upload["upload_id"], "batch_id": batch["batch_id"], "role": "NC",
            "replicate_type": "SINGLE", "replicate_index": 1,
            "condition": {"sample_label": upload["original_filename"]},
        } for upload in imported["files"]], expected_revision=store.get(project_id)["revision"],
    )
    client = TestClient(api.app)
    base = f"/api/v4/experiment-projects/{project_id}"
    gate = client.post(f"{base}/gates", json={
        "name": "All events", "compensation_decision": "NONE",
        "expected_revision": store.get(project_id)["revision"],
        "gate_tree": [{
            "gate_id": "all", "name": "All", "parent_gate_id": "root",
            "kind": "RANGE_1D", "channels": ["FSC-A"],
            "parameters": {"min": -1e12, "max": 1e12},
        }],
    })
    assert gate.status_code == 201, gate.text
    gate_revision_id = gate.json()["data"]["gate_revision_id"]
    ranges = [{"event_start": 0, "event_end": 20}, {"event_start": 20, "event_end": 42}]

    def worker_sample(file_name, event_count=42, removed=0, selected_ranges=None, qc_windows=None):
        return {
            "file_name": file_name, "sample_id": file_name, "event_count": event_count,
            "source_event_count": 42, "events_removed": removed,
            "event_exclusion_ranges": selected_ranges or [],
            "gate_counts": [{
                "gate_id": "all", "gate_name": "All", "gate_path": ["root"],
                "parent_gate": "root", "parent_count": event_count, "count": event_count,
                "absolute_fraction": 1.0, "parent_fraction": 1.0,
            }],
            "qc": {"events_examined": event_count, "events_removed": removed,
                   "findings": [], "time_windows": qc_windows or ranges},
        }

    initial_response = {"status": "COMPLETE", "records": {"samples": [worker_sample("H1.fcs"), worker_sample("H2.fcs")]},
                        "provenance": {"flowkit_version": "1.3.2"}}
    excluded_response = {"status": "COMPLETE", "records": {
        "samples": [worker_sample("H1.fcs", 22, 20, [ranges[0]], [
            {"event_start": 0, "event_end": 10}, {"event_start": 10, "event_end": 22},
        ]), worker_sample("H2.fcs")]},
        "provenance": {"flowkit_version": "1.3.2"}}
    excluded_response_2 = {"status": "COMPLETE", "records": {
        "samples": [worker_sample("H1.fcs", 22, 20, [ranges[0]], [
            {"event_start": 0, "event_end": 10}, {"event_start": 10, "event_end": 22},
        ]), worker_sample("H2.fcs", 22, 20, [ranges[0]], [
            {"event_start": 0, "event_end": 10}, {"event_start": 10, "event_end": 22},
        ])]},
        "provenance": {"flowkit_version": "1.3.2"}}
    excluded_response_3 = {"status": "COMPLETE", "records": {
        "samples": [worker_sample("H1.fcs", 22, 20, [ranges[0]], [
            {"event_start": 0, "event_end": 10}, {"event_start": 10, "event_end": 22},
        ]), worker_sample("H2.fcs", 12, 30, [ranges[0], {"event_start": 20, "event_end": 30}], [
            {"event_start": 0, "event_end": 6}, {"event_start": 6, "event_end": 12},
        ])]},
        "provenance": {"flowkit_version": "1.3.2"}}
    worker_requests = []

    def run_worker(command, **kwargs):
        worker_request = json.loads(kwargs["input"])
        worker_requests.append(worker_request)
        ranges_by_path = worker_request.get("event_exclusion_ranges_by_path", [[], []])
        if not any(ranges_by_path):
            response = initial_response
        elif ranges_by_path == [[ranges[0]], []]:
            response = excluded_response
        elif ranges_by_path == [[ranges[0]], [ranges[0]]]:
            response = excluded_response_2
        elif ranges_by_path == [[ranges[0]], [ranges[0], {"event_start": 20, "event_end": 30}]]:
            response = excluded_response_3
        else:
            raise AssertionError(f"Unexpected exclusion ranges: {ranges_by_path}")
        return SimpleNamespace(returncode=0, stdout=json.dumps(response), stderr="")

    monkeypatch.setattr(api.subprocess, "run", run_worker)
    artifact_ids = [upload["artifact_id"] for upload in imported["files"]]
    initial = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": gate_revision_id, "fcs_artifact_ids": artifact_ids,
        "expected_revision": store.get(project_id)["revision"],
    })
    assert initial.status_code == 201, initial.text
    source = initial.json()["data"]
    sample_id, sample_id_2 = confirmed["sample_ids"]

    invalid = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": gate_revision_id, "fcs_artifact_ids": artifact_ids,
        "source_gate_analysis_id": source["analysis_id"],
        "window_exclusions": [{"sample_id": sample_id, "windows": [
            {"event_start": 0, "event_end": 21},
        ]}], "exclusion_reviewer": "reviewer", "exclusion_reason": "manual review",
        "expected_revision": source["project_revision"],
    })
    assert invalid.status_code == 422
    assert "complete time window" in invalid.json()["detail"]["message"]
    assert len(worker_requests) == 1

    excluded = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": gate_revision_id, "fcs_artifact_ids": artifact_ids,
        "source_gate_analysis_id": source["analysis_id"],
        "window_exclusions": [{"sample_id": sample_id, "windows": [ranges[0]]}],
        "exclusion_reviewer": "reviewer", "exclusion_reason": "unstable acquisition-rate segment",
        "expected_revision": source["project_revision"],
    })
    assert excluded.status_code == 201, excluded.text
    assert worker_requests[1]["event_exclusion_ranges_by_path"] == [[ranges[0]], []]
    result = excluded.json()["data"]
    assert result["samples"][0]["event_count"] == 22
    document = client.get(f"{base}/files/{result['result_artifact']['id']}").json()
    provenance = document["event_exclusions"]
    assert provenance["source_gate_analysis_id"] == source["analysis_id"]
    assert provenance["source_result_sha256"] == source["result_artifact"]["sha256"]
    assert provenance["reviewer"] == "reviewer"
    assert provenance["reason"] == "unstable acquisition-rate segment"
    assert provenance["windows"][0]["events_removed"] == 20
    assert provenance["windows"][0]["events_after_exclusion"] == 22
    assert document["samples"][0]["gate_counts"][0]["parent_count"] == 22
    assert document["samples"][0]["qc"]["events_removed"] == 20

    excluded_again = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": gate_revision_id, "fcs_artifact_ids": artifact_ids,
        "source_gate_analysis_id": result["analysis_id"],
        "window_exclusions": [{"sample_id": sample_id_2, "windows": [ranges[0]]}],
        "exclusion_reviewer": "reviewer 2", "exclusion_reason": "second reviewed segment",
        "expected_revision": result["project_revision"],
    })
    assert excluded_again.status_code == 201, excluded_again.text
    assert worker_requests[2]["event_exclusion_ranges_by_path"] == [[ranges[0]], [ranges[0]]]
    second_result = excluded_again.json()["data"]
    assert [item["event_count"] for item in second_result["samples"]] == [22, 22]
    second_document = client.get(f"{base}/files/{second_result['result_artifact']['id']}").json()
    second_exclusions = {item["sample_id"]: item for item in second_document["event_exclusions"]["windows"]}
    assert second_exclusions[sample_id]["events_removed"] == 20
    assert second_exclusions[sample_id_2]["events_removed"] == 20

    excluded_third = client.post(f"{base}/gate-analyses", json={
        "gate_revision_id": gate_revision_id, "fcs_artifact_ids": artifact_ids,
        "source_gate_analysis_id": second_result["analysis_id"],
        "window_exclusions": [{"sample_id": sample_id_2, "windows": [{"event_start": 0, "event_end": 10}]}],
        "exclusion_reviewer": "reviewer 3", "exclusion_reason": "third reviewed segment",
        "expected_revision": second_result["project_revision"],
    })
    assert excluded_third.status_code == 201, excluded_third.text
    assert worker_requests[3]["event_exclusion_ranges_by_path"] == [
        [ranges[0]], [ranges[0], {"event_start": 20, "event_end": 30}],
    ]
    third_result = excluded_third.json()["data"]
    assert [item["event_count"] for item in third_result["samples"]] == [22, 12]
    assert fcs.read_bytes() == original_bytes
    assert fcs2.read_bytes() == original_bytes_2
