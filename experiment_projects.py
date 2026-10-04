"""Immutable v4 experiment project snapshots and their artifact registry."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import re
import shutil
import sqlite3
import threading
import uuid
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator
from jsonschema import ValidationError as JsonSchemaValidationError

from .experimental_constructs import (
    ConstructError,
    reference_positions,
    validate_antibody_construct,
    validate_antigen_construct,
)
from .experimental_gates import GateError, validate_gate_tree
from .jobs import JobStore, _process_alive
from .util import ClosingSQLiteConnection, read_json, sha256_bytes, sha256_file, write_json


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _discard_unregistered_artifact(path: Path) -> None:
    """Remove a just-published artifact without letting cleanup mask the failure."""
    try:
        if path.is_symlink():
            path.unlink(missing_ok=True)
        elif path.is_file():
            path.chmod(0o644)
            path.unlink(missing_ok=True)
    except OSError:
        pass


def _flowjo_well_id(value: Any) -> str | None:
    """Normalize a FlowJo sample key or registered filename to a plate well ID."""
    text = str(value or "").replace("\\", "/").rsplit("/", 1)[-1].strip().upper()
    if text.endswith(".FCS"):
        text = text[:-4]
    match = re.search(r"(?:^|[^A-Z0-9])([A-P]0?(?:[1-9]|1[0-9]|2[0-4]))$", text)
    if match is None:
        match = re.fullmatch(r"([A-P]0?(?:[1-9]|1[0-9]|2[0-4]))", text)
    if match is None:
        return None
    row, column = re.fullmatch(r"([A-P])0?(\d+)", match.group(1)).groups()
    return f"{row}{int(column)}"


class SnapshotError(ValueError):
    code = "PROJECT_SNAPSHOT_INCOMPLETE"


class IntegrityError(ValueError):
    code = "INTEGRITY_ERROR"


class RevisionConflict(ValueError):
    code = "REVISION_CONFLICT"


class ExperimentStore:
    def __init__(self, root: Path, jobs: JobStore):
        self.root = root.resolve()
        self.jobs = jobs
        self.projects_root = self.root / "experiments"
        self.projects_root.mkdir(parents=True, exist_ok=True)
        self.db_path = self.root / "experiments.sqlite3"
        self._lock = threading.RLock()
        self._initialize()
        with self._connect() as db:
            self._recover_published_project_snapshots(db)
            self._reconcile_registered_project_manifests(db)

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(
            self.db_path, timeout=30, factory=ClosingSQLiteConnection
        )
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys=ON")
        return connection

    def _initialize(self) -> None:
        with self._connect() as db:
            db.execute(
                "CREATE TABLE IF NOT EXISTS schema_migrations (version INTEGER PRIMARY KEY, applied_at TEXT NOT NULL)"
            )
            if db.execute("SELECT 1 FROM schema_migrations WHERE version=1").fetchone():
                self._initialize_constructs(db)
                return
            db.executescript("""
                CREATE TABLE experiment_projects (
                    id TEXT PRIMARY KEY, source_job_id TEXT NOT NULL, name TEXT NOT NULL,
                    reference_sequence TEXT NOT NULL, status TEXT NOT NULL,
                    revision INTEGER NOT NULL, created_at TEXT NOT NULL
                );
                CREATE TABLE project_revisions (
                    project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                    revision INTEGER NOT NULL, reason TEXT NOT NULL,
                    snapshot_json TEXT NOT NULL, created_at TEXT NOT NULL,
                    PRIMARY KEY(project_id, revision)
                );
                CREATE TABLE artifacts (
                    id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                    relative_path TEXT NOT NULL, kind TEXT NOT NULL, sha256 TEXT NOT NULL,
                    size_bytes INTEGER NOT NULL, media_type TEXT NOT NULL,
                    read_only INTEGER NOT NULL, created_at TEXT NOT NULL,
                    UNIQUE(project_id, relative_path)
                );
                CREATE TABLE experiment_tasks (
                    id TEXT PRIMARY KEY, project_id TEXT REFERENCES experiment_projects(id),
                    kind TEXT NOT NULL, status TEXT NOT NULL, progress REAL NOT NULL,
                    input_json TEXT NOT NULL, output_json TEXT, error_code TEXT,
                    created_at TEXT NOT NULL, updated_at TEXT NOT NULL
                );
                CREATE TABLE audit_events (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                    entity_type TEXT NOT NULL, entity_id TEXT NOT NULL,
                    action TEXT NOT NULL, actor TEXT NOT NULL,
                    before_revision INTEGER, after_revision INTEGER,
                    details_json TEXT NOT NULL, created_at TEXT NOT NULL
                );
            """)
            db.execute("INSERT INTO schema_migrations VALUES (1, ?)", (_now(),))
            db.commit()
            self._initialize_constructs(db)

    def _initialize_constructs(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=2").fetchone():
            self._initialize_samples(db)
            return
        db.executescript("""
            CREATE TABLE antigen_constructs (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                name TEXT NOT NULL, source TEXT NOT NULL, sequence TEXT NOT NULL,
                sequence_sha256 TEXT NOT NULL, mutations_json TEXT NOT NULL,
                mapping_status TEXT NOT NULL, risk_flags_json TEXT NOT NULL,
                version INTEGER NOT NULL, project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE TABLE antibody_constructs (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                name TEXT NOT NULL, parent_id TEXT REFERENCES antibody_constructs(id),
                role TEXT NOT NULL, vh_sequence TEXT NOT NULL, vl_sequence TEXT NOT NULL,
                heavy_mutations_json TEXT NOT NULL, light_mutations_json TEXT NOT NULL,
                mutation_class TEXT NOT NULL, version INTEGER NOT NULL,
                project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE TABLE assay_batches (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                name TEXT NOT NULL, date TEXT NOT NULL, assay TEXT NOT NULL,
                cell_line TEXT NOT NULL, dose REAL, dose_unit TEXT,
                positive_gate_id TEXT, independent_batch INTEGER NOT NULL,
                settings_json TEXT NOT NULL, project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE INDEX idx_antigen_project ON antigen_constructs(project_id);
            CREATE INDEX idx_antibody_project ON antibody_constructs(project_id);
            CREATE INDEX idx_batch_project ON assay_batches(project_id);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (2, ?)", (_now(),))
        db.commit()
        self._initialize_samples(db)

    def _initialize_samples(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=3").fetchone():
            self._initialize_flowjo(db)
            return
        db.executescript("""
            CREATE TABLE import_sessions (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                batch_id TEXT REFERENCES assay_batches(id), status TEXT NOT NULL,
                project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE TABLE pending_uploads (
                id TEXT PRIMARY KEY, session_id TEXT NOT NULL REFERENCES import_sessions(id),
                project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                artifact_id TEXT NOT NULL REFERENCES artifacts(id), original_filename TEXT NOT NULL,
                sha256 TEXT NOT NULL, metadata_json TEXT NOT NULL, match_status TEXT NOT NULL,
                created_at TEXT NOT NULL
            );
            CREATE TABLE samples (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                batch_id TEXT NOT NULL REFERENCES assay_batches(id),
                artifact_id TEXT NOT NULL REFERENCES artifacts(id), role TEXT NOT NULL,
                antigen_construct_id TEXT REFERENCES antigen_constructs(id),
                antibody_construct_id TEXT REFERENCES antibody_constructs(id),
                replicate_type TEXT NOT NULL, replicate_index INTEGER NOT NULL,
                condition_json TEXT NOT NULL, source_upload_id TEXT REFERENCES pending_uploads(id),
                created_at TEXT NOT NULL,
                UNIQUE(batch_id, artifact_id, role, replicate_type, replicate_index)
            );
            CREATE INDEX idx_upload_project ON pending_uploads(project_id,session_id);
            CREATE INDEX idx_sample_project ON samples(project_id,batch_id);
            CREATE UNIQUE INDEX idx_project_raw_fcs_hash ON artifacts(project_id,sha256) WHERE kind='raw_fcs';
        """)
        db.execute("INSERT INTO schema_migrations VALUES (3, ?)", (_now(),))
        db.commit()
        self._initialize_flowjo(db)

    def _initialize_flowjo(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=4").fetchone():
            self._initialize_gate_revisions(db)
            return
        db.executescript("""
            CREATE TABLE flowjo_imports (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                workspace_artifact_id TEXT NOT NULL REFERENCES artifacts(id),
                preview_artifact_id TEXT NOT NULL REFERENCES artifacts(id),
                status TEXT NOT NULL, unresolved_json TEXT NOT NULL,
                project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE INDEX idx_flowjo_project ON flowjo_imports(project_id,created_at);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (4, ?)", (_now(),))
        db.commit()
        self._initialize_gate_revisions(db)

    def _initialize_gate_revisions(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=5").fetchone():
            self._initialize_gate_analyses(db)
            return
        db.executescript("""
            CREATE TABLE gate_revisions (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                parent_revision_id TEXT REFERENCES gate_revisions(id), name TEXT NOT NULL,
                version INTEGER NOT NULL, gates_json TEXT NOT NULL,
                compensation_decision TEXT NOT NULL, compensation_artifact_id TEXT REFERENCES artifacts(id),
                project_revision INTEGER NOT NULL, created_at TEXT NOT NULL,
                UNIQUE(project_id,version)
            );
            CREATE INDEX idx_gate_project ON gate_revisions(project_id,version);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (5, ?)", (_now(),))
        db.commit()
        self._initialize_gate_analyses(db)

    def _initialize_gate_analyses(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=6").fetchone():
            self._initialize_qc_reviews(db)
            return
        db.executescript("""
            CREATE TABLE gate_analyses (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                gate_revision_id TEXT NOT NULL REFERENCES gate_revisions(id),
                input_artifact_ids_json TEXT NOT NULL,
                result_artifact_id TEXT NOT NULL REFERENCES artifacts(id),
                status TEXT NOT NULL, project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE INDEX idx_gate_analysis_project ON gate_analyses(project_id,created_at);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (6, ?)", (_now(),))
        db.commit()
        self._initialize_qc_reviews(db)

    def _initialize_qc_reviews(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=7").fetchone():
            self._initialize_qc_review_sources(db)
            return
        db.executescript("""
            CREATE TABLE qc_reviews (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                gate_analysis_id TEXT NOT NULL REFERENCES gate_analyses(id),
                sample_id TEXT NOT NULL REFERENCES samples(id), version INTEGER NOT NULL,
                status TEXT NOT NULL, reason TEXT NOT NULL, reviewer TEXT NOT NULL,
                project_revision INTEGER NOT NULL, created_at TEXT NOT NULL,
                UNIQUE(gate_analysis_id,sample_id,version)
            );
            CREATE INDEX idx_qc_review_project ON qc_reviews(project_id,gate_analysis_id,sample_id,version);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (7, ?)", (_now(),))
        db.commit()
        self._initialize_qc_review_sources(db)

    def _initialize_qc_review_sources(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=8").fetchone():
            self._initialize_effect_analyses(db)
            return
        db.execute("ALTER TABLE qc_reviews RENAME TO qc_reviews_v7")
        db.execute("DROP INDEX IF EXISTS idx_qc_review_project")
        db.executescript("""
            CREATE TABLE qc_reviews_v8 (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                source_type TEXT NOT NULL CHECK(source_type IN ('GATE_ANALYSIS','FLOWJO_IMPORT')),
                source_id TEXT NOT NULL, sample_id TEXT NOT NULL REFERENCES samples(id),
                version INTEGER NOT NULL, status TEXT NOT NULL, reason TEXT NOT NULL,
                reviewer TEXT NOT NULL, project_revision INTEGER NOT NULL, created_at TEXT NOT NULL,
                UNIQUE(source_type,source_id,sample_id,version)
            );
            INSERT INTO qc_reviews_v8
                (id,project_id,source_type,source_id,sample_id,version,status,reason,reviewer,project_revision,created_at)
                SELECT id,project_id,'GATE_ANALYSIS',gate_analysis_id,sample_id,version,status,reason,reviewer,project_revision,created_at
                FROM qc_reviews_v7;
            DROP TABLE qc_reviews_v7;
            ALTER TABLE qc_reviews_v8 RENAME TO qc_reviews;
            CREATE INDEX idx_qc_review_project ON qc_reviews(project_id,source_type,source_id,sample_id,version);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (8, ?)", (_now(),))
        db.commit()
        self._initialize_effect_analyses(db)

    def _initialize_effect_analyses(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=9").fetchone():
            self._initialize_paraplume(db)
            return
        db.executescript("""
            CREATE TABLE effect_analyses (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                gate_analysis_id TEXT NOT NULL REFERENCES gate_analyses(id),
                positive_gate_id TEXT NOT NULL, result_artifact_id TEXT NOT NULL REFERENCES artifacts(id),
                status TEXT NOT NULL, project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE INDEX idx_effect_analysis_project ON effect_analyses(project_id,created_at);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (9, ?)", (_now(),))
        db.commit()
        self._initialize_paraplume(db)

    def _initialize_paraplume(self, db: sqlite3.Connection) -> None:
        if not db.execute("SELECT 1 FROM schema_migrations WHERE version=10").fetchone():
            db.executescript("""
                CREATE TABLE paraplume_cache (
                    cache_key TEXT PRIMARY KEY, result_json TEXT NOT NULL,
                    provenance_json TEXT NOT NULL, created_at TEXT NOT NULL
                );
                CREATE TABLE paraplume_runs (
                    id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                    task_id TEXT NOT NULL REFERENCES experiment_tasks(id),
                    antibody_construct_id TEXT NOT NULL REFERENCES antibody_constructs(id),
                    cache_key TEXT NOT NULL, result_artifact_id TEXT NOT NULL REFERENCES artifacts(id),
                    cache_hit INTEGER NOT NULL, created_at TEXT NOT NULL
                );
                CREATE INDEX idx_paraplume_run_project ON paraplume_runs(project_id,created_at);
                CREATE INDEX idx_paraplume_run_task ON paraplume_runs(task_id);
            """)
            db.execute("INSERT INTO schema_migrations VALUES (10, ?)", (_now(),))
            db.commit()
        columns = {row[1] for row in db.execute("PRAGMA table_info(experiment_tasks)")}
        if "owner_pid" not in columns:
            db.execute("ALTER TABLE experiment_tasks ADD COLUMN owner_pid INTEGER")
            db.commit()
        active = db.execute(
            "SELECT id,owner_pid FROM experiment_tasks WHERE kind='PARAPLUME' AND status IN ('QUEUED','RUNNING')"
        ).fetchall()
        for task in active:
            pid = task["owner_pid"]
            if pid is not None and _process_alive(int(pid)):
                continue
            db.execute(
                """UPDATE experiment_tasks SET status='INTERRUPTED',
                   error_code='TASK_INTERRUPTED_BY_RESTART',output_json=?,updated_at=? WHERE id=?""",
                (
                    json.dumps({"message": "Task was interrupted by a service restart"}),
                    _now(),
                    task["id"],
                ),
            )
        db.commit()
        self._initialize_constraint_exports(db)

    def _initialize_constraint_exports(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=11").fetchone():
            self._initialize_summary_measurements(db)
            return
        db.executescript("""
            CREATE TABLE constraint_exports (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                antibody_construct_id TEXT NOT NULL REFERENCES antibody_constructs(id),
                artifact_id TEXT NOT NULL REFERENCES artifacts(id),
                merge_sources_json TEXT NOT NULL, constraint_count INTEGER NOT NULL,
                status TEXT NOT NULL, project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE INDEX idx_constraint_exports_project ON constraint_exports(project_id,created_at);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (11, ?)", (_now(),))
        db.commit()
        self._initialize_summary_measurements(db)

    def _initialize_summary_measurements(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=12").fetchone():
            self._initialize_flowjo_percentage_imports(db)
            return
        db.executescript("""
            CREATE TABLE summary_measurements (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                assay TEXT NOT NULL CHECK(assay IN ('ANTIGEN_MYC','FAB_MYC')),
                antigen_construct_id TEXT REFERENCES antigen_constructs(id),
                antibody_construct_id TEXT REFERENCES antibody_constructs(id),
                measured_on TEXT NOT NULL, metric TEXT NOT NULL, value REAL NOT NULL,
                unit TEXT NOT NULL, source_label TEXT NOT NULL, notes TEXT NOT NULL,
                reviewer TEXT NOT NULL, historical INTEGER NOT NULL CHECK(historical=1),
                project_revision INTEGER NOT NULL, created_at TEXT NOT NULL,
                CHECK((assay='ANTIGEN_MYC' AND antigen_construct_id IS NOT NULL
                       AND antibody_construct_id IS NULL)
                   OR (assay='FAB_MYC' AND antibody_construct_id IS NOT NULL
                       AND antigen_construct_id IS NULL))
            );
            CREATE INDEX idx_summary_measurement_project
                ON summary_measurements(project_id,measured_on,created_at);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (12, ?)", (_now(),))
        db.commit()
        self._initialize_flowjo_percentage_imports(db)

    def _initialize_flowjo_percentage_imports(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=13").fetchone():
            self._initialize_percentage_qc_sources(db)
            return
        db.executescript("""
            CREATE TABLE flowjo_percentage_imports (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                batch_id TEXT NOT NULL REFERENCES assay_batches(id),
                source_artifact_id TEXT NOT NULL REFERENCES artifacts(id),
                result_artifact_id TEXT NOT NULL REFERENCES artifacts(id),
                flowjo_import_id TEXT REFERENCES flowjo_imports(id),
                status TEXT NOT NULL, project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE INDEX idx_flowjo_percentage_project
                ON flowjo_percentage_imports(project_id,created_at);
            CREATE TABLE flowjo_percentage_effect_analyses (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                import_id TEXT NOT NULL REFERENCES flowjo_percentage_imports(id),
                result_artifact_id TEXT NOT NULL REFERENCES artifacts(id),
                status TEXT NOT NULL, project_revision INTEGER NOT NULL, created_at TEXT NOT NULL
            );
            CREATE INDEX idx_flowjo_percentage_effect_project
                ON flowjo_percentage_effect_analyses(project_id,created_at);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (13, ?)", (_now(),))
        db.commit()
        self._initialize_percentage_qc_sources(db)

    def _initialize_percentage_qc_sources(self, db: sqlite3.Connection) -> None:
        if db.execute("SELECT 1 FROM schema_migrations WHERE version=14").fetchone():
            return
        db.execute("ALTER TABLE qc_reviews RENAME TO qc_reviews_v13")
        db.execute("DROP INDEX IF EXISTS idx_qc_review_project")
        db.executescript("""
            CREATE TABLE qc_reviews_v14 (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL REFERENCES experiment_projects(id),
                source_type TEXT NOT NULL CHECK(source_type IN (
                    'GATE_ANALYSIS','FLOWJO_IMPORT','FLOWJO_PERCENTAGE_IMPORT'
                )),
                source_id TEXT NOT NULL, sample_id TEXT NOT NULL REFERENCES samples(id),
                version INTEGER NOT NULL, status TEXT NOT NULL, reason TEXT NOT NULL,
                reviewer TEXT NOT NULL, project_revision INTEGER NOT NULL, created_at TEXT NOT NULL,
                UNIQUE(source_type,source_id,sample_id,version)
            );
            INSERT INTO qc_reviews_v14
                (id,project_id,source_type,source_id,sample_id,version,status,reason,reviewer,project_revision,created_at)
                SELECT id,project_id,source_type,source_id,sample_id,version,status,reason,reviewer,project_revision,created_at
                FROM qc_reviews_v13;
            DROP TABLE qc_reviews_v13;
            ALTER TABLE qc_reviews_v14 RENAME TO qc_reviews;
            CREATE INDEX idx_qc_review_project ON qc_reviews(project_id,source_type,source_id,sample_id,version);
        """)
        db.execute("INSERT INTO schema_migrations VALUES (14, ?)", (_now(),))
        db.commit()

    def _recover_published_project_snapshots(self, db: sqlite3.Connection) -> None:
        """Restore a published project directory if the process died before DB commit."""
        schema_path = (
            Path(__file__).resolve().parent
            / "schemas"
            / "experiment_project_manifest.schema.json"
        )
        validator = Draft202012Validator(read_json(schema_path))
        required_kinds = {
            "sequence", "structure", "mapping", "sectors", "surface_graph",
            "v4_manifest", "recommendations",
        }
        for folder in sorted(self.projects_root.iterdir()):
            project_id = folder.name
            if (
                folder.is_symlink()
                or not folder.is_dir()
                or len(project_id) != 16
                or any(character not in "0123456789abcdef" for character in project_id)
            ):
                continue
            if db.execute(
                "SELECT 1 FROM experiment_projects WHERE id=?", (project_id,)
            ).fetchone():
                continue
            manifest_path = folder / "project_manifest.json"
            try:
                if manifest_path.is_symlink() or not manifest_path.is_file():
                    continue
                manifest = read_json(manifest_path)
                validator.validate(manifest)
                if manifest["project_id"] != project_id:
                    continue
                artifacts = manifest["artifacts"]
                artifact_ids = [item["id"] for item in artifacts]
                relative_paths = [item["relative_path"] for item in artifacts]
                kinds = [item["kind"] for item in artifacts]
                if (
                    len(set(artifact_ids)) != len(artifact_ids)
                    or len(set(relative_paths)) != len(relative_paths)
                    or not required_kinds.issubset(kinds)
                    or len(set(kinds)) != len(kinds)
                ):
                    continue
                folder_root = folder.resolve()
                for item in artifacts:
                    artifact_path = folder / item["relative_path"]
                    if (
                        artifact_path.is_symlink()
                        or not artifact_path.is_file()
                        or not artifact_path.resolve().is_relative_to(folder_root)
                        or artifact_path.stat().st_size != item["size_bytes"]
                        or sha256_file(artifact_path) != item["sha256"]
                    ):
                        raise IntegrityError("Published snapshot artifact failed recovery verification")
                if any(
                    manifest["frozen_hashes"].get(kind)
                    != next(item["sha256"] for item in artifacts if item["kind"] == kind)
                    for kind in ("structure", "mapping", "sectors", "surface_graph")
                ):
                    continue
            except (OSError, ValueError, TypeError, KeyError, JsonSchemaValidationError):
                # Leave invalid directories intact for manual forensic recovery.
                continue

            db.execute("SAVEPOINT recover_project_snapshot")
            try:
                created_at = manifest["created_at"]
                db.execute(
                    "INSERT INTO experiment_projects VALUES (?,?,?,?,?,?,?)",
                    (
                        project_id, manifest["source_job_id"], manifest["name"],
                        manifest["reference_sequence"], "COMPLETE", 1, created_at,
                    ),
                )
                db.execute(
                    "INSERT INTO project_revisions VALUES (?,?,?,?,?)",
                    (
                        project_id, 1, "INITIAL_SNAPSHOT",
                        json.dumps(manifest, ensure_ascii=False), created_at,
                    ),
                )
                db.executemany(
                    "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                    [
                        (
                            item["id"], project_id, item["relative_path"], item["kind"],
                            item["sha256"], item["size_bytes"], item["media_type"], 1,
                            item["created_at"],
                        )
                        for item in artifacts
                    ],
                )
                db.execute(
                    """INSERT INTO audit_events
                       (project_id,entity_type,entity_id,action,actor,before_revision,
                        after_revision,details_json,created_at)
                       VALUES (?,?,?,?,?,?,?,?,?)""",
                    (
                        project_id, "project", project_id, "CREATE_RECOVERED", "startup_recovery",
                        None, 1,
                        json.dumps({"source_job_id": manifest["source_job_id"]}), created_at,
                    ),
                )
                db.execute("RELEASE SAVEPOINT recover_project_snapshot")
            except Exception:
                db.execute("ROLLBACK TO SAVEPOINT recover_project_snapshot")
                db.execute("RELEASE SAVEPOINT recover_project_snapshot")
                raise

    def _reconcile_registered_project_manifests(self, db: sqlite3.Connection) -> None:
        """Rebuild derived manifest artifact lists from committed registry rows."""
        schema_path = (
            Path(__file__).resolve().parent
            / "schemas"
            / "experiment_project_manifest.schema.json"
        )
        validator = Draft202012Validator(read_json(schema_path))
        immutable_fields = (
            "schema_version", "project_id", "name", "source_job_id",
            "source_v4_schema_version", "source_software", "revision",
            "reference_sequence", "frozen_hashes", "created_at",
        )
        projects = db.execute(
            "SELECT id,revision FROM experiment_projects ORDER BY id"
        ).fetchall()
        for project in projects:
            project_id = project["id"]
            folder = self.projects_root / project_id
            manifest_path = folder / "project_manifest.json"
            if folder.is_symlink() or not folder.is_dir() or manifest_path.is_symlink():
                continue
            initial = db.execute(
                "SELECT snapshot_json FROM project_revisions WHERE project_id=? AND revision=1",
                (project_id,),
            ).fetchone()
            if initial is None:
                continue
            expected = json.loads(initial["snapshot_json"])
            artifacts = db.execute(
                "SELECT * FROM artifacts WHERE project_id=? ORDER BY rowid", (project_id,)
            ).fetchall()
            expected["artifacts"] = [
                {
                    "id": item["id"],
                    "kind": item["kind"],
                    "relative_path": item["relative_path"],
                    "sha256": item["sha256"],
                    "size_bytes": item["size_bytes"],
                    "media_type": item["media_type"],
                    "read_only": bool(item["read_only"]),
                    "created_at": item["created_at"],
                }
                for item in artifacts
            ]
            validator.validate(expected)
            try:
                current = read_json(manifest_path)
            except (OSError, ValueError, TypeError):
                current = None
            if current == expected:
                continue
            if isinstance(current, dict) and any(
                current.get(field) != expected.get(field) for field in immutable_fields
            ):
                # Immutable identity/hash mismatches require integrity review.
                continue

            folder_root = folder.resolve()
            valid_artifacts = True
            for item in artifacts:
                target = folder / item["relative_path"]
                if (
                    target.is_symlink()
                    or not target.is_file()
                    or not target.resolve().is_relative_to(folder_root)
                    or target.stat().st_size != item["size_bytes"]
                    or sha256_file(target) != item["sha256"]
                ):
                    valid_artifacts = False
                    break
            if not valid_artifacts:
                continue

            try:
                observed_hash = sha256_file(manifest_path) if manifest_path.is_file() else None
            except OSError:
                observed_hash = None
            now = _now()
            db.execute(
                """INSERT INTO audit_events
                   (project_id,entity_type,entity_id,action,actor,before_revision,
                    after_revision,details_json,created_at)
                   VALUES (?,?,?,?,?,?,?,?,?)""",
                (
                    project_id, "project_manifest", project_id, "MANIFEST_RECOVERED",
                    "startup_recovery", project["revision"], project["revision"],
                    json.dumps({
                        "observed_manifest_sha256": observed_hash,
                        "restored_artifact_count": len(artifacts),
                    }),
                    now,
                ),
            )
            db.commit()
            write_json(manifest_path, expected)

    @staticmethod
    def _ecd_from_fasta(path: Path) -> str:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
        if not lines or not lines[0].startswith(">"):
            raise SnapshotError("sequences.fasta has no first FASTA record")
        sequence_lines = []
        for line in lines[1:]:
            if line.startswith(">"):
                break
            sequence_lines.append(line.strip())
        sequence = "".join(sequence_lines)
        # The first record of pipeline sequences.fasta is WT_ECD.
        if not sequence or set(sequence.upper()) - set("ACDEFGHIKLMNPQRSTVWY"):
            raise SnapshotError("WT_ECD FASTA sequence is invalid")
        return sequence.upper()

    @staticmethod
    def _ecd_positions_from_fasta(path: Path, sequence: str) -> set[int] | None:
        first_line = path.read_text(encoding="utf-8-sig").splitlines()[0]
        match = re.search(r":(\d+)-(\d+)\s*$", first_line)
        if not match:
            return None
        start, end = map(int, match.groups())
        if start < 1 or end < start or end - start + 1 != len(sequence):
            raise SnapshotError("WT_ECD FASTA coordinate range does not match its sequence")
        return set(range(start, end + 1))

    def _source_files(self, job: dict[str, Any]) -> tuple[dict[str, Path], str, dict[str, Any]]:
        if job["status"] != "COMPLETE" or job["request"].get("design_engine") != "joint_v4_ml":
            raise SnapshotError("Source must be a completed joint_v4_ml job")
        source = Path(job["workdir"]).resolve()
        if not source.is_relative_to(self.jobs.root.resolve()) or source.name != job["id"]:
            raise SnapshotError("Source job directory is outside the job store")
        fixed = {
            "sequence": source / "sequences.fasta",
            "mapping": source / "residue_mapping.csv",
            "sectors": source / "probe_sectors_v4.json",
            "surface_graph": source / "surface_graph_v4.json",
            "v4_manifest": source / "v4_manifest.json",
            "recommendations": source / "designs.json",
        }
        missing = [key for key, path in fixed.items() if not path.is_file() or path.is_symlink()]
        if missing:
            raise SnapshotError(f"Required v4 artifacts missing: {', '.join(missing)}")
        manifest = read_json(fixed["v4_manifest"])
        graph = read_json(fixed["surface_graph"])
        sectors = read_json(fixed["sectors"])
        if manifest.get("schema_version") != "4.0" or manifest.get("engine") != "joint_v4_ml":
            raise SnapshotError("Invalid source v4 manifest")
        if graph.get("schema_version") != "1.0" or graph.get("position_system") != "uniprot_pos":
            raise SnapshotError("Invalid v4 surface graph")
        if (
            not isinstance(graph.get("nodes"), list)
            or not isinstance(graph.get("edges"), list)
            or not isinstance(sectors, list)
        ):
            raise SnapshotError("Invalid graph or sectors format")
        if graph.get("mapping_sha256") != sha256_file(fixed["mapping"]):
            raise SnapshotError("Surface graph mapping hash differs from source mapping")
        recommendations = read_json(fixed["recommendations"])
        if not isinstance(recommendations, list) or any(
            not isinstance(item, dict) or not item.get("id")
            or not isinstance(item.get("mutations"), list)
            or not item.get("mutant_ecd_sequence")
            for item in recommendations
        ):
            raise SnapshotError("Invalid recommended antigen constructs")
        graph_hash = ((manifest.get("provenance") or {}).get("surface_graph") or {}).get(
            "snapshot_sha256"
        )
        if graph_hash != sha256_file(fixed["surface_graph"]):
            raise SnapshotError("V4 manifest surface graph hash differs from file")
        sector_owners = {
            position: sector.get("sector_id")
            for sector in sectors
            for position in sector.get("residue_positions", [])
        }
        graph_owners = {node.get("uniprot_pos"): node.get("sector_id") for node in graph["nodes"]}
        if graph_owners != sector_owners or len(graph_owners) != len(graph["nodes"]):
            raise SnapshotError("Graph node and sector ownership differ")
        structure_hash = graph.get("structure_sha256")
        if not isinstance(structure_hash, str) or len(structure_hash) != 64:
            raise SnapshotError("Surface graph lacks a structure hash")
        structure = next(
            (
                p
                for p in (
                    source / "structure_for_energy.pdb",
                    source / "structure.cif",
                    source / "ecd_structure.cif",
                )
                if p.is_file() and not p.is_symlink() and sha256_file(p) == structure_hash
            ),
            None,
        )
        if structure is None:
            raise SnapshotError("No source structure matches the surface graph hash")
        fixed["structure"] = structure
        sequence = self._ecd_from_fasta(fixed["sequence"])
        sequence_positions = self._ecd_positions_from_fasta(fixed["sequence"], sequence)
        analysis_path = source / "analysis.json"
        if analysis_path.is_file():
            expected = (
                (read_json(analysis_path).get("protein") or {}).get("ecd_sequence") or ""
            ).upper()
            if expected and expected != sequence:
                raise SnapshotError("WT_ECD sequence differs from source analysis")
        run_manifest_path = source / "run_manifest.json"
        if run_manifest_path.is_file() and not run_manifest_path.is_symlink():
            run_manifest = read_json(run_manifest_path)
            declared_files = run_manifest.get("files") or {}
            evidence_files = {
                "residue_reliability": "residue_reliability.csv",
                "sasa_comparison": "sasa_comparison.csv",
                "reliability_qc": "reliability_qc.json",
            }
            declared_evidence = set(evidence_files.values()) & set(declared_files)
            if declared_evidence and declared_evidence != set(evidence_files.values()):
                raise SnapshotError("Source run manifest declares only part of structure evidence")
            if declared_evidence:
                for key, filename in evidence_files.items():
                    path = source / filename
                    if (
                        not path.is_file()
                        or path.is_symlink()
                        or not isinstance(declared_files[filename], str)
                        or sha256_file(path) != declared_files[filename]
                    ):
                        raise SnapshotError(f"Source structure evidence hash differs: {filename}")
                    fixed[key] = path
                self._validate_structure_evidence_files(
                    fixed, fixed["mapping"], sequence_positions=sequence_positions
                )
        return fixed, sequence, manifest

    @staticmethod
    def _validate_structure_evidence_files(
        evidence_files: dict[str, Path],
        mapping_path: Path,
        *,
        sequence_positions: set[int] | None = None,
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
        required = {"residue_reliability", "sasa_comparison", "reliability_qc"}
        if not required.issubset(evidence_files):
            raise SnapshotError("Source structure reliability/SASA evidence is incomplete")
        with mapping_path.open(newline="", encoding="utf-8-sig") as handle:
            mapping_rows = list(csv.DictReader(handle))
        mapped_positions = {
            int(row["uniprot_pos"])
            for row in mapping_rows
            if row.get("uniprot_pos") and row.get("chain_id") and row.get("auth_seq_id")
        }
        if sequence_positions is not None and not mapped_positions.issubset(sequence_positions):
            raise SnapshotError("Frozen mapping contains positions outside the WT_ECD sequence")
        reliability_required = {
            "residue_id", "state", "plddt", "mean_local_pae", "mean_core_pae",
            "quality_flags", "use_for_surface_mesh", "use_for_3d_footprint",
            "use_for_epitope_screening",
        }
        with evidence_files["residue_reliability"].open(
            newline="", encoding="utf-8-sig"
        ) as handle:
            reliability_rows = list(csv.DictReader(handle))
        if not reliability_rows or not reliability_required.issubset(reliability_rows[0]):
            raise SnapshotError("Residue reliability table is missing required columns")
        reliability_positions = []
        reliability_by_position: dict[int, dict[str, str]] = {}
        observed_state_counts: dict[str, int] = {}
        for record in reliability_rows:
            try:
                position = int(record["residue_id"])
            except (KeyError, TypeError, ValueError) as exc:
                raise SnapshotError("Residue reliability table has an invalid residue_id") from exc
            if record["state"] not in {
                "STRUCTURED", "FLEXIBLE_STRUCTURED", "UNCERTAIN", "DISORDERED",
                "UNRESOLVED", "CONDITIONALLY_FOLDED",
            }:
                raise SnapshotError(f"Unsupported residue reliability state at {position}")
            if position in reliability_by_position:
                raise SnapshotError(f"Duplicate residue reliability record at {position}")
            reliability_by_position[position] = record
            observed_state_counts[record["state"]] = observed_state_counts.get(record["state"], 0) + 1
            for column in ("plddt", "mean_local_pae", "mean_core_pae"):
                value = record.get(column, "")
                if value:
                    try:
                        numeric = float(value)
                    except (TypeError, ValueError) as exc:
                        raise SnapshotError(f"Invalid {column} at residue {position}") from exc
                    if not math.isfinite(numeric):
                        raise SnapshotError(f"Non-finite {column} at residue {position}")
            try:
                flags = json.loads(record.get("quality_flags") or "[]")
            except json.JSONDecodeError as exc:
                raise SnapshotError(f"Invalid quality_flags at residue {position}") from exc
            if not isinstance(flags, list) or any(not isinstance(flag, str) for flag in flags):
                raise SnapshotError(f"Invalid quality_flags at residue {position}")
            for column in (
                "use_for_surface_mesh", "use_for_3d_footprint", "use_for_epitope_screening"
            ):
                if str(record.get(column) or "").lower() not in {"true", "false"}:
                    raise SnapshotError(f"Invalid {column} at residue {position}")
            reliability_positions.append(position)
        reliability_position_set = set(reliability_positions)
        if not mapped_positions.issubset(reliability_position_set):
            raise SnapshotError("Residue reliability is missing frozen mapped positions")
        unmapped_reliability_positions = reliability_position_set - mapped_positions
        if (
            sequence_positions is not None
            and not unmapped_reliability_positions.issubset(sequence_positions)
        ):
            raise SnapshotError("Unmapped reliability positions are outside the WT_ECD sequence")
        for position in unmapped_reliability_positions:
            record = reliability_by_position[position]
            if (
                record["state"] != "UNRESOLVED"
                or str(record.get("use_for_surface_mesh") or "").lower() != "false"
                or str(record.get("use_for_3d_footprint") or "").lower() != "false"
            ):
                raise SnapshotError(
                    f"Unmapped reliability position {position} must be UNRESOLVED and excluded from structural geometry"
                )

        sasa_required = {"position", "native_rsasa", "procaliper_rsasa", "absolute_delta", "status"}
        with evidence_files["sasa_comparison"].open(newline="", encoding="utf-8-sig") as handle:
            sasa_rows = list(csv.DictReader(handle))
        if not sasa_rows or not sasa_required.issubset(sasa_rows[0]):
            raise SnapshotError("SASA comparison table is missing required columns")
        sasa_positions = []
        for record in sasa_rows:
            try:
                position = int(record["position"])
            except (KeyError, TypeError, ValueError) as exc:
                raise SnapshotError("SASA comparison table has an invalid position") from exc
            if record["status"] not in {
                "AGREE",
                "SASA_DISAGREEMENT",
                "NOT_COMPARED",
                "COMPARISON_CONTEXT_UNVERIFIED",
                "METHOD_SETTINGS_UNVERIFIED",
            }:
                raise SnapshotError(f"Unsupported SASA comparison status at {position}")
            for column in (
                "native_rsasa",
                "normalized_context_rsasa",
                "procaliper_rsasa",
                "absolute_delta",
                "native_context_absolute_delta",
            ):
                value = record.get(column, "")
                if value:
                    try:
                        numeric = float(value)
                    except (TypeError, ValueError) as exc:
                        raise SnapshotError(f"Invalid {column} at position {position}") from exc
                    if not math.isfinite(numeric) or numeric < 0:
                        raise SnapshotError(f"Invalid {column} at position {position}")
            context_status = record.get("native_context_status", "")
            if context_status and context_status not in {
                "CONTEXT_STABLE",
                "CONTEXT_SENSITIVE",
                "NOT_COMPARED",
            }:
                raise SnapshotError(f"Unsupported native SASA context status at {position}")
            method_status = record.get("comparison_method_status", "")
            if method_status and method_status not in {
                "MATCHED",
                "VERSION_DIFFERS",
                "SETTINGS_DIFFER",
                "PROVENANCE_INCOMPLETE",
                "INPUT_SCOPE_UNVERIFIED",
                "NOT_COMPARED",
            }:
                raise SnapshotError(f"Unsupported SASA method provenance status at {position}")
            for column in (
                "native_context_classification_changed",
                "procaliper_classification_changed",
            ):
                value = record.get(column, "")
                if value and str(value).lower() not in {"true", "false"}:
                    raise SnapshotError(f"Invalid {column} at position {position}")
            sasa_positions.append(position)
        if len(set(sasa_positions)) != len(sasa_positions) or set(sasa_positions) != mapped_positions:
            raise SnapshotError("SASA comparison positions do not match the frozen mapping")

        reliability_qc = read_json(evidence_files["reliability_qc"])
        state_counts = reliability_qc.get("state_counts")
        if (
            reliability_qc.get("classification_version") != "reliability_rules_v1"
            or not isinstance(state_counts, dict)
            or state_counts != observed_state_counts
            or not isinstance(reliability_qc.get("predicted"), bool)
        ):
            raise SnapshotError("Reliability QC summary does not match per-residue classifications")
        return reliability_rows, sasa_rows, reliability_qc

    def create(self, source_job_id: str, name: str) -> dict[str, Any]:
        project_id = uuid.uuid4().hex[:16]
        staging = self.projects_root / f".staging-{project_id}"
        final = self.projects_root / project_id
        try:
            return self._create_snapshot(source_job_id, name, project_id)
        except Exception:
            shutil.rmtree(staging, ignore_errors=True)
            try:
                self._row(project_id)
            except KeyError:
                shutil.rmtree(final, ignore_errors=True)
            except Exception:
                # Keep a published snapshot if database availability is uncertain.
                pass
            raise

    def _create_snapshot(
        self, source_job_id: str, name: str, project_id: str
    ) -> dict[str, Any]:
        job = self.jobs.get(source_job_id)
        sources, sequence, v4_manifest = self._source_files(job)
        final = self.projects_root / project_id
        staging = self.projects_root / f".staging-{project_id}"
        paths = {
            "sequence": "reference/antigen_sequence.fasta",
            "structure": f"reference/structure{sources['structure'].suffix.lower()}",
            "mapping": "reference/residue_mapping.csv",
            "sectors": "reference/probe_sectors_v4.json",
            "surface_graph": "reference/surface_graph_v4.json",
            "v4_manifest": "reference/source_v4_manifest.json",
            "recommendations": "constructs/recommended_antigens.json",
        }
        media = {
            "sequence": "text/plain",
            "structure": "chemical/x-pdb"
            if sources["structure"].suffix == ".pdb"
            else "chemical/x-mmcif",
            "mapping": "text/csv",
            "sectors": "application/json",
            "surface_graph": "application/json",
            "v4_manifest": "application/json",
            "recommendations": "application/json",
        }
        optional_structure_evidence = {
            "residue_reliability": ("reference/residue_reliability.csv", "text/csv"),
            "sasa_comparison": ("reference/sasa_comparison.csv", "text/csv"),
            "reliability_qc": ("reference/reliability_qc.json", "application/json"),
        }
        for kind, (relative, media_type) in optional_structure_evidence.items():
            if kind in sources:
                paths[kind] = relative
                media[kind] = media_type
        now = _now()
        artifacts = []
        with self._lock:
            staging.mkdir()
            for kind, relative in paths.items():
                target = staging / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                temp = target.with_name(target.name + ".tmp")
                shutil.copyfile(sources[kind], temp)
                digest = sha256_file(sources[kind])
                if sha256_file(temp) != digest:
                    raise IntegrityError(f"Snapshot copy failed hash verification: {kind}")
                os.replace(temp, target)
                artifacts.append(
                    {
                        "id": uuid.uuid4().hex,
                        "kind": kind,
                        "relative_path": relative,
                        "sha256": digest,
                        "size_bytes": target.stat().st_size,
                        "media_type": media[kind],
                        "read_only": True,
                        "created_at": now,
                    }
                )
            frozen = {
                key: next(a["sha256"] for a in artifacts if a["kind"] == key)
                for key in ("structure", "mapping", "sectors", "surface_graph")
            }
            snapshot = {
                "schema_version": "1.0",
                "project_id": project_id,
                "name": name,
                "source_job_id": source_job_id,
                "source_v4_schema_version": v4_manifest["schema_version"],
                "source_software": v4_manifest.get("software", {}),
                "revision": 1,
                "reference_sequence": sequence,
                "frozen_hashes": frozen,
                "artifacts": artifacts,
                "created_at": now,
            }
            schema_path = (
                Path(__file__).resolve().parent
                / "schemas"
                / "experiment_project_manifest.schema.json"
            )
            Draft202012Validator(read_json(schema_path)).validate(snapshot)
            write_json(staging / "project_manifest.json", snapshot)
            os.replace(staging, final)
            with self._connect() as db:
                db.execute(
                    "INSERT INTO experiment_projects VALUES (?,?,?,?,?,?,?)",
                    (project_id, source_job_id, name, sequence, "COMPLETE", 1, now),
                )
                db.execute(
                    "INSERT INTO project_revisions VALUES (?,?,?,?,?)",
                    (
                        project_id,
                        1,
                        "INITIAL_SNAPSHOT",
                        json.dumps(snapshot, ensure_ascii=False),
                        now,
                    ),
                )
                db.executemany(
                    "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                    [
                        (
                            a["id"],
                            project_id,
                            a["relative_path"],
                            a["kind"],
                            a["sha256"],
                            a["size_bytes"],
                            a["media_type"],
                            1,
                            now,
                        )
                        for a in artifacts
                    ],
                )
                db.execute(
                    "INSERT INTO audit_events (project_id,entity_type,entity_id,action,actor,before_revision,after_revision,details_json,created_at) VALUES (?,?,?,?,?,?,?,?,?)",
                    (
                        project_id,
                        "project",
                        project_id,
                        "CREATE",
                        "local",
                        None,
                        1,
                        json.dumps({"source_job_id": source_job_id}),
                        now,
                    ),
                )
                db.commit()
        return self.get(project_id)

    def _row(self, project_id: str) -> dict[str, Any]:
        with self._connect() as db:
            row = db.execute(
                "SELECT * FROM experiment_projects WHERE id=?", (project_id,)
            ).fetchone()
        if row is None:
            raise KeyError(project_id)
        return dict(row)

    def get(self, project_id: str) -> dict[str, Any]:
        row = self._row(project_id)
        folder = self.projects_root / project_id
        with self._connect() as db:
            artifacts = [
                dict(a)
                for a in db.execute(
                    "SELECT * FROM artifacts WHERE project_id=? ORDER BY kind", (project_id,)
                )
            ]
            initial_row = db.execute(
                "SELECT snapshot_json FROM project_revisions WHERE project_id=? AND revision=1",
                (project_id,),
            ).fetchone()
        problems = []
        manifest_path = folder / "project_manifest.json"
        try:
            manifest = read_json(manifest_path)
            Draft202012Validator(
                read_json(
                    Path(__file__).resolve().parent
                    / "schemas"
                    / "experiment_project_manifest.schema.json"
                )
            ).validate(manifest)
            initial_manifest = (
                json.loads(initial_row["snapshot_json"]) if initial_row is not None else None
            )
            if initial_manifest is None:
                problems.append("initial project snapshot revision is missing")
            else:
                immutable_fields = (
                    "schema_version", "project_id", "name", "source_job_id",
                    "source_v4_schema_version", "source_software", "revision",
                    "reference_sequence", "frozen_hashes", "created_at",
                )
                if any(
                    manifest.get(field) != initial_manifest.get(field)
                    for field in immutable_fields
                ):
                    problems.append("manifest immutable fields differ from initial snapshot")
            if (
                manifest["project_id"] != project_id
                or manifest["name"] != row["name"]
                or manifest["source_job_id"] != row["source_job_id"]
                or manifest["reference_sequence"] != row["reference_sequence"]
            ):
                problems.append("manifest identity differs from database")
            if {a["id"] for a in manifest["artifacts"]} != {a["id"] for a in artifacts}:
                problems.append("manifest artifact set differs from database")
            for kind, digest in manifest["frozen_hashes"].items():
                registered = next((a for a in artifacts if a["kind"] == kind), None)
                if registered is None or registered["sha256"] != digest:
                    problems.append(f"frozen hash differs from artifact registry: {kind}")
            for item in artifacts:
                declared = next((a for a in manifest["artifacts"] if a["id"] == item["id"]), None)
                if declared is None:
                    if item["kind"] in manifest["frozen_hashes"]:
                        problems.append(f"frozen artifact is absent from manifest: {item['id']}")
                    continue
                if any(
                    declared.get(k) != item[k]
                    for k in (
                        "id", "relative_path", "kind", "sha256", "size_bytes",
                        "media_type", "created_at", "read_only",
                    )
                ):
                    problems.append(f"artifact metadata differs: {item['id']}")
                target = (folder / item["relative_path"]).resolve()
                if not target.is_relative_to(folder.resolve()) or not target.is_file():
                    problems.append(f"artifact missing or unsafe: {item['id']}")
                elif (
                    target.stat().st_size != item["size_bytes"]
                    or sha256_file(target) != item["sha256"]
                ):
                    problems.append(f"artifact hash differs: {item['id']}")
                elif item["kind"] == "compensation_matrix":
                    matrix_document = read_json(target)
                    item["detectors"] = matrix_document.get("detectors", [])
                    item["source_filename"] = matrix_document.get("source_filename", "")
            graph = read_json(folder / "reference/surface_graph_v4.json")
            if (
                graph.get("mapping_sha256") != manifest["frozen_hashes"]["mapping"]
                or graph.get("structure_sha256") != manifest["frozen_hashes"]["structure"]
            ):
                problems.append("graph input hashes differ")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            manifest = None
            problems.append(f"invalid manifest: {exc}")
        status = "INTEGRITY_ERROR" if problems else "COMPLETE"
        if row["status"] != status:
            with self._connect() as db:
                db.execute(
                    "UPDATE experiment_projects SET status=? WHERE id=?", (status, project_id)
                )
                db.commit()
        return {
            "project_id": project_id,
            "name": row["name"],
            "source_job_id": row["source_job_id"],
            "revision": row["revision"],
            "status": status,
            "created_at": row["created_at"],
            "frozen_hashes": manifest.get("frozen_hashes") if manifest else None,
            "artifacts": artifacts,
            "manifest": manifest,
            "integrity_errors": problems,
        }

    def list(self) -> list[dict[str, Any]]:
        with self._connect() as db:
            ids = [
                row[0]
                for row in db.execute("SELECT id FROM experiment_projects ORDER BY created_at DESC")
            ]
        return [self.get(project_id) for project_id in ids]

    def artifact_path(self, project_id: str, artifact_id: str) -> tuple[Path, dict[str, Any]]:
        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        item = next((a for a in project["artifacts"] if a["id"] == artifact_id), None)
        if item is None:
            raise KeyError(artifact_id)
        path = (self.projects_root / project_id / item["relative_path"]).resolve()
        return path, item

    def freeze_source_structure_evidence(
        self, project_id: str, expected_revision: int
    ) -> dict[str, Any]:
        """Append the source job's hash-verified residue reliability and SASA evidence."""
        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        evidence_kinds = {"residue_reliability", "sasa_comparison", "reliability_qc"}
        existing = {item["kind"] for item in project["artifacts"]} & evidence_kinds
        if existing == evidence_kinds:
            return {
                "project_revision": project["revision"],
                "already_frozen": True,
                "artifacts": [item for item in project["artifacts"] if item["kind"] in evidence_kinds],
            }
        if existing:
            raise IntegrityError("Project contains a partial structure evidence snapshot")

        source_job = self.jobs.get(project["source_job_id"])
        sources, _, _ = self._source_files(source_job)
        if not evidence_kinds.issubset(sources):
            raise SnapshotError("Source job has no complete, hash-verified structure evidence set")
        project_artifacts = {item["kind"]: item for item in project["artifacts"]}
        for kind in ("structure", "mapping", "sectors", "surface_graph"):
            if sha256_file(sources[kind]) != project["frozen_hashes"][kind]:
                raise IntegrityError(f"Source job no longer matches the frozen {kind} artifact")
        if sha256_file(sources["v4_manifest"]) != project_artifacts["v4_manifest"]["sha256"]:
            raise IntegrityError("Source v4 manifest no longer matches the frozen project")

        reliability_rows, sasa_rows, reliability_qc = self._validate_structure_evidence_files(
            sources, sources["mapping"]
        )
        artifact_specs = {
            "residue_reliability": ("reference/residue_reliability.csv", "text/csv"),
            "sasa_comparison": ("reference/sasa_comparison.csv", "text/csv"),
            "reliability_qc": ("reference/reliability_qc.json", "application/json"),
        }
        now = _now()
        folder = self.projects_root / project_id
        artifacts = []
        copied_targets: list[Path] = []
        temporary_targets: list[Path] = []
        try:
            for kind, (relative, media_type) in artifact_specs.items():
                source = sources[kind]
                target = folder / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                digest = sha256_file(source)
                if target.exists() or target.is_symlink():
                    if (
                        target.is_symlink()
                        or not target.is_file()
                        or not target.resolve().is_relative_to(folder.resolve())
                        or target.stat().st_size != source.stat().st_size
                        or sha256_file(target) != digest
                    ):
                        raise IntegrityError(
                            f"Existing structure evidence file does not match source: {kind}"
                        )
                else:
                    temporary = target.with_name(target.name + f".{uuid.uuid4().hex}.tmp")
                    temporary_targets.append(temporary)
                    shutil.copyfile(source, temporary)
                    if sha256_file(temporary) != digest:
                        raise IntegrityError(
                            f"Structure evidence copy failed hash verification: {kind}"
                        )
                    os.replace(temporary, target)
                    temporary_targets.remove(temporary)
                    copied_targets.append(target)
                target.chmod(0o444)
                artifacts.append({
                    "id": uuid.uuid4().hex,
                    "project_id": project_id,
                    "relative_path": relative,
                    "kind": kind,
                    "sha256": digest,
                    "size_bytes": target.stat().st_size,
                    "media_type": media_type,
                    "read_only": 1,
                    "created_at": now,
                })
        except Exception:
            for temporary in temporary_targets:
                temporary.unlink(missing_ok=True)
            for target in copied_targets:
                _discard_unregistered_artifact(target)
            raise

        details = {
            "source_job_id": project["source_job_id"],
            "artifact_hashes": {item["kind"]: item["sha256"] for item in artifacts},
            "residue_count": len(reliability_rows),
            "sasa_comparison_status_counts": {
                status: sum(row["status"] == status for row in sasa_rows)
                for status in sorted({row["status"] for row in sasa_rows})
            },
            "sasa_context_status_counts": {
                status: sum(
                    (row.get("native_context_status") or "LEGACY_SCOPE_UNVERIFIED") == status
                    for row in sasa_rows
                )
                for status in sorted({
                    row.get("native_context_status") or "LEGACY_SCOPE_UNVERIFIED"
                    for row in sasa_rows
                })
            },
            "reliability_predicted": reliability_qc["predicted"],
        }
        folder_manifest = folder / "project_manifest.json"
        original_manifest = folder_manifest.read_bytes()
        try:
            with self._lock, self._connect() as db:
                revision = self._append_revision(
                    db, project_id, expected_revision, "structure_evidence", uuid.uuid4().hex, details
                )
                db.executemany(
                    "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                    [(
                        item["id"], item["project_id"], item["relative_path"], item["kind"],
                        item["sha256"], item["size_bytes"], item["media_type"],
                        item["read_only"], item["created_at"],
                    ) for item in artifacts],
                )
                manifest = json.loads(original_manifest)
                manifest["artifacts"].extend({
                    key: item[key]
                    for key in (
                        "id", "kind", "relative_path", "sha256", "size_bytes",
                        "media_type", "created_at",
                    )
                } | {"read_only": True} for item in artifacts)
                write_json(folder_manifest, manifest)
                db.commit()
        except Exception:
            folder_manifest.write_bytes(original_manifest)
            for target in copied_targets:
                _discard_unregistered_artifact(target)
            raise
        return {
            "project_revision": revision,
            "already_frozen": False,
            "artifacts": [{key: item[key] for key in (
                "id", "kind", "sha256", "size_bytes", "media_type", "created_at"
            )} for item in artifacts],
            "summary": details,
        }

    def read_structure_evidence(self, project_id: str) -> dict[str, Any] | None:
        """Read evidence only from registered, integrity-checked project artifacts."""
        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        by_kind = {item["kind"]: item for item in project["artifacts"]}
        required = {"residue_reliability", "sasa_comparison", "reliability_qc"}
        present = required & set(by_kind)
        if not present:
            return None
        if present != required:
            raise IntegrityError("Project contains a partial structure evidence snapshot")
        paths = {kind: self.projects_root / project_id / by_kind[kind]["relative_path"] for kind in required}
        reliability_rows, sasa_rows, reliability_qc = self._validate_structure_evidence_files(
            paths, self.projects_root / project_id / by_kind["mapping"]["relative_path"]
        )
        reliability = [{
            "uniprot_pos": int(row["residue_id"]),
            "state": row["state"],
            "plddt": float(row["plddt"]) if row.get("plddt") else None,
            "mean_local_pae": float(row["mean_local_pae"]) if row.get("mean_local_pae") else None,
            "mean_core_pae": float(row["mean_core_pae"]) if row.get("mean_core_pae") else None,
            "quality_flags": json.loads(row["quality_flags"] or "[]"),
            "use_for_surface_mesh": row["use_for_surface_mesh"].lower() == "true",
            "use_for_3d_footprint": row["use_for_3d_footprint"].lower() == "true",
            "use_for_epitope_screening": row["use_for_epitope_screening"].lower() == "true",
        } for row in reliability_rows]
        sasa = [{
            "uniprot_pos": int(row["position"]),
            "native_rsasa": float(row["native_rsasa"]) if row.get("native_rsasa") else None,
            "normalized_context_rsasa": (
                float(row["normalized_context_rsasa"])
                if row.get("normalized_context_rsasa")
                else None
            ),
            "procaliper_rsasa": float(row["procaliper_rsasa"]) if row.get("procaliper_rsasa") else None,
            "absolute_delta": float(row["absolute_delta"]) if row.get("absolute_delta") else None,
            "native_context_absolute_delta": (
                float(row["native_context_absolute_delta"])
                if row.get("native_context_absolute_delta")
                else None
            ),
            "native_context_status": row.get("native_context_status")
            or "LEGACY_SCOPE_UNVERIFIED",
            "native_context_classification_changed": (
                row.get("native_context_classification_changed", "").lower() == "true"
                if row.get("native_context_classification_changed", "")
                else None
            ),
            "procaliper_classification_changed": (
                row.get("procaliper_classification_changed", "").lower() == "true"
                if row.get("procaliper_classification_changed", "")
                else None
            ),
            "comparison_method_status": row.get("comparison_method_status")
            or "LEGACY_SCOPE_UNVERIFIED",
            "status": row["status"],
        } for row in sasa_rows]
        return {
            "available": True,
            "provenance": {
                kind: {key: by_kind[kind][key] for key in ("id", "sha256")}
                for kind in sorted(required)
            },
            "reliability": {"predicted": reliability_qc["predicted"],
                            "classification_version": reliability_qc["classification_version"],
                            "state_counts": reliability_qc["state_counts"], "records": reliability},
            "sasa": {"records": sasa,
                     "status_counts": {
                         status: sum(row["status"] == status for row in sasa)
                         for status in sorted({row["status"] for row in sasa})
                     },
                     "context_status_counts": {
                         status: sum(row["native_context_status"] == status for row in sasa)
                         for status in sorted({row["native_context_status"] for row in sasa})
                     },
                     "method_status_counts": {
                         status: sum(row["comparison_method_status"] == status for row in sasa)
                         for status in sorted({row["comparison_method_status"] for row in sasa})
                     }},
        }

    def build_support_bundle(self, project_id: str, *, include_originals: bool = False) -> bytes:
        """Export a reproducible project state, audit trail, and verified artifacts."""
        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")

        database_records: dict[str, list[dict[str, Any]]] = {}
        with self._connect() as db:
            tables = [
                row[0]
                for row in db.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name"
                )
            ]
            for table in tables:
                if not table.replace("_", "").isalnum():
                    continue
                columns = {
                    row[1] for row in db.execute(f'PRAGMA table_info("{table}")')
                }
                if table == "experiment_projects":
                    rows = db.execute(
                        f'SELECT * FROM "{table}" WHERE id=? ORDER BY rowid', (project_id,)
                    )
                elif "project_id" in columns:
                    rows = db.execute(
                        f'SELECT * FROM "{table}" WHERE project_id=? ORDER BY rowid', (project_id,)
                    )
                elif table == "schema_migrations":
                    rows = db.execute(f'SELECT * FROM "{table}" ORDER BY version')
                else:
                    continue
                database_records[table] = [dict(row) for row in rows]

        project_root = (self.projects_root / project_id).resolve()
        omitted_originals = []
        included_artifacts = []
        artifact_files: list[tuple[str, Path]] = []
        for artifact in sorted(project["artifacts"], key=lambda item: item["relative_path"]):
            is_original = artifact["kind"] in {
                "raw_fcs", "raw_flowjo", "raw_flowjo_result_table",
            }
            if is_original and not include_originals:
                omitted_originals.append(
                    {"id": artifact["id"], "kind": artifact["kind"],
                     "relative_path": artifact["relative_path"], "sha256": artifact["sha256"]}
                )
                continue
            path = (project_root / artifact["relative_path"]).resolve()
            if not path.is_relative_to(project_root) or not path.is_file():
                raise IntegrityError(f"Artifact is missing or unsafe: {artifact['id']}")
            if path.stat().st_size != artifact["size_bytes"] or sha256_file(path) != artifact["sha256"]:
                raise IntegrityError(f"Artifact integrity check failed: {artifact['id']}")
            archive_path = f"artifacts/{artifact['relative_path']}"
            artifact_files.append((archive_path, path))
            included_artifacts.append(
                {key: artifact[key] for key in (
                    "id", "kind", "relative_path", "sha256", "size_bytes", "media_type", "created_at"
                )}
            )

        bundle_manifest = {
            "schema_version": "1.0",
            "project_id": project_id,
            "project_revision": project["revision"],
            "include_originals": include_originals,
            "included_artifacts": included_artifacts,
            "omitted_originals": omitted_originals,
        }
        documents = {
            "project_manifest.json": project["manifest"],
            "project_state.json": {
                key: project[key]
                for key in ("project_id", "name", "source_job_id", "revision", "status",
                            "created_at", "frozen_hashes", "integrity_errors")
            },
            "database_records.json": database_records,
            "support_bundle_manifest.json": bundle_manifest,
        }
        output = io.BytesIO()
        with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for name, document in sorted(documents.items()):
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                archive.writestr(
                    info,
                    json.dumps(document, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8"),
                )
            for archive_path, source in artifact_files:
                info = zipfile.ZipInfo(archive_path, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                with source.open("rb") as input_file, archive.open(info, "w") as output_file:
                    shutil.copyfileobj(input_file, output_file, length=1024 * 1024)
        return output.getvalue()

    @staticmethod
    def _append_revision(
        db: sqlite3.Connection,
        project_id: str,
        expected: int,
        entity_type: str,
        entity_id: str,
        details: dict[str, Any],
    ) -> int:
        row = db.execute(
            "SELECT revision FROM experiment_projects WHERE id=?", (project_id,)
        ).fetchone()
        if row is None:
            raise KeyError(project_id)
        if row["revision"] != expected:
            raise RevisionConflict(
                f"Expected revision {expected}; current revision is {row['revision']}"
            )
        revision = expected + 1
        now = _now()
        db.execute("UPDATE experiment_projects SET revision=? WHERE id=?", (revision, project_id))
        db.execute(
            "INSERT INTO project_revisions VALUES (?,?,?,?,?)",
            (project_id, revision, f"ADD_{entity_type.upper()}", json.dumps(details), now),
        )
        db.execute(
            """INSERT INTO audit_events
                   (project_id,entity_type,entity_id,action,actor,before_revision,after_revision,details_json,created_at)
                   VALUES (?,?,?,?,?,?,?,?,?)""",
            (
                project_id,
                entity_type,
                entity_id,
                "CREATE",
                "local",
                expected,
                revision,
                json.dumps(details, ensure_ascii=False),
                now,
            ),
        )
        return revision

    def add_antigen_construct(
        self,
        project_id: str,
        *,
        name: str,
        mutations: list[dict[str, Any]] | None,
        full_sequence: str | None,
        expected_revision: int,
    ) -> dict[str, Any]:
        result = self.preview_antigen_construct(
            project_id, mutations=mutations, full_sequence=full_sequence
        )
        construct_id = uuid.uuid4().hex
        now = _now()
        digest = sha256_bytes(result["sequence"].encode("ascii"))
        with self._lock, self._connect() as db:
            revision = self._append_revision(
                db,
                project_id,
                expected_revision,
                "antigen_construct",
                construct_id,
                {"name": name, "sequence_sha256": digest},
            )
            db.execute(
                "INSERT INTO antigen_constructs VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    construct_id,
                    project_id,
                    name,
                    "CUSTOM",
                    result["sequence"],
                    digest,
                    json.dumps(result["mutations"]),
                    result["mapping_status"],
                    json.dumps(result["risk_flags"]),
                    1,
                    revision,
                    now,
                ),
            )
            db.commit()
        return {
            "construct_id": construct_id,
            "name": name,
            "source": "CUSTOM",
            "version": 1,
            "project_revision": revision,
            "sequence_sha256": digest,
            **result,
        }

    def preview_antigen_construct(
        self,
        project_id: str,
        *,
        mutations: list[dict[str, Any]] | None,
        full_sequence: str | None,
    ) -> dict[str, Any]:
        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        folder = self.projects_root / project_id / "reference"
        position_start, mapped = reference_positions(
            folder / "antigen_sequence.fasta", folder / "residue_mapping.csv"
        )
        result = validate_antigen_construct(
            project["manifest"]["reference_sequence"],
            position_start,
            mapped,
            mutations,
            full_sequence,
        )
        return result | {
            "sequence_sha256": sha256_bytes(result["sequence"].encode("ascii")),
            "reference_sequence_sha256": sha256_bytes(
                project["manifest"]["reference_sequence"].encode("ascii")
            ),
        }

    def list_recommended_antigens(self, project_id: str) -> list[dict[str, Any]]:
        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        path = self.projects_root / project_id / "constructs" / "recommended_antigens.json"
        if not path.is_file():
            # Projects created before recommendation snapshots were introduced
            # remain readable, but must not resolve recommendations from a live job.
            return []
        return read_json(path)

    def add_recommended_antigen_construct(
        self, project_id: str, *, recommended_construct_id: str, name: str | None,
        expected_revision: int,
    ) -> dict[str, Any]:
        recommendations = self.list_recommended_antigens(project_id)
        design = next((item for item in recommendations if item["id"] == recommended_construct_id), None)
        if design is None:
            raise KeyError(recommended_construct_id)
        if design.get("status") not in {
            "RECOMMENDED_LOW_PREDICTED_RISK", "RECOMMENDED_MULTI_SITE_PROBE",
            "RECOMMENDED_SINGLE_SITE_PROBE",
        }:
            raise ConstructError("Selected v4 construct is not marked recommended")
        mutations = [
            {"position": item["position"], "wt": item["wt"], "mut": item.get("mutant", item.get("mut"))}
            for item in design["mutations"]
        ]
        project = self.get(project_id)
        folder = self.projects_root / project_id / "reference"
        position_start, mapped = reference_positions(folder / "antigen_sequence.fasta",
                                                     folder / "residue_mapping.csv")
        result = validate_antigen_construct(
            project["manifest"]["reference_sequence"], position_start, mapped, mutations, None
        )
        expected_sequence = "".join(str(design["mutant_ecd_sequence"]).split()).upper()
        if result["sequence"] != expected_sequence:
            raise ConstructError("Recommended mutation list does not match its frozen output sequence")
        flags = list(result["risk_flags"])
        flags.extend(str(item) for item in design.get("warnings", []))
        flags.append(str(design["status"]))
        construct_id = uuid.uuid4().hex
        construct_name = name or str(design.get("name") or design["id"])
        now = _now()
        digest = sha256_bytes(result["sequence"].encode("ascii"))
        with self._lock, self._connect() as db:
            revision = self._append_revision(
                db, project_id, expected_revision, "antigen_construct", construct_id,
                {"name": construct_name, "source": "RECOMMENDED", "candidate_id": design["id"],
                 "sequence_sha256": digest},
            )
            db.execute(
                "INSERT INTO antigen_constructs VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                (construct_id, project_id, construct_name, "RECOMMENDED", result["sequence"], digest,
                 json.dumps(result["mutations"]), result["mapping_status"], json.dumps(flags), 1,
                 revision, now),
            )
            db.commit()
        return {"construct_id": construct_id, "name": construct_name, "source": "RECOMMENDED",
                "recommended_construct_id": design["id"], "version": 1,
                "project_revision": revision, "sequence_sha256": digest,
                **result, "risk_flags": flags}

    def add_antibody_construct(
        self,
        project_id: str,
        *,
        name: str,
        role: str,
        vh_sequence: str,
        vl_sequence: str,
        parent_construct_id: str | None,
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        with self._lock, self._connect() as db:
            parent = None
            if parent_construct_id:
                parent = db.execute(
                    "SELECT * FROM antibody_constructs WHERE id=? AND project_id=?",
                    (parent_construct_id, project_id),
                ).fetchone()
                if parent is None:
                    raise KeyError(parent_construct_id)
            if (role == "PARENT") != (parent is None):
                raise ConstructError("PARENT must have no parent; VARIANT requires a parent")
            result = validate_antibody_construct(
                vh_sequence,
                vl_sequence,
                parent["vh_sequence"] if parent else None,
                parent["vl_sequence"] if parent else None,
            )
            construct_id = uuid.uuid4().hex
            revision = self._append_revision(
                db,
                project_id,
                expected_revision,
                "antibody_construct",
                construct_id,
                {"name": name, "role": role, "parent_id": parent_construct_id},
            )
            now = _now()
            db.execute(
                "INSERT INTO antibody_constructs VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    construct_id,
                    project_id,
                    name,
                    parent_construct_id,
                    role,
                    result["vh_sequence"],
                    result["vl_sequence"],
                    json.dumps(result["heavy_mutations"]),
                    json.dumps(result["light_mutations"]),
                    result["mutation_class"],
                    1,
                    revision,
                    now,
                ),
            )
            db.commit()
        return {
            "construct_id": construct_id,
            "name": name,
            "parent_id": parent_construct_id,
            "role": role,
            "version": 1,
            "project_revision": revision,
            "sequence_hashes": {
                "vh": sha256_bytes(result["vh_sequence"].encode("ascii")),
                "vl": sha256_bytes(result["vl_sequence"].encode("ascii")),
            },
            **result,
        }

    def add_batch(
        self,
        project_id: str,
        *,
        name: str,
        date: str,
        assay: str,
        cell_line: str,
        dose: float | None,
        dose_unit: str | None,
        positive_gate_id: str | None,
        independent_batch: bool,
        settings: dict[str, Any],
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        batch_id = uuid.uuid4().hex
        with self._lock, self._connect() as db:
            # Serialize cross-process writers before they read the revision. A
            # deferred transaction can otherwise fail while upgrading its read lock.
            db.execute("BEGIN IMMEDIATE")
            revision = self._append_revision(
                db,
                project_id,
                expected_revision,
                "assay_batch",
                batch_id,
                {"name": name, "date": date, "assay": assay},
            )
            db.execute(
                "INSERT INTO assay_batches VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    batch_id,
                    project_id,
                    name,
                    date,
                    assay,
                    cell_line,
                    dose,
                    dose_unit,
                    positive_gate_id,
                    int(independent_batch),
                    json.dumps(settings),
                    revision,
                    _now(),
                ),
            )
            db.commit()
        return {
            "batch_id": batch_id,
            "project_revision": revision,
            "name": name,
            "date": date,
            "assay": assay,
            "cell_line": cell_line,
            "dose": dose,
            "dose_unit": dose_unit,
            "positive_gate_id": positive_gate_id,
            "independent_batch": independent_batch,
            "settings": settings,
        }

    def list_entities(self, project_id: str, kind: str) -> list[dict[str, Any]]:
        project = self._row(project_id)
        if kind not in {"antigen_constructs", "antibody_constructs", "assay_batches", "samples"}:
            raise ValueError("Unknown entity table")
        with self._connect() as db:
            query = (
                "SELECT s.*,p.metadata_json AS source_metadata_json FROM samples s "
                "LEFT JOIN pending_uploads p ON p.id=s.source_upload_id "
                "WHERE s.project_id=? ORDER BY s.created_at"
                if kind == "samples"
                else f"SELECT * FROM {kind} WHERE project_id=? ORDER BY created_at"
            )
            rows = [dict(row) for row in db.execute(query, (project_id,))]
            if kind == "assay_batches":
                sample_rows = [
                    dict(row)
                    for row in db.execute(
                        "SELECT batch_id,role,antigen_construct_id,antibody_construct_id,condition_json "
                        "FROM samples WHERE project_id=?",
                        (project_id,),
                    )
                ]
                wt_construct_ids = {
                    row["id"]
                    for row in db.execute(
                        "SELECT id FROM antigen_constructs WHERE project_id=? AND sequence=?",
                        (project_id, project["reference_sequence"]),
                    )
                }
        for row in rows:
            source_metadata = row.pop("source_metadata_json", None)
            if source_metadata is not None:
                row["fcs_metadata"] = json.loads(source_metadata)
            for key in list(row):
                if key.endswith("_json"):
                    row[key.removesuffix("_json")] = json.loads(row.pop(key))
            if kind == "assay_batches":
                warnings = []
                if row["assay"] == "ACTIVATION":
                    batch_samples = [
                        sample for sample in sample_rows if sample["batch_id"] == row["id"]
                    ]
                    activation_samples = [
                        sample for sample in batch_samples if sample["role"] == "ACTIVATION"
                    ]
                    temporary_identity_by_sample = {
                        id(sample): json.loads(sample["condition_json"]).get("temporary_identity")
                        for sample in batch_samples
                    }
                    temporary_wt_present = any(
                        temporary_identity_by_sample[id(sample)]
                        and temporary_identity_by_sample[id(sample)].get("antigen_role") == "REFERENCE"
                        for sample in activation_samples
                    )
                    if any(temporary_identity_by_sample.values()):
                        warnings.append("TEMPORARY_IDENTITIES_UNVERIFIED")
                    if not activation_samples:
                        warnings.append("ACTIVATION_SAMPLE_MISSING")
                    if not any(sample["role"] == "NC" for sample in batch_samples):
                        warnings.append("NC_CONTROL_MISSING")
                    if not wt_construct_ids:
                        warnings.append(
                            "WT_IDENTITY_UNVERIFIED" if temporary_wt_present else "WT_CONSTRUCT_MISSING"
                        )
                    elif not any(
                        sample["antigen_construct_id"] in wt_construct_ids
                        for sample in activation_samples
                    ):
                        warnings.append(
                            "WT_IDENTITY_UNVERIFIED" if temporary_wt_present else "WT_REFERENCE_SAMPLE_MISSING"
                        )
                    missing_links_without_labels = any(
                        (
                            not sample["antigen_construct_id"]
                            or not sample["antibody_construct_id"]
                        )
                        and not temporary_identity_by_sample[id(sample)]
                        for sample in batch_samples
                    )
                    if missing_links_without_labels:
                        warnings.append("CONSTRUCT_LINK_MISSING")
                row["design_warnings"] = warnings
                row["design_status"] = "READY" if not warnings else "INCOMPLETE"
        return rows

    def add_summary_measurement(
        self,
        project_id: str,
        *,
        assay: str,
        construct_id: str,
        measured_on: str,
        metric: str,
        value: float,
        unit: str,
        source_label: str,
        notes: str,
        reviewer: str,
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        if assay not in {"ANTIGEN_MYC", "FAB_MYC"}:
            raise ValueError("Summary-only measurements only support ANTIGEN_MYC or FAB_MYC")
        if metric not in {"POSITIVE_FRACTION", "MEDIAN_FLUORESCENCE", "GEOMETRIC_MEAN", "OTHER"}:
            raise ValueError("Unsupported summary-only metric")
        if not math.isfinite(value) or value < 0:
            raise ValueError("Summary-only value must be finite and non-negative")
        if metric == "POSITIVE_FRACTION" and value > 1:
            raise ValueError("POSITIVE_FRACTION must be between 0 and 1")
        try:
            datetime.strptime(measured_on, "%Y-%m-%d")
        except ValueError as exc:
            raise ValueError("measured_on must use YYYY-MM-DD") from exc
        if not source_label.strip() or not reviewer.strip():
            raise ValueError("Summary-only source label and reviewer are required")
        measurement_id = uuid.uuid4().hex
        now = _now()
        antigen_id = construct_id if assay == "ANTIGEN_MYC" else None
        antibody_id = construct_id if assay == "FAB_MYC" else None
        table = "antigen_constructs" if antigen_id else "antibody_constructs"
        with self._lock, self._connect() as db:
            if not db.execute(
                f"SELECT 1 FROM {table} WHERE id=? AND project_id=?",
                (construct_id, project_id),
            ).fetchone():
                raise KeyError(construct_id)
            revision = self._append_revision(
                db,
                project_id,
                expected_revision,
                "summary_measurement",
                measurement_id,
                {
                    "assay": assay,
                    "construct_id": construct_id,
                    "measured_on": measured_on,
                    "metric": metric,
                    "value": value,
                    "unit": unit.strip(),
                    "source_label": source_label.strip(),
                    "reviewer": reviewer.strip(),
                    "interpretation": "HISTORICAL_SUMMARY_ONLY",
                },
            )
            db.execute(
                "INSERT INTO summary_measurements VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    measurement_id,
                    project_id,
                    assay,
                    antigen_id,
                    antibody_id,
                    measured_on,
                    metric,
                    value,
                    unit.strip(),
                    source_label.strip(),
                    notes.strip(),
                    reviewer.strip(),
                    1,
                    revision,
                    now,
                ),
            )
            db.commit()
        return {
            "measurement_id": measurement_id,
            "project_revision": revision,
            "assay": assay,
            "construct_id": construct_id,
            "measured_on": measured_on,
            "metric": metric,
            "value": value,
            "unit": unit.strip(),
            "source_label": source_label.strip(),
            "notes": notes.strip(),
            "reviewer": reviewer.strip(),
            "historical": True,
            "interpretation": "HISTORICAL_SUMMARY_ONLY_NOT_CONCURRENT_CONTROL",
            "created_at": now,
        }

    def list_summary_measurements(self, project_id: str) -> list[dict[str, Any]]:
        self._row(project_id)
        with self._connect() as db:
            rows = db.execute(
                "SELECT * FROM summary_measurements WHERE project_id=? ORDER BY measured_on,created_at,id",
                (project_id,),
            ).fetchall()
        return [
            dict(row)
            | {
                "construct_id": row["antigen_construct_id"] or row["antibody_construct_id"],
                "historical": bool(row["historical"]),
                "interpretation": "HISTORICAL_SUMMARY_ONLY_NOT_CONCURRENT_CONTROL",
            }
            for row in rows
        ]

    def import_fcs_files(
        self,
        project_id: str,
        files: list[dict[str, Any]],
        *,
        batch_id: str | None,
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        session_id = uuid.uuid4().hex
        now = _now()
        created_artifacts: list[dict[str, Any]] = []
        prepared: list[dict[str, Any]] = []
        seen_hashes: dict[str, dict[str, Any]] = {}
        project_folder = self.projects_root / project_id
        for item in files:
            source = Path(item["path"]).resolve()
            if not source.is_file():
                raise ValueError("Uploaded temporary file is missing")
            digest = sha256_file(source)
            if digest in seen_hashes:
                prepared.append(
                    {
                        "artifact": seen_hashes[digest],
                        "original_filename": item["original_filename"],
                        "metadata": item["metadata"],
                        "sha256": digest,
                    }
                )
                continue
            relative = f"raw/fcs/{digest}.fcs"
            target = project_folder / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            with self._connect() as db:
                existing = db.execute(
                    "SELECT * FROM artifacts WHERE project_id=? AND kind='raw_fcs' AND sha256=?",
                    (project_id, digest),
                ).fetchone()
            def target_matches_uploaded_file() -> bool:
                return (
                    not target.is_symlink()
                    and target.is_file()
                    and target.stat().st_size == source.stat().st_size
                    and sha256_file(target) == digest
                )

            if existing:
                artifact = dict(existing)
                if not target_matches_uploaded_file():
                    raise IntegrityError("Previously registered FCS content is missing or changed")
            else:
                if target.exists() or target.is_symlink():
                    if not target_matches_uploaded_file():
                        raise IntegrityError(
                            "Unregistered FCS path exists and does not match uploaded content"
                        )
                else:
                    temp = target.with_name(target.name + f".{uuid.uuid4().hex}.tmp")
                    try:
                        shutil.copyfile(source, temp)
                        if sha256_file(temp) != digest:
                            raise IntegrityError("FCS copy failed hash verification")
                        try:
                            os.replace(temp, target)
                        except OSError:
                            # A concurrent import may have published the same immutable bytes.
                            if not target_matches_uploaded_file():
                                raise
                        if not target_matches_uploaded_file():
                            raise IntegrityError("Published FCS content failed hash verification")
                    finally:
                        temp.unlink(missing_ok=True)
                target.chmod(0o444)
                artifact = {
                    "id": uuid.uuid4().hex,
                    "project_id": project_id,
                    "relative_path": relative,
                    "kind": "raw_fcs",
                    "sha256": digest,
                    "size_bytes": target.stat().st_size,
                    "media_type": "application/octet-stream",
                    "read_only": 1,
                    "created_at": now,
                }
                created_artifacts.append(artifact)
            seen_hashes[digest] = artifact
            prepared.append(
                {
                    "artifact": artifact,
                    "original_filename": item["original_filename"],
                    "metadata": item["metadata"],
                    "sha256": digest,
                }
            )
        with self._lock, self._connect() as db:
            if (
                batch_id
                and not db.execute(
                    "SELECT 1 FROM assay_batches WHERE id=? AND project_id=?",
                    (batch_id, project_id),
                ).fetchone()
            ):
                raise KeyError(batch_id)
            revision = self._append_revision(
                db,
                project_id,
                expected_revision,
                "fcs_import",
                session_id,
                {"files": len(prepared), "batch_id": batch_id},
            )
            db.execute(
                "INSERT INTO import_sessions VALUES (?,?,?,?,?,?)",
                (session_id, project_id, batch_id, "UNASSIGNED", revision, now),
            )
            for artifact in created_artifacts:
                db.execute(
                    "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                    (
                        artifact["id"],
                        artifact["project_id"],
                        artifact["relative_path"],
                        artifact["kind"],
                        artifact["sha256"],
                        artifact["size_bytes"],
                        artifact["media_type"],
                        artifact["read_only"],
                        artifact["created_at"],
                    ),
                )
            manifest_path = project_folder / "project_manifest.json"
            snapshot = read_json(manifest_path)
            registered_ids = {item["id"] for item in snapshot["artifacts"]}
            snapshot["artifacts"].extend(
                {
                    **{
                        key: artifact[key]
                        for key in (
                            "id",
                            "kind",
                            "relative_path",
                            "sha256",
                            "size_bytes",
                            "media_type",
                            "created_at",
                        )
                    },
                    "read_only": bool(artifact["read_only"]),
                }
                for artifact in created_artifacts
                if artifact["id"] not in registered_ids
            )
            write_json(manifest_path, snapshot)
            uploads = []
            for item in prepared:
                upload_id = uuid.uuid4().hex
                artifact = item["artifact"]
                if "project_id" not in artifact:
                    artifact = db.execute(
                        "SELECT * FROM artifacts WHERE id=?", (artifact["id"],)
                    ).fetchone()
                db.execute(
                    "INSERT INTO pending_uploads VALUES (?,?,?,?,?,?,?,?,?)",
                    (
                        upload_id,
                        session_id,
                        project_id,
                        artifact["id"],
                        item["original_filename"],
                        item["sha256"],
                        json.dumps(item["metadata"], ensure_ascii=False),
                        "UNASSIGNED",
                        now,
                    ),
                )
                uploads.append(
                    {
                        "upload_id": upload_id,
                        "artifact_id": artifact["id"],
                        "original_filename": item["original_filename"],
                        "sha256": item["sha256"],
                        "event_count": item["metadata"]["event_count"],
                        "channels": item["metadata"]["channels"],
                        "metadata": item["metadata"]["metadata"],
                        "compensation_source": item["metadata"]["compensation_source"],
                        "match_status": "UNASSIGNED",
                    }
                )
            db.commit()
        return {"upload_session_id": session_id, "project_revision": revision, "files": uploads}

    def confirm_fcs_matches(
        self,
        project_id: str,
        session_id: str,
        matches: list[dict[str, Any]],
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        with self._lock, self._connect() as db:
            session = db.execute(
                "SELECT * FROM import_sessions WHERE id=? AND project_id=?",
                (session_id, project_id),
            ).fetchone()
            if session is None:
                raise KeyError(session_id)
            uploads = {
                row["id"]: row
                for row in db.execute(
                    "SELECT * FROM pending_uploads WHERE session_id=? AND project_id=?",
                    (session_id, project_id),
                )
            }
            samples = []
            matched_upload_ids = [match["upload_id"] for match in matches]
            if len(matched_upload_ids) != len(set(matched_upload_ids)):
                raise ValueError("An upload can only be matched once")
            for match in matches:
                upload = uploads.get(match["upload_id"])
                if upload is None:
                    raise ValueError("upload_id is not in this session")
                if match["role"] not in {
                    "ACTIVATION",
                    "NC",
                    "ANTIGEN_MYC",
                    "FAB_MYC",
                    "OTHER_CONTROL",
                }:
                    raise ValueError("Unknown sample role")
                batch_id = match.get("batch_id") or session["batch_id"]
                if (
                    not batch_id
                    or not db.execute(
                        "SELECT 1 FROM assay_batches WHERE id=? AND project_id=?",
                        (batch_id, project_id),
                    ).fetchone()
                ):
                    raise ValueError("A valid batch is required for every matched sample")
                for field, table in (
                    ("antigen_construct_id", "antigen_constructs"),
                    ("antibody_construct_id", "antibody_constructs"),
                ):
                    value = match.get(field)
                    if (
                        value
                        and not db.execute(
                            f"SELECT 1 FROM {table} WHERE id=? AND project_id=?",
                            (value, project_id),
                        ).fetchone()
                    ):
                        raise KeyError(value)
                temporary_identity = match.get("temporary_identity")
                if temporary_identity is not None:
                    if (
                        temporary_identity.get("status") != "TEMPORARY_LABEL"
                        or temporary_identity.get("source") != "USER_CONFIRMED_LABELS"
                        or not str(temporary_identity.get("panel_label", "")).strip()
                        or temporary_identity.get("antigen_role")
                        not in {"REFERENCE", "VARIANT", "NONE"}
                        or temporary_identity.get("panel_role")
                        not in {"ANTIBODY", "POSITIVE_LIGAND_CONTROL"}
                        or (
                            temporary_identity.get("antigen_role") == "NONE"
                            and temporary_identity.get("antigen_label") is not None
                        )
                        or (
                            temporary_identity.get("antigen_role") != "NONE"
                            and not str(temporary_identity.get("antigen_label", "")).strip()
                        )
                    ):
                        raise ValueError("Temporary identity labels are incomplete or inconsistent")
                if match["role"] == "ACTIVATION":
                    has_construct_links = bool(
                        match.get("antigen_construct_id")
                        and match.get("antibody_construct_id")
                    )
                    has_temporary_links = bool(
                        temporary_identity
                        and temporary_identity.get("antigen_role") in {"REFERENCE", "VARIANT"}
                        and temporary_identity.get("panel_role") == "ANTIBODY"
                    )
                    if not (has_construct_links or has_temporary_links):
                        raise ValueError(
                            "ACTIVATION samples require antigen and antibody constructs or explicit temporary labels"
                        )
                sample_id = uuid.uuid4().hex
                condition = dict(match.get("condition", {}))
                if temporary_identity is not None:
                    condition["temporary_identity"] = temporary_identity
                samples.append(
                    (
                        sample_id,
                        project_id,
                        batch_id,
                        upload["artifact_id"],
                        match["role"],
                        match.get("antigen_construct_id"),
                        match.get("antibody_construct_id"),
                        match["replicate_type"],
                        match["replicate_index"],
                        json.dumps(condition, ensure_ascii=False),
                        upload["id"],
                    )
                )
            revision = expected_revision
            if samples:
                revision = self._append_revision(
                    db,
                    project_id,
                    expected_revision,
                    "sample_import",
                    session_id,
                    {"sample_count": len(samples)},
                )
                db.executemany(
                    "INSERT INTO samples VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    [(*row, _now()) for row in samples],
                )
                for row in samples:
                    db.execute(
                        "UPDATE pending_uploads SET match_status='MATCHED' WHERE id=?", (row[-1],)
                    )
            else:
                current = db.execute(
                    "SELECT revision FROM experiment_projects WHERE id=?", (project_id,)
                ).fetchone()
                if current["revision"] != expected_revision:
                    raise RevisionConflict(
                        f"Expected revision {expected_revision}; current revision is {current['revision']}"
                    )
            unmatched = [
                upload_id
                for upload_id, row in uploads.items()
                if row["match_status"] != "MATCHED"
                and upload_id not in {m["upload_id"] for m in matches}
            ]
            status = "CONFIRMED" if not unmatched else ("PARTIAL" if samples else "UNASSIGNED")
            db.execute(
                "UPDATE import_sessions SET status=?,project_revision=? WHERE id=?",
                (status, revision, session_id),
            )
            db.commit()
        return {
            "upload_session_id": session_id,
            "sample_ids": [row[0] for row in samples],
            "unmatched_upload_ids": unmatched,
            "project_revision": revision,
            "status": status,
        }

    def get_import_session(self, project_id: str, session_id: str) -> dict[str, Any]:
        self._row(project_id)
        with self._connect() as db:
            session = db.execute(
                "SELECT * FROM import_sessions WHERE id=? AND project_id=?",
                (session_id, project_id),
            ).fetchone()
            if session is None:
                raise KeyError(session_id)
            rows = db.execute(
                "SELECT p.*,a.size_bytes,a.media_type FROM pending_uploads p "
                "JOIN artifacts a ON a.id=p.artifact_id WHERE p.session_id=? ORDER BY p.created_at,p.id",
                (session_id,),
            ).fetchall()
        files = []
        for row in rows:
            item = dict(row)
            item["metadata"] = json.loads(item.pop("metadata_json"))
            files.append(item)
        return {
            "upload_session_id": session_id,
            "status": session["status"],
            "project_revision": session["project_revision"],
            "batch_id": session["batch_id"],
            "files": files,
        }

    def register_flowjo_import(
        self,
        project_id: str,
        workspace_path: Path,
        original_filename: str,
        preview: dict[str, Any],
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        records = preview.get("records") or {}
        unresolved_nodes = records.get("unresolved_nodes") or []
        outer_status = preview.get("status", "COMPLETE")
        records_status = records.get("status", "COMPLETE")
        if not isinstance(unresolved_nodes, list):
            raise ValueError("FlowJo unresolved_nodes must be a list")
        if outer_status not in {"COMPLETE", "WAITING_FOR_USER"}:
            raise ValueError("FlowJo import preview has an invalid outer status")
        if records_status not in {"COMPLETE", "WAITING_FOR_USER"}:
            raise ValueError("FlowJo import preview has an invalid records status")
        import_status = (
            "WAITING_FOR_USER"
            if unresolved_nodes or outer_status == "WAITING_FOR_USER" or records_status == "WAITING_FOR_USER"
            else "COMPLETE"
        )
        if import_status == "WAITING_FOR_USER" and not unresolved_nodes:
            raise ValueError("FlowJo import is waiting but has no unresolved nodes to review")
        import_id = uuid.uuid4().hex
        now = _now()
        workspace_hash = sha256_file(workspace_path)
        folder = self.projects_root / project_id
        workspace_relative = f"raw/flowjo/{workspace_hash}.wsp"
        workspace_target = folder / workspace_relative
        workspace_target.parent.mkdir(parents=True, exist_ok=True)
        with self._connect() as db:
            existing = db.execute(
                "SELECT * FROM artifacts WHERE project_id=? AND kind='raw_flowjo' AND sha256=?",
                (project_id, workspace_hash),
            ).fetchone()
        new_artifacts: list[dict[str, Any]] = []
        if existing:
            workspace_artifact = dict(existing)
            if not workspace_target.is_file() or sha256_file(workspace_target) != workspace_hash:
                raise IntegrityError("Previously registered FlowJo workspace is missing or changed")
        else:
            if workspace_target.is_symlink():
                raise IntegrityError("Unregistered FlowJo workspace path is a symlink")
            if workspace_target.exists():
                if (
                    not workspace_target.is_file()
                    or not workspace_target.resolve().is_relative_to(folder.resolve())
                    or sha256_file(workspace_target) != workspace_hash
                ):
                    raise IntegrityError(
                        "Unregistered FlowJo workspace path exists with different content"
                    )
            else:
                temp = workspace_target.with_name(
                    workspace_target.name + f".{uuid.uuid4().hex}.tmp"
                )
                try:
                    shutil.copyfile(workspace_path, temp)
                    if sha256_file(temp) != workspace_hash:
                        raise IntegrityError("FlowJo workspace copy failed hash verification")
                    try:
                        os.replace(temp, workspace_target)
                    except OSError:
                        # A concurrent import or a retry after a hard crash may
                        # have already published this content-addressed file.
                        if (
                            workspace_target.is_symlink()
                            or not workspace_target.is_file()
                            or sha256_file(workspace_target) != workspace_hash
                        ):
                            raise
                finally:
                    temp.unlink(missing_ok=True)
            workspace_target.chmod(0o444)
            workspace_artifact = {
                "id": uuid.uuid4().hex,
                "project_id": project_id,
                "relative_path": workspace_relative,
                "kind": "raw_flowjo",
                "sha256": workspace_hash,
                "size_bytes": workspace_target.stat().st_size,
                "media_type": "application/xml",
                "read_only": 1,
                "created_at": now,
            }
            new_artifacts.append(workspace_artifact)
        preview_relative = f"analyses/flowjo_imports/{import_id}/preview.json"
        preview_target = folder / preview_relative
        preview_target.parent.mkdir(parents=True, exist_ok=True)
        write_json(preview_target, preview)
        preview_hash = sha256_file(preview_target)
        preview_artifact = {
            "id": uuid.uuid4().hex,
            "project_id": project_id,
            "relative_path": preview_relative,
            "kind": "flowjo_import_preview",
            "sha256": preview_hash,
            "size_bytes": preview_target.stat().st_size,
            "media_type": "application/json",
            "read_only": 1,
            "created_at": now,
        }
        preview_target.chmod(0o444)
        new_artifacts.append(preview_artifact)
        with self._lock, self._connect() as db:
            revision = self._append_revision(
                db,
                project_id,
                expected_revision,
                "flowjo_import",
                import_id,
                {
                    "source": original_filename,
                    "workspace_sha256": workspace_hash,
                    "unresolved_nodes": len(unresolved_nodes),
                },
            )
            for artifact in new_artifacts:
                db.execute(
                    "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                    (
                        artifact["id"],
                        artifact["project_id"],
                        artifact["relative_path"],
                        artifact["kind"],
                        artifact["sha256"],
                        artifact["size_bytes"],
                        artifact["media_type"],
                        artifact["read_only"],
                        artifact["created_at"],
                    ),
                )
            manifest_path = folder / "project_manifest.json"
            snapshot = read_json(manifest_path)
            registered_ids = {item["id"] for item in snapshot["artifacts"]}
            snapshot["artifacts"].extend(
                {
                    **{
                        key: artifact[key]
                        for key in (
                            "id",
                            "kind",
                            "relative_path",
                            "sha256",
                            "size_bytes",
                            "media_type",
                            "created_at",
                        )
                    },
                    "read_only": bool(artifact["read_only"]),
                }
                for artifact in new_artifacts
                if artifact["id"] not in registered_ids
            )
            write_json(manifest_path, snapshot)
            db.execute(
                "INSERT INTO flowjo_imports VALUES (?,?,?,?,?,?,?,?)",
                (
                    import_id,
                    project_id,
                    workspace_artifact["id"],
                    preview_artifact["id"],
                    import_status,
                    json.dumps(unresolved_nodes, ensure_ascii=False),
                    revision,
                    now,
                ),
            )
            db.commit()
        return {
            "flowjo_import_id": import_id,
            "status": import_status,
            "project_revision": revision,
            "workspace_artifact_id": workspace_artifact["id"],
            "preview_artifact_id": preview_artifact["id"],
            "preview": preview,
        }

    def register_flowjo_percentage_import(
        self,
        project_id: str,
        *,
        batch_id: str,
        source_path: Path,
        parsed_table: dict[str, Any],
        compensation_status: str = "UNKNOWN",
        compensation_note: str = "",
        expected_revision: int,
    ) -> dict[str, Any]:
        """Register a read-only, per-well percentage export and its normalized view."""
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        if not source_path.is_file() or source_path.is_symlink():
            raise ValueError("FlowJo result table upload is missing or unsafe")
        observations = parsed_table.get("observations")
        if not isinstance(observations, list) or not observations:
            raise ValueError("FlowJo result table has no per-well observations")
        if compensation_status not in {"NOT_APPLIED", "APPLIED", "UNKNOWN"}:
            raise ValueError("FlowJo table compensation status must be NOT_APPLIED, APPLIED or UNKNOWN")
        if len(compensation_note) > 500:
            raise ValueError("FlowJo table compensation note must be 500 characters or fewer")
        import_id = uuid.uuid4().hex
        now = _now()
        folder = self.projects_root / project_id
        source_hash = sha256_file(source_path)
        suffix = source_path.suffix.lower()
        source_relative = f"raw/flowjo_result_tables/{source_hash}{suffix}"
        source_target = folder / source_relative
        result_relative = f"analyses/flowjo_percentage_imports/{import_id}/results.json"
        result_target = folder / result_relative
        source_created = False
        result_created = False
        source_temp: Path | None = None
        result_temp: Path | None = None
        with self._lock, self._connect() as db:
            batch = db.execute(
                "SELECT id,name,assay FROM assay_batches WHERE id=? AND project_id=?",
                (batch_id, project_id),
            ).fetchone()
            if batch is None:
                raise KeyError(batch_id)
            if batch["assay"] != "ACTIVATION":
                raise ValueError("FlowJo activation percentage tables require an ACTIVATION batch")
            sample_rows = db.execute(
                "SELECT * FROM samples WHERE project_id=? AND batch_id=? ORDER BY replicate_index,id",
                (project_id, batch_id),
            ).fetchall()
            samples_by_well: dict[str, sqlite3.Row] = {}
            for sample in sample_rows:
                condition = json.loads(sample["condition_json"])
                well_id = str(condition.get("well_id", "")).strip().upper()
                if not well_id or well_id in samples_by_well:
                    raise IntegrityError("Batch samples do not have unique registered well IDs")
                samples_by_well[well_id] = sample
            supplied_wells = [str(item.get("well_id", "")).strip().upper() for item in observations]
            if len(supplied_wells) != len(set(supplied_wells)):
                raise ValueError("A well may appear only once in a FlowJo result table")
            unknown_wells = sorted(set(supplied_wells) - set(samples_by_well))
            if unknown_wells:
                raise ValueError("Table wells are not registered in the selected batch: " + ", ".join(unknown_wells))
            if not set(supplied_wells):
                raise ValueError("FlowJo result table has no rows matching the selected batch")
            mapped_observations = []
            for item in observations:
                well_id = str(item["well_id"]).strip().upper()
                fraction = item.get("fraction")
                if not isinstance(fraction, (int, float)) or not math.isfinite(fraction) or not 0 <= fraction <= 1:
                    raise ValueError(f"Invalid percentage value for well {well_id}")
                sample = samples_by_well[well_id]
                condition = json.loads(sample["condition_json"])
                temporary_identity = condition.get("temporary_identity") or {}
                mapped_observations.append({
                    "sample_id": sample["id"],
                    "sample_artifact_id": sample["artifact_id"],
                    "batch_id": batch_id,
                    "well_id": well_id,
                    "sample_role": sample["role"],
                    "condition_label": condition.get("condition_label"),
                    "panel_label": condition.get("panel_label") or temporary_identity.get("panel_label"),
                    "panel_role": condition.get("panel_role") or temporary_identity.get("panel_role"),
                    "antigen_label": temporary_identity.get("antigen_label"),
                    "antigen_role": temporary_identity.get("antigen_role"),
                    "replicate_type": sample["replicate_type"],
                    "replicate_index": sample["replicate_index"],
                    "fraction": float(fraction),
                    "percentage": float(fraction) * 100,
                    "positive_count": None,
                    "parent_count": None,
                    "displayed_value": item.get("raw_value"),
                    "display_precision_decimal_places": item.get("display_precision_decimal_places"),
                    "source_row": item.get("source_row"),
                    "qc_status": "REVIEW_REQUIRED",
                })

            linked_flowjo_import = self._matching_flowjo_import_for_batch(
                db, project_id, {sample["artifact_id"] for sample in sample_rows}
            )
            source_artifact_row = db.execute(
                "SELECT * FROM artifacts WHERE project_id=? AND kind='raw_flowjo_result_table' AND sha256=?",
                (project_id, source_hash),
            ).fetchone()
            if source_artifact_row:
                source_artifact = dict(source_artifact_row)
                source_target = folder / source_artifact["relative_path"]
                if (
                    source_target.is_symlink()
                    or not source_target.is_file()
                    or source_target.stat().st_size != source_artifact["size_bytes"]
                    or sha256_file(source_target) != source_hash
                ):
                    raise IntegrityError("Previously registered FlowJo result table is missing or changed")
            else:
                source_target.parent.mkdir(parents=True, exist_ok=True)
                if source_target.is_symlink():
                    raise IntegrityError("FlowJo result table target is unsafe")
                if source_target.exists():
                    # A prior process may have published the content-addressed source
                    # and exited before SQLite registered it. Reuse only exact bytes.
                    if (
                        not source_target.is_file()
                        or source_target.stat().st_size != source_path.stat().st_size
                        or sha256_file(source_target) != source_hash
                    ):
                        raise IntegrityError(
                            "Unregistered FlowJo result table target differs from the uploaded file"
                        )
                else:
                    source_temp = source_target.with_name(source_target.name + f".{uuid.uuid4().hex}.tmp")
                    shutil.copyfile(source_path, source_temp)
                    if sha256_file(source_temp) != source_hash:
                        raise IntegrityError("FlowJo result table copy failed hash verification")
                    os.replace(source_temp, source_target)
                    source_created = True
                source_target.chmod(0o444)
                media_type = {
                    ".xls": "application/vnd.ms-excel",
                    ".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    ".csv": "text/csv",
                    ".tsv": "text/tab-separated-values",
                    ".txt": "text/plain",
                }.get(suffix, "application/octet-stream")
                source_artifact = {
                    "id": uuid.uuid4().hex, "project_id": project_id,
                    "relative_path": source_relative, "kind": "raw_flowjo_result_table",
                    "sha256": source_hash, "size_bytes": source_target.stat().st_size,
                    "media_type": media_type, "read_only": 1, "created_at": now,
                }

            result_document = {
                "schema_version": "1.0",
                "import_id": import_id,
                "project_id": project_id,
                "batch_id": batch_id,
                "batch_name": batch["name"],
                "source_artifact_id": source_artifact["id"],
                "source_filename": parsed_table.get("source_filename") or source_path.name,
                "source_sha256": source_hash,
                "source_sheet": parsed_table.get("sheet_name", ""),
                "flowjo_import_id": linked_flowjo_import,
                "gate_label": str(parsed_table.get("gate_label", "FlowJo per-well percentage")),
                "metric": parsed_table.get("metric", "USER_SUPPLIED_PERCENTAGE"),
                "measurement_basis": "PERCENTAGE_ONLY",
                "percentage_unit": "percent",
                "denominator_available": False,
                "compensation": {
                    "status": compensation_status,
                    "source": "USER_REPORTED_RESULT_TABLE_CONTEXT",
                    "note": compensation_note.strip(),
                },
                "status": "IMPORTED_PENDING_MANUAL_QC",
                "expected_sample_count": len(sample_rows),
                "imported_sample_count": len(mapped_observations),
                "missing_sample_ids": sorted(
                    sample["id"] for well_id, sample in samples_by_well.items()
                    if well_id not in set(supplied_wells)
                ),
                "observations": mapped_observations,
                "parser_provenance": {
                    "sheet_name": parsed_table.get("sheet_name", ""),
                    "gate_label": parsed_table.get("gate_label", ""),
                    "numeric_scale": parsed_table.get("numeric_scale", "AUTO"),
                    "decimal_separator": parsed_table.get("decimal_separator", "."),
                    "decimal_separator_applies_to_text_only": True,
                    "first_two_columns_used": True,
                    "extra_columns_ignored": True,
                    "raw_percentage_values_preserved": True,
                    "event_counts_not_inferred": True,
                },
                "created_at": now,
            }
            result_target.parent.mkdir(parents=True, exist_ok=True)
            result_temp = result_target.with_name(result_target.name + f".{uuid.uuid4().hex}.tmp")
            write_json(result_temp, result_document)
            result_hash = sha256_file(result_temp)
            result_size = result_temp.stat().st_size
            os.replace(result_temp, result_target)
            result_target.chmod(0o444)
            result_created = True
            result_artifact = {
                "id": uuid.uuid4().hex, "project_id": project_id,
                "relative_path": result_relative, "kind": "flowjo_percentage_result",
                "sha256": result_hash, "size_bytes": result_size,
                "media_type": "application/json", "read_only": 1, "created_at": now,
            }
            manifest_path = folder / "project_manifest.json"
            original_manifest = manifest_path.read_bytes()
            artifacts = [result_artifact]
            if not source_artifact_row:
                artifacts.insert(0, source_artifact)
            try:
                revision = self._append_revision(
                    db, project_id, expected_revision, "flowjo_percentage_import", import_id,
                    {
                        "batch_id": batch_id,
                        "source_sha256": source_hash,
                        "result_sha256": result_hash,
                        "flowjo_import_id": linked_flowjo_import,
                        "gate_label": result_document["gate_label"],
                        "sample_count": len(mapped_observations),
                        "denominator_available": False,
                    },
                )
                if not source_artifact_row:
                    db.execute(
                        "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                        tuple(source_artifact[key] for key in (
                            "id", "project_id", "relative_path", "kind", "sha256",
                            "size_bytes", "media_type", "read_only", "created_at",
                        )),
                    )
                db.execute(
                    "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                    tuple(result_artifact[key] for key in (
                        "id", "project_id", "relative_path", "kind", "sha256",
                        "size_bytes", "media_type", "read_only", "created_at",
                    )),
                )
                manifest = json.loads(original_manifest)
                known_ids = {item["id"] for item in manifest["artifacts"]}
                manifest["artifacts"].extend({
                    key: artifact[key] for key in (
                        "id", "kind", "relative_path", "sha256", "size_bytes", "media_type", "created_at"
                    )
                } | {"read_only": True} for artifact in artifacts if artifact["id"] not in known_ids)
                write_json(manifest_path, manifest)
                db.execute(
                    "INSERT INTO flowjo_percentage_imports VALUES (?,?,?,?,?,?,?,?,?)",
                    (
                        import_id, project_id, batch_id, source_artifact["id"],
                        result_artifact["id"], linked_flowjo_import, "IMPORTED",
                        revision, now,
                    ),
                )
                db.commit()
            except Exception:
                try:
                    manifest_path.write_bytes(original_manifest)
                except OSError:
                    pass
                for path, created in ((result_target, result_created), (source_target, source_created)):
                    if created:
                        _discard_unregistered_artifact(path)
                raise
        return {
            "import_id": import_id,
            "batch_id": batch_id,
            "flowjo_import_id": linked_flowjo_import,
            "status": "IMPORTED_PENDING_MANUAL_QC",
            "project_revision": revision,
            "source_artifact": source_artifact,
            "result_artifact": result_artifact,
            "result": result_document,
        }

    def _matching_flowjo_import_for_batch(
        self, db: sqlite3.Connection, project_id: str, batch_artifact_ids: set[str]
    ) -> str | None:
        matches = []
        rows = db.execute(
            "SELECT id,preview_artifact_id,status FROM flowjo_imports "
            "WHERE project_id=? AND status='COMPLETE' ORDER BY created_at,id",
            (project_id,),
        ).fetchall()
        for row in rows:
            preview_artifact = db.execute(
                "SELECT relative_path FROM artifacts WHERE id=? AND project_id=?",
                (row["preview_artifact_id"], project_id),
            ).fetchone()
            if preview_artifact is None:
                continue
            preview_path = self.projects_root / project_id / preview_artifact["relative_path"]
            if not preview_path.is_file() or sha256_file(preview_path) != db.execute(
                "SELECT sha256 FROM artifacts WHERE id=?", (row["preview_artifact_id"],)
            ).fetchone()["sha256"]:
                continue
            preview = read_json(preview_path)
            input_ids = set(preview.get("flowjo_input_artifact_ids") or [])
            if input_ids == batch_artifact_ids:
                matches.append(row["id"])
        return matches[0] if len(matches) == 1 else None

    def list_flowjo_percentage_imports(self, project_id: str) -> list[dict[str, Any]]:
        self._row(project_id)
        with self._connect() as db:
            rows = db.execute(
                "SELECT * FROM flowjo_percentage_imports WHERE project_id=? ORDER BY created_at,id",
                (project_id,),
            ).fetchall()
        result = []
        for row in rows:
            item = dict(row)
            result_path, result_artifact = self.artifact_path(project_id, item["result_artifact_id"])
            source_path, source_artifact = self.artifact_path(project_id, item["source_artifact_id"])
            item["source_artifact"] = source_artifact
            item["result_artifact"] = result_artifact
            item["result"] = read_json(result_path)
            if not source_path.is_file():
                raise IntegrityError("FlowJo result table source artifact is missing")
            result.append(item)
        return result

    def save_flowjo_percentage_effect_analysis(
        self,
        project_id: str,
        *,
        import_id: str,
        document: dict[str, Any],
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        analysis_id = uuid.uuid4().hex
        artifact_id = uuid.uuid4().hex
        relative = f"analyses/flowjo_percentage_effects/{analysis_id}/effect_results.json"
        target = self.projects_root / project_id / relative
        now = _now()
        result_document = {
            "schema_version": "1.0", "analysis_id": analysis_id,
            "project_id": project_id, "flowjo_percentage_import_id": import_id,
            "status": "COMPLETE", "created_at": now, **document,
        }
        with self._lock, self._connect() as db:
            source = db.execute(
                "SELECT id FROM flowjo_percentage_imports WHERE id=? AND project_id=? AND status='IMPORTED'",
                (import_id, project_id),
            ).fetchone()
            if source is None:
                raise KeyError(import_id)
            target.parent.mkdir(parents=True, exist_ok=True)
            temp = target.with_name(target.name + f".{uuid.uuid4().hex}.tmp")
            manifest_path = self.projects_root / project_id / "project_manifest.json"
            original_manifest = manifest_path.read_bytes()
            target_created = False
            try:
                write_json(temp, result_document)
                digest = sha256_file(temp)
                size_bytes = temp.stat().st_size
                os.replace(temp, target)
                target_created = True
                target.chmod(0o444)
                artifact = {
                    "id": artifact_id, "project_id": project_id, "relative_path": relative,
                    "kind": "flowjo_percentage_effect_analysis", "sha256": digest,
                    "size_bytes": size_bytes, "media_type": "application/json",
                    "read_only": 1, "created_at": now,
                }
                revision = self._append_revision(
                    db, project_id, expected_revision, "flowjo_percentage_effect_analysis",
                    analysis_id,
                    {"import_id": import_id, "result_sha256": digest,
                     "effect_count": len(document.get("effects", []))},
                )
                db.execute(
                    "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                    tuple(artifact[key] for key in (
                        "id", "project_id", "relative_path", "kind", "sha256",
                        "size_bytes", "media_type", "read_only", "created_at",
                    )),
                )
                db.execute(
                    "INSERT INTO flowjo_percentage_effect_analyses VALUES (?,?,?,?,?,?,?)",
                    (analysis_id, project_id, import_id, artifact_id, "COMPLETE", revision, now),
                )
                manifest = json.loads(original_manifest)
                manifest["artifacts"].append({
                    key: artifact[key] for key in (
                        "id", "kind", "relative_path", "sha256", "size_bytes", "media_type", "created_at"
                    )
                } | {"read_only": True})
                write_json(manifest_path, manifest)
                db.commit()
            except Exception:
                try:
                    manifest_path.write_bytes(original_manifest)
                except OSError:
                    pass
                if target_created:
                    _discard_unregistered_artifact(target)
                raise
            finally:
                temp.unlink(missing_ok=True)
        return {
            "analysis_id": analysis_id, "import_id": import_id,
            "project_revision": revision, "status": "COMPLETE",
            "result_artifact": {key: artifact[key] for key in (
                "id", "kind", "sha256", "size_bytes", "media_type"
            )},
        }

    def list_flowjo_percentage_effect_analyses(
        self, project_id: str, import_id: str
    ) -> list[dict[str, Any]]:
        self._row(project_id)
        with self._connect() as db:
            source = db.execute(
                "SELECT id FROM flowjo_percentage_imports WHERE id=? AND project_id=?",
                (import_id, project_id),
            ).fetchone()
            if source is None:
                raise KeyError(import_id)
            rows = db.execute(
                "SELECT * FROM flowjo_percentage_effect_analyses "
                "WHERE project_id=? AND import_id=? ORDER BY created_at,id",
                (project_id, import_id),
            ).fetchall()
        result = []
        for row in rows:
            item = dict(row)
            path, artifact = self.artifact_path(project_id, item["result_artifact_id"])
            item["result_artifact"] = artifact
            item["result"] = read_json(path)
            result.append(item)
        return result

    def register_compensation_matrix(
        self,
        project_id: str,
        matrix_document: dict[str, Any],
        expected_revision: int,
    ) -> dict[str, Any]:
        """Store a validated matrix as an immutable project artifact."""
        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        if matrix_document.get("schema_version") != "1.0":
            raise ValueError("Unsupported compensation matrix schema_version")
        matrix_id = uuid.uuid4().hex
        now = _now()
        folder = self.projects_root / project_id
        payload = json.dumps(matrix_document, ensure_ascii=False, indent=2).encode("utf-8")
        digest = hashlib.sha256(payload).hexdigest()
        size_bytes = len(payload)
        relative = f"compensation/matrices/{digest}.json"
        registered = next((
            item for item in project["artifacts"]
            if item["kind"] == "compensation_matrix" and item["relative_path"] == relative
        ), None)
        if registered is not None:
            if registered["sha256"] != digest or registered["size_bytes"] != size_bytes:
                raise IntegrityError("Registered compensation matrix content differs from its digest")
            return {
                "artifact_id": registered["id"],
                "project_revision": project["revision"],
                "already_registered": True,
                "artifact": {key: registered[key] for key in (
                    "id", "kind", "sha256", "size_bytes", "media_type", "created_at"
                )},
                "detectors": matrix_document["detectors"],
                "source_filename": matrix_document["source_filename"],
            }
        target = folder / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        temp = target.with_name(target.name + f".{uuid.uuid4().hex}.tmp")
        try:
            temp.write_bytes(payload)

            def target_matches_matrix() -> bool:
                return (
                    not target.is_symlink()
                    and target.is_file()
                    and target.stat().st_size == size_bytes
                    and sha256_file(target) == digest
                )

            if target.exists() or target.is_symlink():
                if not target_matches_matrix():
                    raise IntegrityError(
                        "Unregistered compensation matrix path does not match uploaded content"
                    )
            else:
                try:
                    os.replace(temp, target)
                except OSError:
                    # Another import may have published these exact immutable bytes.
                    if not target_matches_matrix():
                        raise
                if not target_matches_matrix():
                    raise IntegrityError("Published compensation matrix failed hash verification")
            target.chmod(0o444)
        finally:
            temp.unlink(missing_ok=True)

        artifact = {
            "id": matrix_id,
            "project_id": project_id,
            "relative_path": relative,
            "kind": "compensation_matrix",
            "sha256": digest,
            "size_bytes": size_bytes,
            "media_type": "application/json",
            "read_only": 1,
            "created_at": now,
        }
        details = {
            "source_filename": matrix_document.get("source_filename"),
            "import_format": matrix_document.get("import_format"),
            "detectors": matrix_document.get("detectors"),
            "detector_count": len(matrix_document.get("detectors", [])),
            "sha256": digest,
        }
        with self._lock, self._connect() as db:
            revision = self._append_revision(
                db, project_id, expected_revision, "compensation_matrix", matrix_id, details
            )
            db.execute(
                "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                (
                    artifact["id"], artifact["project_id"], artifact["relative_path"],
                    artifact["kind"], artifact["sha256"], artifact["size_bytes"],
                    artifact["media_type"], artifact["read_only"], artifact["created_at"],
                ),
            )
            manifest_path = folder / "project_manifest.json"
            manifest = read_json(manifest_path)
            manifest["artifacts"].append({
                key: artifact[key]
                for key in (
                    "id", "kind", "relative_path", "sha256", "size_bytes",
                    "media_type", "created_at",
                )
            } | {"read_only": True})
            write_json(manifest_path, manifest)
            db.commit()
        return {
            "artifact_id": matrix_id,
            "project_revision": revision,
            "already_registered": False,
            "artifact": {key: artifact[key] for key in (
                "id", "kind", "sha256", "size_bytes", "media_type", "created_at"
            )},
            "detectors": matrix_document["detectors"],
            "source_filename": matrix_document["source_filename"],
        }

    def save_gate_revision(
        self,
        project_id: str,
        *,
        name: str,
        gates: list[dict[str, Any]],
        transforms: list[dict[str, Any]] | None = None,
        compensation_decision: str,
        compensation_artifact_id: str | None,
        parent_revision_id: str | None,
        flowjo_rebuild: dict[str, Any] | None = None,
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        with self._lock, self._connect() as db:
            if (
                parent_revision_id
                and not db.execute(
                    "SELECT id FROM gate_revisions WHERE id=? AND project_id=?",
                    (parent_revision_id, project_id),
                ).fetchone()
            ):
                raise KeyError(parent_revision_id)
            if compensation_artifact_id:
                compensation_row = db.execute(
                    "SELECT kind,relative_path,sha256 FROM artifacts WHERE id=? AND project_id=?",
                    (compensation_artifact_id, project_id),
                ).fetchone()
                if compensation_row is None:
                    raise KeyError(compensation_artifact_id)
                if compensation_row["kind"] != "compensation_matrix":
                    raise GateError("External compensation artifact must be a registered compensation matrix")
                matrix_path = self.projects_root / project_id / compensation_row["relative_path"]
                if not matrix_path.is_file() or sha256_file(matrix_path) != compensation_row["sha256"]:
                    raise IntegrityError("Compensation matrix artifact failed integrity verification")
            rebuild_record = None
            if flowjo_rebuild is not None:
                flowjo_import_id = str(flowjo_rebuild.get("flowjo_import_id", ""))
                reviewer = str(flowjo_rebuild.get("reviewer", "")).strip()
                reason = str(flowjo_rebuild.get("reason", "")).strip()
                if not reviewer or not reason:
                    raise ValueError("FlowJo gate reconstruction requires a reviewer and reason")
                flowjo_row = db.execute(
                    "SELECT status,unresolved_json,preview_artifact_id FROM flowjo_imports "
                    "WHERE id=? AND project_id=?",
                    (flowjo_import_id, project_id),
                ).fetchone()
                if flowjo_row is None:
                    raise KeyError(flowjo_import_id)
                unresolved_nodes = json.loads(flowjo_row["unresolved_json"])
                if flowjo_row["status"] != "WAITING_FOR_USER" or not unresolved_nodes:
                    raise ValueError("Only an unresolved FlowJo import can be rebuilt this way")
                preview_artifact = db.execute(
                    "SELECT relative_path,sha256 FROM artifacts WHERE id=? AND project_id=? "
                    "AND kind='flowjo_import_preview'",
                    (flowjo_row["preview_artifact_id"], project_id),
                ).fetchone()
                if preview_artifact is None:
                    raise IntegrityError("FlowJo preview artifact is missing")
                preview_path = self.projects_root / project_id / preview_artifact["relative_path"]
                if not preview_path.is_file() or sha256_file(preview_path) != preview_artifact["sha256"]:
                    raise IntegrityError("FlowJo preview artifact failed integrity verification")
                preview_document = read_json(preview_path)
                input_artifact_ids = preview_document.get("flowjo_input_artifact_ids")
                if (
                    not isinstance(input_artifact_ids, list)
                    or not input_artifact_ids
                    or len(input_artifact_ids) > 128
                    or len(set(input_artifact_ids)) != len(input_artifact_ids)
                ):
                    raise ValueError("Reimport the FlowJo workspace to bind its FCS artifact IDs")
                for artifact_id in input_artifact_ids:
                    artifact_row = db.execute(
                        "SELECT kind FROM artifacts WHERE id=? AND project_id=?",
                        (artifact_id, project_id),
                    ).fetchone()
                    if artifact_row is None or artifact_row["kind"] != "raw_fcs":
                        raise IntegrityError("FlowJo import references an invalid FCS artifact")
                rebuild_record = {
                    "flowjo_import_id": flowjo_import_id,
                    "preview_artifact_id": flowjo_row["preview_artifact_id"],
                    "preview_sha256": preview_artifact["sha256"],
                    "unresolved_node_count": len(unresolved_nodes),
                    "input_artifact_ids": input_artifact_ids,
                    "reviewer": reviewer,
                    "reason": reason,
                }
            channels: set[str] = set()
            for row in db.execute(
                "SELECT metadata_json FROM pending_uploads WHERE project_id=?", (project_id,)
            ):
                metadata = json.loads(row[0])
                channels.update(channel.get("pnn", "") for channel in metadata.get("channels", []))
                channels.update(channel.get("pns", "") for channel in metadata.get("channels", []))
            normalized = validate_gate_tree(
                gates,
                available_channels=channels or None,
                compensation_decision=compensation_decision,
                compensation_artifact_id=compensation_artifact_id,
                transforms=transforms or [],
            )
            version = db.execute(
                "SELECT COALESCE(MAX(version),0)+1 FROM gate_revisions WHERE project_id=?",
                (project_id,),
            ).fetchone()[0]
            revision_id = uuid.uuid4().hex
            project_revision = self._append_revision(
                db,
                project_id,
                expected_revision,
                "gate_revision",
                revision_id,
                {"name": name, "version": version, "parent_revision_id": parent_revision_id,
                 "flowjo_rebuild": rebuild_record},
            )
            now = _now()
            db.execute(
                "INSERT INTO gate_revisions VALUES (?,?,?,?,?,?,?,?,?,?)",
                (
                    revision_id,
                    project_id,
                    parent_revision_id,
                    name,
                    version,
                    json.dumps(
                        {"gate_tree": normalized["gates"], "transforms": normalized["transforms"],
                         "flowjo_rebuild": rebuild_record},
                        ensure_ascii=False,
                    ),
                    compensation_decision,
                    compensation_artifact_id,
                    project_revision,
                    now,
                ),
            )
            db.commit()
        return {
            "gate_revision_id": revision_id,
            "name": name,
            "version": version,
            "parent_revision_id": parent_revision_id,
            "project_revision": project_revision,
            "validated": True,
            "gate_tree": normalized["gates"],
            "transforms": normalized["transforms"],
            "compensation_decision": compensation_decision,
            "compensation_artifact_id": compensation_artifact_id,
            "flowjo_rebuild": rebuild_record,
        }

    def list_gate_revisions(self, project_id: str) -> list[dict[str, Any]]:
        self._row(project_id)
        with self._connect() as db:
            rows = db.execute(
                "SELECT * FROM gate_revisions WHERE project_id=? ORDER BY version", (project_id,)
            ).fetchall()
        output = []
        for row in rows:
            item = dict(row)
            saved = json.loads(item.pop("gates_json"))
            if isinstance(saved, list):
                item["gate_tree"] = saved
                item["transforms"] = []
            else:
                item["gate_tree"] = saved["gate_tree"]
                item["transforms"] = saved.get("transforms", [])
                item["flowjo_rebuild"] = saved.get("flowjo_rebuild")
            output.append(item)
        return output

    def validate_gate_analysis_inputs(
        self, project_id: str, gate_revision_id: str, input_artifact_ids: list[str]
    ) -> None:
        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        if not input_artifact_ids or len(input_artifact_ids) > 128:
            raise ValueError("Select between 1 and 128 FCS artifacts")
        if len(set(input_artifact_ids)) != len(input_artifact_ids):
            raise ValueError("FCS artifact IDs must be unique")
        with self._connect() as db:
            gate_row = db.execute(
                "SELECT gates_json FROM gate_revisions WHERE id=? AND project_id=?",
                (gate_revision_id, project_id),
            ).fetchone()
            if gate_row is None:
                raise KeyError(gate_revision_id)
            saved_gate_revision = json.loads(gate_row["gates_json"])
            flowjo_rebuild = (
                saved_gate_revision.get("flowjo_rebuild")
                if isinstance(saved_gate_revision, dict) else None
            )
            if flowjo_rebuild:
                expected_artifact_ids = flowjo_rebuild.get("input_artifact_ids", [])
                if (len(input_artifact_ids) != len(expected_artifact_ids)
                        or set(input_artifact_ids) != set(expected_artifact_ids)):
                    raise ValueError(
                        "A FlowJo gate rebuild must be analyzed against the same registered FCS artifacts"
                    )
                source_import = db.execute(
                    "SELECT status,preview_artifact_id FROM flowjo_imports "
                    "WHERE id=? AND project_id=?",
                    (flowjo_rebuild.get("flowjo_import_id"), project_id),
                ).fetchone()
                if source_import is None or source_import["status"] != "WAITING_FOR_USER":
                    raise ValueError("The linked FlowJo import is no longer awaiting reconstruction")
                source_preview = db.execute(
                    "SELECT sha256 FROM artifacts WHERE id=? AND project_id=?",
                    (source_import["preview_artifact_id"], project_id),
                ).fetchone()
                if source_preview is None or source_preview["sha256"] != flowjo_rebuild.get("preview_sha256"):
                    raise IntegrityError("The linked FlowJo preview changed after gate reconstruction")
            for artifact_id in input_artifact_ids:
                artifact = db.execute(
                    "SELECT kind FROM artifacts WHERE id=? AND project_id=?",
                    (artifact_id, project_id),
                ).fetchone()
                if artifact is None or artifact["kind"] != "raw_fcs":
                    raise ValueError("Only registered raw FCS artifacts can be analyzed")
                matches = db.execute(
                    "SELECT id FROM samples WHERE project_id=? AND artifact_id=?",
                    (project_id, artifact_id),
                ).fetchall()
                if not matches:
                    raise ValueError(f"SAMPLE_MAPPING_REQUIRED: FCS artifact {artifact_id} is unmatched")
                if len(matches) > 1:
                    raise ValueError(f"FCS artifact {artifact_id} is linked to multiple sample records")

    def save_gate_analysis(
        self,
        project_id: str,
        *,
        gate_revision_id: str,
        input_artifact_ids: list[str],
        worker_result: dict[str, Any],
        expected_revision: int,
    ) -> dict[str, Any]:
        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        if not input_artifact_ids or len(input_artifact_ids) > 128:
            raise ValueError("Select between 1 and 128 FCS artifacts")
        if len(set(input_artifact_ids)) != len(input_artifact_ids):
            raise ValueError("FCS artifact IDs must be unique")
        self.validate_gate_analysis_inputs(project_id, gate_revision_id, input_artifact_ids)
        samples = worker_result.get("samples")
        if not isinstance(samples, list) or len(samples) != len(input_artifact_ids):
            raise ValueError("FlowKit result count does not match the selected FCS files")
        analysis_id = uuid.uuid4().hex
        now = _now()
        folder = self.projects_root / project_id
        relative = f"analyses/gate_runs/{analysis_id}/gate_counts.json"
        target = folder / relative
        artifact_id = uuid.uuid4().hex
        with self._lock, self._connect() as db:
            revision_row = db.execute(
                "SELECT id,version,gates_json,compensation_decision,compensation_artifact_id "
                "FROM gate_revisions WHERE id=? AND project_id=?",
                (gate_revision_id, project_id),
            ).fetchone()
            if revision_row is None:
                raise KeyError(gate_revision_id)
            artifacts = []
            for input_id in input_artifact_ids:
                artifact = db.execute(
                    "SELECT a.id,a.kind,a.sha256,s.id AS sample_id,s.role,s.batch_id,"
                    "s.antigen_construct_id,s.antibody_construct_id,s.replicate_type,"
                    "s.replicate_index,s.condition_json,b.assay,b.date "
                    "FROM artifacts a JOIN samples s ON s.artifact_id=a.id "
                    "JOIN assay_batches b ON b.id=s.batch_id "
                    "WHERE a.id=? AND a.project_id=? AND s.project_id=?",
                    (input_id, project_id, project_id),
                ).fetchone()
                if artifact is None or artifact["kind"] != "raw_fcs":
                    raise ValueError("Only registered raw FCS artifacts can be analyzed")
                artifacts.append(dict(artifact))
            result_samples = []
            for artifact, sample in zip(artifacts, samples, strict=True):
                if not isinstance(sample.get("event_count"), int) or sample["event_count"] < 1:
                    raise ValueError("FlowKit returned an invalid event count")
                result_samples.append({
                    "artifact_id": artifact["id"],
                    "input_sha256": artifact["sha256"],
                    "sample": {
                        "sample_id": artifact["sample_id"], "role": artifact["role"],
                        "batch_id": artifact["batch_id"], "assay": artifact["assay"],
                        "date": artifact["date"],
                        "antigen_construct_id": artifact["antigen_construct_id"],
                        "antibody_construct_id": artifact["antibody_construct_id"],
                        "replicate_type": artifact["replicate_type"],
                        "replicate_index": artifact["replicate_index"],
                        "condition": json.loads(artifact["condition_json"]),
                    },
                    **sample,
                })
            stored_gate_revision = json.loads(revision_row["gates_json"])
            if isinstance(stored_gate_revision, list):
                stored_gate_tree, stored_transforms = stored_gate_revision, []
            else:
                stored_gate_tree = stored_gate_revision["gate_tree"]
                stored_transforms = stored_gate_revision.get("transforms", [])
            flowjo_rebuild = (
                stored_gate_revision.get("flowjo_rebuild")
                if isinstance(stored_gate_revision, dict) else None
            )
            compensation_provenance: dict[str, Any] = {
                "decision": revision_row["compensation_decision"],
                "artifact_id": revision_row["compensation_artifact_id"],
            }
            if revision_row["compensation_artifact_id"]:
                matrix_row = db.execute(
                    "SELECT relative_path,kind,sha256 FROM artifacts WHERE id=? AND project_id=?",
                    (revision_row["compensation_artifact_id"], project_id),
                ).fetchone()
                if matrix_row is None or matrix_row["kind"] != "compensation_matrix":
                    raise IntegrityError("Gate revision references a missing compensation matrix")
                matrix_path = folder / matrix_row["relative_path"]
                if not matrix_path.is_file() or sha256_file(matrix_path) != matrix_row["sha256"]:
                    raise IntegrityError("Compensation matrix artifact failed integrity verification")
                matrix_document = read_json(matrix_path)
                compensation_provenance.update({
                    "artifact_sha256": matrix_row["sha256"],
                    "detectors": matrix_document.get("detectors", []),
                })
            result_document = {
                "schema_version": "1.0",
                "analysis_id": analysis_id,
                "project_id": project_id,
                "gate_revision_id": gate_revision_id,
                "gate_revision_version": revision_row["version"],
                "gate_tree": stored_gate_tree,
                "transforms": stored_transforms,
                "flowjo_rebuild": flowjo_rebuild,
                "compensation": compensation_provenance,
                "status": "COMPLETE",
                "input_artifact_ids": input_artifact_ids,
                "frozen_hashes": project["frozen_hashes"],
                "samples": result_samples,
                "qc_policy": worker_result.get("qc_policy", {}),
                "event_exclusions": worker_result.get("event_exclusions"),
                "provenance": worker_result.get("provenance", {}),
                "created_at": now,
            }
            write_json(target, result_document)
            digest = sha256_file(target)
            size_bytes = target.stat().st_size
            target.chmod(0o444)
            artifact = {
                "id": artifact_id, "project_id": project_id, "relative_path": relative,
                "kind": "gate_analysis", "sha256": digest, "size_bytes": size_bytes,
                "media_type": "application/json", "read_only": 1, "created_at": now,
            }
            try:
                project_revision = self._append_revision(
                    db, project_id, expected_revision, "gate_analysis", analysis_id,
                    {"gate_revision_id": gate_revision_id, "input_artifact_ids": input_artifact_ids,
                     "result_sha256": digest,
                     "event_exclusions": worker_result.get("event_exclusions"),
                     "flowjo_import_id": flowjo_rebuild.get("flowjo_import_id")
                     if flowjo_rebuild else None},
                )
                db.execute(
                    "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                    (artifact_id, project_id, relative, "gate_analysis", digest, size_bytes,
                     "application/json", 1, now),
                )
                db.execute(
                    "INSERT INTO gate_analyses VALUES (?,?,?,?,?,?,?,?)",
                    (analysis_id, project_id, gate_revision_id, json.dumps(input_artifact_ids),
                     artifact_id, "COMPLETE", project_revision, now),
                )
                manifest_path = folder / "project_manifest.json"
                manifest = read_json(manifest_path)
                manifest["artifacts"].append({
                    key: artifact[key]
                    for key in ("id", "kind", "relative_path", "sha256", "size_bytes", "media_type", "created_at")
                } | {"read_only": True})
                write_json(manifest_path, manifest)
                db.commit()
            except Exception:
                _discard_unregistered_artifact(target)
                raise
        return {
            "analysis_id": analysis_id,
            "gate_revision_id": gate_revision_id,
            "project_revision": project_revision,
            "status": "COMPLETE",
            "result_artifact": {key: artifact[key] for key in (
                "id", "kind", "sha256", "size_bytes", "media_type"
            )},
            "samples": result_samples,
        }

    def list_gate_analyses(self, project_id: str) -> list[dict[str, Any]]:
        self._row(project_id)
        with self._connect() as db:
            rows = db.execute(
                "SELECT * FROM gate_analyses WHERE project_id=? ORDER BY created_at,id",
                (project_id,),
            ).fetchall()
        return [{
            **dict(row),
            "input_artifact_ids": json.loads(row["input_artifact_ids_json"]),
        } for row in rows]

    def save_effect_analysis(
        self,
        project_id: str,
        *,
        gate_analysis_id: str,
        positive_gate_id: str,
        document: dict[str, Any],
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        analysis_id = uuid.uuid4().hex
        artifact_id = uuid.uuid4().hex
        relative = f"analyses/effect_runs/{analysis_id}/effect_results.json"
        target = self.projects_root / project_id / relative
        now = _now()
        result_document = {
            "schema_version": "1.0", "analysis_id": analysis_id,
            "project_id": project_id, "gate_analysis_id": gate_analysis_id,
            "positive_gate_id": positive_gate_id, "status": "COMPLETE",
            "created_at": now, **document,
        }
        with self._lock, self._connect() as db:
            source = db.execute(
                "SELECT id FROM gate_analyses WHERE id=? AND project_id=? AND status='COMPLETE'",
                (gate_analysis_id, project_id),
            ).fetchone()
            if source is None:
                raise KeyError(gate_analysis_id)
            target.parent.mkdir(parents=True, exist_ok=True)
            temp = target.with_name(target.name + f".{uuid.uuid4().hex}.tmp")
            write_json(temp, result_document)
            digest = sha256_file(temp)
            size_bytes = temp.stat().st_size
            os.replace(temp, target)
            target.chmod(0o444)
            artifact = {
                "id": artifact_id, "project_id": project_id, "relative_path": relative,
                "kind": "effect_analysis", "sha256": digest, "size_bytes": size_bytes,
                "media_type": "application/json", "read_only": 1, "created_at": now,
            }
            try:
                project_revision = self._append_revision(
                    db, project_id, expected_revision, "effect_analysis", analysis_id,
                    {"gate_analysis_id": gate_analysis_id, "positive_gate_id": positive_gate_id,
                     "result_sha256": digest},
                )
                db.execute(
                    "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                    (artifact_id, project_id, relative, "effect_analysis", digest, size_bytes,
                     "application/json", 1, now),
                )
                db.execute(
                    "INSERT INTO effect_analyses VALUES (?,?,?,?,?,?,?,?)",
                    (analysis_id, project_id, gate_analysis_id, positive_gate_id, artifact_id,
                     "COMPLETE", project_revision, now),
                )
                manifest_path = self.projects_root / project_id / "project_manifest.json"
                manifest = read_json(manifest_path)
                manifest["artifacts"].append({
                    key: artifact[key] for key in (
                        "id", "kind", "relative_path", "sha256", "size_bytes", "media_type", "created_at"
                    )
                } | {"read_only": True})
                write_json(manifest_path, manifest)
                db.commit()
            except Exception:
                _discard_unregistered_artifact(target)
                raise
        return {
            "analysis_id": analysis_id, "gate_analysis_id": gate_analysis_id,
            "positive_gate_id": positive_gate_id, "project_revision": project_revision,
            "status": "COMPLETE", "result_artifact": {
                key: artifact[key] for key in ("id", "kind", "sha256", "size_bytes", "media_type")
            },
        }

    def list_effect_analyses(self, project_id: str) -> list[dict[str, Any]]:
        self._row(project_id)
        with self._connect() as db:
            rows = db.execute(
                "SELECT * FROM effect_analyses WHERE project_id=? ORDER BY created_at,id",
                (project_id,),
            ).fetchall()
        return [dict(row) for row in rows]

    def create_paraplume_task(
        self,
        project_id: str,
        *,
        antibody_construct_ids: list[str],
        model_fingerprint: str,
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        if not antibody_construct_ids or len(antibody_construct_ids) > 20:
            raise ValueError("Select between 1 and 20 antibody constructs")
        if len(set(antibody_construct_ids)) != len(antibody_construct_ids):
            raise ValueError("Antibody construct IDs must be unique")
        task_id = uuid.uuid4().hex
        now = _now()
        with self._lock, self._connect() as db:
            placeholders = ",".join("?" for _ in antibody_construct_ids)
            rows = db.execute(
                f"SELECT * FROM antibody_constructs WHERE project_id=? AND id IN ({placeholders})",
                [project_id, *antibody_construct_ids],
            ).fetchall()
            by_id = {row["id"]: row for row in rows}
            missing = [item for item in antibody_construct_ids if item not in by_id]
            if missing:
                raise KeyError(missing[0])
            inputs = []
            for construct_id in antibody_construct_ids:
                row = by_id[construct_id]
                vh_hash = sha256_bytes(row["vh_sequence"].encode("ascii"))
                vl_hash = sha256_bytes(row["vl_sequence"].encode("ascii"))
                cache_key = sha256_bytes(
                    json.dumps(
                        {
                            "model_fingerprint": model_fingerprint,
                            "vh_sha256": vh_hash,
                            "vl_sha256": vl_hash,
                        },
                        sort_keys=True,
                        separators=(",", ":"),
                    ).encode("utf-8")
                )
                inputs.append(
                    {
                        "antibody_construct_id": construct_id,
                        "name": row["name"],
                        "role": row["role"],
                        "vh_sha256": vh_hash,
                        "vl_sha256": vl_hash,
                        "cache_key": cache_key,
                    }
                )
            revision = self._append_revision(
                db,
                project_id,
                expected_revision,
                "paraplume_task",
                task_id,
                {"construct_ids": antibody_construct_ids, "model_fingerprint": model_fingerprint},
            )
            task_input = {
                "schema_version": "1.0",
                "model_version": "2.0.0",
                "model_size": "large",
                "project_revision": revision,
                "model_fingerprint": model_fingerprint,
                "constructs": inputs,
            }
            db.execute(
                """INSERT INTO experiment_tasks
                   (id,project_id,kind,status,progress,input_json,output_json,error_code,created_at,updated_at,owner_pid)
                   VALUES (?,?,?,'QUEUED',0,?,NULL,NULL,?,?,?)""",
                (task_id, project_id, "PARAPLUME", json.dumps(task_input), now, now, os.getpid()),
            )
            db.commit()
        return self.get_experiment_task(task_id)

    def get_experiment_task(self, task_id: str) -> dict[str, Any]:
        with self._connect() as db:
            row = db.execute("SELECT * FROM experiment_tasks WHERE id=?", (task_id,)).fetchone()
        if row is None:
            raise KeyError(task_id)
        item = dict(row)
        item.pop("owner_pid", None)
        for key in ("input_json", "output_json"):
            value = item.pop(key)
            item[key.removesuffix("_json")] = json.loads(value) if value else None
        return item

    def update_experiment_task(
        self,
        task_id: str,
        *,
        status: str,
        progress: float,
        output: dict[str, Any] | None = None,
        error_code: str | None = None,
    ) -> dict[str, Any]:
        if status not in {"QUEUED", "RUNNING", "COMPLETE", "FAILED", "INTERRUPTED"}:
            raise ValueError("Unsupported experiment task status")
        with self._lock, self._connect() as db:
            cursor = db.execute(
                """UPDATE experiment_tasks SET status=?,progress=?,output_json=?,error_code=?,updated_at=?
                   WHERE id=?""",
                (
                    status,
                    min(1.0, max(0.0, float(progress))),
                    json.dumps(output, ensure_ascii=False) if output is not None else None,
                    error_code,
                    _now(),
                    task_id,
                ),
            )
            if cursor.rowcount != 1:
                raise KeyError(task_id)
            db.commit()
        return self.get_experiment_task(task_id)

    def complete_paraplume_task(
        self,
        task_id: str,
        *,
        results_by_cache_key: dict[str, dict[str, Any]],
        provenance_by_cache_key: dict[str, dict[str, Any]],
        cache_hits: set[str],
    ) -> dict[str, Any]:
        task = self.get_experiment_task(task_id)
        if task["kind"] != "PARAPLUME" or task["status"] != "RUNNING":
            raise ValueError("Task is not a running Paraplume task")
        project_id = task["project_id"]
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        with self._lock:
            with self._connect() as db:
                constructs_by_id = {
                    row["id"]: dict(row)
                    for row in db.execute(
                        "SELECT * FROM antibody_constructs WHERE project_id=?", (project_id,)
                    )
                }
            records = []
            for item in task["input"]["constructs"]:
                cached = results_by_cache_key.get(item["cache_key"])
                if not isinstance(cached, dict):
                    raise ValueError("A requested Paraplume cache key has no result")
                construct = constructs_by_id.get(item["antibody_construct_id"])
                if construct is None:
                    raise KeyError(item["antibody_construct_id"])
                if (
                    cached.get("vh_length") != len(construct["vh_sequence"])
                    or cached.get("vl_length") != len(construct["vl_sequence"])
                ):
                    raise ValueError("Cached Paraplume residue counts do not match construct sequence")
                vh_scores = cached.get("vh_probabilities")
                vl_scores = cached.get("vl_probabilities")
                if not isinstance(vh_scores, list) or not isinstance(vl_scores, list):
                    raise ValueError("Paraplume probability arrays are missing")
                if len(vh_scores) != cached["vh_length"] or len(vl_scores) != cached["vl_length"]:
                    raise ValueError("Paraplume probability arrays have invalid lengths")
                if any(
                    not math.isfinite(float(value)) or not 0 <= float(value) <= 1
                    for value in vh_scores + vl_scores
                ):
                    raise ValueError("Paraplume probabilities must be finite values in [0, 1]")
                records.append(
                    {
                        "antibody_construct_id": item["antibody_construct_id"],
                        "name": construct["name"],
                        "role": construct["role"],
                        "sequence_hashes": {"vh": item["vh_sha256"], "vl": item["vl_sha256"]},
                        "vh_length": cached["vh_length"],
                        "vl_length": cached["vl_length"],
                        "vh_probabilities": cached["vh_probabilities"],
                        "vl_probabilities": cached["vl_probabilities"],
                        "prediction_status": "PREDICTED_UNREVIEWED",
                        "cache_hit": item["cache_key"] in cache_hits,
                        "provenance": provenance_by_cache_key[item["cache_key"]],
                    }
                )
            now = _now()
            artifact_id = uuid.uuid4().hex
            relative = f"analyses/paraplume_runs/{task_id}/paratope_predictions.json"
            folder = self.projects_root / project_id
            target = folder / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            document = {
                "schema_version": "1.0",
                "task_id": task_id,
                "project_id": project_id,
                "status": "COMPLETE",
                "model_version": "2.0.0",
                "model_size": "large",
                "prediction_scope": "antibody-side paratope probabilities; not antigen epitope evidence",
                "records": records,
                "created_at": now,
            }
            temp = target.with_name(target.name + f".{uuid.uuid4().hex}.tmp")
            write_json(temp, document)
            digest = sha256_file(temp)
            size_bytes = temp.stat().st_size
            os.replace(temp, target)
            target.chmod(0o444)
            artifact = {
                "id": artifact_id,
                "project_id": project_id,
                "relative_path": relative,
                "kind": "paraplume_prediction",
                "sha256": digest,
                "size_bytes": size_bytes,
                "media_type": "application/json",
                "read_only": 1,
                "created_at": now,
            }
            try:
                with self._connect() as db:
                    current = db.execute(
                        "SELECT revision FROM experiment_projects WHERE id=?", (project_id,)
                    ).fetchone()
                    if current is None:
                        raise KeyError(project_id)
                    revision = self._append_revision(
                        db,
                        project_id,
                        current["revision"],
                        "paraplume_result",
                        task_id,
                        {"task_id": task_id, "result_sha256": digest},
                    )
                    db.execute(
                        "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)",
                        (
                            artifact_id,
                            project_id,
                            relative,
                            artifact["kind"],
                            digest,
                            size_bytes,
                            artifact["media_type"],
                            1,
                            now,
                        ),
                    )
                    for item, record in zip(task["input"]["constructs"], records, strict=True):
                        cache_key = item["cache_key"]
                        cache_result = {
                            key: record[key]
                            for key in (
                                "vh_length",
                                "vl_length",
                                "vh_probabilities",
                                "vl_probabilities",
                            )
                        }
                        provenance = provenance_by_cache_key[cache_key]
                        db.execute(
                            """INSERT OR IGNORE INTO paraplume_cache
                               (cache_key,result_json,provenance_json,created_at) VALUES (?,?,?,?)""",
                            (
                                cache_key,
                                json.dumps(cache_result),
                                json.dumps(provenance),
                                now,
                            ),
                        )
                        db.execute(
                            """INSERT INTO paraplume_runs
                               (id,project_id,task_id,antibody_construct_id,cache_key,
                                result_artifact_id,cache_hit,created_at)
                               VALUES (?,?,?,?,?,?,?,?)""",
                            (
                                uuid.uuid4().hex,
                                project_id,
                                task_id,
                                item["antibody_construct_id"],
                                cache_key,
                                artifact_id,
                                int(item["cache_key"] in cache_hits),
                                now,
                            ),
                        )
                    manifest_path = folder / "project_manifest.json"
                    manifest = read_json(manifest_path)
                    manifest["artifacts"].append(
                        {
                            key: artifact[key]
                            for key in (
                                "id",
                                "kind",
                                "relative_path",
                                "sha256",
                                "size_bytes",
                                "media_type",
                                "created_at",
                            )
                        }
                        | {"read_only": True}
                    )
                    write_json(manifest_path, manifest)
                    output = {
                        "artifact_id": artifact_id,
                        "sha256": digest,
                        "size_bytes": size_bytes,
                        "project_revision": revision,
                        "record_count": len(records),
                        "records": [
                            {
                                "antibody_construct_id": record["antibody_construct_id"],
                                "prediction_status": record["prediction_status"],
                                "cache_hit": record["cache_hit"],
                            }
                            for record in records
                        ],
                    }
                    db.execute(
                        """UPDATE experiment_tasks SET status='COMPLETE',progress=1,
                           output_json=?,error_code=NULL,updated_at=? WHERE id=?""",
                        (json.dumps(output, ensure_ascii=False), now, task_id),
                    )
                    db.commit()
            except Exception:
                _discard_unregistered_artifact(target)
                raise
        return {"task_id": task_id, "artifact": artifact, "records": records}

    def list_paraplume_runs(self, project_id: str) -> list[dict[str, Any]]:
        self._row(project_id)
        with self._connect() as db:
            rows = db.execute(
                """SELECT r.*,a.relative_path,a.sha256,a.size_bytes
                   FROM paraplume_runs r JOIN artifacts a ON a.id=r.result_artifact_id
                   WHERE r.project_id=? ORDER BY r.created_at,r.id""",
                (project_id,),
            ).fetchall()
        folder = self.projects_root / project_id
        latest: dict[str, dict[str, Any]] = {}
        for row in rows:
            item = dict(row)
            path = (folder / item["relative_path"]).resolve()
            if not path.is_relative_to(folder.resolve()) or not path.is_file():
                raise IntegrityError("Paraplume result artifact is missing or unsafe")
            if path.stat().st_size != item["size_bytes"] or sha256_file(path) != item["sha256"]:
                raise IntegrityError("Paraplume result artifact hash differs")
            document = read_json(path)
            record = next(
                (
                    value
                    for value in document["records"]
                    if value["antibody_construct_id"] == item["antibody_construct_id"]
                ),
                None,
            )
            if record is None:
                raise IntegrityError("Paraplume artifact does not include its registered construct")
            latest[item["antibody_construct_id"]] = record | {
                "run_id": item["id"],
                "task_id": item["task_id"],
                "artifact_id": item["result_artifact_id"],
                "artifact_sha256": item["sha256"],
                "created_at": item["created_at"],
            }
        return list(latest.values())

    def create_constraint_export(
        self,
        project_id: str,
        *,
        antibody_construct_id: str,
        experimental_ir: list[dict[str, Any]],
        predicted_ir: list[dict[str, Any]],
        manual_rpr: list[dict[str, Any]],
        merge_sources: list[str],
        reviewer: str,
        expected_revision: int,
    ) -> dict[str, Any]:
        from .experimental_grasp import (
            ConstraintError,
            FastaChain,
            build_constraint_package,
            select_experimental_ir,
        )

        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        if manual_rpr:
            raise ConstraintError(
                "RPR_COMPLEX_STRUCTURE_REQUIRED: no validated antigen-antibody complex and residue mapping are registered",
                code="RPR_COMPLEX_STRUCTURE_REQUIRED",
            )
        if not reviewer.strip():
            raise ValueError("A reviewer is required")
        with self._connect() as db:
            construct = db.execute(
                "SELECT * FROM antibody_constructs WHERE project_id=? AND id=?",
                (project_id, antibody_construct_id),
            ).fetchone()
        if construct is None:
            raise KeyError(antibody_construct_id)
        # Enforce evidence links against immutable source records before serialization.
        with self._connect() as db:
            effects = {
                row["id"]: dict(row)
                for row in db.execute(
                    """SELECT e.id,e.result_artifact_id,a.relative_path,a.sha256
                       FROM effect_analyses e JOIN artifacts a ON a.id=e.result_artifact_id
                       WHERE e.project_id=? AND e.status='COMPLETE'""",
                    (project_id,),
                )
            }
        runs = {item["run_id"]: item for item in self.list_paraplume_runs(project_id)}
        for row in experimental_ir:
            effect = effects.get(row.get("evidence_id"))
            if effect is None:
                raise ValueError("EXPERIMENTAL_IR_EVIDENCE_REQUIRED: evidence_id must be a completed effect analysis")
            probe_id = row.get("details", {}).get("probe_id")
            if not probe_id:
                raise ValueError("EXPERIMENTAL_IR_PROBE_REQUIRED: details.probe_id must identify a source probe")
            effect_document = read_json(self.projects_root / project_id / effect["relative_path"])
            mapping_csv = self.projects_root / project_id / "reference" / "residue_mapping.csv"
            with mapping_csv.open(newline="", encoding="utf-8-sig") as handle:
                mapping_rows = list(csv.DictReader(handle))
            expected_constraint, trace = select_experimental_ir(
                effect_document,
                probe_id,
                mapping_rows,
                project["manifest"]["reference_sequence"],
            )
            if row.get("constraint") != expected_constraint:
                raise ValueError(
                    f"EXPERIMENTAL_IR_FASTA_MAPPING_MISMATCH: probe {probe_id} maps to {expected_constraint}"
                )
            row["details"] = dict(row.get("details", {})) | {
                **trace,
                "effect_artifact_sha256": effect["sha256"],
            }
            row["reviewer"] = reviewer
        for row in predicted_ir:
            run = runs.get(row.get("evidence_id"))
            if run is None or run["antibody_construct_id"] != antibody_construct_id:
                raise ValueError("PREDICTED_IR_EVIDENCE_REQUIRED: run must belong to the selected antibody")
            try:
                chain = int(row["constraint"].split("-")[0])
                position = int(row["constraint"].split("-")[1])
            except (KeyError, ValueError, IndexError) as exc:
                raise ConstraintError("Predicted IR must use chain-position-amino-acid syntax") from exc
            if chain not in {2, 3}:
                raise ValueError("Paraplume evidence only supports VH or VL chains")
            scores = run["vh_probabilities"] if chain == 2 else run["vl_probabilities"]
            if position < 1 or position > len(scores):
                raise ValueError("Predicted IR position is outside the selected Paraplume result")
            row["details"] = dict(row.get("details", {})) | {
                "paraplume_probability": scores[position - 1],
                "paraplume_run_id": run["run_id"],
                "prediction_status": "PREDICTED_UNREVIEWED",
            }
            row["reviewer"] = reviewer
        if not experimental_ir and not predicted_ir and not manual_rpr:
            raise ValueError("Select at least one reviewed constraint")
        artifact_id = uuid.uuid4().hex
        now = _now()
        relative = f"analyses/constraint_exports/{artifact_id}.zip"
        folder = self.projects_root / project_id
        target = folder / relative
        package, export_manifest = build_constraint_package(
            chains=[
                FastaChain(1, "ECD", project["manifest"]["reference_sequence"]),
                FastaChain(2, "VH", construct["vh_sequence"]),
                FastaChain(3, "VL", construct["vl_sequence"]),
            ],
            experimental_ir=experimental_ir,
            predicted_ir=predicted_ir,
            manual_rpr=manual_rpr,
            merge_sources=merge_sources,
            manifest_context={
                "limitations": [
                    "IR-only export: no validated antigen-antibody complex structure is registered; RPR is unavailable.",
                    "Predicted IR remains PREDICTED_UNREVIEWED and is not antigen epitope evidence.",
                ],
                "provenance": {
                    "project_id": project_id,
                    "project_revision": expected_revision,
                    "antibody_construct_id": antibody_construct_id,
                    "reviewer": reviewer,
                },
            },
        )
        target.parent.mkdir(parents=True, exist_ok=True)
        temp = target.with_name(target.name + f".{uuid.uuid4().hex}.tmp")
        temp.write_bytes(package)
        digest = sha256_file(temp)
        size_bytes = temp.stat().st_size
        os.replace(temp, target)
        target.chmod(0o444)
        artifact = {
            "id": artifact_id, "project_id": project_id, "relative_path": relative,
            "kind": "grasp_constraint_export", "sha256": digest,
            "size_bytes": size_bytes, "media_type": "application/zip", "read_only": 1,
            "created_at": now,
        }
        try:
            with self._lock, self._connect() as db:
                revision = self._append_revision(
                    db, project_id, expected_revision, "constraint_export", artifact_id,
                    {"artifact_sha256": digest, "constraint_count": export_manifest["constraint_count"],
                     "antibody_construct_id": antibody_construct_id},
                )
                db.execute("INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?)", (
                    artifact_id, project_id, relative, artifact["kind"], digest,
                    size_bytes, artifact["media_type"], 1, now,
                ))
                db.execute("""INSERT INTO constraint_exports
                    VALUES (?,?,?,?,?,?,?,?,?)""", (
                    artifact_id, project_id, antibody_construct_id, artifact_id,
                    json.dumps(merge_sources), export_manifest["constraint_count"],
                    export_manifest["status"], revision, now,
                ))
                manifest_path = folder / "project_manifest.json"
                manifest = read_json(manifest_path)
                manifest["artifacts"].append({
                    key: artifact[key] for key in (
                        "id", "kind", "relative_path", "sha256", "size_bytes", "media_type", "created_at"
                    )
                } | {"read_only": True})
                write_json(manifest_path, manifest)
                db.commit()
        except Exception:
            _discard_unregistered_artifact(target)
            raise
        return {
            "export_id": artifact_id, "artifact": artifact, "manifest": export_manifest,
            "project_revision": revision,
            "download_url": f"/api/v4/experiment-projects/{project_id}/files/{artifact_id}",
        }

    def list_constraint_exports(self, project_id: str) -> list[dict[str, Any]]:
        self._row(project_id)
        with self._connect() as db:
            rows = db.execute(
                """SELECT e.*,a.relative_path,a.sha256,a.size_bytes FROM constraint_exports e
                   JOIN artifacts a ON a.id=e.artifact_id WHERE e.project_id=? ORDER BY e.created_at,e.id""",
                (project_id,),
            ).fetchall()
        result = []
        for row in rows:
            item = dict(row)
            item["merge_sources"] = json.loads(item.pop("merge_sources_json"))
            target = (self.projects_root / project_id / item["relative_path"]).resolve()
            if not target.is_relative_to((self.projects_root / project_id).resolve()) or not target.is_file():
                raise IntegrityError("GRASP constraint export is missing or unsafe")
            if target.stat().st_size != item["size_bytes"] or sha256_file(target) != item["sha256"]:
                raise IntegrityError("GRASP constraint export hash differs")
            result.append(item)
        return result

    def get_paraplume_cache(self, cache_key: str) -> dict[str, Any] | None:
        with self._connect() as db:
            row = db.execute(
                "SELECT result_json,provenance_json FROM paraplume_cache WHERE cache_key=?",
                (cache_key,),
            ).fetchone()
        if row is None:
            return None
        return {"result": json.loads(row["result_json"]), "provenance": json.loads(row["provenance_json"])}

    def paraplume_construct_inputs(self, task_id: str) -> list[dict[str, str]]:
        task = self.get_experiment_task(task_id)
        if task["kind"] != "PARAPLUME":
            raise ValueError("Task is not a Paraplume task")
        ids = [item["antibody_construct_id"] for item in task["input"]["constructs"]]
        placeholders = ",".join("?" for _ in ids)
        with self._connect() as db:
            rows = db.execute(
                f"SELECT * FROM antibody_constructs WHERE project_id=? AND id IN ({placeholders})",
                [task["project_id"], *ids],
            ).fetchall()
        by_id = {row["id"]: row for row in rows}
        if len(by_id) != len(ids):
            raise KeyError("An antibody construct was removed before prediction")
        expected = {
            item["antibody_construct_id"]: item for item in task["input"]["constructs"]
        }
        result = []
        for item in task["input"]["constructs"]:
            row = by_id[item["antibody_construct_id"]]
            vh = row["vh_sequence"]
            vl = row["vl_sequence"]
            if (
                sha256_bytes(vh.encode("ascii")) != expected[item["antibody_construct_id"]]["vh_sha256"]
                or sha256_bytes(vl.encode("ascii")) != expected[item["antibody_construct_id"]]["vl_sha256"]
            ):
                raise IntegrityError("Antibody sequence changed after prediction task creation")
            result.append({"construct_id": item["antibody_construct_id"], "vh_sequence": vh, "vl_sequence": vl})
        return result

    def create_qc_review(
        self,
        project_id: str,
        *,
        source_type: str,
        source_id: str,
        sample_id: str,
        status: str,
        reason: str,
        reviewer: str,
        expected_revision: int,
    ) -> dict[str, Any]:
        if self.get(project_id)["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        allowed = {"PASS", "PASS_WITH_OVERRIDE", "HOLD", "EXCLUDE"}
        if status not in allowed:
            raise ValueError(f"Unsupported QC status: {status}")
        if status in {"PASS_WITH_OVERRIDE", "HOLD", "EXCLUDE"} and not reason.strip():
            raise ValueError(f"A reason is required for {status}")
        if not reviewer.strip():
            raise ValueError("QC reviewer is required")
        if source_type not in {"GATE_ANALYSIS", "FLOWJO_IMPORT", "FLOWJO_PERCENTAGE_IMPORT"}:
            raise ValueError(f"Unsupported QC source type: {source_type}")
        with self._lock, self._connect() as db:
            sample = db.execute(
                "SELECT artifact_id,batch_id,condition_json FROM samples WHERE id=? AND project_id=?",
                (sample_id, project_id),
            ).fetchone()
            if sample is None:
                raise KeyError(sample_id)
            if source_type == "GATE_ANALYSIS":
                source = db.execute(
                    "SELECT input_artifact_ids_json,result_artifact_id FROM gate_analyses "
                    "WHERE id=? AND project_id=? AND status='COMPLETE'",
                    (source_id, project_id),
                ).fetchone()
                if source is None:
                    raise KeyError(source_id)
                if sample["artifact_id"] not in json.loads(source["input_artifact_ids_json"]):
                    raise ValueError("Sample is not part of this gate analysis")
                artifact = db.execute(
                    "SELECT relative_path FROM artifacts WHERE id=? AND project_id=?",
                    (source["result_artifact_id"], project_id),
                ).fetchone()
                if artifact is None:
                    raise IntegrityError("Gate analysis artifact is missing")
                analysis_doc = read_json(self.projects_root / project_id / artifact["relative_path"])
                analysis_sample = next((item for item in analysis_doc.get("samples", [])
                                        if item.get("sample", {}).get("sample_id") == sample_id), None)
                review_flags = (analysis_sample or {}).get("qc", {}).get("findings", [])
            elif source_type == "FLOWJO_IMPORT":
                source = db.execute(
                    "SELECT preview_artifact_id,status FROM flowjo_imports WHERE id=? AND project_id=?",
                    (source_id, project_id),
                ).fetchone()
                if source is None:
                    raise KeyError(source_id)
                if source["status"] != "COMPLETE":
                    raise ValueError("FlowJo import is unresolved and cannot receive a QC review")
                preview_row = db.execute(
                    "SELECT relative_path FROM artifacts WHERE id=? AND project_id=?",
                    (source["preview_artifact_id"], project_id),
                ).fetchone()
                if preview_row is None:
                    raise IntegrityError("FlowJo preview artifact is missing")
                preview = read_json(self.projects_root / project_id / preview_row["relative_path"])
                condition = json.loads(sample["condition_json"])
                sample_well = _flowjo_well_id(condition.get("well_id"))
                if sample_well is None:
                    sample_well = _flowjo_well_id(condition.get("sample_label"))
                preview_samples = preview.get("records", {}).get("samples", [])
                matched_preview_samples = [
                    item for item in preview_samples
                    if _flowjo_well_id(item.get("sample_id")) == sample_well
                ]
                if sample_well is None or len(matched_preview_samples) != 1:
                    raise ValueError("Sample is not represented in this FlowJo import")
                preview_sample = matched_preview_samples[0]
                review_flags = (preview_sample or {}).get("qc", {}).get("findings", [])
            else:
                source = db.execute(
                    "SELECT batch_id,result_artifact_id,status FROM flowjo_percentage_imports "
                    "WHERE id=? AND project_id=?",
                    (source_id, project_id),
                ).fetchone()
                if source is None:
                    raise KeyError(source_id)
                if source["status"] != "IMPORTED":
                    raise ValueError("Percentage result import is not available for QC review")
                if source["batch_id"] != sample["batch_id"]:
                    raise ValueError("Sample is not part of this percentage result table batch")
                artifact = db.execute(
                    "SELECT relative_path,sha256,size_bytes FROM artifacts "
                    "WHERE id=? AND project_id=? AND kind='flowjo_percentage_result'",
                    (source["result_artifact_id"], project_id),
                ).fetchone()
                if artifact is None:
                    raise IntegrityError("Percentage result artifact is missing")
                result_path = self.projects_root / project_id / artifact["relative_path"]
                if (
                    result_path.is_symlink()
                    or not result_path.is_file()
                    or result_path.stat().st_size != artifact["size_bytes"]
                    or sha256_file(result_path) != artifact["sha256"]
                ):
                    raise IntegrityError("Percentage result artifact failed integrity verification")
                result_document = read_json(result_path)
                matched_observations = [
                    item for item in result_document.get("observations", [])
                    if item.get("sample_id") == sample_id
                ]
                if len(matched_observations) != 1:
                    raise ValueError("Sample is not uniquely represented in this percentage result table")
                review_flags = []
            if status == "PASS" and any(item.get("severity") == "REVIEW" for item in review_flags) and not reason.strip():
                raise ValueError("A reason is required when passing a sample with automated QC review flags")
            version = db.execute(
                "SELECT COALESCE(MAX(version),0)+1 FROM qc_reviews "
                "WHERE source_type=? AND source_id=? AND sample_id=?",
                (source_type, source_id, sample_id),
            ).fetchone()[0]
            review_id = uuid.uuid4().hex
            project_revision = self._append_revision(
                db, project_id, expected_revision, "qc_review", review_id,
                {"source_type": source_type, "source_id": source_id, "sample_id": sample_id,
                 "version": version, "status": status, "reason": reason, "reviewer": reviewer},
            )
            created_at = _now()
            db.execute(
                "INSERT INTO qc_reviews VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (review_id, project_id, source_type, source_id, sample_id, version, status,
                 reason, reviewer, project_revision, created_at),
            )
            db.commit()
        return {
            "qc_review_id": review_id, "source_type": source_type, "source_id": source_id,
            "sample_id": sample_id, "version": version, "status": status,
            "reason": reason, "reviewer": reviewer, "project_revision": project_revision,
            "created_at": created_at,
        }

    def build_flowjo_percentage_qc_worksheet(
        self, project_id: str, import_id: str
    ) -> bytes:
        from .flowjo_qc_worksheet import build_qc_review_worksheet

        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        imported = next(
            (item for item in self.list_flowjo_percentage_imports(project_id)
             if item["id"] == import_id),
            None,
        )
        if imported is None:
            raise KeyError(import_id)
        result = imported["result"]
        if imported["status"] != "IMPORTED":
            raise ValueError("Percentage result import is not available for QC review")
        latest_reviews: dict[str, dict[str, Any]] = {}
        for review in self.list_qc_reviews(
            project_id, "FLOWJO_PERCENTAGE_IMPORT", import_id
        ):
            old = latest_reviews.get(review["sample_id"])
            if old is None or review["version"] > old["version"]:
                latest_reviews[review["sample_id"]] = review

        source_sha256 = imported["source_artifact"]["sha256"]
        rows = []
        for observation in result.get("observations", []):
            previous = latest_reviews.get(observation["sample_id"], {})
            rows.append({
                "schema_version": "1.0",
                "project_id": project_id,
                "export_project_revision": project["revision"],
                "import_id": import_id,
                "batch_id": imported["batch_id"],
                "source_sha256": source_sha256,
                "sample_id": observation["sample_id"],
                "sample_artifact_id": observation.get("sample_artifact_id"),
                "well_id": observation.get("well_id"),
                "panel_label": observation.get("panel_label"),
                "condition_label": observation.get("condition_label"),
                "antigen_label": observation.get("antigen_label"),
                "replicate_type": observation.get("replicate_type"),
                "replicate_index": observation.get("replicate_index"),
                "displayed_value": observation.get("displayed_value"),
                "percentage": format(float(observation["percentage"]), ".15g"),
                "current_qc_review_id": previous.get("id"),
                "current_status": previous.get("status"),
                "current_reviewer": previous.get("reviewer"),
                "current_reason": previous.get("reason"),
                "decision": "",
                "reviewer": "",
                "reason": "",
            })
        return build_qc_review_worksheet(rows)

    def import_flowjo_percentage_qc_worksheet(
        self, project_id: str, import_id: str, payload: bytes
    ) -> dict[str, Any]:
        from .flowjo_qc_worksheet import (
            QC_WORKSHEET_STATUSES,
            parse_qc_review_worksheet,
        )

        project = self.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        rows = parse_qc_review_worksheet(payload)
        worksheet_sha256 = sha256_bytes(payload)

        with self._lock, self._connect() as db:
            # Lock the revision and QC ledger before validating the exported snapshot.
            db.execute("BEGIN IMMEDIATE")
            project_row = db.execute(
                "SELECT status,revision FROM experiment_projects WHERE id=?", (project_id,)
            ).fetchone()
            if project_row is None:
                raise KeyError(project_id)
            if project_row["status"] != "COMPLETE":
                raise IntegrityError("Project snapshot failed integrity verification")
            current_revision = project_row["revision"]
            imported = db.execute(
                "SELECT i.batch_id,i.result_artifact_id,i.source_artifact_id,i.status,"
                "result.relative_path AS result_path,result.sha256 AS result_sha256,"
                "result.size_bytes AS result_size,source.sha256 AS source_sha256 "
                "FROM flowjo_percentage_imports i "
                "JOIN artifacts result ON result.id=i.result_artifact_id AND result.project_id=i.project_id "
                "JOIN artifacts source ON source.id=i.source_artifact_id AND source.project_id=i.project_id "
                "WHERE i.id=? AND i.project_id=?",
                (import_id, project_id),
            ).fetchone()
            if imported is None:
                raise KeyError(import_id)
            if imported["status"] != "IMPORTED":
                raise ValueError("Percentage result import is not available for QC review")
            folder = self.projects_root / project_id
            result_path = folder / imported["result_path"]
            if (
                result_path.is_symlink()
                or not result_path.is_file()
                or result_path.stat().st_size != imported["result_size"]
                or sha256_file(result_path) != imported["result_sha256"]
            ):
                raise IntegrityError("Percentage result artifact failed integrity verification")
            result_document = read_json(result_path)
            source_artifact = db.execute(
                "SELECT relative_path,size_bytes,kind FROM artifacts WHERE id=? AND project_id=?",
                (imported["source_artifact_id"], project_id),
            ).fetchone()
            if source_artifact is None or source_artifact["kind"] != "raw_flowjo_result_table":
                raise IntegrityError("Percentage result source artifact is missing")
            source_path = folder / source_artifact["relative_path"]
            if (
                source_path.is_symlink()
                or not source_path.is_file()
                or source_path.stat().st_size != source_artifact["size_bytes"]
                or sha256_file(source_path) != imported["source_sha256"]
                or result_document.get("source_sha256") != imported["source_sha256"]
            ):
                raise IntegrityError("Percentage result source artifact failed integrity verification")

            # A retry after a lost HTTP response is idempotent for the exact uploaded bytes.
            for audit in db.execute(
                "SELECT details_json FROM audit_events WHERE project_id=? AND entity_type='qc_review'",
                (project_id,),
            ):
                details = json.loads(audit["details_json"])
                if (
                    details.get("qc_source_type") == "FLOWJO_PERCENTAGE_IMPORT"
                    and details.get("qc_source_id") == import_id
                    and details.get("worksheet_sha256") == worksheet_sha256
                ):
                    return {
                        "import_id": import_id,
                        "worksheet_sha256": worksheet_sha256,
                        "applied_count": 0,
                        "unchanged_count": len(rows),
                        "project_revision": current_revision,
                        "already_applied": True,
                    }

            observations = result_document.get("observations") or []
            by_sample = {item.get("sample_id"): item for item in observations}
            if len(by_sample) != len(observations) or not observations:
                raise IntegrityError("Percentage result observations are missing or duplicated")
            by_row = {row["sample_id"].strip(): row for row in rows}
            if len(by_row) != len(rows) or set(by_row) != set(by_sample):
                raise ValueError("QC worksheet must contain each imported observation exactly once")
            export_revisions = {row["export_project_revision"] for row in rows}
            if (
                len(export_revisions) != 1
                or not next(iter(export_revisions)).isdigit()
                or int(next(iter(export_revisions))) < 1
            ):
                raise ValueError("QC worksheet export revision metadata is invalid")
            export_revision = int(next(iter(export_revisions)))
            source_sha256 = imported["source_sha256"]
            batch_id = imported["batch_id"]

            review_rows = db.execute(
                "SELECT * FROM qc_reviews WHERE project_id=? AND source_type=? AND source_id=? "
                "ORDER BY sample_id,version",
                (project_id, "FLOWJO_PERCENTAGE_IMPORT", import_id),
            ).fetchall()
            latest_reviews: dict[str, dict[str, Any]] = {}
            for review in review_rows:
                old = latest_reviews.get(review["sample_id"])
                if old is None or review["version"] > old["version"]:
                    latest_reviews[review["sample_id"]] = dict(review)

            sample_ids = list(by_sample)
            placeholders = ",".join("?" for _ in sample_ids)
            sample_rows = db.execute(
                f"SELECT id,artifact_id,batch_id FROM samples WHERE project_id=? AND id IN ({placeholders})",
                (project_id, *sample_ids),
            ).fetchall()
            sample_by_id = {item["id"]: dict(item) for item in sample_rows}
            if set(sample_by_id) != set(by_sample):
                raise IntegrityError("Percentage result references missing samples")

            decisions: list[tuple[dict[str, Any], str, str, str, int]] = []
            unchanged_count = 0
            for sample_id, observation in by_sample.items():
                row = by_row[sample_id]
                sample = sample_by_id[sample_id]
                if sample["batch_id"] != batch_id or sample["artifact_id"] != observation.get("sample_artifact_id"):
                    raise IntegrityError("Percentage result sample mapping changed")
                previous = latest_reviews.get(sample_id, {})
                immutable_expected = {
                    "schema_version": "1.0",
                    "project_id": project_id,
                    "import_id": import_id,
                    "batch_id": batch_id,
                    "source_sha256": source_sha256,
                    "sample_id": sample_id,
                    "sample_artifact_id": sample["artifact_id"],
                    "well_id": observation.get("well_id"),
                    "panel_label": observation.get("panel_label"),
                    "condition_label": observation.get("condition_label"),
                    "antigen_label": observation.get("antigen_label"),
                    "replicate_type": observation.get("replicate_type"),
                    "replicate_index": observation.get("replicate_index"),
                    "displayed_value": observation.get("displayed_value"),
                    "percentage": format(float(observation["percentage"]), ".15g"),
                    "current_qc_review_id": previous.get("id"),
                    "current_status": previous.get("status"),
                    "current_reviewer": previous.get("reviewer"),
                    "current_reason": previous.get("reason"),
                }
                for key, expected in immutable_expected.items():
                    if row[key] != ("" if expected is None else str(expected)):
                        raise ValueError(f"QC worksheet source field was changed: {key} for {sample_id}")
                decision = row["decision"].strip()
                reviewer = row["reviewer"].strip()
                reason = row["reason"].strip()
                if not decision:
                    if reviewer or reason:
                        raise ValueError(f"Reviewer or reason is present without a decision for {sample_id}")
                    unchanged_count += 1
                    continue
                if decision not in QC_WORKSHEET_STATUSES:
                    raise ValueError(f"Unsupported QC decision for {sample_id}: {decision}")
                if not reviewer:
                    raise ValueError(f"A reviewer is required for {sample_id}")
                if len(reviewer) > 120 or len(reason) > 1000:
                    raise ValueError(f"Reviewer or reason is too long for {sample_id}")
                if decision in {"PASS_WITH_OVERRIDE", "HOLD", "EXCLUDE"} and not reason:
                    raise ValueError(f"A reason is required for {decision} ({sample_id})")
                if previous and (
                    decision == previous["status"]
                    and reviewer == previous["reviewer"]
                    and reason == previous["reason"]
                ):
                    unchanged_count += 1
                    continue
                decisions.append((observation, decision, reviewer, reason, previous.get("version", 0) + 1))

            revision = current_revision
            now = _now()
            for observation, decision, reviewer, reason, version in decisions:
                sample_id = observation["sample_id"]
                review_id = uuid.uuid4().hex
                revision = self._append_revision(
                    db, project_id, revision, "qc_review", review_id,
                    {
                        "source_type": "FLOWJO_PERCENTAGE_IMPORT",
                        "source_id": import_id,
                        "source_sha256": source_sha256,
                        "sample_id": sample_id,
                        "well_id": observation.get("well_id"),
                        "version": version,
                        "status": decision,
                        "reason": reason,
                        "reviewer": reviewer,
                        "worksheet_sha256": worksheet_sha256,
                        "worksheet_export_revision": export_revision,
                        "reviewed_via": "XLSX_WORKSHEET",
                    },
                )
                db.execute(
                    "INSERT INTO qc_reviews VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    (review_id, project_id, "FLOWJO_PERCENTAGE_IMPORT", import_id,
                     sample_id, version, decision, reason, reviewer, revision, now),
                )
            db.commit()
        return {
            "import_id": import_id,
            "worksheet_sha256": worksheet_sha256,
            "applied_count": len(decisions),
            "unchanged_count": unchanged_count,
            "project_revision": revision,
            "already_applied": False,
        }

    def list_qc_reviews(
        self, project_id: str, source_type: str | None = None, source_id: str | None = None
    ) -> list[dict[str, Any]]:
        self._row(project_id)
        with self._connect() as db:
            if source_type and source_id:
                rows = db.execute(
                    "SELECT * FROM qc_reviews WHERE project_id=? AND source_type=? AND source_id=? "
                    "ORDER BY sample_id,version",
                    (project_id, source_type, source_id),
                ).fetchall()
            else:
                rows = db.execute(
                    "SELECT * FROM qc_reviews WHERE project_id=? ORDER BY created_at,id",
                    (project_id,),
                ).fetchall()
        return [dict(row) for row in rows]
