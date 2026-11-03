"""Run or inspect the one-case ordinary site-first v2 acceptance job.

Run from the repository parent with the project venv. Model execution only
starts with the explicit ``--execute`` option.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "workspace_data"
ACCEPTANCE_DIR = DATA / "site_first_acceptance_20261027"
SOURCE_REQUEST = DATA / "2b522bedd0f5" / "request.json"
SOURCE_STRUCTURE = DATA / "2b522bedd0f5" / "uploads" / "structure.pdb"
SOURCE_RUN = DATA / "f606506c8e4a"
EXPECTED_STRUCTURE_SHA256 = "8ff23af24604861cec56365dc835068f97fa714222b8313a71291da2dc7851e4"
EXPECTED_PROTOCOL = "ordinary_site_first_v2"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    temporary.replace(path)


def preflight(source_request: Path = SOURCE_REQUEST, structure: Path = SOURCE_STRUCTURE,
              acceptance_dir: Path = ACCEPTANCE_DIR) -> dict:
    request = json.loads(Path(source_request).read_text(encoding="utf-8"))
    digest = sha256(Path(structure))
    if digest != EXPECTED_STRUCTURE_SHA256:
        raise ValueError(f"AIB234 source structure hash mismatch: {digest}")
    for field in ("input_type", "identifier", "chain", "analysis_target", "design_engine", "ml_mode"):
        if request.get(field) != {
            "input_type": "structure", "identifier": "IL7Ra_binder_AF2_35_AIB234.pdb", "chain": "A",
            "analysis_target": "binder", "design_engine": "joint_v4_ml", "ml_mode": "strict",
        }[field]:
            raise ValueError(f"Saved AIB234 request has unexpected {field}: {request.get(field)!r}")
    if request.get("confidence_path"):
        raise ValueError("AIB234 saved input is expected to have no confidence file")
    historical_trace = SOURCE_RUN / "execution_trace.jsonl"
    if not historical_trace.is_file():
        raise FileNotFoundError(f"Historical AIB234 trace is missing: {historical_trace}")
    state_path = Path(acceptance_dir) / "state.json"
    report_path = Path(acceptance_dir) / "report.json"
    return {
        "status": "PREFLIGHT_PASS", "mode": "read_only", "job_created": False,
        "structure_path": str(Path(structure).resolve()), "structure_sha256": digest,
        "source_request_path": str(Path(source_request).resolve()), "source_run": str(SOURCE_RUN),
        "protocol": EXPECTED_PROTOCOL, "analysis_target": "binder", "chain": "A",
        "confidence_path": None, "execution_already_started": state_path.exists() or report_path.exists(),
        "historical_trace_path": str(historical_trace.resolve()),
        "historical_trace_sha256": sha256(historical_trace),
        "output_directory": str(Path(acceptance_dir).resolve()),
    }


def request_for(source: dict, structure: Path, workdir: Path):
    if str(ROOT.parent) not in sys.path:
        sys.path.insert(0, str(ROOT.parent))
    from ecd_epitope_scanner.pipeline import PipelineRequest

    values = {key: source.get(key) for key in (
        "input_type", "identifier", "chain", "strict_af3", "allow_partial", "minimum_regions",
        "maximum_regions", "assay_format", "assay_profile", "design_engine", "ml_mode",
        "surface_strategy", "max_sectors", "target_sector_residues", "maximum_sector_geodesic_diameter_a",
        "maximum_sector_ca_diameter_a", "max_primary_constructs", "minimum_constructs", "target_constructs",
        "maximum_constructs", "adaptive_budget", "allow_spurs_only_final", "biological_context_mode",
        "optimization_profile", "candidate_ids", "calibration_repeat_schedule", "export_mode",
        "model_assets_root", "case_time_budget_seconds", "case_byte_budget", "analysis_target",
    ) if key in source}
    values.update(structure_path=str(structure.resolve()), confidence_path=None,
        output_dir=str(workdir.resolve()), metadata={
            "joint_v4_ml": {
                "execution_protocol": {"name": EXPECTED_PROTOCOL},
                "runtime": {"cache_dir": str((workdir / "site_first_score_cache").resolve())},
            },
            "benchmark_cache_condition": "cold",
        })
    return PipelineRequest(**values)


def assert_no_active_jobs(db_path: Path) -> None:
    if not Path(db_path).exists():
        return
    with sqlite3.connect(db_path) as connection:
        active = connection.execute("SELECT id,status FROM jobs WHERE status IN ('RUNNING','QUEUED')").fetchall()
    if active:
        raise RuntimeError(f"Active jobs prevent site-first acceptance: {active}")


def _cache_paths(workdir: Path) -> dict[str, Path]:
    workdir = Path(workdir).resolve()
    repair_work_root = workdir / "fx_shared" / "work"
    # engine.predict_site_first_foldx passes repair_work_root explicitly; the worker
    # places its SHA-keyed repair cache beside that directory, not beside engine_work_root.
    return {"score_cache": workdir / "site_first_score_cache",
            "repair_work_root": repair_work_root,
            "foldx_repair_cache": repair_work_root.parent / "foldx_repair_cache"}


def _cache_inventory(path: Path) -> dict:
    path = Path(path)
    files = [item for item in path.rglob("*") if item.is_file()] if path.exists() else []
    return {"path": str(path.resolve()), "exists": path.exists(), "file_count": len(files),
            "bytes": sum(item.stat().st_size for item in files)}


def _cache_evidence(workdir: Path) -> dict:
    paths = _cache_paths(workdir)
    inventory = {name: _cache_inventory(path) for name, path in paths.items() if name != "repair_work_root"}
    inventory["repair_work_root"] = _cache_inventory(paths["repair_work_root"])
    inventory["cold_verified"] = not inventory["score_cache"]["exists"] and not inventory["foldx_repair_cache"]["exists"]
    inventory["foldx_repair_cache_resolution"] = (
        "repair_work_root.parent / foldx_repair_cache; worker receives repair_work_root=output_dir/fx_shared/work"
    )
    return inventory


def _model_report(workdir: Path) -> dict:
    def read(name):
        path = workdir / name
        return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    ledger, performance = read("ordinary_candidate_ledger.json"), read("site_first_performance.json")
    foldx = ledger.get("foldx", performance.get("foldx", {}))
    return {
        "site": ledger.get("counts", performance.get("counts", {})), "candidate_timings": performance.get("candidate_timings", {}),
        "spurs": ledger.get("spurs", performance.get("spurs", {})),
        "foldx": foldx,
        "repair": {key: value for key, value in foldx.items() if "repair" in key.lower()},
        "worker_metrics": ledger.get("worker_metrics", performance.get("worker_metrics", {})),
        "peak_process_memory_bytes": performance.get("peak_process_memory_bytes"),
        "peak_memory_scope": performance.get("peak_memory_scope"),
        "disk_bytes_before_receipt": performance.get("disk_bytes_before_receipt"),
    }


def _region_principles(workdir: Path, analysis: dict) -> list[dict]:
    path = workdir / "sector_recommendations_v4.json"
    rows = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else []
    principles = ("epitope_support", "mechanism_pareto", "sector_coverage")
    by_region = {}
    for primary in rows if isinstance(rows, list) else []:
        if not isinstance(primary, dict):
            continue
        region = str(primary.get("sector_id", primary.get("region_id", "UNKNOWN")))
        entry = by_region.setdefault(region, {key: {"status": "NOT_REPORTED", "gap": None}
                                              for key in principles})
        unavailable = {item.get("principle"): item for item in primary.get("unavailable_principles", [])
                       if isinstance(item, dict) and item.get("principle")}
        cards = [primary, *(primary.get("alternatives") or [])]
        for card in cards:
            if not isinstance(card, dict):
                continue
            for principle in card.get("principles", []):
                if principle in entry:
                    entry[principle] = {"status": card.get("status", card.get("qualification_status")),
                                        "formal_eligible": card.get("formal_eligible"),
                                        "candidate_id": card.get("candidate_id"), "gap": None}
        for principle, gap in unavailable.items():
            if principle in entry and entry[principle]["status"] == "NOT_REPORTED":
                entry[principle] = {"status": gap.get("state", "UNAVAILABLE"),
                                    "candidate_id": gap.get("candidate_id"),
                                    "gap": gap.get("reason", gap.get("state"))}
    regions = analysis.get("regions", [])
    ids = []
    for index, row in enumerate(regions):
        if not isinstance(row, dict):
            continue
        region = str(row.get("sector_id", row.get("region_id", row.get("id", index))))
        if region.isdecimal() and f"VS{int(region):03d}" in by_region:
            region = f"VS{int(region):03d}"
        ids.append(region)
    ids.extend(region for region in by_region if region not in ids)
    return [{"region": region, "principles": by_region.get(
                region, {key: {"status": "NOT_REPORTED", "gap": None} for key in principles}),
             "missing_principles": [p for p in principles
                                      if str(by_region.get(region, {}).get(p, {}).get("status", "NOT_REPORTED"))
                                      .upper() in {"NOT_REPORTED", "TECHNICAL_FAILURE"}]}
            for region in ids]


def _independent_geometry_evidence(workdir: Path, verifier=None) -> dict:
    manifest = json.loads((workdir / "run_manifest.json").read_text(encoding="utf-8"))
    configured = manifest.get("parameters", {}).get("joint_v4_ml", {}).get("sectoring", {})
    limits = (configured.get("maximum_geodesic_diameter_a"),
              configured.get("maximum_euclidean_diameter_a"))
    if limits != (30.0, 22.0):
        raise ValueError(f"Independent geometry verifier is pinned to 30A/22A; run limits are {limits!r}")
    if verifier is None:
        scripts_dir = str(ROOT / "scripts")
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        source = DATA / "af3_binder_batch_20261017" / "verify_af3_batch.py"
        import importlib.util
        spec = importlib.util.spec_from_file_location("af3_single_geometry_verifier", source)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Could not load independent geometry verifier: {source}")
        verifier = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(verifier)
    evidence = verifier._geometry_evidence(Path(workdir))
    return {**evidence, "status": "PASS", "connectivity": "independently verified",
            "maximum_geodesic_diameter_a": limits[0], "maximum_ca_diameter_a": limits[1],
            "verifier": "workspace_data/af3_binder_batch_20261017/verify_af3_batch.py::_geometry_evidence"}


def _write_report(directory: Path, report: dict) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    atomic_json(directory / "report.json", report)
    comparison = report.get("historical_reference", {})
    markdown = ["# AIB234 site-first v2 single-case acceptance", "",
        f"Status: `{report.get('status')}`; job: `{report.get('job_id') or 'not created'}`.",
        f"Protocol: `{report.get('protocol')}`. Runtime comparison: **{comparison.get('runtime_conclusion', 'not established')}**.", "",
        "## Evidence", "", f"- Input SHA256: `{report.get('input_sha256')}`",
        f"- Source digest checkpoint/run/end: `{report.get('source_digest_checkpoint')}` / `{report.get('source_digest_run')}` / `{report.get('source_digest_end')}`",
        f"- Elapsed wall seconds: `{report.get('elapsed_seconds')}`",
        f"- Validation: `{(report.get('validation') or {}).get('status', 'NOT_RUN')}`", "",
        f"- Independent display geometry: `{(report.get('geometry_validation') or {}).get('status', 'NOT_RUN')}`", "",
        f"- Cold-cache preflight: `{(report.get('cache_evidence') or {}).get('before_run', {}).get('cold_verified', 'NOT_RUN')}`", "",
        "## Historical AIB234 reference", "",
        f"Historical trace: `{comparison.get('trace_path')}`; reported legacy runtime: {comparison.get('legacy_minutes')} minutes.",
        f"Cache condition: {comparison.get('cache_condition')}; protocol match: {comparison.get('protocol_match')}.",
        "50% reduction: **not determined**. Historical runtime is context only and is not a matched benchmark.", "",
        "## Model and export evidence", "", "```json", json.dumps(report.get("metrics", {}), ensure_ascii=False, indent=2), "```", "",
        "## Region principle coverage", "", "```json", json.dumps(report.get("region_principles", []), ensure_ascii=False, indent=2), "```", "",
        "## Errors", "", "```json", json.dumps(report.get("errors", []), ensure_ascii=False, indent=2), "```", ""]
    (directory / "report.md").write_text("\n".join(markdown), encoding="utf-8")


def execute(*, acceptance_dir: Path = ACCEPTANCE_DIR, source_request: Path = SOURCE_REQUEST,
            structure: Path = SOURCE_STRUCTURE, store_factory=None, runner_factory=None,
            validate_output=None, source_digest=None, data_dir: Path = DATA,
            geometry_check=None) -> dict:
    plan = preflight(source_request, structure, acceptance_dir)
    if plan["execution_already_started"]:
        raise FileExistsError(f"Single-case execution already has a state or report: {acceptance_dir}")
    acceptance_dir = Path(acceptance_dir)
    acceptance_dir.mkdir(parents=True, exist_ok=False)
    state_path = acceptance_dir / "state.json"
    state = {"status": "STARTING", "protocol": EXPECTED_PROTOCOL, "started_at_utc": datetime.now(timezone.utc).isoformat(),
             "input_path": plan["structure_path"], "input_sha256": plan["structure_sha256"], "job_id": None,
             "latest_event": None}
    atomic_json(state_path, state)
    report = {"status": "NEEDS_REVIEW", "protocol": EXPECTED_PROTOCOL, "input_path": plan["structure_path"],
              "input_sha256": plan["structure_sha256"], "historical_reference": {
                  "trace_path": plan["historical_trace_path"],
                  "trace_sha256": plan["historical_trace_sha256"], "legacy_minutes": 31.8,
                  "cache_condition": "UNKNOWN", "protocol_match": False,
                  "runtime_conclusion": "NOT_DETERMINED; historical runtime is not a matched benchmark"},
              "errors": [], "validation": None}
    job_id = uuid.uuid4().hex[:12]
    report["job_id"] = job_id
    workdir = Path(data_dir) / job_id
    runner = None
    store = None
    elapsed = None
    start_digest = end_digest = None
    run_digest = None
    original_record = None
    try:
        assert_no_active_jobs(Path(data_dir) / "jobs.sqlite3")
        if workdir.exists():
            raise FileExistsError(f"Refusing to reuse job directory: {workdir}")
        upload = workdir / "uploads"
        upload.mkdir(parents=True, exist_ok=False)
        local_structure = upload / "structure.pdb"
        local_structure.write_bytes(Path(structure).read_bytes())
        if sha256(local_structure) != EXPECTED_STRUCTURE_SHA256:
            raise ValueError("Copied AIB234 input hash mismatch")
        before_cache = _cache_evidence(workdir)
        report["cache_evidence"] = {"before_run": before_cache, "after_run": None}
        state["cache_evidence"] = before_cache
        atomic_json(state_path, state)
        if not before_cache["cold_verified"]:
            raise FileExistsError("New job score cache or FoldX repair cache already exists")
        if source_digest is None:
            if str(ROOT.parent) not in sys.path:
                sys.path.insert(0, str(ROOT.parent))
            from ecd_epitope_scanner.joint_design_v4.run_audit import source_snapshot
            source_digest = lambda: source_snapshot(ROOT)["source_files_sha256"]
        start_digest = source_digest()
        source = json.loads(Path(source_request).read_text(encoding="utf-8"))
        request = request_for(source, local_structure, workdir)
        if store_factory is None:
            if str(ROOT.parent) not in sys.path:
                sys.path.insert(0, str(ROOT.parent))
            from ecd_epitope_scanner.jobs import JobRunner, JobStore
            store_factory = lambda: JobStore(DATA)
            runner_factory = runner_factory or JobRunner
        assert_no_active_jobs(Path(data_dir) / "jobs.sqlite3")
        store = store_factory()
        job = store.create(request, job_id=job_id)
        state.update(status="RUNNING", job_id=job_id, workdir=job["workdir"], source_digest_start=start_digest)
        atomic_json(state_path, state)
        runner = runner_factory(store)

        if original_record is None:
            if str(ROOT.parent) not in sys.path:
                sys.path.insert(0, str(ROOT.parent))
            from ecd_epitope_scanner.joint_design_v4.run_audit import ExecutionTrace
            original_record = ExecutionTrace.record
            def record_with_state(trace, event, stage, fraction, message, **details):
                row = original_record(trace, event, stage, fraction, message, **details)
                if Path(trace.path).resolve().parent == Path(job["workdir"]).resolve():
                    state["latest_event"] = row
                    atomic_json(state_path, state)
                return row
            ExecutionTrace.record = record_with_state
        started = time.monotonic()
        runner._run(job_id)
        elapsed = time.monotonic() - started
        saved_job = store.get(job_id)
        report.update(job_status=saved_job["status"], job_error=saved_job.get("error"),
                      workdir=saved_job["workdir"], elapsed_seconds=elapsed)
        workdir = Path(saved_job["workdir"])
        snapshot_path = workdir / "run_source_snapshot_v1.json"
        if snapshot_path.is_file():
            run_digest = json.loads(snapshot_path.read_text(encoding="utf-8")).get("source_files_sha256")
        else:
            report["errors"].append("Run source snapshot is missing")
        end_digest = source_digest()
        report["source_digest_checkpoint"] = start_digest
        report["source_digest_run"] = run_digest
        report["source_digest_end"] = end_digest
        if run_digest != start_digest or end_digest != start_digest:
            report["errors"].append("Production source digest differs from the pre-run checkpoint or changed during the run")
        if saved_job["status"] != "COMPLETE":
            report["errors"].append(f"Job ended {saved_job['status']}: {saved_job.get('error')}")
        else:
            analysis_path = workdir / "analysis.json"
            analysis = json.loads(analysis_path.read_text(encoding="utf-8")) if analysis_path.is_file() else {}
            if validate_output is None:
                if str(ROOT.parent) not in sys.path:
                    sys.path.insert(0, str(ROOT.parent))
                sys.path.insert(0, str(ROOT / "scripts"))
                from run_ordinary_adaptive_rerun import validate_case_output
                validate_output = validate_case_output
            validation = validate_output(workdir, analysis, expected_protocol=EXPECTED_PROTOCOL)
            try:
                if geometry_check is None:
                    geometry_check = _independent_geometry_evidence
                geometry = geometry_check(workdir)
                if geometry.get("status") != "PASS":
                    raise ValueError(f"Independent geometry verifier returned {geometry.get('status')!r}")
            except Exception as geometry_exc:
                geometry = {"status": "NEEDS_REVIEW", "error": f"{type(geometry_exc).__name__}: {geometry_exc}"}
                validation.setdefault("errors", []).append("Independent connectivity/diameter verification failed: " + geometry["error"])
                validation["status"] = "NEEDS_REVIEW"
            validation.setdefault("checks", {})["independent_geometry"] = geometry
            report["geometry_validation"] = geometry
            recs = json.loads((workdir / "sector_recommendations_v4.json").read_text(encoding="utf-8")) if (workdir / "sector_recommendations_v4.json").is_file() else []
            formal_cards = [row for row in recs if isinstance(row, dict) and row.get("candidate_id")
                            and row.get("formal_eligible") is True and row.get("status") == "PASSED"]
            validation.setdefault("checks", {})["formal_cards_present"] = bool(formal_cards)
            if not formal_cards:
                validation.setdefault("errors", []).append("No formal recommendation cards were produced")
                validation["status"] = "NEEDS_REVIEW"
            report["validation"] = validation
            if validation.get("status") == "PASS" and not report["errors"]:
                report["status"] = "PASS"
            else:
                report["errors"].extend(validation.get("errors", []))
        report["metrics"] = _model_report(workdir)
        report.setdefault("cache_evidence", {"before_run": before_cache if 'before_cache' in locals() else None,
                                               "after_run": None})
        if 'before_cache' in locals():
            report["cache_evidence"]["after_run"] = _cache_evidence(workdir)
        report["region_principles"] = _region_principles(workdir, analysis if saved_job["status"] == "COMPLETE" else {})
        report["source_snapshot_path"] = str(workdir / "run_source_snapshot_v1.json")
        report["performance_receipt_path"] = str(workdir / "site_first_performance.json")
    except Exception as exc:
        report["status"] = "NEEDS_REVIEW"
        report["errors"].append(f"{type(exc).__name__}: {exc}")
        if store is not None:
            try:
                current = store.get(job_id)
                if current["status"] in {"QUEUED", "RUNNING"}:
                    store.update(job_id, status="FAILED", stage="FAILED",
                                 message="Single-case acceptance failed before terminal state",
                                 error=report["errors"][-1])
            except Exception as status_exc:
                report["errors"].append(f"Could not finalize job state: {status_exc}")
        if elapsed is not None:
            report["elapsed_seconds"] = elapsed
    finally:
        if original_record is not None:
            ExecutionTrace.record = original_record
        if runner is not None:
            runner.executor.shutdown(wait=False, cancel_futures=True)
        report.setdefault("source_digest_start", start_digest)
        report.setdefault("source_digest_checkpoint", start_digest)
        report.setdefault("source_digest_run", run_digest)
        report.setdefault("source_digest_end", end_digest)
        report.setdefault("elapsed_seconds", elapsed)
        state.update(status=report["status"], finished_at_utc=datetime.now(timezone.utc).isoformat(),
                     job_id=report.get("job_id"), report_path=str(acceptance_dir / "report.json"))
        atomic_json(state_path, state)
        _write_report(acceptance_dir, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true", help="Read-only input and duplicate-run check")
    mode.add_argument("--execute", action="store_true", help="Run one real AIB234 v2 job synchronously")
    args = parser.parse_args()
    if args.preflight:
        print(json.dumps(preflight(), indent=2, ensure_ascii=False))
    else:
        report = execute()
        print(json.dumps({"status": report["status"], "job_id": report.get("job_id"),
                          "report": str(ACCEPTANCE_DIR / "report.json")}, indent=2))
        if report["status"] != "PASS":
            raise SystemExit(1)


if __name__ == "__main__":
    main()
