"""Prepare or explicitly execute a cold legacy baseline for the AIB234 case.

Run from the repository parent with the project venv. The default command is
read-only; the archived pipeline runs only with ``--execute``.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
ACCEPTANCE_DIR = ROOT / "workspace_data/site_first_acceptance_20261027"
BASELINE_ROOT = ACCEPTANCE_DIR / "baseline"
SNAPSHOT_ROOT = BASELINE_ROOT / "baseline_snapshot"
SNAPSHOT_PACKAGE = SNAPSHOT_ROOT / "ecd_epitope_scanner"
REPORT_JSON = ACCEPTANCE_DIR / "baseline_report.json"
REPORT_MD = ACCEPTANCE_DIR / "baseline_report.md"
SOURCE_REQUEST = ROOT / "workspace_data/2b522bedd0f5/request.json"
SOURCE_STRUCTURE = ROOT / "workspace_data/2b522bedd0f5/uploads/structure.pdb"
EXPECTED_INPUT_SHA256 = "8ff23af24604861cec56365dc835068f97fa714222b8313a71291da2dc7851e4"
BASELINE_COMMIT = "1fe65072af0fdc80b9b3e79e24f2dd1b1813405"
FOLDX_RELATIVE = Path("vendor/foldx/foldx_1_20270131.exe")


def sha256_file(path: Path) -> str:
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


def resolve_commit(root: Path = ROOT, commit: str = BASELINE_COMMIT) -> str:
    resolved = subprocess.run(["git", "rev-parse", "--verify", f"{commit}^{{commit}}"], cwd=root,
                              check=True, capture_output=True, text=True).stdout.strip()
    if not resolved.startswith(commit):
        raise ValueError(f"Resolved baseline commit does not match requested revision: {resolved}")
    return resolved


def safe_extract_archive(archive: bytes, destination: Path) -> Path:
    """Extract a git tar after validating every member and refusing links."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    root = destination.resolve()
    seen = set()
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:*") as bundle:
        for member in bundle.getmembers():
            name = PurePosixPath(member.name)
            if name.is_absolute() or ".." in name.parts or not name.parts or ":" in name.parts[0]:
                raise ValueError(f"Unsafe archive member: {member.name}")
            target = destination.joinpath(*name.parts)
            if not target.resolve().is_relative_to(root):
                raise ValueError(f"Archive member escapes destination: {member.name}")
            normalized = name.as_posix().rstrip("/")
            if normalized in seen:
                raise ValueError(f"Duplicate archive member: {member.name}")
            seen.add(normalized)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                source = bundle.extractfile(member)
                if source is None:
                    raise ValueError(f"Unreadable archive member: {member.name}")
                with source, target.open("xb") as output:
                    while block := source.read(1024 * 1024):
                        output.write(block)
            else:
                raise ValueError(f"Archive links and special files are not accepted: {member.name}")
    package = destination / "ecd_epitope_scanner"
    if not (package / "pipeline.py").is_file():
        raise ValueError("Archive did not contain the expected ecd_epitope_scanner package")
    return package


def model_identity_paths(root: Path, foldx_binary: Path) -> dict[str, Path]:
    return {
        "foldx_binary": foldx_binary,
        "spurs_single_checkpoint": root / "models/spurs/spurs/checkpoints/best.ckpt",
        "spurs_multi_checkpoint": root / "models/spurs/spurs_multi/checkpoints/best.ckpt",
        "spurs_esm_weights": root / "models/spurs/esm2_t33_650M_UR50D.pt",
        "discotope_model_archive": root / "vendor/discotope3/models.zip",
    }


def validate_model_paths(root: Path = ROOT, foldx_binary: Path | None = None) -> dict[str, str]:
    binary = Path(foldx_binary or (root / FOLDX_RELATIVE)).resolve()
    paths = model_identity_paths(Path(root).resolve(), binary)
    missing = [f"{name}: {path}" for name, path in paths.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required legacy model identity files are missing: " + "; ".join(missing))
    return {name: str(path.resolve()) for name, path in paths.items()}


def preflight(*, root: Path = ROOT, acceptance_dir: Path = ACCEPTANCE_DIR,
              source_request: Path = SOURCE_REQUEST, structure: Path = SOURCE_STRUCTURE,
              foldx_binary: Path | None = None, model_paths_override: dict[str, Path] | None = None,
              commit_resolver=resolve_commit, max_foldx_root_length: int = 100) -> dict:
    root, acceptance_dir = Path(root).resolve(), Path(acceptance_dir).resolve()
    request = json.loads(Path(source_request).read_text(encoding="utf-8"))
    input_hash = sha256_file(Path(structure))
    if input_hash != EXPECTED_INPUT_SHA256:
        raise ValueError(f"AIB234 input hash mismatch: {input_hash}")
    expected = {"input_type": "structure", "identifier": "IL7Ra_binder_AF2_35_AIB234.pdb",
                "chain": "A", "analysis_target": "binder", "design_engine": "joint_v4_ml",
                "ml_mode": "strict"}
    for key, value in expected.items():
        if request.get(key) != value:
            raise ValueError(f"Saved AIB234 request has unexpected {key}: {request.get(key)!r}")
    if request.get("confidence_path"):
        raise ValueError("AIB234 legacy replay must use the saved no-confidence input")
    output_paths = (acceptance_dir / "baseline", acceptance_dir / "baseline_report.json",
                    acceptance_dir / "baseline_report.md")
    existing = [str(path) for path in output_paths if path.exists()]
    if existing:
        raise FileExistsError("Refusing to reuse existing legacy baseline output: " + ", ".join(existing))
    binary = Path(foldx_binary or root / FOLDX_RELATIVE).resolve()
    if model_paths_override is None:
        identities = validate_model_paths(root, binary)
    else:
        missing = [str(path) for path in model_paths_override.values() if not Path(path).is_file()]
        if missing:
            raise FileNotFoundError("Required legacy model identity files are missing: " + "; ".join(missing))
        identities = {name: str(Path(path).resolve()) for name, path in model_paths_override.items()}
        if identities.get("foldx_binary") != str(binary):
            raise ValueError("Model identity override must identify the configured FoldX binary")
    resolved = commit_resolver(root)
    short_root = acceptance_dir / "baseline" / "fx"
    if len(str(short_root.resolve())) >= max_foldx_root_length:
        raise ValueError(f"FoldX short work root must be under 100 characters: {short_root}")
    return {"status": "PREFLIGHT_PASS", "execution_started": False, "baseline_commit": resolved,
            "input_path": str(Path(structure).resolve()), "input_sha256": input_hash,
            "source_request": str(Path(source_request).resolve()), "protocol": "legacy",
            "model_identity_paths": identities, "model_assets_root": str(root),
            "output_root": str(acceptance_dir / "baseline"), "short_foldx_work_root": str(short_root.resolve()),
            "outputs_available": not any(path.exists() for path in output_paths)}


def build_legacy_request(source: dict, *, structure: Path, output: Path, baseline_root: Path,
                         model_assets_root: Path, foldx_binary: Path) -> dict:
    request = dict(source)
    metadata = dict(request.get("metadata") or {})
    metadata.pop("execution_queue", None)
    metadata.pop("execution_protocol", None)
    joint = dict(metadata.get("joint_v4_ml") or {})
    joint.pop("execution_queue", None)
    joint.pop("execution_protocol", None)
    runtime = dict(joint.get("runtime") or {})
    score_cache = baseline_root / "score_cache"
    runtime["cache_dir"] = str(score_cache.resolve())
    joint["runtime"] = runtime
    foldx = dict(joint.get("foldx_confirmation") or {})
    foldx.update(binary=str(foldx_binary.resolve()),
                 engine_work_root=str((baseline_root / "fx").resolve()), short_binary_alias=True)
    joint["foldx_confirmation"] = foldx
    metadata["joint_v4_ml"] = joint
    request.update(structure_path=str(structure.resolve()), confidence_path=None,
                   output_dir=str(output.resolve()), model_assets_root=str(model_assets_root.resolve()),
                   metadata=metadata)
    return request


CHILD_CODE = r'''import json, sys
from pathlib import Path
snapshot_parent, request_path = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(snapshot_parent))
from ecd_epitope_scanner.pipeline import PipelineRequest, run_pipeline
data = json.loads(request_path.read_text(encoding="utf-8"))
request = PipelineRequest(**data)
def progress(*args, **kwargs):
    print(json.dumps({"event": "progress", "args": args, "kwargs": kwargs}, default=str), flush=True)
try:
    print(json.dumps({"event": "pipeline_start"}), flush=True)
    result = run_pipeline(request, progress)
    print(json.dumps({"event": "result", "result_keys": sorted(result) if isinstance(result, dict) else [],
                      "result_type": type(result).__name__}), flush=True)
except BaseException as exc:
    print(json.dumps({"event": "error", "error_type": type(exc).__name__, "error": str(exc)}), flush=True)
    raise
'''


def validate_outputs(output: Path) -> dict:
    analysis_path = output / "analysis.json"
    zip_path = output / "ecd_epitope_design.zip"
    card_path = output / "card_evidence_completeness.json"
    checks = {"analysis": False, "zip": False, "cards": False}
    errors = []
    analysis = {}
    if analysis_path.is_file():
        try:
            analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
            checks["analysis"] = isinstance(analysis, dict) and bool(analysis)
            if not checks["analysis"]:
                errors.append("analysis.json is empty or not an object")
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"analysis.json invalid: {exc}")
    else:
        errors.append("analysis.json is missing")
    if zip_path.is_file():
        try:
            with zipfile.ZipFile(zip_path) as archive:
                damaged = archive.testzip()
                names = set(archive.namelist())
            checks["zip"] = damaged is None and "analysis.json" in names
            if not checks["zip"]:
                errors.append(f"Design ZIP failed integrity/content check: {damaged or 'analysis.json missing'}")
        except (OSError, zipfile.BadZipFile) as exc:
            errors.append(f"Design ZIP invalid: {exc}")
    else:
        errors.append("ecd_epitope_design.zip is missing")
    if card_path.is_file():
        try:
            cards = json.loads(card_path.read_text(encoding="utf-8"))
            checks["cards"] = cards.get("status") == "PASS" and int(cards.get("card_count", 0)) > 0
            if not checks["cards"]:
                errors.append("Card evidence receipt is incomplete or contains no cards")
        except (OSError, json.JSONDecodeError, AttributeError, TypeError, ValueError) as exc:
            errors.append(f"Card evidence receipt invalid: {exc}")
    else:
        errors.append("card_evidence_completeness.json is missing")
    normalized = output / "structure_for_energy.pdb"
    checks["normalized_structure"] = normalized.is_file()
    if not checks["normalized_structure"]:
        errors.append("structure_for_energy.pdb is missing")
    source_identity_path = output / "run_source_identity_v1.json"
    worker_identity = {"status": "NOT_DETERMINED", "reason": "run_source_identity_v1.json is missing",
                       "path": str(source_identity_path)}
    if source_identity_path.is_file():
        try:
            identity = json.loads(source_identity_path.read_text(encoding="utf-8"))
            workers = identity.get("worker_model_identities", {})
            missing = []
            for name in ("spurs_single", "spurs_multi", "discotope"):
                provenance = workers.get(name) if isinstance(workers, dict) else None
                if not isinstance(provenance, dict) or not provenance.get("actual_device") or not provenance.get("weights_sha256"):
                    missing.append(name)
                elif provenance.get("actual_device") == "cuda" and not provenance.get("gpu_name"):
                    missing.append(name + ":gpu_name")
            worker_identity = {"status": "PASS" if not missing else "NOT_DETERMINED",
                "path": str(source_identity_path), "missing_identities": missing,
                "worker_model_identities": workers,
                "normalized_structure_sha256": identity.get("normalized_structure_sha256")}
        except (OSError, json.JSONDecodeError, AttributeError) as exc:
            worker_identity = {"status": "NOT_DETERMINED", "path": str(source_identity_path), "reason": str(exc)}
    joint = analysis.get("joint_v4_ml")
    joint = joint if isinstance(joint, dict) else {}
    return {"status": "PASS" if all(checks.values()) and not errors else "NEEDS_REVIEW",
            "checks": checks, "errors": errors, "analysis_path": str(analysis_path),
            "zip_path": str(zip_path), "card_evidence_path": str(card_path),
            "normalized_structure": {"path": str(normalized), "sha256": sha256_file(normalized)}
                if normalized.is_file() else None,
            "worker_model_identity": worker_identity,
            "protocol": analysis.get("execution_protocol", joint.get("execution_protocol"))}


def legacy_config_identity(snapshot_package: Path) -> dict:
    config = snapshot_package / "config.py"
    text = config.read_text(encoding="utf-8")
    import re
    match = re.search(r'"procaliper"\s*:\s*\{.*?"python"\s*:\s*r?["\']([^"\']+)', text, re.S)
    procaliper_python = Path(match.group(1)) if match else None
    return {"path": str(config.resolve()), "sha256": sha256_file(config),
            "procaliper_python": str(procaliper_python) if procaliper_python else None,
            "procaliper_python_exists": bool(procaliper_python and procaliper_python.is_file()),
            "procaliper_python_sha256": sha256_file(procaliper_python)
                if procaliper_python and procaliper_python.is_file() else None}


def write_report(acceptance_dir: Path, report: dict) -> None:
    atomic_json(acceptance_dir / "baseline_report.json", report)
    lines = ["# AIB234 legacy cold baseline", "", f"Status: `{report.get('status')}`.",
             f"Commit: `{report.get('baseline_commit')}`.",
             f"Pipeline wall seconds: `{report.get('elapsed_pipeline_seconds')}`; total including preparation: `{report.get('elapsed_total_seconds')}`.",
             f"Input SHA256: `{report.get('input_sha256')}`.",
             f"Normalized structure: `{(report.get('validation') or {}).get('normalized_structure')}`.",
             f"Worker identity comparison: `{(report.get('validation') or {}).get('worker_model_identity', {}).get('status', 'NOT_DETERMINED')}`.",
             f"Cold cache evidence: `{report.get('cold_cache_evidence')}`.", "",
             "A runtime reduction is not determined by this run alone. The controller must compare both runs' GPU, model weights, FoldX binary, cache condition, and scientific acceptance evidence.", "",
             "## Validation", "", "```json", json.dumps(report.get("validation"), indent=2, ensure_ascii=False), "```", "",
             "## Errors", "", "```json", json.dumps(report.get("errors", []), indent=2, ensure_ascii=False), "```", ""]
    (acceptance_dir / "baseline_report.md").write_text("\n".join(lines), encoding="utf-8")


def _fingerprints(paths: dict[str, str]) -> dict:
    return {name: {"path": path, "sha256": sha256_file(Path(path)), "size_bytes": Path(path).stat().st_size}
            for name, path in paths.items()}


def execute(*, root: Path = ROOT, acceptance_dir: Path = ACCEPTANCE_DIR,
            source_request: Path = SOURCE_REQUEST, structure: Path = SOURCE_STRUCTURE,
            foldx_binary: Path | None = None, archive_provider=None, popen_factory=None,
            model_paths_override: dict[str, Path] | None = None, commit_resolver=resolve_commit,
            max_foldx_root_length: int = 100) -> dict:
    root, acceptance_dir = Path(root).resolve(), Path(acceptance_dir).resolve()
    report = {"status": "NEEDS_REVIEW", "baseline_commit": None, "input_sha256": None,
              "started_at_utc": datetime.now(timezone.utc).isoformat(), "elapsed_pipeline_seconds": None,
              "elapsed_total_seconds": None,
              "errors": [], "validation": None, "runtime_comparison": "NOT_DETERMINED"}
    started = time.monotonic()
    try:
        plan = preflight(root=root, acceptance_dir=acceptance_dir, source_request=source_request,
                         structure=structure, foldx_binary=foldx_binary,
                         model_paths_override=model_paths_override, commit_resolver=commit_resolver,
                         max_foldx_root_length=max_foldx_root_length)
        commit = plan["baseline_commit"]
        baseline = acceptance_dir / "baseline"
        baseline.mkdir(parents=True, exist_ok=False)
        report.update(baseline_commit=commit, input_path=plan["input_path"], input_sha256=plan["input_sha256"])
        binary = Path(plan["model_identity_paths"]["foldx_binary"])
        identities = plan["model_identity_paths"]
        report["model_identity"] = _fingerprints(identities)
        if archive_provider is None:
            archive = subprocess.run(["git", "archive", "--format=tar", "--prefix=ecd_epitope_scanner/", commit],
                                     cwd=root, check=True, capture_output=True).stdout
        else:
            archive = archive_provider(commit)
        package = safe_extract_archive(archive, baseline / "baseline_snapshot")
        report["source_snapshot"] = {"commit": commit, "archive_sha256": hashlib.sha256(archive).hexdigest(),
                                     "package_path": str(package)}
        report["legacy_config_identity"] = legacy_config_identity(package)
        uploads = baseline / "uploads"
        uploads.mkdir()
        local_structure = uploads / "structure.pdb"
        local_structure.write_bytes(Path(structure).read_bytes())
        if sha256_file(local_structure) != EXPECTED_INPUT_SHA256:
            raise ValueError("Copied baseline input hash mismatch")
        output = baseline / "legacy_output"
        output.mkdir()
        score_cache = baseline / "score_cache"
        score_cache.mkdir()
        engine_root = baseline / "fx"
        engine_root.mkdir()
        direct_repair_cache = baseline / "foldx_repair_cache"
        short_repair_cache = engine_root / "foldx_repair_cache"
        if direct_repair_cache.exists() or short_repair_cache.exists() or any(score_cache.iterdir()):
            raise FileExistsError("Cold score/RepairPDB cache condition could not be established")
        request_data = build_legacy_request(json.loads(Path(source_request).read_text(encoding="utf-8")),
            structure=local_structure, output=output, baseline_root=baseline,
            model_assets_root=root, foldx_binary=binary)
        request_path = baseline / "request.json"
        atomic_json(request_path, request_data)
        output_length = len(str(output.resolve()))
        direct_execution_root = os.name != "nt" or output_length + 90 <= 220
        report["input_copy_path"] = str(local_structure)
        report["score_cache"] = {"path": str(score_cache), "entries_at_start": 0}
        report["cold_cache_evidence"] = {
            "repair_cache_possible_paths": [str(direct_repair_cache), str(short_repair_cache)],
            "all_possible_repair_cache_paths_absent_before_run": (
                not direct_repair_cache.exists() and not short_repair_cache.exists()),
            "execution_root_selection": {"platform": os.name, "output_path_length": output_length,
                "direct_output_root_selected": direct_execution_root,
                "rule": "old worker selects output when non-Windows or len(output)+90 <= 220; otherwise engine_work_root"},
            "score_cache_path": str(score_cache), "score_cache_entries_at_start": 0}
        event_path = baseline / "events.jsonl"
        child_args = [sys.executable, "-c", CHILD_CODE, str(package.parent), str(request_path)]
        process_factory = popen_factory or subprocess.Popen
        pipeline_started = pipeline_finished = None
        child = process_factory(child_args, cwd=str(root.parent), stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, bufsize=1)
        with event_path.open("w", encoding="utf-8") as events:
            for raw in child.stdout:
                raw = raw.rstrip("\r\n")
                try:
                    event = json.loads(raw)
                except json.JSONDecodeError:
                    event = {"event": "stdout", "text": raw}
                if event.get("event") == "pipeline_start":
                    pipeline_started = time.monotonic()
                elif event.get("event") in {"result", "error"}:
                    pipeline_finished = time.monotonic()
                events.write(json.dumps(event, ensure_ascii=False, default=str) + "\n")
                if event.get("event") == "error":
                    report["errors"].append(f"{event.get('error_type')}: {event.get('error')}")
                if event.get("event") == "result":
                    report["pipeline_result_keys"] = event.get("result_keys", [])
                events.flush()
        return_code = child.wait()
        if pipeline_started is not None:
            report["elapsed_pipeline_seconds"] = (pipeline_finished or time.monotonic()) - pipeline_started
        report["child_return_code"] = return_code
        if return_code != 0:
            report["errors"].append(f"Archived pipeline exited with code {return_code}")
        report["validation"] = validate_outputs(output)
        if report["validation"]["status"] != "PASS":
            report["errors"].extend(report["validation"]["errors"])
        report["normalized_structure"] = report["validation"].get("normalized_structure")
        report["elapsed_total_seconds"] = time.monotonic() - started
        report["elapsed_seconds"] = report["elapsed_total_seconds"]
        report["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        report["status"] = "COMPLETE" if return_code == 0 and report["validation"]["status"] == "PASS" else "NEEDS_REVIEW"
    except Exception as exc:
        report["errors"].append(f"{type(exc).__name__}: {exc}")
        report["elapsed_total_seconds"] = time.monotonic() - started
        report["elapsed_seconds"] = report["elapsed_total_seconds"]
        report["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    if not (acceptance_dir / "baseline_report.json").exists():
        write_report(acceptance_dir, report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true", help="Read-only paths, hashes, and duplicate-run check")
    mode.add_argument("--execute", action="store_true", help="Run the archived AIB234 legacy pipeline")
    args = parser.parse_args()
    result = preflight() if args.preflight else execute()
    print(json.dumps({"status": result.get("status"), "report": str(REPORT_JSON)}, indent=2, ensure_ascii=False))
    if args.execute and result.get("status") != "COMPLETE":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
