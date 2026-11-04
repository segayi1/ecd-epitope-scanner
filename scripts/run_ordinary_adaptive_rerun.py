"""Prepare or run the frozen 13-case ordinary adaptive batch.

Run from the repository parent with its venv Python. Execution is deliberately
opt-in and requires a PASS receipt tied to this exact manifest hash.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sqlite3
import sys
import time
import uuid
import zipfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "workspace_data/ordinary_adaptive_rerun_20261026_manifest.json"
STATE = ROOT / "workspace_data/ordinary_adaptive_rerun_20261026_state.json"
OUT = ROOT / "workspace_data/ordinary_adaptive_rerun_20261026_report"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def production_source_sha256() -> str:
    root = Path(__file__).resolve().parents[1]
    if str(root.parent) not in sys.path:
        sys.path.insert(0, str(root.parent))
    from ecd_epitope_scanner.joint_design_v4.run_audit import source_snapshot
    return source_snapshot(root)["source_files_sha256"]


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)


def preflight(manifest: dict, root: Path = ROOT) -> list[dict]:
    cases = manifest.get("rerun_cases", [])
    if manifest.get("case_count") != 13 or len(cases) != 13:
        raise ValueError("Manifest must contain exactly 13 rerun cases")
    if len({case["case_id"] for case in cases}) != len(cases):
        raise ValueError("Rerun case IDs must be unique")
    staged = []
    for case in cases:
        structure = Path(case["structure_path"])
        if not structure.is_file() or sha256(structure) != case["structure_sha256"]:
            raise ValueError(f"Frozen structure missing or hash changed: {structure}")
        confidence = Path(case["confidence_path"]) if case.get("confidence_path") else None
        if case["input_type"] == "af3_upload":
            if not confidence or not confidence.is_file() or sha256(confidence) != case["confidence_sha256"]:
                raise ValueError(f"Frozen AF3 confidence missing or hash changed: {confidence}")
        if case.get("metadata", {}).get("disable_external_engines"):
            raise ValueError(f"External engines may not be disabled: {case['case_id']}")
        staged.append(case)
    return staged


def request_for(case: dict, structure: Path, confidence: Path | None):
    from ecd_epitope_scanner.pipeline import PipelineRequest

    old = case.get("request_contract", {})
    source_request = json.loads(Path(case["source_request"]).read_text(encoding="utf-8")) if case.get("source_request") else {}
    binder = case.get("cohort", "").startswith("prior_binder")
    values = {
        "input_type": case["input_type"], "identifier": case.get("identifier") or case.get("name") or case["case_id"],
        "accession": case.get("accession"), "chain": case["chain"],
        "structure_path": str(structure), "confidence_path": str(confidence) if confidence else None,
        "strict_af3": True, "allow_partial": False,
        "minimum_regions": old.get("minimum_regions", 5), "maximum_regions": old.get("maximum_regions", 12),
        "assay_format": "unspecified" if binder else old.get("assay_format", "cell_surface"),
        "design_engine": old.get("design_engine", "joint_v4_ml"),
        "ml_mode": "strict", "surface_strategy": "legacy", "max_sectors": old.get("max_sectors", 12),
        # Let current chain-length defaults choose the display target; old saved
        # requests may carry the superseded 8-residue short-chain target.
        "target_sector_residues": None,
        "optimization_profile": "legacy", "export_mode": "deliverable",
        "analysis_target": "binder" if binder else "antigen",
        "biological_context_mode": "provided_structure",
        "sequence": source_request.get("sequence"),
        "ecd_start": source_request.get("ecd_start"), "ecd_end": source_request.get("ecd_end"),
        # Fresh metadata prevents old prediction_results, execution_queue,
        # disabled-engine flags, and external cache paths from crossing runs.
        "metadata": {"ordinary_adaptive_rerun_case": case["case_id"]},
    }
    return PipelineRequest(**values)


def assert_no_active_jobs(db_path: Path) -> None:
    if not db_path.exists():
        return
    with sqlite3.connect(db_path) as connection:
        rows = connection.execute("SELECT id,status FROM jobs WHERE status IN ('RUNNING','QUEUED')").fetchall()
    if rows:
        raise RuntimeError(f"Active jobs prevent batch start: {rows}")


def execute(manifest: dict, receipt_path: Path, state_path: Path = STATE, resume: bool = False) -> dict:
    expected = hashlib.sha256(MANIFEST.read_bytes()).hexdigest()
    runner_sha256 = sha256(Path(__file__).resolve())
    source_sha256 = production_source_sha256()
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if (receipt.get("status") != "PASS" or receipt.get("manifest_sha256") != expected
            or receipt.get("runner_sha256") != runner_sha256
            or receipt.get("production_source_sha256") != source_sha256):
        raise ValueError("Acceptance receipt must be PASS and match the exact manifest, runner, and production source SHA256")
    cases = preflight(manifest)
    sys.path.insert(0, str(ROOT.parent))
    from ecd_epitope_scanner.jobs import JobRunner, JobStore

    db_path = ROOT / "workspace_data" / "jobs.sqlite3"
    assert_no_active_jobs(db_path)
    store = JobStore(ROOT / "workspace_data")
    if state_path.exists() and not resume:
        raise FileExistsError(f"Run state already exists; pass --resume to continue: {state_path}")
    if resume and not state_path.exists():
        raise FileNotFoundError(f"Cannot resume without saved run state: {state_path}")
    state = json.loads(state_path.read_text(encoding="utf-8")) if resume else {
        "manifest_sha256": expected, "runner_sha256": runner_sha256,
        "production_source_sha256": source_sha256, "jobs": {}, "status": "RUNNING"}
    if (state.get("manifest_sha256") != expected or state.get("runner_sha256") != runner_sha256
            or state.get("production_source_sha256") != source_sha256):
        raise ValueError("State belongs to a different frozen manifest, runner, or production source")
    runner = JobRunner(store)
    try:
        state["status"] = "RUNNING"
        state.pop("error", None)
        atomic_json(state_path, state)
        for case in cases:
            if production_source_sha256() != source_sha256:
                raise RuntimeError("Production source changed during batch; stopping before the next case")
            saved = state["jobs"].get(case["case_id"])
            if saved:
                job = store.get(saved["job_id"])
                if job["status"] in {"COMPLETE", "FAILED", "WAITING_FOR_AF3_UPLOAD"}:
                    saved.update(job_status=job["status"], error=job.get("error"))
                    validation = capture_post_validation(case, job, saved)
                    saved["post_validation"] = validation
                    saved["status"] = validation["status"] if job["status"] == "COMPLETE" and validation else job["status"]
                    atomic_json(state_path, state)
                    publish_reports(state, manifest, state_path)
                    if job["status"] != "COMPLETE" or not validation or validation["status"] not in {"PASS", "DIAGNOSTIC_COMPLETE"}:
                        state.update(status="HALTED_FOR_REVIEW", blocked_case_id=case["case_id"])
                        atomic_json(state_path, state)
                        publish_reports(state, manifest, state_path)
                        return state
                    continue
                raise RuntimeError(f"Existing job is nonterminal; inspect before resume: {saved['job_id']} {job['status']}")
            # Recheck immediately before snapshotting each input.
            original = Path(case["structure_path"])
            if sha256(original) != case["structure_sha256"]:
                raise ValueError(f"Frozen source changed before case: {original}")
            job_id = uuid.uuid4().hex[:12]
            upload = ROOT / "workspace_data" / job_id / "uploads"
            upload.mkdir(parents=True, exist_ok=True)
            local_structure = upload / original.name
            local_structure.write_bytes(original.read_bytes())
            if sha256(local_structure) != case["structure_sha256"]:
                raise ValueError(f"Snapshotted structure hash mismatch: {local_structure}")
            local_confidence = None
            if case.get("confidence_path"):
                source_confidence = Path(case["confidence_path"])
                if sha256(source_confidence) != case["confidence_sha256"]:
                    raise ValueError(f"Frozen confidence changed before case: {source_confidence}")
                local_confidence = upload / source_confidence.name
                local_confidence.write_bytes(source_confidence.read_bytes())
                if sha256(local_confidence) != case["confidence_sha256"]:
                    raise ValueError(f"Snapshotted confidence hash mismatch: {local_confidence}")
            request = request_for(case, local_structure, local_confidence)
            job = store.create(request, job_id=job_id)
            saved = {"job_id": job["id"], "status": "QUEUED", "source_hashes": {
                "structure": case["structure_sha256"], "confidence": case.get("confidence_sha256")}}
            state["jobs"][case["case_id"]] = saved
            atomic_json(state_path, state)
            started = time.monotonic()
            runner._run(job["id"])
            job = store.get(job["id"])
            saved.update(job_status=job["status"], elapsed_seconds=time.monotonic() - started,
                         error=job.get("error"), workdir=job["workdir"])
            validation = capture_post_validation(case, job, saved)
            saved["post_validation"] = validation
            saved["status"] = validation["status"] if job["status"] == "COMPLETE" and validation else job["status"]
            atomic_json(state_path, state)
            publish_reports(state, manifest, state_path)
            if job["status"] != "COMPLETE" or not validation or validation["status"] not in {"PASS", "DIAGNOSTIC_COMPLETE"}:
                state.update(status="HALTED_FOR_REVIEW", blocked_case_id=case["case_id"])
                atomic_json(state_path, state)
                publish_reports(state, manifest, state_path)
                return state
        terminal = all(v.get("job_status", v.get("status")) in {"COMPLETE", "FAILED", "WAITING_FOR_AF3_UPLOAD"}
                       for v in state["jobs"].values())
        validated = all(v.get("post_validation", {}).get("status") in {"PASS", "DIAGNOSTIC_COMPLETE"}
                        for v in state["jobs"].values() if v.get("job_status", v.get("status")) == "COMPLETE")
        state["status"] = ("COMPLETE" if terminal and validated and
            all(v.get("job_status", v.get("status")) == "COMPLETE" for v in state["jobs"].values())
            else "FINISHED_WITH_REVIEW" if terminal and not validated else "FINISHED_WITH_FAILURES")
        atomic_json(state_path, state)
    except BaseException as exc:
        state.update(status="INTERRUPTED", error=f"{type(exc).__name__}: {exc}")
        atomic_json(state_path, state)
        publish_reports(state, manifest, state_path)
        raise
    finally:
        runner.executor.shutdown(wait=False, cancel_futures=True)
    publish_reports(state, manifest, state_path)
    return state


def write_case_evidence(case: dict, job: dict, saved: dict) -> None:
    workdir = Path(job["workdir"])
    analysis = workdir / "analysis.json"
    if not analysis.exists():
        if job.get("status") == "COMPLETE":
            return {"status": "NEEDS_REVIEW", "errors": ["COMPLETE job is missing analysis.json"]}
        return None
    data = json.loads(analysis.read_text(encoding="utf-8"))
    designs = data.get("designs", [])
    regions = data.get("regions", [])
    saved["result_evidence"] = {"analysis": str(analysis), "status": job["status"],
        "design_count": len(designs) if isinstance(designs, list) else None,
        "region_count": len(regions) if isinstance(regions, list) else None,
        "recommendation_status": (data.get("design_summary") or {}).get("recommendation_status"),
        "protocol": case.get("request_contract", {})}
    for name in ("panel_candidates.json", "sector_recommendations_v4.json", "model_evidence_diagnostics_v1.json",
                 "ordinary_candidate_ledger.json", "card_evidence_completeness.json", "run_manifest.json",
                 "sector_recommendations_v4.csv", "ecd_epitope_design.zip"):
        path = workdir / name
        if path.is_file():
            saved["result_evidence"][name.removesuffix(".json")] = str(path)
    return validate_case_output(workdir, data)


def capture_post_validation(case: dict, job: dict, saved: dict) -> dict | None:
    try:
        return write_case_evidence(case, job, saved)
    except Exception as exc:
        return {"status": "NEEDS_REVIEW", "checks": {},
                "technical_failures": {"post_validation": True},
                "errors": [f"post-validation raised {type(exc).__name__}: {exc}"]}


def _panel_hypothesis_signature_errors(workdir: Path) -> list[str]:
    paths = {name: workdir / name for name in (
        "selected_panel_v4.json", "panel_identifiability_v41.json", "panel_candidates.json")}
    missing = [name for name, path in paths.items() if not path.is_file()]
    if missing:
        return [f"site-first panel hypothesis signature artifact missing: {name}" for name in missing]
    try:
        panel, identifiability, candidate_rows = (
            json.loads(paths[name].read_text(encoding="utf-8"))
            for name in ("selected_panel_v4.json", "panel_identifiability_v41.json", "panel_candidates.json"))
    except Exception as exc:
        return [f"site-first panel hypothesis signature artifact invalid: {type(exc).__name__}: {exc}"]

    selected_ids = panel.get("selected_candidate_ids") if isinstance(panel, dict) else None
    ident_ids = identifiability.get("selected_candidate_ids") if isinstance(identifiability, dict) else None
    if not isinstance(selected_ids, list) or not selected_ids or not all(isinstance(cid, str) for cid in selected_ids):
        return ["site-first panel selected candidate IDs are missing or invalid"]
    errors = []
    if ident_ids != selected_ids:
        errors.append("site-first panel selected IDs differ from panel identifiability IDs")
    signatures = identifiability.get("hypothesis_signatures") if isinstance(identifiability, dict) else None
    panel_signatures = panel.get("hypothesis_signatures") if isinstance(panel, dict) else None
    if not isinstance(signatures, dict) or not signatures:
        return errors + ["site-first panel hypothesis signatures are missing or unknown"]
    if panel_signatures != signatures:
        errors.append("site-first panel signatures differ between selected_panel_v4 and panel_identifiability_v41")
    if not isinstance(candidate_rows, list):
        return errors + ["site-first panel candidate evidence is missing or invalid"]
    by_id = {}
    duplicates = set()
    for row in candidate_rows:
        cid = row.get("candidate_id") if isinstance(row, dict) else None
        if not isinstance(cid, str):
            continue
        if cid in by_id:
            duplicates.add(cid)
        by_id[cid] = row
    if duplicates:
        errors.append("site-first panel candidate IDs are duplicated: " + ", ".join(sorted(duplicates)[:3]))

    unknown = []
    mismatches = []
    for hypothesis_id, values in signatures.items():
        if not isinstance(values, list) or len(values) != len(selected_ids):
            errors.append(f"site-first signature vector has invalid length for {hypothesis_id}")
            continue
        for index, candidate_id in enumerate(selected_ids):
            expected = values[index]
            if (not isinstance(expected, (int, float)) or isinstance(expected, bool)
                    or not math.isfinite(expected)):
                unknown.append(f"{hypothesis_id}/{candidate_id}: panel signature value")
                continue
            candidate = by_id.get(candidate_id)
            if candidate is None:
                unknown.append(f"{hypothesis_id}/{candidate_id}: selected candidate")
                continue
            maps = [(key, candidate.get(key)) for key in ("hypothesis_perturbation", "hypothesis_response")
                    if isinstance(candidate.get(key), dict) and candidate.get(key)]
            if not maps:
                unknown.append(f"{hypothesis_id}/{candidate_id}: candidate response")
                continue
            for key, response in maps:
                observed = response.get(hypothesis_id)
                if (not isinstance(observed, (int, float)) or isinstance(observed, bool)
                        or not math.isfinite(observed)):
                    unknown.append(f"{hypothesis_id}/{candidate_id}: {key}")
                elif not math.isclose(float(expected), float(observed), rel_tol=1e-9, abs_tol=1e-8):
                    mismatches.append(f"{hypothesis_id}/{candidate_id}/{key}: {expected} != {observed}")
    if unknown:
        errors.append(f"site-first panel signature evidence is unknown ({len(unknown)} values; first: {unknown[0]})")
    if mismatches:
        errors.append(f"site-first panel signatures disagree with selected candidate responses ({len(mismatches)} values; first: {mismatches[0]})")
    return errors


def _site_first_display_errors(workdir: Path) -> list[str]:
    names = ("selected_panel_v4.json", "panel_candidates.json", "sector_recommendations_v4.json")
    try:
        panel, candidates, recommendations = (json.loads((workdir / name).read_text(encoding="utf-8"))
                                               for name in names)
        if not isinstance(panel, dict) or not isinstance(candidates, list) or not isinstance(recommendations, list):
            raise ValueError("unexpected artifact shape")
    except (OSError, ValueError, TypeError) as exc:
        return [f"site-first display artifacts missing or invalid: {exc}"]
    errors = []
    by_id = {row.get("candidate_id"): row for row in candidates if isinstance(row, dict)}

    def confirmed_repeat_pass(audit):
        if str(ROOT.parent) not in sys.path:
            sys.path.insert(0, str(ROOT.parent))
        from ecd_epitope_scanner.joint_design_v4.site_first_display import (
            confirmed_repeat_pass as shared_confirmed_repeat_pass)

        return shared_confirmed_repeat_pass(audit)

    def repeat_audit_errors(audit, candidate_id):
        values = audit.get("foldx_repeat_slot_values_kcal_mol")
        confirmed = confirmed_repeat_pass(audit)
        result = []
        if (audit.get("sampling_decision_diagnostic_status") != "COMPLETE"
                or (not confirmed and (audit.get("sampling_decision_sensitive") is not False
                    or audit.get("recommendation_risk_label") != "NO_WITHIN_RUN_GATE_CROSSING"))
                or audit.get("cross_execution_conflict") is True):
            result.append(f"site-first display candidate needs repeat evidence review: {candidate_id}")
        if (not isinstance(values, list) or len(values) not in (3, 10)
                or any(not isinstance(value, (int, float)) or isinstance(value, bool)
                       or not math.isfinite(value) for value in values)):
            result.append(f"site-first display candidate has invalid repeat evidence: {candidate_id}")
        elif any(value > 2 for value in values) and not confirmed:
            result.append(f"site-first display candidate crosses the FoldX gate: {candidate_id}")
        return result

    def confirmation_evidence_matches(audit, diagnostic):
        keys = ("repeat_confirmation_policy", "repeat_confirmation_status",
                "sampling_decision_diagnostic_status", "sampling_decision_sensitive",
                "recommendation_risk_label", "foldx_repeat_slots",
                "foldx_repeat_slot_values_kcal_mol", "foldx_repeat_gate_kcal_mol",
                "cross_execution_conflict")
        return all(audit.get(key) == diagnostic.get(key) for key in keys)

    def candidate_errors(candidate_id):
        row = by_id.get(candidate_id)
        if not row:
            return [f"site-first display candidate missing from raw evidence: {candidate_id}"]
        result = []
        impact = row.get("impact_evidence") or {}
        audit = impact.get("stability_repeat_audit") or {}
        values = audit.get("foldx_repeat_slot_values_kcal_mol")
        if row.get("formal_eligible") is not True or row.get("display_eligible") is False:
            result.append(f"site-first display candidate is not eligible: {candidate_id}")
        elif row.get("display_eligible") is not True:
            result.append(f"site-first display qualification unknown: {candidate_id}")
        if row.get("qualification") == "REVIEW_REQUIRED" or any(
                mutation.get("review_required") is True for mutation in row.get("mutations", [])):
            result.append(f"site-first display candidate needs structural review: {candidate_id}")
        if impact.get("policy") == "regional_mechanism_pareto_v1" and impact.get("missing_axes"):
            result.append(f"site-first display candidate has incomplete mechanism evidence: {candidate_id}")
        result.extend(repeat_audit_errors(audit, candidate_id))
        return result

    for candidate_id in panel.get("selected_candidate_ids", []):
        errors.extend(candidate_errors(candidate_id))
    if panel.get("sector_recommendations") != recommendations:
        errors.append("site-first panel and sector recommendation artifacts disagree")
    for primary in recommendations:
        if not isinstance(primary, dict):
            errors.append("site-first sector recommendation is invalid")
            continue
        for card in [primary, *(primary.get("alternatives") or []),
                     *(primary.get("all_formal_local_candidates") or [])]:
            if not isinstance(card, dict) or not card.get("candidate_id"):
                continue
            candidate_id = card["candidate_id"]
            errors.extend(candidate_errors(candidate_id))
            diagnostic = card.get("stability_repeat_diagnostic") or {}
            raw_audit = (by_id.get(candidate_id, {}).get("impact_evidence") or {}).get("stability_repeat_audit") or {}
            confirmed = confirmed_repeat_pass(raw_audit)
            if confirmed and (not confirmed_repeat_pass(diagnostic)
                              or not confirmation_evidence_matches(raw_audit, diagnostic)):
                errors.append(f"site-first card repeat confirmation disagrees with raw evidence: {candidate_id}")
            if (card.get("formal_eligible") is not True or card.get("status") != "PASSED"
                    or repeat_audit_errors(diagnostic, candidate_id)):
                errors.append(f"site-first card needs review or lacks formal qualification: {candidate_id}")
    return list(dict.fromkeys(errors))


def regional_principle_coverage(ledger: dict) -> dict:
    queues = ledger.get("queues")
    if (not isinstance(queues, list) or not queues or any(
            not isinstance(row, dict) or not row.get("sector_id") or not row.get("principle")
            for row in queues)):
        return {"status": "UNKNOWN", "gaps": []}
    terminal = {"candidate_failures", "budget_limited", "technical_failure", "missing_evidence", "exhausted"}
    gaps = [{"sector_id": row.get("sector_id"), "principle": row.get("principle"), "state": row.get("state")}
            for row in queues if row.get("state") in terminal]
    other = [row for row in queues if row.get("state") not in terminal | {"complete", "not_applicable"}]
    return {"status": "INCOMPLETE" if gaps else "UNKNOWN" if other else "COMPLETE", "gaps": gaps}


def validate_case_output(workdir: Path, analysis: dict | None = None,
                         expected_protocol: str = "ordinary_three_principle_adaptive_v1",
                         spurs_request_limit: int = 1000) -> dict:
    workdir = Path(workdir)
    errors, checks = [], {}
    analysis = analysis or json.loads((workdir / "analysis.json").read_text(encoding="utf-8"))
    protocol = analysis.get("execution_protocol") or {}
    protocol_name = protocol.get("name", protocol.get("policy")) if isinstance(protocol, dict) else protocol
    checks["ordinary_protocol"] = protocol_name == expected_protocol
    if not checks["ordinary_protocol"]:
        errors.append(f"unexpected execution protocol: {protocol_name!r}; expected {expected_protocol!r}")
    if protocol_name == "ordinary_site_first_v2":
        panel_signature_errors = _panel_hypothesis_signature_errors(workdir)
        checks["panel_hypothesis_signatures"] = not panel_signature_errors
        errors.extend(panel_signature_errors)
        display_errors = _site_first_display_errors(workdir)
        checks["formal_display_eligibility"] = not display_errors
        errors.extend(display_errors)
    ledger_path = workdir / "ordinary_candidate_ledger.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8")) if ledger_path.is_file() else {}
    spurs, foldx = ledger.get("spurs") or {}, ledger.get("foldx") or {}
    spurs_requests = spurs.get("unique_requests")
    foldx_slots = foldx.get("attempted_slots")
    checks["spurs_budget"] = isinstance(spurs_requests, int) and 0 <= spurs_requests <= spurs_request_limit
    checks["foldx_budget"] = isinstance(foldx_slots, int) and 0 <= foldx_slots <= 300
    if not checks["spurs_budget"]: errors.append(f"invalid or exceeded SPURS request budget: {spurs_requests!r}")
    if not checks["foldx_budget"]: errors.append(f"invalid or exceeded FoldX slot budget: {foldx_slots!r}")
    technical = {"spurs": bool(spurs.get("technical_failure")), "foldx": bool(foldx.get("technical_failure"))}
    checks["no_technical_failure"] = not any(technical.values())
    if not checks["no_technical_failure"]:
        errors.append(f"model technical failures: {technical}")
    passed_foldx = [row for row in ledger.get("foldx_results", []) if row.get("status") == "PASS"]
    for row in passed_foldx:
        slots, values = row.get("slots"), row.get("scores")
        if (not isinstance(slots, list) or len(slots) not in (3, 10) or len(set(slots)) != len(slots)
                or not isinstance(values, list) or len(values) != len(slots)
                or not all(isinstance(value, (int, float)) and math.isfinite(value) for value in values)
                or not isinstance(row.get("median"), (int, float)) or not math.isfinite(row["median"])
                or not math.isclose(row["median"], __import__("statistics").median(values), rel_tol=1e-9, abs_tol=1e-9)
                or row["median"] > 2.0 or slots != list(range(len(slots)))
                or (len(slots) == 3 and any(value > 2.0 for value in values))
                or (len(slots) == 10 and sum(value > 2.0 for value in values[:3]) != 1)):
            errors.append(f"invalid PASS FoldX evidence for {row.get('candidate_id')}: slots/values/median/gate")
            break
        pairs = row.get("paired_models") or [pair for stage in (row.get("audit") or {}).get("stages", [])
                                             for record in stage.get("records", [])
                                             for pair in record.get("paired_models", [])]
        if len(pairs) != len(slots):
            errors.append(f"FoldX PASS pair count mismatch for {row.get('candidate_id')}: {len(pairs)} != {len(slots)}")
            break
        try:
            from ecd_epitope_scanner.surface_perturbation.paired import validate_pair
            for pair in pairs:
                validate_pair(pair)
        except Exception as exc:
            errors.append(f"paired model source/hash validation failed for {row.get('candidate_id')}: {exc}")
            break
    checks["passed_foldx_evidence"] = not any("FoldX PASS" in error or "PASS FoldX" in error or "paired model" in error for error in errors)
    checks["pass_foldx_count"] = len(passed_foldx)

    recommendation_path, panel_path, csv_path = (workdir / name for name in (
        "sector_recommendations_v4.json", "panel_candidates.json", "sector_recommendations_v4.csv"))
    recommendations = json.loads(recommendation_path.read_text(encoding="utf-8")) if recommendation_path.is_file() else []
    cards = [row for recommendation in recommendations if isinstance(recommendation, dict)
             for row in [recommendation, *(recommendation.get("alternatives") or [])]
             if row.get("candidate_id")]
    receipt_path = workdir / "card_evidence_completeness.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8")) if receipt_path.is_file() else None
    if cards and (not receipt or receipt.get("status") != "PASS"):
        errors.append("recommendation cards lack PASS card_evidence_completeness.json")
    if receipt:
        for name, digest in (receipt.get("inputs_sha256") or {}).items():
            source = workdir / name
            if not source.is_file() or sha256(source) != digest:
                errors.append(f"card completeness source hash mismatch: {name}")
                break

    manifest_path, zip_path = workdir / "run_manifest.json", workdir / "ecd_epitope_design.zip"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    hashes = manifest.get("files") or {}
    key_files = ("analysis.json", "panel_candidates.json", "sector_recommendations_v4.json",
                 "sector_recommendations_v4.csv", "ordinary_candidate_ledger.json")
    analysis_files = analysis.get("files") or {}
    for name in key_files:
        path = workdir / name
        if (not path.is_file() or name not in hashes or sha256(path) != hashes.get(name)
                or analysis_files.get(name) != name):
            errors.append(f"run_manifest source hash mismatch or missing artifact: {name}")
            break
    if csv_path.is_file() and recommendation_path.is_file():
        try:
            previous_limit = csv.field_size_limit()
            csv.field_size_limit(min(max(previous_limit, csv_path.stat().st_size), 16 * 1024 * 1024))
            try:
                with csv_path.open(encoding="utf-8-sig", newline="") as stream:
                    csv_ids = sorted((row.get("sector_id", ""), row.get("candidate_id", ""))
                                     for row in csv.DictReader(stream))
            finally:
                csv.field_size_limit(previous_limit)
            json_ids = sorted((str(row.get("sector_id", "")), str(row.get("candidate_id") or ""))
                for primary in recommendations if isinstance(primary, dict)
                for row in [primary, *(primary.get("alternatives") or [])])
            if csv_ids != json_ids:
                errors.append("sector recommendation JSON/CSV membership mismatch")
        except Exception as exc:
            errors.append(f"sector recommendation CSV validation failed: {exc}")
    if not zip_path.is_file() or not manifest_path.is_file():
        errors.append("scientific ZIP or run_manifest.json missing")
    else:
        try:
            with zipfile.ZipFile(zip_path) as archive:
                names = set(archive.namelist())
                if any(Path(name).suffix.lower() in {".sqlite", ".sqlite3", ".db"} for name in names):
                    errors.append("scientific ZIP contains a SQLite database")
                for name in (*key_files, "run_manifest.json"):
                    if name not in names:
                        errors.append(f"scientific ZIP missing {name}")
                        break
                    if name != "run_manifest.json" and hashlib.sha256(archive.read(name)).hexdigest() != hashes.get(name):
                        errors.append(f"scientific ZIP source hash mismatch: {name}")
                        break
                if "run_manifest.json" in names and archive.read("run_manifest.json") != manifest_path.read_bytes():
                    errors.append("scientific ZIP run_manifest differs from on-disk manifest")
        except (OSError, zipfile.BadZipFile) as exc:
            errors.append(f"scientific ZIP validation failed: {exc}")
    checks["export_consistency"] = not any("hash mismatch" in error or "ZIP" in error or "JSON/CSV" in error for error in errors)
    checks["card_evidence"] = not cards or bool(receipt and receipt.get("status") == "PASS")
    if errors:
        status = "NEEDS_REVIEW"
    elif not cards:
        status = "DIAGNOSTIC_COMPLETE"
    else:
        status = "PASS"
    return {"status": status, "checks": checks, "technical_failures": technical, "errors": errors,
            "regional_principle_coverage": regional_principle_coverage(ledger),
            "recommendation_card_count": len(cards), "foldx_pass_count": len(passed_foldx)}


def mutation_signature(value: Any) -> str | None:
    mutations = value.get("mutations") if isinstance(value, dict) else None
    if not isinstance(mutations, list) or not mutations:
        return None
    parsed = []
    for mutation in mutations:
        if isinstance(mutation, str):
            match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", mutation.strip().upper())
            if match:
                parsed.append(f"{match[1]}{int(match[2])}{match[3]}")
        elif isinstance(mutation, dict):
            notation = mutation.get("notation")
            if notation:
                match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", str(notation).strip().upper())
                if match:
                    parsed.append(f"{match[1]}{int(match[2])}{match[3]}")
            else:
                position = mutation.get("position", mutation.get("residue_position"))
                wt = mutation.get("wt", mutation.get("wild_type"))
                mutant = mutation.get("mutant", mutation.get("mutant_aa"))
                if position is not None and wt and mutant:
                    parsed.append(f"{str(wt).upper()}{int(position)}{str(mutant).upper()}")
    return "|".join(sorted(parsed)) if len(parsed) == len(mutations) else None


def _candidate_rows(analysis: dict, workdir: Path) -> tuple[list[dict], str]:
    panel = workdir / "panel_candidates.json"
    if panel.is_file():
        data = json.loads(panel.read_text(encoding="utf-8"))
        rows = data if isinstance(data, list) else data.get("candidates", data.get("panel_candidates", []))
        if isinstance(rows, list):
            return rows, "panel_candidates"
    return [row for row in analysis.get("designs", []) if isinstance(row, dict)], "displayed_designs"


def _region_members(rows: Any) -> list[set[int]]:
    result = []
    if not isinstance(rows, list):
        return result
    for row in rows:
        if not isinstance(row, dict):
            continue
        members = row.get("residues", row.get("members", row.get("positions", [])))
        parsed = set()
        for value in members if isinstance(members, list) else []:
            if isinstance(value, dict):
                value = value.get("position", value.get("residue_position", value.get("uniprot_pos")))
            try:
                parsed.add(int(value))
            except (TypeError, ValueError):
                continue
        if parsed:
            result.append(parsed)
    return result


def _chain_sequence(path: Path, chain_id: str) -> str | None:
    try:
        from ecd_epitope_scanner.structure_analysis import chain_sequence as get_sequence, load_structure
        parsed = load_structure(path)
        chain = next((chain for chain in parsed.model.get_chains() if str(chain.id) == str(chain_id)), None)
        return get_sequence(chain) if chain is not None else None
    except Exception:
        return None


def cross_compare_case(case: dict, state_job: dict) -> dict:
    baseline = case.get("baseline_evidence")
    result = state_job.get("result_evidence") or {}
    report = {"case_id": case["case_id"], "source_run_id": case["source_run_id"],
              "baseline_status": "MISSING", "mapping_status": "UNKNOWN", "shared_signatures": [],
              "new_reported_signatures": [], "legacy_scored_not_in_rerun": [], "same_combo_scores": [],
              "region_best_jaccard": None, "recommendation_comparison": None,
              "budget_and_sensitivity": None}
    if not baseline:
        report["reason"] = "No frozen historical baseline evidence is declared for this case."
        return report
    authority = baseline.get("authority_path")
    if authority and (not Path(authority).is_file() or sha256(Path(authority)) != baseline.get("authority_sha256")):
        report["baseline_status"] = "BASELINE_AUTHORITY_HASH_MISMATCH"
        report["reason"] = f"Historical acceptance evidence missing or changed: {authority}"
        return report
    report["legacy_summary"] = baseline.get("legacy_summary")
    old_scores, old_seq = {}, None
    for entry in baseline.get("paired_manifests", []):
        path = Path(entry["path"])
        if not path.is_file() or sha256(path) != entry["sha256"]:
            report["baseline_status"] = "BASELINE_HASH_MISMATCH"
            report["reason"] = f"Historical paired manifest missing or changed: {path}"
            return report
        data = json.loads(path.read_text(encoding="utf-8"))
        old_seq = data.get("execution_identity", {}).get("sequence", old_seq)
        for record in data.get("records", []):
            signature = mutation_signature(record)
            if not signature:
                continue
            values = record.get("ddg_runs_kcal_mol") or [row.get("ddg_kcal_mol") for row in record.get("paired_models", [])]
            values = [float(value) for value in values if value is not None]
            old_scores.setdefault(signature, []).extend(values)
    if not old_seq:
        old_dir = Path(case.get("source_run_dir", ""))
        old_structure = next((path for path in (old_dir / "structure.cif", old_dir / "structure.pdb") if path.is_file()), None)
        if old_structure:
            old_seq = _chain_sequence(old_structure, case["chain"])
    if old_scores:
        report["baseline_status"] = "PAIRED_FOLDX_EVIDENCE_AVAILABLE"
    elif baseline.get("legacy_summary"):
        report["baseline_status"] = "LEGACY_AGGREGATE_ONLY"
    new_path = Path(result["analysis"]) if result.get("analysis") else None
    if not new_path or not new_path.is_file():
        report["reason"] = "Rerun analysis is not available yet."
        return report
    analysis = json.loads(new_path.read_text(encoding="utf-8"))
    current_seq = _chain_sequence(Path(case["structure_path"]), case["chain"])
    if not current_seq:
        current_seq = (analysis.get("request") or {}).get("sequence") or analysis.get("sequence")
    if not current_seq and isinstance(analysis.get("protein"), dict):
        current_seq = analysis["protein"].get("sequence")
    if old_seq and current_seq:
        report["mapping_status"] = "SEQUENCE_IDENTICAL" if old_seq == current_seq else "SEQUENCE_MISMATCH"
    else:
        report["mapping_status"] = "SEQUENCE_UNAVAILABLE"
    workdir = new_path.parent
    candidates, candidate_source = _candidate_rows(analysis, workdir)
    new_scores = {}
    for candidate in candidates:
        signature = mutation_signature(candidate)
        if not signature:
            continue
        stability = candidate.get("stability") or {}
        new_scores[signature] = {"spurs": stability.get("spurs_multi_ddg_kcal_mol", candidate.get("spurs_multi_ddg_kcal_mol")),
            "foldx": stability.get("foldx_ddg_kcal_mol", candidate.get("foldx_ddg_kcal_mol")),
            "foldx_repeats": ((candidate.get("impact_evidence") or {}).get("stability_repeat_audit") or {}).get("foldx_repeat_values_kcal_mol", [])}
    report["candidate_source"] = candidate_source
    if report["mapping_status"] == "SEQUENCE_IDENTICAL":
        if old_scores:
            old_set, new_set = set(old_scores), set(new_scores)
            report["shared_signatures"] = sorted(old_set & new_set)
            report["new_reported_signatures"] = sorted(new_set - old_set)
            report["legacy_scored_not_in_rerun"] = sorted(old_set - new_set)
            report["same_combo_scores"] = [{"signature": signature, "legacy_foldx_slots": old_scores[signature],
                "legacy_foldx_median": __import__("statistics").median(old_scores[signature]),
                "rerun_foldx": new_scores[signature]["foldx"], "rerun_spurs": new_scores[signature]["spurs"],
                "rerun_foldx_repeats": new_scores[signature]["foldx_repeats"]}
                for signature in report["shared_signatures"]]
            report["signature_comparison_status"] = "COMPARED_TO_SCORED_LEGACY_COMBINATIONS"
        else:
            report["signature_comparison_status"] = "UNKNOWN_LEGACY_PER_COMBINATION_SCORES_UNAVAILABLE"
    else:
        report["reason"] = "Mutation-position comparisons require identical sequence provenance; legacy scoring signatures are not mapped across different or unknown sequences."
    regions = analysis.get("regions", [])
    report["rerun_region_members"] = [sorted(group) for group in _region_members(regions)]
    baseline_region_path = baseline.get("sector_recommendations_path") or baseline.get("analysis_path")
    baseline_regions = []
    if baseline_region_path and Path(baseline_region_path).is_file():
        base_data = json.loads(Path(baseline_region_path).read_text(encoding="utf-8"))
        baseline_regions = _region_members(base_data.get("regions", base_data.get("recommendations", [])))
    if baseline_regions and report["rerun_region_members"]:
        pairs = [(len(a & b) / len(a | b), i, j) for i, a in enumerate(baseline_regions)
                 for j, b in enumerate(report["rerun_region_members"]) if a | b]
        report["best_region_jaccard"] = max(pairs)[0] if pairs else None
        report["region_jaccard_status"] = "COMPUTED_FROM_SAVED_MEMBERSHIPS"
    else:
        report["best_region_jaccard"] = None
        report["region_jaccard_status"] = "UNKNOWN: legacy region membership unavailable"
    legacy_summary = baseline.get("legacy_summary", {})
    sector_path = result.get("sector_recommendations_v4")
    sector_data = json.loads(Path(sector_path).read_text(encoding="utf-8")) if sector_path and Path(sector_path).is_file() else None
    diagnostics_path = result.get("model_evidence_diagnostics_v1")
    diagnostics = json.loads(Path(diagnostics_path).read_text(encoding="utf-8")) if diagnostics_path and Path(diagnostics_path).is_file() else {}
    diagnostic_rows = diagnostics if isinstance(diagnostics, list) else diagnostics.get("candidates", diagnostics.get("records", []))
    report["recommendation_comparison"] = {"legacy_formal_count": legacy_summary.get("formal"),
        "legacy_main_panel_count": legacy_summary.get("main_panel"),
        "rerun_recommendation_status": (analysis.get("design_summary") or {}).get("recommendation_status"),
        "rerun_principle_summary": sector_data,
        "legacy_principle_winners": "UNKNOWN: historical sector recommendation file is unavailable"}
    report["budget_and_sensitivity"] = {"legacy_spurs_multi": legacy_summary.get("spurs_multi"),
        "legacy_foldx_unique": legacy_summary.get("foldx_unique"), "legacy_foldx_slots": legacy_summary.get("base_slots"),
        "legacy_repeat_protocol": "5 slots",
        "rerun_protocol": case.get("request_contract", {}),
        "rerun_ten_slot_sensitivity": [entry for entry in diagnostic_rows
            if entry.get("foldx_repeat_count") == 10 or entry.get("sampling_decision_sensitive") is not None]}
    return report


def publish_reports(state: dict, manifest: dict, state_path: Path) -> None:
    """Keep report failures separate from scientific case status."""
    try:
        write_reports(state, manifest)
    except Exception as exc:
        state["report_error"] = f"{type(exc).__name__}: {exc}"
    else:
        state.pop("report_error", None)
    atomic_json(state_path, state)


def write_reports(state: dict, manifest: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    state.setdefault("status", "RUNNING")
    comparisons = []
    rows = []
    post_statuses = [(job.get("post_validation") or {}).get("status") for job in state["jobs"].values()]
    acceptance_gate = ("PASS" if state["status"] == "COMPLETE" and len(state["jobs"]) == len(manifest["rerun_cases"]) and all(
        status in {"PASS", "DIAGNOSTIC_COMPLETE"} for status in post_statuses)
        else "NEEDS_REVIEW" if state["status"] in {"FINISHED_WITH_REVIEW", "FINISHED_WITH_FAILURES",
                                                  "HALTED_FOR_REVIEW", "INTERRUPTED", "STOPPED_BY_USER"}
        or "NEEDS_REVIEW" in post_statuses else "PENDING")
    for case in manifest["rerun_cases"]:
        job = state["jobs"].get(case["case_id"], {})
        source_dir = Path(case["source_run_dir"]) if case.get("source_run_dir") else Path()
        old_analysis = next((p for p in (source_dir / "analysis.json", source_dir / "results" / "analysis.json") if p.is_file()), None)
        result = job.get("result_evidence") or {}
        new_analysis = Path(result["analysis"]) if result.get("analysis") else None
        comparison = cross_compare_case(case, job)
        comparison.update(rerun_status=job.get("status", "NOT_RUN"),
            baseline_analysis=str(old_analysis) if old_analysis else None,
            rerun_analysis=str(new_analysis) if new_analysis and new_analysis.is_file() else None,
            baseline_provenance=case.get("provenance_status", case.get("provenance")))
        comparisons.append(comparison)
        rows.append({"case_id": case["case_id"], "name": case["name"], "source_run_id": case["source_run_id"],
                     "job_id": job.get("job_id", ""), "status": job.get("status", "NOT_RUN"),
                     "elapsed_seconds": job.get("elapsed_seconds", ""), "region_count": result.get("region_count", ""),
                     "design_count": result.get("design_count", ""),
                     "recommendation_status": result.get("recommendation_status", ""), "error": job.get("error", ""),
                     "job_status": job.get("job_status", job.get("status", "NOT_RUN")),
                     "post_validation": json.dumps(job.get("post_validation")),
                     "baseline_status": comparison["baseline_status"], "mapping_status": comparison["mapping_status"],
                     "shared_signatures": json.dumps(comparison["shared_signatures"]),
                     "new_reported_signatures": json.dumps(comparison["new_reported_signatures"]),
                     "legacy_scored_not_in_rerun": json.dumps(comparison["legacy_scored_not_in_rerun"]),
                     "same_combo_scores": json.dumps(comparison["same_combo_scores"]),
                     "best_region_jaccard": comparison.get("best_region_jaccard", ""),
                     "recommendation_comparison": json.dumps(comparison.get("recommendation_comparison")),
                     "budget_and_sensitivity": json.dumps(comparison.get("budget_and_sensitivity"))})
    atomic_json(OUT / "report.json", {"status": state["status"], "manifest_sha256": state["manifest_sha256"],
                "error": state.get("error", state.get("stop_reason")), "blocked_case_id": state.get("blocked_case_id"),
                "protocol": manifest["protocol"], "cases": state["jobs"], "cross_comparison": comparisons,
                "baseline_provenance": manifest["comparison_baselines"],
                "comparison_dimensions": manifest["comparison_dimensions"],
                "acceptance_gate": acceptance_gate})
    with (OUT / "report.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys()); writer.writeheader(); writer.writerows(rows)
    comparison_text = []
    for item in comparisons:
        comparison_text.extend((
            f"### {item['case_id']} vs {item['source_run_id']}",
            f"Baseline: `{item['baseline_status']}`; residue mapping: `{item['mapping_status']}`; region best Jaccard: `{item.get('best_region_jaccard', 'UNKNOWN')}` ({item.get('region_jaccard_status', 'not computed')}).",
            f"Shared signatures: `{json.dumps(item['shared_signatures'])}`",
            f"Newly reported signatures: `{json.dumps(item['new_reported_signatures'])}`",
            f"Legacy scored signatures absent from rerun report: `{json.dumps(item['legacy_scored_not_in_rerun'])}`",
            f"Same-combination SPURS/FoldX evidence: `{json.dumps(item['same_combo_scores'], ensure_ascii=False)}`",
            f"Recommendations: `{json.dumps(item.get('recommendation_comparison'), ensure_ascii=False)}`",
            f"Budgets and repeat sensitivity: `{json.dumps(item.get('budget_and_sensitivity'), ensure_ascii=False)}`", ""))
    protocol = manifest["protocol"]
    protocol_name = protocol.get("engine", protocol.get("name")) if isinstance(protocol, dict) else protocol
    (OUT / "report.md").write_text("# Ordinary adaptive rerun\n\n"
        f"Status: `{state['status']}`. Cases: {len(rows)}. Protocol: `{protocol_name}`.\n\n"
        f"Acceptance gate: `{acceptance_gate}`. Independent checks cover scientific protocol, budgets, paired evidence, cards, and export consistency. PD-1 and CTLA4 source differences are recorded in report.json; scores from different protocols must not be described as improvements.\n\n"
        "Regional principle coverage is evaluated separately from the technical acceptance gate:\n\n"
        + "\n".join(f"- {case_id}: {json.dumps((job.get('post_validation') or {}).get('regional_principle_coverage', {'status': 'UNKNOWN'}), ensure_ascii=False)}" for case_id, job in state['jobs'].items())
        + "\n\n## Cross-comparison\n\n"
        + "\n".join(comparison_text)
        + "\nSignature differences describe available reports and model evidence; they do not infer validation or efficacy. Historical region membership and per-combination SPURS remain unavailable where the baseline lacks them.\n\n"
        + "\n".join(f"- {r['case_id']}: {r['status']} ({r['job_id'] or 'no job'})" for r in rows) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="Verify frozen inputs and print cases without creating jobs")
    parser.add_argument("--receipt", type=Path, help="Acceptance receipt with status PASS and matching manifest_sha256")
    parser.add_argument("--resume", action="store_true", help="Resume from the saved state; terminal jobs are never retried")
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    cases = preflight(manifest)
    if args.dry_run:
        print(json.dumps({"status": "PREFLIGHT_PASS", "cases": [c["case_id"] for c in cases]}, indent=2))
        return
    if not args.receipt:
        raise SystemExit("Execution requires --receipt with PASS for this manifest SHA256")
    execute(manifest, args.receipt, resume=args.resume)


if __name__ == "__main__":
    main()
