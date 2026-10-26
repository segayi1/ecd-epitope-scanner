"""Prepare or run bounded FoldX completion for displayed V4 recommendations."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import sys
import uuid
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def make_plan(root: Path, audit_path: Path):
    audit = read_json(audit_path)
    jobs = []
    for audited in audit["jobs"]:
        job_id = audited["job_id"]
        source = root / "workspace_data" / job_id
        expected_hashes = audited["source_files_sha256"]
        for name, expected in expected_hashes.items():
            path = source / name
            if not path.is_file() or sha256(path) != expected:
                raise ValueError(f"source evidence changed: {job_id}/{name}")
        candidates = {item["candidate_id"]: item for item in read_json(source / "panel_candidates.json")}
        pairs_manifest = read_json(source / "paired_model_manifest.json")
        source_hashes = {**expected_hashes, **{
            name: sha256(source / name) for name in
            ("paired_model_manifest.json", "v4_manifest.json", "structure_for_energy.pdb")
        }}
        by_signature = {
            tuple(sorted(row.get("mutations", []))): row
            for row in pairs_manifest.get("records", [])
        }
        unique = {}
        for displayed in audited["pending_display_rows"]:
            candidate_id = displayed["candidate"]
            candidate = candidates[candidate_id]
            notation = sorted(item["notation"] for item in candidate["mutations"])
            identity = tuple(notation)
            evidence = candidate.get("impact_evidence", {}).get("stability_repeat_audit", {})
            slots = sorted(set(int(slot) for slot in evidence.get("foldx_repeat_slots", [])))
            if slots not in ([], [0, 1, 2]):
                raise ValueError(f"unexpected existing FoldX slots for {job_id}/{candidate_id}: {slots}")
            previous = unique.get(candidate_id)
            if previous and previous["existing_slots"] != slots:
                raise ValueError(f"conflicting repeated display evidence: {job_id}/{candidate_id}")
            if candidate_id not in unique:
                unique[candidate_id] = {
                    "candidate_id": candidate_id,
                    "mutations": notation,
                    "display_rows": [],
                    "existing_slots": slots,
                    "missing_slots": [slot for slot in range(5) if slot not in slots],
                    "existing_pair_record": by_signature.get(identity),
                    "spurs_multi_ddg_kcal_mol": candidate.get("stability", {}).get("spurs_multi_ddg_kcal_mol"),
                    "foldx_ddg_kcal_mol": candidate.get("stability", {}).get("foldx_ddg_kcal_mol"),
                    "source_files": source_hashes,
                }
            unique[candidate_id]["display_rows"].append({
                "sector": displayed["sector"], "role": displayed["role"],
                "status": displayed["status"], "failed_gates": displayed["failed_gates"],
            })
        candidate_rows = list(unique.values())
        if any(not row["missing_slots"] for row in candidate_rows):
            raise ValueError(f"audit includes already complete candidate in {job_id}")
        if any(row["existing_slots"] and row["existing_pair_record"] is None for row in candidate_rows):
            raise ValueError(f"original FoldX pair record missing for {job_id}")
        jobs.append({"job_id": job_id, "source_dir": str(source), "source_files_sha256": source_hashes,
                     "candidates": candidate_rows,
                     "candidate_count": len(candidate_rows),
                     "pair_slots": sum(len(row["missing_slots"]) for row in candidate_rows)})
    return {
        "schema_version": "1.0", "protocol": "displayed_recommendation_foldx_completion_v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_audit": str(audit_path), "source_audit_sha256": sha256(audit_path),
        "source_manifest": audit.get("source_manifest", str(root / "workspace_data/p0_p3_upgrade_20261016/rerun_manifest.json")),
        "new_model_slots_only": True, "original_jobs_and_evidence_are_read_only": True,
        "budget": {"new_unique_combinations": sum(row["candidate_count"] for row in jobs),
                   "new_pair_slots": sum(row["pair_slots"] for row in jobs)},
        "jobs": jobs,
    }


def append_event(path: Path, event: dict):
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(event, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def execute_job(root: Path, output_root: Path, job: dict, *, execution_id: str,
                ledger_path: Path, reserved: dict, budget: dict | None = None):
    from ecd_epitope_scanner.joint_design_v4.engine import _adapter, _prediction

    source = Path(job["source_dir"])
    for name, expected in job["source_files_sha256"].items():
        if sha256(source / name) != expected:
            raise ValueError(f"source evidence changed before execution: {job['job_id']}/{name}")
    v4 = read_json(source / "v4_manifest.json")
    analysis = read_json(source / "analysis.json")
    params = v4["parameters"]["foldx_confirmation"]
    structure = source / "structure_for_energy.pdb"
    if not structure.is_file():
        raise FileNotFoundError(structure)
    sequence = analysis["protein"]["ecd_sequence"]
    chain = analysis["structure_selection"]["selected_chain"]
    job_output = output_root / job["job_id"]
    job_output.mkdir(parents=True, exist_ok=False)
    grouped = defaultdict(list)
    for candidate in job["candidates"]:
        slots = candidate["missing_slots"]
        grouped[(slots[0], len(slots))].append(candidate)
    records, groups = [], []
    for (slot_start, run_count), candidates in sorted(grouped.items()):
        candidate_ids = [row["candidate_id"] for row in candidates]
        mutations = [row["mutations"] for row in candidates]
        slots = list(range(slot_start, slot_start + run_count))
        cap = budget if budget is not None else {"new_unique_combinations": 20, "new_pair_slots": 43}
        if (reserved["unique"] + len(candidates) > cap["new_unique_combinations"]
                or reserved["slots"] + len(candidates) * run_count > cap["new_pair_slots"]):
            raise RuntimeError("supplement budget exceeds the authorized combination / slot limits")
        reserved["unique"] += len(candidates)
        reserved["slots"] += len(candidates) * run_count
        append_event(ledger_path, {"event": "GROUP_RESERVED", "execution_id": execution_id,
            "job_id": job["job_id"], "candidate_ids": candidate_ids, "mutations": mutations,
            "sampling_slots": slots, "unique_combinations": len(candidates),
            "pair_slots": len(candidates) * run_count,
            "timestamp_utc": datetime.now(timezone.utc).isoformat()})
        settings = copy.deepcopy(params)
        settings.update(number_of_runs=run_count, sampling_slot_start=slot_start,
                        sampling_protocol="time_rotabase_v3", sampling_role="production",
                        sampling_cohort="ordinary_v4_sequential_v1", emit_paired_models=True)
        adapter = _adapter("foldx", settings, job_output / "cache")
        group_dir = job_output / f"slots_{slot_start}_{slot_start + run_count - 1}"
        try:
            result = _prediction(adapter, structure, chain, sequence, group_dir, "pdb", settings, None,
                "foldx", mutations)
            expected = {tuple(sorted(items)) for items in mutations}
            observed = {tuple(sorted(row.get("mutations", []))) for row in result.records}
            observed_slots = {tuple(sorted(row.get("mutations", []))): sorted(
                int(pair["sampling_slot"]) for pair in row.get("paired_models", [])) for row in result.records}
            if result.status != "success" or result.errors or observed != expected or any(
                    observed_slots.get(identity) != slots for identity in expected):
                raise RuntimeError(f"FoldX response incomplete or mismatched: {result.errors}")
            group_file = job_output / f"group_slots_{slot_start}_{slot_start + run_count - 1}.json"
            group_file.write_text(json.dumps({"records": result.records, "provenance": result.provenance}, indent=2), encoding="utf-8")
            append_event(ledger_path, {"event": "GROUP_COMPLETE", "execution_id": execution_id,
                "job_id": job["job_id"], "candidate_ids": candidate_ids, "sampling_slots": slots,
                "record_count": len(result.records), "raw_records_path": str(group_file),
                "raw_records_sha256": sha256(group_file), "timestamp_utc": datetime.now(timezone.utc).isoformat()})
        except BaseException as exc:
            append_event(ledger_path, {"event": "GROUP_FAILED", "execution_id": execution_id,
                "job_id": job["job_id"], "candidate_ids": candidate_ids, "sampling_slots": slots,
                "error": f"{type(exc).__name__}: {exc}", "retry": "NOT_AUTOMATIC",
                "timestamp_utc": datetime.now(timezone.utc).isoformat()})
            raise
        groups.append({"sampling_slot_start": slot_start, "number_of_runs": run_count,
                       "candidate_ids": candidate_ids,
                       "status": result.status, "errors": result.errors,
                       "provenance": result.provenance})
        records.extend(result.records)
    (job_output / "supplemental_records.json").write_text(
        json.dumps({"records": records, "groups": groups}, indent=2), encoding="utf-8")
    return {"job_id": job["job_id"], "output_dir": str(job_output), "groups": groups,
            "record_count": len(records), "status": "RAW_EVIDENCE_WRITTEN_NOT_REEXPORTED"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--audit", type=Path, default=ROOT / "workspace_data/recommendation_completion_20261017/prechange_display_audit.json")
    parser.add_argument("--output", type=Path, default=ROOT / "workspace_data/recommendation_completion_20261017/supplement_run_01")
    parser.add_argument("--execute", action="store_true", help="run only missing FoldX slots")
    args = parser.parse_args()
    plan = make_plan(args.root, args.audit)
    args.output.mkdir(parents=True, exist_ok=False)
    plan_path = args.output / "supplement_plan.json"
    plan_path.write_text(json.dumps(plan, indent=2), encoding="utf-8")
    if args.execute:
        from ecd_epitope_scanner.surface_perturbation.power_guard import prevent_sleep
        execution_id = str(uuid.uuid4())
        ledger = args.output / "execution_ledger.jsonl"
        script_hash = sha256(Path(__file__).resolve())
        reserved = {"unique": 0, "slots": 0}
        plan["execution"] = {"execution_id": execution_id, "script_sha256": script_hash,
            "status": "RUNNING", "pid": os.getpid(), "started_at_utc": datetime.now(timezone.utc).isoformat(),
            "ledger": str(ledger), "reserved_unique_combinations": 0, "reserved_pair_slots": 0}
        plan_path.write_text(json.dumps(plan, indent=2), encoding="utf-8")
        append_event(ledger, {"event": "EXECUTION_STARTED", "execution_id": execution_id,
            "script_sha256": script_hash, "pid": os.getpid(),
            "timestamp_utc": datetime.now(timezone.utc).isoformat()})
        receipts = []
        try:
            with prevent_sleep() as power:
                plan["execution"]["power_guard"] = power
                plan_path.write_text(json.dumps(plan, indent=2), encoding="utf-8")
                for job in plan["jobs"]:
                    receipts.append(execute_job(args.root, args.output / "raw", job,
                        execution_id=execution_id, ledger_path=ledger, reserved=reserved))
            plan["execution"].update(status="RAW_EVIDENCE_COMPLETE", finished_at_utc=datetime.now(timezone.utc).isoformat(),
                jobs=receipts, reserved_unique_combinations=reserved["unique"], reserved_pair_slots=reserved["slots"])
            append_event(ledger, {"event": "EXECUTION_COMPLETE", "execution_id": execution_id,
                "reserved_unique_combinations": reserved["unique"], "reserved_pair_slots": reserved["slots"],
                "timestamp_utc": datetime.now(timezone.utc).isoformat()})
        except BaseException as exc:
            plan["execution"].update(status="FAILED", finished_at_utc=datetime.now(timezone.utc).isoformat(),
                error=f"{type(exc).__name__}: {exc}", reserved_unique_combinations=reserved["unique"],
                reserved_pair_slots=reserved["slots"], retry="NOT_AUTOMATIC")
            append_event(ledger, {"event": "EXECUTION_FAILED", "execution_id": execution_id,
                "error": f"{type(exc).__name__}: {exc}", "reserved_unique_combinations": reserved["unique"],
                "reserved_pair_slots": reserved["slots"], "retry": "NOT_AUTOMATIC",
                "timestamp_utc": datetime.now(timezone.utc).isoformat()})
            plan_path.write_text(json.dumps(plan, indent=2), encoding="utf-8")
            raise
        plan_path.write_text(json.dumps(plan, indent=2), encoding="utf-8")
    print(json.dumps({"plan": str(plan_path), "candidate_count": plan["budget"]["new_unique_combinations"],
                      "new_pair_slots": plan["budget"]["new_pair_slots"],
                      "execute": args.execute}, indent=2))


if __name__ == "__main__":
    main()
