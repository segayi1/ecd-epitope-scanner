"""Rebuild independent V4 delivery bundles from saved scores plus a supplement run."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def signature(row):
    mutations = row.get("mutations", row.get("mutation_list", []))
    if isinstance(mutations, str):
        mutations = mutations.split(",")
    return tuple(sorted(str(item).strip() for item in mutations))


def merge_pairs(existing, additional):
    """Return a five-slot row while retaining all original paired records."""
    first = {signature(row): dict(row) for row in existing}
    for row in additional:
        key = signature(row)
        old = first.get(key)
        pairs = list(old.get("paired_models", [])) if old else []
        pairs.extend(row.get("paired_models", []))
        pairs.sort(key=lambda pair: int(pair["sampling_slot"]))
        slots = [int(pair["sampling_slot"]) for pair in pairs]
        values = [float(pair["ddg_kcal_mol"]) for pair in pairs]
        identities = [(pair.get("batch_path"), pair.get("group"), pair.get("run")) for pair in pairs]
        if slots != [0, 1, 2, 3, 4] or len(set(identities)) != 5 or not all(math.isfinite(value) for value in values):
            raise ValueError(f"incomplete, duplicate or non-finite paired evidence: {key}")
        merged = dict(old or {})
        merged.update(mutations=list(key), paired_models=pairs,
                      ddg_runs_kcal_mol=values, ddg=statistics.median(values), number_of_runs=5)
        first[key] = merged
    return first


def repeat_audit(pairs, spurs, foldx, spurs_parameters, foldx_parameters):
    from ecd_epitope_scanner.joint_design_v4.stability_repeat_diagnostics import foldx_repeat_gate_diagnostic
    slots = [int(pair["sampling_slot"]) for pair in pairs]
    values = [float(pair["ddg_kcal_mol"]) for pair in pairs]
    complete = slots == [0, 1, 2, 3, 4] and all(math.isfinite(value) for value in values)
    direction = ("DESTABILIZING_CONCORDANT" if (spurs > 0) == (foldx > 0) and spurs > 0 else
                 "NONDESTABILIZING_CONCORDANT" if (spurs > 0) == (foldx > 0) else "DISCORDANT")
    threshold = float(spurs_parameters.get("absolute_gap_diagnostic_kcal_mol",
        spurs_parameters.get("disagreement_gate_kcal_mol", 1.5)))
    content = [(pair.get("wt_sha256"), pair.get("mut_sha256")) for pair in pairs]
    duplicate_count = len(content) - len(set(content)) if all(all(item) for item in content) else None
    return {
        "foldx_pair_identities": [[pair.get("batch_path"), pair.get("group"), pair.get("run")] for pair in pairs],
        "foldx_duplicate_content_count": duplicate_count, "foldx_repeat_count": len(values),
        "foldx_repeat_slots": slots, "foldx_repeat_values_kcal_mol": values,
        "foldx_repeat_slot_values_kcal_mol": values,
        "foldx_phase1_median_kcal_mol": statistics.median([v for s, v in zip(slots, values) if s in (0, 1, 2)]),
        "foldx_phase2_median_kcal_mol": statistics.median([v for s, v in zip(slots, values) if s in (3, 4)]),
        "foldx_evidence_complete": complete, "direction_consistency": direction,
        "absolute_gap_kcal_mol": abs(spurs - foldx), "absolute_gap_threshold_kcal_mol": threshold,
        "absolute_gap_is_diagnostic_only": True,
        **foldx_repeat_gate_diagnostic(values, foldx_parameters.get("combination_gate_kcal_mol", 2.0), slots=slots),
    }


def update_target(candidate, merged, spurs_parameters, foldx_parameters):
    key = tuple(sorted(item["notation"] for item in candidate["mutations"]))
    row = merged[key]
    pairs = row["paired_models"]
    values = [float(pair["ddg_kcal_mol"]) for pair in pairs]
    median = statistics.median(values)
    spurs = float(candidate["stability"]["spurs_multi_ddg_kcal_mol"])
    directional = spurs_parameters.get("stability_threshold_policy") == "ordinary_v4_directional_v1"
    spurs_pass = spurs <= float(spurs_parameters.get("acceptable_max_ddg", 1.0) if directional else spurs_parameters.get("preferred_max_ddg", 0.5))
    spurs_review = spurs <= float(spurs_parameters.get("review_max_ddg", 2.0) if directional else spurs_parameters.get("review_max_ddg", 1.0))
    foldx_pass = median <= float(foldx_parameters.get("combination_gate_kcal_mol", 2.0))
    foldx_review = median <= float(foldx_parameters.get("review_max_kcal_mol", 3.0) if directional else foldx_parameters.get("combination_gate_kcal_mol", 2.0))
    gap = abs(spurs - median)
    if not directional and gap > float(spurs_parameters.get("absolute_gap_diagnostic_kcal_mol", spurs_parameters.get("disagreement_gate_kcal_mol", 1.5))):
        status = "MODEL_DISAGREEMENT"
    elif spurs_pass and foldx_pass:
        status = "CONSENSUS_PASS"
    elif spurs_review and foldx_review:
        status = "CONSENSUS_REVIEW"
    else:
        status = "COMBINATION_REJECT"
    audit = repeat_audit(pairs, spurs, median, spurs_parameters, foldx_parameters)
    candidate["stability"].update(foldx_ddg_kcal_mol=median, disagreement=gap, consensus_status=status)
    candidate["stability"].setdefault("provenance", {})["repeat_audit"] = audit
    candidate.setdefault("impact_evidence", {})["stability_repeat_audit"] = audit
    candidate["recommendation_risk_label"] = audit["recommendation_risk_label"]
    candidate["recommendation_risk_explanation"] = (
        "FoldX重复值跨越配置门限；中位数资格按现行合同保留，但结果对采样敏感，建议复核。"
        if audit["sampling_decision_sensitive"] is True else
        "FoldX重复证据缺失或含无效值，采样敏感性未知。"
        if audit["sampling_decision_diagnostic_status"] != "COMPLETE" else
        "本批FoldX重复值未跨越配置门限；请结合门限两侧槽数与中位数查看。"
    )
    reasons = set(candidate.get("formal_exclusion_reasons", [])) - {
        "FOLDX_REPEAT_EVIDENCE_INCOMPLETE", "COMBINATION_FOLDX_MISSING", "SPURS_ONLY_REVIEW_REQUIRED",
        "SPURS_ONLY_REVIEW_BAND", "FOLDX_ONLY_LEGACY", "NOT_SCORED",
    }
    if status != "CONSENSUS_PASS":
        reasons.add(status)
    candidate["formal_exclusion_reasons"] = sorted(reasons)
    candidate["formal_eligible"] = not reasons
    return row


def validate_completed_supplement(raw_root, plan):
    """Refuse to consume a partial run, even when some group files exist."""
    execution = plan.get("execution") or {}
    if execution.get("status") != "RAW_EVIDENCE_COMPLETE":
        raise ValueError("supplement execution is not RAW_EVIDENCE_COMPLETE")
    execution_id = execution.get("execution_id")
    ledger_path = raw_root / "execution_ledger.jsonl"
    if not execution_id or not ledger_path.is_file():
        raise ValueError("supplement execution identity or durable ledger is missing")
    events = [json.loads(line) for line in ledger_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    events = [event for event in events if event.get("execution_id") == execution_id]
    if any(event.get("event") in {"GROUP_FAILED", "EXECUTION_FAILED"} for event in events):
        raise ValueError("supplement ledger contains a failed group or execution")
    reserved = [event for event in events if event.get("event") == "GROUP_RESERVED"]
    completed = [event for event in events if event.get("event") == "GROUP_COMPLETE"]
    unique_total = sum(int(event.get("unique_combinations", 0)) for event in reserved)
    slot_total = sum(int(event.get("pair_slots", 0)) for event in reserved)
    if (unique_total, slot_total) != (int(plan["budget"]["new_unique_combinations"]),
                                       int(plan["budget"]["new_pair_slots"])):
        raise ValueError("supplement durable reservations do not match the plan budget")
    if (unique_total, slot_total) != (int(execution.get("reserved_unique_combinations", -1)),
                                       int(execution.get("reserved_pair_slots", -1))):
        raise ValueError("supplement execution receipt does not match durable reservations")
    key = lambda event: (event.get("job_id"), tuple(sorted(event.get("candidate_ids", []))),
                         tuple(int(slot) for slot in event.get("sampling_slots", [])))
    reserved_keys = [key(event) for event in reserved]
    completed_keys = [key(event) for event in completed]
    if len(reserved_keys) != len(set(reserved_keys)) or sorted(reserved_keys) != sorted(completed_keys):
        raise ValueError("not every reserved FoldX group has exactly one matching completion event")
    if not any(event.get("event") == "EXECUTION_COMPLETE" for event in events):
        raise ValueError("supplement ledger has no EXECUTION_COMPLETE event")
    completed_paths = set()
    project_root = Path(__file__).resolve().parents[1]
    for event in completed:
        path = Path(event.get("raw_records_path", ""))
        if not path.is_absolute():
            path = project_root / path
        path = path.resolve()
        if not path.is_relative_to(raw_root.resolve()) or not path.is_file():
            raise ValueError("completed group record is missing or escapes the supplement directory")
        if digest(path) != event.get("raw_records_sha256"):
            raise ValueError(f"completed group record hash mismatch: {path}")
        completed_paths.add(path)
    observed_paths = {path.resolve() for path in raw_root.glob("raw/*/group_slots_*.json")}
    if completed_paths != observed_paths:
        raise ValueError("raw group files do not exactly match ledger completions")
    return {"execution_id": execution_id, "completed_groups": len(completed),
            "reserved_unique_combinations": unique_total, "reserved_pair_slots": slot_total}


def write_card_evidence_receipt(output):
    """Bind the receipt to final output inputs after analysis paths are fixed."""
    analysis_path = output / "analysis.json"
    analysis = read(analysis_path)
    files = analysis.setdefault("files", {})
    files["paired_model_manifest.json"] = "paired_model_manifest.json"
    files["card_evidence_completeness.json"] = "card_evidence_completeness.json"
    write(analysis_path, analysis)
    input_names = ("analysis.json", "panel_candidates.json", "sector_recommendations_v4.json")
    input_hashes = {name: digest(output / name) for name in input_names}
    from ecd_epitope_scanner.joint_design_v4.impact_completeness import card_evidence_completeness
    receipt = card_evidence_completeness(read(output / "sector_recommendations_v4.json"), input_hashes)
    receipt_path = output / "card_evidence_completeness.json"
    write(receipt_path, receipt)
    run_manifest = read(output / "run_manifest.json")
    manifest_files = run_manifest.setdefault("files", {})
    manifest_files["paired_model_manifest.json"] = digest(output / "paired_model_manifest.json")
    manifest_files[receipt_path.name] = digest(receipt_path)
    write(output / "run_manifest.json", run_manifest)
    if digest(analysis_path) != input_hashes["analysis.json"]:
        raise AssertionError("analysis changed after completeness receipt hashes were captured")
    return receipt


def copy_refinement_manifests(source, output):
    """Keep active refinement manifests in replayed deliveries and their bundle inventory."""
    copied = []
    for path in sorted((source / "foldx_refinement").glob("**/paired_model_manifest.json")):
        relative = path.relative_to(source)
        if {"prior_attempts", "paired_manifest_history"}.intersection(relative.parts):
            continue
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and digest(target) != digest(path):
            raise ValueError("Existing refinement evidence differs: " + relative.as_posix())
        if not target.exists():
            shutil.copy2(path, target)
        copied.append(relative.as_posix())
    if copied:
        analysis_path = output / "analysis.json"
        analysis = read(analysis_path)
        analysis.setdefault("files", {}).update({name: name for name in copied})
        write(analysis_path, analysis)
        manifest_path = output / "run_manifest.json"
        manifest = read(manifest_path)
        manifest.setdefault("files", {}).update({name: digest(output / name) for name in copied})
        write(manifest_path, manifest)
    return copied


def prepare_delivery(source, raw_job, output, audit_job, raw_plan):
    if output.exists():
        raise FileExistsError(output)
    expected_hashes = {**audit_job["source_files_sha256"], **raw_plan["source_files_sha256"]}
    for name, expected in expected_hashes.items():
        if digest(source / name) != expected:
            raise ValueError(f"original source changed: {source / name}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staged = output.parent / (output.name + "_staging")
    if staged.exists():
        raise FileExistsError(staged)
    shutil.copytree(source, staged)
    candidates = read(staged / "panel_candidates.json")
    by_id = {row["candidate_id"]: row for row in candidates}
    old_manifest = read(staged / "paired_model_manifest.json")
    existing = {signature(row): dict(row) for row in old_manifest.get("records", [])}
    supplementary = []
    for path in sorted((raw_job).glob("group_slots_*.json")):
        supplementary.extend(read(path).get("records", []))
    candidate_plan = {row["candidate_id"]: row for row in raw_plan["candidates"]}
    expected_ids = set(candidate_plan)
    observed_ids = {tuple(sorted(row["mutations"])) for row in supplementary}
    if observed_ids != {tuple(row["mutations"]) for row in raw_plan["candidates"]}:
        raise ValueError(f"raw FoldX identities do not match the supplement plan for {raw_plan['job_id']}")
    merged = merge_pairs(list(existing.values()), supplementary)
    v4 = read(staged / "v4_manifest.json")
    spurs_parameters = v4["parameters"]["spurs_multi"]
    foldx_parameters = v4["parameters"]["foldx_confirmation"]
    for candidate_id in expected_ids:
        candidate = by_id[candidate_id]
        row = update_target(candidate, merged, spurs_parameters, foldx_parameters)
        existing[signature(row)] = row
    old_manifest["records"] = sorted(existing.values(), key=lambda row: signature(row))
    old_manifest.setdefault("supplemental_runs", []).append({
        "raw_output_dir": str(raw_job), "plan_job_id": raw_plan["job_id"],
        "candidate_ids": sorted(expected_ids), "new_pair_slots": raw_plan["pair_slots"],
        "records_sha256": {path.name: digest(path) for path in sorted(raw_job.glob("group_slots_*.json"))},
    })
    write(staged / "panel_candidates.json", candidates)
    write(staged / "paired_model_manifest.json", old_manifest)
    replay = ROOT / "scripts" / "replay_v4_selection.py"
    subprocess.run([sys.executable, str(replay), str(staged), str(output)], check=True, cwd=ROOT.parent)
    write(output / "paired_model_manifest.json", old_manifest)
    copy_refinement_manifests(staged, output)
    write_card_evidence_receipt(output)
    from ecd_epitope_scanner.exporting import finalize_scientific_bundle
    finalize_scientific_bundle(output)
    for name, expected in expected_hashes.items():
        if digest(source / name) != expected:
            raise ValueError(f"original source changed during replay: {source / name}")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output", "--output-root", dest="output_root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, default=ROOT / "workspace_data")
    args = parser.parse_args()
    plan = read(args.raw_root / "supplement_plan.json")
    audit = read(ROOT / "workspace_data/recommendation_completion_20261017/prechange_display_audit.json")
    audits = {row["job_id"]: row for row in audit["jobs"]}
    plans = {row["job_id"]: row for row in plan["jobs"]}
    if set(plans) != set(audits) or plan["budget"]["new_unique_combinations"] != 20 or plan["budget"]["new_pair_slots"] != 43:
        parser.error("Supplement plan does not match the audited 20 candidate / 43 slot contract")
    try:
        validate_completed_supplement(args.raw_root, plan)
    except (OSError, KeyError, TypeError, ValueError) as exc:
        parser.error(str(exc))
    if args.output_root.exists():
        parser.error("Use a new output root; existing deliveries are preserved")
    args.output_root.mkdir(parents=True)
    for job_id in sorted(plans):
        source = args.source_root / job_id
        prepare_delivery(source, args.raw_root / "raw" / job_id,
            args.output_root / job_id, audits[job_id], plans[job_id])
    print(json.dumps({"status": "REEXPORTED", "output_root": str(args.output_root),
        "jobs": sorted(plans), "additional_pair_slots": 43}, indent=2))


if __name__ == "__main__":
    main()
