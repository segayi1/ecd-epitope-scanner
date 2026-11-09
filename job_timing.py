"""Remaining wall time from comparable completed runs, not progress percentages."""
from __future__ import annotations

import datetime as dt
import json
import math
from pathlib import Path
from statistics import median


def timing_key(request):
    return tuple(request.get(key, default) for key, default in (
        ("input_type", "uniprot"), ("analysis_target", "antigen"), ("design_engine", "joint_v4_ml"),
        ("ml_mode", "strict"), ("surface_strategy", "legacy"), ("optimization_profile", "legacy"),
        ("execution_protocol", "legacy"))) + (
        bool(request.get("af3_manifest_path")), bool((request.get("metadata") or {}).get("disable_external_engines")))


def read_timing(workdir):
    path = Path(workdir) / "execution_trace.jsonl"
    try:
        if path.stat().st_size > 2_000_000:
            return None
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError):
        return None
    current = None
    for line in lines:
        try:
            row = json.loads(line)
            if not isinstance(row.get("execution_id"), str) or not row["execution_id"]:
                continue
            timestamp = dt.datetime.fromisoformat(row["timestamp_utc"])
            elapsed = float(row["elapsed_monotonic_seconds"])
            if timestamp.tzinfo is None or not math.isfinite(elapsed) or elapsed < 0:
                continue
            if row["event"] == "PIPELINE_STARTED":
                current = {"execution_id": row["execution_id"], "started_at": timestamp,
                           "stage": None, "stage_starts": {}, "elapsed_seconds": elapsed}
            if not current or row.get("execution_id") != current["execution_id"]:
                continue
            stage = row.get("stage")
            if stage and stage != current["stage"]:
                current["stage_starts"][stage] = elapsed
                current.update(stage=stage, stage_started_at=timestamp)
            current["elapsed_seconds"] = elapsed
        except (ValueError, TypeError, KeyError, AttributeError):
            continue
    return current


def remaining_time(current, references, now):
    stage = current.get("stage") if current else None
    remaining = []
    seen = {current["execution_id"]} if current else set()
    for reference in references:
        if reference["execution_id"] in seen or stage not in reference["stage_starts"]:
            continue
        seen.add(reference["execution_id"])
        value = reference["total_seconds"] - reference["stage_starts"][stage]
        if math.isfinite(value) and value > 0:
            remaining.append(value)
    result = {"status": "INSUFFICIENT_HISTORY", "stage": stage, "reference_count": len(remaining),
              "remaining_seconds": None, "range_seconds": None, "basis": "comparable_completed_stage_wall_time"}
    if not stage or len(remaining) < 3:
        return result
    elapsed = max(0., (now - current["stage_started_at"]).total_seconds())
    if elapsed >= max(remaining):
        return {**result, "status": "OVERRUN"}
    return {**result, "status": "ESTIMATED", "remaining_seconds": math.ceil(max(0., median(remaining) - elapsed) / 30) * 30,
            "range_seconds": [math.floor(max(0., min(remaining) - elapsed) / 30) * 30,
                              math.ceil((max(remaining) - elapsed) / 30) * 30]}
