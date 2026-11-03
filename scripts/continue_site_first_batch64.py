"""One-shot, fail-closed continuation after the authorized batch64 source boundary."""
from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_site_first_13_cases as launch

batch = launch.batch
BASE = ROOT / "workspace_data/site_first_13_cases_20261029"
OUT = ROOT / "workspace_data/site_first_13_cases_20261029_batch64"
ACTIVATION = ROOT / "workspace_data/site_first_spurs_batch64_20261029/activation.json"
BOUNDARY_ERROR = "RuntimeError: Production source changed during batch; stopping before the next case"


def continuation_state(state, activation, source_sha256):
    if state.get("status") != "INTERRUPTED" or state.get("error") != BOUNDARY_ERROR:
        raise ValueError("Original batch did not stop at the authorized source boundary")
    if state.get("production_source_sha256") != activation["old_source_sha256"]:
        raise ValueError("Original batch source does not match activation")
    if source_sha256 != activation["new_source_sha256"]:
        raise ValueError("Source changed again after batch64 activation")
    jobs = state.get("jobs", {})
    if activation["current_case_id"] not in jobs:
        raise ValueError("Boundary case is missing")
    for saved in jobs.values():
        if (saved.get("job_status") != "COMPLETE" or
                saved.get("post_validation", {}).get("status") not in {"PASS", "DIAGNOSTIC_COMPLETE"}):
            raise ValueError("A case is incomplete or failed acceptance; continuation refused")
    if len(jobs) != activation["completed_before_boundary"] + 1:
        raise ValueError("Case boundary moved; continuation refused")
    result = deepcopy(state)
    result.pop("error", None)
    result.update(status="RUNNING", production_source_sha256=source_sha256)
    result["source_transition"] = {
        "authorization": "User requested batch64 starting with the next case",
        "previous_state": str(BASE / "state.json"),
        "old_source_sha256": activation["old_source_sha256"],
        "new_source_sha256": source_sha256,
        "old_scheduler_job_ids": [job["job_id"] for job in jobs.values()],
        "new_refill_policy": "head_miss_round_robin_v1",
    }
    return result


def main():
    activation = json.loads(ACTIVATION.read_text(encoding="utf-8"))
    state = json.loads((BASE / "state.json").read_text(encoding="utf-8"))
    source = batch.production_source_sha256()
    next_state = continuation_state(state, activation, source)
    old_receipt = json.loads((BASE / "launcher_receipt.json").read_text(encoding="utf-8"))
    manifest = json.loads(batch.MANIFEST.read_text(encoding="utf-8"))
    batch.preflight(manifest)
    batch.assert_no_active_jobs(ROOT / "workspace_data/jobs.sqlite3")
    if (old_receipt["runner_sha256"] != batch.sha256(Path(batch.__file__)) or
            old_receipt["manifest_sha256"] != batch.sha256(batch.MANIFEST)):
        raise ValueError("Runner or frozen inputs changed")
    OUT.mkdir(exist_ok=False)
    batch.atomic_json(OUT / "previous_state.json", state)
    receipt = dict(old_receipt, production_source_sha256=source,
                   basis="73 regression checks and explicit next-case activation",
                   continuation_sha256=batch.sha256(Path(__file__)),
                   previous_receipt=str(BASE / "launcher_receipt.json"))
    batch.atomic_json(OUT / "launcher_receipt.json", receipt)
    batch.atomic_json(OUT / "state.json", next_state)
    batch.request_for = launch.request_for
    batch.validate_case_output = launch.validate
    batch.OUT = OUT / "report"
    # Report actual protocol; source_transition preserves the two scheduling segments.
    manifest["protocol"] = "ordinary_site_first_v2"
    batch.execute(manifest, OUT / "launcher_receipt.json", state_path=OUT / "state.json", resume=True)


if __name__ == "__main__":
    main()
