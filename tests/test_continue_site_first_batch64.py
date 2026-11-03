import importlib.util
from copy import deepcopy
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "continue_site_first_batch64", Path(__file__).resolve().parents[1] / "scripts/continue_site_first_batch64.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def evidence():
    state = {"status": "INTERRUPTED", "error": module.BOUNDARY_ERROR,
             "production_source_sha256": "old", "jobs": {
                 "current": {"job_id": "job1", "job_status": "COMPLETE", "post_validation": {"status": "PASS"}}}}
    activation = {"old_source_sha256": "old", "new_source_sha256": "new",
                  "current_case_id": "current", "completed_before_boundary": 0}
    return state, activation


def test_authorized_boundary_keeps_old_job_evidence_and_source_lineage():
    state, activation = evidence()
    before = deepcopy(state)
    result = module.continuation_state(state, activation, "new")
    assert state == before
    assert result["jobs"] == before["jobs"]
    assert result["source_transition"]["old_scheduler_job_ids"] == ["job1"]
    assert result["production_source_sha256"] == "new"
    assert "error" not in result


@pytest.mark.parametrize("change", ["failure", "running", "wrong_stop", "moved_boundary", "old_hash"])
def test_refuses_unrelated_failure_or_changed_boundary(change):
    state, activation = evidence()
    if change == "failure":
        state["jobs"]["current"]["post_validation"]["status"] = "NEEDS_REVIEW"
    elif change == "running":
        state["jobs"]["current"]["job_status"] = "RUNNING"
    elif change == "wrong_stop":
        state["error"] = "ToolError: failure"
    elif change == "moved_boundary":
        activation["completed_before_boundary"] = 1
    else:
        state["production_source_sha256"] = "unexpected"
    with pytest.raises(ValueError):
        module.continuation_state(state, activation, "new")


def test_refuses_another_source_change_after_activation():
    state, activation = evidence()
    with pytest.raises(ValueError, match="changed again"):
        module.continuation_state(state, activation, "newer")
