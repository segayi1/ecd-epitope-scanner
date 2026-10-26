from pathlib import Path

from ecd_epitope_scanner.scripts.supplement_displayed_recommendations import make_plan


ROOT = Path(__file__).resolve().parents[1]


def test_supplement_plan_contains_only_missing_display_slots():
    plan = make_plan(ROOT, ROOT / "workspace_data/recommendation_completion_20261017/prechange_display_audit.json")
    assert plan["budget"] == {"new_unique_combinations": 20, "new_pair_slots": 43}
    candidates = [candidate for job in plan["jobs"] for candidate in job["candidates"]]
    assert len({(job["job_id"], candidate["candidate_id"])
                for job in plan["jobs"] for candidate in job["candidates"]}) == 20
    assert all(candidate["missing_slots"] == [3, 4] for candidate in candidates
               if candidate["existing_slots"] == [0, 1, 2])
    assert sum(len(candidate["missing_slots"]) for candidate in candidates) == 43
    assert all(candidate["existing_pair_record"] is None for candidate in candidates
               if not candidate["existing_slots"])


def test_supplement_rejects_zero_budget_before_worker(tmp_path):
    import pytest
    from ecd_epitope_scanner.scripts.supplement_displayed_recommendations import execute_job
    plan = make_plan(ROOT, ROOT / "workspace_data/recommendation_completion_20261017/prechange_display_audit.json")
    with pytest.raises(RuntimeError, match="authorized combination / slot limits"):
        execute_job(ROOT, tmp_path, plan["jobs"][0], execution_id="budget-test",
                    ledger_path=tmp_path / "ledger.jsonl", reserved={"unique": 0, "slots": 0},
                    budget={"new_unique_combinations": 0, "new_pair_slots": 0})
    assert not (tmp_path / "ledger.jsonl").exists()
