import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "write_p0_p3_comparison_figures.py"
spec = importlib.util.spec_from_file_location("p0figs", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_terminal_audit_label_accepts_report_status_and_keeps_running_pending():
    complete = {"status": "COMPLETE_ARTIFACTS_PRESENT", "audit": {"status": "AUDIT_COMPLETE"}}
    running = {"status": "RUNNING", "audit": {"status": "PASS"}}
    assert module.pending_tasks({"tasks": [complete]}) == []
    assert module.pending_tasks({"tasks": [complete, running]}) == [running]


def test_pending_cases_remain_blank_and_old_regions_are_parsed():
    report = {"tasks": [{
        "protein": "X", "status": "RUNNING", "new_job_id": "j1",
        "new_region_sizes_surface_nodes": None,
        "audit": {"status": "PENDING_REQUIRED_OUTPUTS"},
    }]}
    rows = [{"binder": "X", "job_status": "COMPLETE", "foldx_candidates": "60",
             "delivery_partition": json.dumps({"sizes": [7, 8, 9]})}]
    src = module.build_sources(report, rows)
    assert [r["surface_residues"] for r in src["regions"] if r["series"] == "原始流程/旧分区"] == [7, 8, 9]
    new = [r for r in src["regions"] if r["series"] == "新流程/新分区"]
    assert len(new) == 1 and new[0]["status"] == "PENDING" and new[0]["surface_residues"] == ""
    budget = src["budgets"][0]
    assert budget["old_attempt_slots_60x5"] == 300
    assert budget["new_total_slots_attempted"] is None


def test_initial_three_slot_count_and_missing_bias_are_distinct():
    report = {"tasks": [{
        "protein": "Y", "status": "COMPLETE_ARTIFACTS_PRESENT", "new_job_id": "j2",
        "new_sequential_foldx_budget": {"unique_combinations": 10, "initial_slots": 30,
                                         "supplemental_slots": 4, "total_slots": 34},
        "audit": {"status": "PASS", "counts": {"candidate_denominator": 20,
                "foldx_complete_five": 2, "formal_eligible": 2, "selected_panel_count": 1},
                  "selection_diversity_and_chemistry": {"main_panel": {
                      "unique_chemical_identities": 1, "unique_site_sets": 1,
                      "native_interface_mutation_fraction": None}}}
    }]}
    src = module.build_sources(report, [])
    by_stage = {r["stage"]: r["count"] for r in src["funnels"]}
    assert by_stage["scheduled initial 3 slots"] == 10
    assert by_stage["confirmed 5 slots"] == 2
    assert by_stage["main panel"] == 1
    assert src["biases"][-1]["native_interface_fraction"] is None
