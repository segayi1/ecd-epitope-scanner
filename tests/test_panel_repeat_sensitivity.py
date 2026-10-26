import sys
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_panel_repeat_sensitivity import _choose_panel, _status


def test_audit_module_imports_project_package_when_launched_outside_repository(tmp_path):
    scripts = ROOT / "scripts"
    code = (
        "import sys; "
        f"sys.path.insert(0, {str(scripts)!r}); "
        "import audit_panel_repeat_sensitivity; "
        "from ecd_epitope_scanner.joint_design_v4.models import V4Sector; "
        "assert V4Sector is not None"
    )
    result = subprocess.run([sys.executable, "-c", code], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def _row(candidate_id, notation, spurs=0.1):
    import re
    from ecd_epitope_scanner.joint_design_v4.models import ResidueKey
    mutations = []
    for item in notation:
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", item)
        wt, pos, mutant = match.groups()
        mutations.append({"residue": ResidueKey("0", "A", int(pos), uniprot_pos=int(pos), one_letter=wt).__dict__,
            "mutant_aa": mutant, "notation": item, "chemistry_class": "other",
            "spurs_ddg_kcal_mol": 0.0, "score": 0.2, "eligible": True})
    return {"candidate_id": candidate_id, "sector_id": "R1", "mutations": mutations,
        "connected": True, "order": 2, "coverage_by_scenario": {"geodesic_0": 1.0},
        "covered_positions_by_scenario": {"geodesic_0": [1, 2]},
        "stability": {"mutations": list(notation), "spurs_multi_ddg_kcal_mol": spurs,
            "foldx_ddg_kcal_mol": 1.0, "order": 2, "consensus_status": "CONSENSUS_PASS"},
        "formal_eligible": True, "formal_exclusion_reasons": [], "risk_components": {},
        "recommendation_class": "STABILITY_FALLBACK", "impact_evidence": {}}


def test_fixed_completed_pool_selects_different_candidate_when_three_to_five_median_changes_rank():
    from ecd_epitope_scanner.joint_design_v4.models import V4Sector
    rows = [_row("c1", ["A1V", "A2L"]), _row("c2", ["A1L", "A2V"])]
    sectors = [V4Sector("R1", [1, 2], 1, 100.0, 100.0, 0.0)]
    params = {"stability_threshold_policy": "ordinary_v4_directional_v1", "acceptable_max_ddg": 1.0,
        "review_max_ddg": 2.0, "foldx_confirmation": {"combination_gate_kcal_mol": 2.0, "review_max_kcal_mol": 3.0},
        "panel_optimizer": {"higher_order_max_primary": 2, "risk_penalty_weights": {"stability": 1.0}},
        "max_primary_constructs": 1}
    budget = {"minimum_constructs": 1, "target_constructs": 1, "maximum_constructs": 1, "adaptive_budget": False}
    three, _, eligible_three = _choose_panel(rows,
        {("A1V", "A2L"): 1.8, ("A1L", "A2V"): 0.6}, set(), params, sectors, {1: 100, 2: 100}, [], budget)
    five, _, eligible_five = _choose_panel(rows,
        {("A1V", "A2L"): 0.4, ("A1L", "A2V"): 1.7}, set(), params, sectors, {1: 100, 2: 100}, [], budget)
    assert eligible_three == eligible_five == 2
    assert three != five
    assert _status(0.1, 1.7, params) == "PASS_PROXY_GATES"


def test_saved_nested_v4_stability_policy_and_thresholds_are_read_from_run_manifest():
    params = {
        "spurs_multi": {
            "stability_threshold_policy": "ordinary_v4_directional_v1",
            "acceptable_max_ddg": 1.0,
            "review_max_ddg": 2.0,
        },
        "foldx_confirmation": {"combination_gate_kcal_mol": 2.0, "review_max_kcal_mol": 3.0},
    }
    assert _status(1.0, 2.0, params) == "PASS_PROXY_GATES"
    assert _status(1.1, 2.0, params) == "REVIEW_PROXY_BAND"
    assert _status(2.1, 2.0, params) == "OUTSIDE_PROXY_BANDS"
