"""Coverage-first ranking and review isolation, including a saved CTLA4 replay."""
import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import pytest

from ecd_epitope_scanner.joint_design_v4.engine import (
    _deconvolution_presets, _legacy_designs, _sector_recommendation_designs,
)
from ecd_epitope_scanner.joint_design_v4.identifiability import greedy_panel
from ecd_epitope_scanner.joint_design_v4.models import (
    CombinationStabilityEvidence, InferenceHypothesis, ResidueKey,
    V4Candidate, V4MutationOption, V4Panel, V4Sector,
)
from ecd_epitope_scanner.joint_design_v4.optimizer import select_panel
from ecd_epitope_scanner.joint_design_v4.sector_recommendations import build_sector_recommendations
from ecd_epitope_scanner.models import ProteinContext


def option(position, mutant="A"):
    return V4MutationOption(
        ResidueKey("0", "A", position, uniprot_pos=position, one_letter="D"),
        mutant, f"D{position}{mutant}", "charge", 0.0, 1.0, True,
    )


def candidate(identifier, positions, response=(0.1, 0.1), mutant="A"):
    return V4Candidate(
        identifier, "VS001", [option(p, mutant) for p in positions], order=len(positions),
        formal_eligible=True, stability=CombinationStabilityEvidence(tuple(), 0.0, 0.0, len(positions)),
        hypothesis_response={"H1": response[0], "H2": response[1]},
        covered_positions_by_scenario={"geodesic_0": list(positions)},
    )


def hypotheses():
    return [SimpleNamespace(hypothesis_id=identifier, discotope_prior=0.5) for identifier in ("H1", "H2")]


def budget(minimum=1, maximum=3):
    return dict(minimum_constructs=minimum, target_constructs=maximum,
                maximum_constructs=maximum, adaptive_budget=True)


def sector(identifier, positions):
    return V4Sector(identifier, list(positions), 1, 100.0, 100.0, 0.0)


def test_surface_area_gain_precedes_information_gain():
    small_information = candidate("small_high_information", [1, 2], (1.0, 0.0))
    large_cover = candidate("large_low_information", [3, 4], (0.05, 0.05))
    selected, _, decision = greedy_panel(
        [small_information, large_cover], hypotheses(), budget(maximum=1), 2,
        lambda *args: 0.0, surface_areas={1: 1, 2: 1, 3: 20, 4: 20},
    )
    assert selected == [large_cover]
    assert decision.budget_curve[0]["coverage_gain_area_a2"] == 40
    assert decision.budget_curve[0]["coverage_gain"] == pytest.approx(40 / 42)


def test_positive_cover_remains_selectable_when_net_objective_is_negative():
    item = candidate("high_risk", [1, 2])
    selected, _, decision = greedy_panel(
        [item], hypotheses(), budget(), 2, lambda *args: 100.0,
        surface_areas={1: 5, 2: 5},
    )
    assert selected == [item]
    assert decision.budget_curve[0]["net_gain"] < 0
    assert decision.budget_curve[0]["coverage_gain"] == 1


def test_identical_responses_at_different_direct_sites_are_not_deduplicated():
    left, right = candidate("left", [1, 2]), candidate("right", [3, 4])
    selected, _, decision = greedy_panel(
        [left, right], hypotheses(), budget(), 2, lambda *args: 0,
        surface_areas=dict.fromkeys(range(1, 5), 1),
    )
    assert {item.candidate_id for item in selected} == {"left", "right"}
    assert [row["coverage_gain"] for row in decision.budget_curve] == [0.5, 0.5]


def test_minimum_budget_can_fill_unique_legal_combinations_before_zero_cover_stop():
    candidates = [candidate("a", [1, 2]), candidate("b", [1, 2], mutant="W"),
                  candidate("c", [1, 2], mutant="I")]
    selected, _, decision = greedy_panel(
        candidates, hypotheses(), budget(minimum=2), 2, lambda *args: 0,
        surface_areas={1: 5, 2: 5},
    )
    assert len(selected) == 2 and decision.minimum_target_met
    assert decision.budget_curve[1]["coverage_gain"] == 0
    assert decision.stop_reason == "NO_POSITIVE_COVERAGE_GAIN"
    one, _, fixed = greedy_panel(candidates, hypotheses(), budget(maximum=1), 2,
                                 lambda *args: 0, surface_areas={1: 5, 2: 5})
    assert len(one) == 1 and fixed.stop_reason == "MAXIMUM_CONSTRUCT_BUDGET_REACHED"


def test_review_exploratory_disconnected_and_higher_order_caps_remain_isolated():
    valid = candidate("valid", [1, 2])
    review = candidate("review", [3, 4])
    review.qualification = "REVIEW_REQUIRED"
    mutation_review = candidate("single_review", [5, 6])
    mutation_review.mutations[0].review_required = True
    exploratory = candidate("exploratory", [7, 8])
    exploratory.search_role = "EXPLORATORY"
    four = candidate("four", [9, 10, 11, 12])
    disconnected = candidate("disconnected", [13, 14])
    disconnected.connected = False
    triple = candidate("triple", [15, 16, 17])
    selected, _, _ = greedy_panel(
        [review, mutation_review, exploratory, four, disconnected, triple, valid],
        hypotheses(), budget(), 0, lambda *args: 0,
        surface_areas={p: 1 if p <= 2 else 100 for p in range(1, 18)},
    )
    assert selected == [valid]


def test_select_panel_uses_only_display_sector_surface_as_denominator():
    inside = candidate("inside", [1, 2])
    outside = candidate("outside", [99, 100])
    panel, selected = select_panel(
        [outside, inside], sectors=[sector("VS001", [1, 2])],
        surface_areas={1: 5, 2: 5, 99: 1000, 100: 1000}, hypotheses=hypotheses(), panel_budget=budget(),
    )
    assert selected == [inside]
    assert panel.budget_decision.budget_curve[0]["coverage_gain"] == 1
    assert panel.coverage_denominator_area_a2 == 10
    assert panel.achieved_area_coverage["geodesic_0"] == 1


def test_region_design_projection_reuses_primary_and_does_not_promote_review():
    shared = candidate("primary", [1, 2])
    extra = candidate("additional_formal", [3, 4])
    review = candidate("review", [5, 6])
    review.formal_eligible = False
    review.formal_exclusion_reasons = ["SPURS_ONLY_REVIEW_REQUIRED"]
    review.stability.foldx_ddg_kcal_mol = None
    sectors = [sector("VS001", [1]), sector("VS002", [2]), sector("VS003", [3, 4]), sector("VS004", [5, 6])]
    protein = ProteinContext("TEST", "test", "test", "DDDDDD", 1, 6, 6, "test")
    primary = _legacy_designs([shared], sectors, protein, None, "FORMAL_RECOMMENDATION")
    primary_before = asdict(primary[0])
    panel = V4Panel(["primary"], {"STAGE1_PROBE": ["primary"]}, {}, [], "test")
    panel.sector_recommendations = build_sector_recommendations(
        [review, extra, shared], sectors, dict.fromkeys(range(1, 7), 10), {}, selected=[shared],
    )
    appended = _sector_recommendation_designs(panel, [review, extra, shared], primary, sectors, protein, None)
    rows = panel.sector_recommendations
    assert rows[0]["design_id"] is None and rows[1]["design_id"] is None
    assert rows[0]["alternatives"][0]["design_id"] == rows[1]["alternatives"][0]["design_id"] == primary[0].id
    assert len(appended) == 2
    assert asdict(primary[0]) == primary_before
    assert panel.selected_candidate_ids == ["primary"]
    formal, unqualified = appended
    assert formal.design_role == "SECTOR_RECOMMENDATION"
    assert formal.status == "PASSED_SECTOR_RECOMMENDATION"
    assert all(m.allowed for m in formal.mutations)
    assert unqualified.design_role == "SECTOR_REVIEW_CANDIDATE"
    assert unqualified.status == "NEAR_THRESHOLD_NOT_QUALIFIED"
    assert not any(m.allowed for m in unqualified.mutations)
    assert unqualified.foldx_diagnostics["failed_gates"] == review.formal_exclusion_reasons
    stage2 = _deconvolution_presets(primary + appended, [], protein, None)
    assert stage2 and {design.parent_design_id for design in stage2} == {primary[0].id}
    assert review.formal_eligible is False and review.stability.foldx_ddg_kcal_mol is None


def test_saved_ctla4_coverage_panel_and_region_recommendations_preserve_scientific_scores():
    run = Path(__file__).resolve().parents[1] / "workspace_data" / "5068b62450ec"
    if not run.exists():
        pytest.skip("Saved CTLA4 run is not available in this checkout")
    names = ["panel_candidates.json", "probe_sectors_v4.json", "surface_graph_v41.json",
             "inference_hypotheses_v411.json", "run_manifest.json"]
    content = {name: (run / name).read_bytes() for name in names}
    hashes = {name: hashlib.sha256(data).hexdigest() for name, data in content.items()}
    candidates = []
    for row in json.loads(content["panel_candidates.json"]):
        row["mutations"] = [V4MutationOption(**{**option_row, "residue": ResidueKey(**option_row["residue"])})
                            for option_row in row["mutations"]]
        if row["stability"]:
            row["stability"] = CombinationStabilityEvidence(**row["stability"])
        candidates.append(V4Candidate(**row))
    sectors = [V4Sector(**row) for row in json.loads(content["probe_sectors_v4.json"])]
    graph = json.loads(content["surface_graph_v41.json"])
    areas = {node["uniprot_pos"]: node["surface_area_a2"] for node in graph["nodes"]}
    saved_hypotheses = [InferenceHypothesis(**row) for row in json.loads(content["inference_hypotheses_v411.json"])]
    params = json.loads(content["run_manifest.json"])["parameters"]["joint_v4_ml"]
    scores_before = {
        item.candidate_id: (asdict(item.stability) if item.stability else None,
                            [asdict(mutation) for mutation in item.mutations], item.formal_eligible,
                            list(item.formal_exclusion_reasons), dict(item.hypothesis_response))
        for item in candidates
    }
    panel, selected = select_panel(
        [item for item in candidates if item.formal_eligible], sectors=sectors,
        surface_areas=areas, hypotheses=saved_hypotheses,
        panel_budget=budget(minimum=4, maximum=8),
        risk_penalty_weights=params["panel_optimizer"].get("risk_penalty_weights"),
    )
    rows = build_sector_recommendations(candidates, sectors, areas, params, selected)
    assert 4 <= len(selected) <= 8
    assert all(item.formal_eligible for item in selected)
    assert panel.budget_decision.minimum_target_met
    assert panel.budget_decision.coverage_basis == "geodesic_0_eligible_surface_area"
    assert len(rows) == 9
    assert [row["status"] for row in rows] == ["PASSED", "NO_COMBINATION", "NO_COMBINATION",
        "PASSED", "PASSED", "PASSED", "PASSED", "NO_COMBINATION", "NEAR_THRESHOLD"]
    assert all(row["candidate_id"] or row["alternatives"] or row["unavailable_principles"] for row in rows)
    assert all(item["formal_eligible"] for row in rows for item in row["alternatives"])
    assert all(row["ownership"] == "EXCLUSIVE" for row in rows if row["candidate_id"])
    for item in candidates:
        assert scores_before[item.candidate_id] == (
            asdict(item.stability) if item.stability else None,
            [asdict(mutation) for mutation in item.mutations], item.formal_eligible,
            list(item.formal_exclusion_reasons), dict(item.hypothesis_response),
        )
    assert hashes == {name: hashlib.sha256((run / name).read_bytes()).hexdigest() for name in names}
