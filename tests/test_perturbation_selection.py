"""Coverage-first ranking and review isolation, including a saved CTLA4 replay."""
import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import pytest

from ecd_epitope_scanner.joint_design_v4.engine import (
    _deconvolution_presets, _legacy_designs,
    _sector_recommendation_designs,
)
from ecd_epitope_scanner.joint_design_v4.identifiability import greedy_panel
from ecd_epitope_scanner.joint_design_v4.regional_perturbation import PERTURBATION_SELECTOR_POLICY
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


def selector_candidate(identifier, sector_id, positions, prior=.5, exposure=.8):
    item = candidate(identifier, positions)
    item.sector_id = sector_id
    for mutation in item.mutations:
        mutation.epitope_evidence_status = "SCORED"
        mutation.epitope_prior = prior
        mutation.chemical_descriptors = {
            "charge_delta": 1., "volume_delta_a3": 50.,
            "hydropathy_delta": 1., "donor_delta": 1.,
            "acceptor_delta": 1., "aromatic_delta": 0.,
        }
        mutation.impact_evidence = {
            "sidechain_exposure": {"value": exposure},
            "selector_policy": PERTURBATION_SELECTOR_POLICY,
        }
    item.impact_evidence = {"selector_policy": PERTURBATION_SELECTOR_POLICY}
    return item


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


def test_selector_fills_each_sector_before_global_coverage_and_uses_epitope_weighted_utility():
    low = selector_candidate("a_low", "A", [1, 2], prior=.1)
    high = selector_candidate("a_high", "A", [1, 3], prior=.9)
    other = selector_candidate("b_other", "B", [4, 5], prior=.9)
    selected, _, decision = greedy_panel(
        [low, high, other], hypotheses(), budget(minimum=2, maximum=2), 0,
        lambda *args: 0., surface_areas={1: 1, 2: 100, 3: 1, 4: 1, 5: 1},
    )
    assert [item.candidate_id for item in selected] == ["a_high", "b_other"]
    assert all(row["sector_first"] for row in decision.budget_curve[:2])
    assert decision.selection_priority.startswith("sector_first_epitope_weighted")


def test_selector_covers_every_sector_before_a_twelfth_supplement():
    candidates = []
    expected = {f"VS{index:03d}" for index in range(1, 12)}
    for index in range(11):
        position = index * 2 + 1
        candidates.extend([
            selector_candidate(f"s{index}_primary", f"VS{index + 1:03d}", [position, position + 1], prior=.8),
            selector_candidate(f"s{index}_alternate", f"VS{index + 1:03d}", [position + 1, position + 2], prior=.4),
        ])
    selected, _, decision = greedy_panel(
        candidates, hypotheses(), budget(minimum=4, maximum=12), 0,
        lambda *args: 0., surface_areas=dict.fromkeys(range(1, 25), 1.),
    )
    assert {item.sector_id for item in selected[:11]} == expected
    assert all(row["sector_first"] for row in decision.budget_curve[:11])
    assert len(selected) == 12 and not decision.budget_curve[11]["sector_first"]


def test_selector_does_not_stop_before_available_positive_coverage():
    first = selector_candidate("first", "A", [1, 2], prior=1.)
    no_gain = selector_candidate("no_gain", "A", [1, 3], prior=.9)
    remaining = selector_candidate("remaining", "A", [4, 5], prior=.1)
    first.covered_positions_by_scenario = {"geodesic_0": [1, 2, 3]}
    no_gain.covered_positions_by_scenario = {"geodesic_0": [1, 2, 3]}
    selected, _, decision = greedy_panel(
        [first, no_gain, remaining], hypotheses(), budget(minimum=1, maximum=3), 12,
        lambda *args: 0., surface_areas=dict.fromkeys(range(1, 6), 1.),
    )
    assert [item.candidate_id for item in selected] == ["first", "remaining"]
    assert decision.budget_curve[1]["coverage_gain_area_a2"] == 2.


def test_new_panel_field_gain_precedes_triple_parent_margin(monkeypatch):
    from ecd_epitope_scanner.joint_design_v4 import spatial_support
    monkeypatch.setattr(spatial_support, "annotate_pointwise_perturbation_fields", lambda *args: None)
    first = selector_candidate("first", "A", [1, 2])
    duplicate = selector_candidate("duplicate", "A", [1, 3, 5])
    complementary = selector_candidate("complementary", "A", [4, 5])
    for item, field in ((first, {"1": 4.}), (duplicate, {"1": 2.}),
                        (complementary, {"4": 1.})):
        item.impact_evidence["surface_perturbation_field"] = {
            "physical_field": field, "epitope_weighted_field": field}
    duplicate.risk_components["triple_enhancement_delta"] = 100.
    selected, _, decision = greedy_panel(
        [first, duplicate, complementary], hypotheses(), budget(minimum=1, maximum=2), 12,
        lambda *args: 0., surface_areas=dict.fromkeys(range(1, 6), 1.))
    assert [item.candidate_id for item in selected] == ["first", "complementary"]
    assert decision.budget_curve[1]["physical_field_gain_area_a2"] == 1.


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
