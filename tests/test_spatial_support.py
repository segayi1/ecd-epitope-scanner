import json
from pathlib import Path

import pytest

from ecd_epitope_scanner.joint_design_v4.models import CombinationStabilityEvidence, ResidueKey, V4Candidate, V4MutationOption, V4Sector
from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates, hypothesis_balanced_shortlist
from ecd_epitope_scanner.joint_design_v4.engine import _foldx_shortlist
from ecd_epitope_scanner.joint_design_v4.identifiability import greedy_panel
from ecd_epitope_scanner.joint_design_v4.regional_perturbation import PERTURBATION_SELECTOR_POLICY
from ecd_epitope_scanner.joint_design_v4.sector_recommendations import build_sector_recommendations
from ecd_epitope_scanner.joint_design_v4.spatial_support import (
    annotate_pointwise_perturbation_fields, annotate_spatial_support,
    spatial_shortlist, site_signature,
)


def option(position, mutant="A", risk=0.0):
    return V4MutationOption(ResidueKey("0", "A", position, uniprot_pos=position, one_letter="D"),
                            mutant, f"D{position}{mutant}", "charge", risk, 1.0, True)


def candidate(identifier, positions, risk=0.0, mutant="A"):
    return V4Candidate(identifier, "S", [option(p, mutant) for p in positions], order=len(positions),
                       structural_risk=risk, formal_eligible=True,
                       stability=CombinationStabilityEvidence((), .1, .1, len(positions), consensus_status="CONSENSUS_PASS"))


def chain():
    graph = {p: {} for p in range(1, 6)}
    for p in range(1, 5):
        graph[p][p + 1] = graph[p + 1][p] = 6.0
    return graph


def test_containment_zero_gain_even_with_nonidentical_neighborhoods():
    graph = {1: {2: 1, 3: 9}, 2: {1: 1}, 3: {1: 9}}
    pair = candidate("pair", [1, 2])
    annotate_spatial_support([pair], graph, {p: 10 for p in graph})
    nine = pair.geometry_support["scenarios"]["geodesic_9"]
    assert nine["marginal_area_a2"] == 0
    assert nine["mean_pair_overlap_fraction"] == pytest.approx(2 / 3)


def test_pointwise_field_keeps_physical_channel_separate_from_epitope_weight():
    item = candidate("field", [2, 3])
    for mutation in item.mutations:
        mutation.chemical_descriptors = {
            "charge_delta": 1., "volume_delta_a3": 50., "hydropathy_delta": 1.,
            "donor_delta": 1., "acceptor_delta": 1., "aromatic_delta": 0.,
        }
        mutation.impact_evidence = {"sidechain_exposure": {"value": .8}}
        mutation.epitope_evidence_status = "SCORED"
        mutation.epitope_prior = .9
    graph = chain()
    areas = dict.fromkeys(graph, 10.)
    annotate_spatial_support([item], graph, areas)
    annotate_pointwise_perturbation_fields([item], areas, {1: .9, 5: .1})
    field = item.impact_evidence["surface_perturbation_field"]
    assert field["status"] == "SCORED"
    assert field["epitope_status"] == "PARTIAL"
    assert field["physical_field"]["1"] > 0
    assert field["physical_field"]["5"] > 0
    assert field["physical_field"]["1"] == pytest.approx(
        field["epitope_weighted_field"]["1"] / 1.4)
    assert field["epitope_weighted_field"]["5"] == pytest.approx(field["physical_field"]["5"] * .6)
    assert field["epitope_weighted_field"]["2"] == pytest.approx(field["physical_field"]["2"])
    assert field["epitope_mapping_basis"] == "target_surface_residue_prior"


def test_selector_uses_incremental_point_field_area_when_utility_ties():
    left, right = candidate("left_field", [1, 2]), candidate("right_field", [4, 5])
    for item in (left, right):
        item.impact_evidence["selector_policy"] = PERTURBATION_SELECTOR_POLICY
        for mutation in item.mutations:
            mutation.chemical_descriptors = {
                "charge_delta": 1., "volume_delta_a3": 50., "hydropathy_delta": 1.,
                "donor_delta": 1., "acceptor_delta": 1., "aromatic_delta": 0.,
            }
            mutation.impact_evidence = {"sidechain_exposure": {"value": .8}}
            mutation.epitope_evidence_status = "SCORED"
            mutation.epitope_prior = .5
    graph = chain()
    areas = {1: 1., 2: 1., 3: 1., 4: 50., 5: 50.}
    annotate_spatial_support([left, right], graph, areas)
    selected, _, decision = greedy_panel(
        [left, right], [], {"minimum_constructs": 1, "target_constructs": 1,
        "maximum_constructs": 1, "adaptive_budget": False}, 0, lambda *args: 0.,
        surface_areas=areas)
    assert selected == [right]
    assert decision.budget_curve[0]["epitope_weighted_field_gain_area_a2"] > 0


def test_triple_gain_compares_best_exact_pair_not_single():
    triple = candidate("triple", [1, 2, 3])
    annotate_spatial_support([triple], {1: {2: 10}, 2: {1: 10, 3: 10}, 3: {2: 10}}, {1: 10, 2: 10, 3: 10})
    scenario = triple.geometry_support["scenarios"]["geodesic_12"]
    assert scenario["union_area_a2"] == scenario["best_lower_order_area_a2"] == 30
    assert scenario["marginal_area_a2"] == 0


def test_pareto_pool_preserves_coverage_risk_and_site_diversity():
    pool = [candidate("near1", [2, 3]), candidate("near2", [2, 3], mutant="W"),
            candidate("wide", [1, 5], risk=0.8)]
    annotate_spatial_support(pool, chain(), dict.fromkeys(chain(), 10))
    chosen = spatial_shortlist(pool, 2)
    assert chosen[0].candidate_id == "wide"
    assert len({site_signature(c) for c in chosen}) == 2
    assert spatial_shortlist(pool, 0) == []
    assert [c.candidate_id for c in spatial_shortlist(list(reversed(pool)), 2)] == [c.candidate_id for c in chosen]


def test_generation_and_foldx_keep_complementary_pairs_before_variants():
    graph = chain()
    areas = dict.fromkeys(graph, 10)
    sector = V4Sector("S", list(graph), 1, 50, 50, 0)
    options = [option(p, aa, 0.4 if p in {1, 5} else 0.0) for p in graph for aa in ("A", "W", "E")]
    params = {"orders": [2], "site_pool_per_sector": 5, "beam_width": 3, "shortlist_per_order": 3}
    pool = generate_candidates([sector], options, graph, params, surface_areas=areas)
    assert len(pool) == 3 and len({site_signature(c) for c in pool}) == 3
    assert (1, 5) in {site_signature(c) for c in pool}
    again = generate_candidates([sector], list(reversed(options)), graph, params, surface_areas=areas)
    assert [c.candidate_id for c in again] == [c.candidate_id for c in pool]
    assert all(c.geometry_support and not c.formal_eligible for c in pool)
    chosen = _foldx_shortlist(pool, 2)
    assert len(chosen) == 2 and len({site_signature(c) for c in chosen}) == 2
    assert chosen[0].geometry_support["mean_marginal_area_a2"] > 0


def test_site_cap_uses_geometric_diversity_before_enumeration():
    graph = chain()
    sector = V4Sector("S", list(graph), 1, 50, 50, 0)
    pool = generate_candidates([sector], [option(p) for p in graph], graph,
                               {"orders": [2], "site_pool_per_sector": 2}, surface_areas=dict.fromkeys(graph, 10))
    assert len(pool) == 1
    assert max(site_signature(pool[0])) - min(site_signature(pool[0])) >= 2


def test_foldx_hypothesis_quota_preserves_unique_site_sets_under_budget():
    from ecd_epitope_scanner.joint_design_v4.models import V4SurfaceHypothesis
    pool = [candidate("near1", [2, 3]), candidate("near2", [2, 3], mutant="W"), candidate("wide", [1, 5], 0.8)]
    for item in pool:
        item.targeted_hypotheses = ["H"]
    annotate_spatial_support(pool, chain(), dict.fromkeys(chain(), 10))
    hypotheses = [V4SurfaceHypothesis("H", dict.fromkeys(chain(), 10))]
    chosen = _foldx_shortlist(pool, 2, hypotheses=hypotheses)
    assert len(chosen) == 2 and len({site_signature(c) for c in chosen}) == 2
    assert chosen[0].candidate_id == "wide"
    diagnostics = {}
    selected = hypothesis_balanced_shortlist(pool, hypotheses, {"orders": [2], "global_shortlist_per_order": 2}, diagnostics)
    assert len(selected) == 2 and len({site_signature(c) for c in selected}) == 2


def test_area_based_global_pool_independent_of_display_partition_ownership():
    from ecd_epitope_scanner.joint_design_v4.models import V4SurfaceHypothesis
    graph = chain()
    areas = dict.fromkeys(graph, 10)
    hypotheses = [V4SurfaceHypothesis("small", {1: 10, 2: 10, 3: 10}),
                  V4SurfaceHypothesis("large", dict.fromkeys(graph, 10))]
    one = [V4Sector("A", list(graph), 1, 50, 50, 0)]
    split = [V4Sector("B", [1, 2], 1, 20, 20, 0), V4Sector("C", [3, 4, 5], 1, 30, 30, 0)]
    params = {"orders": [2], "global_shortlist_per_order": 3, "sector_resource_guarantees": False}
    options = [option(p, aa) for p in graph for aa in ("A", "W")]
    first = generate_candidates(one, options, graph, params, hypotheses, surface_areas=areas)
    second = generate_candidates(split, list(reversed(options)), graph, params, list(reversed(hypotheses)), surface_areas=areas)
    assert [c.candidate_id for c in first] == [c.candidate_id for c in second]
    assert all(c.geometry_support["target_scope_id"] == "full_eligible_surface" for c in first)


def test_regional_complementarity_never_overrides_formal_gate():
    pool = [candidate("near", [2, 3]), candidate("wide", [1, 5], risk=0.8)]
    annotate_spatial_support(pool, chain(), dict.fromkeys(chain(), 10))
    sector = V4Sector("S", list(chain()), 1, 50, 50, 0)
    rows = build_sector_recommendations(pool, [sector], dict.fromkeys(chain(), 10), {})
    assert rows[0]["candidate_id"] == "wide" and rows[0]["geometry_support"]
    pool[1].formal_eligible = False
    pool[1].formal_exclusion_reasons = ["COMBINATION_REJECT"]
    assert build_sector_recommendations(pool, [sector], dict.fromkeys(chain(), 10), {})[0]["candidate_id"] == "near"


def test_whole_region_coverage_precedes_gains_outside_region():
    pool = [candidate("near", [2, 3]), candidate("wide", [1, 5])]
    areas = dict.fromkeys(chain(), 10)
    annotate_spatial_support(pool, chain(), areas)
    assert pool[1].geometry_support["mean_marginal_area_a2"] > pool[0].geometry_support["mean_marginal_area_a2"]
    sector = V4Sector("S", [1, 2, 3], 1, 30, 30, 0)
    other = V4Sector("T", [4, 5], 1, 20, 20, 0)
    row = build_sector_recommendations(pool, [sector, other], areas, {})[0]
    assert row["candidate_id"] == "near"
    assert row["regional_geometry_support"]["mean_coverage_fraction"] == 1
    wide = build_sector_recommendations(pool[1:], [sector, other], areas, {})[0]
    assert wide["status"] == "NO_COMBINATION"
    cross = next(row for row in wide["alternatives"] if row["ownership"] == "CROSS_REGION")
    assert cross["regional_geometry_support"]["scenarios"]["geodesic_9"]["uncovered_positions"] == [3]
    assert cross["regional_geometry_support"]["scenarios"]["geodesic_9"]["coverage_fraction"] == pytest.approx(2 / 3)


def test_zero_region_area_and_missing_support_do_not_claim_full_coverage():
    pool = [candidate("c", [1, 2])]
    annotate_spatial_support(pool, chain(), dict.fromkeys(chain(), 0))
    sector = V4Sector("S", [1, 2], 1, 0, 0, 0)
    row = build_sector_recommendations(pool, [sector], dict.fromkeys(chain(), 0), {})[0]
    assert row["regional_geometry_support"]["mean_coverage_fraction"] is None
    pool[0].geometry_support = {}
    assert build_sector_recommendations(pool, [sector], {1: 1, 2: 1}, {})[0]["regional_geometry_support"] == {}


def test_invalid_area_rejected_and_missing_area_is_not_fake_zero():
    pool = [candidate("c", [1, 2])]
    annotate_spatial_support(pool, chain(), None)
    assert pool[0].geometry_support == {}
    with pytest.raises(ValueError, match="finite nonnegative"):
        annotate_spatial_support(pool, chain(), {1: float("nan")})


def test_saved_ctla4_metrics_match_audited_values_without_changing_scores():
    root = Path(__file__).resolve().parents[1] / "workspace_data" / "5068b62450ec"
    if not root.exists():
        pytest.skip("Saved CTLA4 run not available")
    raw = (root / "panel_candidates.json").read_bytes()
    graph_data = json.loads((root / "surface_graph_v4.json").read_text(encoding="utf-8"))
    areas = {node["uniprot_pos"]: node["surface_area_a2"] for node in graph_data["nodes"]}
    graph = {p: {} for p in areas}
    for edge in graph_data["edges"]:
        graph[edge["source"]][edge["target"]] = graph[edge["target"]][edge["source"]] = edge["distance_a"]
    pool = [candidate("K65A66", [65, 66]), candidate("L98D99", [98, 99]), candidate("T82E83", [82, 83])]
    annotate_spatial_support(pool, graph, areas)
    assert [c.geometry_support["scenarios"]["geodesic_12"]["marginal_area_a2"] for c in pool] == pytest.approx([135.40264, 0, 138.98406], abs=1e-4)
    assert all(row["marginal_area_a2"] == pytest.approx(0) for row in pool[1].geometry_support["scenarios"].values())
    assert (root / "panel_candidates.json").read_bytes() == raw
