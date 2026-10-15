from copy import deepcopy
from types import SimpleNamespace

from ecd_epitope_scanner.joint_design_v4.spatial_support import annotate_higher_order_regional_gain


def sector(identifier, positions):
    return SimpleNamespace(sector_id=identifier, residue_positions=positions)


def candidate(identifier, positions, order, covered=None, connected=True):
    mutations = [SimpleNamespace(residue=SimpleNamespace(uniprot_pos=position)) for position in positions]
    scenarios = None if covered is None else {
        f"geodesic_{radius}": {"covered_positions": values}
        for radius, values in zip((9, 12, 15), covered)
    }
    return SimpleNamespace(
        candidate_id=identifier,
        mutations=mutations,
        order=order,
        connected=connected,
        geometry_support={} if scenarios is None else {"scenarios": scenarios},
        risk_components={"existing": "preserved"},
    )


def setup_candidates():
    regions = [sector("VS001", [1, 2, 3, 4]), sector("VS002", [5, 6, 7, 8]),
               sector("VS003", [9, 10, 11, 12])]
    pair = candidate("pair-z", [1, 5], 2, [[1, 2], [1, 2, 3], [1, 2, 3, 5, 6]])
    pair_tie = candidate("pair-a", [2, 6], 2, [[1, 2], [1, 2, 3], [1, 2, 3, 5, 6]])
    triple = candidate("triple", [1, 5, 2], 3, [[1], [1, 5], [1, 2, 5, 6]])
    outside = candidate("outside", [9, 10], 2, [[9, 10, 11, 12]] * 3)
    return regions, [pair, pair_tie, triple, outside]


def test_uses_union_of_mutation_regions_and_stable_best_lower_order_baseline():
    regions, candidates = setup_candidates()
    original_geometry = deepcopy(candidates[2].geometry_support)
    annotate_higher_order_regional_gain(candidates, regions, dict.fromkeys(range(1, 13), 1.0))

    gain = candidates[2].risk_components
    # VS001 + VS002 form the eight-position scope; each mean is area weighted.
    assert gain["regional_gain_region_ids"] == ["VS001", "VS002"]
    assert gain["regional_candidate_coverage"] == (1 / 8 + 2 / 8 + 4 / 8) / 3
    assert gain["regional_lower_order_coverage"] == (2 / 8 + 3 / 8 + 5 / 8) / 3
    assert gain["regional_higher_order_gain"] == -1 / 8
    assert gain["regional_gain_baseline_candidate_id"] == "pair-a"
    assert gain["regional_gain_status"] == "COMPUTED"
    assert gain["regional_gain_interpretation"] == "counterfactual_geometry_only_not_candidate_qualification"
    assert candidates[2].geometry_support == original_geometry


def test_ordering_and_unrelated_region_candidates_do_not_change_gain():
    regions, candidates = setup_candidates()
    areas = dict.fromkeys(range(1, 13), 1.0)
    annotate_higher_order_regional_gain(candidates, regions, areas)
    expected = deepcopy(candidates[2].risk_components)

    regions2, candidates2 = setup_candidates()
    unrelated = candidate("unrelated-large", [9, 11, 12], 3, [[9, 10, 11, 12]] * 3)
    annotate_higher_order_regional_gain(list(reversed(candidates2 + [unrelated])), list(reversed(regions2)), areas)
    actual = next(item for item in candidates2 if item.candidate_id == "triple").risk_components
    assert actual == expected


def test_missing_baseline_geometry_and_area_are_not_treated_as_zero():
    region = [sector("VS001", [1, 2, 3])]
    triple = candidate("triple", [1, 2, 3], 3, [[1], [1, 2], [1, 2, 3]])
    annotate_higher_order_regional_gain([triple], region, {1: 1.0, 2: 1.0, 3: 1.0})
    assert triple.risk_components["regional_gain_status"] == "NO_LOWER_ORDER_BASELINE"
    assert triple.risk_components["regional_candidate_coverage"] == 2 / 3
    assert triple.risk_components["regional_higher_order_gain"] is None
    assert triple.risk_components["regional_lower_order_coverage"] is None

    incomplete_pair = candidate("pair", [1, 2], 2, None)
    annotate_higher_order_regional_gain([triple, incomplete_pair], region, {1: 1.0, 2: 1.0, 3: 1.0})
    assert triple.risk_components["regional_gain_status"] == "MISSING_LOWER_ORDER_GEOMETRY"
    assert triple.risk_components["regional_higher_order_gain"] is None

    no_geometry = candidate("no-geometry", [1, 2, 3], 3, None)
    annotate_higher_order_regional_gain([no_geometry], region, {1: 1.0, 2: 1.0, 3: 1.0})
    assert no_geometry.risk_components["regional_gain_status"] == "MISSING_CANDIDATE_GEOMETRY"
    assert no_geometry.risk_components["regional_candidate_coverage"] is None

    missing_area = candidate("missing-area", [1, 2, 3], 3, [[1], [1, 2], [1, 2, 3]])
    annotate_higher_order_regional_gain([missing_area], region, {1: 1.0, 2: None, 3: 1.0})
    assert missing_area.risk_components["regional_gain_status"] == "MISSING_AREA"
    assert missing_area.risk_components["regional_higher_order_gain"] is None


def test_unowned_mutation_position_prevents_partial_scope_scoring():
    triple = candidate("triple", [1, 2, 99], 3, [[1, 2], [1, 2], [1, 2]])
    annotate_higher_order_regional_gain([triple], [sector("VS001", [1, 2, 3])], {1: 1, 2: 1, 3: 1})
    assert triple.risk_components["regional_gain_status"] == "MISSING_REGION_OWNERSHIP"
    assert triple.risk_components["regional_candidate_coverage"] is None
