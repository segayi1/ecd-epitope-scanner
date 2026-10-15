from itertools import combinations
import json

from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates, hypothesis_balanced_shortlist, _distances
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, V4Candidate, V4MutationOption, V4Sector, V4SurfaceHypothesis
from ecd_epitope_scanner.joint_design_v4.spatial_support import candidate_region_membership, site_signature, regional_support


def option(p, aa="A"):
    return V4MutationOption(ResidueKey("0", "A", p, uniprot_pos=p, one_letter="D"), aa,
                            f"D{p}{aa}", "charge", 0.0, 1.0, True)


def fixture():
    graph = {p: {q: 3.0 for q in range(1, 9) if p != q} for p in range(1, 9)}
    sectors = [V4Sector("VS001", [1, 2, 3, 4], 1, 40, 40, 0),
               V4Sector("VS002", [5, 6, 7, 8], 1, 40, 40, 0)]
    return graph, sectors, dict.fromkeys(graph, 10.0)


def test_membership_uses_actual_positions_and_unmapped_is_not_exclusive():
    _, sectors, _ = fixture()
    c = V4Candidate("x", "wrong", [option(1), option(2)], order=2)
    assert candidate_region_membership(c, sectors)["exclusive_region_id"] == "VS001"
    c.mutations = [option(1), option(6)]
    assert candidate_region_membership(c, sectors)["assignment"] == "CROSS_REGION"
    c.mutations = [option(1), option(99)]
    assert candidate_region_membership(c, sectors)["exclusive_region_id"] is None


def test_supplement_finds_regions_absent_from_hypotheses_and_order_stable():
    graph, sectors, areas = fixture()
    hypotheses = [V4SurfaceHypothesis("H1", {1: 10, 2: 10}), V4SurfaceHypothesis("H2", {1: 10, 2: 10, 3: 10})]
    params = {"orders": [2, 3], "global_shortlist_per_order": 1}
    audit = {}
    pool = generate_candidates(sectors, [option(p, aa) for p in graph for aa in ("A", "W")],
                               graph, params, hypotheses, audit, areas)
    for sector in sectors:
        exclusive = [c for c in pool if candidate_region_membership(c, sectors)["exclusive_region_id"] == sector.sector_id]
        assert len({site_signature(c) for c in exclusive if c.order == 2}) >= 4
        assert len({site_signature(c) for c in exclusive if c.order == 3}) >= 2
    assert all(row["supplemental_window"] for row in audit["sector_resources"])
    assert json.loads(json.dumps(audit)) == audit
    again = generate_candidates(list(reversed(sectors)), [option(p, aa) for p in reversed(graph) for aa in ("W", "A")],
                                graph, params, list(reversed(hypotheses)), surface_areas=areas)
    assert [c.candidate_id for c in pool] == [c.candidate_id for c in again]
    assert all(row["exclusive_base_expansion"] <= (8 if row["order"] == 2 else 4) for row in audit["global_order_pools"])


def test_insufficient_sites_reports_deficit_and_disabled_path_unchanged():
    graph, sectors, areas = fixture()
    audit = {}
    pool = generate_candidates(sectors, [option(1), option(2)], graph, {"orders": [2, 3]}, diagnostics=audit, surface_areas=areas)
    assert len(pool) == 1
    assert audit["sector_resources"][1]["quota_deficits"] == {"2": 4, "3": 2}
    legacy = generate_candidates(sectors, [option(1), option(2)], graph,
                                 {"orders": [2, 3], "sector_resource_guarantees": False}, surface_areas=areas)
    assert len(legacy) == 1 and legacy[0].candidate_id.startswith("V4_VS001")


def test_coverage_trigger_not_a_gate_or_unbounded_capacity():
    graph, sectors, areas = fixture()
    # Sparse geometric neighborhoods leave uncovered area even with enough site sets.
    graph = {p: {q: 20.0 for q in graph if p != q} for p in graph}
    audit = {}
    pool = generate_candidates(sectors, [option(p) for p in graph], graph,
                               {"orders": [2], "global_shortlist_per_order": 1}, diagnostics=audit, surface_areas=areas)
    assert all(row["initial_exclusive_site_sets"]["2"] >= 4 and row["supplemental_window"] for row in audit["sector_resources"])
    assert len(pool) <= 9
    assert all(not c.formal_eligible for c in pool)
    row = regional_support(pool[0], sectors[0], areas, graph=graph)
    assert row["scenarios"]["geodesic_9"]["target_uncovered_positions"]


def test_same_sites_variants_do_not_satisfy_base_quota():
    graph, sectors, areas = fixture()
    pool = [V4Candidate(f"p{i}", "old", [option(1, aa), option(2, aa)], order=2)
            for i, aa in enumerate(("A", "W", "E", "F"))]
    pool.extend(V4Candidate(f"q{p}{q}", "old", [option(p), option(q)], order=2)
                for p, q in combinations([1, 2, 3, 4], 2) if (p, q) != (1, 2))
    distances = {p: _distances(graph, p, set(graph)) for p in graph}
    chosen = hypothesis_balanced_shortlist(pool, [], {"orders": [2], "global_shortlist_per_order": 1},
                                          sectors=sectors, surface_areas=areas, distances=distances)
    assert len({site_signature(c) for c in chosen}) >= 4
