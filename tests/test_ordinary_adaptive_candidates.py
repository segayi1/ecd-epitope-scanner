import json
from types import SimpleNamespace as NS

from ecd_epitope_scanner.joint_design_v4 import ordinary_adaptive
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, V4Sector
from ecd_epitope_scanner.joint_design_v4.ordinary_adaptive import (
    PRINCIPLES,
    build_candidate_queues,
    strict_triple_upgrade,
    triple_upgrade_for_principle,
)
from ecd_epitope_scanner.joint_design_v4.regional_perturbation import AXES, candidate_impact


def _single(position, mutant, *, supported=False, exposure=.5, descriptor=1.0):
    residue = ResidueKey("m", "A", position, uniprot_pos=position, one_letter="A")
    return NS(residue=residue, mutant_aa=mutant, spurs_ddg_kcal_mol=.1,
              eligible=True, review_required=False, hard_exclusions=[], soft_warnings=[],
              chemistry_class="test", stability_band="PASS", chemical_descriptors={
                  **{axis: descriptor for axis in AXES}}, impact_evidence={
                      "sidechain_exposure": {"value": exposure, "status": "COMPUTED"}},
              epitope_prior=.95 if supported else .1, epitope_supported=supported)


def _inputs(count=5, supported_positions=()):
    sequence = "A" * count
    singles = [_single(position, mutant, supported=position in supported_positions,
                       descriptor=float(position))
               for position in range(1, count + 1) for mutant in ("D", "E")]
    graph = {position: {other: float(abs(position - other)) for other in range(1, count + 1)
                        if other != position} for position in range(1, count + 1)}
    sector = V4Sector("r1", list(range(1, count + 1)), 1, 100.0, 100.0, .5)
    areas = {position: 10.0 for position in range(1, count + 1)}
    return singles, [sector], graph, areas, sequence


def test_enumerates_past_legacy_double_triple_and_parent_caps():
    store = build_candidate_queues(*_inputs())
    rows = store.db.execute("SELECT metadata FROM candidates").fetchall()
    counts = {order: sum(ordinary_adaptive._decode_metadata(row[0])["order"] == order for row in rows) for order in (2, 3)}
    assert counts[2] > 20
    assert counts[3] > 10
    assert store.status("r1", "mechanism_pareto")["pending"] == counts[2] + counts[3]
    store.close()


def test_principle_orders_are_independent_and_have_no_discotope_dependency():
    store = build_candidate_queues(*_inputs(supported_positions=(4, 5)))
    heads = {principle: store.head("r1", principle, 1)[0] for principle in PRINCIPLES}
    assert heads["epitope_support"].epitope_supported_sites > 0
    assert heads["sector_coverage"].geometry_support["target_mean_coverage_fraction"] is not None
    assert heads["mechanism_pareto"].impact_evidence["regional_impact"]["ranking_version"]
    # No epitope/DiscoTope score is read by the mechanism-ranked queue.
    assert heads["mechanism_pareto"].impact_evidence["ranking"]["reference_size"] > 0
    store.close()


def test_same_candidate_merges_all_principle_labels():
    store = build_candidate_queues(*_inputs())
    candidate = store.head("r1", "epitope_support", 1)[0]
    assert set(store.labels(candidate.candidate_id)) == set(PRINCIPLES)
    assert candidate.impact_evidence["ordinary_principles"] == list(PRINCIPLES)
    store.close()


def test_materialized_candidates_restore_triple_gate_metadata():
    store = build_candidate_queues(*_inputs(count=4))
    triple = next(candidate for candidate in store.materialize_region("r1")
                  if candidate.order == 3)
    assert "regional_higher_order_gain" in triple.risk_components
    assert triple.impact_evidence["regional_impact"]["ranking_version"]
    assert triple.impact_evidence["ranking"]["reference_size"] > 0
    store.close()


def test_materialization_reconstructs_shared_reference_profiles_without_losing_evidence(monkeypatch):
    from copy import deepcopy

    expected = {}
    original_add_region = ordinary_adaptive.CandidateQueue.add_region

    def capture_before_persist(store, sector_id, candidates, *args, **kwargs):
        expected.update({(sector_id, candidate.candidate_id): {
            "impact_evidence": deepcopy(candidate.impact_evidence),
            "risk_components": deepcopy(candidate.risk_components),
            "geometry_support": deepcopy(candidate.geometry_support),
            "coverage_by_scenario": deepcopy(candidate.coverage_by_scenario),
            "covered_positions_by_scenario": deepcopy(candidate.covered_positions_by_scenario),
        } for candidate in candidates})
        result = original_add_region(store, sector_id, candidates, *args, **kwargs)
        for candidate in candidates:
            key = (sector_id, candidate.candidate_id.split("|", 1)[1])
            row = store.db.execute("SELECT metadata FROM candidates WHERE candidate_id=?",
                                   (candidate.candidate_id,)).fetchone()
            persisted = ordinary_adaptive._decode_metadata(row["metadata"])["impact_evidence"]
            expected[key]["persisted_impact_evidence"] = persisted
        return result

    monkeypatch.setattr(ordinary_adaptive.CandidateQueue, "add_region", capture_before_persist)
    store = build_candidate_queues(*_inputs(count=4))
    candidate = store.head("r1", "mechanism_pareto", 1)[0]
    before = expected[("r1", candidate.candidate_id.split("|", 1)[1])]
    row = store.db.execute("SELECT metadata FROM candidates WHERE candidate_id=?",
                           (candidate.candidate_id,)).fetchone()
    assert isinstance(row["metadata"], bytes)
    metadata = ordinary_adaptive._decode_metadata(row["metadata"])
    assert "regional_reference_profiles" not in metadata["impact_evidence"]
    assert "geometry_support" not in metadata
    assert metadata["regional_reference_profile_keys"] == list(
        before["impact_evidence"]["regional_reference_profiles"])
    expected_persisted = {key: value for key, value in before["impact_evidence"].items()
                          if key != "regional_reference_profiles"}
    expected_persisted["ranking"] = dict(expected_persisted["ranking"])
    if metadata["ranking_reference_profile_key"] is not None:
        expected_persisted["ranking"].pop("reference_values")
        assert "reference_values" not in metadata["impact_evidence"]["ranking"]
    expected_persisted["ordinary_principles"] = metadata["impact_evidence"]["ordinary_principles"]
    assert before["persisted_impact_evidence"] == json.loads(json.dumps(expected_persisted))
    restored = store.materialize(candidate.candidate_id)
    restored_core = {key: value for key, value in restored.impact_evidence.items()
                     if key not in {"ordinary_principles", "ordinary_queue"}}
    assert restored_core == json.loads(json.dumps(before["impact_evidence"]))
    assert restored.risk_components == before["risk_components"]
    for field in ("geometry_support", "coverage_by_scenario", "covered_positions_by_scenario"):
        assert getattr(restored, field) == json.loads(json.dumps(before[field]))
    store.close()


def test_geometry_coverage_is_cached_per_site_set_without_shared_mutables(monkeypatch):
    original = ordinary_adaptive._candidate_coverage
    computed = []

    def capture(candidate, *args):
        result = original(candidate, *args)
        computed.append((candidate.order,
                         tuple(sorted(option.residue.uniprot_pos for option in candidate.mutations)),
                         candidate))
        return result

    monkeypatch.setattr(ordinary_adaptive, "_candidate_coverage", capture)
    store = build_candidate_queues(*_inputs(count=4))
    all_candidates = store.materialize_region("r1")
    duplicate_sets = {}
    for candidate in all_candidates:
        key = (candidate.order, tuple(sorted(option.residue.uniprot_pos
                                              for option in candidate.mutations)))
        duplicate_sets.setdefault(key, []).append(candidate)
    assert len(computed) == len(duplicate_sets)
    pair = next(group for group in duplicate_sets.values() if len(group) > 1)
    assert pair[0].geometry_support is not pair[1].geometry_support
    assert pair[0].geometry_support["scenarios"] is not pair[1].geometry_support["scenarios"]
    store.close()


def test_triple_upgrade_is_independent_by_principle():
    store = build_candidate_queues(*_inputs(supported_positions=(4, 5)))
    triple = next(c for c in store.materialize_region("r1") if c.order == 3
                  and [o.residue.uniprot_pos for o in c.mutations] == [3, 4, 5])
    double = next(
        c for c in store.materialize_region("r1")
        if c.order == 2
        and [o.residue.uniprot_pos for o in c.mutations] == [1, 2]
        and {o.mutant_aa for o in c.mutations} == {"D"}
    )
    for option in triple.mutations:
        option.epitope_supported = True
    for option in double.mutations:
        option.epitope_supported = option.residue.uniprot_pos == 1
    triple.geometry_support["target_mean_coverage_fraction"] = .60
    double.geometry_support["target_mean_coverage_fraction"] = .50
    assert strict_triple_upgrade(triple, double, "sector_coverage")
    assert strict_triple_upgrade(triple, double, "epitope_support")
    triple.geometry_support["target_mean_coverage_fraction"] = .54
    assert not strict_triple_upgrade(triple, double, "sector_coverage")
    assert strict_triple_upgrade(triple, double, "epitope_support")
    assert (
        candidate_impact(triple)["mean_absolute_delta"]["volume_delta_a3"]
        > candidate_impact(double)["mean_absolute_delta"]["volume_delta_a3"]
    )
    assert triple_upgrade_for_principle(triple, double, "mechanism_pareto")
    assert strict_triple_upgrade(triple, double)
    store.db.execute("DELETE FROM queue WHERE sector_id='r1'")
    store.db.execute("DELETE FROM candidates")
    rankings = {principle: [double, triple] for principle in PRINCIPLES}
    store.add_region("r1", [double, triple], rankings, {
        "epitope_support": set(), "mechanism_pareto": {triple.candidate_id}, "sector_coverage": set()})
    assert store.head("r1", "mechanism_pareto", 1)[0].candidate_id == triple.candidate_id
    assert store.head("r1", "sector_coverage", 1)[0].candidate_id == double.candidate_id
    store.close()


def test_queue_completion_and_budget_status_are_explicit():
    store = build_candidate_queues(*_inputs(count=3))
    head = store.head("r1", "epitope_support", 1)[0]
    store.complete("r1", "epitope_support", head.candidate_id)
    store._refresh_state("r1", "epitope_support")
    assert store.status("r1", "epitope_support")["status"] == "complete"
    store.mark_budget_limited("r1", "mechanism_pareto", 64)
    assert store.status("r1", "mechanism_pareto") == {
        "status": "budget_limited", "budget_limit": 64}
    store.close()


def test_binder_epitope_lane_is_recorded_not_scheduled():
    store = build_candidate_queues(*_inputs(count=3), epitope_applicable=False)
    assert store.status("r1", "epitope_support") == {"status": "not_applicable"}
    assert store.status("r1", "mechanism_pareto")["status"] == "ready"
    store.close()
