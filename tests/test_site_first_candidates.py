import json
from types import SimpleNamespace

from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, V4Sector
from ecd_epitope_scanner.joint_design_v4.site_first import (
    POLICY, _exact_parent_mutations, _triple_precedes_parents, build_site_first_queues,
)


def _pool(positions, n=7):
    rows = []
    aas = "DEFHIKL"
    for pos in positions:
        for i in range(n):
            aa = aas[i]
            rows.append(SimpleNamespace(
                residue=ResidueKey("1", "A", pos, uniprot_pos=pos, one_letter="A"),
                mutant_aa=aa, eligible=True, review_required=False, hard_exclusions=[],
                soft_warnings=[], spurs_ddg_kcal_mol=i / 10, chemistry_class="test",
                stability_band="PASS", epitope_prior=i / 10, epitope_supported=i % 2 == 0,
                epitope_evidence_status="SCORED",
                impact_evidence={"sidechain_exposure": {"value": .5 + i / 100}},
                chemical_descriptors={"charge_delta": i % 2, "volume_delta_a3": i,
                    "hydropathy_delta": i / 2, "donor_delta": i % 3,
                    "acceptor_delta": i % 4, "aromatic_delta": i % 2},
            ))
    return rows


def _inputs(n=7):
    graph = {1: {2: 8.0, 3: 16.0}, 2: {1: 8.0, 3: 8.0}, 3: {1: 16.0, 2: 8.0}}
    sector = V4Sector("s", [1, 2, 3], 1, 30, 30, .5)
    return _pool([1, 2, 3], n), [sector], graph, {1: 10., 2: 10., 3: 10.}, "AAAAAA"


def test_site_count_and_geometry_are_independent_of_replacement_pool_size():
    first = build_site_first_queues(*_inputs(3))
    larger = build_site_first_queues(*_inputs(7))
    a, b = first.ledger(), larger.ledger()
    assert a["site_set_count"] == b["site_set_count"]
    assert a["geometry_evaluations"] == b["geometry_evaluations"]
    assert b["geometry_evaluations"] == b["site_set_count"]
    assert b["expanded_candidate_count"] == 0
    assert b["initial_candidate_count"] == 0
    assert b["timings"]["site_preparation_seconds"] >= 0
    json.dumps(b)
    assert b["policy"] == POLICY
    first.close(); larger.close()


def test_one_lane_head_never_crosses_site_set_and_expands_only_after_failures():
    store = build_site_first_queues(*_inputs())
    sector, principle = "s", "epitope_support"
    head = store.head(sector, principle, limit=99)
    assert len(head) <= 9
    group_id = head[0].impact_evidence["site_set_id"]
    assert all(item.impact_evidence["site_set_id"] == group_id for item in head)
    before = store.ledger()["expanded_candidate_count"]
    for candidate in list(head):
        store.record_outcome(sector, principle, candidate.candidate_id, "SCIENTIFIC_FAILURE")
    expanded = store.ledger()["expanded_candidate_count"]
    assert 0 < expanded <= 27
    assert store.ledger()["initial_candidate_count"] == 9
    assert expanded == 27  # 6 x 6 minus the already evaluated 3 x 3.
    expanded_head = store.head(sector, principle, 999)
    assert all(c.impact_evidence["site_set_id"] == group_id for c in expanded_head)
    for candidate in expanded_head:
        store.record_outcome(sector, principle, candidate.candidate_id, "SCIENTIFIC_FAILURE")
    deep_head = store.head(sector, principle, 999)
    assert len(deep_head) == 13  # A 7 x 7 pass over the prior 6 x 6 failures.
    assert all(c.impact_evidence["representative_tier"] == "deep" for c in deep_head)
    assert store.ledger()["deep_candidate_count"] == 13
    store.close()


def test_prefetch_reads_hidden_groups_without_advancing_or_expanding_queue():
    store = build_site_first_queues(*_inputs())
    head = store.head("s", "sector_coverage", 1)[0]
    before = store.ledger()
    groups = list(store.prefetch_groups("s", "sector_coverage"))
    assert len(groups) >= 2
    assert groups[0][0].candidate_id == head.candidate_id
    assert len({group[0].impact_evidence["site_set_id"] for group in groups}) == len(groups)
    assert all(c.impact_evidence["representative_tier"] == "initial" for group in groups for c in group)
    after = store.ledger()
    assert after["expanded_candidate_count"] == before["expanded_candidate_count"]
    assert after["initial_candidate_count"] == before["initial_candidate_count"]
    assert store.head("s", "sector_coverage", 1)[0].candidate_id == head.candidate_id
    store.close()


def test_joint_sequence_gate_and_deterministic_queue_order():
    data = _inputs()
    one = build_site_first_queues(*data)
    two = build_site_first_queues(*data)
    def ids(store):
        return [c.candidate_id for c in store.head("s", "mechanism_pareto", 4)]
    assert ids(one) == ids(two)
    for lane in one.queue_states():
        assert lane in two.queue_states()
    one.close(); two.close()


def test_triples_without_five_percent_gain_are_excluded():
    # Complete small triangle: every pair and triple covers the same area.
    graph = {1: {2: 1., 3: 1.}, 2: {1: 1., 3: 1.}, 3: {1: 1., 2: 1.}}
    sector = V4Sector("s", [1, 2, 3], 1, 30, 30, .5)
    store = build_site_first_queues(_pool([1, 2, 3]), [sector], graph,
                                   {1: 10., 2: 10., 3: 10.}, "AAAAAA")
    # All offered site sets remain pairs; no weak-parent fallback can admit a triple.
    assert all(len(row["site_set_id"].split("|", 1)[1].split(",")) == 2
               for row in store.ledger()["queues"] if row.get("site_set_id"))
    assert not any(c.order == 3 for c in store.head("s", "sector_coverage", 999))
    triple = next(group for group in store._lanes[("s", "sector_coverage")]["groups"] if group["order"] == 3)
    assert triple["triple_gain_gate_passed"] is False
    store.close()


def test_triple_upgrade_compares_with_strongest_exact_parent():
    from types import SimpleNamespace
    def mutations(order):
        return [SimpleNamespace(epitope_evidence_status="SCORED", epitope_prior=.5) for _ in range(order)]
    weak = SimpleNamespace(order=2, epitope_supported_sites=1, epitope_score=.2, candidate_id="weak", mutations=mutations(2))
    strong = SimpleNamespace(order=2, epitope_supported_sites=2, epitope_score=2., candidate_id="strong", mutations=mutations(2))
    triple = SimpleNamespace(order=3, epitope_supported_sites=3, epitope_score=3., candidate_id="triple", mutations=mutations(3))
    assert not _triple_precedes_parents(triple, [weak, strong], "epitope_support", .06)
    strong.epitope_score = 0.
    assert _triple_precedes_parents(triple, [weak, strong], "epitope_support", .06)


def test_triple_parents_preserve_its_exact_mutant_amino_acids():
    from types import SimpleNamespace
    triple = SimpleNamespace(mutations=[
        SimpleNamespace(residue=SimpleNamespace(uniprot_pos=1), mutant_aa="D"),
        SimpleNamespace(residue=SimpleNamespace(uniprot_pos=2), mutant_aa="K"),
        SimpleNamespace(residue=SimpleNamespace(uniprot_pos=3), mutant_aa="Y"),
    ])
    parents = _exact_parent_mutations(triple)
    signatures = {tuple(sorted((item.residue.uniprot_pos, item.mutant_aa) for item in parent))
                  for parent in parents}
    assert signatures == {
        ((1, "D"), (2, "K")), ((1, "D"), (3, "Y")), ((2, "K"), (3, "Y")),
    }


def test_joint_sequence_gate_rejects_sequon_created_by_two_mutations():
    rows = _pool([2, 3], n=1)
    rows[0].residue = ResidueKey("1", "A", 2, uniprot_pos=2, one_letter="A")
    rows[0].mutant_aa = "D"
    rows[1].residue = ResidueKey("1", "A", 3, uniprot_pos=3, one_letter="A")
    rows[1].mutant_aa = "S"
    sector = V4Sector("seq", [2, 3], 1, 20, 20, .5)
    graph = {2: {3: 1.}, 3: {2: 1.}}
    store = build_site_first_queues(rows, [sector], graph, {2: 10., 3: 10.}, "NAA")
    assert store.head("seq", "epitope_support", 99) == []
    assert store.ledger()["expanded_candidate_count"] <= 1
    store.close()


def test_epitope_priority_changes_first_site_when_geometry_is_tied():
    rows = _pool([1, 2, 3, 4])
    for row in rows:
        supported = row.residue.uniprot_pos in {1, 2}
        row.epitope_supported = supported
        row.epitope_prior = 1.0 if supported else 0.0
    graph = {p: {q: 1.0 for q in (1, 2, 3, 4) if q != p} for p in (1, 2, 3, 4)}
    sector = V4Sector("tie", [1, 2, 3, 4], 1, 40, 40, .5)
    store = build_site_first_queues(rows, [sector], graph, {p: 10. for p in graph}, "AAAAAA")
    first = store.head("tie", "epitope_support", 1)[0]
    assert first.impact_evidence["site_set_id"] == "tie|1,2"
    store.close()


def test_site_with_only_mixed_sequence_legal_combinations_is_retained(monkeypatch):
    from ecd_epitope_scanner.joint_design_v4 import coverage_profile
    monkeypatch.setattr(coverage_profile, "sequence_reasons", lambda sequence, mutations:
                        ["fixture_joint_gate"] if len(mutations) > 1
                        and len({mutation[2] for mutation in mutations}) == 1 else [])
    rows = _pool([1, 2])
    sector = V4Sector("mixed", [1, 2], 1, 20, 20, .5)
    graph = {1: {2: 1.}, 2: {1: 1.}}
    store = build_site_first_queues(rows, [sector], graph, None, "AAAA")
    group = store._lanes[("mixed", "epitope_support")]["groups"][0]
    assert group["witness_count"] == 0
    candidates = store.head("mixed", "epitope_support", 99)
    assert candidates and all(len({option.mutant_aa for option in c.mutations}) == 2 for c in candidates)
    store.close()


def test_epitope_site_order_uses_native_evidence_without_role_witnesses(monkeypatch):
    from ecd_epitope_scanner.joint_design_v4 import coverage_profile
    monkeypatch.setattr(coverage_profile, "sequence_reasons", lambda sequence, mutations:
                        ["fixture_joint_gate"] if len(mutations) > 1
                        and len({mutation[2] for mutation in mutations}) == 1 else [])
    rows = _pool([1, 2, 3, 4])
    for row in rows:
        row.epitope_supported = row.residue.uniprot_pos in {1, 2}
        row.epitope_prior = 0.9 if row.residue.uniprot_pos in {1, 2} else 0.1
    graph = {p: {q: 1.0 for q in (1, 2, 3, 4) if p != q} for p in (1, 2, 3, 4)}
    sector = V4Sector("roleless", [1, 2, 3, 4], 1, 40, 40, .5)
    store = build_site_first_queues(rows, [sector], graph, {p: 10. for p in graph}, "AAAAAA")
    first = store._lanes[("roleless", "epitope_support")]["groups"][0]
    assert first["site_set_id"] == "roleless|1,2"
    assert first["witness_count"] == 0
    store.close()


def test_missing_epitope_prior_stays_unknown_and_ranks_after_known_zero():
    rows = _pool([1, 2, 3, 4])
    for row in rows:
        if row.residue.uniprot_pos in {1, 2}:
            row.epitope_prior = 0.0
            row.epitope_supported = False
        else:
            row.epitope_prior = None
            row.epitope_supported = False
            row.epitope_evidence_status = "MISSING"
    graph = {p: {q: 1.0 for q in (1, 2, 3, 4) if p != q} for p in (1, 2, 3, 4)}
    sector = V4Sector("missing", [1, 2, 3, 4], 1, 40, 40, .5)
    store = build_site_first_queues(rows, [sector], graph, {p: 10. for p in graph}, "AAAAAA")
    groups = store._lanes[("missing", "epitope_support")]["groups"]
    assert groups[0]["site_set_id"] == "missing|1,2"
    unknown = next(group for group in groups if group["positions"] == (3, 4))
    assert unknown["site_priority"] is None
    assert any(gap.get("site_set_id") == "missing|3,4" and gap["status"] == "epitope_prior_missing"
               for gap in store.ledger()["uncompleted_gaps"])
    store.close()


def test_missing_epitope_evidence_cannot_win_principle_queue():
    rows = _pool([1, 2])
    for row in rows:
        row.epitope_prior = None
        row.epitope_evidence_status = "MISSING"
    graph = {1: {2: 1.}, 2: {1: 1.}}
    sector = V4Sector("missing", [1, 2], 1, 20, 20, .5)
    store = build_site_first_queues(rows, [sector], graph, {1: 10., 2: 10.}, "AAAA")
    assert store.head("missing", "epitope_support", 4) == []
    assert store.status("missing", "epitope_support")["status"] == "missing_evidence"
    assert store.head("missing", "sector_coverage", 1)
    store.close()


def test_materialization_keeps_source_record_numbering_and_geometry_evidence():
    rows, sectors, graph, areas, sequence = _inputs()
    records = {p: SimpleNamespace(uniprot_pos=p, aa="A", auth_seq_id=p + 100,
                insertion_code="", label_seq_id=p, sidechain_exposure={"value": .7})
               for p in (1, 2, 3)}
    store = build_site_first_queues(rows, sectors, graph, areas, sequence, records=records)
    candidate = store.head("s", "sector_coverage", 1)[0]
    assert candidate.impact_evidence["protocol"] == POLICY
    assert candidate.coverage_by_scenario
    assert candidate.impact_evidence["regional_geometry_support"] == candidate.geometry_support
    assert all(option.residue.auth_seq_id == option.residue.uniprot_pos + 100
               for option in candidate.mutations)
    store.close()


def test_terminal_queue_gaps_include_missing_candidates_and_budget_once():
    from ecd_epitope_scanner.joint_design_v4.site_first import SiteFirstQueue
    store = SiteFirstQueue()
    store._lanes = {
        ('VS006', 'mechanism_pareto'): {'groups': [], 'group_index': 0, 'state': 'candidate_failures'},
        ('VS004', 'sector_coverage'): {'groups': [{'items': None, 'history': [], 'site_set_id': 'x'}],
                                      'group_index': 0, 'state': 'budget_limited', 'budget_limit': 300},
        ('VS005', 'epitope_support'): {'groups': [], 'group_index': 0, 'state': 'not_applicable'},
    }
    first = store.ledger()['uncompleted_gaps']
    assert first == [
        {'sector_id': 'VS004', 'principle': 'sector_coverage', 'status': 'budget_limited', 'budget_limit': 300},
        {'sector_id': 'VS006', 'principle': 'mechanism_pareto', 'status': 'no_legal_candidates'},
    ]
    assert store.ledger()['uncompleted_gaps'] == first
    store.close()
