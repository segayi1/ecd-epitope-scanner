from ecd_epitope_scanner.joint_design_v4.engine import _foldx_shortlist
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, V4Candidate, V4MutationOption, V4Sector


def candidate(identifier, positions, mutant="A"):
    options = [V4MutationOption(ResidueKey("0", "A", p, uniprot_pos=p, one_letter="D"),
        mutant, f"D{p}{mutant}", "charge", 0.0, 1.0, True) for p in positions]
    return V4Candidate(identifier, "OLD_WINDOW", options, order=len(options))


def fixture():
    graph = {p: {q: 6.0 for q in range(1, 9) if abs(p-q) == 1} for p in range(1, 9)}
    sectors = [V4Sector("A", [1, 2, 3, 4], 0, 40, 40, 0), V4Sector("B", [5, 6, 7, 8], 0, 40, 40, 0)]
    pool = [candidate("a", [1, 2]), candidate("a_variant", [1, 2], "W"), candidate("a2", [3, 4]),
            candidate("b", [5, 6]), candidate("b2", [7, 8]), candidate("cross", [4, 5])]
    return graph, sectors, pool


def test_feedback_search_is_bounded_and_skips_acceptable_regions(monkeypatch):
    from types import SimpleNamespace
    from ecd_epitope_scanner.joint_design_v4 import candidate_search
    graph, sectors, pool = fixture()
    known = pool[0]
    known.stability = SimpleNamespace(consensus_status="COMBINATION_REJECT")
    accepted = pool[3]
    accepted.stability = SimpleNamespace(consensus_status="SPURS_ONLY_REVIEW_REQUIRED")
    available = [candidate(f"new{i}", [1, 3], chr(65+i)) for i in range(15)]
    called = []
    def generate(regions, options, graph, parameters, **kwargs):
        called.append(regions[0].sector_id)
        assert parameters["orders"] == [2]
        return [known] + available
    monkeypatch.setattr(candidate_search, "generate_candidates", generate)
    output, audit = candidate_search.supplemental_exclusive_doubles(
        [known, accepted], sectors, [], graph, {}, dict.fromkeys(graph, 10))
    assert called == ["A"]
    assert len(output) == 12
    assert known not in output
    assert audit[0]["rounds"] == 1
    assert audit[0]["maximum_new_candidates"] == 12


def test_two_exclusive_groups_per_region_before_variants_and_cross():
    graph, sectors, pool = fixture()
    audit = {}
    selected = _foldx_shortlist(pool, 4, sectors=sectors, surface_areas=dict.fromkeys(graph, 10), graph=graph, diagnostics=audit)
    assert {c.candidate_id for c in selected} == {"a", "a2", "b", "b2"}
    assert all(row["scheduled_exclusive_groups"] == 2 for row in audit["foldx_region_resources"])
    assert [c.candidate_id for c in selected] == [c.candidate_id for c in _foldx_shortlist(
        list(reversed(pool)), 4, sectors=list(reversed(sectors)), surface_areas=dict.fromkeys(graph, 10), graph=graph)]


def test_first_exclusive_confirmation_preserves_a_double_baseline():
    graph, sectors, pool = fixture()
    pool.extend([candidate("a_triple", [1, 2, 4]), candidate("b_triple", [5, 6, 8])])
    selected = _foldx_shortlist(pool, 2, sectors=sectors,
        surface_areas=dict.fromkeys(graph, 10), graph=graph)
    assert len(selected) == 2
    assert all(c.order == 2 for c in selected)
    assert [c.candidate_id for c in selected] == [c.candidate_id for c in _foldx_shortlist(
        list(reversed(pool)), 2, sectors=list(reversed(sectors)),
        surface_areas=dict.fromkeys(graph, 10), graph=graph)]


def test_first_round_under_budget_and_cross_fallback():
    graph, sectors, pool = fixture()
    selected = _foldx_shortlist(pool, 2, sectors=sectors, surface_areas=dict.fromkeys(graph, 10), graph=graph)
    assert len(selected) == 2 and any(c.candidate_id.startswith("a") for c in selected) and any(c.candidate_id.startswith("b") for c in selected)
    audit = {}
    selected = _foldx_shortlist([pool[0], pool[-1]], 42, sectors=sectors,
        surface_areas=dict.fromkeys(graph, 10), graph=graph, diagnostics=audit)
    assert len(selected) == 2
    b = next(row for row in audit["foldx_region_resources"] if row["sector_id"] == "B")
    assert b["scheduled_exclusive_groups"] == 0 and b["scheduled_associated_groups"] == 1


def test_exploratory_and_review_do_not_consume_confirmation():
    graph, sectors, pool = fixture()
    pool[0].qualification = "REVIEW_REQUIRED"
    pool[2].search_role = "EXPLORATORY"
    selected = _foldx_shortlist(pool, 42, sectors=sectors, surface_areas=dict.fromkeys(graph, 10), graph=graph)
    assert not {"a", "a2"}.intersection(c.candidate_id for c in selected)


def test_progress_complete_cannot_publish_empty_result(tmp_path, monkeypatch):
    import ecd_epitope_scanner.jobs as jobs
    store = jobs.JobStore(tmp_path)
    job = store.create(jobs.PipelineRequest(input_type="uniprot", identifier="P16410"))
    def run(request, progress):
        progress("COMPLETE", 1.0, "Finished files")
        visible = store.get(job["id"])
        assert visible["status"] == "RUNNING" and visible["result"] is None
        return {"design_summary": {"recommendation_status": "FORMAL_RECOMMENDATION"}}
    monkeypatch.setattr(jobs, "run_pipeline", run)
    runner = jobs.JobRunner(store)
    try:
        runner._run(job["id"])
        assert store.get(job["id"])["status"] == "COMPLETE"
        assert store.get(job["id"])["result"] is not None
    finally:
        runner.executor.shutdown()
