from types import SimpleNamespace

from ecd_epitope_scanner.pipeline import _epitope_prior_used_for_panel_selection


def _result(*candidates, priority="sector_first_epitope_weighted_perturbation_then_designable_area"):
    return SimpleNamespace(
        candidates=list(candidates),
        panel=SimpleNamespace(
            budget_decision=SimpleNamespace(selection_priority=priority),
        ),
    )


def _candidate(*, scored=True, display_eligible=True, order=1):
    option = SimpleNamespace(
        epitope_evidence_status="SCORED" if scored else "MISSING",
        epitope_prior=0.8 if scored else None,
        review_required=False,
    )
    return SimpleNamespace(
        mutations=[option], display_eligible=display_eligible,
        qualification="FORMAL_POOL", connected=True, order=order,
        search_role="PRIMARY_POOL",
    )


def test_epitope_flag_uses_eligible_sorted_candidates_even_when_not_selected():
    result = _result(_candidate(scored=True), _candidate(scored=False))

    assert _epitope_prior_used_for_panel_selection(result) is True


def test_epitope_flag_requires_scored_prior_in_sorted_candidate_pool():
    result = _result(_candidate(scored=True, display_eligible=False))

    assert _epitope_prior_used_for_panel_selection(result) is False


def test_epitope_flag_requires_epitope_selector_provenance():
    result = _result(
        _candidate(scored=True),
        priority="direct_coverage_then_geometric_complementarity_then_information_then_risk",
    )

    assert _epitope_prior_used_for_panel_selection(result) is False


def test_epitope_flag_recognizes_target_surface_prior():
    candidate = _candidate(scored=False)
    candidate.impact_evidence = {"surface_perturbation_field": {
        "status": "SCORED", "epitope_status": "PARTIAL"}}
    assert _epitope_prior_used_for_panel_selection(_result(candidate)) is True
