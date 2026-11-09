from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, V4Candidate, V4MutationOption
from ecd_epitope_scanner.joint_design_v4.site_first_display import (
    apply_display_qualification, apply_site_first_repeat_review,
)


def make_candidate():
    mutation = V4MutationOption(
        ResidueKey("0", "A", 1, uniprot_pos=1, one_letter="D"),
        "A", "D1A", "charge", 0.1, 1.0, True,
    )
    return V4Candidate("c1", "VS001", [mutation], order=1, formal_eligible=True)


def test_unknown_cross_execution_state_is_reviewed_without_fabricating_conflict():
    candidate = make_candidate()
    audit = {
        "recommendation_risk_label": "SAMPLING_DIAGNOSTIC_UNKNOWN",
        "sampling_decision_diagnostic_status": "COMPLETE",
        "sampling_decision_sensitive": False,
        "cross_execution_conflict": None,
        "cross_execution_evidence_status": "NOT_ASSESSED",
    }
    reasons = apply_display_qualification(candidate, audit)
    assert not candidate.display_eligible
    assert reasons == ["SAMPLING_DIAGNOSTIC_UNKNOWN"]
    assert audit["cross_execution_conflict"] is None


def test_confirmed_cross_execution_conflict_requires_display_review():
    candidate = make_candidate()
    audit = {
        "recommendation_risk_label": "NO_WITHIN_RUN_GATE_CROSSING",
        "sampling_decision_diagnostic_status": "COMPLETE",
        "sampling_decision_sensitive": False,
        "cross_execution_conflict": True,
        "cross_execution_evidence_status": "CONFIRMED_CONFLICT",
    }
    reasons = apply_display_qualification(candidate, audit)
    assert reasons == ["CROSS_EXECUTION_CONFLICT_REVIEW_REQUIRED"]
    assert not candidate.display_eligible
    assert candidate.formal_eligible


def test_unassessed_cross_execution_state_stays_unknown():
    candidate = make_candidate()
    audit = {
        "recommendation_risk_label": "NO_WITHIN_RUN_GATE_CROSSING",
        "sampling_decision_diagnostic_status": "COMPLETE",
        "sampling_decision_sensitive": False,
        "cross_execution_conflict": None,
        "cross_execution_evidence_status": "NOT_ASSESSED",
    }
    assert apply_display_qualification(candidate, audit) == []
    assert candidate.display_eligible
    assert audit["cross_execution_conflict"] is None
    assert audit["cross_execution_evidence_status"] == "NOT_ASSESSED"


def test_missing_repeat_audit_fails_closed_and_confirmed_conflict_survives_recheck():
    candidate = make_candidate()
    assert apply_display_qualification(candidate) == ["REPEAT_EVIDENCE_INCOMPLETE"]
    assert not candidate.display_eligible

    candidate = make_candidate()
    candidate.impact_evidence["stability_repeat_audit"] = {
        "cross_execution_conflict": True,
        "cross_execution_evidence_status": "CONFIRMED_CONFLICT",
    }
    result = {"scores": [0.5] * 10, "slots": list(range(10)), "repeat_audit": {
        "cross_execution_conflict": True,
        "cross_execution_evidence_status": "CONFIRMED_CONFLICT",
    }}
    assert apply_site_first_repeat_review(candidate, result) == "DISPLAY_REVIEW"
    assert result["repeat_audit"]["sampling_decision_diagnostic_status"] == "COMPLETE"
    assert result["repeat_audit"]["sampling_decision_sensitive"] is False
    assert result["repeat_audit"]["cross_execution_conflict"] is True
    assert "CROSS_EXECUTION_CONFLICT_REVIEW_REQUIRED" in result["display_review_reasons"]


def test_structural_and_mechanism_review_are_independent_display_exclusions():
    candidate = make_candidate()
    candidate.mutations[0].review_required = True
    candidate.impact_evidence = {
        "policy": "regional_mechanism_pareto_v1",
        "missing_axes": ["aromatic_delta"],
        "stability_repeat_audit": {
            "recommendation_risk_label": "NO_WITHIN_RUN_GATE_CROSSING",
            "sampling_decision_diagnostic_status": "COMPLETE",
            "sampling_decision_sensitive": False,
        },
    }
    reasons = apply_display_qualification(candidate)
    assert reasons == ["MECHANISM_EVIDENCE_LIMITED", "STRUCTURAL_CONFIDENCE_REVIEW_REQUIRED"]
    assert not candidate.display_eligible
    assert candidate.formal_eligible  # display review does not rewrite the scientific flag
