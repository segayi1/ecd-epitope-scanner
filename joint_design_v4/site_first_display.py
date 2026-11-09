"""Display-only qualification for ordinary site-first recommendations.

These checks never change the physical FoldX result or scientific eligibility.
They only keep items that explicitly need review out of formal panels/cards.
"""

import math
import statistics

REVIEW_RISK_LABELS = {
    "SAMPLING_SENSITIVE_REVIEW",
    "REPEAT_EVIDENCE_INCOMPLETE",
    "SAMPLING_DIAGNOSTIC_UNKNOWN",
}


def confirmed_repeat_pass(audit):
    """Recognize a complete, explicitly versioned 3+7 confirmation."""
    values = audit.get("foldx_repeat_slot_values_kcal_mol")
    return (
        audit.get("repeat_confirmation_policy") == "ten_slot_median_v1"
        and audit.get("repeat_confirmation_status") == "CONFIRMED_PASS"
        and audit.get("recommendation_risk_label") == "TEN_SLOT_MEDIAN_CONFIRMED"
        and audit.get("sampling_decision_diagnostic_status") == "COMPLETE"
        and audit.get("sampling_decision_sensitive") is True
        and audit.get("cross_execution_conflict") is not True
        and audit.get("foldx_repeat_gate_kcal_mol") == 2.0
        and audit.get("foldx_repeat_slots") == list(range(10))
        and not any(isinstance(slot, bool) for slot in audit["foldx_repeat_slots"])
        and isinstance(values, list) and len(values) == 10
        and all(isinstance(value, (int, float)) and not isinstance(value, bool)
                and math.isfinite(value) for value in values)
        and sum(value > 2.0 for value in values[:3]) == 1
        and statistics.median(values) <= 2.0
    )


def known_display_review_reasons(candidate):
    reasons = []
    if candidate.qualification == "REVIEW_REQUIRED" or any(
        mutation.review_required for mutation in candidate.mutations
    ):
        reasons.append("STRUCTURAL_CONFIDENCE_REVIEW_REQUIRED")
    impact = candidate.impact_evidence or {}
    if impact.get("policy") == "regional_mechanism_pareto_v1" and impact.get("missing_axes"):
        reasons.append("MECHANISM_EVIDENCE_LIMITED")
    return reasons


def display_review_reasons(candidate, repeat_audit=None):
    audit = repeat_audit
    if audit is None:
        audit = (candidate.impact_evidence or {}).get("stability_repeat_audit", {})
    audit = audit or {}
    reasons = known_display_review_reasons(candidate)
    confirmed = confirmed_repeat_pass(audit)
    risk_label = audit.get("recommendation_risk_label")
    if risk_label in REVIEW_RISK_LABELS and not confirmed:
        reasons.append(risk_label)
    if (audit.get("sampling_decision_diagnostic_status") != "COMPLETE"
            or (audit.get("sampling_decision_sensitive") is not False and not confirmed)):
        reasons.append("REPEAT_EVIDENCE_INCOMPLETE")
    if audit.get("cross_execution_conflict") is True:
        reasons.append("CROSS_EXECUTION_CONFLICT_REVIEW_REQUIRED")
    return sorted(set(reasons))


def apply_display_qualification(candidate, repeat_audit=None):
    reasons = display_review_reasons(candidate, repeat_audit)
    candidate.display_exclusion_reasons = reasons
    candidate.display_eligible = not reasons
    return reasons


def apply_site_first_repeat_review(candidate, result, *, confirmation_policy=None):
    """Attach the actual repeat diagnostic and display verdict to a PASS result."""
    from .stability_repeat_diagnostics import foldx_repeat_gate_diagnostic

    scores = list(result.get("scores") or [])
    expected_count = len(scores) if len(scores) in {3, 10} else 10
    audit = dict(result.get("repeat_audit") or {})
    diagnostic = foldx_repeat_gate_diagnostic(
        scores, 2.0, expected_count=expected_count, slots=result.get("slots"))
    audit.update({
        **diagnostic,
        "cross_execution_conflict": audit.get("cross_execution_conflict"),
        "cross_execution_evidence_status": audit.get(
            "cross_execution_evidence_status", "NOT_ASSESSED"),
    })
    if confirmation_policy == "ten_slot_median_v1":
        audit.update(repeat_confirmation_policy=confirmation_policy,
                     repeat_confirmation_status="NOT_CONFIRMED")
        if len(scores) == 10 and result.get("status") == "PASS":
            audit["repeat_confirmation_status"] = "CONFIRMED_PASS"
            audit["recommendation_risk_label"] = "TEN_SLOT_MEDIAN_CONFIRMED"
            if not confirmed_repeat_pass(audit):
                audit["repeat_confirmation_status"] = "NOT_CONFIRMED"
                audit["recommendation_risk_label"] = diagnostic["recommendation_risk_label"]
    result["repeat_audit"] = audit
    candidate.impact_evidence["stability_repeat_audit"] = audit
    reasons = display_review_reasons(candidate, audit)
    if reasons:
        result["display_review_reasons"] = reasons
        return "DISPLAY_REVIEW"
    return None
