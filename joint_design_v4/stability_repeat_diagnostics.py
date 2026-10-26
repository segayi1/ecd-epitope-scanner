"""Deterministic diagnostics for FoldX repeat spread against the configured gate."""

import math


def foldx_repeat_gate_diagnostic(values, gate_kcal_mol, *, expected_count=5, slots=None):
    parsed = []
    invalid_count = 0
    for value in values or []:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            invalid_count += 1
            continue
        if math.isfinite(numeric):
            parsed.append(numeric)
        else:
            invalid_count += 1
    try:
        gate = float(gate_kcal_mol)
    except (TypeError, ValueError):
        gate = math.nan
    gate_valid = math.isfinite(gate)
    above = sum(value > gate for value in parsed) if gate_valid else None
    at_or_below = sum(value <= gate for value in parsed) if gate_valid else None
    slots_complete = slots is None or list(slots) == list(range(expected_count))
    complete = len(parsed) == expected_count and invalid_count == 0 and gate_valid and slots_complete
    sensitive = bool(above and at_or_below) if complete else None
    return {
        "foldx_repeat_gate_kcal_mol": gate if gate_valid else None,
        "foldx_repeat_min_kcal_mol": min(parsed) if parsed else None,
        "foldx_repeat_max_kcal_mol": max(parsed) if parsed else None,
        "foldx_repeat_above_gate_count": above,
        "foldx_repeat_at_or_below_gate_count": at_or_below,
        "foldx_repeat_invalid_value_count": invalid_count,
        "sampling_decision_diagnostic_status": "COMPLETE" if complete else "INCOMPLETE",
        "sampling_decision_sensitive": sensitive,
        "recommendation_risk_label": ("REPEAT_EVIDENCE_INCOMPLETE" if not complete else
                                       "SAMPLING_SENSITIVE_REVIEW" if sensitive else
                                       "NO_WITHIN_RUN_GATE_CROSSING"),
        "cross_execution_conflict": None,
        "cross_execution_evidence_status": "NOT_ASSESSED",
    }
