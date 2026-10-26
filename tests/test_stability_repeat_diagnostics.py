import pytest

from ecd_epitope_scanner.joint_design_v4.stability_repeat_diagnostics import foldx_repeat_gate_diagnostic


def test_legacy_ctla4_five_slots_cross_gate_without_changing_median_contract():
    values = [2.99488, 0.276247, 3.4539, 1.21002, 0.710114]
    diagnostic = foldx_repeat_gate_diagnostic(values, 2.0)

    assert diagnostic == {
        "foldx_repeat_gate_kcal_mol": 2.0,
        "foldx_repeat_min_kcal_mol": 0.276247,
        "foldx_repeat_max_kcal_mol": 3.4539,
        "foldx_repeat_above_gate_count": 2,
        "foldx_repeat_at_or_below_gate_count": 3,
        "foldx_repeat_invalid_value_count": 0,
        "sampling_decision_diagnostic_status": "COMPLETE",
        "sampling_decision_sensitive": True,
        "recommendation_risk_label": "SAMPLING_SENSITIVE_REVIEW",
        "cross_execution_conflict": None,
        "cross_execution_evidence_status": "NOT_ASSESSED",
    }
    # The existing qualification rule consumes the median; this diagnosis does not impose a five-slot gate.
    assert sorted(values)[2] <= diagnostic["foldx_repeat_gate_kcal_mol"]


def test_engine_audit_reports_crossing_but_keeps_median_consensus_pass():
    from ecd_epitope_scanner.joint_design_v4.engine import _annotate_stability
    from ecd_epitope_scanner.joint_design_v4.models import (
        PredictionResult, ResidueKey, V4Candidate, V4MutationOption,
    )

    mutations = [
        V4MutationOption(ResidueKey("0", "A", pos, uniprot_pos=pos, one_letter="A"), aa,
                         f"A{pos}{aa}", "test", .1, 1.0, True)
        for pos, aa in ((1, "D"), (2, "V"))
    ]
    candidate = V4Candidate("ctla4", "region", mutations, order=2, structural_risk=0)
    values = [2.99488, 0.276247, 3.4539, 1.21002, 0.710114]
    foldx = PredictionResult("success", [{
        "mutations": ["A1D", "A2V"], "ddg": 1.21002,
        "paired_models": [{"sampling_slot": i, "ddg_kcal_mol": value} for i, value in enumerate(values)],
    }])
    multi = PredictionResult("success", [{"mutations": ["A1D", "A2V"], "ddg": .2}])

    _annotate_stability([candidate], [], multi, foldx, {}, {"combination_gate_kcal_mol": 2.0})

    audit = candidate.impact_evidence["stability_repeat_audit"]
    assert candidate.stability.consensus_status == "CONSENSUS_PASS"
    assert audit["foldx_repeat_above_gate_count"] == 2
    assert audit["sampling_decision_sensitive"] is True
    assert audit["foldx_repeat_slot_values_kcal_mol"] == values


@pytest.mark.parametrize("values", [[0.1, 0.2, 1.0, 1.5, 2.0], [-2.0, 0.0, 0.5, 1.0, 1.9]])
def test_all_five_slots_at_or_below_gate_are_not_sampling_sensitive(values):
    diagnostic = foldx_repeat_gate_diagnostic(values, 2.0)
    assert diagnostic["sampling_decision_diagnostic_status"] == "COMPLETE"
    assert diagnostic["sampling_decision_sensitive"] is False
    assert diagnostic["foldx_repeat_above_gate_count"] == 0


def test_all_five_slots_above_gate_are_not_crossing_the_gate():
    diagnostic = foldx_repeat_gate_diagnostic([2.01, 2.2, 2.5, 3.0, 4.0], 2.0)
    assert diagnostic["sampling_decision_sensitive"] is False
    assert diagnostic["foldx_repeat_above_gate_count"] == 5
    assert diagnostic["foldx_repeat_at_or_below_gate_count"] == 0


@pytest.mark.parametrize("values", [[1.0, 1.2, 1.3, 1.4], [1.0, 1.2, 1.3, 1.4, float("nan")],
                                     [1.0, 1.2, 1.3, 1.4, None]])
def test_missing_or_invalid_slot_never_claims_stable_sampling(values):
    diagnostic = foldx_repeat_gate_diagnostic(values, 2.0)
    assert diagnostic["sampling_decision_diagnostic_status"] == "INCOMPLETE"
    assert diagnostic["sampling_decision_sensitive"] is None


def test_configured_gate_is_used_for_crossing_decision():
    diagnostic = foldx_repeat_gate_diagnostic([1.0, 1.2, 1.3, 1.4, 1.8], 1.5)
    assert diagnostic["foldx_repeat_gate_kcal_mol"] == 1.5
    assert diagnostic["foldx_repeat_above_gate_count"] == 1
    assert diagnostic["sampling_decision_sensitive"] is True


def test_wrong_slot_identity_does_not_claim_complete_repeat_diagnostic():
    diagnostic = foldx_repeat_gate_diagnostic([1.0] * 5, 2.0, slots=[0, 1, 2, 3, 8])
    assert diagnostic["sampling_decision_diagnostic_status"] == "INCOMPLETE"
    assert diagnostic["sampling_decision_sensitive"] is None
