from types import SimpleNamespace as NS

import pytest

from ecd_epitope_scanner.joint_design_v4.models import PredictionResult
from ecd_epitope_scanner.joint_design_v4.ordinary_adaptive_validation import AdaptiveValidation


def candidate(identifier, *mutations):
    return NS(candidate_id=identifier, mutations=[NS(notation=value) for value in mutations])


def spurs(rows, status="success", errors=None):
    return PredictionResult(status, records=rows, errors=errors or [])


def foldx(mutations, values, slot_start=0, *, status="success", attempted=None):
    return PredictionResult(status, records=[{
        "mutations": list(mutations),
        "paired_models": [
            {"sampling_slot": slot_start + index, "ddg_kcal_mol": value}
            for index, value in enumerate(values)
        ],
    }], provenance={"attempted_slots": attempted if attempted is not None else len(values)})


@pytest.mark.parametrize(("score", "expected"), [
    (1.0, "PASS"), (1.00001, "REVIEW"), (2.0, "REVIEW"), (2.00001, "OVER_LIMIT"),
])
def test_spurs_exact_boundary_bands(score, expected):
    item = candidate("c", "A1D")
    result = AdaptiveValidation().run_spurs({"q": [item]}, lambda batch: spurs([
        {"mutations": ["A1D"], "ddg": score}
    ]))
    assert result["results"][("A1D",)]["status"] == expected
    assert ("A1D",) in result["passed"] if expected == "PASS" else ("A1D",) not in result["passed"]


@pytest.mark.parametrize("row", [
    {"mutations": ["A1D"], "ddg": float("nan")},
    {"mutations": ["A1D"], "ddg": "not-a-number"},
    {"mutations": ["A2D"], "ddg": 0.2},
])
def test_spurs_invalid_score_or_identity_is_technical_and_cannot_qualify(row):
    item = candidate("c", "A1D")
    result = AdaptiveValidation().run_spurs({"q": [item]}, lambda batch: spurs([row]))
    assert not result["passed"]
    assert result["ledger"]["technical_failure"] is True
    assert result["results"][("A1D",)]["status"] in {"INVALID_SCORE", "TECHNICAL_FAILURE"}


def test_spurs_partial_keeps_valid_score_but_missing_or_failed_request_never_qualifies():
    valid, missing = candidate("valid", "A1D"), candidate("missing", "A2D")
    result = AdaptiveValidation().run_spurs({"valid-lane": [valid], "missing-lane": [missing]},
        lambda batch: spurs([
            {"mutations": ["A1D"], "ddg": 1.0}
        ], status="partial"))
    assert result["results"][("A1D",)]["status"] == "PASS"
    assert result["results"][("A2D",)]["status"] == "MISSING_RESULT"
    assert set(result["passed"]) == {("A1D",)}

    failed = AdaptiveValidation().run_spurs({"q": [valid]}, lambda batch: spurs([
        {"mutations": ["A1D"], "ddg": 0.2}
    ], status="failed", errors=[{"code": "worker"}]))
    assert failed["results"][("A1D",)]["score"] == 0.2
    assert failed["results"][("A1D",)]["status"] == "BLOCKED_UNSCORED"
    assert not failed["passed"]


def test_foldx_exact_two_passes_without_append():
    item = candidate("c", "A1D", "A2V")
    calls = []

    def predict(c, *, slot_start, runs):
        calls.append((slot_start, runs))
        return foldx([m.notation for m in c.mutations], [2.0] * runs, slot_start)

    validation = AdaptiveValidation()
    validation.spurs_results[("A1D", "A2V")] = {"status": "PASS"}
    result = validation.run_foldx({"q": [item]}, {("A1D", "A2V"): item}, predict)
    assert calls == [(0, 3)]
    assert result["results"][("A1D", "A2V")]["status"] == "PASS"
    assert result["results"][("A1D", "A2V")]["slots"] == [0, 1, 2]


@pytest.mark.parametrize("response", [
    lambda identity: foldx(["A99G"], [1.0, 1.0, 1.0]),
    lambda identity: PredictionResult("success", records=[{
        "mutations": list(identity), "paired_models": [
            {"sampling_slot": 0, "ddg_kcal_mol": 1.0},
            {"sampling_slot": 0, "ddg_kcal_mol": 1.1},
            {"sampling_slot": 2, "ddg_kcal_mol": 1.0},
        ],
    }], provenance={"attempted_slots": 3}),
    lambda identity: PredictionResult("success", records=[{
        "mutations": list(identity), "paired_models": [
            {"sampling_slot": 0, "ddg_kcal_mol": float("nan")},
            {"sampling_slot": 1, "ddg_kcal_mol": 1.0},
            {"sampling_slot": 2, "ddg_kcal_mol": 1.0},
        ],
    }], provenance={"attempted_slots": 3}),
])
def test_foldx_duplicate_identity_or_nonfinite_slot_rejects_without_append(response):
    item = candidate("c", "A1D", "A2V")
    calls = []

    def predict(c, *, slot_start, runs):
        calls.append((slot_start, runs))
        if slot_start:
            pytest.fail("invalid first-stage data must not trigger an append")
        identity = [m.notation for m in c.mutations]
        return response(identity)

    validation = AdaptiveValidation()
    validation.spurs_results[("A1D", "A2V")] = {"status": "PASS"}
    result = validation.run_foldx({"q": [item]}, {("A1D", "A2V"): item}, predict)
    assert calls == [(0, 3)]
    assert result["results"][("A1D", "A2V")]["status"] == "TECHNICAL_FAILURE"


def test_foldx_reserves_and_releases_slots_across_candidates_with_300_cap():
    candidates = [candidate(f"c{i}", f"A{i + 1}D", f"A{i + 101}V") for i in range(31)]
    passed = {tuple(sorted(m.notation for m in item.mutations)): item for item in candidates}
    validation = AdaptiveValidation()
    validation.spurs_results.update({signature: {"status": "PASS"} for signature in passed})
    calls = []

    def predict(item, *, slot_start, runs):
        calls.append((item.candidate_id, slot_start, runs))
        values = [2.1, 1.0, 1.0] if slot_start == 0 else [1.0] * runs
        return foldx([m.notation for m in item.mutations], values, slot_start)

    queues = {f"q{i}": [item] for i, item in enumerate(candidates)}
    result = validation.run_foldx(queues, passed, predict, slot_limit=1000)
    assert result["ledger"]["slot_limit"] == 300
    assert len(calls) == 60
    assert calls[:2] == [("c0", 0, 3), ("c0", 3, 7)]
    assert calls[-2:] == [("c29", 0, 3), ("c29", 3, 7)]
    assert result["ledger"]["attempted_slots"] == 300
    assert result["ledger"]["valid_slots"] == 300
    assert result["ledger"]["reserved_slots"] == 0
    assert result["ledger"]["released_slots"] == 0
    final_signature = tuple(sorted(("A31D", "A131V")))
    assert result["results"][final_signature]["status"] == "BUDGET_INSUFFICIENT"
