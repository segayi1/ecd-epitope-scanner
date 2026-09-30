"""Explain an empty panel without confusing tool failure with scientific rejection."""

from __future__ import annotations

from collections import Counter

from .models import PredictionResult, V4Candidate


def recommendation_diagnostics(
    candidates: list[V4Candidate], selected: list[V4Candidate],
    workers: dict[str, PredictionResult], *, mutation_options: list,
) -> dict:
    skipped_codes = {"NO_CANDIDATES_TO_SCORE", "NO_CONFIRMABLE_CANDIDATES"}
    skipped = {
        name: result.errors for name, result in workers.items()
        if result.errors and all(error.get("code") in skipped_codes for error in result.errors)
    }
    failures = {
        name: result.errors for name, result in workers.items()
        if name not in skipped and (result.status != "success" or result.errors)
    }
    eligible = [item for item in candidates if item.formal_eligible]
    if selected:
        outcome = "PANEL_SELECTED"
    elif "spurs_single" in failures:
        outcome = "SINGLE_STABILITY_EVIDENCE_UNAVAILABLE"
    elif not candidates:
        outcome = "NO_GEOMETRICALLY_FEASIBLE_COMBINATIONS"
    elif "spurs_multi" in failures:
        outcome = "COMBINATION_STABILITY_EVIDENCE_UNAVAILABLE"
    elif "foldx" in failures:
        outcome = "FOLDX_CONFIRMATION_UNAVAILABLE"
    elif not eligible:
        outcome = "ALL_CANDIDATES_FAILED_FORMAL_GATES"
    else:
        outcome = "NO_INFORMATIVE_PANEL_WITHIN_CONSTRAINTS"
    return {
        "outcome": outcome,
        "mutation_option_count": len(mutation_options),
        "eligible_mutation_option_count": sum(item.eligible for item in mutation_options),
        "candidate_count": len(candidates),
        "unique_combination_count": len({
            tuple(sorted(option.notation for option in item.mutations))
            for item in candidates
        }),
        "formal_eligible_candidate_count": len(eligible),
        "selected_construct_count": len(selected),
        "formal_exclusion_counts": dict(sorted(Counter(
            reason for item in candidates for reason in item.formal_exclusion_reasons
        ).items())),
        "worker_failures": failures,
        "workers_not_run": skipped,
    }
