"""Build independent, auditable risk components for V4.1 candidates."""

from __future__ import annotations

from typing import Any

from .models import V4Candidate, V4SurfaceHypothesis


def _response_weighted_mean(
    candidate: V4Candidate,
    hypotheses: list[V4SurfaceHypothesis],
    key: str,
) -> float | None:
    numerator = denominator = 0.0
    for hypothesis in hypotheses:
        value = hypothesis.structural_uncertainty.get(key)
        response = candidate.hypothesis_perturbation or candidate.hypothesis_response
        weight = max(0.0, float(response.get(hypothesis.hypothesis_id, 0.0)))
        if value is None or weight <= 0:
            continue
        numerator += float(value) * weight
        denominator += weight
    return numerator / denominator if denominator else None


def annotate_candidate_risks(
    candidates: list[V4Candidate],
    records: list[Any],
    hypotheses: list[V4SurfaceHypothesis],
) -> None:
    """Attach raw risk evidence without collapsing it to one probability.

    Hard mutation exclusions remain governed by the existing evidence gate.
    These values are separate planning evidence for tie-breaking and reporting.
    """
    records_by_position = {item.uniprot_pos: item for item in records}
    for candidate in candidates:
        mutations = candidate.mutations
        count = max(1, len(mutations))
        combination = candidate.stability
        single_positive_sum = sum(
            max(0.0, float(option.spurs_ddg_kcal_mol))
            for option in mutations if option.spurs_ddg_kcal_mol is not None
        )
        shielded = sum(bool(option.glycan_risk) for option in mutations)
        prohibited_positions = {
            position
            for hypothesis in hypotheses
            for position in hypothesis.glycan_evidence.get("mutation_prohibition_positions", [])
        }
        prohibited_mutations = sum(
            option.residue.uniprot_pos in prohibited_positions for option in mutations
        )
        disulfide = sum(
            bool(getattr(records_by_position.get(option.residue.uniprot_pos), "disulfide", False))
            for option in mutations
        )
        interface = 0
        for option in mutations:
            record = records_by_position.get(option.residue.uniprot_pos)
            flags = {str(flag).lower() for flag in getattr(record, "risk_flags", [])}
            interface += bool(
                getattr(record, "interface", False)
                or getattr(record, "domain_interface", False)
                or bool(flags & {"interface", "domain-interface", "domain interface"})
            )
        candidate.risk_components = {
            **candidate.risk_components,
            "single_spurs_positive_sum_kcal_mol": round(single_positive_sum, 6),
            "spurs_multi_ddg_kcal_mol": (
                combination.spurs_multi_ddg_kcal_mol if combination else None
            ),
            "foldx_ddg_kcal_mol": combination.foldx_ddg_kcal_mol if combination else None,
            "glycan_shielded_mutation_fraction": round(shielded / count, 6),
            "glycan_prohibition_mutation_fraction": round(prohibited_mutations / count, 6),
            "hypothesis_near_glycan_fraction": None,
            "disulfide_mutation_fraction": round(disulfide / count, 6),
            "native_interface_mutation_fraction": round(interface / count, 6),
            "hypothesis_low_confidence_surface_fraction": _response_weighted_mean(
                candidate, hypotheses, "low_confidence_surface_fraction"
            ),
            "response_weighted_median_local_pae_a": _response_weighted_mean(
                candidate, hypotheses, "median_local_pae_a"
            ),
        }
        glycan_values = [
            hypothesis.glycan_evidence.get("near_resolved_glycan_area_fraction")
            for hypothesis in hypotheses
            if hypothesis.glycan_evidence.get("near_resolved_glycan_area_fraction") is not None
            and (candidate.hypothesis_perturbation or candidate.hypothesis_response).get(hypothesis.hypothesis_id, 0.0) > 0
        ]
        glycan_weights = [
            (candidate.hypothesis_perturbation or candidate.hypothesis_response)[hypothesis.hypothesis_id]
            for hypothesis in hypotheses
            if hypothesis.glycan_evidence.get("near_resolved_glycan_area_fraction") is not None
            and (candidate.hypothesis_perturbation or candidate.hypothesis_response).get(hypothesis.hypothesis_id, 0.0) > 0
        ]
        if glycan_weights:
            candidate.risk_components["hypothesis_near_glycan_fraction"] = round(
                sum(float(value) * weight for value, weight in zip(glycan_values, glycan_weights))
                / sum(glycan_weights),
                6,
            )
        else:
            candidate.risk_components["hypothesis_near_glycan_fraction"] = None
