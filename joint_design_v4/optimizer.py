from __future__ import annotations

from collections import defaultdict

from .models import V4Candidate, V4Panel, V4Sector, V4SurfaceHypothesis


def _candidate_rank(candidate: V4Candidate) -> tuple:
    stability = candidate.stability
    spurs = (
        stability.spurs_multi_ddg_kcal_mol
        if stability and stability.spurs_multi_ddg_kcal_mol is not None
        else float("inf")
    )
    foldx = (
        stability.foldx_ddg_kcal_mol
        if stability and stability.foldx_ddg_kcal_mol is not None
        else float("inf")
    )
    return (
        candidate.order > 2,
        max(spurs, foldx),
        candidate.structural_risk,
        -min(candidate.coverage_by_scenario.values(), default=0.0),
        -sum(candidate.coverage_by_scenario.values()),
        candidate.order,
        -candidate.epitope_score,
        candidate.candidate_id,
    )


def _position_signature(candidate: V4Candidate) -> tuple[int, ...]:
    return tuple(sorted(option.residue.uniprot_pos or -1 for option in candidate.mutations))


def _construct_signature(candidate: V4Candidate) -> tuple[tuple[int, str], ...]:
    return tuple(sorted(
        (option.residue.uniprot_pos or -1, option.mutant_aa)
        for option in candidate.mutations
    ))


def _risk_penalty(
    candidate: V4Candidate,
    weights: dict[str, float] | None = None,
) -> float:
    """Convert raw risk evidence to a small, explicit tie-break penalty."""
    components = candidate.risk_components
    positive_stability = [
        max(0.0, float(components[key]))
        for key in (
            "single_spurs_positive_sum_kcal_mol",
            "spurs_multi_ddg_kcal_mol",
            "foldx_ddg_kcal_mol",
        )
        if components.get(key) is not None
    ]
    stability = min(1.0, max(positive_stability, default=max(0.0, candidate.structural_risk)) / 5.0)
    glycan = max(
        float(components.get("glycan_shielded_mutation_fraction") or 0.0),
        float(components.get("hypothesis_near_glycan_fraction") or 0.0),
        float(components.get("glycan_prohibition_mutation_fraction") or 0.0),
    )
    pae = components.get("response_weighted_median_local_pae_a")
    uncertainty = max(
        float(components.get("hypothesis_low_confidence_surface_fraction") or 0.0),
        min(1.0, max(0.0, float(pae)) / 30.0) if pae is not None else 0.0,
    )
    risk_values = {
        "stability": stability,
        "glycan": glycan,
        "disulfide": float(components.get("disulfide_mutation_fraction") or 0.0),
        "native_interface": float(components.get("native_interface_mutation_fraction") or 0.0),
        "structure_uncertainty": uncertainty,
    }
    penalties = {
        "stability": 0.04,
        "glycan": 0.05,
        "disulfide": 0.08,
        "native_interface": 0.03,
        "structure_uncertainty": 0.04,
    }
    penalties.update({key: max(0.0, float(value)) for key, value in (weights or {}).items()})
    return sum(risk_values[key] * penalties.get(key, 0.0) for key in risk_values)


def _diverse_take(
    pool: list[V4Candidate], count: int, already_selected: list[V4Candidate]
) -> list[V4Candidate]:
    if count <= 0:
        return []
    ranked = sorted(pool, key=_candidate_rank)
    used_candidates = {item.candidate_id for item in already_selected}
    used_positions = {_position_signature(item) for item in already_selected}
    chosen: list[V4Candidate] = []
    for candidate in ranked:
        signature = _position_signature(candidate)
        if candidate.candidate_id in used_candidates or signature in used_positions:
            continue
        chosen.append(candidate)
        used_candidates.add(candidate.candidate_id)
        used_positions.add(signature)
        if len(chosen) >= count:
            return chosen
    for candidate in ranked:
        if candidate.candidate_id in used_candidates:
            continue
        chosen.append(candidate)
        used_candidates.add(candidate.candidate_id)
        if len(chosen) >= count:
            break
    return chosen


def select_panel(
    candidates: list[V4Candidate],
    max_primary: int = 7,
    higher_order_max: int = 2,
    recommendations_per_sector: int = 6,
    epitope_priority_per_sector: int = 3,
    stability_fallback_per_sector: int = 3,
    sectors: list[V4Sector] | None = None,
    surface_areas: dict[int, float] | None = None,
    hypotheses: list[V4SurfaceHypothesis] | None = None,
    risk_penalty_weights: dict[str, float] | None = None,
    panel_budget: dict | None = None,
    stop_parameters: dict | None = None,
) -> tuple[V4Panel, list[V4Candidate]]:
    return _select_panel(candidates, max_primary, higher_order_max, recommendations_per_sector,
                         epitope_priority_per_sector, stability_fallback_per_sector,
                         sectors, surface_areas, hypotheses, risk_penalty_weights,
                         panel_budget, stop_parameters)


def _select_panel(
    candidates, max_primary=7, higher_order_max=2, recommendations_per_sector=6,
    epitope_priority_per_sector=3, stability_fallback_per_sector=3,
    sectors=None, surface_areas=None, hypotheses=None, risk_penalty_weights=None,
    panel_budget=None, stop_parameters=None,
):
    candidates = [c for c in candidates if c.qualification != "REVIEW_REQUIRED" and not any(m.review_required for m in c.mutations)]
    # Legacy geometry spread and V4.1 objectives share review isolation.
    by_sector: dict[str, list[V4Candidate]] = defaultdict(list)
    for candidate in candidates:
        if candidate.connected:
            by_sector[candidate.sector_id].append(candidate)

    sector_series: dict[str, list[V4Candidate]] = {}
    per_sector = max(1, int(recommendations_per_sector))
    for sector_id, sector_candidates in sorted(by_sector.items()):
        primary_candidates = [
            item for item in sector_candidates
            if item.order < 4 and item.search_role != "EXPLORATORY"
        ]
        prioritized = [
            item
            for item in primary_candidates
            if item.recommendation_class == "EPITOPE_PRIORITIZED"
        ]
        fallback = [
            item
            for item in primary_candidates
            if item.recommendation_class == "STABILITY_FALLBACK"
        ]
        series: list[V4Candidate] = []
        series.extend(
            _diverse_take(prioritized, min(per_sector, epitope_priority_per_sector), series)
        )
        series.extend(
            _diverse_take(
                fallback,
                min(per_sector - len(series), stability_fallback_per_sector),
                series,
            )
        )
        if len(series) < per_sector:
            series.extend(_diverse_take(primary_candidates, per_sector - len(series), series))
        sector_series[sector_id] = series[:per_sector]

    selected: list[V4Candidate] = []
    mean_response = mean_separation = None
    prior_weighted_response = None
    marginal_information_gain: dict[str, float] = {}
    indistinguishable = None
    minimum_distance = None
    prior_weighted_separation = None
    unresolved_groups: list[list[str]] = []
    uncovered_hypotheses: list[str] = []
    hypothesis_signatures: dict[str, list[float]] = {}
    budget_decision = None
    nearby_pairs = []
    identifiability_diagnostics = {}
    if panel_budget is not None or surface_areas is not None or hypotheses:
        from .identifiability import greedy_panel
        effective_budget = panel_budget or dict(
            minimum_constructs=min(4, max_primary), target_constructs=max_primary,
            maximum_constructs=max_primary, adaptive_budget=False)
        eligible_areas = ({p: surface_areas.get(p, 0.)
                          for sector in sectors for p in sector.residue_positions}
                         if surface_areas is not None and sectors is not None else surface_areas)
        selected, metrics, budget_decision = greedy_panel(
            candidates, hypotheses or [], effective_budget, max(0, higher_order_max),
            _risk_penalty, risk_penalty_weights, stop_parameters, eligible_areas,
        )
        mean_response = metrics["coverage"]
        prior_weighted_response = metrics["weighted_coverage"]
        mean_separation = metrics["mean_separation"]
        minimum_distance = metrics["minimum_distance"]
        prior_weighted_separation = metrics["weighted_separation"]
        nearby_pairs = metrics["near_pairs"]
        identifiability_diagnostics = metrics["identifiability_diagnostics"]
        indistinguishable = len(nearby_pairs)
        unresolved_groups = metrics["groups"]
        uncovered_hypotheses = metrics["uncovered"]
        hypothesis_signatures = metrics["signatures"]
        marginal_information_gain = {row["candidate_id"]: row["information_gain"] for row in budget_decision.budget_curve}
    high_order_count = 0
    offset = 0
    while not hypotheses and panel_budget is None and budget_decision is None and len(selected) < max(1, int(max_primary)):
        added = False
        for sector_id in sorted(sector_series):
            series = sector_series[sector_id]
            if offset >= len(series):
                continue
            candidate = series[offset]
            if candidate.order > 2:
                if high_order_count >= max(0, int(higher_order_max)):
                    continue
                high_order_count += 1
            selected.append(candidate)
            added = True
            if len(selected) >= max(1, int(max_primary)):
                break
        if not added:
            break
        offset += 1

    roles: dict[str, list[str]] = {"PRIMARY": []}
    for candidate in selected:
        candidate.role = candidate.recommendation_class
        roles["PRIMARY"].append(candidate.candidate_id)
        roles.setdefault(candidate.recommendation_class, []).append(candidate.candidate_id)
    denominator = (
        {position for sector in sectors for position in sector.residue_positions}
        if sectors is not None
        else {
            position
            for candidate in candidates
            for positions in candidate.covered_positions_by_scenario.values()
            for position in positions
        }
    )
    coverage_positions: dict[str, set[int]] = defaultdict(set)
    for candidate in selected:
        for scenario, positions in candidate.covered_positions_by_scenario.items():
            coverage_positions[scenario].update(set(positions) & denominator)
    coverage = {
        scenario: round(len(positions) / len(denominator), 6) if denominator else 0.0
        for scenario, positions in sorted(coverage_positions.items())
    }
    areas = surface_areas or {}
    total_area = sum(max(0.0, areas.get(position, 0.0)) for position in denominator)
    area_coverage = {
        scenario: round(
            sum(max(0.0, areas.get(position, 0.0)) for position in positions) / total_area,
            6,
        ) if total_area else 0.0
        for scenario, positions in sorted(coverage_positions.items())
    }
    selected_sectors = {item.sector_id for item in selected}
    uncovered = [
        {"sector_id": sector_id, "reason": "NO_SELECTED_FEASIBLE_CANDIDATE"}
        for sector_id in sorted(by_sector)
        if sector_id not in selected_sectors
    ]
    panel = V4Panel(
        selected_candidate_ids=[item.candidate_id for item in selected],
        candidates_by_role=roles,
        achieved_coverage=coverage,
        uncovered_components=uncovered,
        solver_status=("COVERAGE_FIRST_GREEDY" if budget_decision else
                       "FOOTPRINT_IDENTIFIABILITY_GREEDY" if hypotheses else "SECTOR_SERIES_DETERMINISTIC"),
        coverage_denominator_residues=len(denominator),
        covered_residues_by_scenario={
            scenario: len(positions) for scenario, positions in sorted(coverage_positions.items())
        },
        coverage_denominator_area_a2=round(total_area, 3),
        achieved_area_coverage=area_coverage,
        mean_hypothesis_response=round(mean_response, 6) if mean_response is not None else None,
        mean_pairwise_separation=round(mean_separation, 6) if mean_separation is not None else None,
        indistinguishable_hypothesis_pairs=indistinguishable,
        minimum_signature_distance=round(minimum_distance, 6) if minimum_distance is not None else None,
        unresolved_hypothesis_groups=unresolved_groups,
        uncovered_hypotheses=uncovered_hypotheses,
        hypothesis_signatures=hypothesis_signatures,
        prior_weighted_mean_pairwise_separation=(
            round(prior_weighted_separation, 6) if prior_weighted_separation is not None else None
        ),
        prior_weighted_hypothesis_response=(
            round(prior_weighted_response, 6) if prior_weighted_response is not None else None
        ),
        marginal_information_gain=marginal_information_gain,
        budget_decision=budget_decision,
        nearby_hypothesis_pairs=nearby_pairs,
        identifiability_diagnostics=identifiability_diagnostics,
    )
    return panel, selected
