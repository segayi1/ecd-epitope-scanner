from __future__ import annotations

from collections import defaultdict

from .models import V4Candidate, V4Panel, V4Sector


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
) -> tuple[V4Panel, list[V4Candidate]]:
    """Select a stability-gated recommendation series for every sector."""
    by_sector: dict[str, list[V4Candidate]] = defaultdict(list)
    for candidate in candidates:
        if candidate.connected:
            by_sector[candidate.sector_id].append(candidate)

    sector_series: dict[str, list[V4Candidate]] = {}
    per_sector = max(1, int(recommendations_per_sector))
    for sector_id, sector_candidates in sorted(by_sector.items()):
        prioritized = [
            item
            for item in sector_candidates
            if item.recommendation_class == "EPITOPE_PRIORITIZED"
        ]
        fallback = [
            item
            for item in sector_candidates
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
            series.extend(_diverse_take(sector_candidates, per_sector - len(series), series))
        sector_series[sector_id] = series[:per_sector]

    selected: list[V4Candidate] = []
    high_order_count = 0
    offset = 0
    while len(selected) < max(1, int(max_primary)):
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
        solver_status="SECTOR_SERIES_DETERMINISTIC",
        coverage_denominator_residues=len(denominator),
        covered_residues_by_scenario={
            scenario: len(positions) for scenario, positions in sorted(coverage_positions.items())
        },
        coverage_denominator_area_a2=round(total_area, 3),
        achieved_area_coverage=area_coverage,
    )
    return panel, selected
