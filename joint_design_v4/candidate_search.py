from __future__ import annotations

from itertools import combinations, product
from typing import Any

from .models import V4Candidate, V4MutationOption, V4Sector


def _distances(
    graph: dict[int, dict[int, float]], start: int, allowed: set[int]
) -> dict[int, float]:
    import heapq

    result = {start: 0.0}
    queue = [(0.0, start)]
    while queue:
        distance, current = heapq.heappop(queue)
        if distance != result.get(current):
            continue
        for neighbor, step in graph.get(current, {}).items():
            if neighbor not in allowed:
                continue
            proposed = distance + step
            if proposed < result.get(neighbor, float("inf")):
                result[neighbor] = proposed
                heapq.heappush(queue, (proposed, neighbor))
    return result


def _connected(positions: set[int], graph: dict[int, dict[int, float]], allowed: set[int]) -> bool:
    if not positions:
        return False
    seen = {min(positions)}
    changed = True
    while changed:
        changed = False
        for left in tuple(seen):
            reachable = _distances(graph, left, allowed)
            newly_seen = positions.intersection(reachable).difference(seen)
            if newly_seen:
                seen.update(newly_seen)
                changed = True
    return seen == positions


def generate_candidates(
    sectors: list[V4Sector],
    options: list[V4MutationOption],
    graph: dict[int, dict[int, float]],
    parameters: dict[str, Any],
) -> list[V4Candidate]:
    by_position: dict[int, list[V4MutationOption]] = {}
    for option in options:
        by_position.setdefault(option.residue.uniprot_pos or -1, []).append(option)
    orders = [int(item) for item in parameters.get("orders", [2, 3, 4])]
    maximum_span = float(parameters.get("maximum_geodesic_span_a", 24.0))
    site_cap = int(parameters.get("site_pool_per_sector", 12))
    per_sector = int(parameters.get("shortlist_per_sector", 30))
    epitope_site_fraction = min(
        1.0, max(0.0, float(parameters.get("epitope_site_pool_fraction", 0.5)))
    )
    fallback_candidate_fraction = min(
        1.0, max(0.0, float(parameters.get("fallback_candidate_fraction", 0.5)))
    )
    output = []
    for sector in sectors:
        positions = [position for position in sector.residue_positions if position in by_position]
        positions.sort(key=lambda p: (-max(option.score for option in by_position[p]), p))
        supported_positions = [
            position
            for position in positions
            if any(option.epitope_supported for option in by_position[position])
        ]
        fallback_positions = [position for position in positions if position not in supported_positions]
        if supported_positions and fallback_positions:
            supported_cap = min(
                len(supported_positions), max(1, round(site_cap * epitope_site_fraction))
            )
            fallback_cap = min(len(fallback_positions), max(1, site_cap - supported_cap))
            positions = supported_positions[:supported_cap] + fallback_positions[:fallback_cap]
            remaining = [
                position
                for position in supported_positions[supported_cap:] + fallback_positions[fallback_cap:]
                if position not in positions
            ]
            positions.extend(remaining[: max(0, site_cap - len(positions))])
        else:
            positions = positions[:site_cap]
        candidates: list[V4Candidate] = []
        for order in orders:
            if len(positions) < order:
                continue
            for site_combo in combinations(positions, order):
                if not _connected(set(site_combo), graph, set(sector.residue_positions)):
                    continue
                distances = [
                    _distances(graph, left, set(sector.residue_positions)) for left in site_combo
                ]
                span = max(
                    (
                        distances[i].get(right, float("inf"))
                        for i, left in enumerate(site_combo)
                        for right in site_combo[i + 1 :]
                    ),
                    default=0.0,
                )
                if span > maximum_span:
                    continue
                pools = [
                    sorted(by_position[position], key=lambda item: (-item.score, item.notation))[:3]
                    for position in site_combo
                ]
                for selected in product(*pools):
                    if len({item.mutant_aa for item in selected}) == 0:
                        continue
                    risk_values = [
                        item.spurs_ddg_kcal_mol
                        for item in selected
                        if item.spurs_ddg_kcal_mol is not None
                    ]
                    candidate = V4Candidate(
                        candidate_id=f"V4_{sector.sector_id}_{order}_{len(candidates) + 1:04d}",
                        sector_id=sector.sector_id,
                        mutations=list(selected),
                        order=order,
                        epitope_score=round(sum(item.epitope_prior for item in selected), 6),
                        structural_risk=round(sum(max(0.0, value) for value in risk_values), 6),
                        connected=True,
                        warnings=(
                            []
                            if order == 2
                            else [
                                f"HIGHER_ORDER_EXTRAPOLATION_{'MEDIUM' if order == 3 else 'HIGH'}"
                            ]
                        ),
                        epitope_supported_sites=sum(
                            1 for item in selected if item.epitope_supported
                        ),
                        recommendation_class=(
                            "EPITOPE_PRIORITIZED"
                            if any(item.epitope_supported for item in selected)
                            else "STABILITY_FALLBACK"
                        ),
                    )
                    candidate.coverage_by_scenario, candidate.covered_positions_by_scenario = _coverage(
                        candidate, sector, graph
                    )
                    candidates.append(candidate)
            if len(candidates) >= int(parameters.get("beam_width", 300)):
                break
        candidates.sort(
            key=lambda item: (
                item.structural_risk,
                -item.epitope_score,
                item.order,
                item.candidate_id,
            )
        )
        prioritized = [
            item for item in candidates if item.recommendation_class == "EPITOPE_PRIORITIZED"
        ]
        fallback = [
            item for item in candidates if item.recommendation_class == "STABILITY_FALLBACK"
        ]
        if prioritized and fallback:
            fallback_cap = min(
                len(fallback), max(1, round(per_sector * fallback_candidate_fraction))
            )
            priority_cap = min(len(prioritized), max(1, per_sector - fallback_cap))
            shortlisted = prioritized[:priority_cap] + fallback[:fallback_cap]
            remaining = [
                item
                for item in prioritized[priority_cap:] + fallback[fallback_cap:]
                if item not in shortlisted
            ]
            shortlisted.extend(remaining[: max(0, per_sector - len(shortlisted))])
        else:
            shortlisted = candidates[:per_sector]
        output.extend(shortlisted[:per_sector])
    return sorted(output, key=lambda item: (item.sector_id, item.order, item.candidate_id))


def _coverage(
    candidate: V4Candidate, sector: V4Sector, graph: dict[int, dict[int, float]]
) -> tuple[dict[str, float], dict[str, list[int]]]:
    positions = set(sector.residue_positions)
    covered: dict[str, float] = {}
    covered_positions: dict[str, list[int]] = {}
    for radius in (0.0, 9.0, 12.0, 15.0):
        hit = set()
        for mutation in candidate.mutations:
            distance = _distances(graph, mutation.residue.uniprot_pos or -1, positions)
            hit |= {position for position, value in distance.items() if value <= radius}
        scenario = f"geodesic_{int(radius)}"
        covered[scenario] = round(len(hit) / max(1, len(positions)), 6)
        covered_positions[scenario] = sorted(hit)
    return covered, covered_positions
