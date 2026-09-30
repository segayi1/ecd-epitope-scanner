from __future__ import annotations

from itertools import combinations, product
from typing import Any

from .models import V4Candidate, V4MutationOption, V4Sector, V4SurfaceHypothesis


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
    hypotheses: list[V4SurfaceHypothesis] | None = None,
) -> list[V4Candidate]:
    by_position: dict[int, list[V4MutationOption]] = {}
    for option in options:
        by_position.setdefault(option.residue.uniprot_pos or -1, []).append(option)
    orders = list(dict.fromkeys(
        [int(item) for item in parameters.get("orders", [2, 3])]
        + [int(item) for item in parameters.get("exploratory_orders", [])]
    ))
    maximum_span = float(parameters.get("maximum_geodesic_span_a", 24.0))
    site_cap = int(parameters.get("site_pool_per_sector", 12))
    per_sector = int(parameters.get("shortlist_per_sector", 30))
    epitope_site_fraction = min(
        1.0, max(0.0, float(parameters.get("epitope_site_pool_fraction", 0.5)))
    )
    fallback_candidate_fraction = min(
        1.0, max(0.0, float(parameters.get("fallback_candidate_fraction", 0.5)))
    )
    display_owner = {
        position: sector.sector_id
        for sector in sectors for position in sector.residue_positions
    }
    if hypotheses:
        windows = []
        for hypothesis in hypotheses:
            covered = [position for position in hypothesis.residue_weights if position in by_position and position in graph]
            all_positions = sorted(position for position in hypothesis.residue_weights if position in graph)
            if not covered or not all_positions:
                continue
            owner_counts: dict[str, int] = {}
            for position in all_positions:
                owner = display_owner.get(position)
                if owner:
                    owner_counts[owner] = owner_counts.get(owner, 0) + 1
            owner = min(owner_counts, key=lambda item: (-owner_counts[item], item)) if owner_counts else "UNASSIGNED"
            window = V4Sector(
                sector_id=owner, residue_positions=all_positions, component_id=None,
                designable_area_a2=sum(hypothesis.residue_weights[p] for p in covered),
                total_area_a2=sum(hypothesis.residue_weights[p] for p in all_positions),
                epitope_score=hypothesis.discotope_prior or 0.0,
            )
            windows.append((window, hypothesis.hypothesis_id))
    else:
        windows = [(sector, None) for sector in sectors]
    output = []
    for sector, hypothesis_id in windows:
        positions = [position for position in sector.residue_positions if position in by_position]
        positions.sort(key=lambda p: (
            -max(option.accessibility_score for option in by_position[p]),
            -max(option.chemical_disruption_score for option in by_position[p]),
            min(option.spurs_ddg_kcal_mol if option.spurs_ddg_kcal_mol is not None else float("inf") for option in by_position[p]),
            p,
        ))
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
        order_shortlists: list[V4Candidate] = []
        if hypotheses:
            positions = positions[:min(site_cap, int(parameters.get("site_pool_per_hypothesis", 8)))]
        for order in orders:
            candidates = []
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
                    sorted(by_position[position], key=lambda item: (
                        -item.chemical_disruption_score,
                        item.spurs_ddg_kcal_mol if item.spurs_ddg_kcal_mol is not None else float("inf"),
                        -item.accessibility_score,
                        item.notation,
                    ))[:3]
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
                        candidate_id=(
                            f"V4_{sector.sector_id}_{hypothesis_id}_{order}_{len(candidates) + 1:04d}"
                            if hypothesis_id else
                            f"V4_{sector.sector_id}_{order}_{len(candidates) + 1:04d}"
                        ),
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
                        search_role="EXPLORATORY" if order >= 4 else "PRIMARY_POOL",
                        targeted_hypotheses=[hypothesis_id] if hypothesis_id else [],
                    )
                    candidate.coverage_by_scenario, candidate.covered_positions_by_scenario = _coverage(
                        candidate, sector, graph
                    )
                    candidates.append(candidate)
            candidates.sort(
                key=lambda item: (
                    item.structural_risk,
                    -item.epitope_score,
                    item.order,
                    item.candidate_id,
                )
            )
            beam_by_order = parameters.get("beam_width_by_order", {})
            beam_cap = beam_by_order.get(str(order), beam_by_order.get(order, parameters.get("beam_width", 300)))
            candidates = candidates[: max(1, int(beam_cap))]
            prioritized = [
                item for item in candidates if item.recommendation_class == "EPITOPE_PRIORITIZED"
            ]
            fallback = [
                item for item in candidates if item.recommendation_class == "STABILITY_FALLBACK"
            ]
            shortlist_by_order = parameters.get("shortlist_by_order", {})
            default_order_cap = parameters.get(
                "shortlist_per_order", max(1, per_sector // max(1, len(orders)))
            )
            order_cap = max(
                1,
                int(shortlist_by_order.get(str(order), shortlist_by_order.get(order, default_order_cap))),
            )
            if prioritized and fallback:
                fallback_cap = min(
                    len(fallback), max(1, round(order_cap * fallback_candidate_fraction))
                )
                fallback_cap = min(fallback_cap, max(1, order_cap - 1))
                priority_cap = min(len(prioritized), max(1, order_cap - fallback_cap))
                shortlisted = prioritized[:priority_cap] + fallback[:fallback_cap]
                remaining = [
                    item
                    for item in prioritized[priority_cap:] + fallback[fallback_cap:]
                    if item not in shortlisted
                ]
                shortlisted.extend(remaining[: max(0, order_cap - len(shortlisted))])
            else:
                shortlisted = candidates[:order_cap]
            order_shortlists.extend(shortlisted[:order_cap])
        output.extend(order_shortlists)
    if hypotheses:
        grouped: dict[tuple[str, int], list[V4Candidate]] = {}
        for candidate in output:
            grouped.setdefault((candidate.sector_id, candidate.order), []).append(candidate)
        pooled = []
        for (_sector_id, order), items in sorted(grouped.items()):
            by_signature: dict[tuple[int, ...], V4Candidate] = {}
            for candidate in sorted(items, key=lambda item: (
                item.structural_risk, -item.epitope_score, item.candidate_id
            )):
                signature = tuple(sorted(option.residue.uniprot_pos or -1 for option in candidate.mutations))
                existing = by_signature.get(signature)
                if existing is None:
                    by_signature[signature] = candidate
                elif candidate.targeted_hypotheses:
                    existing.targeted_hypotheses = sorted(set(existing.targeted_hypotheses + candidate.targeted_hypotheses))
            shortlist_by_order = parameters.get("shortlist_by_order", {})
            cap = max(1, int(shortlist_by_order.get(str(order), shortlist_by_order.get(order, default_order_cap))))
            ranked = sorted(by_signature.values(), key=lambda item: (
                item.structural_risk, -item.epitope_score, item.candidate_id
            ))
            prioritized = [item for item in ranked if item.recommendation_class == "EPITOPE_PRIORITIZED"]
            fallback = [item for item in ranked if item.recommendation_class == "STABILITY_FALLBACK"]
            if prioritized and fallback:
                fcap = min(len(fallback), max(1, round(cap * fallback_candidate_fraction)))
                selected = prioritized[:max(1, cap - fcap)] + fallback[:fcap]
                selected_ids = {item.candidate_id for item in selected}
                selected.extend(item for item in ranked if item.candidate_id not in selected_ids and len(selected) < cap)
            else:
                selected = ranked[:cap]
            pooled.extend(selected[:cap])
        output = pooled
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
