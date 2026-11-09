from __future__ import annotations

import heapq

from .models import MicroFootprint, ProbeCandidate, SurfaceResidueUnit


def _shortest_paths(graph: dict[int, dict[int, float]], start: int) -> dict[int, float]:
    distances = {start: 0.0}
    queue: list[tuple[float, int]] = [(0.0, start)]
    while queue:
        distance, current = heapq.heappop(queue)
        if distance != distances.get(current):
            continue
        for neighbor, step in sorted(graph.get(current, {}).items()):
            proposed = distance + step
            if proposed < distances.get(neighbor, float("inf")):
                distances[neighbor] = proposed
                heapq.heappush(queue, (proposed, neighbor))
    return distances


def attach_coverage(
    candidates: list[ProbeCandidate],
    units: list[SurfaceResidueUnit],
    graph: dict[int, dict[int, float]],
    footprints: list[MicroFootprint],
    distance_scenarios_a: dict[str, float] | None = None,
) -> list[ProbeCandidate]:
    positions = {unit.position for unit in units}
    areas = {unit.position: unit.mapped_area_a2 for unit in units}
    designable = {unit.position for unit in units if unit.designable}
    id_to_position = {unit.residue: unit.position for unit in units}
    footprint_position_sets = [
        {id_to_position[item] for item in footprint.residue_ids if item in id_to_position}
        for footprint in footprints
    ]
    designable_total = max(1.0, sum(areas[position] for position in designable))
    total_surface = max(1.0, sum(areas.values()))
    for candidate in candidates:
        mutations = {item.position for item in candidate.mutations}
        direct = mutations & positions
        distances: dict[int, float] = {}
        for mutation in mutations:
            distances.update({position: min(value, distances.get(position, float("inf"))) for position, value in _shortest_paths(graph, mutation).items()})
        scenario_distances = distance_scenarios_a or {
            "geodesic_tight": 9.0,
            "geodesic_medium": 12.0,
            "geodesic_wide": 15.0,
        }
        for scenario, cutoff in scenario_distances.items():
            covered = {position for position, distance in distances.items() if distance <= cutoff and position in designable}
            candidate.coverage_units_by_scenario[scenario] = sorted(covered)
            candidate.coverage_by_scenario[scenario] = {
                "designable_area_a2": round(sum(areas[position] for position in covered), 3),
                "fraction_of_designable_surface": round(sum(areas[position] for position in covered) / designable_total, 4),
                "fraction_of_total_surface": round(sum(areas[position] for position in covered) / total_surface, 4),
            }
        candidate.coverage_units_by_scenario["direct_site"] = sorted(direct & designable)
        candidate.coverage_by_scenario["direct_site"] = {
            "designable_area_a2": round(sum(areas[position] for position in direct & designable), 3),
            "fraction_of_designable_surface": round(sum(areas[position] for position in direct & designable) / designable_total, 4),
            "fraction_of_total_surface": round(sum(areas[position] for position in direct & designable) / total_surface, 4),
        }
        for name, minimum in (("footprint_hit_1", 1), ("footprint_hit_2", 2)):
            hit_count = 0
            covered_by_footprints: set[int] = set()
            for footprint_positions in footprint_position_sets:
                hits = len(mutations & footprint_positions)
                if hits >= minimum:
                    hit_count += 1
                    covered_by_footprints.update(footprint_positions & designable)
            hit_area = sum(areas.get(position, 0.0) for position in covered_by_footprints)
            candidate.coverage_by_scenario[name] = {
                "footprint_count": float(hit_count),
                "designable_area_a2": round(hit_area, 3),
                "fraction_of_designable_surface": round(hit_area / designable_total, 4),
                "fraction_of_total_surface": round(hit_area / total_surface, 4),
            }
            candidate.coverage_units_by_scenario[name] = sorted(covered_by_footprints)
    return candidates
