from __future__ import annotations

import heapq
from collections.abc import Iterable

from .models import MicroFootprint, SurfaceResidueUnit


def _distances(graph: dict[int, dict[int, float]], start: int) -> dict[int, float]:
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


def sample_centers(
    positions: Iterable[int],
    graph: dict[int, dict[int, float]],
    areas: dict[int, float],
    max_centers: int,
) -> list[int]:
    positions = sorted(set(positions))
    if not positions or max_centers <= 0:
        return []
    selected = [max(positions, key=lambda position: (areas.get(position, 0.0), -position))]
    distances_cache = {selected[0]: _distances(graph, selected[0])}
    while len(selected) < min(max_centers, len(positions)):
        best = None
        for position in positions:
            if position in selected:
                continue
            distances = [distances_cache[center].get(position, float("inf")) for center in selected]
            key = (min(distances), areas.get(position, 0.0), -position)
            if best is None or key > best[0]:
                best = (key, position)
        if best is None:
            break
        selected.append(best[1])
        distances_cache[best[1]] = _distances(graph, best[1])
    return selected


def build_microfootprints(
    units: list[SurfaceResidueUnit],
    graph: dict[int, dict[int, float]],
    areas: dict[int, float],
    target_areas: list[float],
    max_centers: int,
) -> list[MicroFootprint]:
    positions = [unit.position for unit in units]
    centers = sample_centers(positions, graph, areas, max_centers)
    by_position = {unit.position: unit for unit in units}
    footprints: list[MicroFootprint] = []
    for center in centers:
        distances = _distances(graph, center)
        ordered = sorted(positions, key=lambda position: (distances.get(position, float("inf")), -areas.get(position, 0.0), position))
        for index, target in enumerate(sorted(float(value) for value in target_areas), start=1):
            selected: list[int] = []
            total = 0.0
            for position in ordered:
                if position not in distances:
                    continue
                selected.append(position)
                total += areas.get(position, 1.0)
                if total >= target:
                    break
            if len(selected) < 2:
                continue
            flags = []
            if total < target:
                flags.append("boundary_truncated")
            footprints.append(
                MicroFootprint(
                    footprint_id=f"MF{len(footprints) + 1:04d}",
                    center=by_position[center].residue,
                    scale_id=f"area_{int(target)}",
                    target_area_a2=target,
                    residue_ids=[by_position[position].residue for position in selected],
                    realized_area_a2=round(total, 3),
                    boundary_truncated=total < target,
                    quality_flags=flags,
                )
            )
    return footprints
