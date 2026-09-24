from __future__ import annotations

import math
from collections import deque
from typing import Any

import numpy as np

from .models import EpitopeEvidence, SurfaceEvidence, V4Sector


def _components(graph: dict[int, dict[int, float]]) -> list[set[int]]:
    remaining = set(graph)
    output = []
    while remaining:
        start = min(remaining)
        group = {start}
        queue = deque([start])
        remaining.remove(start)
        while queue:
            current = queue.popleft()
            for neighbor in graph.get(current, {}):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    group.add(neighbor)
                    queue.append(neighbor)
        output.append(group)
    return output


def _distance(
    graph: dict[int, dict[int, float]], start: int, allowed: set[int]
) -> dict[int, float]:
    import heapq

    distances = {start: 0.0}
    queue = [(0.0, start)]
    while queue:
        distance, current = heapq.heappop(queue)
        if distance != distances.get(current):
            continue
        for neighbor, step in graph.get(current, {}).items():
            if neighbor not in allowed:
                continue
            proposed = distance + step
            if proposed < distances.get(neighbor, float("inf")):
                distances[neighbor] = proposed
                heapq.heappush(queue, (proposed, neighbor))
    return distances


def build_epitope_sectors(
    surface: list[SurfaceEvidence],
    epitope: list[EpitopeEvidence],
    graph: dict[int, dict[int, float]],
    records: dict[int, Any],
    max_sectors: int = 7,
    parameters: dict[str, Any] | None = None,
) -> list[V4Sector]:
    params = parameters or {}
    min_area = float(params.get("minimum_sector_area_a2", 50.0))
    min_residues = max(1, int(params.get("minimum_sector_residues", 4)))
    max_geodesic = float(params.get("maximum_geodesic_diameter_a", 35.0))
    max_euclidean = float(params.get("maximum_euclidean_diameter_a", 30.0))
    by_position = {item.residue.uniprot_pos: item for item in surface}
    epi = {item.residue.uniprot_pos: item for item in epitope}
    components = _components(graph)
    component_by_position = {
        position: index for index, component in enumerate(components) for position in component
    }
    distances = {
        position: _distance(graph, position, components[component_by_position[position]])
        for position in graph
    }

    def prior(position: int) -> float:
        value = epi.get(position)
        return float(value.surface_percentile or 0.0) if value else 0.0

    # Allocate at least one center to every physical surface component. Within
    # the remaining budget use deterministic farthest-point sampling. This is
    # a graph k-center approximation: it controls sector diameter directly and
    # avoids allowing a one-residue epitope island to become its own sector.
    centers = [
        max(
            component,
            key=lambda position: (
                prior(position),
                by_position.get(position).mapped_surface_area_a2
                if position in by_position
                else 0.0,
                -position,
            ),
        )
        for component in components
        if component
    ]
    target_count = min(len(graph), max(len(centers), max(1, max_sectors)))
    while len(centers) < target_count:
        candidates = [position for position in graph if position not in centers]
        center = max(
            candidates,
            key=lambda position: (
                min(
                    distances[position].get(existing, float("inf"))
                    for existing in centers
                    if component_by_position[existing] == component_by_position[position]
                ),
                prior(position),
                by_position.get(position).mapped_surface_area_a2
                if position in by_position
                else 0.0,
                -position,
            ),
        )
        centers.append(center)

    # A few deterministic minimax-medoid iterations improve balance without
    # changing the requested sector count or introducing random seeds.
    groups_by_center: dict[int, set[int]] = {}
    for _ in range(8):
        groups_by_center = {center: set() for center in centers}
        for position in sorted(graph):
            owner = min(
                (
                    center
                    for center in centers
                    if component_by_position[center] == component_by_position[position]
                ),
                key=lambda center: (distances[center].get(position, float("inf")), center),
            )
            groups_by_center[owner].add(position)
        new_centers = []
        for center in centers:
            group = groups_by_center[center]
            new_centers.append(
                min(
                    group,
                    key=lambda position: (
                        max(distances[position].get(other, float("inf")) for other in group),
                        sum(distances[position].get(other, float("inf")) for other in group),
                        -prior(position),
                        position,
                    ),
                )
            )
        if new_centers == centers:
            break
        centers = new_centers
    groups_by_center = {center: set() for center in centers}
    for position in sorted(graph):
        owner = min(
            (
                center
                for center in centers
                if component_by_position[center] == component_by_position[position]
            ),
            key=lambda center: (distances[center].get(position, float("inf")), center),
        )
        groups_by_center[owner].add(position)
    groups = [set(group) for group in groups_by_center.values() if group]

    def adjacent(left: set[int], right: set[int]) -> bool:
        return any(neighbor in right for position in left for neighbor in graph.get(position, {}))

    # Merge graph-adjacent undersized groups before enforcing the requested
    # maximum. A one-residue epitope island is evidence attached to a surface
    # neighborhood, not a useful standalone experimental sector.
    groups.sort(
        key=lambda group: (
            -sum(by_position.get(p).mapped_surface_area_a2 for p in group if p in by_position),
            min(group),
        )
    )
    while True:
        mergeable = [
            index
            for index, group in enumerate(groups)
            if any(adjacent(group, other) for other in groups if other is not group)
        ]
        if not mergeable:
            break
        undersized = [
            index
            for index in mergeable
            if len(groups[index]) < min_residues
            or sum(
                by_position.get(p).mapped_surface_area_a2
                for p in groups[index]
                if p in by_position
            )
            < min_area
        ]
        if len(groups) <= max(1, max_sectors) and not undersized:
            break
        pool = undersized or mergeable
        tiny_index = min(
            pool,
            key=lambda i: (
                sum(
                    by_position.get(p).mapped_surface_area_a2 for p in groups[i] if p in by_position
                ),
                i,
            ),
        )
        tiny = groups.pop(tiny_index)
        adjacent_targets = [i for i, group in enumerate(groups) if adjacent(tiny, group)]
        target_index = min(
            adjacent_targets,
            key=lambda i: (
                -sum(
                    1
                    for position in tiny
                    for neighbor in graph.get(position, {})
                    if neighbor in groups[i]
                ),
                min(
                    _distance(graph, position, tiny | groups[i]).get(neighbor, float("inf"))
                    for position in tiny
                    for neighbor in groups[i]
                ),
                -sum(by_position[p].mapped_surface_area_a2 for p in groups[i] if p in by_position),
                min(groups[i]),
            ),
        )
        groups[target_index] |= tiny
    output = []
    for index, group in enumerate(sorted(groups, key=lambda value: min(value)), start=1):
        area = sum(by_position[p].mapped_surface_area_a2 for p in group if p in by_position)
        designable = sum(
            by_position[p].mapped_surface_area_a2
            for p in group
            if p in by_position and p in records and not records[p].hard_exclusions
        )
        scores = [
            epi[p].surface_percentile
            for p in group
            if p in epi and epi[p].surface_percentile is not None
        ]
        islands = sorted({epi[p].island_id for p in group if p in epi and epi[p].island_id})
        flags = [] if islands else ["LOW_EPITOPE_PRIOR"]
        if area < min_area:
            flags.append("LOW_SURFACE_AREA")
        # Diameter is a molecular-surface geodesic and therefore uses the full
        # surface graph; artificial sector boundaries must not lengthen it.
        group_distances = {position: distances[position] for position in group}
        induced_reachable = _distance(graph, min(group), group)
        geodesic_diameter = max(
            (
                group_distances[left].get(right, float("inf"))
                for left in group
                for right in group
            ),
            default=0.0,
        )
        coordinates = [
            np.asarray(records[position].ca, dtype=float)
            for position in group
            if position in records and records[position].ca is not None
        ]
        euclidean_diameter = max(
            (
                float(np.linalg.norm(left - right))
                for left in coordinates
                for right in coordinates
            ),
            default=0.0,
        )
        if len(induced_reachable) != len(group):
            flags.append("DISCONNECTED_SECTOR")
        if not math.isfinite(geodesic_diameter):
            flags.append("GEODESIC_DISTANCE_UNAVAILABLE")
        elif geodesic_diameter > max_geodesic:
            flags.append("GEODESIC_DIAMETER_EXCEEDED")
        if euclidean_diameter > max_euclidean:
            flags.append("EUCLIDEAN_DIAMETER_EXCEEDED")
        center = min(
            group,
            key=lambda position: (
                max(group_distances[position].get(other, float("inf")) for other in group),
                sum(group_distances[position].get(other, float("inf")) for other in group),
                -prior(position),
                position,
            ),
        )
        output.append(
            V4Sector(
                sector_id=f"VS{index:03d}",
                residue_positions=sorted(group),
                component_id=next(
                    (i + 1 for i, component in enumerate(components) if group & component), None
                ),
                designable_area_a2=round(designable, 3),
                total_area_a2=round(area, 3),
                epitope_score=round(sum(scores) / len(scores), 4) if scores else 0.0,
                center_position=center,
                geodesic_diameter_a=(
                    round(geodesic_diameter, 3) if math.isfinite(geodesic_diameter) else None
                ),
                euclidean_diameter_a=round(euclidean_diameter, 3),
                island_ids=islands,
                sector_type="EPITOPE_PRIOR" if islands else "COVERAGE_ONLY",
                quality_flags=flags,
            )
        )
    return output
