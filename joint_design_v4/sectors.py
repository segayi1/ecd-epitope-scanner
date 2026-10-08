from __future__ import annotations

import math
from collections import deque
from typing import Any

import numpy as np

from .models import EpitopeEvidence, SurfaceEvidence, V4Sector


PARTITION_SEARCH_LIMIT = 100_000


class PartitionConstraintError(ValueError):
    def __init__(self, message: str, diagnostics: dict[str, Any]):
        super().__init__(message)
        self.diagnostics = diagnostics


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
    diagnostics: dict[str, Any] | None = None,
) -> list[V4Sector]:
    params = parameters or {}
    min_area = float(params.get("minimum_sector_area_a2", 50.0))
    min_residues = max(1, int(params.get("minimum_sector_residues", 4)))
    max_geodesic = float(params.get("maximum_geodesic_diameter_a", 26.0))
    max_euclidean = float(params.get("maximum_euclidean_diameter_a", 18.0))
    if not all(math.isfinite(v) and v > 0 for v in (min_area, max_geodesic, max_euclidean)):
        raise ValueError("Display region area and diameter limits must be positive and finite")
    audit = diagnostics if diagnostics is not None else {}
    audit.clear()
    audit.update({
        "version": "connected_bounded_v3",
        "status": "IN_PROGRESS",
        "requested_region_count": max(1, max_sectors),
        "count_is_soft_target": True,
        "minimum_sector_residues": min_residues,
        "minimum_sector_area_a2": min_area,
        "maximum_geodesic_diameter_a": max_geodesic,
        "maximum_euclidean_diameter_a": max_euclidean,
        "search_limit": PARTITION_SEARCH_LIMIT,
        "search_states": 0,
        "search_scope": "local_boundary_reassignment; failure is not a global infeasibility proof",
        "repairs": [],
    })
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
    # a graph k-center approximation; the later split/repair and final checks
    # enforce physical extent and minimum size independently of epitope islands.
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
    coordinates_by_position = {
        p: np.asarray(records[p].ca, dtype=float)
        for p in graph if p in records and records[p].ca is not None
    }

    def bounded(group: set[int]) -> bool:
        return all(
            distances[left].get(right, float("inf")) <= max_geodesic
            and (
                left not in coordinates_by_position or right not in coordinates_by_position
                or float(np.linalg.norm(coordinates_by_position[left] - coordinates_by_position[right])) <= max_euclidean
            )
            for left in group for right in group if left < right
        )

    # The requested count is a display target. Connectivity and physical extent
    # are invariants: never join a tail and core merely to keep seven colors.
    # Growing through the induced graph also repairs disconnected Voronoi cells.
    groups = []
    for original in groups_by_center.values():
        remaining = set(original)
        while remaining:
            seed = min(remaining)
            group = {seed}
            remaining.remove(seed)
            while True:
                frontier = {
                    neighbor for position in group for neighbor in graph.get(position, {})
                    if neighbor in remaining
                }
                feasible = [p for p in frontier if bounded(group | {p})]
                if not feasible:
                    break
                chosen = min(feasible, key=lambda p: (distances[seed].get(p, float("inf")), p))
                group.add(chosen)
                remaining.remove(chosen)
            groups.append(group)

    def adjacent(left: set[int], right: set[int]) -> bool:
        return any(neighbor in right for position in left for neighbor in graph.get(position, {}))

    def area(group: set[int]) -> float:
        return sum(by_position[p].mapped_surface_area_a2 for p in sorted(group) if p in by_position)

    def too_small(group: set[int]) -> bool:
        return len(group) < min_residues or area(group) < min_area

    def fail(reason: str) -> None:
        audit.update({
            "status": "FAILED", "reason": reason,
            "draft_regions": [
                {"residue_positions": sorted(g), "residue_count": len(g),
                 "total_area_a2": area(g), "undersized": too_small(g)}
                for g in sorted(groups, key=min)
            ],
        })
        detail = ("搜索预算耗尽；不能据此判定无可行分区" if reason == "SEARCH_BUDGET_EXHAUSTED"
                  else "无法满足显示分区的最小残基数、面积、连通性或直径约束")
        raise PartitionConstraintError(f"{reason}: {detail}", audit)

    def reassign(
        proposed: list[set[int]], pending: set[int], source_region: int,
    ) -> tuple[list[set[int]], list[dict[str, int]]] | None:
        # ponytail: bounded local DFS, not a global partition solver; fail explicitly
        # at the ceiling, extend the neighborhood only if real inputs require it.
        stack = [(proposed, pending, [])]
        while stack:
            if audit["search_states"] >= PARTITION_SEARCH_LIMIT:
                fail("SEARCH_BUDGET_EXHAUSTED")
            audit["search_states"] += 1
            current, remaining, moves = stack.pop()
            if not remaining:
                return current, moves
            choices = []
            for position in sorted(remaining):
                for index, group in enumerate(current):
                    if adjacent({position}, group) and bounded(group | {position}):
                        choices.append((
                            min(distances[position].get(p, math.inf) for p in group),
                            min(group), position, index,
                        ))
            for _, destination, position, index in reversed(sorted(choices)):
                updated = [set(g) for g in current]
                updated[index].add(position)
                stack.append((updated, remaining - {position}, moves + [{
                    "position": position, "from_region": source_region,
                    "to_region": destination,
                }]))
        return None

    def rebuild_boundary(tiny: set[int], target: set[int]) -> set[int] | None:
        # If full absorption strands the far side of a target, rebuild its
        # boundary into two valid regions rather than certifying a tiny region.
        stack = [set(tiny)]
        visited = set()
        while stack:
            retained = stack.pop()
            identity = frozenset(retained)
            if identity in visited:
                continue
            if audit["search_states"] >= PARTITION_SEARCH_LIMIT:
                fail("SEARCH_BUDGET_EXHAUSTED")
            audit["search_states"] += 1
            visited.add(identity)
            remaining = target - retained
            if (not too_small(retained) and not too_small(remaining)
                    and set(_distance(graph, min(remaining), remaining)) == remaining):
                return retained
            if len(remaining) <= min_residues:
                continue
            frontier = {q for p in retained for q in graph[p] if q in remaining}
            candidates = sorted(
                (min(distances[p].get(q, math.inf) for p in tiny), q)
                for q in frontier if bounded(retained | {q})
            )
            stack.extend(retained | {q} for _, q in reversed(candidates))
        return None

    # Preserve the initial bounded whole-region consolidation; repair its
    # remaining undersized groups by adjusting boundaries below.
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
            if any(adjacent(group, other) and bounded(group | other)
                   for other in groups if other is not group)
        ]
        if not mergeable:
            break
        undersized = [index for index in mergeable if too_small(groups[index])]
        if len(groups) <= max(1, max_sectors) and not undersized:
            break
        tiny_index = min(
            undersized or mergeable,
            key=lambda i: (
                sum(
                    by_position.get(p).mapped_surface_area_a2 for p in groups[i] if p in by_position
                ),
                i,
            ),
        )
        tiny = groups.pop(tiny_index)
        adjacent_targets = [i for i, group in enumerate(groups)
                            if adjacent(tiny, group) and bounded(tiny | group)]
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

    while any(too_small(g) for g in groups):
        repaired = False
        for source_index in sorted(
            (i for i, g in enumerate(groups) if too_small(g)),
            key=lambda i: (len(groups[i]), area(groups[i]), min(groups[i])),
        ):
            tiny = groups[source_index]
            targets = sorted(
                (min(distances[p].get(q, math.inf) for p in tiny for q in group), min(group), i)
                for i, group in enumerate(groups)
                if i != source_index and adjacent(tiny, group)
            )
            for _, _, target_index in targets:
                target = groups[target_index]
                allowed = tiny | {p for p in target if bounded(tiny | {p})}
                retained = set(_distance(graph, min(tiny), allowed))
                if not tiny <= retained or too_small(retained) or not bounded(retained):
                    continue
                displaced = target - retained
                proposed = [retained if i == target_index else set(g)
                            for i, g in enumerate(groups) if i != source_index]
                solution = (proposed, []) if not displaced else reassign(proposed, displaced, min(target))
                if solution is None:
                    continue
                updated, moves = solution
                final_target_index = target_index - int(source_index < target_index)
                audit["repairs"].append({
                    "repair_type": "ABSORPTION",
                    "source_positions": sorted(tiny),
                    "target_positions_before": sorted(target),
                    "target_positions_after": sorted(updated[final_target_index]),
                    "boundary_moves": moves,
                })
                groups = updated
                repaired = True
                break
            if not repaired:
                for _, _, target_index in targets:
                    target = groups[target_index]
                    retained = rebuild_boundary(tiny, target)
                    if retained is None:
                        continue
                    remaining = target - retained
                    audit["repairs"].append({
                        "repair_type": "BOUNDARY_REBUILD",
                        "source_positions": sorted(tiny),
                        "target_positions_before": sorted(target),
                        "target_positions_after": sorted(retained),
                        "rebuilt_remainder_positions": sorted(remaining),
                        "boundary_moves": [
                            {"position": p, "from_region": min(target), "to_region": min(tiny)}
                            for p in sorted(retained - tiny)
                        ],
                    })
                    groups[source_index] = retained
                    groups[target_index] = remaining
                    repaired = True
                    break
            if repaired:
                break
        if not repaired:
            fail("UNSATISFIED_PARTITION_CONSTRAINTS")

    assigned = [p for g in groups for p in g]
    if (len(assigned) != len(set(assigned)) or set(assigned) != set(graph)
            or any(too_small(g) or not bounded(g) or set(_distance(graph, min(g), g)) != g
                   for g in groups)):
        fail("UNSATISFIED_PARTITION_CONSTRAINTS")
    audit.update({"status": "COMPLETE", "realized_region_count": len(groups)})
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
        if len(group) < min_residues:
            flags.append("SMALL_DISPLAY_REGION")
        if len(groups) > max(1, max_sectors):
            flags.append("DISPLAY_COUNT_EXPANDED_FOR_GEOMETRY")
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
        if len(induced_reachable) != len(group) or not bounded(group):
            raise ValueError("Display region violates connectivity or diameter constraints")
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
