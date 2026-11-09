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


def _build_legacy_sectors(
    surface: list[SurfaceEvidence],
    epitope: list[EpitopeEvidence],
    graph: dict[int, dict[int, float]],
    records: dict[int, Any],
    max_sectors: int = 7,
    parameters: dict[str, Any] | None = None,
    diagnostics: dict[str, Any] | None = None,
    _budget: dict[str, int] | None = None,
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

    def spend() -> None:
        states = _budget['states'] if _budget is not None else audit['search_states']
        if states >= PARTITION_SEARCH_LIMIT:
            fail('SEARCH_BUDGET_EXHAUSTED')
        audit['search_states'] += 1
        if _budget is not None:
            _budget['states'] += 1

    def reassign(
        proposed: list[set[int]], pending: set[int], source_region: int,
    ) -> tuple[list[set[int]], list[dict[str, int]]] | None:
        # ponytail: bounded local DFS, not a global partition solver; fail explicitly
        # at the ceiling, extend the neighborhood only if real inputs require it.
        stack = [(proposed, pending, [])]
        visited = set()
        while stack:
            current, remaining, moves = stack.pop()
            identity = (tuple(sorted(tuple(sorted(g)) for g in current)), tuple(sorted(remaining)))
            if identity in visited:
                continue
            visited.add(identity)
            spend()
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
            spend()
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
    return _sector_outputs(groups, by_position, epi, graph, records, distances,
                           max_sectors, min_residues, min_area, max_geodesic, max_euclidean)


def _sector_outputs(groups, by_position, epi, graph, records, distances,
                    max_sectors, min_residues, min_area, max_geodesic, max_euclidean):
    components = _components(graph)
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
        if (len(induced_reachable) != len(group) or geodesic_diameter > max_geodesic
                or euclidean_diameter > max_euclidean):
            raise ValueError("Display region violates connectivity or diameter constraints")
        center = min(
            group,
            key=lambda position: (
                max(group_distances[position].get(other, float("inf")) for other in group),
                sum(group_distances[position].get(other, float("inf")) for other in group),
                -float(epi[position].surface_percentile or 0.0) if position in epi else 0.0,
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


def build_epitope_sectors(
    surface: list[SurfaceEvidence],
    epitope: list[EpitopeEvidence],
    graph: dict[int, dict[int, float]],
    records: dict[int, Any],
    max_sectors: int = 7,
    parameters: dict[str, Any] | None = None,
    diagnostics: dict[str, Any] | None = None,
) -> list[V4Sector]:
    params = dict(parameters or {})
    policy = params.get('policy', 'balanced_connected_v4')
    if policy not in {'balanced_connected_v4', 'connected_bounded_v3'}:
        raise ValueError(f'Unknown display partition policy: {policy}')
    legacy = policy == 'connected_bounded_v3'
    maximum_geo = float(params.setdefault('maximum_geodesic_diameter_a', 26.0 if legacy else 30.0))
    maximum_ca = float(params.setdefault('maximum_euclidean_diameter_a', 18.0 if legacy else 22.0))
    minimum_area = float(params.get('minimum_sector_area_a2', 50.0))
    minimum_size = int(params.get('minimum_sector_residues', 4))
    target = int(params.get('target_sector_residues', 12))
    preferred_min = int(params.get('preferred_minimum_sector_residues', 10))
    preferred_max = int(params.get('preferred_maximum_sector_residues', 15))
    preferred_count_min = params.get('preferred_minimum_sector_count')
    preferred_count_max = params.get('preferred_maximum_sector_count')
    has_count_range = preferred_count_min is not None or preferred_count_max is not None
    if has_count_range:
        if (isinstance(preferred_count_min, bool) or not isinstance(preferred_count_min, int)
                or isinstance(preferred_count_max, bool) or not isinstance(preferred_count_max, int)):
            raise ValueError('Preferred sector-count bounds must be positive integers') from None
    if (not all(math.isfinite(v) and v > 0 for v in (maximum_geo, maximum_ca, minimum_area))
            or not 1 <= minimum_size <= preferred_min <= target <= preferred_max
            or (has_count_range and not 1 <= preferred_count_min <= preferred_count_max)):
        raise ValueError('Display size preferences and positive finite diameter/area limits are required')
    audit = diagnostics if diagnostics is not None else {}
    audit.clear()
    audit.update({'version': policy, 'status': 'IN_PROGRESS',
        'requested_region_count': max(1, max_sectors), 'count_is_soft_target': True,
        'minimum_sector_residues': minimum_size, 'minimum_sector_area_a2': minimum_area,
        'maximum_geodesic_diameter_a': maximum_geo, 'maximum_euclidean_diameter_a': maximum_ca,
        'target_sector_residues': target, 'preferred_minimum_sector_residues': preferred_min,
        'preferred_maximum_sector_residues': preferred_max, 'size_is_soft_target': True,
        'diameter_basis': 'full_eligible_residue_graph_shortest_paths',
        'connectivity_basis': 'induced_eligible_residue_graph_with_contracted_surface_transit',
        'search_limit': PARTITION_SEARCH_LIMIT, 'search_states': 0, 'repairs': []})
    if has_count_range:
        audit.update(preferred_minimum_sector_count=preferred_count_min,
            preferred_maximum_sector_count=preferred_count_max, count_range_is_soft=True)

    by_position = {row.residue.uniprot_pos: row for row in surface}
    invalid = []
    if len(by_position) != len(surface):
        invalid.append('DUPLICATE_SURFACE_EVIDENCE')
    coordinates = {}
    for p in sorted(graph):
        record = records.get(p)
        try:
            ca = np.asarray(getattr(record, 'ca', None), dtype=float)
        except (ValueError, TypeError):
            ca = np.array([])
        if ca.shape != (3,) or not np.isfinite(ca).all():
            invalid.append(f'INVALID_CA_COORDINATE:{p}')
        else:
            coordinates[p] = ca
        row = by_position.get(p)
        if row is None or not math.isfinite(row.mapped_surface_area_a2) or row.mapped_surface_area_a2 < 0:
            invalid.append(f'INVALID_SURFACE_AREA:{p}')
        for q, value in graph[p].items():
            if (q not in graph or not math.isfinite(value) or value <= 0
                    or graph.get(q, {}).get(p) != value):
                invalid.append(f'INVALID_SURFACE_GRAPH_EDGE:{p}:{q}')
    if invalid:
        audit.update(status='FAILED', reason='INVALID_PARTITION_INPUT', input_errors=invalid)
        raise PartitionConstraintError('INVALID_PARTITION_INPUT: 显示分区缺少有效坐标、面积或表面图证据', audit)
    if legacy:
        result = _build_legacy_sectors(surface, epitope, graph, records, max_sectors, params, audit)
    else:
        components = _components(graph)
        distances = {p: _distance(graph, p, component) for component in components for p in component}
        compatible = {p: {q for q in graph if distances[p].get(q, math.inf) <= maximum_geo
                           and np.linalg.norm(coordinates[p] - coordinates[q]) <= maximum_ca} for p in graph}
        budget = {'states': 0}
        visited = set()
        valid_cache = {}

        class BudgetExhausted(RuntimeError):
            pass

        def identity(groups):
            return tuple(sorted(tuple(sorted(g)) for g in groups))

        def spend(state):
            if state in visited:
                return False
            if budget['states'] >= PARTITION_SEARCH_LIMIT:
                raise BudgetExhausted
            budget['states'] += 1
            visited.add(state)
            return True

        def valid(group):
            key = frozenset(group)
            if key not in valid_cache:
                valid_cache[key] = (len(group) >= minimum_size
                    and sum(by_position[p].mapped_surface_area_a2 for p in sorted(group)) >= minimum_area
                    and all(group <= compatible[p] for p in group)
                    and set(_distance(graph, min(group), group)) == group)
            return valid_cache[key]

        def quality(groups):
            sizes = [len(g) for g in groups]
            size_penalty = sum(max(preferred_min - n, 0) ** 2 + max(n - preferred_max, 0) ** 2
                               for n in sizes)
            if has_count_range:
                count_penalty = (max(preferred_count_min - len(groups), 0) ** 2
                                 + max(len(groups) - preferred_count_max, 0) ** 2)
                return (size_penalty, count_penalty, len(groups), sum((n - target) ** 2 for n in sizes))
            return (size_penalty, len(groups), sum((n - target) ** 2 for n in sizes))

        def adjacent(left, right):
            return any(q in right for p in left for q in graph[p])

        initial_count = max(len(components), math.ceil(len(graph) / target)) if graph else 1
        starts = list(dict.fromkeys(max(len(components), min(len(graph) or 1, count))
            for count in (initial_count, initial_count - 1, initial_count + 1, max(1, max_sectors))))
        audit.update(preferred_region_count=initial_count, initializations=[],
            search_scope='deterministic_multistart_and_local_boundary_balance; not a global optimality proof')
        candidates = []
        incumbent = None
        exhausted = False
        for count in starts:
            attempt = {}
            try:
                candidate = _build_legacy_sectors(surface, epitope, graph, records, count, params, attempt, budget)
                groups = [set(row.residue_positions) for row in candidate]
                attempt['initial_sizes'] = [len(g) for g in groups]
                candidates.append((groups, attempt['repairs'], count))
                option = (quality(groups), identity(groups), groups, attempt['repairs'], count)
                if incumbent is None or option[:2] < incumbent[:2]:
                    incumbent = option
            except PartitionConstraintError:
                if attempt.get('reason') == 'SEARCH_BUDGET_EXHAUSTED':
                    exhausted = True
            audit['initializations'].append(attempt)
            if exhausted:
                break

        def optimize(groups, history, count):
            nonlocal incumbent
            while True:
                current_quality = quality(groups)
                choice = None

                def consider(updated, operation):
                    nonlocal choice, incumbent
                    signature = identity(updated)
                    if not spend(('partition', signature)) or quality(updated) >= current_quality:
                        return
                    if not all(valid(g) for g in updated):
                        return
                    option = (quality(updated), signature, updated, history + [operation], count)
                    if choice is None or option[:2] < choice[:2]:
                        choice = option
                    if incumbent is None or option[:2] < incumbent[:2]:
                        incumbent = option

                pairs = sorted(((i, j) for i, a in enumerate(groups) for j, b in enumerate(groups)
                    if i < j and adjacent(a, b)), key=lambda pair: (
                        min(len(groups[pair[0]]), len(groups[pair[1]])),
                        min(groups[pair[0]]), min(groups[pair[1]])))
                for i, j in pairs:
                    merged = [set(g) for k, g in enumerate(groups) if k not in {i, j}] + [groups[i] | groups[j]]
                    consider(merged, {'repair_type': 'BALANCED_MERGE',
                        'source_positions': sorted(groups[i]), 'target_positions_before': sorted(groups[j]),
                        'target_positions_after': sorted(groups[i] | groups[j]), 'boundary_moves': []})
                    for donor, recipient in ((i, j), (j, i)):
                        boundary = sorted(p for p in groups[donor] if adjacent({p}, groups[recipient]))
                        for p in boundary:
                            updated = [set(g) for g in groups]
                            updated[donor].remove(p)
                            updated[recipient].add(p)
                            consider(updated, {'repair_type': 'BALANCED_TRANSFER',
                                'boundary_moves': [{'position': p, 'from_region': min(groups[donor]),
                                                    'to_region': min(groups[recipient])}]})
                    left_boundary = sorted(p for p in groups[i] if adjacent({p}, groups[j]))
                    right_boundary = sorted(p for p in groups[j] if adjacent({p}, groups[i]))
                    for p in left_boundary:
                        for q in right_boundary:
                            updated = [set(g) for g in groups]
                            updated[i] = groups[i] - {p} | {q}
                            updated[j] = groups[j] - {q} | {p}
                            consider(updated, {'repair_type': 'BALANCED_SWAP', 'boundary_moves': [
                                {'position': p, 'from_region': min(groups[i]), 'to_region': min(groups[j])},
                                {'position': q, 'from_region': min(groups[j]), 'to_region': min(groups[i])}]})
                if choice is None:
                    for i, j in pairs:
                        recipient, donor = min((i, j), key=lambda k: (len(groups[k]), min(groups[k]))), max(
                            (i, j), key=lambda k: (len(groups[k]), min(groups[k])))
                        if recipient == donor or len(groups[recipient]) >= preferred_min:
                            continue
                        union = groups[i] | groups[j]
                        stack = [set(groups[recipient])]
                        pair_seen = set()
                        while stack:
                            retained = stack.pop()
                            if frozenset(retained) in pair_seen:
                                continue
                            pair_seen.add(frozenset(retained))
                            spend(('rebuild', identity(groups), i, j, tuple(sorted(retained))))
                            remaining = union - retained
                            if valid(retained) and valid(remaining):
                                updated = [set(g) for g in groups]
                                updated[recipient], updated[donor] = retained, remaining
                                consider(updated, {'repair_type': 'BALANCED_BOUNDARY_REBUILD',
                                    'boundary_moves': [{'position': p, 'from_region': min(groups[donor]),
                                                        'to_region': min(groups[recipient])}
                                                       for p in sorted(retained - groups[recipient])]})
                            if len(retained) >= preferred_max or len(remaining) <= minimum_size:
                                continue
                            frontier = {q for p in retained for q in graph[p] if q in remaining}
                            choices = sorted((min(distances[p].get(q, math.inf) for p in groups[recipient]), q)
                                for q in frontier if retained | {q} <= compatible[q]
                                and all(retained | {q} <= compatible[p] for p in retained))
                            stack.extend(retained | {q} for _, q in reversed(choices))
                if choice is None:
                    return
                _, _, groups, history, _ = choice

        if incumbent is not None and not exhausted:
            try:
                for groups, history, count in sorted(candidates, key=lambda row: (quality(row[0]), identity(row[0]))):
                    optimize(groups, history, count)
            except BudgetExhausted:
                exhausted = True
        if incumbent is None:
            reason = 'SEARCH_BUDGET_EXHAUSTED' if exhausted else 'UNSATISFIED_PARTITION_CONSTRAINTS'
            drafts = [attempt['draft_regions'] for attempt in audit['initializations'] if attempt.get('draft_regions')]
            draft = min(drafts, key=lambda rows: (
                sum(row['undersized'] for row in rows), len(rows),
                tuple(tuple(row['residue_positions']) for row in rows))) if drafts else []
            violations = sorted({rule for row in draft for rule, failed in [
                ('MINIMUM_RESIDUES_NOT_MET', row['residue_count'] < minimum_size),
                ('MINIMUM_SURFACE_AREA_NOT_MET', row['total_area_a2'] < minimum_area)] if failed})
            audit.update(status='FAILED', reason=reason, search_states=budget['states'],
                         draft_regions=draft, violations=violations)
            raise PartitionConstraintError(f'{reason}: 未找到满足硬约束的显示分区；不能据此证明全局无解', audit)
        objective, _, groups, repairs, chosen_count = incumbent
        assigned = [p for g in groups for p in g]
        if len(assigned) != len(set(assigned)) or set(assigned) != set(graph) or not all(valid(g) for g in groups):
            audit.update(status='FAILED', reason='UNSATISFIED_PARTITION_CONSTRAINTS')
            raise PartitionConstraintError('UNSATISFIED_PARTITION_CONSTRAINTS: 最终显示分区未通过硬验收', audit)
        audit.update(status='COMPLETE', realized_region_count=len(groups), search_states=budget['states'],
            chosen_initial_region_count=chosen_count, repairs=repairs, objective=list(objective),
            objective_order=(['squared_preferred_size_range_deviation', 'squared_preferred_count_range_deviation',
                              'region_count', 'squared_target_size_deviation'] if has_count_range else
                             ['squared_preferred_range_deviation', 'region_count', 'squared_target_size_deviation']),
            optimization_status='SEARCH_BUDGET_EXHAUSTED' if exhausted else 'LOCAL_SEARCH_COMPLETE')
        result = _sector_outputs(groups, by_position, {r.residue.uniprot_pos: r for r in epitope},
            graph, records, distances, max_sectors, minimum_size, minimum_area, maximum_geo, maximum_ca)

    distances = {p: _distance(graph, p, component) for component in _components(graph) for p in component}
    audit['region_metrics'] = []
    audit['size_exceptions'] = []
    components_by_position = {p: component for component in _components(graph) for p in component}
    for row in result:
        group = set(row.residue_positions)
        induced = {p: _distance(graph, p, group) for p in group}
        audit['region_metrics'].append({'sector_id': row.sector_id, 'residue_count': len(group),
            'full_graph_geodesic_diameter_a': max(distances[p][q] for p in group for q in group),
            'induced_graph_geodesic_diameter_a': max(induced[p][q] for p in group for q in group),
            'ca_diameter_a': max(float(np.linalg.norm(coordinates[p] - coordinates[q])) for p in group for q in group)})
        if not legacy and not preferred_min <= len(group) <= preferred_max:
            reason = ('SEARCH_BUDGET_EXHAUSTED' if audit['optimization_status'] == 'SEARCH_BUDGET_EXHAUSTED'
                else 'COMPONENT_SIZE_BELOW_PREFERRED_RANGE' if len(components_by_position[min(group)]) < preferred_min
                else 'LOCAL_BALANCE_LIMIT')
            description = ('平衡搜索预算耗尽，保留已满足硬约束的方案；未证明无法进一步改善。'
                if reason == 'SEARCH_BUDGET_EXHAUSTED' else '该独立表面连通分量小于偏好大小，不能跨越不连通表面归并。'
                if reason == 'COMPONENT_SIZE_BELOW_PREFERRED_RANGE' else '有界局部平衡未找到更优方案，仍偏离大小偏好；不代表全局无解。')
            audit['size_exceptions'].append({'sector_id': row.sector_id, 'residue_count': len(group),
                                            'reason': reason, 'description': description})
            row.quality_flags.append('BELOW_PREFERRED_DISPLAY_SIZE' if len(group) < preferred_min
                                     else 'ABOVE_PREFERRED_DISPLAY_SIZE')
    if has_count_range:
        realized_count = len(result)
        audit['count_preference_met'] = preferred_count_min <= realized_count <= preferred_count_max
        if not audit['count_preference_met']:
            audit['count_range_exception'] = {
                'realized_region_count': realized_count,
                'reason': 'SOFT_COUNT_RANGE_EXCEEDED',
                'description': '连通性、覆盖或硬几何约束优先于分区数量偏好。',
            }
    return result
