from __future__ import annotations

import heapq

from .models import MicroFootprint, ProbeSector, SurfaceResidueUnit


def build_probe_sectors(
    units: list[SurfaceResidueUnit],
    footprints: list[MicroFootprint],
    max_sectors: int,
    graph: dict[int, dict[int, float]] | None = None,
    maximum_geodesic_diameter_a: float = 35.0,
    maximum_euclidean_diameter_a: float = 30.0,
) -> list[ProbeSector]:
    """Merge overlapping footprints into compact, connected surface sectors.

    Geodesic connectivity alone is insufficient after the surface point graph
    is collapsed to residue nodes: a short/incorrect bridge can otherwise make
    two visibly distant patches one sector. Every proposed union must be
    connected, geodesically bounded, component-consistent, and bounded in
    three-dimensional C-alpha span.
    """
    if not footprints:
        return []
    by_residue = {unit.residue: unit for unit in units}
    by_position = {unit.position: unit for unit in units}
    sets = [
        {by_residue[item].position for item in footprint.residue_ids if item in by_residue}
        for footprint in footprints
    ]
    graph = graph or {position: {} for position in by_residue}

    def distances(start: int, allowed: set[int]) -> dict[int, float]:
        output = {start: 0.0}
        queue = [(0.0, start)]
        while queue:
            distance, current = heapq.heappop(queue)
            if distance != output.get(current):
                continue
            for neighbor, step in graph.get(current, {}).items():
                if neighbor not in allowed:
                    continue
                proposed = distance + step
                if proposed < output.get(neighbor, float("inf")):
                    output[neighbor] = proposed
                    heapq.heappush(queue, (proposed, neighbor))
        return output

    def valid_union(residues: set[int]) -> bool:
        if not residues:
            return False
        component_ids = {
            by_position[position].component_id
            for position in residues
            if position in by_position
        }
        if len(component_ids) > 1:
            return False
        start = min(residues)
        paths = distances(start, residues)
        if len(paths) != len(residues):
            return False
        sample = sorted(residues)
        if not all(
            max(distances(position, residues).values(), default=0.0) <= maximum_geodesic_diameter_a
            for position in sample
        ):
            return False
        coordinates = [
            by_position[position].coordinate_xyz
            for position in sample
            if position in by_position and by_position[position].coordinate_xyz is not None
        ]
        if len(coordinates) == len(sample):
            for index, left in enumerate(coordinates):
                for right in coordinates[index + 1 :]:
                    squared = sum((a - b) ** 2 for a, b in zip(left, right))
                    if squared ** 0.5 > maximum_euclidean_diameter_a:
                        return False
        return True

    # Constrained agglomeration avoids transitive overlap chains collapsing an
    # entire protein face into a single sector.
    groups: list[list[int]] = []
    for footprint_index, footprint_set in enumerate(sets):
        best: tuple[float, int] | None = None
        for group_index, group in enumerate(groups):
            group_set = set().union(*(sets[item] for item in group))
            overlap = len(footprint_set & group_set) / max(1, min(len(footprint_set), len(group_set)))
            if overlap < 0.35 or not valid_union(footprint_set | group_set):
                continue
            if best is None or overlap > best[0]:
                best = (overlap, group_index)
        if best is None:
            groups.append([footprint_index])
        else:
            groups[best[1]].append(footprint_index)
    sectors: list[ProbeSector] = []
    for index, group in enumerate(groups, start=1):
        residue_ids = sorted(set().union(*(sets[item] for item in group)))
        total = sum(by_residue[item].mapped_area_a2 for item in residue_ids if item in by_residue)
        designable = sum(by_residue[item].mapped_area_a2 for item in residue_ids if item in by_residue and by_residue[item].designable)
        sectors.append(
            ProbeSector(
                sector_id=f"PS{index:03d}",
                residue_ids=residue_ids,
                footprint_ids=[footprints[item].footprint_id for item in group],
                total_area_a2=round(total, 3),
                designable_area_a2=round(designable, 3),
                compactness=round(sum(len(sets[item]) for item in group) / max(1, len(residue_ids)), 3),
                quality_flags=(
                    ["single_surface_component"] if len(group) == len(footprints) else []
                ) + (["sector_count_exceeds_soft_cap"] if len(groups) > max_sectors else []),
            )
        )
    sectors.sort(key=lambda item: (-item.designable_area_a2, item.sector_id))
    return sectors
