from __future__ import annotations

import heapq
from collections import deque
from typing import Any

import numpy as np

from ..models import ResidueRecord
from .models import SurfaceResidueUnit
from .structure_audit import residue_id


def _component_count(graph: dict[int, dict[int, float]]) -> int:
    remaining = set(graph)
    count = 0
    while remaining:
        count += 1
        queue = deque([remaining.pop()])
        while queue:
            current = queue.popleft()
            for neighbor in graph[current]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    queue.append(neighbor)
    return count


def _contract_surface_transit_nodes(
    full_graph: dict[int, dict[int, float]],
    eligible_positions: set[int],
    maximum_path_a: float,
) -> dict[int, dict[int, float]]:
    """Collapse non-design surface residues into bounded geodesic edges.

    Accessibility thresholds define which residues can be sector members; they
    must not punch holes in the physical surface used to measure connectivity.
    Search stops at the first eligible residue on each branch so the resulting
    edge remains a local surface step rather than a shortcut across a domain.
    """
    output = {position: {} for position in sorted(eligible_positions)}
    for start in sorted(eligible_positions):
        distances = {start: 0.0}
        queue = [(0.0, start)]
        while queue:
            distance, current = heapq.heappop(queue)
            if distance != distances.get(current) or distance > maximum_path_a:
                continue
            if current != start and current in eligible_positions:
                output[start][current] = min(distance, output[start].get(current, float("inf")))
                output[current][start] = min(distance, output[current].get(start, float("inf")))
                continue
            for neighbor, step in full_graph.get(current, {}).items():
                proposed = distance + step
                if proposed <= maximum_path_a and proposed < distances.get(neighbor, float("inf")):
                    distances[neighbor] = proposed
                    heapq.heappush(queue, (proposed, neighbor))
    return output


def build_surface_units(
    records: list[ResidueRecord],
    rsasa_threshold: float = 0.25,
    neighbor_cutoff_a: float = 10.0,
    max_neighbors: int = 8,
    surface_geometry: Any | None = None,
    graph_provenance: dict | None = None,
) -> tuple[list[SurfaceResidueUnit], dict[int, dict[int, float]], dict[int, float]]:
    """Build residue-level surface units and a CA surface graph.

    Prefer the explicit solvent-surface point graph produced by Surface V2.
    The record-only fallback applies normal/clearance guards and is labelled as
    a proxy; production pipeline runs provide the explicit point graph.
    """
    threshold_eligible = [
        item for item in records
        if item.ca is not None and item.coordinate_eligible and item.surface_accessible and item.rsasa >= rsasa_threshold
    ]
    geometry_owner_positions = (
        {int(owner) for owner in surface_geometry.owner_residue}
        if surface_geometry is not None
        else None
    )
    # A reliability-trimmed SurfaceGeometry is authoritative. Residues omitted
    # from that geometry must not be reintroduced as isolated nodes merely
    # because their SASA was calculated before trimming.
    eligible = [
        item
        for item in threshold_eligible
        if geometry_owner_positions is None or item.uniprot_pos in geometry_owner_positions
    ]
    eligible.sort(key=lambda item: (item.chain_id, item.auth_seq_id, item.insertion_code, item.uniprot_pos))
    eligible_positions = {item.uniprot_pos for item in eligible}
    areas = {
        item.uniprot_pos: max(float(item.sasa), float(item.sidechain_sasa), 1.0)
        for item in eligible
    }
    graph: dict[int, dict[int, float]] = {item.uniprot_pos: {} for item in eligible}
    ca_by_position = {
        item.uniprot_pos: np.asarray(item.ca, dtype=float) for item in records if item.ca is not None
    }
    if surface_geometry is not None:
        surface_areas: dict[int, float] = {}
        for owner, area in zip(surface_geometry.owner_residue, surface_geometry.vertex_area):
            position = int(owner)
            if position in eligible_positions:
                surface_areas[position] = surface_areas.get(position, 0.0) + float(area)
        areas.update({position: max(value, 1.0) for position, value in surface_areas.items()})
        surface_owners = {
            int(owner) for owner in surface_geometry.owner_residue if int(owner) in ca_by_position
        }
        full_graph: dict[int, dict[int, float]] = {position: {} for position in surface_owners}
        for edge, _length in zip(surface_geometry.edges, surface_geometry.edge_length):
            left = int(surface_geometry.owner_residue[int(edge[0])])
            right = int(surface_geometry.owner_residue[int(edge[1])])
            if left == right or left not in full_graph or right not in full_graph:
                continue
            # Surface-point adjacency proves the edge does not tunnel through
            # excluded molecular volume; CA distance provides a residue-scale
            # path length without collapsing within-residue travel to ~0 A.
            value = float(np.linalg.norm(ca_by_position[left] - ca_by_position[right]))
            # A short point-cloud edge can join two owner residues that are
            # far apart in the folded chain (for example across a near-contact
            # between separate surface segments).  Do not collapse such a
            # bridge into one residue-level surface graph edge.
            if value > neighbor_cutoff_a:
                continue
            full_graph[left][right] = min(value, full_graph[left].get(right, float("inf")))
            full_graph[right][left] = min(value, full_graph[right].get(left, float("inf")))

        # Sparse solvent point sampling can miss an owner transition at the
        # peptide bond. Consecutive observed residues are a safe local bridge
        # when both contribute solvent-surface points and their C-alpha atoms
        # retain normal peptide geometry.
        sequence_bridges = 0
        ordered = sorted(
            (item for item in records if item.uniprot_pos in surface_owners),
            key=lambda item: (item.chain_id, item.auth_seq_id, item.insertion_code),
        )
        for left_record, right_record in zip(ordered, ordered[1:]):
            if left_record.chain_id != right_record.chain_id:
                continue
            if right_record.uniprot_pos - left_record.uniprot_pos != 1:
                continue
            left, right = left_record.uniprot_pos, right_record.uniprot_pos
            value = float(np.linalg.norm(ca_by_position[left] - ca_by_position[right]))
            if value > 4.5 or value <= 1e-8:
                continue
            if right not in full_graph[left]:
                sequence_bridges += 1
            full_graph[left][right] = min(value, full_graph[left].get(right, float("inf")))
            full_graph[right][left] = min(value, full_graph[right].get(left, float("inf")))

        direct_eligible_edges = sum(
            1
            for left in eligible_positions
            for right in full_graph.get(left, {})
            if right in eligible_positions and left < right
        )
        graph = _contract_surface_transit_nodes(
            full_graph,
            eligible_positions,
            maximum_path_a=max(4.5, neighbor_cutoff_a * 1.5),
        )
        contracted_edges = sum(len(neighbors) for neighbors in graph.values()) // 2
        if graph_provenance is not None:
            graph_provenance.update({
                "backend": str(surface_geometry.backend),
                "basis": "explicit_solvent_surface_point_graph",
                "surface_vertex_count": int(len(surface_geometry.vertices)),
                "surface_edge_count": int(len(surface_geometry.edges)),
                "distance_unit": "angstrom",
                "through_core_guard": True,
                "eligible_residue_count": len(eligible_positions),
                "reliability_omitted_eligible_count": len(threshold_eligible) - len(eligible),
                "transit_surface_residue_count": len(surface_owners - eligible_positions),
                "sequence_bridge_count": sequence_bridges,
                "direct_eligible_edge_count": direct_eligible_edges,
                "contracted_edge_count": contracted_edges,
                "component_count": _component_count(graph),
                "transit_path_limit_a": max(4.5, neighbor_cutoff_a * 1.5),
            })
    else:
        all_ca = [(item.uniprot_pos, np.asarray(item.ca, dtype=float)) for item in records if item.ca is not None]
        for index, item in enumerate(eligible):
            point = np.asarray(item.ca, dtype=float)
            distances: list[tuple[float, int]] = []
            for other in eligible[index + 1 :]:
                other_point = np.asarray(other.ca, dtype=float)
                delta = other_point - point
                distance = float(np.linalg.norm(delta))
                if distance > neighbor_cutoff_a or distance <= 1e-8:
                    continue
                direction = delta / distance
                left_normal = np.asarray(item.surface_normal, dtype=float) if item.surface_normal is not None else None
                right_normal = np.asarray(other.surface_normal, dtype=float) if other.surface_normal is not None else None
                if left_normal is not None and right_normal is not None:
                    if float(np.dot(left_normal, right_normal)) < 0.0:
                        continue
                    if abs(float(np.dot(left_normal, direction))) > 0.75:
                        continue
                    if abs(float(np.dot(right_normal, direction))) > 0.75:
                        continue
                blocked = False
                for fraction in (0.25, 0.5, 0.75):
                    sample = point + fraction * delta
                    if any(
                        position not in {item.uniprot_pos, other.uniprot_pos}
                        and float(np.linalg.norm(sample - coordinate)) < 2.8
                        for position, coordinate in all_ca
                    ):
                        blocked = True
                        break
                if not blocked:
                    distances.append((distance, other.uniprot_pos))
            for distance, position in sorted(distances)[:max_neighbors]:
                graph[item.uniprot_pos][position] = distance
                graph[position][item.uniprot_pos] = distance
        if graph_provenance is not None:
            graph_provenance.update({
                "backend": "residue_visibility_fallback",
                "basis": "ca_distance_with_normal_and_clearance_guards",
                "distance_unit": "angstrom",
                "through_core_guard": True,
                "warning": "No explicit solvent-surface geometry was supplied.",
            })
    components: dict[int, int] = {}
    component_id = 0
    for start in sorted(graph):
        if start in components:
            continue
        component_id += 1
        queue = deque([start])
        components[start] = component_id
        while queue:
            current = queue.popleft()
            for neighbor in sorted(graph[current]):
                if neighbor not in components:
                    components[neighbor] = component_id
                    queue.append(neighbor)
    units: list[SurfaceResidueUnit] = []
    for item in eligible:
        reasons = list(item.hard_exclusions)
        if item.disulfide or item.aa == "C":
            reasons.append("disulfide-cysteine" if item.disulfide else "cysteine-default-protection")
        if item.glycosylation_site or item.glycosylation_sequon:
            reasons.append("glycosylation-sequon")
        if item.interface:
            reasons.append("native-interface")
        if item.domain_interface:
            reasons.append("domain-interface")
        if item.flexible_terminal:
            reasons.append("flexible-terminal")
        if item.salt_bridge:
            reasons.append("salt-bridge-core")
        if item.aa in {"G", "P"}:
            reasons.append("backbone-conformation-sensitive")
        designable = not reasons
        units.append(
            SurfaceResidueUnit(
                residue=residue_id(item), position=item.uniprot_pos,
                rasa=float(item.rsasa), mapped_area_a2=areas[item.uniprot_pos],
                designable=designable, exclusion_reasons=sorted(set(reasons)),
                component_id=components[item.uniprot_pos],
                accessibility_confidence=min(1.0, max(0.0, float(item.rsasa))),
                coordinate_xyz=(
                    tuple(float(value) for value in item.ca)
                    if item.ca is not None else None
                ),
            )
        )
    return units, graph, areas


def unit_records(units: list[SurfaceResidueUnit], records: list[ResidueRecord]) -> dict[int, ResidueRecord]:
    wanted = {unit.position for unit in units}
    return {record.uniprot_pos: record for record in records if record.uniprot_pos in wanted}
