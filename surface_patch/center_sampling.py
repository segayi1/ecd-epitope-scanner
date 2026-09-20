from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix

from .geodesic import distances_from
from .models import SurfaceGeometry


def eligible_center_vertices(
    surface: SurfaceGeometry,
    exposed_residues: set[int],
    minimum_component_area: float,
    minimum_confidence: float = 70.0,
) -> np.ndarray:
    component_area = {
        int(component): float(surface.vertex_area[surface.component_id == component].sum())
        for component in np.unique(surface.component_id)
    }
    output = []
    for vertex_id, residue_id in enumerate(surface.owner_residue):
        confidence = surface.confidence[vertex_id]
        if int(residue_id) not in exposed_residues:
            continue
        if component_area[int(surface.component_id[vertex_id])] < minimum_component_area:
            continue
        if np.isfinite(confidence) and confidence < minimum_confidence:
            continue
        output.append(vertex_id)
    return np.asarray(output, dtype=np.int64)


def sample_surface_centers(
    surface: SurfaceGeometry,
    graph: csr_matrix,
    eligible_vertices: np.ndarray,
    center_spacing: float = 7.0,
    coverage_radius: float = 7.0,
) -> list[int]:
    """Deterministic geodesic farthest-point sampling, component by component."""
    if not len(eligible_vertices):
        return []
    residue_area: dict[int, float] = {}
    for residue_id in np.unique(surface.owner_residue):
        residue_area[int(residue_id)] = float(
            surface.vertex_area[surface.owner_residue == residue_id].sum()
        )
    centers: list[int] = []
    for component in sorted(set(int(surface.component_id[item]) for item in eligible_vertices)):
        candidates = np.asarray(
            [item for item in eligible_vertices if int(surface.component_id[item]) == component],
            dtype=np.int64,
        )
        first = min(
            candidates,
            key=lambda item: (
                -residue_area[int(surface.owner_residue[item])],
                int(surface.canonical_residue_id[item]),
                int(item),
            ),
        )
        component_centers = [int(first)]
        nearest = distances_from(graph, int(first))[candidates]
        while len(candidates):
            worst_index = max(
                range(len(candidates)),
                key=lambda index: (float(nearest[index]), -int(candidates[index])),
            )
            worst_distance = float(nearest[worst_index])
            if worst_distance <= coverage_radius:
                break
            next_center = int(candidates[worst_index])
            # With disconnected labels this normally cannot happen, but keep
            # spacing as an explicit invariant rather than an assumption.
            if worst_distance + 1e-9 < center_spacing:
                break
            component_centers.append(next_center)
            nearest = np.minimum(nearest, distances_from(graph, next_center)[candidates])
        centers.extend(component_centers)
    return centers
