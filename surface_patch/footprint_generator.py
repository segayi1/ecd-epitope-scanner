from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix

from .geodesic import distances_from
from .models import Footprint, SurfaceGeometry


def _raw_footprint(
    surface: SurfaceGeometry,
    center: int,
    vertex_ids: np.ndarray,
    mode: str,
    label: str,
    distances: np.ndarray,
    radius: float | None = None,
    target_area: float | None = None,
) -> Footprint:
    ids = [int(item) for item in vertex_ids]
    points = surface.vertices[vertex_ids]
    weights = surface.vertex_area[vertex_ids]
    centroid = np.average(points, axis=0, weights=weights) if len(points) else surface.vertices[center]
    normal = np.average(surface.normals[vertex_ids], axis=0, weights=weights) if len(points) else surface.normals[center]
    norm = float(np.linalg.norm(normal))
    unit_normal = normal / norm if norm else np.zeros(3)
    dispersion = 1.0 - float(np.average(surface.normals[vertex_ids] @ unit_normal, weights=weights)) if len(points) else 0.0
    finite_distances = distances[vertex_ids]
    return Footprint(
        footprint_id=label,
        center_vertex_id=center,
        center_residue_id=int(surface.owner_residue[center]),
        generation_mode=mode,
        radius=radius,
        target_area=target_area,
        vertex_ids=ids,
        total_surface_area=float(weights.sum()),
        geodesic_radius=float(finite_distances.max(initial=0.0)),
        centroid_xyz=tuple(float(item) for item in centroid),
        mean_normal=tuple(float(item) for item in unit_normal),
        normal_dispersion=round(dispersion, 6),
        surface_component_id=int(surface.component_id[center]),
    )


def generate_footprints(
    surface: SurfaceGeometry,
    graph: csr_matrix,
    centers: list[int],
    radii: tuple[float, ...] = (9.0, 12.0, 15.0),
    adaptive_target_area: float = 600.0,
    adaptive_min_area: float = 350.0,
    adaptive_max_area: float = 1000.0,
    adaptive_min_radius: float = 8.0,
    adaptive_max_radius: float = 16.0,
) -> list[Footprint]:
    output: list[Footprint] = []
    for center_number, center in enumerate(centers, start=1):
        distances = distances_from(graph, center)
        component = surface.component_id == surface.component_id[center]
        for radius in radii:
            members = np.flatnonzero(component & (distances <= radius + 1e-9))
            output.append(_raw_footprint(
                surface, center, members, "fixed_radius",
                f"C{center_number:03d}_R{int(round(radius)):02d}", distances, radius=radius,
            ))

        candidates = np.flatnonzero(component & np.isfinite(distances) & (distances <= adaptive_max_radius))
        ordered = sorted(candidates, key=lambda item: (float(distances[item]), int(item)))
        selected: list[int] = []
        area = 0.0
        for vertex_id in ordered:
            distance = float(distances[vertex_id])
            contribution = float(surface.vertex_area[vertex_id])
            must_include = distance <= adaptive_min_radius or area < adaptive_min_area
            if not must_include and area >= adaptive_target_area:
                break
            if not must_include and area + contribution > adaptive_max_area:
                break
            selected.append(int(vertex_id))
            area += contribution
        output.append(_raw_footprint(
            surface, center, np.asarray(selected, dtype=np.int64), "adaptive_area",
            f"C{center_number:03d}_A{int(round(adaptive_target_area)):04d}", distances,
            target_area=adaptive_target_area,
        ))
    return output
