from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix

from .geodesic import geodesic_diameter
from .models import Footprint, SurfaceGeometry


def evaluate_footprint_quality(
    footprints: list[Footprint],
    surface: SurfaceGeometry,
    graph: csr_matrix,
    min_residue_count: int = 8,
    max_residue_count: int = 30,
    geodesic_span_warning: float = 32.0,
    low_confidence_threshold: float = 70.0,
) -> list[Footprint]:
    for footprint in footprints:
        flags: list[str] = []
        component_ids = {int(surface.component_id[item]) for item in footprint.vertex_ids}
        if len(component_ids) != 1:
            flags.append("ERROR_DISCONNECTED_SURFACE")
        footprint.geodesic_diameter = round(geodesic_diameter(graph, footprint.vertex_ids), 3)
        if footprint.geodesic_diameter > geodesic_span_warning:
            flags.append("WARNING_GEODESIC_SPAN")
        if len(footprint.residues) < min_residue_count:
            flags.append("WARNING_FEW_RESIDUES")
        elif len(footprint.residues) > max_residue_count:
            flags.append("WARNING_MANY_RESIDUES")
        vertex_ids = np.asarray(footprint.vertex_ids, dtype=np.int64)
        confidence = surface.confidence[vertex_ids]
        known = np.isfinite(confidence)
        if known.any():
            weights = surface.vertex_area[vertex_ids]
            low_area = float(weights[known & (confidence < low_confidence_threshold)].sum())
            known_area = float(weights[known].sum())
            footprint.low_confidence_surface_fraction = low_area / known_area if known_area else 0.0
            if footprint.low_confidence_surface_fraction > 0.25:
                flags.append("WARNING_LOW_CONFIDENCE")
        if not footprint.mutable_residues:
            flags.append("WARNING_NO_MUTABLE_RESIDUES")
        footprint.quality_flags = flags
    return footprints


def global_surface_coverage(
    footprints: list[Footprint], surface: SurfaceGeometry, eligible_vertices: np.ndarray
) -> float:
    if not len(eligible_vertices):
        return 0.0
    covered = set()
    for footprint in footprints:
        if footprint.duplicate_of is None:
            covered.update(footprint.vertex_ids)
    eligible = set(int(item) for item in eligible_vertices)
    numerator = float(surface.vertex_area[list(eligible & covered)].sum()) if eligible & covered else 0.0
    denominator = float(surface.vertex_area[list(eligible)].sum())
    return numerator / denominator if denominator else 0.0
