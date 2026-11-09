from __future__ import annotations

import math

import numpy as np
from scipy.sparse import csr_matrix

from .footprint_dedup import weighted_jaccard
from .geodesic import distances_from
from .models import Footprint, FootprintEdge


def build_footprint_graph(
    footprints: list[Footprint],
    graph: csr_matrix,
    overlap_threshold: float = 0.20,
    center_distance_threshold: float = 18.0,
) -> list[FootprintEdge]:
    output: list[FootprintEdge] = []
    distance_cache: dict[int, np.ndarray] = {}
    for index, left in enumerate(footprints):
        distance_cache.setdefault(left.center_vertex_id, distances_from(graph, left.center_vertex_id))
        for right in footprints[index + 1 :]:
            jaccard = weighted_jaccard(left, right)
            distance = float(distance_cache[left.center_vertex_id][right.center_vertex_id])
            if jaccard < overlap_threshold and distance > center_distance_threshold:
                continue
            shared = sum(
                min(left.residue_surface_area.get(item, 0.0), right.residue_surface_area.get(item, 0.0))
                for item in set(left.residue_surface_area) & set(right.residue_surface_area)
            )
            euclidean = float(np.linalg.norm(np.asarray(left.centroid_xyz) - np.asarray(right.centroid_xyz)))
            output.append(FootprintEdge(
                footprint_a=left.footprint_id,
                footprint_b=right.footprint_id,
                weighted_jaccard=round(jaccard, 6),
                shared_surface_area=round(shared, 6),
                center_geodesic_distance=round(distance, 6) if math.isfinite(distance) else math.inf,
                center_euclidean_distance=round(euclidean, 6),
            ))
    return output
