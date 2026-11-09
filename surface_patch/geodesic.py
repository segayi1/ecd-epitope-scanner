from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra


def distances_from(graph: csr_matrix, sources: int | list[int] | np.ndarray) -> np.ndarray:
    return np.asarray(dijkstra(graph, directed=False, indices=sources), dtype=float)


def geodesic_diameter(graph: csr_matrix, vertices: list[int]) -> float:
    """Deterministic two-sweep estimate of induced-subgraph diameter.

    Exact all-pairs Dijkstra is prohibitive for hundreds of overlapping point
    cloud footprints.  Two sweeps are exact on trees and a stable lower bound
    on meshed surfaces, which is sufficient for the warning-only QC metric.
    """
    if len(vertices) < 2:
        return 0.0
    indices = np.asarray(sorted(set(vertices)), dtype=np.int64)
    subgraph = graph[indices][:, indices]
    first = np.asarray(dijkstra(subgraph, directed=False, indices=0), dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(first))
    if not len(finite_indices):
        return 0.0
    farthest = int(max(finite_indices, key=lambda item: (float(first[item]), -int(item))))
    second = np.asarray(dijkstra(subgraph, directed=False, indices=farthest), dtype=float)
    finite = second[np.isfinite(second)]
    return float(finite.max(initial=0.0))
