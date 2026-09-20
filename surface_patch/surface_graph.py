from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components as scipy_connected_components
from scipy.spatial import cKDTree

from .models import SurfaceGeometry


def _segment_clear(
    left: np.ndarray,
    right: np.ndarray,
    ignored_atoms: set[int],
    atom_xyz: np.ndarray,
    atom_radii: np.ndarray,
    atom_tree: cKDTree,
) -> bool:
    # Interior samples prevent a graph edge from tunnelling through excluded
    # molecular volume. Endpoint owner atoms are ignored because a chord
    # between two points on the same expanded sphere naturally lies inside it.
    for fraction in (0.25, 0.5, 0.75):
        point = left + fraction * (right - left)
        for atom_id in atom_tree.query_ball_point(point, float(atom_radii.max(initial=0.0))):
            if atom_id in ignored_atoms:
                continue
            if np.linalg.norm(point - atom_xyz[atom_id]) < atom_radii[atom_id] - 1e-5:
                return False
    return True


def build_point_cloud_edges(
    vertices: np.ndarray,
    normals: np.ndarray,
    owner_atom: np.ndarray,
    atom_xyz: np.ndarray,
    atom_radii: np.ndarray,
    max_distance: float = 2.5,
    normal_cosine_min: float = 0.3,
) -> tuple[np.ndarray, np.ndarray]:
    if not len(vertices):
        return np.empty((0, 2), dtype=np.int64), np.empty(0, dtype=float)
    point_tree = cKDTree(vertices)
    atom_tree = cKDTree(atom_xyz)
    output: list[tuple[int, int]] = []
    lengths: list[float] = []
    for left, right in sorted(point_tree.query_pairs(max_distance)):
        if float(np.dot(normals[left], normals[right])) < normal_cosine_min:
            continue
        if not _segment_clear(
            vertices[left], vertices[right],
            {int(owner_atom[left]), int(owner_atom[right])},
            atom_xyz, atom_radii, atom_tree,
        ):
            continue
        output.append((left, right))
        lengths.append(float(np.linalg.norm(vertices[left] - vertices[right])))
    return np.asarray(output, dtype=np.int64).reshape((-1, 2)), np.asarray(lengths, dtype=float)


def connected_components(vertex_count: int, edges: np.ndarray) -> np.ndarray:
    if vertex_count == 0:
        return np.empty(0, dtype=np.int64)
    graph = graph_csr(vertex_count, edges, np.ones(len(edges), dtype=float))
    _, labels = scipy_connected_components(graph, directed=False, return_labels=True)
    return labels.astype(np.int64)


def graph_csr(vertex_count: int, edges: np.ndarray, edge_length: np.ndarray) -> csr_matrix:
    if not len(edges):
        return csr_matrix((vertex_count, vertex_count), dtype=float)
    rows = np.concatenate((edges[:, 0], edges[:, 1]))
    cols = np.concatenate((edges[:, 1], edges[:, 0]))
    values = np.concatenate((edge_length, edge_length))
    return coo_matrix((values, (rows, cols)), shape=(vertex_count, vertex_count)).tocsr()


def build_surface_graph(surface: SurfaceGeometry) -> csr_matrix:
    return graph_csr(len(surface.vertices), surface.edges, surface.edge_length)
