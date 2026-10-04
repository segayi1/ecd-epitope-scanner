"""Triangle distance and parity occupancy in a fixed coordinate system.

The implicit atom-ball clearance function is deliberately NOT used as SES SDF.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree


def triangle_distance(point, triangles):
    """Exact unsigned point-to-triangle distance, including edges."""
    a, b, c = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    ab, ac = b - a, c - a
    normal = np.cross(ab, ac)
    nn = np.einsum("ij,ij->i", normal, normal)
    displacement = point - a
    height = np.einsum("ij,ij->i", displacement, normal)
    projection = point - (height / np.maximum(nn, 1e-30))[:, None] * normal
    v = projection - a
    d00 = np.einsum("ij,ij->i", ab, ab)
    d01 = np.einsum("ij,ij->i", ab, ac)
    d11 = np.einsum("ij,ij->i", ac, ac)
    d20 = np.einsum("ij,ij->i", v, ab)
    d21 = np.einsum("ij,ij->i", v, ac)
    denominator = np.maximum(d00 * d11 - d01 * d01, 1e-30)
    u = (d11 * d20 - d01 * d21) / denominator
    w = (d00 * d21 - d01 * d20) / denominator
    inside = (u >= -1e-12) & (w >= -1e-12) & (u + w <= 1 + 1e-12) & (nn > 1e-24)
    distance2 = np.where(inside, height * height / np.maximum(nn, 1e-30), np.inf)
    for left, right in ((a, b), (b, c), (c, a)):
        edge = right - left
        t = np.clip(np.einsum("ij,ij->i", point - left, edge) /
                    np.maximum(np.einsum("ij,ij->i", edge, edge), 1e-30), 0, 1)
        delta = point - (left + t[:, None] * edge)
        distance2 = np.minimum(distance2, np.einsum("ij,ij->i", delta, delta))
    return np.sqrt(distance2)


@dataclass
class Mesh:
    vertices: np.ndarray
    faces: np.ndarray

    def __post_init__(self):
        self.vertices = np.asarray(self.vertices, dtype=float)
        self.faces = np.asarray(self.faces, dtype=int)
        if self.vertices.ndim != 2 or self.vertices.shape[1] != 3 or not np.isfinite(self.vertices).all():
            raise ValueError("invalid mesh coordinates")
        if self.faces.ndim != 2 or self.faces.shape[1] != 3 or not len(self.faces):
            raise ValueError("invalid triangles")
        if self.faces.min() < 0 or self.faces.max() >= len(self.vertices):
            raise ValueError("triangle index out of bounds")
        self.triangles = self.vertices[self.faces]
        cross = np.cross(self.triangles[:, 1] - self.triangles[:, 0],
                         self.triangles[:, 2] - self.triangles[:, 0])
        length = np.linalg.norm(cross, axis=1)
        self.area = length / 2
        self.normals = cross / np.maximum(length[:, None], 1e-30)
        self.centers = self.triangles.mean(axis=1)
        self.radii = np.linalg.norm(self.triangles - self.centers[:, None], axis=2).max(axis=1)
        self.max_radius = self.radii.max()
        self.tree = cKDTree(self.centers)
        edges = np.concatenate([self.faces[:, [0, 1]], self.faces[:, [1, 2]], self.faces[:, [2, 0]]])
        unique_edges, inverse, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True,return_inverse=True)
        orientation=np.bincount(inverse,weights=np.where(edges[:,0]<edges[:,1],1.,-1.))
        self.closed = bool(np.all(counts == 2) and np.all(orientation==0) and np.all(self.area > 1e-16))
        self.euler_characteristic = int(len(np.unique(self.faces))-len(unique_edges)+len(self.faces))
        # Irrational-ish ray avoids repeated alignments with mesher axes.
        ray = np.array([1., .371390676, .529150262]); ray /= np.linalg.norm(ray)
        tangent = np.cross(ray, [0., 0., 1.]); tangent /= np.linalg.norm(tangent)
        basis = np.stack([ray, tangent, np.cross(ray, tangent)], axis=1)
        self.ray_basis = basis
        self.ray_triangles = self.triangles @ basis
        self.ray_centers = self.centers @ basis
        self.ray_tree = cKDTree(self.ray_centers[:, 1:])

    def signed_distance(self, points, *, return_projection=False):
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        output = np.full(len(points), np.nan)
        valid = np.zeros(len(points), dtype=bool)
        normals = np.full((len(points),3),np.nan)
        if not self.closed:
            return (output,valid,normals) if return_projection else (output, valid)
        _, nearest = self.tree.query(points)
        rotated = points @ self.ray_basis
        for index, point in enumerate(points):
            first = int(nearest[index])
            upper = float(triangle_distance(point, self.triangles[first:first + 1])[0])
            candidates = self.tree.query_ball_point(point, upper + self.max_radius + 1e-10)
            distances = triangle_distance(point, self.triangles[candidates])
            distance = float(distances.min())
            tied=np.asarray(candidates)[np.abs(distances-distance)<1e-8]
            projected=self.normals[tied]
            if len(projected) and np.all(projected@projected[0]>.5):
                normals[index]=projected[0]
            if distance < 1e-9:
                output[index], valid[index] = 0., True
                continue
            q = rotated[index]
            ids = self.ray_tree.query_ball_point(q[1:], self.max_radius + 1e-10)
            tri = self.ray_triangles[ids]
            a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
            # Barycentric projection in the ray's transverse plane.
            e1, e2, p = b[:, 1:] - a[:, 1:], c[:, 1:] - a[:, 1:], q[1:] - a[:, 1:]
            det = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
            safe = np.abs(det) > 1e-14
            u = np.divide(p[:, 0] * e2[:, 1] - p[:, 1] * e2[:, 0], det,
                          out=np.zeros(len(det)), where=safe)
            v = np.divide(e1[:, 0] * p[:, 1] - e1[:, 1] * p[:, 0], det,
                          out=np.zeros(len(det)), where=safe)
            hit = safe & (u >= 0) & (v >= 0) & (u + v <= 1)
            t = a[:, 0] + u * (b[:, 0] - a[:, 0]) + v * (c[:, 0] - a[:, 0]) - q[0]
            edge_hit = hit & (t > 1e-9) & ((u < 1e-9) | (v < 1e-9) | (1-u-v < 1e-9))
            if edge_hit.any():  # Don't invent an inside/outside sign at an ambiguous ray.
                continue
            inside = int(np.sum(hit & (t > 1e-9))) % 2
            output[index], valid[index] = distance * (-1 if inside else 1), True
        normals[~valid]=np.nan
        return (output,valid,normals) if return_projection else (output, valid)


def align_pair(reference_atoms, wt_atoms, mut_atoms, core_keys):
    """One proper rigid transform for the entire pair; no local mutant fit."""
    keys = sorted(set(core_keys) & reference_atoms.keys() & wt_atoms.keys())
    if len(keys) < 6:
        raise ValueError("REFERENCE_UNRELIABLE: fewer than six core atoms")
    target = np.array([reference_atoms[k] for k in keys])
    moving = np.array([wt_atoms[k] for k in keys])
    if np.linalg.matrix_rank(moving - moving.mean(axis=0)) < 2:
        raise ValueError("REFERENCE_UNRELIABLE: degenerate core")
    u, _, vt = np.linalg.svd((moving-moving.mean(axis=0)).T @ (target-target.mean(axis=0)))
    correction = np.diag([1., 1., np.linalg.det(u @ vt)])
    rotation = u @ correction @ vt
    translation = target.mean(axis=0) - moving.mean(axis=0) @ rotation
    rmsd = float(np.sqrt(np.mean(np.sum((moving @ rotation + translation-target)**2, axis=1))))
    return rotation, translation, rmsd
