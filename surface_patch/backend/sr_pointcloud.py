from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from ...structure_analysis import TIEN_DSSP_RADII, _fibonacci_sphere
from ..models import SurfaceGeometry
from ..surface_graph import build_point_cloud_edges, connected_components


def _radius(atom: Any) -> float:
    element = str(getattr(atom, "element", "") or "").strip().upper()
    return float(TIEN_DSSP_RADII.get(element, 1.70))


class SRPointCloudBackend:
    """Build an explicit solvent-accessible point cloud without changing SASA.

    The production residue SASA remains Bio.PDB's pinned 960-point calculation.
    This independent deterministic sampling exists only to supply geometry.
    """

    def build(
        self,
        structure: Any,
        chain_id: str,
        ecd_residues: dict[int, Any],
        residue_confidence: dict[int, float | None],
        config: dict[str, Any],
        *,
        excluded_residue_objects: set[int] | None = None,
    ) -> SurfaceGeometry:
        probe = float(config.get("probe_radius", 1.4))
        point_count = int(config.get("sampling_points", 96))
        sphere = _fibonacci_sphere(point_count)
        all_atoms = [
            atom for atom in structure.get_atoms()
            if str(getattr(atom, "element", "") or "").strip().upper() != "H"
            and id(atom.get_parent()) not in (excluded_residue_objects or set())
        ]
        atom_xyz = np.asarray([atom.coord for atom in all_atoms], dtype=float).reshape((-1, 3))
        atom_radii = np.asarray([_radius(atom) + probe for atom in all_atoms], dtype=float)
        atom_index = {id(atom): index for index, atom in enumerate(all_atoms)}
        atom_tree = cKDTree(atom_xyz)

        vertices: list[np.ndarray] = []
        normals: list[np.ndarray] = []
        areas: list[float] = []
        owners: list[int] = []
        owner_atoms: list[int] = []
        chains: list[str] = []
        confidences: list[float] = []
        max_radius = float(atom_radii.max(initial=0.0))

        for residue_id, residue in sorted(ecd_residues.items()):
            for atom in residue:
                if str(getattr(atom, "element", "") or "").strip().upper() == "H":
                    continue
                global_atom_id = atom_index.get(id(atom))
                if global_atom_id is None:
                    continue
                center = atom_xyz[global_atom_id]
                radius = atom_radii[global_atom_id]
                local_area = 4.0 * math.pi * radius * radius / point_count
                for normal in sphere:
                    point = center + radius * normal
                    blocked = False
                    for other_id in atom_tree.query_ball_point(point, max_radius):
                        if other_id == global_atom_id:
                            continue
                        if np.linalg.norm(point - atom_xyz[other_id]) < atom_radii[other_id] - 1e-6:
                            blocked = True
                            break
                    if blocked:
                        continue
                    vertices.append(point)
                    normals.append(normal)
                    areas.append(local_area)
                    owners.append(int(residue_id))
                    owner_atoms.append(global_atom_id)
                    chains.append(chain_id)
                    value = residue_confidence.get(residue_id)
                    confidences.append(float(value) if value is not None else np.nan)

        xyz = np.asarray(vertices, dtype=float).reshape((-1, 3))
        normal_array = np.asarray(normals, dtype=float).reshape((-1, 3))
        owner_atom_array = np.asarray(owner_atoms, dtype=np.int64)
        edges, lengths = build_point_cloud_edges(
            xyz,
            normal_array,
            owner_atom_array,
            atom_xyz,
            atom_radii,
            max_distance=float(config.get("point_graph_max_distance", 2.5)),
            normal_cosine_min=float(config.get("normal_cosine_min", 0.3)),
        )
        components = connected_components(len(xyz), edges)
        geometry = SurfaceGeometry(
            vertices=xyz,
            normals=normal_array,
            vertex_area=np.asarray(areas, dtype=float),
            edges=edges,
            edge_length=lengths,
            owner_residue=np.asarray(owners, dtype=np.int64),
            owner_atom=owner_atom_array,
            chain_id=np.asarray(chains, dtype=object),
            canonical_residue_id=np.asarray(owners, dtype=np.int64),
            is_ecd=np.ones(len(xyz), dtype=bool),
            confidence=np.asarray(confidences, dtype=float),
            component_id=components,
            backend="sr_pointcloud",
            metadata={
                "probe_radius": probe,
                "sampling_points_per_atom": point_count,
                "point_graph_max_distance": float(config.get("point_graph_max_distance", 2.5)),
                "normal_cosine_min": float(config.get("normal_cosine_min", 0.3)),
                "sasa_values_unchanged": True,
            },
        )
        geometry.validate()
        return geometry
