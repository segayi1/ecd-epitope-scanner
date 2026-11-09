from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(slots=True)
class SurfaceGeometry:
    vertices: np.ndarray
    normals: np.ndarray
    vertex_area: np.ndarray
    edges: np.ndarray
    edge_length: np.ndarray
    owner_residue: np.ndarray
    owner_atom: np.ndarray
    chain_id: np.ndarray
    canonical_residue_id: np.ndarray
    is_ecd: np.ndarray
    confidence: np.ndarray
    component_id: np.ndarray
    faces: np.ndarray | None = None
    backend: str = "unknown"
    metadata: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        count = len(self.vertices)
        for name in (
            "normals", "vertex_area", "owner_residue", "owner_atom", "chain_id",
            "canonical_residue_id", "is_ecd", "confidence", "component_id",
        ):
            if len(getattr(self, name)) != count:
                raise ValueError(f"SurfaceGeometry.{name} length does not match vertices")
        if self.edges.ndim != 2 or self.edges.shape[1:] != (2,):
            raise ValueError("SurfaceGeometry.edges must have shape (N, 2)")
        if len(self.edge_length) != len(self.edges):
            raise ValueError("SurfaceGeometry.edge_length length does not match edges")


@dataclass(slots=True)
class Footprint:
    footprint_id: str
    center_vertex_id: int
    center_residue_id: int | None
    generation_mode: str
    radius: float | None
    target_area: float | None
    vertex_ids: list[int]
    residues: list[int] = field(default_factory=list)
    residue_surface_area: dict[int, float] = field(default_factory=dict)
    residue_surface_fraction: dict[int, float] = field(default_factory=dict)
    total_surface_area: float = 0.0
    geodesic_radius: float = 0.0
    geodesic_diameter: float = 0.0
    centroid_xyz: tuple[float, float, float] = (0.0, 0.0, 0.0)
    mean_normal: tuple[float, float, float] = (0.0, 0.0, 0.0)
    normal_dispersion: float = 0.0
    mean_rsasa: float = 0.0
    mean_plddt: float | None = None
    minimum_plddt: float | None = None
    low_confidence_surface_fraction: float = 0.0
    structural_domains: list[str] = field(default_factory=list)
    mutable_residues: list[int] = field(default_factory=list)
    protected_residues: list[int] = field(default_factory=list)
    quality_flags: list[str] = field(default_factory=list)
    region_id: str | None = None
    duplicate_of: str | None = None
    surface_component_id: int | None = None
    epitope_prior: float = 0.0


@dataclass(slots=True)
class FootprintEdge:
    footprint_a: str
    footprint_b: str
    weighted_jaccard: float
    shared_surface_area: float
    center_geodesic_distance: float
    center_euclidean_distance: float


@dataclass(slots=True)
class SurfaceRegion:
    region_id: str
    footprint_ids: list[str]
    representative_footprint_id: str
    residues_union: list[int]
    surface_area_union: float
    centroid_xyz: tuple[float, float, float]
    domain_ids: list[str]
    quality_flags: list[str] = field(default_factory=list)
    display_residues: list[int] = field(default_factory=list)
    primary_footprint_ids: list[str] = field(default_factory=list)
    geodesic_extent: float = 0.0
