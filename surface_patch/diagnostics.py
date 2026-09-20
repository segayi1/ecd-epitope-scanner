from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from ..util import write_json
from .models import Footprint, FootprintEdge, SurfaceGeometry, SurfaceRegion


def _footprint_dict(item: Footprint) -> dict[str, Any]:
    return asdict(item)


def write_surface_v2_outputs(
    output_dir: Path,
    surface: SurfaceGeometry,
    raw_footprints: list[Footprint],
    footprints: list[Footprint],
    duplicate_footprints: list[Footprint],
    edges: list[FootprintEdge],
    regions: list[SurfaceRegion],
    qc: dict[str, Any],
) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    # CSV is intentionally written without pandas so exports remain available
    # in minimal CLI installations.
    import csv

    with (output_dir / "surface_vertices.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["vertex_id", "x", "y", "z", "nx", "ny", "nz", "area", "owner_residue", "component_id", "confidence"])
        for index, point in enumerate(surface.vertices):
            writer.writerow([
                index, *point.tolist(), *surface.normals[index].tolist(),
                float(surface.vertex_area[index]), int(surface.owner_residue[index]),
                int(surface.component_id[index]),
                "" if not np.isfinite(surface.confidence[index]) else float(surface.confidence[index]),
            ])
    np.savez_compressed(
        output_dir / "surface_edges.npz",
        edges=surface.edges,
        edge_length=surface.edge_length,
    )
    write_json(output_dir / "raw_footprints.json", [_footprint_dict(item) for item in raw_footprints])
    write_json(
        output_dir / "footprints.json",
        [_footprint_dict(item) for item in footprints + duplicate_footprints],
    )
    with (output_dir / "footprints.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        fields = [
            "footprint_id", "generation_mode", "center_residue_id", "radius",
            "target_area", "total_surface_area", "n_residues", "n_mutable_residues",
            "geodesic_radius", "geodesic_diameter", "mean_rsasa", "mean_plddt",
            "normal_dispersion", "region_id", "quality_status", "duplicate_of",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in footprints + duplicate_footprints:
            writer.writerow({
                "footprint_id": item.footprint_id,
                "generation_mode": item.generation_mode,
                "center_residue_id": item.center_residue_id,
                "radius": item.radius,
                "target_area": item.target_area,
                "total_surface_area": item.total_surface_area,
                "n_residues": len(item.residues),
                "n_mutable_residues": len(item.mutable_residues),
                "geodesic_radius": item.geodesic_radius,
                "geodesic_diameter": item.geodesic_diameter,
                "mean_rsasa": item.mean_rsasa,
                "mean_plddt": item.mean_plddt,
                "normal_dispersion": item.normal_dispersion,
                "region_id": item.region_id,
                "quality_status": ";".join(item.quality_flags) or "PASS",
                "duplicate_of": item.duplicate_of,
            })
    with (output_dir / "footprint_membership.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["footprint_id", "residue_id", "surface_area_in_footprint", "surface_fraction_in_footprint", "is_mutable", "is_protected"])
        for item in footprints + duplicate_footprints:
            for residue_id in item.residues:
                writer.writerow([
                    item.footprint_id, residue_id,
                    item.residue_surface_area.get(residue_id, 0.0),
                    item.residue_surface_fraction.get(residue_id, 0.0),
                    residue_id in item.mutable_residues,
                    residue_id in item.protected_residues,
                ])
    with (output_dir / "footprint_adjacency.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        fields = list(asdict(edges[0])) if edges else [
            "footprint_a", "footprint_b", "weighted_jaccard", "shared_surface_area",
            "center_geodesic_distance", "center_euclidean_distance",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(asdict(item) for item in edges)
    write_json(output_dir / "surface_regions.json", [asdict(item) for item in regions])
    with (output_dir / "surface_region_membership.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["region_id", "footprint_id", "is_representative"])
        for region in regions:
            for footprint_id in region.footprint_ids:
                writer.writerow([region.region_id, footprint_id, footprint_id == region.representative_footprint_id])
    write_json(output_dir / "surface_v2_qc.json", qc)
    return [
        "surface_vertices.csv", "surface_edges.npz", "raw_footprints.json", "footprints.json", "footprints.csv",
        "footprint_membership.csv", "footprint_adjacency.csv", "surface_regions.json",
        "surface_region_membership.csv", "surface_v2_qc.json",
    ]
