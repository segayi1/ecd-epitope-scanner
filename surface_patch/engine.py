from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..models import ProteinContext, RegionRecord, ResidueRecord
from .backend.sr_pointcloud import SRPointCloudBackend
from .center_sampling import eligible_center_vertices, sample_surface_centers
from .diagnostics import write_surface_v2_outputs
from .footprint_dedup import deduplicate_footprints
from .footprint_generator import generate_footprints
from .footprint_graph import build_footprint_graph
from .footprint_qc import evaluate_footprint_quality, global_surface_coverage
from .footprint_residue_mapping import map_footprints_to_residues
from .legacy_adapter import footprint_regions_to_legacy
from .models import Footprint, FootprintEdge, SurfaceGeometry, SurfaceRegion
from .region_builder import build_surface_regions
from .surface_graph import build_surface_graph


@dataclass(slots=True)
class SurfaceV2Result:
    legacy_regions: list[RegionRecord]
    summary: dict[str, Any]
    surface: SurfaceGeometry
    raw_footprints: list[Footprint]
    footprints: list[Footprint]
    duplicate_footprints: list[Footprint]
    footprint_edges: list[FootprintEdge]
    surface_regions: list[SurfaceRegion]
    output_files: list[str]


def run_surface_patch_v2(
    structure: Any,
    chain_id: str,
    ecd_residues: dict[int, Any],
    records: list[ResidueRecord],
    protein: ProteinContext,
    output_dir: Path | None = None,
    confidence: dict[str, Any] | None = None,
    predicted: bool = False,
    parameters: dict[str, Any] | None = None,
    reliability: dict[str, Any] | None = None,
) -> SurfaceV2Result:
    params = parameters or {}
    geometry_config = {
        "probe_radius": 1.4,
        "sampling_points": 48,
        "point_graph_max_distance": 2.5,
        "normal_cosine_min": 0.3,
    } | params.get("surface_geometry", {})
    by_pos = {item.uniprot_pos: item for item in records}
    backend = SRPointCloudBackend()
    excluded_objects = set()
    reliability_rows = {}
    diagnostic_excluded_positions = {}
    if reliability and reliability['qc']['enabled']:
        reliability_rows = {r['residue_id']: r for r in reliability['residues']}
        excluded = set(reliability['qc']['excluded_positions'])
        if reliability['qc']['parameters']['surface_mode'] == 'strict':
            excluded_objects = {id(r) for p, r in ecd_residues.items() if p in excluded}
        # Omit boundary samples but retain their atoms as occluders. This is a
        # conservative whole-residue cap guard, not a reconstructed cut mesh.
        omitted = excluded | set(reliability['qc']['boundary_positions']) | {
            p for p, row in reliability_rows.items()
            if not row.get('use_for_3d_footprint', True)
        }
        diagnostic_excluded_positions = {
            str(position): {
                "status": "OMITTED_BEFORE_POINT_SAMPLING",
                "observed_reliability_state": reliability_rows.get(position, {}).get("state"),
                "observed_exclusion": position in excluded,
                "observed_boundary": position in set(reliability['qc']['boundary_positions']),
            }
            for position in sorted(omitted)
        }
        ecd_residues = {p: r for p, r in ecd_residues.items() if p not in omitted}
    geometry_config["diagnostic_excluded_positions"] = diagnostic_excluded_positions
    surface = backend.build(
        structure,
        chain_id,
        ecd_residues,
        {position: by_pos[position].plddt for position in by_pos},
        geometry_config,
        excluded_residue_objects=excluded_objects,
    )
    graph = build_surface_graph(surface)
    center_config = {
        "center_spacing": 7.0,
        "coverage_radius": 7.0,
        "minimum_component_area": 25.0,
    } | params.get("footprint_centers", {})
    exposed = {item.uniprot_pos for item in records if item.surface_accessible}
    eligible = eligible_center_vertices(
        surface,
        exposed,
        float(center_config["minimum_component_area"]),
    )
    centers = sample_surface_centers(
        surface,
        graph,
        eligible,
        float(center_config["center_spacing"]),
        float(center_config["coverage_radius"]),
    )
    footprint_config = params.get("footprints", {})
    raw = generate_footprints(
        surface,
        graph,
        centers,
        tuple(float(item) for item in footprint_config.get("radii", (9.0, 12.0, 15.0))),
        float(footprint_config.get("target_area", 600.0)),
        float(footprint_config.get("min_area", 350.0)),
        float(footprint_config.get("max_area", 1000.0)),
        float(footprint_config.get("min_radius", 8.0)),
        float(footprint_config.get("max_radius", 16.0)),
    )
    inclusion = params.get("residue_inclusion", {})
    map_footprints_to_residues(
        raw,
        surface,
        records,
        float(inclusion.get("min_surface_area", 8.0)),
        float(inclusion.get("min_surface_fraction", 0.20)),
    )
    evaluate_footprint_quality(raw, surface, graph)
    for footprint in raw:
        for state in ('UNCERTAIN', 'FLEXIBLE_STRUCTURED', 'CONDITIONALLY_FOLDED'):
            if any(reliability_rows.get(p, {}).get('state') == state for p in footprint.residues):
                footprint.quality_flags.append('CONTAINS_' + state)
    unique, duplicates = deduplicate_footprints(
        raw, float(params.get("deduplication", {}).get("threshold", 0.75))
    )
    edges = build_footprint_graph(unique, graph)
    regions = build_surface_regions(
        unique, surface=surface, graph=graph, primary_candidates=raw,
        parameters=params.get("regions", {}),
    )
    legacy = footprint_regions_to_legacy(
        regions, raw, records, protein, confidence, predicted
    )
    coverage = global_surface_coverage(unique, surface, eligible)
    component_areas = {
        str(int(component)): round(float(surface.vertex_area[surface.component_id == component].sum()), 3)
        for component in np.unique(surface.component_id)
    }
    qc = {
        "surface_backend": surface.backend,
        "reliability": reliability['qc'] if reliability else None,
        "surface_vertices": len(surface.vertices),
        "surface_edges": len(surface.edges),
        "surface_components": len(component_areas),
        "geometry_owner_diagnostics": surface.metadata.get("owner_diagnostics", {}),
        "surface_component_area_A2": component_areas,
        "eligible_center_vertices": len(eligible),
        "center_count": len(centers),
        "raw_footprint_count": len(raw),
        "deduplicated_footprint_count": len(unique),
        "duplicate_footprint_count": len(duplicates),
        "surface_region_count": len(regions),
        "region_clustering": {
            "method": "primary-scale constrained Louvain", "seed": 0,
            "parameters": params.get("regions", {}),
            "primary_footprint_count": sum(len(r.primary_footprint_ids) for r in regions),
            "maximum_center_geodesic_extent": max((r.geodesic_extent for r in regions), default=0),
            "display_residue_count": sum(len(r.display_residues) for r in regions),
            "singleton_regions": sum(len(r.primary_footprint_ids) == 1 for r in regions),
        },
        "eligible_surface_coverage": round(coverage, 6),
        "minimum_global_surface_coverage": 0.95,
        "coverage_pass": coverage >= 0.95,
        "parameters": {
            "surface_geometry": geometry_config,
            "footprint_centers": center_config,
            "fixed_radii_A": list(footprint_config.get("radii", (9.0, 12.0, 15.0))),
            "adaptive_target_area_A2": float(footprint_config.get("target_area", 600.0)),
            "residue_inclusion": {
                "min_surface_area_A2": float(inclusion.get("min_surface_area", 8.0)),
                "min_surface_fraction": float(inclusion.get("min_surface_fraction", 0.20)),
            },
            "dedup_weighted_jaccard": float(params.get("deduplication", {}).get("threshold", 0.75)),
        },
    }
    output_files = []
    if output_dir is not None:
        output_files = write_surface_v2_outputs(
            output_dir, surface, raw, unique, duplicates, edges, regions, qc
        )
    coverage_plan = [
        {
            "patch_id": item.id,
            "tier": item.tier,
            "epitope_prior": item.epitope_prior,
            "required_probes": item.required_probe_coverage,
            "footprint_residues": len(item.footprint_residues),
            "envelope_residues": len(item.envelope_residues),
        }
        for item in legacy
    ]
    summary = qc | {
        "partition_mode": "surface_v2",
        "partition_basis": "explicit solvent surface point graph",
        "method": "Fab-footprint Surface V2: overlapping antibody-scale windows from sparse surface geodesics",
        "backend": surface.backend,
        "coverage": round(coverage, 6),
        "actual_regions": len(legacy),
        "patch_radius_a": None,
        "coverage_plan": coverage_plan,
        "warnings": ([] if coverage >= 0.95 else ["Surface-v2 eligible-area coverage is below 95%."]),
        "output_files": output_files,
    }
    return SurfaceV2Result(
        legacy, summary, surface, raw, unique, duplicates, edges, regions, output_files
    )
