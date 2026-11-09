from __future__ import annotations

from collections import defaultdict

import numpy as np

from ..models import ResidueRecord
from .models import Footprint, SurfaceGeometry


def map_footprints_to_residues(
    footprints: list[Footprint],
    surface: SurfaceGeometry,
    records: list[ResidueRecord],
    min_surface_area: float = 8.0,
    min_surface_fraction: float = 0.20,
) -> list[Footprint]:
    by_pos = {item.uniprot_pos: item for item in records}
    totals: dict[int, float] = defaultdict(float)
    for residue_id, area in zip(surface.owner_residue, surface.vertex_area):
        totals[int(residue_id)] += float(area)
    for footprint in footprints:
        local: dict[int, float] = defaultdict(float)
        for vertex_id in footprint.vertex_ids:
            local[int(surface.owner_residue[vertex_id])] += float(surface.vertex_area[vertex_id])
        footprint.residue_surface_area = {key: round(value, 6) for key, value in sorted(local.items())}
        footprint.residue_surface_fraction = {
            key: round(value / totals[key], 6) if totals[key] else 0.0
            for key, value in sorted(local.items())
        }
        footprint.residues = [
            key for key in sorted(local)
            if local[key] >= min_surface_area or footprint.residue_surface_fraction[key] >= min_surface_fraction
        ]
        footprint.mutable_residues = [
            key for key in footprint.residues if key in by_pos and by_pos[key].folding_screen_pass
        ]
        footprint.protected_residues = [key for key in footprint.residues if key not in footprint.mutable_residues]
        selected_records = [by_pos[key] for key in footprint.residues if key in by_pos]
        footprint.mean_rsasa = round(
            float(np.mean([item.rsasa for item in selected_records])) if selected_records else 0.0, 4
        )
        plddt = [item.plddt for item in selected_records if item.plddt is not None]
        footprint.mean_plddt = round(float(np.mean(plddt)), 2) if plddt else None
        footprint.minimum_plddt = round(float(min(plddt)), 2) if plddt else None
        footprint.structural_domains = sorted({item.domain for item in selected_records})
        # Keep biological scoring out of geometry; this score only preserves
        # the current downstream ordering signal for the compatibility layer.
        footprint.epitope_prior = round(
            float(np.mean([min(1.0, item.rsasa / 0.5) for item in selected_records]))
            if selected_records else 0.0, 4
        )
    return footprints
