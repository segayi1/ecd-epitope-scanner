from __future__ import annotations

from collections import Counter
from statistics import mean
from typing import Any

from .. import config
from ..models import ProteinContext, RegionRecord, ResidueRecord
from ..structure_analysis import median_region_pae
from .models import Footprint, SurfaceRegion


def _region_color(index: int) -> str:
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#6F4E7C", "#666666", "#8C564B", "#BCBD22"]
    if index <= len(palette):
        return palette[index - 1]
    import colorsys
    rgb = colorsys.hsv_to_rgb((index * .61803398875) % 1, .65, .8)
    return "#" + "".join(f"{round(c*255):02X}" for c in rgb)


def footprint_to_legacy_region(
    footprint_or_region: Footprint | SurfaceRegion,
    footprints: list[Footprint],
    records: list[ResidueRecord],
    protein: ProteinContext,
    confidence: dict[str, Any] | None = None,
    predicted: bool = False,
) -> RegionRecord:
    """Adapt one V2 footprint/region for callers that still consume patches."""
    if isinstance(footprint_or_region, Footprint):
        footprint = footprint_or_region
        region = SurfaceRegion(
            region_id=footprint.region_id or "SR001",
            footprint_ids=[footprint.footprint_id],
            representative_footprint_id=footprint.footprint_id,
            residues_union=list(footprint.residues),
            surface_area_union=footprint.total_surface_area,
            centroid_xyz=footprint.centroid_xyz,
            domain_ids=list(footprint.structural_domains),
            quality_flags=list(footprint.quality_flags),
        )
        available = footprints if any(item.footprint_id == footprint.footprint_id for item in footprints) else footprints + [footprint]
    else:
        region = footprint_or_region
        available = footprints
    adapted = footprint_regions_to_legacy(
        [region], available, records, protein, confidence, predicted
    )
    if not adapted:
        raise ValueError("Cannot adapt an empty SurfaceRegion")
    return adapted[0]


def footprint_regions_to_legacy(
    regions: list[SurfaceRegion],
    footprints: list[Footprint],
    records: list[ResidueRecord],
    protein: ProteinContext,
    confidence: dict[str, Any] | None,
    predicted: bool,
) -> list[RegionRecord]:
    by_footprint = {item.footprint_id: item for item in footprints}
    by_pos = {item.uniprot_pos: item for item in records}
    output: list[RegionRecord] = []
    for index, surface_region in enumerate(regions, start=1):
        representative = by_footprint[surface_region.representative_footprint_id]
        items = [by_pos[item] for item in representative.residues if item in by_pos]
        mutable = [item for item in representative.mutable_residues if item in by_pos]
        protected = [
            item for item in surface_region.residues_union
            if item in by_pos and not by_pos[item].folding_screen_pass
        ]
        domain = Counter(item.domain for item in items).most_common(1)[0][0] if items else "ECD"
        pae = median_region_pae(representative.residues, records, confidence) if predicted else None
        secondary = Counter(item.secondary_structure for item in items)
        structural_character = (
            "beta-sheet face" if secondary["E"] > max(secondary["H"], secondary["C"])
            else "helical face" if secondary["H"] > secondary["C"] else "mixed surface"
        )
        region = RegionRecord(
            id=index,
            name=f"{domain} / surface region {index}",
            domain=domain,
            membrane_position="geometry-independent",
            structural_character=structural_character,
            residues=mutable,
            surface_area=round(surface_region.surface_area_union, 2),
            geodesic_diameter=round(representative.geodesic_diameter, 2),
            interface_fraction=round(mean(item.interface for item in items), 3) if items else 0.0,
            pae_median=pae,
            auto_design_allowed=bool(mutable) and (not predicted or pae is None or pae <= 10.0),
            rationale=(
                f"surface-v2 geodesic region; representative={representative.footprint_id}; "
                f"overlapping footprints={len(surface_region.footprint_ids)}; "
                f"surface area={surface_region.surface_area_union:.1f} A2"
            ),
            color=_region_color(index),
            domain_id=domain,
            epitope_prior=representative.epitope_prior,
            prior_components={
                "surface_exposure": representative.mean_rsasa,
                "surface_area": round(representative.total_surface_area, 3),
                "normal_coherence": round(1.0 - representative.normal_dispersion, 4),
                "protected_fraction": round(
                    len(representative.protected_residues) / max(1, len(representative.residues)), 4
                ),
            },
            tier="MEDIUM",
            required_probe_coverage=1,
            partition_mode="surface_v2",
            partition_basis="surface_point_graph_geodesic",
            envelope_residues=surface_region.display_residues or surface_region.residues_union,
            mutable_residues=mutable,
            protected_surface_residues=protected,
            # The legacy designer reads footprint_residues directly. Its input
            # must be the mutable subset; full geometry stays in Footprint.
            footprint_residues=mutable,
            candidate_footprints=[by_footprint[item].residues for item in surface_region.footprint_ids],
            footprint_center=representative.center_residue_id,
            surface_connectivity="connected",
            normal_coherence=round(1.0 - representative.normal_dispersion, 4),
            quality_flags=representative.quality_flags,
        )
        output.append(region)
    ordered = sorted(output, key=lambda item: (-(item.epitope_prior or 0.0), item.id))
    total = len(ordered)
    for rank, item in enumerate(ordered):
        item.tier = "HIGH" if rank < max(1, round(total * config.TIER_HIGH_FRACTION)) else (
            "MEDIUM" if rank < max(1, round(total * (config.TIER_HIGH_FRACTION + config.TIER_MEDIUM_FRACTION))) else "LOW"
        )
    return output
