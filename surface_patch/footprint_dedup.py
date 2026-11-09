from __future__ import annotations

from .models import Footprint


def weighted_jaccard(left: Footprint, right: Footprint) -> float:
    residues = set(left.residue_surface_area) | set(right.residue_surface_area)
    numerator = sum(
        min(left.residue_surface_area.get(item, 0.0), right.residue_surface_area.get(item, 0.0))
        for item in residues
    )
    denominator = sum(
        max(left.residue_surface_area.get(item, 0.0), right.residue_surface_area.get(item, 0.0))
        for item in residues
    )
    return numerator / denominator if denominator else 0.0


def _representative_key(footprint: Footprint) -> tuple[float, float, float, int, str]:
    confidence = footprint.mean_plddt if footprint.mean_plddt is not None else 100.0
    medoid_offset = footprint.geodesic_radius
    return (
        -confidence,
        -footprint.total_surface_area,
        medoid_offset,
        len(footprint.quality_flags),
        footprint.footprint_id,
    )


def deduplicate_footprints(
    footprints: list[Footprint], threshold: float = 0.75
) -> tuple[list[Footprint], list[Footprint]]:
    """Greedy deterministic deduplication with retained audit records."""
    kept: list[Footprint] = []
    removed: list[Footprint] = []
    for candidate in sorted(footprints, key=_representative_key):
        duplicate = next(
            (item for item in kept if weighted_jaccard(candidate, item) >= threshold),
            None,
        )
        if duplicate is None:
            kept.append(candidate)
        else:
            candidate.duplicate_of = duplicate.footprint_id
            removed.append(candidate)
    kept.sort(key=lambda item: item.footprint_id)
    removed.sort(key=lambda item: item.footprint_id)
    return kept, removed
