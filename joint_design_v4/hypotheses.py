"""Overlapping Surface V2 footprints as testable surface hypotheses.

Responses are geometry and chemistry proxies, not predicted binding changes.
"""

from __future__ import annotations

import heapq
import math
from copy import deepcopy
from typing import Any

from ..surface_patch.footprint_dedup import deduplicate_footprints
from .models import V4Candidate, V4SurfaceHypothesis


def hypotheses_from_surface_v2(
    surface_result: Any,
    records: list[Any] | None = None,
    epitope_evidence: list[Any] | None = None,
) -> list[V4SurfaceHypothesis]:
    if surface_result is None:
        return []
    records_by_position = {item.uniprot_pos: item for item in (records or [])}
    prior_by_position = {
        item.residue.uniprot_pos: item.surface_percentile
        for item in (epitope_evidence or [])
        if item.surface_percentile is not None
    }
    scales: dict[int, dict[int, set[int]]] = {}
    for footprint in getattr(surface_result, "raw_footprints", []):
        radius = getattr(footprint, "radius", None)
        if radius is not None:
            scales.setdefault(footprint.center_vertex_id, {})[int(round(radius))] = set(footprint.residues)
    hypotheses = []
    raw_footprints = list(getattr(surface_result, "raw_footprints", []))
    primary_scale = [
        item for item in raw_footprints
        if item.generation_mode == "adaptive_area"
        or (item.generation_mode == "fixed_radius" and item.radius is not None and abs(float(item.radius) - 12.0) < 1e-6)
    ]
    if primary_scale:
        primary_scale = [deepcopy(item) for item in primary_scale]
        candidate_footprints, _ = deduplicate_footprints(primary_scale)
    else:
        candidate_footprints = list(surface_result.footprints)
    for footprint in candidate_footprints:
        if "ERROR_DISCONNECTED_SURFACE" in footprint.quality_flags:
            continue
        weights = {
            int(position): max(0.0, float(area))
            for position, area in footprint.residue_surface_area.items()
            if float(area) > 0 and int(position) in set(footprint.residues)
        }
        if not weights:
            continue
        total = sum(weights.values())
        membership = {
            position: float(getattr(footprint, "residue_surface_fraction", {}).get(position, 0.0))
            for position in weights
        }
        center_scales = scales.get(getattr(footprint, "center_vertex_id", -1), {})
        core = center_scales.get(9, set()) & set(weights)
        if not core:
            core = {position for position, fraction in membership.items() if fraction >= 0.5}
        boundary = (center_scales.get(15, set()) | set(weights)) - core
        glycan_area = sum(
            weight for position, weight in weights.items()
            if getattr(records_by_position.get(position), "glycan_shielded", False)
        )
        prior_area = sum(
            weight * float(prior_by_position[position])
            for position, weight in weights.items() if position in prior_by_position
        )
        scored_area = sum(weight for position, weight in weights.items() if position in prior_by_position)
        hypotheses.append(V4SurfaceHypothesis(
            hypothesis_id=footprint.footprint_id,
            residue_weights=weights,
            region_id=footprint.region_id,
            source_footprint_id=footprint.footprint_id,
            quality_flags=list(footprint.quality_flags),
            center_position=getattr(footprint, "center_residue_id", None),
            generation_mode=getattr(footprint, "generation_mode", None),
            geodesic_radius_a=getattr(footprint, "geodesic_radius", None),
            geodesic_diameter_a=getattr(footprint, "geodesic_diameter", None),
            realized_area_a2=getattr(footprint, "total_surface_area", None),
            target_area_a2=getattr(footprint, "target_area", None),
            center_xyz=getattr(footprint, "centroid_xyz", None),
            mean_normal=getattr(footprint, "mean_normal", None),
            mutable_positions=sorted(set(getattr(footprint, "mutable_residues", [])) & set(weights)),
            protected_positions=sorted(set(getattr(footprint, "protected_residues", [])) & set(weights)),
            normal_dispersion=getattr(footprint, "normal_dispersion", None),
            mean_plddt=getattr(footprint, "mean_plddt", None),
            minimum_plddt=getattr(footprint, "minimum_plddt", None),
            low_confidence_surface_fraction=getattr(footprint, "low_confidence_surface_fraction", None),
            surface_component_id=getattr(footprint, "surface_component_id", None),
            core_positions=sorted(core),
            boundary_positions=sorted(boundary),
            membership_fraction=membership,
            discotope_prior=round(prior_area / scored_area, 6) if scored_area else None,
            glycan_fraction=round(glycan_area / total, 6),
        ))
    return sorted(hypotheses, key=lambda item: item.hypothesis_id)


def _distances(graph: dict[int, dict[int, float]], source: int, cutoff: float) -> dict[int, float]:
    distances = {source: 0.0}
    queue = [(0.0, source)]
    while queue:
        distance, position = heapq.heappop(queue)
        if distance != distances[position]:
            continue
        for neighbor, step in graph.get(position, {}).items():
            new_distance = distance + step
            if new_distance <= cutoff and new_distance < distances.get(neighbor, float("inf")):
                distances[neighbor] = new_distance
                heapq.heappush(queue, (new_distance, neighbor))
    return distances


def attach_hypothesis_responses(
    candidates: list[V4Candidate],
    hypotheses: list[V4SurfaceHypothesis],
    graph: dict[int, dict[int, float]],
    radius_a: float = 15.0,
    decay_length_a: float = 6.0,
) -> None:
    """Map each probe to each footprint using surface geodesic reach.

    Each mutation contributes a surface-area weighted geodesic overlap times
    its explicit chemical disruption strength. Contributions combine by a
    noisy-OR rule. No experimental effect size or probability is implied.
    """
    if radius_a <= 0 or decay_length_a <= 0:
        raise ValueError("Hypothesis response radius and decay length must be positive")
    cache: dict[int, dict[int, float]] = {}
    for candidate in candidates:
        positions = [option.residue.uniprot_pos for option in candidate.mutations]
        for position in positions:
            if position is not None and position not in cache:
                cache[position] = _distances(graph, position, radius_a)
        response = {}
        for hypothesis in hypotheses:
            total = sum(hypothesis.residue_weights.values())
            if total <= 0:
                continue
            survival = 1.0
            for option in candidate.mutations:
                position = option.residue.uniprot_pos
                if position is None:
                    continue
                overlap = sum(
                    weight * math.exp(-cache[position][residue] / decay_length_a)
                    for residue, weight in hypothesis.residue_weights.items()
                    if residue in cache[position]
                ) / total
                chemical = option.chemical_disruption_score or (
                    0.8 if option.chemistry_class in {"charge_reversal", "hotspot_removal"}
                    else 0.5
                )
                survival *= 1.0 - min(1.0, max(0.0, overlap * chemical))
            response[hypothesis.hypothesis_id] = round(1.0 - survival, 6)
        candidate.hypothesis_response = response
