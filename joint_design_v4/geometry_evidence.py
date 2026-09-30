"""Attach explicit membrane-orientation and PAE evidence to hypotheses."""

from __future__ import annotations

import math
from typing import Any

from ..partitioning import _membrane_accessibility
from ..structure_analysis import median_region_pae
from .models import V4SurfaceHypothesis


def _unit(vector: tuple[float, ...] | list[float] | None) -> tuple[float, float, float] | None:
    if vector is None or len(vector) != 3:
        return None
    values = tuple(float(value) for value in vector)
    norm = math.sqrt(sum(value * value for value in values))
    if norm < 1e-8:
        return None
    return tuple(value / norm for value in values)


def _weighted_mean(values: dict[int, float], weights: dict[int, float]) -> float | None:
    matched = {position: value for position, value in values.items() if position in weights}
    denominator = sum(weights[position] for position in matched)
    if denominator <= 0:
        return None
    return sum(value * weights[position] for position, value in matched.items()) / denominator


def annotate_hypothesis_geometry(
    hypotheses: list[V4SurfaceHypothesis],
    records: list[Any],
    *,
    membrane_anchor: tuple[float, float, float] | None = None,
    membrane_normal: tuple[float, float, float] | None = None,
    assay_orientation_geometry: dict[str, Any] | None = None,
    confidence: dict[str, Any] | None = None,
    predicted: bool = False,
) -> None:
    """Attach structure-coordinate proxies without inventing membrane reach.

    A transmembrane-derived plane supports signed and absolute plane distances.
    The assay C-terminal-axis fallback supports only relative projections; it
    must not be converted into a membrane accessibility or Fab-reach score.
    """
    record_by_position = {item.uniprot_pos: item for item in records}
    plane_normal = _unit(membrane_normal)
    has_plane = membrane_anchor is not None and plane_normal is not None
    if has_plane:
        scores, plane_distances, source = _membrane_accessibility(
            records, membrane_anchor, membrane_normal
        )
        source = "mapped_tm_helix_axis_surrogate"
    else:
        scores, plane_distances, source = {}, {}, "UNAVAILABLE"

    orientation = assay_orientation_geometry or {}
    relative_anchor = _unit(orientation.get("extracellular_axis"))
    relative_origin_raw = orientation.get("anchor_xyz")
    relative_origin = (
        tuple(float(value) for value in relative_origin_raw)
        if isinstance(relative_origin_raw, (list, tuple)) and len(relative_origin_raw) == 3
        else None
    )
    has_relative_axis = (
        not has_plane
        and orientation.get("status") == "relative_axis_only"
        and relative_anchor is not None
        and relative_origin is not None
    )

    for hypothesis in hypotheses:
        weights = hypothesis.residue_weights
        positions = list(weights)
        ca_points = {
            position: record_by_position[position].ca
            for position in positions
            if position in record_by_position and record_by_position[position].ca is not None
        }
        absolute_distance = _weighted_mean(
            {position: plane_distances[position] for position in positions if position in plane_distances},
            weights,
        )
        signed_distance = None
        relative_projection = None
        if has_plane and membrane_anchor is not None and plane_normal is not None:
            anchor = tuple(float(value) for value in membrane_anchor)
            signed_by_position = {
                position: sum(
                    (float(point[axis]) - anchor[axis]) * plane_normal[axis]
                    for axis in range(3)
                )
                for position, point in ca_points.items()
            }
            signed_distance = _weighted_mean(signed_by_position, weights)
        elif has_relative_axis and relative_origin is not None and relative_anchor is not None:
            projection_by_position = {
                position: sum(
                    (float(point[axis]) - relative_origin[axis]) * relative_anchor[axis]
                    for axis in range(3)
                )
                for position, point in ca_points.items()
            }
            relative_projection = _weighted_mean(projection_by_position, weights)

        orientation_vector = _unit(hypothesis.mean_normal)
        reference_axis = plane_normal if has_plane else relative_anchor if has_relative_axis else None
        normal_cosine = (
            sum(orientation_vector[index] * reference_axis[index] for index in range(3))
            if orientation_vector is not None and reference_axis is not None
            else None
        )
        hypothesis.membrane_geometry = {
            "orientation_source": source if has_plane else (
                str(orientation.get("axis_source", "C-terminal anchor surrogate"))
                if has_relative_axis else "UNAVAILABLE"
            ),
            "orientation_confidence": (
                "predicted_structure_tm_axis_surrogate" if predicted else "experimental_structure_tm_axis_surrogate"
            ) if has_plane else (
                "relative_axis_surrogate_only" if has_relative_axis else "unavailable"
            ),
            "weighted_plane_distance_a": round(absolute_distance, 3) if absolute_distance is not None else None,
            "weighted_distance_a": round(absolute_distance, 3) if absolute_distance is not None else None,
            "weighted_signed_plane_distance_a": round(signed_distance, 3) if signed_distance is not None else None,
            "relative_axis_projection_from_anchor_a": (
                round(relative_projection, 3) if relative_projection is not None else None
            ),
            "surface_normal_reference_axis_cosine": (
                round(max(-1.0, min(1.0, normal_cosine)), 6) if normal_cosine is not None else None
            ),
            "weighted_membrane_position_factor_proxy": (
                round(_weighted_mean(scores, weights), 6)
                if has_plane and _weighted_mean(scores, weights) is not None else None
            ),
            # Retain the old key for V4.0 readers, but mark its limited meaning.
            "weighted_accessibility_proxy": (
                round(_weighted_mean(scores, weights), 6)
                if has_plane and _weighted_mean(scores, weights) is not None else None
            ),
            "accessibility_proxy_interpretation": "distance-based membrane-position prior; not molecular or assay accessibility",
            "membrane_plane_surrogate_available": bool(has_plane),
            "physical_membrane_plane_inferred": False,
            "flags": (
                ["TM_HELIX_AXIS_USED_AS_MEMBRANE_PLANE_SURROGATE", "INTERMEMBRANE_DISTANCE_UNKNOWN", "FAB_SPECIFIC_REACH_UNKNOWN"]
                if has_plane else
                ["RELATIVE_AXIS_ONLY", "NO_MEMBRANE_PLANE_INFERRED", "INTERMEMBRANE_DISTANCE_UNKNOWN", "FAB_SPECIFIC_REACH_UNKNOWN"]
                if has_relative_axis else
                ["MEMBRANE_ORIENTATION_UNAVAILABLE", "INTERMEMBRANE_DISTANCE_UNKNOWN", "FAB_SPECIFIC_REACH_UNKNOWN"]
            ),
        }
        local_pae = median_region_pae(positions, records, confidence) if predicted else None
        hypothesis.structural_uncertainty = {
            "source": "predicted_structure_PAE" if predicted else "experimental_or_unclassified_structure",
            "median_local_pae_a": local_pae,
            "pae_available": local_pae is not None,
            "mean_plddt": hypothesis.mean_plddt,
            "minimum_plddt": hypothesis.minimum_plddt,
            "low_confidence_surface_fraction": hypothesis.low_confidence_surface_fraction,
            "confidence_note": (
                "Experimental B factors are not interpreted as pLDDT."
                if not predicted else "pLDDT and PAE remain structural confidence evidence, not epitope evidence."
            ),
        }
