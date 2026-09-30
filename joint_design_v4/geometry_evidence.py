"""Attach existing membrane and PAE evidence to footprint hypotheses."""

from __future__ import annotations

from typing import Any

from ..partitioning import _membrane_accessibility
from ..structure_analysis import median_region_pae
from .models import V4SurfaceHypothesis


def annotate_hypothesis_geometry(
    hypotheses: list[V4SurfaceHypothesis],
    records: list[Any],
    *,
    membrane_anchor: tuple[float, float, float] | None = None,
    membrane_normal: tuple[float, float, float] | None = None,
    confidence: dict[str, Any] | None = None,
    predicted: bool = False,
) -> None:
    scores, distances, _ = _membrane_accessibility(records, membrane_anchor, membrane_normal)
    for hypothesis in hypotheses:
        positions = list(hypothesis.residue_weights)
        weights = hypothesis.residue_weights
        area = sum(weights.values())
        known_distance_area = sum(weights[position] for position in positions if position in distances)
        hypothesis.membrane_geometry = {
            "orientation_source": (
                "ECD_anchor_axis_surrogate" if distances else "UNAVAILABLE"
            ),
            "weighted_distance_a": (
                round(sum(weights[position] * distances[position] for position in positions if position in distances) / known_distance_area, 3)
                if known_distance_area else None
            ),
            "weighted_accessibility_proxy": (
                round(sum(weights[position] * scores.get(position, 1.0) for position in positions) / area, 6)
                if area else None
            ),
            "flags": ["INTERMEMBRANE_DISTANCE_UNKNOWN", "FAB_SPECIFIC_REACH_UNKNOWN"],
        }
        hypothesis.structural_uncertainty = {
            "source": "predicted_structure_PAE" if predicted else "experimental_or_unclassified_structure",
            "median_local_pae_a": (
                median_region_pae(positions, records, confidence) if predicted else None
            ),
            "pae_available": bool(predicted and confidence and (confidence.get("pae") or confidence.get("predicted_aligned_error"))),
        }
