"""Shadow-only fusion of independently validated surface-field mechanisms."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components, dijkstra

from .chemistry import CHANNELS
from .directions import VECTOR_GROUPS, consistent_four
from .null_store import GEOMETRY_CHANNELS

PROFILE = "field_fusion_v1_shadow"
PREPARATION_PROTOCOL = "heavy_atom_preserving_v1"
PRODUCTION_SLOTS = 5
NULL_BLOCKS_PER_ROLE = 5
VECTOR_COMPONENTS = {name for group in VECTOR_GROUPS.values() for name in group}
CHEMISTRY_SCALARS = tuple(name for name in CHANNELS if name not in VECTOR_COMPONENTS)
SCALAR_CHANNELS = (*tuple(name for name in GEOMETRY_CHANNELS if name not in VECTOR_COMPONENTS),
                   *CHEMISTRY_SCALARS, "delta_phi")


def protocol_id(*, reference_sha256: str, source_sha256: str, pdb2pqr_sha256: str,
                apbs_sha256: str | None, query_offset_a: float = 1.5) -> str:
    payload = {
        "profile": PROFILE, "reference": reference_sha256, "source": source_sha256,
        "pdb2pqr": pdb2pqr_sha256, "apbs": apbs_sha256,
        "preparation": PREPARATION_PROTOCOL, "force_field": "PARSE",
        "pH": 7.4, "temperature_K": 298.15, "salt_molar": 0.15,
        "query_offset_a": query_offset_a, "fusion": "shared-index-maxstat-4of5-v1",
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def pair_null_blocks(values: np.ndarray) -> np.ndarray:
    """Create five disjoint paired differences from ten frozen WT slots."""
    values = np.asarray(values, dtype=float)
    if values.shape[0] != 10:
        raise ValueError("each estimate/holdout role must contain exactly ten frozen WT slots")
    return values[0::2] - values[1::2]


def validate_null_roles(train: list[dict], holdout: list[dict], production: list[dict]) -> dict:
    """Validate disjoint run identities while reporting repeated coordinate content."""
    identities = {}; content = {}
    for role, rows in (("null_train", train), ("null_holdout", holdout), ("production", production)):
        if role != "production" and len(rows) != 10:
            raise ValueError(f"{role} requires exactly ten frozen WT slots")
        identities[role] = [tuple(row.get("wt_execution_id",
            (row.get("batch_path", f"{role}:{i}"), row.get("group", 0), row.get("run", i))))
            for i, row in enumerate(rows)]
        if len(set(identities[role])) != len(identities[role]):
            raise ValueError(f"{role} reuses a FoldX run identity")
        content[role] = [row.get("coordinate_sha256") or row["wt_sha256"] for row in rows]
    if set(identities["null_train"]) & set(identities["null_holdout"]):
        raise ValueError("estimate and holdout FoldX run identities overlap")
    production_ids = set(identities["production"])
    if production_ids & (set(identities["null_train"]) | set(identities["null_holdout"])):
        raise ValueError("production and null FoldX run identities overlap")
    return {role: {"slots": len(identities[role]), "unique_run_identities": len(set(identities[role])),
                   "unique_coordinate_models": len(set(content[role]))}
            for role in identities}


def _quorum_scalar(values: np.ndarray, threshold: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    above = values > threshold
    below = values < -threshold
    return (above.sum(axis=0) >= 4) | (below.sum(axis=0) >= 4)


def _channel_stack(fields: dict[str, np.ndarray], names: tuple[str, ...], vector_groups: dict[str, tuple[str, ...]]):
    rows = []
    for name in names:
        if name not in fields:
            continue
        value = np.asarray(fields[name], float)
        if value.ndim != 2:
            raise ValueError(f"{name} must have shape repeat x surface-point")
        rows.append((name, value))
    vectors = []
    for label, components in vector_groups.items():
        if not all(name in fields for name in components):
            continue
        value = np.stack([np.asarray(fields[name], float) for name in components], axis=-1)
        if value.ndim != 3:
            raise ValueError(f"{label} must have shape repeat x surface-point x component")
        vectors.append((label, components, value))
    return rows, vectors


def fuse_scope(production: dict[str, np.ndarray], estimate: dict[str, np.ndarray],
               heldout: dict[str, np.ndarray], numerical_error: dict[str, np.ndarray],
               area: np.ndarray, *, scalar_names: tuple[str, ...] = SCALAR_CHANNELS,
               vector_groups: dict[str, tuple[str, ...]] | None = None,
               protocol: str, minimum_vector_cosine: float = 0.0) -> dict:
    """Return lower/unknown masks for one predeclared mechanism scope.

    All channel maxima for each null block retain the same block index. A
    candidate channel must itself meet the four-of-five directional quorum;
    maxima calibrate multiplicity but never let channels lend votes to each other.
    """
    area = np.asarray(area, float)
    if area.ndim != 1 or not np.isfinite(area).all() or np.any(area < 0):
        raise ValueError("surface area weights must be finite, nonnegative and one-dimensional")
    vectors_declared = VECTOR_GROUPS if vector_groups is None else vector_groups
    p_scalars, p_vectors = _channel_stack(production, scalar_names, vectors_declared)
    e_scalars, e_vectors = _channel_stack(estimate, scalar_names, vectors_declared)
    h_scalars, h_vectors = _channel_stack(heldout, scalar_names, vectors_declared)
    expected_scalars = set(scalar_names)
    if any(any(name not in fields for name in expected_scalars)
           for fields in (production, estimate, heldout)):
        raise ValueError("a predeclared scalar channel is missing from a role")
    if ({x[0] for x in p_scalars} != expected_scalars or {x[0] for x in e_scalars} != expected_scalars
            or {x[0] for x in h_scalars} != expected_scalars):
        raise ValueError("production, estimate and holdout scalar channel scopes differ")
    if {x[0] for x in p_scalars} != {x[0] for x in e_scalars} or {x[0] for x in p_scalars} != {x[0] for x in h_scalars}:
        raise ValueError("production, estimate and holdout scalar channel names differ")
    if {x[0] for x in p_vectors} != {x[0] for x in e_vectors} or {x[0] for x in p_vectors} != {x[0] for x in h_vectors}:
        raise ValueError("production, estimate and holdout vector scopes differ")
    for _, components in vectors_declared.items():
        if any(any(name not in fields for name in components)
               for fields in (production, estimate, heldout)):
            raise ValueError("a predeclared vector channel is missing from a role")
    n = len(area)
    all_rows = [("scalar", name, values) for name, values in p_scalars]
    all_rows += [("vector", label, values) for label, _, values in p_vectors]
    estimate_by = {name: value for name, value in e_scalars}
    estimate_by.update({label: np.linalg.norm(value, axis=2) for label, _, value in e_vectors})
    held_by = {name: value for name, value in h_scalars}
    held_by.update({label: np.linalg.norm(value, axis=2) for label, _, value in h_vectors})
    production_by = {name: value for name, value in p_scalars}
    production_by.update({label: np.linalg.norm(value, axis=2) for label, _, value in p_vectors})
    error_by = {}
    for kind, name, _ in all_rows:
        if kind == "scalar":
            error_by[name] = np.asarray(numerical_error.get(name, np.full(n, np.nan)), float)
        else:
            components = next(parts for label, parts, _ in p_vectors if label == name)
            component_errors = [np.asarray(numerical_error.get(c, np.full(n, np.nan)), float) for c in components]
            error_by[name] = np.linalg.norm(np.stack(component_errors, axis=0), axis=0)
        if error_by[name].shape != (n,):
            raise ValueError(f"numerical error shape mismatch for {name}")

    scales = {}; calibration_valid = {}; production_valid = {}; heldout_valid = {}
    null_complete = {}; error_valid = {}
    for kind, name, values in all_rows:
        null = estimate_by[name]
        prod_magnitude = production_by[name]
        if prod_magnitude.shape != (PRODUCTION_SLOTS, n) or null.shape != (NULL_BLOCKS_PER_ROLE, n) or held_by[name].shape != null.shape:
            raise ValueError(f"fixed five-slot evidence shape mismatch for {name}")
        null_magnitude = np.abs(null)
        q95 = np.quantile(null_magnitude, .95, axis=0)
        err = error_by[name]
        scale = np.maximum(q95, err)
        scale = np.maximum(scale, np.finfo(float).tiny)
        # Five finite logical WT-difference slots and a matching numerical envelope are required;
        # content duplicates are audited separately and are not counted as new independent samples.
        null_complete[name] = np.isfinite(null).all(axis=0)
        error_valid[name] = np.isfinite(err) & (err >= 0)
        calibration_ok = null_complete[name] & error_valid[name]
        production_ok = np.isfinite(prod_magnitude).sum(axis=0) >= 4
        heldout_ok = np.isfinite(held_by[name]).all(axis=0)
        scales[name] = scale
        calibration_valid[name] = calibration_ok
        production_valid[name] = production_ok
        heldout_valid[name] = heldout_ok

    if not all_rows:
        empty = np.zeros(n, bool)
        return {"protocol_id": protocol, "status": "UNKNOWN", "positive": empty,
                "known": empty, "unknown": np.ones(n, bool), "area_a2": 0.0,
                "channels": {}, "reason": "no calibrated channels",
                "unknown_reason_masks": {"no_calibrated_channels": np.ones(n, bool)},
                "channel_unknown_reasons": {}}

    # Calibrate the maximum over the channels that have complete null and error
    # evidence at each point. A missing channel remains unknown, but cannot erase
    # a positive established by another calibrated channel.
    names = [name for _, name, _ in all_rows]
    standardized = np.stack([np.abs(estimate_by[name]) / scales[name] for name in names])
    calibratable = np.stack([calibration_valid[name] for name in names])
    null_max = np.full((NULL_BLOCKS_PER_ROLE, n), np.nan)
    joint_threshold = np.full(n, np.nan)
    for point in range(n):
        active = calibratable[:, point]
        if active.any():
            values = standardized[active, :, point]
            null_max[:, point] = np.max(values, axis=0)
            # Do not let an all-zero finite null erase the empirical single-channel floor.
            joint_threshold[point] = max(float(np.quantile(null_max[:, point], .95)), 1.0)
    channel_positive = {}; channel_unknown = {}; channel_audits = {}; channel_unknown_reasons = {}
    for kind, name, values in all_rows:
        limit = joint_threshold * scales[name]
        if kind == "scalar":
            positive = _quorum_scalar(values, limit)
            false_holdout = _quorum_scalar(held_by[name], limit)
        else:
            vectors = next(value for label, _, value in p_vectors if label == name)
            held_vectors = next(value for label, _, value in h_vectors if label == name)
            eligible = np.linalg.norm(vectors, axis=2) > limit
            held_eligible = np.linalg.norm(held_vectors, axis=2) > limit
            positive = consistent_four(vectors, eligible, minimum_vector_cosine)
            false_holdout = consistent_four(held_vectors, held_eligible, minimum_vector_cosine)
        available = (calibration_valid[name] & production_valid[name] & heldout_valid[name]
                     & np.isfinite(joint_threshold))
        channel_positive[name] = positive & available & ~false_holdout
        channel_unknown[name] = ~available | false_holdout
        channel_unknown_reasons[name] = {
            "estimate_null_incomplete": ~null_complete[name],
            "numerical_error_unavailable": ~error_valid[name],
            "production_repeats_incomplete": ~production_valid[name],
            "holdout_null_incomplete": ~heldout_valid[name],
            "holdout_quorum_failure": false_holdout,
            "joint_threshold_unavailable": ~np.isfinite(joint_threshold),
        }
        channel_audits[name] = {
            "positive_area_a2": float(area[channel_positive[name]].sum()),
            "unknown_area_a2": float(area[channel_unknown[name]].sum()),
            "holdout_quorum_failures": int(false_holdout.sum()),
            "estimate_blocks": NULL_BLOCKS_PER_ROLE,
            "holdout_blocks": NULL_BLOCKS_PER_ROLE,
            "noise_scale_median": float(np.nanmedian(scales[name][available])) if available.any() else None,
            "unknown_reason_area_a2": {reason: float(area[mask].sum())
                                       for reason, mask in channel_unknown_reasons[name].items()},
        }
    positive = np.logical_or.reduce(list(channel_positive.values()))
    all_channels_known = np.logical_and.reduce([~value for value in channel_unknown.values()])
    known = positive | all_channels_known
    reason_names = tuple(next(iter(channel_unknown_reasons.values())))
    unknown_reason_masks = {
        reason: np.logical_or.reduce([channel_unknown_reasons[name][reason] for name in channel_unknown_reasons])
        & ~positive for reason in reason_names
    }
    return {
        "protocol_id": protocol, "status": "RESOLVED" if known.all() else "PARTIAL",
        "positive": positive, "known": known, "unknown": ~known,
        "joint_threshold": joint_threshold, "channel_positive": channel_positive,
        "channel_unknown": channel_unknown, "channel_scales": scales,
        "channel_unknown_reasons": channel_unknown_reasons,
        "unknown_reason_masks": unknown_reason_masks,
        "channel_negative": {name: ~channel_positive[name] & ~channel_unknown[name]
                             for name in channel_positive},
        "channel_audits": channel_audits, "area_a2": float(area[positive].sum()),
        "unknown_area_a2": float(area[~known].sum()),
        "channels": [name for _, name, _ in all_rows],
        "interpretation": "computational lower bound; finite empirical null and holdout, not global error control",
    }


def geometry_summary(mask: np.ndarray, area: np.ndarray, graph, seed_mask: np.ndarray | None = None) -> dict:
    mask = np.asarray(mask, bool); area = np.asarray(area, float)
    if mask.shape != area.shape or graph.shape != (len(area), len(area)):
        raise ValueError("mask, area and surface graph dimensions differ")
    count, labels = connected_components(graph, directed=False)
    components = []
    for label in range(count):
        component = mask & (labels == label)
        if component.any(): components.append(float(area[component].sum()))
    positive_ids = np.flatnonzero(mask)
    distance = None
    if seed_mask is not None:
        seeds = np.flatnonzero(np.asarray(seed_mask, bool))
        if seeds.size and positive_ids.size:
            distances = np.asarray(dijkstra(graph, indices=seeds, min_only=True), float).reshape(-1)
            distance = distances[positive_ids]
    area_sum = float(area[mask].sum())
    result = {"positive_area_a2": area_sum, "equivalent_area_radius_a": float(np.sqrt(area_sum / np.pi)),
              "positive_samples": int(mask.sum()), "connected_component_areas_a2": sorted(components, reverse=True),
              "largest_component_a2": max(components, default=0.0),
              "interpretation": "equivalent-area radius only; the positive region remains the saved irregular mask"}
    if distance is not None:
        finite = np.isfinite(distance)
        ids = positive_ids[finite]; dist = distance[finite]; weights = area[ids]
        if len(dist) and weights.sum() > 0:
            order = np.argsort(dist); cumulative = np.cumsum(weights[order]) / weights.sum()
            result["geodesic_distance_area_quantiles_a"] = {
                str(q): float(dist[order[np.searchsorted(cumulative, q, side="left")]]) for q in (.5, .9, .95)}
        result["unreachable_positive_samples"] = int((~finite).sum())
    return result
