"""Per-channel empirical resolution; no biological confidence claims."""
from __future__ import annotations

import numpy as np


def calibrate_channel(null_blocks, numerical_error, protocol_id, channel, *, matched=True,
                      required_mask=None, area=None):
    blocks = np.asarray(null_blocks, dtype=float)
    error = np.asarray(numerical_error, dtype=float)
    if blocks.ndim != 2 or error.shape != blocks.shape[1:] or len(blocks) < 4:
        return {"status": "UNRESOLVED", "reason": "insufficient matched null blocks", "channel": channel}
    valid = np.isfinite(blocks).all(axis=0) & np.isfinite(error) & (error >= 0)
    magnitude = np.abs(blocks)
    q95 = np.quantile(magnitude, .95, axis=0)
    q99 = np.quantile(magnitude, .99, axis=0)
    threshold = np.maximum(q95, error)
    # Explicit convergence diagnostic; zeros cannot waive unknown numerical error.
    half = max(2, len(blocks)//2)
    previous = np.quantile(magnitude[:half], .95, axis=0)
    stable = np.abs(previous-q95) <= np.maximum(error, .1 * np.maximum(q95, 1e-12))
    valid &= stable & bool(matched)
    required = np.ones(error.shape, bool) if required_mask is None else np.asarray(required_mask, bool)
    weights = np.ones(error.shape) if area is None else np.asarray(area, float)
    if required.shape != error.shape or weights.shape != error.shape:
        raise ValueError('calibration scope/reference shape mismatch')
    if not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError('invalid calibration area weights')
    denominator = float(weights[required].sum())
    fraction = float(weights[required & valid].sum()) / denominator if denominator > 0 else None
    scope_converged = bool(required.any() and denominator > 0 and valid[required].all())
    return {"status": "COMPLETE" if valid.all() else "PARTIAL", "channel": channel,
            "protocol_id": protocol_id, "threshold": threshold, "q95": q95, "q99": q99,
            "numerical_error": error, "valid": valid, "null_blocks": blocks,
            "requested_scope_converged": scope_converged,
            "required_mask": required, "required_area": denominator,
            "calibrated_local_area_fraction": fraction,
            "unresolved_local_area_fraction": 1. - fraction if fraction is not None else None,
            "block_count": len(blocks), "interpretation": "empirical computational noise; not false-positive control"}


def resolved_changes(deltas, calibration, protocol_id, *, directional=True):
    values = np.asarray(deltas, dtype=float)
    if values.ndim != 2:
        raise ValueError("expected run x canonical sample channel")
    unknown = np.ones(values.shape[1], dtype=bool)
    if calibration.get("status") not in {"COMPLETE", "PARTIAL"} or calibration.get("protocol_id") != protocol_id or len(values) < 5:
        return np.zeros(values.shape[1], bool), unknown
    threshold = np.asarray(calibration["threshold"])
    if threshold.shape != values.shape[1:]:
        raise ValueError("calibration reference shape mismatch")
    values = values[:5]
    finite = np.isfinite(values)
    valid = np.asarray(calibration["valid"])
    known = valid & finite.all(axis=0)
    positive = np.sum(values > threshold, axis=0) >= 4
    negative = np.sum(values < -threshold, axis=0) >= 4
    hits = (positive | negative) if directional else (np.sum(np.abs(values) > threshold, axis=0) >= 4)
    if calibration.get('positive_lower_bound'):
        missing = (~finite).sum(axis=0)
        if directional:
            possible = ((np.sum(values > threshold, axis=0)+missing >= 4)
                        | (np.sum(values < -threshold, axis=0)+missing >= 4))
        else:
            possible = np.sum(np.abs(values)>threshold,axis=0)+missing >= 4
        known = valid & (hits | ~possible)
    return hits & known, ~known
