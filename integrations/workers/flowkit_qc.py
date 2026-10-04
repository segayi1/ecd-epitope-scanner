"""Event-level acquisition QC that does not remove or resample events."""

from __future__ import annotations

import math
from statistics import median
from typing import Any

DEFAULT_QC_POLICY = {
    "minimum_gate_events": 500,
    "time_window_count": 20,
    "time_mad_multiplier": 5.0,
    "boundary_fraction_threshold": 0.01,
}


def calculate_sample_qc(
    events: Any,
    channel_names: list[str],
    metadata: dict[str, Any],
    gate_counts: list[dict[str, Any]],
    *,
    policy: dict[str, Any] | None = None,
    events_removed: int = 0,
) -> dict[str, Any]:
    """Flag event-rate outliers, detector rails and low-count gates for review."""
    selected = {**DEFAULT_QC_POLICY, **(policy or {})}
    event_count = int(events.shape[0])
    findings: list[dict[str, Any]] = []
    windows: list[dict[str, Any]] = []

    time_indices = [index for index, name in enumerate(channel_names) if name.casefold() == "time"]
    if not time_indices:
        findings.append({
            "code": "TIME_QC_UNAVAILABLE", "severity": "INFO",
            "message": "No Time channel; acquisition-rate stability was not assessed",
        })
    elif event_count < selected["time_window_count"] * 20:
        findings.append({
            "code": "TIME_QC_INSUFFICIENT_EVENTS", "severity": "INFO",
            "message": "Too few events for the configured event-rate windows",
        })
    else:
        time_index = time_indices[0]
        window_size = math.ceil(event_count / selected["time_window_count"])
        rates = []
        for start in range(0, event_count, window_size):
            end = min(event_count, start + window_size)
            values = events[start:end, time_index]
            first, last = float(values[0]), float(values[-1])
            span = last - first
            monotonic = all(float(values[index]) >= float(values[index - 1])
                            for index in range(1, len(values)))
            rate = (end - start - 1) / span if span > 0 and monotonic else None
            windows.append({
                "event_start": start, "event_end": end, "time_span": span,
                "time_monotonic": monotonic, "event_rate": rate,
            })
            rates.append(rate)
        valid_rates = [rate for rate in rates if rate is not None and math.isfinite(rate)]
        if len(valid_rates) != len(rates):
            findings.append({
                "code": "TIME_AXIS_NONINCREASING", "severity": "REVIEW",
                "message": "At least one acquisition window has a nonincreasing Time axis",
            })
        if len(valid_rates) >= 3:
            center = float(median(valid_rates))
            mad = float(median([abs(rate - center) for rate in valid_rates]))
            robust_scale = 1.4826 * mad
            threshold = selected["time_mad_multiplier"] * robust_scale
            outliers = []
            for window in windows:
                rate = window["event_rate"]
                is_outlier = rate is not None and abs(rate - center) > threshold
                window["rate_outlier"] = is_outlier
                if is_outlier:
                    outliers.append({"event_start": window["event_start"], "event_end": window["event_end"]})
            if outliers:
                findings.append({
                    "code": "TIME_EVENT_RATE_OUTLIER", "severity": "REVIEW",
                    "message": "Event rate differs from the sample median by the configured MAD threshold",
                    "windows": outliers,
                })
        for window in windows:
            window.setdefault("rate_outlier", False)

    boundaries = []
    for index, channel in enumerate(channel_names):
        if channel.casefold() in {"event", "time"} or "flag" in channel.casefold():
            continue
        values = events[:, index]
        nonfinite_count = sum(not math.isfinite(float(value)) for value in values)
        if nonfinite_count:
            findings.append({
                "code": "NONFINITE_CHANNEL_VALUES", "severity": "REVIEW",
                "channel": channel, "count": nonfinite_count,
                "message": "Channel contains non-finite event values",
            })
        try:
            detector_range = float(metadata[f"p{index + 1}r"])
        except (KeyError, TypeError, ValueError):
            continue
        if not math.isfinite(detector_range) or detector_range <= 1:
            continue
        lower_fraction = sum(value <= 0 for value in values) / event_count if event_count else 0.0
        upper_fraction = sum(value >= detector_range - 1 for value in values) / event_count if event_count else 0.0
        boundaries.append({
            "channel": channel, "detector_range": detector_range,
            "fraction_at_lower_rail": lower_fraction,
            "fraction_at_upper_rail": upper_fraction,
        })
        if lower_fraction > selected["boundary_fraction_threshold"]:
            findings.append({
                "code": "DETECTOR_LOWER_BOUNDARY_EVENTS", "severity": "INFO",
                "channel": channel, "fraction_at_boundary": lower_fraction,
                "message": "Events accumulate at the detector's lower range",
            })
        if upper_fraction > selected["boundary_fraction_threshold"]:
            findings.append({
                "code": "DETECTOR_UPPER_SATURATION", "severity": "REVIEW",
                "channel": channel, "fraction_at_boundary": upper_fraction,
                "message": "Events accumulate at the detector's upper range",
            })

    for gate in gate_counts:
        if gate["parent_gate"] != "root" and gate["count"] < selected["minimum_gate_events"]:
            findings.append({
                "code": "LOW_GATE_EVENT_COUNT", "severity": "REVIEW",
                "gate_id": gate.get("gate_id", gate.get("gate_name")), "count": gate["count"],
                "threshold": selected["minimum_gate_events"],
                "message": "Gate count is below the configured minimum for review",
            })
    return {
        "status": "REVIEW_REQUIRED" if any(item["severity"] == "REVIEW" for item in findings)
        else "NO_AUTOMATED_FLAGS",
        "policy": selected,
        "findings": findings,
        "time_windows": windows,
        "detector_boundaries": boundaries,
        "events_examined": event_count,
        "events_removed": events_removed,
    }
