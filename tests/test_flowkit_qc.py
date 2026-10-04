from __future__ import annotations

import numpy as np
from ecd_epitope_scanner.integrations.workers.flowkit_qc import calculate_sample_qc


def test_qc_flags_rate_boundary_and_small_gates_without_removing_events():
    events = np.zeros((1000, 2), dtype=float)
    events[:, 0] = np.arange(1000)
    events[800:, 0] = 799
    events[100, 1] = np.nan
    events[980:, 1] = 99.5
    result = calculate_sample_qc(
        events,
        ["Time", "RL1-A"],
        {"p2r": "100"},
        [{"gate_id": "positive", "parent_gate": "live", "count": 30}],
        policy={
            "minimum_gate_events": 50,
            "time_window_count": 5,
            "time_mad_multiplier": 3.0,
            "boundary_fraction_threshold": 0.01,
        },
    )

    codes = {finding["code"] for finding in result["findings"]}
    assert "TIME_AXIS_NONINCREASING" in codes
    assert "DETECTOR_LOWER_BOUNDARY_EVENTS" in codes
    assert "DETECTOR_UPPER_SATURATION" in codes
    assert "NONFINITE_CHANNEL_VALUES" in codes
    assert "LOW_GATE_EVENT_COUNT" in codes
    assert result["status"] == "REVIEW_REQUIRED"
    assert result["events_examined"] == 1000
    assert result["events_removed"] == 0
    assert len(result["time_windows"]) == 5


def test_normal_time_axis_has_no_event_rate_outlier():
    events = np.zeros((1000, 1), dtype=float)
    events[:, 0] = np.arange(1000)
    result = calculate_sample_qc(
        events, ["Time"], {}, [],
        policy={
            "minimum_gate_events": 10,
            "time_window_count": 5,
            "time_mad_multiplier": 5.0,
            "boundary_fraction_threshold": 0.01,
        },
    )
    assert result["status"] == "NO_AUTOMATED_FLAGS"
    assert not result["findings"]
    assert result["events_removed"] == 0
