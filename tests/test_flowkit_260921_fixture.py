"""Optional live fixture: user-provided FlowJo 10 Myc experiment, never copied to repo."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FLOWKIT_PYTHON = ROOT / ".tool_envs" / "flowkit" / "Scripts" / "python.exe"
FIXTURE = Path(os.environ.get("FACS_260921_DIR", r"D:\MFL lab\FACS\260921"))


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file() or not FIXTURE.is_dir(),
    reason="local FCS/FlowKit fixture unavailable",
)
def test_real_flowjo_myc_cached_count_mismatches_are_informational_with_consistent_method():
    files = sorted(FIXTURE.glob("*.fcs"))
    assert len(files) == 4
    request = {
        "schema_version": "1.0",
        "request_id": "260921-myc",
        "operation": "apply_wsp_gates",
        "path": str(FIXTURE / "21-Sep-2026.wsp"),
        "fcs_paths": [str(path) for path in files],
    }
    process = subprocess.run(
        [str(FLOWKIT_PYTHON), str(ROOT / "integrations/workers/flowkit_worker.py")],
        input=json.dumps(request),
        capture_output=True,
        text=True,
        timeout=45,
        check=True,
    )
    response = json.loads(process.stdout)
    assert response["status"] == "COMPLETE", response["errors"]
    records = response["records"]
    assert records["status"] == "COMPLETE"
    assert len(records["gate_counts"]) == 8
    assert records["unresolved_nodes"] == []
    assert records["gate_method_consistency"]["status"] == "CONSISTENT"
    assert records["gate_method_consistency"]["samples_checked"] == 4
    assert records["gate_method_consistency"]["method_variants"] == 1
    assert records["gate_count_verification"] == {
        "status": "MISMATCH",
        "compared_gate_count": 8,
        "matched_gate_count": 0,
        "mismatched_gate_count": 8,
        "unmatched_gate_count": 0,
    }
    assert all(
        sample["qc"]["events_examined"] == sample["event_count"]
        and sample["qc"]["events_removed"] == 0
        for sample in records["samples"]
    )
    diagnostics = {
        (row["sample_id"], row["gate_name"]): row
        for row in records["gate_count_diagnostics"]
    }
    assert [
        diagnostics[(f"H{i}.fcs", "Live")]["flowkit_recomputed_count"]
        for i in range(1, 5)
    ] == [3686, 12616, 21737, 20219]
    assert [
        diagnostics[(f"H{i}.fcs", "Live")]["flowjo_cached_count"]
        for i in range(1, 5)
    ] == [3700, 12668, 21915, 20387]
    positive_gate = "Alexa Fluor    647-A, SSC-A subset"
    assert [
        diagnostics[(f"H{i}.fcs", positive_gate)]["flowkit_recomputed_count"]
        for i in range(1, 5)
    ] == [3582, 12059, 21530, 237]
    assert [
        diagnostics[(f"H{i}.fcs", positive_gate)]["flowjo_cached_count"]
        for i in range(1, 5)
    ] == [3595, 12110, 21709, 246]
    assert any(
        row["gate_name"] == "Live" and row["source_gate_type"] == "EllipsoidGate"
        for row in records["gate_count_diagnostics"]
    )


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file() or not FIXTURE.is_dir(),
    reason="local FCS/FlowKit fixture unavailable",
)
def test_internal_gate_revision_executes_parent_child_on_all_events():
    fcs = next(path for path in FIXTURE.glob("*H1.fcs"))
    request = {
        "schema_version": "1.0",
        "request_id": "260921-internal-gates",
        "operation": "apply_gate_revision",
        "paths": [str(fcs)],
        "gate_revision": {
            "compensation_decision": "NONE",
            "gate_tree": [
                {
                    "gate_id": "live",
                    "parent_gate_id": "root",
                    "name": "Live",
                    "kind": "POLYGON_2D",
                    "channels": ["FSC-A", "SSC-A"],
                    "parameters": {"vertices": [
                        [-1e9, -1e9], [1e9, -1e9], [1e9, 1e9], [-1e9, 1e9]
                    ]},
                },
                {
                    "gate_id": "positive",
                    "parent_gate_id": "live",
                    "name": "Positive",
                    "kind": "RANGE_1D",
                    "channels": ["RL1-A"],
                    "parameters": {"min": -1e9, "max": 1e9},
                },
            ],
        },
    }
    process = subprocess.run(
        [str(FLOWKIT_PYTHON), str(ROOT / "integrations/workers/flowkit_worker.py")],
        input=json.dumps(request), capture_output=True, text=True, timeout=45, check=True,
    )
    response = json.loads(process.stdout)
    assert response["status"] == "COMPLETE", response["errors"]
    sample = response["records"]["samples"][0]
    assert sample["event_count"] == 13863
    assert [row["count"] for row in sample["gate_counts"]] == [13863, 13863]
    assert [row["parent_fraction"] for row in sample["gate_counts"]] == [1.0, 1.0]
    assert sample["qc"]["events_examined"] == 13863
    assert sample["qc"]["events_removed"] == 0
    assert len(sample["qc"]["time_windows"]) == 20


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file() or not FIXTURE.is_dir(),
    reason="local FCS/FlowKit fixture unavailable",
)
def test_fcs_embedded_compensation_excludes_nonfluorescent_channels():
    fcs = next(path for path in FIXTURE.glob("*H1.fcs"))
    request = {
        "schema_version": "1.0", "request_id": "260921-fcs-compensation",
        "operation": "apply_gate_revision", "paths": [str(fcs)],
        "gate_revision": {
            "compensation_decision": "FCS_EMBEDDED",
            "gate_tree": [{
                "gate_id": "cells", "parent_gate_id": "root", "name": "Cells",
                "kind": "RECTANGLE_2D", "channels": ["FSC-A", "SSC-A"],
                "parameters": {"bounds": [
                    {"min": 0, "max": 500000}, {"min": 0, "max": 500000}
                ]},
            }],
        },
    }
    process = subprocess.run(
        [str(FLOWKIT_PYTHON), str(ROOT / "integrations/workers/flowkit_worker.py")],
        input=json.dumps(request), capture_output=True, text=True, timeout=45, check=True,
    )
    response = json.loads(process.stdout)
    assert response["status"] == "COMPLETE", response["errors"]
    sample = response["records"]["samples"][0]
    assert sample["event_count"] == 13863
    assert sample["gate_counts"][0]["count"] == 7024


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file() or not FIXTURE.is_dir(),
    reason="local FCS/FlowKit fixture unavailable",
)
def test_external_compensation_matrix_is_applied_by_flowkit_worker():
    fcs = next(path for path in FIXTURE.glob("*H1.fcs"))
    gate = {
        "gate_id": "window", "parent_gate_id": "root", "name": "RL1 window",
        "kind": "RANGE_1D", "channels": ["RL1-A"],
        "parameters": {"min": 50, "max": 500},
    }
    request = {
        "schema_version": "1.0", "request_id": "260921-external-compensation",
        "operation": "apply_gate_revision", "paths": [str(fcs)],
        "gate_revision": {"compensation_decision": "EXTERNAL_MATRIX", "gate_tree": [gate]},
        # This synthetic spill matrix checks execution only; its values are not assay evidence.
        "compensation_matrix": {
            "schema_version": "1.0", "detectors": ["RL1-A", "YL2-A"],
            "fluorochromes": ["Alexa Fluor 647-A", "mCherry-A"],
            "matrix": [[1.0, 0.12], [0.08, 1.0]],
        },
    }
    process = subprocess.run(
        [str(FLOWKIT_PYTHON), str(ROOT / "integrations/workers/flowkit_worker.py")],
        input=json.dumps(request), capture_output=True, text=True, timeout=45, check=True,
    )
    response = json.loads(process.stdout)
    assert response["status"] == "COMPLETE", response["errors"]
    sample = response["records"]["samples"][0]
    assert sample["event_count"] == 13863
    assert sample["gate_counts"][0]["count"] == 151

    without_compensation = {**request, "gate_revision": {
        "compensation_decision": "NONE", "gate_tree": [gate]
    }, "compensation_matrix": None}
    plain_process = subprocess.run(
        [str(FLOWKIT_PYTHON), str(ROOT / "integrations/workers/flowkit_worker.py")],
        input=json.dumps(without_compensation), capture_output=True, text=True,
        timeout=45, check=True,
    )
    plain_response = json.loads(plain_process.stdout)
    assert plain_response["status"] == "COMPLETE", plain_response["errors"]
    assert plain_response["records"]["samples"][0]["gate_counts"][0]["count"] == 150


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file() or not FIXTURE.is_dir(),
    reason="local FCS/FlowKit fixture unavailable",
)
def test_versioned_logicle_transform_is_applied_before_gate_boundaries():
    fcs = next(path for path in FIXTURE.glob("*H1.fcs"))
    request = {
        "schema_version": "1.0", "request_id": "260921-logicle-gate",
        "operation": "apply_gate_revision", "paths": [str(fcs)],
        "gate_revision": {
            "compensation_decision": "NONE",
            "transforms": [{
                "transform_id": "fl1-logicle-v1", "kind": "LOGICLE",
                "parameters": {
                    "param_t": 262144, "param_w": 0.5, "param_m": 4.5, "param_a": 0,
                },
            }],
            "gate_tree": [{
                "gate_id": "all", "parent_gate_id": "root", "name": "All transformed",
                "kind": "RECTANGLE_2D", "channels": ["FSC-A", "RL1-A"],
                "transform_ids": {"RL1-A": "fl1-logicle-v1"},
                "parameters": {"bounds": [
                    {"min": -1e9, "max": 1e9}, {"min": -100, "max": 100}
                ]},
            }],
        },
    }
    process = subprocess.run(
        [str(FLOWKIT_PYTHON), str(ROOT / "integrations/workers/flowkit_worker.py")],
        input=json.dumps(request), capture_output=True, text=True, timeout=45, check=True,
    )
    response = json.loads(process.stdout)
    assert response["status"] == "COMPLETE", response["errors"]
    sample = response["records"]["samples"][0]
    assert sample["event_count"] == 13863
    assert sample["gate_counts"][0]["count"] == 13863


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file() or not FIXTURE.is_dir(),
    reason="local FCS/FlowKit fixture unavailable",
)
def test_real_fcs_scatter_sampling_is_repeatable_and_read_only():
    fcs = next(path for path in FIXTURE.glob("*H1.fcs"))
    original_sha256 = hashlib.sha256(fcs.read_bytes()).hexdigest()
    request = {
        "schema_version": "1.0", "request_id": "260921-sample-plot",
        "operation": "sample_events", "path": str(fcs),
        "x_channel": "FSC-A", "y_channel": "SSC-A", "maximum_points": 257,
    }
    responses = []
    for _ in range(2):
        process = subprocess.run(
            [str(FLOWKIT_PYTHON), str(ROOT / "integrations/workers/flowkit_worker.py")],
            input=json.dumps(request), capture_output=True, text=True, timeout=45, check=True,
        )
        response = json.loads(process.stdout)
        assert response["status"] == "COMPLETE", response["errors"]
        responses.append(response["records"])

    first, second = responses
    assert first["display_only"] is True
    assert first["source"] == "RAW_UNTRANSFORMED"
    assert first["sampling_method"] == "DETERMINISTIC_EVEN_INDEX"
    assert first["event_count"] == 13863
    assert first["finite_event_count"] == 13863
    assert first["sampled_event_count"] == 257
    assert len(first["points"]) == 257
    assert all(len(point) == 2 and all(map(lambda value: isinstance(value, (float, int)), point))
               for point in first["points"])
    assert first["domain"]["quantiles"] == [0.01, 0.99]
    assert first["points"] == second["points"]
    assert first["domain"] == second["domain"]
    assert hashlib.sha256(fcs.read_bytes()).hexdigest() == original_sha256


@pytest.mark.skipif(
    not FLOWKIT_PYTHON.is_file() or not FIXTURE.is_dir(),
    reason="local FCS/FlowKit fixture unavailable",
)
def test_manual_time_window_exclusion_recounts_without_mutating_real_fcs():
    fcs = next(path for path in FIXTURE.glob("*H1.fcs"))
    original_sha256 = hashlib.sha256(fcs.read_bytes()).hexdigest()
    request = {
        "schema_version": "1.0", "request_id": "260921-manual-time-window-exclusion",
        "operation": "apply_gate_revision", "paths": [str(fcs)],
        "gate_revision": {
            "compensation_decision": "NONE",
            "gate_tree": [{
                "gate_id": "all", "parent_gate_id": "root", "name": "All events",
                "kind": "RANGE_1D", "channels": ["FSC-A"],
                "parameters": {"min": -1e12, "max": 1e12},
            }],
        },
        "qc_policy": {"minimum_gate_events": 0, "time_window_count": 20,
                      "time_mad_multiplier": 5.0, "boundary_fraction_threshold": 0.01},
    }

    def run_worker(payload):
        process = subprocess.run(
            [str(FLOWKIT_PYTHON), str(ROOT / "integrations/workers/flowkit_worker.py")],
            input=json.dumps(payload), capture_output=True, text=True, timeout=45, check=True,
        )
        response = json.loads(process.stdout)
        assert response["status"] == "COMPLETE", response["errors"]
        return response["records"]["samples"][0]

    original = run_worker(request)
    windows = original["qc"]["time_windows"]
    if not windows:
        pytest.skip("H1 fixture has no usable Time-channel windows")
    selected = windows[0]
    selected_range = {"event_start": selected["event_start"], "event_end": selected["event_end"]}
    removed = selected_range["event_end"] - selected_range["event_start"]
    assert 0 < removed < original["event_count"]

    excluded = run_worker({**request, "event_exclusion_ranges_by_path": [[selected_range]]})
    assert excluded["source_event_count"] == original["event_count"]
    assert excluded["events_removed"] == removed
    assert excluded["event_count"] == original["event_count"] - removed
    assert excluded["event_exclusion_ranges"] == [selected_range]
    assert excluded["gate_counts"][0]["count"] == excluded["event_count"]
    assert excluded["gate_counts"][0]["parent_count"] == excluded["event_count"]
    assert excluded["qc"]["events_examined"] == excluded["event_count"]
    assert excluded["qc"]["events_removed"] == removed
    assert hashlib.sha256(fcs.read_bytes()).hexdigest() == original_sha256
