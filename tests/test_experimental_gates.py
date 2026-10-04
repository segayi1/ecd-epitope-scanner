from __future__ import annotations

import pytest
from ecd_epitope_scanner.experimental_gates import (
    GateError,
    parse_compensation_matrix,
    validate_gate_tree,
)


def test_gate_tree_orders_parents_and_validates_channels_and_compensation():
    child = {
        "gate_id": "positive",
        "parent_gate_id": "live",
        "kind": "RANGE_1D",
        "channels": ["RL1-A"],
        "parameters": {"min": 100.0, "max": 500.0},
    }
    parent = {
        "gate_id": "live",
        "kind": "POLYGON_2D",
        "channels": ["FSC-A", "SSC-A"],
        "parameters": {"vertices": [[0, 0], [1, 0], [1, 1], [0, 1]]},
    }
    result = validate_gate_tree(
        [child, parent],
        available_channels={"FSC-A", "SSC-A", "RL1-A"},
        compensation_decision="NONE",
    )
    assert [item["gate_id"] for item in result["gates"]] == ["live", "positive"]
    assert result["gates"][0]["parameters"]["vertices"][0] == [0.0, 0.0]
    with pytest.raises(GateError, match="UNKNOWN_COMPENSATION_SOURCE"):
        validate_gate_tree([parent])
    with pytest.raises(GateError, match="absent from uploaded FCS"):
        validate_gate_tree(
            [parent], available_channels={"FSC-A", "RL1-A"}, compensation_decision="NONE"
        )


@pytest.mark.parametrize(
    "bad",
    [
        {"min": 5, "max": 1},
        {"min": None, "max": None},
        {"min": float("nan"), "max": 1},
    ],
)
def test_invalid_range_is_rejected(bad):
    gate = {"gate_id": "g", "kind": "RANGE_1D", "channels": ["A"], "parameters": bad}
    with pytest.raises(GateError):
        validate_gate_tree([gate], compensation_decision="NONE")


def test_unknown_parent_and_cycle_are_rejected():
    a = {
        "gate_id": "a",
        "parent_gate_id": "b",
        "kind": "RANGE_1D",
        "channels": ["A"],
        "parameters": {"min": 0, "max": 1},
    }
    b = {
        "gate_id": "b",
        "parent_gate_id": "a",
        "kind": "RANGE_1D",
        "channels": ["A"],
        "parameters": {"min": 0, "max": 1},
    }
    with pytest.raises(GateError, match="cycle"):
        validate_gate_tree([a, b], compensation_decision="NONE")


def test_polygon_geometry_and_unavailable_processing_metadata_are_rejected():
    gate = {
        "gate_id": "poly", "kind": "POLYGON_2D", "channels": ["FSC-A", "SSC-A"],
        "parameters": {"vertices": [[0, 0], [4, 4], [0, 4], [4, 0]]},
    }
    with pytest.raises(GateError, match="zero area|self-intersects"):
        validate_gate_tree([gate], compensation_decision="NONE")

    gate["parameters"]["vertices"] = [[0, 0], [2, 0], [2, 2], [0, 2], [0, 0]]
    with pytest.raises(GateError, match="repeats a vertex"):
        validate_gate_tree([gate], compensation_decision="NONE")

    gate["parameters"]["vertices"] = [[0, 0], [2, 0], [2, 2], [0, 2]]
    gate["transform_id"] = "logicle"
    with pytest.raises(GateError, match="undefined.*transform"):
        validate_gate_tree([gate], compensation_decision="NONE")


def test_channel_validation_handles_non_string_values_as_client_errors():
    gate = {
        "gate_id": "bad", "kind": "RECTANGLE_2D", "channels": [["FSC-A"], "SSC-A"],
        "parameters": {"bounds": [{"min": 0, "max": 1}, {"min": 0, "max": 1}]},
    }
    with pytest.raises(GateError, match="empty channel"):
        validate_gate_tree([gate], compensation_decision="NONE")


def test_transform_definitions_are_versioned_and_gate_references_are_resolved():
    gate = {
        "gate_id": "positive", "kind": "RECTANGLE_2D", "channels": ["FSC-A", "FL1-A"],
        "transform_ids": {"FL1-A": "logicle-v1"},
        "parameters": {"bounds": [{"min": 0, "max": 100}, {"min": 0.25, "max": 1.0}]},
    }
    transform = {
        "transform_id": "logicle-v1", "kind": "LOGICLE",
        "parameters": {"param_t": 262144, "param_w": 0.5, "param_m": 4.5, "param_a": 0},
    }
    result = validate_gate_tree([gate], compensation_decision="NONE", transforms=[transform])
    assert result["gates"][0]["transform_ids"] == {"FL1-A": "logicle-v1"}
    assert result["transforms"][0]["kind"] == "LOGICLE"
    with pytest.raises(GateError, match="undefined.*transform"):
        validate_gate_tree([gate], compensation_decision="NONE")
    transform["parameters"]["param_t"] = float("inf")
    with pytest.raises(GateError, match="non-finite"):
        validate_gate_tree([gate], compensation_decision="NONE", transforms=[transform])


def test_compensation_matrix_json_and_labelled_csv_are_normalized():
    document = parse_compensation_matrix(
        b'{"schema_version":"1.0","detectors":["RL1-A","YL2-A"],'
        b'"matrix":[[1,0.12],[0.08,1]]}',
        "spill.json",
    )
    assert document["detectors"] == ["RL1-A", "YL2-A"]
    assert document["fluorochromes"] == ["", ""]
    assert document["matrix"] == [[1.0, 0.12], [0.08, 1.0]]
    csv_document = parse_compensation_matrix(
        b",RL1-A,YL2-A\nRL1-A,1,0.12\nYL2-A,0.08,1\n",
        "spill.csv",
    )
    assert csv_document["matrix"] == document["matrix"]
    assert csv_document["import_format"] == "LABELLED_CSV"


@pytest.mark.parametrize(
    ("content", "filename", "message"),
    [
        (b'{"detectors":["A","A"],"matrix":[[1,0],[0,1]]}', "dup.json", "unique"),
        (b'{"detectors":["A","B"],"matrix":[[1,1],[1,1]]}', "singular.json", "unstable"),
        (b',A,B\nA,1,0\nB,0\n', "ragged.csv", "square"),
        (b'{"detectors":["A"],"matrix":[[NaN]]}', "nan.json", "Invalid JSON number"),
    ],
)
def test_invalid_compensation_matrices_are_rejected(content, filename, message):
    with pytest.raises(GateError, match=message):
        parse_compensation_matrix(content, filename)


def test_gate_external_compensation_reference_requires_registered_matrix_id():
    gate = {
        "gate_id": "positive", "kind": "RANGE_1D", "channels": ["RL1-A"],
        "parameters": {"min": 0, "max": 10},
    }
    with pytest.raises(GateError, match="requires a registered matrix artifact"):
        validate_gate_tree([gate], compensation_decision="EXTERNAL_MATRIX")
    result = validate_gate_tree(
        [gate], compensation_decision="EXTERNAL_MATRIX", compensation_artifact_id="matrix-id"
    )
    assert result["compensation_artifact_id"] == "matrix-id"
    gate["compensation_id"] = "external"
    assert validate_gate_tree(
        [gate], compensation_decision="EXTERNAL_MATRIX", compensation_artifact_id="matrix-id"
    )["gates"][0]["compensation_id"] == "external"
    with pytest.raises(GateError, match="conflicts with EXTERNAL_MATRIX"):
        gate["compensation_id"] = "FCS"
        validate_gate_tree(
            [gate], compensation_decision="EXTERNAL_MATRIX", compensation_artifact_id="matrix-id"
        )
