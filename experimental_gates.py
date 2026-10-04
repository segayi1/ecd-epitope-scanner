"""Validation for immutable, versioned Flow cytometry gate trees."""

from __future__ import annotations

import csv
import io
import json
import math
import re
from typing import Any

import numpy as np


class GateError(ValueError):
    pass


_ID = re.compile(r"^[A-Za-z0-9_. -]{1,120}$")
_KINDS = {"RANGE_1D", "RECTANGLE_2D", "POLYGON_2D"}
_COMPENSATION = {"NONE", "FCS_EMBEDDED", "FLOWJO_WORKSPACE", "EXTERNAL_MATRIX"}
_TRANSFORMS = {
    "ASINH": ("param_t", "param_m", "param_a"),
    "HYPERLOG": ("param_t", "param_w", "param_m", "param_a"),
    "LINEAR": ("param_t", "param_a"),
    "LOG": ("param_t", "param_m"),
    "LOGICLE": ("param_t", "param_w", "param_m", "param_a"),
    "WSP_BIEXP": ("negative", "width", "positive", "max_value"),
    "WSP_LOG": ("offset", "decades"),
}


def parse_compensation_matrix(data: bytes, filename: str) -> dict[str, Any]:
    """Parse a labelled square CSV or JSON matrix into a bounded canonical record."""
    if len(data) > 1024 * 1024:
        raise GateError("Compensation matrix upload exceeds 1 MiB")
    safe_filename = str(filename).replace("\\", "/").split("/")[-1].strip()
    if not safe_filename or len(safe_filename) > 240:
        raise GateError("Compensation matrix filename is invalid")
    suffix = safe_filename.rsplit(".", 1)[-1].lower() if "." in safe_filename else ""
    try:
        content = data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise GateError("Compensation matrix must be UTF-8 CSV or JSON") from exc
    if suffix == "json":
        try:
            document = json.loads(content, parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"Invalid JSON number {value}")
            ))
        except (json.JSONDecodeError, ValueError) as exc:
            raise GateError(f"Invalid compensation matrix JSON: {exc}") from exc
        if not isinstance(document, dict):
            raise GateError("Compensation matrix JSON must be an object")
        if document.get("schema_version", "1.0") != "1.0":
            raise GateError("Unsupported compensation matrix schema_version")
        detectors = document.get("detectors")
        values = document.get("matrix")
        fluorochromes = document.get("fluorochromes")
        import_format = "JSON"
    elif suffix == "csv":
        try:
            rows = list(csv.reader(io.StringIO(content, newline="")))
        except csv.Error as exc:
            raise GateError(f"Invalid compensation matrix CSV: {exc}") from exc
        rows = [[cell.strip() for cell in row] for row in rows]
        while rows and not any(rows[-1]):
            rows.pop()
        if len(rows) < 2 or len(rows[0]) < 2 or rows[0][0].lower() not in {
            "", "detector", "channel", "parameter"
        }:
            raise GateError(
                "Labelled CSV must have a blank/Detector first header, detector columns, and detector row labels"
            )
        detectors = rows[0][1:]
        row_labels = [row[0] for row in rows[1:]]
        try:
            values = [[float(cell) for cell in row[1:]] for row in rows[1:]]
        except ValueError as exc:
            raise GateError("Compensation matrix CSV contains a non-numeric cell") from exc
        if row_labels != detectors:
            raise GateError("Compensation matrix CSV row labels must match column labels in the same order")
        if any(len(row) != len(detectors) + 1 for row in rows[1:]):
            raise GateError("Compensation matrix CSV must be square with no extra columns")
        fluorochromes = None
        import_format = "LABELLED_CSV"
    else:
        raise GateError("Compensation matrix must use .csv or .json")

    if not isinstance(detectors, list) or not 1 <= len(detectors) <= 128:
        raise GateError("Compensation matrix must contain 1 to 128 detector labels")
    if any(not isinstance(item, str) or not item.strip() or len(item.strip()) > 120 for item in detectors):
        raise GateError("Compensation matrix detector labels must be non-empty strings up to 120 characters")
    detectors = [item.strip() for item in detectors]
    if len(set(detectors)) != len(detectors):
        raise GateError("Compensation matrix detector labels must be unique")
    if (
        not isinstance(values, list)
        or len(values) != len(detectors)
        or any(not isinstance(row, list) or len(row) != len(detectors) for row in values)
    ):
        raise GateError("Compensation matrix values must form a square matrix matching detector labels")
    matrix: list[list[float]] = []
    for row in values:
        normalized_row = []
        for value in row:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise GateError("Compensation matrix values must be finite numbers")
            number = float(value)
            if not math.isfinite(number) or abs(number) > 1_000_000:
                raise GateError("Compensation matrix values must be finite and within ±1,000,000")
            normalized_row.append(number)
        matrix.append(normalized_row)
    if fluorochromes is None:
        fluorochromes = [""] * len(detectors)
    if (
        not isinstance(fluorochromes, list)
        or len(fluorochromes) != len(detectors)
        or any(not isinstance(item, str) or len(item) > 120 for item in fluorochromes)
    ):
        raise GateError("fluorochromes must contain one text label per detector")
    try:
        condition = float(np.linalg.cond(np.asarray(matrix, dtype=float)))
    except np.linalg.LinAlgError as exc:
        raise GateError("Compensation matrix is singular") from exc
    if not math.isfinite(condition) or condition > 1e12:
        raise GateError("Compensation matrix is singular or numerically unstable (condition number > 1e12)")
    return {
        "schema_version": "1.0",
        "source_filename": safe_filename,
        "import_format": import_format,
        "detectors": detectors,
        "fluorochromes": fluorochromes,
        "matrix": matrix,
    }


def _normalize_transforms(raw: Any) -> list[dict[str, Any]]:
    if not isinstance(raw, list):
        raise GateError("transforms must be a list of versioned definitions")
    normalized = []
    seen: set[str] = set()
    for transform in raw:
        if not isinstance(transform, dict):
            raise GateError("Each transform definition must be an object")
        transform_id = transform.get("transform_id")
        kind = transform.get("kind")
        parameters = transform.get("parameters")
        if not isinstance(transform_id, str) or not _ID.fullmatch(transform_id):
            raise GateError("Invalid transform_id")
        if transform_id in seen:
            raise GateError(f"Duplicate transform_id: {transform_id}")
        if not isinstance(kind, str) or kind not in _TRANSFORMS or not isinstance(parameters, dict):
            raise GateError(f"Unsupported transform kind or parameters for {transform_id}")
        required = set(_TRANSFORMS[kind])
        if set(parameters) != required:
            raise GateError(f"Transform {transform_id} requires parameters {sorted(required)}")
        values = {}
        for key, value in parameters.items():
            if not isinstance(value, (int, float)) or not math.isfinite(value):
                raise GateError(f"Transform {transform_id} has a non-finite parameter")
            values[key] = float(value)
        if kind in {"ASINH", "HYPERLOG", "LINEAR", "LOG", "LOGICLE"}:
            if values["param_t"] <= 0 or (
                "param_m" in values and values["param_m"] <= 0
            ):
                raise GateError(f"Transform {transform_id} has invalid scale parameters")
            if "param_w" in values and values["param_w"] < 0:
                raise GateError(f"Transform {transform_id} has negative width")
            if "param_a" in values and values["param_a"] < 0:
                raise GateError(f"Transform {transform_id} has negative additional scale")
        if kind == "WSP_BIEXP" and (
            values["positive"] <= 0 or values["max_value"] <= 0
        ):
            raise GateError(f"Transform {transform_id} requires positive and max_value > 0")
        if kind == "WSP_LOG" and any(value <= 0 for value in values.values()):
            raise GateError(f"Transform {transform_id} requires positive parameters")
        normalized.append({"transform_id": transform_id, "kind": kind, "parameters": values})
        seen.add(transform_id)
    return normalized


def _orientation(a: list[float], b: list[float], c: list[float]) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _on_segment(a: list[float], b: list[float], point: list[float]) -> bool:
    return (
        min(a[0], b[0]) <= point[0] <= max(a[0], b[0])
        and min(a[1], b[1]) <= point[1] <= max(a[1], b[1])
    )


def _segments_intersect(
    a: list[float], b: list[float], c: list[float], d: list[float]
) -> bool:
    ab_c = _orientation(a, b, c)
    ab_d = _orientation(a, b, d)
    cd_a = _orientation(c, d, a)
    cd_b = _orientation(c, d, b)
    if (ab_c > 0) != (ab_d > 0) and (cd_a > 0) != (cd_b > 0):
        return True
    return (
        (ab_c == 0 and _on_segment(a, b, c))
        or (ab_d == 0 and _on_segment(a, b, d))
        or (cd_a == 0 and _on_segment(c, d, a))
        or (cd_b == 0 and _on_segment(c, d, b))
    )


def _validate_polygon(vertices: list[list[float]], gate_id: str) -> None:
    if len(set(map(tuple, vertices))) != len(vertices):
        raise GateError(f"Polygon gate {gate_id} repeats a vertex")
    area_twice = sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))
    )
    if math.isclose(area_twice, 0.0, abs_tol=1e-12):
        raise GateError(f"Polygon gate {gate_id} has zero area")
    edge_count = len(vertices)
    for first in range(edge_count):
        a, b = vertices[first], vertices[(first + 1) % edge_count]
        for second in range(first + 1, edge_count):
            if second == first + 1 or (first == 0 and second == edge_count - 1):
                continue
            c, d = vertices[second], vertices[(second + 1) % edge_count]
            if _segments_intersect(a, b, c, d):
                raise GateError(f"Polygon gate {gate_id} self-intersects")


def validate_gate_tree(
    gates: list[dict[str, Any]],
    *,
    available_channels: set[str] | None = None,
    compensation_decision: str = "UNKNOWN",
    compensation_artifact_id: str | None = None,
    transforms: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    if compensation_decision not in _COMPENSATION:
        raise GateError(
            "UNKNOWN_COMPENSATION_SOURCE: select the compensation source before saving gates"
        )
    if compensation_decision == "EXTERNAL_MATRIX" and not compensation_artifact_id:
        raise GateError("External compensation requires a registered matrix artifact")
    if compensation_decision != "EXTERNAL_MATRIX" and compensation_artifact_id:
        raise GateError("A compensation artifact is only valid with EXTERNAL_MATRIX")
    normalized_transforms = _normalize_transforms(transforms or [])
    transform_ids = {item["transform_id"] for item in normalized_transforms}
    if not gates:
        raise GateError("Gate tree must contain at least one gate")
    by_id: dict[str, dict[str, Any]] = {}
    for raw in gates:
        gate_id = str(raw.get("gate_id", ""))
        parent = str(raw.get("parent_gate_id") or "root")
        name = str(raw.get("name", gate_id)).strip()
        kind = raw.get("kind")
        channels = raw.get("channels")
        params = raw.get("parameters")
        if not _ID.fullmatch(gate_id) or gate_id == "root" or gate_id in by_id:
            raise GateError(f"Invalid or duplicate gate_id: {gate_id!r}")
        if not _ID.fullmatch(parent):
            raise GateError(f"Invalid parent gate id for {gate_id}")
        if not name or not _ID.fullmatch(name) or name in {".", ".."}:
            raise GateError(f"Invalid gate name for {gate_id}")
        if (
            not isinstance(kind, str) or kind not in _KINDS
            or not isinstance(channels, list) or not isinstance(params, dict)
        ):
            raise GateError(f"Invalid gate schema for {gate_id}")
        expected_dims = 1 if kind == "RANGE_1D" else 2
        if len(channels) != expected_dims:
            raise GateError(f"{kind} requires {expected_dims} channels")
        if any(not isinstance(channel, str) or not channel for channel in channels):
            raise GateError(f"Gate {gate_id} has an empty channel")
        if len(set(channels)) != expected_dims:
            raise GateError(f"{kind} requires {expected_dims} distinct channels")
        if available_channels is not None and not set(channels).issubset(available_channels):
            raise GateError(f"Gate {gate_id} references a channel absent from uploaded FCS")
        normalized = {
            "gate_id": gate_id,
            "parent_gate_id": parent,
            "name": name,
            "kind": kind,
            "channels": channels,
            "parameters": params,
            "transform_id": raw.get("transform_id"),
            "compensation_id": raw.get("compensation_id"),
        }
        transform_id = raw.get("transform_id")
        if transform_id is not None:
            if not isinstance(transform_id, str) or not _ID.fullmatch(transform_id):
                raise GateError(f"Gate {gate_id} has an invalid transform_id")
            if transform_id not in transform_ids:
                raise GateError(f"Gate {gate_id} references an undefined transform")
        transform_ids_by_channel = raw.get("transform_ids", {})
        if not isinstance(transform_ids_by_channel, dict):
            raise GateError(f"Gate {gate_id} transform_ids must map channels to transform IDs")
        if transform_id is not None and transform_ids_by_channel:
            raise GateError(f"Gate {gate_id} cannot mix transform_id and transform_ids")
        for channel, channel_transform_id in transform_ids_by_channel.items():
            if channel not in channels:
                raise GateError(f"Gate {gate_id} transform channel is not used by the gate")
            if (
                not isinstance(channel_transform_id, str)
                or not _ID.fullmatch(channel_transform_id)
                or channel_transform_id not in transform_ids
            ):
                raise GateError(f"Gate {gate_id} references an undefined channel transform")
        normalized["transform_ids"] = transform_ids_by_channel
        compensation_id = raw.get("compensation_id")
        if compensation_id is not None and (
            not isinstance(compensation_id, str)
            or compensation_id not in {"FCS", "external", "uncompensated"}
        ):
            raise GateError(f"Unsupported gate compensation reference: {compensation_id!r}")
        if compensation_decision == "NONE" and compensation_id == "FCS":
            raise GateError("Gate compensation reference conflicts with compensation_decision NONE")
        if compensation_decision != "EXTERNAL_MATRIX" and compensation_id == "external":
            raise GateError("Gate external compensation reference conflicts with compensation_decision")
        if compensation_decision == "EXTERNAL_MATRIX" and compensation_id == "FCS":
            raise GateError("Gate FCS compensation reference conflicts with EXTERNAL_MATRIX")
        if kind in {"RANGE_1D", "RECTANGLE_2D"}:
            bounds = params.get("bounds") if kind == "RECTANGLE_2D" else [params]
            if not isinstance(bounds, list) or len(bounds) != expected_dims:
                raise GateError(f"{kind} requires bounds for every dimension")
            normalized_bounds = []
            for item in bounds:
                if (
                    not isinstance(item, dict)
                    or item.get("min") is None
                    and item.get("max") is None
                ):
                    raise GateError(f"Gate {gate_id} has an unbounded dimension")
                lower = item.get("min")
                upper = item.get("max")
                for value in (lower, upper):
                    if value is not None and (
                        not isinstance(value, (int, float)) or not math.isfinite(value)
                    ):
                        raise GateError(f"Gate {gate_id} has a non-finite boundary")
                if lower is not None and upper is not None and lower >= upper:
                    raise GateError(f"Gate {gate_id} has min >= max")
                normalized_bounds.append({"min": lower, "max": upper})
            normalized["parameters"] = (
                {"bounds": normalized_bounds} if kind == "RECTANGLE_2D" else normalized_bounds[0]
            )
        elif kind == "POLYGON_2D":
            vertices = params.get("vertices")
            if (
                not isinstance(vertices, list)
                or len(vertices) < 3
                or any(
                    not isinstance(point, (list, tuple))
                    or len(point) != 2
                    or any(
                        not isinstance(value, (int, float)) or not math.isfinite(value)
                        for value in point
                    )
                    for point in vertices
                )
            ):
                raise GateError(f"Polygon gate {gate_id} needs at least three finite 2D vertices")
            normalized["parameters"] = {
                "vertices": [[float(v) for v in point] for point in vertices]
            }
            _validate_polygon(normalized["parameters"]["vertices"], gate_id)
        by_id[gate_id] = normalized
    for gate_id, gate in by_id.items():
        parent = gate["parent_gate_id"]
        if parent != "root" and parent not in by_id:
            raise GateError(f"Unknown parent gate {parent} for {gate_id}")
        seen = {gate_id}
        while parent != "root":
            if parent in seen:
                raise GateError("Gate tree contains a parent cycle")
            seen.add(parent)
            parent = by_id[parent]["parent_gate_id"]
    ordered: list[dict[str, Any]] = []
    remaining = dict(by_id)
    while remaining:
        ready = [
            key
            for key, gate in remaining.items()
            if gate["parent_gate_id"] == "root"
            or any(item["gate_id"] == gate["parent_gate_id"] for item in ordered)
        ]
        if not ready:
            raise GateError("Gate tree could not be ordered")
        for key in ready:
            ordered.append(remaining.pop(key))
    return {
        "gates": ordered,
        "transforms": normalized_transforms,
        "compensation_decision": compensation_decision,
        "compensation_artifact_id": compensation_artifact_id,
    }
