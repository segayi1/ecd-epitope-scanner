"""Isolated FlowKit 1.3.2 worker; JSON in and JSON out.

The parent service resolves artifact IDs to trusted local paths before invoking
this worker. Paths in this process protocol must never be accepted from HTTP.
"""

from __future__ import annotations

import json
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import flowkit
import numpy as np
from flowkit_qc import DEFAULT_QC_POLICY, calculate_sample_qc


def inspect_fcs(path: Path) -> dict[str, Any]:
    sample = flowkit.Sample(str(path))
    metadata = sample.get_metadata()
    selected = {
        key: value
        for key, value in metadata.items()
        if key.lower()
        in {
            "date",
            "btim",
            "etim",
            "cyt",
            "cytsn",
            "tot",
            "par",
            "spillover",
            "spill",
            "spil",
            "fil",
            "timestep",
        }
    }
    compensation_keys = [key for key in metadata if key.lower() in {"spillover", "spill", "spil"}]
    return {
        "file_name": path.name,
        "sample_id": sample.id,
        "event_count": int(sample.event_count),
        "channels": [
            {"number": index + 1, "pnn": pnn, "pns": pns}
            for index, (pnn, pns) in enumerate(zip(sample.pnn_labels, sample.pns_labels))
        ],
        "metadata": selected,
        "compensation_source": "FCS_MATRIX_AVAILABLE"
        if compensation_keys
        else "UNKNOWN_COMPENSATION_SOURCE",
        "compensation_keys": compensation_keys,
    }


def sample_events(
    path: Path, x_channel: str, y_channel: str, maximum_points: int = 3000
) -> dict[str, Any]:
    """Return a deterministic display-only subset; analysis always uses full events."""
    if maximum_points < 100 or maximum_points > 10_000:
        raise ValueError("maximum_points must be between 100 and 10000")
    if x_channel == y_channel:
        raise ValueError("Scatter axes must use different channels")
    sample = flowkit.Sample(str(path))
    pnn = [str(item) for item in sample.pnn_labels]
    pns = [str(item) if item is not None else "" for item in sample.pns_labels]

    def resolve(channel: str) -> int:
        direct = [index for index, label in enumerate(pnn) if label == channel]
        secondary = [index for index, label in enumerate(pns) if label == channel]
        matches = direct or secondary
        if len(matches) != 1:
            raise ValueError(f"Channel {channel!r} is missing or ambiguous")
        return matches[0]

    x_index, y_index = resolve(x_channel), resolve(y_channel)
    events = np.asarray(sample.get_events(source="raw", subsample=False), dtype=float)
    xy = events[:, [x_index, y_index]]
    finite = np.isfinite(xy).all(axis=1)
    valid = xy[finite]
    if valid.size == 0:
        raise ValueError("Selected channels have no finite paired events")
    count = min(maximum_points, len(valid))
    indices = np.linspace(0, len(valid) - 1, count, dtype=int)
    selected = valid[indices]
    quantiles = np.quantile(valid, [0.01, 0.99], axis=0)
    return {
        "file_name": path.name,
        "sample_id": sample.id,
        "event_count": int(sample.event_count),
        "finite_event_count": int(len(valid)),
        "omitted_nonfinite_count": int(len(xy) - len(valid)),
        "sampled_event_count": int(len(selected)),
        "maximum_points": maximum_points,
        "sampling_method": "DETERMINISTIC_EVEN_INDEX",
        "display_only": True,
        "source": "RAW_UNTRANSFORMED",
        "x_channel": x_channel,
        "y_channel": y_channel,
        "domain": {
            "x": [float(quantiles[0, 0]), float(quantiles[1, 0])],
            "y": [float(quantiles[0, 1]), float(quantiles[1, 1])],
            "quantiles": [0.01, 0.99],
        },
        "points": selected.tolist(),
    }


WSP_SAMPLE_KEYWORDS = ("fil", "wellid", "platename", "smno")
WSP_GATING_NAMESPACE = "http://www.isac-net.org/std/Gating-ML/v2.0/gating"
WSP_DATA_TYPE_NAMESPACE = "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes"
WSP_TRANSFORM_NAMESPACE = "http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
WSP_TRANSFORM_REQUIRED_ATTRIBUTES = {
    "linear": ("minRange", "maxRange"),
    "log": ("offset", "decades"),
    "logicle": ("T", "W", "M", "A"),
    "biex": ("neg", "width", "length", "maxRange", "pos"),
    "fasinh": ("T", "A", "M"),
}
WSP_CONVERTIBLE_GATE_TYPES = {
    "RectangleGate", "PolygonGate", "EllipsoidGate", "BooleanGate",
}


def _local_xml_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _wsp_cached_gate_counts(path: Path) -> list[dict[str, Any]]:
    """Read FlowJo's saved counts and original gate types for replay verification."""
    root = ET.parse(path).getroot()
    records = []
    for sample_el in root.iter():
        if _local_xml_name(sample_el.tag) != "Sample":
            continue
        sample_node = next(
            (child for child in list(sample_el) if _local_xml_name(child.tag) == "SampleNode"),
            None,
        )
        if sample_node is None:
            continue
        keywords = {
            str(child.attrib.get("name", "")).lstrip("$").casefold():
            str(child.attrib.get("value", ""))
            for keyword_list in list(sample_el)
            if _local_xml_name(keyword_list.tag) == "Keywords"
            for child in list(keyword_list)
            if _local_xml_name(child.tag) == "Keyword"
        }
        counts: dict[tuple[tuple[str, ...], str], dict[str, Any]] = {}

        def collect_populations(container: ET.Element, parent_path: tuple[str, ...]) -> None:
            for population in list(container):
                if _local_xml_name(population.tag) != "Population":
                    continue
                gate_name = str(population.attrib.get("name", ""))
                gate_element = next(
                    (child for child in list(population) if _local_xml_name(child.tag) == "Gate"),
                    None,
                )
                gate_element_type = None
                if gate_element is not None:
                    gate_child = next(iter(list(gate_element)), None)
                    if gate_child is not None:
                        gate_element_type = _local_xml_name(gate_child.tag)
                count_value = population.attrib.get("count")
                if gate_name and count_value is not None:
                    try:
                        cached_count = int(float(count_value))
                    except (TypeError, ValueError):
                        cached_count = None
                    if cached_count is not None:
                        counts[(parent_path, gate_name)] = {
                            "count": cached_count,
                            "source_gate_type": gate_element_type,
                        }
                gate_path = (*parent_path, gate_name)
                for child in list(population):
                    if _local_xml_name(child.tag) == "Subpopulations":
                        collect_populations(child, gate_path)

        for child in list(sample_node):
            if _local_xml_name(child.tag) == "Subpopulations":
                collect_populations(child, ("root",))
        workspace_sample_id = f"__wsp_sample_{sample_node.attrib.get('sampleID', '')}_{len(records)}"
        records.append(
            {
                "sample_name": str(sample_node.attrib.get("name", "")),
                "sample_id": str(sample_node.attrib.get("sampleID", "")),
                "workspace_sample_id": workspace_sample_id,
                "keywords": keywords,
                "counts": counts,
            }
        )
    return records


def _wsp_transformed_ellipsoid_nodes(
    path: Path, cached_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Find WSP ellipsoids FlowKit 1.3.2 cannot parse due to dimension transforms."""
    root = ET.parse(path).getroot()
    nodes = []
    sample_index = 0

    for sample_el in root.iter():
        if _local_xml_name(sample_el.tag) != "Sample":
            continue
        sample_node = next(
            (child for child in list(sample_el) if _local_xml_name(child.tag) == "SampleNode"),
            None,
        )
        if sample_node is None:
            continue
        cached_record = (
            cached_records[sample_index] if sample_index < len(cached_records) else None
        )
        sample_index += 1
        if cached_record is None:
            continue

        def collect_populations(container: ET.Element, parent_path: tuple[str, ...]) -> None:
            for population in list(container):
                if _local_xml_name(population.tag) != "Population":
                    continue
                gate_name = str(population.attrib.get("name", ""))
                gate_element = next(
                    (child for child in list(population) if _local_xml_name(child.tag) == "Gate"),
                    None,
                )
                gate_child = next(iter(list(gate_element)), None) if gate_element is not None else None
                if gate_child is not None and _local_xml_name(gate_child.tag) == "EllipsoidGate":
                    has_transform = False
                    channels = []
                    for dimension in gate_child.iter():
                        if _local_xml_name(dimension.tag) != "dimension":
                            continue
                        if f"{{{WSP_GATING_NAMESPACE}}}transformation-ref" not in dimension.attrib:
                            continue
                        has_transform = True
                        for dimension_child in dimension.iter():
                            if _local_xml_name(dimension_child.tag) not in {
                                "fcs-dimension", "new-dimension"
                            }:
                                continue
                            channel = dimension_child.attrib.get(
                                f"{{{WSP_DATA_TYPE_NAMESPACE}}}name", ""
                            ).strip()
                            if channel:
                                channels.append(channel)
                    if has_transform:
                        nodes.append({
                            "workspace_sample_id": cached_record["workspace_sample_id"],
                            "name": gate_name,
                            "parent_path": list(parent_path),
                            "source_gate_type": "EllipsoidGate",
                            "kind": None,
                            "channels": channels,
                            "status": "UNRESOLVED_FLOWJO_NODE",
                            "reason_code": "FLOWJO_ELLIPSOID_TRANSFORM_UNSUPPORTED",
                        })
                gate_path = (*parent_path, gate_name) if gate_name else parent_path
                for child in list(population):
                    if _local_xml_name(child.tag) == "Subpopulations":
                        collect_populations(child, gate_path)

        for child in list(sample_node):
            if _local_xml_name(child.tag) == "Subpopulations":
                collect_populations(child, ("root",))
    return nodes


def _wsp_unsupported_gate_nodes(
    path: Path, cached_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Find WSP population gates FlowKit cannot convert into executable gates."""
    root = ET.parse(path).getroot()
    nodes = []
    sample_index = 0

    for sample_el in root.iter():
        if _local_xml_name(sample_el.tag) != "Sample":
            continue
        sample_node = next(
            (child for child in list(sample_el) if _local_xml_name(child.tag) == "SampleNode"),
            None,
        )
        if sample_node is None:
            continue
        cached_record = (
            cached_records[sample_index] if sample_index < len(cached_records) else None
        )
        sample_index += 1
        if cached_record is None:
            continue

        def collect(container: ET.Element, parent_path: tuple[str, ...]) -> None:
            for child in list(container):
                child_type = _local_xml_name(child.tag)
                if child_type == "Population":
                    gate_name = str(child.attrib.get("name", ""))
                    gate_el = next(
                        (item for item in list(child) if _local_xml_name(item.tag) == "Gate"),
                        None,
                    )
                    gate_children = list(gate_el) if gate_el is not None else []
                    gate_type = _local_xml_name(gate_children[0].tag) if len(gate_children) == 1 else ""
                    if gate_type not in WSP_CONVERTIBLE_GATE_TYPES:
                        channels = []
                        for dimension in gate_el.iter() if gate_el is not None else ():
                            if _local_xml_name(dimension.tag) != "dimension":
                                continue
                            for dimension_child in dimension.iter():
                                if _local_xml_name(dimension_child.tag) not in {
                                    "fcs-dimension", "new-dimension"
                                }:
                                    continue
                                channel = dimension_child.attrib.get(
                                    f"{{{WSP_DATA_TYPE_NAMESPACE}}}name", ""
                                ).strip()
                                if channel:
                                    channels.append(channel)
                        nodes.append({
                            "workspace_sample_id": cached_record["workspace_sample_id"],
                            "name": gate_name or "Unnamed gate",
                            "parent_path": list(parent_path),
                            "source_gate_type": gate_type or "INVALID_GATE_STRUCTURE",
                            "kind": None,
                            "channels": channels,
                            "status": "UNRESOLVED_FLOWJO_NODE",
                            "reason_code": (
                                "FLOWJO_GATE_TYPE_UNSUPPORTED"
                                if gate_type else "FLOWJO_GATE_STRUCTURE_INVALID"
                            ),
                        })
                    gate_path = (*parent_path, gate_name) if gate_name else parent_path
                    for nested in list(child):
                        if _local_xml_name(nested.tag) == "Subpopulations":
                            collect(nested, gate_path)
                elif child_type in {"AndNode", "OrNode", "NotNode"}:
                    gate_name = str(child.attrib.get("name", ""))
                    gate_path = (*parent_path, gate_name) if gate_name else parent_path
                    for nested in list(child):
                        if _local_xml_name(nested.tag) == "Subpopulations":
                            collect(nested, gate_path)
                elif child_type == "Subpopulations":
                    collect(child, parent_path)

        for child in list(sample_node):
            if _local_xml_name(child.tag) == "Subpopulations":
                collect(child, ("root",))
    return nodes


def _wsp_transform_nodes(
    path: Path, cached_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Find FlowJo WSP transforms unsupported or malformed for FlowKit 1.3.2."""
    root = ET.parse(path).getroot()
    nodes = []
    sample_index = 0

    for sample_el in root.iter():
        if _local_xml_name(sample_el.tag) != "Sample":
            continue
        sample_node = next(
            (child for child in list(sample_el) if _local_xml_name(child.tag) == "SampleNode"),
            None,
        )
        if sample_node is None:
            continue
        cached_record = (
            cached_records[sample_index] if sample_index < len(cached_records) else None
        )
        sample_index += 1
        if cached_record is None:
            continue
        transforms_el = next(
            (child for child in list(sample_el) if _local_xml_name(child.tag) == "Transformations"),
            None,
        )
        if transforms_el is None:
            nodes.append({
                "workspace_sample_id": cached_record["workspace_sample_id"],
                "name": "Channel transformations",
                "parent_path": ["root"],
                "source_gate_type": "FlowJoTransform",
                "kind": None,
                "channels": [],
                "transform_type": "UNSPECIFIED",
                "status": "UNRESOLVED_FLOWJO_NODE",
                "reason_code": "FLOWJO_TRANSFORM_BLOCK_MISSING",
            })
            continue
        for transform in list(transforms_el):
            transform_type = _local_xml_name(transform.tag).casefold()
            parameter = next(
                (item for item in transform.iter() if _local_xml_name(item.tag) == "parameter"),
                None,
            )
            channel = (
                parameter.attrib.get(f"{{{WSP_DATA_TYPE_NAMESPACE}}}name", "").strip()
                if parameter is not None else ""
            )
            required_attributes = WSP_TRANSFORM_REQUIRED_ATTRIBUTES.get(transform_type, ())
            attributes = {
                name: str(
                    transform.attrib.get(f"{{{WSP_TRANSFORM_NAMESPACE}}}{name}", "") or ""
                ).strip()
                for name in required_attributes
            }
            reason_code = None
            if transform_type not in WSP_TRANSFORM_REQUIRED_ATTRIBUTES:
                reason_code = "FLOWJO_TRANSFORM_UNSUPPORTED"
            elif not channel:
                reason_code = "FLOWJO_TRANSFORM_PARAMETERS_INVALID"
            else:
                required = WSP_TRANSFORM_REQUIRED_ATTRIBUTES[transform_type]
                if any(not attributes.get(name) for name in required):
                    reason_code = "FLOWJO_TRANSFORM_PARAMETERS_INVALID"
                elif transform_type == "biex" and attributes["length"] != "256":
                    reason_code = "FLOWJO_TRANSFORM_PARAMETER_UNSUPPORTED"
                else:
                    numeric_attributes = [
                        name for name in required
                        if not (transform_type == "biex" and name == "length")
                    ]
                    try:
                        if any(
                            not bool(np.isfinite(float(attributes[name])))
                            for name in numeric_attributes
                        ):
                            reason_code = "FLOWJO_TRANSFORM_PARAMETERS_INVALID"
                    except (TypeError, ValueError):
                        reason_code = "FLOWJO_TRANSFORM_PARAMETERS_INVALID"
            if reason_code is None:
                continue
            nodes.append({
                "workspace_sample_id": cached_record["workspace_sample_id"],
                "name": f"Channel transform: {channel or 'unknown channel'}",
                "parent_path": ["root"],
                "source_gate_type": "FlowJoTransform",
                "kind": None,
                "channels": [channel] if channel else [],
                "transform_type": transform_type,
                "status": "UNRESOLVED_FLOWJO_NODE",
                "reason_code": reason_code,
            })
    return nodes


def _flowjo_unsupported_gate_result(
    path: Path,
    cached_records: list[dict[str, Any]],
    input_sample_ids: dict[str, str],
) -> dict[str, Any] | None:
    nodes = (
        _wsp_transformed_ellipsoid_nodes(path, cached_records)
        + _wsp_unsupported_gate_nodes(path, cached_records)
        + _wsp_transform_nodes(path, cached_records)
    )
    if not nodes:
        return None
    for node in nodes:
        workspace_sample_id = node.pop("workspace_sample_id")
        node["sample_id"] = input_sample_ids.get(workspace_sample_id, workspace_sample_id)
    return {
        "file_name": path.name,
        "groups": [],
        "samples": [],
        "unresolved_nodes": nodes,
        "gate_counts": [],
        "status": "WAITING_FOR_USER",
    }


def _match_wsp_sample(metadata: dict[str, Any], records: list[dict[str, Any]]) -> dict[str, Any] | None:
    values = {
        key: str(metadata.get(key, "")).strip().casefold()
        for key in WSP_SAMPLE_KEYWORDS
    }
    values = {key: "" if value in {"", "na", "none"} else value for key, value in values.items()}
    if not values["fil"]:
        return None
    matches = []
    for record in records:
        cached = record["keywords"]
        shared_keys = [
            key for key in WSP_SAMPLE_KEYWORDS
            if values[key] and str(cached.get(key, "")).strip().casefold() not in {"", "na", "none"}
        ]
        if "fil" not in shared_keys:
            continue
        if all(values[key] == str(cached.get(key, "")).strip().casefold() for key in shared_keys):
            matches.append(record)
    return matches[0] if len(matches) == 1 else None


def _stable_method_value(value: Any) -> Any:
    """Make FlowKit transform/compensation settings comparable across samples."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, np.generic):
        return _stable_method_value(value.item())
    if isinstance(value, np.ndarray):
        return _stable_method_value(value.tolist())
    if isinstance(value, dict):
        return {
            str(key): _stable_method_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_stable_method_value(item) for item in value]
    attributes = getattr(value, "__dict__", None)
    if isinstance(attributes, dict):
        return {
            str(key): _stable_method_value(item)
            for key, item in sorted(attributes.items())
            if not str(key).startswith("_")
        } or {"type": type(value).__name__}
    return {"type": type(value).__name__}


def _gate_membership_semantics(gate: Any) -> dict[str, Any]:
    """Return gate settings that change which events belong to a population."""
    semantics: dict[str, Any] = {}
    if hasattr(gate, "use_complement"):
        semantics["use_complement"] = bool(gate.use_complement)
    if type(gate).__name__ == "BooleanGate":
        semantics["boolean_operation"] = str(getattr(gate, "type", "")).casefold()
        semantics["gate_references"] = _stable_method_value(
            getattr(gate, "gate_refs", [])
        )
    return semantics


def _verify_wsp_gate_method(
    workspace: flowkit.Workspace, cached_records: list[dict[str, Any]],
    input_sample_ids: dict[str, str],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Check that a single imported batch uses one gate strategy and processing setup.

    Per-sample gate coordinates are intentionally excluded: manual boundaries may
    differ while the gate hierarchy, gate types, channels, membership semantics,
    transforms, and compensation method remain consistent within the batch.
    """
    signatures: dict[str, list[str]] = {}
    unmatched = []
    for sample_id in workspace.get_sample_ids():
        sample = workspace.get_sample(sample_id)
        record = _match_wsp_sample(sample.get_metadata(), cached_records)
        if record is None:
            unmatched.append({
                "sample_id": input_sample_ids.get(sample_id, sample_id),
                "name": "WSP sample mapping",
                "parent_path": ["root"],
                "source_gate_type": "WSP_SAMPLE",
                "kind": None,
                "status": "WSP_SAMPLE_UNMATCHED",
                "reason_code": "WSP_SAMPLE_UNMATCHED",
            })
            continue

        gates = []
        for gate_name, gate_path in workspace.get_gate_ids(sample_id):
            gate_path_tuple = tuple(str(item) for item in gate_path)
            cached_gate = record["counts"].get((gate_path_tuple, str(gate_name)))
            gate = workspace.get_gate(sample_id, gate_name, gate_path)
            gates.append({
                "name": str(gate_name),
                "parent_path": list(gate_path_tuple),
                "source_gate_type": cached_gate.get("source_gate_type") if cached_gate else None,
                "flowkit_gate_type": type(gate).__name__,
                "channels": [
                    str(getattr(dimension, "id", ""))
                    for dimension in (getattr(gate, "dimensions", None) or [])
                ],
                "membership_semantics": _gate_membership_semantics(gate),
            })
        gate_channels = {channel for gate in gates for channel in gate["channels"]}
        transforms = {
            str(channel): {
                "type": type(transform).__name__,
                "parameters": _stable_method_value(transform),
            }
            for channel, transform in sorted(workspace.get_transforms(sample_id).items())
            if str(channel) in gate_channels
        }
        compensation = workspace.get_comp_matrix(sample_id)
        signature = {
            "gates": sorted(gates, key=lambda item: (item["parent_path"], item["name"])),
            "transforms": transforms,
            "wsp_compensation": _stable_method_value(compensation),
        }
        signature_key = json.dumps(signature, sort_keys=True, separators=(",", ":"))
        signatures.setdefault(signature_key, []).append(sample_id)

    status = "INCOMPLETE" if unmatched else "CONSISTENT" if len(signatures) <= 1 else "INCONSISTENT"
    unresolved = list(unmatched)
    if len(signatures) > 1:
        unresolved.append({
            "name": "Within-batch gate method",
            "parent_path": ["root"],
            "source_gate_type": "BATCH_METHOD_COMPARISON",
            "kind": None,
            "status": "FLOWJO_GATE_METHOD_INCONSISTENT",
            "reason_code": "FLOWJO_GATE_METHOD_INCONSISTENT",
            "method_variants": len(signatures),
            "samples_checked": sum(len(sample_ids) for sample_ids in signatures.values()),
            "examples": [
                {"sample_ids": sample_ids[:3]}
                for sample_ids in [
                    [input_sample_ids.get(item, item) for item in values]
                    for values in list(signatures.values())[:3]
                ]
            ],
        })
    return ({
        "status": status,
        "samples_checked": sum(len(sample_ids) for sample_ids in signatures.values()),
        "method_variants": len(signatures),
        "comparison_scope": [
            "gate_hierarchy", "source_gate_type", "flowkit_gate_type", "channels",
            "gate_membership_semantics", "gate_channel_transforms", "wsp_compensation",
        ],
        "gate_coordinates_compared": False,
        "unmatched_sample_count": len(unmatched),
    }, unresolved)


def _verify_wsp_counts(
    workspace: flowkit.Workspace, report: Any, cached_records: list[dict[str, Any]],
    input_sample_ids: dict[str, str],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    cached_total = sum(len(record["counts"]) for record in cached_records)
    if cached_total == 0:
        return ({"status": "UNAVAILABLE", "compared_gate_count": 0,
                 "matched_gate_count": 0, "mismatched_gate_count": 0,
                 "unmatched_gate_count": 0}, [], [])

    matched_samples: dict[str, dict[str, Any] | None] = {}
    for sample_id in workspace.get_sample_ids():
        matched_samples[sample_id] = _match_wsp_sample(
            workspace.get_sample(sample_id).get_metadata(), cached_records
        )

    diagnostics = []
    unmatched = []
    mismatches = []
    for row in report.itertuples(index=False):
        sample_id = str(row.sample_id)
        parent_path = tuple(str(item) for item in row.gate_path)
        gate_name = str(row.gate_name)
        record = matched_samples.get(sample_id)
        cached_gate = record["counts"].get((parent_path, gate_name)) if record else None
        if cached_gate is None:
            unmatched.append({
                "sample_id": input_sample_ids.get(sample_id, sample_id),
                "name": gate_name,
                "parent_path": list(parent_path),
                "source_gate_type": "WSP_CACHED_COUNT",
                "kind": None,
                "status": "WSP_CACHED_COUNT_UNMATCHED",
            })
            continue
        cached_count = int(cached_gate["count"])
        recomputed_count = int(row.count)
        diagnostic = {
            "sample_id": input_sample_ids.get(sample_id, sample_id),
            "gate_name": gate_name,
            "gate_path": list(parent_path),
            "source_gate_type": cached_gate["source_gate_type"],
            "flowjo_cached_count": cached_count,
            "flowkit_recomputed_count": recomputed_count,
            "difference": recomputed_count - cached_count,
            "matches": recomputed_count == cached_count,
        }
        diagnostics.append(diagnostic)
        if not diagnostic["matches"]:
            mismatches.append(diagnostic)

    unresolved = list(unmatched)
    by_gate: dict[tuple[tuple[str, ...], str, str], list[dict[str, Any]]] = {}
    for item in mismatches:
        key = (tuple(item["gate_path"]), item["gate_name"], str(item["source_gate_type"]))
        by_gate.setdefault(key, []).append(item)
    for (parent_path, gate_name, source_type), items in by_gate.items():
        unresolved.append({
            "name": gate_name,
            "parent_path": list(parent_path),
            "source_gate_type": source_type,
            "kind": None,
            "status": "FLOWJO_CACHED_COUNT_MISMATCH",
            "reason_code": "FLOWJO_CACHED_COUNT_MISMATCH",
            "mismatched_samples": len(items),
            "samples_checked": len({item["sample_id"] for item in diagnostics
                                     if item["gate_name"] == gate_name
                                     and item["gate_path"] == list(parent_path)}),
            "examples": [
                {
                    "sample_id": item["sample_id"],
                    "flowjo_cached_count": item["flowjo_cached_count"],
                    "flowkit_recomputed_count": item["flowkit_recomputed_count"],
                }
                for item in items[:3]
            ],
        })

    verification_status = "MISMATCH" if mismatches else "MATCHED"
    if unmatched:
        verification_status = "INCOMPLETE"
    return ({
        "status": verification_status,
        "compared_gate_count": len(diagnostics),
        "matched_gate_count": len(diagnostics) - len(mismatches),
        "mismatched_gate_count": len(mismatches),
        "unmatched_gate_count": len(unmatched),
    }, diagnostics, unresolved)


def _write_rekeyed_wsp(path: Path, destination: Path) -> None:
    """Give duplicate FlowJo sample names unique internal keys in a temp copy."""
    tree = ET.parse(path)
    sample_index = 0
    for sample_el in tree.getroot().iter():
        if _local_xml_name(sample_el.tag) != "Sample":
            continue
        sample_node = next(
            (child for child in list(sample_el) if _local_xml_name(child.tag) == "SampleNode"),
            None,
        )
        if sample_node is None:
            continue
        sample_id = str(sample_node.attrib.get("sampleID", ""))
        sample_node.set("name", f"__wsp_sample_{sample_id}_{sample_index}")
        sample_index += 1
    tree.write(destination, encoding="utf-8", xml_declaration=True)


def import_wsp(
    path: Path, fcs_paths: list[Path], *, apply_gates: bool = False,
    qc_policy: dict[str, Any] | None = None,
) -> dict[str, Any]:
    # Preserve input IDs for output, but hold Sample objects so they can be
    # paired to unique WSP records even when several batches reuse B3.fcs.
    sample_paths: dict[str, list[Path]] = {}
    loaded_samples: list[tuple[Path, flowkit.Sample, str]] = []
    for fcs_path in fcs_paths:
        sample = flowkit.Sample(str(fcs_path))
        sample_id = str(sample.id)
        sample_paths.setdefault(str(sample_id), []).append(fcs_path)
        loaded_samples.append((fcs_path, sample, sample_id))
    duplicate_ids = {
        sample_id: paths for sample_id, paths in sample_paths.items() if len(paths) > 1
    }
    if duplicate_ids:
        examples = ", ".join(list(duplicate_ids)[:8])
        raise ValueError(
            "FCS inputs contain duplicate FlowJo sample IDs, which FlowKit cannot "
            "distinguish in one workspace import. Import one sample group at a time. "
            f"Found {len(duplicate_ids)} duplicated IDs (examples: {examples})."
        )

    cached_wsp_counts = _wsp_cached_gate_counts(path)
    input_sample_ids: dict[str, str] = {}
    workspace_ids_by_input: dict[str, str] = {}
    matched_wsp_records: dict[str, dict[str, Any]] = {}
    for _, sample, input_sample_id in loaded_samples:
        record = _match_wsp_sample(sample.get_metadata(), cached_wsp_counts)
        if record is None:
            raise ValueError(
                f"FCS sample {input_sample_id!r} does not uniquely match a WSP sample "
                "by $FIL/$WELLID/$PLATENAME/$SMNO metadata"
            )
        workspace_sample_id = str(record["workspace_sample_id"])
        if workspace_sample_id in input_sample_ids:
            raise ValueError("Multiple FCS files matched the same FlowJo workspace sample")
        input_sample_ids[workspace_sample_id] = input_sample_id
        workspace_ids_by_input[input_sample_id] = workspace_sample_id
        matched_wsp_records[workspace_sample_id] = record
        sample.id = workspace_sample_id

    unsupported_gate_result = _flowjo_unsupported_gate_result(
        path, cached_wsp_counts, input_sample_ids
    )
    if unsupported_gate_result is not None:
        return unsupported_gate_result

    with tempfile.TemporaryDirectory(prefix="ecd-flowjo-") as temp_dir:
        temporary_wsp_path = Path(temp_dir) / "workspace.wsp"
        _write_rekeyed_wsp(path, temporary_wsp_path)
        workspace = flowkit.Workspace(
            str(temporary_wsp_path),
            fcs_samples=[sample for _, sample, _ in loaded_samples],
        )
    samples = []
    supported = {
        "PolygonGate": "POLYGON_2D",
        "RectangleGate": "RECTANGLE_2D",
        "BooleanGate": "BOOLEAN",
    }
    unresolved = []
    for sample_id in workspace.get_sample_ids():
        input_sample_id = input_sample_ids[sample_id]
        matched_record = matched_wsp_records[sample_id]
        gates = []
        for gate_name, gate_path in workspace.get_gate_ids(sample_id):
            gate = workspace.get_gate(sample_id, gate_name, gate_path)
            gate_type = type(gate).__name__
            cached_gate = matched_record["counts"].get(
                (tuple(str(item) for item in gate_path), str(gate_name))
            )
            source_gate_type = (
                cached_gate.get("source_gate_type") if cached_gate else gate_type
            )
            mapped = supported.get(gate_type)
            if gate_type == "RectangleGate" and len(getattr(gate, "dimensions", [])) == 1:
                mapped = "RANGE_1D"
            channels = [
                getattr(dim, "id", None)
                for dim in (getattr(gate, "dimensions", None) or [])
            ]
            record = {
                "name": gate_name,
                "parent_path": list(gate_path),
                "source_gate_type": source_gate_type,
                "flowkit_gate_type": gate_type,
                "kind": mapped,
                "channels": channels,
                "membership_semantics": _gate_membership_semantics(gate),
                "status": "RECOGNIZED" if mapped else "UNRESOLVED_FLOWJO_NODE",
            }
            if source_gate_type == "EllipsoidGate" and gate_type == "PolygonGate":
                record["mapping_note"] = "FLOWKIT_APPROXIMATES_FLOWJO_ELLIPSOID_AS_POLYGON"
            gates.append(record)
            if mapped is None:
                unresolved.append({"sample_id": input_sample_id, **record})
        matrix = workspace.get_comp_matrix(sample_id)
        samples.append(
            {
                "sample_id": input_sample_id,
                "gates": gates,
                "wsp_compensation": "AVAILABLE" if matrix is not None else "NONE",
                "transforms": {
                    channel: type(transform).__name__
                    for channel, transform in workspace.get_transforms(sample_id).items()
                },
            }
        )
    result: dict[str, Any] = {
        "file_name": path.name,
        "groups": workspace.get_sample_groups(),
        "samples": samples,
        "unresolved_nodes": unresolved,
    }
    if apply_gates:
        if unresolved:
            result["gate_counts"] = []
            result["status"] = "WAITING_FOR_USER"
        else:
            # Disable multiprocessing: Windows spawn may recursively import a
            # worker launched from a pipe, and event counts need no pool.
            workspace.analyze_samples(use_mp=False)
            report = workspace.get_analysis_report()
            recomputed_gate_counts = [
                {
                    "sample_id": input_sample_ids.get(str(row.sample_id), str(row.sample_id)),
                    "gate_name": str(row.gate_name),
                    "gate_path": list(row.gate_path),
                    "parent_gate": str(row.parent),
                    "count": int(row.count),
                    "absolute_fraction": float(row.absolute_percent) / 100,
                    "parent_fraction": float(row.relative_percent) / 100,
                }
                for row in report.itertuples(index=False)
            ]
            verification, diagnostics, cache_unresolved = _verify_wsp_counts(
                workspace, report, cached_wsp_counts, input_sample_ids
            )
            result["gate_count_verification"] = verification
            result["gate_count_diagnostics"] = diagnostics
            method_consistency, method_unresolved = _verify_wsp_gate_method(
                workspace, cached_wsp_counts, input_sample_ids
            )
            result["gate_method_consistency"] = method_consistency
            if method_consistency["status"] != "CONSISTENT":
                result["unresolved_nodes"].extend(method_unresolved)
            if verification["status"] == "INCOMPLETE":
                result["unresolved_nodes"].extend(cache_unresolved)
            if method_consistency["status"] != "CONSISTENT" or verification["status"] == "INCOMPLETE":
                result["gate_counts"] = []
                result["status"] = "WAITING_FOR_USER"
            else:
                result["gate_counts"] = recomputed_gate_counts
                result["status"] = "COMPLETE"
            for sample_record in samples:
                sample_id = sample_record["sample_id"]
                flow_sample = workspace.get_sample(workspace_ids_by_input[sample_id])
                sample_counts = [
                    row for row in result["gate_counts"] if row["sample_id"] == sample_id
                ]
                sample_record["event_count"] = int(flow_sample.event_count)
                sample_record["qc"] = calculate_sample_qc(
                    flow_sample.get_events(source="raw", subsample=False),
                    list(flow_sample.pnn_labels), flow_sample.get_metadata(),
                    sample_counts, policy=qc_policy,
                )
            result["qc_policy"] = {**DEFAULT_QC_POLICY, **(qc_policy or {})}
    return result


def apply_gate_revision(
    path: Path, gate_tree: list[dict[str, Any]], compensation_decision: str,
    qc_policy: dict[str, Any] | None = None,
    transforms: list[dict[str, Any]] | None = None,
    compensation_matrix: dict[str, Any] | None = None,
    event_exclusion_ranges: list[dict[str, int]] | None = None,
) -> dict[str, Any]:
    """Execute the project's explicit rectangular/polygon gate tree on all events."""
    if compensation_decision not in {"NONE", "FCS_EMBEDDED", "EXTERNAL_MATRIX"}:
        raise ValueError(
            "This gate executor supports NONE, FCS_EMBEDDED, or EXTERNAL_MATRIX compensation; "
            "replay FlowJo gates for workspace compensation"
        )
    external_matrix = None
    if compensation_decision == "EXTERNAL_MATRIX":
        if not isinstance(compensation_matrix, dict):
            raise ValueError("EXTERNAL_MATRIX was selected but the registered matrix document is missing")
        detectors = compensation_matrix.get("detectors")
        values = compensation_matrix.get("matrix")
        fluorochromes = compensation_matrix.get("fluorochromes")
        if (
            compensation_matrix.get("schema_version") != "1.0"
            or not isinstance(detectors, list)
            or not isinstance(values, list)
            or not isinstance(fluorochromes, list)
        ):
            raise ValueError("External compensation matrix document is invalid")
        external_matrix = flowkit.Matrix(
            np.asarray(values, dtype=float), detectors, fluorochromes
        )
    elif compensation_matrix is not None:
        raise ValueError("A compensation matrix was supplied without EXTERNAL_MATRIX selection")
    if compensation_decision == "FCS_EMBEDDED":
        metadata_sample = flowkit.Sample(str(path))
        metadata = metadata_sample.get_metadata()
        spill_value = next(
            (value for key, value in metadata.items()
             if key.lower() in {"spillover", "spill", "spil"}),
            None,
        )
        if not isinstance(spill_value, str):
            raise ValueError("FCS_EMBEDDED was selected but the FCS file has no compensation matrix")
        spill_fields = [item.strip() for item in spill_value.split(",")]
        try:
            spill_count = int(spill_fields[0])
            compensated_channels = set(spill_fields[1 : spill_count + 1])
        except (ValueError, IndexError) as exc:
            raise ValueError("FCS compensation matrix channel list is invalid") from exc
        if spill_count < 1 or len(compensated_channels) != spill_count:
            raise ValueError("FCS compensation matrix has invalid or duplicate channel labels")
        null_channels = [
            channel for channel in metadata_sample.pnn_labels
            if channel not in compensated_channels
        ]
        sample = flowkit.Sample(str(path), null_channel_list=null_channels)
    else:
        sample = flowkit.Sample(str(path))
    strategy = flowkit.GatingStrategy()
    if external_matrix is not None:
        strategy.add_comp_matrix("external", external_matrix)
    transform_types = {
        "ASINH": flowkit.transforms.AsinhTransform,
        "HYPERLOG": flowkit.transforms.HyperlogTransform,
        "LINEAR": flowkit.transforms.LinearTransform,
        "LOG": flowkit.transforms.LogTransform,
        "LOGICLE": flowkit.transforms.LogicleTransform,
        "WSP_BIEXP": flowkit.transforms.WSPBiexTransform,
        "WSP_LOG": flowkit.transforms.WSPLogTransform,
    }
    transform_ids = set()
    for definition in transforms or []:
        transform_id = definition["transform_id"]
        kind = definition["kind"]
        if transform_id in transform_ids or kind not in transform_types:
            raise ValueError(f"Invalid or duplicate transform definition {transform_id!r}")
        strategy.add_transform(
            transform_id, transform_types[kind](**definition["parameters"])
        )
        transform_ids.add(transform_id)
    by_id = {str(gate["gate_id"]): gate for gate in gate_tree}
    for gate in gate_tree:
        channel_transform_ids = gate.get("transform_ids", {})
        default_transform_ref = gate.get("transform_id")
        def transform_ref(channel: str) -> str | None:
            reference = channel_transform_ids.get(channel, default_transform_ref)
            if reference is not None and reference not in transform_ids:
                raise ValueError(f"Gate {gate['gate_id']} references an undefined transform")
            return reference
        channels = gate["channels"]
        compensation_ref = gate.get("compensation_id")
        if compensation_ref is None:
            compensation_ref = (
                "FCS" if compensation_decision == "FCS_EMBEDDED"
                else "external" if compensation_decision == "EXTERNAL_MATRIX"
                else "uncompensated"
            )
        if compensation_ref not in {"FCS", "external", "uncompensated"}:
            raise ValueError(f"Unsupported compensation reference {compensation_ref!r}")
        if compensation_ref == "external" and external_matrix is None:
            raise ValueError("Gate references external compensation but no matrix was selected")
        if compensation_ref == "FCS" and compensation_decision != "FCS_EMBEDDED":
            raise ValueError("Gate references FCS compensation without FCS_EMBEDDED selection")
        parameters = gate["parameters"]
        dimensions = []
        if gate["kind"] == "RANGE_1D":
            bound = parameters
            dimensions.append(
                flowkit.Dimension(
                    channels[0], compensation_ref=compensation_ref,
                    transformation_ref=transform_ref(channels[0]),
                    range_min=bound.get("min"), range_max=bound.get("max"),
                )
            )
            flow_gate = flowkit.gates.RectangleGate(gate["name"], dimensions)
        elif gate["kind"] == "RECTANGLE_2D":
            dimensions = [
                flowkit.Dimension(
                    channel, compensation_ref=compensation_ref,
                    transformation_ref=transform_ref(channel),
                    range_min=bound.get("min"), range_max=bound.get("max"),
                )
                for channel, bound in zip(channels, parameters["bounds"], strict=True)
            ]
            flow_gate = flowkit.gates.RectangleGate(gate["name"], dimensions)
        elif gate["kind"] == "POLYGON_2D":
            dimensions = [
                flowkit.Dimension(
                    channel, compensation_ref=compensation_ref,
                    transformation_ref=transform_ref(channel),
                )
                for channel in channels
            ]
            flow_gate = flowkit.gates.PolygonGate(
                gate["name"], dimensions, np.asarray(parameters["vertices"], dtype=float)
            )
        else:
            raise ValueError(f"Unsupported gate kind {gate['kind']!r}")

        parent_id = gate["parent_gate_id"]
        parent_names = []
        while parent_id != "root":
            parent = by_id[parent_id]
            parent_names.insert(0, parent["name"])
            parent_id = parent["parent_gate_id"]
        strategy.add_gate(flow_gate, tuple(["root", *parent_names]))

    results = strategy.gate_sample(sample, cache_events=True)
    events = np.asarray(sample.get_events(source="raw", subsample=False))
    source_event_count = int(sample.event_count)
    excluded_mask = np.zeros(source_event_count, dtype=bool)
    normalized_ranges = []
    for item in event_exclusion_ranges or []:
        start, end = item.get("event_start"), item.get("event_end")
        if (not isinstance(start, int) or isinstance(start, bool)
                or not isinstance(end, int) or isinstance(end, bool)
                or start < 0 or end <= start or end > source_event_count):
            raise ValueError("Manual event exclusion range is outside the source FCS event indices")
        if excluded_mask[start:end].any():
            raise ValueError("Manual event exclusion ranges must not overlap")
        excluded_mask[start:end] = True
        normalized_ranges.append({"event_start": start, "event_end": end})
    events_removed = int(excluded_mask.sum())
    if events_removed >= source_event_count:
        raise ValueError("Manual event exclusions must leave at least one source event")
    keep_mask = ~excluded_mask
    analyzed_events = events[keep_mask]
    rows = []
    for gate in gate_tree:
        parent_id = gate["parent_gate_id"]
        parent_names = []
        while parent_id != "root":
            parent = by_id[parent_id]
            parent_names.insert(0, parent["name"])
            parent_id = parent["parent_gate_id"]
        gate_path = tuple(["root", *parent_names])
        membership = results.get_gate_membership(gate["name"], gate_path=gate_path)
        count = int(np.count_nonzero(membership & keep_mask))
        parent_membership = (
            results.get_gate_membership(
                parent_names[-1], gate_path=tuple(["root", *parent_names[:-1]])
            ) if parent_names else np.ones(source_event_count, dtype=bool)
        )
        parent_count = int(np.count_nonzero(parent_membership & keep_mask))
        analyzed_event_count = source_event_count - events_removed
        rows.append({
            "gate_id": gate["gate_id"], "gate_name": gate["name"],
            "gate_path": list(gate_path), "parent_gate": parent_names[-1] if parent_names else "root",
            "parent_count": parent_count, "count": count,
            "absolute_fraction": count / analyzed_event_count,
            "parent_fraction": count / parent_count if parent_count else None,
        })
    qc = calculate_sample_qc(
        analyzed_events, list(sample.pnn_labels), sample.get_metadata(), rows,
        policy=qc_policy, events_removed=events_removed,
    )
    if events_removed:
        qc["findings"].append({
            "code": "MANUAL_EVENT_WINDOWS_EXCLUDED", "severity": "INFO",
            "windows": normalized_ranges, "count": events_removed,
            "message": "Reviewer-selected event windows were excluded from this immutable gate analysis",
        })
    return {
        "file_name": path.name, "sample_id": sample.id,
        "event_count": source_event_count - events_removed,
        "source_event_count": source_event_count,
        "events_removed": events_removed,
        "event_exclusion_ranges": normalized_ranges,
        "gate_counts": rows, "qc": qc,
    }


def run(request: dict[str, Any]) -> dict[str, Any]:
    started = time.monotonic()
    operation = request.get("operation")
    if request.get("schema_version") != "1.0":
        raise ValueError("Unsupported request schema_version")
    if operation == "inspect_fcs":
        records = inspect_fcs(Path(request["path"]))
    elif operation == "sample_events":
        records = sample_events(
            Path(request["path"]),
            str(request["x_channel"]),
            str(request["y_channel"]),
            int(request.get("maximum_points", 3000)),
        )
    elif operation in {"import_wsp", "apply_wsp_gates"}:
        records = import_wsp(
            Path(request["path"]),
            [Path(item) for item in request["fcs_paths"]],
            apply_gates=operation == "apply_wsp_gates",
            qc_policy=request.get("qc_policy"),
        )
    elif operation == "apply_gate_revision":
        paths = request.get("paths") or [request["path"]]
        exclusions_by_path = request.get("event_exclusion_ranges_by_path") or [[] for _ in paths]
        if len(exclusions_by_path) != len(paths):
            raise ValueError("Event exclusion selections must match the analyzed FCS files")
        records = {
            "samples": [
                apply_gate_revision(
                    Path(path), request["gate_revision"]["gate_tree"],
                    request["gate_revision"]["compensation_decision"],
                    request.get("qc_policy"),
                    request["gate_revision"].get("transforms", []),
                    request.get("compensation_matrix"),
                    exclusions,
                )
                for path, exclusions in zip(paths, exclusions_by_path, strict=True)
            ]
        }
    else:
        raise ValueError(f"Unsupported operation: {operation}")
    return {
        "schema_version": "1.0",
        "request_id": request.get("request_id"),
        "status": records.get("status", "COMPLETE"),
        "records": records,
        "warnings": [],
        "errors": [],
        "provenance": {"flowkit_version": flowkit.__version__, "worker": "flowkit_worker"},
        "runtime_seconds": round(time.monotonic() - started, 3),
    }


def main() -> None:
    try:
        response = run(json.load(sys.stdin))
    except Exception as exc:
        response = {
            "schema_version": "1.0",
            "status": "FAILED",
            "records": None,
            "warnings": [],
            "errors": [{"code": type(exc).__name__, "message": str(exc)}],
            "provenance": {"flowkit_version": flowkit.__version__},
        }
    print(json.dumps(response, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()
