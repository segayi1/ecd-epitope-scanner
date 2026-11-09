"""Freeze the residue graph used for V4 sectoring and later experiment views."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..util import sha256_file, write_json
from .models import V4Sector


def write_surface_graph_snapshot(
    path: Path,
    graph: dict[int, dict[int, float]],
    areas: dict[int, float],
    sectors: list[V4Sector] | None,
    *,
    mapping_path: Path,
    structure_path: Path | None,
    provenance: dict[str, Any],
    include_display_ownership: bool = True,
) -> str:
    """Write the exact weighted graph, with stable UniProt positions as node IDs."""
    if sectors is None and include_display_ownership:
        raise ValueError("Display ownership requires verified V4 sectors")
    owners = {
        position: sector.sector_id
        for sector in sectors or []
        for position in sector.residue_positions
    }
    if len(owners) != sum(len(sector.residue_positions) for sector in sectors or []):
        raise ValueError("A surface residue belongs to multiple V4 sectors")
    if sectors is not None and set(owners) != set(graph):
        raise ValueError("V4 sectors do not cover the surface graph exactly")
    edges = []
    for left in sorted(graph):
        for right, distance in sorted(graph[left].items()):
            if right not in graph or graph[right].get(left) != distance:
                raise ValueError("Surface graph edge is not symmetric")
            if left < right:
                edges.append({"source": left, "target": right, "distance_a": distance})
    payload = {
        "schema_version": "1.0" if include_display_ownership else "1.1",
        "position_system": "uniprot_pos",
        "mapping_sha256": sha256_file(mapping_path),
        "structure_sha256": (
            sha256_file(structure_path) if structure_path is not None and structure_path.is_file()
            else None
        ),
        "graph_provenance": provenance,
        "nodes": [
            {
                "uniprot_pos": position,
                "surface_area_a2": areas[position],
                **({"sector_id": owners[position]} if include_display_ownership else {}),
            }
            for position in sorted(graph)
        ],
        "edges": edges,
    }
    write_json(path, payload)
    return sha256_file(path)
