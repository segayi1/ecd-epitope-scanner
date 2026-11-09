"""Replay region aggregation from saved geometry, without running design engines.

Usage: python -m ecd_epitope_scanner.surface_patch.replay_regions RUN_DIR OUTPUT_DIR
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from ..config import SURFACE_V2_PARAMETERS
from ..util import write_json
from .legacy_adapter import _region_color
from .models import Footprint, SurfaceGeometry
from .region_builder import build_surface_regions
from .surface_graph import build_surface_graph


def replay(source: Path, output: Path) -> dict:
    if source.resolve() == output.resolve():
        raise ValueError("Replay output must be separate from the original run")
    with (source / "surface_vertices.csv").open(encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    n = len(rows)
    with np.load(source / "surface_edges.npz") as edges:
        surface = SurfaceGeometry(
            vertices=np.array([[float(r[k]) for k in ("x", "y", "z")] for r in rows]),
            normals=np.array([[float(r[k]) for k in ("nx", "ny", "nz")] for r in rows]),
            vertex_area=np.array([float(r["area"]) for r in rows]),
            edges=edges["edges"],
            edge_length=edges["edge_length"],
            owner_residue=np.array([int(r["owner_residue"]) for r in rows]),
            owner_atom=np.full(n, -1),
            chain_id=np.full(n, ""),
            canonical_residue_id=np.array([int(r["owner_residue"]) for r in rows]),
            is_ecd=np.ones(n, dtype=bool),
            confidence=np.full(n, np.nan),
            component_id=np.array([int(r["component_id"]) for r in rows]),
        )
    raw = [
        Footprint(**r)
        for r in json.loads(
            (source / "raw_footprints.json").read_text(encoding="utf-8")
        )
    ]
    for f in raw:
        f.residue_surface_area = {int(k): v for k, v in f.residue_surface_area.items()}
        f.residue_surface_fraction = {
            int(k): v for k, v in f.residue_surface_fraction.items()
        }
    unique = [f for f in raw if f.duplicate_of is None]
    snapshots = []
    graph = build_surface_graph(surface)
    for _ in range(3):
        regions = build_surface_regions(
            unique,
            surface=surface,
            graph=graph,
            primary_candidates=raw,
            parameters=SURFACE_V2_PARAMETERS["regions"],
        )
        snapshots.append([asdict(r) for r in regions])
    assert snapshots[0] == snapshots[1] == snapshots[2]
    display = [p for r in regions for p in r.display_residues]
    assert len(display) == len(set(display))
    lookup = {f.footprint_id: f for f in raw}
    report = dict(
        source=str(source.resolve()),
        region_count=len(regions),
        retained_unique_footprints=len(unique),
        unique_display_residues=len(display),
        deterministic_three_runs=True,
        parameters=SURFACE_V2_PARAMETERS["regions"],
        regions=[],
    )
    for i, r in enumerate(regions, 1):
        f = lookup[r.representative_footprint_id]
        report["regions"].append(
            dict(
                region_id=r.region_id,
                primary_count=len(r.primary_footprint_ids),
                display_residue_count=len(r.display_residues),
                color=_region_color(i),
                center_residue=f.center_residue_id,
                mutable_count=len(f.mutable_residues),
                center_geodesic_extent=r.geodesic_extent,
                domains=r.domain_ids,
            )
        )
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "region_replay_report.json", report)
    write_json(output / "surface_regions.json", snapshots[0])
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(replay(args.source, args.output), indent=2))
