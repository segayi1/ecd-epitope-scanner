from __future__ import annotations

import math
from dataclasses import asdict

import numpy as np
import pytest
from ecd_epitope_scanner.models import ResidueRecord
from ecd_epitope_scanner.surface_patch.center_sampling import sample_surface_centers
from ecd_epitope_scanner.surface_patch.footprint_dedup import deduplicate_footprints
from ecd_epitope_scanner.surface_patch.footprint_generator import generate_footprints
from ecd_epitope_scanner.surface_patch.footprint_residue_mapping import map_footprints_to_residues
from ecd_epitope_scanner.surface_patch.models import SurfaceGeometry
from ecd_epitope_scanner.surface_patch.region_builder import build_surface_regions
from ecd_epitope_scanner.surface_patch.surface_graph import (
    build_surface_graph,
    connected_components,
)
from scipy.sparse.csgraph import dijkstra


def test_region_communities_keep_opposite_components_and_display_ownership():
    vertices = np.array([(float(x), float(y), z) for z in (0., 4.) for x in range(15) for y in range(10)])
    edges = []
    for offset in (0, 150):
        for x in range(15):
            for y in range(10):
                i = offset + x * 10 + y
                if x < 14:
                    edges.append((i, i + 10))
                if y < 9:
                    edges.append((i, i + 1))
    surface = _geometry(vertices, edges, np.tile((0., 0., 1.), (300, 1)))
    graph = build_surface_graph(surface)
    raw = generate_footprints(surface, graph, [0, 79, 149, 150, 229, 299])
    records = [ResidueRecord(i, 'A', 'A', i, folding_screen_pass=True) for i in range(1,301)]
    map_footprints_to_residues(raw, surface, records)
    before = [(f.footprint_id, list(f.vertex_ids), list(f.residues)) for f in raw]
    snapshots = []
    for _ in range(3):
        regions = build_surface_regions(raw, surface=surface, graph=graph, primary_candidates=raw)
        snapshots.append([asdict(r) for r in regions])
        display = [p for r in regions for p in r.display_residues]
        assert len(display) == len(set(display)) == 300
        assert sum(len(r.primary_footprint_ids) for r in regions) == 6
        by_id = {f.footprint_id:f for f in raw}
        for r in regions:
            assert len({surface.component_id[by_id[f].center_vertex_id] for f in r.footprint_ids}) == 1
            assert r.geodesic_extent <= 60
    assert snapshots[0] == snapshots[1] == snapshots[2]
    assert before == [(f.footprint_id, f.vertex_ids, f.residues) for f in raw]


def _geometry(vertices: np.ndarray, edges: list[tuple[int, int]], normals: np.ndarray) -> SurfaceGeometry:
    edge_array = np.asarray(edges, dtype=np.int64).reshape((-1, 2))
    lengths = np.asarray(
        [np.linalg.norm(vertices[left] - vertices[right]) for left, right in edge_array],
        dtype=float,
    )
    components = connected_components(len(vertices), edge_array)
    ids = np.arange(1, len(vertices) + 1, dtype=np.int64)
    return SurfaceGeometry(
        vertices=vertices,
        normals=normals,
        vertex_area=np.ones(len(vertices), dtype=float),
        edges=edge_array,
        edge_length=lengths,
        owner_residue=ids,
        owner_atom=ids,
        chain_id=np.asarray(["A"] * len(vertices), dtype=object),
        canonical_residue_id=ids,
        is_ecd=np.ones(len(vertices), dtype=bool),
        confidence=np.full(len(vertices), 90.0),
        component_id=components,
        backend="synthetic_mesh",
    )


def test_parallel_surfaces_have_no_cross_face_shortcut() -> None:
    top = np.asarray([(float(x), 0.0, 0.0) for x in range(6)])
    bottom = np.asarray([(float(x), 0.0, 4.0) for x in range(6)])
    vertices = np.vstack((top, bottom))
    edges = [(index, index + 1) for index in range(5)] + [(index + 6, index + 7) for index in range(5)]
    normals = np.vstack((np.tile((0.0, 0.0, 1.0), (6, 1)), np.tile((0.0, 0.0, -1.0), (6, 1))))
    surface = _geometry(vertices, edges, normals)
    graph = build_surface_graph(surface)
    assert surface.component_id[0] != surface.component_id[6]
    assert math.isinf(float(dijkstra(graph, indices=0)[6]))


def test_spherical_surface_geodesic_exceeds_chord_and_matches_arc() -> None:
    radius = 10.0
    angles = np.linspace(0.0, math.pi / 2.0, 33)
    vertices = np.column_stack((radius * np.cos(angles), radius * np.sin(angles), np.zeros(len(angles))))
    edges = [(index, index + 1) for index in range(len(vertices) - 1)]
    normals = vertices / radius
    surface = _geometry(vertices, edges, normals)
    graph = build_surface_graph(surface)
    distance = float(dijkstra(graph, indices=0)[-1])
    chord = float(np.linalg.norm(vertices[0] - vertices[-1]))
    assert distance > chord
    assert distance == pytest.approx(radius * math.pi / 2.0, rel=0.01)


def test_center_and_footprint_generation_is_deterministic() -> None:
    vertices = np.asarray([(float(x), 0.0, 0.0) for x in range(25)])
    edges = [(index, index + 1) for index in range(24)]
    surface = _geometry(vertices, edges, np.tile((0.0, 1.0, 0.0), (25, 1)))
    graph = build_surface_graph(surface)
    eligible = np.arange(25, dtype=np.int64)
    first = sample_surface_centers(surface, graph, eligible, 7.0, 7.0)
    second = sample_surface_centers(surface, graph, eligible, 7.0, 7.0)
    assert first == second
    left = generate_footprints(surface, graph, first)
    right = generate_footprints(surface, graph, second)
    assert [(item.footprint_id, item.vertex_ids) for item in left] == [
        (item.footprint_id, item.vertex_ids) for item in right
    ]


def test_concave_surface_follows_path_instead_of_crossing_cavity() -> None:
    vertices = np.asarray([
        (0.0, 3.0, 0.0), (0.0, 2.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
        (3.0, 0.0, 0.0), (3.0, 1.0, 0.0), (3.0, 2.0, 0.0), (3.0, 3.0, 0.0),
    ])
    edges = [(index, index + 1) for index in range(len(vertices) - 1)]
    surface = _geometry(vertices, edges, np.tile((0.0, 0.0, 1.0), (len(vertices), 1)))
    distance = float(dijkstra(build_surface_graph(surface), indices=0)[-1])
    assert distance == pytest.approx(9.0)
    assert distance > float(np.linalg.norm(vertices[0] - vertices[-1]))


def test_protected_residue_stays_in_footprint_but_not_mutable() -> None:
    vertices = np.asarray([(float(x), 0.0, 0.0) for x in range(4)])
    surface = _geometry(
        vertices,
        [(0, 1), (1, 2), (2, 3)],
        np.tile((0.0, 1.0, 0.0), (4, 1)),
    )
    graph = build_surface_graph(surface)
    footprints = generate_footprints(surface, graph, [0], radii=(9.0,))
    records = [
        ResidueRecord(
            uniprot_pos=index,
            aa="A",
            chain_id="A",
            auth_seq_id=index,
            ca=(float(index - 1), 0.0, 0.0),
            surface_accessible=True,
            folding_screen_pass=index != 2,
        )
        for index in range(1, 5)
    ]
    map_footprints_to_residues(
        footprints, surface, records, min_surface_area=0.0, min_surface_fraction=0.0
    )
    assert 2 in footprints[0].residues
    assert 2 in footprints[0].protected_residues
    assert 2 not in footprints[0].mutable_residues


def test_weighted_jaccard_dedup_is_deterministic() -> None:
    vertices = np.asarray([(float(x), 0.0, 0.0) for x in range(8)])
    surface = _geometry(
        vertices,
        [(index, index + 1) for index in range(7)],
        np.tile((0.0, 1.0, 0.0), (8, 1)),
    )
    graph = build_surface_graph(surface)
    footprints = generate_footprints(surface, graph, [0], radii=(6.0, 7.0))
    records = [
        ResidueRecord(index, "A", "A", index, folding_screen_pass=True)
        for index in range(1, 9)
    ]
    map_footprints_to_residues(
        footprints, surface, records, min_surface_area=0.0, min_surface_fraction=0.0
    )
    first_kept, first_removed = deduplicate_footprints(footprints, threshold=0.75)
    for item in footprints:
        item.duplicate_of = None
    second_kept, second_removed = deduplicate_footprints(footprints, threshold=0.75)
    assert [item.footprint_id for item in first_kept] == [item.footprint_id for item in second_kept]
    assert [item.footprint_id for item in first_removed] == [item.footprint_id for item in second_removed]
