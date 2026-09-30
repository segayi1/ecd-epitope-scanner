import json

import pytest
from ecd_epitope_scanner.joint_design_v4.models import V4Sector
from ecd_epitope_scanner.joint_design_v4.surface_snapshot import write_surface_graph_snapshot
from ecd_epitope_scanner.util import sha256_file


def _sector(sector_id, positions):
    return V4Sector(sector_id, positions, 1, 20.0, 20.0, 0.0)


def test_surface_snapshot_freezes_weighted_edges_and_inputs(tmp_path):
    mapping = tmp_path / "residue_mapping.csv"
    structure = tmp_path / "structure.pdb"
    mapping.write_text("position,aa\n10,A\n11,D\n", encoding="utf-8")
    structure.write_text("MODEL        1\nENDMDL\n", encoding="utf-8")
    path = tmp_path / "surface_graph_v4.json"
    digest = write_surface_graph_snapshot(
        path,
        {11: {10: 3.8}, 10: {11: 3.8}},
        {10: 12.0, 11: 8.0},
        [_sector("S1", [10, 11])],
        mapping_path=mapping,
        structure_path=structure,
        provenance={"basis": "explicit_solvent_surface_point_graph"},
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert digest == sha256_file(path)
    assert payload["mapping_sha256"] == sha256_file(mapping)
    assert payload["structure_sha256"] == sha256_file(structure)
    assert payload["nodes"] == [
        {"uniprot_pos": 10, "surface_area_a2": 12.0, "sector_id": "S1"},
        {"uniprot_pos": 11, "surface_area_a2": 8.0, "sector_id": "S1"},
    ]
    assert payload["edges"] == [{"source": 10, "target": 11, "distance_a": 3.8}]


def test_surface_snapshot_rejects_uncovered_nodes(tmp_path):
    mapping = tmp_path / "mapping.csv"
    mapping.write_text("mapping", encoding="utf-8")
    with pytest.raises(ValueError, match="cover"):
        write_surface_graph_snapshot(
            tmp_path / "graph.json",
            {10: {}, 11: {}},
            {10: 1.0, 11: 1.0},
            [_sector("S1", [10])],
            mapping_path=mapping,
            structure_path=None,
            provenance={},
        )
