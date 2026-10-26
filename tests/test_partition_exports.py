"""Engine partition gates and exports, without external scientific calls."""

import json
from types import SimpleNamespace

import pytest

from ecd_epitope_scanner.joint_design_v4 import engine
from ecd_epitope_scanner.joint_design_v4.sectors import PartitionConstraintError
from ecd_epitope_scanner.joint_design_v4.surface_snapshot import write_surface_graph_snapshot
from ecd_epitope_scanner.models import ProteinContext, ResidueRecord


def _inputs(count):
    sequence = "DEKRFWYAILMS"[:count]
    records = [
        ResidueRecord(
            uniprot_pos=position, aa=aa, chain_id="A", auth_seq_id=position,
            ca=(position * 3.4, 0.0, 0.0), sasa=100.0,
            sidechain_sasa=40.0, rsasa=0.5, surface_accessible=True,
        )
        for position, aa in enumerate(sequence, 1)
    ]
    protein = ProteinContext(
        accession="TEST", protein_name="test", organism="test", sequence=sequence,
        ecd_start=1, ecd_end=count, membrane_anchor_position=count, topology="test",
    )
    return records, protein


def _read(output_dir, name):
    return json.loads((output_dir / name).read_text(encoding="utf-8"))


@pytest.mark.parametrize("surface_v2", [False, True])
def test_partition_success_exports_share_membership_and_diagnostics(tmp_path, surface_v2):
    records, protein = _inputs(12)
    result = engine.run_joint_v4_ml(
        records, protein, tmp_path, disable_external_engines=True,
        parameters={"mode": "strict", "use_hypothesis_hierarchy": False},
        surface_v2_result=(
            SimpleNamespace(surface_regions=[], footprints=[], raw_footprints=[])
            if surface_v2 else None
        ),
    )
    diagnostics = _read(tmp_path, "partition_diagnostics_v4.json")
    assert diagnostics["status"] == "COMPLETE"
    assert diagnostics["version"] == "balanced_connected_v4"
    assert diagnostics["minimum_sector_residues"] == 4
    assert diagnostics["minimum_sector_area_a2"] == 50.0
    assert diagnostics["count_is_soft_target"] is True
    assert diagnostics["target_sector_residues"] == 8
    assert diagnostics["preferred_minimum_sector_count"] == 4
    assert diagnostics["preferred_maximum_sector_count"] == 12
    sectors = _read(tmp_path, "probe_sectors_v4.json")
    expected = {
        sector.sector_id: sector.residue_positions for sector in result.sectors
    }
    assert {row["sector_id"]: row["residue_positions"] for row in sectors} == expected
    assert [region.residues for region in result.legacy_regions] == list(expected.values())
    owners = {position: key for key, positions in expected.items() for position in positions}
    assert set(owners) == set(range(1, 13))
    assert all(len(row["residue_positions"]) >= 4 and row["total_area_a2"] >= 50 for row in sectors)
    snapshot = _read(tmp_path, "surface_graph_v4.json")
    assert {node["uniprot_pos"]: node["sector_id"] for node in snapshot["nodes"]} == owners
    assert snapshot["graph_provenance"]["display_partition_policy"] == diagnostics
    assert result.provenance["surface_graph"]["display_partition_policy"] == diagnostics
    for name in ("v4_manifest.json", "v4_manifest_v41.json"):
        assert _read(tmp_path, name)["provenance"]["surface_graph"]["display_partition_policy"] == diagnostics
    assert "partition_diagnostics_v4.json" in _read(tmp_path, "v4_manifest_v41.json")["v41_artifacts"]
    if surface_v2:
        display = _read(tmp_path, "surface_display_regions_v41.json")
        assert display["regions"] == sectors
        assert display["partition_policy"] == diagnostics
        physical = _read(tmp_path, "surface_graph_v41.json")
        assert physical["graph_provenance"]["display_partition_policy"] == diagnostics
        assert {node["uniprot_pos"] for node in physical["nodes"]} == set(owners)
        assert all("sector_id" not in node for node in physical["nodes"])


def test_partition_failure_saves_unowned_graph_and_stops_design(tmp_path, monkeypatch):
    records, protein = _inputs(3)

    def must_not_generate(*args, **kwargs):
        pytest.fail("Candidate generation must not run after partition failure")

    monkeypatch.setattr(engine, "generate_candidates", must_not_generate)
    with pytest.raises(PartitionConstraintError) as caught:
        engine.run_joint_v4_ml(
            records, protein, tmp_path, disable_external_engines=True,
            parameters={"mode": "strict"},
        )
    diagnostics = _read(tmp_path, "partition_diagnostics_v4.json")
    assert diagnostics == caught.value.diagnostics
    assert diagnostics["status"] == "FAILED"
    assert diagnostics["draft_regions"]
    snapshot = _read(tmp_path, "surface_graph_v4.json")
    assert {node["uniprot_pos"] for node in snapshot["nodes"]} == {1, 2, 3}
    assert snapshot["edges"]
    assert all("sector_id" not in node for node in snapshot["nodes"])
    assert snapshot["graph_provenance"]["display_partition_policy"] == diagnostics
    assert not (tmp_path / "probe_sectors_v4.json").exists()
    assert not (tmp_path / "selected_panel_v4.json").exists()
    assert not (tmp_path / "v4_manifest.json").exists()


def test_snapshot_requires_explicit_unowned_mode_for_failed_partition(tmp_path):
    mapping = tmp_path / "mapping.csv"
    mapping.write_text("uniprot_pos\n1\n", encoding="utf-8")
    output = tmp_path / "graph.json"
    arguments = dict(mapping_path=mapping, structure_path=None, provenance={})
    with pytest.raises(ValueError, match="requires verified"):
        write_surface_graph_snapshot(output, {1: {}}, {1: 100.0}, None, **arguments)
    with pytest.raises(ValueError, match="cover the surface graph exactly"):
        write_surface_graph_snapshot(
            output, {1: {}}, {1: 100.0}, [], include_display_ownership=False, **arguments
        )
    assert not output.exists()
    write_surface_graph_snapshot(
        output, {1: {}}, {1: 100.0}, None, include_display_ownership=False, **arguments
    )
    assert _read(tmp_path, "graph.json")["nodes"] == [
        {"uniprot_pos": 1, "surface_area_a2": 100.0}
    ]
