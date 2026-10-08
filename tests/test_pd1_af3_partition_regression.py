from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest
import numpy as np
from Bio.PDB import Atom, Chain, Model, Residue, MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from ecd_epitope_scanner.joint_design.surface_units import build_surface_units
from ecd_epitope_scanner.joint_design_v4.engine import _legacy_regions
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, SurfaceEvidence
from ecd_epitope_scanner.joint_design_v4.sectors import _distance, build_epitope_sectors
from ecd_epitope_scanner.models import ProteinContext, ResidueRecord
from ecd_epitope_scanner.scripts.replay_joint_v41_run import replay_request_data
from ecd_epitope_scanner.structure_reliability import classify_reliability
from ecd_epitope_scanner.surface_patch import run_surface_patch_v2


@pytest.fixture
def pd1():
    data = json.loads((Path(__file__).parent / "fixtures/pd1_af3_partition.json").read_text(encoding="utf-8"))
    records = [ResidueRecord(**r) for r in data["records"]]
    protein = ProteinContext("Q15116", "PD1", "PD1", data["sequence"], 25, 170, 171, "single-pass")
    return data, records, protein


def surface_rows(graph, areas):
    return [SurfaceEvidence(ResidueKey("0", "A", p, uniprot_pos=p), mapped_surface_area_a2=areas[p]) for p in graph]


def assert_partition(sectors, graph, geo=26.0, euclidean=18.0):
    memberships = [p for s in sectors for p in s.residue_positions]
    assert sorted(memberships) == sorted(graph)
    assert len(memberships) == len(set(memberships))
    for s in sectors:
        positions = set(s.residue_positions)
        assert len(positions) >= 4
        assert s.total_area_a2 >= 50.0
        assert set(_distance(graph, min(positions), positions)) == positions
        assert s.geodesic_diameter_a <= geo + 0.001
        assert s.euclidean_diameter_a <= euclidean + 0.001
        assert not {"DISCONNECTED_SECTOR", "GEODESIC_DIAMETER_EXCEEDED", "EUCLIDEAN_DIAMETER_EXCEEDED"} & set(s.quality_flags)


def test_original_pd1_graph_is_repaired_before_output(pd1):
    data, records, _ = pd1
    graph = {r["uniprot_pos"]: {} for r in data["graph_nodes"]}
    areas = {r["uniprot_pos"]: r["surface_area_a2"] for r in data["graph_nodes"]}
    for edge in data["graph_edges"]:
        graph[edge["source"]][edge["target"]] = graph[edge["target"]][edge["source"]] = edge["distance_a"]
    sectors = build_epitope_sectors(surface_rows(graph, areas), [], graph, {r.uniprot_pos: r for r in records})
    assert_partition(sectors, graph)
    assert len(sectors) > 7
    regions = _legacy_regions(sectors, records)
    assert len({r.color for r in regions}) == len(sectors)
    ownership = {p: r.id for r in regions for p in r.residues}
    assert all(r.region_id == ownership.get(r.uniprot_pos) for r in records)
    assert all(set(r.mutable_residues) <= set(r.residues) for r in regions)
    protected = {r.uniprot_pos for r in records if r.hard_exclusions}
    assert all(not set(r.mutable_residues) & protected for r in regions)


def test_af3_uncertain_tail_is_geometry_ineligible_without_claiming_disorder(pd1):
    data, records, protein = pd1
    result = classify_reliability(records, protein, predicted=True, confidence=data["confidence"])
    omitted = set(result["qc"]["uncertain_terminal_geometry_positions"])
    assert set(range(25, 32)) | set(range(150, 171)) <= omitted
    assert not result["qc"]["excluded_positions"]
    assert not result["idr_candidates"]
    rows = {r["residue_id"]: r for r in result["residues"]}
    assert all(rows[p]["state"] == "UNCERTAIN" and not rows[p]["use_for_3d_footprint"] for p in omitted)
    assert all(rows[p]["use_for_3d_footprint"] for p in [58, 59, 60, 61, 86, 87, 88, 89, 90, 91])
    # No deletion from absent PAE or reinterpretation of experimental B factors.
    no_pae = classify_reliability(records, protein, predicted=True)
    assert not no_pae["qc"]["uncertain_terminal_geometry_positions"]
    experimental = classify_reliability(records, protein, predicted=False, confidence=data["confidence"])
    assert all(r["use_for_3d_footprint"] for r in experimental["residues"])


def test_count_target_cannot_force_long_chain_into_one_region():
    records = {p: ResidueRecord(p, "A", "A", p, ca=(p * 3.8, 0, 0)) for p in range(1, 31)}
    graph = {p: {} for p in records}
    for p in range(1, 30):
        graph[p][p + 1] = graph[p + 1][p] = 3.8
    sectors = build_epitope_sectors(surface_rows(graph, {p: 80 for p in graph}), [], graph, records, max_sectors=1)
    assert len(sectors) > 1
    assert_partition(sectors, graph)


def test_geometry_mask_reaches_surface_graph_and_footprints(pd1):
    data, records, protein = pd1
    reliability = classify_reliability(records, protein, predicted=True, confidence=data["confidence"])
    model = Model.Model(0)
    chain = Chain.Chain("A")
    model.add(chain)
    mapping = {}
    for record in records:
        residue = Residue.Residue((" ", record.auth_seq_id, " "), "ALA", " ")
        residue.add(Atom.Atom("CA", np.asarray(record.ca), record.plddt, 1.0, " ", " CA ", record.auth_seq_id, element="C"))
        chain.add(residue)
        mapping[record.uniprot_pos] = residue
    result = run_surface_patch_v2(model, "A", mapping, records, protein, reliability=reliability)
    invalid = set(reliability["qc"]["footprint_excluded_positions"]) | set(reliability["qc"]["boundary_positions"])
    assert not set(result.surface.owner_residue) & invalid
    assert all(not set(f.residues) & invalid for f in result.raw_footprints)
    _, graph, areas = build_surface_units(records, surface_geometry=result.surface)
    assert not set(graph) & invalid
    sectors = build_epitope_sectors(surface_rows(graph, areas), [], graph, {r.uniprot_pos: r for r in records})
    assert_partition(sectors, graph)


def test_af3_replay_retains_model_type_and_confidence(tmp_path):
    source = tmp_path / "source"
    upload = source / "af3_upload"
    upload.mkdir(parents=True)
    (upload / "official_af3_model.cif").write_text("model")
    confidence = upload / "official_af3_confidences.json"
    confidence.write_text("{}")
    analysis = {"protein": {"accession": "Q15116"}, "source": {}, "structure_selection": {}}
    request = replay_request_data(source, tmp_path / "output", analysis, {"input_type": "af3_upload", "strict_af3": True})
    assert request["input_type"] == "af3_upload"
    assert request["strict_af3"]
    assert Path(request["confidence_path"]) == confidence
    confidence.unlink()
    with pytest.raises(ValueError, match="confidence"):
        replay_request_data(source, tmp_path / "output", analysis, {"input_type": "af3_upload"})


def test_frontend_queries_exported_author_ids_not_input_af3_labels(tmp_path):
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required to evaluate the real frontend selector")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    model.add(chain)
    for pos in (25, 26):
        residue = Residue.Residue((" ", pos, " "), "ALA", " ")
        residue.add(Atom.Atom("CA", np.array([float(pos), 0, 0]), 90, 1, " ", " CA ", pos, element="C"))
        chain.add(residue)
    writer = MMCIFIO()
    writer.set_structure(model)
    output = tmp_path / "ecd.cif"
    writer.save(str(output))
    exported = MMCIF2Dict(str(output))
    assert exported["_atom_site.auth_seq_id"] == ["25", "26"]
    assert exported["_atom_site.label_seq_id"] == ["1", "2"]
    app = Path(__file__).parents[1] / "web/app.js"
    script = """
const fs=require('node:fs'),vm=require('node:vm');
const source=fs.readFileSync(process.argv[1],'utf8');
const fn=source.slice(source.indexOf('function residueQuery('),source.indexOf('let molstarReady'));
const context={};vm.createContext(context);vm.runInContext(fn,context);
console.log(JSON.stringify([
context.residueQuery({chain_id:'A',auth_seq_id:25,label_chain_id:'A',label_seq_id:25}),
context.residueQuery({chain_id:'B',auth_seq_id:26,insertion_code:'A',label_chain_id:'X',label_seq_id:90})
]));
"""
    result = subprocess.run([node, "-e", script, str(app)], capture_output=True, text=True, check=True)
    queries = json.loads(result.stdout)
    assert queries == [
        {"auth_asym_id": "A", "auth_seq_id": 25},
        {"auth_asym_id": "B", "auth_seq_id": 26, "pdbx_PDB_ins_code": "A"},
    ]
