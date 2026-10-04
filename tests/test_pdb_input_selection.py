from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from Bio.PDB import Atom, Chain, Model, Residue
from ecd_epitope_scanner.models import FeatureRange, ProteinContext, StructureCandidate
from ecd_epitope_scanner.pipeline import _apply_admission_gate
from ecd_epitope_scanner.structure_analysis import (
    ChainMapping,
    ParsedStructure,
    _interface_residue_pairs,
    _parse_chain_descriptions,
    analyze_residues,
    best_chain_mapping,
    calculate_sasa,
    infer_antigen_chain,
)


def _add_protein_chain(model, chain_id: str, residue_name: str, length: int, x: float):
    chain = Chain.Chain(chain_id)
    model.add(chain)
    residues = []
    for index in range(1, length + 1):
        residue = Residue.Residue((" ", index, " "), residue_name, " ")
        residue.add(
            Atom.Atom(
                "CA", np.asarray([x, float(index), 0.0]), 90.0, 1.0, " ", " CA ",
                index, element="C",
            )
        )
        chain.add(residue)
        residues.append(residue)
    return chain, residues


def test_sasa_uses_selected_protein_chain_without_crystal_solvent_or_partner(tmp_path):
    model = Model.Model(0)
    target, residues = _add_protein_chain(model, "A", "ALA", 1, 0.0)
    _add_protein_chain(model, "H", "ALA", 1, 0.2)
    water = Residue.Residue(("W", 2, " "), "HOH", " ")
    water.add(Atom.Atom("O", np.asarray([0.1, 1.0, 0.0]), 20.0, 1.0, " ", " O  ", 3, element="O"))
    target.add(water)
    parsed = SimpleNamespace(model=model)

    sasa, method = calculate_sasa(parsed, tmp_path, "A")

    assert sasa[("A", 1, "")][0] > 120.0
    assert ("A", 2, "") not in sasa
    assert ("H", 1, "") not in sasa
    assert not hasattr(residues[0]["CA"], "sasa")
    assert "isolated selected protein chain" in method


def test_antigen_chain_is_selected_from_fab_interface_and_entity_annotations():
    model = Model.Model(0)
    _add_protein_chain(model, "H", "TYR", 14, 0.0)
    antigen, _ = _add_protein_chain(model, "A", "ALA", 10, 3.0)
    _add_protein_chain(model, "X", "GLY", 16, 50.0)
    parsed = SimpleNamespace(
        model=model,
        chain_descriptions={
            "H": "Fab heavy chain",
            "A": "CTLA4 antigen",
            "X": "unrelated enzyme",
        },
    )

    assert _interface_residue_pairs(antigen, model["H"]) > 0
    selection = infer_antigen_chain(parsed)

    assert selection["chain_id"] == "A"
    assert selection["antibody_chains"] == ["H"]
    assert selection["method"] == "entity annotation and antibody interface"
    assert selection["confidence"] == "high"
    assert not selection["warnings"]


def test_mmcif_entity_descriptions_are_mapped_to_author_chain_ids(tmp_path):
    path = tmp_path / "complex.cif"
    path.write_text(
        """data_test
loop_
_entity.id
_entity.pdbx_description
1 'Fab heavy chain'
2 'CTLA4 antigen'
loop_
_atom_site.auth_asym_id
_atom_site.label_entity_id
H 1
A 2
""",
        encoding="utf-8",
    )

    descriptions = _parse_chain_descriptions(path, "mmcif")

    assert descriptions == {"H": "Fab heavy chain", "A": "CTLA4 antigen"}


def test_explicit_chain_overrides_antigen_chain_heuristic():
    model = Model.Model(0)
    _add_protein_chain(model, "H", "TYR", 14, 0.0)
    _add_protein_chain(model, "A", "ALA", 10, 3.0)
    parsed = SimpleNamespace(
        model=model,
        chain_descriptions={"H": "Fab heavy chain", "A": "antigen"},
    )

    selection = infer_antigen_chain(parsed, "H")

    assert selection["chain_id"] == "H"
    assert selection["method"] == "user specified"
    assert selection["confidence"] == "explicit"


def test_sequence_mapping_prefers_target_identity_over_a_longer_unrelated_chain():
    model = Model.Model(0)
    _add_protein_chain(model, "H", "TYR", 20, 0.0)
    _, target_residues = _add_protein_chain(model, "A", "ALA", 12, 30.0)
    parsed = SimpleNamespace(model=model)
    protein = ProteinContext(
        accession="TEST",
        protein_name="Test antigen",
        organism="test",
        sequence="A" * 12,
        ecd_start=1,
        ecd_end=12,
        membrane_anchor_position=12,
        topology="test",
    )

    mapping = best_chain_mapping(parsed, protein)

    assert mapping.chain_id == "A"
    assert mapping.identity == 1.0
    assert mapping.observed_ecd_coverage == 1.0
    assert len(target_residues) == 12


def test_failed_admission_thresholds_are_recorded_without_raising():
    protein = ProteinContext(
        accession="TEST",
        protein_name="Test antigen",
        organism="test",
        sequence="A" * 20,
        ecd_start=1,
        ecd_end=20,
        membrane_anchor_position=20,
        topology="test",
        features=[FeatureRange("Domain", 1, 20, "major domain")],
    )
    candidate = StructureCandidate(
        identifier="test",
        category="USER-PROVIDED",
        provider="user upload",
        model_url="local",
    )
    mapping = ChainMapping(
        chain_id="A",
        residues=[],
        query_to_uniprot={},
        uniprot_to_residue={index: object() for index in range(1, 5)},
        identity=0.5,
        observed_ecd_coverage=0.2,
        max_internal_gap=10,
        engineered_mismatches=0,
    )

    _apply_admission_gate(candidate, protein, mapping)

    assert not candidate.eligible
    assert candidate.rejection_reasons == [
        "observed ECD coverage < 95%",
        "internal coordinate gap > 8 aa",
        "aligned sequence identity < 95%",
        "one or more major ECD domains have < 90% observed coverage",
    ]


def test_engineered_structure_positions_are_omitted_from_residue_scoring(tmp_path):
    model = Model.Model(0)
    chain = Chain.Chain("A")
    model.add(chain)
    residues = []
    for index in range(1, 13):
        residue = Residue.Residue(
            (" ", index, " "), "VAL" if index == 5 else "ALA", " "
        )
        x = float(index * 3.8)
        atoms = (
            ("N", [x - 1.2, 0.0, 0.0], "N"),
            ("CA", [x, 0.0, 0.0], "C"),
            ("C", [x + 1.2, 0.0, 0.0], "C"),
            ("O", [x + 1.8, 0.8, 0.0], "O"),
            ("CB", [x, 1.5, 0.0], "C"),
        )
        for atom_index, (name, xyz, element) in enumerate(atoms):
            residue.add(
                Atom.Atom(
                    name, np.asarray(xyz), 90.0, 1.0, " ", f"{name:>4}",
                    (index - 1) * len(atoms) + atom_index, element=element,
                )
            )
        chain.add(residue)
        residues.append(residue)
    parsed = ParsedStructure(
        path=tmp_path / "source.pdb",
        format="pdb",
        structure=None,
        model=model,
        raw_text="",
    )
    protein = ProteinContext(
        accession="TEST",
        protein_name="Test antigen",
        organism="test",
        sequence="A" * 12,
        ecd_start=1,
        ecd_end=12,
        membrane_anchor_position=12,
        topology="test",
    )
    mapping = ChainMapping(
        chain_id="A",
        residues=residues,
        query_to_uniprot={index: index + 1 for index in range(12)},
        uniprot_to_residue={index + 1: residue for index, residue in enumerate(residues)},
        identity=11 / 12,
        observed_ecd_coverage=1.0,
        max_internal_gap=0,
        engineered_mismatches=1,
        mismatch_positions=[5],
    )

    records, methods = analyze_residues(
        parsed, mapping, protein, tmp_path, ignored_positions={5}
    )

    assert len(records) == 11
    assert 5 not in {item.uniprot_pos for item in records}
    assert methods["mapped_ecd_residues"] == 11
