from pathlib import Path

import numpy as np
from Bio.PDB import Atom, Chain, Model, Residue, Structure

from ecd_epitope_scanner.models import ProteinContext
from ecd_epitope_scanner.structure_analysis import (
    ChainMapping,
    ParsedStructure,
    analyze_residues,
    calculate_sidechain_exposure,
)


def _atom(name, coord, serial):
    return Atom.Atom(
        name,
        np.asarray(coord, dtype=float),
        90.0,
        1.0,
        " ",
        f"{name:>4}",
        serial,
        element=name[0],
    )


def _parsed(residues):
    structure = Structure.Structure("test")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    for residue in residues:
        chain.add(residue)
    model.add(chain)
    structure.add(model)
    return ParsedStructure(Path("test.pdb"), "pdb", structure, model, "")


def _ala(seq, cb):
    residue = Residue.Residue((" ", seq, " "), "ALA", " ")
    residue.add(_atom("CB", cb, seq))
    return residue


def test_single_exposed_sidechain_has_full_relative_exposure():
    parsed = _parsed([_ala(1, (0, 0, 0))])

    result = calculate_sidechain_exposure(parsed, "A")[("A", 1, "")]

    assert result["status"] == "KNOWN"
    assert result["value"] == 1.0
    assert result["numerator_a2"] == result["denominator_a2"] > 0


def test_neighbor_buries_sidechain_and_does_not_change_coordinates():
    residues = [_ala(1, (0, 0, 0)), _ala(2, (3, 0, 0))]
    parsed = _parsed(residues)
    before = [[atom.coord.copy() for atom in residue] for residue in residues]

    result = calculate_sidechain_exposure(parsed, "A")[("A", 1, "")]

    assert result["status"] == "KNOWN"
    assert 0 < result["value"] < 1
    assert all(
        np.array_equal(atom.coord, coord)
        for residue, coords in zip(residues, before, strict=True)
        for atom, coord in zip(residue, coords, strict=True)
    )


def test_incomplete_sidechain_is_unknown_with_no_imputed_value():
    residue = Residue.Residue((" ", 1, " "), "ASP", " ")
    residue.add(_atom("CB", (0, 0, 0), 1))
    result = calculate_sidechain_exposure(_parsed([residue]), "A")[("A", 1, "")]

    assert result["status"] == "UNKNOWN_INCOMPLETE_SIDECHAIN"
    assert result["value"] is None
    assert result["numerator_a2"] is None
    assert result["denominator_a2"] is None


def test_invalid_reused_native_sasa_remains_unknown():
    parsed = _parsed([_ala(1, (0, 0, 0))])
    for invalid in (float("nan"), float("inf"), -1.0):
        result = calculate_sidechain_exposure(parsed, "A", {("A", 1, ""): (1.0, invalid)})[("A", 1, "")]
        assert result["value"] is None and result["status"].startswith("UNKNOWN")


def test_isolated_residue_copy_has_no_parent_structure(monkeypatch):
    from Bio.PDB.SASA import ShrakeRupley
    original = ShrakeRupley.compute
    observed = []

    def compute(calculator, entity, *args, **kwargs):
        if entity.level == "R":
            observed.append(entity.get_parent())
        return original(calculator, entity, *args, **kwargs)

    monkeypatch.setattr(ShrakeRupley, "compute", compute)
    calculate_sidechain_exposure(_parsed([_ala(1, (0, 0, 0)), _ala(2, (3, 0, 0))]), "A")
    assert observed == [None, None]


def test_analyze_residues_attaches_sidechain_exposure(tmp_path):
    residue = _ala(1, (0, 1.5, 0))
    for serial, (name, coord) in enumerate(
        (("N", (-1, 0, 0)), ("CA", (0, 0, 0)), ("C", (1, 0, 0)), ("O", (1, -1, 0))),
        start=2,
    ):
        residue.add(_atom(name, coord, serial))
    parsed = _parsed([residue])
    protein = ProteinContext("TEST", "test", "test", "A", 1, 1, 2, "single-pass")
    mapping = ChainMapping("A", [residue], {0: 1}, {1: residue}, 1.0, 1.0, 0, 0)

    records, _ = analyze_residues(parsed, mapping, protein, tmp_path)

    assert records[0].sidechain_exposure["status"] == "KNOWN"
    assert records[0].sidechain_exposure["value"] is not None
