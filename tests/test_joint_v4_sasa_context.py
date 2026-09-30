from __future__ import annotations

import hashlib

import numpy as np
from Bio import __version__ as biopython_version
from Bio.PDB import PDBIO, Atom, Chain, Model, Residue, Structure
from ecd_epitope_scanner.joint_design_v4.contexts import (
    normalized_prediction_surface_rows,
    residue_keys,
    write_normalized_prediction_chain,
)
from ecd_epitope_scanner.joint_design_v4.engine import _validate_procaliper_input_scope
from ecd_epitope_scanner.joint_design_v4.evidence import surface_evidence
from ecd_epitope_scanner.joint_design_v4.models import PredictionResult
from ecd_epitope_scanner.models import ResidueRecord
from ecd_epitope_scanner.structure_analysis import TIEN_DSSP_RADII
from ecd_epitope_scanner.util import sha256_file


def _records() -> list[ResidueRecord]:
    return [
        ResidueRecord(
            uniprot_pos=position,
            aa="A",
            chain_id="A",
            auth_seq_id=auth_position,
            ca=(float(auth_position * 3.2), 0.0, 0.0),
            sasa=100.0,
            sidechain_sasa=35.0,
            rsasa=0.5,
            surface_accessible=True,
        )
        for position, auth_position in zip((10, 11, 12), (2, 3, 4), strict=True)
    ]


def _write_source_pdb(path) -> None:
    structure = Structure.Structure("source")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    model.add(chain)
    structure.add(model)
    serial = 1
    for position in range(1, 5):
        x = float(position * 3.2)
        residue = Residue.Residue((" ", position, " "), "ALA", " ")
        for name, xyz, element in (
            ("CA", [x, 0.0, 0.0], "C"),
            ("CB", [x, 1.5, 0.0], "C"),
        ):
            residue.add(
                Atom.Atom(
                    name,
                    np.asarray(xyz, dtype=float),
                    90.0,
                    1.0,
                    " ",
                    f"{name:>4}",
                    serial,
                    element=element,
                )
            )
            serial += 1
        chain.add(residue)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path))


def test_normalized_prediction_sasa_is_measured_on_the_worker_input(tmp_path):
    source = tmp_path / "source.pdb"
    normalized = tmp_path / "prediction_context.pdb"
    records = _records()
    _write_source_pdb(source)

    assert write_normalized_prediction_chain(source, normalized, records, "A") == normalized
    rows, method = normalized_prediction_surface_rows(normalized, "A", records)

    assert set(rows) == {10, 11, 12}
    assert all(value["sasa"] > 0 and value["rsasa"] > 0 for value in rows.values())
    assert method["scope"].startswith("normalized mapped canonical ECD")
    assert method["biopython_version"] == biopython_version


def test_sasa_method_comparison_separates_normalized_scope_from_source_context():
    records = _records()[:1]
    records[0].auth_seq_id = 2
    normalized_rows = {10: {"sasa": 20.0, "rsasa": 0.20}}
    native_method = {
        "algorithm": "Bio.PDB Shrake-Rupley",
        "biopython_version": biopython_version,
        "probe_radius_a": 1.4,
        "sphere_points": 960,
        "radii_reference": "Tien et al. 2013/DSSP",
        "atom_radii_a": dict(TIEN_DSSP_RADII),
        "scope": "normalized mapped canonical ECD; no partners or engineered mismatches",
    }
    procaliper_provenance = {
        "algorithm": "Bio.PDB ShrakeRupley via ProCaliper",
        "biopython_version": biopython_version,
        "probe_radius_a": 1.4,
        "sphere_points": 960,
        "rsasa_reference": "Tien et al. 2013 maximum ASA",
        "atom_radii_a": dict(TIEN_DSSP_RADII),
    }

    evidence = surface_evidence(
        records,
        residue_keys(records),
        [{"position": 10, "sasa": 20.0, "rsasa": 0.20}],
        normalized_context_rows=normalized_rows,
        normalized_context_method=native_method,
        procaliper_provenance=procaliper_provenance,
    )[0]

    assert evidence.sasa_agreement_status == "AGREE"
    assert evidence.comparison_method_status == "MATCHED"
    assert evidence.normalized_context_rsasa == 0.20
    assert evidence.native_context_status == "CONTEXT_SENSITIVE"
    assert evidence.native_context_classification_changed is True
    assert evidence.procaliper_classification_changed is False


def test_sasa_comparison_does_not_call_different_sampling_settings_agreement():
    record = _records()[0]
    method = {
        "algorithm": "Bio.PDB Shrake-Rupley",
        "biopython_version": biopython_version,
        "probe_radius_a": 1.4,
        "sphere_points": 960,
        "radii_reference": "Tien et al. 2013/DSSP",
        "atom_radii_a": dict(TIEN_DSSP_RADII),
        "scope": "normalized mapped canonical ECD",
    }
    provenance = {
        "algorithm": "Bio.PDB ShrakeRupley via ProCaliper",
        "biopython_version": biopython_version,
        "probe_radius_a": 1.4,
        "sphere_points": 480,
        "rsasa_reference": "Tien et al. 2013 maximum ASA",
        "atom_radii_a": dict(TIEN_DSSP_RADII),
    }

    evidence = surface_evidence(
        [record],
        residue_keys([record]),
        [{"position": 10, "sasa": 20.0, "rsasa": 0.20}],
        normalized_context_rows={10: {"sasa": 20.0, "rsasa": 0.20}},
        normalized_context_method=method,
        procaliper_provenance=provenance,
    )[0]

    assert evidence.comparison_method_status == "SETTINGS_DIFFER"
    assert evidence.sasa_agreement_status == "METHOD_SETTINGS_UNVERIFIED"


def test_missing_procaliper_radii_remain_diagnostic():
    record = _records()[0]
    method = {
        "algorithm": "Bio.PDB Shrake-Rupley", "biopython_version": biopython_version,
        "probe_radius_a": 1.4, "sphere_points": 960,
        "radii_reference": "Tien et al. 2013/DSSP",
        "atom_radii_a": dict(TIEN_DSSP_RADII),
    }
    provenance = {
        "algorithm": "Bio.PDB ShrakeRupley via ProCaliper",
        "biopython_version": biopython_version,
        "probe_radius_a": 1.4, "sphere_points": 960,
        "rsasa_reference": "Tien et al. 2013 maximum ASA",
    }
    evidence = surface_evidence(
        [record], residue_keys([record]),
        [{"position": 10, "sasa": 20.0, "rsasa": 0.20}],
        normalized_context_rows={10: {"sasa": 20.0, "rsasa": 0.20}},
        normalized_context_method=method,
        procaliper_provenance=provenance,
    )[0]
    assert evidence.comparison_method_status == "PROVENANCE_INCOMPLETE"
    assert evidence.sasa_agreement_status == "METHOD_SETTINGS_UNVERIFIED"


def test_biopython_version_difference_is_visible_but_values_remain_comparable():
    record = _records()[0]
    method = {
        "algorithm": "Bio.PDB Shrake-Rupley",
        "biopython_version": biopython_version,
        "probe_radius_a": 1.4,
        "sphere_points": 960,
        "radii_reference": "Tien et al. 2013/DSSP",
        "atom_radii_a": dict(TIEN_DSSP_RADII),
        "scope": "normalized mapped canonical ECD",
    }
    provenance = {
        "algorithm": "Bio.PDB ShrakeRupley via ProCaliper",
        "biopython_version": "different-installed-version",
        "probe_radius_a": 1.4,
        "sphere_points": 960,
        "rsasa_reference": "Tien et al. 2013 maximum ASA",
        "atom_radii_a": dict(TIEN_DSSP_RADII),
    }
    evidence = surface_evidence(
        [record],
        residue_keys([record]),
        [{"position": 10, "sasa": 20.0, "rsasa": 0.20}],
        normalized_context_rows={10: {"sasa": 20.0, "rsasa": 0.20}},
        normalized_context_method=method,
        procaliper_provenance=provenance,
    )[0]

    assert evidence.comparison_method_status == "VERSION_DIFFERS"
    assert evidence.sasa_agreement_status == "AGREE"


def test_procaliper_scope_validation_checks_exact_structure_chain_and_sequence(tmp_path):
    normalized = tmp_path / "prediction_context.pdb"
    normalized.write_text("normalized input\n", encoding="ascii")
    sequence = "ACD"
    provenance = {
        "input_structure_sha256": sha256_file(normalized),
        "chain_id": "A",
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "residue_count": len(sequence),
    }
    result = PredictionResult(status="success", records=[{}], provenance=provenance)

    assert _validate_procaliper_input_scope(result, normalized, "A", sequence) is result

    mismatched = PredictionResult(
        status="success",
        records=[{}],
        provenance={**provenance, "input_structure_sha256": "0" * 64},
    )
    failed = _validate_procaliper_input_scope(mismatched, normalized, "A", sequence)
    assert failed.status == "failed"
    assert failed.errors[0]["code"] == "PROCALIPER_INPUT_SCOPE_UNVERIFIED"
    assert failed.provenance["scope_validation_fields"] == ["input_structure_sha256"]
