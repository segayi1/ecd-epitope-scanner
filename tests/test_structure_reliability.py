import json
from dataclasses import asdict

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, Residue
from ecd_epitope_scanner.models import ProteinContext, ResidueRecord
from ecd_epitope_scanner.structure_reliability import (
    classify_reliability,
    export_reliability,
    import_iupred,
)
from ecd_epitope_scanner.surface_patch import run_surface_patch_v2


def fixture(n=40):
    protein = ProteinContext(
        "TEST", "test", "test", "A" * n, 1, n, n + 1, "single-pass"
    )
    records = [
        ResidueRecord(
            i,
            "A",
            "A",
            i,
            ca=(float(i * 10), 0.0, 0.0),
            plddt=90.0 if i <= 10 else 30.0,
            token_index=i - 1,
            surface_accessible=True,
            folding_screen_pass=True,
        )
        for i in range(1, n + 1)
    ]
    pae = np.full((n, n), 25.0)
    np.fill_diagonal(pae, 0.0)
    return protein, records, {"pae": pae.tolist()}


def classify(protein, records, confidence=None, **kwargs):
    return classify_reliability(
        records,
        protein,
        predicted=True,
        confidence=confidence,
        chain_length=len(records),
        single_chain=True,
        **kwargs,
    )


def test_long_tail_multievidence_and_records_unchanged(tmp_path):
    p, records, confidence = fixture()
    before = [asdict(r) for r in records]
    for r in records:
        r.ca = np.asarray(r.ca, dtype=np.float32)
    result = classify(p, records, confidence)
    json.dumps(result)
    assert result["qc"]["excluded_positions"] == list(range(14, 41))
    assert result["qc"]["boundary_positions"] == [12, 13]
    covered = {
        i for c in result["idr_candidates"] for i in range(c["start"], c["end"] + 1)
    }
    assert covered == set(range(14, 41))
    assert all(
        not c["run_foldx"] and c["sasa"] is None for c in result["idr_candidates"]
    )
    assert all(
        tuple(r.ca) == old["ca"] and r.sasa == old["sasa"]
        for r, old in zip(records, before)
    )
    export_reliability(result, tmp_path)
    assert len(list(tmp_path.iterdir())) == 5


def test_missing_evidence_and_experimental_do_not_delete():
    p, records, confidence = fixture()
    assert not classify(p, records)["qc"]["excluded_positions"]
    result = classify_reliability(records, p, predicted=False, confidence=confidence)
    assert result["qc"]["state_counts"] == {"STRUCTURED": 40}


def test_embedded_loop_and_mobile_domain():
    p, records, confidence = fixture()
    for r in records:
        r.plddt = 30.0 if 15 <= r.uniprot_pos <= 20 else 90.0
    result = classify(p, records, confidence)
    assert result["qc"]["state_counts"] == {"STRUCTURED": 34, "FLEXIBLE_STRUCTURED": 6}
    assert not result["qc"]["excluded_positions"]


def test_unmapped_pae_is_ignored_and_missing_coordinates_retained_as_unresolved():
    p, records, confidence = fixture()
    result = classify_reliability(
        records, p, predicted=True, confidence=confidence, single_chain=False
    )
    assert result["qc"]["pae_status"] == "unmapped_pae_ignored"
    assert not result["qc"]["excluded_positions"]
    result = classify(p, records[:-1])
    assert result["residues"][-1]["state"] == "UNRESOLVED"


def test_iupred_validation(tmp_path):
    path = tmp_path / "iupred.txt"
    path.write_text("# IUPred2A\n1 A 0.8\n2 C 0.7\n")
    scores, provenance = import_iupred(path, "AC")
    assert scores == {1: 0.8, 2: 0.7} and provenance["sha256"]
    with pytest.raises(ValueError):
        import_iupred(path, "AA")
    with pytest.raises(ValueError):
        import_iupred(path, "ACA")


def model_for(records, shift=0.0):
    model = Model.Model(0)
    chain = Chain.Chain("A")
    model.add(chain)
    mapping = {}
    for r in records:
        residue = Residue.Residue((" ", r.uniprot_pos, " "), "ALA", " ")
        xyz = np.array(r.ca) + (
            np.array([0.0, shift, 0.0]) if r.uniprot_pos >= 14 else 0.0
        )
        residue.add(
            Atom.Atom("CA", xyz, 90.0, 1.0, " ", " CA ", r.uniprot_pos, element="C")
        )
        chain.add(residue)
        mapping[r.uniprot_pos] = residue
    return model, mapping


def test_strict_tail_invariance_boundary_guard_and_disabled_compatibility():
    p, records, confidence = fixture()
    reliability = classify(p, records, confidence)
    snapshots = []
    for shift in (0.0, 100.0):
        model, mapping = model_for(records, shift)
        result = run_surface_patch_v2(
            model, "A", mapping, records, p, reliability=reliability
        )
        snapshots.append(result.surface.vertices)
        assert not set(result.surface.owner_residue) & set(range(12, 41))
    np.testing.assert_array_equal(*snapshots)
    model, mapping = model_for(records)
    off = classify(p, records, confidence, parameters={"enabled": False})
    a = run_surface_patch_v2(model, "A", mapping, records, p, reliability=off)
    b = run_surface_patch_v2(model, "A", mapping, records, p)
    np.testing.assert_array_equal(a.surface.vertices, b.surface.vertices)


def test_all_disordered_empty_surface():
    p, records, _ = fixture(20)
    for r in records:
        r.plddt = 30.0
    reliability = classify(p, records, sequence_scores={i: 0.9 for i in range(1, 21)})
    model, mapping = model_for(records)
    result = run_surface_patch_v2(
        model, "A", mapping, records, p, reliability=reliability
    )
    assert len(result.surface.vertices) == 0 and not result.legacy_regions


def test_conditional_requires_sequence_support_and_preserves_protected_windows():
    p, records, confidence = fixture()
    p.sequence = "ACGPNAST" + "A" * 32
    for r in records:
        r.plddt = 90.0
        r.secondary_structure = "H"
    assert not classify(p, records, confidence)["idr_candidates"]
    result = classify(
        p, records, confidence, sequence_scores={i: 0.9 for i in range(1, 41)}
    )
    assert result["qc"]["state_counts"] == {"CONDITIONALLY_FOLDED": 40}
    first = result["idr_candidates"][0]
    assert {2, 3, 4, 5, 6, 7} <= set(first["protected_positions"])
    assert not set(first["review_positions"]) & set(first["protected_positions"])
    assert not first["run_foldx"]


def test_occlusion_sensitivity_retains_tail_as_blocker():
    p, records, confidence = fixture()
    strict = classify(p, records, confidence)
    sensitivity = classify(
        p, records, confidence, parameters={"surface_mode": "occlusion_sensitivity"}
    )
    model, mapping = model_for(records)
    # Move one excluded tail atom onto a retained surface shell.
    mapping[14]["CA"].coord = mapping[1]["CA"].coord + np.array([3.0, 0.0, 0.0])
    a = run_surface_patch_v2(model, "A", mapping, records, p, reliability=strict)
    b = run_surface_patch_v2(model, "A", mapping, records, p, reliability=sensitivity)
    assert len(a.surface.vertices) > len(b.surface.vertices)
    assert 14 not in b.surface.owner_residue
