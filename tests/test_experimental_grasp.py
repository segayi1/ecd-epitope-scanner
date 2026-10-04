from __future__ import annotations

import io
import zipfile

import pytest
from ecd_epitope_scanner.experimental_grasp import (
    ConstraintError,
    FastaChain,
    build_constraint_package,
    format_ir,
    format_rpr,
    parse_constraint,
    select_experimental_ir,
    validate_constraint_against_fasta,
)


def test_grasp_ir_rpr_syntax_and_fasta_validation():
    ir = parse_constraint("1-2-C")
    assert ir == {"kind": "IR", "constraint": "1-2-C"}
    validate_constraint_against_fasta(ir, {1: "ACD", 2: "EFG", 3: "HIK"})
    rpr = format_rpr(1, 2, "C", 2, 3, "G", 8.0, 0.05)
    assert rpr == "1-2-C, 2-3-G, 8.000, 0.05"
    parsed = parse_constraint(rpr)
    validate_constraint_against_fasta(parsed, {1: "ACD", 2: "EFG", 3: "HIK"})
    with pytest.raises(ConstraintError, match="FASTA has D"):
        validate_constraint_against_fasta(ir, {1: "ADD", 2: "EFG", 3: "HIK"})
    with pytest.raises(ConstraintError):
        format_ir(1, 0, "A")
    with pytest.raises(ConstraintError, match="GRASP bins"):
        format_rpr(1, 1, "A", 2, 1, "E", 3.9)


def test_package_has_fixed_chain_order_evidence_and_deterministic_zip():
    kwargs = dict(
        chains=[FastaChain(1, "ECD", "ACD"), FastaChain(2, "VH", "EFG"), FastaChain(3, "VL", "HIK")],
        experimental_ir=[],
        predicted_ir=[{
            "constraint": "2-2-F", "evidence_id": "run-1", "reviewer": "reviewer",
            "rationale": "Manually selected after review", "details": {"score": 0.8},
        }],
        manual_rpr=[],
        merge_sources=["predicted_ir"],
        manifest_context={"limitations": ["RPR unavailable"], "provenance": {"project_id": "p1"}},
    )
    archive_a, manifest = build_constraint_package(**kwargs)
    archive_b, _ = build_constraint_package(**kwargs)
    assert archive_a == archive_b
    assert manifest["status"] == "IR_ONLY"
    assert manifest["rpr_count"] == 0
    with zipfile.ZipFile(io.BytesIO(archive_a)) as bundle:
        assert "1=ECD" not in bundle.namelist()
        assert bundle.read("constraints.fasta").decode().startswith(">1|ECD\nACD\n>2|VH")
        assert bundle.read("grasp_constraints.txt").decode() == "2-2-F\n"
        assert bundle.read("predicted_ir.tsv").decode().splitlines()[1].startswith("2-2-F\tpredicted_ir\trun-1")


def test_package_rejects_wrong_chain_sequence_and_unreviewed_row():
    common = dict(
        chains=[FastaChain(1, "ECD", "ACD"), FastaChain(2, "VH", "EFG"), FastaChain(3, "VL", "HIK")],
        experimental_ir=[], predicted_ir=[], manual_rpr=[], merge_sources=[], manifest_context={},
    )
    with pytest.raises(ConstraintError, match="FASTA has F"):
        build_constraint_package(**(common | {"experimental_ir": [{
            "constraint": "2-2-G", "reviewer": "r", "rationale": "checked"
        }]}))
    with pytest.raises(ConstraintError, match="reviewer and rationale"):
        build_constraint_package(**(common | {"predicted_ir": [{
            "constraint": "2-2-F", "reviewer": "", "rationale": ""
        }]}))


def test_experimental_ir_requires_single_site_and_uses_frozen_mapping():
    effect = {"surface_projection": {"probes": [{
        "probe_id": "single-20", "attribution": "SINGLE_SITE", "mutated_positions": [20],
        "decline_pct": 42.5, "valid_decline_batch_count": 2,
    }]}}
    mapping = [
        {"uniprot_pos": "10", "normalized_position": "1", "one_letter": "A"},
        {"uniprot_pos": "20", "normalized_position": "3", "one_letter": "D"},
    ]
    restraint, details = select_experimental_ir(effect, "single-20", mapping, "ACD")
    assert restraint == "1-3-D"
    assert details["uniprot_position"] == 20
    assert details["decline_pct"] == 42.5
    effect["surface_projection"]["probes"][0]["attribution"] = "SAME_SECTOR_COMBINATION"
    with pytest.raises(ConstraintError, match="SINGLE_SITE_REQUIRED"):
        select_experimental_ir(effect, "single-20", mapping, "ACD")


def test_reverse_rpr_pair_is_duplicate():
    row = {"reviewer": "reviewer", "rationale": "validated pair"}
    with pytest.raises(ConstraintError, match="duplicate constraint"):
        build_constraint_package(
            chains=[FastaChain(1, "ECD", "ACD"), FastaChain(2, "VH", "EFG"), FastaChain(3, "VL", "HIK")],
            experimental_ir=[], predicted_ir=[],
            manual_rpr=[
                row | {"constraint": "1-1-A, 2-1-E, 8.0, 0.05"},
                row | {"constraint": "2-1-E, 1-1-A, 9.0, 0.05"},
            ],
            merge_sources=["manual_rpr"], manifest_context={},
        )
