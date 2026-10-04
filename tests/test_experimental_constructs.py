from __future__ import annotations

import pytest
from ecd_epitope_scanner.experimental_constructs import (
    ConstructError,
    validate_antibody_construct,
    validate_antigen_construct,
)


def test_antigen_substitutions_validate_wt_duplicates_and_partial_mapping():
    result = validate_antigen_construct(
        "ACDE",
        19,
        {19, 20, 21},
        [{"position": 20, "wt": "C", "mut": "A"}, {"position": 22, "wt": "E", "mut": "D"}],
    )
    assert result["sequence"] == "AADD"
    assert result["mapping_status"] == "PARTIAL_OR_UNMAPPED"
    assert [item["mapped"] for item in result["mutations"]] == [True, False]
    with pytest.raises(ConstructError, match="Duplicate"):
        validate_antigen_construct(
            "ACDE",
            19,
            {20},
            [{"position": 20, "wt": "C", "mut": "A"}, {"position": 20, "wt": "C", "mut": "D"}],
        )
    with pytest.raises(ConstructError, match="WT residue mismatch"):
        validate_antigen_construct("ACDE", 19, {20}, [{"position": 20, "wt": "A", "mut": "D"}])


def test_full_sequence_reports_substitution_and_indel():
    changed = validate_antigen_construct("ACDE", 19, {19, 20, 21, 22}, full_sequence="ACDF")
    assert changed["mutations"] == [
        {"kind": "SUBSTITUTION", "position": 22, "wt": "E", "mut": "F", "mapped": True}
    ]
    inserted = validate_antigen_construct("ACDE", 19, {19, 20, 21, 22}, full_sequence="ACDGE")
    assert any(item["kind"] == "INSERTION" for item in inserted["mutations"])


def test_full_reference_sequence_can_be_registered_as_wild_type_control():
    result = validate_antigen_construct("ACD", 10, {10, 11, 12}, full_sequence="ACD")
    assert result["sequence"] == "ACD"
    assert result["mutations"] == []
    assert result["mapping_status"] == "MAPPED"


def test_antibody_parent_and_variant_track_chain_positions():
    parent = validate_antibody_construct("ACD", "EFG")
    assert parent["mutation_class"] == "PARENT"
    variant = validate_antibody_construct("AAD", "EFG", "ACD", "EFG")
    assert variant["mutation_class"] == "SINGLE"
    assert variant["heavy_mutations"] == [{"position": 2, "wt": "C", "mut": "A"}]
    with pytest.raises(ConstructError, match="lengths"):
        validate_antibody_construct("AA", "EFG", "ACD", "EFG")
