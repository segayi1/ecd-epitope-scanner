"""Sequence validation for experiment constructs; independent of HTTP and storage."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

from Bio.Align import PairwiseAligner

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


class ConstructError(ValueError):
    pass


def reference_positions(fasta_path: Path, mapping_path: Path) -> tuple[int, set[int]]:
    header = fasta_path.read_text(encoding="utf-8-sig").splitlines()[0]
    match = re.search(r":(\d+)-(\d+)\s*$", header)
    with mapping_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    mapped = {
        int(row["uniprot_pos"])
        for row in rows
        if row.get("uniprot_pos") and row.get("chain_id") and row.get("auth_seq_id")
    }
    if match:
        return int(match.group(1)), mapped
    first = next(
        (int(row["uniprot_pos"]) for row in rows if row.get("normalized_position") == "1"), None
    )
    if first is None:
        raise ConstructError("Cannot determine reference position numbering")
    return first, mapped


def validate_antigen_construct(
    reference: str,
    position_start: int,
    mapped_positions: set[int],
    mutations: list[dict[str, Any]] | None = None,
    full_sequence: str | None = None,
) -> dict[str, Any]:
    reference = reference.upper()
    if not reference or set(reference) - AMINO_ACIDS:
        raise ConstructError("Reference sequence contains invalid residues")
    if bool(mutations) == bool(full_sequence):
        raise ConstructError("Provide either mutations or full_sequence")
    differences: list[dict[str, Any]] = []
    if mutations:
        seen: set[int] = set()
        sequence = list(reference)
        for item in mutations:
            position = int(item["position"])
            wt = str(item["wt"]).upper()
            mutant = str(item["mut"]).upper()
            if position in seen:
                raise ConstructError(f"Duplicate mutation position {position}")
            seen.add(position)
            offset = position - position_start
            if offset < 0 or offset >= len(reference):
                raise ConstructError(f"Mutation position {position} is outside the reference")
            if len(wt) != 1 or wt not in AMINO_ACIDS or wt != reference[offset]:
                raise ConstructError(
                    f"WT residue mismatch at {position}: expected {reference[offset]}"
                )
            if len(mutant) != 1 or mutant not in AMINO_ACIDS or mutant == wt:
                raise ConstructError(f"Invalid replacement at {position}")
            sequence[offset] = mutant
            differences.append(
                {
                    "kind": "SUBSTITUTION",
                    "position": position,
                    "wt": wt,
                    "mut": mutant,
                    "mapped": position in mapped_positions,
                }
            )
        variant = "".join(sequence)
    else:
        variant = "".join(str(full_sequence).split()).upper()
        if not variant or set(variant) - AMINO_ACIDS:
            raise ConstructError("Full sequence contains invalid residues")
        aligner = PairwiseAligner(
            mode="global",
            match_score=2,
            mismatch_score=-1,
            open_gap_score=-5,
            extend_gap_score=-0.5,
        )
        alignment = aligner.align(reference, variant)[0]
        ref_coordinates, var_coordinates = alignment.coordinates
        for index in range(len(ref_coordinates) - 1):
            r0, r1 = int(ref_coordinates[index]), int(ref_coordinates[index + 1])
            v0, v1 = int(var_coordinates[index]), int(var_coordinates[index + 1])
            if r1 > r0 and v1 > v0:
                for offset in range(r1 - r0):
                    if reference[r0 + offset] != variant[v0 + offset]:
                        position = position_start + r0 + offset
                        differences.append(
                            {
                                "kind": "SUBSTITUTION",
                                "position": position,
                                "wt": reference[r0 + offset],
                                "mut": variant[v0 + offset],
                                "mapped": position in mapped_positions,
                            }
                        )
            elif r1 > r0:
                for offset in range(r0, r1):
                    position = position_start + offset
                    differences.append(
                        {
                            "kind": "DELETION",
                            "position": position,
                            "wt": reference[offset],
                            "mut": None,
                            "mapped": position in mapped_positions,
                        }
                    )
            elif v1 > v0:
                differences.append(
                    {
                        "kind": "INSERTION",
                        "after_position": position_start + r0 - 1,
                        "residues": variant[v0:v1],
                        "mapped": False,
                    }
                )
    if not differences and full_sequence and variant == reference:
        return {
            "sequence": variant,
            "mutations": [],
            "mapping_status": "MAPPED",
            "risk_flags": [],
        }
    if not differences:
        raise ConstructError("Construct is identical to the reference")
    flags = []
    if any(not item["mapped"] for item in differences):
        flags.append("STRUCTURE_UNMAPPED")
    if any(item["kind"] != "SUBSTITUTION" for item in differences):
        flags.append("INSERTION_OR_DELETION")
    return {
        "sequence": variant,
        "mutations": differences,
        "mapping_status": "PARTIAL_OR_UNMAPPED" if "STRUCTURE_UNMAPPED" in flags else "MAPPED",
        "risk_flags": flags,
    }


def validate_antibody_construct(
    vh: str,
    vl: str,
    parent_vh: str | None = None,
    parent_vl: str | None = None,
) -> dict[str, Any]:
    vh, vl = vh.strip().upper(), vl.strip().upper()
    if not vh or not vl or set(vh + vl) - AMINO_ACIDS:
        raise ConstructError("VH and VL must contain standard amino acids")
    if parent_vh is None and parent_vl is None:
        return {
            "vh_sequence": vh,
            "vl_sequence": vl,
            "heavy_mutations": [],
            "light_mutations": [],
            "mutation_class": "PARENT",
        }
    if (
        parent_vh is None
        or parent_vl is None
        or len(vh) != len(parent_vh)
        or len(vl) != len(parent_vl)
    ):
        raise ConstructError("Variant VH/VL lengths must match its parent")
    heavy = [
        {"position": index, "wt": before, "mut": after}
        for index, (before, after) in enumerate(zip(parent_vh, vh), 1)
        if before != after
    ]
    light = [
        {"position": index, "wt": before, "mut": after}
        for index, (before, after) in enumerate(zip(parent_vl, vl), 1)
        if before != after
    ]
    count = len(heavy) + len(light)
    if count == 0:
        raise ConstructError("Variant is identical to parent")
    return {
        "vh_sequence": vh,
        "vl_sequence": vl,
        "heavy_mutations": heavy,
        "light_mutations": light,
        "mutation_class": "SINGLE" if count == 1 else "COMBINATION",
    }
