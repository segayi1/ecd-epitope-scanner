from __future__ import annotations

import copy
import csv
from pathlib import Path

from Bio import __version__ as biopython_version
from Bio.PDB import PDBIO, Chain, MMCIFParser, Model, PDBParser, Structure
from Bio.PDB.Polypeptide import protein_letters_3to1

from ..models import ProteinContext, ResidueRecord
from ..structure_analysis import MAX_ASA, TIEN_DSSP_RADII, calculate_sasa, load_structure
from .models import ResidueKey


def residue_keys(records: list[ResidueRecord], model_id: str = "0") -> dict[int, ResidueKey]:
    seen: set[tuple[str, int, str]] = set()
    result: dict[int, ResidueKey] = {}
    for item in records:
        key = (item.chain_id, item.auth_seq_id, item.insertion_code or "")
        if key in seen or item.uniprot_pos in result:
            raise ValueError("CHAIN_MAPPING_AMBIGUOUS")
        seen.add(key)
        result[item.uniprot_pos] = ResidueKey(
            model_id=model_id,
            chain_id=item.chain_id,
            auth_seq_id=item.auth_seq_id,
            insertion_code=item.insertion_code or "",
            label_seq_id=item.label_seq_id,
            uniprot_pos=item.uniprot_pos,
            one_letter=item.aa,
        )
    return result


def validate_mapping(records: list[ResidueRecord], protein: ProteinContext) -> list[dict]:
    errors = []
    for item in records:
        if not item.aa or item.aa not in "ACDEFGHIKLMNPQRSTVWY":
            errors.append({"code": "WT_RESIDUE_MISMATCH", "position": item.uniprot_pos})
        if protein.sequence[item.uniprot_pos - 1 : item.uniprot_pos] != item.aa:
            errors.append({"code": "WT_RESIDUE_MISMATCH", "position": item.uniprot_pos})
    return errors


def write_mapping_sidecar(
    path: Path, records: list[ResidueRecord], keys: dict[int, ResidueKey]
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "normalized_position",
                "uniprot_pos",
                "chain_id",
                "auth_seq_id",
                "insertion_code",
                "label_seq_id",
                "one_letter",
            ],
        )
        writer.writeheader()
        for normalized_position, item in enumerate(records, 1):
            key = keys[item.uniprot_pos]
            writer.writerow(
                {
                    "normalized_position": normalized_position,
                    "uniprot_pos": item.uniprot_pos,
                    "chain_id": key.chain_id,
                    "auth_seq_id": key.auth_seq_id,
                    "insertion_code": key.insertion_code,
                    "label_seq_id": key.label_seq_id or "",
                    "one_letter": key.one_letter,
                }
            )
    return path


def context_manifest(protein: ProteinContext, records: list[ResidueRecord], chain_id: str) -> dict:
    return {
        "prediction_context": {
            "chain_id": chain_id,
            "residue_count": len(records),
            "mode": "target_chain",
        },
        "biological_context": {"mode": "provided_structure", "assembly_status": "unknown"},
        "ecd_range": [protein.ecd_start, protein.ecd_end],
        "mapping_status": "validated",
    }


def write_normalized_prediction_chain(
    structure_path: Path | None,
    output_path: Path,
    records: list[ResidueRecord],
    chain_id: str,
) -> Path | None:
    """Write a deterministic single-chain, continuously numbered PDB.

    This is the prediction context.  The original input remains the biological
    context and is never silently replaced by this normalized file.
    """
    if structure_path is None or not structure_path.is_file():
        return None
    parser = (
        MMCIFParser(QUIET=True)
        if structure_path.suffix.lower() in {".cif", ".mmcif"}
        else PDBParser(QUIET=True)
    )
    try:
        source = parser.get_structure("source", str(structure_path))
        source_chain = next(model[chain_id] for model in source if chain_id in model)
    except (KeyError, StopIteration, OSError, ValueError) as exc:
        raise ValueError("CHAIN_MAPPING_AMBIGUOUS") from exc
    source_by_id = {
        (int(residue.id[1]), str(residue.id[2]).strip()): residue
        for residue in source_chain
        if residue.id[0] == " "
    }
    target = Chain.Chain(chain_id)
    for item in records:
        residue = source_by_id.get((item.auth_seq_id, item.insertion_code or ""))
        if residue is None:
            raise ValueError(f"RESIDUE_COUNT_MISMATCH:{item.uniprot_pos}")
        observed = protein_letters_3to1.get(str(residue.resname).upper(), "X")
        if observed != item.aa:
            raise ValueError(f"WT_RESIDUE_MISMATCH:{item.uniprot_pos}")
        clone = copy.deepcopy(residue)
        clone.detach_parent()
        clone.id = (clone.id[0], len(target.child_list) + 1, " ")
        target.add(clone)
    if len(target.child_list) != len(records):
        raise ValueError("RESIDUE_COUNT_MISMATCH")
    structure = Structure.Structure("normalized_prediction")
    model = Model.Model(0)
    model.add(target)
    structure.add(model)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))
    return output_path


def normalized_prediction_surface_rows(
    structure_path: Path,
    chain_id: str,
    records: list[ResidueRecord],
) -> tuple[dict[int, dict[str, float]], dict[str, object]]:
    """Measure SASA on the exact normalized single-chain prediction input.

    The pipeline's native SASA retains all standard residues from the selected
    source chain, including engineered substitutions and residues outside the
    mapped ECD. Prediction workers receive the normalized canonical ECD chain.
    This second reference separates structural-context changes from differences
    between the native and external SASA implementations.
    """
    parsed = load_structure(structure_path)
    if chain_id not in parsed.model:
        raise ValueError("NORMALIZED_PREDICTION_CHAIN_MISSING")
    observed = "".join(
        protein_letters_3to1.get(residue.resname.upper(), "X")
        for residue in parsed.model[chain_id]
        if residue.id[0] == " "
    )
    expected = "".join(item.aa for item in records)
    if observed != expected:
        raise ValueError("NORMALIZED_PREDICTION_SEQUENCE_MISMATCH")

    measured, native_method = calculate_sasa(parsed, structure_path.parent, chain_id)
    if len(measured) != len(records):
        raise ValueError("NORMALIZED_PREDICTION_RESIDUE_COUNT_MISMATCH")
    rows: dict[int, dict[str, float]] = {}
    for normalized_position, item in enumerate(records, start=1):
        total_sasa = measured.get((chain_id, normalized_position, ""))
        if total_sasa is None:
            raise ValueError(f"NORMALIZED_PREDICTION_RESIDUE_MISSING:{normalized_position}")
        absolute_sasa = float(total_sasa[0])
        rows[item.uniprot_pos] = {
            "sasa": absolute_sasa,
            "rsasa": absolute_sasa / MAX_ASA.get(item.aa, 200.0),
        }
    return rows, {
        "algorithm": "Bio.PDB Shrake-Rupley",
        "biopython_version": biopython_version,
        "probe_radius_a": 1.4,
        "sphere_points": 960,
        "radii_reference": "Tien et al. 2013/DSSP",
        "atom_radii_a": dict(TIEN_DSSP_RADII),
        "scope": "normalized mapped canonical ECD; no partners, engineered mismatches, or unmapped flanks",
        "method_description": native_method,
    }
