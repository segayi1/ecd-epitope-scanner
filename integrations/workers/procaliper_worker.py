from __future__ import annotations

import hashlib
import importlib.metadata
import json
import sys
import tempfile
import time
import traceback
from collections import Counter
from pathlib import Path

from Bio import __version__ as biopython_version
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley

MAX_ASA = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "Q": 225.0,
    "E": 223.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}


def _residue_sequence(path: Path, chain_id: str) -> str:
    structure = PDBParser(QUIET=True).get_structure("prediction", str(path))
    model = next(structure.get_models())
    if chain_id not in model:
        raise RuntimeError(f"ProCaliper prediction chain {chain_id!r} was not found")
    letters = []
    for residue in model[chain_id]:
        if residue.id[0] != " " or not is_aa(residue, standard=True):
            continue
        name = residue.resname.strip().upper()
        if name not in protein_letters_3to1:
            raise RuntimeError(f"Unsupported residue {name} in ProCaliper input")
        letters.append(protein_letters_3to1[name])
    return "".join(letters)


class _ProteinChainSelect(Select):
    """Keep only standard amino-acid residues from one target chain."""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_model(self, model) -> bool:
        return model.id == 0

    def accept_chain(self, chain) -> bool:
        return chain.id == self.chain_id

    def accept_residue(self, residue) -> bool:
        return residue.id[0] == " " and is_aa(residue, standard=True)


def _write_protein_only_structure(
    source_path: Path, output_path: Path, chain_id: str
) -> dict[str, object]:
    """Strip solvent, ions, ligands and partner chains before SASA calculation."""
    structure = PDBParser(QUIET=True).get_structure("sasa_source", str(source_path))
    model = next(structure.get_models(), None)
    if model is None or chain_id not in model:
        raise RuntimeError(f"Target protein chain {chain_id!r} was not found")

    protein_residues = [
        residue
        for residue in model[chain_id]
        if residue.id[0] == " " and is_aa(residue, standard=True)
    ]
    if not protein_residues:
        raise RuntimeError(f"Target chain {chain_id!r} has no standard protein residues")

    partner_protein_chains = [
        chain
        for chain in model
        if chain.id != chain_id
        and any(is_aa(residue, standard=True) for residue in chain)
    ]
    nonstandard_residues = [
        residue
        for chain in model
        for residue in chain
        if not (residue.id[0] == " " and is_aa(residue, standard=True))
    ]
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), _ProteinChainSelect(chain_id))
    return {
        "retained_protein_residue_count": len(protein_residues),
        "excluded_nonstandard_or_hetero_residue_count": len(nonstandard_residues),
        "excluded_nonstandard_or_hetero_residue_names": dict(
            Counter(residue.resname.strip().upper() for residue in nonstandard_residues)
        ),
        "excluded_partner_protein_chain_count": len(partner_protein_chains),
    }


def main() -> int:
    started = time.perf_counter()
    try:
        request = json.load(sys.stdin)
        from procaliper.protein_structure.sasa import (
            N_POINTS,
            PROBE_RADIUS,
            calculate_sasa,
        )
        default_radii = ShrakeRupley(
            probe_radius=PROBE_RADIUS, n_points=N_POINTS, radii_dict=None
        ).radii_dict
        atom_radii = {
            element: float(default_radii[element]) for element in ("C", "N", "O", "S")
        }

        source_structure_path = Path(request["structure_path"])
        chain_id = str(request["chain_id"])
        requested_sequence = str(request["sequence"])
        with tempfile.TemporaryDirectory(prefix="procaliper_protein_only_") as temp_dir:
            structure_path = Path(temp_dir) / "target_protein_chain.pdb"
            calculation_scope = _write_protein_only_structure(
                source_structure_path, structure_path, chain_id
            )
            sequence = _residue_sequence(structure_path, chain_id)
            if sequence != requested_sequence:
                raise RuntimeError(
                    "ProCaliper input sequence does not match the normalized prediction chain"
                )
            sasa = calculate_sasa(str(structure_path))["all_sasa_value"]
            if len(sasa) != len(sequence):
                raise RuntimeError(
                    f"ProCaliper returned {len(sasa)} residues for sequence length {len(sequence)}"
                )
            calculation_structure_sha256 = hashlib.sha256(
                structure_path.read_bytes()
            ).hexdigest()
        records = [
            {
                "position": position,
                "wt": aa,
                "sasa": float(value),
                "rsasa": float(value) / MAX_ASA[aa],
            }
            for position, (aa, value) in enumerate(zip(sequence, sasa, strict=True), 1)
        ]
        print(
            json.dumps(
                {
                    "status": "success",
                    "records": records,
                    "warnings": [],
                    "errors": [],
                    "provenance": {
                        "worker": "procaliper_worker",
                        "package": "procaliper",
                        "version": importlib.metadata.version("procaliper"),
                        "device": "cpu",
                        "algorithm": "Bio.PDB ShrakeRupley via ProCaliper",
                        "scientific_role": "SASA_CROSSCHECK_ONLY",
                        "capabilities": ["residue_sasa", "relative_sasa"],
                        "not_computed": ["finite_probe_accessibility", "steric_probe_reach", "residue_distance_network"],
                        "biopython_version": biopython_version,
                        "probe_radius_a": PROBE_RADIUS,
                        "sphere_points": N_POINTS,
                        "rsasa_reference": "Tien et al. 2013 maximum ASA",
                        "atom_radii_a": atom_radii,
                        "radii_reference": "Bio.PDB default Shrake-Rupley radii",
                        "calculation_scope": "single_target_standard_protein_chain_only",
                        **calculation_scope,
                        "source_structure_sha256": hashlib.sha256(
                            source_structure_path.read_bytes()
                        ).hexdigest(),
                        "input_structure_sha256": calculation_structure_sha256,
                        "chain_id": chain_id,
                        "sequence_sha256": hashlib.sha256(
                            requested_sequence.encode("ascii")
                        ).hexdigest(),
                        "residue_count": len(sequence),
                    },
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
        return 0
    except ModuleNotFoundError as exc:
        if exc.name != 'procaliper':
            raise
        print(json.dumps({"status": "failed", "records": [], "warnings": [],
            "errors": [{"code": "TOOL_NOT_CONFIGURED", "message": "ProCaliper is not installed"}],
            "provenance": {"worker": "procaliper_worker"},
            "runtime_seconds": time.perf_counter() - started}))
        return 0
    except Exception as exc:
        traceback.print_exc(file=sys.stderr)
        print(
            json.dumps(
                {
                    "status": "failed",
                    "records": [],
                    "warnings": [],
                    "errors": [{"code": "TOOL_CRASHED", "message": str(exc)}],
                    "provenance": {"worker": "procaliper_worker"},
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
