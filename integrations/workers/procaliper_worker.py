from __future__ import annotations

import importlib.metadata
import json
import sys
import time
import traceback
from pathlib import Path

from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

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


def main() -> int:
    started = time.perf_counter()
    try:
        request = json.load(sys.stdin)
        from procaliper.protein_structure.sasa import (
            N_POINTS,
            PROBE_RADIUS,
            calculate_sasa,
        )

        structure_path = Path(request["structure_path"])
        sequence = _residue_sequence(structure_path, str(request["chain_id"]))
        if sequence != request["sequence"]:
            raise RuntimeError(
                "ProCaliper input sequence does not match the normalized prediction chain"
            )
        sasa = calculate_sasa(str(structure_path))["all_sasa_value"]
        if len(sasa) != len(sequence):
            raise RuntimeError(
                f"ProCaliper returned {len(sasa)} residues for sequence length {len(sequence)}"
            )
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
                        "probe_radius_a": PROBE_RADIUS,
                        "sphere_points": N_POINTS,
                        "rsasa_reference": "Tien et al. 2013 maximum ASA",
                    },
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
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
