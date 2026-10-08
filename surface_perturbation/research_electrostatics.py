"""Input preparation checks for the shadow APBS protocol."""
from __future__ import annotations

import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np

from ..util import sha256_file, write_json
from .field_fusion_research import PREPARATION_PROTOCOL
from .tools import fingerprint, read_pqr


def _heavy_atom_name(name: str) -> bool:
    return not name.strip().lstrip("0123456789").upper().startswith("H")


def validate_preserved_pqr(atom_map: dict, residue_ids: dict, pqr: dict) -> dict:
    """Require exact heavy-atom identity/coordinates and complete finite PQR parameters."""
    expected = {}
    for key, value in atom_map.items():
        if not _heavy_atom_name(key[3]):
            continue
        pqr_key = (key[0], int(residue_ids[key[:3]]), key[3])
        expected[pqr_key] = value
    actual = {key: value for key, value in pqr.items() if _heavy_atom_name(key[2])}
    missing = sorted(set(expected) - set(actual))
    added = sorted(set(actual) - set(expected))
    mismatched_residues = sorted(key for key in set(expected) & set(actual)
                                 if expected[key][1] != actual[key]["residue"])
    displacements = [math.dist(expected[key][0], actual[key]["xyz"])
                     for key in set(expected) & set(actual)]
    invalid_parameters = sorted(key for key, value in pqr.items()
        if not all(math.isfinite(float(value[field])) for field in ("charge", "radius"))
        or float(value["radius"]) < 0)
    heavy_radii = [float(value["radius"]) for key, value in actual.items()]
    status = ("MATCHED" if not missing and not added and not mismatched_residues and
              not invalid_parameters and displacements and max(displacements) <= .002 and
              all(radius > 0 for radius in heavy_radii) else "UNKNOWN_PREPARATION_MISMATCH")
    return {"status": status, "expected_heavy_atoms": len(expected), "observed_heavy_atoms": len(actual),
            "missing_heavy_atoms": [list(key) for key in missing],
            "added_heavy_atoms": [list(key) for key in added],
            "residue_name_mismatches": [list(key) for key in mismatched_residues],
            "invalid_parameter_atoms": [list(key) for key in invalid_parameters],
            "max_heavy_coordinate_displacement_a": max(displacements, default=None),
            "pH_assignment": "PDB2PQR PARSE at pH 7.4; explicit hydrogen optimization and debumping disabled"}


def prepare_preserved_model(atom_map: dict, residue_ids: dict, tools, folder: Path) -> dict:
    """Run/cached PDB2PQR preparation without moving existing heavy atoms."""
    from .electrostatics import write_atoms

    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    model_digest = hashlib.sha256()
    for key in sorted(atom_map):
        xyz, residue, element = atom_map[key]
        model_digest.update(repr(key).encode()); model_digest.update(repr(residue_ids[key[:3]]).encode())
        model_digest.update(np.asarray(xyz, float).tobytes())
        model_digest.update(str(residue).encode()); model_digest.update(str(element).encode())
    model_sha = model_digest.hexdigest()
    pdb_path = folder / f"source_{model_sha}.pdb"
    if not pdb_path.exists():
        write_atoms(pdb_path, atom_map, residue_ids)
    source_sha = sha256_file(pdb_path)
    pqr_path, pqr = tools.prepare_pqr(pdb_path, protocol=PREPARATION_PROTOCOL)
    audit = validate_preserved_pqr(atom_map, residue_ids, pqr)
    pqr_sha = sha256_file(pqr_path)
    portable = folder / "prepared_pqr" / f"{pqr_sha}.pqr"
    portable.parent.mkdir(parents=True, exist_ok=True)
    if portable.exists():
        if sha256_file(portable) != pqr_sha:
            raise ValueError("prepared PQR content-addressed object is corrupt")
    else:
        shutil.copyfile(pqr_path, portable)
    contract = {"preparation_protocol": PREPARATION_PROTOCOL, "source_atom_map_sha256": model_sha,
                "source_pdb_sha256": source_sha, "pdb2pqr_sha256": tools.identity(tools.pqr),
                "prepared_pqr_sha256": pqr_sha, "force_field": "PARSE", "pH": 7.4,
                "optimize_hydrogens": False, "debump": False}
    write_json(folder / f"preparation_{model_sha}.json", {"contract": contract, "audit": audit})
    return {"path": portable, "cache_path": Path(pqr_path), "atoms": pqr, "audit": audit,
            "contract": contract}


def parameter_signature(pqr: dict) -> str:
    rows = [(key, value["residue"], round(float(value["charge"]), 8), round(float(value["radius"]), 8))
            for key, value in sorted(pqr.items())]
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()


def validate_pair_parameters(wt: dict, mut: dict, mutated_author: set[tuple[str, int]]) -> dict:
    from .electrostatics import protonation_audit
    audit = protonation_audit(wt, mut, mutated_author)
    audit["status"] = "MATCHED" if not audit["nonmutated_parameter_changes"] else "UNKNOWN_NONMUTATED_PARAMETER_CHANGE"
    return audit
