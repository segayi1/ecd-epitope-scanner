import hashlib
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "cross_inputs", Path(__file__).resolve().parents[1] / "scripts/prepare_cross_protein_validation.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_author_identity_and_frozen_corruption(tmp_path):
    keys = ("pdbx_PDB_model_num", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
            "label_asym_id", "label_seq_id", "label_comp_id")
    values = [("1", "1"), ("Z", "Z"), ("42", "42"), ("?", "A"),
              ("B", "B"), ("7", "8"), ("ALA", "GLY")]
    rows = module.residue_rows({"_atom_site." + k: v for k, v in zip(keys, values)})
    assert len(rows) == 2
    assert {row[3] for row in rows} == {"?", "A"}
    assert all(row[1] == "Z" and row[4] == "B" for row in rows)
    assert len({r[0] for r in module.COHORT if r[3] == "primary"}) == 8
    assert {r[1] for r in module.COHORT if r[0] == "MBP"} == {"holdout"}
    path = tmp_path / "input.cif"
    path.write_bytes(b"original")
    entries = [{"role": r[3]} for r in module.COHORT]
    manifest = {"entries": entries, "files": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()}}
    assert module.verify(tmp_path, manifest)["input_hashes"] == "PASS"
    path.write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="corrupted"):
        module.verify(tmp_path, manifest)
