"""Freeze cross-protein research inputs; no model solving or canonical mapping."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

COHORT = [
    ("CTLA4", "development", "1AH1", "primary"),
    ("ubiquitin", "development", "1UBQ", "primary"),
    ("PD1", "confirmation", "3RRQ", "primary"),
    ("PDL1", "confirmation", "4Z18", "primary"),
    ("B7H6", "confirmation", "3PV7", "primary"),
    ("hen_lysozyme", "confirmation", "1LYZ", "primary"),
    ("MBP", "holdout", "1OMP", "primary"),
    ("T4_lysozyme", "holdout", "2LZM", "primary"),
    ("ubiquitin", "development", "1D3Z", "NMR_stress"),
    ("PD1_PDL1", "confirmation", "4ZQK", "complex_engineered_stress"),
    ("B7H6", "confirmation", "3PV6", "complex_stress"),
    ("MBP", "holdout", "1ANF", "ligand_stress"),
]


def residue_rows(data):
    keys = ("pdbx_PDB_model_num", "auth_asym_id", "auth_seq_id",
            "pdbx_PDB_ins_code", "label_asym_id", "label_seq_id", "label_comp_id")
    columns = [data["_atom_site." + key] for key in keys]
    if len({len(col) for col in columns}) != 1:
        raise ValueError("inconsistent atom identity columns")
    # Preserve author identifiers and insertion codes, including non-polymers.
    return sorted(set(zip(*columns)))


def protein_chains(data):
    types = data.get('_entity_poly.type', [])
    chains = data.get('_entity_poly.pdbx_strand_id', [])
    if len(types) != len(chains):
        raise ValueError('polymer type/author chain metadata mismatch')
    return sorted({chain.strip() for kind, group in zip(types, chains)
                   if kind in ('polypeptide(L)', 'polypeptide(D)')
                   for chain in group.split(',') if chain.strip()})


def verify(root, manifest):
    for name, expected in manifest["files"].items():
        actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError("frozen input corrupted: " + name)
    assert len(manifest["entries"]) == 12
    assert len([e for e in manifest["entries"] if e["role"] == "primary"]) == 8
    return {"input_hashes": "PASS", "files": len(manifest["files"]),
            "scientific_comparisons": "NOT_RUN", "canonical_mapping": "REVIEW_REQUIRED"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "input_manifest.json"
    if manifest_path.exists():
        print(json.dumps(verify(root, json.loads(manifest_path.read_text(encoding="utf-8")))))
        return
    entries = []
    for protein, split, pdb, role in COHORT:
        source = f"https://files.rcsb.org/download/{pdb}.cif"
        target = root / f"{pdb}.cif"
        if not target.exists():
            with urlopen(source, timeout=60) as response:
                content = response.read()
            partial = target.with_suffix(".download")
            partial.write_bytes(content)
            partial.replace(target)
        data = MMCIF2Dict(str(target))
        if data["_entry.id"][0].upper() != pdb:
            raise ValueError("wrong deposited identity: " + pdb)
        rows = residue_rows(data)
        with (root / f"{pdb}_residues.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("model", "author_chain", "author_residue", "insertion_code",
                             "label_chain", "label_residue", "residue_name"))
            writer.writerows(rows)
        metadata = {key: value for key, value in data.items() if key.startswith(
            ("_entity.", "_entity_poly.", "_struct.", "_exptl.", "_struct_ref.", "_struct_ref_seq."))}
        (root / f"{pdb}_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        entries.append(dict(protein=protein, split=split, pdb=pdb, role=role, source=source,
                            models=len({row[0] for row in rows}),
                            author_chains=sorted({row[1] for row in rows}),
                            protein_author_chains=protein_chains(data),
                            method=data.get("_exptl.method", []),
                            status="INPUT_PARSED_MAPPING_AND_CONTEXT_REVIEW_REQUIRED"))
        print(pdb + ": input parsed", flush=True)
    files = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
             for path in sorted(root.iterdir()) if path.suffix in {".cif", ".csv", ".json"}}
    manifest = dict(schema="cross_protein_inputs_v1", entries=entries, files=files,
                    preparation_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    scientific_status="NOT_RUN", scope="deposited inputs, not validated recommendations")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(verify(root, manifest)))


if __name__ == "__main__":
    main()
