from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import tempfile
from collections import Counter
from pathlib import Path

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "HID": "H",
    "HIE": "H", "HIP": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
    "TRP": "W", "TYR": "Y", "VAL": "V",
}
MUTATION = re.compile(r"([A-Z])(\d+)([A-Z])")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ca_sequence(path: Path) -> dict[int, str]:
    residues: dict[int, str] = {}
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith("ATOM  ") and line[12:16].strip() == "CA" and line[21] == "A":
            aa = AA3.get(line[17:20].strip())
            if aa is None:
                raise ValueError(f"Unknown chain A residue {line[17:20].strip()} in {path}")
            residues[int(line[22:26])] = aa
    if not residues:
        raise ValueError(f"No chain A C-alpha atoms in {path}")
    return residues


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
        temp = Path(stream.name)
    os.replace(temp, path)


def store_structure(source: Path, expected_hash: str, out_dir: Path) -> str:
    if not source.is_file() or sha256(source) != expected_hash:
        raise ValueError(f"Source structure identity mismatch: {source}")
    relative = f"structures/{expected_hash}.pdb"
    target = out_dir / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if sha256(target) != expected_hash:
            raise ValueError(f"Existing web structure has wrong content: {target}")
    else:
        with tempfile.NamedTemporaryFile("wb", dir=target.parent, delete=False) as stream:
            temp = Path(stream.name)
        shutil.copyfile(source, temp)
        if sha256(temp) != expected_hash:
            temp.unlink(missing_ok=True)
            raise ValueError(f"Copied structure failed hash verification: {source}")
        os.replace(temp, target)
    return relative


def main() -> None:
    parser = argparse.ArgumentParser(description="Export hash-verified CTLA4 production structures for the web viewer.")
    parser.add_argument("--study-root", type=Path)
    parser.add_argument("--mapping-json", type=Path)
    parser.add_argument("--out-dir", type=Path)
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    study = args.study_root or repo.parent / "ctla4_field_fusion_20261006"
    mapping_path = args.mapping_json or repo.parent / "cross_protein_validation_20261006" / "mapping" / "canonical_residue_maps.json"
    out_dir = args.out_dir or repo / "web" / "data" / "ctla4_six_mutants"
    source_manifest = study / "study_manifest.json"
    report = study / "deliverables" / "CTLA4_field_fusion_completion_20261007.md"
    manifest = json.loads(source_manifest.read_text(encoding="utf-8"))
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    if manifest.get("research_profile") != "field_fusion_v1_shadow" or len(manifest["candidate_ids"]) != 6:
        raise ValueError("The frozen six-candidate CTLA4 field-fusion manifest was not found.")

    entity = next(
        item for item in mapping["structures"]["1AH1"]["entity_mappings"]
        if item.get("target") == "CTLA4" and "A" in item.get("chains", [])
    )
    if mapping.get("references", {}).get("CTLA4", {}).get("source", "").split("/")[-1].split(".")[0] != "P16410":
        raise ValueError("Canonical map is not frozen to P16410.")
    reference = {
        int(item["reference_position"]): item["reference_aa"]
        for item in entity["residues"]
        if item.get("reference_position") is not None and item.get("reference_aa")
    }
    canonical_start = 37
    local_offset = canonical_start - 1
    first_production = next(
        row for rows in manifest["pair_manifest"].values() for row in rows
        if row.get("sampling_role") == "production"
    )
    baseline_path = Path(first_production["batch_path"]) / "input.pdb"
    baseline = ca_sequence(baseline_path)
    expected_reference = "".join(reference[position] for position in range(canonical_start, 162))
    baseline_sequence = "".join(baseline[position] for position in range(1, len(baseline) + 1))
    if baseline_sequence != expected_reference:
        raise ValueError("Normalized FoldX input does not exactly match the frozen P16410 reference segment.")
    if sha256(baseline_path) != first_production["input_sha256"]:
        raise ValueError("Normalized FoldX input hash does not match the pair manifest.")

    extension_id = manifest.get("candidate_extension", {}).get("candidate_id")
    candidates = []
    for index, candidate_id in enumerate(manifest["candidate_ids"]):
        label = manifest["candidate_labels"][candidate_id]
        display_label = label.split(" (")[0]
        label_mutations = []
        for token in display_label.split("/"):
            match = MUTATION.fullmatch(token)
            if not match:
                raise ValueError(f"Invalid canonical mutation label: {token}")
            wt, position, mutant = match.group(1), int(match.group(2)), match.group(3)
            local_resi = position - local_offset
            if reference.get(position) != wt or baseline.get(local_resi) != wt:
                raise ValueError(f"Canonical/local WT residue mismatch for {token}.")
            label_mutations.append({
                "canonical": token, "from": wt, "position": position,
                "to": mutant, "local_resi": local_resi,
            })
        expected_local = {item["local_resi"]: (item["from"], item["to"]) for item in label_mutations}
        rows = sorted(
            (row for row in manifest["pair_manifest"][candidate_id] if row.get("sampling_role") == "production"),
            key=lambda row: int(row["sampling_slot"]),
        )
        role_counts = Counter(row.get("sampling_role") for row in manifest["pair_manifest"][candidate_id])
        if len(rows) != 5 or [int(row["sampling_slot"]) for row in rows] != list(range(5)):
            raise ValueError(f"{label}: expected exactly production slots 0-4.")
        if role_counts != Counter({"production": 5, "null_train": 10, "null_holdout": 10}):
            raise ValueError(f"{label}: production/null role counts do not match the frozen contract.")

        models = []
        wt_hashes, mut_hashes, ddgs = set(), set(), []
        protocols = set()
        for row in rows:
            wt_path, mut_path = Path(row["wt_path"]), Path(row["mut_path"])
            wt_hash, mut_hash = row["wt_sha256"], row["mut_sha256"]
            wt_file = store_structure(wt_path, wt_hash, out_dir)
            mut_file = store_structure(mut_path, mut_hash, out_dir)
            wt_hashes.add(wt_hash)
            mut_hashes.add(mut_hash)
            wt_seq, mut_seq = ca_sequence(wt_path), ca_sequence(mut_path)
            if wt_seq != baseline:
                raise ValueError(f"{label}: WT chain A sequence changed in slot {row['sampling_slot']}.")
            differences = {
                position for position in set(wt_seq) | set(mut_seq)
                if wt_seq.get(position) != mut_seq.get(position)
            }
            if differences != set(expected_local):
                raise ValueError(f"{label}: unexpected sequence changes in slot {row['sampling_slot']}.")
            for position, (wt_aa, mut_aa) in expected_local.items():
                if wt_seq.get(position) != wt_aa or mut_seq.get(position) != mut_aa:
                    raise ValueError(f"{label}: mutation identity mismatch at local residue {position}.")
            row_mutations = {}
            for token in row["mutations"]:
                match = MUTATION.fullmatch(token)
                if not match:
                    raise ValueError(f"Invalid FoldX mutation record: {token}")
                row_mutations[int(match.group(2))] = (match.group(1), match.group(3))
            if row_mutations != expected_local:
                raise ValueError(f"{label}: candidate label does not match pair-manifest mutations.")
            ddg = float(row["ddg_kcal_mol"])
            if not math.isfinite(ddg):
                raise ValueError(f"{label}: non-finite FoldX score.")
            ddgs.append(ddg)
            protocols.add(row.get("sampling_protocol"))
            models.append({
                "slot": int(row["sampling_slot"]) + 1,
                "foldx_slot": int(row["sampling_slot"]),
                "wt_file": wt_file,
                "mut_file": mut_file,
                "wt_sha256": wt_hash,
                "mut_sha256": mut_hash,
                "foldx_ddg_kcal_mol": ddg,
            })
        if len(protocols) != 1:
            raise ValueError(f"{label}: mixed sampling protocols.")
        candidates.append({
            "id": candidate_id,
            "label": label,
            "mutations": label_mutations,
            "role": "coverage_extension" if candidate_id == extension_id else "original_five_panel",
            "production_models": models,
            "unique_wt_files": len(wt_hashes),
            "unique_mutant_files": len(mut_hashes),
            "foldx_ddg_median_kcal_mol": statistics.median(ddgs),
            "foldx_ddg_range_kcal_mol": [min(ddgs), max(ddgs)],
            "sampling_protocol": next(iter(protocols)),
        })

    payload = {
        "schema_version": "ctla4_structure_viewer_v1",
        "study_id": manifest["study_id"],
        "research_profile": manifest["research_profile"],
        "protein": {
            "name": "CTLA-4 extracellular domain",
            "accession": "P16410",
            "structure": "PDB 1AH1",
            "chain": "A",
            "reference_segment": [37, 161],
            "local_to_canonical_offset": local_offset,
            "sequence_identity_verified": True,
        },
        "source_manifest_sha256": sha256(source_manifest),
        "canonical_mapping_sha256": sha256(mapping_path),
        "completion_report": {
            "file": "CTLA4_field_fusion_completion_20261007.md",
            "sha256": sha256(report),
        },
        "interpretation": {
            "role": "research structure display; does not grant or change formal eligibility",
            "replicates": "FoldX side-chain sampling slots, not independent biological or backbone replicates",
            "unknowns": "Research field-fusion report retains unresolved surface regions as unknown",
        },
        "candidates": candidates,
    }
    if report.is_file():
        report_target = out_dir / "CTLA4_field_fusion_completion_20261007.md"
        if not report_target.exists():
            shutil.copyfile(report, report_target)
        if sha256(report_target) != payload["completion_report"]["sha256"]:
            raise ValueError("Copied completion report does not match its frozen hash.")
    atomic_json(out_dir / "manifest.json", payload)
    print(json.dumps({
        "output": str(out_dir),
        "candidate_count": len(candidates),
        "production_pair_count": sum(len(item["production_models"]) for item in candidates),
        "unique_wt_contents": {item["label"]: item["unique_wt_files"] for item in candidates},
        "unique_mutant_contents": {item["label"]: item["unique_mutant_files"] for item in candidates},
        "manifest_sha256": sha256(out_dir / "manifest.json"),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
