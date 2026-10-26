#!/usr/bin/env python
"""Independent integrity and mapping check for the frozen benchmark preflight."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "workspace_data/p0_p3_upgrade_20261016/benchmark"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve_repo_path(value: str) -> Path:
    path = (ROOT / value).resolve()
    if path != ROOT and ROOT not in path.parents:
        raise ValueError(f"Path escapes repository: {value}")
    return path


def verify() -> dict:
    manifest = json.loads((BASE / "candidate_manifest.json").read_text(encoding="utf-8"))
    summary, candidates = manifest["summary"], manifest["candidates"]
    checks = {}
    for page in summary["sources"][0]["pages"]:
        path = resolve_repo_path(page["path"])
        checks[f"fireprot_page_{page['offset']}"] = path.is_file() and digest(path) == page["sha256"]
    for source in summary["sources"][1:]:
        path = resolve_repo_path(source["path"])
        checks[f"{source['name']}_sha"] = path.is_file() and digest(path) == source["sha256"]
    for reference in summary["reference_documents"]:
        if reference.get("path"):
            path = resolve_repo_path(reference["path"])
            checks[f"reference_{reference['name']}"] = path.is_file() and digest(path) == reference["sha256"]
        for sequence in reference.get("records", []):
            path = resolve_repo_path(sequence["path"])
            checks[f"fireprotdb_sequence_{sequence['path']}"] = path.is_file() and digest(path) == sequence["sha256"]

    train = {}
    for name in ("S2648", "VB1432"):
        path = BASE / "raw" / f"ptmul_train_{name.lower()}.csv"
        with path.open(newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        train[name] = rows

    identities = set()
    mapped = 0
    for i, candidate in enumerate(candidates):
        if candidate.get("source") == "ProtDDG-Bench PTMUL (ProTherm-derived)":
            for name in train:
                cluster = candidate["cluster"]
                expected_cluster = any(row["CLUID"] == cluster for row in train[name])
                if candidate["train_cluster_overlap"][name] != expected_cluster:
                    raise AssertionError(f"PTMUL {name} cluster-overlap drift: {candidate['source_variant_id']}")
                expected_components = [m for m in candidate["mutations_author"] if any(row["CLUID"] == cluster and row["MUT"] == m for row in train[name])]
                if candidate["train_single_component_overlap"][name] != expected_components:
                    raise AssertionError(f"PTMUL {name} mutation-overlap drift: {candidate['source_variant_id']}")
        if not candidate.get("structural_mapping_verified"):
            continue
        source_path = resolve_repo_path(candidate["normalized_pdb"])
        raw_path = BASE / "raw" / "pdb" / f"{candidate['pdb']}.pdb"
        if digest(source_path) != candidate["normalized_pdb_sha256"] or digest(raw_path) != candidate["source_pdb_sha256"]:
            raise AssertionError(f"Structure hash mismatch: {candidate['source_variant_id']}")
        checks[f"normalized_structure_{i}"] = True
        checks[f"source_pdb_{i}"] = True
        structure = PDBParser(QUIET=True).get_structure("bench", str(source_path))[0]
        if list(structure.child_dict) != ["A"]:
            raise AssertionError(f"Expected only normalized chain A: {candidate['source_variant_id']}")
        sequence = "".join(seq1(r.resname) for r in structure["A"] if r.id[0] == " ")
        if sequence != candidate["sequence"]:
            raise AssertionError(f"Normalized sequence mismatch: {candidate['source_variant_id']}")
        for mutation in candidate["mutations_normalized"]:
            m = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", mutation)
            if not m or not 1 <= int(m.group(2)) <= len(sequence) or sequence[int(m.group(2)) - 1] != m.group(1):
                raise AssertionError(f"Normalized mutation mismatch: {candidate['source_variant_id']} {mutation}")
        if candidate["order"] == 1:
            if candidate["experimental_ddg_kcal_mol"] != candidate["experimental_ddg_source_value"]:
                raise AssertionError("FireProtDB sign changed unexpectedly")
            if not candidate.get("sequence_source_verified"):
                raise AssertionError(f"FireProtDB source sequence mismatch: {candidate['source_variant_id']}")
        else:
            if candidate["experimental_ddg_kcal_mol"] != -candidate["experimental_ddg_source_value"]:
                raise AssertionError("PTMUL sign conversion mismatch")
        identity = (candidate["source"], str(candidate["source_variant_id"]))
        if identity in identities:
            raise AssertionError(f"Duplicate source row: {identity}")
        identities.add(identity)
        mapped += 1

    counts = {order: sum(c["order"] == order for c in candidates) for order in (1, 2, 3)}
    if any(counts[o] > 20 for o in counts) or len(candidates) > 60 or mapped != len(candidates):
        raise AssertionError(f"Sample cap or mapping failed: {counts}; mapped={mapped}")
    if not all(checks.values()):
        raise AssertionError(f"Source file hashes failed: {[k for k,v in checks.items() if not v]}")
    result = {
        "status": "PASS_DATA_PREFLIGHT_ONLY",
        "prediction_status": summary["prediction_status"],
        "candidate_counts_by_order": counts,
        "mapped_candidates": mapped,
        "foldx_slots_used": summary["foldx_slot_budget"]["used"],
        "foldx_max_slots": summary["foldx_slot_budget"]["max_pair_slots"],
        "source_hash_checks": len(checks),
        "source_hash_failures": 0,
        "checks": checks,
    }
    (BASE / "verification.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, ensure_ascii=False))
