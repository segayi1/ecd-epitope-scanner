"""Freeze canonical residue maps for the 2026 cross-protein structure cohort."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.request
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

COHORT = {
    "1AH1": {"1": ("CTLA4", "P16410")},
    "1UBQ": {"1": ("ubiquitin_mature", None)},
    "3RRQ": {"1": ("PD1", "Q15116")},
    "4Z18": {"1": ("PDL1", "Q9NZQ7")},
    "3PV7": {"1": ("B7H6", "Q68D85")},
    "1LYZ": {"1": ("hen_lysozyme", "P00698")},
    "1OMP": {"1": ("MBP", "P0AEX9")},
    "2LZM": {"1": ("T4_lysozyme", "P00720")},
    "1D3Z": {"1": ("ubiquitin_mature", None)},
    "4ZQK": {"1": ("PDL1", "Q9NZQ7"), "2": ("PD1", "Q15116")},
    "3PV6": {"1": ("B7H6", "Q68D85"), "2": ("NKp30_stress", "O14931")},
    "1ANF": {"1": ("MBP", "P0AEX9")},
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _one_letter(comp: str) -> str:
    return protein_letters_3to1_extended.get(comp.upper(), "X")


def _polymer_sequences(cif: dict) -> dict[str, str]:
    entity_ids = cif.get("_entity_poly.entity_id", [])
    sequences = cif.get("_entity_poly.pdbx_seq_one_letter_code_can", [])
    if len(entity_ids) != len(sequences):
        raise ValueError("entity polymer sequence columns do not agree")
    return {str(entity): re.sub(r"\s+", "", str(seq)).upper()
            for entity, seq in zip(entity_ids, sequences)}


def align_polymer(reference: str, polymer: str) -> tuple[dict[int, int], dict]:
    """Map zero-based polymer positions to one-based reference positions."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -8.0
    aligner.extend_gap_score = -0.5
    aligner.left_insertion_score = 0.0
    aligner.right_insertion_score = 0.0
    alignment = aligner.align(reference, polymer)[0]
    mapped: dict[int, int] = {}
    for (r0, r1), (p0, p1) in zip(*alignment.aligned):
        for ref_index, polymer_index in zip(range(r0, r1), range(p0, p1)):
            if reference[ref_index] != "X" and polymer[polymer_index] != "X":
                mapped[polymer_index] = ref_index + 1
    paired = [(polymer[i], reference[pos - 1]) for i, pos in mapped.items()]
    matches = sum(left == right for left, right in paired)
    return mapped, {"score": float(alignment.score), "mapped_positions": len(mapped),
                    "sequence_identity": matches / len(paired) if paired else None,
                    "reference_positions": sorted(mapped.values())}


def _fetch_uniprot(accession: str, source_dir: Path) -> tuple[str, dict]:
    path = source_dir / f"{accession}.fasta"
    if not path.is_file() or path.stat().st_size == 0:
        temp = path.with_suffix(".fasta.part")
        for attempt in range(3):
            request = urllib.request.Request(
                f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
                headers={"User-Agent": "ecd-epitope-scanner-cross-protein-map/1.0",
                         "Accept": "text/x-fasta"},
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                content = response.read()
            if content.startswith(b">"):
                temp.write_bytes(content); temp.replace(path); break
            if attempt == 2:
                raise ValueError(f"UniProt returned no FASTA content: {accession}")
            time.sleep(attempt + 1)
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"invalid UniProt FASTA: {accession}")
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
        raise ValueError(f"invalid UniProt sequence: {accession}")
    return sequence, {"source": f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
                      "fasta_sha256": sha256(path), "length": len(sequence)}


def build(input_root: Path, output_root: Path) -> dict:
    input_root = input_root.resolve(); output_root = output_root.resolve()
    input_manifest = json.loads((input_root / "input_manifest.json").read_text(encoding="utf-8"))
    for name, expected in input_manifest["files"].items():
        if sha256(input_root / name) != expected:
            raise ValueError(f"frozen input hash mismatch: {name}")
    source_dir = output_root / "sequence_sources"
    source_dir.mkdir(parents=True, exist_ok=True)
    references: dict[str, tuple[str, dict]] = {}
    polyubiquitin, ubiquitin_source = _fetch_uniprot("P0CG48", source_dir)
    ubiquitin = polyubiquitin[:76]
    ubiquitin_cif = MMCIF2Dict(str(input_root / "1UBQ.cif"))
    if _polymer_sequences(ubiquitin_cif)["1"] != ubiquitin:
        raise ValueError("UniProt P0CG48 residues 1-76 do not match the frozen mature ubiquitin structure")
    references["ubiquitin_mature"] = (ubiquitin, {**ubiquitin_source,
        "reference_region_1based": [1, 76],
        "reference_sequence_sha256": hashlib.sha256(ubiquitin.encode()).hexdigest(),
        "interpretation": "mature ubiquitin domain, mapped to P0CG48 residues 1-76"})
    accession_for_target = {target: accession for pdb_map in COHORT.values()
                            for target, accession in pdb_map.values() if accession}
    for target, accession in sorted(accession_for_target.items()):
        references[target] = _fetch_uniprot(accession, source_dir)

    observed: dict[str, dict[tuple[str, str, str, str], set[str]]] = {}
    for pdb in COHORT:
        with (input_root / f"{pdb}_residues.csv").open(newline="", encoding="utf-8") as stream:
            rows = csv.DictReader(stream)
            observed[pdb] = {}
            for row in rows:
                key = (row["author_chain"], row["author_residue"],
                       "" if row["insertion_code"] in {".", "?"} else row["insertion_code"],
                       row["residue_name"].upper())
                observed[pdb].setdefault(key, set()).add(row["model"])

    results = {}
    for pdb, target_entities in COHORT.items():
        cif = MMCIF2Dict(str(input_root / f"{pdb}.cif"))
        polymer = _polymer_sequences(cif)
        scheme_keys = ("asym_id", "entity_id", "seq_id", "mon_id", "auth_seq_num",
                       "pdb_strand_id", "pdb_ins_code")
        columns = [cif[f"_pdbx_poly_seq_scheme.{key}"] for key in scheme_keys]
        if len({len(column) for column in columns}) != 1:
            raise ValueError(f"{pdb}: inconsistent polymer residue mapping columns")
        chain_entities = {}
        for entity, chains in zip(cif.get("_entity_poly.entity_id", []),
                                  cif.get("_entity_poly.pdbx_strand_id", [])):
            for chain in str(chains).split(","):
                chain_entities[chain.strip()] = str(entity)
        entity_to_ref = {}
        for ref_id, db_name, accession, entity in zip(
                cif.get("_struct_ref.id", []), cif.get("_struct_ref.db_name", []),
                cif.get("_struct_ref.pdbx_db_accession", []), cif.get("_struct_ref.entity_id", [])):
            entity_to_ref.setdefault(str(entity), []).append({"id": str(ref_id), "db": str(db_name),
                                                               "accession": str(accession)})
        annotation_conflicts = {}
        conflict_keys = ("pdbx_pdb_strand_id", "pdbx_auth_seq_num", "pdbx_pdb_ins_code",
                         "mon_id", "pdbx_seq_db_accession_code", "pdbx_seq_db_seq_num", "db_mon_id", "details")
        conflict_columns = [cif.get(f"_struct_ref_seq_dif.{key}", []) for key in conflict_keys]
        if conflict_columns and len({len(column) for column in conflict_columns}) != 1:
            raise ValueError(f"{pdb}: inconsistent reference-conflict annotation columns")
        for chain, auth_seq, ins, model_aa, accession, ref_pos, ref_aa, detail in zip(*conflict_columns):
            key = (str(chain), str(auth_seq), "" if ins in {".", "?"} else str(ins))
            annotation_conflicts[key] = {"deposited_model_aa": _one_letter(str(model_aa)),
                "deposition_reference_accession": str(accession),
                "deposition_reference_position": str(ref_pos),
                "deposition_reference_aa": _one_letter(str(ref_aa)), "detail": str(detail)}
        entity_maps = {}; entity_rows = {}
        for entity, (target, accession) in target_entities.items():
            sequence, source = references[target]
            entity_sequence = polymer[entity]
            mapped, alignment = align_polymer(sequence, entity_sequence)
            entity_maps[entity] = (target, sequence, entity_sequence, mapped, alignment, source)
            entity_rows[entity] = []
        for asym, entity, seq_id, mon_id, auth_seq, auth_chain, ins_code in zip(*columns):
            entity = str(entity)
            if entity not in entity_maps: continue
            target, reference, chain_sequence, mapping, alignment, _source = entity_maps[entity]
            pindex = int(seq_id) - 1
            ins_code = "" if ins_code in {".", "?"} else str(ins_code)
            p_aa = _one_letter(str(mon_id))
            position = mapping.get(pindex)
            observed_models = observed[pdb].get(
                (str(auth_chain), str(auth_seq), ins_code, str(mon_id).upper()), set())
            conflict = annotation_conflicts.get((str(auth_chain), str(auth_seq), ins_code))
            if conflict and not str(conflict["deposition_reference_position"]).isdigit():
                position = None
            ref_aa = reference[position - 1] if position else None
            sequence_mismatch = bool(position and p_aa != ref_aa)
            entity_rows[entity].append({"label_chain": str(asym), "author_chain": str(auth_chain),
                "author_residue": str(auth_seq), "insertion_code": ins_code,
                "polymer_position": pindex + 1, "reference_position": position,
                "pdb_aa": p_aa, "reference_aa": ref_aa,
                "resolved_in_any_model": bool(observed_models), "resolved_model_count": len(observed_models),
                "engineered_or_sequence_mismatch": sequence_mismatch,
                "deposition_reference_conflict": conflict,
                "candidate_freeze_exclusion": sequence_mismatch or conflict is not None,
                "unmapped_polymer_residue": position is None})
        entity_reports = []
        for entity, (target, sequence, chain_sequence, mapped, alignment, source) in entity_maps.items():
            rows = entity_rows[entity]
            mismatches = [{"author_chain": r["author_chain"], "author_residue": r["author_residue"],
                           "insertion_code": r["insertion_code"], "reference_position": r["reference_position"],
                           "pdb_aa": r["pdb_aa"], "reference_aa": r["reference_aa"],
                           "resolved_in_any_model": r["resolved_in_any_model"],
                           "deposition_reference_conflict": r["deposition_reference_conflict"]}
                          for r in rows if r["engineered_or_sequence_mismatch"]]
            annotation_discrepancies = [r for r in rows if r["deposition_reference_conflict"]]
            unmapped = [r for r in rows if r["unmapped_polymer_residue"]]
            entity_reports.append({"entity_id": entity, "target": target,
                "deposition_struct_ref": entity_to_ref.get(entity, []),
                "reference_sequence": source, "polymer_sequence_length": len(chain_sequence),
                "alignment": alignment, "chains": sorted({r["author_chain"] for r in rows}),
                "resolved_residue_count": sum(r["resolved_in_any_model"] for r in rows),
                "residue_count": len(rows), "sequence_mismatches": mismatches,
                "deposition_annotation_discrepancies": annotation_discrepancies,
                "candidate_freeze_excluded_positions": sorted({r["reference_position"] for r in rows
                    if r["candidate_freeze_exclusion"] and r["reference_position"] is not None}),
                "unmapped_polymer_residue_count": len(unmapped),
                "unmapped_polymer_residues": unmapped, "residues": rows})
        model_ids = sorted({model for model_set in observed[pdb].values() for model in model_set})
        results[pdb] = {"input_cif_sha256": sha256(input_root / f"{pdb}.cif"),
                        "coordinate_model_count": len(model_ids), "coordinate_model_ids": model_ids,
                        "entity_mappings": entity_reports}

    output_root.mkdir(parents=True, exist_ok=True)
    result_path = output_root / "canonical_residue_maps.json"
    result = {"schema": "cross_protein_canonical_maps_v1", "input_manifest_sha256": sha256(input_root / "input_manifest.json"),
              "mapping_method": "global pairwise sequence alignment of frozen mmCIF entity sequence to the specified reference; author residue IDs and insertion codes preserved",
              "reference_accession_policy": "use current UniProt accessions; historical/inactive mmCIF accessions are retained as deposition references but are not treated as canonical coordinates",
              "stress_mapping_note": "1UBQ and 1D3Z map to the mature ubiquitin domain P0CG48 residues 1-76; frozen legacy mmCIF records P62988, now inactive, while RCSB currently links both entities to P0CG48. MBP uses current P0AEX9 where frozen structures carry inactive P02928.",
              "references": {target: source for target, (_sequence, source) in references.items()},
              "structures": results}
    result["mapping_sha256"] = hashlib.sha256(json.dumps(results, sort_keys=True,
        separators=(",", ":")).encode()).hexdigest()
    result_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    manifest = {"schema": "cross_protein_canonical_map_manifest_v1", "input_manifest_sha256": result["input_manifest_sha256"],
                "canonical_residue_maps_sha256": sha256(result_path),
                "mapping_script_sha256": sha256(Path(__file__).resolve()),
                "sequence_sources": {p.name: sha256(p) for p in sorted(source_dir.glob("*.fasta"))},
                "mapping_status": "REVIEW_REQUIRED_FOR_CANDIDATE_FREEZE"}
    (output_root / "canonical_mapping_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path(r"E:\user\New project\cross_protein_validation_20261006\inputs"))
    parser.add_argument("--output", type=Path, default=Path(r"E:\user\New project\cross_protein_validation_20261006\mapping"))
    args = parser.parse_args()
    result = build(args.inputs, args.output)
    print(json.dumps({"schema": result["schema"], "structures": len(result["structures"]),
                      "mapping_sha256": result["mapping_sha256"]}))


if __name__ == "__main__":
    main()
