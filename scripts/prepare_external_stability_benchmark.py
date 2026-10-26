#!/usr/bin/env python
"""Fetch and identity-check a small, score-blind experimental ddG benchmark."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from Bio.PDB import PDBIO, PDBParser, Select
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "workspace_data/p0_p3_upgrade_20261016/benchmark"
PTMUL_URL = "https://raw.githubusercontent.com/protddg-bench/protddg-bench/master/PTMUL/ptmul-5fold.csv"
PTMUL_TRAIN_URLS = {
    "S2648": "https://raw.githubusercontent.com/protddg-bench/protddg-bench/master/PTMUL/train-s2648-test-ptmul.csv",
    "VB1432": "https://raw.githubusercontent.com/protddg-bench/protddg-bench/master/PTMUL/train-vb1432-test-ptmul.csv",
}
FIREPROT_API = "https://loschmidt.chemi.muni.cz/fireprotdb/api/search"
FIREPROT_QUERY = {
    "tree": {
        "operator": "AND",
        "children": [
            {"variable": "MUTATED_POSITION", "operator": "IS_EMPTY", "value": False},
            {"variable": "DDG", "operator": "IS_EMPTY", "value": False},
            {"variable": "DATASET", "operator": "EQUAL_TO", "value": "ProTherm"},
        ],
    }
}
AA = "ACDEFGHIKLMNPQRSTVWY"


def fetch(url: str, path: Path) -> bytes:
    if path.is_file() and path.stat().st_size:
        return path.read_bytes()
    request = Request(url, headers={"User-Agent": "ECD-stability-benchmark/1.0"})
    with urlopen(request, timeout=90) as response:
        data = response.read()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return data


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fireprot_pages() -> list[dict]:
    pages = []
    encoded = json.dumps(FIREPROT_QUERY, separators=(",", ":"))
    # Evenly spaced pages reduce the database's large-protein ordering bias.
    # Selection below is by protein identity and source order, never by score.
    for offset in range(0, 5001, 250):
        url = FIREPROT_API + "?" + urlencode({"format": "json", "limit": 50, "offset": offset, "query": encoded})
        path = OUT / "raw" / f"fireprotdb_protherm_offset_{offset:05d}.json"
        raw = fetch(url, path)
        payload = json.loads(raw)
        pages.append({"offset": offset, "url": url, "path": str(path.relative_to(ROOT)), "sha256": sha256(raw), "entries": payload})
    return pages


def fetch_reference_data() -> tuple[list[dict], dict[str, bytes]]:
    pages = fireprot_pages()
    refs = {}
    for name, url in {"PTMUL": PTMUL_URL, **{f"PTMUL_train_{k}": v for k, v in PTMUL_TRAIN_URLS.items()}}.items():
        path = OUT / "raw" / f"{name.lower()}.csv"
        refs[name] = fetch(url, path)
    for name, url in {
        "fireprotdb_help.html": "https://loschmidt.chemi.muni.cz/fireprotdb/help/",
        "fireprotdb_api_docs.html": "https://loschmidt.chemi.muni.cz/fireprotdb/api-docs/",
        "ptmul_source_article.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6606456/fullTextXML",
    }.items():
        refs[name] = fetch(url, OUT / "raw" / name)
    return pages, refs


def _finite_ddg(experiment: dict) -> float | None:
    vals = [m.get("numValue") for m in experiment.get("measurements", []) if m.get("type") == "DDG"]
    vals = [float(v) for v in vals if isinstance(v, (int, float))]
    return vals[0] if len(vals) == 1 else None


def fireprot_candidates(pages: list[dict]) -> tuple[list[dict], dict]:
    candidates, excluded = [], {"non_mutant": 0, "not_single_substitution": 0, "missing_or_ambiguous_ddg": 0, "no_protein_accession": 0, "no_coordinate_mapping_annotation": 0, "duplicate_protein": 0}
    seen_proteins = set()
    seen_variants = set()
    for page in pages:
        for entry in page["entries"]:
            mutant = entry.get("mutant")
            if not mutant:
                excluded["non_mutant"] += 1
                continue
            if len(mutant.get("substitutions") or []) != 1 or mutant.get("insertions") or mutant.get("deletions"):
                excluded["not_single_substitution"] += 1
                continue
            experiment = mutant.get("experiment") or {}
            ddg = _finite_ddg(experiment)
            if ddg is None:
                excluded["missing_or_ambiguous_ddg"] += 1
                continue
            links = (mutant.get("sourceSequence") or {}).get("proteinLinks") or []
            proteins = [link.get("protein") or {} for link in links]
            accession = next((ref.get("accession") for protein in proteins for ref in protein.get("references", []) if ref.get("type") == "UNIPROTKB"), None)
            protein_name = next((p.get("name") for p in proteins if p.get("name")), None)
            if not accession:
                excluded["no_protein_accession"] += 1
                continue
            substitution = mutant["substitutions"][0]
            annotation = next((a.get("strValue") for a in experiment.get("annotations", []) if a.get("type") == "_PDB_CHAIN_MUTATION" and a.get("strValue")), None)
            match = re.fullmatch(r"([0-9A-Za-z]{4})_([A-Za-z0-9]):([A-Z])([0-9]+)([A-Z])", annotation or "")
            if not match:
                excluded["no_coordinate_mapping_annotation"] += 1
                continue
            pdb, chain, wt, author_pos, aa_to = match.groups()
            if (wt, int(author_pos), aa_to) != (substitution["sourceAa"], int(substitution["position"]), substitution["targetAa"]):
                excluded["no_coordinate_mapping_annotation"] += 1
                continue
            variant_key = (accession, wt, int(author_pos), aa_to)
            if variant_key in seen_variants:
                continue
            if accession in seen_proteins:
                excluded["duplicate_protein"] += 1
                continue
            seen_variants.add(variant_key)
            seen_proteins.add(accession)
            annotations = {a.get("type"): (a.get("numValue") if a.get("numValue") is not None else a.get("strValue")) for a in experiment.get("annotations", [])}
            structure = next((s for s in mutant.get("structures", []) if s.get("wwpdb", "").upper() == pdb.upper() and any(r.get("chainName") == chain and int(r.get("structPosition", -99999)) == int(author_pos) and int(r.get("seqPosition", -99999)) == int(substitution["position"]) and r.get("aminoAcid") == wt for r in s.get("residues", []))), None)
            candidates.append({
                "order": 1, "source": "FireProtDB 2.0 / ProTherm", "source_variant_id": mutant.get("id"), "experiment_id": experiment.get("id"),
                "uniprot": accession, "protein": protein_name, "pdb": pdb.upper(), "chain": chain,
                "source_sequence_id": (mutant.get("sourceSequence") or {}).get("id"), "sequence_position": int(substitution["position"]),
                "mutations_author": [f"{wt}{author_pos}{aa_to}"], "experimental_ddg_source_value": ddg, "experimental_ddg_kcal_mol": ddg,
                "source_sign_convention": "FireProtDB: negative stabilizing, positive destabilizing",
                "comparison_sign_convention": "negative stabilizing, positive destabilizing",
                "publication": experiment.get("publication"),
                "ph": annotations.get("PH"), "temperature_c": annotations.get("EXP_TEMPERATURE"),
                "method": annotations.get("METHOD"), "measure": annotations.get("MEASURE"),
                "mapping_metadata_match": structure is not None, "map_annotation": annotation,
                "raw_page_offset": page["offset"], "selection_rule": "first eligible unique UniProt in evenly spaced source pages; no prediction-score selection",
            })
            if len(candidates) == 20:
                return candidates, excluded
    return candidates, excluded


def _parse_mutations(text: str) -> list[tuple[str, int, str]]:
    rows = []
    for item in text.split(";"):
        m = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", item.strip())
        if not m:
            raise ValueError(f"Unrecognized PTMUL mutation: {item}")
        rows.append((m.group(1), int(m.group(2)), m.group(3)))
    return rows


def ptmul_candidates(ptmul_bytes: bytes, train_bytes: dict[str, bytes]) -> tuple[list[dict], dict]:
    rows = list(csv.DictReader(ptmul_bytes.decode("utf-8-sig").splitlines()))
    train_clusters, train_mutations = {}, {}
    for name, data in train_bytes.items():
        train_rows = list(csv.DictReader(data.decode("utf-8-sig").splitlines()))
        train_clusters[name] = {row["CLUID"] for row in train_rows if row.get("CLUID")}
        train_mutations[name] = {}
        for train_row in train_rows:
            if train_row.get("CLUID") and train_row.get("MUT"):
                train_mutations[name].setdefault(train_row["CLUID"], set()).add(train_row["MUT"])
    selected, seen_cluster_order = [], {2: set(), 3: set()}
    skipped = {"order_not_2_or_3": 0, "duplicate_cluster": 0, "duplicate_variant": 0, "invalid_label": 0}
    seen_variants = set()
    for row in rows:
        try:
            order, value = int(row["NMUT"]), float(row["DDG"])
        except (ValueError, KeyError):
            skipped["invalid_label"] += 1
            continue
        if order not in (2, 3):
            skipped["order_not_2_or_3"] += 1
            continue
        if sum(c["order"] == order for c in selected) >= 20:
            continue
        cluster = row["CLUID"]
        if cluster in seen_cluster_order[order]:
            skipped["duplicate_cluster"] += 1
            continue
        variant_key = (cluster, row["MUTS"])
        if variant_key in seen_variants:
            skipped["duplicate_variant"] += 1
            continue
        try:
            mutations = _parse_mutations(row["MUTS"])
        except ValueError:
            skipped["invalid_label"] += 1
            continue
        if len(mutations) != order or not all(w in AA and m in AA for w, _, m in mutations):
            skipped["invalid_label"] += 1
            continue
        selected.append({
            "order": order, "source": "ProtDDG-Bench PTMUL (ProTherm-derived)", "source_variant_id": f"{row['CLUID']}:{row['MUTS']}",
            "pdb_chain_id": row["PDB"], "pdb": row["PDB"][:4].upper(), "chain": row["PDB"][4:], "cluster": cluster,
            "mutations_author": [f"{w}{pos}{m}" for w, pos, m in mutations], "experimental_ddg_source_value": value,
            "experimental_ddg_kcal_mol": -value,
            "source_sign_convention": "DDGun/PTMUL source paper classifies stabilizing as ΔΔG >= 0 and destabilizing as ΔΔG < 0",
            "comparison_sign_convention": "negative stabilizing, positive destabilizing; PTMUL value sign inverted for comparison",
            "sign_conversion": "comparison_ddg = -source_ddg",
            "ph": None, "temperature_c": None, "method": None, "measure": None,
            "conditions_status": "not supplied in PTMUL table; source paper / ProTherm experiment lookup required",
            "train_cluster_overlap": {name: cluster in clusters for name, clusters in train_clusters.items()},
            "train_single_component_overlap": {name: [f"{w}{pos}{m}" for w, pos, m in mutations if f"{w}{pos}{m}" in train_mutations[name].get(cluster, set())] for name in train_mutations},
            "selection_rule": "first row per mutation order and CLUID in source order; no prediction-score selection",
        })
        seen_cluster_order[order].add(cluster)
        seen_variants.add(variant_key)
        if sum(c["order"] == order for c in selected) >= 20 and all(sum(c["order"] == o for c in selected) >= 20 for o in (2, 3)):
            break
    for order in (2, 3):
        count = sum(c["order"] == order for c in selected)
        if count < 20:
            skipped[f"order_{order}_available_unique_clusters"] = count
    return selected, skipped


class _ChainSelect(Select):
    def __init__(self, chain): self.chain = chain
    def accept_chain(self, chain): return int(chain.id == self.chain)
    def accept_model(self, model): return int(model.id == 0)
    def accept_residue(self, residue): return int(residue.id[0] == " " and residue.resname.upper() in {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"})
    def accept_atom(self, atom): return int(atom.element not in {"H", "D"})


def normalize_candidate(candidate: dict, source_dir: Path) -> tuple[dict, str | None]:
    pdb = candidate["pdb"]
    chain = candidate["chain"]
    source = source_dir / f"{pdb}.pdb"
    if not source.exists():
        fetch(f"https://files.rcsb.org/download/{pdb}.pdb", source)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb, str(source))
    if not structure or not len(structure[0]) or chain not in structure[0]:
        return candidate, "structure_or_chain_missing"
    residues = [r for r in structure[0][chain] if r.id[0] == " " and r.resname.upper() in {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}]
    residue_map = {}
    sequence = []
    for index, residue in enumerate(residues, 1):
        aa = seq1(residue.resname, custom_map={"MSE": "M"})
        sequence.append(aa)
        key = (int(residue.id[1]), str(residue.id[2]).strip())
        if key in residue_map:
            return candidate, "ambiguous_author_residue"
        residue_map[key] = (index, aa)
    normalized_mutations = []
    for notation in candidate["mutations_author"]:
        match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", notation)
        wt, pos, mutant = match.groups()
        mapped = residue_map.get((int(pos), ""))
        if not mapped or mapped[1] != wt:
            return candidate, f"author_mapping_mismatch:{notation}"
        normalized_mutations.append(f"{wt}{mapped[0]}{mutant}")
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(candidate["source_variant_id"]))
    out = OUT / "normalized_structures_final" / f"{safe_id}.pdb"
    out.parent.mkdir(parents=True, exist_ok=True)
    from copy import deepcopy
    normalized = deepcopy(structure)
    model = normalized[0]
    for chain_obj in list(model):
        if chain_obj.id != chain:
            model.detach_child(chain_obj.id)
    chain_obj = model[chain]
    chain_obj.id = "A"
    standard_residues = [r for r in chain_obj if r.id[0] == " " and r.resname.upper() in {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}]
    for residue in list(chain_obj):
        chain_obj.detach_child(residue.id)
    for index, residue in enumerate(standard_residues, 1):
        residue.id = (" ", index, " ")
        chain_obj.add(residue)
    writer = PDBIO()
    writer.set_structure(normalized)
    writer.save(str(out), _ChainSelect("A"))
    roundtrip = parser.get_structure("normalized", str(out))[0]["A"]
    roundtrip_sequence = "".join(seq1(r.resname) for r in roundtrip if r.id[0] == " ")
    if roundtrip_sequence != "".join(sequence):
        return candidate, "normalized_sequence_roundtrip_mismatch"
    for notation in normalized_mutations:
        match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", notation)
        if not match or roundtrip_sequence[int(match.group(2)) - 1] != match.group(1):
            return candidate, f"normalized_wildtype_mismatch:{notation}"
    candidate["source_pdb_sha256"] = sha256(source.read_bytes())
    candidate["normalized_pdb"] = str(out.relative_to(ROOT))
    candidate["normalized_pdb_sha256"] = sha256(out.read_bytes())
    candidate["sequence"] = "".join(sequence)
    candidate["mutations_normalized"] = normalized_mutations
    candidate["structural_mapping_verified"] = True
    candidate["source_pdb_url"] = f"https://files.rcsb.org/download/{pdb}.pdb"
    return candidate, None


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    pages, refs = fetch_reference_data()
    fireprot, fireprot_excluded = fireprot_candidates(pages)
    train = {name: refs[f"PTMUL_train_{name}"] for name in PTMUL_TRAIN_URLS}
    ptmul, ptmul_excluded = ptmul_candidates(refs["PTMUL"], train)
    raw_sources = [
        {"name": "FireProtDB ProTherm API sampled pages", "url": FIREPROT_API, "query": FIREPROT_QUERY, "pages": [{k: p[k] for k in ("offset", "url", "path", "sha256")} for p in pages]},
        {"name": "PTMUL", "url": PTMUL_URL, "path": "workspace_data/p0_p3_upgrade_20261016/benchmark/raw/ptmul.csv", "sha256": sha256(refs["PTMUL"])},
        *[{"name": f"PTMUL train {name}", "url": PTMUL_TRAIN_URLS[name], "path": f"workspace_data/p0_p3_upgrade_20261016/benchmark/raw/ptmul_train_{name.lower()}.csv", "sha256": sha256(train[name])} for name in PTMUL_TRAIN_URLS],
    ]
    candidates, map_failures = [], []
    fireprot_sequences = {}
    for candidate in fireprot:
        sequence_id = candidate["source_sequence_id"]
        url = f"https://loschmidt.chemi.muni.cz/fireprotdb/api/sequences/{sequence_id}/sequence"
        path = OUT / "raw" / f"fireprotdb_sequence_{sequence_id}.json"
        raw = fetch(url, path)
        fireprot_sequences[sequence_id] = {"url": url, "path": str(path.relative_to(ROOT)), "sha256": sha256(raw), "sequence": json.loads(raw)}
        pos = candidate["sequence_position"]
        candidate["sequence_source_verified"] = len(fireprot_sequences[sequence_id]["sequence"]) >= pos and fireprot_sequences[sequence_id]["sequence"][pos - 1] == candidate["mutations_author"][0][0]
    for candidate in fireprot + ptmul:
        try:
            candidate, error = normalize_candidate(candidate, OUT / "raw" / "pdb")
        except Exception as exc:
            error = f"structure_fetch_or_parse:{type(exc).__name__}:{exc}"
        if error:
            candidate["structural_mapping_verified"] = False
            candidate["mapping_failure"] = error
            map_failures.append(candidate["source_variant_id"])
        candidates.append(candidate)
    per_order = {order: [c for c in candidates if c["order"] == order and c.get("structural_mapping_verified")] for order in (1, 2, 3)}
    pilot, cursor = [], {order: 0 for order in (1, 2, 3)}
    while len(pilot) < 10 and any(cursor[o] < len(per_order[o]) for o in per_order):
        for order in (1, 2, 3):
            if cursor[order] < len(per_order[order]) and len(pilot) < 10:
                candidate = per_order[order][cursor[order]]
                candidate["timing_pilot"] = True
                pilot.append(candidate)
                cursor[order] += 1
    pilot_ids = {str(c["source_variant_id"]) for c in pilot}
    for candidate in candidates:
        candidate.setdefault("timing_pilot", str(candidate.get("source_variant_id")) in pilot_ids)
    for candidate in fireprot:
        match_ids = {}
        for name, data in train.items():
            train_rows = list(csv.DictReader(data.decode("utf-8-sig").splitlines()))
            match_ids[name] = any(row.get("PDB", "").upper() == f"{candidate['pdb']}{candidate['chain']}".upper() and row.get("MUT") == candidate["mutations_author"][0] for row in train_rows)
        candidate["training_exact_variant_overlap"] = match_ids
    eligible = [c for c in candidates if c.get("structural_mapping_verified")]
    summary = {
        "schema": "external_stability_benchmark_preflight_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "prediction_status": "NOT_RUN_COORDINATED_WITH_PRODUCTION_MODEL_WORK",
        "foldx_slot_budget": {"max_variants": 60, "max_runs_per_variant": 5, "max_pair_slots": 300, "used": 0},
        "selection_caps": {"single": 20, "double": 20, "triple": 20},
        "sources": raw_sources,
        "reference_documents": [
            {"name": "FireProtDB 2.0 API documentation", "url": "https://loschmidt.chemi.muni.cz/fireprotdb/api-docs/", "path": "workspace_data/p0_p3_upgrade_20261016/benchmark/raw/fireprotdb_api_docs.html", "sha256": sha256(refs["fireprotdb_api_docs.html"])},
            {"name": "FireProtDB stability help and sign definition", "url": "https://loschmidt.chemi.muni.cz/fireprotdb/help/", "path": "workspace_data/p0_p3_upgrade_20261016/benchmark/raw/fireprotdb_help.html", "sha256": sha256(refs["fireprotdb_help.html"])},
            {"name": "DDGun/PTMUL source paper PMID 31266447", "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6606456/fullTextXML", "path": "workspace_data/p0_p3_upgrade_20261016/benchmark/raw/ptmul_source_article.xml", "sha256": sha256(refs["ptmul_source_article.xml"])},
            {"name": "FireProtDB per-sequence endpoints", "records": list(fireprot_sequences.values())},
        ],
        "fireprotdb_exclusions": fireprot_excluded,
        "ptmul_exclusions": ptmul_excluded,
        "candidate_counts": {"single_selected": len(fireprot), "double_selected": sum(c["order"] == 2 for c in ptmul), "triple_selected": sum(c["order"] == 3 for c in ptmul), "all_structurally_mapped": len(eligible), "mapping_failures": len(map_failures), "fireprot_sequence_mismatches": sum(not c.get("sequence_source_verified", False) for c in fireprot)},
        "scientific_limits": [
            "Selection is capped by order and score-blind; uneven API/source coverage and family dependence remain.",
            "FireProtDB entries report experimental DDG with conditions and author mapping where recorded; missing conditions or coordinate metadata are kept as gaps.",
            "PTMUL table lacks per-row pH, temperature, and assay details; supplementary/source-paper lookup is still needed.",
            "PTMUL cluster overlap with S2648/VB1432 training candidates is recorded; no claim of broad independent validation is allowed.",
            "SPURS training-set overlap is unknown until its training records are compared; model execution and predictions have not been run.",
            "No ΔTm or binding ΔΔG values enter the folding stability endpoint.",
        ],
    }
    (OUT / "candidate_manifest.json").write_text(json.dumps({"summary": summary, "candidates": candidates}, indent=2, ensure_ascii=False), encoding="utf-8")
    (OUT / "benchmark_acceptance.json").write_text(json.dumps({"status": "DATA_PREFLIGHT_ONLY", "summary": summary["candidate_counts"], "prediction_status": summary["prediction_status"], "foldx_slots_used": 0, "mapping_failures": map_failures}, indent=2, ensure_ascii=False), encoding="utf-8")
    with (OUT / "candidates.csv").open("w", newline="", encoding="utf-8-sig") as f:
        if candidates:
            keys = sorted(set().union(*(c.keys() for c in candidates)))
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            for row in candidates:
                writer.writerow({k: json.dumps(row[k], ensure_ascii=False) if isinstance(row.get(k), (list, dict)) else row.get(k) for k in keys})
    lines = [
        "# 外部实测折叠稳定性基准：数据预审",
        "",
        "状态：**数据预审完成；预测未运行**。本清单用于后续限额FoldX/SPURS验证，不构成生产阈值校准。",
        "",
        "## 数据与映射",
        "",
        f"FireProtDB ProTherm：在API结果偏移0–5000、每250行取一页的固定分页中，score-blind规则选出{len(fireprot)}条单点记录；每个UniProt最多1条。PTMUL：原始914条中取{sum(c['order']==2 for c in ptmul)}个双点和{sum(c['order']==3 for c in ptmul)}个三点，每阶每个官方CLUID最多1条。",
        f"结构映射与WT字母核对通过{len(eligible)}/{len(candidates)}条；失败{len(map_failures)}条。FireProtDB源序列核对失败{summary['candidate_counts']['fireprot_sequence_mismatches']}条。归一化结构只保留目标蛋白链、连续重编号，保留作者编号到新编号映射及源/输出SHA256。",
        "",
        "| 阶数 | 已选 | 实验条件缺失/不全 | S2648训练簇重叠 | VB1432训练簇重叠 |",
        "|---|---:|---:|---:|---:|",
    ]
    for order, label in ((1, "单点"), (2, "双点"), (3, "三点")):
        group = [c for c in candidates if c["order"] == order]
        unknown = sum(c.get("ph") is None or c.get("temperature_c") is None for c in group)
        overlap_s = sum(c.get("train_cluster_overlap", {}).get("S2648", False) for c in group)
        overlap_v = sum(c.get("train_cluster_overlap", {}).get("VB1432", False) for c in group)
        lines.append(f"| {label} | {len(group)} | {unknown} | {overlap_s} | {overlap_v} |")
    lines += [
        "",
        "PTMUL每条记录还保留了其变体中各单点是否出现在相应训练文件的诊断。当前抽样与S2648/VB1432训练集存在同簇重叠，因此这些样本不能整体称为蛋白簇独立测试集；SPURS训练集重叠也尚未确认。FireProtDB 8个单点虽来自不同UniProt ID，但未完成序列同源聚类，故独立性仍属有限。",
        "",
        "## 符号与边界",
        "",
        "FireProtDB帮助页定义ΔΔG负值稳定化、正值去稳定化。DDGun/PTMUL源论文表格把ΔΔG≥0归为稳定化、ΔΔG<0归为去稳定化，方向相反；候选清单保留源值，并把PTMUL值取反到比较字段。ΔTm、结合ΔΔG没有混入折叠ΔΔG。PTMUL表没有逐实验pH、温度及检测条件，必须列为缺口。",
        "",
        "## 预留计算",
        "",
        f"当前可映射样本{len(eligible)}条：按阶数最多{len(per_order[1])}/{len(per_order[2])}/{len(per_order[3])}（单/双/三），完整5次FoldX理论上限{len(eligible)*5}个WT/Mut配对槽，低于300上限。首10条计时试跑按4单、3双、3三轮转选择（50槽）；其余38条若执行，再用190槽。所有预测仍待与生产任务错峰协调，当前实际消耗0槽。",
        "SPURS当前预测尚未运行；不得据此计算MAE/RMSE/偏移/排序或评判准入线。待产生分数后应按蛋白簇、阶数报告误差和符号，并将未知训练重叠单列。",
        "",
        "## 来源",
        "",
        "- [FireProtDB 2.0 API](https://loschmidt.chemi.muni.cz/fireprotdb/api-docs/) 与 [符号/字段帮助](https://loschmidt.chemi.muni.cz/fireprotdb/help/)；仅筛实验DDG、ProTherm数据集的单替换。",
        "- [PTMUL数据页](https://protddg-bench.github.io/ptmul/) 与 [公开文件/训练划分](https://github.com/protddg-bench/protddg-bench/tree/master/PTMUL)。PTMUL来源论文：[DDGun (PMID 31266447)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6606456/)。",
        "- 坐标结构从 [RCSB PDB](https://www.rcsb.org/) 原始PDB下载，来源文件和归一化文件hash见`candidate_manifest.json`。",
    ]
    (OUT / "benchmark_data_review_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary["candidate_counts"], ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
