"""Paired, geometry-only comparison of archived AF3 model_0 and AF3 selection.

Run from the repository parent as documented in docs/project_memory.md.
"""
from __future__ import annotations

import csv
import json
import math
import platform
import sys
import uuid
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))

from ecd_epitope_scanner import config  # noqa: E402
from ecd_epitope_scanner.af3_preprocessing import (  # noqa: E402
    import_af3_directory,
    select_af3_model,
)
from ecd_epitope_scanner.models import ProteinContext  # noqa: E402
from ecd_epitope_scanner.structure_analysis import (  # noqa: E402
    analyze_residues,
    load_structure,
    map_chain_to_uniprot,
)
from ecd_epitope_scanner.structure_reliability import classify_reliability  # noqa: E402
from ecd_epitope_scanner.surface_patch import run_surface_patch_v2  # noqa: E402
from ecd_epitope_scanner.util import sha256_file, write_json  # noqa: E402

MANIFEST = ROOT / "workspace_data" / "af3_binder_batch_20261017" / "batch_manifest.json"
PROTOCOL_FILES = ["af3_preprocessing.py", "structure_analysis.py", "structure_reliability.py",
                  "config.py", "models.py", "util.py", "compare_af3_preprocessing.py"]


def jaccard(left: set, right: set) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def adjusted_rand_index(left: list[list[int]], right: list[list[int]]) -> float:
    """ARI over the shared universe; missing elements are ignored by caller."""
    left_map = {item: i for i, group in enumerate(left) for item in group}
    right_map = {item: i for i, group in enumerate(right) for item in group}
    universe = sorted(left_map.keys() & right_map.keys())
    n = len(universe)
    if n < 2:
        return 1.0
    table = Counter((left_map[item], right_map[item]) for item in universe)
    rows, cols = Counter(), Counter()
    for (i, j), count in table.items():
        rows[i] += count
        cols[j] += count
    def choose2(values):
        return values * (values - 1) / 2
    sum_cells = sum(choose2(v) for v in table.values())
    sum_rows = sum(choose2(v) for v in rows.values())
    sum_cols = sum(choose2(v) for v in cols.values())
    total = choose2(n)
    expected = sum_rows * sum_cols / total
    maximum = (sum_rows + sum_cols) / 2
    return 1.0 if math.isclose(maximum, expected) else (sum_cells - expected) / (maximum - expected)


def _is_partition(groups: list[list[int]]) -> bool:
    members = [item for group in groups for item in group]
    return len(members) == len(set(members))


def kabsch_rmsd(reference, mobile) -> float | None:
    """RMSD after optimal rigid alignment; arrays contain corresponding XYZ rows."""
    a, b = np.asarray(reference, dtype=float), np.asarray(mobile, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3 or len(a) < 3:
        return None
    a0, b0 = a - a.mean(axis=0), b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(b0.T @ a0)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    delta = b0 @ rotation - a0
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _summary(source: Path, model_index: int) -> dict:
    stem = source.glob("*_summary_confidences_*.json")
    for path in stem:
        if path.name.endswith(f"_{model_index}.json"):
            return _read_json(path)
    return {}


def _metadata(source: Path) -> tuple[dict, dict]:
    request = next(iter(source.glob("*_job_request.json")), None)
    if request is None:
        return {}, {}
    data = _read_json(request)
    if isinstance(data, list):
        data = data[0] if data else {}
    if not isinstance(data, dict):
        return {}, {}
    return data, {"file": request.name, "sha256": sha256_file(request)}


def _source_inventory(source: Path) -> list[dict]:
    return [
        {"relative_path": str(path.relative_to(source)), "size_bytes": path.stat().st_size,
         "sha256": sha256_file(path)}
        for path in sorted(source.rglob("*")) if path.is_file()
    ]


def _protocol_fingerprints() -> dict:
    geometry_sources = sorted(path.relative_to(ROOT).as_posix() for path in (ROOT / "surface_patch").rglob("*.py"))
    return {name: sha256_file(ROOT / name) for name in PROTOCOL_FILES + geometry_sources}


def _confidence_for(source: Path, index: int) -> Path:
    paths = [p for p in source.glob(f"*_full_data_{index}.json") if p.is_file()]
    if len(paths) != 1:
        raise RuntimeError(f"Expected one confidence file for model {index} in {source}, got {len(paths)}")
    return paths[0]


def _metric_structure(structure_path: Path, confidence_path: Path, protein: ProteinContext,
                      chain: str, workdir: Path, *, build_surface: bool = True) -> dict:
    parsed = load_structure(structure_path)
    mapping = map_chain_to_uniprot(parsed, chain, protein)
    confidence = _read_json(confidence_path)
    workdir.mkdir(parents=True, exist_ok=True)
    records, methods = analyze_residues(
        parsed, mapping, protein, workdir, predicted=True, confidence=confidence,
    )
    reliability = classify_reliability(
        records, protein, predicted=True, confidence=confidence,
        chain_length=len(mapping.residues), single_chain=len(list(parsed.model.get_chains())) == 1,
    )
    # Use exactly the current pipeline's surface and reliability parameters.
    surface = run_surface_patch_v2(
        parsed.model, chain, mapping.uniprot_to_residue, records, protein,
        output_dir=None, confidence=confidence, predicted=True,
        parameters=config.SURFACE_V2_PARAMETERS, reliability=reliability,
    ) if build_surface else None
    by_position = {r.uniprot_pos: r for r in records}
    ca = {r.uniprot_pos: np.asarray(r.ca, dtype=float) for r in records if r.ca is not None}
    pl = {r.uniprot_pos: r.plddt for r in records}
    sasa = {r.uniprot_pos: float(r.sasa) for r in records}
    rsasa = {r.uniprot_pos: float(r.rsasa) for r in records}
    contacts = set()
    positions = sorted(ca)
    for i, a in enumerate(positions):
        for b in positions[i + 1:]:
            if b - a > 2 and np.linalg.norm(ca[a] - ca[b]) <= 8.0:
                contacts.add((a, b))
    regions = []
    for region in surface.legacy_regions if surface else []:
        members = sorted(set(region.envelope_residues or region.residues))
        regions.append({"region_id": region.id, "members": members,
                        "surface_area": float(region.surface_area),
                        "diameter_a": float(region.geodesic_diameter)})
    excluded = set(reliability["qc"].get("excluded_positions", []))
    accessible = {p for p, record in by_position.items() if record.surface_accessible}
    eligible = {p for p in accessible if p not in excluded and p in mapping.uniprot_to_residue}
    return {
        "parsed": parsed, "mapping": mapping, "records": records,
        "ca": ca, "plddt": pl, "sasa": sasa, "rsasa": rsasa, "contacts": contacts,
        "accessible": accessible, "eligible": eligible, "excluded": excluded,
        "regions": regions, "surface_summary": surface.summary if surface else None,
        "mean_plddt": float(np.mean([value for value in pl.values() if value is not None])) if any(value is not None for value in pl.values()) else None,
        "analysis_methods": methods, "reliability_qc": reliability["qc"],
    }


def _compare_geometry(old: dict, new: dict) -> dict:
    common = sorted(old["ca"].keys() & new["ca"].keys())
    high = [p for p in common if old["plddt"].get(p) is not None and new["plddt"].get(p) is not None
            and old["plddt"][p] >= 70 and new["plddt"][p] >= 70]
    global_rmsd = kabsch_rmsd([old["ca"][p] for p in common], [new["ca"][p] for p in common])
    core_rmsd = kabsch_rmsd([old["ca"][p] for p in high], [new["ca"][p] for p in high])
    # Align on the shared high-confidence core before reporting residue motion.
    aa = np.asarray([old["ca"][p] for p in high])
    bb = np.asarray([new["ca"][p] for p in high])
    if len(high) >= 3:
        a0, b0 = aa - aa.mean(0), bb - bb.mean(0)
        u, _, vt = np.linalg.svd(b0.T @ a0)
        correction = np.eye(3)
        correction[-1, -1] = np.linalg.det(u @ vt)
        rotation = u @ correction @ vt
        aligned = (np.asarray([new["ca"][p] for p in common]) - bb.mean(0)) @ rotation + aa.mean(0)
        displacements = {p: float(np.linalg.norm(aligned[i] - old["ca"][p])) for i, p in enumerate(common)}
    else:
        displacements = {p: None for p in common}
    old_regions = [x["members"] for x in old["regions"]]
    new_regions = [x["members"] for x in new["regions"]]
    exposure_positions = sorted(old["rsasa"].keys() & new["rsasa"].keys())
    old_membership = {p for group in old_regions for p in group}
    new_membership = {p for group in new_regions for p in group}
    return {
        "matched_ca_count": len(common), "high_confidence_core_count": len(high),
        "global_ca_rmsd_a": global_rmsd, "matched_high_confidence_core_rmsd_a": core_rmsd,
        "core_aligned_residue_displacement_a": {str(k): v for k, v in displacements.items()},
        "per_residue_plddt_delta": {str(p): (new["plddt"][p] - old["plddt"][p])
                                    for p in sorted(old["plddt"].keys() & new["plddt"].keys())
                                    if old["plddt"][p] is not None and new["plddt"][p] is not None},
        "per_residue_sasa_delta_a2": {str(p): new["sasa"][p] - old["sasa"][p] for p in exposure_positions},
        "per_residue_rsasa_delta": {str(p): new["rsasa"][p] - old["rsasa"][p] for p in exposure_positions},
        "mean_rsasa_old": float(np.mean([old["rsasa"][p] for p in exposure_positions])) if exposure_positions else None,
        "mean_rsasa_new": float(np.mean([new["rsasa"][p] for p in exposure_positions])) if exposure_positions else None,
        "mean_sasa_delta_a2": float(np.mean([new["sasa"][p] - old["sasa"][p] for p in exposure_positions])) if exposure_positions else None,
        "displacement_median_a": float(np.median([v for v in displacements.values() if v is not None])) if any(v is not None for v in displacements.values()) else None,
        "displacement_p90_a": float(np.percentile([v for v in displacements.values() if v is not None], 90)) if any(v is not None for v in displacements.values()) else None,
        "old_contact_count": len(old["contacts"]), "new_contact_count": len(new["contacts"]),
        "contact_jaccard": jaccard(old["contacts"], new["contacts"]),
        "contacts_gained": len(new["contacts"] - old["contacts"]),
        "contacts_lost": len(old["contacts"] - new["contacts"]),
        "old_accessible_count": len(old["accessible"]), "new_accessible_count": len(new["accessible"]),
        "accessible_gain": sorted(new["accessible"] - old["accessible"]),
        "accessible_loss": sorted(old["accessible"] - new["accessible"]),
        "old_eligible_count": len(old["eligible"]), "new_eligible_count": len(new["eligible"]),
        "eligible_gain": sorted(new["eligible"] - old["eligible"]),
        "eligible_loss": sorted(old["eligible"] - new["eligible"]),
        "old_excluded_count": len(old["excluded"]), "new_excluded_count": len(new["excluded"]),
        "excluded_gain": sorted(new["excluded"] - old["excluded"]),
        "excluded_loss": sorted(old["excluded"] - new["excluded"]),
        "old_region_count": len(old_regions), "new_region_count": len(new_regions),
        "region_membership_jaccard": jaccard(old_membership, new_membership),
        "region_membership_ari": adjusted_rand_index(old_regions, new_regions)
                                   if _is_partition(old_regions) and _is_partition(new_regions) else None,
        "region_best_match_jaccard_mean": float(np.mean(
            [max((jaccard(set(a), set(b)) for b in new_regions), default=0.0) for a in old_regions]
            + [max((jaccard(set(b), set(a)) for a in old_regions), default=0.0) for b in new_regions]
        )) if old_regions or new_regions else 1.0,
        "old_regions": old["regions"], "new_regions": new["regions"],
        "old_surface_summary": old["surface_summary"], "new_surface_summary": new["surface_summary"],
    }


def _read_table(path: Path, delimiter=",") -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def _write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run() -> Path:
    protocol_start = _protocol_fingerprints()
    batch = _read_json(MANIFEST)
    jobs = {item["binder"]: item for item in batch["jobs"]}
    preflight = {item["binder"]: item for item in batch["preflight_items"]}
    tag = datetime.now(UTC).strftime("%Y%m%d") + "_" + uuid.uuid4().hex[:8]
    out = ROOT / "test_output" / f"af3_preprocessing_comparison_{tag}"
    out.mkdir(parents=True)
    cases, compact = [], []
    for binder in sorted(jobs):
        check = preflight[binder]
        source = Path(check["source_directory"])
        if not source.is_dir():
            raise FileNotFoundError(source)
        protein = ProteinContext(
            accession=binder, protein_name=binder, organism="AF3 binder input",
            sequence=check["sequence"], ecd_start=1, ecd_end=len(check["sequence"]),
            membrane_anchor_position=len(check["sequence"]), topology="binder chain",
        )
        case_dir = out / binder
        case_dir.mkdir()
        before = _source_inventory(source)
        # Archive and select through the new public preprocessing implementation.
        af3_manifest = import_af3_directory(source, case_dir / "frozen")
        selected = select_af3_model(af3_manifest, protein, chain=check["chain"],
                                    output_dir=case_dir / "selection")
        new_struct, new_conf_path = selected.structure_path, selected.confidence_path
        selection_report, report_path = selected.report, selected.report_path
        selected_row = next(row for row in selection_report["models"]
                            if row["model_id"] == selection_report["selected_model_id"])
        frozen_rows = _read_json(af3_manifest)["models"]
        frozen_selected = next(row for row in frozen_rows if row["id"] == selection_report["selected_model_id"])
        # Anchor the old condition to the manifest's exact historical model_0 bytes.
        old_struct = source / check["selected_structure"]["name"]
        old_conf_path = _confidence_for(source, 0)
        old_summary = _summary(source, 0)
        model0_item = next(item for item in check["all_model_hashes"] if item["model"] == 0)
        if sha256_file(old_struct) != model0_item["structure_sha256"] or sha256_file(old_conf_path) != model0_item["confidence_sha256"]:
            raise RuntimeError(f"{binder}: historical model_0 source hashes differ from batch manifest")
        old_metrics = _metric_structure(old_struct, old_conf_path, protein, check["chain"],
                                        case_dir / "old_model0")
        new_metrics = _metric_structure(new_struct, new_conf_path, protein, check["chain"],
                                        case_dir / "new_selected")
        geometry = _compare_geometry(old_metrics, new_metrics)
        old_source_rows = [row for row in frozen_rows
                           if row["model_index"] == 0 and row["sha256"]["structure"] == model0_item["structure_sha256"]]
        if not old_source_rows:
            raise RuntimeError(f"{binder}: original model_0 was not represented in the frozen source manifest")
        old_source_row = old_source_rows[0]
        old_report_model = next((row for row in selection_report["models"]
                                 if row["model_id"] == old_source_row["id"]), None)
        old_selection = {"model_id": old_source_row["id"], "model_index": 0,
                         "job_id": old_source_row["job_id"], "seed_id": old_source_row["seed_id"],
                         "sample_id": old_source_row["sample_id"],
                         "source_name": old_source_row["source_name"],
                         "source_sha256": sha256_file(old_struct), "confidence_sha256": sha256_file(old_conf_path),
                         "ptm": old_summary.get("ptm"), "ranking_score": old_summary.get("ranking_score"),
                         "mean_plddt": old_metrics["mean_plddt"],
                         "consensus_score": old_report_model.get("consensus_score") if old_report_model else None}
        new_selection = {"model_id": selected_row["model_id"], "model_index": frozen_selected["model_index"],
                         "job_id": frozen_selected["job_id"], "seed_id": frozen_selected["seed_id"],
                         "sample_id": frozen_selected["sample_id"],
                         "source_sha256": selection_report["selected_sha256"]["structure"],
                         "confidence_sha256": selection_report["selected_sha256"]["confidence"],
                         "ptm": selected_row.get("ptm"), "ranking_score": selected_row.get("ranking_score"),
                         "mean_plddt": new_metrics["mean_plddt"],
                         "consensus_score": selected_row.get("consensus_score")}
        _, request_hash = _metadata(source)
        after = _source_inventory(source)
        if before != after:
            raise RuntimeError(f"{binder}: source directory changed during comparison")
        result = {
            "binder": binder, "old": old_selection, "new": new_selection,
            "geometry_reused_from_previous": None,
            "selection_report": selection_report, "selection_report_path": str(report_path),
            "frozen_manifest": str(af3_manifest), "original_metadata": request_hash,
            "source_files_unchanged": True, "source_inventory": before,
            "mapping": {
                "old_identity": old_metrics["mapping"].identity,
                "new_identity": new_metrics["mapping"].identity,
                "old_coverage": old_metrics["mapping"].observed_ecd_coverage,
                "new_coverage": new_metrics["mapping"].observed_ecd_coverage,
                "old_mismatches": old_metrics["mapping"].mismatch_positions,
                "new_mismatches": new_metrics["mapping"].mismatch_positions,
            },
            "geometry": geometry,
            "confidence_delta": {
                "ptm": (new_selection["ptm"] - old_selection["ptm"])
                        if new_selection["ptm"] is not None and old_selection["ptm"] is not None else None,
                "mean_plddt": (new_selection["mean_plddt"] - old_selection["mean_plddt"])
                        if new_selection["mean_plddt"] is not None and old_selection["mean_plddt"] is not None else None,
                "ranking_score": (new_selection["ranking_score"] - old_selection["ranking_score"])
                        if new_selection["ranking_score"] is not None and old_selection["ranking_score"] is not None else None,
                "consensus_score": (new_selection["consensus_score"] - old_selection["consensus_score"])
                        if new_selection["consensus_score"] is not None and old_selection["consensus_score"] is not None else None,
            },
            "geometry_metrics_methods": {
                "surface_parameters": config.SURFACE_V2_PARAMETERS,
                "reliability_parameters": "classify_reliability defaults",
                "accessibility": "analyze_residues native Shrake-Rupley SASA; current default thresholds",
                "contacts": "mapped C-alpha pairs <=8 A and sequence separation >2",
            },
        }
        write_json(case_dir / "case_audit.json", result)
        cases.append(result)
        compact.append({
            "binder": binder, "old_model": "model_0", "new_model": new_selection["model_id"],
            "selection_changed": new_selection["source_sha256"] != old_selection["source_sha256"],
            "old_ptm": old_selection["ptm"], "new_ptm": new_selection["ptm"],
            "delta_ptm": result["confidence_delta"]["ptm"],
            "old_mean_plddt": old_selection["mean_plddt"], "new_mean_plddt": new_selection["mean_plddt"],
            "delta_mean_plddt": result["confidence_delta"]["mean_plddt"],
            "old_ranking_score": old_selection["ranking_score"], "new_ranking_score": new_selection["ranking_score"],
            "delta_ranking_score": result["confidence_delta"]["ranking_score"],
            "delta_consensus_score": result["confidence_delta"]["consensus_score"],
            "delta_mean_rsasa": (geometry["mean_rsasa_new"] - geometry["mean_rsasa_old"])
                               if geometry["mean_rsasa_new"] is not None and geometry["mean_rsasa_old"] is not None else None,
            "matched_ca": geometry["matched_ca_count"], "highconf_core_ca": geometry["high_confidence_core_count"],
            "global_rmsd_a": geometry["global_ca_rmsd_a"], "highconf_core_rmsd_a": geometry["matched_high_confidence_core_rmsd_a"],
            "median_ca_displacement_a": geometry["displacement_median_a"], "p90_ca_displacement_a": geometry["displacement_p90_a"],
            "contact_jaccard": geometry["contact_jaccard"], "contacts_gained": geometry["contacts_gained"], "contacts_lost": geometry["contacts_lost"],
            "old_eligible": geometry["old_eligible_count"], "new_eligible": geometry["new_eligible_count"],
            "eligible_gain": len(geometry["eligible_gain"]), "eligible_loss": len(geometry["eligible_loss"]),
            "old_regions": geometry["old_region_count"], "new_regions": geometry["new_region_count"],
            "region_membership_jaccard": geometry["region_membership_jaccard"], "region_membership_ari": geometry["region_membership_ari"],
            "region_best_match_jaccard_mean": geometry["region_best_match_jaccard_mean"],
        })
    protocol = protocol_start
    protocol["source_files_unchanged_during_run"] = protocol_start == _protocol_fingerprints()
    protocol["python"] = sys.version
    protocol["platform"] = platform.platform()
    protocol["command"] = "direct canonical-mapping + current analyze_residues/classify_reliability/run_surface_patch_v2; no pipeline mutation scoring"
    result = {
        "schema": "af3_preprocessing_comparison_v1", "created_at_utc": datetime.now(UTC).isoformat(),
        "source_manifest": str(MANIFEST), "source_manifest_sha256": sha256_file(MANIFEST),
        "protocol_fingerprints": protocol, "case_count": len(cases),
        "changed_selection_count": sum(row["selection_changed"] for row in compact),
        "cases": cases,
        "interpretation": "Geometry and confidence diagnostics only. Historical mutation scores were not re-used for the newly selected conformer. No experimental truth, new SPURS/FoldX scores, or validated mutation recommendations are established.",
    }
    if not protocol["source_files_unchanged_during_run"]:
        raise RuntimeError("Protocol source changed during comparison; rerun for a stable SHA-bound result")
    write_json(out / "comparison.json", result)
    _write_csv(out / "comparison.csv", compact)
    write_json(out / "protocol_fingerprints.json", protocol)
    _write_report(out / "report_zh.md", compact, result)
    print(out)
    return out


def _write_report(path: Path, rows: list[dict], result: dict):
    header = "| Binder | 新选构象 | ΔpTM | Δranking | Δ平均pLDDT | Δ共识分数 | 核心RMSD Å | 中位位移 Å | 接触Jaccard | 合格暴露位点旧→新 | 区域旧→新 | 成员Jaccard / ARI / 双向最佳匹配 |"
    lines = [
        "# AF3 model_0 与预处理选模：五例配对几何比较",
        "",
        "本比较对五个binder分别冻结全部原始AF3模型与各自confidence文件，再按当前预处理流程选择实际来源构象。每例以原有model_0为旧条件；两边均按相同序列映射、当前SASA/可靠性门与Surface V2几何参数重算。以下变化说明输入构象如何改变当前几何描述。",
        "",
        f"五例中有 **{result['changed_selection_count']}/5** 例实际选择了不同来源模型。",
        "",
        header,
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['binder']} | {row['new_model']} | {_fmt(row['delta_ptm'])} | {_fmt(row['delta_ranking_score'])} | "
            f"{_fmt(row['delta_mean_plddt'])} | {_fmt(row['delta_consensus_score'])} | "
            f"{_fmt(row['highconf_core_rmsd_a'])} | {_fmt(row['median_ca_displacement_a'])} | {_fmt(row['contact_jaccard'])} | "
            f"{row['old_eligible']}→{row['new_eligible']} | {row['old_regions']}→{row['new_regions']} | "
            f"{_fmt(row['region_membership_jaccard'])} / {_fmt(row['region_membership_ari'])} / "
            f"{_fmt(row['region_best_match_jaccard_mean'])} |"
        )
    lines.extend([
        "",
        "核心RMSD按两边均有pLDDT≥70的共同Cα对齐计算；中位位移基于该核心对齐。逐案JSON保留全局Cα RMSD及逐残基位移、pLDDT、SASA和相对SASA差值。接触为映射后Cα距离≤8 Å且序列间隔>2。区域成员指标把当前几何区域投影到规范残基编号后比较；它不代表抗体结合界面。合格暴露位点是原生SASA可及且未被当前可靠性严格模式排除的残基集合。",
        "",
        "区域ARI仅在两侧成员各自构成互斥分区时计算；若区域有成员重叠，ARI记为NA，并用旧→新与新→旧的最佳成员Jaccard均值作匹配摘要。共识分数仅在两构象都进入可比较集合并具有分数时报告差值。",
        "",
        "旧模型与新选模型使用各自同索引confidence；pTM、平均pLDDT和ranking score均来自各自原始summary/full-data元数据。源目录在比较前后逐文件SHA及大小核对，任何变化都会中止。每例`case_audit.json`保留逐残基位移、增减的接触/位点、各区域成员和选模诊断。",
        "",
        "本结果仅说明预处理选模对结构几何、暴露和区域分区的影响。它不证明相对实验结构或真实表位更准确；没有为新构象调用SPURS或FoldX，也没有把旧构象的突变分数移作新构象证据，因此不构成验证过的突变建议。",
        "",
        "逐案数据见同目录`comparison.json`和`comparison.csv`。",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _fmt(value):
    return "NA" if value is None else f"{value:.3f}"


if __name__ == "__main__":
    run()
