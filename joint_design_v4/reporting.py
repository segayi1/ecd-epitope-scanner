from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator

from ..util import sha256_bytes, sha256_file, write_json
from .models import JointV4Result, PredictionResult


def _rows(value: Any) -> Any:
    if is_dataclass(value):
        return {key: _rows(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): _rows(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_rows(item) for item in value]
    return value


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_gzip_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _code_snapshot() -> dict[str, Any]:
    root = Path(__file__).resolve().parents[1]
    paths = list(root.glob("*.py"))
    for directory in ("joint_design_v4", "joint_design", "integrations", "surface_patch"):
        paths.extend((root / directory).rglob("*.py"))
    digest = hashlib.sha256()
    for path in sorted(set(paths)):
        if path.is_file():
            digest.update(path.relative_to(root).as_posix().encode("utf-8"))
            digest.update(path.read_bytes())
    snapshot: dict[str, Any] = {"source_code_sha256": digest.hexdigest()}
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, timeout=5
        )
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=all"],
            cwd=root, capture_output=True, text=True, timeout=5,
        )
        snapshot["git_commit"] = revision.stdout.strip() if revision.returncode == 0 else None
        snapshot["git_dirty"] = bool(status.stdout.strip()) if status.returncode == 0 else None
    except (OSError, subprocess.TimeoutExpired):
        snapshot.update({"git_commit": None, "git_dirty": None})
    return snapshot


def write_v4_outputs(
    output_dir: Path,
    result: JointV4Result,
    *,
    parameters: dict,
    contexts: dict,
    source_provenance: dict | None = None,
    worker_results: dict[str, PredictionResult] | None = None,
    tool_commands: dict[str, list[str]] | None = None,
    normalized_structure_path: Path | None = None,
) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "tool_inventory.json", _rows(result.tool_inventory))
    write_json(output_dir / "structure_contexts.json", contexts)
    write_json(output_dir / "probe_sectors_v4.json", _rows(result.sectors))
    write_json(output_dir / "panel_candidates.json", _rows(result.candidates))
    write_json(output_dir / "selected_panel_v4.json", _rows(result.panel))
    surface_rows = [_rows(item) for item in result.surface_evidence]
    _write_csv(
        output_dir / "surface_evidence.csv",
        surface_rows,
        [
            "residue",
            "procaliper_sasa_a2",
            "procaliper_rsasa",
            "native_sasa_a2",
            "native_rsasa",
            "sidechain_sasa_a2",
            "mapped_surface_area_a2",
            "surface_component_id",
            "sasa_agreement_status",
        ],
    )
    _write_csv(
        output_dir / "sasa_comparison.csv",
        [
            {
                "position": item.residue.uniprot_pos,
                "native_rsasa": item.native_rsasa,
                "procaliper_rsasa": item.procaliper_rsasa,
                "absolute_delta": (
                    abs(item.native_rsasa - item.procaliper_rsasa)
                    if item.native_rsasa is not None and item.procaliper_rsasa is not None
                    else None
                ),
                "status": item.sasa_agreement_status,
            }
            for item in result.surface_evidence
        ],
        ["position", "native_rsasa", "procaliper_rsasa", "absolute_delta", "status"],
    )
    _write_csv(
        output_dir / "discotope_residues.csv",
        [_rows(item) for item in result.epitope_evidence],
        [
            "residue",
            "raw_score",
            "calibrated_score",
            "surface_percentile",
            "local_peak",
            "predicted_epitope_label",
            "island_id",
            "status",
        ],
    )
    _write_csv(
        output_dir / "single_mutation_options.csv",
        [_rows(item) for item in result.mutation_options],
        [
            "residue",
            "mutant_aa",
            "spurs_ddg_kcal_mol",
            "stability_band",
            "chemistry_class",
            "eligible",
        ],
    )
    _write_gzip_csv(
        output_dir / "spurs_single_matrix.csv.gz",
        [
            {
                "position": item.residue.uniprot_pos,
                "wt": item.residue.one_letter,
                "mutant_aa": item.mutant_aa,
                "ddg": item.spurs_ddg_kcal_mol,
                "band": item.stability_band,
            }
            for item in result.mutation_options
        ],
        ["position", "wt", "mutant_aa", "ddg", "band"],
    )
    islands: dict[str, list[int]] = {}
    for item in result.epitope_evidence:
        if item.island_id:
            islands.setdefault(item.island_id, []).append(item.residue.uniprot_pos or 0)
    write_json(output_dir / "epitope_islands.json", islands)
    _write_gzip_csv(
        output_dir / "candidate_combinations.csv.gz",
        [
            {
                "candidate_id": item.candidate_id,
                "sector_id": item.sector_id,
                "order": item.order,
                "mutations": ";".join(option.notation for option in item.mutations),
                "epitope_score": item.epitope_score,
                "epitope_supported_sites": item.epitope_supported_sites,
                "recommendation_class": item.recommendation_class,
                "structural_risk": item.structural_risk,
                "role": item.role,
                "formal_eligible": item.formal_eligible,
                "formal_exclusion_reasons": ";".join(item.formal_exclusion_reasons),
                "stability": _rows(item.stability),
            }
            for item in result.candidates
        ],
        [
            "candidate_id",
            "sector_id",
            "order",
            "mutations",
            "epitope_score",
            "epitope_supported_sites",
            "recommendation_class",
            "structural_risk",
            "role",
            "formal_eligible",
            "formal_exclusion_reasons",
            "stability",
        ],
    )
    _write_csv(
        output_dir / "combination_stability.csv",
        [
            {
                "candidate_id": item.candidate_id,
                "mutations": ";".join(option.notation for option in item.mutations),
                **(_rows(item.stability) or {}),
            }
            for item in result.candidates
        ],
        [
            "candidate_id",
            "mutations",
            "spurs_multi_ddg_kcal_mol",
            "foldx_ddg_kcal_mol",
            "order",
            "epistasis_estimate",
            "extrapolation_level",
            "consensus_status",
            "disagreement",
        ],
    )
    manifest = {
            "schema_version": "4.0",
            "engine": "joint_v4_ml",
            "recommendation_status": result.recommendation_status,
            "tools": _rows(result.tool_inventory),
            "parameters": parameters,
            "contexts": contexts,
            "provenance": source_provenance or {},
            "software": _code_snapshot(),
            "normalized_structure_sha256": (
                sha256_file(normalized_structure_path)
                if normalized_structure_path is not None and normalized_structure_path.is_file()
                else None
            ),
            "tool_commands": tool_commands or {},
            "worker_runs": {
                name: {
                    "status": worker.status,
                    "record_count": len(worker.records),
                    "warnings": worker.warnings,
                    "errors": worker.errors,
                    "provenance": worker.provenance,
                    "runtime_seconds": worker.runtime_seconds,
                }
                for name, worker in (worker_results or {}).items()
            },
            "warnings": result.warnings,
            "errors": result.errors,
            "cache_key": sha256_bytes(json.dumps(parameters, sort_keys=True).encode()),
    }
    schema_path = Path(__file__).resolve().parents[1] / "schemas" / "joint_v4_manifest.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    Draft202012Validator(schema).validate(manifest)
    write_json(output_dir / "v4_manifest.json", manifest)
    report = [
        "# Joint V4-ML scientific report",
        "",
        f"- Recommendation status: `{result.recommendation_status}`",
        f"- Sectors: {len(result.sectors)}; candidates: {len(result.candidates)}",
        f"- Recommended series: {len(result.panel.selected_candidate_ids)} combinations",
        f"- Global surface coverage denominator: {result.panel.coverage_denominator_residues} mapped residues",
        f"- Global union coverage by geodesic scenario: {result.panel.achieved_coverage}",
        f"- Epitope-prioritized: {len(result.panel.candidates_by_role.get('EPITOPE_PRIORITIZED', []))}",
        f"- Stability fallback: {len(result.panel.candidates_by_role.get('STABILITY_FALLBACK', []))}",
        "",
        "## Interpretation boundaries",
        "",
        "DiscoTope scores are soft epitope-prior evidence, not assertions of true antibody contacts.",
        "SPURS/FoldX results are predicted stability evidence and do not ensure expression, folding, or cell-surface display.",
        "Preview outputs are not formal recommendations.",
        "",
        "## Warnings",
        "",
        *[f"- {item}" for item in result.warnings],
    ]
    (output_dir / "v4_scientific_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return {path.name: path.name for path in output_dir.iterdir() if path.is_file()}
