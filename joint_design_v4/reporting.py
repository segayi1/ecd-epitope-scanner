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

from ..util import portable_command, sha256_bytes, sha256_file, write_json
from .models import JointV4Result, PredictionResult


def _rows(value: Any) -> Any:
    if is_dataclass(value):
        return {key: _rows(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): _rows(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_rows(item) for item in value]
    return value


def _portable_tool_commands(commands: dict[str, list[str]]) -> dict[str, list[str]]:
    return {name: portable_command(command) for name, command in commands.items()}


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


def _code_snapshot(effective_parameters: dict | None = None) -> dict[str, Any]:
    root = Path(__file__).resolve().parents[1]
    paths = list(root.glob("*.py"))
    for directory in ("joint_design_v4", "joint_design", "integrations", "surface_patch", "surface_perturbation"):
        paths.extend((root / directory).rglob("*.py"))
    digest = hashlib.sha256()
    source_files = []
    for path in sorted(set(paths)):
        if path.is_file():
            content = path.read_bytes()
            relative = path.relative_to(root).as_posix()
            digest.update(relative.encode("utf-8"))
            digest.update(content)
            source_files.append({"path": relative, "sha256": hashlib.sha256(content).hexdigest(),
                                 "size_bytes": len(content)})
    lock_names = {"pyproject.toml", "requirements.txt", "requirements.lock", "uv.lock", "poetry.lock",
                  "Pipfile.lock", "environment.yml", "environment.yaml"}
    dependency_files = [{"path": path.name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
                        for path in sorted(root.iterdir()) if path.is_file() and path.name in lock_names]
    snapshot: dict[str, Any] = {"source_code_sha256": digest.hexdigest(),
                                "source_files_sha256": sha256_bytes(json.dumps(source_files, sort_keys=True,
                                    separators=(",", ":")).encode("utf-8")),
                                "source_files": source_files,
                                "dependency_files": dependency_files,
                                "effective_parameters_sha256": sha256_bytes(json.dumps(
                                    effective_parameters or {}, sort_keys=True, separators=(",", ":"),
                                    default=str).encode("utf-8"))}
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
        snapshot["git_dirty_paths"] = [line[3:] for line in status.stdout.splitlines() if len(line) >= 4]
    except (OSError, subprocess.TimeoutExpired):
        snapshot.update({"git_commit": None, "git_dirty": None})
    return snapshot


def _model_evidence_diagnostics(result: JointV4Result, worker_results: dict[str, PredictionResult] | None) -> dict[str, Any]:
    foldx_records = (worker_results or {}).get("foldx")
    foldx_by_mutations = {}
    if foldx_records:
        for row in foldx_records.records:
            mutations = row.get("mutations", [])
            if isinstance(mutations, (list, tuple)):
                foldx_by_mutations[tuple(sorted(map(str, mutations)))] = row
    details, grouped = [], {}
    for candidate in result.candidates:
        stability = candidate.stability
        spurs = stability.spurs_multi_ddg_kcal_mol if stability else None
        foldx = stability.foldx_ddg_kcal_mol if stability else None
        window = ("MISSING" if spurs is None else "LE_0_5" if spurs <= .5 else
                  "GT_0_5_LE_1" if spurs <= 1 else "GT_1_LE_2" if spurs <= 2 else "GT_2")
        row = foldx_by_mutations.get(tuple(sorted(option.notation for option in candidate.mutations)), {})
        pairs = row.get("paired_models", [])
        repeat_audit = candidate.impact_evidence.get("stability_repeat_audit", {})
        values = repeat_audit.get("foldx_repeat_values_kcal_mol", row.get("ddg_runs_kcal_mol")) or []
        identities = [(item.get("batch_path"), item.get("group"), item.get("run")) for item in pairs]
        detail = {"candidate_id": candidate.candidate_id, "order": candidate.order,
                  "spurs_window": window, "spurs_multi_ddg_kcal_mol": spurs,
                  "foldx_ddg_kcal_mol": foldx, "formal_eligible": bool(candidate.formal_eligible),
                  "foldx_repeat_count": len(values), "foldx_repeat_values_kcal_mol": values,
                  "signed_model_difference_spurs_minus_foldx": spurs-foldx if spurs is not None and foldx is not None else None,
                  "absolute_model_difference": abs(spurs-foldx) if spurs is not None and foldx is not None else None,
                  "foldx_repeat_median": __import__("statistics").median(values) if values else None,
                  "foldx_repeat_range": max(values)-min(values) if values else None,
                  "foldx_repeat_sample_sd": __import__("statistics").stdev(values) if len(values)>1 else None,
                  "foldx_repeat_gate_kcal_mol": repeat_audit.get("foldx_repeat_gate_kcal_mol"),
                  "foldx_repeat_min_kcal_mol": repeat_audit.get("foldx_repeat_min_kcal_mol"),
                  "foldx_repeat_max_kcal_mol": repeat_audit.get("foldx_repeat_max_kcal_mol"),
                  "foldx_repeat_above_gate_count": repeat_audit.get("foldx_repeat_above_gate_count"),
                  "foldx_repeat_at_or_below_gate_count": repeat_audit.get("foldx_repeat_at_or_below_gate_count"),
                  "sampling_decision_diagnostic_status": repeat_audit.get("sampling_decision_diagnostic_status"),
                  "sampling_decision_sensitive": repeat_audit.get("sampling_decision_sensitive"),
                  "recommendation_risk_label": repeat_audit.get("recommendation_risk_label"),
                  "cross_execution_conflict": repeat_audit.get("cross_execution_conflict"),
                  "cross_execution_evidence_status": repeat_audit.get("cross_execution_evidence_status", "NOT_ASSESSED"),
                  "repeat_audit": repeat_audit,
                  "scheduling": candidate.impact_evidence.get("foldx_scheduling", {}),
                  "foldx_pair_identity_count": len({tuple(p) for p in repeat_audit.get("foldx_pair_identities", [])}) if repeat_audit.get("foldx_pair_identities") else len(set(identities)),
                  "foldx_duplicate_content_count": repeat_audit.get("foldx_duplicate_content_count"),
                  "foldx_missing_reason": None if foldx is not None else (
                      "NOT_IN_WORKER_RECORDS" if row else "NOT_SCHEDULED_OR_NO_WORKER_EVIDENCE")}
        details.append(detail)
        key = f"order_{candidate.order}/{window}"
        bucket = grouped.setdefault(key, {"order": candidate.order, "spurs_window": window,
            "candidate_count": 0, "spurs_scored_count": 0, "foldx_scored_count": 0,
            "formal_candidate_count": 0, "repeat_observation_count": 0,
            "candidate_denominator": "all candidates in this saved run and order/window"})
        bucket["candidate_count"] += 1
        bucket["spurs_scored_count"] += spurs is not None
        bucket["foldx_scored_count"] += foldx is not None
        bucket["formal_candidate_count"] += bool(candidate.formal_eligible)
        bucket["repeat_observation_count"] += len(values)
    return {"schema_version": "model_evidence_diagnostics_v1",
            "by_order_and_spurs_window": list(grouped.values()), "candidates": details,
            "foldx_worker_runtime_seconds": foldx_records.runtime_seconds if foldx_records else None,
            "interpretation": "Saved prediction evidence and run repetitions only; repeat spread is not an experimental confidence interval."}


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
    write_json(output_dir / "execution_protocol.json", result.execution_protocol)
    write_json(output_dir / "ordinary_candidate_ledger.json",
               result.provenance.get("global_candidate_pool", {}).get(
                   "adaptive_candidate_ledger", {}))
    write_json(output_dir / "tool_inventory.json", _rows(result.tool_inventory))
    write_json(output_dir / "structure_contexts.json", contexts)
    write_json(output_dir / "probe_sectors_v4.json", _rows(result.sectors))
    write_json(output_dir / "panel_candidates.json", _rows(result.candidates))
    write_json(output_dir / "selected_panel_v4.json", _rows(result.panel))
    sector_recommendations = _rows(getattr(result.panel, "sector_recommendations", []))
    write_json(output_dir / "sector_recommendations_v4.json", sector_recommendations)
    surface_threshold = parameters.get("surface", {}).get("rsasa_threshold")
    surface_gate = {
        "policy": "native_total_rsasa_only",
        "native_total_rsasa_threshold": surface_threshold,
        "sidechain_sasa_role": "descriptive_only_not_surface_gate",
        "geometry_eligibility": "eligible design nodes require a retained solvent-surface geometry owner",
        "interpretation": "This records the actual V4 surface accessibility rule; other coordinate, reliability and mutation gates remain separate.",
    }
    software = _code_snapshot(parameters)
    write_json(output_dir / "run_source_identity_v1.json", {
        "schema_version": "run_source_identity_v1", "software": software,
        "effective_parameters": parameters, "surface_gate": surface_gate,
        "normalized_structure_sha256": (sha256_file(normalized_structure_path)
            if normalized_structure_path is not None and normalized_structure_path.is_file() else None),
        "worker_model_identities": {name: worker.provenance for name, worker in (worker_results or {}).items()},
        "interpretation": "Identity captured during result export; it cannot reconstruct identities omitted by a historical run.",
    })
    write_json(output_dir / "model_evidence_diagnostics_v1.json",
               _model_evidence_diagnostics(result, worker_results))
    recommendation_rows = [item for row in sector_recommendations for item in [row, *row.get("alternatives", [])]]
    _write_csv(output_dir / "sector_recommendations_v4.csv", [
        {key: json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list)) else value
         for key, value in row.items()} for row in recommendation_rows
    ], ["sector_id", "candidate_id", "status", "formal_eligible", "failed_gates", "gate_details",
        "gate_distance", "coverage_gain", "coverage_gain_area_a2", "geometry_support", "regional_geometry_support", "mutation_positions_in_sector",
        "mutation_notations", "selected_in_panel", "design_id", "recommendation_role", "ownership", "principles",
        "recommendation_reasons", "recommendation_risk_label", "recommendation_risk_explanation", "stability_repeat_diagnostic", "gate_checks", "site_details", "stability_evidence", "pair_distances", "provenance", "primary_comparison", "validation_status",
        "recommendation_policy", "impact_evidence", "sensitivity", "resource_audit", "unavailable_principles", "direct_area_a2", "local_12_covered_area_a2",
        "local_12_coverage_fraction", "local_12_uncovered_area_a2", "local_12_uncovered_positions"])
    if parameters.get("optimization_profile") in {"coverage_scan_v2_preview", "coverage_scan_v3_preview"}:
        menu_rows = [_rows(item) for item in result.mutation_options]
        write_json(output_dir / "mutation_menu_v2.json", {
            "schema_version": "coverage_menu_2.0", "optimization_profile": parameters["optimization_profile"],
            "interpretation": "Chemical scheduling descriptors; not measured surface perturbation or binding evidence",
            "rows": menu_rows})
        _write_csv(output_dir / "mutation_menu_v2.csv", [
            {**row, "residue": json.dumps(row["residue"], sort_keys=True),
             "chemical_descriptors": json.dumps(row["chemical_descriptors"], sort_keys=True),
             "hard_exclusions": json.dumps(row["hard_exclusions"])} for row in menu_rows],
            ["residue", "mutant_aa", "spurs_ddg_kcal_mol", "stability_band", "eligible",
             "hard_exclusions", "chemical_descriptors", "pareto_layer", "menu_selected", "selection_reason"])
    extension_names = [
        "qualification_audit_v411.json", "hypothesis_hierarchy_v411.json",
        "global_candidate_pool_v411.json", "review_candidates_v411.json",
        "panel_budget_decision_v411.json", "panel_budget_curve_v411.csv",
        "raw_hypotheses_v411.json", "footprint_families_v411.json",
        "inference_hypotheses_v411.json", "surface_complexity_v411.json",
        "retrospective_benchmark_contract_v411.json",
        "raw_surface_footprints_v411.json",
        "sector_recommendations_v4.json", "sector_recommendations_v4.csv",
        "run_source_identity_v1.json", "model_evidence_diagnostics_v1.json",
    ]
    if parameters.get("optimization_profile") in {"coverage_scan_v2_preview", "coverage_scan_v3_preview"}:
        extension_names.extend(["mutation_menu_v2.json", "mutation_menu_v2.csv"])
    write_json(output_dir / extension_names[0], result.provenance.get("qualification_audit", {}))
    write_json(output_dir / extension_names[1], {
        "interpretation": "Raw -> family -> inference; engineering membership weights are not probabilities.",
        "parameters": parameters.get("hypothesis_hierarchy", {}),
        "families": _rows(result.footprint_families),
        "inference_hypotheses": _rows(result.surface_hypotheses),
        "raw_to_family": {fid: f.family_id for f in result.footprint_families for fid in f.member_footprint_ids},
        "unassigned_raw_footprint_ids": sorted(set(result.provenance.get("hierarchy_raw_footprint_ids", [])) - {fid for f in result.footprint_families for fid in f.member_footprint_ids}),
        "family_to_inference": {fid: h.hypothesis_id for h in result.surface_hypotheses for fid in getattr(h, "family_ids", [])},
        "resolution_by_level": result.panel.hierarchy_resolution,
    })
    write_json(output_dir / extension_names[2], {"formal": result.provenance.get("global_candidate_pool", {}), "review": result.provenance.get("review_candidate_pool", {})})
    write_json(output_dir / extension_names[3], _rows([c for c in result.candidates if c.qualification == "REVIEW_REQUIRED"]))
    write_json(output_dir / extension_names[4], _rows(result.panel.budget_decision))
    _write_csv(output_dir / extension_names[5], result.panel.budget_decision.budget_curve if result.panel.budget_decision else [],
               ["construct_count", "candidate_id", "objective", "net_gain", "information_gain", "relative_information_gain", "coverage_gain", "mean_separation_gain", "minimum_separation_gain", "unresolved_pair_reduction", "risk_delta", "redundancy_delta"])
    write_json(output_dir / extension_names[6], _rows(result.raw_hypotheses))
    write_json(output_dir / extension_names[7], _rows(result.footprint_families))
    write_json(output_dir / extension_names[8], _rows(result.surface_hypotheses))
    write_json(output_dir / extension_names[9], _rows(result.complexity))
    from .benchmark import benchmark_contract
    write_json(output_dir / extension_names[10], benchmark_contract())
    write_json(output_dir / "surface_hypotheses_v4.json", _rows(result.surface_hypotheses))
    write_json(
        output_dir / "epitope_hypotheses_v41.json",
        {"schema_version": "4.1", "interpretation": "Candidate antibody-compatible surface footprints, not established epitopes.",
         "hypotheses": _rows(result.surface_hypotheses)},
    )
    membership_rows = [
        {
            "hypothesis_id": hypothesis.hypothesis_id,
            "uniprot_pos": position,
            "surface_area_in_hypothesis_a2": hypothesis.residue_weights.get(position, 0.0),
            "membership_fraction": hypothesis.membership_fraction.get(position, 0.0),
            "core_or_boundary": (
                "CORE" if position in hypothesis.core_positions else
                "BOUNDARY" if position in hypothesis.boundary_positions else "NOMINAL"
            ),
        }
        for hypothesis in result.surface_hypotheses
        for position in sorted(set(hypothesis.residue_weights) | set(hypothesis.boundary_positions))
    ]
    _write_csv(
        output_dir / "hypothesis_membership_v41.csv", membership_rows,
        ["hypothesis_id", "uniprot_pos", "surface_area_in_hypothesis_a2",
         "membership_fraction", "core_or_boundary"],
    )
    write_json(
        output_dir / "hypothesis_probe_matrix_v4.json",
        {
            "interpretation": "Geometry and chemistry perturbation proxy; not a binding-change prediction or epitope probability.",
            "model": result.provenance.get("hypothesis_response_model", {}),
            "hypothesis_ids": [item.hypothesis_id for item in result.surface_hypotheses],
            "candidate_responses": {
                item.candidate_id: item.hypothesis_perturbation or item.hypothesis_response
                for item in result.candidates
            },
        },
    )
    _write_csv(
        output_dir / "probe_hypothesis_matrix_v41.csv",
        [
            {"candidate_id": item.candidate_id, "hypothesis_id": hypothesis.hypothesis_id,
             "perturbation_proxy": (item.hypothesis_perturbation or item.hypothesis_response).get(hypothesis.hypothesis_id, 0.0)}
            for item in result.candidates for hypothesis in result.surface_hypotheses
        ],
        ["candidate_id", "hypothesis_id", "perturbation_proxy"],
    )
    write_json(
        output_dir / "panel_identifiability_v41.json",
        {
            "schema_version": "4.1",
            "interpretation": "Geometry and chemistry perturbation signatures; separation is not experimental binding evidence.",
            "selected_candidate_ids": result.panel.selected_candidate_ids,
            "uniform_hypothesis_response": result.panel.mean_hypothesis_response,
            "prior_weighted_hypothesis_response": result.panel.prior_weighted_hypothesis_response,
            "mean_pairwise_separation": result.panel.mean_pairwise_separation,
            "prior_weighted_mean_pairwise_separation": result.panel.prior_weighted_mean_pairwise_separation,
            "minimum_signature_distance": result.panel.minimum_signature_distance,
            "marginal_information_gain": result.panel.marginal_information_gain,
            "objective": {
                "coverage": "equal blend of uniform and DiscoTope-prior weighted hypothesis coverage",
                "separation": "equal blend of uniform and prior-weighted pairwise signature separation",
                "minimum_separation": "additional weight on the closest hypothesis pair",
                "penalties": "separate stability, glycan, disulfide, interface, uncertainty, and probe-redundancy tie breaks",
            },
            "unresolved_hypothesis_groups": result.panel.unresolved_hypothesis_groups,
            "uncovered_hypotheses": result.panel.uncovered_hypotheses,
            "hypothesis_signatures": result.panel.hypothesis_signatures,
            "identifiability_diagnostics": result.panel.identifiability_diagnostics,
            "interface_comparison": result.panel.interface_comparison,
            "surface_gate": surface_gate,
            "model_evidence_diagnostics_file": "model_evidence_diagnostics_v1.json",
        },
    )
    exposure_by_position = {
        item.residue.uniprot_pos: item.impact_evidence.get("sidechain_exposure", {})
        for item in result.mutation_options
        if item.impact_evidence.get("sidechain_exposure")
    }
    surface_rows = [
        {**_rows(item), **{
            f"sidechain_exposure_{key}": exposure_by_position.get(item.residue.uniprot_pos, {}).get(key)
            for key in ("value", "numerator_a2", "denominator_a2", "status", "method")
        }}
        for item in result.surface_evidence
    ]
    _write_csv(
        output_dir / "surface_evidence.csv",
        surface_rows,
        [
            "residue",
            "procaliper_sasa_a2",
            "procaliper_rsasa",
            "native_sasa_a2",
            "native_rsasa",
            "normalized_context_sasa_a2",
            "normalized_context_rsasa",
            "native_context_absolute_delta",
            "native_context_status",
            "native_context_classification_changed",
            "procaliper_classification_changed",
            "comparison_method_status",
            "sidechain_sasa_a2",
            "mapped_surface_area_a2",
            "surface_component_id",
            "sasa_agreement_status",
            "sidechain_exposure_value",
            "sidechain_exposure_numerator_a2",
            "sidechain_exposure_denominator_a2",
            "sidechain_exposure_status",
            "sidechain_exposure_method",
        ],
    )
    _write_csv(
        output_dir / "sasa_comparison.csv",
        [
            {
                "position": item.residue.uniprot_pos,
                "native_rsasa": item.native_rsasa,
                "normalized_context_rsasa": item.normalized_context_rsasa,
                "procaliper_rsasa": item.procaliper_rsasa,
                "absolute_delta": (
                    abs(item.normalized_context_rsasa - item.procaliper_rsasa)
                    if item.normalized_context_rsasa is not None
                    and item.procaliper_rsasa is not None
                    else None
                ),
                "native_context_absolute_delta": item.native_context_absolute_delta,
                "native_context_status": item.native_context_status,
                "native_context_classification_changed": item.native_context_classification_changed,
                "procaliper_classification_changed": item.procaliper_classification_changed,
                "comparison_method_status": item.comparison_method_status,
                "status": item.sasa_agreement_status,
            }
            for item in result.surface_evidence
        ],
        [
            "position",
            "native_rsasa",
            "normalized_context_rsasa",
            "procaliper_rsasa",
            "absolute_delta",
            "native_context_absolute_delta",
            "native_context_status",
            "native_context_classification_changed",
            "procaliper_classification_changed",
            "comparison_method_status",
            "status",
        ],
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
        [{**_rows(item),
          "chemical_descriptors": json.dumps(item.chemical_descriptors, ensure_ascii=False, sort_keys=True),
          "impact_evidence": json.dumps(item.impact_evidence, ensure_ascii=False, sort_keys=True)}
         for item in result.mutation_options],
        [
            "residue",
            "mutant_aa",
            "spurs_ddg_kcal_mol",
            "stability_band",
            "chemistry_class",
            "chemical_descriptors",
            "impact_evidence",
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
                "anchor_region_id": item.anchor_region_id or item.sector_id,
                "order": item.order,
                "mutations": ";".join(option.notation for option in item.mutations),
                "targeted_hypotheses": ";".join(item.targeted_hypotheses),
                "hypothesis_perturbation": json.dumps(item.hypothesis_perturbation or item.hypothesis_response, sort_keys=True),
                "risk_components": json.dumps(item.risk_components, sort_keys=True),
                "geometry_support": json.dumps(item.geometry_support, sort_keys=True),
                "epitope_score": item.epitope_score,
                "epitope_supported_sites": item.epitope_supported_sites,
                "recommendation_class": item.recommendation_class,
                "structural_risk": item.structural_risk,
                "role": item.role,
                "triple_role": item.triple_role,
                "parent_candidate_ids": ";".join(item.parent_candidate_ids),
                "control_for_candidate_ids": ";".join(item.control_for_candidate_ids),
                "formal_eligible": item.formal_eligible,
                "formal_exclusion_reasons": ";".join(item.formal_exclusion_reasons),
                "stability": _rows(item.stability),
            }
            for item in result.candidates
        ],
        [
            "candidate_id",
            "sector_id",
            "anchor_region_id",
            "order",
            "mutations",
            "targeted_hypotheses",
            "hypothesis_perturbation",
            "risk_components",
            "geometry_support",
            "epitope_score",
            "epitope_supported_sites",
            "recommendation_class",
            "structural_risk",
            "role",
            "triple_role",
            "parent_candidate_ids",
            "control_for_candidate_ids",
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
            "software": software,
            "surface_gate": surface_gate,
            "normalized_structure_sha256": (
                sha256_file(normalized_structure_path)
                if normalized_structure_path is not None and normalized_structure_path.is_file()
                else None
            ),
            "tool_commands": _portable_tool_commands(tool_commands or {}),
            "tool_command_paths_redacted": True,
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
    diagnostics = result.provenance.get("recommendation_diagnostics", {})
    write_json(output_dir / "recommendation_diagnostics_v41.json", diagnostics)
    v41_manifest = {
        **manifest,
        "schema_version": "4.1",
        "base_schema_version": "4.0",
        "compatibility_manifest": "v4_manifest.json",
        "hypothesis_count": len(result.surface_hypotheses),
        "display_region_count": len(result.sectors),
        "candidate_order_counts": {
            str(order): sum(item.order == order for item in result.candidates)
            for order in sorted({item.order for item in result.candidates})
        },
        "primary_candidate_ids": list(result.panel.selected_candidate_ids),
        "panel_identifiability": _rows(result.panel),
        "semantic_extensions": [
            "independent_exposure_coordinate_mutation_qualification",
            "footprint_family_inference_hierarchy",
            "global_hypothesis_balanced_candidate_pool",
            "adaptive_budget_v3_profile",
            "overlapping_surface_hypotheses",
            "probe_hypothesis_perturbation_matrix",
            "identifiability_greedy_panel",
            "independent_evidence_and_risk_components",
            "four_site_candidates_exploratory_only",
        ],
        "v41_artifacts": [name for name in [
            *extension_names,
            "surface_graph_v41.json",
            "surface_display_regions_v41.json",
            "partition_diagnostics_v4.json",
            "epitope_hypotheses_v41.json",
            "hypothesis_membership_v41.csv",
            "probe_hypothesis_matrix_v41.csv",
            "panel_identifiability_v41.json",
            "recommendation_diagnostics_v41.json",
        ] if (output_dir / name).is_file()],
    }
    v41_schema_path = Path(__file__).resolve().parents[1] / "schemas" / "joint_v41_manifest.schema.json"
    v41_schema = json.loads(v41_schema_path.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(v41_schema)
    Draft202012Validator(v41_schema).validate(v41_manifest)
    write_json(output_dir / "v4_manifest_v41.json", v41_manifest)
    foldx_status = result.provenance.get("foldx_evidence_status", {})
    candidates_by_id = {item.candidate_id: item for item in result.candidates}
    selected_rows = [
        candidates_by_id[candidate_id]
        for candidate_id in result.panel.selected_candidate_ids
        if candidate_id in candidates_by_id
    ]
    selected_probe_table = [
        "| Candidate | Display anchor | Order | Mutations | SPURS multi ΔΔG | FoldX ΔΔG | Risk components |",
        "|---|---|---:|---|---:|---:|---|",
        *[
            "| " + " | ".join([
                item.candidate_id,
                item.anchor_region_id or item.sector_id,
                str(item.order),
                ", ".join(option.notation for option in item.mutations),
                str(item.stability.spurs_multi_ddg_kcal_mol if item.stability else None),
                str(item.stability.foldx_ddg_kcal_mol if item.stability else None),
                json.dumps(item.risk_components, sort_keys=True),
            ]) + " |"
            for item in selected_rows
        ],
    ]
    hypotheses_with_direct_glycan = sum(
        bool(item.glycan_evidence.get("mutation_prohibition_positions"))
        for item in result.surface_hypotheses
    )
    hypotheses_with_near_glycan = sum(
        bool(item.glycan_evidence.get("near_resolved_glycan_positions"))
        for item in result.surface_hypotheses
    )
    hypotheses_with_pae = sum(bool(item.structural_uncertainty.get("pae_available")) for item in result.surface_hypotheses)
    hypotheses_with_tm_axis = sum(bool(item.membrane_geometry.get("membrane_plane_surrogate_available")) for item in result.surface_hypotheses)
    hypotheses_with_relative_axis = sum(
        item.membrane_geometry.get("orientation_confidence") == "relative_axis_surrogate_only"
        for item in result.surface_hypotheses
    )
    mechanism_policy = "regional_mechanism_pareto_v1"
    mechanism_recommendations = [
        row for row in recommendation_rows if row.get("recommendation_policy") == mechanism_policy
    ]
    recommendation_basis_note = (
        "- Exclusive primary suggestions are selected as regional engineering representatives."
        if mechanism_recommendations else
        "- Exclusive primary suggestions are selected from mutations contained within the display region."
    )
    ranking_note = (
        "- Regional mechanism policy uses bounded regional coverage as context alongside mutation mechanism evidence. Geodesic neighborhoods at 9/12/15 Å, overlap and marginal area are geometry summaries, not measured perturbation or binding evidence. Missing geometry remains unknown."
        if mechanism_recommendations else
        "- Within-sector ranking prioritizes the fraction of the actual display sector covered by geodesic neighborhoods across 9/12/15 Å, then added neighborhood area, overlap and direct mutation-site coverage. This requires eligible surface areas and saved neighborhood positions. Neighborhood coverage is a geometry proxy, not measured perturbation or binding evidence, and cannot guarantee complete biological perturbation; residue-count proxies are not areas. Missing geometry remains unknown."
    )
    mechanism_notes = ([
        "- Regional mechanism recommendations select an engineering representative for testing; they do not predict the maximum antibody-affinity loss.",
        "- Epitope tendency is retained as an annotation/context prior, not an independent mutation-effect axis or antibody-binding probability. Chemical deltas are uncalibrated mechanistic descriptors; missing axes remain unknown.",
        "- Selection sensitivity reports whether predefined variants change the selected representative; it is not a probability or confidence estimate.",
        *(["- Ranking mechanism_channels_v2 uses maximum normalized axis utility within charge/polarity and size/aromaticity groups, with exposure as a separate group and leximin across groups. This is an explicit engineering preference, not an unbiased optimum or affinity model.",
           "- Signed net changes, aromatic substitution counts and WT C-alpha proximity are review-only diagnostics; zero net change does not imply physical cancellation, salt bridges or mutant clashes. REVIEW-band single-mutation losses are audited without relaxing stability gates."]
          if any(row.get("impact_evidence", {}).get("ranking_version") == "mechanism_channels_v2" for row in mechanism_recommendations) else []),
    ] if mechanism_recommendations else [])
    mechanism_columns = ([
        "Policy", "Minimum relative side-chain exposure (%)", "Direct mutated area (Å²)",
        "12 Å local covered area (Å²)", "12 Å local coverage fraction", "12 Å uncovered area / positions",
        "Impact evidence", "Selection sensitivity",
    ] if mechanism_recommendations else [])
    recommendation_headers = [
        "Sector", "Suggestion role and principles", "Candidate", "Gate status", "Mutations", "Sector coverage gain",
        "Selected in primary panel", "Failed gates and evidence", "Geometry proxy (9/12/15 Å)",
        "Regional proxy coverage and uncovered residues", *mechanism_columns,
    ]
    report = [
        "# Joint V4-ML scientific report",
        "",
        f"- Recommendation status: `{result.recommendation_status}`",
        f"- Sectors: {len(result.sectors)}; candidates: {len(result.candidates)}",
        f"- Selected panel: {len(result.panel.selected_candidate_ids)} combinations",
        "- Per-sector suggestions are separate from the formal primary panel; near-threshold suggestions do not count toward formal coverage or the primary construct budget.",
        recommendation_basis_note,
        "- Per-sector coverage_gain uses geodesic_0: mapped surface area of mutated sites within the sector, without neighborhood expansion, divided by the sector surface area; coverage_gain_area_a2 records its area numerator. Each region is evaluated independently, without subtracting mutations displayed in earlier regions.",
        ranking_note,
        *mechanism_notes,
        "",
        "| " + " | ".join(recommendation_headers) + " |",
        "|" + "|".join(["---"] * len(recommendation_headers)) + "|",
        *["| " + " | ".join(str(value).replace("|", "\\|").replace("\n", " ") for value in [
            row["sector_id"], (row.get("recommendation_role", "PRIMARY") + ": " + ", ".join(row.get("principles", []))), row.get("candidate_id") or "—", row["status"],
            ", ".join(row.get("mutation_notations", [])), row.get("coverage_gain", 0),
            row.get("selected_in_panel", False),
            json.dumps({"failed_gates": row.get("failed_gates", []), "gate_details": row.get("gate_details", [])}, ensure_ascii=False),
            json.dumps(row.get("geometry_support", {}), ensure_ascii=False),
            json.dumps(row.get("regional_geometry_support", {}), ensure_ascii=False),
            *([row.get("recommendation_policy", "—"),
               (f"{100 * row['impact_evidence']['exposure_min']:.1f}%"
                if isinstance(row.get("impact_evidence", {}).get("exposure_min"), (int, float)) else "—"),
               row.get("direct_area_a2", "—"), row.get("local_12_covered_area_a2", "—"),
               row.get("local_12_coverage_fraction", "—"),
               json.dumps({"area_a2": row.get("local_12_uncovered_area_a2"),
                           "positions": row.get("local_12_uncovered_positions")}, ensure_ascii=False),
               json.dumps(row.get("impact_evidence", {}), ensure_ascii=False),
               json.dumps(row.get("sensitivity", {}), ensure_ascii=False)] if mechanism_recommendations else []),
        ]) + " |" for row in recommendation_rows],
        "",
        f"- Budget decision: {_rows(result.panel.budget_decision)}",
        f"- Surface complexity: {_rows(result.complexity)}",
        "- Near-distance connected groups are pending-resolution summaries, not pairwise equivalence classes.",
        "- Single-site accessibility and relative solvent accessibility use native total rSASA. For the regional mechanism policy, exposure is the side-chain SASA ratio (numerator / denominator), shown as relative side-chain exposure; it is not the same quantity as total rSASA. Confidence review is distinct from folding damage.",
        f"- Overlapping surface hypotheses: {len(result.surface_hypotheses)}",
        f"- Hypothesis selection method: {result.panel.solver_status}",
        f"- Mean and minimum hypothesis signature separation: {result.panel.mean_pairwise_separation}, {result.panel.minimum_signature_distance}",
        f"- DiscoTope-prior weighted mean separation (reported alongside uniform): {result.panel.prior_weighted_mean_pairwise_separation}",
        f"- Uniform and DiscoTope-prior weighted hypothesis coverage: {result.panel.mean_hypothesis_response}, {result.panel.prior_weighted_hypothesis_response}",
        f"- Marginal information gain by selected construct: {result.panel.marginal_information_gain}",
        f"- Unresolved hypothesis groups: {len(result.panel.unresolved_hypothesis_groups)}",
        *[
            f"  - {', '.join(group)}"
            for group in result.panel.unresolved_hypothesis_groups
        ],
        f"- Uncovered hypotheses: {len(result.panel.uncovered_hypotheses)}",
        f"- Four site exploratory candidates: {sum(item.order >= 4 for item in result.candidates)}; these are excluded from the primary panel",
        f"- FoldX run successful: {foldx_status.get('run_success')}; candidates scored: {foldx_status.get('scored_candidates')}/{foldx_status.get('total_candidates')}; selected panel complete: {foldx_status.get('selected_panel_complete')}",
        f"- Global surface coverage denominator: {result.panel.coverage_denominator_residues} mapped residues",
        f"- Global union coverage by geodesic scenario: {result.panel.achieved_coverage}",
        f"- Hypotheses with direct glycosylation/sequon mutation prohibitions: {hypotheses_with_direct_glycan}",
        f"- Hypotheses near resolved glycans (accessibility uncertainty): {hypotheses_with_near_glycan}",
        f"- Hypotheses with local PAE: {hypotheses_with_pae}",
        f"- Hypotheses with a TM-helix plane surrogate / relative anchor axis: {hypotheses_with_tm_axis} / {hypotheses_with_relative_axis}",
        f"- Epitope-prioritized: {len(result.panel.candidates_by_role.get('EPITOPE_PRIORITIZED', []))}",
        f"- Stability fallback: {len(result.panel.candidates_by_role.get('STABILITY_FALLBACK', []))}",
        "",
        "## Selected probe evidence",
        "",
        *selected_probe_table,
        "",
        "## Recommendation diagnostics",
        "",
        f"- Outcome: {diagnostics.get('outcome', 'UNASSESSED')}",
        f"- Unique combinations / candidate records: {diagnostics.get('unique_combination_count')} / {diagnostics.get('candidate_count')}",
        f"- Formal eligible candidates: {diagnostics.get('formal_eligible_candidate_count')}",
        f"- Formal exclusion counts: {json.dumps(diagnostics.get('formal_exclusion_counts', {}), sort_keys=True)}",
        *[
            f"- {name}: {error.get('code')}: {error.get('message')}"
            for name, errors in diagnostics.get("worker_failures", {}).items()
            for error in errors
        ],
        *[
            f"- {name} was not run: {error.get('message')}"
            for name, errors in diagnostics.get("workers_not_run", {}).items()
            for error in errors
        ],
        "",
        "## Interpretation boundaries",
        "",
        "DiscoTope scores are soft epitope-prior evidence, not assertions of true antibody contacts.",
        "SPURS/FoldX results are predicted stability evidence and do not ensure expression, folding, or cell-surface display.",
        "Preview outputs are not formal recommendations.",
        "Hypothesis responses are surface-geodesic and chemistry proxies, not predicted binding changes.",
        "Hypotheses are candidate antibody-compatible footprints, not experimentally established epitopes.",
        "SASA defines solvent exposure; ProCaliper is an independently audited residue-SASA cross-check in this worker and does not compute finite-probe steric reach.",
        "Direct glycosylation sites and sequons prohibit mutation; nearby resolved glycans are reported as accessibility uncertainty.",
        "Membrane-position distances use a mapped transmembrane helix axis as a plane surrogate; they are not assay-accessibility predictions. A C-terminal anchor axis is relative geometry only; intermembrane distance and Fab reach remain unknown.",
        "Risk components retain stability, glycan, native-interface, disulfide, and structural-uncertainty evidence separately; optimizer penalties are tie-break heuristics, not probabilities.",
        "",
        "## Warnings",
        "",
        *[f"- {item}" for item in result.warnings],
    ]
    physical = result.provenance.get("surface_perturbation")
    if physical and physical.get('panel_information'):
        report.extend(['', 'Geometric information and robustness (not measured binding-response codes):',
                      '```json',json.dumps(physical['panel_information'],indent=2),'```'])
    if physical:
        report.extend(["", "## Physical surface preview", "", f"Status: {physical['status']}",
                       f"Combination target: {physical.get('goal_status','NOT_EVALUATED')}; search: {physical.get('search_status','NOT_EVALUATED')}",
                       "Hard uses resolved local channels on native SES; APBS never changes the hard objective.",
                       "Nested K curves are not independently optimal panels or antibody detection probabilities.",
                       f"Verified constructs: {physical.get('actual_constructs',0)}",
                       "Full raw channels, unknown masks and matched references: surface_perturbation_evidence.zip",
                       "", "```json", json.dumps(physical.get('coverage_curve',[]),indent=2), "```"])
    (output_dir / "v4_scientific_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return {path.name: path.name for path in output_dir.iterdir() if path.is_file()}
def externalize_checkpoint_pair_diagnostics(analysis, output_dir):
    """Keep full proxy diagnostics in their hashed artifact, once per checkpoint."""
    from ..util import sha256_file
    if not any(analysis.get(name, {}).get('selected_panel', {}).get('hierarchy_resolution')
               for name in ('partitioning', 'joint_v4_ml')):
        return
    artifact=Path(output_dir)/'hypothesis_hierarchy_v411.json'
    if not artifact.is_file():raise ValueError('Missing full hierarchy diagnostics artifact')
    saved=json.loads(artifact.read_text(encoding='utf-8'))['resolution_by_level']
    digest=sha256_file(artifact)
    for container in (analysis.get('partitioning',{}),analysis.get('joint_v4_ml',{})):
        panel=container.get('selected_panel',{})
        hierarchy=panel.get('hierarchy_resolution')
        if hierarchy is None:continue
        if hierarchy!=saved:raise ValueError('Checkpoint hierarchy diagnostics identity mismatch')
        projected=dict(hierarchy)
        for level,details in hierarchy.items():
            if not isinstance(details,dict):continue
            projected[level]=dict(details)
            for name in ('near_pairs','pairwise_distances'):
                if name not in details:continue
                projected[level].pop(name)
                projected[level][name+'_reference']={'file':artifact.name,'sha256':digest,
                    'json_path':['resolution_by_level',level,name],'count':len(details[name])}
        panel['hierarchy_resolution']=projected
    analysis.setdefault('files',{})[artifact.name]=artifact.name


def project_display_analysis(
    analysis: dict[str, Any], *, compact_bulk_identity: bool = True
) -> dict[str, Any]:
    """Project saved results for both viewers without mutating scientific evidence."""
    if not analysis.get("joint_v4_ml"):
        return dict(analysis)

    def reference(file, path, value):
        return {"file": file, "json_path": path, "count": len(value)}

    projected = dict(analysis)
    v4 = projected["joint_v4_ml"] = dict(analysis["joint_v4_ml"])
    for name, artifact, id_field in (
        ("candidates", "panel_candidates.json", "candidate_id"),
        ("surface_hypotheses", "surface_hypotheses_v4.json", "hypothesis_id"),
    ):
        if name in v4:
            rows = v4.pop(name)
            ref = reference(artifact, [], rows)
            ref["ids"] = [row.get(id_field) for row in rows
                          if isinstance(row, dict) and row.get(id_field) is not None]
            v4[name + "_reference"] = ref
    panel = v4["selected_panel"] = dict(v4.get("selected_panel", {}))
    if "nearby_hypothesis_pairs" in panel:
        name = "nearby_hypothesis_pairs"
        panel[name + "_reference"] = reference("selected_panel_v4.json", [name], panel.pop(name))
    hierarchy = panel.get("hierarchy_resolution")
    if isinstance(hierarchy, dict):
        hierarchy = panel["hierarchy_resolution"] = dict(hierarchy)
        for level, details in list(hierarchy.items()):
            if not isinstance(details, dict):
                continue
            details = hierarchy[level] = dict(details)
            for name in ("near_pairs", "pairwise_distances"):
                if name in details:
                    details[name + "_reference"] = reference(
                        "hypothesis_hierarchy_v411.json", ["resolution_by_level", level, name], details.pop(name))
    for container_name, field, path in (
        ("partitioning", "selected_panel", ["joint_v4_ml", "selected_panel"]),
        ("design_summary", "sector_recommendations", ["joint_v4_ml", "selected_panel", "sector_recommendations"]),
        ("design_summary", "epitope_mapping_plan", ["epitope_mapping_plan"]),
    ):
        if field in projected.get(container_name, {}):
            container = projected[container_name] = dict(projected[container_name])
            container.pop(field)
            container[field + "_reference"] = {"display_json_path": path, "file": "analysis.json"}

    # Per-file hashes are audit data; keep the aggregate fingerprint and model identity visible.
    bulk_identity_fields = {"source_files_sha256", "model_files_sha256", "git_tracked_modifications"}

    def compact_identity(value, path=()):
        if isinstance(value, list):
            return [compact_identity(item, (*path, index)) for index, item in enumerate(value)]
        if not isinstance(value, dict):
            return value
        compact = {}
        for name, item in value.items():
            if name in bulk_identity_fields and isinstance(item, (dict, list)):
                compact[name + "_reference"] = {
                    "file": "analysis.json", "json_path": [*path, name],
                    "count": len(item),
                }
            else:
                compact[name] = compact_identity(item, (*path, name))
        return compact

    af3 = (projected.get("methods") or {}).get("af3_preprocessing")
    if af3:
        projected["methods"] = dict(projected["methods"])
        compact = projected["methods"]["af3_preprocessing"] = dict(af3)
        artifact = (projected.get("source", {}).get("af3_preprocessing") or {}).get("report", "af3_selection_report.json")
        for name in ("Rcore", "Rdiag", "Ecore", "diagnostic_edges", "pair_rmsd",
                     "per_residue_diagnostics", "domain_pae", "minority_contacts"):
            if name in compact:
                compact[name + "_reference"] = reference(artifact, [name], compact.pop(name))
    return compact_identity(projected) if compact_bulk_identity else projected

