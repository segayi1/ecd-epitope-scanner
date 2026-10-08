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
    for directory in ("joint_design_v4", "joint_design", "integrations", "surface_patch", "surface_perturbation"):
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
        },
    )
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
                "anchor_region_id": item.anchor_region_id or item.sector_id,
                "order": item.order,
                "mutations": ";".join(option.notation for option in item.mutations),
                "targeted_hypotheses": ";".join(item.targeted_hypotheses),
                "hypothesis_perturbation": json.dumps(item.hypothesis_perturbation or item.hypothesis_response, sort_keys=True),
                "risk_components": json.dumps(item.risk_components, sort_keys=True),
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
    report = [
        "# Joint V4-ML scientific report",
        "",
        f"- Recommendation status: `{result.recommendation_status}`",
        f"- Sectors: {len(result.sectors)}; candidates: {len(result.candidates)}",
        f"- Selected panel: {len(result.panel.selected_candidate_ids)} combinations",
        f"- Budget decision: {_rows(result.panel.budget_decision)}",
        f"- Surface complexity: {_rows(result.complexity)}",
        "- Near-distance connected groups are pending-resolution summaries, not pairwise equivalence classes.",
        "- Exposure uses native total rSASA. Side-chain SASA is descriptive; confidence review is distinct from folding damage.",
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

