"""Replay V4 selection and region suggestions from saved scores, without workers."""
from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import shutil
import sys
from dataclasses import fields
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ecd_epitope_scanner.exporting import finalize_scientific_bundle, write_scientific_outputs
from ecd_epitope_scanner.joint_design_v4.engine import (
    _deconvolution_presets, _legacy_designs, _sector_recommendation_designs,
)
from ecd_epitope_scanner.joint_design_v4.models import (
    CombinationStabilityEvidence, EpitopeEvidence, FootprintFamily, InferenceHypothesis,
    JointV4Result, ResidueKey, SingleMutationEvidence, SurfaceComplexitySummary,
    SurfaceEvidence, ToolInventory, V4Candidate, V4MutationOption, V4Sector, V4SurfaceHypothesis,
)
from ecd_epitope_scanner.joint_design_v4.optimizer import select_panel
from ecd_epitope_scanner.joint_design_v4.reporting import write_v4_outputs
from ecd_epitope_scanner.joint_design_v4.sector_recommendations import build_sector_recommendations, read_structure_evidence
from ecd_epitope_scanner.joint_design_v4.spatial_support import annotate_spatial_support
from ecd_epitope_scanner.models import ProteinContext, RegionRecord, ResidueRecord, to_dict
from ecd_epitope_scanner.util import sha256_file, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_run", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--version-identification", type=Path,
                        help="Verified current installation audit bound to this frozen scoring manifest")
    args = parser.parse_args()
    source, output = args.source_run.resolve(), args.output_dir.resolve()
    if source == output or source in output.parents or output in source.parents or output.exists():
        parser.error("Use a new output directory separate from historical runs")

    def read(name):
        return json.loads((source / name).read_text(encoding="utf-8"))

    def load(cls, row):
        data = {field.name: row[field.name] for field in fields(cls) if field.name in row}
        for name in ("residue_weights", "membership_fraction", "membership"):
            if name in data:
                data[name] = {int(position): value for position, value in data[name].items()}
        return cls(**data)

    source_names = ("analysis.json", "run_manifest.json", "panel_candidates.json", "probe_sectors_v4.json",
                    "surface_graph_v41.json", "inference_hypotheses_v411.json", "structure.cif",
                    "residue_mapping.csv", "single_mutation_options.csv")
    source_names += tuple(name for name in ("v4_manifest.json", "structure_contexts.json", "prediction_context.pdb", "structure_for_energy.pdb")
                          if (source / name).is_file())
    hashes = {name: sha256_file(source / name) for name in source_names}
    analysis = read("analysis.json")
    records = [load(ResidueRecord, row) for row in analysis["residues"]]
    epitope = [EpitopeEvidence(**{**row, "residue": ResidueKey(**row["residue"])})
               for row in analysis.get("joint_v4_ml", {}).get("epitope_evidence", [])]
    manifest = read("run_manifest.json")
    v4_manifest = read("v4_manifest.json") if (source / "v4_manifest.json").is_file() else {}
    models = {name: {"provenance": dict(row.get("provenance", {})),
                     "tool_inventory": dict(v4_manifest.get("tools", {}).get(name,
                         v4_manifest.get("tools", {}).get(name.split("_")[0], {}))),
                     "status": row.get("status")}
              for name, row in v4_manifest.get("worker_runs", {}).items()}
    version_identification = None
    if args.version_identification:
        version_identification = json.loads(args.version_identification.read_text(encoding="utf-8"))
        if version_identification.get("source_manifest_sha256") != sha256_file(source / "v4_manifest.json"):
            parser.error("Version audit does not match the frozen scoring manifest")
        for name, identity in version_identification["tools"].items():
            frozen = identity["frozen_run"]
            original = models[name]["provenance"]
            if not frozen["all_current_weight_files_match"] or frozen["weights_sha256"] != original.get("weights_sha256"):
                parser.error(f"Version audit weights do not match {name}")
            supplement = {"basis": "current_environment_metadata_supplement",
                          "original_run_version_recorded": bool(original.get("version")),
                          "weights_match_frozen_run": True,
                          "historical_limitation": version_identification["historical_limitation"]}
            model_revision = ("DiscoTope-" + identity["model_release_version"]
                              if name == "discotope" and identity.get("model_release_version") else None)
            original.update(version=identity["package_version"], package_name=identity["package_name"],
                            git_revision=identity["git_revision_current_checkout"],
                            git_tracked_modifications=identity["git_tracked_modifications"],
                            model_revision=model_revision, version_identification=supplement)
            models[name]["tool_inventory"].update(version=identity["package_version"],
                git_revision=identity["git_revision_current_checkout"], model_revision=model_revision)
    evidence_source = {"models": models, "source_files_sha256": hashes,
                       "normalized_structure_sha256": v4_manifest.get("normalized_structure_sha256"),
                       "source_run": str(source), "structure_evidence": read_structure_evidence(source / "structure.cif",
                            predicted=bool(analysis.get("source", {}).get("explicit_af3") or
                                           "alphafold" in str(analysis.get("source", {}).get("kind", "")).lower()))}
    evidence_source["resource_audits"] = {row["sector_id"]: row["resource_audit"] for row in
        analysis.get("joint_v4_ml", {}).get("selected_panel", {}).get("sector_recommendations", []) if "resource_audit" in row}
    params = manifest["parameters"]["joint_v4_ml"]
    if analysis.get("design_summary", {}).get("surface_perturbation"):
        parser.error("Physical panels require their original physical replay workflow")
    sectors = [V4Sector(**row) for row in read("probe_sectors_v4.json")]
    candidates = []
    for row in read("panel_candidates.json"):
        row["mutations"] = [V4MutationOption(**{**item, "residue": ResidueKey(**item["residue"])}) for item in row["mutations"]]
        row["stability"] = CombinationStabilityEvidence(**row["stability"]) if row["stability"] else None
        if row["stability"]:
            row["stability"].provenance = {**row["stability"].provenance, "models": models,
                "normalized_structure_sha256": evidence_source["normalized_structure_sha256"]}
        row["role"] = "UNASSIGNED"
        missing_models = []
        for field, reason in (("spurs_multi_ddg_kcal_mol", "COMBINATION_SPURS_MISSING"),
                              ("foldx_ddg_kcal_mol", "COMBINATION_FOLDX_MISSING")):
            value = getattr(row["stability"], field, None)
            if value is None or not math.isfinite(float(value)):
                missing_models.append(reason)
        if missing_models:
            row["formal_eligible"] = False
            row["formal_exclusion_reasons"] = sorted(set(row["formal_exclusion_reasons"] + missing_models))
        candidates.append(V4Candidate(**row))
    hypotheses = [load(InferenceHypothesis, row) for row in read("inference_hypotheses_v411.json")]
    areas = {node["uniprot_pos"]: node["surface_area_a2"] for node in read("surface_graph_v41.json")["nodes"]}
    graph_data = read("surface_graph_v41.json")
    graph = {position: {} for position in areas}
    for edge in graph_data["edges"]:
        graph[edge["source"]][edge["target"]] = edge["distance_a"]
        graph[edge["target"]][edge["source"]] = edge["distance_a"]
    annotate_spatial_support(candidates, graph, areas)
    budget = read("panel_budget_decision_v411.json")
    budget = {key: budget[key] for key in ("minimum_constructs", "target_constructs", "maximum_constructs", "adaptive_budget")}
    panel, selected = select_panel([item for item in candidates if item.formal_eligible], sectors=sectors,
        surface_areas=areas, hypotheses=hypotheses, panel_budget=budget,
        higher_order_max=params["panel_optimizer"].get("higher_order_max_primary", 2),
        risk_penalty_weights=params["panel_optimizer"].get("risk_penalty_weights"))
    panel.sector_recommendations = build_sector_recommendations(candidates, sectors, areas, params, selected,
        graph=graph, records=records, source_provenance=evidence_source, epitope_evidence=epitope)
    prior_recommendations = {row["sector_id"]: row for row in
        analysis.get("joint_v4_ml", {}).get("selected_panel", {}).get("sector_recommendations", [])}
    for row in panel.sector_recommendations:
        prior = prior_recommendations.get(row["sector_id"], {})
        if "resource_audit" in prior:
            row["resource_audit"] = prior["resource_audit"]
    selected_positions = {m.residue.uniprot_pos for candidate in selected for m in candidate.mutations}
    panel.uncovered_components = [{"sector_id": row["sector_id"], "reason":
        "NO_SELECTED_FEASIBLE_CANDIDATE" if row["formal_eligible"] else "NO_FEASIBLE_CANDIDATE"}
        for row, sector in zip(panel.sector_recommendations, sorted(sectors, key=lambda s: s.sector_id))
        if not selected_positions.intersection(sector.residue_positions)]
    old_protein = analysis["protein"]
    entry = read("uniprot.json")
    protein = ProteinContext(old_protein["accession"], old_protein["protein_name"], old_protein["organism"],
        entry["sequence"]["value"], old_protein["ecd_start"], old_protein["ecd_end"],
        old_protein["membrane_anchor_position"], old_protein["topology"])
    regions = [load(RegionRecord, row) for row in analysis["regions"]]
    primary = _legacy_designs(selected, sectors, protein, None, analysis["joint_v4_ml"]["recommendation_status"])
    controls = _deconvolution_presets(primary, records, protein, None)
    alternatives = _sector_recommendation_designs(panel, candidates, primary, sectors, protein, None)
    designs = primary + controls + alternatives
    single = []
    with (source / "single_mutation_options.csv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            single.append(SingleMutationEvidence(ResidueKey(**ast.literal_eval(row["residue"])), row["mutant_aa"],
                float(row["spurs_ddg_kcal_mol"]) if row["spurs_ddg_kcal_mol"] else None,
                row["stability_band"], row["chemistry_class"], eligible=row["eligible"] == "True"))
    provenance = {"source_run": str(source), "source_files_sha256": hashes,
        "scope": "Saved-score selection replay with geometric support recomputed from frozen graph; no surface reconstruction or workers rerun",
        "new_scientific_calls": 0}
    v4 = analysis["joint_v4_ml"]
    diagnostics = {**v4["recommendation_diagnostics"], "selected_construct_count": len(selected)}
    result = JointV4Result("4.0", v4["recommendation_status"],
        {name: ToolInventory(**row) for name, row in v4["tool_inventory"].items()},
        [SurfaceEvidence(**{**row, "residue": ResidueKey(**row["residue"])}) for row in v4["surface_evidence"]],
        [EpitopeEvidence(**{**row, "residue": ResidueKey(**row["residue"])}) for row in v4["epitope_evidence"]],
        single, sectors, candidates, panel, surface_hypotheses=hypotheses,
        legacy_regions=regions, legacy_designs=designs, warnings=v4["warnings"] + [provenance["scope"]],
        provenance={"selection_replay": provenance, "recommendation_diagnostics": diagnostics},
        raw_hypotheses=[load(V4SurfaceHypothesis, row) for row in read("raw_hypotheses_v411.json")],
        footprint_families=[load(FootprintFamily, row) for row in read("footprint_families_v411.json")],
        complexity=SurfaceComplexitySummary(**read("surface_complexity_v411.json")))
    if version_identification:
        for name in version_identification["tools"]:
            for field in ("version", "git_revision", "model_revision"):
                setattr(result.tool_inventory[name], field, models[name]["tool_inventory"][field])
    output.mkdir(parents=True)
    updated_params = {**params, "panel_optimizer": {**params["panel_optimizer"],
        "selection_priority": panel.budget_decision.selection_priority, "coverage_basis": panel.budget_decision.coverage_basis}}
    write_v4_outputs(output, result, parameters=updated_params, contexts=read("structure_contexts.json"), source_provenance=result.provenance)
    replay_manifest = json.loads((output / "v4_manifest.json").read_text(encoding="utf-8"))
    replay_manifest.update(worker_runs=v4_manifest.get("worker_runs", {}),
        tool_commands=v4_manifest.get("tool_commands", {}),
        normalized_structure_sha256=v4_manifest.get("normalized_structure_sha256"),
        worker_runs_origin="saved_source_scores_no_new_worker_calls")
    if version_identification:
        replay_manifest["version_identification"] = version_identification
        write_json(output / "model_version_identification.json", version_identification)
    write_json(output / "v4_manifest.json", replay_manifest)
    for name in ("structure.cif", "ecd_structure.cif", "prediction_context.pdb", "structure_for_energy.pdb",
                 "uniprot.json", "residue_mapping.csv", "surface_graph_v4.json", "surface_graph_v41.json",
                 "surface_display_regions_v41.json", "partition_diagnostics_v4.json", "native_surface_rows.json",
                 "raw_surface_footprints_v411.json", "single_mutation_options.csv", "spurs_single_matrix.csv.gz"):
        if (source / name).is_file():
            shutil.copy2(source / name, output / name)
    analysis["request"]["output_dir"] = str(output)
    analysis["parameters"]["joint_v4_ml"] = updated_params
    analysis["selection_replay"] = provenance
    analysis["designs"] = to_dict(designs)
    analysis["partitioning"].update(primary_construct_count=len(selected), selected_panel=to_dict(panel),
                                  stage1_selection_basis=panel.solver_status.lower())
    analysis["joint_v4_ml"].update(candidates=to_dict(candidates), selected_panel=to_dict(panel), recommendation_diagnostics=diagnostics,
                                 tool_inventory=to_dict(result.tool_inventory))
    if version_identification:
        analysis["model_version_identification"] = version_identification
    summary = analysis["design_summary"]
    summary.update(selected_candidate_ids=panel.selected_candidate_ids, panel_budget_decision=to_dict(panel.budget_decision),
        recommendation_diagnostics=diagnostics, stage2_preset_count=len(controls),
        sector_recommendations=panel.sector_recommendations, sector_recommendation_count=len(sectors),
        sector_review_count=sum(row["status"] == "NEAR_THRESHOLD" for row in panel.sector_recommendations),
        unscored_saturation_substitution_count=sum(len(options) for design in primary for options in design.saturation_scan.values()),
        recommended_by_class={key: len(value) for key, value in panel.candidates_by_role.items() if key != "PRIMARY"})
    plan = summary["epitope_mapping_plan"]
    plan["stage1_probes"] = [{"probe_id": design.id, "patch_id": design.region_id, "candidate_id": candidate.candidate_id,
        "mutations": [m.notation for m in design.mutations], "targeted_hypotheses": candidate.targeted_hypotheses,
        "hypothesis_response_proxy": candidate.hypothesis_response, "status": design.status} for design, candidate in zip(primary, selected)]
    plan["coverage_matrix"] = [{"probe_id": design.id, "patch_id": design.region_id, "scenario": scenario,
        "predicted_coverage": coverage, "evidence_type": "GEOMETRIC_MODEL_ONLY"}
        for design, candidate in zip(primary, selected) for scenario, coverage in candidate.coverage_by_scenario.items()]
    plan["stage2_presets"] = [{"design_id": item.id, "parent_design_id": item.parent_design_id,
        "role": item.design_role, "mutations": [choice.notation for choice in item.mutations],
        "conditional_on_parent_probe_hit": True} for item in controls]
    plan["blind_spots"] = panel.uncovered_components
    plan["panel_identifiability"].update(uniform_hypothesis_response=panel.mean_hypothesis_response,
        prior_weighted_hypothesis_response=panel.prior_weighted_hypothesis_response,
        marginal_information_gain=panel.marginal_information_gain, mean_pairwise_separation=panel.mean_pairwise_separation,
        minimum_signature_distance=panel.minimum_signature_distance,
        prior_weighted_mean_pairwise_separation=panel.prior_weighted_mean_pairwise_separation,
        unresolved_hypothesis_groups=panel.unresolved_hypothesis_groups, uncovered_hypotheses=panel.uncovered_hypotheses)
    analysis["epitope_mapping_plan"] = plan
    analysis["summary"]["recommended_constructs"] = len(primary)
    analysis["files"] = {path.name: path.name for path in output.iterdir() if path.is_file()}
    analysis["files"].update(write_scientific_outputs(output, analysis, protein, records, regions, designs, defer_archive=True))
    write_json(output / "analysis.json", analysis)
    finalize_scientific_bundle(output)
    assert hashes == {name: sha256_file(source / name) for name in source_names}
    assert all(candidate.formal_eligible for candidate in selected)
    assert len(panel.sector_recommendations) == len(sectors)
    assert json.loads((output / "v4_manifest.json").read_text(encoding="utf-8"))["worker_runs"] == v4_manifest.get("worker_runs", {})
    assert all(row.get("resource_audit") == prior_recommendations.get(row["sector_id"], {}).get("resource_audit")
               for row in panel.sector_recommendations)
    print(json.dumps({"output_dir": str(output), "selected_primary_count": len(selected),
        "sector_statuses": {row["sector_id"]: row["status"] for row in panel.sector_recommendations},
        "direct_area_coverage": panel.achieved_area_coverage.get("geodesic_0"),
        "new_scientific_calls": 0, "historical_hashes_unchanged": True}, ensure_ascii=False))


if __name__ == "__main__":
    main()
