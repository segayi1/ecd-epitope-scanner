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
    _deconvolution_presets, _legacy_designs, _legacy_regions, _sector_recommendation_designs,
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
from ecd_epitope_scanner.foldx_source_evidence import copy_foldx_source_evidence




def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_run", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--repartition-small-protein", action="store_true",
                        help="Repartition a chain <=150 aa at about 8 surface residues per region; reuse only saved candidate scores")
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
    source_names += tuple(name for name in ("v4_manifest.json", "structure_contexts.json", "prediction_context.pdb", "structure_for_energy.pdb",
                                           "source_confidence.json", "source_input_identity.json", "native_input.cif", "native_input.pdb")
                          if (source / name).is_file())
    af3_names = tuple(path.relative_to(source).as_posix()
        for path in sorted((source / "af3_inputs").rglob("*")) if path.is_file())
    af3_names += tuple(path.name for path in sorted(source.glob("af3_preprocessing*")) if path.is_file())
    saved_analysis = json.loads((source / "analysis.json").read_text(encoding="utf-8"))
    af3_names += tuple((saved_analysis.get("source", {}).get("af3_preprocessing") or {}).get("files", []))
    af3_names = tuple(dict.fromkeys(af3_names))
    source_names += af3_names
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
    if args.repartition_small_protein and len(records) > 150:
        parser.error("Small-protein repartitioning requires an analyzed chain <=150 aa")
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
    partition_diagnostics = None
    if args.repartition_small_protein:
        from ecd_epitope_scanner.config import display_sectoring_parameters
        from ecd_epitope_scanner.joint_design_v4.sectors import build_epitope_sectors
        size_keys = {"target_sector_residues", "preferred_minimum_sector_residues",
                     "preferred_maximum_sector_residues", "preferred_minimum_sector_count",
                     "preferred_maximum_sector_count"}
        params = {**params, "sectoring": display_sectoring_parameters(len(records), {
            key: value for key, value in params["sectoring"].items() if key not in size_keys})}
        partition_diagnostics = {}
        sectors = build_epitope_sectors([
            SurfaceEvidence(**{**row, "residue": ResidueKey(**row["residue"])})
            for row in analysis["joint_v4_ml"]["surface_evidence"]], epitope, graph,
            {row.uniprot_pos: row for row in records}, params.get("max_sectors", 7),
            params["sectoring"], diagnostics=partition_diagnostics)
        # Frozen generation anchors and score identities refer to the original partition.
        evidence_source["resource_audits"] = {}
    annotate_spatial_support(candidates, graph, areas)
    single = []
    with (source / "single_mutation_options.csv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            single.append(SingleMutationEvidence(ResidueKey(**ast.literal_eval(row["residue"])), row["mutant_aa"],
                float(row["spurs_ddg_kcal_mol"]) if row["spurs_ddg_kcal_mol"] else None,
                row["stability_band"], row["chemistry_class"], eligible=row["eligible"] == "True"))
    from ecd_epitope_scanner.joint_design_v4.regional_perturbation import (
        POLICY, RANKING_VERSION, annotate_local_interactions, audit_single_thresholds)
    regional = params.get("regional_recommendation_policy") == POLICY
    if regional:
        from ecd_epitope_scanner.joint_design_v4.evidence import build_single_evidence
        rebuilt = build_single_evidence(records, {item.residue.uniprot_pos: item.residue for item in single},
            [{"position": item.residue.uniprot_pos, "mutant_aa": item.mutant_aa,
              "ddg": item.spurs_ddg_kcal_mol} for item in single], params["spurs_single"],
            surface_by_position={row["residue"]["uniprot_pos"]: SurfaceEvidence(**{
                **row, "residue": ResidueKey(**row["residue"])}) for row in analysis["joint_v4_ml"]["surface_evidence"]})
        assert [(item.residue.uniprot_pos, item.mutant_aa, item.eligible) for item in rebuilt] == [
            (item.residue.uniprot_pos, item.mutant_aa, item.eligible) for item in single], "Single eligibility drift"
        record_map = {record.uniprot_pos: record for record in records}
        for item in rebuilt:
            item.impact_evidence = {"policy": POLICY, "sidechain_exposure":
                record_map[item.residue.uniprot_pos].sidechain_exposure,
                "chemical_descriptors": item.chemical_descriptors}
        single = rebuilt
        annotate_local_interactions(candidates, records)
        threshold_audits = audit_single_thresholds(single, sectors, records, graph, params["spurs_single"])
    budget = read("panel_budget_decision_v411.json")
    budget = {key: budget[key] for key in ("minimum_constructs", "target_constructs", "maximum_constructs", "adaptive_budget")}
    panel, selected = select_panel([item for item in candidates if item.formal_eligible], sectors=sectors,
        surface_areas=areas, hypotheses=hypotheses, panel_budget=budget,
        higher_order_max=params["panel_optimizer"].get("higher_order_max_primary", 2),
        risk_penalty_weights=params["panel_optimizer"].get("risk_penalty_weights"))
    panel.sector_recommendations = build_sector_recommendations(candidates, sectors, areas, params, selected,
        graph=graph, records=records, source_provenance=evidence_source, epitope_evidence=epitope)
    prior_recommendations = {} if args.repartition_small_protein else {row["sector_id"]: row for row in
        analysis.get("joint_v4_ml", {}).get("selected_panel", {}).get("sector_recommendations", [])}
    for row in panel.sector_recommendations:
        prior = prior_recommendations.get(row["sector_id"], {})
        if "resource_audit" in prior:
            row["resource_audit"] = dict(prior["resource_audit"])
        if regional:
            row.setdefault("resource_audit", {})["single_threshold_audit"] = threshold_audits[row["sector_id"]]
        if args.repartition_small_protein:
            row.setdefault("resource_audit", {}).update(candidate_pool_scope="saved_source_pool_after_repartition",
                new_candidate_generation=False, new_model_calls=0)
    selected_positions = {m.residue.uniprot_pos for candidate in selected for m in candidate.mutations}
    panel.uncovered_components = [{"sector_id": row["sector_id"], "reason":
        "NO_SELECTED_FEASIBLE_CANDIDATE" if row["formal_eligible"] else "NO_FEASIBLE_CANDIDATE"}
        for row, sector in zip(panel.sector_recommendations, sorted(sectors, key=lambda s: s.sector_id))
        if not selected_positions.intersection(sector.residue_positions)]
    old_protein = analysis["protein"]
    if (source / "uniprot.json").is_file():
        sequence = read("uniprot.json")["sequence"]["value"]
    else:
        sequence = old_protein["ecd_sequence"]
        assert old_protein["ecd_start"] == 1 and len(sequence) == old_protein["sequence_length"], "Missing full sequence evidence"
    protein = ProteinContext(old_protein["accession"], old_protein["protein_name"], old_protein["organism"],
        sequence, old_protein["ecd_start"], old_protein["ecd_end"],
        old_protein["membrane_anchor_position"], old_protein["topology"])
    regions = (_legacy_regions(sectors, records) if args.repartition_small_protein else
               [load(RegionRecord, row) for row in analysis["regions"]])
    primary = _legacy_designs(selected, sectors, protein, None, analysis["joint_v4_ml"]["recommendation_status"])
    controls = _deconvolution_presets(primary, records, protein, None)
    alternatives = _sector_recommendation_designs(panel, candidates, primary, sectors, protein, None)
    designs = primary + controls + alternatives
    if analysis.get("analysis_target") == "binder":
        for design in designs:
            design.rationale = design.rationale.replace("DiscoTope-prioritized, stability-gated V4 candidate.",
                "Binder surface candidate with saved stability evidence.")
    provenance = {"source_run": str(source), "source_files_sha256": hashes,
        "scope": "Saved-score selection replay with geometric support recomputed from frozen graph; no surface reconstruction or workers rerun",
        "new_scientific_calls": 0}
    if regional:
        provenance.update(ranking_version=RANKING_VERSION,
            candidate_population="frozen source pool; no full candidate regeneration or new confirmation scheduling")
    if args.repartition_small_protein:
        provenance.update(scope="Small-protein display repartition and saved-score selection replay; no new candidate generation or model calls",
                          original_sector_count=len(analysis["regions"]), repartitioned_sector_count=len(sectors),
                          candidate_population="frozen source pool; a new region may have no scored exclusive combination")
    v4 = analysis["joint_v4_ml"]
    diagnostics = {**v4["recommendation_diagnostics"], "selected_construct_count": len(selected)}
    result = JointV4Result("4.0", v4["recommendation_status"],
        {name: ToolInventory(**row) for name, row in v4["tool_inventory"].items()},
        [SurfaceEvidence(**{**row, "residue": ResidueKey(**row["residue"])}) for row in v4["surface_evidence"]],
        [EpitopeEvidence(**{**row, "residue": ResidueKey(**row["residue"])}) for row in v4["epitope_evidence"]],
        single, sectors, candidates, panel, surface_hypotheses=hypotheses,
        legacy_regions=regions, legacy_designs=designs, warnings=v4["warnings"] + [provenance["scope"]],
        provenance={**v4_manifest.get("provenance", {}), "selection_replay": provenance,
                    "recommendation_diagnostics": diagnostics},
        raw_hypotheses=[load(V4SurfaceHypothesis, row) for row in read("raw_hypotheses_v411.json")],
        footprint_families=[load(FootprintFamily, row) for row in read("footprint_families_v411.json")],
        complexity=SurfaceComplexitySummary(**read("surface_complexity_v411.json")))
    if version_identification:
        for name in version_identification["tools"]:
            for field in ("version", "git_revision", "model_revision"):
                setattr(result.tool_inventory[name], field, models[name]["tool_inventory"][field])
    output.mkdir(parents=True)
    portable_evidence = copy_foldx_source_evidence(source, output, v4_manifest)
    if portable_evidence.get("status") == "UNKNOWN":
        provenance["portable_source_evidence"] = portable_evidence
        portable_evidence = None
    else:
        provenance["portable_source_evidence"] = {key: portable_evidence[key] for key in
            ("manifest_path", "manifest_sha256", "artifact_paths")}
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
    for name in ("source_confidence.json", "source_input_identity.json", "native_input.cif", "native_input.pdb",
                 "structure.cif", "ecd_structure.cif", "prediction_context.pdb", "structure_for_energy.pdb",
                 "uniprot.json", "residue_mapping.csv", "surface_graph_v4.json", "surface_graph_v41.json",
                 "surface_display_regions_v41.json", "partition_diagnostics_v4.json", "native_surface_rows.json",
                 "raw_surface_footprints_v411.json", "spurs_single_matrix.csv.gz"):
        if (source / name).is_file():
            shutil.copy2(source / name, output / name)
    for name in af3_names:
        original, destination = source / name, output / name
        if not original.resolve().is_relative_to(source):
            raise ValueError("AF3 replay artifact escapes the frozen source run")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, destination)
    if args.repartition_small_protein:
        from ecd_epitope_scanner.joint_design_v4.surface_snapshot import write_surface_graph_snapshot
        graph_provenance = {**graph_data["graph_provenance"], "display_partition_policy": partition_diagnostics}
        graph_provenance.pop("snapshot_sha256", None)
        write_json(output / "partition_diagnostics_v4.json", partition_diagnostics)
        for name, ownership in (("surface_graph_v4.json", True), ("surface_graph_v41.json", False)):
            snapshot_hash = write_surface_graph_snapshot(output / name, graph, areas, sectors,
                mapping_path=output / "residue_mapping.csv", structure_path=output / "structure_for_energy.pdb",
                provenance=graph_provenance, include_display_ownership=ownership)
            if ownership:
                graph_provenance.update(snapshot_file=name, snapshot_sha256=snapshot_hash)
        display = read("surface_display_regions_v41.json")
        display.update(regions=to_dict(sectors), partition_policy=partition_diagnostics)
        write_json(output / "surface_display_regions_v41.json", display)
        result.provenance["surface_graph"] = graph_provenance
        replay_manifest["provenance"] = result.provenance
        write_json(output / "v4_manifest.json", replay_manifest)
        v41 = json.loads((output / "v4_manifest_v41.json").read_text(encoding="utf-8"))
        v41["provenance"] = result.provenance
        write_json(output / "v4_manifest_v41.json", v41)
        analysis["regions"] = to_dict(regions)
        analysis["residues"] = to_dict(records)
        analysis["joint_v4_ml"]["sectors"] = to_dict(sectors)
        analysis["summary"]["regions"] = len(sectors)
        analysis["partitioning"].update(sector_count=len(sectors), actual_regions=len(sectors),
            display_partition_policy=partition_diagnostics)
    analysis["request"]["output_dir"] = str(output)
    if args.repartition_small_protein:
        analysis["request"]["target_sector_residues"] = 8
        metadata = analysis["request"].setdefault("metadata", {})
        metadata.setdefault("joint_v4_ml", {})["sectoring"] = updated_params["sectoring"]
    analysis["parameters"]["joint_v4_ml"] = updated_params
    analysis["selection_replay"] = provenance
    analysis["designs"] = to_dict(designs)
    analysis["partitioning"].update(primary_construct_count=len(selected), selected_panel=to_dict(panel),
                                  stage1_selection_basis=panel.solver_status.lower(),
                                  stage1_selected_region_ids=sorted({item.region_id for item in primary if item.region_id}))
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
    plan["sectors"] = to_dict(sectors)
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
    if portable_evidence:
        analysis["files"].update({name: name for name in portable_evidence["artifact_paths"]})
        run_manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
        run_manifest.setdefault("files", {}).update({name: "" for name in portable_evidence["artifact_paths"]})
        write_json(output / "run_manifest.json", run_manifest)
    write_json(output / "analysis.json", analysis)
    finalize_scientific_bundle(output)
    assert hashes == {name: sha256_file(source / name) for name in source_names}
    assert all(candidate.formal_eligible for candidate in selected)
    assert len(panel.sector_recommendations) == len(sectors)
    assert json.loads((output / "v4_manifest.json").read_text(encoding="utf-8"))["worker_runs"] == v4_manifest.get("worker_runs", {})
    score_identity = lambda rows: [(row["candidate_id"],
        [(m["residue"]["uniprot_pos"], m["mutant_aa"]) for m in row["mutations"]],
        (row.get("stability") or {}).get("spurs_multi_ddg_kcal_mol"),
        (row.get("stability") or {}).get("foldx_ddg_kcal_mol")) for row in rows]
    assert score_identity(to_dict(candidates)) == score_identity(read("panel_candidates.json")), "Frozen score or mutation identity drift"
    assert all(all(row.get("resource_audit", {}).get(key) == value for key, value in
        prior_recommendations.get(row["sector_id"], {}).get("resource_audit", {}).items())
        for row in panel.sector_recommendations)
    print(json.dumps({"output_dir": str(output), "selected_primary_count": len(selected),
        "sector_statuses": {row["sector_id"]: row["status"] for row in panel.sector_recommendations},
        "direct_area_coverage": panel.achieved_area_coverage.get("geodesic_0"),
        "new_scientific_calls": 0, "historical_hashes_unchanged": True}, ensure_ascii=False))


if __name__ == "__main__":
    main()
