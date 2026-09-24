from __future__ import annotations

from pathlib import Path
from typing import Any

from .. import config
from ..models import DesignRecord, MutationChoice, ProteinContext, RegionRecord, ResidueRecord
from ..sequences import apply_mutations, coding_fragment
from .candidate_search import generate_probe_candidates
from .coverage import attach_coverage
from .microfootprints import build_microfootprints
from .models import JointV3Result, ProbeCandidate
from .optimizer import select_probe_panel
from .reporting import write_joint_v3_outputs
from .sectors import build_probe_sectors
from .stability import annotate_stability_evidence
from .structure_audit import audit_structure
from .surface_units import build_surface_units, unit_records
from .topology import build_topology_context


def _color(index: int) -> str:
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#6F4E7C", "#666666"]
    return palette[(index - 1) % len(palette)]


def _legacy_regions(sectors: list[Any], units: list[Any], records: list[ResidueRecord]) -> list[RegionRecord]:
    by_position = {item.uniprot_pos: item for item in records}
    by_unit = {item.position: item for item in units}
    output: list[RegionRecord] = []
    for index, sector in enumerate(sectors, start=1):
        present = [by_position[position] for position in sector.residue_ids if position in by_position]
        mutable = [position for position in sector.residue_ids if position in by_unit and by_unit[position].designable]
        domain = present[0].domain if present else "ECD"
        secondary = {item.secondary_structure for item in present}
        character = "beta-sheet face" if "E" in secondary else "helical face" if "H" in secondary else "mixed surface"
        output.append(
            RegionRecord(
                id=index, name=f"{domain} / probe sector {index}", domain=domain,
                membrane_position="topology-only; membrane plane unknown",
                structural_character=character, residues=mutable,
                surface_area=round(sector.total_area_a2, 2), geodesic_diameter=0.0,
                interface_fraction=round(sum(item.interface for item in present) / max(1, len(present)), 3),
                pae_median=None, auto_design_allowed=bool(mutable),
                rationale=("Joint V3 connected coarse localization sector; this is a ProbeSector, "
                           "not a measured antibody footprint."), color=_color(index),
                domain_id=domain, epitope_prior=None,
                prior_components={"designable_area_a2": sector.designable_area_a2, "total_area_a2": sector.total_area_a2},
                tier="MEDIUM", assay_accessible=True, predicted_accessibility="unknown",
                required_probe_coverage=1, partition_mode="joint_v3", partition_basis="surface_residue_geodesic_graph",
                envelope_residues=list(sector.residue_ids), mutable_residues=mutable,
                protected_surface_residues=[position for position in sector.residue_ids if position not in mutable],
                footprint_residues=mutable, candidate_footprints=[], surface_connectivity="connected",
                quality_flags=list(sector.quality_flags), legacy_projection=True,
            )
        )
    return output


def _legacy_designs(
    selected: list[ProbeCandidate],
    protein: ProteinContext,
    original_cds: str | None,
    records: list[ResidueRecord],
    sector_region_ids: dict[str, int] | None = None,
) -> list[DesignRecord]:
    del records
    designs: list[DesignRecord] = []
    for index, candidate in enumerate(selected, start=1):
        has_combination_energy = (
            candidate.evidence_level == "EXTERNAL_COMBINATION"
            and candidate.stability_value is not None
        )
        choices = [
            MutationChoice(
                position=item.position, wt=item.wild_type, mutant=item.mutant_aa,
                stability_score=round(1.0 - item.rule_risk, 4),
                perturbation_score=round(item.chemical_perturbation_strength, 4),
                stability_percentile=round(1.0 - item.rule_risk, 4), allowed=True,
                reasons=(
                    ["Combination-level FoldX evidence is recorded on the design; no single-site ΔΔG was assigned."]
                    if has_combination_energy else
                    ["dimensionless rule proxy; not an energy value"]
                ),
                stability_class="proxy", menu_class=item.perturbation_class,
                disruption_score=item.chemical_perturbation_strength,
            )
            for item in candidate.mutations
        ]
        mutant_sequence = apply_mutations(protein.ecd_sequence, choices, protein.ecd_start)
        mutant_cds, cds_source = coding_fragment(protein, mutant_sequence, original_cds)
        designs.append(
            DesignRecord(
                id=f"V3_{index:02d}", name=f"Joint V3 probe {index}",
                region_id=(sector_region_ids or {}).get(candidate.sector_id),
                mutations=choices, status="RECOMMENDED_MULTI_SITE_PROBE",
                stability_engine=(
                    "FoldX 5.1 (combination-level)"
                    if has_combination_energy else
                    "builtin rule proxy (dimensionless)"
                ),
                combined_stability_score=(
                    round(float(candidate.stability_value), 4)
                    if has_combination_energy else None
                ),
                perturbation_score=round(sum(item.perturbation_score for item in choices), 4),
                rationale=(f"Selected by deterministic joint panel optimization for sector {candidate.sector_id}; "
                           "no anchor mutation was forced."),
                warnings=list(candidate.warnings) + [f"Coding fragment source: {cds_source}."],
                mutant_ecd_sequence=mutant_sequence, mutant_coding_sequence=mutant_cds,
                score_origin="foldx" if has_combination_energy else "proxy",
                foldx_diagnostics=dict(candidate.stability_provenance),
                design_role="STAGE1_PROBE",
                disruption_score=round(sum(item.disruption_score for item in choices), 4),
                anchor_notation=None,
            )
        )
    return designs


def run_joint_v3(
    records: list[ResidueRecord],
    protein: ProteinContext,
    output_dir: Path,
    original_cds: str | None = None,
    chain_id: str | None = None,
    tm_anchor_xyz: tuple[float, float, float] | None = None,
    tm_axis_xyz: tuple[float, float, float] | None = None,
    parameters: dict[str, Any] | None = None,
    surface_geometry: Any | None = None,
    structure: Any | None = None,
    source_provenance: dict[str, Any] | None = None,
    structure_pdb: Path | None = None,
    engine_inventory: Any | None = None,
) -> JointV3Result:
    params = parameters or config.JOINT_V3_PARAMETERS
    selected_chain = chain_id or (records[0].chain_id if records else "")
    audit = audit_structure(records, protein, selected_chain, structure=structure)
    topology = build_topology_context(records, protein, tm_anchor_xyz, tm_axis_xyz)
    graph_params = params.get("surface_graph", {})
    graph_provenance: dict[str, Any] = {}
    units, graph, areas = build_surface_units(
        records,
        rsasa_threshold=float(params.get("surface_rsasa_threshold", config.RSASA_THRESHOLD)),
        neighbor_cutoff_a=float(graph_params.get("neighbor_cutoff_a", 10.0)),
        max_neighbors=int(graph_params.get("max_neighbors", 8)),
        surface_geometry=surface_geometry,
        graph_provenance=graph_provenance,
    )
    target_areas = [float(item) for item in params.get("footprint_target_areas_a2", [600, 800, 1000])]
    budget = int(params.get("max_probes", 6))
    footprints = build_microfootprints(
        units,
        graph,
        areas,
        target_areas,
        max_centers=int(params.get("max_footprint_centers", max(24, budget * 6))),
    )
    sectors = build_probe_sectors(
        units,
        footprints,
        max_sectors=max(1, budget * 3),
        graph=graph,
        maximum_geodesic_diameter_a=float(params.get("sector_maximum_geodesic_diameter_a", 35.0)),
        maximum_euclidean_diameter_a=float(
            params.get("sector_maximum_euclidean_diameter_a", 30.0)
        ),
    )
    record_map = unit_records(units, records)
    search = params.get("candidate_search", {})
    candidates, exclusions = generate_probe_candidates(
        sectors, list(record_map.values()), protein, graph,
        {unit.position for unit in units if unit.designable},
        site_pool_cap=int(search.get("site_pool_cap_per_sector", 12)),
        options_per_site=int(search.get("mutation_options_cap_per_site", 3)),
        beam_width=int(search.get("beam_width", 200)),
        shortlist_per_sector=int(search.get("shortlist_per_sector", 20)),
        min_mutations=int(params.get("mutations_per_probe", [2, 4])[0]),
        max_mutations=int(params.get("mutations_per_probe", [2, 4])[1]),
        maximum_geodesic_span_a=float(search.get("maximum_geodesic_span_a", 24.0)),
    )
    candidates = attach_coverage(
        candidates,
        units,
        graph,
        footprints,
        distance_scenarios_a={
            str(key): float(value)
            for key, value in params.get("coverage_distances_a", {
                "geodesic_tight": 9.0,
                "geodesic_medium": 12.0,
                "geodesic_wide": 15.0,
            }).items()
        },
    )
    solver_params = params.get("solver", {})
    coverage_scenarios = list(params.get("coverage_scenarios", []))
    proxy_panel, proxy_selected = select_probe_panel(
        candidates,
        units,
        budget,
        float(solver_params.get("time_limit_seconds", 60)),
        float(solver_params.get("mip_relative_gap", 0.01)),
        coverage_scenarios=coverage_scenarios,
    )
    priority_candidate_ids = [item.candidate_id for item in proxy_selected]
    for alternative in proxy_panel.alternative_panels:
        priority_candidate_ids.extend(alternative.get("candidate_ids", []))
    priority_candidate_ids = list(dict.fromkeys(priority_candidate_ids))
    stability_summary = annotate_stability_evidence(
        candidates,
        inventory=engine_inventory,
        structure_pdb=structure_pdb,
        workdir=output_dir,
        parameters=params.get("stability", {}),
        priority_candidate_ids=priority_candidate_ids,
    )
    # Stability is a secondary ranking signal after geometric coverage. Keep
    # deterministic ordering and put measured combination energies before
    # failed/proxy candidates when coverage ties.
    candidates.sort(key=lambda item: (
        item.stability_status not in {"PASS", "ABOVE_GATE"},
        item.stability_value is None,
        item.stability_value if item.stability_value is not None else float("inf"),
        item.structural_risk,
        item.candidate_id,
    ))
    selection_candidates = candidates
    stability_params = params.get("stability", {})
    require_pass = bool(stability_params.get("require_pass_for_final_panel", True))
    if stability_summary.get("scored_count") and require_pass:
        selection_candidates = [
            candidate for candidate in candidates
            if candidate.stability_value is not None
            and candidate.stability_status == "PASS"
        ]
    elif stability_summary.get("scored_count"):
        selection_candidates = [
            candidate for candidate in candidates
            if candidate.stability_value is not None
            and candidate.stability_status in {"PASS", "ABOVE_GATE"}
        ]
    panel, selected = select_probe_panel(
        selection_candidates,
        units,
        budget,
        float(solver_params.get("time_limit_seconds", 60)),
        float(solver_params.get("mip_relative_gap", 0.01)),
        coverage_scenarios=coverage_scenarios,
    )
    stability_summary.update({
        "proxy_preselection_candidate_ids": [item.candidate_id for item in proxy_selected],
        "final_selection_pool": (
            "foldx_pass_candidates_only"
            if stability_summary.get("scored_count") and require_pass else
            "foldx_scored_candidates_only"
            if stability_summary.get("scored_count") else
            "all_proxy_candidates"
        ),
        "require_pass_for_final_panel": require_pass,
        "pass_candidate_count": sum(item.stability_status == "PASS" for item in candidates),
        "final_selection_candidate_count": len(selection_candidates),
        "final_selected_candidate_ids": [item.candidate_id for item in selected],
        "final_selected_all_scored": all(
            item.stability_value is not None for item in selected
        ) if selected else False,
    })
    designable_positions = {unit.position for unit in units if unit.designable}
    if any(
        mutation.position not in designable_positions
        for candidate in selected
        for mutation in candidate.mutations
    ):
        raise RuntimeError("Joint V3 selected a mutation outside the designable surface mask")
    selected_sector_ids = {candidate.sector_id for candidate in selected}
    projected_sectors = [
        sector for sector in sectors if sector.sector_id in selected_sector_ids
    ]
    legacy_regions = _legacy_regions(projected_sectors, units, records)
    sector_region_ids = {
        sector.sector_id: index for index, sector in enumerate(projected_sectors, start=1)
    }
    designs = _legacy_designs(
        selected, protein, original_cds, records,
        sector_region_ids,
    )
    partition_summary = {
        "engine": "joint_v3", "schema_version": "3.0", "surface_unit_count": len(units),
        "designable_surface_units": sum(item.designable for item in units),
        "total_surface_area_a2": round(sum(item.mapped_area_a2 for item in units), 3),
        "designable_surface_area_a2": round(sum(item.mapped_area_a2 for item in units if item.designable), 3),
        "microfootprint_count": len(footprints), "probe_sector_count": len(sectors),
        "coverage_plan": [
            {"patch_id": sector_region_ids[sector.sector_id], "sector_id": sector.sector_id, "probe_budget": budget}
            for sector in projected_sectors
        ],
        "all_probe_sector_count": len(sectors),
        "legacy_projected_sector_count": len(projected_sectors),
        "stage1_selected_region_ids": [item.id for item in legacy_regions],
        "stage1_target_constructs": budget, "stage1_max_constructs": budget,
        "solver_status": panel.solver_status,
        "surface_graph": graph_provenance,
        "warnings": ["ProbeSector is a coarse localization area, not a true antibody footprint."]
        + (["No membrane plane was inferred from the isolated ECD."] if topology.mode == "topology_only" else []),
        "candidate_exclusions": exclusions,
        "stability": stability_summary,
    }
    if not candidates:
        partition_summary["status"] = "NO_FEASIBLE_PROBE"
        partition_summary["warnings"].append(
            "No connected 2-site candidate survived protected-site and structure-risk rules; no design was fabricated."
        )
    else:
        partition_summary["status"] = "COMPLETE"
    design_summary = {
        "engine": "joint_v3", "schema_version": "3.0", "final_library_size": len(designs),
        "selected_panel": {"panel_id": panel.panel_id, "candidate_ids": panel.selected_candidate_ids, "achieved_coverage": panel.achieved_coverage},
        "stability_evidence": (
            "EXTERNAL_COMBINATION; FoldX BuildModel combination ΔΔG in kcal/mol"
            if stability_summary.get("scored_count") else
            "RULE_ONLY; dimensionless structural-risk proxy; no kcal/mol values"
        ),
        "anchor_mutation_policy": "not forced",
        "foldx_diagnostics": stability_summary,
        "coverage_selection_scenarios": [
            item for item in params.get("coverage_scenarios", []) if item != "direct_site"
        ],
    }
    if not designs:
        design_summary["no_feasible_design_reason"] = (
            "The configured budget and protected-site rules produced no selectable 2-4 site probe."
        )
    if stability_summary.get("scored_count") and require_pass:
        design_summary["stability_gate_policy"] = (
            "Final panel accepts only FoldX PASS combinations; ABOVE_GATE candidates are excluded."
        )
        if len(selected) < budget:
            design_summary["stability_gate_shortfall"] = {
                "requested_panel_size": budget,
                "selected_pass_candidates": len(selected),
                "pass_candidate_count": stability_summary["pass_candidate_count"],
                "reason": "insufficient_combination_stability_pass_candidates",
            }
    result = JointV3Result(
        audit=audit, topology=topology, surface_units=units, microfootprints=footprints,
        sectors=sectors, candidates=candidates, panel=panel, legacy_regions=legacy_regions,
        designs=designs, partition_summary=partition_summary, design_summary=design_summary,
        candidate_exclusions=exclusions, graph_provenance=graph_provenance,
        source_provenance=dict(source_provenance or {}),
    )
    write_joint_v3_outputs(output_dir, result, params)
    return result
