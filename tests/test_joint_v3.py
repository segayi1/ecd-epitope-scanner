from pathlib import Path

import numpy as np
from ecd_epitope_scanner.joint_design import engine as v3_engine
from ecd_epitope_scanner.joint_design import run_joint_v3
from ecd_epitope_scanner.joint_design import stability as v3_stability
from ecd_epitope_scanner.joint_design.coverage import attach_coverage
from ecd_epitope_scanner.joint_design.engine import _legacy_designs
from ecd_epitope_scanner.joint_design.models import (
    MicroFootprint,
    MutationOption,
    ProbeCandidate,
    ResidueId,
    SurfaceResidueUnit,
)
from ecd_epitope_scanner.joint_design.optimizer import select_probe_panel
from ecd_epitope_scanner.joint_design.sectors import build_probe_sectors
from ecd_epitope_scanner.joint_design.surface_units import build_surface_units
from ecd_epitope_scanner.models import ProteinContext, ResidueRecord
from ecd_epitope_scanner.stability import EngineInventory
from ecd_epitope_scanner.surface_patch.models import SurfaceGeometry


def _records(count: int = 18) -> list[ResidueRecord]:
    sequence = "DEKRFWYAILMSTNQ" * 2
    return [
        ResidueRecord(
            uniprot_pos=index + 1,
            aa=sequence[index],
            chain_id="A",
            auth_seq_id=index + 1,
            residue_name="ALA" if sequence[index] == "A" else "GLY",
            ca=(float(index * 3.7), 0.0, 0.0),
            sasa=100.0,
            sidechain_sasa=40.0,
            rsasa=0.5,
            surface_accessible=True,
            folding_screen_pass=True,
        )
        for index in range(count)
    ]


def test_joint_v3_is_deterministic_and_does_not_force_anchor(tmp_path: Path) -> None:
    records = _records()
    protein = ProteinContext(
        accession="TEST",
        protein_name="synthetic ECD",
        organism="test",
        sequence="".join(item.aa for item in records),
        ecd_start=1,
        ecd_end=len(records),
        membrane_anchor_position=len(records),
        topology="C-terminal anchor",
    )
    parameters = {
        "schema_version": "3.0",
        "max_probes": 2,
        "mutations_per_probe": [2, 3],
        "footprint_target_areas_a2": [160.0],
        "candidate_search": {
            "site_pool_cap_per_sector": 8,
            "mutation_options_cap_per_site": 2,
            "beam_width": 20,
            "shortlist_per_sector": 8,
        },
        "surface_graph": {"neighbor_cutoff_a": 10.0, "max_neighbors": 4},
        "solver": {"time_limit_seconds": 1},
    }
    first = run_joint_v3(records, protein, tmp_path / "first", parameters=parameters)
    second = run_joint_v3(records, protein, tmp_path / "second", parameters=parameters)

    assert first.panel.selected_candidate_ids == second.panel.selected_candidate_ids
    assert all(2 <= len(item.mutations) <= 3 for item in first.designs)
    assert all(item.anchor_notation is None for item in first.designs)
    assert first.topology.mode == "topology_only"
    assert first.topology.membrane_plane is None
    assert (tmp_path / "first" / "scientific_report.md").is_file()


def test_joint_v3_keeps_protected_surface_in_denominator(tmp_path: Path) -> None:
    records = _records()
    records[3].disulfide = True
    records[3].hard_exclusions = ["disulfide-cysteine"]
    protein = ProteinContext(
        accession="TEST",
        protein_name="synthetic ECD",
        organism="test",
        sequence="".join(item.aa for item in records),
        ecd_start=1,
        ecd_end=len(records),
        membrane_anchor_position=len(records),
        topology="C-terminal anchor",
    )
    result = run_joint_v3(
        records,
        protein,
        tmp_path,
        parameters={"max_probes": 1, "mutations_per_probe": [2, 2], "footprint_target_areas_a2": [160]},
    )
    assert result.partition_summary["surface_unit_count"] >= result.partition_summary["designable_surface_units"]
    protected = next(item for item in result.surface_units if item.position == 4)
    assert not protected.designable
    assert "disulfide-cysteine" in protected.exclusion_reasons


def test_joint_v3_candidates_cover_all_configured_sizes_without_duplicates(tmp_path: Path) -> None:
    records = _records(24)
    records[4].interface = True
    protein = ProteinContext(
        accession="TEST", protein_name="synthetic ECD", organism="test",
        sequence="".join(item.aa for item in records), ecd_start=1, ecd_end=len(records),
        membrane_anchor_position=len(records), topology="C-terminal anchor",
    )
    result = run_joint_v3(
        records, protein, tmp_path,
        parameters={
            "max_probes": 2,
            "mutations_per_probe": [2, 4],
            "footprint_target_areas_a2": [600],
            "max_footprint_centers": 12,
            "candidate_search": {
                "site_pool_cap_per_sector": 10,
                "mutation_options_cap_per_site": 2,
                "beam_width": 100,
                "shortlist_per_sector": 30,
                "maximum_geodesic_span_a": 24,
            },
            "surface_graph": {"neighbor_cutoff_a": 10.0, "max_neighbors": 6},
            "coverage_scenarios": ["geodesic_tight", "geodesic_medium", "geodesic_wide"],
            "solver": {"time_limit_seconds": 2},
        },
    )
    sizes = {len(item.mutations) for item in result.candidates}
    assert sizes == {2, 3, 4}
    signatures = [
        tuple(sorted((mutation.position, mutation.mutant_aa) for mutation in item.mutations))
        for item in result.candidates
    ]
    assert len(signatures) == len(set(signatures))
    assert all(mutation.position != 5 for item in result.candidates for mutation in item.mutations)
    assert all(item.connected for item in result.candidates)


def test_residue_fallback_graph_rejects_opposing_surface_normals() -> None:
    records = _records(2)
    records[0].ca = (0.0, 0.0, 0.0)
    records[1].ca = (4.0, 0.0, 0.0)
    records[0].surface_normal = (1.0, 0.0, 0.0)
    records[1].surface_normal = (-1.0, 0.0, 0.0)
    _, graph, _ = build_surface_units(records, neighbor_cutoff_a=10.0)
    assert graph[1] == {}
    assert graph[2] == {}


def test_explicit_surface_graph_contracts_non_design_transit_residue() -> None:
    records = _records(4)
    records[1].surface_accessible = False
    records[1].rsasa = 0.1
    geometry = SurfaceGeometry(
        vertices=np.asarray([(0.0, 1.0, 0.0), (3.7, 1.0, 0.0), (7.4, 1.0, 0.0)]),
        normals=np.asarray([(0.0, 1.0, 0.0)] * 3),
        vertex_area=np.asarray([20.0] * 3),
        edges=np.asarray([(0, 1), (1, 2)], dtype=np.int64),
        edge_length=np.asarray([3.7, 3.7]),
        owner_residue=np.asarray([1, 2, 3], dtype=np.int64),
        owner_atom=np.asarray([0, 1, 2], dtype=np.int64),
        chain_id=np.asarray(["A"] * 3, dtype=object),
        canonical_residue_id=np.asarray([1, 2, 3], dtype=np.int64),
        is_ecd=np.asarray([True] * 3),
        confidence=np.asarray([100.0] * 3),
        component_id=np.asarray([1] * 3, dtype=np.int64),
    )
    provenance = {}
    _, graph, _ = build_surface_units(records, surface_geometry=geometry, graph_provenance=provenance)
    assert graph == {1: {3: 7.4}, 3: {1: 7.4}}
    assert provenance["transit_surface_residue_count"] == 1
    assert provenance["reliability_omitted_eligible_count"] == 1
    assert provenance["component_count"] == 1


def test_probe_sectors_do_not_merge_far_three_dimensional_segments() -> None:
    residues = [ResidueId("0", "A", index, "", "A") for index in (1, 2, 3)]
    units = [
        SurfaceResidueUnit(
            residue,
            index,
            0.5,
            100.0,
            True,
            [],
            coordinate_xyz=(float((index - 1) * 5), 0.0, 0.0),
        )
        for index, residue in enumerate(residues, start=1)
    ]
    footprints = [
        MicroFootprint("F1", residues[0], "small", 100.0, residues[:2], 100.0),
        MicroFootprint("F2", residues[1], "small", 100.0, residues[1:], 100.0),
    ]
    graph = {1: {2: 5.0}, 2: {1: 5.0, 3: 5.0}, 3: {2: 5.0}}
    sectors = build_probe_sectors(
        units,
        footprints,
        max_sectors=4,
        graph=graph,
        maximum_geodesic_diameter_a=20.0,
        maximum_euclidean_diameter_a=6.0,
    )
    assert len(sectors) == 2
    assert {tuple(item.residue_ids) for item in sectors} == {(1, 2), (2, 3)}


def test_footprint_coverage_uses_unique_units_and_designable_denominator() -> None:
    residue_ids = [ResidueId("0", "A", index, "", "A") for index in (1, 2, 3)]
    units = [
        SurfaceResidueUnit(residue, index, 0.5, 100.0, True, [])
        for index, residue in enumerate(residue_ids, start=1)
    ]
    mutation = MutationOption(
        residue_ids[0], 1, "A", "S", "polar_substitution", 0.1, 0.5,
        None, None, "PROXY_ONLY",
    )
    candidate = ProbeCandidate("C1", "S1", [mutation], True)
    footprints = [
        MicroFootprint("F1", residue_ids[0], "a", 200, residue_ids[:2], 200),
        MicroFootprint("F2", residue_ids[0], "b", 300, residue_ids, 300),
    ]
    graph = {1: {2: 4.0}, 2: {1: 4.0, 3: 4.0}, 3: {2: 4.0}}
    attach_coverage([candidate], units, graph, footprints)
    assert candidate.coverage_units_by_scenario["footprint_hit_1"] == [1, 2, 3]
    assert candidate.coverage_by_scenario["footprint_hit_1"]["designable_area_a2"] == 300.0
    assert candidate.coverage_by_scenario["footprint_hit_1"]["fraction_of_designable_surface"] == 1.0


def test_milp_uses_scenario_specific_coverage_variables() -> None:
    residues = [ResidueId("0", "A", index, "", "A") for index in (1, 2, 3)]
    units = [
        SurfaceResidueUnit(residue, index, 0.5, area, True, [])
        for index, (residue, area) in enumerate(zip(residues, (100.0, 100.0, 60.0)), start=1)
    ]
    def mutation(index: int) -> MutationOption:
        return MutationOption(
            residues[index - 1], index, "A", "S", "test", 0.1, 0.5,
            None, None, "PROXY_ONLY",
        )
    cross_scenario = ProbeCandidate(
        "A", "S1", [mutation(1)], True,
        coverage_units_by_scenario={"tight": [1], "wide": [2], "direct_site": [1]},
    )
    intersection_only = ProbeCandidate(
        "B", "S2", [mutation(3)], True,
        coverage_units_by_scenario={"tight": [3], "wide": [3], "direct_site": [3]},
    )
    panel, selected = select_probe_panel(
        [cross_scenario, intersection_only], units, 1,
        timeout_seconds=2, coverage_scenarios=["tight", "wide"],
    )
    assert [item.candidate_id for item in selected] == ["A"]
    assert panel.achieved_coverage["tight"] == 0.3846
    assert panel.achieved_coverage["wide"] == 0.3846


def test_v3_foldx_adapter_scores_complete_combinations_and_records_provenance(tmp_path: Path, monkeypatch) -> None:
    residue = ResidueId("0", "A", 12, "", "K")
    mutation = MutationOption(
        residue, 12, "K", "A", "charge_removal", 0.1, 0.7,
        None, None, "PROXY_ONLY",
    )
    candidate = ProbeCandidate("PC_1", "PS001", [mutation], True)
    repaired = tmp_path / "repaired.pdb"
    repaired.write_text("ATOM\n", encoding="ascii")
    monkeypatch.setattr(v3_stability, "repair_with_foldx", lambda *args, **kwargs: repaired)
    captured = {}

    def fake_scores(binary, pdb_path, groups, workdir, **kwargs):
        captured["groups"] = groups
        captured["kwargs"] = kwargs
        return [1.25]

    monkeypatch.setattr(v3_stability, "score_many_with_foldx", fake_scores)
    summary = v3_stability.annotate_stability_evidence(
        [candidate],
        inventory=EngineInventory(True, None, "foldx.exe"),
        structure_pdb=tmp_path / "structure.pdb",
        workdir=tmp_path,
        parameters={"foldx_runs": 3, "combination_gate_kcal_mol": 2.0},
    )
    assert summary["scored_count"] == 1
    assert candidate.evidence_level == "EXTERNAL_COMBINATION"
    assert candidate.stability_value == 1.25
    assert candidate.stability_unit == "kcal/mol"
    assert candidate.stability_status == "PASS"
    assert captured["groups"] == [[("K", 12, "A", "A")]]
    assert captured["kwargs"]["runs"] == 3
    assert candidate.stability_provenance["combination_level"] is True
    assert candidate.stability_provenance["single_site_sum_used"] is False


def test_v3_foldx_cap_prioritizes_proxy_panel_candidates(tmp_path: Path, monkeypatch) -> None:
    candidates = []
    for index, sector in enumerate(("PS001", "PS001", "PS002"), start=1):
        residue = ResidueId("0", "A", index, "", "K")
        mutation = MutationOption(
            residue, index, "K", "A", "charge_removal", 0.1, 0.7,
            None, None, "PROXY_ONLY",
        )
        candidates.append(ProbeCandidate(f"PC_{index}", sector, [mutation], True))
    repaired = tmp_path / "repaired.pdb"
    repaired.write_text("ATOM\n", encoding="ascii")
    monkeypatch.setattr(v3_stability, "repair_with_foldx", lambda *args, **kwargs: repaired)
    captured = {}

    def fake_scores(binary, pdb_path, groups, workdir, **kwargs):
        captured["positions"] = [group[0][1] for group in groups]
        return [-0.2] * len(groups)

    monkeypatch.setattr(v3_stability, "score_many_with_foldx", fake_scores)
    summary = v3_stability.annotate_stability_evidence(
        candidates,
        inventory=EngineInventory(True, None, "foldx.exe"),
        structure_pdb=tmp_path / "structure.pdb",
        workdir=tmp_path,
        parameters={"foldx_max_candidates": 2},
        priority_candidate_ids=["PC_3"],
    )
    assert 3 in captured["positions"]
    assert candidates[2].stability_value == -0.2
    assert summary["priority_submitted_count"] == 1
    assert summary["deferred_count"] == 1
    assert summary["status"] == "COMPLETE_CAPPED"


def test_v3_legacy_design_exposes_combination_foldx_energy() -> None:
    records = _records(3)
    protein = ProteinContext(
        accession="TEST",
        protein_name="synthetic ECD",
        organism="test",
        sequence="".join(item.aa for item in records),
        ecd_start=1,
        ecd_end=3,
        membrane_anchor_position=3,
        topology="C-terminal anchor",
    )
    residue = ResidueId("0", "A", 1, "", records[0].aa)
    mutation = MutationOption(
        residue, 1, records[0].aa, "A", "charge_removal", 0.1, 0.7,
        None, None, "PROXY_ONLY",
    )
    candidate = ProbeCandidate(
        "PC_1", "PS001", [mutation], True,
        evidence_level="EXTERNAL_COMBINATION",
        stability_value=-0.45,
        stability_unit="kcal/mol",
        stability_status="PASS",
        stability_provenance={"engine": "FoldX", "combination_level": True},
    )
    design = _legacy_designs([candidate], protein, None, records, {"PS001": 1})[0]
    assert design.combined_stability_score == -0.45
    assert design.score_origin == "foldx"
    assert "FoldX" in design.stability_engine
    assert design.foldx_diagnostics["combination_level"] is True
    assert all(choice.foldx_ddg is None for choice in design.mutations)


def test_v3_final_panel_is_selected_from_scored_candidates(tmp_path: Path, monkeypatch) -> None:
    records = _records()
    protein = ProteinContext(
        accession="TEST",
        protein_name="synthetic ECD",
        organism="test",
        sequence="".join(item.aa for item in records),
        ecd_start=1,
        ecd_end=len(records),
        membrane_anchor_position=len(records),
        topology="C-terminal anchor",
    )
    captured = {}

    def fake_stability(candidates, **kwargs):
        priority = list(kwargs["priority_candidate_ids"])
        captured["priority"] = priority
        selected_ids = set(priority[:2])
        for candidate in candidates:
            if candidate.candidate_id in selected_ids:
                is_pass = candidate.candidate_id == priority[0]
                candidate.stability_value = -0.25 if is_pass else 3.0
                candidate.stability_unit = "kcal/mol"
                candidate.stability_status = "PASS" if is_pass else "ABOVE_GATE"
                candidate.evidence_level = "EXTERNAL_COMBINATION"
                candidate.stability_provenance = {
                    "engine": "FoldX", "combination_level": True,
                }
            else:
                candidate.stability_status = "NOT_RUN_CAP"
        return {
            "enabled": True,
            "engine": "FoldX",
            "status": "COMPLETE_CAPPED",
            "scored_count": len(selected_ids),
            "all_failed": False,
            "batch_constructs": len(selected_ids),
        }

    monkeypatch.setattr(v3_engine, "annotate_stability_evidence", fake_stability)
    result = run_joint_v3(
        records,
        protein,
        tmp_path,
        parameters={
            "max_probes": 2,
            "mutations_per_probe": [2, 3],
            "footprint_target_areas_a2": [160.0],
            "candidate_search": {
                "site_pool_cap_per_sector": 8,
                "mutation_options_cap_per_site": 2,
                "beam_width": 20,
                "shortlist_per_sector": 8,
            },
            "surface_graph": {"neighbor_cutoff_a": 10.0, "max_neighbors": 4},
            "solver": {"time_limit_seconds": 1},
        },
    )
    assert captured["priority"]
    assert result.designs
    assert all(design.score_origin == "foldx" for design in result.designs)
    assert all(design.combined_stability_score == -0.25 for design in result.designs)
    assert result.design_summary["foldx_diagnostics"]["final_selected_all_scored"] is True
    assert result.design_summary["foldx_diagnostics"]["final_selection_pool"] == "foldx_pass_candidates_only"
    assert result.design_summary["stability_gate_policy"]
    assert result.design_summary["stability_gate_shortfall"]["selected_pass_candidates"] == 1
