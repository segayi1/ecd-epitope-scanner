import json
import zipfile
from pathlib import Path

import pytest
from ecd_epitope_scanner.pipeline import PipelineRequest, run_pipeline

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def fixture_result(tmp_path_factory):
    lines = (ROOT / "examples" / "pdl1_run" / "canonical.fasta").read_text(
        encoding="utf-8"
    ).splitlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    output = tmp_path_factory.mktemp("pipeline")
    request = PipelineRequest(
        input_type="structure",
        identifier="pdl1_fixture",
        structure_path=str(ROOT / "examples" / "pdl1_run" / "structure.pdb"),
        sequence=sequence,
        ecd_start=19,
        ecd_end=238,
        allow_partial=True,
        ml_mode="preview",
        output_dir=str(output),
        metadata={"disable_external_engines": True},
    )
    return run_pipeline(request), output


@pytest.mark.integration
def test_full_pipeline_outputs_auditable_preview(fixture_result):
    result, output = fixture_result
    assert result["summary"]["mapped_ecd_residues"] == 220
    # Display count is soft and independent of the construct ceiling. Measure
    # connectivity and diameter instead of forcing unsafe region merges.
    sectors=result['joint_v4_ml']['sectors']
    graph_data=json.loads((output/'surface_graph_v4.json').read_text(encoding='utf-8'))
    graph={r['uniprot_pos']:{} for r in graph_data['nodes']}
    for edge in graph_data['edges']:
        graph[edge['source']][edge['target']]=edge['distance_a']
        graph[edge['target']][edge['source']]=edge['distance_a']
    from ecd_epitope_scanner.joint_design_v4.sectors import _distance
    policy=graph_data['graph_provenance']['display_partition_policy']
    positions=[p for sector in sectors for p in sector['residue_positions']]
    assert len(result['regions'])==len(sectors)>0
    assert sorted(positions)==sorted(graph) and len(positions)==len(set(positions))
    for sector in sectors:
        owned=set(sector['residue_positions'])
        assert len(owned) >= policy['minimum_sector_residues']
        assert sector['total_area_a2'] >= policy['minimum_sector_area_a2']
        assert set(_distance(graph,min(owned),owned))==owned
        assert sector['geodesic_diameter_a']<=policy['maximum_geodesic_diameter_a']+.001
        assert sector['euclidean_diameter_a']<=policy['maximum_euclidean_diameter_a']+.001
    panel = result["partitioning"]["selected_panel"]
    assert result["partitioning"]["coverage"] == 1.0
    assert result["partitioning"]["coverage_basis"] == "eligible_reliable_surface_graph_nodes"
    assert panel["coverage_denominator_residues"] > 0
    assert panel["achieved_coverage"] == {}
    assert result["partitioning"]["stage1_max_constructs"] == result['parameters']['stage1_construct_budget']['maximum'] <= 8
    diagnostics = result["joint_v4_ml"]["recommendation_diagnostics"]
    assert diagnostics["outcome"] == "SINGLE_STABILITY_EVIDENCE_UNAVAILABLE"
    assert diagnostics["candidate_count"] > 0
    assert diagnostics["formal_eligible_candidate_count"] == 0
    assert result["designs"] == []
    assert all(
        "COMBINATION_SPURS_MISSING" in candidate["formal_exclusion_reasons"]
        and "COMBINATION_FOLDX_MISSING" in candidate["formal_exclusion_reasons"]
        for candidate in result["joint_v4_ml"]["candidates"]
    )
    assert result["joint_v4_ml"]["recommendation_status"] == (
        "PREVIEW_NOT_FOR_FORMAL_RECOMMENDATION"
    )

    assert (output / "ecd_epitope_design.zip").is_file()
    with zipfile.ZipFile(output / "ecd_epitope_design.zip") as archive:
        for name in ("partition_diagnostics_v4.json", "surface_graph_v4.json",
                     "surface_display_regions_v41.json", "probe_sectors_v4.json",
                     "sector_recommendations_v4.json", "sector_recommendations_v4.csv",
                     "execution_protocol.json"):
            assert name in result["files"]
            assert archive.read(name) == (output / name).read_bytes()
    suggestions = result["joint_v4_ml"]["selected_panel"]["sector_recommendations"]
    assert {row["sector_id"] for row in suggestions} == {sector["sector_id"] for sector in sectors}
    assert all(row["candidate_id"] is None and not row["formal_eligible"] for row in suggestions)
    assert suggestions == json.loads((output / "sector_recommendations_v4.json").read_text(encoding="utf-8"))
    assert all(row["failed_gates"] for row in suggestions if row["status"] == "NEAR_THRESHOLD")
    assert all(candidate["geometry_support"]["area_basis"] == "eligible_surface_area_a2"
               for candidate in result["joint_v4_ml"]["candidates"])
    # Preview lacks single-model evidence; it must not invent mechanism recommendations or JSON infinities.
    assert not any(candidate.get("impact_evidence", {}).get("policy") == "regional_mechanism_pareto_v1"
                   for candidate in result["joint_v4_ml"]["candidates"])
    def reject_nonfinite(value):
        raise AssertionError(f"Non-finite candidate export: {value}")
    json.loads((output / "panel_candidates.json").read_text(encoding="utf-8"), parse_constant=reject_nonfinite)
    assert all(set(row["geometry_support"]["scenarios"]) == {"geodesic_9", "geodesic_12", "geodesic_15"}
               for row in suggestions if row["candidate_id"])
    assert (output / "run_manifest.json").is_file()
    reliability = result["structure_reliability"]
    assert reliability["qc"]["state_counts"] == {"STRUCTURED": 220}
    manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    for name in ("residue_reliability.csv", "structural_segments.csv", "surface_mask.csv", "idr_candidates.csv", "reliability_qc.json"):
        assert name in manifest["files"] and name in result["files"]
    assert (output / "methods_report.md").is_file()
    offline = output / "offline_report.html"
    assert offline.is_file()
    offline_text = offline.read_text(encoding="utf-8")
    assert "__DATA__" not in offline_text
    assert "__PDB__" not in offline_text
    assert "file://" not in offline_text
    assert (output / "panel_coverage_matrix.csv").is_file()
    assert not (output / "full_length_protein_sequences.fasta").exists()


@pytest.mark.integration
def test_epitope_mapping_plan_reports_missing_evidence(fixture_result):
    result, output = fixture_result
    plan = result["epitope_mapping_plan"]
    assert plan["stage1_probes"] == []
    assert isinstance(plan["blind_spots"], list)
    assert "Stage 1" in plan["assay_note"]
    assert (
        result["joint_v4_ml"]["recommendation_diagnostics"]["outcome"]
        == "SINGLE_STABILITY_EVIDENCE_UNAVAILABLE"
    )


@pytest.mark.integration
def test_patch_priors_and_tiers_are_present(fixture_result):
    result, _ = fixture_result
    assert result["partitioning"]["method"].startswith("Joint V4-ML")
    for region in result["regions"]:
        assert region["epitope_prior"] is not None
        assert region["tier"] in {"HIGH", "MEDIUM", "LOW"}
        assert region["prior_components"]
        assert region["required_probe_coverage"] in {1, 2}
    tiers = {region["tier"] for region in result["regions"]}
    has_epitope_scores = any(
        item.get("raw_score") is not None
        for item in result["joint_v4_ml"]["epitope_evidence"]
    )
    if not has_epitope_scores:
        assert tiers == {"LOW"}, "missing external epitope scores must not create a HIGH-priority tier"


@pytest.mark.integration
def test_saturation_scan_is_not_fabricated_without_scoring(fixture_result):
    result, _ = fixture_result
    assert result["design_summary"]["saturation_scan_status"] == "ENUMERATED_NOT_SCORED"
    assert result["design_summary"]["unscored_saturation_substitution_count"] == 0
    diagnostics = result["joint_v4_ml"]["recommendation_diagnostics"]
    assert diagnostics["mutation_option_count"] > 0
    assert diagnostics["eligible_mutation_option_count"] > 0
    assert result["designs"] == []
