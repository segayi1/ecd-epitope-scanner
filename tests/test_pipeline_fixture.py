import json
from pathlib import Path

import pytest
from ecd_epitope_scanner.pipeline import PipelineRequest, run_pipeline
from ecd_epitope_scanner.sequences import translate_cds

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def fixture_result(tmp_path_factory):
    old = json.loads((ROOT / "examples" / "pdl1_run" / "analysis.json").read_text(encoding="utf-8"))
    output = tmp_path_factory.mktemp("pipeline")
    request = PipelineRequest(
        input_type="structure",
        identifier="pdl1_fixture",
        structure_path=str(ROOT / "examples" / "pdl1_run" / "structure.pdb"),
        sequence=old["sequence"],
        ecd_start=19,
        ecd_end=238,
        allow_partial=True,
        ml_mode="preview",
        output_dir=str(output),
        metadata={"disable_external_engines": True},
    )
    return run_pipeline(request), output


@pytest.mark.integration
def test_full_pipeline_outputs_meaningful_regional_constructs(fixture_result):
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
        assert set(_distance(graph,min(owned),owned))==owned
        assert sector['geodesic_diameter_a']<=policy['maximum_geodesic_diameter_a']+.001
        assert sector['euclidean_diameter_a']<=policy['maximum_euclidean_diameter_a']+.001
    panel = result["partitioning"]["selected_panel"]
    assert panel["coverage_denominator_residues"] > 0
    assert panel["achieved_coverage"]
    assert all(0.0 <= value <= 1.0 for value in panel["achieved_coverage"].values())
    assert result["partitioning"]["stage1_max_constructs"] == result['parameters']['stage1_construct_budget']['maximum'] <= 8
    region_positions = {region["id"]: set(region["residues"]) for region in result["regions"]}
    protected = {
        residue["uniprot_pos"]
        for residue in result["residues"]
        if residue["hard_exclusions"]
    }

    primaries = [design for design in result["designs"] if design["design_role"] == "STAGE1_PROBE"]
    assert 0 < len(primaries) <= result["partitioning"]["stage1_max_constructs"]
    for design in primaries:
        assert design["region_id"] in region_positions
        assert design["mutations"], "every testable patch must yield a probe or an explicit refusal"
        if design["status"] == "RECOMMENDED_LOW_PREDICTED_RISK":
            assert 2 <= len(design["mutations"]) <= 6
            assert design["anchor_notation"], "every Stage-1 probe must carry an anchor mutation"
            assert len({item["position"] for item in design["mutations"]}) == len(design["mutations"])
            for mutation in design["mutations"]:
                assert mutation["wt"] != mutation["mutant"]
                assert mutation["position"] in region_positions[design["region_id"]]
                assert mutation["position"] not in protected
            assert translate_cds(design["mutant_coding_sequence"]) == design["mutant_ecd_sequence"]

    stage2 = [design for design in result["designs"] if design["design_role"].startswith("STAGE2")]
    assert stage2, "Stage-2 deconvolution presets must be generated alongside the panel"
    for design in stage2:
        assert design["region_id"] is None
        assert design["parent_design_id"]
        assert translate_cds(design["mutant_coding_sequence"]) == design["mutant_ecd_sequence"]

    assert (output / "ecd_epitope_design.zip").is_file()
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
    assert (output / "full_length_protein_sequences.fasta").is_file()


@pytest.mark.integration
def test_epitope_mapping_plan_is_complete(fixture_result):
    result, output = fixture_result
    plan = result["epitope_mapping_plan"]
    assert plan["stage1_probes"], "the plan must list every Stage-1 probe"
    assert len(plan["coverage_matrix"]) >= len(plan["stage1_probes"])
    assert isinstance(plan["blind_spots"], list)
    assert "Stage 1" in plan["assay_note"]
    patch_ids = {region["id"] for region in result["regions"]}
    for probe in plan["stage1_probes"]:
        assert probe["patch_id"] in patch_ids


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
def test_saturation_remains_full_but_lives_on_probes(fixture_result):
    result, _ = fixture_result
    assert result["design_summary"]["saturation_scan_status"] == "ENUMERATED_NOT_SCORED"
    assert result["design_summary"]["unscored_saturation_substitution_count"] > 100
    probes = [design for design in result["designs"] if design["design_role"] == "STAGE1_PROBE"]
    with_scan = [design for design in probes if design["saturation_scan"]]
    assert with_scan
    for design in with_scan:
        assert all(len(choices) == 19 for choices in design["saturation_scan"].values())
        for choices in design["saturation_scan"].values():
            classes = {choice["menu_class"] for choice in choices}
            assert classes == {"unscored_saturation_option"}
            assert all(choice["status"] == "NOT_SCORED" for choice in choices)
