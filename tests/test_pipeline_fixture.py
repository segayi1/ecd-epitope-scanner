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
        output_dir=str(output),
        metadata={"disable_external_engines": True},
    )
    return run_pipeline(request), output


@pytest.mark.integration
def test_full_pipeline_outputs_meaningful_regional_constructs(fixture_result):
    result, output = fixture_result
    assert result["summary"]["mapped_ecd_residues"] == 220
    # Patch count is data-driven but bounded by the panel budget (maximum_regions).
    assert 1 <= len(result["regions"]) <= 8
    assert result["partitioning"]["coverage"] >= 0.85
    region_positions = {region["id"]: set(region["residues"]) for region in result["regions"]}
    protected = {
        residue["uniprot_pos"]
        for residue in result["residues"]
        if residue["hard_exclusions"]
    }

    primaries = [design for design in result["designs"] if design["design_role"] == "STAGE1_PROBE"]
    assert len(primaries) == len(result["regions"])
    for region, design in zip(result["regions"], primaries):
        assert design["region_id"] == region["id"]
        assert design["mutations"], "every testable patch must yield a probe or an explicit refusal"
        if design["status"] == "RECOMMENDED_LOW_PREDICTED_RISK":
            assert 2 <= len(design["mutations"]) <= 6
            assert design["anchor_notation"], "every Stage-1 probe must carry an anchor mutation"
            assert len({item["position"] for item in design["mutations"]}) == len(design["mutations"])
            for mutation in design["mutations"]:
                assert mutation["wt"] != mutation["mutant"]
                assert mutation["position"] in region_positions[region["id"]]
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
    assert result["partitioning"]["method"].startswith("Fab-footprint")
    for region in result["regions"]:
        assert region["epitope_prior"] is not None
        assert region["tier"] in {"HIGH", "MEDIUM", "LOW"}
        assert region["prior_components"]
        assert region["required_probe_coverage"] in {1, 2}
    tiers = {region["tier"] for region in result["regions"]}
    assert "HIGH" in tiers


@pytest.mark.integration
def test_saturation_remains_full_but_lives_on_probes(fixture_result):
    result, _ = fixture_result
    assert result["design_summary"]["in_silico_substitutions_evaluated"] > 100
    probes = [design for design in result["designs"] if design["design_role"] == "STAGE1_PROBE"]
    with_scan = [design for design in probes if design["saturation_scan"]]
    assert with_scan
    for design in with_scan:
        assert all(len(choices) == 19 for choices in design["saturation_scan"].values())
        for choices in design["saturation_scan"].values():
            classes = {choice["menu_class"] for choice in choices}
            assert classes, "every substitution carries its directed menu class"
