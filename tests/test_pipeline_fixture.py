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
    assert 5 <= len(result["regions"]) <= 8
    assert result["partitioning"]["coverage"] >= 0.9
    assert len(result["designs"]) == len(result["regions"])
    region_positions = {region["id"]: set(region["residues"]) for region in result["regions"]}
    protected = {
        residue["uniprot_pos"]
        for residue in result["residues"]
        if residue["hard_exclusions"]
    }
    for design in result["designs"]:
        assert design["status"] == "RECOMMENDED_LOW_PREDICTED_RISK"
        assert 2 <= len(design["mutations"]) <= 5
        assert len({item["position"] for item in design["mutations"]}) == len(design["mutations"])
        for mutation in design["mutations"]:
            assert mutation["wt"] != mutation["mutant"]
            assert mutation["position"] in region_positions[design["region_id"]]
            assert mutation["position"] not in protected
        assert translate_cds(design["mutant_coding_sequence"]) == design["mutant_ecd_sequence"]
    assert (output / "ecd_epitope_design.zip").is_file()
    assert (output / "run_manifest.json").is_file()
    assert (output / "methods_report.md").is_file()


@pytest.mark.integration
def test_saturation_is_internal_but_final_library_is_one_per_region(fixture_result):
    result, _ = fixture_result
    assert result["design_summary"]["in_silico_substitutions_evaluated"] > 100
    assert result["design_summary"]["final_library_size"] == len(result["regions"])
    for design in result["designs"]:
        assert design["saturation_scan"]
        assert all(len(choices) == 19 for choices in design["saturation_scan"].values())
