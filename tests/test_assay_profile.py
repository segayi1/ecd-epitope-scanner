from pathlib import Path

import numpy as np
from ecd_epitope_scanner.assay_profiles import (
    DEFAULT_ASSAY_PROFILE_ID,
    get_assay_profile,
)
from ecd_epitope_scanner.pipeline import (
    PipelineRequest,
    _pdb_only_context,
    _profile_orientation_geometry,
    run_pipeline,
)
from ecd_epitope_scanner.structure_analysis import best_chain_mapping, load_structure

ROOT = Path(__file__).resolve().parents[1]
PDL1_PDB = ROOT / "examples" / "pdl1_run" / "structure.pdb"


def test_fixed_profile_encodes_correct_cell_sides_and_backbone_hashes():
    profile = get_assay_profile(DEFAULT_ASSAY_PROFILE_ID)
    assert profile["experimental_roles"] == {
        "antibody_cell": "Jurkat",
        "antigen_cell": "K562",
        "readout_cell": "Jurkat",
        "readout": "Notch/Gal4-VP64 downstream reporter response",
        "interaction": "two-cell membrane-membrane contact",
    }
    assert profile["antigen_display"]["anchor_side"] == "C_terminal"
    assert profile["fab_display"]["paratope_orientation"] == "extracellular"
    assert profile["design_policy"]["stage1_maximum_constructs"] == 6
    assert profile["design_policy"]["combination_maximum_sites"] == 4
    assert profile["antigen_display"]["sha256"] == (
        "CD314E6BE64E431D21EF491D9A18E021FABCEB6C1A0EFF9D03118A7E02B576BA"
    )


def test_pdb_only_context_uses_observed_longest_chain_and_c_terminal_anchor():
    parsed = load_structure(PDL1_PDB)
    request = PipelineRequest(input_type="structure", structure_path=str(PDL1_PDB))
    profile = get_assay_profile(DEFAULT_ASSAY_PROFILE_ID)
    protein, chain_id = _pdb_only_context(parsed, request, profile)
    mapping = best_chain_mapping(parsed, protein, [chain_id])

    assert protein.ecd_start == 1
    assert protein.ecd_end == len(protein.sequence)
    assert protein.membrane_anchor_position == protein.ecd_end
    assert mapping.identity == 1.0
    assert mapping.observed_ecd_coverage == 1.0
    assert protein.raw["sequence_source"] == "observed coordinate residues"

    geometry = _profile_orientation_geometry(mapping, protein)
    assert geometry["status"] == "relative_axis_only"
    assert geometry["physical_membrane_plane_inferred"] is False
    assert np.isclose(np.linalg.norm(geometry["extracellular_axis"]), 1.0, atol=1e-5)
    assert geometry["anchor_residues"][-1] == protein.ecd_end


def test_complete_pdb_only_profile_run_needs_no_sequence_or_accession(tmp_path):
    result = run_pipeline(PipelineRequest(
        input_type="structure",
        identifier="pdb_only_profile",
        structure_path=str(PDL1_PDB),
        assay_profile=DEFAULT_ASSAY_PROFILE_ID,
        ml_mode="preview",
        output_dir=str(tmp_path),
        metadata={"disable_external_engines": True},
    ))

    probes = [item for item in result["designs"] if item["design_role"] == "STAGE1_PROBE"]
    assert result["protein"]["accession"].startswith("PDB_ONLY_CHAIN_")
    assert result["assay_profile"]["experimental_roles"]["antibody_cell"] == "Jurkat"
    assert result["assay_profile"]["experimental_roles"]["antigen_cell"] == "K562"
    assert result["assay_geometry"]["physical_membrane_plane_inferred"] is False
    assert 1 <= len(probes) <= 6
    assert result["partitioning"]["stage1_max_constructs"] == 6
    assert all(
        not item["available"] and item["failure_reason"] == "disabled by request"
        for item in result["design_summary"]["tool_inventory"].values()
    )
    assert result["partitioning"]["epitope_prior_used_for_panel_selection"] is False
    assert all(len(item["mutations"]) <= 4 for item in probes)
    assert all(region["predicted_accessibility"] == "unknown" for region in result["regions"])
    assert (tmp_path / "fixed_assay_profile.json").is_file()
    assert "fixed_assay_profile.json" in result["files"]
