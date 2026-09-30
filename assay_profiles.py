from __future__ import annotations

from copy import deepcopy
from typing import Any

# This profile is intentionally encoded in source rather than inferred anew for
# each run.  The hashes make the biological assumption auditable: changing a
# backbone requires a new profile/version and a deliberate review.
JURKAT_FAB_K562_ANTIGEN_V1 = {
    "id": "jurkat_fab_k562_antigen_v1",
    "version": 1,
    "display_name": "Jurkat membrane-Fab x K562 membrane-antigen",
    "experimental_roles": {
        "antibody_cell": "Jurkat",
        "antigen_cell": "K562",
        "readout_cell": "Jurkat",
        "readout": "Notch/Gal4-VP64 downstream reporter response",
        "interaction": "two-cell membrane-membrane contact",
    },
    "antigen_display": {
        "plasmid": "ZP18-pHR_myc-CTLA4-PDGFR-FUS-7A-T2A-puro.dna",
        "sha256": "CD314E6BE64E431D21EF491D9A18E021FABCEB6C1A0EFF9D03118A7E02B576BA",
        "cell": "K562",
        "architecture_n_to_c": [
            "signal peptide",
            "Myc tag",
            "antigen ECD insert",
            "PDGFRB extracellular stalk/TM",
            "PDGFRB intracellular segment",
            "FUS-7A",
        ],
        "extracellular_terminus": "N",
        "membrane_proximal_ecd_terminus": "C",
        "anchor_side": "C_terminal",
        "snapgene_feature_evidence": {
            "coordinate_system": "0-based half-open nucleotide intervals",
            "leader": [8296, 8359],
            "myc": [8359, 8389],
            "antigen_insert_in_supplied_backbone": [8389, 8770],
            "pdgfrb_anchor_region": [8770, 9034],
            "fus_7a": [9034, 9676],
        },
        "expression_control": (
            "N-terminal Myc reports tag accessibility/surface display; it does not by itself "
            "establish native ECD folding."
        ),
    },
    "fab_display": {
        "cell": "Jurkat",
        "light_chain_plasmid": "PYW320_pHR_pGK_flag_mIgK_ALZ_AscI_T2A-neo.dna",
        "light_chain_sha256": "C129CBA84F8880210E6044FB309B7F945E5727C0485976944376E2290F7573E1",
        "light_chain_architecture_n_to_c": [
            "signal peptide", "FLAG", "VL cloning site", "mouse IgK constant", "ALZ"
        ],
        "light_chain_snapgene_feature_evidence": {
            "coordinate_system": "0-based half-open nucleotide intervals",
            "leader": [8276, 8336],
            "flag": [8336, 8360],
            "mouse_igk_constant": [8369, 8690],
            "alz": [8702, 8792],
        },
        "heavy_chain_plasmid": (
            "PYW734_pHR-EV-mIgG1-CH1-BLZ-CD8a-hnotch1TMD-notch2-JMD-T2A-puro.dna"
        ),
        "heavy_chain_sha256": "87A35B4FC87C9ACD80FEEFC24D9D7A6DE2D630F3765B12C64B13B26008A7AB6F",
        "heavy_chain_architecture_n_to_c": [
            "signal peptide",
            "Myc",
            "VH cloning site",
            "mouse IgG1 CH1",
            "BLZ",
            "CD8A stalk",
            "Notch1 TM",
            "Notch2 JMD",
            "Gal4 DBD",
            "VP64",
        ],
        "heavy_chain_snapgene_feature_evidence": {
            "coordinate_system": "0-based half-open nucleotide intervals",
            "leader": [8974, 9043],
            "myc": [9043, 9073],
            "mouse_igg1_ch1": [9076, 9367],
            "blz": [9379, 9469],
            "cd8a_stalk": [9469, 9550],
            "notch1_tm": [9550, 9616],
            "notch2_jmd": [9616, 9634],
            "gal4_dbd": [9634, 10075],
            "vp64": [10111, 10261],
        },
        "paratope_orientation": "extracellular",
        "membrane_anchor": "heavy_chain_C_terminal_Notch1_TM",
        "assembly": "extracellular ALZ-BLZ heterodimerization",
        "backbone_scope": (
            "The supplied backbone files do not contain the experiment-specific VH/VL inserts; "
            "the fixed profile therefore encodes display topology, not antibody-specific reach."
        ),
    },
    "design_policy": {
        "stage1_target_constructs": 6,
        "stage1_maximum_constructs": 6,
        "combination_minimum_sites": 2,
        "combination_maximum_sites": 4,
        "selection_basis": "geometry_only_surface_spread",
        "membrane_proximal_penalty": False,
        "epitope_prior_used_for_panel_selection": False,
    },
    "orientation_model": {
        "known": [
            "antigen ECD C terminus is membrane-proximal on K562",
            "Fab heavy-chain C terminus is membrane-proximal on Jurkat",
            "both binding surfaces occupy the extracellular intercellular space",
        ],
        "not_inferable_from_ecd_pdb": [
            "physical membrane plane",
            "intermembrane distance",
            "receptor tilt and rotational ensemble",
            "steric reach of a particular Fab-antigen pair",
        ],
        "pdb_only_surrogate": (
            "The last three observed C-terminal ECD CA atoms define an anchor point; the vector "
            "from that point to the ECD centroid defines only a relative extracellular axis."
        ),
    },
}

LEGACY_ASSAY_PROFILE_ID = JURKAT_FAB_K562_ANTIGEN_V1["id"]

# New V4.1 runs use an information-oriented panel objective. Keep V1 available
# for replay and preserve its original geometry-spread semantics.
JURKAT_FAB_K562_ANTIGEN_V2 = deepcopy(JURKAT_FAB_K562_ANTIGEN_V1)
JURKAT_FAB_K562_ANTIGEN_V2.update({
    "id": "jurkat_fab_k562_antigen_v2",
    "version": 2,
    "display_name": "Jurkat membrane-Fab x K562 membrane-antigen (hypothesis panel)",
})
JURKAT_FAB_K562_ANTIGEN_V2["design_policy"].update({
    "selection_basis": "hypothesis_identifiability_with_stability_constraints",
    "epitope_prior_role": "soft_hypothesis_prior",
    "geometry_role": "annotation_not_reporter_prediction",
    "reporter_is_affinity": False,
    "surface_tag_proves_native_fold": False,
})

ASSAY_PROFILES = {
    JURKAT_FAB_K562_ANTIGEN_V1["id"]: JURKAT_FAB_K562_ANTIGEN_V1,
    JURKAT_FAB_K562_ANTIGEN_V2["id"]: JURKAT_FAB_K562_ANTIGEN_V2,
}
DEFAULT_ASSAY_PROFILE_ID = JURKAT_FAB_K562_ANTIGEN_V2["id"]
SUPPORTED_ASSAY_PROFILE_IDS = tuple(ASSAY_PROFILES)


def get_assay_profile(profile_id: str | None) -> dict[str, Any] | None:
    """Return an isolated profile payload, or None for generic workflows."""
    if profile_id is None:
        return None
    if profile_id not in ASSAY_PROFILES:
        raise ValueError(f"Unknown assay profile: {profile_id}")
    return deepcopy(ASSAY_PROFILES[profile_id])
