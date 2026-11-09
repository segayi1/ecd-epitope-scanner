from __future__ import annotations

from copy import deepcopy
from typing import Any

# Display topology assumptions are explicit; users must verify their own constructs.
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
        "expression_control": (
            "N-terminal Myc reports tag accessibility/surface display; it does not by itself "
            "establish native ECD folding."
        ),
    },
    "fab_display": {
        "cell": "Jurkat",
        "light_chain_architecture_n_to_c": [
            "signal peptide", "FLAG", "VL cloning site", "mouse IgK constant", "ALZ"
        ],
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
        "paratope_orientation": "extracellular",
        "membrane_anchor": "heavy_chain_C_terminal_Notch1_TM",
        "assembly": "extracellular ALZ-BLZ heterodimerization",
        "backbone_scope": (
            "This profile describes display topology assumptions. Verify the topology and "
            "expression controls against your own constructs; antibody-specific reach is unknown."
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

JURKAT_FAB_K562_ANTIGEN_V3 = deepcopy(JURKAT_FAB_K562_ANTIGEN_V2)
JURKAT_FAB_K562_ANTIGEN_V3.update({
    "id": "jurkat_fab_k562_antigen_v3", "version": 3,
    "display_name": "Jurkat membrane-Fab x K562 membrane-antigen (adaptive 4-8)",
})
JURKAT_FAB_K562_ANTIGEN_V3["design_policy"].update({
    "stage1_minimum_constructs": 4, "stage1_target_constructs": 7,
    "stage1_maximum_constructs": 8, "adaptive_budget": True,
    "selection_basis": "hierarchical_hypothesis_identifiability_with_adaptive_budget",
    "minimum_is_target_not_forced_fill": True,
    "minimum_marginal_information_gain": 0.02,
    "minimum_relative_information_gain": 0.05,
})

ASSAY_PROFILES = {
    JURKAT_FAB_K562_ANTIGEN_V1["id"]: JURKAT_FAB_K562_ANTIGEN_V1,
    JURKAT_FAB_K562_ANTIGEN_V2["id"]: JURKAT_FAB_K562_ANTIGEN_V2,
    JURKAT_FAB_K562_ANTIGEN_V3["id"]: JURKAT_FAB_K562_ANTIGEN_V3,
}
DEFAULT_ASSAY_PROFILE_ID = JURKAT_FAB_K562_ANTIGEN_V3["id"]
SUPPORTED_ASSAY_PROFILE_IDS = tuple(ASSAY_PROFILES)


def get_assay_profile(profile_id: str | None) -> dict[str, Any] | None:
    """Return an isolated profile payload, or None for generic workflows."""
    if profile_id is None:
        return None
    if profile_id not in ASSAY_PROFILES:
        raise ValueError(f"Unknown assay profile: {profile_id}")
    return deepcopy(ASSAY_PROFILES[profile_id])
