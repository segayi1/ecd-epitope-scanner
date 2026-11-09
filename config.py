"""Central scientific thresholds for the epitope-mapping workflow.

Every number here is cited in docs/改进方案.md and echoed into
run_manifest.json / methods_report.md so that results stay auditable.
Changing a value changes what the output means; keep it in one place.
"""

from __future__ import annotations

import os
import sys
from copy import deepcopy
from typing import Any

# --- Structure eligibility gate (UniProt/PDB candidate audit) ---------------
MIN_ECD_COVERAGE = 0.95
MAX_INTERNAL_GAP_AA = 8
MIN_SEQUENCE_IDENTITY = 0.95
MIN_MAJOR_DOMAIN_COVERAGE = 0.90

# --- Surface screen ----------------------------------------------------------
RSASA_THRESHOLD = 0.25
SIDECHAIN_SASA_THRESHOLD_A2 = 20.0

# --- Patch partitioning (Fab-footprint scale) --------------------------------
PATCH_ENGINE_VERSION = "joint_v4_ml"

# These are versioned engineering defaults; they are not calibrated
# probabilities and must be reported in the
# run manifest. External model workers require separately installed dependencies.
JOINT_V4_ML_PARAMETERS = {
    "schema_version": "4.0",
    "mode": "strict",
    "max_sectors": 7,
    "max_primary_constructs": 7,
    "regional_recommendation_policy": "regional_mechanism_pareto_v1",
    "execution_protocol": {
        "name": "ordinary_site_first_v2",
        "spurs_batch_size": 64,
        "spurs_unique_request_limit": 1000,
        "foldx_slot_limit": 300,
        "foldx_initial_slots": 3,
        "foldx_append_slots": 7,
    },
    "use_hypothesis_hierarchy": True,
    "hypothesis_hierarchy": {
        "weighted_jaccard": 0.75, "maximum_center_geodesic_a": 6.0,
        "maximum_normal_angle_deg": 60.0, "maximum_nominal_geodesic_a": 30.0,
        "maximum_nominal_ca_a": 24.0,
        "core_weight": 1.0, "nominal_weight": 0.7, "extended_weight": 0.3,
    },
    "surface": {
        "rsasa_threshold": RSASA_THRESHOLD,
        "sidechain_sasa_threshold_a2": SIDECHAIN_SASA_THRESHOLD_A2,
        "minimum_mapped_area_a2": 8.0,
        "borderline_rsasa_threshold": 0.20,
    },
    "procaliper": {
        "enabled": True,
        "command": [],
        "backend": "local_worker",
        "python": sys.executable,
        "version": "0.2.1",
        "device": "cpu",
        "timeout_seconds": 180,
    },
    "spurs_single": {
        "model_assets_root": os.environ.get("ECD_MODEL_ASSETS_ROOT"),
        "enabled": True,
        "command": [],
        "backend": "local_worker",
        "device": "auto",
        "timeout_seconds": 900,
        "cuda_oom_retries": 1,
        "allow_cpu_fallback_on_oom": True,
        "preferred_max_ddg": 0.5,
        "acceptable_max_ddg": 1.0,
        "review_max_ddg": 2.0,
        "max_options_per_site": 3,
    },
    "discotope": {
        "model_assets_root": os.environ.get("ECD_MODEL_ASSETS_ROOT"),
        "enabled": True,
        "command": [],
        "backend": "local_worker",
        "device": "auto",
        "timeout_seconds": 900,
        "cuda_oom_retries": 1,
        "allow_cpu_fallback_on_oom": True,
        "seed_percentile": 0.80,
        "local_peak_percentile": 0.65,
        "neighborhood_radii_a": [9.0, 12.0, 15.0],
    },
    "sectoring": {
        "policy": "balanced_connected_v4",
        "target_sector_residues": 12,
        "preferred_minimum_sector_residues": 10,
        "preferred_maximum_sector_residues": 15,
        "maximum_geodesic_diameter_a": 30.0,
        "maximum_euclidean_diameter_a": 22.0,
        "minimum_sector_area_a2": 50.0,
        "minimum_sector_residues": 4,
    },
    "candidate_search": {
        "site_pool_per_sector": 12,
        "site_pool_per_hypothesis": 8,
        "options_per_site": 3,
        "beam_width": 300,
        "beam_width_by_order": {},
        "shortlist_per_sector": 30,
        "shortlist_per_order": 10,
        "shortlist_by_order": {},
        "global_shortlist_per_order": 60,
        "global_shortlist_by_order": {},
        "orders": [2, 3],
        "primary_orders": [2, 3],
        "exploratory_orders": [],
        "maximum_combination_order": 3,
        "foldx_local_budget": 60,
        "foldx_global_budget": 20,
        "foldx_local_exploration_budget": 16,
        "foldx_local_exploration_spurs_min_exclusive": 0.5,
        "foldx_local_exploration_spurs_max_inclusive": 1.0,
        "foldx_local_exploration_mechanism_percentile": 0.25,
        "foldx_first_pass_runs": 3,
        "foldx_second_pass_runs": 2,
        "foldx_second_pass_max_candidates": 30,
        "maximum_geodesic_span_a": 24.0,
        "epitope_priority_weight": 1.0,
        "epitope_site_pool_fraction": 0.5,
        "fallback_candidate_fraction": 0.5,
        "sector_resource_guarantees": True,
        "exclusive_groups_by_order": {"2": 4, "3": 2},
        "sector_coverage_search_trigger": 0.9,
    },
    "spurs_multi": {
        "model_assets_root": os.environ.get("ECD_MODEL_ASSETS_ROOT"),
        "enabled": True,
        "command": [],
        "backend": "local_worker",
        "device": "auto",
        "timeout_seconds": 1800,
        "cuda_oom_retries": 1,
        "allow_cpu_fallback_on_oom": True,
        "preferred_max_ddg": 0.5,
        "acceptable_max_ddg": 1.0,
        "review_max_ddg": 2.0,
        "absolute_gap_diagnostic_kcal_mol": 1.5,
    },
    "foldx_confirmation": {
        "enabled": True,
        "command": [],
        "backend": "local_worker",
        "device": "cpu",
        "binary": os.environ.get("ECD_FOLDX_BINARY", ""),
        "version": "5.1",
        "timeout_seconds": 7200,
        "number_of_runs": 5,
        "max_workers": 8,
        "combination_gate_kcal_mol": 2.0,
        "review_max_kcal_mol": 3.0,
        "max_candidates": 80,
    },
    "panel_optimizer": {
        "selection_priority": "direct_surface_coverage_then_geometric_complementarity_then_information_then_risk",
        "coverage_basis": "geodesic_0_eligible_surface_area",
        "solver": "deterministic_greedy",
        "time_limit_seconds": 60,
        "higher_order_max_primary": 2,
        "minimum_coverage_gain": 0.05,
        "hypothesis_response_radius_a": 15.0,
        "hypothesis_response_decay_a": 6.0,
        "risk_penalty_weights": {
            "stability": 0.04,
            "glycan": 0.05,
            "disulfide": 0.08,
            "native_interface": 0.03,
            "structure_uncertainty": 0.04,
            "redundancy": 0.03,
        },
        "recommendations_per_sector": 6,
        "epitope_priority_per_sector": 3,
        "stability_fallback_per_sector": 3,
    },
    "runtime": {"device": "auto", "cache_dir": "", "license_id": "research_only"},
}


def display_sectoring_parameters(
    chain_length: int,
    parameters: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return length-aware display-sector preferences without mutating defaults.

    Small analyzed chains use a smaller preferred display footprint. Explicit
    sectoring values override those defaults; hard geometry and minimum gates
    remain the same for every chain length.
    """
    if isinstance(chain_length, bool) or not isinstance(chain_length, int) or chain_length < 0:
        raise ValueError("chain_length must be a non-negative integer")

    result = deepcopy(JOINT_V4_ML_PARAMETERS["sectoring"])
    result.update(preferred_minimum_sector_count=4, preferred_maximum_sector_count=12)
    if chain_length < 120:
        result.update({
            "target_sector_residues": 7,
            "preferred_minimum_sector_residues": 6,
            "preferred_maximum_sector_residues": 8,
        })
    elif chain_length <= 150:
        result.update({
            "target_sector_residues": 8,
            "preferred_minimum_sector_residues": 6,
            "preferred_maximum_sector_residues": 10,
        })

    overrides = dict(parameters or {})
    result.update(overrides)

    # An explicit request to keep the historical 10+ residue footprint also
    # restores its established 10-15 preference interval unless the caller
    # supplied its own bounds.
    if "target_sector_residues" in overrides and int(overrides["target_sector_residues"]) >= 10:
        if "preferred_minimum_sector_residues" not in overrides:
            result["preferred_minimum_sector_residues"] = 10
        if "preferred_maximum_sector_residues" not in overrides:
            result["preferred_maximum_sector_residues"] = 15

    # A caller specifying an 8- or 9-residue target without explicit size
    # bounds should get a feasible preference interval around that target.
    if result.get("target_sector_residues") in {7, 8, 9}:
        target = int(result["target_sector_residues"])
        if "preferred_minimum_sector_residues" not in overrides:
            result["preferred_minimum_sector_residues"] = min(
                int(result["preferred_minimum_sector_residues"]), target
            )
        if "preferred_maximum_sector_residues" not in overrides:
            result["preferred_maximum_sector_residues"] = max(
                int(result["preferred_maximum_sector_residues"]), target
            )

    return result

JOINT_V4_ERROR_CODES = [
    "TOOL_NOT_CONFIGURED",
    "TOOL_ENVIRONMENT_INVALID",
    "MODEL_WEIGHTS_MISSING",
    "MODEL_REVISION_MISMATCH",
    "TOOL_TIMEOUT",
    "TOOL_CRASHED",
    "OUTPUT_MISSING",
    "OUTPUT_SCHEMA_INVALID",
    "RESIDUE_COUNT_MISMATCH",
    "WT_RESIDUE_MISMATCH",
    "CHAIN_MAPPING_AMBIGUOUS",
    "UNSUPPORTED_STRUCTURE",
    "LICENSE_RESTRICTED",
    "NO_SURFACE_CANDIDATES",
    "NO_STABLE_SINGLE_OPTIONS",
    "NO_FEASIBLE_COMBINATIONS",
    "NO_CONSENSUS_STABLE_COMBINATION",
]

# Joint V3 is deliberately opt-in until the benchmark set freezes its
# thresholds.  These are auditable startup parameters, not scientific
# acceptance cut-offs or calibrated epitope probabilities.
JOINT_V3_PARAMETERS = {
    "schema_version": "3.0",
    "max_probes": 6,
    "mutations_per_probe": [2, 4],
    "footprint_target_areas_a2": [600.0, 800.0, 1000.0],
    "max_footprint_centers": 64,
    "sector_maximum_geodesic_diameter_a": 35.0,
    "sector_maximum_euclidean_diameter_a": 30.0,
    "surface_rsasa_threshold": RSASA_THRESHOLD if "RSASA_THRESHOLD" in globals() else 0.25,
    "coverage_scenarios": [
        "direct_site",
        "geodesic_tight",
        "geodesic_medium",
        "geodesic_wide",
        "footprint_hit_1",
        "footprint_hit_2",
    ],
    "coverage_distances_a": {
        "geodesic_tight": 9.0,
        "geodesic_medium": 12.0,
        "geodesic_wide": 15.0,
    },
    "candidate_search": {
        "site_pool_cap_per_sector": 12,
        "mutation_options_cap_per_site": 3,
        "beam_width": 200,
        "shortlist_per_sector": 20,
        "maximum_geodesic_span_a": 24.0,
    },
    "solver": {
        "backend": "scipy_milp",
        "time_limit_seconds": 60,
        "mip_relative_gap": 0.01,
        "deterministic_seed": 0,
    },
    "stability": {
        "enabled": True,
        "foldx_runs": 5,
        "foldx_max_workers": 1,
        "foldx_max_candidates": 81,
        "combination_gate_kcal_mol": 2.0,
        "require_pass_for_final_panel": True,
    },
    "surface_graph": {"neighbor_cutoff_a": 10.0, "max_neighbors": 8},
}

# Surface V2 parameters. These affect geometry/footprint enumeration, not the
# production Shrake-Rupley residue SASA thresholds above.
SURFACE_V2_PARAMETERS = {
    "surface_geometry": {
        "backend": "sr_pointcloud",
        "probe_radius": 1.4,
        "sampling_points": 48,
        "point_graph_max_distance": 2.5,
        "normal_cosine_min": 0.3,
    },
    "footprint_centers": {
        "center_spacing": 7.0,
        "coverage_radius": 7.0,
        "minimum_component_area": 25.0,
    },
    "footprints": {
        "radii": [9.0, 12.0, 15.0],
        "target_area": 600.0,
        "min_area": 350.0,
        "max_area": 1000.0,
        "min_radius": 8.0,
        "max_radius": 16.0,
    },
    "residue_inclusion": {"min_surface_area": 8.0, "min_surface_fraction": 0.20},
    "deduplication": {"threshold": 0.75},
    "regions": {
        "resolution": 0.65,
        "maximum_center_geodesic": 24.0,
        "maximum_geodesic_diameter": 60.0,
        "domain_boundary_penalty": 0.15,
        "minimum_surface_area": 700.0,
        "minimum_shared_surface_area": 40.0,
        "minimum_weighted_jaccard": 0.12,
    },
}

# An antibody footprint buries ~600-900 A2 over 15-25 residues
# (Kringelum 2012).  13 A geodesic radius on the C-alpha surface graph
# approximates that scale on the targets this tool is used for.
PATCH_GEODESIC_RADIUS_A = 13.0
# Surface-geodesic distance wraps around thin sheets: a residue can be 12 A
# along the surface yet 40 A away in 3D (the "patch spans both faces" bug).
# Membership therefore also requires this EUCLIDEAN distance to the seed.
PATCH_EUCLIDEAN_CAP_A = 18.0
PATCH_MIN_RESIDUES = 8
# A patch is an antibody footprint, not an arbitrary connected component.
# These are pairwise hard limits checked after assignment/any merge.
PATCH_MAX_EUCLIDEAN_DIAMETER_A = 26.0
# Non-adjacent residues whose side chains point in strongly opposing
# directions are generally opposite faces of a beta-sandwich, not one surface.
PATCH_OPPOSING_FACE_COSINE_MAX = -0.80
PATCH_MIN_SITE_SPACING_A = 5.0  # key sites inside one patch
STAGE2_PERTURBATION_RADIUS_A = 10.0  # coverage-matrix neighbourhood

# --- Compact single-domain surface-face rescue --------------------------------
# These thresholds are deliberately separate from PATCH_*: they only apply
# after the normal footprint workflow has failed its own quality checks for a
# compact, single-domain ECD.  They must never retune legacy jobs.
COMPACT_ECD_MAX_LENGTH_AA = 180
COMPACT_DOMINANT_DOMAIN_MIN_FRACTION = 0.70
COMPACT_MIN_ENVELOPE_RESIDUES = 30
COMPACT_MIN_CORE_CONFIDENCE = 70.0
COMPACT_MUTABLE_ENVELOPE_RATIO_TRIGGER = 0.65
COMPACT_MAX_SURFACE_REGIONS = 3
COMPACT_GRAPH_CONTACT_A = 10.5
COMPACT_MIN_NORMAL_COSINE = -0.20
COMPACT_FACE_SPLIT_DIAMETER_A = 28.0
COMPACT_ENVELOPE_RSASA = 0.08
COMPACT_SURFACE_NORMAL_POINTS = 96
COMPACT_SURFACE_NORMAL_MIN_COHERENCE = 0.65
COMPACT_MAX_FOOTPRINTS_PER_FACE = 3
COMPACT_KEY_SITES_PER_FACE = 12
COMPACT_BETA_STRAND_CONTACT_A = 6.5
COMPACT_BETA_STRAND_MIN_CONTACTS = 3
COMPACT_BETA_SHEET_MIN_STRANDS = 2
COMPACT_BETA_SHEET_MIN_COVERAGE = 0.60
COMPACT_FOOTPRINT_GEODESIC_RADIUS_A = 13.0
COMPACT_FOOTPRINT_EUCLIDEAN_CAP_A = 18.0
COMPACT_COMBO_MAX_SPREAD_A = 15.0

# Stage-1 panel budget scales with ECD size but never exceeds 7 patches:
# each patch costs FoldX time linearly and a 7-construct panel already
# covers a typical single-pass ECD surface.
ADAPTIVE_REGIONS_DIVISOR = 22  # ~one patch per 22 ECD residues
ADAPTIVE_REGIONS_CAP = 7
ADAPTIVE_REGIONS_FLOOR = 3
SMALL_ECD_LENGTH_AA = 100
SMALL_ECD_SURFACE_RESIDUES = 40
SMALL_ECD_FINE_MAP_WINDOW_RESIDUES = 3


def adaptive_maximum_regions(ecd_length: int) -> int:
    """Panel budget from ECD size: ~1 patch per 22 residues, capped at 7.

    150 aa -> 7, 220 aa -> 7, 100 aa -> 5, 66 aa -> 3.  Explicit user
    values are still honoured when *smaller* (pipeline takes the min).
    """
    scaled = round(ecd_length / ADAPTIVE_REGIONS_DIVISOR)
    return max(ADAPTIVE_REGIONS_FLOOR, min(ADAPTIVE_REGIONS_CAP, scaled))


# Epitope-prior weights (sum of positives = 1.0; penalties subtracted).
PRIOR_W_EXPOSURE = 0.30  # mean rSASA of the patch
PRIOR_W_PROTRUSION = 0.20  # mean low-contact score
PRIOR_W_LOOP = 0.15  # loop/turn secondary-structure fraction
PRIOR_W_VARIABILITY = 0.15  # 1 - conservation (when available)
PRIOR_W_INTERFACE = 0.05  # native-interface fraction (small bonus)
PRIOR_PENALTY_GLYCAN = 0.10  # glycan-shielded fraction
# Cell-surface accessibility is now continuous.  The legacy value remains
# only for reading old manifests; new runs calculate a geometry-informed
# multiplicative penalty from MEMBRANE_ACCESSIBILITY_* below.
PRIOR_PENALTY_MEMBRANE_PROXIMAL_CELL = 0.50
MEMBRANE_ACCESSIBILITY_SCALE_A = 25.0
MEMBRANE_ACCESSIBILITY_MIN_FACTOR = 0.35
MEMBRANE_ACCESSIBILITY_REDUCED_THRESHOLD = 0.55
MEMBRANE_ACCESSIBILITY_GLYCAN_PENALTY = 0.20

# Tier assignment by prior rank among patches of one run.
TIER_HIGH_FRACTION = 0.4
TIER_MEDIUM_FRACTION = 0.35  # next band; remainder = LOW

# --- Mutation policy ----------------------------------------------------------
# WT classes drive the directed substitution menu.
BULKY_WT = set("WYFRKMIL")
SMALL_WT = set("GASTDN")
CHARGED_WT = set("KRHDE")
HOTSPOT_WT = set("WYF")  # Bogan & Thorn 1998 enrichment
ANCHOR_MUTANTS = {"E", "D"}  # charge-reversal targets from K/R
# A one-residue "hotspot removal" (W/Y/F -> A) is an anchor as well.

# Deterministic chemistry gates (stability proxy).  Hard protections
# (Cys/P/G introduction, native P/G, sequon, PTM, disulfide) are unchanged.
RISK_GATE_BURIED = 1.45  # rSASA < EDGE_RSASA boundary
RISK_GATE_EXPOSED = 1.90  # rSASA >= EXPOSED_RSASA boundary
EDGE_RSASA = 0.30
EXPOSED_RSASA = 0.35

# FoldX single-point gates by exposure.  Literature RMSE for FoldX ddG is
# 0.85-1.1 kcal/mol; the old flat +0.46 gate sat at ~1 sigma even for fully
# exposed residues whose folding risk is intrinsically low.
FOLDX_GATE_EXPOSED = 0.80  # rSASA >= EXPOSED_RSASA
FOLDX_GATE_EDGE = 0.60  # EDGE_RSASA <= rSASA < EXPOSED_RSASA
FOLDX_GATE_DEFAULT = 0.46  # tighter elsewhere (conservative)
FOLDX_COMBINATION_GATE = 2.0

# --- Flexible terminal segments ------------------------------------------------
# Long N/C-terminal runs with no observed secondary structure are the least
# trustworthy part of any structure (crystal truncation artefacts, one frame
# of an ensemble, low pLDDT in models) and are poor mapping instruments: the
# tail's conformation adapts, so a retained antibody binding cannot be read
# as "epitope not here".  They are excluded from patch design by default and
# reported separately as linear-epitope candidates (alanine-scan blocks).
EXCLUDE_FLEXIBLE_TERMINAL = True
FLEXIBLE_TERMINAL_MIN_RESIDUES = 6
# If fewer than this fraction of residues carry H/E annotation, the structure
# has no usable SSE record (e.g. a bare ATOM-only PDB): skip tail flagging.
FLEXIBLE_TERMINAL_MIN_SSE_FRACTION = 0.2

# Combination enumeration.
MAX_CANDIDATES_PER_SITE = 2
MAX_COMBINATIONS_PER_REGION = 512
ENGINEERED_NEIGHBORHOOD_RADIUS_A = 6.0
KEY_SITES_PER_PATCH = 5
MIN_COMBO_SITES = 3
MAX_COMBO_SITES = 6
MIN_COMBO_SPACING_A = 5.5
# A probe must fit inside one antibody footprint (~600-900 A2): max pairwise
# C-alpha distance between its mutations.  Without this cap the enumerator
# happily pairs residues on opposite faces of a domain.
MAX_COMBO_SPREAD_A = 22.0
DISRUPTION_LAMBDA_RISK = 0.5  # single-site menu ranking only; combo selection is Pareto
PARETO_KNEE_RISK_ANCHOR = 1.0  # knee = max D/(risk+anchor); anchor ~ 1x FoldX noise
MAX_PAIRWISE_CONFIRMS = 2  # stage-2 confirmatory pairs per probe

# --- FoldX execution -----------------------------------------------------------
FOLDX_SINGLE_SCREEN_RUNS = 1  # gate check only; median unnecessary
FOLDX_COMBINATION_RUNS = 5  # probes/alternates: robust median
FOLDX_STAGE2_RUNS = 1  # presets are screened; 5-run median on demand
FOLDX_MAX_WORKERS = 8  # parallel BuildModel chunk processes

# Disruption-score weights (dimensionless instrument-power score).
D_W_EXPOSURE = 0.9  # min(1.2, rsasa*2)
D_W_VOLUME_GAIN = 1.0  # max(0, dV)/85 per unit
D_W_CHARGE_FLIP = 0.8  # |dq| >= 1 (full reversal)
D_W_PARTIAL_CHARGE = 0.4  # 0 < |dq| < 1
D_W_HOTSPOT_REMOVAL = 0.7  # W/Y/F -> A
D_W_AROMATIC_INSERT = 0.5  # -> W/F/Y from non-aromatic

# --- Conservation --------------------------------------------------------------
CONSERVATION_MAX_ORTHOLOGS = 8
CONSERVATION_MIN_ORTHOLOGS = 3

# --- Legacy region semantics kept for backwards-compatible tooling -------------
LEGACY_REGION_COLORS = [
    "#087E8B",
    "#E4572E",
    "#6A4C93",
    "#2E7D32",
    "#C47F00",
    "#1565C0",
    "#AD1457",
    "#5D4037",
    "#00838F",
    "#8E24AA",
    "#6D4C41",
    "#558B2F",
]


def foldx_single_gate(rsasa: float) -> float:
    if rsasa >= EXPOSED_RSASA:
        return FOLDX_GATE_EXPOSED
    if rsasa >= EDGE_RSASA:
        return FOLDX_GATE_EDGE
    return FOLDX_GATE_DEFAULT


def risk_gate(rsasa: float) -> float:
    if rsasa >= EXPOSED_RSASA:
        return RISK_GATE_EXPOSED
    if rsasa >= EDGE_RSASA:
        return (RISK_GATE_BURIED + RISK_GATE_EXPOSED) / 2.0
    return RISK_GATE_BURIED
