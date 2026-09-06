"""Central scientific thresholds for the epitope-mapping workflow.

Every number here is cited in docs/改进方案.md and echoed into
run_manifest.json / methods_report.md so that results stay auditable.
Changing a value changes what the output means; keep it in one place.
"""
from __future__ import annotations

# --- Structure eligibility gate (UniProt/PDB candidate audit) ---------------
MIN_ECD_COVERAGE = 0.95
MAX_INTERNAL_GAP_AA = 8
MIN_SEQUENCE_IDENTITY = 0.95
MIN_MAJOR_DOMAIN_COVERAGE = 0.90

# --- Surface screen ----------------------------------------------------------
RSASA_THRESHOLD = 0.25
SIDECHAIN_SASA_THRESHOLD_A2 = 20.0

# --- Patch partitioning (Fab-footprint scale) --------------------------------
# An antibody footprint buries ~600-900 A2 over 15-25 residues
# (Kringelum 2012).  13 A geodesic radius on the C-alpha surface graph
# approximates that scale on the targets this tool is used for.
PATCH_GEODESIC_RADIUS_A = 13.0
# Surface-geodesic distance wraps around thin sheets: a residue can be 12 A
# along the surface yet 40 A away in 3D (the "patch spans both faces" bug).
# Membership therefore also requires this EUCLIDEAN distance to the seed.
PATCH_EUCLIDEAN_CAP_A = 18.0
PATCH_MIN_RESIDUES = 8
PATCH_MIN_SITE_SPACING_A = 5.0          # key sites inside one patch
STAGE2_PERTURBATION_RADIUS_A = 10.0     # coverage-matrix neighbourhood

# Stage-1 panel budget scales with ECD size but never exceeds 7 patches:
# each patch costs FoldX time linearly and a 7-construct panel already
# covers a typical single-pass ECD surface.
ADAPTIVE_REGIONS_DIVISOR = 22           # ~one patch per 22 ECD residues
ADAPTIVE_REGIONS_CAP = 7
ADAPTIVE_REGIONS_FLOOR = 3


def adaptive_maximum_regions(ecd_length: int) -> int:
    """Panel budget from ECD size: ~1 patch per 22 residues, capped at 7.

    150 aa -> 7, 220 aa -> 7, 100 aa -> 5, 66 aa -> 3.  Explicit user
    values are still honoured when *smaller* (pipeline takes the min).
    """
    scaled = round(ecd_length / ADAPTIVE_REGIONS_DIVISOR)
    return max(ADAPTIVE_REGIONS_FLOOR, min(ADAPTIVE_REGIONS_CAP, scaled))

# Epitope-prior weights (sum of positives = 1.0; penalties subtracted).
PRIOR_W_EXPOSURE = 0.30                 # mean rSASA of the patch
PRIOR_W_PROTRUSION = 0.20               # mean low-contact score
PRIOR_W_LOOP = 0.15                     # loop/turn secondary-structure fraction
PRIOR_W_VARIABILITY = 0.15              # 1 - conservation (when available)
PRIOR_W_INTERFACE = 0.05                # native-interface fraction (small bonus)
PRIOR_PENALTY_GLYCAN = 0.10             # glycan-shielded fraction
PRIOR_PENALTY_MEMRANE_PROXIMAL_CELL = 0.25  # cell_surface assay only
MEMBRANE_PROXIMAL_PRIOR_SCALE = 0.5     # multiplicative fallback for cell assay

# Tier assignment by prior rank among patches of one run.
TIER_HIGH_FRACTION = 0.4
TIER_MEDIUM_FRACTION = 0.35            # next band; remainder = LOW

# --- Mutation policy ----------------------------------------------------------
# WT classes drive the directed substitution menu.
BULKY_WT = set("WYFRKMIL")
SMALL_WT = set("GASTDN")
CHARGED_WT = set("KRHDE")
HOTSPOT_WT = set("WYF")                 # Bogan & Thorn 1998 enrichment
ANCHOR_MUTANTS = {"E", "D"}             # charge-reversal targets from K/R
# A one-residue "hotspot removal" (W/Y/F -> A) is an anchor as well.

# Deterministic chemistry gates (stability proxy).  Hard protections
# (Cys/P/G introduction, native P/G, sequon, PTM, disulfide) are unchanged.
RISK_GATE_BURIED = 1.45                 # rSASA < EDGE_RSASA boundary
RISK_GATE_EXPOSED = 1.90                # rSASA >= EXPOSED_RSASA boundary
EDGE_RSASA = 0.30
EXPOSED_RSASA = 0.35

# FoldX single-point gates by exposure.  Literature RMSE for FoldX ddG is
# 0.85-1.1 kcal/mol; the old flat +0.46 gate sat at ~1 sigma even for fully
# exposed residues whose folding risk is intrinsically low.
FOLDX_GATE_EXPOSED = 0.80               # rSASA >= EXPOSED_RSASA
FOLDX_GATE_EDGE = 0.60                  # EDGE_RSASA <= rSASA < EXPOSED_RSASA
FOLDX_GATE_DEFAULT = 0.46               # tighter elsewhere (conservative)
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
MAX_CANDIDATES_PER_SITE = 3
KEY_SITES_PER_PATCH = 5
MIN_COMBO_SITES = 3
MAX_COMBO_SITES = 6
MIN_COMBO_SPACING_A = 5.5
# A probe must fit inside one antibody footprint (~600-900 A2): max pairwise
# C-alpha distance between its mutations.  Without this cap the enumerator
# happily pairs residues on opposite faces of a domain.
MAX_COMBO_SPREAD_A = 22.0
DISRUPTION_LAMBDA_RISK = 0.5            # single-site menu ranking only; combo selection is Pareto
PARETO_KNEE_RISK_ANCHOR = 1.0           # knee = max D/(risk+anchor); anchor ~ 1x FoldX noise
MAX_PAIRWISE_CONFIRMS = 2               # stage-2 confirmatory pairs per probe

# --- FoldX execution -----------------------------------------------------------
FOLDX_SINGLE_SCREEN_RUNS = 1            # gate check only; median unnecessary
FOLDX_COMBINATION_RUNS = 5              # probes/alternates: robust median
FOLDX_STAGE2_RUNS = 1                   # presets are screened; 5-run median on demand
FOLDX_MAX_WORKERS = 8                   # parallel BuildModel chunk processes

# Disruption-score weights (dimensionless instrument-power score).
D_W_EXPOSURE = 0.9                      # min(1.2, rsasa*2)
D_W_VOLUME_GAIN = 1.0                   # max(0, dV)/85 per unit
D_W_CHARGE_FLIP = 0.8                   # |dq| >= 1 (full reversal)
D_W_PARTIAL_CHARGE = 0.4                # 0 < |dq| < 1
D_W_HOTSPOT_REMOVAL = 0.7               # W/Y/F -> A
D_W_AROMATIC_INSERT = 0.5               # -> W/F/Y from non-aromatic

# --- Conservation --------------------------------------------------------------
CONSERVATION_MAX_ORTHOLOGS = 8
CONSERVATION_MIN_ORTHOLOGS = 3

# --- Legacy region semantics kept for backwards-compatible tooling -------------
LEGACY_REGION_COLORS = [
    "#087E8B", "#E4572E", "#6A4C93", "#2E7D32",
    "#C47F00", "#1565C0", "#AD1457", "#5D4037",
    "#00838F", "#8E24AA", "#6D4C41", "#558B2F",
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
