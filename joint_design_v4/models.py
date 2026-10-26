from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True, order=True)
class ResidueKey:
    model_id: str
    chain_id: str
    auth_seq_id: int
    insertion_code: str = ""
    label_seq_id: int | None = None
    uniprot_pos: int | None = None
    one_letter: str = "X"

    @property
    def stable_id(self) -> str:
        insertion = self.insertion_code or ""
        return f"{self.model_id}:{self.chain_id}:{self.auth_seq_id}{insertion}"


@dataclass
class ToolInventory:
    name: str
    available: bool
    executable: str | None = None
    environment_type: str = "local_process"
    version: str | None = None
    git_revision: str | None = None
    model_revision: str | None = None
    weights_sha256: str | None = None
    device: str | None = None
    license_id: str | None = None
    failure_reason: str | None = None


@dataclass
class PredictionRequest:
    schema_version: str
    request_id: str
    structure_path: str
    chain_id: str
    structure_type: str = "unknown"
    sequence: str = ""
    mutations: list[list[str]] = field(default_factory=list)
    output_dir: str = ""
    device: str = "auto"
    parameters: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class PredictionResult:
    status: str
    records: list[dict[str, Any]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[dict[str, Any]] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    runtime_seconds: float = 0.0

    @property
    def ok(self) -> bool:
        return self.status in {"success", "partial"}


@dataclass
class SurfaceEvidence:
    residue: ResidueKey
    procaliper_sasa_a2: float | None = None
    procaliper_rsasa: float | None = None
    native_sasa_a2: float | None = None
    native_rsasa: float | None = None
    normalized_context_sasa_a2: float | None = None
    normalized_context_rsasa: float | None = None
    native_context_absolute_delta: float | None = None
    native_context_status: str = "NOT_COMPARED"
    native_context_classification_changed: bool | None = None
    procaliper_classification_changed: bool | None = None
    comparison_method_status: str = "NOT_COMPARED"
    sidechain_sasa_a2: float | None = None
    mapped_surface_area_a2: float = 0.0
    surface_component_id: int | None = None
    sasa_agreement_status: str = "NOT_COMPARED"
    provenance: dict[str, Any] = field(default_factory=dict)


@dataclass
class EpitopeEvidence:
    residue: ResidueKey
    raw_score: float | None = None
    calibrated_score: float | None = None
    surface_percentile: float | None = None
    local_peak: bool = False
    predicted_epitope_label: bool | None = None
    neighborhood_mean_by_radius: dict[str, float] = field(default_factory=dict)
    island_id: str | None = None
    status: str = "MISSING"
    provenance: dict[str, Any] = field(default_factory=dict)


@dataclass
class SingleMutationEvidence:
    residue: ResidueKey
    mutant_aa: str
    spurs_ddg_kcal_mol: float | None = None
    stability_band: str = "NOT_SCORED"
    chemistry_class: str = "unknown"
    hard_exclusions: list[str] = field(default_factory=list)
    soft_warnings: list[str] = field(default_factory=list)
    eligible: bool = False
    provenance: dict[str, Any] = field(default_factory=dict)
    review_required: bool = False
    chemical_descriptors: dict[str, Any] = field(default_factory=dict)
    pareto_layer: int | None = None
    menu_selected: bool = False
    selection_reason: str | None = None
    impact_evidence: dict[str, Any] = field(default_factory=dict)


@dataclass
class CombinationStabilityEvidence:
    mutations: tuple[str, ...]
    spurs_multi_ddg_kcal_mol: float | None = None
    foldx_ddg_kcal_mol: float | None = None
    order: int = 0
    epistasis_estimate: float | None = None
    extrapolation_level: str = "UNKNOWN"
    consensus_status: str = "NOT_SCORED"
    disagreement: float | None = None
    provenance: dict[str, Any] = field(default_factory=dict)


@dataclass
class V4MutationOption:
    residue: ResidueKey
    mutant_aa: str
    notation: str
    chemistry_class: str
    spurs_ddg_kcal_mol: float | None
    score: float  # compatibility exposure x chemistry proxy; components are authoritative
    eligible: bool
    hard_exclusions: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    epitope_prior: float = 0.0
    epitope_supported: bool = False
    accessibility_score: float = 0.0
    chemical_disruption_score: float = 0.0
    stability_band: str = "NOT_SCORED"
    glycan_risk: list[str] = field(default_factory=list)
    structural_risk_flags: list[str] = field(default_factory=list)
    review_required: bool = False

    impact_evidence: dict[str, Any] = field(default_factory=dict)


@dataclass
class V4Candidate:
    candidate_id: str
    sector_id: str
    mutations: list[V4MutationOption]
    connected: bool = True
    order: int = 0
    epitope_score: float = 0.0
    structural_risk: float = 0.0
    coverage_by_scenario: dict[str, float] = field(default_factory=dict)
    covered_positions_by_scenario: dict[str, list[int]] = field(default_factory=dict)
    stability: CombinationStabilityEvidence | None = None
    formal_eligible: bool = False
    formal_exclusion_reasons: list[str] = field(default_factory=list)
    role: str = "UNASSIGNED"
    warnings: list[str] = field(default_factory=list)
    epitope_supported_sites: int = 0
    recommendation_class: str = "STABILITY_FALLBACK"
    hypothesis_response: dict[str, float] = field(default_factory=dict)
    # Explicit V4.1 names; legacy sector/response fields remain serialized for
    # consumers of the V4.0 compatibility projection.
    anchor_region_id: str | None = None
    hypothesis_perturbation: dict[str, float] = field(default_factory=dict)
    risk_components: dict[str, float | None] = field(default_factory=dict)
    search_role: str = "PRIMARY_POOL"
    targeted_hypotheses: list[str] = field(default_factory=list)
    qualification: str = "FORMAL_POOL"
    triple_role: str | None = None
    parent_candidate_ids: list[str] = field(default_factory=list)
    control_for_candidate_ids: list[str] = field(default_factory=list)
    coupling_class: str = 'NOT_EVALUATED'
    positive_lower_bound_qualified: bool = False
    footprint_characterization_status: str = 'NOT_EVALUATED'
    geometry_support: dict[str, Any] = field(default_factory=dict)
    impact_evidence: dict[str, Any] = field(default_factory=dict)


@dataclass
class V4Sector:
    sector_id: str
    residue_positions: list[int]
    component_id: int | None
    designable_area_a2: float
    total_area_a2: float
    epitope_score: float
    center_position: int | None = None
    geodesic_diameter_a: float | None = None
    euclidean_diameter_a: float | None = None
    island_ids: list[str] = field(default_factory=list)
    sector_type: str = "EPITOPE_PRIOR"
    quality_flags: list[str] = field(default_factory=list)
    semantic_type: str = "display_region"


@dataclass
class V4SurfaceHypothesis:
    hypothesis_id: str
    residue_weights: dict[int, float]
    region_id: str | None = None
    source_footprint_id: str | None = None
    quality_flags: list[str] = field(default_factory=list)
    center_position: int | None = None
    generation_mode: str | None = None
    geodesic_radius_a: float | None = None
    geodesic_diameter_a: float | None = None
    realized_area_a2: float | None = None
    target_area_a2: float | None = None
    center_xyz: tuple[float, float, float] | None = None
    mean_normal: tuple[float, float, float] | None = None
    mutable_positions: list[int] = field(default_factory=list)
    protected_positions: list[int] = field(default_factory=list)
    normal_dispersion: float | None = None
    mean_plddt: float | None = None
    minimum_plddt: float | None = None
    low_confidence_surface_fraction: float | None = None
    surface_component_id: int | None = None
    core_positions: list[int] = field(default_factory=list)
    boundary_positions: list[int] = field(default_factory=list)
    membership_fraction: dict[int, float] = field(default_factory=dict)
    discotope_prior: float | None = None
    glycan_fraction: float | None = None
    glycan_evidence: dict[str, Any] = field(default_factory=dict)
    membrane_geometry: dict[str, Any] = field(default_factory=dict)
    structural_uncertainty: dict[str, Any] = field(default_factory=dict)


@dataclass
class V4Panel:
    selected_candidate_ids: list[str]
    candidates_by_role: dict[str, list[str]]
    achieved_coverage: dict[str, float]
    uncovered_components: list[dict[str, Any]]
    solver_status: str
    coverage_denominator_residues: int = 0
    covered_residues_by_scenario: dict[str, int] = field(default_factory=dict)
    coverage_denominator_area_a2: float = 0.0
    achieved_area_coverage: dict[str, float] = field(default_factory=dict)
    mean_hypothesis_response: float | None = None
    mean_pairwise_separation: float | None = None
    indistinguishable_hypothesis_pairs: int | None = None
    minimum_signature_distance: float | None = None
    unresolved_hypothesis_groups: list[list[str]] = field(default_factory=list)
    uncovered_hypotheses: list[str] = field(default_factory=list)
    hypothesis_signatures: dict[str, list[float]] = field(default_factory=dict)
    prior_weighted_mean_pairwise_separation: float | None = None
    prior_weighted_hypothesis_response: float | None = None
    marginal_information_gain: dict[str, float] = field(default_factory=dict)
    budget_decision: PanelBudgetDecision | None = None
    nearby_hypothesis_pairs: list[dict[str, Any]] = field(default_factory=list)
    identifiability_diagnostics: dict[str, Any] = field(default_factory=dict)
    interface_comparison: dict[str, Any] = field(default_factory=dict)
    hierarchy_resolution: dict[str, Any] = field(default_factory=dict)
    sector_recommendations: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class FootprintFamily:
    family_id: str
    center_vertex_id: int
    member_footprint_ids: list[str]
    core_positions: list[int]
    nominal_positions: list[int]
    extended_positions: list[int]
    residue_weights: dict[int, float]
    membership: dict[int, float]
    representative_footprint_id: str
    center_position: int | None = None
    component_id: int | None = None
    domains: list[str] = field(default_factory=list)
    reliability_states: list[str] = field(default_factory=list)
    mean_normal: tuple[float, float, float] | None = None
    structural_uncertainty: dict[str, Any] = field(default_factory=dict)


@dataclass
class InferenceHypothesis(V4SurfaceHypothesis):
    family_ids: list[str] = field(default_factory=list)
    nominal_positions: list[int] = field(default_factory=list)
    source_footprint_ids: list[str] = field(default_factory=list)


@dataclass
class SurfaceComplexitySummary:
    eligible_residues: int
    eligible_area_a2: float
    display_regions: int
    raw_footprints: int
    footprint_families: int
    inference_hypotheses: int
    surface_components: int
    structural_domains: int


@dataclass
class PanelBudgetDecision:
    minimum_constructs: int
    target_constructs: int
    maximum_constructs: int
    adaptive_budget: bool
    selected_constructs: int = 0
    stop_reason: str = "NO_FEASIBLE_CANDIDATES"
    budget_curve: list[dict[str, Any]] = field(default_factory=list)
    minimum_target_met: bool = False
    normalization: str = "fixed_sqrt_8"
    selection_priority: str = "coverage_then_information_then_risk"
    coverage_basis: str = "hypothesis_response_proxy"


@dataclass
class JointV4Result:
    schema_version: str
    recommendation_status: str
    tool_inventory: dict[str, ToolInventory]
    surface_evidence: list[SurfaceEvidence]
    epitope_evidence: list[EpitopeEvidence]
    mutation_options: list[SingleMutationEvidence]
    sectors: list[V4Sector]
    candidates: list[V4Candidate]
    panel: V4Panel
    surface_hypotheses: list[V4SurfaceHypothesis] = field(default_factory=list)
    legacy_regions: list[Any] = field(default_factory=list)
    legacy_designs: list[Any] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[dict[str, Any]] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    raw_hypotheses: list[V4SurfaceHypothesis] = field(default_factory=list)
    footprint_families: list[FootprintFamily] = field(default_factory=list)
    complexity: SurfaceComplexitySummary | None = None
