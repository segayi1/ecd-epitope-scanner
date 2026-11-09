from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


@dataclass(frozen=True, order=True)
class ResidueId:
    model_id: str
    chain_id: str
    auth_seq_id: int
    insertion_code: str
    one_letter: str

    @property
    def label(self) -> str:
        return f"{self.chain_id}:{self.one_letter}{self.auth_seq_id}{self.insertion_code}"


@dataclass
class StructureAudit:
    selected_model: str
    selected_chains: list[str]
    chain_decision: str
    residue_mapping_status: str
    missing_segments: list[dict] = field(default_factory=list)
    nonstandard_residues: list[dict] = field(default_factory=list)
    duplicate_chain_groups: list[list[str]] = field(default_factory=list)
    quality_warnings: list[str] = field(default_factory=list)
    fatal_errors: list[str] = field(default_factory=list)
    model_count: int = 1
    available_chains: list[str] = field(default_factory=list)
    low_information_regions: list[dict] = field(default_factory=list)
    alternate_location_count: int = 0
    partial_occupancy_atom_count: int = 0


@dataclass
class TopologyContext:
    mode: Literal["topology_only", "tm_geometry"]
    anchor_residue: ResidueId | None
    anchor_point_xyz: tuple[float, float, float] | None
    relative_axis_xyz: tuple[float, float, float] | None
    membrane_plane: dict | None
    uncertainty: Literal["low", "medium", "high", "not_estimable"]
    assumptions: list[str] = field(default_factory=list)


@dataclass
class SurfaceResidueUnit:
    residue: ResidueId
    position: int
    rasa: float
    mapped_area_a2: float
    designable: bool
    exclusion_reasons: list[str]
    component_id: int = 0
    accessibility_confidence: float = 1.0
    # Optional C-alpha coordinate used to enforce compact 3-D sectors.
    coordinate_xyz: tuple[float, float, float] | None = None


@dataclass
class MicroFootprint:
    footprint_id: str
    center: ResidueId
    scale_id: str
    target_area_a2: float
    residue_ids: list[ResidueId]
    realized_area_a2: float
    boundary_truncated: bool = False
    quality_flags: list[str] = field(default_factory=list)


@dataclass
class ProbeSector:
    sector_id: str
    # Canonical mapped positions are used for optimization; each source
    # footprint still retains full ResidueId values for provenance.
    residue_ids: list[int]
    footprint_ids: list[str]
    total_area_a2: float
    designable_area_a2: float
    compactness: float
    quality_flags: list[str] = field(default_factory=list)


@dataclass
class MutationOption:
    residue: ResidueId
    position: int
    wild_type: str
    mutant_aa: str
    perturbation_class: str
    rule_risk: float
    chemical_perturbation_strength: float
    stability_value: float | None
    stability_unit: str | None
    evidence_level: Literal["EXTERNAL_ENGINE", "PROXY_ONLY"]
    rejection_reasons: list[str] = field(default_factory=list)

    @property
    def notation(self) -> str:
        return f"{self.wild_type}{self.position}{self.mutant_aa}"


@dataclass
class ProbeCandidate:
    candidate_id: str
    sector_id: str
    mutations: list[MutationOption]
    connected: bool
    coverage_by_scenario: dict[str, dict[str, float]] = field(default_factory=dict)
    coverage_units_by_scenario: dict[str, list[int]] = field(default_factory=dict)
    structural_risk: float = 0.0
    evidence_level: str = "PROXY_ONLY"
    warnings: list[str] = field(default_factory=list)
    stability_value: float | None = None
    stability_unit: str | None = None
    stability_status: str = "NOT_RUN"
    stability_provenance: dict = field(default_factory=dict)


@dataclass
class ProbePanel:
    panel_id: str
    selected_candidate_ids: list[str]
    objective_vector: list[float]
    achieved_coverage: dict[str, float]
    uncovered_components: list[dict]
    solver_status: str
    optimality_gap: float | None
    alternative_panel_ids: list[str] = field(default_factory=list)
    alternative_panels: list[dict] = field(default_factory=list)


@dataclass
class JointV3Result:
    audit: StructureAudit
    topology: TopologyContext
    surface_units: list[SurfaceResidueUnit]
    microfootprints: list[MicroFootprint]
    sectors: list[ProbeSector]
    candidates: list[ProbeCandidate]
    panel: ProbePanel
    legacy_regions: list[object]
    designs: list[object]
    partition_summary: dict
    design_summary: dict
    candidate_exclusions: dict[str, list[str]] = field(default_factory=dict)
    graph_provenance: dict = field(default_factory=dict)
    source_provenance: dict = field(default_factory=dict)
