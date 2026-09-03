from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


TERMINAL_STATUSES = {"COMPLETE", "FAILED", "WAITING_FOR_AF3_UPLOAD"}


@dataclass(slots=True)
class FeatureRange:
    type: str
    start: int
    end: int
    description: str = ""
    evidence: str = ""

    def contains(self, position: int) -> bool:
        return self.start <= position <= self.end


@dataclass(slots=True)
class ProteinContext:
    accession: str
    protein_name: str
    organism: str
    sequence: str
    ecd_start: int
    ecd_end: int
    membrane_anchor_position: int
    topology: str
    reviewed: bool = False
    features: list[FeatureRange] = field(default_factory=list)
    pdb_ids: list[str] = field(default_factory=list)
    source_url: str = ""
    raw: dict[str, Any] = field(default_factory=dict, repr=False)

    @property
    def ecd_sequence(self) -> str:
        return self.sequence[self.ecd_start - 1 : self.ecd_end]

    @property
    def ecd_length(self) -> int:
        return self.ecd_end - self.ecd_start + 1

    def feature_at(self, position: int, *types: str) -> list[FeatureRange]:
        wanted = {item.lower() for item in types}
        return [
            item
            for item in self.features
            if item.contains(position) and (not wanted or item.type.lower() in wanted)
        ]


@dataclass(slots=True)
class StructureCandidate:
    identifier: str
    category: str
    provider: str
    model_url: str
    page_url: str = ""
    model_format: str = "MMCIF"
    method: str = ""
    resolution: float | None = None
    sequence_identity: float | None = None
    uniprot_start: int | None = None
    uniprot_end: int | None = None
    reported_coverage: float | None = None
    chain_ids: list[str] = field(default_factory=list)
    created: str = ""
    explicit_af3: bool = False
    observed_ecd_coverage: float | None = None
    max_internal_gap: int | None = None
    engineered_mismatches: int | None = None
    selected_chain: str | None = None
    eligible: bool = False
    rejection_reasons: list[str] = field(default_factory=list)
    rank: int | None = None
    local_path: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class ResidueRecord:
    uniprot_pos: int
    aa: str
    chain_id: str
    auth_seq_id: int
    insertion_code: str = ""
    label_chain_id: str = ""
    label_seq_id: int | None = None
    residue_name: str = ""
    ca: tuple[float, float, float] | None = None
    sasa: float = 0.0
    sidechain_sasa: float = 0.0
    rsasa: float = 0.0
    contacts: int = 0
    plddt: float | None = None
    secondary_structure: str = "C"
    domain: str = "ECD"
    interface: bool = False
    interface_partners: list[str] = field(default_factory=list)
    glycan_shielded: bool = False
    glycosylation_site: bool = False
    glycosylation_sequon: bool = False
    disulfide: bool = False
    ptm: bool = False
    salt_bridge: bool = False
    domain_interface: bool = False
    surface_accessible: bool = False
    folding_screen_pass: bool = False
    hard_exclusions: list[str] = field(default_factory=list)
    risk_flags: list[str] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)
    region_id: int | None = None
    token_index: int | None = None


@dataclass(slots=True)
class RegionRecord:
    id: int
    name: str
    domain: str
    membrane_position: str
    structural_character: str
    residues: list[int]
    surface_area: float
    geodesic_diameter: float
    interface_fraction: float
    pae_median: float | None
    auto_design_allowed: bool
    rationale: str
    color: str


@dataclass(slots=True)
class MutationChoice:
    position: int
    wt: str
    mutant: str
    stability_score: float
    perturbation_score: float
    stability_percentile: float
    allowed: bool
    reasons: list[str] = field(default_factory=list)
    foldx_ddg: float | None = None
    stability_class: str = "proxy"

    @property
    def notation(self) -> str:
        return f"{self.wt}{self.position}{self.mutant}"


@dataclass(slots=True)
class DesignRecord:
    id: str
    name: str
    region_id: int | None
    mutations: list[MutationChoice]
    status: str
    stability_engine: str
    combined_stability_score: float | None
    perturbation_score: float
    rationale: str
    warnings: list[str]
    mutant_ecd_sequence: str
    mutant_coding_sequence: str
    saturation_scan: dict[str, list[dict[str, Any]]] = field(default_factory=dict)
    custom: bool = False


def to_dict(value: Any) -> Any:
    if hasattr(value, "__dataclass_fields__"):
        return {key: to_dict(item) for key, item in asdict(value).items()}
    if isinstance(value, list):
        return [to_dict(item) for item in value]
    if isinstance(value, dict):
        return {key: to_dict(item) for key, item in value.items()}
    return value
