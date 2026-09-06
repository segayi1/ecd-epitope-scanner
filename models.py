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
    # Fraction of aligned reviewed orthologs carrying the same residue
    # (1.0 = fully conserved).  None when conservation was not computed.
    conservation_score: float | None = None
    # Terminal run with no observed secondary structure (docs: linear-epitope
    # candidate, excluded from patch design by default).
    flexible_terminal: bool = False


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
    # --- epitope-mapping patch extensions (backwards compatible) ---
    domain_id: str = "ECD"                  # dominant structural domain label
    epitope_prior: float | None = None      # dimensionless 0..~1
    prior_components: dict[str, float] = field(default_factory=dict)
    tier: str = "MEDIUM"                    # HIGH | MEDIUM | LOW
    assay_accessible: bool = True           # IgG reachability for the assay format
    required_probe_coverage: int = 1        # primary + alternates demanded by plan


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
    # --- epitope-mapping extensions ---
    menu_class: str = "other"
    # charge_reversal | hotspot_to_ala | bulky_to_ala | trp_insert |
    # conservative | other
    disruption_score: float = 0.0           # instrument-power score
    esm_llr: float | None = None            # PLM zero-shot log-likelihood ratio

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
    # Provenance of `combined_stability_score`.  MUST stay consistent with
    # `stability_engine`: a "FoldX 5.1" engine label is only legal when
    # score_origin == "foldx" (i.e. the number really is kcal/mol).
    score_origin: str = "proxy"          # foldx | evoef2 | proxy
    foldx_diagnostics: dict[str, Any] = field(default_factory=dict)
    # --- epitope-mapping extensions ---
    design_role: str = "STAGE1_PROBE"
    # STAGE1_PROBE | ALTERNATE_PROBE | STAGE2_SINGLE | STAGE2_PAIRWISE
    parent_design_id: str | None = None  # links stage-2 designs to their probe
    disruption_score: float | None = None    # probe instrument power (sum of D)
    anchor_notation: str | None = None       # the anchor mutation, if any


def to_dict(value: Any) -> Any:
    if hasattr(value, "__dataclass_fields__"):
        data = {key: to_dict(item) for key, item in asdict(value).items()}
        # `notation` is a derived property (WT+position+mutant); include it in
        # every serialized choice so JSON/CSV consumers never re-derive it.
        if isinstance(value, MutationChoice):
            data["notation"] = value.notation
        return data
    if isinstance(value, list):
        return [to_dict(item) for item in value]
    if isinstance(value, dict):
        return {key: to_dict(item) for key, item in value.items()}
    return value
