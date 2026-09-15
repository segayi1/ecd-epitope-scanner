"""Unit tests for the epitope-mapping redesign: patch partitioning,
directed mutation menu, anchor policy, exhaustive enumeration and the
conservation adapter.  Everything here is synthetic and offline.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from ecd_epitope_scanner import config
from ecd_epitope_scanner.conservation import fetch_conservation
from ecd_epitope_scanner.design import (
    _combo_is_valid,
    _enumerate_combinations,
    saturation_scan,
    select_key_sites,
)
from ecd_epitope_scanner.models import FeatureRange, ProteinContext, RegionRecord, ResidueRecord
from ecd_epitope_scanner.partitioning import partition_surface
from ecd_epitope_scanner.stability import evaluate_substitution, is_anchor
from ecd_epitope_scanner.util import PipelineError


def make_record(
    position: int,
    aa: str,
    ca: tuple[float, float, float],
    rsasa: float = 0.5,
    contacts: int = 2,
) -> ResidueRecord:
    return ResidueRecord(
        uniprot_pos=position,
        aa=aa,
        chain_id="A",
        auth_seq_id=position,
        ca=ca,
        sasa=100.0,
        sidechain_sasa=40.0,
        rsasa=rsasa,
        contacts=contacts,
        secondary_structure="C",
        surface_accessible=True,
        folding_screen_pass=True,
    )


def make_protein(sequence: str) -> ProteinContext:
    return ProteinContext(
        accession="TEST",
        protein_name="Test protein",
        organism="test",
        sequence=sequence,
        ecd_start=1,
        ecd_end=len(sequence),
        membrane_anchor_position=len(sequence),
        topology="test",
    )


# --------------------------------------------------------------------------
# Directed mutation menu and anchors
# --------------------------------------------------------------------------

def test_menu_classes_follow_wt_chemistry():
    protein = make_protein("KWS")
    k = make_record(1, "K", (0.0, 0.0, 0.0))
    w = make_record(2, "W", (6.0, 0.0, 0.0), rsasa=0.4)
    s = make_record(3, "S", (12.0, 0.0, 0.0))

    assert evaluate_substitution(k, "E", protein).menu_class == "charge_reversal"
    assert evaluate_substitution(w, "A", protein).menu_class == "hotspot_to_ala"
    assert evaluate_substitution(s, "W", protein).menu_class == "trp_insert"
    assert is_anchor(k, "E")
    assert is_anchor(w, "A")
    assert not is_anchor(s, "W")  # Trp insertion is strong but not an anchor


def test_exposed_residues_get_a_relaxed_proxy_gate():
    protein = make_protein("KKK")
    exposed = make_record(1, "K", (0.0, 0.0, 0.0), rsasa=0.6)
    edge = make_record(2, "K", (6.0, 0.0, 0.0), rsasa=0.32)
    buried = make_record(3, "K", (12.0, 0.0, 0.0), rsasa=0.05)
    assert config.risk_gate(exposed.rsasa) == config.RISK_GATE_EXPOSED
    assert config.risk_gate(edge.rsasa) > config.RISK_GATE_BURIED
    assert config.risk_gate(buried.rsasa) == config.RISK_GATE_BURIED
    assert config.foldx_single_gate(0.6) == config.FOLDX_GATE_EXPOSED
    assert config.foldx_single_gate(0.3) == config.FOLDX_GATE_EDGE
    assert config.foldx_single_gate(0.1) == config.FOLDX_GATE_DEFAULT


def test_saturation_scan_covers_19_substitutions_with_menu_classes():
    protein = make_protein("KKKKK")
    record = make_record(2, "K", (0.0, 0.0, 0.0))
    choices = saturation_scan(record, protein)
    assert len(choices) == 19
    assert all(choice.menu_class for choice in choices)
    assert any(choice.menu_class == "charge_reversal" for choice in choices)


# --------------------------------------------------------------------------
# Combination enumeration
# --------------------------------------------------------------------------

def _site_setup() -> tuple[dict[int, ResidueRecord], dict[int, list]]:
    protein = make_protein("KKKKKK")
    records = {
        1: make_record(1, "K", (0.0, 0.0, 0.0)),
        2: make_record(2, "K", (10.0, 0.0, 0.0)),
        3: make_record(3, "K", (20.0, 0.0, 0.0)),
    }
    scans = {position: saturation_scan(record, protein) for position, record in records.items()}
    return records, scans


def test_enumeration_requires_anchor_and_spacing():
    records, scans = _site_setup()
    front, downgrades = _enumerate_combinations(
        [records[1], records[2], records[3]], scans, records, measured=False
    )
    assert front, "charge-reversal anchors are available, a full probe must exist"
    best = front[0][2]
    assert {choice.position for choice in best} == {1, 2, 3}
    assert any(is_anchor(records[choice.position], choice.mutant) for choice in best)
    assert not downgrades

    valid, reason = _combo_is_valid(best, records, measured=False, require_anchor=True)
    assert valid


def test_enumeration_downgrades_gracefully_without_anchors():
    records, scans = _site_setup()
    # Strip every anchor-capable substitution from the menus.
    for position, choices in scans.items():
        for choice in choices:
            if is_anchor(records[position], choice.mutant):
                choice.allowed = False
    ranked, downgrades = _enumerate_combinations(
        [records[1], records[2], records[3]], scans, records, measured=False
    )
    ranked = [mutations for _, _, mutations in ranked]
    if ranked:
        assert any("anchor" in note for note in downgrades)
    # Either way the workflow must not raise; empty means NO_LOW_RISK_COMBINATION.


def test_combo_rejects_close_spacing():
    records, scans = _site_setup()
    records[2] = make_record(2, "K", (3.0, 0.0, 0.0))  # 3 A from site 1
    scans[2] = saturation_scan(records[2], make_protein("KKK"))
    picks = [scans[1][0], scans[2][0], scans[3][0]]
    valid, reason = _combo_is_valid(picks, records, measured=False, require_anchor=False)
    assert not valid and reason == "spacing"


def test_combo_rejects_footprint_overflow():
    """Mutations spread across opposite faces of a domain can never perturb
    one paratope together -- max pairwise spread caps the probe."""
    protein = make_protein("KKKK")
    records = {
        1: make_record(1, "K", (0.0, 0.0, 0.0)),
        2: make_record(2, "K", (10.0, 0.0, 0.0)),
        3: make_record(3, "K", (20.0, 0.0, 0.0)),
        4: make_record(4, "K", (30.0, 0.0, 0.0)),
    }
    scans = {position: saturation_scan(record, protein) for position, record in records.items()}
    picks = [scans[position][0] for position in (1, 2, 3, 4)]
    valid, reason = _combo_is_valid(picks, records, measured=False, require_anchor=False)
    assert not valid and reason == "spread"
    valid, _ = _combo_is_valid(picks[:3], records, measured=False, require_anchor=False)
    assert valid, "3 sites within one footprint (spread 20 A) must stay valid"


def test_assignment_never_stretches_a_patch_when_no_local_seed_exists():
    """Surface-geodesic distance wraps around thin sheets: a residue 13 A
    along the graph but 50 A through space must not join the far seed."""
    import networkx as nx

    from ecd_epitope_scanner.partitioning import _assign_to_patches

    records = [
        make_record(1, "K", (0.0, 0.0, 0.0)),
        make_record(2, "K", (8.0, 0.0, 0.0)),
        make_record(3, "K", (50.0, 0.0, 0.0)),
        make_record(4, "K", (72.0, 0.0, 0.0)),
        make_record(5, "K", (80.0, 0.0, 0.0)),
    ]
    graph = nx.Graph()
    for record in records:
        graph.add_node(record.uniprot_pos)
    graph.add_edge(1, 2, weight=8.0, euclidean=8.0)
    graph.add_edge(2, 3, weight=5.0, euclidean=42.0)   # geodesically close, spatially far
    graph.add_edge(3, 4, weight=8.0, euclidean=22.0)
    graph.add_edge(4, 5, weight=8.0, euclidean=8.0)
    assignments = _assign_to_patches(graph, records, [1, 4])
    assert sorted(assignments[1]) == [1, 2]
    assert sorted(assignments[4]) == [4, 5]
    assert assignments[3] == [3], "a stray residue must seed its own patch rather than stretch another patch"


def test_key_sites_respect_minimum_spacing():
    records = [
        make_record(1, "K", (0.0, 0.0, 0.0)),
        make_record(2, "K", (3.0, 0.0, 0.0)),
        make_record(3, "K", (6.0, 0.0, 0.0)),
        make_record(4, "K", (20.0, 0.0, 0.0)),
    ]
    selected = select_key_sites(records, limit=4)
    positions = [record.uniprot_pos for record in selected]
    # The spaced sites must be chosen first; the legacy fill-to-limit fallback
    # may only append extras after them.
    assert positions[:3] == [1, 3, 4]


# --------------------------------------------------------------------------
# Patch partitioning
# --------------------------------------------------------------------------

def _cluster(offset: float, start: int, count: int, aa: str = "K") -> list[ResidueRecord]:
    return [
        make_record(
            start + index,
            aa,
            (offset + 2.0 * index, 0.0, 0.0),
        )
        for index in range(count)
    ]


def test_patches_cover_a_compact_surface_with_one_patch():
    protein = make_protein("K" * 9)
    records = _cluster(0.0, 1, 9)
    regions, summary = partition_surface(records, protein, maximum_regions=8)
    assert summary["coverage"] == 1.0
    assert len(regions) == 1
    assert regions[0].epitope_prior is not None
    assert regions[0].tier in {"HIGH", "MEDIUM", "LOW"}
    assert regions[0].residues == list(range(1, 10))


def test_distant_surfaces_form_separate_patches():
    protein = make_protein("K" * 18)
    records = _cluster(0.0, 1, 9) + _cluster(120.0, 10, 9)
    regions, summary = partition_surface(records, protein, maximum_regions=8)
    assert len(regions) == 2
    assert summary["coverage"] == 1.0
    assert {tuple(region.residues) for region in regions} == {
        tuple(range(1, 10)),
        tuple(range(10, 19)),
    }


def test_compact_single_domain_rescue_uses_continuous_envelope_not_sparse_mutants():
    """Protected surface nodes must preserve compact-domain geometry.

    The legacy footprint route intentionally sees many local patches on this
    circular surface.  The rescue router is allowed to activate only because
    the ECD is a compact, dominant single domain and the old result is
    pathological; it returns at most three continuous surface faces.
    """
    import math

    protein = make_protein("K" * 60)
    records = []
    for position in range(1, 61):
        theta = 2.0 * math.pi * (position - 1) / 60.0
        record = make_record(position, "K", (20.0 * math.cos(theta), 20.0 * math.sin(theta), 0.0))
        record.domain = "Ig-like V-type"
        record.plddt = 90.0
        record.surface_normal = (math.cos(theta), math.sin(theta), 0.0)
        # Every fifth residue is protected but remains part of the surface
        # envelope.  Removing it before partitioning would create islands.
        if position % 5 == 0:
            record.folding_screen_pass = False
        records.append(record)
    regions, summary = partition_surface(records, protein, maximum_regions=7)
    assert summary["partition_mode"] == "compact_surface_face"
    assert 1 <= len(regions) <= config.COMPACT_MAX_SURFACE_REGIONS
    assert all(region.surface_connectivity == "connected" for region in regions)
    assert sum(len(region.envelope_residues) for region in regions) == len(records)
    assert any(region.protected_surface_residues for region in regions)
    assert all(
        set(region.footprint_residues).issubset(set(region.mutable_residues))
        for region in regions
    )


def test_parallel_beta_strands_with_one_surface_normal_stay_one_face():
    protein = make_protein("K" * 36)
    records = []
    position = 1
    for strand_y in (-4.0, 0.0, 4.0):
        for x in range(12):
            record = make_record(position, "K", (float(x * 3), strand_y, 0.0))
            record.domain = "Ig-like V-type"
            record.plddt = 90.0
            record.surface_normal = (0.0, 0.0, 1.0)
            records.append(record)
            position += 1
    # Call the V2 implementation directly: this geometry represents one side
    # of a beta sheet and must not be split into middle-versus-edge strands.
    from ecd_epitope_scanner.partitioning import _partition_compact_surface_face

    output = _partition_compact_surface_face(
        records, protein, None, False, "cell_surface", None, None, records,
        ["synthetic beta-sheet regression"],
    )
    assert output is not None
    regions, summary = output
    assert summary["selected_face_count"] == 1
    assert len(regions) == 1
    assert set(regions[0].envelope_residues) == set(range(1, 37))


def test_beta_sheet_backbone_topology_keeps_parallel_strands_together():
    from ecd_epitope_scanner.partitioning import _beta_sheet_assignments

    protein = make_protein("K" * 64)
    segment_starts = [1, 11, 21, 31, 41, 51]
    protein.features = [
        FeatureRange("Beta strand", start, start + 3) for start in segment_starts
    ]
    records = []
    for sheet_index, starts in enumerate((segment_starts[:3], segment_starts[3:])):
        for strand_index, start in enumerate(starts):
            for offset in range(4):
                record = make_record(
                    start + offset,
                    "K",
                    (float(offset * 3.5), float(strand_index * 4.5), float(sheet_index * 12.0)),
                )
                record.surface_normal = (0.0, 0.0, 1.0 if sheet_index == 0 else -1.0)
                records.append(record)
    assignments = _beta_sheet_assignments(records, protein, records)
    assert assignments is not None
    groups = [set(members) for members in assignments.values()]
    assert len(groups) == 2
    assert {1, 2, 3, 4, 11, 12, 13, 14, 21, 22, 23, 24} in groups
    assert {31, 32, 33, 34, 41, 42, 43, 44, 51, 52, 53, 54} in groups


def test_compact_face_emits_safe_single_site_stage1_fallback(tmp_path):
    from ecd_epitope_scanner.design import design_regions
    from ecd_epitope_scanner.stability import EngineInventory

    protein = make_protein("K")
    record = make_record(1, "K", (0.0, 0.0, 0.0))
    region = RegionRecord(
        id=1, name="surface zone", domain="Ig", membrane_position="distal",
        structural_character="beta-sheet", residues=[1], surface_area=40.0,
        geodesic_diameter=0.0, interface_fraction=0.0, pae_median=None,
        auto_design_allowed=True, rationale="test", color="#000",
        partition_mode="compact_surface_face", mutable_residues=[1],
        footprint_residues=[1], candidate_footprints=[[1]],
    )
    designs, summary = design_regions(
        [region], [record], protein, EngineInventory(True, None, None), None, tmp_path
    )
    primary = next(item for item in designs if item.region_id == 1)
    assert primary.status == "RECOMMENDED_SINGLE_SITE_PROBE"
    assert len(primary.mutations) == 1
    assert summary["epitope_mapping_plan"]["design_downgrades"]


def test_compact_probe_has_a_stricter_spread_gate():
    records, _ = _site_setup()
    mutations = [
        evaluate_substitution(records[1], "E", make_protein("K" * 4)),
        evaluate_substitution(records[3], "E", make_protein("K" * 4)),
    ]
    # The synthetic endpoints are 20 A apart: legal in the legacy 22 A
    # footprint but not in the 15 A compact-domain probe footprint.
    assert _combo_is_valid(mutations, records, measured=False, require_anchor=False)[0]
    assert not _combo_is_valid(
        mutations, records, measured=False, require_anchor=False,
        max_combo_spread=config.COMPACT_COMBO_MAX_SPREAD_A,
    )[0]


def test_opposing_beta_sandwich_faces_are_not_connected_by_ca_proximity():
    from ecd_epitope_scanner.partitioning import _build_graph

    left = make_record(1, "K", (0.0, 0.0, 0.0))
    right = make_record(5, "K", (5.0, 0.0, 0.0))
    left.sidechain_direction = (1.0, 0.0, 0.0)
    right.sidechain_direction = (-1.0, 0.0, 0.0)
    graph = _build_graph([left, right])
    assert not graph.has_edge(1, 5), "opposing faces must not become one surface patch"


def test_patch_hard_diameter_limits_are_never_exceeded():
    protein = make_protein("K" * 12)
    records = [make_record(index, "K", (float(index * 6), 0.0, 0.0)) for index in range(1, 13)]
    regions, _ = partition_surface(records, protein, maximum_regions=3)
    by_position = {record.uniprot_pos: record for record in records}
    for region in regions:
        assert all(
            ((by_position[left].ca[0] - by_position[right].ca[0]) ** 2) ** 0.5
            <= config.PATCH_MAX_EUCLIDEAN_DIAMETER_A
            for left in region.residues for right in region.residues
        )


def test_cell_surface_assay_downweights_membrane_proximal_patches():
    protein = make_protein("K" * 18)
    records = _cluster(0.0, 1, 9) + _cluster(120.0, 10, 9)
    # One patch sits at the membrane anchor (last residue of the sequence).
    protein.membrane_anchor_position = 18
    regions_purified, purified_summary = partition_surface(
        records, protein, maximum_regions=8, assay_format="purified",
        membrane_anchor=(130.0, 0.0, 0.0), membrane_normal=(-1.0, 0.0, 0.0),
    )
    regions_cell, cell_summary = partition_surface(
        records, protein, maximum_regions=8, assay_format="cell_surface",
        membrane_anchor=(130.0, 0.0, 0.0), membrane_normal=(-1.0, 0.0, 0.0),
    )
    purified_priors = {row["patch_id"]: row["epitope_prior"] for row in purified_summary["coverage_plan"]}
    cell_priors = {row["patch_id"]: row["epitope_prior"] for row in cell_summary["coverage_plan"]}
    # The membrane-proximal patch must never gain priority under a
    # cell-surface assay.
    proximal_cell = min(regions_cell, key=lambda r: min(r.residues))
    proximal_purified = min(regions_purified, key=lambda r: min(r.residues))
    if proximal_cell.membrane_position == "membrane-proximal":
        assert cell_priors[proximal_cell.id] <= purified_priors[proximal_purified.id]


# --------------------------------------------------------------------------
# Conservation adapter
# --------------------------------------------------------------------------

class _StubClient:
    def __init__(self, payload):
        self._payload = payload

    def json(self, url: str):
        return self._payload


def test_conservation_scores_from_reviewed_orthologs(tmp_path: Path):
    protein = make_protein("AAAA")
    protein.raw = {"genes": [{"geneName": {"value": "TEST1"}}]}
    payload = {
        "results": [
            {"primaryAccession": f"P{i}", "sequence": {"value": seq}}
            for i, seq in enumerate(["AAAA", "AAAA", "AAAD"], start=1)
        ]
    }
    scores = fetch_conservation(protein, _StubClient(payload), tmp_path)
    assert scores is not None
    assert scores[1] == 1.0
    assert scores[4] == pytest.approx(2 / 3, abs=1e-3)


def test_conservation_requires_three_orthologs(tmp_path: Path):
    protein = make_protein("AAAA")
    protein.raw = {"genes": [{"geneName": {"value": "TEST1"}}]}
    payload = {"results": [{"primaryAccession": "P1", "sequence": {"value": "AAAA"}}]}
    assert fetch_conservation(protein, _StubClient(payload), tmp_path) is None
    assert "insufficient" in (tmp_path / "conservation.json").read_text(encoding="utf-8")


def test_conservation_without_gene_returns_none(tmp_path: Path):
    protein = make_protein("AAAA")
    assert fetch_conservation(protein, _StubClient({"results": []}), tmp_path) is None


def test_conservation_ignores_malformed_records_without_keyerror(tmp_path: Path):
    protein = make_protein("AAAA")
    protein.raw = {"genes": [{"geneName": {"value": "TEST1"}}]}
    payload = {
        "results": [{"primaryAccession": "BAD", "sequence": {}}]
        + [
            {"primaryAccession": f"P{i}", "sequence": {"value": "AAAA"}}
            for i in range(1, 4)
        ]
    }
    assert fetch_conservation(protein, _StubClient(payload), tmp_path) is not None
    assert '"status": "computed"' in (tmp_path / "conservation.json").read_text(encoding="utf-8")


def test_gene_symbol_resolution_is_human_and_unambiguous():
    from ecd_epitope_scanner.retrieval import resolve_uniprot_identifier

    class Resolver:
        def json(self, url: str):
            if "PECAM1.json" in url:
                raise PipelineError("not an accession")
            return {
                "results": [
                    {
                        "primaryAccession": "P16284",
                        "organism": {"taxonId": 9606},
                    }
                ]
            }

    assert resolve_uniprot_identifier("PECAM1", Resolver()) == ("P16284", "gene_symbol")


# --------------------------------------------------------------------------
# Performance budget: adaptive regions, batch parallelism, stage-2 screening
# --------------------------------------------------------------------------

def test_adaptive_maximum_regions_scales_with_ecd_size():
    from ecd_epitope_scanner.config import adaptive_maximum_regions

    assert adaptive_maximum_regions(150) == 7
    assert adaptive_maximum_regions(220) == 7      # cap: panel budget < 8
    assert adaptive_maximum_regions(100) == 5
    assert adaptive_maximum_regions(66) == 3       # floor
    assert adaptive_maximum_regions(1000) == 7


def test_pareto_front_and_knee_selection():
    """Front keeps every non-dominated trade-off; knee picks the balance point."""
    from ecd_epitope_scanner.design import _pareto_front, _pareto_knee

    conservative = (0.2, 8.0, "conservative")
    mid = (0.5, 12.0, "mid")
    aggressive = (1.8, 18.0, "aggressive")
    dominated = (0.6, 6.0, "dominated")
    front = _pareto_front([dominated, mid, conservative, aggressive])
    assert [item[2] for item in front] == ["conservative", "mid", "aggressive"], (
        "front must be ordered conservative -> aggressive and drop dominated points"
    )
    knee = _pareto_knee([(risk, disruption, marker) for risk, disruption, marker in front])
    assert knee == "mid", "knee = max D/(risk+1.0) must select the balance point"
    assert _pareto_knee([(0.2, 8.0, "solo")]) == "solo", "single-point front degenerates gracefully"


def test_enumeration_returns_pareto_front_ordered_by_risk():
    records, scans = _site_setup()
    front, downgrades = _enumerate_combinations(
        [records[1], records[2], records[3]], scans, records, measured=False
    )
    assert front, "anchor-capable menus must yield a front"
    risks = [risk for risk, _, _ in front]
    assert risks == sorted(risks), "front must be ordered conservative -> aggressive"
    assert all(
        any(is_anchor(records[choice.position], choice.mutant) for choice in combo)
        for _, _, combo in front
    )


# --------------------------------------------------------------------------
# Flexible terminal segments (no observed secondary structure)
# --------------------------------------------------------------------------

def _record_with_sse(position, ca, sse):
    record = make_record(position, "K", ca)
    record.secondary_structure = sse
    return record


def test_flexible_terminal_detection_flags_unstructured_ends_only():
    from ecd_epitope_scanner.structure_analysis import flag_flexible_terminal_segments

    records = (
        [_record_with_sse(i, (6.0 * (i - 1), 0.0, 0.0), "C") for i in range(1, 8)]    # 1-7: C-tail (7 res)
        + [_record_with_sse(i, (6.0 * (i - 1), 5.0, 0.0), "E") for i in range(8, 13)]  # 8-12: beta core
    )
    ranges = flag_flexible_terminal_segments(records)
    assert ranges == [(1, 7)]
    flagged = {record.uniprot_pos for record in records if record.flexible_terminal}
    assert flagged == set(range(1, 8))
    assert all("flexible-terminal" in record.risk_flags for record in records if record.flexible_terminal)
    # Short structured termini must never be flagged (walk stops at first SSE).
    assert not any(record.flexible_terminal for record in records if record.uniprot_pos >= 8)


def test_flexible_terminal_residues_are_excluded_from_patches():
    protein = make_protein("K" * 12)
    records = (
        [_record_with_sse(i, (6.0 * (i - 1), 0.0, 0.0), "C") for i in range(1, 8)]
        + [_record_with_sse(i, (6.0 * (i - 1), 5.0, 0.0), "E") for i in range(8, 13)]
    )
    from ecd_epitope_scanner.structure_analysis import flag_flexible_terminal_segments

    flag_flexible_terminal_segments(records)
    for record in records:
        record.surface_accessible = True
        record.folding_screen_pass = True
    regions, summary = partition_surface(records, protein, maximum_regions=6)
    in_patches = {position for region in regions for position in region.residues}
    assert in_patches == set(range(8, 13)), "flagged tail residues must not enter any patch"
    assert summary["excluded_flexible_terminal_residues"] == 7
    assert any("flexible terminal" in warning for warning in summary["warnings"])


def test_membrane_accessibility_is_continuous_not_binary_exclusion():
    protein = make_protein("K" * 12)
    protein.membrane_anchor_position = 12
    records = [
        make_record(index, "K", (float(index * 5), 0.0, 0.0))
        for index in range(1, 13)
    ]
    regions, summary = partition_surface(
        records, protein, maximum_regions=3, assay_format="cell_surface",
        membrane_anchor=(60.0, 0.0, 0.0), membrane_normal=(-1.0, 0.0, 0.0),
    )
    assert summary["membrane_accessibility_model"] != "unknown"
    assert all(region.membrane_accessibility_score is not None for region in regions)
    assert all(region.predicted_accessibility in {"normal", "reduced"} for region in regions)


def test_engineered_neighborhood_is_deprioritised_for_key_site_selection():
    regular = make_record(1, "K", (0.0, 0.0, 0.0))
    mismatch = make_record(2, "K", (20.0, 0.0, 0.0))
    mismatch.engineered_neighborhood = True
    selected = select_key_sites([mismatch, regular], limit=1)
    assert selected == [regular]


def test_split_groups_makes_near_equal_chunks():
    from ecd_epitope_scanner.stability import _split_groups

    chunks = _split_groups(list(range(7)), 3)
    assert [len(chunk) for chunk in chunks] == [3, 2, 2]
    flattened = [item for chunk in chunks for item in chunk]
    assert flattened == list(range(7)), "order must be preserved for re-aggregation"


def test_stage2_presets_are_screened_at_one_run(tmp_path, monkeypatch):
    """Panel designs get the 5-run median; Stage-2 presets only runs=1."""
    import ecd_epitope_scanner.design as design_module
    from ecd_epitope_scanner.design import design_regions
    from ecd_epitope_scanner.models import RegionRecord
    from ecd_epitope_scanner.stability import EngineInventory

    calls: list[tuple[int, int]] = []

    def fake_score(binary, pdb, groups, workdir, runs=5, already_repaired=False, **kwargs):
        calls.append((len(groups), runs))
        return [-0.5] * len(groups)

    monkeypatch.setattr(design_module, "repair_with_foldx", lambda *a, **k: tmp_path / "repaired.pdb")
    monkeypatch.setattr(design_module, "score_many_with_foldx", fake_score)

    protein = make_protein("KKKKK")
    records = [
        make_record(1, "K", (0.0, 0.0, 0.0)),
        make_record(2, "K", (10.0, 0.0, 0.0)),
        make_record(3, "K", (20.0, 0.0, 0.0)),
        make_record(4, "K", (30.0, 0.0, 0.0)),
        make_record(5, "K", (40.0, 0.0, 0.0)),
    ]
    region = RegionRecord(
        id=1, name="patch 1", domain="ECD", membrane_position="middle ECD",
        structural_character="loop-rich surface", residues=[1, 2, 3, 4, 5],
        surface_area=100.0, geodesic_diameter=40.0, interface_fraction=0.0,
        pae_median=None, auto_design_allowed=True, rationale="test", color="#000000",
    )
    inventory = EngineInventory(proxy=True, evoef2_path=None, foldx_path="foldx.exe")
    designs, summary = design_regions([region], records, protein, inventory, "structure.pdb", tmp_path)

    # single-point screen (4 sites x top-2, runs=1) -> panel batch (1 probe, runs=5)
    # -> stage-2 batch (4 singles + <=2 pairs, runs=1)
    assert calls[0] == (config.MAX_CANDIDATES_PER_SITE * config.KEY_SITES_PER_PATCH, config.FOLDX_SINGLE_SCREEN_RUNS)
    assert (1, config.FOLDX_COMBINATION_RUNS) in calls
    stage2_calls = [call for call in calls if call != calls[0] and call != (1, config.FOLDX_COMBINATION_RUNS)]
    assert stage2_calls and all(runs == config.FOLDX_STAGE2_RUNS for _, runs in stage2_calls)

    for design in designs:
        if design.design_role in {"STAGE2_SINGLE", "STAGE2_PAIRWISE"} and design.score_origin == "foldx":
            assert any("runs=1" in warning for warning in design.warnings)
