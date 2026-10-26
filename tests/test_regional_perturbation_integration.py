"""Integration checks for regional mechanism ranking and its worker boundary."""

from copy import deepcopy
import pytest

from ecd_epitope_scanner.joint_design_v4.engine import (
    _annotate_stability,
    _foldx_shortlist,
    _prediction,
)
from ecd_epitope_scanner.joint_design_v4.models import (
    CombinationStabilityEvidence,
    PredictionResult,
    ResidueKey,
    V4Candidate,
    V4MutationOption,
    V4Sector,
)
from ecd_epitope_scanner.joint_design_v4.regional_perturbation import (
    POLICY,
    dominates,
)
from ecd_epitope_scanner.joint_design_v4.sector_recommendations import build_sector_recommendations


AXES = (
    "charge_delta",
    "volume_delta_a3",
    "hydropathy_delta",
    "donor_delta",
    "acceptor_delta",
    "aromatic_delta",
)


def test_default_engine_policy_scores_unique_normalized_requests_and_exports(tmp_path, monkeypatch):
    import re
    from ecd_epitope_scanner.joint_design_v4 import engine
    from ecd_epitope_scanner.models import ProteinContext, ResidueRecord

    sequence = "DEKRFY"
    records = [ResidueRecord(uniprot_pos=37 + index, aa=aa, chain_id="A", auth_seq_id=37 + index,
        ca=(float(index % 3) * 4, float(index // 3) * 4, 0.0), sasa=100, sidechain_sasa=60,
        rsasa=.6, surface_accessible=True, sidechain_exposure={"value": .8, "status": "KNOWN",
            "numerator_a2": 60, "denominator_a2": 75}) for index, aa in enumerate(sequence)]
    protein = ProteinContext("TEST", "test", "test", "A" * 36 + sequence, 37, 42, 43, "test")
    supplied = {
        "spurs": PredictionResult("success", records=[{"position": index, "wt": aa,
            "scores": {mutant: .1 for mutant in "ACDEFGHIKLMNPQRSTVWY" if mutant != aa}}
            for index, aa in enumerate(sequence, 1)]),
        "discotope": PredictionResult("success", records=[{"position": index, "score": 0}
            for index in range(1, len(sequence) + 1)]),
        "procaliper": PredictionResult("failed"),
    }
    captured = {}
    original = engine._prediction

    def predict(adapter, structure, chain, prediction_sequence, output, structure_type,
                parameters, prediction_results, key=None, mutations=None):
        if key in {"spurs_multi", "foldx"}:
            captured.setdefault(key if key != "foldx" or not parameters.get("sampling_slot_start") else str(output), []).extend(mutations or [])
            return PredictionResult("success", records=[{"mutations": group, "ddg": .1,
                **({"paired_models": [{"sampling_slot": slot, "ddg_kcal_mol": .1,
                    "batch_path": str(output), "group": index, "run": slot}
                    for slot in range(int(parameters.get("sampling_slot_start", 0)),
                        int(parameters.get("sampling_slot_start", 0))+int(parameters.get("number_of_runs", 3)))]} if key == "foldx" else {})}
                for index, group in enumerate(mutations or [])], provenance={"test_request_driven": True})
        return original(adapter, structure, chain, prediction_sequence, output, structure_type,
                        parameters, prediction_results, key, mutations)

    monkeypatch.setattr(engine, "_prediction", predict)
    result = engine.run_joint_v4_ml(records, protein, tmp_path,
        parameters={"mode": "strict", "use_hypothesis_hierarchy": False,
                    "candidate_search": {"exploratory_orders": []}}, prediction_results=supplied)
    assert result.provenance["sector_recommendation_policy"]["version"] == POLICY
    assert result.provenance["sector_recommendation_policy"]["ranking_version"] == "mechanism_channels_v2"
    assert captured["spurs_multi"] and 0 < len(captured["foldx"]) <= 80
    for requests in captured.values():
        signatures = [tuple(sorted(group)) for group in requests]
        assert len(signatures) == len(set(signatures))
        assert all(1 <= int(re.search(r"\d+", notation).group()) <= len(sequence)
                   for group in requests for notation in group)
    assert all(candidate.covered_positions_by_scenario for candidate in result.candidates)
    assert result.panel.selected_candidate_ids
    assert all(row["formal_eligible"] and row["impact_evidence"]["policy"] == POLICY
               for row in result.panel.sector_recommendations)
    assert (tmp_path / "sector_recommendations_v4.csv").is_file()
    import csv, json
    rows = list(csv.DictReader((tmp_path / "sector_recommendations_v4.csv").open(encoding="utf-8")))
    assert any(json.loads(row["impact_evidence"]).get("local_interaction_diagnostics") for row in rows)
    assert all(row["resource_audit"]["single_threshold_audit"]["review_only"] for row in result.panel.sector_recommendations)


def _option(position, mutant="A", exposure=0.7, magnitude=1.0, epitope=False):
    residue = ResidueKey("0", "A", position, uniprot_pos=position, one_letter="D")
    descriptors = {axis: magnitude for axis in AXES}
    return V4MutationOption(
        residue, mutant, f"D{position}{mutant}", "chemistry", 0.0, 1.0, True,
        epitope_prior=0.95 if epitope else 0.0,
        epitope_supported=epitope,
        impact_evidence={
            "sidechain_exposure": {"value": exposure, "status": "KNOWN"},
            "chemical_descriptors": descriptors,
        },
    )


def _candidate(identifier, positions, *, order=None, exposure=0.7, magnitude=1.0,
               formal=True, epitope=False, selection_rank=0):
    mutations = [_option(position, exposure=exposure, magnitude=magnitude, epitope=epitope)
                 for position in positions]
    order = order or len(mutations)
    return V4Candidate(
        identifier,
        "VS001",
        mutations,
        order=order,
        epitope_score=0.95 if epitope else 0.0,
        stability=CombinationStabilityEvidence(
            tuple(option.notation for option in mutations), 0.1, 0.1, order,
            consensus_status="CONSENSUS_PASS",
        ),
        formal_eligible=formal,
        impact_evidence={"policy": POLICY, "selection_rank": selection_rank},
    )


def _regional_build(candidates, *, minimum_gain=0.05):
    sector = V4Sector("VS001", [1, 2, 3], 1, 30, 30, 0)
    graph = {position: {neighbor: 6.0 for neighbor in (position - 1, position + 1)
                        if 1 <= neighbor <= 3} for position in range(1, 4)}
    return build_sector_recommendations(
        candidates,
        [sector],
        {position: 10.0 for position in range(1, 4)},
        {"regional_recommendation_policy": POLICY,
         "panel_optimizer": {"minimum_coverage_gain": minimum_gain}},
        graph=graph,
    )


def test_cross_order_triple_requires_known_gain_and_raw_axis_dominance():
    double = _candidate("double", [1, 2], exposure=0.4, magnitude=1.0)
    triple = _candidate("triple", [1, 2, 3], exposure=0.9, magnitude=2.0)
    triple.risk_components.update(
        regional_gain_basis="mean_regional_area_coverage_9_12_15",
        regional_higher_order_gain=0.08,
    )

    # Cross-order comparison is valid when both candidates have the same known axes.
    from ecd_epitope_scanner.joint_design_v4.regional_perturbation import candidate_impact
    assert dominates(candidate_impact(triple), candidate_impact(double))
    assert _regional_build([double, triple])[0]["candidate_id"] == "triple"

    # A triple with unknown gain cannot replace a qualified double.
    triple.risk_components["regional_higher_order_gain"] = None
    assert _regional_build([double, triple])[0]["candidate_id"] == "double"

    # Known but subthreshold gain also cannot replace it.
    triple.risk_components["regional_higher_order_gain"] = 0.049
    assert _regional_build([double, triple])[0]["candidate_id"] == "double"


@pytest.mark.parametrize("cap", [4, 8])
def test_twelve_regions_keep_all_qualified_principles_outside_compact_panel(cap):
    from ecd_epitope_scanner.joint_design_v4.optimizer import select_panel
    from ecd_epitope_scanner.joint_design_v4.spatial_support import annotate_spatial_support
    sectors, pool, graph, areas = [], [], {}, {}
    for index in range(12):
        positions = list(range(index * 4 + 1, index * 4 + 5))
        sectors.append(V4Sector(f"VS{index + 1:03d}", positions, 1, 40, 40, 0))
        for position in positions:
            graph[position] = {neighbor: 6.0 for neighbor in (position - 1, position + 1) if neighbor in positions}
            areas[position] = 10.0
        strong = _candidate(f"{index:02d}_strong", positions[:2], magnitude=2.0, exposure=.9)
        low_risk = _candidate(f"{index:02d}_risk", positions[-2:], magnitude=.5, exposure=.6)
        low_risk.stability.spurs_multi_ddg_kcal_mol = -.2
        low_risk.stability.foldx_ddg_kcal_mol = -.2
        pool.extend([strong, low_risk])
    annotate_spatial_support(pool, graph, areas)
    _, selected = select_panel(pool, cap, sectors=sectors, surface_areas=areas,
        panel_budget={"minimum_constructs": 4, "target_constructs": cap, "maximum_constructs": cap, "adaptive_budget": False})
    assert len(selected) == cap
    rows = build_sector_recommendations(pool, sectors, areas,
        {"regional_recommendation_policy": POLICY}, selected, graph=graph)
    assert len(rows) == 12
    expected = {"mechanism_pareto", "sector_coverage", "lower_structural_risk", "charge_extreme",
                "volume_extreme", "hydropathy_extreme", "hbond_extreme"}
    for row in rows:
        items = [row, *row["alternatives"]]
        assert all(item["formal_eligible"] for item in items)
        assert {principle for item in items for principle in item["principles"]} >= expected
        assert len({item["candidate_id"] for item in items}) == len(items)
        assert len(row["all_formal_local_candidates"]) == 2
    assert any(not item["selected_in_panel"] for row in rows for item in [row, *row["alternatives"]])


def test_failed_principles_are_gaps_but_remain_available_for_confirmation_planning():
    failed = _candidate("failed", [1, 2], formal=False)
    failed.formal_exclusion_reasons = ["CONSENSUS_REVIEW"]
    row = _regional_build([failed])[0]
    assert row["status"] == "NEAR_THRESHOLD" and not row["formal_eligible"]
    assert not row["alternatives"] and not row["all_formal_local_candidates"]
    assert {item["principle"] for item in row["unavailable_principles"]} >= {
        "sector_coverage", "lower_structural_risk", "charge_extreme", "volume_extreme", "hydropathy_extreme", "hbond_extreme"}
    assert all(item["reason"] == "NO_QUALIFIED_PRINCIPLE_CANDIDATE" for item in row["unavailable_principles"])
    planned = build_sector_recommendations([failed], [V4Sector("VS001", [1, 2, 3], 1, 30, 30, 0)],
        {1: 10., 2: 10., 3: 10.}, {"regional_recommendation_policy": POLICY}, confirmation_planning=True)[0]
    assert {"lower_structural_risk", "charge_extreme", "volume_extreme"} <= set(planned["principles"])
    assert not planned["formal_eligible"] and not failed.formal_eligible


def test_missing_axes_are_not_imputed_for_cross_order_comparisons():
    known = {axis: 1.0 for axis in AXES}
    partial = dict(known)
    partial["aromatic_delta"] = None
    high = {"order": 3, "exposure_min": 0.9,
            "mean_absolute_delta": known, "missing_mask": [], "complete": True}
    low_partial = {"order": 2, "exposure_min": 0.5,
                   "mean_absolute_delta": partial, "missing_mask": ["aromatic_delta"], "complete": False}
    assert not dominates(high, low_partial)

    # Even matching masks do not turn unknown values into zero evidence.
    same_mask_left = {**high, "missing_mask": ["aromatic_delta"],
                      "mean_absolute_delta": partial}
    same_mask_right = {**same_mask_left, "exposure_min": 0.4,
                       "mean_absolute_delta": {**partial, "charge_delta": 0.5}}
    assert dominates(same_mask_left, same_mask_right)


def test_epitope_signal_changes_do_not_change_mechanism_primary_or_foldx_region_order():
    candidates = [
        _candidate("a_best", [1, 2], exposure=0.9, magnitude=2.0, selection_rank=0),
        _candidate("a_next", [2, 3], exposure=0.7, magnitude=1.5, selection_rank=1),
        _candidate("b_best", [4, 5], exposure=0.9, magnitude=2.0, selection_rank=2),
        _candidate("b_next", [5, 6], exposure=0.7, magnitude=1.5, selection_rank=3),
    ]
    sectors = [V4Sector("VS001", [1, 2, 3], 1, 30, 30, 0),
               V4Sector("VS002", [4, 5, 6], 1, 30, 30, 0)]
    graph = {position: {neighbor: 6.0 for neighbor in (position - 1, position + 1)
                        if 1 <= neighbor <= 6} for position in range(1, 7)}
    areas = {position: 10.0 for position in range(1, 7)}

    def primary_ids(pool):
        rows = build_sector_recommendations(
            pool, sectors, areas, {"regional_recommendation_policy": POLICY}, graph=graph)
        return [row["candidate_id"] for row in rows]

    before_primary = primary_ids(candidates)
    before_foldx = _foldx_shortlist(candidates, 4, sectors=sectors,
                                    surface_areas=areas, graph=graph, regional_policy=True)
    before_foldx = [item.candidate_id for item in before_foldx]

    perturbed = deepcopy(candidates)
    for item in perturbed:
        item.epitope_score = 1.0 - item.epitope_score
        for option in item.mutations:
            option.epitope_supported = not option.epitope_supported
            option.epitope_prior = 1.0 - option.epitope_prior
    after_primary = primary_ids(perturbed)
    after_foldx = _foldx_shortlist(perturbed, 4, sectors=sectors,
                                   surface_areas=areas, graph=graph, regional_policy=True)
    after_foldx = [item.candidate_id for item in after_foldx]

    assert before_primary == after_primary
    assert before_primary == ["a_best", "b_best"]
    assert before_foldx == after_foldx
    assert before_foldx[:2] == ["a_best", "b_best"]


def test_regional_foldx_shortlist_is_balanced_and_respects_default_sixty_budget():
    sectors = []
    candidates = []
    graph = {}
    areas = {}
    for region_index in range(25):
        start = region_index * 4 + 1
        positions = [start, start + 1, start + 2]
        sector_id = f"VS{region_index + 1:03d}"
        sectors.append(V4Sector(sector_id, positions, 1, 30, 30, 0))
        for position in positions:
            graph[position] = {neighbor: 6.0 for neighbor in (position - 1, position + 1)
                               if neighbor in positions}
            areas[position] = 10.0
        for offset, pair in enumerate(((positions[0], positions[1]),
                                       (positions[1], positions[2]),
                                       (positions[0], positions[2]))):
            candidate = _candidate(f"{sector_id}_{offset}", pair,
                                   exposure=0.8 - 0.1 * offset,
                                   magnitude=1.0 + offset,
                                   selection_rank=offset)
            candidate.sector_id = sector_id
            candidates.append(candidate)

    selected = _foldx_shortlist(candidates, 60, sectors=sectors,
                                surface_areas=areas, graph=graph, regional_policy=True)
    assert len(selected) == 60
    counts = {sector.sector_id: 0 for sector in sectors}
    for candidate in selected:
        counts[candidate.sector_id] += 1
    assert min(counts.values()) >= 2
    assert max(counts.values()) <= 3


def test_normalized_worker_mutations_are_deduplicated_and_mapped_back_to_canonical_positions(tmp_path):
    class CaptureAdapter:
        name = "spurs_multi"
        timeout_seconds = 30

        def __init__(self):
            self.request = None

        def checked_result(self, request, result):
            self.request = request
            return result

    adapter = CaptureAdapter()
    result = PredictionResult(status="success", records=[])
    returned = _prediction(
        adapter, tmp_path / "normalized.pdb", "A", "DE", tmp_path,
        "experimental", {}, {"spurs_multi": result}, "spurs_multi",
        [["D1A", "E2A"], ["E2A", "D1A"], ["D1A", "E2A"]],
    )
    assert returned is result
    assert adapter.request.mutations == [["D1A", "E2A"]]

    candidate = _candidate("canonical", [101, 105])
    _annotate_stability(
        [candidate], [],
        PredictionResult(status="success", records=[{"mutations": ["D1A", "D2A"], "ddg": 0.2}]),
        PredictionResult(status="success", records=[{"mutations": ["D1A", "D2A"], "ddg": 0.3}]),
        {}, {}, position_map={1: 101, 2: 105},
    )
    assert candidate.stability.mutations == ("D101A", "D105A")
    assert candidate.stability.spurs_multi_ddg_kcal_mol == 0.2
    assert candidate.stability.foldx_ddg_kcal_mol == 0.3
