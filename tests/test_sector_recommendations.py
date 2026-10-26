import json
from pathlib import Path

import pytest

from ecd_epitope_scanner.joint_design_v4.models import (
    CombinationStabilityEvidence, ResidueKey, V4Candidate, V4MutationOption, V4Sector,
)
from ecd_epitope_scanner.joint_design_v4.sector_recommendations import build_sector_recommendations


def option(position, native="D", mutant="A"):
    return V4MutationOption(
        ResidueKey("0", "A", position, uniprot_pos=position, one_letter=native),
        mutant, f"{native}{position}{mutant}", "charge", 0.0, 1.0, True,
    )


def candidate(identifier, positions, *, passed=True, spurs=0.0, foldx=0.0, reasons=(), anchor="VS999"):
    return V4Candidate(
        identifier, anchor, [option(p) for p in positions], order=len(positions),
        stability=CombinationStabilityEvidence(tuple(), spurs, foldx, len(positions)),
        formal_eligible=passed, formal_exclusion_reasons=list(reasons),
        covered_positions_by_scenario={"geodesic_0": list(positions)},
    )


def sector(identifier, positions):
    return V4Sector(identifier, list(positions), 1, 100.0, 100.0, 0.0)


def test_actual_mutation_membership_and_coverage_precede_risk_and_anchor():
    low_risk = candidate("low", [1, 2], anchor="VS002")
    high_cover = candidate("high", [1, 3], anchor="VS002")
    high_cover.structural_risk = 100.0
    rows = build_sector_recommendations(
        [low_risk, high_cover], [sector("VS001", [1, 2, 3])], {1: 5, 2: 1, 3: 20}, {}, [high_cover],
    )
    assert rows[0]["candidate_id"] == "high"
    assert rows[0]["status"] == "PASSED"
    assert rows[0]["coverage_gain_area_a2"] == 25
    assert rows[0]["coverage_gain"] == pytest.approx(25 / 26)
    assert rows[0]["mutation_positions_in_sector"] == [1, 3]
    assert rows[0]["selected_in_panel"]


def test_nearest_threshold_uses_configured_threshold_and_eligibility_stays_false():
    far = candidate("far", [1, 2], passed=False, spurs=0.8, reasons=["COMBINATION_REJECT"])
    near = candidate("near", [1, 2], passed=False, spurs=0.65, reasons=["COMBINATION_REJECT"])
    rows = build_sector_recommendations([far, near], [sector("VS001", [1, 2])], {1: 10, 2: 10},
                                        {"spurs_multi": {"preferred_max_ddg": 0.6}})
    row = rows[0]
    assert row["candidate_id"] == "near"
    assert row["status"] == "NEAR_THRESHOLD"
    assert row["gate_distance"]["normalized_excess"] == pytest.approx(0.05 / 0.6)
    assert row["gate_details"][0]["observed"] == 0.65
    assert row["gate_details"][0]["threshold"] == 0.6
    assert not near.formal_eligible


def test_missing_foldx_is_one_gap_even_with_two_reasons_and_unknown_gain_stays_unknown():
    item = candidate("triple", [1, 2, 3], passed=False, foldx=None,
                     reasons=["SPURS_ONLY_REVIEW_REQUIRED", "HIGHER_ORDER_REQUIRES_FOLDX"])
    row = build_sector_recommendations([item], [sector("VS001", [1, 2, 3])], {1: 1, 2: 1, 3: 1}, {})[0]
    assert row["gate_distance"] == {"unmet_gate_count": 1, "normalized_excess": 0.0, "missing_evidence_count": 1}
    assert row["gate_details"][0]["observed"] is None
    assert set(row["failed_gates"]) == set(item.formal_exclusion_reasons)
    item.formal_exclusion_reasons.append("INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN")
    row = build_sector_recommendations([item], [sector("VS001", [1, 2, 3])], {1: 1, 2: 1, 3: 1}, {})[0]
    gain = next(gate for gate in row["gate_details"] if gate["gate"] == "higher_order_gain")
    assert gain["observed"] is None and gain["excess"] is None and gain["evidence_missing"]


def test_equal_failed_gate_count_prefers_measured_small_excess_to_more_missing_evidence():
    measured = candidate("measured", [1, 2], passed=False, spurs=0.51, foldx=None,
                         reasons=["COMBINATION_REJECT"])
    unknown = candidate("unknown", [1, 2, 3], passed=False, spurs=0.1, foldx=None,
                        reasons=["SPURS_ONLY_REVIEW_REQUIRED", "HIGHER_ORDER_REQUIRES_FOLDX",
                                 "INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN"])
    rows = build_sector_recommendations([unknown, measured], [sector("VS001", [1, 2, 3])],
                                        {1: 1, 2: 1, 3: 100}, {})
    row = rows[0]
    assert row["candidate_id"] == "measured"
    assert row["gate_distance"]["unmet_gate_count"] == 2
    assert row["gate_distance"]["missing_evidence_count"] == 1
    assert row["gate_distance"]["normalized_excess"] == pytest.approx(0.01 / 0.5)
    assert row["ranking_key"][:3] == [2, 1, row["gate_distance"]["normalized_excess"]]


@pytest.mark.parametrize("unsafe", ["native", "mutant", "hard", "disconnected", "order"])
def test_hard_protection_and_invalid_geometry_never_become_suggestions(unsafe):
    item = candidate("bad", [1, 2])
    if unsafe == "native":
        item.mutations[0] = option(1, native="C")
    elif unsafe == "mutant":
        item.mutations[0] = option(1, mutant="G")
    elif unsafe == "hard":
        item.mutations[0].hard_exclusions = ["DISULFIDE_PROTECTED"]
    elif unsafe == "disconnected":
        item.connected = False
    else:
        item.order = 4
    row = build_sector_recommendations([item], [sector("VS001", [1, 2])], {1: 1, 2: 1}, {})[0]
    assert row["status"] == "NO_COMBINATION"
    assert row["candidate_id"] is None and row["failed_gates"]


def test_cross_region_reuse_is_alternative_and_metrics_are_independent():
    shared = candidate("shared", [1, 2])
    review = candidate("review", [3, 4], passed=False, foldx=None, reasons=["SPURS_ONLY_REVIEW_REQUIRED"])
    formal = candidate("formal", [5, 6])
    formal.covered_positions_by_scenario["geodesic_0"] = [3, 4, 5, 6]
    rows = build_sector_recommendations([review, formal, shared], [
        sector("VS001", [1]), sector("VS002", [2]), sector("VS003", [3]), sector("VS004", [4]),
        sector("VS005", [5]), sector("VS006", [6]),
    ], dict.fromkeys(range(1, 7), 1), {})
    assert all(row["candidate_id"] is None for row in rows)
    assert [[item["candidate_id"] for item in row["alternatives"]] for row in rows] == [
        ["shared"], ["shared"], [], [], ["formal"], ["formal"]]
    assert rows[0]["alternatives"][0]["coverage_gain"] == 1
    assert rows[1]["alternatives"][0]["coverage_gain"] == 1
    gap = next(item for item in rows[2]["unavailable_principles"] if item["principle"] == "cross_region_coverage")
    assert gap["candidate_id"] == "review" and gap["gate_details"]
    assert gap["reason"] == "NO_QUALIFIED_PRINCIPLE_CANDIDATE"
    assert rows[4]["alternatives"][0]["coverage_gain_area_a2"] == 1


def test_permutations_do_not_change_assignment_or_ranking():
    candidates = [candidate("b", [1, 2]), candidate("a", [1, 2]), candidate("c", [3, 4])]
    sectors = [sector("VS002", [3, 4]), sector("VS001", [1, 2])]
    areas = dict.fromkeys(range(1, 5), 1)
    expected = build_sector_recommendations(candidates, sectors, areas, {})
    assert expected == build_sector_recommendations(list(reversed(candidates)), list(reversed(sectors)), areas, {})
    assert expected[0]["candidate_id"] == "a"


def test_latest_ctla4_replay_uses_real_membership_and_keeps_gate_rejections():
    run = Path(__file__).resolve().parents[1] / "workspace_data" / "5068b62450ec"
    if not run.exists():
        pytest.skip("Saved CTLA4 run is not available in this checkout")
    sectors = [V4Sector(**row) for row in json.loads((run / "probe_sectors_v4.json").read_text(encoding="utf-8"))]
    candidates = []
    for row in json.loads((run / "panel_candidates.json").read_text(encoding="utf-8")):
        row["mutations"] = [V4MutationOption(**{**option_row, "residue": ResidueKey(**option_row["residue"])})
                            for option_row in row["mutations"]]
        if row["stability"]:
            row["stability"] = CombinationStabilityEvidence(**row["stability"])
        candidates.append(V4Candidate(**row))
    graph = json.loads((run / "surface_graph_v41.json").read_text(encoding="utf-8"))
    areas = {node["uniprot_pos"]: node["surface_area_a2"] for node in graph["nodes"]}
    manifest = json.loads((run / "run_manifest.json").read_text(encoding="utf-8"))
    rows = build_sector_recommendations(candidates, sectors, areas, manifest["parameters"]["joint_v4_ml"])
    assert len(rows) == 9
    assert [row["status"] for row in rows] == ["PASSED", "NO_COMBINATION", "NO_COMBINATION",
        "PASSED", "PASSED", "PASSED", "PASSED", "NO_COMBINATION", "NEAR_THRESHOLD"]
    assert rows[5]["status"] == "PASSED"  # VS006's old display anchor falsely suggested no formal pool.
    by_id = {candidate.candidate_id: candidate for candidate in candidates}
    for row, display_sector in zip(rows, sectors):
        if row["candidate_id"] is not None:
            assert {option.residue.uniprot_pos for option in by_id[row["candidate_id"]].mutations} <= set(display_sector.residue_positions)
        for item in [row, *row["alternatives"]]:
            if item["candidate_id"] is not None:
                assert set(item["mutation_positions_in_sector"]) & set(display_sector.residue_positions)
                assert item["formal_eligible"] == by_id[item["candidate_id"]].formal_eligible
    assert rows == build_sector_recommendations(list(reversed(candidates)), list(reversed(sectors)), areas,
                                                manifest["parameters"]["joint_v4_ml"])
