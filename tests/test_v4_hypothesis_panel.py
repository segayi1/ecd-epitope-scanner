from types import SimpleNamespace

from ecd_epitope_scanner.joint_design_v4.candidate_search import generate_candidates
from ecd_epitope_scanner.joint_design_v4.hypotheses import (
    attach_hypothesis_responses,
    hypotheses_from_surface_v2,
)
from ecd_epitope_scanner.joint_design_v4.models import (
    ResidueKey,
    V4Candidate,
    V4MutationOption,
    V4Sector,
    V4SurfaceHypothesis,
)
from ecd_epitope_scanner.joint_design_v4.optimizer import select_panel


def _option(position):
    return V4MutationOption(
        residue=ResidueKey("0", "A", position, uniprot_pos=position, one_letter="D"),
        mutant_aa="A", notation=f"D{position}A", chemistry_class="charge_removal",
        spurs_ddg_kcal_mol=0.0, score=1.0, eligible=True,
    )


def test_surface_v2_footprints_remain_overlapping_and_response_crosses_sector():
    footprints = [
        SimpleNamespace(footprint_id="F1", residue_surface_area={1: 10, 2: 10},
                        region_id="R1", quality_flags=[]),
        SimpleNamespace(footprint_id="F2", residue_surface_area={2: 10, 3: 10},
                        region_id="R2", quality_flags=[]),
    ]
    hypotheses = hypotheses_from_surface_v2(SimpleNamespace(footprints=footprints))
    assert set(hypotheses[0].residue_weights) & set(hypotheses[1].residue_weights) == {2}
    candidate = V4Candidate("C1", "S1", [_option(2)], order=1)
    graph = {1: {2: 3.0}, 2: {1: 3.0, 3: 3.0}, 3: {2: 3.0}}
    attach_hypothesis_responses([candidate], hypotheses, graph)
    assert candidate.hypothesis_response["F1"] > 0
    assert candidate.hypothesis_response["F2"] > 0


def test_panel_prefers_probe_that_distinguishes_hypotheses():
    hypotheses = [V4SurfaceHypothesis("F1", {1: 1}), V4SurfaceHypothesis("F2", {2: 1})]
    broad = V4Candidate("A", "S1", [_option(1)], order=1,
                        hypothesis_response={"F1": 0.8, "F2": 0.8})
    specific = V4Candidate("B", "S2", [_option(2)], order=1,
                           hypothesis_response={"F1": 0.7, "F2": 0.0})
    panel, selected = select_panel([broad, specific], max_primary=1, hypotheses=hypotheses)
    assert [item.candidate_id for item in selected] == ["B"]
    assert panel.solver_status == "FOOTPRINT_IDENTIFIABILITY_GREEDY"
    assert panel.mean_pairwise_separation == 0.7


def test_candidate_search_keeps_each_mutation_order_when_beam_is_small():
    graph = {position: {} for position in range(1, 7)}
    for position in range(1, 6):
        graph[position][position + 1] = 2.0
        graph[position + 1][position] = 2.0
    sector = V4Sector("S1", list(graph), 1, 60.0, 60.0, 0.0)
    candidates = generate_candidates(
        [sector], [_option(position) for position in graph], graph,
        {"orders": [2, 3], "beam_width": 1, "shortlist_per_order": 1},
    )
    assert {item.order for item in candidates} == {2, 3}
