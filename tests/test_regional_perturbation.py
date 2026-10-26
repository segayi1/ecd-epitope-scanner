from copy import deepcopy
from types import SimpleNamespace as NS

from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, V4Candidate, V4Sector
from ecd_epitope_scanner.joint_design_v4.regional_perturbation import (
    AXES, POLICY, candidate_impact, dominates, generate_regional_candidates,
    rank_impacts, annotate_local_interactions, audit_single_thresholds, RANKING_VERSION,
)
from ecd_epitope_scanner.joint_design_v4.evidence import SingleMutationEvidence
from ecd_epitope_scanner.joint_design_v4.coverage_profile import sequence_reasons


def _option(position, wt, mutant, exposure=0.5, descriptors=None, ddg=0.1):
    residue = ResidueKey("m", "A", position, uniprot_pos=position, one_letter=wt)
    return NS(residue=residue, mutant_aa=mutant, spurs_ddg_kcal_mol=ddg,
              chemical_descriptors=descriptors, impact_evidence={
                  "sidechain_exposure": {"value": exposure, "status": "COMPUTED"}})


def _candidate(identifier, options):
    return V4Candidate(candidate_id=identifier, sector_id="r", mutations=options,
                       order=len(options), structural_risk=0.0)


def test_impact_uses_six_axes_absolute_means_and_exposure_minimum():
    first = _option(1, "A", "D", exposure=.8)
    second = _option(2, "A", "V", exposure=.3)
    candidate = _candidate("pair", [first, second])
    impact = candidate_impact(candidate)
    assert impact["policy"] == POLICY
    assert impact["exposure_min"] == .3
    assert impact["exposure_mean"] == .55
    assert set(impact["axes"]) == set(AXES)
    assert impact["charge_delta"] == .5
    assert impact["complete"] is True


def test_unknown_histidine_is_masked_and_not_filled_with_zero():
    known = _candidate("known", [_option(1, "A", "D", exposure=.5)])
    unknown = _candidate("unknown", [_option(1, "A", "H", exposure=.5)])
    unknown.residue = None
    impact = candidate_impact(unknown)
    assert {"charge_delta", "donor_delta", "acceptor_delta"} <= set(impact["missing_axes"])
    assert impact["charge_delta"] is None
    assert impact["complete"] is False
    ranked = rank_impacts([unknown, known])
    assert ranked[0].candidate_id == "known"


def test_pareto_dominance_and_ordinal_utility_choose_high_axes():
    strong_axes = {axis: 1.0 for axis in AXES}
    weak_axes = {axis: 0.0 for axis in AXES}
    strong = _candidate("strong", [_option(1, "A", "D", exposure=.9, descriptors=strong_axes)])
    weak = _candidate("weak", [_option(1, "A", "S", exposure=.3, descriptors=weak_axes)])
    assert dominates(strong, weak)
    assert not dominates(weak, strong)
    ranked = rank_impacts([weak, strong])
    assert ranked[0].candidate_id == "strong"
    utilities = ranked[0].impact_evidence["ranking"]["mechanism_group_utilities"]
    assert all(value > 0.5 for value in utilities.values())


def test_reference_profile_is_frozen_when_candidates_are_filtered():
    candidates = []
    for index, exposure in enumerate((.2, .5, .9)):
        descriptors = {axis: float(index + 1) for axis in AXES}
        candidates.append(_candidate(str(index), [
            _option(index + 1, "A", "D", exposure=exposure, descriptors=descriptors)]))
    initially_ranked = rank_impacts(candidates)
    reference = deepcopy(initially_ranked[0].impact_evidence["regional_reference_profiles"])
    expected_order = [item.candidate_id for item in initially_ranked if item.candidate_id != "1"]
    retained = [deepcopy(item) for item in initially_ranked if item.candidate_id != "1"]
    reranked = rank_impacts(retained)
    assert [item.candidate_id for item in reranked] == expected_order
    for item in reranked:
        assert item.impact_evidence["regional_reference_profiles"] == reference
        assert item.impact_evidence["ranking"]["reference_size"] == 3


def test_rank_and_pareto_layers_are_stable_to_candidate_input_permutation():
    vector_high = {axis: 2.0 for axis in AXES}
    vector_low = {axis: 1.0 for axis in AXES}
    first = [_candidate("same_a", [_option(1, "A", "D", .8, vector_high)]),
             _candidate("same_b", [_option(2, "A", "D", .8, vector_high)]),
             _candidate("low", [_option(3, "A", "D", .2, vector_low)])]
    left = rank_impacts(deepcopy(first))
    right = rank_impacts(deepcopy(list(reversed(first))))
    assert [item.candidate_id for item in left] == [item.candidate_id for item in right]
    left_layers = {item.candidate_id: item.impact_evidence["ranking"]["pareto_layer"] for item in left}
    right_layers = {item.candidate_id: item.impact_evidence["ranking"]["pareto_layer"] for item in right}
    assert left_layers == right_layers


def _singles(sequence, positions, mutants="ACDEFGHIKLMNPQRSTVWY"):
    rows = []
    for position in positions:
        wt = sequence[position - 1]
        residue = ResidueKey("m", "A", position, uniprot_pos=position, one_letter=wt)
        for mutant in mutants:
            if mutant == wt:
                continue
            rows.append(SingleMutationEvidence(
                residue=residue, mutant_aa=mutant, spurs_ddg_kcal_mol=.2,
                chemistry_class="nominal", eligible=True,
                chemical_descriptors=None,
                impact_evidence={"policy": POLICY, "sidechain_exposure": {
                    "value": .4 + (position % 3) * .1, "status": "COMPUTED"}},
            ))
    return rows


def test_regional_generation_uses_all_sites_and_all_replacements_beyond_legacy_caps():
    sequence = "A" * 15
    positions = list(range(1, 15))
    singles = _singles(sequence, positions)
    graph = {position: {} for position in positions}
    for position in positions[:-1]:
        graph[position][position + 1] = 10.0
        graph[position + 1][position] = 10.0
    sector = V4Sector("r", positions, 1, 100.0, 100.0, 0.0)
    candidates, audit = generate_regional_candidates(
        singles, [sector], graph, {"maximum_geodesic_span_a": 24.0}, sequence,
        surface_areas={position: 10.0 for position in positions})
    region = audit["regions"][0]
    assert region["eligible_site_count"] == 14  # beyond the normal 12-site pool
    assert region["eligible_replacement_count"] > 3 * len(positions)
    assert region["raw_double_count"] > 3
    assert region["retained_double_count"] <= 20
    assert region["retained_triple_count"] <= 10
    assert {candidate.order for candidate in candidates} == {2, 3}
    assert all(candidate.impact_evidence["policy"] == POLICY for candidate in candidates)


def test_regional_cap_retains_site_thirteen_and_replacement_outside_legacy_top_three():
    sequence = "A" * 13
    positions = list(range(1, 14))
    rows = []
    for position in positions:
        wt = "A"
        residue = ResidueKey("m", "A", position, uniprot_pos=position, one_letter=wt)
        for index, mutant in enumerate(("D", "E", "F", "W", "Y")):
            magnitude = 10.0 if mutant == "W" else float(index + 1)
            rows.append(SingleMutationEvidence(
                residue=residue, mutant_aa=mutant,
                spurs_ddg_kcal_mol=(-.5 + index * .01) if mutant in {"D", "E", "F"} else .4 + index * .01,
                chemistry_class="nominal", eligible=True,
                chemical_descriptors={axis: magnitude for axis in AXES},
                impact_evidence={"policy": POLICY, "sidechain_exposure": {
                    "value": .99 if position in {12, 13} else .2, "status": "COMPUTED"}},
            ))
    graph = {position: {other: 1.0 for other in positions if other != position}
             for position in positions}
    sector = V4Sector("r", positions, 1, 100.0, 100.0, 0.0)
    candidates, audit = generate_regional_candidates(
        rows, [sector], graph, {"maximum_geodesic_span_a": 24.0}, sequence)
    doubles = [item for item in candidates if item.order == 2]
    triples = [item for item in candidates if item.order == 3]
    assert any(any(option.residue.uniprot_pos == 13 and option.mutant_aa == "W"
                   for option in item.mutations) for item in doubles)
    assert any(any(option.residue.uniprot_pos == 13 and option.mutant_aa == "W"
                   for option in item.mutations) for item in triples)
    legacy_three = {"D", "E", "F"}
    assert any(any(option.residue.uniprot_pos == 13 and option.mutant_aa not in legacy_three
                   for option in item.mutations) for item in doubles)
    assert audit["regions"][0]["eligible_site_count"] == 13
    assert audit["regions"][0]["double_lowest_risk_representative_id"] in {item.candidate_id for item in doubles}
    assert audit["regions"][0]["triple_lowest_risk_representative_id"] in {item.candidate_id for item in triples}


def test_region_generation_enforces_24a_geodesic_span_and_joint_sequon_gate():
    sequence = "AAAAAA"
    positions = [1, 2, 3, 4, 5, 6]
    graph = {position: {} for position in positions}
    for position in positions[:-1]:
        graph[position][position + 1] = 10.0
        graph[position + 1][position] = 10.0
    sector = V4Sector("r", positions, 1, 50.0, 50.0, 0.0)
    rows = _singles(sequence, positions, mutants="ANST")
    # Only the pair A1N + A3S jointly creates N-X-S. Neither single does.
    candidates, audit = generate_regional_candidates(
        rows, [sector], graph, {"maximum_geodesic_span_a": 24.0}, sequence)
    region = audit["regions"][0]
    assert region["double_span_rejected_count"] > 0
    assert region["double_sequence_rejected_count"] > 0
    assert all(not sequence_reasons(sequence, [
        (option.residue.uniprot_pos, option.residue.one_letter, option.mutant_aa)
        for option in candidate.mutations]) for candidate in candidates)
    # No retained pair/triple can span more than two 10 A graph edges.
    assert all(max((abs(a.residue.uniprot_pos - b.residue.uniprot_pos) * 10
                    for i, a in enumerate(candidate.mutations)
                    for b in candidate.mutations[i + 1:]), default=0) <= 20
               for candidate in candidates)


def test_unscored_preview_singles_do_not_create_mechanism_recommendations():
    rows = _singles("AAA", [1, 2, 3], mutants="DW")
    for row in rows:
        row.spurs_ddg_kcal_mol = None
    graph = {1: {2: 5}, 2: {1: 5, 3: 5}, 3: {2: 5}}
    sector = V4Sector("r", [1, 2, 3], 1, 30, 30, 0)
    candidates, audit = generate_regional_candidates(rows, [sector], graph, {}, "AAA")
    assert candidates == [] and audit["eligible_sites"] == 0
    assert audit["regions"][0]["double_site_set_deficit"] == 4


def test_charge_reversal_does_not_need_aromaticity_to_fill_mechanism_groups():
    reversal = _candidate("R", [_option(1, "D", "R", .8)])
    aromatic = _candidate("Y", [_option(1, "D", "Y", .8)])
    weak = _candidate("weak", [_option(1, "D", "E", .2)])
    ranked = rank_impacts([aromatic, weak, reversal])
    assert ranked[0].candidate_id == "R"
    assert reversal.impact_evidence["ranking_version"] == RANKING_VERSION
    groups = reversal.impact_evidence["ranking"]["mechanism_group_utilities"]
    assert groups == {"exposure": 1.0, "charge_polarity": 1.0, "size_aromatic": 1.0}
    assert reversal.impact_evidence["axes"]["aromatic_delta"] == 0
    assert not dominates(reversal, aromatic) and not dominates(aromatic, reversal)


def test_net_zero_and_nearby_opposite_charges_are_diagnostics_not_penalties():
    opposed = _candidate("opposed", [_option(1, "E", "K"), _option(2, "R", "D")])
    same_direction = _candidate("same", [_option(1, "E", "K"), _option(2, "D", "R")])
    before = [item.candidate_id for item in rank_impacts([opposed, same_direction])]
    annotate_local_interactions([opposed, same_direction], [NS(uniprot_pos=1, ca=(0, 0, 0)),
                                                          NS(uniprot_pos=2, ca=(8, 0, 0))])
    signed = candidate_impact(opposed)["signed_change_summary"]["charge_delta"]
    assert signed == {"status": "KNOWN", "net_change": 0.0, "absolute_sum": 4.0, "opposition_fraction": 1.0}
    assert candidate_impact(same_direction)["signed_change_summary"]["charge_delta"]["opposition_fraction"] == 0
    assert opposed.impact_evidence["local_interaction_diagnostics"]["pairs"][0]["nearby_opposed_charge"] is True
    assert opposed.impact_evidence["local_interaction_diagnostics"]["pairs"][0]["nearby_opposed_hydropathy"] is True
    assert same_direction.impact_evidence["local_interaction_diagnostics"]["pairs"][0]["nearby_opposed_charge"] is False
    assert [item.candidate_id for item in rank_impacts([opposed, same_direction])] == before
    assert candidate_impact(opposed)["charge_delta"] == 2.0


def test_unknown_signed_axes_and_missing_ca_are_not_zero_or_contacts():
    candidate = _candidate("unknown", [_option(1, "H", "K"), _option(2, "R", "D")])
    signed = candidate_impact(candidate)["signed_change_summary"]["charge_delta"]
    assert signed["status"] == "UNKNOWN" and signed["net_change"] is None
    assert signed["opposition_fraction"] is None and signed["absolute_sum"] is None
    annotate_local_interactions([candidate], [NS(uniprot_pos=1, ca=None), NS(uniprot_pos=2, ca=(8, 0, 0))])
    pair = candidate.impact_evidence["local_interaction_diagnostics"]["pairs"][0]
    assert pair["wt_ca_distance_a"] is None and pair["nearby_opposed_charge"] is None
    zero = _candidate("zero", [_option(1, "A", "V")])
    assert candidate_impact(zero)["signed_change_summary"]["charge_delta"]["status"] == "NO_CHANGE"
    assert candidate_impact(zero)["signed_change_summary"]["charge_delta"]["opposition_fraction"] is None


def test_high_exposure_review_band_is_audited_without_threshold_relaxation():
    singles = _singles("AA", [1, 2], mutants="DW")
    for row in singles:
        row.spurs_ddg_kcal_mol = .7
        row.eligible = False
        row.stability_band = "REVIEW"
        row.impact_evidence["sidechain_exposure"]["value"] = .9 if row.residue.uniprot_pos == 1 else None
    singles[-1].hard_exclusions = ["test-protection"]
    graph = {1: {2: 5}, 2: {1: 5}}
    sector = V4Sector("r", [1, 2], 1, 30, 30, 0)
    records = [NS(uniprot_pos=1, rsasa=.3), NS(uniprot_pos=2, rsasa=.8)]
    audit = audit_single_thresholds(singles, [sector], records, graph, {})["r"]
    assert audit["review_band_count"] == 4 and audit["otherwise_eligible_count"] == 3
    assert audit["high_rsasa_count"] == 1 and audit["high_sidechain_exposure_count"] == 2
    assert audit["sidechain_exposure_unknown_count"] == 1
    candidates, _ = generate_regional_candidates(singles, [sector], graph, {}, "AA")
    assert candidates == []
