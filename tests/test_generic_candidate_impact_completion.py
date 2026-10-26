import sys
from pathlib import Path
from types import SimpleNamespace as NS

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ecd_epitope_scanner.joint_design_v4.impact_completeness import AXES, _complete_display_ranges, card_evidence_completeness
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, V4MutationOption
from ecd_epitope_scanner.joint_design_v4.regional_perturbation import complete_candidate_impact


def test_generic_cross_region_impact_is_computed_without_mechanism_policy():
    options = [
        V4MutationOption(ResidueKey("1", "A", 1, uniprot_pos=1, one_letter="E"), "A", "E1A", "charge_change", .2, .5, True),
        V4MutationOption(ResidueKey("1", "A", 2, uniprot_pos=2, one_letter="K"), "L", "K2L", "charge_change", .3, .5, True),
    ]
    candidate = NS(candidate_id="cross", mutations=options, impact_evidence={})
    records = [
        NS(uniprot_pos=1, ca=(0.0, 0.0, 0.0), sidechain_exposure={"value": .7, "numerator_a2": 7.0,
           "denominator_a2": 10.0, "status": "KNOWN", "method": "fixture"}),
        NS(uniprot_pos=2, ca=(8.0, 0.0, 0.0), sidechain_exposure={"value": .6, "numerator_a2": 6.0,
           "denominator_a2": 10.0, "status": "KNOWN", "method": "fixture"}),
    ]

    complete_candidate_impact([candidate], records)

    evidence = candidate.impact_evidence
    assert "policy" not in evidence
    assert evidence["complete"] is True
    assert set(evidence["axes"]) == set(AXES)
    assert all(evidence["axes"][axis] is not None for axis in AXES)
    assert [row["notation"] for row in evidence["mutation_site_details"]] == ["E1A", "K2L"]
    pair = evidence["local_interaction_diagnostics"]["pairs"][0]
    assert pair["wt_ca_distance_a"] == 8.0
    assert pair["nearby_opposed_charge"] is True
    assert pair["nearby_opposed_hydropathy"] is False


def test_histidine_microstates_add_display_ranges_without_overwriting_unknown_ranking_axes():
    options = [
        V4MutationOption(ResidueKey("1", "A", 1, uniprot_pos=1, one_letter="E"), "H", "E1H", "charge_change", .2, .5, True),
        V4MutationOption(ResidueKey("1", "A", 2, uniprot_pos=2, one_letter="K"), "H", "K2H", "charge_change", .3, .5, True),
    ]
    candidate = NS(candidate_id="his-cross", mutations=options, impact_evidence={})
    records = [NS(uniprot_pos=1, ca=(0.0, 0.0, 0.0), sidechain_exposure={"value": .7, "numerator_a2": 7.0,
        "denominator_a2": 10.0, "status": "KNOWN", "method": "fixture"}),
        NS(uniprot_pos=2, ca=(8.0, 0.0, 0.0), sidechain_exposure={"value": .6, "numerator_a2": 6.0,
        "denominator_a2": 10.0, "status": "KNOWN", "method": "fixture"})]

    complete_candidate_impact([candidate], records)

    evidence = candidate.impact_evidence
    evidence["stability_repeat_audit"] = {"foldx_repeat_slots": [0, 1, 2, 3, 4],
        "foldx_repeat_values_kcal_mol": [1.0] * 5}
    assert evidence["axes"]["charge_delta"] is None
    assert evidence["axes"]["donor_delta"] is None
    assert evidence["axes"]["acceptor_delta"] is None
    ranges = evidence["display_chemistry_ranges"]
    assert _complete_display_ranges(ranges, [option.notation for option in options]), ranges
    assert ranges["scenario_count"] == 4
    assert len(ranges["microstates"]) == 2
    assert all(all(value is not None for value in scenario["mean_absolute_delta"].values())
               for scenario in ranges["scenarios"])
    assert all(row["nearby_opposed_charge"] is not None for scenario in ranges["scenarios"]
               for row in scenario["local_interaction_diagnostics"]["pairs"])

    card = {"sector_id": "S1", "candidate_id": candidate.candidate_id,
        "mutation_notations": [option.notation for option in options], "site_details": [
            {"notation": option.notation, "native_aa": option.residue.one_letter, "mutant_aa": option.mutant_aa,
             "chemical_descriptors": option.impact_evidence["chemical_descriptors"],
             "sidechain_exposure": option.impact_evidence["sidechain_exposure"]} for option in options],
        "impact_evidence": evidence, "regional_geometry_support": {"scenarios": {
            f"geodesic_{radius}": {"coverage_fraction": .5, "uncovered_positions": []} for radius in (9, 12, 15)}},
        "pair_distances": [{"positions": [1, 2], "geodesic_distance_a": 8.0}],
        "stability_evidence": {"spurs_multi_ddg_kcal_mol": .2, "foldx_ddg_kcal_mol": 1.0}}
    receipt = card_evidence_completeness([card], {"analysis.json": "a" * 64,
        "panel_candidates.json": "b" * 64, "sector_recommendations_v4.json": "c" * 64})
    assert receipt["status"] == "PASS", receipt["cards"][0]["missing_fields"]


def _complete_card():
    descriptors = {axis: 1.0 for axis in AXES}
    site = {"notation": "E1A", "native_aa": "E", "mutant_aa": "A",
            "chemical_descriptors": descriptors,
            "sidechain_exposure": {"value": .5, "numerator_a2": 5.0, "denominator_a2": 10.0,
                                   "status": "KNOWN", "method": "fixture"}}
    return {"sector_id": "S1", "candidate_id": "c1", "recommendation_role": "ALTERNATIVE",
            "mutation_notations": ["E1A"], "site_details": [site],
            "impact_evidence": {"axes": descriptors,
              "signed_change_summary": {axis: {"net_change": 1.0, "absolute_sum": 1.0,
                  "opposition_fraction": 0.0} for axis in ("charge_delta", "hydropathy_delta")},
              "local_interaction_diagnostics": {"pairs": []},
              "stability_repeat_audit": {"foldx_repeat_slots": [0, 1, 2, 3, 4],
                  "foldx_repeat_values_kcal_mol": [1.0] * 5}},
            "regional_geometry_support": {"scenarios": {f"geodesic_{r}":
              {"coverage_fraction": .5, "uncovered_positions": []} for r in (9, 12, 15)}},
            "pair_distances": [],
            "stability_evidence": {"spurs_multi_ddg_kcal_mol": .2, "foldx_ddg_kcal_mol": 1.0}}


def test_card_completeness_receipt_checks_all_axes_and_source_hashes():
    card = _complete_card()
    alternative = {**_complete_card(), "candidate_id": "c2", "recommendation_role": "ALTERNATIVE"}
    formal = {**_complete_card(), "candidate_id": "c3"}
    card["alternatives"] = [alternative]
    card["all_formal_local_candidates"] = [formal]
    receipt = card_evidence_completeness([card], {"analysis.json": "a" * 64,
        "panel_candidates.json": "b" * 64, "sector_recommendations_v4.json": "c" * 64})
    assert receipt["status"] == "PASS"
    assert receipt["card_count"] == 3
    assert {row["role"] for row in receipt["cards"]} == {"PRIMARY", "ALTERNATIVE", "ALL_FORMAL_LOCAL"}
    assert {row["parent_sector_id"] for row in receipt["cards"]} == {"S1"}
    assert set(receipt["inputs_sha256"]) == {"analysis.json", "panel_candidates.json", "sector_recommendations_v4.json"}

    del card["impact_evidence"]["axes"]["aromatic_delta"]
    incomplete = card_evidence_completeness([card], receipt["inputs_sha256"])
    assert incomplete["status"] == "INCOMPLETE"
    assert "impact_evidence.axes.aromatic_delta" in incomplete["cards"][0]["missing_fields"]


def test_card_publication_receipt_hashes_final_analysis_and_preserves_failure(tmp_path):
    import json
    from ecd_epitope_scanner.joint_design_v4.impact_completeness import write_card_acceptance
    from ecd_epitope_scanner.util import sha256_file
    for name in ('panel_candidates.json', 'sector_recommendations_v4.json'):
        (tmp_path / name).write_text('[]', encoding='utf-8')
    analysis = {'design_summary': {'sector_recommendations': [_complete_card()]}}
    receipt = write_card_acceptance(tmp_path, analysis)
    assert receipt['status'] == 'PASS'
    assert receipt['inputs_sha256']['analysis.json'] == sha256_file(tmp_path / 'analysis.json')
    assert json.loads((tmp_path / 'analysis.json').read_text()) == analysis
    del analysis['design_summary']['sector_recommendations'][0]['stability_evidence']['foldx_ddg_kcal_mol']
    receipt = write_card_acceptance(tmp_path, analysis)
    assert receipt['status'] == 'INCOMPLETE'
    assert analysis['design_summary']['display_evidence_acceptance']['incomplete_card_count'] == 1
    assert (tmp_path / 'card_evidence_completeness.json').is_file()
