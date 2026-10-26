"""Regression for missing-evidence strata and full sensitivity objectives."""
from copy import deepcopy
from types import SimpleNamespace as NS

from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, V4Candidate


from ecd_epitope_scanner.joint_design_v4 import regional_perturbation as regional


def _candidate(identifier, exposures, *, value=1.0, missing=()):
    mutations = []
    for position, exposure in enumerate(exposures, 1):
        descriptors = {axis: value for axis in regional.AXES}
        for axis in missing:
            descriptors[axis] = None
        residue = ResidueKey("m", "A", position, uniprot_pos=position, one_letter="A")
        option = NS(
            residue=residue,
            mutant_aa="D",
            spurs_ddg_kcal_mol=0.1,
            chemical_descriptors=descriptors,
            impact_evidence={"sidechain_exposure": {"value": exposure}},
        )
        mutations.append(option)
    return V4Candidate(identifier, "region", mutations, order=len(mutations))


def test_complete_and_missing_mask_strata_precede_incomparable_group_vectors():
    exposure_missing_best = _candidate("A_exposure_missing_best", [None], value=1.0)
    exposure_missing_low = _candidate("A_exposure_missing_low", [None], value=0.2)
    chemical_missing = ("charge_delta", "donor_delta", "acceptor_delta")
    three_axes_missing_best = _candidate("B_three_axes_missing_best", [0.9], value=1.0, missing=chemical_missing)
    three_axes_missing_low = _candidate("B_three_axes_missing_low", [0.1], value=0.2, missing=chemical_missing)

    ranked = regional.rank_impacts([
        three_axes_missing_best, exposure_missing_low,
        exposure_missing_best, three_axes_missing_low,
    ], assess_sensitivity=False)
    assert ranked[0].candidate_id == "A_exposure_missing_best"
    assert ranked[0].impact_evidence["regional_impact"]["missing_mask"] == ["sidechain_exposure"]
    assert ranked[2].candidate_id == "B_three_axes_missing_best"
    assert ranked[0].impact_evidence["sensitivity"]["status"] == "NOT_EVALUATED_GENERATION_STAGE"


def test_mean_exposure_sensitivity_recomputes_its_pareto_front():
    low_min_high_mean = _candidate("A_mean_wins", [0.2, 1.0])
    middle_min_mean = _candidate("B_min_wins", [0.3, 0.3])

    ranked = regional.rank_impacts([low_min_high_mean, middle_min_mean])
    sensitivity = ranked[0].impact_evidence["sensitivity"]
    assert ranked[0].candidate_id == "B_min_wins"
    assert sensitivity["selection_winners"]["ordinal"] == "B_min_wins"
    assert sensitivity["selection_winners"]["mean_exposure"] == "A_mean_wins"
    assert sensitivity["sensitive"] is True


def test_leave_one_exposure_out_recomputes_pareto_front():
    lower_exposure = _candidate("A_exposure_omitted_wins", [0.2])
    higher_exposure = _candidate("B_default_wins", [0.9])

    ranked = regional.rank_impacts([lower_exposure, higher_exposure])
    sensitivity = ranked[0].impact_evidence["sensitivity"]
    assert ranked[0].candidate_id == "B_default_wins"
    assert sensitivity["selection_winners"]["leave_one_out:exposure"] == "A_exposure_omitted_wins"
    assert sensitivity["sensitive"] is True


def test_raw_reference_profiles_survive_filtering_and_generation_skips_sensitivity():
    original = [
        _candidate("low", [0.1], value=0.1),
        _candidate("mid", [0.5], value=0.5),
        _candidate("high", [0.9], value=0.9),
    ]
    first = regional.rank_impacts(original, assess_sensitivity=False)
    reference_profile = deepcopy(first[0].impact_evidence["regional_reference_profiles"])
    retained = [deepcopy(item) for item in first if item.candidate_id != "mid"]
    reranked = regional.rank_impacts(retained, assess_sensitivity=True)
    assert [item.candidate_id for item in reranked] == ["high", "low"]
    for item in reranked:
        assert item.impact_evidence["regional_reference_profiles"] == reference_profile
        assert item.impact_evidence["ranking"]["reference_size"] == 3


def test_aromaticity_ablation_recomputes_front_without_imposing_aromatic_requirement():
    a = _candidate("A_no_aromatic_change", [0.7])
    b = _candidate("B_aromatic_change", [0.7])
    a.mutations[0].chemical_descriptors["aromatic_delta"] = 0
    b.mutations[0].chemical_descriptors["aromatic_delta"] = 1
    ranked = regional.rank_impacts([a, b])
    winners = ranked[0].impact_evidence["sensitivity"]["selection_winners"]
    assert winners["ordinal"] == "B_aromatic_change"
    assert winners["without_aromaticity"] == "A_no_aromatic_change"


def test_single_candidate_is_not_claimed_robust():
    candidate = _candidate("only", [0.7])
    regional.rank_impacts([candidate])
    sensitivity = candidate.impact_evidence["sensitivity"]
    assert sensitivity["status"] == "INSUFFICIENT_COMPARATORS"
    assert sensitivity["robust"] is None and sensitivity["sensitive"] is None
