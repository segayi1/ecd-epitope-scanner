from ecd_epitope_scanner.scripts.freeze_cross_protein_candidates import (
    STRATA,
    _rows_by_author_id,
    select_stratified,
)


def _candidate(identifier, stratum, pair, band, mechanisms, area=100.0, risk=0.0):
    return {
        "candidate_id": identifier,
        "coupling_class": stratum,
        "site_pair": pair,
        "exposure_profile": band,
        "mechanisms": mechanisms,
        "proxy_support_area_a2": area,
        "positive_single_spurs_risk": risk,
    }


def test_candidate_strata_are_deterministic_and_do_not_repeat_site_pairs():
    pool = []
    for index, stratum in enumerate(STRATA):
        pool.extend(
            [
                _candidate(
                    f"{index}a",
                    stratum,
                    [index * 10 + 1, index * 10 + 2],
                    ["EXPOSED", "EXPOSED"],
                    ["volume"],
                    100.0,
                ),
                _candidate(
                    f"{index}a_alt",
                    stratum,
                    [index * 10 + 1, index * 10 + 2],
                    ["EXPOSED", "EXPOSED"],
                    ["charge"],
                    500.0,
                ),
                _candidate(
                    f"{index}b",
                    stratum,
                    [index * 10 + 3, index * 10 + 4],
                    ["BOUNDARY", "EXPOSED"],
                    ["charge"],
                    200.0,
                ),
            ]
        )
    first = select_stratified(pool)
    second = select_stratified(list(reversed(pool)))
    assert {
        key: [x["candidate_id"] for x in value["selected"]] for key, value in first.items()
    } == {key: [x["candidate_id"] for x in value["selected"]] for key, value in second.items()}
    for value in first.values():
        assert value["status"] == "FROZEN"
        assert len(value["selected"]) == 2
        assert len({tuple(x["site_pair"]) for x in value["selected"]}) == 2


def test_underfilled_stratum_is_empty_and_does_not_borrow_another_class():
    pool = [
        _candidate("one", STRATA[0], [1, 2], ["EXPOSED", "EXPOSED"], ["charge"]),
        _candidate("two", STRATA[1], [3, 4], ["EXPOSED", "EXPOSED"], ["volume"]),
        _candidate("three", STRATA[1], [5, 6], ["BOUNDARY", "EXPOSED"], ["aromatic"]),
    ]
    result = select_stratified(pool)
    assert result[STRATA[0]]["status"] == "EMPTY"
    assert result[STRATA[0]]["selected"] == []
    assert len(result[STRATA[1]]["selected"]) == 2
    assert result[STRATA[2]]["status"] == "EMPTY"


def test_unassigned_author_residue_ids_remain_unmapped():
    rows = [
        {"author_chain": "A", "author_residue": "12", "insertion_code": ""},
        {"author_chain": "A", "author_residue": "?", "insertion_code": "?"},
        {"author_chain": "B", "author_residue": "13", "insertion_code": ""},
    ]
    assert list(_rows_by_author_id(rows, "A")) == [(12, "")]
