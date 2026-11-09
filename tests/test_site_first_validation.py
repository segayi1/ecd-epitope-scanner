from types import SimpleNamespace as NS
from copy import deepcopy
import pytest

from ecd_epitope_scanner.joint_design_v4.models import PredictionResult
from ecd_epitope_scanner.joint_design_v4.models import (
    CombinationStabilityEvidence, ResidueKey, V4Candidate, V4MutationOption, V4Sector,
)
from ecd_epitope_scanner.joint_design_v4.site_first_display import (
    apply_display_qualification, apply_site_first_repeat_review, confirmed_repeat_pass,
)
from ecd_epitope_scanner.joint_design_v4.site_first_validation import (
    PROTOCOL, SiteFirstValidation,
)


def candidate(identifier, mutation):
    return NS(candidate_id=identifier, mutations=[NS(notation=mutation)])


@pytest.mark.parametrize("case", ["pass", "failed_median", "missing_slot", "duplicate_slot",
                                  "invalid_value", "cross_execution", "known_review", "legacy"])
def test_ten_slot_confirmation_requires_complete_authorized_evidence(case):
    item = NS(qualification="ELIGIBLE", mutations=[], impact_evidence={})
    values = [2.5, 1., 1.] + [1.] * 7
    if case == "failed_median":
        values[3:] = [3.] * 7
    if case == "missing_slot":
        values.pop()
    if case == "invalid_value":
        values[-1] = float("nan")
    result = {"status": "PASS", "scores": values, "slots": list(range(len(values))),
              "repeat_audit": {"foldx_repeat_slot_values_kcal_mol": values,
                               "foldx_repeat_slots": list(range(len(values)))}}
    if case == "duplicate_slot":
        result["slots"][-1] = 8
        result["repeat_audit"]["foldx_repeat_slots"][-1] = 8
    if case == "cross_execution":
        result["repeat_audit"]["cross_execution_conflict"] = True
    if case == "known_review":
        item.impact_evidence.update(policy="regional_mechanism_pareto_v1", missing_axes=["charge_delta"])
    verdict = apply_site_first_repeat_review(item, result,
        confirmation_policy=None if case == "legacy" else "ten_slot_median_v1")
    assert (verdict is None) == (case == "pass")
    if case == "pass":
        audit = result["repeat_audit"]
        assert confirmed_repeat_pass(audit)
        assert audit["recommendation_risk_label"] == "TEN_SLOT_MEDIAN_CONFIRMED"
        assert audit["sampling_decision_sensitive"] is True
        assert audit["foldx_repeat_above_gate_count"] == 1
        for key, value in [("sampling_decision_sensitive", False),
                           ("recommendation_risk_label", "SAMPLING_SENSITIVE_REVIEW"),
                           ("foldx_repeat_slots", [False, 1, 2, 3, 4, 5, 6, 7, 8, 9])]:
            assert not confirmed_repeat_pass({**audit, key: value})


def test_ten_slot_confirmation_keeps_shared_priority_heads_without_extra_physical_work():
    heads = [candidate(f"{sector}|priority", "A1D") for sector in ("north", "south")]
    for item in heads:
        item.qualification, item.impact_evidence = "ELIGIBLE", {}
        item.mutations[0].review_required = False
    store = Store({(sector, "coverage"): [item, candidate(f"{sector}|backup", "A2D")]
                   for sector, item in zip(("north", "south"), heads)})
    calls = []

    def foldx(batch, *, slot_start, runs):
        calls.append((tuple(c.candidate_id for c in batch), slot_start, runs))
        return PredictionResult("success", records=[pair_row(c,
            [2.5, 1., 1.] if slot_start == 0 else [1.] * runs, slot_start) for c in batch])

    saved = SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx,
        foldx_recheck=lambda c, r: apply_site_first_repeat_review(
            c, r, confirmation_policy="ten_slot_median_v1"), foldx_slot_limit=10)
    assert [(start, runs) for _, start, runs in calls] == [(0, 3), (3, 7)]
    assert store.outcomes == []
    assert store.completed == [((sector, "coverage"), f"{sector}|priority") for sector in ("north", "south")]
    assert saved["foldx"]["ledger"]["attempted_slots"] == 10
    assert saved["foldx"]["ledger"]["reserved_slots"] == 0
    gates = saved["foldx"]["results"][("A1D",)]["candidate_rechecks"]
    assert all(confirmed_repeat_pass(g["repeat_audit"]) for g in gates.values())
    assert all(g["display_status"] == "DISPLAY_ELIGIBLE" for g in gates.values())


class Store:
    def __init__(self, queues):
        self.queues = {lane: list(items) for lane, items in queues.items()}
        self.states = {lane: "ready" for lane in queues}
        self.outcomes, self.completed = [], []

    def queue_states(self):
        return list(self.queues)

    def status(self, *lane):
        return {"status": self.states[tuple(lane)]}

    def head(self, *lane, limit=4):
        return self.queues[tuple(lane)][:limit]

    def record_outcome(self, *args):
        lane, candidate_id, outcome = args[:2], args[2], args[3]
        assert self.queues[tuple(lane)][0].candidate_id == candidate_id
        self.outcomes.append((tuple(lane), self.queues[tuple(lane)].pop(0).candidate_id, outcome))

    def complete(self, *args):
        lane, candidate_id = args[:2], args[2]
        assert self.queues[tuple(lane)][0].candidate_id == candidate_id
        self.completed.append((tuple(lane), self.queues[tuple(lane)].pop(0).candidate_id))
        self.states[tuple(lane)] = "complete"

    def mark_queue_state(self, *args):
        self.states[tuple(args[:2])] = args[2]

    def mark_all_technical_failure(self):
        for lane, state in self.states.items():
            if state == "ready":
                self.states[lane] = "technical_failure"

    def mark_all_budget_limited(self, limit):
        for lane, state in self.states.items():
            if state == "ready":
                self.states[lane] = "budget_limited"

    def ledger(self, _validation):
        return {"outcomes": self.outcomes, "completed": self.completed}


class CopyingStore(Store):
    def head(self, *lane, limit=4):
        return deepcopy(self.queues[tuple(lane)][:limit])


def test_ranked_prefetch_shares_signatures_and_consumes_only_formal_heads():
    class PrefetchStore(Store):
        def prefetch_groups(self, *lane):
            for item in self.queues[tuple(lane)]:
                yield [item]

    store = PrefetchStore({("s1", "coverage"): [candidate("a", "A1D"), candidate("later", "A2D")],
                           ("s2", "coverage"): [candidate("shared", "A1D"), candidate("other", "A3D")]})
    requests = []
    def predict(batch):
        requests.append([c.mutations[0].notation for c in batch])
        return spurs(batch)
    def foldx(batch, *, slot_start, runs):
        return PredictionResult("success", records=[pair_row(c, [1.] * runs, slot_start) for c in batch])
    result = SiteFirstValidation().run_indexed(store, predict, foldx, spurs_request_limit=3,
        spurs_prefetch_policy="ranked_window_v1", spurs_prefetch_window=3)
    assert requests == [["A1D", "A2D", "A3D"]]
    assert len(store.completed) == 2 and not store.outcomes
    assert result["spurs"]["ledger"]["prefetched_unused"] == 2
    assert result["foldx"]["ledger"]["attempted_slots"] == 3


@pytest.mark.parametrize("limit,expected", [(10000, 1024), (70, 70)])
def test_ranked_window_is_bounded_and_split_into_64(limit, expected):
    class PrefetchStore(Store):
        def prefetch_groups(self, *lane):
            for item in self.queues[tuple(lane)]:
                yield [item]
    items = [candidate(str(i), f"A{i}D") for i in range(1, 1101)]
    store = PrefetchStore({("s", "coverage"): items})
    calls, progress, events = [], [], []
    def predict(batch):
        calls.append(len(batch))
        return spurs(batch)
    def foldx(batch, *, slot_start, runs):
        return PredictionResult("success", records=[pair_row(c, [1.] * runs, slot_start) for c in batch])
    result = SiteFirstValidation().run_indexed(store, predict, foldx, spurs_request_limit=limit,
        spurs_prefetch_policy="ranked_window_v1", spurs_progress=lambda done, total: progress.append((done, total)),
        spurs_audit=lambda event, details: events.append((event, details.copy())))
    assert sum(calls) == expected and max(calls) == 64
    assert progress[-1] == (expected, expected)
    assert events[0][0] == "window_started" and events[-1][0] == "window_finished"
    assert sum(event == "chunk" for event, _ in events) == len(calls)
    assert result["spurs"]["ledger"]["prefetched_unused"] == expected - 1


def spurs(batch, scores=None):
    scores = scores or {}
    return PredictionResult("success", records=[
        {"mutations": [c.mutations[0].notation], "ddg": scores.get(c.candidate_id, .2)}
        for c in batch
    ], provenance={"cache_hit": True})


def pair_row(c, values, start=0):
    return {"mutations": [m.notation for m in c.mutations], "paired_models": [
        {"sampling_slot": start + i, "ddg_kcal_mol": value, "source_path": f"s{start+i}"}
        for i, value in enumerate(values)
    ]}


def test_foldx_batches_multiple_heads_and_matches_unordered_rows():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("s1", "coverage"): [a], ("s2", "coverage"): [b]})
    calls = []

    def foldx(batch, *, slot_start, runs):
        calls.append((list(c.candidate_id for c in batch), slot_start, runs))
        return PredictionResult("success", records=[pair_row(c, [1.] * runs, slot_start)
                                                     for c in reversed(batch)],
                                provenance={"cache_hit": True})

    result = SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx)
    assert calls == [(["a", "b"], 0, 3)]
    assert len(store.completed) == 2
    assert result["foldx"]["results"][("A1D",)]["scores"] == [1., 1., 1.]
    assert result["foldx"]["results"][("A1D",)]["repeat_audit"]["execution_protocol"] == PROTOCOL
    assert result["spurs"]["ledger"]["batches"][0]["provenance"]["cache_hit"] is True
    assert result["foldx"]["results"][("A1D",)]["provenance"][0]["cache_hit"] is True


def test_foldx_prioritizes_an_unattempted_site_when_budget_allows():
    shared1, shared2, unseen = candidate("shared1", "A1D"), candidate("shared2", "A2D"), candidate("new", "C3D")
    shared1.impact_evidence = shared2.impact_evidence = {"site_set_id": "same"}
    unseen.impact_evidence = {"site_set_id": "different"}
    store = Store({("s1", "coverage"): [shared1],
                   ("s2", "coverage"): [shared2],
                   ("s3", "coverage"): [unseen]})
    calls = []

    def foldx(batch, *, slot_start, runs):
        calls.append([item.candidate_id for item in batch])
        return PredictionResult("success", records=[
            pair_row(item, ([1.] * 3 if slot_start == 0 else [1.] * runs), slot_start)
            for item in batch
        ])

    result = SiteFirstValidation().run_indexed(
        store, lambda batch: spurs(batch), foldx, foldx_slot_limit=30, max_workers=2
    )
    assert calls == [["shared1", "new"], ["shared2"]]
    assert result["foldx"]["ledger"]["unseen_site_candidates"] == 2
    assert result["foldx"]["ledger"]["attempted_site_keys"] == ["different", "same"]
    assert result["foldx"]["ledger"]["slot_limit"] == 30


def test_mixed_three_of_ten_rules_and_releases_reserved_budget():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("s1", "coverage"): [a], ("s2", "coverage"): [b]})
    calls = []

    def foldx(batch, *, slot_start, runs):
        calls.append((slot_start, runs, [c.candidate_id for c in batch]))
        rows = []
        for c in batch:
            vals = ([1., 1.2, 2.1] if c is a else [2.1, 2.2, .5]) if slot_start == 0 else [1.] * 7
            rows.append(pair_row(c, vals, slot_start))
        return PredictionResult("success", records=rows)

    result = SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx)
    assert calls == [(0, 3, ["a", "b"]), (3, 7, ["a"])]
    assert store.completed == [(('s1', 'coverage'), 'a')]
    assert store.outcomes == [(('s2', 'coverage'), 'b', 'FOLDX_SCIENTIFIC_FAILURE')]
    assert result["foldx"]["results"][("A1D",)]["status"] == "PASS"
    assert result["foldx"]["results"][("A2D",)]["status"] == "SCIENTIFIC_FAILURE"


def test_partial_foldx_keeps_good_row_and_stops_failed_lane():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("a", "coverage"): [a], ("b", "coverage"): [b]})

    def foldx(batch, *, slot_start, runs):
        good = pair_row(a, [1.] * runs, slot_start)
        bad = {"mutations": ["A2D"], "status": "failed", "errors": [{"code": "timeout"}]}
        return PredictionResult("partial", records=[bad, good], errors=[{"code": "one_failed"}])

    result = SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx)
    assert result["foldx"]["results"][("A1D",)]["status"] == "PASS"
    assert result["foldx"]["results"][("A2D",)]["status"] == "TECHNICAL_FAILURE"
    assert store.completed == [(('a', 'coverage'), 'a')]
    assert result["foldx"]["ledger"]["technical_failure"] is True


def test_unknown_foldx_signature_stops_dispatch_without_invalidating_good_rows():
    a, b, later = candidate("a", "A1D"), candidate("b", "A2D"), candidate("later", "A3D")
    store = Store({("a", "coverage"): [a, later], ("b", "coverage"): [b]})
    calls = []

    def foldx(batch, *, slot_start, runs):
        calls.append((slot_start, runs))
        rows = [pair_row(item, [1.] * runs, slot_start) for item in batch]
        rows.append({"mutations": ["A99D"], "paired_models": []})
        return PredictionResult("partial", records=rows)

    result = SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx)
    assert calls == [(0, 3)]
    assert all(result["foldx"]["results"][signature]["status"] == "PASS"
               for signature in (("A1D",), ("A2D",)))
    assert store.completed == [(('a', 'coverage'), 'a'), (('b', 'coverage'), 'b')]
    assert ("A3D",) not in result["foldx"]["results"]
    assert result["foldx"]["ledger"]["attempted_slots"] == 6
    assert result["foldx"]["ledger"]["technical_failure"] is True
    assert any(error["code"] == "FOLDX_UNKNOWN_SIGNATURE"
               for batch in result["foldx"]["ledger"]["batches"]
               for stage in batch["stages"] for error in stage["errors"])


def test_duplicate_foldx_signature_fails_only_that_candidate_and_keeps_both_rows():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("a", "coverage"): [a], ("b", "coverage"): [b]})

    def foldx(batch, *, slot_start, runs):
        rows = [pair_row(item, [1.] * runs, slot_start) for item in batch]
        duplicated = pair_row(b, [1.] * runs, slot_start)
        duplicated["paired_models"][0]["ddg_kcal_mol"] = 9.0
        rows.append(duplicated)
        return PredictionResult("partial", records=rows)

    result = SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx)
    assert result["foldx"]["results"][("A1D",)]["status"] == "PASS"
    bad = result["foldx"]["results"][("A2D",)]
    assert bad["status"] == "TECHNICAL_FAILURE"
    assert len(bad["paired_models"]) == 6
    assert len(bad["audit"]["stages"][0]["records"]) == 2
    assert result["foldx"]["ledger"]["attempted_slots"] == 6
    assert any(error["code"] == "FOLDX_DUPLICATE_SIGNATURE"
               for batch in result["foldx"]["ledger"]["batches"]
               for stage in batch["stages"] for error in stage["errors"])


def test_foldx_row_without_identity_stops_dispatch_and_keeps_good_rows():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("a", "coverage"): [a], ("b", "coverage"): [b]})

    def foldx(batch, *, slot_start, runs):
        rows = [pair_row(item, [1.] * runs, slot_start) for item in batch]
        rows.append({"paired_models": [{"sampling_slot": 0, "ddg_kcal_mol": 100.}]})
        return PredictionResult("partial", records=rows)

    result = SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx)
    assert all(result["foldx"]["results"][signature]["status"] == "PASS"
               for signature in (("A1D",), ("A2D",)))
    assert store.completed == [(('a', 'coverage'), 'a'), (('b', 'coverage'), 'b')]
    assert result["foldx"]["ledger"]["attempted_slots"] == 6
    errors = [error for batch in result["foldx"]["ledger"]["batches"]
              for stage in batch["stages"] for error in stage["errors"]]
    assert any(error["code"] == "FOLDX_MISSING_SIGNATURE" and "row" in error for error in errors)


def test_spurs_prefetch_is_limited_to_four_and_duplicate_heads_score_once():
    items = [candidate(f"c{i}", f"A{i+1}D") for i in range(6)]
    store = Store({("a", "coverage"): items, ("b", "coverage"): [items[0]]})
    seen = []

    def predict(batch):
        seen.extend(c.candidate_id for c in batch)
        return spurs(batch, {c.candidate_id: 2.5 for c in batch})

    SiteFirstValidation().run_indexed(store, predict, lambda *a, **k: None)
    assert set(seen) == {"c0", "c1", "c2", "c3", "c4", "c5"}
    assert len(seen) == len(set(seen)) == 6
    assert store.outcomes[0][1] == "c0"  # c4/c5 prefetch never displaced queue order.


def test_spurs_recheck_failure_advances_to_next_candidate():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("s", "coverage"): [a, b]})
    seen = []

    def recheck(item, result):
        if item is a:
            return False
        return True

    def foldx(batch, *, slot_start, runs):
        seen.extend(c.candidate_id for c in batch)
        return PredictionResult("success", records=[pair_row(c, [1.] * runs, slot_start) for c in batch])

    SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx, spurs_recheck=recheck)
    assert seen == ["b"]
    assert store.outcomes[0][2] == "SPURS_SCIENTIFIC_FAILURE"
    assert store.completed[-1][1] == "b"


def test_foldx_science_gate_failure_advances_lane_to_next_candidate():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("s", "coverage"): [a, b]})
    seen = []

    def foldx(batch, *, slot_start, runs):
        seen.extend(c.candidate_id for c in batch)
        return PredictionResult("success", records=[pair_row(c, [1.] * runs, slot_start) for c in batch])

    def recheck(item, _result):
        return item is not a

    SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx, foldx_recheck=recheck)
    assert seen == ["a", "b"]
    assert store.outcomes[0][2] == "FOLDX_SCIENTIFIC_FAILURE"
    assert store.completed[-1][1] == "b"


def test_same_sector_candidate_id_uses_canonical_object_across_lane_copies():
    a, b = candidate("north|same", "A1D"), candidate("north|same", "A1D")
    a.region = b.region = "north"
    a.impact_evidence = b.impact_evidence = {}
    store = CopyingStore({("north", "epitope_support"): [a],
                          ("north", "sector_coverage"): [b]})
    rechecked = []

    def foldx(batch, *, slot_start, runs):
        assert len(batch) == 1
        return PredictionResult("success", records=[pair_row(batch[0], [1.] * runs, slot_start)])

    def recheck(item, result):
        rechecked.append(item)
        item.stability = "checked"
        item.impact_evidence["stability_repeat_audit"] = result["repeat_audit"]

    validation = SiteFirstValidation()
    result = validation.run_indexed(store, lambda batch: spurs(batch), foldx,
                                    foldx_recheck=recheck)
    assert rechecked == [validation.materialized_candidates["north|same"]]
    assert validation.materialized_candidates["north|same"].stability == "checked"
    assert store.completed == [(('north', 'epitope_support'), 'north|same'),
                               (('north', 'sector_coverage'), 'north|same')]
    assert result["foldx"]["results"][("A1D",)]["candidate_rechecks"]["north|same"]["status"] == "PASS"


def test_same_signature_different_regions_recheck_formal_gate_per_candidate():
    north = candidate("north|same", "A1D")
    south = candidate("south|same", "A1D")
    north.region, south.region = "north", "south"
    north.impact_evidence = south.impact_evidence = {}
    store = CopyingStore({("north", "epitope_support"): [north],
                          ("south", "epitope_support"): [south]})
    foldx_calls, rechecked = [], []

    def foldx(batch, *, slot_start, runs):
        foldx_calls.append([item.candidate_id for item in batch])
        return PredictionResult("success", records=[pair_row(batch[0], [1.] * runs, slot_start)])

    def recheck(item, result):
        rechecked.append(item.candidate_id)
        item.stability = f"checked:{item.region}"
        item.impact_evidence["stability_repeat_audit"] = result["repeat_audit"]
        return item.region != "north"

    validation = SiteFirstValidation()
    result = validation.run_indexed(store, lambda batch: spurs(batch), foldx,
                                    foldx_recheck=recheck)
    assert foldx_calls == [["north|same"]]  # Physical result is reused by mutation signature.
    assert rechecked == ["north|same", "south|same"]
    assert store.outcomes == [(('north', 'epitope_support'), 'north|same', 'FOLDX_SCIENTIFIC_FAILURE')]
    assert store.completed == [(('south', 'epitope_support'), 'south|same')]
    assert validation.materialized_candidates["north|same"].stability == "checked:north"
    assert validation.materialized_candidates["south|same"].stability == "checked:south"
    gates = result["foldx"]["results"][("A1D",)]["candidate_rechecks"]
    assert gates["north|same"]["status"] == "SCIENTIFIC_FAILURE"
    assert gates["south|same"]["status"] == "PASS"


def test_display_review_advances_lane_and_reuses_same_ten_slot_foldx_result():
    review = candidate("sensitive", "A1D")
    replacement = candidate("replacement", "A1D")
    store = Store({("s", "mechanism_pareto"): [review, replacement]})
    foldx_calls = []

    def foldx(batch, *, slot_start, runs):
        foldx_calls.append((slot_start, runs))
        # Median passes the unchanged 2 kcal/mol gate, while repeats span it.
        values = [2.4, 1.0, 1.0] if slot_start == 0 else [1.0] * runs
        return PredictionResult("success", records=[pair_row(batch[0], values, slot_start)])

    def recheck(item, result):
        if item.candidate_id == "sensitive":
            result["repeat_audit"].update(
                recommendation_risk_label="SAMPLING_SENSITIVE_REVIEW",
                sampling_decision_sensitive=True,
            )
            result["display_review_reasons"] = ["SAMPLING_SENSITIVE_REVIEW"]
            return "DISPLAY_REVIEW"

    validation = SiteFirstValidation()
    result = validation.run_indexed(store, lambda batch: spurs(batch), foldx,
                                    foldx_recheck=recheck)
    physical = result["foldx"]["results"][("A1D",)]
    gates = physical["candidate_rechecks"]
    assert physical["status"] == "PASS" and physical["median"] == 1.0
    assert physical["slots"] == list(range(10))
    assert foldx_calls == [(0, 3), (3, 7)]
    assert gates["sensitive"]["status"] == "PASS"
    assert gates["sensitive"]["display_status"] == "DISPLAY_REVIEW"
    assert gates["replacement"]["status"] == "PASS"
    assert gates["replacement"]["display_status"] == "DISPLAY_ELIGIBLE"
    assert store.outcomes == [(("s", "mechanism_pareto"), "sensitive", "FOLDX_DISPLAY_REVIEW")]
    assert store.completed == [(("s", "mechanism_pareto"), "replacement")]


def test_real_site_first_repeat_review_rejects_shared_cross_region_signature_then_advances():
    from ecd_epitope_scanner.joint_design_v4.optimizer import select_panel
    from ecd_epitope_scanner.joint_design_v4.sector_recommendations import build_sector_recommendations

    def scored(identifier, sector_id, notation):
        native, position, mutant = notation[0], int(notation[1:-1]), notation[-1]
        return V4Candidate(
            identifier, sector_id,
            [V4MutationOption(ResidueKey("0", "A", position, uniprot_pos=position,
                                         one_letter=native), mutant, notation, "charge",
                              0.1, 1.0, True)],
            order=1, formal_eligible=True,
            stability=CombinationStabilityEvidence((notation,), 0.1, 1.0, 1),
            impact_evidence={},
        )

    north_shared = scored("north|shared", "north", "A1D")
    south_shared = scored("south|shared", "south", "A1D")
    north_backup = scored("north|backup", "north", "A2D")
    south_backup = scored("south|backup", "south", "A2D")
    store = Store({("north", "mechanism_pareto"): [north_shared, north_backup],
                   ("south", "mechanism_pareto"): [south_shared, south_backup]})
    physical_signatures = []

    def foldx(batch, *, slot_start, runs):
        signature = tuple(sorted(m.notation for m in batch[0].mutations))
        physical_signatures.append((signature, slot_start, runs))
        values = ([2.5, 1.0, 1.0] if signature == ("A1D",) and slot_start == 0
                  else [1.0] * runs)
        return PredictionResult("success", records=[pair_row(batch[0], values, slot_start)])

    def recheck(item, result):
        verdict = apply_site_first_repeat_review(item, result)
        item.stability = CombinationStabilityEvidence(
            tuple(sorted(m.notation for m in item.mutations)), 0.1, result["median"],
            item.order, consensus_status="CONSENSUS_PASS",
            provenance={"repeat_audit": result["repeat_audit"]})
        return verdict

    validation = SiteFirstValidation()
    saved = validation.run_indexed(
        store, lambda batch: spurs(batch), foldx, foldx_recheck=recheck,
        foldx_slot_limit=40)
    physical = saved["foldx"]["results"][("A1D",)]
    gates = physical["candidate_rechecks"]
    assert physical["status"] == "PASS" and physical["slots"] == list(range(10))
    assert gates["north|shared"]["display_status"] == "DISPLAY_REVIEW"
    assert gates["south|shared"]["display_status"] == "DISPLAY_REVIEW"
    assert gates["north|shared"]["repeat_audit"]["sampling_decision_sensitive"] is True
    assert gates["south|shared"]["repeat_audit"]["sampling_decision_sensitive"] is True
    assert store.outcomes == [
        (("north", "mechanism_pareto"), "north|shared", "FOLDX_DISPLAY_REVIEW"),
        (("south", "mechanism_pareto"), "south|shared", "FOLDX_DISPLAY_REVIEW"),
    ]
    assert store.completed == [
        (("north", "mechanism_pareto"), "north|backup"),
        (("south", "mechanism_pareto"), "south|backup"),
    ]
    # The shared physical signature ran once (3 + 7 slots); only the next
    # distinct signature was scheduled after the two independent rechecks.
    assert physical_signatures == [
        (("A1D",), 0, 3), (("A1D",), 3, 7),
        (("A2D",), 0, 3),
    ]

    candidates = list(validation.materialized_candidates.values())
    for item in candidates:
        apply_display_qualification(item, item.impact_evidence.get("stability_repeat_audit"))
    assert not north_shared.display_eligible and not south_shared.display_eligible
    assert north_backup.display_eligible and south_backup.display_eligible
    sectors = [V4Sector("north", [1, 2], 1, 20, 20, 0),
               V4Sector("south", [1, 2], 2, 20, 20, 0)]
    _panel, selected = select_panel(candidates, max_primary=4)
    assert {item.candidate_id for item in selected} == {"north|backup", "south|backup"}
    recommendations = build_sector_recommendations(
        candidates, sectors, {1: 10, 2: 10}, {}, selected)
    for row in recommendations:
        serialized_ids = {row.get("candidate_id"), *[
            item.get("candidate_id") for item in row.get("alternatives", [])]}
        assert not {"north|shared", "south|shared"} & serialized_ids


def test_display_review_without_replacement_leaves_queue_gap_not_a_completion():
    only = candidate("sensitive", "A1D")
    store = Store({("s", "sector_coverage"): [only]})

    def foldx(batch, *, slot_start, runs):
        values = [2.4, 1.0, 1.0] if slot_start == 0 else [1.0] * runs
        return PredictionResult("success", records=[pair_row(batch[0], values, slot_start)])

    def recheck(_item, result):
        result["display_review_reasons"] = ["SAMPLING_DIAGNOSTIC_UNKNOWN"]
        return "DISPLAY_REVIEW"

    result = SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx,
                                                foldx_recheck=recheck)
    assert result["foldx"]["results"][("A1D",)]["status"] == "PASS"
    assert store.completed == []
    assert store.outcomes == [(("s", "sector_coverage"), "sensitive", "FOLDX_DISPLAY_REVIEW")]
    assert result["foldx"]["results"][("A1D",)]["candidate_rechecks"]["sensitive"]["display_status"] == "DISPLAY_REVIEW"


def test_late_lane_reuses_foldx_signature_without_charging_budget_twice():
    north_head = candidate("north|shared", "A1D")
    south_top = candidate("south|top", "A2D")
    south_late = candidate("south|shared", "A1D")
    north_head.region, south_top.region, south_late.region = "north", "south", "south"
    north_head.impact_evidence = south_top.impact_evidence = south_late.impact_evidence = {}
    store = CopyingStore({("north", "epitope_support"): [north_head],
                          ("south", "epitope_support"): [south_top, south_late]})
    foldx_calls, rechecked = [], []

    def predict(batch):
        return PredictionResult("success", records=[
            {"mutations": [item.mutations[0].notation],
             "ddg": 1.5 if item.candidate_id == south_top.candidate_id else .2}
            for item in batch
        ])

    def foldx(batch, *, slot_start, runs):
        foldx_calls.append([item.candidate_id for item in batch])
        return PredictionResult("success", records=[
            pair_row(item, [1.] * runs, slot_start) for item in batch
        ])

    def recheck(item, result):
        rechecked.append(item.candidate_id)
        item.stability = f"checked:{item.region}"
        item.impact_evidence["stability_repeat_audit"] = result["repeat_audit"]
        return item.region != "south"

    validation = SiteFirstValidation()
    result = validation.run_indexed(store, predict, foldx, foldx_slot_limit=20,
                                    foldx_recheck=recheck)
    assert foldx_calls == [["north|shared"]]
    assert rechecked == ["north|shared", "south|shared"]
    assert result["foldx"]["ledger"]["attempted_slots"] == 3
    assert result["foldx"]["ledger"]["reserved_budget"] == 10
    assert store.completed == [(('north', 'epitope_support'), 'north|shared')]
    assert store.outcomes == [(('south', 'epitope_support'), 'south|top', 'SPURS_REVIEW'),
                              (('south', 'epitope_support'), 'south|shared', 'FOLDX_SCIENTIFIC_FAILURE')]
    assert validation.materialized_candidates["north|shared"].stability == "checked:north"
    assert validation.materialized_candidates["south|shared"].stability == "checked:south"
    gates = result["foldx"]["results"][("A1D",)]["candidate_rechecks"]
    assert gates["north|shared"]["status"] == "PASS"
    assert gates["south|shared"]["status"] == "SCIENTIFIC_FAILURE"


def test_spurs_budget_limited_lane_does_not_block_other_lane_cached_head():
    a1, a2 = candidate("a1", "A1D"), candidate("a2", "A2D")
    b1 = candidate("b1", "A3D")
    store = Store({("a", "coverage"): [a1, a2], ("b", "coverage"): [b1]})
    scored, foldx_calls = [], []

    def predict(batch):
        scored.extend(item.candidate_id for item in batch)
        return spurs(batch)

    def foldx(batch, *, slot_start, runs):
        foldx_calls.extend(c.candidate_id for c in batch)
        return PredictionResult("success", records=[
            pair_row(c, [2.1, 2.2, .5] if c is a1 else [1.] * runs, slot_start)
            for c in batch])

    result = SiteFirstValidation().run_indexed(store, predict, foldx,
        spurs_request_limit=2, foldx_slot_limit=20)
    assert set(scored) == {"a1", "b1"}
    assert foldx_calls == ["a1", "b1"]
    assert result["spurs"]["ledger"]["unique_requests"] == 2
    assert result["foldx"]["ledger"]["attempted_slots"] == 6
    assert store.completed == [(('b', 'coverage'), 'b1')]
    assert result["lane_status"][("a", "coverage")] == "STOPPED_BY_BUDGET"


def test_less_than_ten_foldx_slots_dispatches_nothing():
    a = candidate("a", "A1D")
    store = Store({("s", "coverage"): [a]})
    calls = []
    result = SiteFirstValidation().run_indexed(
        store, lambda batch: spurs(batch), lambda *a, **k: calls.append(1), foldx_slot_limit=9)
    assert not calls
    assert result["foldx"]["ledger"]["attempted_slots"] == 0
    assert store.states[("s", "coverage")] == "budget_limited"


def test_spurs_unique_budget_counts_dispatched_failures_without_retry():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("s1", "coverage"): [a], ("s2", "coverage"): [b]})
    calls = []

    def predict(batch):
        calls.append([item.candidate_id for item in batch])
        return PredictionResult("failed", errors=[{"code": "worker_error"}])

    result = SiteFirstValidation().run_indexed(store, predict, lambda *a, **k: None,
                                                spurs_request_limit=1)
    assert sum(map(len, calls)) == result["spurs"]["ledger"]["unique_requests"] == 1
    assert result["spurs"]["ledger"]["failed_requests"] == 1
    assert result["foldx"]["ledger"]["attempted_slots"] == 0


def test_spurs_partial_failure_preserves_good_row_and_stops_before_foldx():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("s1", "coverage"): [a], ("s2", "coverage"): [b]})
    foldx_calls = []

    def predict(batch):
        return PredictionResult("partial", records=[
            {"mutations": ["A1D"], "ddg": .2},
            {"mutations": ["A2D"], "status": "failed", "errors": [{"code": "timeout"}]},
        ], errors=[{"code": "one_candidate_failed"}])

    result = SiteFirstValidation().run_indexed(
        store, predict, lambda *args, **kwargs: foldx_calls.append(args))
    assert not foldx_calls
    assert result["spurs"]["results"][("A1D",)]["status"] == "PASS"
    assert result["spurs"]["results"][("A2D",)]["status"] == "TECHNICAL_FAILURE"
    assert result["foldx"]["ledger"]["attempted_slots"] == 0
    assert all(value == "STOPPED_TECHNICAL_FAILURE" for value in result["lane_status"].values())


def test_spurs_prefetch_batches_are_sequential_and_capped_at_sixty_four():
    queues = {(f"s{i}", "coverage"): [candidate(f"c{i}-{j}", f"A{(i*4+j)+1}D")
                                        for j in range(4)] for i in range(17)}
    store = Store(queues)
    batch_sizes, concurrent = [], 0

    def predict(batch):
        nonlocal concurrent
        concurrent += 1
        assert concurrent == 1
        batch_sizes.append(len(batch))
        result = spurs(batch, {item.candidate_id: 2.5 for item in batch})
        concurrent -= 1
        return result

    result = SiteFirstValidation().run_indexed(store, predict, lambda *a, **k: None)
    assert batch_sizes[:2] == [64, 4]
    assert max(batch_sizes) <= 64
    assert sum(batch_sizes) == result["spurs"]["ledger"]["unique_requests"] == 68


def test_released_reservation_allows_the_next_head_to_reserve_ten():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    store = Store({("s1", "coverage"): [a], ("s2", "coverage"): [b]})
    calls = []

    def foldx(batch, *, slot_start, runs):
        calls.append([c.candidate_id for c in batch])
        return PredictionResult("success", records=[pair_row(c, [1.] * runs, slot_start) for c in batch])

    result = SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx,
                                                foldx_slot_limit=13)
    assert calls == [["a"], ["b"]]
    assert result["foldx"]["ledger"]["attempted_slots"] == 6
    assert result["foldx"]["ledger"]["released_slots"] == 14



def test_refill_64_consumes_cached_replacements_before_another_worker_call():
    items = [candidate(str(i), f"A{i + 1}D") for i in range(130)]
    store = Store({("s", "coverage"): items})
    calls, physical = [], []

    def predict(batch):
        calls.append([c.candidate_id for c in batch])
        return spurs(batch, {str(i): 1.5 for i in range(129)})

    def foldx(batch, *, slot_start, runs):
        physical.extend(c.candidate_id for c in batch)
        return PredictionResult("success", records=[pair_row(c, [1.] * runs, slot_start) for c in batch])

    result = SiteFirstValidation().run_indexed(store, predict, foldx)
    assert [len(batch) for batch in calls] == [64, 64, 2]
    assert [identifier for batch in calls for identifier in batch] == [str(i) for i in range(130)]
    assert physical == ["129"]
    assert [row[1] for row in store.outcomes] == [str(i) for i in range(129)]
    assert result["spurs"]["ledger"]["unique_requests"] == 130
    assert result["foldx"]["ledger"]["attempted_slots"] == 3


def test_refill_deduplicates_shared_mutations_and_prioritizes_other_lane_head():
    first = [candidate(str(i), f"A{i + 1}D") for i in range(70)]
    other = [candidate("shared", "A1D"), candidate("other", "A200D")]
    store = Store({("a", "coverage"): first, ("b", "coverage"): other})
    calls = []

    def predict(batch):
        calls.append([c.mutations[0].notation for c in batch])
        return spurs(batch, {c.candidate_id: 1.5 for c in batch})

    result = SiteFirstValidation().run_indexed(store, predict, lambda *a, **k: None)
    flat = [s for batch in calls for s in batch]
    assert [len(batch) for batch in calls] == [64, 7]
    assert len(flat) == len(set(flat)) == 71
    assert calls[0][:3] == ["A1D", "A2D", "A200D"]
    assert result["foldx"]["ledger"]["attempted_slots"] == 0


def test_cached_head_completion_does_not_refill_unneeded_tail():
    items = [candidate(str(i), f"A{i + 1}D") for i in range(130)]
    store = Store({("s", "coverage"): items})
    calls = []

    def predict(batch):
        calls.append(len(batch))
        return spurs(batch, {"0": 1.5})

    def foldx(batch, *, slot_start, runs):
        return PredictionResult("success", records=[pair_row(c, [1.] * runs, slot_start) for c in batch])

    result = SiteFirstValidation().run_indexed(store, predict, foldx)
    assert calls == [64]
    assert store.completed == [(("s", "coverage"), "1")]
    assert result["spurs"]["ledger"]["unique_requests"] == 64


def test_refill_partial_last_budget_batch_never_dispatches_1001st_request():
    items = [candidate(str(i), f"A{i + 1}D") for i in range(1020)]
    store = Store({("s", "coverage"): items})
    calls = []

    def predict(batch):
        calls.append(len(batch))
        return spurs(batch, {c.candidate_id: 1.5 for c in batch})

    result = SiteFirstValidation().run_indexed(store, predict, lambda *a, **k: None)
    assert calls == [64] * 15 + [40]
    assert result["spurs"]["ledger"]["unique_requests"] == 1000
    assert result["lane_status"][("s", "coverage")] == "STOPPED_BY_BUDGET"


def test_more_than_one_batch_of_lane_heads_remains_pending_until_scored():
    store = Store({(str(i), "coverage"): [candidate(str(i), f"A{i + 1}D")] for i in range(65)})
    calls = []

    def predict(batch):
        calls.append(len(batch))
        return spurs(batch)

    def foldx(batch, *, slot_start, runs):
        return PredictionResult("success", records=[pair_row(c, [1.] * runs, slot_start) for c in batch])

    result = SiteFirstValidation().run_indexed(store, predict, foldx)
    assert calls == [64, 1]
    assert len(store.completed) == 65
    assert set(result["lane_status"].values()) == {"TARGET_REACHED"}
    assert result["spurs"]["ledger"]["refill_policy"] == "head_miss_round_robin_v1"


def _prefilter_candidate(identifier, mutation, missing=False):
    row = candidate(identifier, mutation)
    row.qualification = 'FORMAL_POOL'
    row.impact_evidence = {'policy': 'regional_mechanism_pareto_v1',
                           'missing_axes': ['charge_delta', 'donor_delta', 'acceptor_delta'] if missing else []}
    for option in row.mutations:
        option.review_required = False
    return row


def test_known_review_prefilter_preserves_budget_and_replenishes():
    bad = _prefilter_candidate('histidine_review', 'H1A', True)
    good = _prefilter_candidate('replacement', 'A2D')
    calls = []
    def foldx(batch, *, slot_start, runs):
        calls.append([c.candidate_id for c in batch])
        return PredictionResult('success', records=[pair_row(c, [.1] * runs, slot_start) for c in batch])
    scored = []
    def predict(batch):
        scored.extend(c.candidate_id for c in batch)
        return spurs(batch)
    store = Store({('s', 'sector_coverage'): [bad, good]})
    result = SiteFirstValidation().run_indexed(store, predict, foldx, foldx_slot_limit=10,
                                              display_prefilter='known_review_v1')
    assert scored == ['replacement'] and calls == [['replacement']]
    assert store.completed == [(('s', 'sector_coverage'), 'replacement')]
    assert result['foldx']['ledger']['attempted_slots'] == 3
    review = result['spurs']['ledger']['display_prefilter']['candidates']['histidine_review']
    assert review['scientific_status'] == 'NOT_ASSESSED'
    assert review['physical_scoring_requested'] is False
    assert store.outcomes[0][2] == 'DISPLAY_REVIEW_BEFORE_SCORING'


def test_prefilter_does_not_exclude_same_mutations_in_another_context():
    bad = _prefilter_candidate('bad_context', 'H1A', True)
    good = _prefilter_candidate('good_context', 'H1A')
    store = Store({('s1', 'sector_coverage'): [bad], ('s2', 'mechanism_pareto'): [good]})
    calls=[]
    def foldx(batch, *, slot_start, runs):
        calls.append([c.candidate_id for c in batch])
        return PredictionResult('success', records=[pair_row(c, [.1] * runs, slot_start) for c in batch])
    result=SiteFirstValidation().run_indexed(store, lambda batch: spurs(batch), foldx,
                                            display_prefilter='known_review_v1')
    assert calls == [['good_context']]
    assert store.completed == [(('s2', 'mechanism_pareto'), 'good_context')]
    assert len(result['foldx']['results']) == 1
