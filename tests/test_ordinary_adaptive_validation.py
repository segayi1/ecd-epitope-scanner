from types import SimpleNamespace as NS

import pytest
from ecd_epitope_scanner.joint_design_v4.models import PredictionResult
from ecd_epitope_scanner.joint_design_v4.ordinary_adaptive_validation import AdaptiveValidation


def candidate(identifier, *mutations):
    return NS(candidate_id=identifier, mutations=[NS(notation=value) for value in mutations])


def spurs_result(rows, status="success", errors=None):
    return PredictionResult(status, records=rows, errors=errors or [])


def foldx_result(mutations, values, slot_start=0, status="success", attempted=None, errors=None):
    return PredictionResult(status, records=[{
        "mutations": list(mutations),
        "paired_models": [
            {"sampling_slot": slot_start + index, "ddg_kcal_mol": value,
             "batch_path": f"batch_{slot_start}", "group": 1,
             "run": index + 1, "source_path": f"source_{slot_start}_{index + 1}"}
            for index, value in enumerate(values)
        ],
    }], errors=errors or [], provenance={"attempted_slots": attempted if attempted is not None else len(values),
                                           "model_revision": "fixture-r1"})


class FakeIndex:
    def __init__(self, queues):
        self.queues = {lane: list(items) for lane, items in queues.items()}
        self.states = {lane: "ready" for lane in queues}
        self.outcomes = []
        self.completed = []

    def queue_states(self):
        return list(self.queues)

    def status(self, *lane):
        return {"status": self.states[tuple(lane)]}

    def queue_candidate_id(self, *lane):
        items = self.queues[tuple(lane)]
        return items[0].candidate_id if items else None

    def materialize(self, candidate_id):
        return next(item for items in self.queues.values() for item in items
                    if item.candidate_id == candidate_id)

    def record_outcome(self, *args):
        lane, candidate_id, outcome = args[:2], args[2], args[3]
        items = self.queues[tuple(lane)]
        assert items[0].candidate_id == candidate_id
        self.outcomes.append((tuple(lane), items.pop(0).candidate_id, outcome))

    def complete(self, *args):
        lane, candidate_id = args[:2], args[2]
        assert self.queues[tuple(lane)][0].candidate_id == candidate_id
        self.completed.append((tuple(lane), self.queues[tuple(lane)].pop(0).candidate_id))
        self.states[tuple(lane)] = "complete"

    def mark_queue_state(self, *args):
        lane, state = args[:2], args[2]
        self.states[tuple(lane)] = state

    def mark_all_technical_failure(self):
        for lane, state in self.states.items():
            if state == "ready":
                self.states[lane] = "technical_failure"

    def mark_all_budget_limited(self, limit):
        for lane, state in self.states.items():
            if state == "ready":
                self.states[lane] = "budget_limited"

    def ledger(self, validation):
        return {"outcomes": self.outcomes, "completed": self.completed}


def test_spurs_batches_at_64_and_deduplicates_shared_queue_head():
    rows = [candidate(f"c{i}", f"A{i + 1}D") for i in range(65)]
    queues = {("r", f"p{i}"): [item] for i, item in enumerate(rows)}
    queues[("shared-a", "epitope")] = [rows[0]]
    queues[("shared-b", "coverage")] = [rows[0]]
    calls = []

    def predict(batch):
        calls.append([c.candidate_id for c in batch])
        return spurs_result([{"mutations": [c.mutations[0].notation], "ddg": 2.1} for c in batch])

    result = AdaptiveValidation().run_spurs(queues, predict)
    assert [len(batch) for batch in calls] == [64, 1]
    assert result["ledger"]["unique_requests"] == 65
    assert all(status == "QUEUE_EXHAUSTED" for status in result["lane_status"].values())


def test_spurs_classifies_bands_and_does_not_skip_missing_queue_head():
    top = candidate("top", "A1D")
    lower = candidate("lower", "A2D")
    calls = []

    def predict(batch):
        calls.append([item.candidate_id for item in batch])
        return spurs_result([{"mutations": ["A1D"], "ddg": 1.4}]) if batch[0] is top else spurs_result([])

    result = AdaptiveValidation().run_spurs({"q": [top, lower]}, predict)
    assert result["results"][("A1D",)]["status"] == "REVIEW"
    assert result["results"][("A2D",)]["status"] == "MISSING_RESULT"
    assert result["lane_status"]["q"] == "BLOCKED_UNSCORED"
    assert calls == [["top"], ["lower"]]


def test_spurs_technical_failure_preserves_score_without_qualifying_candidate():
    a, b = candidate("a", "A1D"), candidate("b", "A2D")
    result = AdaptiveValidation().run_spurs({"q": [a, b]}, lambda batch: spurs_result(
        [{"mutations": ["A1D"], "ddg": .2}], status="failed", errors=[{"code": "worker"}]))
    assert not result["passed"]
    assert result["ledger"]["unique_requests"] == result["ledger"]["failed_requests"] == 1
    assert result["ledger"]["stop_reason"] == "SPURS_TECHNICAL_FAILURE"
    assert result["results"][("A1D",)]["status"] == "BLOCKED_UNSCORED"
    assert result["results"][("A1D",)]["score"] == .2


def test_indexed_scheduler_interleaves_spurs_and_foldx_before_advancing_passed_head():
    a1, a2 = candidate("a1", "A1D"), candidate("a2", "A2D")
    b1, b2 = candidate("b1", "A3D"), candidate("b2", "A4D")
    store = FakeIndex({("a", "epitope_support"): [a1, a2],
                       ("b", "mechanism_pareto"): [b1, b2]})
    events = []

    def spurs(batch):
        events.append(("spurs", [item.candidate_id for item in batch]))
        return spurs_result([{"mutations": [item.mutations[0].notation],
                              "ddg": 0.2 if item is a1 or item is b2 else 1.5}
                             for item in batch])

    def foldx(item, *, slot_start, runs):
        events.append(("foldx", item.candidate_id))
        values = [1.0] * runs if item is a1 or item is b2 else [2.1, 2.2, 2.3]
        return foldx_result([x.notation for x in item.mutations], values, slot_start)

    result = AdaptiveValidation().run_indexed(store, spurs, foldx)
    assert events == [("spurs", ["a1", "b1"]), ("foldx", "a1"),
                      ("spurs", ["b2"]), ("foldx", "b2")]
    assert store.completed == [(('a', 'epitope_support'), "a1"),
                               (('b', 'mechanism_pareto'), "b2")]
    assert result["ledger"]["outcomes"] == [(("b", "mechanism_pareto"), "b1", "SPURS_REVIEW")]


def test_indexed_shared_result_waits_for_lane_head_and_is_scored_once():
    shared = candidate("shared", "A1D")
    top = candidate("top", "A2D")
    store = FakeIndex({("a", "epitope_support"): [shared],
                       ("b", "mechanism_pareto"): [top, shared]})
    calls = []

    def spurs(batch):
        calls.append(("spurs", [item.candidate_id for item in batch]))
        return spurs_result([{"mutations": [item.mutations[0].notation],
                              "ddg": 0.2 if item is shared else 1.5}
                             for item in batch])

    def foldx(item, *, slot_start, runs):
        calls.append(("foldx", item.candidate_id))
        return foldx_result([x.notation for x in item.mutations], [1.0] * runs, slot_start)

    AdaptiveValidation().run_indexed(store, spurs, foldx)
    assert calls == [("spurs", ["shared", "top"]), ("foldx", "shared")]
    assert store.completed == [(('a', 'epitope_support'), "shared"),
                               (('b', 'mechanism_pareto'), "shared")]
    assert store.outcomes == [(('b', 'mechanism_pareto'), "top", "SPURS_REVIEW")]


def test_indexed_spurs_limit_allows_exact_budget_across_rounds():
    first, second = candidate("first", "A1D"), candidate("second", "A2D")
    store = FakeIndex({("q", "sector_coverage"): [first, second]})
    spurs_calls = []
    foldx_calls = []

    def spurs(batch):
        spurs_calls.append(batch[0].candidate_id)
        return spurs_result([{"mutations": [batch[0].mutations[0].notation],
                              "ddg": 1.5 if batch[0] is first else 0.2}])

    def foldx(item, *, slot_start, runs):
        foldx_calls.append(item.candidate_id)
        return foldx_result([option.notation for option in item.mutations],
                            [0.5] * runs, slot_start)

    validation = AdaptiveValidation()
    result = validation.run_indexed(store, spurs, foldx, spurs_request_limit=2)
    assert spurs_calls == ["first", "second"]
    assert foldx_calls == ["second"]
    assert result["spurs"]["ledger"]["unique_requests"] == 2
    assert store.completed == [(("q", "sector_coverage"), "second")]


def test_indexed_technical_failure_does_not_advance_or_call_foldx():
    first, later = candidate("first", "A1D"), candidate("later", "A2D")
    store = FakeIndex({("q", "sector_coverage"): [first, later]})
    calls = []

    def spurs(batch):
        calls.append("spurs")
        return spurs_result([{"mutations": ["A1D"], "ddg": .2}], status="failed")

    def foldx(*args, **kwargs):
        calls.append("foldx")
        pytest.fail("technical SPURS failure must stop downstream scoring")

    result = AdaptiveValidation().run_indexed(store, spurs, foldx)
    assert calls == ["spurs"]
    assert store.outcomes == []
    assert store.states[("q", "sector_coverage")] == "technical_failure"
    assert result["lane_status"][("q", "sector_coverage")] == "STOPPED_TECHNICAL_FAILURE"


def test_spurs_request_limit_counts_dispatched_failures_without_refund():
    items = [candidate(str(i), f"A{i + 1}D") for i in range(3)]
    calls = []

    def predict(batch):
        calls.append(batch)
        return spurs_result([], status="failed")

    result = AdaptiveValidation().run_spurs({str(i): [item] for i, item in enumerate(items)},
        predict, request_limit=2)
    assert [len(batch) for batch in calls] == [2]
    assert result["ledger"]["unique_requests"] == result["ledger"]["failed_requests"] == 2
    assert result["ledger"]["stop_reason"] == "SPURS_TECHNICAL_FAILURE"


def test_spurs_configuration_cannot_exceed_task_or_batch_caps():
    items = [candidate(str(i), f"A{i + 1}D") for i in range(1001)]
    calls = []

    def predict(batch):
        calls.append(len(batch))
        return spurs_result([{"mutations": [item.mutations[0].notation], "ddg": 3.0}
                             for item in batch])

    result = AdaptiveValidation().run_spurs(
        {str(i): [item] for i, item in enumerate(items)}, predict,
        request_limit=2000, batch_size=1000)
    assert result["ledger"]["unique_requests"] == sum(calls) == 1000
    assert max(calls) == 64
    assert len(result["results"]) == 1000


def test_foldx_three_slots_pass_and_shared_result_is_not_repeated_or_skipped():
    first, lower = candidate("first", "A1D", "A2V"), candidate("lower", "A3D")
    calls = []

    def predict(item, *, slot_start, runs):
        calls.append((item.candidate_id, slot_start, runs))
        return foldx_result([x.notation for x in item.mutations], [1.0] * runs, slot_start)

    result = AdaptiveValidation().run_foldx({"e": [first, lower], "m": [first]}, {("A1D", "A2V"): first}, predict)
    assert calls == [("first", 0, 3)]
    assert result["lane_status"] == {"e": "TARGET_REACHED", "m": "TARGET_REACHED"}
    assert result["ledger"]["attempted_slots"] == result["ledger"]["valid_slots"] == 3
    assert result["ledger"]["released_slots"] == 7
    assert result["results"][("A1D", "A2V")]["repeat_audit"]["execution_protocol"] == (
        "ordinary_three_principle_adaptive_v1")


def test_foldx_two_high_slots_are_scientific_failure_without_append():
    item = candidate("c", "A1D", "A2V")
    calls = []

    def predict(c, *, slot_start, runs):
        calls.append((slot_start, runs))
        return foldx_result([x.notation for x in c.mutations], [2.1, 2.2, .5], slot_start)

    result = AdaptiveValidation().run_foldx({"q": [item]}, {("A1D", "A2V"): item}, predict)
    assert calls == [(0, 3)]
    assert result["results"][("A1D", "A2V")]["status"] == "SCIENTIFIC_FAILURE"
    assert result["ledger"]["released_slots"] == 7


@pytest.mark.parametrize("tail,expected", [([2.1] + [1.0] * 6, "PASS"), ([2.1] * 7, "SCIENTIFIC_FAILURE")])
def test_foldx_mixed_first_pass_runs_one_seven_slot_append_and_uses_all_ten(tail, expected):
    item = candidate("c", "A1D", "A2V")
    calls = []

    def predict(c, *, slot_start, runs):
        calls.append((slot_start, runs))
        values = [1.0, 1.1, 2.1] if slot_start == 0 else tail
        return foldx_result([x.notation for x in c.mutations], values, slot_start)

    result = AdaptiveValidation().run_foldx({"q": [item]}, {("A1D", "A2V"): item}, predict)
    assert calls == [(0, 3), (3, 7)]
    assert result["results"][("A1D", "A2V")]["status"] == expected
    assert len(result["results"][("A1D", "A2V")]["scores"]) == 10
    assert result["ledger"]["attempted_slots"] == result["ledger"]["valid_slots"] == 10
    stored = result["results"][("A1D", "A2V")]
    assert [pair["sampling_slot"] for pair in stored["paired_models"]] == list(range(10))
    assert stored["paired_models"][0]["batch_path"] == "batch_0"
    assert stored["paired_models"][-1]["source_path"] == "source_3_7"
    assert stored["provenance"] == [
        {"attempted_slots": 3, "model_revision": "fixture-r1"},
        {"attempted_slots": 7, "model_revision": "fixture-r1"},
    ]


def test_foldx_reservation_requires_ten_remaining_slots_and_stops_without_science_rejection():
    item = candidate("c", "A1D", "A2V")
    calls = []
    result = AdaptiveValidation().run_foldx({"q": [item]}, {("A1D", "A2V"): item},
        lambda *args, **kwargs: calls.append(1), slot_limit=9)
    assert calls == []
    assert result["results"][("A1D", "A2V")]["status"] == "BUDGET_INSUFFICIENT"
    assert result["lane_status"]["q"] == "STOPPED_BY_BUDGET"
    assert result["ledger"]["attempted_slots"] == 0
    assert result["ledger"]["budget_exhausted"] is True


def test_foldx_does_not_skip_higher_ranked_spurs_unverified_head():
    uncertain, lower = candidate("uncertain", "A1D", "A2V"), candidate("lower", "A3D", "A4V")
    validation = AdaptiveValidation()
    validation.spurs_results[("A1D", "A2V")] = {"status": "MISSING_RESULT"}
    calls = []
    result = validation.run_foldx({"q": [uncertain, lower]}, {("A3D", "A4V"): lower},
        lambda c, **kwargs: calls.append(c))
    assert calls == []
    assert result["lane_status"]["q"] == "BLOCKED_UNSCORED"


def test_foldx_technical_failure_preserves_valid_partial_slots_and_stops_later_candidates():
    first, later = candidate("first", "A1D", "A2V"), candidate("later", "A3D", "A4V")
    calls = []

    def predict(item, *, slot_start, runs):
        calls.append(item.candidate_id)
        return foldx_result([x.notation for x in item.mutations], [1.0], slot_start,
                            status="failed", attempted=2, errors=[{"code": "timeout"}])

    result = AdaptiveValidation().run_foldx({"q": [first, later]}, {("A1D", "A2V"): first,
        ("A3D", "A4V"): later}, predict)
    assert calls == ["first"]
    assert result["results"][("A1D", "A2V")]["scores"] == [1.0]
    assert result["results"][("A1D", "A2V")]["status"] == "TECHNICAL_FAILURE"
    assert result["ledger"]["attempted_slots"] == 2
    assert result["ledger"]["valid_slots"] == 1
    assert result["ledger"]["failed_slots"] == 1
    assert result["lane_status"]["q"] == "STOPPED_TECHNICAL_FAILURE"


def test_foldx_budget_gap_is_not_reported_as_candidate_failure():
    item = candidate("c", "A1D", "A2V")
    validation = AdaptiveValidation()
    validation.spurs_results[("A1D", "A2V")] = {"status": "PASS"}
    result = validation.run_foldx(
        {"q": [item]}, {("A1D", "A2V"): item},
        lambda *args, **kwargs: pytest.fail("budget gap must not call FoldX"),
        slot_limit=9,
    )
    assert result["results"][("A1D", "A2V")]["status"] == "BUDGET_INSUFFICIENT"
    assert result["lane_status"]["q"] == "STOPPED_BY_BUDGET"
