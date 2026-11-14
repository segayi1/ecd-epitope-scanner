"""Pure adaptive SPURS/FoldX scheduling for ordinary V4 validation.

Queue keys identify a (region, principle) lane. Each lane advances only from
its current head. Model callables return ``PredictionResult`` and can be
replaced with deterministic fixtures in tests.
"""
from __future__ import annotations

import math
import statistics
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

from .models import PredictionResult

SPURS_BATCH_SIZE = 64
SPURS_REQUEST_LIMIT = 1000
FOLDX_ATTEMPT_LIMIT = 1000


def _mutations(candidate: Any) -> tuple[str, ...]:
    values = candidate.get("mutations", ()) if isinstance(candidate, Mapping) else getattr(candidate, "mutations", ())
    result = []
    for value in values:
        if isinstance(value, str):
            result.append(value.strip().upper())
        else:
            result.append(str(getattr(value, "notation", value)).strip().upper())
    if not result:
        raise ValueError("candidate has no mutation identity")
    return tuple(sorted(result))


def _candidate_id(candidate: Any) -> str:
    value = candidate.get("candidate_id") if isinstance(candidate, Mapping) else getattr(candidate, "candidate_id", None)
    return str(value if value is not None else ",".join(_mutations(candidate)))


def _records(result: Any) -> list[dict[str, Any]]:
    if isinstance(result, PredictionResult):
        return result.records
    if isinstance(result, Mapping):
        return list(result.get("records", []))
    return list(getattr(result, "records", []))


def _status(result: Any) -> str:
    if isinstance(result, Mapping):
        return str(result.get("status", "failed")).lower()
    return str(getattr(result, "status", "failed")).lower()


def _errors(result: Any) -> list[Any]:
    if isinstance(result, Mapping):
        return list(result.get("errors", []))
    return list(getattr(result, "errors", []))


def _score(row: Mapping[str, Any]) -> tuple[float | None, str]:
    if not any(key in row for key in ("ddg", "ddg_kcal_mol", "spurs_multi_ddg_kcal_mol")):
        return None, "MISSING_RESULT"
    value = row.get("ddg", row.get("ddg_kcal_mol", row.get("spurs_multi_ddg_kcal_mol")))
    try:
        score = float(value)
    except (TypeError, ValueError, OverflowError):
        return None, "INVALID_SCORE"
    return (score, "SCORED") if math.isfinite(score) else (None, "INVALID_SCORE")


def _row_signature(row: Mapping[str, Any]) -> tuple[str, ...] | None:
    values = row.get("mutations", row.get("mutation_list"))
    if isinstance(values, str):
        values = values.split(",")
    if not isinstance(values, (list, tuple)) or not values:
        return None
    return tuple(sorted(str(value).strip().upper() for value in values))


def _queues(queues: Mapping[Any, Iterable[Any]]) -> dict[Any, list[Any]]:
    return {key: list(values) for key, values in queues.items()}


@dataclass
class AdaptiveValidation:
    """State and audit ledger shared by the two model stages."""

    spurs_results: dict[tuple[str, ...], dict[str, Any]] = field(default_factory=dict)
    foldx_results: dict[tuple[str, ...], dict[str, Any]] = field(default_factory=dict)
    materialized_candidates: dict[str, Any] = field(default_factory=dict)
    spurs_ledger: dict[str, Any] = field(default_factory=lambda: {
        "policy": "ordinary_three_principle_adaptive_v1", "unique_requests": 0,
        "failed_requests": 0, "technical_failure": False, "batches": [],
        "stop_reason": None,
    })
    foldx_ledger: dict[str, Any] = field(default_factory=lambda: {
        "policy": "ordinary_three_principle_adaptive_v1", "slot_limit": FOLDX_ATTEMPT_LIMIT,
        "reserved_slots": 0, "attempted_slots": 0, "valid_slots": 0,
        "failed_slots": 0, "released_slots": 0, "technical_failure": False,
        "budget_exhausted": False, "batches": [], "stop_reason": None,
    })

    def run_spurs(self, queues: Mapping[Any, Iterable[Any]], predict: Callable,
                  recheck: Callable | None = None, *, request_limit: int = SPURS_REQUEST_LIMIT,
                  batch_size: int = SPURS_BATCH_SIZE) -> dict[str, Any]:
        """Score lane heads in deduplicated batches; only scores <=1 pass."""
        if request_limit < 0 or batch_size < 1:
            raise ValueError("invalid SPURS limits")
        request_limit = min(request_limit, SPURS_REQUEST_LIMIT)
        batch_size = min(batch_size, SPURS_BATCH_SIZE)
        lanes = _queues(queues)
        offsets = {key: 0 for key in lanes}
        finished: dict[Any, str] = {}
        ledger = self.spurs_ledger
        ledger["stop_reason"] = None

        while len(finished) < len(lanes) and not ledger["technical_failure"]:
            heads: dict[tuple[str, ...], Any] = {}
            for lane, items in lanes.items():
                if lane in finished:
                    continue
                while offsets[lane] < len(items):
                    candidate = items[offsets[lane]]
                    signature = _mutations(candidate)
                    prior = self.spurs_results.get(signature)
                    if prior is None:
                        heads.setdefault(signature, candidate)
                        break
                    if prior["status"] == "PASS":
                        finished[lane] = "TARGET_REACHED"
                        break
                    if prior["status"] in {"MISSING_RESULT", "INVALID_SCORE", "TECHNICAL_FAILURE"}:
                        finished[lane] = "BLOCKED_UNSCORED"
                        break
                    offsets[lane] += 1
                else:
                    finished[lane] = "QUEUE_EXHAUSTED"
            if len(finished) == len(lanes):
                break

            remaining = request_limit - ledger["unique_requests"]
            if remaining <= 0:
                ledger["stop_reason"] = "SPURS_UNIQUE_REQUEST_LIMIT"
                for lane in lanes:
                    if lane not in finished:
                        finished[lane] = "BUDGET_EXHAUSTED"
                break
            batch_signatures = list(heads)[:min(batch_size, remaining)]
            if not batch_signatures:
                ledger["stop_reason"] = "SPURS_NO_DISPATCHABLE_HEAD"
                break
            batch = [heads[signature] for signature in batch_signatures]
            ledger["unique_requests"] += len(batch_signatures)
            try:
                response = predict(batch)
            except Exception as exc:  # injected callable boundary is a worker boundary
                response = PredictionResult("failed", errors=[{"code": "CALL_RAISED", "message": str(exc)}])
            indexed: dict[tuple[str, ...], dict[str, Any]] = {}
            invalid_identity = False
            expected = set(batch_signatures)
            for row in _records(response):
                signature = _row_signature(row)
                if signature not in expected or signature in indexed:
                    invalid_identity = True
                    continue
                indexed[signature] = row
            failed_status = _status(response) not in {"success", "partial"}
            batch_audit = {"mutations": [list(signature) for signature in batch_signatures],
                           "status": _status(response), "errors": _errors(response),
                           "invalid_identity": invalid_identity}
            ledger["batches"].append(batch_audit)
            missing_signatures = [
                signature for signature in batch_signatures if signature not in indexed
            ]
            if failed_status:
                ledger["failed_requests"] += len(batch_signatures)
                ledger["technical_failure"] = True
                ledger["stop_reason"] = "SPURS_TECHNICAL_FAILURE"
            elif invalid_identity:
                ledger["failed_requests"] += len(batch_signatures)
                ledger["technical_failure"] = True
                ledger["stop_reason"] = "SPURS_INVALID_RESPONSE_IDENTITY"
            for signature in batch_signatures:
                row = indexed.get(signature)
                score, score_status = _score(row) if row is not None else (None, "MISSING_RESULT")
                if score is None:
                    if score_status == "INVALID_SCORE":
                        status = "TECHNICAL_FAILURE"
                        ledger["failed_requests"] += 1
                        ledger["technical_failure"] = True
                        ledger["stop_reason"] = "SPURS_INVALID_SCORE"
                    elif failed_status or invalid_identity:
                        status = "TECHNICAL_FAILURE"
                    else:
                        status = score_status
                        if score_status == "MISSING_RESULT":
                            ledger["failed_requests"] += 1
                elif failed_status or invalid_identity:
                    # Preserve a usable numeric observation for audit, but a
                    # technically failed request cannot qualify a candidate.
                    status = "BLOCKED_UNSCORED"
                elif score <= 1:
                    status = "PASS"
                elif score <= 2:
                    status = "REVIEW"
                else:
                    status = "OVER_LIMIT"
                self.spurs_results[signature] = {
                    "status": status, "score": score,
                    "candidate_id": _candidate_id(heads[signature]),
                    "mutations": list(signature), "row": dict(row) if row else None,
                }
                if recheck is not None and status == "PASS":
                    recheck(heads[signature], self.spurs_results[signature])
            # A partial batch may omit a requested head. Continue with other
            # complete rows, but never reinterpret the omission as a score.
            if missing_signatures and not (
                    failed_status or invalid_identity
                    or any(self.spurs_results[signature]["status"] == "TECHNICAL_FAILURE"
                           for signature in batch_signatures)):
                for lane in lanes:
                    if (lane not in finished and offsets[lane] < len(lanes[lane])
                            and _mutations(lanes[lane][offsets[lane]]) in missing_signatures):
                        finished[lane] = "BLOCKED_UNSCORED"
            elif ledger["technical_failure"]:
                for lane in lanes:
                    if (lane not in finished and offsets[lane] < len(lanes[lane])
                            and _mutations(lanes[lane][offsets[lane]]) in expected):
                        finished[lane] = "STOPPED_TECHNICAL_FAILURE"

        passed = {}
        for lane, items in lanes.items():
            for candidate in items:
                result = self.spurs_results.get(_mutations(candidate))
                if result and result["status"] == "PASS":
                    passed[_mutations(candidate)] = candidate
        for lane in lanes:
            if lane not in finished:
                finished[lane] = "STOPPED_TECHNICAL_FAILURE" if ledger["technical_failure"] else "STOPPED"
        return {"passed": passed, "lane_status": finished, "ledger": ledger,
                "results": self.spurs_results}

    def run_indexed(self, store, predict: Callable, foldx_predict: Callable,
                    *, spurs_request_limit: int = SPURS_REQUEST_LIMIT,
                    foldx_slot_limit: int = FOLDX_ATTEMPT_LIMIT,
                    spurs_recheck: Callable | None = None,
                    foldx_recheck: Callable | None = None) -> dict[str, Any]:
        """Run both stages over SQLite queues without materializing the pool.

        Every queue remains blocked at its own unverified head. A shared model
        result can satisfy several queues only after each reaches that item.
        """
        spurs_request_limit = min(spurs_request_limit, SPURS_REQUEST_LIMIT)
        foldx_slot_limit = min(foldx_slot_limit, FOLDX_ATTEMPT_LIMIT)
        lanes = [
            lane for lane in store.queue_states()
            if store.status(*lane)["status"] != "not_applicable"
        ]
        lane_status = {}
        foldx = {"passed": {}, "lane_status": {}, "ledger": self.foldx_ledger,
                 "results": self.foldx_results}
        while len(lane_status) < len(lanes):
            heads = {}
            candidates_by_lane = {}
            for lane in lanes:
                if lane in lane_status:
                    continue
                candidate_id = store.queue_candidate_id(*lane)
                if candidate_id is None:
                    lane_status[lane] = "QUEUE_EXHAUSTED"
                    store.mark_queue_state(*lane, "candidate_failures")
                    continue
                candidate = store.materialize(candidate_id)
                self.materialized_candidates[candidate.candidate_id] = candidate
                signature = _mutations(candidate)
                candidates_by_lane[lane] = (candidate_id, candidate, signature)
                heads.setdefault(lane, [candidate])

            if not candidates_by_lane:
                break
            spurs = self.run_spurs(
                heads, predict, spurs_recheck, request_limit=spurs_request_limit,
            )
            if self.spurs_ledger["technical_failure"]:
                store.mark_all_technical_failure()
                for lane in candidates_by_lane:
                    lane_status[lane] = "STOPPED_TECHNICAL_FAILURE"
                break

            passed_heads = {}
            passed_candidates = {}
            unrequested = []
            for lane, (candidate_id, candidate, signature) in candidates_by_lane.items():
                result = self.spurs_results.get(signature)
                if result is None:
                    unrequested.append(lane)
                elif result["status"] in {
                        "MISSING_RESULT", "INVALID_SCORE", "TECHNICAL_FAILURE", "BLOCKED_UNSCORED"}:
                    store.record_outcome(*lane, candidate_id,
                                         "SPURS_" + result["status"])
                    lane_status[lane] = "BLOCKED_UNSCORED"
                elif result["status"] == "PASS":
                    passed_heads[lane] = [candidate]
                    passed_candidates[signature] = candidate
                else:
                    store.record_outcome(*lane, candidate_id, "SPURS_" + result["status"])

            if unrequested:
                store.mark_all_budget_limited(spurs_request_limit)
                lane_status.update({lane: "STOPPED_BY_BUDGET" for lane in unrequested})
                break

            if not passed_heads:
                continue
            foldx = self.run_foldx(
                passed_heads, passed_candidates, foldx_predict, foldx_recheck,
                slot_limit=foldx_slot_limit,
            )
            for lane, (candidate_id, candidate, signature) in candidates_by_lane.items():
                if lane not in passed_heads:
                    continue
                status = foldx["lane_status"].get(lane)
                result = self.foldx_results.get(signature)
                if status == "TARGET_REACHED":
                    store.complete(*lane, candidate_id)
                    lane_status[lane] = "TARGET_REACHED"
                elif status == "QUEUE_EXHAUSTED":
                    store.record_outcome(*lane, candidate_id,
                                         "FOLDX_" + str((result or {}).get("status", "SCIENTIFIC_FAILURE")))
                elif status in {"STOPPED_BY_BUDGET", "BLOCKED_UNSCORED"}:
                    store.mark_queue_state(
                        *lane,
                        "technical_failure" if self.foldx_ledger["technical_failure"] else "budget_limited",
                        foldx_slot_limit if not self.foldx_ledger["technical_failure"] else None,
                    )
                    lane_status[lane] = status
                elif status == "STOPPED_TECHNICAL_FAILURE":
                    store.mark_queue_state(*lane, "technical_failure")
                    lane_status[lane] = status

            if self.foldx_ledger["technical_failure"]:
                store.mark_all_technical_failure()
                for lane in candidates_by_lane:
                    lane_status.setdefault(lane, "STOPPED_TECHNICAL_FAILURE")
                break
            if self.foldx_ledger["budget_exhausted"]:
                store.mark_all_budget_limited(foldx_slot_limit)
                for lane in candidates_by_lane:
                    lane_status.setdefault(lane, "STOPPED_BY_BUDGET")
                break

        if not self.spurs_ledger["stop_reason"]:
            self.spurs_ledger["stop_reason"] = "QUEUES_EXHAUSTED"
        if not foldx["passed"] and self.foldx_ledger["attempted_slots"] == 0:
            self.foldx_ledger["stop_reason"] = "NO_SPURS_PASSED_CANDIDATES"
        result = {"spurs": {"ledger": self.spurs_ledger, "results": self.spurs_results},
                  "foldx": foldx, "lane_status": lane_status}
        result["ledger"] = store.ledger(self)
        return result

    def run_foldx(self, queues: Mapping[Any, Iterable[Any]], spurs_passed: Mapping,
                  predict: Callable, recheck: Callable | None = None, *,
                  slot_limit: int = FOLDX_ATTEMPT_LIMIT) -> dict[str, Any]:
        """Evaluate lane heads using 3 then optionally 7 physical-pair slots."""
        if slot_limit < 0:
            raise ValueError("invalid FoldX slot limit")
        lanes = _queues(queues)
        offsets = {key: 0 for key in lanes}
        finished: dict[Any, str] = {}
        ledger = self.foldx_ledger
        ledger["slot_limit"] = min(FOLDX_ATTEMPT_LIMIT, slot_limit)
        ledger["stop_reason"] = None
        ledger["budget_exhausted"] = False
        passed_signatures = {_mutations(value) for value in spurs_passed.values()}
        lane_sequence = list(lanes)
        lane_cursor = 0

        while len(finished) < len(lanes) and not ledger["technical_failure"]:
            heads: dict[tuple[str, ...], Any] = {}
            ordered_lanes = lane_sequence[lane_cursor:] + lane_sequence[:lane_cursor]
            for lane in ordered_lanes:
                items = lanes[lane]
                if lane in finished:
                    continue
                while offsets[lane] < len(items):
                    candidate = items[offsets[lane]]
                    signature = _mutations(candidate)
                    if signature not in passed_signatures:
                        spurs_status = self.spurs_results.get(signature, {}).get("status")
                        if spurs_status in {None, "MISSING_RESULT", "INVALID_SCORE",
                                            "TECHNICAL_FAILURE", "BLOCKED_UNSCORED"}:
                            finished[lane] = "BLOCKED_UNSCORED"
                            break
                        offsets[lane] += 1
                        continue
                    prior = self.foldx_results.get(signature)
                    if prior is None:
                        heads.setdefault(signature, candidate)
                        break
                    if prior["status"] == "PASS":
                        finished[lane] = "TARGET_REACHED"
                        break
                    if prior["status"] == "BUDGET_INSUFFICIENT":
                        finished[lane] = "STOPPED_BY_BUDGET"
                        break
                    if prior["status"] in {"MISSING_RESULT", "INVALID_SLOTS", "TECHNICAL_FAILURE"}:
                        finished[lane] = "BLOCKED_UNSCORED"
                        break
                    offsets[lane] += 1
                else:
                    finished[lane] = "QUEUE_EXHAUSTED"
            if len(finished) == len(lanes):
                break
            if not heads:
                continue
            signature = next(iter(heads))
            candidate = heads[signature]
            selected_lane = next(lane for lane in ordered_lanes if lane not in finished and offsets[lane] < len(lanes[lane])
                                 and _mutations(lanes[lane][offsets[lane]]) == signature)
            lane_cursor = (lane_sequence.index(selected_lane) + 1) % len(lane_sequence)
            remaining = ledger["slot_limit"] - ledger["attempted_slots"] - ledger["reserved_slots"]
            if remaining < 10:
                ledger["budget_exhausted"] = True
                ledger["stop_reason"] = "FOLDX_BUDGET_INSUFFICIENT"
                self.foldx_results[signature] = {"status": "BUDGET_INSUFFICIENT", "mutations": list(signature),
                                                 "candidate_id": _candidate_id(candidate), "scores": []}
                continue

            ledger["reserved_slots"] += 10
            audit = {"mutations": list(signature), "candidate_id": _candidate_id(candidate),
                     "reserved_slots": 10, "attempted_slots": 0, "valid_slots": 0,
                     "failed_slots": 0, "released_slots": 0, "stages": []}
            all_scores: list[tuple[int, float]] = []
            try:
                response = predict(candidate, slot_start=0, runs=3)
            except Exception as exc:
                response = PredictionResult("failed", errors=[{"code": "CALL_RAISED", "message": str(exc)}])
            stage_scores, stage_attempted, stage_failed, stage_error = self._foldx_stage(response, signature, 0, 3)
            audit["stages"].append({"slot_start": 0, "requested_slots": 3, "status": _status(response),
                                    "scores": stage_scores, "records": _records(response),
                                    "provenance": dict(getattr(response, "provenance", {}) or {}),
                                    "errors": _errors(response)})
            self._account_stage(audit, ledger, stage_attempted, len(stage_scores), stage_failed)
            all_scores.extend(stage_scores)
            if stage_error:
                audit["released_slots"] = 10 - audit["attempted_slots"]
                ledger["reserved_slots"] -= 10
                ledger["released_slots"] += audit["released_slots"]
                self.foldx_results[signature] = self._set_foldx_result(signature, candidate, "TECHNICAL_FAILURE", all_scores, audit)
                ledger["technical_failure"] = True
                ledger["stop_reason"] = "FOLDX_TECHNICAL_FAILURE"
                ledger["batches"].append(audit)
                break

            if len(all_scores) != 3:
                audit["released_slots"] = 10 - audit["attempted_slots"]
                ledger["reserved_slots"] -= 10
                ledger["released_slots"] += audit["released_slots"]
                self.foldx_results[signature] = self._set_foldx_result(signature, candidate, "MISSING_RESULT", all_scores, audit)
                ledger["batches"].append(audit)
                for lane, items in lanes.items():
                    if lane not in finished and offsets[lane] < len(items) and _mutations(items[offsets[lane]]) == signature:
                        finished[lane] = "BLOCKED_UNSCORED"
                continue

            over = sum(score > 2 for _, score in all_scores)
            if over == 0:
                status = "PASS"
            elif over >= 2:
                status = "SCIENTIFIC_FAILURE"
            else:
                try:
                    response = predict(candidate, slot_start=3, runs=7)
                except Exception as exc:
                    response = PredictionResult("failed", errors=[{"code": "CALL_RAISED", "message": str(exc)}])
                extra_scores, extra_attempted, extra_failed, stage_error = self._foldx_stage(response, signature, 3, 7)
                audit["stages"].append({"slot_start": 3, "requested_slots": 7, "status": _status(response),
                                        "scores": extra_scores, "records": _records(response),
                                        "provenance": dict(getattr(response, "provenance", {}) or {}),
                                        "errors": _errors(response)})
                self._account_stage(audit, ledger, extra_attempted, len(extra_scores), extra_failed)
                all_scores.extend(extra_scores)
                if stage_error or len(all_scores) != 10:
                    status = "TECHNICAL_FAILURE" if stage_error else "MISSING_RESULT"
                    ledger["technical_failure"] = stage_error
                    if stage_error:
                        ledger["stop_reason"] = "FOLDX_TECHNICAL_FAILURE"
                else:
                    status = "PASS" if statistics.median(score for _, score in all_scores) <= 2 else "SCIENTIFIC_FAILURE"
            audit["released_slots"] = max(0, 10 - audit["attempted_slots"])
            ledger["reserved_slots"] -= 10
            ledger["released_slots"] += audit["released_slots"]
            self.foldx_results[signature] = self._set_foldx_result(signature, candidate, status, all_scores, audit)
            ledger["batches"].append(audit)
            if recheck is not None:
                recheck(candidate, self.foldx_results[signature])
            if status == "TECHNICAL_FAILURE":
                ledger["technical_failure"] = True
                ledger["stop_reason"] = "FOLDX_TECHNICAL_FAILURE"
                break

        for lane in lanes:
            if lane not in finished:
                finished[lane] = "STOPPED_TECHNICAL_FAILURE" if ledger["technical_failure"] else "STOPPED_BY_BUDGET"
        return {"passed": {signature: result for signature, result in self.foldx_results.items()
                            if result["status"] == "PASS"},
                "lane_status": finished, "ledger": ledger, "results": self.foldx_results}

    @staticmethod
    def _foldx_stage(response: Any, signature: tuple[str, ...], slot_start: int, runs: int):
        rows = [row for row in _records(response) if _row_signature(row) == signature]
        unexpected = any(_row_signature(row) not in {signature, None} for row in _records(response))
        failed_status = _status(response) not in {"success", "partial"}
        attempted = min(runs, int(getattr(response, "provenance", {}).get("attempted_slots", runs)))
        if len(rows) != 1:
            return [], attempted, runs, failed_status or unexpected
        paired = rows[0].get("paired_models")
        if not isinstance(paired, list):
            return [], attempted, runs, True
        slots: set[int] = set()
        scores = []
        invalid = False
        for pair in paired:
            try:
                slot = int(pair["sampling_slot"])
                value = float(pair["ddg_kcal_mol"])
            except (KeyError, TypeError, ValueError, OverflowError):
                invalid = True
                continue
            if slot < slot_start or slot >= slot_start + runs or slot in slots or not math.isfinite(value):
                invalid = True
                continue
            slots.add(slot)
            scores.append((slot, value))
        if invalid:
            return sorted(scores), attempted, max(1, attempted - len(scores)), True
        return sorted(scores), attempted, max(0, attempted - len(scores)), failed_status or unexpected

    @staticmethod
    def _account_stage(audit: dict[str, Any], ledger: dict[str, Any], attempted: int,
                       valid: int, failed: int) -> None:
        audit["attempted_slots"] += attempted
        audit["valid_slots"] += valid
        audit["failed_slots"] += failed
        ledger["attempted_slots"] += attempted
        ledger["valid_slots"] += valid
        ledger["failed_slots"] += failed

    @staticmethod
    def _set_foldx_result(signature, candidate, status, scores, audit):
        values = [value for _, value in scores]
        repeat_audit = {
            "execution_protocol": "ordinary_three_principle_adaptive_v1",
            "foldx_repeat_count": len(values),
            "foldx_repeat_slots": [slot for slot, _ in scores],
            "foldx_repeat_values_kcal_mol": values,
            "foldx_repeat_slot_values_kcal_mol": values,
            "foldx_evidence_complete": status == "PASS",
            "sampling_decision_diagnostic_status": (
                "COMPLETE_WITH_UNCERTAINTY" if status == "PASS" and len(values) == 10
                else "COMPLETE" if status == "PASS" else "INCOMPLETE"
            ),
        }
        records = [row for stage in audit["stages"] for row in stage.get("records", [])]
        paired_models = [pair for row in records for pair in row.get("paired_models", [])]
        return {"status": status, "candidate_id": _candidate_id(candidate),
                "mutations": list(signature), "scores": values,
                "median": statistics.median(values) if values else None,
                "slots": [slot for slot, _ in scores], "audit": audit,
                "repeat_audit": repeat_audit, "paired_models": paired_models,
                "records": records,
                "provenance": [stage.get("provenance", {}) for stage in audit["stages"]]}
