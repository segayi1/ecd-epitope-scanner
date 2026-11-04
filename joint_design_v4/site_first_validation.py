"""Site-first batched scheduler for the ordinary V4 validation protocol."""
from __future__ import annotations

import statistics
from collections.abc import Mapping
from copy import deepcopy
from time import perf_counter
from typing import Any, Callable

from .models import PredictionResult
from .ordinary_adaptive_validation import (
    AdaptiveValidation, FOLDX_ATTEMPT_LIMIT, SPURS_BATCH_SIZE,
    SPURS_REQUEST_LIMIT, _candidate_id, _mutations, _records, _row_signature,
    _score, _status, _errors,
)

PROTOCOL = "ordinary_site_first_v2"


def _candidate_site_key(candidate):
    evidence = getattr(candidate, "impact_evidence", {}) or {}
    site_set_id = evidence.get("site_set_id")
    if site_set_id:
        return str(site_set_id)
    positions = sorted(
        mutation.residue.uniprot_pos
        for mutation in (getattr(candidate, "mutations", ()) or ())
        if getattr(getattr(mutation, "residue", None), "uniprot_pos", None) is not None
    )
    return "positions:" + ",".join(map(str, positions)) if positions else (
        "candidate:" + str(getattr(candidate, "candidate_id", "unknown"))
    )


def _candidate_positions(candidate):
    return {mutation.residue.uniprot_pos
            for mutation in (getattr(candidate, "mutations", ()) or ())
            if getattr(getattr(mutation, "residue", None), "uniprot_pos", None) is not None}


class SiteFirstValidation(AdaptiveValidation):
    """Batch validation while preserving each lane's current site-first head."""

    def _init_ledgers(self, spurs_limit, foldx_limit):
        # The dataclass base has no post-init hook; initialize protocol-specific fields per run.
        self.spurs_ledger.update(policy=PROTOCOL, batch_size=SPURS_BATCH_SIZE,
                                 refill_policy="head_miss_round_robin_v1",
                                 unique_requests=0, failed_requests=0,
                                 technical_failure=False, batches=[], stop_reason=None,
                                 timings={"seconds": 0.0}, model_batches=0, valid=0,
                                 failed=0, reserved_budget=0, request_limit=spurs_limit)
        self.foldx_ledger.update(policy=PROTOCOL, slot_limit=foldx_limit,
                                 reserved_slots=0, attempted_slots=0, valid_slots=0,
                                 failed_slots=0, released_slots=0, technical_failure=False,
                                 budget_exhausted=False, batches=[], stop_reason=None,
                                 timings={"first_pass_seconds": 0.0, "append_seconds": 0.0},
                                 model_batches=0, reserved_budget=0,
                                 unseen_site_priority="distinct_site_first_v1",
                                 unseen_site_reserve_candidates=1,
                                 unseen_site_candidates=0, attempted_site_keys=[])

    @staticmethod
    def _recheck(callback, candidate, result):
        if callback is None:
            return True
        verdict = callback(candidate, result)
        if verdict is False:
            result["status"] = "SCIENTIFIC_FAILURE"
            result["recheck"] = "SCIENTIFIC_FAILURE"
            return False
        return verdict

    def _candidate_gate(self, result, candidate, callback):
        candidate_id = candidate.candidate_id
        gates = result.setdefault("candidate_rechecks", {})
        if candidate_id not in gates:
            checked = deepcopy(result)
            if checked["status"] == "PASS":
                verdict = self._recheck(callback, candidate, checked)
                if verdict == "DISPLAY_REVIEW":
                    checked["display_status"] = "DISPLAY_REVIEW"
            gates[candidate_id] = {
                "status": checked["status"],
                "display_status": checked.get("display_status", "DISPLAY_ELIGIBLE"),
                "display_review_reasons": list(checked.get("display_review_reasons") or []),
                "formal_exclusion_reasons": list(checked.get("formal_exclusion_reasons") or []),
                "repeat_audit": deepcopy(checked.get("repeat_audit")),
                "provenance": deepcopy(checked.get("provenance", [])),
            }
        gate = gates[candidate_id]
        return "DISPLAY_REVIEW" if gate.get("display_status") == "DISPLAY_REVIEW" else gate["status"]

    def _score_spurs_batch(self, candidates, predict, recheck, progress=None, audit=None):
        ledger = self.spurs_ledger
        unseen = {}
        for candidate in candidates:
            signature = _mutations(candidate)
            unseen.setdefault(signature, candidate)
        signatures = list(unseen)
        if not signatures:
            return
        initial_requests = ledger["unique_requests"]
        batches = [signatures[i:i + SPURS_BATCH_SIZE]
                   for i in range(0, len(signatures), SPURS_BATCH_SIZE)]
        for batch_sigs in batches:
            started = perf_counter()
            batch = [unseen[sig] for sig in batch_sigs]
            ledger["unique_requests"] += len(batch_sigs)
            ledger["reserved_budget"] += len(batch_sigs)
            try:
                response = predict(batch)
            except Exception as exc:
                response = PredictionResult("failed", errors=[{"code": "CALL_RAISED", "message": str(exc)}])
            rows = {}
            invalid = False
            expected = set(batch_sigs)
            for row in _records(response):
                sig = _row_signature(row)
                if sig not in expected or sig in rows:
                    invalid = True
                else:
                    rows[sig] = row
            failed_status = _status(response) not in {"success", "partial"}
            row_failures = {sig for sig, row in rows.items()
                            if str(row.get("status", "success")).lower() in {"failed", "error"}}
            missing_rows = expected.difference(rows)
            global_failure = failed_status or invalid
            provenance = getattr(response, "provenance", {}) or {}
            if isinstance(response, Mapping):
                provenance = response.get("provenance", provenance) or {}
            batch_audit = {"mutations": [list(s) for s in batch_sigs],
                           "status": _status(response), "errors": _errors(response),
                           "invalid_identity": invalid, "provenance": dict(provenance)}
            ledger["batches"].append(batch_audit)
            if audit is not None:
                audit("chunk", batch_audit)
            ledger["model_batches"] += 1
            if global_failure:
                ledger["technical_failure"] = True
                ledger["failed_requests"] += len(batch_sigs)
                ledger["failed"] += len(batch_sigs)
                ledger["stop_reason"] = "SPURS_TECHNICAL_FAILURE"
            elif _errors(response) or row_failures:
                ledger["technical_failure"] = True
                ledger["stop_reason"] = "SPURS_TECHNICAL_FAILURE"
            elif missing_rows:
                ledger["technical_failure"] = True
                ledger["stop_reason"] = "SPURS_MISSING_RESULT"
            for sig in batch_sigs:
                row = rows.get(sig)
                score, state = _score(row) if row else (None, "MISSING_RESULT")
                row_failed = sig in row_failures
                if row_failed:
                    score, state = None, "TECHNICAL_FAILURE"
                if score is None:
                    status = "TECHNICAL_FAILURE" if global_failure or state in {"TECHNICAL_FAILURE", "INVALID_SCORE"} else state
                    if not global_failure:
                        ledger["failed_requests"] += 1
                        ledger["failed"] += 1
                    if state == "INVALID_SCORE":
                        ledger["technical_failure"] = True
                        ledger["stop_reason"] = "SPURS_INVALID_SCORE"
                elif global_failure:
                    status = "BLOCKED_UNSCORED"
                elif row_failed:
                    status = "TECHNICAL_FAILURE"
                else:
                    status = "PASS" if score <= 1 else "REVIEW" if score <= 2 else "OVER_LIMIT"
                    ledger["valid"] += 1
                value = {"status": status, "score": score,
                         "candidate_id": _candidate_id(unseen[sig]),
                         "mutations": list(sig), "row": dict(row) if row else None}
                self.spurs_results[sig] = value
                if status == "PASS":
                    self._recheck(recheck, unseen[sig], value)
            ledger["timings"]["seconds"] += perf_counter() - started
            if progress is not None:
                progress(min(len(signatures), ledger["unique_requests"] - initial_requests), len(signatures))
            if ledger["technical_failure"]:
                break

    def _foldx_batch_stage(self, candidates, slot_start, runs, predict, max_workers):
        """Run one multi-candidate FoldX batch and map unordered rows by signature."""
        ledger = self.foldx_ledger
        started = perf_counter()
        try:
            response = predict(candidates, slot_start=slot_start, runs=runs)
        except Exception as exc:
            response = PredictionResult("failed", errors=[{"code": "CALL_RAISED", "message": str(exc)}])
        elapsed = perf_counter() - started
        timer_key = "first_pass_seconds" if slot_start == 0 else "append_seconds"
        ledger["timings"][timer_key] += elapsed
        ledger["model_batches"] += 1
        expected = {_mutations(c): c for c in candidates}
        indexed = {}
        rows_by_signature = {}
        duplicate_signatures = set()
        identity_errors = []
        for row in _records(response):
            sig = _row_signature(row)
            if sig is None:
                identity_errors.append({"code": "FOLDX_MISSING_SIGNATURE", "row": dict(row)})
            elif sig not in expected:
                identity_errors.append({"code": "FOLDX_UNKNOWN_SIGNATURE", "mutations": list(sig), "row": dict(row)})
            elif sig in indexed or sig in duplicate_signatures:
                duplicate_signatures.add(sig)
                indexed.pop(sig, None)
                rows_by_signature.setdefault(sig, []).append(row)
                identity_errors.append({"code": "FOLDX_DUPLICATE_SIGNATURE", "mutations": list(sig), "row": dict(row)})
            else:
                indexed[sig] = [row]
                rows_by_signature[sig] = [row]
        failed_status = _status(response) not in {"success", "partial"}
        errors = [*_errors(response), *identity_errors]
        results = {}
        for sig, candidate in expected.items():
            rows = indexed.get(sig, [])
            duplicate = sig in duplicate_signatures
            row = rows[0] if rows else None
            row_failed = (duplicate or (row is not None and
                          str(row.get("status", "success")).lower() in {"failed", "error"}))
            item_response = PredictionResult(
                "failed" if failed_status or row_failed else "success",
                records=rows_by_signature.get(sig, []) if duplicate else rows, errors=errors)
            provenance = (row or {}).get("provenance", {}) or {}
            if not provenance and len(candidates) == 1:
                provenance = getattr(response, "provenance", {}) or {}
                if isinstance(response, dict):
                    provenance = response.get("provenance", {}) or {}
            item_response.provenance = {"attempted_slots": min(runs, int(provenance.get("attempted_slots", runs)))}
            scores, attempted, failed, error = self._foldx_stage(item_response, sig, slot_start, runs)
            if row is None and not duplicate:
                error = True
            audit_rows = rows_by_signature.get(sig, []) if duplicate else rows
            results[sig] = (audit_rows, scores, attempted, failed,
                            error or failed_status or row_failed)
        stop_dispatch = failed_status or bool(errors) or any(item[4] for item in results.values())
        provenance = getattr(response, "provenance", {}) or {}
        if isinstance(response, Mapping):
            provenance = response.get("provenance", provenance) or {}
        return results, _status(response), errors, bool(identity_errors), stop_dispatch, dict(provenance)

    def _record_foldx_stage(self, candidates, audits, slot_start, runs, predict, max_workers):
        batch_results, status, errors, invalid, stop_dispatch, provenance = self._foldx_batch_stage(
            candidates, slot_start, runs, predict, max_workers)
        ledger = self.foldx_ledger
        if stop_dispatch:
            ledger["technical_failure"] = True
            ledger["stop_reason"] = "FOLDX_TECHNICAL_FAILURE"
        for candidate in candidates:
            sig = _mutations(candidate)
            row, scores, attempted, failed, technical = batch_results[sig]
            audit = audits[sig]
            records = row if isinstance(row, list) else [row] if row else []
            stage_provenance = {**provenance, "attempted_slots": attempted}
            audit["stages"].append({"slot_start": slot_start, "requested_slots": runs,
                                    "status": status, "scores": scores, "records": records,
                                    "provenance": stage_provenance, "errors": errors})
            self._account_stage(audit, ledger, attempted, len(scores), failed)
            if technical:
                ledger["technical_failure"] = True
                ledger["stop_reason"] = "FOLDX_TECHNICAL_FAILURE"
            yield sig, (row, scores, technical)

    def run_indexed(self, store, predict: Callable, foldx_predict: Callable, *,
                    spurs_request_limit=SPURS_REQUEST_LIMIT,
                    foldx_slot_limit=FOLDX_ATTEMPT_LIMIT, spurs_recheck=None,
                    foldx_recheck=None, max_workers=8, display_prefilter=None,
                    spurs_prefetch_policy=None, spurs_prefetch_window=1024, spurs_progress=None,
                    spurs_audit=None):
        if display_prefilter not in {None, "known_review_v1"}:
            raise ValueError("Unknown site-first display prefilter")
        if spurs_prefetch_policy not in {None, "ranked_window_v1"}:
            raise ValueError("Unknown SPURS prefetch policy")
        ranked = spurs_prefetch_policy == "ranked_window_v1"
        if ranked and not 1 <= int(spurs_prefetch_window) <= 1024:
            raise ValueError("SPURS prefetch window must be between 1 and 1024")
        if ranked and not hasattr(store, "prefetch_groups"):
            raise ValueError("Ranked prefetch requires site-first queue support")
        from .site_first_display import known_display_review_reasons
        spurs_limit = min(max(0, int(spurs_request_limit)), 10000 if ranked else SPURS_REQUEST_LIMIT)
        foldx_limit = min(max(0, int(foldx_slot_limit)), FOLDX_ATTEMPT_LIMIT)
        workers = min(8, max(1, int(max_workers)))
        self._init_ledgers(spurs_limit, foldx_limit)
        self.spurs_ledger.update(prefetch_policy=spurs_prefetch_policy,
                                 prefetch_window=int(spurs_prefetch_window) if ranked else None,
                                 windows=[])
        consumed = set()
        self.spurs_ledger["display_prefilter"] = {"policy": display_prefilter, "candidates": {}}
        lanes = [tuple(lane) for lane in store.queue_states()
                 if store.status(*lane)["status"] != "not_applicable"]
        lane_status = {}
        lane_heads = {}
        explored_positions = set()
        self.foldx_ledger["exploration_policy"] = "uncovered_residue_before_retry_v1"
        foldx = {"passed": {}, "lane_status": {}, "ledger": self.foldx_ledger,
                 "results": self.foldx_results}

        while len(lane_status) < len(lanes):
            if hasattr(store, "prioritize_uncovered"):
                store.prioritize_uncovered(explored_positions)
            lane_heads.clear()
            for lane in lanes:
                if lane in lane_status:
                    continue
                candidates = store.head(*lane, limit=1)
                if not candidates:
                    lane_status[lane] = "QUEUE_EXHAUSTED"
                    store.mark_queue_state(*lane, "candidate_failures")
                    continue
                current = candidates[0]
                if display_prefilter:
                    reasons = known_display_review_reasons(current)
                    if reasons:
                        self.materialized_candidates.setdefault(current.candidate_id, current)
                        current.display_eligible = False
                        current.display_exclusion_reasons = reasons
                        self.spurs_ledger["display_prefilter"]["candidates"][current.candidate_id] = {
                            "reasons": reasons, "missing_axes": list((current.impact_evidence or {}).get("missing_axes", [])),
                            "scientific_status": "NOT_ASSESSED", "physical_scoring_requested": False}
                        store.record_outcome(*lane, current.candidate_id, "DISPLAY_REVIEW_BEFORE_SCORING")
                        continue
                lane_heads[lane] = [self.materialized_candidates.setdefault(
                    candidate.candidate_id, candidate) for candidate in candidates]
            remaining_spurs = spurs_limit - self.spurs_ledger["unique_requests"]
            if remaining_spurs > 0 and any(
                _mutations(candidates[0]) not in self.spurs_results
                for candidates in lane_heads.values()
            ):
                # Refill only on a cache miss; round-robin keeps every lane's head ahead of tails.
                # Preserve SPURS capacity for hidden site groups before prefetching retries.
                unseen = {}
                offered = 0
                batch_limit = min(int(spurs_prefetch_window) if ranked else SPURS_BATCH_SIZE, remaining_spurs)
                def offer(candidate):
                    nonlocal offered
                    if display_prefilter and known_display_review_reasons(candidate):
                        return
                    offered += 1
                    sig = _mutations(candidate)
                    if sig not in self.spurs_results and len(unseen) < batch_limit:
                        unseen.setdefault(sig, candidate)

                for candidates in lane_heads.values():
                    offer(candidates[0])
                if ranked:
                    iterators = [iter(store.prefetch_groups(*lane)) for lane in lane_heads]
                    tails = []
                    while iterators and len(unseen) < batch_limit:
                        active = []
                        for iterator in iterators:
                            group = next(iterator, None)
                            if group:
                                offer(group[0])
                                tails.append(iter(group[1:]))
                                active.append(iterator)
                        iterators = active
                    while tails and len(unseen) < batch_limit:
                        active = []
                        for iterator in tails:
                            candidate = next(iterator, None)
                            if candidate is not None:
                                offer(candidate)
                                active.append(iterator)
                        tails = active
                else:
                    exploring = (hasattr(store, "needs_exploration")
                                 and store.needs_exploration(explored_positions))
                    frontier_limit = 1 if exploring else SPURS_BATCH_SIZE
                    frontiers = [store.head(*lane, limit=frontier_limit) for lane in lane_heads]
                    for offset in range(max(map(len, frontiers), default=0)):
                        for frontier in frontiers:
                            if offset < len(frontier):
                                offer(frontier[offset])
                batch = list(unseen.values())
                before = self.spurs_ledger["unique_requests"]
                batch_start = len(self.spurs_ledger["batches"])
                window = {"candidates_considered": offered, "unique_candidates": len(batch),
                          "duplicates_or_already_scored": offered - len(batch),
                          "budget_before": remaining_spurs, "status": "running"}
                self.spurs_ledger["windows"].append(window)
                if ranked and spurs_audit is not None:
                    spurs_audit("window_started", dict(window))
                self._score_spurs_batch(batch, predict, spurs_recheck, spurs_progress if ranked else None,
                                        spurs_audit if ranked else None)
                audits = self.spurs_ledger["batches"][batch_start:]
                window.update(requests=self.spurs_ledger["unique_requests"] - before,
                              cache_hit_batches=sum(bool(a.get("provenance", {}).get("cache_hit")) for a in audits),
                              budget_after=spurs_limit - self.spurs_ledger["unique_requests"],
                              status="failed" if self.spurs_ledger["technical_failure"] else "completed")
                if ranked and spurs_audit is not None:
                    spurs_audit("window_finished", dict(window))
            passed_heads = {}
            by_signature = {}
            blocked = []
            for lane, candidates in lane_heads.items():
                current = candidates[0]
                sig = _mutations(current)
                outcome = self.spurs_results.get(sig)
                if outcome is not None:
                    consumed.add(sig)
                if outcome is None:
                    blocked.append(lane)
                elif outcome["status"] == "PASS":
                    passed_heads.setdefault(sig, current)
                    by_signature.setdefault(sig, []).append((lane, current))
                elif outcome["status"] in {"TECHNICAL_FAILURE", "MISSING_RESULT", "INVALID_SCORE", "BLOCKED_UNSCORED"}:
                    store.mark_queue_state(*lane, "technical_failure" if outcome["status"] in {"TECHNICAL_FAILURE", "BLOCKED_UNSCORED"} else "candidate_failures")
                    lane_status[lane] = "BLOCKED_UNSCORED"
                else:
                    explored_positions.update(_candidate_positions(current))
                    store.record_outcome(*lane, current.candidate_id, "SPURS_" + outcome["status"])
            if self.spurs_ledger["technical_failure"]:
                store.mark_all_technical_failure()
                lane_status.update({lane: "STOPPED_TECHNICAL_FAILURE" for lane in lane_heads})
                break
            if blocked:
                if self.spurs_ledger["unique_requests"] >= spurs_limit:
                    for lane in blocked:
                        store.mark_queue_state(*lane, "budget_limited", spurs_limit)
                        lane_status[lane] = "STOPPED_BY_BUDGET"

            if not passed_heads:
                continue
            cached = [sig for sig in passed_heads if sig in self.foldx_results]
            for sig in cached:
                physical = self.foldx_results[sig]
                for lane, candidate in by_signature[sig]:
                    explored_positions.update(_candidate_positions(candidate))
                    self.foldx_ledger["attempted_site_keys"] = sorted(
                        set(self.foldx_ledger["attempted_site_keys"])
                        | {_candidate_site_key(candidate)}
                    )
                    status = self._candidate_gate(physical, candidate, foldx_recheck)
                    if status == "PASS":
                        store.complete(*lane, candidate.candidate_id)
                        lane_status[lane] = "TARGET_REACHED"
                    elif status == "DISPLAY_REVIEW":
                        store.record_outcome(*lane, candidate.candidate_id, "FOLDX_DISPLAY_REVIEW")
                    elif status == "SCIENTIFIC_FAILURE":
                        store.record_outcome(*lane, candidate.candidate_id, "FOLDX_SCIENTIFIC_FAILURE")
                    elif status in {"TECHNICAL_FAILURE", "MISSING_RESULT"}:
                        store.mark_queue_state(*lane, "technical_failure")
                        lane_status[lane] = "BLOCKED_UNSCORED"
                del passed_heads[sig]
                del by_signature[sig]
            if self.foldx_ledger["technical_failure"]:
                store.mark_all_technical_failure()
                for lane in lane_heads:
                    lane_status.setdefault(lane, "STOPPED_TECHNICAL_FAILURE")
                break
            if not passed_heads:
                continue
            remaining = foldx_limit - self.foldx_ledger["attempted_slots"] - self.foldx_ledger["reserved_slots"]
            slots = min(workers, max(0, remaining // 10))
            attempted_sites = set(self.foldx_ledger["attempted_site_keys"])
            preferred, deferred = [], []
            # Fill each batch with untried site sets before revisiting a lane.
            for item in passed_heads.items():
                key = _candidate_site_key(item[1])
                if key not in attempted_sites:
                    attempted_sites.add(key)
                    preferred.append(item)
                else:
                    deferred.append(item)
            selected = (preferred + deferred)[:slots]
            if not selected:
                self.foldx_ledger["budget_exhausted"] = True
                self.foldx_ledger["stop_reason"] = "FOLDX_BUDGET_INSUFFICIENT"
                for sig, members in by_signature.items():
                    for lane, candidate in members:
                        store.mark_queue_state(*lane, "budget_limited", foldx_limit)
                        lane_status[lane] = "STOPPED_BY_BUDGET"
                break
            selected_candidates = [c for _, c in selected]
            for candidate in selected_candidates:
                explored_positions.update(_candidate_positions(candidate))
            self.foldx_ledger["explored_positions"] = sorted(explored_positions)
            dispatched_sites = {_candidate_site_key(c) for c in selected_candidates}
            new_sites = dispatched_sites.difference(self.foldx_ledger["attempted_site_keys"])
            self.foldx_ledger["unseen_site_candidates"] += len(new_sites)
            self.foldx_ledger["attempted_site_keys"] = sorted(
                set(self.foldx_ledger["attempted_site_keys"]) | dispatched_sites
            )
            audits = {}
            for sig, candidate in selected:
                self.foldx_ledger["reserved_slots"] += 10
                self.foldx_ledger["reserved_budget"] += 10
                audits[sig] = {"mutations": list(sig), "candidate_id": _candidate_id(candidate),
                               "reserved_slots": 10, "attempted_slots": 0, "valid_slots": 0,
                               "failed_slots": 0, "released_slots": 0, "stages": []}
            first = dict(self._record_foldx_stage(selected_candidates, audits, 0, 3,
                                                   foldx_predict, workers))
            append = []
            scores_by_sig = {}
            for sig, candidate in selected:
                _, scores, technical = first[sig]
                scores_by_sig[sig] = list(scores)
                if len(scores) != 3 and not technical:
                    technical = True
                    self.foldx_ledger["technical_failure"] = True
                    self.foldx_ledger["stop_reason"] = "FOLDX_MISSING_RESULT"
                if technical:
                    status = "TECHNICAL_FAILURE" if self.foldx_ledger["technical_failure"] else "MISSING_RESULT"
                    self.foldx_results[sig] = self._set_foldx_result(sig, candidate, status, scores, audits[sig])
                elif sum(value > 2 for _, value in scores) == 0:
                    self.foldx_results[sig] = self._set_foldx_result(sig, candidate, "PASS", scores, audits[sig])
                elif sum(value > 2 for _, value in scores) >= 2:
                    self.foldx_results[sig] = self._set_foldx_result(sig, candidate, "SCIENTIFIC_FAILURE", scores, audits[sig])
                else:
                    append.append((sig, candidate))
            if append and not self.foldx_ledger["technical_failure"]:
                append_result = dict(self._record_foldx_stage([c for _, c in append], audits, 3, 7,
                                                               foldx_predict, workers))
                for sig, candidate in append:
                    _, scores, technical = append_result[sig]
                    scores_by_sig[sig].extend(scores)
                    status = "TECHNICAL_FAILURE" if technical else (
                        "MISSING_RESULT" if len(scores_by_sig[sig]) != 10 else
                        "PASS" if statistics.median(v for _, v in scores_by_sig[sig]) <= 2 else "SCIENTIFIC_FAILURE")
                    self.foldx_results[sig] = self._set_foldx_result(sig, candidate, status,
                                                                     scores_by_sig[sig], audits[sig])
            elif append:
                for sig, candidate in append:
                    self.foldx_results[sig] = self._set_foldx_result(
                        sig, candidate, "TECHNICAL_FAILURE", scores_by_sig[sig], audits[sig])
            for sig, candidate in selected:
                result = self.foldx_results[sig]
                result["repeat_audit"]["execution_protocol"] = PROTOCOL
                audit = audits[sig]
                self.foldx_ledger["reserved_slots"] -= 10
                audit["released_slots"] = max(0, 10 - audit["attempted_slots"])
                self.foldx_ledger["released_slots"] += audit["released_slots"]
                self.foldx_ledger["batches"].append(audit)
                for lane, lane_candidate in by_signature[sig]:
                    lane_status_result = self._candidate_gate(result, lane_candidate, foldx_recheck)
                    if lane_status_result == "PASS":
                        store.complete(*lane, lane_candidate.candidate_id)
                        lane_status[lane] = "TARGET_REACHED"
                    elif lane_status_result == "DISPLAY_REVIEW":
                        store.record_outcome(*lane, lane_candidate.candidate_id, "FOLDX_DISPLAY_REVIEW")
                    elif lane_status_result == "SCIENTIFIC_FAILURE":
                        store.record_outcome(*lane, lane_candidate.candidate_id, "FOLDX_SCIENTIFIC_FAILURE")
                    elif lane_status_result in {"TECHNICAL_FAILURE", "MISSING_RESULT"}:
                        store.mark_queue_state(*lane, "technical_failure")
                        lane_status[lane] = "BLOCKED_UNSCORED"
            if self.foldx_ledger["technical_failure"]:
                store.mark_all_technical_failure()
                for lane in lane_heads:
                    lane_status.setdefault(lane, "STOPPED_TECHNICAL_FAILURE")
                break

        for lane in lanes:
            lane_status.setdefault(lane, "STOPPED")
        if not self.spurs_ledger["stop_reason"]:
            self.spurs_ledger["stop_reason"] = "QUEUES_EXHAUSTED"
        if not self.foldx_ledger["stop_reason"]:
            self.foldx_ledger["stop_reason"] = "QUEUES_EXHAUSTED"
        self.spurs_ledger["consumed_signatures"] = len(consumed)
        self.spurs_ledger["prefetched_unused"] = len(set(self.spurs_results) - consumed)
        foldx["passed"] = {sig: result for sig, result in self.foldx_results.items()
                            if result["status"] == "PASS"}
        foldx["lane_status"] = dict(lane_status)
        result = {"spurs": {"ledger": self.spurs_ledger, "results": self.spurs_results},
                  "foldx": foldx, "lane_status": lane_status}
        result["ledger"] = store.ledger(self)
        return result

# Case boundary: user requested the remaining cases use the next tested runtime.
