"""Task-local candidate index and independent ordinary V4 principle queues.

The SQLite database stores only candidate identity, ordering, labels, and queue
state. Mutation evidence remains in the caller's already computed single-site
pool and is materialized only for candidates that are returned from a queue.
"""
from __future__ import annotations

import json
import math
import sqlite3
import zlib
from copy import deepcopy
from dataclasses import replace
from itertools import combinations, product
from pathlib import Path
from types import SimpleNamespace

from .candidate_search import _coverage_from_distances, _distances
from .regional_perturbation import (
    _candidate,
    _option,
    _position,
    _span_ok,
    _valid_sequence,
    candidate_impact,
    complete_candidate_impact,
    dominates,
    rank_impacts,
)

PRINCIPLES = ("epitope_support", "mechanism_pareto", "sector_coverage")
POLICY = "ordinary_three_principle_adaptive_v1"


def _decode_metadata(value):
    return json.loads(zlib.decompress(value) if isinstance(value, bytes) else value)


class CandidateQueue:
    """SQLite queue index backed by the already computed single-site evidence."""

    def __init__(self, db_path=":memory:"):
        self.db_path = str(db_path)
        if self.db_path != ":memory:":
            Path(self.db_path).parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.db_path)
        self.db.row_factory = sqlite3.Row
        self.db.executescript("""
            CREATE TABLE IF NOT EXISTS candidates (
                candidate_id TEXT PRIMARY KEY, sector_id TEXT NOT NULL,
                signature TEXT NOT NULL, labels TEXT NOT NULL, metadata TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS queue (
                sector_id TEXT NOT NULL, principle TEXT NOT NULL,
                ordinal INTEGER NOT NULL, candidate_id TEXT NOT NULL,
                state TEXT NOT NULL DEFAULT 'pending',
                PRIMARY KEY (sector_id, principle, candidate_id),
                UNIQUE (sector_id, principle, ordinal)
            );
            CREATE TABLE IF NOT EXISTS queue_state (
                sector_id TEXT NOT NULL, principle TEXT NOT NULL,
                state TEXT NOT NULL DEFAULT 'ready', budget_limit INTEGER,
                PRIMARY KEY (sector_id, principle)
            );
        """)
        # Only single-site evidence and region definitions stay resident. Full
        # combination objects are materialized from the index on demand.
        self._options_by_position = {}
        self._sectors = {}
        self._reference_profiles_by_sector = {}
        self._geometry_by_sector = {}

    def close(self):
        self.db.close()

    def head(self, sector_id, principle, limit=1):
        if principle not in PRINCIPLES:
            raise ValueError(f"Unknown principle: {principle}")
        if limit < 0:
            raise ValueError("Queue head limit must be nonnegative")
        rows = self.db.execute(
            "SELECT c.candidate_id,c.signature,c.metadata FROM queue q "
            "JOIN candidates c USING(candidate_id) "
            "WHERE q.sector_id=? AND q.principle=? AND q.state='pending' "
            "ORDER BY q.ordinal LIMIT ?", (str(sector_id), principle, int(limit)))
        return [self._materialize(row) for row in rows]

    def complete(self, sector_id, principle, candidate_id):
        cursor = self.db.execute(
            "UPDATE queue SET state='complete' WHERE sector_id=? AND principle=? "
            "AND candidate_id=? AND state='pending'", (str(sector_id), principle, candidate_id))
        if cursor.rowcount != 1:
            raise KeyError(f"No pending queue item: {sector_id}/{principle}/{candidate_id}")
        self.db.execute(
            "UPDATE queue_state SET state='complete',budget_limit=NULL "
            "WHERE sector_id=? AND principle=?", (str(sector_id), principle))
        self.db.commit()

    def record_outcome(self, sector_id, principle, candidate_id, outcome):
        cursor = self.db.execute(
            "UPDATE queue SET state=? WHERE sector_id=? AND principle=? "
            "AND candidate_id=? AND state='pending'",
            (json.dumps({"outcome": str(outcome)}, sort_keys=True), str(sector_id),
             principle, str(candidate_id)))
        if cursor.rowcount != 1:
            raise KeyError(f"No pending queue item: {sector_id}/{principle}/{candidate_id}")
        self._refresh_state(sector_id, principle)
        self.db.commit()

    def mark_budget_limited(self, sector_id, principle, limit):
        self.db.execute("INSERT INTO queue_state(sector_id,principle,state,budget_limit) VALUES(?,?,?,?) "
                        "ON CONFLICT(sector_id,principle) DO UPDATE SET state='budget_limited',budget_limit=excluded.budget_limit",
                        (str(sector_id), principle, "budget_limited", int(limit)))
        self.db.commit()

    def mark_queue_state(self, sector_id, principle, state, limit=None):
        if state not in {
                "candidate_failures", "budget_limited", "technical_failure",
                "not_applicable"}:
            raise ValueError(f"Unknown queue state: {state}")
        current = self.db.execute(
            "SELECT state FROM queue_state WHERE sector_id=? AND principle=?",
            (str(sector_id), principle)).fetchone()
        if current is None or current["state"] in {
                "complete", "candidate_failures", "not_applicable"}:
            return
        self.db.execute(
            "UPDATE queue_state SET state=?,budget_limit=? "
            "WHERE sector_id=? AND principle=?",
            (state, int(limit) if limit is not None else None,
             str(sector_id), principle))
        self.db.commit()

    def mark_all_budget_limited(self, limit):
        self.db.execute(
            "UPDATE queue_state SET state='budget_limited',budget_limit=? "
            "WHERE state='ready'", (int(limit),))
        self.db.commit()

    def mark_all_technical_failure(self):
        self.db.execute(
            "UPDATE queue_state SET state='technical_failure' WHERE state='ready'")
        self.db.commit()

    def status(self, sector_id, principle):
        row = self.db.execute("SELECT state,budget_limit FROM queue_state WHERE sector_id=? AND principle=?",
                              (str(sector_id), principle)).fetchone()
        if row and row["state"] != "ready":
            result = {"status": row["state"]}
            if row["state"] == "budget_limited":
                result["budget_limit"] = row["budget_limit"]
            return result
        pending = self.db.execute("SELECT count(*) FROM queue WHERE sector_id=? AND principle=? AND state='pending'",
                                  (str(sector_id), principle)).fetchone()[0]
        if pending:
            return {"status": "ready", "pending": pending}
        exists = self.db.execute("SELECT count(*) FROM queue WHERE sector_id=? AND principle=?",
                                 (str(sector_id), principle)).fetchone()[0]
        return {"status": "exhausted" if exists else "empty", "pending": 0}

    def labels(self, candidate_id):
        row = self.db.execute("SELECT labels FROM candidates WHERE candidate_id=?", (candidate_id,)).fetchone()
        return json.loads(row[0]) if row else []

    def signatures(self):
        return [
            tuple(str(item).strip().upper() for item in json.loads(row[0]))
            for row in self.db.execute("SELECT signature FROM candidates ORDER BY candidate_id")
        ]

    def queue_candidate_ids(self, sector_id, principle):
        rows = self.db.execute(
            "SELECT candidate_id FROM queue WHERE sector_id=? AND principle=? "
            "AND state='pending' ORDER BY ordinal", (str(sector_id), principle))
        return [row[0] for row in rows]

    def queue_states(self):
        return [
            (row["sector_id"], row["principle"]) for row in self.db.execute(
                "SELECT sector_id,principle FROM queue_state ORDER BY sector_id,principle")
        ]

    def queue_candidate_id(self, sector_id, principle, *, offset=0):
        row = self.db.execute(
            "SELECT candidate_id FROM queue WHERE sector_id=? AND principle=? "
            "AND state='pending' ORDER BY ordinal LIMIT 1 OFFSET ?",
            (str(sector_id), principle, int(offset))).fetchone()
        return row[0] if row else None

    def materialize(self, candidate_id):
        row = self.db.execute(
            "SELECT candidate_id,signature,metadata FROM candidates WHERE candidate_id=?",
            (str(candidate_id),)).fetchone()
        if row is None:
            raise KeyError(f"Missing indexed candidate: {candidate_id}")
        return self._materialize(row)

    def ledger(self, validation=None):
        principles = list(PRINCIPLES)
        rows = []
        for sector_id, principle in self.db.execute(
                "SELECT sector_id,principle FROM queue_state ORDER BY sector_id,principle"):
            status = self.status(sector_id, principle)
            candidate_ids = self.queue_candidate_ids(sector_id, principle)
            rows.append({
                "sector_id": sector_id,
                "principle": principle,
                "state": status["status"],
                "pending_count": len(candidate_ids),
                "next_candidate_id": candidate_ids[0] if candidate_ids else None,
                "completed": [
                    {"candidate_id": row["candidate_id"],
                     "outcome": self._outcome(row["state"])}
                    for row in self.db.execute(
                        "SELECT candidate_id,state FROM queue WHERE sector_id=? "
                        "AND principle=? AND state!='pending' ORDER BY ordinal",
                        (sector_id, principle))
                ],
            })
        counts = {
            "legal_candidates": self.db.execute("SELECT count(*) FROM candidates").fetchone()[0],
            "queues": len(rows),
            "ready_queues": sum(row["state"] == "ready" for row in rows),
            "exhausted_queues": sum(row["state"] == "exhausted" for row in rows),
            "budget_limited_queues": sum(row["state"] == "budget_limited" for row in rows),
        }
        payload = {
            "schema_version": "ordinary_adaptive_candidate_ledger_v1",
            "policy": "ordinary_three_principle_adaptive_v1",
            "principles": principles,
            "counts": counts,
            "queues": rows,
            "spurs": dict(validation.spurs_ledger) if validation is not None else None,
            "foldx": dict(validation.foldx_ledger) if validation is not None else None,
            "status_counts": {
                "no_legal_candidates": sum(
                    self.db.execute(
                        "SELECT count(*) FROM queue WHERE sector_id=? AND principle=?",
                        (sector_id, principle)).fetchone()[0] == 0
                    for sector_id, principle in self.queue_states()),
                "candidate_failures": sum(
                    any(token in str(row[0]) for token in (
                        "REVIEW", "OVER_LIMIT", "SCIENTIFIC_FAILURE"))
                    for row in self.db.execute(
                        "SELECT state FROM queue WHERE state!='pending'")),
                "not_yet_scored": self.db.execute(
                    "SELECT count(*) FROM queue WHERE state='pending'").fetchone()[0],
                "technical_failures": self.db.execute(
                    "SELECT count(*) FROM queue_state WHERE state='technical_failure'"
                ).fetchone()[0],
            },
            "ledger_file": Path(self.db_path).name if self.db_path != ":memory:" else None,
        }
        if validation is not None:
            payload["evaluated_candidates"] = [
                {
                    "candidate_id": result.get("candidate_id"),
                    "mutations": list(signature),
                    **{key: value for key, value in result.items()
                       if key not in {"candidate_id", "mutations", "row", "audit"}},
                }
                for signature, result in sorted(validation.spurs_results.items())
            ]
            payload["foldx_results"] = [
                {
                    "candidate_id": result.get("candidate_id"),
                    "mutations": list(signature),
                    **{key: value for key, value in result.items()
                       if key not in {"candidate_id", "mutations", "audit"}},
                }
                for signature, result in sorted(validation.foldx_results.items())
            ]
        return payload

    def write_ledger(self, path, validation=None):
        from ..util import write_json
        payload = self.ledger(validation)
        write_json(Path(path), payload)
        return payload

    def materialize_region(self, sector_id):
        """Rebuild one region's candidates without retaining the whole task pool."""
        rows = self.db.execute(
            "SELECT candidate_id,signature,metadata FROM candidates "
            "WHERE sector_id=? ORDER BY candidate_id", (str(sector_id),))
        return [self._materialize(row) for row in rows]

    @staticmethod
    def _outcome(state):
        try:
            value = json.loads(state)
        except (TypeError, ValueError):
            return {"outcome": str(state)}
        return value if isinstance(value, dict) else {"outcome": value}

    def _materialize(self, row):
        signature = tuple(tuple(item) for item in json.loads(row["signature"]))
        sector_row = self.db.execute(
            "SELECT sector_id FROM candidates WHERE candidate_id=?",
            (row["candidate_id"],)).fetchone()
        sector = self._sectors.get(sector_row["sector_id"]) if sector_row else None
        if sector is None:
            raise KeyError(f"Missing source region for candidate: {row['candidate_id']}")
        selected = []
        for position, _wt, mutant in signature:
            option = self._options_by_position.get((int(position), str(mutant)))
            if option is None:
                raise KeyError(f"Missing source single evidence: {position}{mutant}")
            selected.append(option)
        candidate = _candidate(sector, selected, len(selected))
        metadata = _decode_metadata(row["metadata"])
        metadata.update(deepcopy(self._geometry_by_sector[sector_row["sector_id"]][
            tuple(metadata["positions"])]))
        candidate.candidate_id = row["candidate_id"]
        candidate.epitope_supported_sites = sum(bool(option.epitope_supported) for option in selected)
        candidate.epitope_score = sum(float(option.epitope_prior or 0.0) for option in selected)
        candidate.coverage_by_scenario = dict(metadata.get("coverage_by_scenario") or {})
        candidate.covered_positions_by_scenario = {
            key: list(value) for key, value in (metadata.get("covered_positions_by_scenario") or {}).items()}
        candidate.geometry_support = dict(metadata.get("geometry_support") or {})
        candidate.risk_components = dict(metadata.get("risk_components") or {})
        candidate.impact_evidence = dict(metadata.get("impact_evidence") or {})
        profile_keys = metadata.get("regional_reference_profile_keys") or ()
        if profile_keys and sector_row["sector_id"] in self._reference_profiles_by_sector:
            shared_profiles = self._reference_profiles_by_sector[sector_row["sector_id"]]
            candidate.impact_evidence["regional_reference_profiles"] = deepcopy(
                {key: shared_profiles[key] for key in profile_keys})
            reference_key = metadata.get("ranking_reference_profile_key")
            if reference_key is not None:
                candidate.impact_evidence["ranking"]["reference_values"] = deepcopy(
                    shared_profiles[reference_key]["values"])
        candidate.impact_evidence["ordinary_principles"] = json.loads(
            self.db.execute(
                "SELECT labels FROM candidates WHERE candidate_id=?",
                (row["candidate_id"],)).fetchone()[0])
        candidate.impact_evidence["ordinary_queue"] = metadata
        return candidate

    def _refresh_state(self, sector_id, principle):
        current = self.db.execute(
            "SELECT state FROM queue_state WHERE sector_id=? AND principle=?",
            (str(sector_id), principle)).fetchone()
        if current and current[0] == "complete":
            return
        pending = self.db.execute("SELECT count(*) FROM queue WHERE sector_id=? AND principle=? AND state='pending'",
                                  (str(sector_id), principle)).fetchone()[0]
        state = "ready" if pending else "exhausted"
        self.db.execute("INSERT INTO queue_state(sector_id,principle,state) VALUES(?,?,?) "
                        "ON CONFLICT(sector_id,principle) DO UPDATE SET state=excluded.state",
                        (str(sector_id), principle, state))

    def add_region(self, sector_id, candidates, rankings, upgrades, *, epitope_applicable=True):
        """Persist independently ordered per-principle queues for one region."""
        identity = {candidate.candidate_id: f"{sector_id}|{candidate.candidate_id}"
                    for candidate in candidates}
        identities = set(identity.values())
        if len(identities) != len(identity):
            raise ValueError("Candidate identities are not unique within the region")
        original_ids = {id(candidate): candidate.candidate_id for candidate in candidates}
        ranked_ids = {name: {item.candidate_id for item in rankings[name]}
                      for name in PRINCIPLES}
        reference_profiles = {}
        for candidate in candidates:
            reference_profiles.update(
                candidate.impact_evidence.get("regional_reference_profiles") or {})
        if reference_profiles:
            self._reference_profiles_by_sector[str(sector_id)] = deepcopy(reference_profiles)
        geometry = self._geometry_by_sector.setdefault(str(sector_id), {})
        for candidate in candidates:
            signature = _mutation_signature(candidate)
            labels = [name for name in PRINCIPLES if candidate.candidate_id in ranked_ids[name]]
            candidate.impact_evidence["ordinary_principles"] = labels
            original_id = candidate.candidate_id
            impact_evidence = dict(candidate.impact_evidence)
            impact_evidence.pop("regional_reference_profiles", None)
            positions = sorted({_position(o) for o in candidate.mutations})
            if tuple(positions) not in geometry:
                geometry[tuple(positions)] = json.loads(json.dumps({
                    "coverage_by_scenario": candidate.coverage_by_scenario,
                    "covered_positions_by_scenario": candidate.covered_positions_by_scenario,
                    "geometry_support": candidate.geometry_support,
                }))
            ranking = dict(impact_evidence.get("ranking") or {})
            profile_key = f"{ranking.get('order')}|{','.join(ranking.get('missing_mask') or [])}"
            reference_key = None
            if (profile_key in reference_profiles
                    and ranking.get("reference_values") is reference_profiles[profile_key]["values"]):
                ranking.pop("reference_values")
                reference_key = profile_key
                impact_evidence["ranking"] = ranking
            metadata = {"order": candidate.order,
                        "positions": positions,
                        "substitutions": [option.mutant_aa for option in candidate.mutations],
                        "risk_components": candidate.risk_components,
                        "impact_evidence": impact_evidence,
                        "regional_reference_profile_keys": list(
                            candidate.impact_evidence.get("regional_reference_profiles") or {}),
                        "ranking_reference_profile_key": reference_key,
                        "triple_upgrade": {name: original_id in upgrades.get(name, set())
                                           for name in PRINCIPLES}}
            self.db.execute("INSERT OR REPLACE INTO candidates VALUES(?,?,?,?,?)",
                            (identity[original_id], str(sector_id), json.dumps(signature),
                             json.dumps(labels), zlib.compress(json.dumps(
                                 metadata, sort_keys=True, separators=(",", ":")).encode("utf-8"), 1)))
            candidate.candidate_id = identity[original_id]
        for principle, ordered in rankings.items():
            principle_upgrades = {f"{sector_id}|{value}" for value in upgrades.get(principle, set())}
            is_upgrade = [c for c in ordered if identity[original_ids[id(c)]] in principle_upgrades]
            doubles = [c for c in ordered if c.order == 2]
            triples = [c for c in ordered if c.order == 3
                       and identity[original_ids[id(c)]] not in principle_upgrades]
            sequence = [*is_upgrade, *doubles, *triples]
            for index, candidate in enumerate(sequence):
                self.db.execute("INSERT OR REPLACE INTO queue(sector_id,principle,ordinal,candidate_id,state) "
                                "VALUES(?,?,?,?, 'pending')",
                                (str(sector_id), principle, index,
                                 identity[original_ids[id(candidate)]]))
            state = (
                "not_applicable"
                if principle == "epitope_support" and not epitope_applicable
                else "ready" if sequence else "exhausted"
            )
            self.db.execute("INSERT OR REPLACE INTO queue_state VALUES(?,?,?,NULL)",
                            (str(sector_id), principle, state))
        self.db.commit()

    def mark_all_budget_limited(self, limit):
        """Persist budget exhaustion for lanes that still need scoring."""
        self.db.execute(
            "UPDATE queue_state SET state='budget_limited',budget_limit=? "
            "WHERE state='ready'", (int(limit),))
        self.db.commit()

    def mark_all_technical_failure(self):
        self.db.execute("UPDATE queue_state SET state='technical_failure' WHERE state='ready'")
        self.db.commit()


def _mutation_signature(candidate):
    return tuple(sorted((_position(o), o.residue.one_letter, o.mutant_aa) for o in candidate.mutations))


def queue_identity(candidate):
    """Return the notation signature used by model workers and validation."""
    return tuple(sorted(option.notation.strip().upper() for option in candidate.mutations))


def _candidate_signature(candidate):
    return tuple(sorted(option.notation for option in candidate.mutations))


def _epitope_key(candidate):
    return (-sum(bool(o.epitope_supported) for o in candidate.mutations) / candidate.order,
            -sum(float(o.epitope_prior or 0.0) for o in candidate.mutations) / candidate.order)


def triple_upgrade_for_principle(triple, best_double, principle):
    """Apply the explicit upgrade rule owned by one independent principle."""
    if triple.order != 3 or best_double.order != 2:
        return False
    if principle == "mechanism_pareto":
        return dominates(candidate_impact(triple), candidate_impact(best_double))
    if principle == "sector_coverage":
        triple_coverage = float(triple.geometry_support.get("target_mean_coverage_fraction") or 0.0)
        double_coverage = float(best_double.geometry_support.get("target_mean_coverage_fraction") or 0.0)
        return triple_coverage >= double_coverage + .05
    if principle == "epitope_support":
        return _epitope_key(triple) < _epitope_key(best_double)
    raise ValueError(f"Unknown principle: {principle}")


def strict_triple_upgrade(triple, double, principle=None):
    """Return a per-principle upgrade decision.

    With no principle, report whether any independent principle upgrades the
    triple. Callers building queues must pass the principle explicitly.
    """
    if principle is not None:
        return triple_upgrade_for_principle(triple, double, principle)
    return any(triple_upgrade_for_principle(triple, double, name) for name in PRINCIPLES)


def _candidate_coverage(candidate, sector, graph, areas, distance_cache, global_distance_cache):
    positions = sorted({_position(option) for option in candidate.mutations})
    for position in positions:
        key = (sector.sector_id, position)
        if key not in distance_cache:
            distance_cache[key] = _distances(graph, position, set(sector.residue_positions))
    distances = [distance_cache[(sector.sector_id, position)] for position in positions]
    candidate.coverage_by_scenario, candidate.covered_positions_by_scenario = _coverage_from_distances(
        distances, set(sector.residue_positions))
    from .spatial_support import support_metrics
    global_distances = {}
    for position in positions:
        if position not in global_distance_cache:
            global_distance_cache[position] = _distances(graph, position, set(graph))
        global_distances[position] = global_distance_cache[position]
    candidate.geometry_support = support_metrics(positions, global_distances, areas,
                                                 sector.residue_positions, sector.sector_id)
    return candidate.geometry_support


def build_candidate_queues(single, sectors, graph, areas, canonical_sequence, db_path=":memory:",
                            maximum_span=24.0, records=None, epitope_applicable=True):
    """Enumerate all legal same-region doubles/triples and build three queues.

    `single` must contain already-qualified, scored single-site evidence. No
    DiscoTope value is read; epitope ordering uses the supplied support flag and
    prior on each mutation option. When epitope evidence does not apply, the
    lane is recorded but excluded from execution.
    """
    store = CandidateQueue(db_path)
    options_by_position = {}
    options_by_key = {}
    records_by_position = (
        dict(records) if isinstance(records, dict)
        else {record.uniprot_pos: record for record in (records or [])}
    )
    for item in single:
        position = getattr(item.residue, "uniprot_pos", None)
        score = getattr(item, "spurs_ddg_kcal_mol", None)
        if (not getattr(item, "eligible", False) or getattr(item, "review_required", False)
                or item.hard_exclusions or score is None or not math.isfinite(float(score))
                or position not in graph):
            continue
        option = _option(item)
        if records_by_position:
            source_record = records_by_position.get(position)
            if source_record is not None:
                option = replace(option, residue=replace(
                    option.residue,
                    one_letter=source_record.aa,
                    auth_seq_id=source_record.auth_seq_id,
                    insertion_code=source_record.insertion_code or "",
                    label_seq_id=getattr(source_record, "label_seq_id", None),
                ))
        option.epitope_prior = float(getattr(item, "epitope_prior", 0.0) or 0.0)
        option.epitope_supported = bool(getattr(item, "epitope_supported", False))
        options_by_position.setdefault(position, []).append(option)
        options_by_key[(position, option.mutant_aa)] = option
    store._options_by_position = options_by_key
    distances, global_distances = {}, {}
    for sector in sorted(sectors, key=lambda value: value.sector_id):
        store._sectors[sector.sector_id] = sector
        positions = sorted(set(sector.residue_positions) & set(options_by_position))
        candidates_by_signature = {}
        coverage_cache = {}
        pair_span_ok = {
            (left, right): _span_ok((left, right), graph, set(sector.residue_positions), maximum_span)[0]
            for left, right in combinations(positions, 2)
        }
        for order in (2, 3):
            for site_set in combinations(positions, order):
                if not all(pair_span_ok.get((left, right), False)
                           for left, right in combinations(site_set, 2)):
                    continue
                for selected in product(*(options_by_position[p] for p in site_set)):
                    if _valid_sequence(canonical_sequence, selected):
                        continue
                    candidate = _candidate(sector, selected, order)
                    candidate.epitope_supported_sites = sum(bool(option.epitope_supported) for option in selected)
                    candidate.epitope_score = sum(float(option.epitope_prior or 0.0) for option in selected)
                    signature = _mutation_signature(candidate)
                    candidates_by_signature.setdefault(signature, candidate)
        region_candidates = list(candidates_by_signature.values())
        for candidate in region_candidates:
            geometry_key = (candidate.order, tuple(_position(option)
                                                   for option in candidate.mutations))
            cached_geometry = coverage_cache.get(geometry_key)
            if cached_geometry is None:
                _candidate_coverage(candidate, sector, graph, areas, distances, global_distances)
                coverage_cache[geometry_key] = (
                    candidate.coverage_by_scenario,
                    candidate.covered_positions_by_scenario,
                    candidate.geometry_support,
                )
            else:
                (candidate.coverage_by_scenario,
                 candidate.covered_positions_by_scenario,
                 candidate.geometry_support) = cached_geometry
        from .spatial_support import annotate_higher_order_regional_gain
        geometry_representatives = {}
        for candidate in region_candidates:
            key = (candidate.order, tuple(_position(option) for option in candidate.mutations))
            previous = geometry_representatives.get(key)
            if previous is None or candidate.candidate_id < previous.candidate_id:
                geometry_representatives[key] = candidate
        annotate_higher_order_regional_gain(
            list(geometry_representatives.values()), [sector], areas)
        for candidate in region_candidates:
            if candidate.order > 2:
                key = (candidate.order, tuple(_position(option) for option in candidate.mutations))
                candidate.risk_components = dict(geometry_representatives[key].risk_components)
        ranked_mechanism = {order: rank_impacts([c for c in region_candidates if c.order == order])
                            for order in (2, 3)}
        mechanism = [*ranked_mechanism[2], *ranked_mechanism[3]]
        epitope = sorted(region_candidates, key=lambda c: (
            -sum(bool(o.epitope_supported) for o in c.mutations) / c.order,
            -sum(float(o.epitope_prior or 0.0) for o in c.mutations) / c.order,
            c.order != 2, c.candidate_id))
        coverage = sorted(region_candidates, key=lambda c: (
            -float(c.geometry_support.get("target_mean_coverage_fraction") or 0.0),
            -float(c.geometry_support.get("mean_marginal_area_a2") or 0.0),
            float(c.geometry_support.get("mean_overlap_fraction") or 0.0),
            -float(c.geometry_support.get("direct_area_a2") or 0.0),
            c.order != 2, c.candidate_id))
        rankings = {
            "epitope_support": epitope,
            "mechanism_pareto": mechanism,
            "sector_coverage": coverage,
        }
        if not epitope_applicable:
            rankings["epitope_support"] = []
        doubles_by_signature = {
            _candidate_signature(candidate): candidate
            for candidate in region_candidates if candidate.order == 2
        }
        upgrades = {principle: set() for principle in PRINCIPLES}
        for principle, ordered in rankings.items():
            best_double = next((c for c in ordered if c.order == 2), None)
            if best_double is None:
                continue
            for triple in (c for c in ordered if c.order == 3):
                parents = [
                    doubles_by_signature.get(tuple(sorted(
                        option.notation for option in candidate.mutations if option is not omitted)))
                    for omitted in triple.mutations
                ]
                comparable = [candidate for candidate in parents if candidate is not None]
                exact_parent = min(
                    comparable,
                    key=lambda candidate: (
                        candidate.geometry_support.get("target_mean_coverage_fraction") or 0.0,
                        candidate.candidate_id,
                    ),
                    default=None,
                )
                if exact_parent is None:
                    continue
                gain = triple.risk_components.get("regional_higher_order_gain")
                if gain is None or not math.isfinite(gain) or gain < .05:
                    continue
                if strict_triple_upgrade(triple, exact_parent, principle):
                    upgrades[principle].add(triple.candidate_id)
        source_records = list(records_by_position.values())
        if not source_records:
            source_records = []
            for position, position_options in options_by_position.items():
                option = position_options[0]
                source_records.append(SimpleNamespace(
                    uniprot_pos=position,
                    aa=option.residue.one_letter,
                    sidechain_exposure=option.impact_evidence.get("sidechain_exposure", {}),
                ))
        complete_candidate_impact(region_candidates, source_records)
        for candidate in region_candidates:
            candidate.impact_evidence["regional_impact"] = candidate_impact(candidate)
        store.add_region(
            sector.sector_id, region_candidates, rankings, upgrades,
            epitope_applicable=epitope_applicable,
        )
    return store
