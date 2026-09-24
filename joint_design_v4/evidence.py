from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from typing import Any

from .. import config
from ..models import ResidueRecord
from .models import (
    EpitopeEvidence,
    ResidueKey,
    SingleMutationEvidence,
    SurfaceEvidence,
    V4MutationOption,
)


def _record(row: dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in row:
            return row[name]
    return None


def _position(row: dict[str, Any]) -> int | None:
    value = _record(row, "uniprot_pos", "position", "residue_index", "index", "pos")
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def surface_evidence(
    records: list[ResidueRecord],
    keys: dict[int, ResidueKey],
    procaliper_rows: list[dict[str, Any]] | None = None,
    areas: dict[int, float] | None = None,
    components: dict[int, int] | None = None,
) -> list[SurfaceEvidence]:
    external = {
        _position(row): row for row in (procaliper_rows or []) if _position(row) is not None
    }
    result = []
    for item in records:
        row = external.get(item.uniprot_pos, {})
        pro_sasa = _record(row, "sasa", "residue_sasa", "sasa_a2")
        pro_rsasa = _record(row, "rsasa", "rSASA", "relative_sasa")
        try:
            pro_sasa = float(pro_sasa) if pro_sasa is not None else None
            pro_rsasa = float(pro_rsasa) if pro_rsasa is not None else None
        except (TypeError, ValueError):
            pro_sasa = pro_rsasa = None
        status = "NOT_COMPARED"
        if pro_rsasa is not None:
            delta = abs(pro_rsasa - item.rsasa)
            status = (
                "AGREE"
                if delta <= 0.05 or delta / max(item.rsasa, 0.01) <= 0.20
                else "SASA_DISAGREEMENT"
            )
        result.append(
            SurfaceEvidence(
                residue=keys[item.uniprot_pos],
                procaliper_sasa_a2=pro_sasa,
                procaliper_rsasa=pro_rsasa,
                native_sasa_a2=float(item.sasa),
                native_rsasa=float(item.rsasa),
                sidechain_sasa_a2=float(item.sidechain_sasa),
                mapped_surface_area_a2=float(
                    (areas or {}).get(item.uniprot_pos, max(item.sasa, item.sidechain_sasa, 1.0))
                ),
                surface_component_id=(components or {}).get(item.uniprot_pos),
                sasa_agreement_status=status,
                provenance={"native_source": "Shrake-Rupley", "procaliper_present": bool(row)},
            )
        )
    return result


def epitope_evidence(
    records: list[ResidueRecord],
    keys: dict[int, ResidueKey],
    discotope_rows: list[dict[str, Any]] | None = None,
    surface_graph: dict[int, dict[int, float]] | None = None,
) -> list[EpitopeEvidence]:
    external = {_position(row): row for row in (discotope_rows or []) if _position(row) is not None}
    raw: dict[int, float | None] = {}
    calibrated: dict[int, float | None] = {}
    score_for_rank: dict[int, float | None] = {}
    for item in records:
        row = external.get(item.uniprot_pos, {})
        value = _record(row, "raw_score", "score", "propensity")
        calibrated_value = _record(row, "calibrated_score")
        try:
            raw[item.uniprot_pos] = float(value) if value is not None else None
        except (TypeError, ValueError):
            raw[item.uniprot_pos] = None
        try:
            calibrated[item.uniprot_pos] = (
                float(calibrated_value) if calibrated_value is not None else None
            )
        except (TypeError, ValueError):
            calibrated[item.uniprot_pos] = None
        score_for_rank[item.uniprot_pos] = (
            calibrated[item.uniprot_pos]
            if calibrated[item.uniprot_pos] is not None
            else raw[item.uniprot_pos]
        )
    surface_positions = set(surface_graph) if surface_graph is not None else set(score_for_rank)
    finite = sorted(
        score_for_rank[position]
        for position in surface_positions
        if score_for_rank.get(position) is not None
    )
    result = []
    for item in records:
        position = item.uniprot_pos
        value = score_for_rank[position]
        percentile = (
            (bisect_left(finite, value) + bisect_right(finite, value) + 1) / (2 * len(finite))
            if value is not None and position in surface_positions and finite
            else None
        )
        row = external.get(item.uniprot_pos, {})
        neighbor_scores = [
            score_for_rank[neighbor]
            for neighbor in (surface_graph or {}).get(position, {})
            if score_for_rank.get(neighbor) is not None
        ]
        local_peak = bool(_record(row, "local_peak")) or bool(
            value is not None and neighbor_scores and value > max(neighbor_scores)
        )
        result.append(
            EpitopeEvidence(
                residue=keys[item.uniprot_pos],
                raw_score=raw[position],
                calibrated_score=calibrated[position],
                surface_percentile=percentile,
                local_peak=local_peak,
                predicted_epitope_label=(bool(row["epitope"]) if "epitope" in row else None),
                status="PREDICTED" if value is not None else "MISSING",
                provenance={
                    "tool": "DiscoTope-3.0" if row else "none",
                    "rank_population": "mapped_surface_graph",
                    "local_peak_source": "worker" if "local_peak" in row else "surface_graph",
                },
            )
        )
    return result


def _chemistry(wt: str, mutant: str) -> str:
    if wt in "KR" and mutant in "ED":
        return "charge_reversal"
    if wt in "WYF" and mutant == "A":
        return "hotspot_removal"
    if wt in config.BULKY_WT and mutant == "A":
        return "volume_reduction"
    if wt in config.SMALL_WT and mutant == "W":
        return "aromatic_insertion"
    if wt in config.CHARGED_WT and mutant not in config.CHARGED_WT:
        return "charge_change"
    return "other"


def _spurs_rows(rows: list[dict[str, Any]] | None) -> dict[tuple[int, str], float]:
    result = {}
    for row in rows or []:
        position = _position(row)
        mutant = str(_record(row, "mutant_aa", "mutant", "aa") or "")
        value = _record(row, "ddg", "spurs_ddg_kcal_mol", "score")
        if position is None or len(mutant) != 1 or value is None:
            scores = row.get("scores")
            if position is not None and isinstance(scores, dict):
                for aa, score in scores.items():
                    try:
                        result[(position, aa)] = float(score)
                    except (TypeError, ValueError):
                        pass
            continue
        try:
            result[(position, mutant)] = float(value)
        except (TypeError, ValueError):
            pass
    return result


def stability_band(value: float | None, params: dict[str, Any]) -> str:
    if value is None:
        return "NOT_SCORED"
    if value <= float(params.get("preferred_max_ddg", 0.0)):
        return "PREFERRED"
    if value <= float(params.get("acceptable_max_ddg", 0.5)):
        return "ACCEPTABLE"
    if value <= float(params.get("review_max_ddg", 1.0)):
        return "REVIEW"
    return "REJECT"


def build_single_evidence(
    records: list[ResidueRecord],
    keys: dict[int, ResidueKey],
    spurs_rows: list[dict[str, Any]] | None,
    params: dict[str, Any],
    *,
    allow_unscored: bool = False,
    surface_by_position: dict[int, SurfaceEvidence] | None = None,
) -> list[SingleMutationEvidence]:
    scores = _spurs_rows(spurs_rows)
    result = []
    for item in records:
        for mutant in "ACDEFGHIKLMNPQRSTVWY":
            if mutant == item.aa:
                continue
            hard = list(item.hard_exclusions)
            if item.disulfide or item.aa == "C":
                hard.append("disulfide-cysteine")
            if item.glycosylation_site or item.glycosylation_sequon:
                hard.append("glycosylation-sequon")
            if item.aa in "GP" or mutant in "GP":
                hard.append("backbone-conformation-sensitive")
            value = scores.get((item.uniprot_pos, mutant))
            band = stability_band(value, params)
            if value is not None and band == "REJECT":
                hard.append("spurs-reject")
            result.append(
                SingleMutationEvidence(
                    residue=keys[item.uniprot_pos],
                    mutant_aa=mutant,
                    spurs_ddg_kcal_mol=value,
                    stability_band=band,
                    chemistry_class=_chemistry(item.aa, mutant),
                    hard_exclusions=sorted(set(hard)),
                    soft_warnings=(
                        ["SASA_DISAGREEMENT"]
                        if (surface_by_position or {}).get(item.uniprot_pos) is not None
                        and surface_by_position[item.uniprot_pos].sasa_agreement_status == "SASA_DISAGREEMENT"
                        else []
                    ),
                    eligible=not hard
                    and (
                        band in {"PREFERRED", "ACCEPTABLE"}
                        or (allow_unscored and band == "NOT_SCORED")
                    ),
                    provenance={"spurs_present": value is not None},
                )
            )
    return result


def compressed_options(
    records: list[ResidueRecord],
    single: list[SingleMutationEvidence],
    max_options: int = 3,
    epitope_rows: list[EpitopeEvidence] | None = None,
    epitope_priority_weight: float = 1.0,
    epitope_supported_percentile: float = 0.80,
) -> list[V4MutationOption]:
    by_pos: dict[int, list[SingleMutationEvidence]] = defaultdict(list)
    for item in single:
        if item.eligible:
            by_pos[item.residue.uniprot_pos or -1].append(item)
    epitope_by_position = {
        item.residue.uniprot_pos: item
        for item in (epitope_rows or [])
        if item.residue.uniprot_pos is not None
    }
    result = []
    for record in records:
        options = by_pos.get(record.uniprot_pos, [])
        options.sort(
            key=lambda item: (
                item.stability_band not in {"PREFERRED", "ACCEPTABLE"},
                item.spurs_ddg_kcal_mol if item.spurs_ddg_kcal_mol is not None else 0.0,
                item.mutant_aa,
            )
        )
        selected: list[SingleMutationEvidence] = []
        classes: set[str] = set()
        for item in options:
            if item.chemistry_class not in classes or len(selected) < max_options:
                selected.append(item)
                classes.add(item.chemistry_class)
            if len(selected) >= max_options:
                break
        for item in selected:
            accessibility = min(1.0, max(0.0, record.rsasa))
            epitope = epitope_by_position.get(record.uniprot_pos)
            epitope_prior = (
                min(1.0, max(0.0, epitope.surface_percentile))
                if epitope is not None and epitope.surface_percentile is not None
                else 0.0
            )
            epitope_supported = bool(
                epitope is not None
                and (
                    epitope.local_peak
                    or (
                        epitope.surface_percentile is not None
                        and epitope.surface_percentile >= epitope_supported_percentile
                    )
                )
            )
            perturbation = (
                0.6 if item.chemistry_class in {"charge_reversal", "hotspot_removal"} else 0.3
            )
            risk = item.spurs_ddg_kcal_mol if item.spurs_ddg_kcal_mol is not None else 0.0
            result.append(
                V4MutationOption(
                    residue=item.residue,
                    mutant_aa=item.mutant_aa,
                    notation=f"{record.aa}{record.uniprot_pos}{item.mutant_aa}",
                    chemistry_class=item.chemistry_class,
                    spurs_ddg_kcal_mol=item.spurs_ddg_kcal_mol,
                    score=round(
                        accessibility
                        + perturbation
                        + epitope_priority_weight * epitope_prior
                        - max(0.0, risk) * 0.25,
                        6,
                    ),
                    eligible=True,
                    hard_exclusions=list(item.hard_exclusions),
                    warnings=list(item.soft_warnings),
                    epitope_prior=round(epitope_prior, 6),
                    epitope_supported=epitope_supported,
                )
            )
    return result


def attach_epitope_islands(
    evidence: list[EpitopeEvidence],
    graph: dict[int, dict[int, float]],
    seed_percentile: float = 0.80,
    local_peak_percentile: float = 0.65,
    max_distance: float = 12.0,
) -> None:
    seeds = {
        item.residue.uniprot_pos
        for item in evidence
        if item.residue.uniprot_pos is not None
        and (
            (item.surface_percentile is not None and item.surface_percentile >= seed_percentile)
            or (
                item.local_peak
                and item.surface_percentile is not None
                and item.surface_percentile >= local_peak_percentile
            )
        )
    }
    islands: list[set[int]] = []
    for seed in sorted(seeds):
        # DiscoTope scores the complete normalized chain, while the surface
        # graph intentionally contains only mapped surface residues. A high
        # score on a buried/non-graph residue must not create an empty island
        # that later fails on min(empty).
        if seed not in graph:
            continue
        group = {position for position in graph if _shortest(graph, seed, position) <= max_distance}
        if not group:
            continue
        merged = [island for island in islands if island & group]
        if merged:
            group |= set().union(*merged)
            islands = [island for island in islands if island not in merged]
        islands.append(group)
    for index, island in enumerate(sorted(islands, key=lambda item: min(item)), start=1):
        for item in evidence:
            if item.residue.uniprot_pos in island:
                item.island_id = f"EI{index:03d}"


def _shortest(graph: dict[int, dict[int, float]], start: int, end: int) -> float:
    if start == end:
        return 0.0
    distances = {start: 0.0}
    queue = [(0.0, start)]
    import heapq

    while queue:
        distance, current = heapq.heappop(queue)
        if current == end:
            return distance
        if distance != distances.get(current):
            continue
        for neighbor, step in graph.get(current, {}).items():
            proposed = distance + step
            if proposed < distances.get(neighbor, math.inf):
                distances[neighbor] = proposed
                heapq.heappush(queue, (proposed, neighbor))
    return math.inf
