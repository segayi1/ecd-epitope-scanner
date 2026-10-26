from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from typing import Any

from .. import config
from ..models import ResidueRecord
from ..stability import AUTO_EXCLUDED_MUTANTS
from .coverage_profile import descriptors
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


def _sasa_values_close(left: float, right: float) -> bool:
    delta = abs(left - right)
    return delta <= 0.05 or delta / max(left, 0.01) <= 0.20


def _procaliper_method_status(
    procaliper_rows: list[dict[str, Any]],
    procaliper_provenance: dict[str, Any] | None,
    normalized_context_method: dict[str, Any] | None,
) -> str:
    if not procaliper_rows:
        return "NOT_COMPARED"
    if normalized_context_method is None:
        return "INPUT_SCOPE_UNVERIFIED"
    provenance = procaliper_provenance or {}
    required = {
        "algorithm",
        "biopython_version",
        "probe_radius_a",
        "sphere_points",
        "rsasa_reference",
        "atom_radii_a",
    }
    if not required.issubset(provenance):
        return "PROVENANCE_INCOMPLETE"
    try:
        settings_match = (
            "ShrakeRupley" in str(provenance["algorithm"])
            and float(provenance["probe_radius_a"])
            == float(normalized_context_method["probe_radius_a"])
            and int(provenance["sphere_points"])
            == int(normalized_context_method["sphere_points"])
            and all(
                float(provenance["atom_radii_a"][element])
                == float(normalized_context_method["atom_radii_a"][element])
                for element in ("C", "N", "O", "S")
            )
            and "Tien" in str(provenance["rsasa_reference"])
            and "Tien" in str(normalized_context_method["radii_reference"])
        )
    except (KeyError, TypeError, ValueError):
        return "PROVENANCE_INCOMPLETE"
    if not settings_match:
        return "SETTINGS_DIFFER"
    return (
        "MATCHED"
        if str(provenance["biopython_version"])
        == str(normalized_context_method["biopython_version"])
        else "VERSION_DIFFERS"
    )


def surface_evidence(
    records: list[ResidueRecord],
    keys: dict[int, ResidueKey],
    procaliper_rows: list[dict[str, Any]] | None = None,
    areas: dict[int, float] | None = None,
    components: dict[int, int] | None = None,
    normalized_context_rows: dict[int, dict[str, float]] | None = None,
    normalized_context_method: dict[str, Any] | None = None,
    procaliper_provenance: dict[str, Any] | None = None,
    rsasa_threshold: float = config.RSASA_THRESHOLD,
) -> list[SurfaceEvidence]:
    external = {_position(row): row for row in (procaliper_rows or []) if _position(row) is not None}
    comparison_method_status = _procaliper_method_status(
        procaliper_rows or [], procaliper_provenance, normalized_context_method
    )
    result = []
    for item in records:
        row = external.get(item.uniprot_pos, {})
        normalized_row = (normalized_context_rows or {}).get(item.uniprot_pos, {})
        pro_sasa = _record(row, "sasa", "residue_sasa", "sasa_a2")
        pro_rsasa = _record(row, "rsasa", "rSASA", "relative_sasa")
        normalized_sasa = normalized_row.get("sasa")
        normalized_rsasa = normalized_row.get("rsasa")
        try:
            pro_sasa = float(pro_sasa) if pro_sasa is not None else None
            pro_rsasa = float(pro_rsasa) if pro_rsasa is not None else None
            normalized_sasa = float(normalized_sasa) if normalized_sasa is not None else None
            normalized_rsasa = float(normalized_rsasa) if normalized_rsasa is not None else None
        except (TypeError, ValueError):
            pro_sasa = pro_rsasa = normalized_sasa = normalized_rsasa = None
        status = "NOT_COMPARED"
        if pro_rsasa is not None:
            if normalized_rsasa is None:
                status = "COMPARISON_CONTEXT_UNVERIFIED"
            elif comparison_method_status not in {"MATCHED", "VERSION_DIFFERS"}:
                status = "METHOD_SETTINGS_UNVERIFIED"
            else:
                same_class = (normalized_rsasa >= rsasa_threshold) == (
                    pro_rsasa >= rsasa_threshold
                )
                status = (
                    "AGREE"
                    if same_class and _sasa_values_close(normalized_rsasa, pro_rsasa)
                    else "SASA_DISAGREEMENT"
                )
        native_context_delta = (
            abs(float(item.rsasa) - normalized_rsasa)
            if normalized_rsasa is not None
            else None
        )
        native_context_classification_changed = (
            (float(item.rsasa) >= rsasa_threshold) != (normalized_rsasa >= rsasa_threshold)
            if normalized_rsasa is not None
            else None
        )
        native_context_status = "NOT_COMPARED"
        if native_context_delta is not None and normalized_rsasa is not None:
            native_context_status = (
                "CONTEXT_STABLE"
                if not native_context_classification_changed
                and _sasa_values_close(float(item.rsasa), normalized_rsasa)
                else "CONTEXT_SENSITIVE"
            )
        procaliper_classification_changed = (
            (pro_rsasa >= rsasa_threshold) != (normalized_rsasa >= rsasa_threshold)
            if pro_rsasa is not None and normalized_rsasa is not None
            else None
        )
        result.append(
            SurfaceEvidence(
                residue=keys[item.uniprot_pos],
                procaliper_sasa_a2=pro_sasa,
                procaliper_rsasa=pro_rsasa,
                native_sasa_a2=float(item.sasa),
                native_rsasa=float(item.rsasa),
                normalized_context_sasa_a2=normalized_sasa,
                normalized_context_rsasa=normalized_rsasa,
                native_context_absolute_delta=native_context_delta,
                native_context_status=native_context_status,
                native_context_classification_changed=native_context_classification_changed,
                procaliper_classification_changed=procaliper_classification_changed,
                comparison_method_status=comparison_method_status,
                sidechain_sasa_a2=float(item.sidechain_sasa),
                mapped_surface_area_a2=float(
                    (areas or {}).get(item.uniprot_pos, max(item.sasa, item.sidechain_sasa, 1.0))
                ),
                surface_component_id=(components or {}).get(item.uniprot_pos),
                sasa_agreement_status=status,
                provenance={
                    "native_scope": "selected source protein chain; standard residues retained, partners excluded",
                    "normalized_context_scope": (
                        normalized_context_method.get("scope")
                        if normalized_context_method
                        else None
                    ),
                    "procaliper_present": bool(row),
                    "procaliper_input_structure_sha256": (
                        (procaliper_provenance or {}).get("input_structure_sha256")
                    ),
                    "comparison_method_status": comparison_method_status,
                },
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
    if (wt in "KR" and mutant in "ED") or (wt in "ED" and mutant in "KR"):
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
    if value is None or not math.isfinite(value):
        return "NOT_SCORED"
    if value <= float(params.get("preferred_max_ddg", 0.5)):
        return "PREFERRED"
    if value <= float(params.get("acceptable_max_ddg", 1.0)):
        return "ACCEPTABLE"
    if value <= float(params.get("review_max_ddg", 2.0)):
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
    sequence: str | None = None,
    coverage_profile: bool = False,
) -> list[SingleMutationEvidence]:
    scores = _spurs_rows(spurs_rows)
    result = []
    for item in records:
        surface_item = (surface_by_position or {}).get(item.uniprot_pos)
        surface_warnings = []
        if .20 <= item.rsasa < .25:
            surface_warnings.append("RSASA_BOUNDARY_020_025")
        if surface_item is not None:
            if surface_item.sasa_agreement_status == "SASA_DISAGREEMENT":
                surface_warnings.append("SASA_DISAGREEMENT")
            elif surface_item.sasa_agreement_status == "METHOD_SETTINGS_UNVERIFIED":
                surface_warnings.append("SASA_METHOD_SETTINGS_UNVERIFIED")
            if surface_item.native_context_status == "CONTEXT_SENSITIVE":
                surface_warnings.append("SASA_CONTEXT_SENSITIVE")
        for mutant in "ACDEFGHIKLMNPQRSTVWY":
            if mutant == item.aa:
                continue
            hard = list(item.hard_exclusions)
            if coverage_profile and surface_item is not None and surface_item.native_context_classification_changed is True:
                hard.append('native-normalized-sasa-threshold-conflict')
            if not item.coordinate_eligible:
                hard.append("geometry-ineligible")
            if item.disulfide or item.aa == "C":
                hard.append("disulfide-cysteine")
            if item.glycosylation_site or item.glycosylation_sequon:
                hard.append("glycosylation-sequon")
            if mutant in AUTO_EXCLUDED_MUTANTS:
                hard.append("auto-excluded-mutant")
            if item.aa in "GP" or mutant in "GP":
                hard.append("backbone-conformation-sensitive")
            value = scores.get((item.uniprot_pos, mutant))
            if value is not None and not math.isfinite(value):
                hard.append("nonfinite-spurs-score")
                value = None
            if coverage_profile:
                from .coverage_profile import sequence_reasons
                if sequence is None:
                    hard.append("sequence-unavailable")
                else:
                    hard.extend(sequence_reasons(sequence, [(item.uniprot_pos, item.aa, mutant)]))
                if value is None:
                    hard.append("missing-spurs-score")
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
                    soft_warnings=sorted(set(surface_warnings + item.review_reasons)),
                    review_required=bool(item.review_reasons),
                    eligible=not hard
                    and (
                        band in {"PREFERRED", "ACCEPTABLE"}
                        or (allow_unscored and band == "NOT_SCORED")
                    ),
                    provenance={"spurs_present": value is not None},
                    chemical_descriptors=descriptors(item.aa, mutant),
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
    coverage_profile: bool = False,
) -> list[V4MutationOption]:
    # Retained for call compatibility; DiscoTope is no longer folded into a
    # mutation quality score. It stays as separate hypothesis-prior evidence.
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
        if coverage_profile:
            from .coverage_profile import select_menu
            selected = select_menu(options, max_options)
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
            result.append(
                V4MutationOption(
                    residue=item.residue,
                    mutant_aa=item.mutant_aa,
                    notation=f"{record.aa}{record.uniprot_pos}{item.mutant_aa}",
                    chemistry_class=item.chemistry_class,
                    spurs_ddg_kcal_mol=item.spurs_ddg_kcal_mol,
                    # Compatibility field only: bounded exposure × chemistry
                    # proxy. Candidate ranking consumes the component fields.
                    score=round(accessibility * perturbation, 6),
                    eligible=True,
                    review_required=item.review_required,
                    hard_exclusions=list(item.hard_exclusions),
                    warnings=list(item.soft_warnings),
                    accessibility_score=round(accessibility, 6),
                    chemical_disruption_score=perturbation,
                    stability_band=item.stability_band,
                    glycan_risk=(["GLYCAN_SHIELDED"] if record.glycan_shielded else []),
                    structural_risk_flags=list(item.soft_warnings),
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
