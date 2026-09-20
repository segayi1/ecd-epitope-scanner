"""Combinatorial mutation design for epitope-mapping probes.

Design goal (docs/改进方案.md): each Stage-1 probe must be strong enough that
a binding assay can produce an interpretable *negative* for its patch, while
staying fold-safe.  Substitutions are directed by the WT residue class
(Bogan & Thorn 1998 hotspot chemistry), every probe carries an anchor
mutation, and combinations are found by exhaustive enumeration over the
FoldX/proxy-cleared candidate lists (<= 81 combos) instead of a per-site
greedy pick.  Stage-2 deconvolution presets (singles + confirmatory pairs)
are generated up front so a Stage-1 hit converts to a residue-level answer
without another design cycle.
"""
from __future__ import annotations

import itertools
import os
from typing import Any

import numpy as np

from . import config
from .models import (
    DesignRecord,
    MutationChoice,
    ProteinContext,
    RegionRecord,
    ResidueRecord,
    to_dict,
)
from .sequences import apply_mutations, coding_fragment
from .stability import (
    AMINO_ACIDS,
    EngineInventory,
    assign_percentiles,
    evaluate_substitution,
    is_anchor,
    repair_with_foldx,
    score_many_with_foldx,
    score_with_evoef2,
)

# Kept for backwards compatibility with external tooling importing this name.
FOLDX_SINGLE_GATE = config.FOLDX_GATE_DEFAULT
FOLDX_COMBINATION_GATE = config.FOLDX_COMBINATION_GATE


def _distance(left: ResidueRecord, right: ResidueRecord) -> float:
    if left.ca is None or right.ca is None:
        return float("inf")
    return float(np.linalg.norm(np.asarray(left.ca) - np.asarray(right.ca)))


def _site_priority(record: ResidueRecord) -> float:
    score = record.rsasa * 2.0 + min(1.2, record.sidechain_sasa / 100.0)
    score += 0.6 / (1 + record.contacts / 6.0)
    # Interface exposure can look attractive geometrically, but a binding loss
    # caused by disrupted native assembly is hard to interpret in a
    # cell-surface assay.  Keep these sites available, but prefer a non-
    # interface alternative whenever one exists.
    score -= 0.30 if record.interface else 0.0
    score -= 0.35 if record.salt_bridge else 0.0
    score -= 0.25 if record.domain_interface else 0.0
    score -= 0.35 if record.glycan_shielded else 0.0
    score -= 0.45 if record.aa in "GP" else 0.0
    score -= 0.60 if record.engineered_neighborhood else 0.0
    if record.conservation_score is not None:
        # Variable (low-conservation) sites are likelier antigenic and cheaper
        # to mutate functionally; fully conserved sites are downweighted.
        score += 0.5 * (1.0 - record.conservation_score)
    return score


def select_key_sites(items: list[ResidueRecord], limit: int = 8) -> list[ResidueRecord]:
    ranked = sorted(items, key=lambda item: (-_site_priority(item), item.uniprot_pos))
    selected: list[ResidueRecord] = []
    for record in ranked:
        if all(_distance(record, existing) >= config.PATCH_MIN_SITE_SPACING_A for existing in selected):
            selected.append(record)
        if len(selected) >= limit:
            break
    if len(selected) < min(limit, len(ranked)):
        for record in ranked:
            if record not in selected:
                selected.append(record)
            if len(selected) >= limit:
                break
    return selected


def saturation_scan(
    record: ResidueRecord, protein: ProteinContext
) -> list[MutationChoice]:
    choices = [
        evaluate_substitution(record, mutant, protein)
        for mutant in AMINO_ACIDS
        if mutant != record.aa
    ]
    assign_percentiles(choices)
    return choices


def _ensure_anchor_site(
    items: list[ResidueRecord],
    key_sites: list[ResidueRecord],
    scans: dict[int, list[MutationChoice]],
    protein: ProteinContext,
) -> tuple[list[ResidueRecord], dict[int, list[MutationChoice]]]:
    """Guarantee the anchor policy has material to work with.

    If no key site can carry an anchor substitution (no K/R/E/D/W/Y/F among
    them), admit the best anchor-capable neighbour from the remaining patch
    members -- spacing permitting.  Without this, membrane-proximal patches
    dominated by V/L/H residues would systematically lose anchor protection
    and be downgraded to weaker probes.
    """
    from .stability import is_anchor as _is_anchor

    def has_anchor(record: ResidueRecord) -> bool:
        scan = scans.get(record.uniprot_pos)
        return bool(scan) and any(_is_anchor(record, choice.mutant) for choice in scan if choice.allowed)

    if any(has_anchor(record) for record in key_sites):
        return key_sites, scans
    selected = {record.uniprot_pos for record in key_sites}
    candidates = [record for record in items if record.uniprot_pos not in selected]
    candidates.sort(key=lambda record: (-_site_priority(record), record.uniprot_pos))
    for record in candidates:
        if any(
            _distance(record, chosen) < config.PATCH_MIN_SITE_SPACING_A
            for chosen in key_sites
        ):
            continue
        scan = saturation_scan(record, protein)
        if not any(_is_anchor(record, choice.mutant) for choice in scan if choice.allowed):
            continue
        return key_sites + [record], {**scans, record.uniprot_pos: scan}
    return key_sites, scans


def _probe_utility(choice: MutationChoice) -> float:
    """Mapping utility: disruption power minus a stability-risk tariff."""
    return choice.disruption_score - config.DISRUPTION_LAMBDA_RISK * choice.stability_score


def _best_choice(choices: list[MutationChoice]) -> MutationChoice | None:
    allowed = [choice for choice in choices if choice.allowed and choice.stability_percentile <= 0.9]
    if not allowed:
        return None
    return max(allowed, key=lambda choice: (_probe_utility(choice), -choice.stability_score, choice.mutant))


def _foldx_gate_for(choice: MutationChoice, by_position: dict[int, ResidueRecord]) -> float:
    record = by_position[choice.position]
    return config.foldx_single_gate(record.rsasa)


def _foldx_filter_single_choices(
    key_sites: list[ResidueRecord],
    scans: dict[int, list[MutationChoice]],
    records_by_position: dict[int, ResidueRecord],
    binary: str,
    repaired_pdb: Any,
    workdir: Any,
) -> bool:
    """Gate the per-site candidate menus with measured FoldX single ΔΔG.

    Only the top `MAX_CANDIDATES_PER_SITE` probe-ranked substitutions per key
    site enter the screen; unscored substitutions are explicitly disabled so a
    proxy-only choice can never reach the final combination.
    """
    jobs: list[list[tuple[str, int, str, str]]] = []
    job_choices: list[MutationChoice] = []
    for record in key_sites:
        options = [choice for choice in scans[record.uniprot_pos] if choice.allowed]
        options.sort(key=lambda choice: (-_probe_utility(choice), choice.mutant))
        scored_positions = {choice.mutant for choice in options[: config.MAX_CANDIDATES_PER_SITE]}
        for choice in scans[record.uniprot_pos]:
            if choice.mutant not in scored_positions:
                choice.allowed = False
                choice.stability_class = "not_scored"
                if "not evaluated by FoldX single-point screen" not in choice.reasons:
                    choice.reasons.append("not evaluated by FoldX single-point screen")
        for choice in options[: config.MAX_CANDIDATES_PER_SITE]:
            jobs.append(
                [
                    (
                        choice.wt,
                        records_by_position[choice.position].auth_seq_id,
                        records_by_position[choice.position].chain_id,
                        choice.mutant,
                    )
                ]
            )
            job_choices.append(choice)
    if not jobs:
        return False
    scores = score_many_with_foldx(
        binary,
        repaired_pdb,
        jobs,
        workdir,
        runs=config.FOLDX_SINGLE_SCREEN_RUNS,
        already_repaired=True,
    )
    valid = 0
    for choice, score in zip(job_choices, scores):
        if score is None:
            choice.allowed = False
            choice.stability_class = "foldx_failed"
            choice.reasons.append("FoldX single-point score unavailable")
            continue
        valid += 1
        gate = _foldx_gate_for(choice, records_by_position)
        choice.foldx_ddg = score
        choice.stability_score = score
        if score < -gate:
            choice.stability_class = "stabilizing"
        elif score <= gate:
            choice.stability_class = "neutral"
        else:
            choice.stability_class = "destabilizing"
        choice.allowed = choice.allowed and score <= gate
        choice.reasons.append(
            f"FoldX single-point ΔΔG={score:.3f} kcal/mol (gate +{gate:.2f} by exposure)"
        )
        if score > gate:
            choice.reasons.append(f"exceeds the +{gate:.2f} kcal/mol exposure-linked single-point gate")
    for choices in scans.values():
        assign_percentiles(choices)
    return valid > 0


def _combination_candidates(
    key_sites: list[ResidueRecord],
    scans: dict[int, list[MutationChoice]],
) -> dict[int, list[MutationChoice]]:
    """Per-site candidate menu: FoldX/proxy-cleared, probe-ranked.

    One slot is reserved for the best available anchor when the utility top-N
    would otherwise crowd it out: the anchor policy only means something if an
    anchor actually survives into the enumeration pool.
    """
    from .stability import is_anchor as _is_anchor

    menu: dict[int, list[MutationChoice]] = {}
    for record in key_sites:
        options = [choice for choice in scans[record.uniprot_pos] if choice.allowed]
        options.sort(key=lambda choice: (-_probe_utility(choice), choice.mutant))
        top = options[: config.MAX_CANDIDATES_PER_SITE]
        if not any(_is_anchor(record, choice.mutant) for choice in top):
            anchor = next(
                (choice for choice in options[config.MAX_CANDIDATES_PER_SITE :] if _is_anchor(record, choice.mutant)),
                None,
            )
            if anchor is not None:
                top = top[: config.MAX_CANDIDATES_PER_SITE - 1] + [anchor]
        menu[record.uniprot_pos] = top
    return menu


def _combo_is_valid(
    mutations: list[MutationChoice],
    by_position: dict[int, ResidueRecord],
    measured: bool,
    require_anchor: bool,
    max_combo_spread: float | None = None,
) -> tuple[bool, str]:
    max_combo_spread = max_combo_spread or config.MAX_COMBO_SPREAD_A
    records = [by_position[item.position] for item in mutations]
    for index, left in enumerate(records):
        for right in records[index + 1 :]:
            distance = _distance(left, right)
            if distance < config.MIN_COMBO_SPACING_A:
                return False, "spacing"
            # Probe must fit inside one antibody footprint: mutations spread
            # across opposite faces cannot perturb one paratope together.
            if distance > max_combo_spread:
                return False, "spread"
    if require_anchor and not any(is_anchor(record, choice.mutant) for record, choice in zip(records, mutations)):
        return False, "anchor"
    total_risk = sum(choice.stability_score for choice in mutations)
    if measured:
        if total_risk > config.FOLDX_COMBINATION_GATE:
            return False, "combination gate"
    else:
        if total_risk > max(2.0, len(mutations) * 1.05):
            return False, "proxy cap"
    return True, ""


def _pareto_front(
    scored: list[tuple[float, float, list[MutationChoice]]],
) -> list[tuple[float, float, list[MutationChoice]]]:
    """Skyline over (minimise stability risk, maximise disruption power).

    Returns the non-dominated combinations ordered from the conservative end
    (lowest risk) to the aggressive end (highest disruption).  Unlike a linear
    utility `D - lambda*risk`, the front keeps every trade-off visible and
    defers the weight choice to the experiment (which probe's folding control
    actually passed).
    """
    ordered = sorted(scored, key=lambda item: (item[0], -item[1]))
    front: list[tuple[float, float, list[MutationChoice]]] = []
    best_disruption = float("-inf")
    for risk, disruption, mutations in ordered:
        if disruption > best_disruption:
            front.append((risk, disruption, mutations))
            best_disruption = disruption
    return front


def _pareto_knee(
    front: list[tuple[float, float, list[MutationChoice]]],
) -> list[MutationChoice]:
    """Primary probe = knee of the front: max disruption per unit risk.

    The +1.0 kcal/mol anchor is the documented FoldX single-mutation noise
    scale (docs/改进方案.md): it stops tiny-risk single-ish points from being
    ranked above genuinely strong probes, without restoring a tunable lambda.
    """
    return max(
        front,
        key=lambda item: (item[1] / (item[0] + config.PARETO_KNEE_RISK_ANCHOR), -item[0]),
    )[2]


def _enumerate_combinations(
    key_sites: list[ResidueRecord],
    scans: dict[int, list[MutationChoice]],
    by_position: dict[int, ResidueRecord],
    measured: bool,
    max_combo_spread: float | None = None,
    minimum_sites: int | None = None,
    maximum_sites: int | None = None,
) -> tuple[list[tuple[float, float, list[MutationChoice]]], list[str]]:
    """Exhaustive search over the cleared per-site menus.

    Returns the Pareto front as (risk, disruption, mutations) tuples ordered
    conservative -> aggressive, plus any policy downgrades the search had to
    take (anchor requirement relaxed, 2-site fallback).
    """
    menu = _combination_candidates(key_sites, scans)
    active_sites = [record for record in key_sites if menu.get(record.uniprot_pos)]
    downgrades: list[str] = []
    if len(active_sites) < 2:
        return [], downgrades

    policy_minimum = config.MIN_COMBO_SITES if minimum_sites is None else max(2, minimum_sites)
    policy_maximum = config.MAX_COMBO_SITES if maximum_sites is None else max(2, maximum_sites)
    policy_minimum = min(policy_minimum, policy_maximum)

    def search(require_anchor: bool, min_sites: int) -> list[tuple[float, float, list[MutationChoice]]]:
        scored: list[tuple[float, float, list[MutationChoice]]] = []
        considered = 0
        maximum = min(policy_maximum, len(active_sites))
        for size in range(maximum, min_sites - 1, -1):
            for site_subset in itertools.combinations(active_sites, size):
                pools = [menu[record.uniprot_pos] for record in site_subset]
                for picks in itertools.product(*pools):
                    considered += 1
                    if considered > config.MAX_COMBINATIONS_PER_REGION:
                        return _pareto_front(scored)
                    mutations = list(picks)
                    valid, _ = _combo_is_valid(
                        mutations, by_position, measured, require_anchor, max_combo_spread
                    )
                    if not valid:
                        continue
                    scored.append(
                        (
                            round(sum(choice.stability_score for choice in mutations), 4),
                            round(sum(choice.disruption_score for choice in mutations), 4),
                            mutations,
                        )
                    )
        return _pareto_front(scored)

    front = search(require_anchor=True, min_sites=policy_minimum)
    if front:
        return front, downgrades
    front = search(require_anchor=False, min_sites=policy_minimum)
    if front:
        downgrades.append(
            "No combination with an anchor mutation passed the stability gates; the anchor requirement was relaxed."
        )
        return front, downgrades
    # 2-site fallback: still prefer an anchor-carrying pair before giving up
    # on the anchor policy entirely.
    front = search(require_anchor=True, min_sites=2)
    if front:
        downgrades.append(
            "No >=3-site combination passed; a 2-site probe with an anchor mutation was accepted after full enumeration."
        )
        return front, downgrades
    front = search(require_anchor=False, min_sites=2)
    if front:
        downgrades.append(
            "No >=3-site combination passed; a 2-site probe was accepted after full enumeration."
        )
        return front, downgrades
    return [], downgrades


def _combo_score(mutations: list[MutationChoice]) -> float:
    return round(sum(_probe_utility(choice) for choice in mutations), 4)


class _CodingCache:
    """Per-job cache so each distinct mutant sequence is optimized once."""

    def __init__(self, protein: ProteinContext, original_cds: str | None):
        self.protein = protein
        self.original_cds = original_cds
        self.cache: dict[str, tuple[str, str]] = {}

    def fragment(self, mutant_sequence: str) -> tuple[str, str]:
        if mutant_sequence not in self.cache:
            self.cache[mutant_sequence] = coding_fragment(
                self.protein, mutant_sequence, self.original_cds
            )
        return self.cache[mutant_sequence]


def _make_design(
    design_id: str,
    name: str,
    region_id: int | None,
    mutations: list[MutationChoice],
    status: str,
    inventory: EngineInventory,
    protein: ProteinContext,
    coding: _CodingCache,
    rationale: str,
    warnings: list[str],
    role: str,
    parent_id: str | None = None,
    saturation: dict[int, list[MutationChoice]] | None = None,
) -> DesignRecord:
    mutant_sequence = apply_mutations(protein.ecd_sequence, mutations, protein.ecd_start)
    mutant_dna, dna_source = coding.fragment(mutant_sequence)
    if f"Coding fragment source: {dna_source}." not in warnings:
        warnings = warnings + [f"Coding fragment source: {dna_source}."]
    design = DesignRecord(
        id=design_id,
        name=name,
        region_id=region_id,
        mutations=mutations,
        status=status,
        stability_engine=inventory.primary_name,
        combined_stability_score=round(sum(item.stability_score for item in mutations), 4),
        perturbation_score=round(sum(item.perturbation_score for item in mutations), 4),
        rationale=rationale,
        warnings=warnings,
        mutant_ecd_sequence=mutant_sequence,
        mutant_coding_sequence=mutant_dna,
        score_origin="proxy",
        design_role=role,
        parent_design_id=parent_id,
        disruption_score=round(sum(item.disruption_score for item in mutations), 4),
        anchor_notation=next(
            (choice.notation for choice in mutations if choice.menu_class in {"charge_reversal", "hotspot_to_ala", "bulky_to_ala"}),
            None,
        ),
    )
    if saturation:
        design.saturation_scan = {
            str(position): [to_dict(choice) | {"notation": choice.notation} for choice in values]
            for position, values in saturation.items()
        }
    return design


def _coverage_matrix(
    designs: list[DesignRecord],
    regions: list[RegionRecord],
    by_position: dict[int, ResidueRecord],
) -> list[dict[str, Any]]:
    """Fraction of each patch within perturbation reach of each probe."""
    rows: list[dict[str, Any]] = []
    for design in designs:
        if not design.mutations:
            continue
        mutation_ca = [
            np.asarray(by_position[choice.position].ca)
            for choice in design.mutations
            if choice.position in by_position and by_position[choice.position].ca is not None
        ]
        for region in regions:
            reach = 0
            considered = 0
            for position in (region.envelope_residues or region.residues):
                record = by_position.get(position)
                if record is None or record.ca is None:
                    continue
                considered += 1
                point = np.asarray(record.ca)
                if any(float(np.linalg.norm(point - origin)) <= config.STAGE2_PERTURBATION_RADIUS_A for origin in mutation_ca):
                    reach += 1
            rows.append(
                {
                    "design_id": design.id,
                    "patch_id": region.id,
                    "perturbed_fraction": round(reach / considered, 4) if considered else 0.0,
                }
            )
    return rows


def _apply_plm_evidence(
    scans: dict[int, list[MutationChoice]],
    plm_llrs: dict[tuple[int, str], float] | None,
) -> None:
    """Attach PLM zero-shot evidence and apply the consensus veto.

    A strongly deleterious PLM LLR (< -3.0 nats) vetoes the substitution even
    when the structural proxy passed: the two evidence lines are independent,
    which is the point of a consensus gate.
    """
    if not plm_llrs:
        return
    for choices in scans.values():
        for choice in choices:
            llr = plm_llrs.get((choice.position, choice.mutant))
            if llr is None:
                continue
            choice.esm_llr = llr
            choice.reasons.append(f"ESM LLR={llr:.2f}")
            if llr < -3.0:
                choice.allowed = False
                choice.reasons.append("consensus veto: PLM strongly deleterious")


def _annotate_design_risks(
    designs: list[DesignRecord], by_position: dict[int, ResidueRecord]
) -> None:
    """Preserve interpretation risks next to every exported construct."""
    for design in designs:
        items = [by_position[item.position] for item in design.mutations if item.position in by_position]
        if any(item.interface for item in items):
            design.warnings.append(
                "Native interface mutation present: interpret binding loss only with surface-expression and assembly/function controls."
            )
        if any(item.engineered_neighborhood for item in items):
            design.warnings.append(
                "Mutation lies in an engineered-structure mismatch neighbourhood; target-sequence structural review is required."
            )
            if design.design_role == "STAGE1_PROBE" and design.status == "RECOMMENDED_LOW_PREDICTED_RISK":
                design.status = "REVIEW_REQUIRED"


def design_regions(
    regions: list[RegionRecord],
    records: list[ResidueRecord],
    protein: ProteinContext,
    inventory: EngineInventory,
    structure_pdb: Any | None,
    workdir: Any,
    original_cds: str | None = None,
    minimum_combo_sites: int | None = None,
    maximum_combo_sites: int | None = None,
) -> tuple[list[DesignRecord], dict[str, Any]]:
    by_position = {record.uniprot_pos: record for record in records}
    designs: list[DesignRecord] = []
    all_scans = 0
    foldx_repaired = None
    if inventory.foldx_path and structure_pdb:
        foldx_repaired = repair_with_foldx(inventory.foldx_path, structure_pdb, workdir)
    coding = _CodingCache(protein, original_cds)

    plm_llrs: dict[tuple[int, str], float] | None = None
    plm_notes: list[str] = []
    if os.environ.get("ECD_ESM_MODEL", "").strip():
        from .esm_scores import open_scorer

        scorer, reason = open_scorer()
        if scorer is None:
            plm_notes.append(reason)
        else:
            try:
                plm_llrs = scorer.substitution_llrs(protein.ecd_sequence)
                plm_notes.append(
                    f"PLM zero-shot evidence enabled ({scorer.model_name}); consensus veto at LLR < -3.0."
                )
            except Exception as exc:
                plm_notes.append(f"PLM scoring failed ({type(exc).__name__}); continuing without it.")

    probes: dict[int, str] = {}
    downgrades: list[str] = []
    pareto_fronts: dict[str, list[dict[str, Any]]] = {}

    for region in regions:
        # `residues` remains the legacy-compatible active footprint.  Compact
        # surface-face regions additionally expose all mutable face positions,
        # but Stage-1 must stay inside their local footprint.
        active_positions = (
            region.mutable_residues
            if region.partition_mode == "compact_surface_face" and region.mutable_residues
            else region.footprint_residues or region.residues
        )
        items = [by_position[position] for position in active_positions if position in by_position]
        max_combo_spread = (
            config.COMPACT_COMBO_MAX_SPREAD_A
            if region.partition_mode == "compact_surface_face"
            else config.MAX_COMBO_SPREAD_A
        )
        if not region.auto_design_allowed:
            designs.append(
                _make_design(
                    f"REGION_{region.id:02d}",
                    region.name,
                    region.id,
                    [],
                    "PAE_REVIEW_REQUIRED",
                    inventory,
                    protein,
                    coding,
                    "Predicted relative geometry failed the median intra-region PAE <= 10 A gate.",
                    ["No automatic probe was generated; this patch is a declared blind spot."],
                    role="STAGE1_PROBE",
                )
            )
            probes[region.id] = f"REGION_{region.id:02d}"
            continue

        key_site_limit = (
            config.COMPACT_KEY_SITES_PER_FACE
            if region.partition_mode == "compact_surface_face"
            else config.KEY_SITES_PER_PATCH
        )
        key_sites = select_key_sites(items, limit=key_site_limit)
        scans = {record.uniprot_pos: saturation_scan(record, protein) for record in key_sites}
        key_sites, scans = _ensure_anchor_site(items, key_sites, scans, protein)
        all_scans += sum(len(values) for values in scans.values())
        _apply_plm_evidence(scans, plm_llrs)
        foldx_single_ok = False
        if inventory.foldx_path and foldx_repaired:
            foldx_single_ok = _foldx_filter_single_choices(
                key_sites, scans, by_position, inventory.foldx_path, foldx_repaired, workdir
            )
        measured = bool(inventory.foldx_path and foldx_repaired)
        if region.partition_mode == "compact_surface_face" and region.candidate_footprints:
            local_fronts: list[tuple[float, float, list[MutationChoice]]] = []
            downgrade_notes = []
            for footprint in region.candidate_footprints:
                footprint_set = set(footprint)
                local_sites = [
                    record for record in key_sites
                    if record.uniprot_pos in footprint_set
                ]
                local_scans = {
                    position: values for position, values in scans.items()
                    if position in footprint_set
                }
                candidate_front, candidate_notes = _enumerate_combinations(
                    local_sites, local_scans, by_position, measured, max_combo_spread,
                    minimum_combo_sites, maximum_combo_sites,
                )
                local_fronts.extend(candidate_front)
                downgrade_notes.extend(candidate_notes)
            front = _pareto_front(local_fronts)
        else:
            front, downgrade_notes = _enumerate_combinations(
                key_sites, scans, by_position, measured, max_combo_spread,
                minimum_combo_sites, maximum_combo_sites,
            )
        for note in downgrade_notes:
            downgrades.append(f"Patch {region.id}: {note}")

        if not front:
            allowed_singles = sorted(
                (choice for values in scans.values() for choice in values if choice.allowed),
                key=lambda choice: (-_probe_utility(choice), choice.stability_score, choice.position, choice.mutant),
            )
            if region.partition_mode == "compact_surface_face" and allowed_singles:
                best_single = allowed_singles[0]
                primary_id = f"REGION_{region.id:02d}"
                designs.append(
                    _make_design(
                        primary_id,
                        region.name,
                        region.id,
                        [best_single],
                        "RECOMMENDED_SINGLE_SITE_PROBE",
                        inventory,
                        protein,
                        coding,
                        "No spatially valid multi-site probe passed; selected the strongest FoldX-cleared local single-site perturbation.",
                        [
                            "Single-site Stage-1 fallback: a negative result cannot exclude the full surface zone.",
                            "Pair with another independent single-site construct when broader regional coverage is required.",
                        ],
                        role="STAGE1_PROBE",
                        saturation=scans,
                    )
                )
                matching = next(
                    (
                        footprint for footprint in region.candidate_footprints
                        if best_single.position in footprint
                    ),
                    [best_single.position],
                )
                if matching not in region.candidate_footprints:
                    region.candidate_footprints.append(matching)
                region.footprint_residues = matching
                region.residues = matching
                region.footprint_center = best_single.position
                probes[region.id] = primary_id
                downgrades.append(
                    f"Surface zone {region.id}: multi-site probe unavailable; emitted a FoldX-cleared single-site Stage-1 probe."
                )
                continue
            status = (
                "FOLDX_SINGLE_FILTER_FAILED"
                if inventory.foldx_path and not foldx_single_ok
                else "NO_LOW_RISK_COMBINATION"
            )
            designs.append(
                _make_design(
                    f"REGION_{region.id:02d}",
                    region.name,
                    region.id,
                    [],
                    status,
                    inventory,
                    protein,
                    coding,
                    "No 2-6-site combination passed the deterministic chemistry and structural-risk constraints.",
                    ["Manual review is required; the workflow did not fabricate a fallback design.",
                     "This patch is a declared blind spot in the Stage-1 plan."],
                    role="STAGE1_PROBE",
                    saturation=scans,
                )
            )
            probes[region.id] = f"REGION_{region.id:02d}"
            continue

        # Pareto selection: the front spans the conservative->aggressive
        # trade-off; the knee point is the primary probe, the aggressive end
        # becomes its backup so the folding-control fallback is pre-designed.
        pareto_fronts[str(region.id)] = [
            {
                "mutations": [choice.notation for choice in mutations],
                "risk": risk,
                "disruption": disruption,
            }
            for risk, disruption, mutations in front
        ]
        best = _pareto_knee(front)
        if region.partition_mode == "compact_surface_face":
            best_positions = {choice.position for choice in best}
            matching = next(
                (
                    footprint for footprint in region.candidate_footprints
                    if best_positions.issubset(set(footprint))
                ),
                sorted(best_positions),
            )
            region.footprint_residues = matching
            region.residues = matching
            region.footprint_center = next(iter(sorted(best_positions)), region.footprint_center)
        alternate = front[-1][2] if len(front) > 1 and front[-1][2] is not best else None

        primary_id = f"REGION_{region.id:02d}"
        designs.append(
            _make_design(
                primary_id,
                region.name,
                region.id,
                best,
                "RECOMMENDED_LOW_PREDICTED_RISK",
                inventory,
                protein,
                coding,
                (
                    f"Stage-1 probe: knee of the exhaustive Pareto front ({len(front)} non-dominated combinations; "
                    f"{len(best)} substitutions, spacing >= {config.MIN_COMBO_SPACING_A:g} A, "
                    f"spread <= {max_combo_spread:g} A). Alternates cover the aggressive end."
                ),
                ["Computational low-risk classification is not experimental proof of preserved folding."],
                role="STAGE1_PROBE",
                saturation=scans,
            )
        )
        probes[region.id] = primary_id

        if alternate is not None:
            designs.append(
                _make_design(
                    f"ALT_{region.id:02d}",
                    f"{region.name} (alternate probe)",
                    None,
                    alternate,
                    "ALTERNATE_PROBE",
                    inventory,
                    protein,
                    coding,
                    "Aggressive end of the Pareto front: maximal disruption at higher stability risk.",
                    ["Run when the primary probe fails its folding/expression controls."],
                    role="ALTERNATE_PROBE",
                    parent_id=primary_id,
                )
            )

        # Stage-2 deconvolution: one construct per probe mutation.
        for choice in best:
            designs.append(
                _make_design(
                    f"SING_{region.id:02d}_{choice.position}",
                    f"{region.name} single {choice.notation}",
                    None,
                    [choice],
                    "STAGE2_SINGLE",
                    inventory,
                    protein,
                    coding,
                    "Stage-2 deconvolution preset: isolates one probe mutation to pinpoint the contact residue.",
                    ["Order and assay only after the Stage-1 probe loses binding."],
                    role="STAGE2_SINGLE",
                    parent_id=primary_id,
                )
            )
        # Confirmatory pairs: for epitopes where only the combination loses
        # binding, two-site subsets localise the interacting pair.
        pair_candidates = []
        for left, right in itertools.combinations(best, 2):
            left_record = by_position[left.position]
            right_record = by_position[right.position]
            if _distance(left_record, right_record) >= config.MIN_COMBO_SPACING_A:
                pair_candidates.append((left, right))
        pair_candidates.sort(key=lambda pair: -(pair[0].disruption_score + pair[1].disruption_score))
        for left, right in pair_candidates[: config.MAX_PAIRWISE_CONFIRMS]:
            designs.append(
                _make_design(
                    f"PAIR_{region.id:02d}_{left.position}_{right.position}",
                    f"{region.name} pair {left.notation}+{right.notation}",
                    None,
                    [left, right],
                    "STAGE2_PAIRWISE",
                    inventory,
                    protein,
                    coding,
                    "Confirmatory two-site subset for combination-only binding loss (redundant-epitope case).",
                    ["Order only when Stage-2 singles individually do not reproduce the loss."],
                    role="STAGE2_PAIRWISE",
                    parent_id=primary_id,
                )
            )

    _annotate_design_risks(designs, by_position)
    foldx_diagnostics: dict[str, Any] = {}
    scored_designs = [design for design in designs if design.mutations]
    if inventory.foldx_path and structure_pdb and foldx_repaired and scored_designs:
        # Two independent batches: the ordering panel gets a robust 5-run
        # median, while Stage-2 presets are only screened (runs=1) -- they are
        # ordered selectively after a Stage-1 hit, and a 5-run median can be
        # requested then.  This cuts the dominant FoldX cost roughly in half.
        panel_designs = [
            design for design in scored_designs
            if design.design_role in {"STAGE1_PROBE", "ALTERNATE_PROBE"}
        ]
        stage2_designs = [
            design for design in scored_designs
            if design.design_role in {"STAGE2_SINGLE", "STAGE2_PAIRWISE"}
        ]
        # Stage-2 singles are identical to substitutions already measured in
        # the FoldX one-run screen.  Reuse that audited measurement instead
        # of launching the same BuildModel job again.
        reused_single_designs = [
            design
            for design in stage2_designs
            if design.design_role == "STAGE2_SINGLE"
            and len(design.mutations) == 1
            and design.mutations[0].foldx_ddg is not None
        ]
        stage2_to_score = [design for design in stage2_designs if design not in reused_single_designs]

        def _groups(selected: list[DesignRecord]) -> list[list[tuple[str, int, str, str]]]:
            return [
                [
                    (
                        mutation.wt,
                        by_position[mutation.position].auth_seq_id,
                        by_position[mutation.position].chain_id,
                        mutation.mutant,
                    )
                    for mutation in design.mutations
                ]
                for design in selected
            ]

        panel_scores = score_many_with_foldx(
            inventory.foldx_path,
            foldx_repaired,
            _groups(panel_designs),
            workdir,
            runs=config.FOLDX_COMBINATION_RUNS,
            already_repaired=True,
        )
        stage2_scores = score_many_with_foldx(
            inventory.foldx_path,
            foldx_repaired,
            _groups(stage2_to_score),
            workdir,
            runs=config.FOLDX_STAGE2_RUNS,
            already_repaired=True,
        )
        reused_scores = [design.mutations[0].foldx_ddg for design in reused_single_designs]
        all_scores = list(panel_scores) + list(reused_scores) + list(stage2_scores)
        n_failed = sum(1 for item in all_scores if item is None)
        total = len(panel_designs) + len(stage2_designs)
        foldx_diagnostics = {
            "batch_constructs": total,
            "panel_constructs": len(panel_designs),
            "panel_runs": config.FOLDX_COMBINATION_RUNS,
            "stage2_constructs": len(stage2_designs),
            "stage2_runs": config.FOLDX_STAGE2_RUNS,
            "batch_scores_returned": total - n_failed,
            "batch_scores_failed": n_failed,
            "stage2_single_scores_reused": len(reused_single_designs),
            "single_gate_kcal_mol": {
                "exposed": config.FOLDX_GATE_EXPOSED,
                "edge": config.FOLDX_GATE_EDGE,
                "default": config.FOLDX_GATE_DEFAULT,
            },
            "combination_gate_kcal_mol": config.FOLDX_COMBINATION_GATE,
            "all_failed": bool(total) and n_failed == total,
        }
        panel_ids = {id(design) for design in panel_designs}
        scored_pairs = (
            list(zip(panel_designs, panel_scores))
            + list(zip(reused_single_designs, reused_scores))
            + list(zip(stage2_to_score, stage2_scores))
        )
        for design, foldx_score in scored_pairs:
            if foldx_score is None:
                # score_origin stays "proxy": the number on the record is the
                # proxy, never kcal/mol.
                design.status = "FOLDX_REVIEW_REQUIRED"
                design.score_origin = "proxy"
                design.stability_engine = "builtin structural-risk proxy (FoldX failed)"
                design.warnings.append(
                    "FoldX BuildModel output could not be parsed or completed too few runs; no stability conclusion was assigned."
                )
                continue
            single = len(design.mutations) == 1
            gate = (
                config.foldx_single_gate(
                    by_position[design.mutations[0].position].rsasa
                )
                if single
                else config.FOLDX_COMBINATION_GATE
            )
            design.combined_stability_score = foldx_score
            design.score_origin = "foldx"
            design.stability_engine = "FoldX 5.1 BuildModel median delta-delta-G"
            runs_used = (
                config.FOLDX_COMBINATION_RUNS
                if id(design) in panel_ids
                else config.FOLDX_STAGE2_RUNS
            )
            design.warnings.append(
                f"FoldX BuildModel delta-delta-G: {foldx_score:.3f} kcal/mol over {runs_used} run(s)."
            )
            if design.design_role in {"STAGE2_SINGLE", "STAGE2_PAIRWISE"} and runs_used < config.FOLDX_COMBINATION_RUNS:
                design.warnings.append(
                    "Stage-2 preset was screened at runs=1; request a 5-run median before ordering decisions."
                )
            if foldx_score > gate:
                design.status = "FOLDX_REVIEW_REQUIRED"
                design.warnings.append(
                    f"FoldX delta-delta-G exceeds the configured +{gate:.2f} kcal/mol "
                    f"{'single-point (exposure-linked)' if single else 'combination'} gate."
                )
    elif inventory.foldx_path and structure_pdb and not foldx_repaired:
        for design in scored_designs:
            design.status = "FOLDX_REVIEW_REQUIRED"
            design.warnings.append("FoldX was detected but RepairPDB did not produce a reusable repaired structure.")

    # EvoEF2 optional refinement for Stage-1 probes (unchanged policy).
    if inventory.evoef2_path and structure_pdb:
        for design in designs:
            if design.design_role != "STAGE1_PROBE" or not design.mutations:
                continue
            external_mutations = [
                (
                    mutation.wt,
                    by_position[mutation.position].auth_seq_id,
                    by_position[mutation.position].chain_id,
                    mutation.mutant,
                )
                for mutation in design.mutations
            ]
            evoef_score = score_with_evoef2(
                inventory.evoef2_path, structure_pdb, external_mutations, workdir
            )
            if evoef_score is not None:
                design.combined_stability_score = evoef_score
                design.score_origin = "evoef2"
                design.stability_engine = "EvoEF2 mutant-minus-WT energy"
            else:
                design.warnings.append("EvoEF2 was detected but its output could not be parsed; proxy score retained.")

    if not inventory.evoef2_path and not inventory.foldx_path:
        for design in designs:
            if design.mutations:
                design.warnings.append(
                    "Neither FoldX nor EvoEF2 is configured; the reported score is a dimensionless within-target structural-risk proxy, not kcal/mol."
                )

    # --- epitope mapping experiment plan ------------------------------------
    matrix = _coverage_matrix(designs, regions, by_position)
    blind_spots = [
        {"patch_id": region.id, "reason": "PAE gate" if not region.auto_design_allowed else "no combination passed the stability gates"}
        for region in regions
        if not any(design.region_id == region.id and design.mutations for design in designs)
    ]
    plan = {
        "assay_note": (
            "Stage 1: order every STAGE1_PROBE (+ ALTERNATE_PROBE for HIGH patches), assay binding of all antibodies "
            "against each construct with expression/folding controls. A loss pattern localises the epitope to the "
            "perturbed patch; identical loss patterns across antibodies are binning evidence. "
            "Stage 2: for hits, order the parent probe's STAGE2_SINGLE presets; add STAGE2_PAIRWISE presets when "
            "singles do not reproduce the loss (redundant epitope). Confirm by restoring the WT residue (back-mutation)."
        ),
        "selection_policy": (
            "Pareto front over (disruption D, stability risk) from exhaustive enumeration; primary probe = knee "
            "(max D per unit risk, 1.0 kcal/mol noise anchor), alternate = aggressive end of the same front."
        ),
        "pareto_fronts": pareto_fronts,
        "stage1_probes": [
            {
                "design_id": design.id,
                "patch_id": design.region_id,
                "mutations": [choice.notation for choice in design.mutations],
                "anchor": design.anchor_notation,
                "disruption_score": design.disruption_score,
                "status": design.status,
            }
            for design in designs
            if design.design_role == "STAGE1_PROBE"
        ],
        "coverage_matrix": matrix,
        "blind_spots": blind_spots,
        "stage2_presets": [
            {
                "design_id": design.id,
                "parent_design_id": design.parent_design_id,
                "role": design.design_role,
                "mutations": [choice.notation for choice in design.mutations],
            }
            for design in designs
            if design.design_role in {"STAGE2_SINGLE", "STAGE2_PAIRWISE", "ALTERNATE_PROBE"}
        ],
        "design_downgrades": downgrades,
    }

    return designs, {
        "regions_considered": len(regions),
        "recommended_designs": sum(
            item.status in {
                "RECOMMENDED_LOW_PREDICTED_RISK",
                "RECOMMENDED_MULTI_SITE_PROBE",
                "RECOMMENDED_SINGLE_SITE_PROBE",
            }
            for item in designs
        ),
        "in_silico_substitutions_evaluated": all_scans,
        "final_library_size": len(designs),
        "foldx_diagnostics": foldx_diagnostics,
        "epitope_mapping_plan": plan,
        "plm_notes": plm_notes,
        "score_interpretation": (
            "FoldX 5.1 single-point and combination ΔΔG in kcal/mol when available (single-point gate linked to side-chain exposure); "
            "otherwise a dimensionless within-target proxy"
        ),
    }
