from __future__ import annotations

import itertools
from statistics import mean
from typing import Any

import numpy as np

from .models import DesignRecord, MutationChoice, ProteinContext, RegionRecord, ResidueRecord, to_dict
from .sequences import apply_mutations, coding_fragment
from .stability import (
    AMINO_ACIDS,
    EngineInventory,
    assign_percentiles,
    evaluate_substitution,
    repair_with_foldx,
    score_many_with_foldx,
    score_with_foldx,
    score_with_evoef2,
)


def _distance(left: ResidueRecord, right: ResidueRecord) -> float:
    if left.ca is None or right.ca is None:
        return float("inf")
    return float(np.linalg.norm(np.asarray(left.ca) - np.asarray(right.ca)))


def _site_priority(record: ResidueRecord) -> float:
    score = record.rsasa * 2.0 + min(1.2, record.sidechain_sasa / 100.0)
    score += 0.6 / (1 + record.contacts / 6.0)
    score += 0.25 if record.interface else 0.0
    score -= 0.35 if record.salt_bridge else 0.0
    score -= 0.25 if record.domain_interface else 0.0
    score -= 0.35 if record.glycan_shielded else 0.0
    score -= 0.45 if record.aa in "GP" else 0.0
    return score


def select_key_sites(items: list[ResidueRecord], limit: int = 8) -> list[ResidueRecord]:
    ranked = sorted(items, key=lambda item: (-_site_priority(item), item.uniprot_pos))
    selected: list[ResidueRecord] = []
    for record in ranked:
        if all(_distance(record, existing) >= 5.0 for existing in selected):
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


def _choice_utility(choice: MutationChoice) -> float:
    return choice.perturbation_score - 0.72 * choice.stability_score - 0.35 * choice.stability_percentile


def _best_choice(choices: list[MutationChoice]) -> MutationChoice | None:
    allowed = [choice for choice in choices if choice.allowed and choice.stability_percentile <= 0.7]
    if not allowed:
        return None
    return max(allowed, key=lambda choice: (_choice_utility(choice), -choice.stability_score, choice.mutant))


def _foldx_filter_single_choices(
    key_sites: list[ResidueRecord],
    scans: dict[int, list[MutationChoice]],
    records_by_position: dict[int, ResidueRecord],
    binary: str,
    repaired_pdb: Any,
    workdir: Any,
) -> bool:
    """Use FoldX PositionScan-equivalent BuildModel singles as the combination gate.

    The proxy saturation scan remains visible for audit, but only the top three
    proxy-ranked substitutions per key site are admitted to this FoldX screen.
    Unscored substitutions are explicitly disallowed from the final combination.
    """
    jobs: list[list[tuple[str, int, str, str]]] = []
    job_choices: list[MutationChoice] = []
    for record in key_sites:
        options = [choice for choice in scans[record.uniprot_pos] if choice.allowed]
        options.sort(key=lambda choice: (-_choice_utility(choice), choice.mutant))
        scored_positions = {choice.mutant for choice in options[:3]}
        for choice in scans[record.uniprot_pos]:
            if choice.mutant not in scored_positions:
                choice.allowed = False
                choice.stability_class = "not_scored"
                if "not evaluated by FoldX single-point screen" not in choice.reasons:
                    choice.reasons.append("not evaluated by FoldX single-point screen")
        for choice in options[:3]:
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
        runs=1,
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
        choice.foldx_ddg = score
        choice.stability_score = score
        if score < -0.46:
            choice.stability_class = "stabilizing"
        elif score <= 0.46:
            choice.stability_class = "neutral"
        else:
            choice.stability_class = "destabilizing"
        choice.allowed = choice.allowed and score <= 0.46
        choice.reasons.append(f"FoldX single-point ΔΔG={score:.3f} kcal/mol")
        if score > 0.46:
            choice.reasons.append("exceeds +0.46 kcal/mol single-point gate")
    for choices in scans.values():
        assign_percentiles(choices)
    return valid > 0


def _combined_proxy(
    mutations: list[MutationChoice], by_position: dict[int, ResidueRecord]
) -> float:
    risk = sum(item.stability_score for item in mutations)
    for left, right in itertools.combinations(mutations, 2):
        distance = _distance(by_position[left.position], by_position[right.position])
        if distance < 5.0:
            risk += 0.55
        elif distance < 7.0:
            risk += 0.2
        left_charge = 1 if left.mutant in "KR" else -1 if left.mutant in "DE" else 0
        right_charge = 1 if right.mutant in "KR" else -1 if right.mutant in "DE" else 0
        if distance < 8.0 and left_charge and left_charge == right_charge:
            risk += 0.15
    return round(risk, 4)


def _choose_combination(
    key_sites: list[ResidueRecord],
    scans: dict[int, list[MutationChoice]],
) -> tuple[list[MutationChoice], float | None]:
    by_position = {item.uniprot_pos: item for item in key_sites}
    options = [(record, _best_choice(scans[record.uniprot_pos])) for record in key_sites]
    options = [(record, choice) for record, choice in options if choice is not None]
    if len(options) < 2:
        return [], None
    options.sort(
        key=lambda pair: (
            -(_site_priority(pair[0]) + _choice_utility(pair[1])),  # type: ignore[arg-type]
            pair[0].uniprot_pos,
        )
    )
    maximum = min(5, len(options))
    for count in range(maximum, 1, -1):
        chosen: list[tuple[ResidueRecord, MutationChoice]] = []
        for record, choice in options:
            if not chosen or min(_distance(record, previous[0]) for previous in chosen) >= 5.5:
                chosen.append((record, choice))  # type: ignore[arg-type]
            if len(chosen) >= count:
                break
        if len(chosen) < count:
            chosen = options[:count]  # type: ignore[assignment]
        mutations = [choice for _, choice in chosen]
        proxy = _combined_proxy(mutations, by_position)
        if proxy / count <= 1.15 and proxy <= max(2.0, count * 1.05):
            return mutations, proxy
    return [], None


def design_regions(
    regions: list[RegionRecord],
    records: list[ResidueRecord],
    protein: ProteinContext,
    inventory: EngineInventory,
    structure_pdb: Any | None,
    workdir: Any,
    original_cds: str | None = None,
) -> tuple[list[DesignRecord], dict[str, Any]]:
    by_position = {record.uniprot_pos: record for record in records}
    designs: list[DesignRecord] = []
    all_scans = 0
    foldx_repaired = None
    if inventory.foldx_path and structure_pdb:
        foldx_repaired = repair_with_foldx(inventory.foldx_path, structure_pdb, workdir)
    for region in regions:
        items = [by_position[position] for position in region.residues if position in by_position]
        if not region.auto_design_allowed:
            designs.append(
                DesignRecord(
                    id=f"REGION_{region.id:02d}",
                    name=region.name,
                    region_id=region.id,
                    mutations=[],
                    status="PAE_REVIEW_REQUIRED",
                    stability_engine=inventory.primary_name,
                    combined_stability_score=None,
                    perturbation_score=0.0,
                    rationale="Predicted relative geometry failed the median intra-region PAE <= 10 A gate.",
                    warnings=["No automatic combination was generated."],
                    mutant_ecd_sequence=protein.ecd_sequence,
                    mutant_coding_sequence=coding_fragment(protein, protein.ecd_sequence, original_cds)[0],
                )
            )
            continue
        key_sites = select_key_sites(items, limit=4)
        scans = {record.uniprot_pos: saturation_scan(record, protein) for record in key_sites}
        all_scans += sum(len(values) for values in scans.values())
        foldx_single_ok = False
        if inventory.foldx_path and foldx_repaired:
            foldx_single_ok = _foldx_filter_single_choices(
                key_sites,
                scans,
                by_position,
                inventory.foldx_path,
                foldx_repaired,
                workdir,
            )
        mutations, proxy_score = _choose_combination(key_sites, scans)
        if not mutations:
            status = "FOLDX_SINGLE_FILTER_FAILED" if inventory.foldx_path and not foldx_single_ok else "NO_LOW_RISK_COMBINATION"
            designs.append(
                DesignRecord(
                    id=f"REGION_{region.id:02d}",
                    name=region.name,
                    region_id=region.id,
                    mutations=[],
                    status=status,
                    stability_engine=inventory.primary_name,
                    combined_stability_score=None,
                    perturbation_score=0.0,
                    rationale="No 2-5-site combination passed the deterministic chemistry and structural-risk constraints.",
                    warnings=["Manual review is required; the workflow did not fabricate a fallback design."],
                    mutant_ecd_sequence=protein.ecd_sequence,
                    mutant_coding_sequence=coding_fragment(protein, protein.ecd_sequence, original_cds)[0],
                    saturation_scan={
                        str(position): [to_dict(choice) | {"notation": choice.notation} for choice in values]
                        for position, values in scans.items()
                    },
                )
            )
            continue

        mutant_sequence = apply_mutations(protein.ecd_sequence, mutations, protein.ecd_start)
        mutant_dna, dna_source = coding_fragment(protein, mutant_sequence, original_cds)
        engine_name = inventory.primary_name
        combined_score = proxy_score
        design_status = "RECOMMENDED_LOW_PREDICTED_RISK"
        warnings = [
            "Computational low-risk classification is not experimental proof of preserved folding.",
            f"Coding fragment source: {dna_source}.",
        ]
        if inventory.evoef2_path and structure_pdb:
            external_mutations = [
                (
                    mutation.wt,
                    by_position[mutation.position].auth_seq_id,
                    by_position[mutation.position].chain_id,
                    mutation.mutant,
                )
                for mutation in mutations
            ]
            evoef_score = score_with_evoef2(
                inventory.evoef2_path,
                structure_pdb,
                external_mutations,
                workdir,
            )
            if evoef_score is not None:
                combined_score = evoef_score
                engine_name = "EvoEF2 mutant-minus-WT energy"
            else:
                warnings.append("EvoEF2 was detected but its output could not be parsed; proxy score retained.")
        elif inventory.foldx_path and structure_pdb and not foldx_repaired:
            warnings.append("FoldX was detected but RepairPDB did not produce a reusable repaired structure; no FoldX score was assigned.")
        elif not inventory.evoef2_path and not inventory.foldx_path:
            warnings.append(
                "Neither FoldX nor EvoEF2 is configured; the reported score is a dimensionless within-target structural-risk proxy, not kcal/mol."
            )
        designs.append(
            DesignRecord(
                id=f"REGION_{region.id:02d}",
                name=region.name,
                region_id=region.id,
                mutations=mutations,
                status=design_status,
                stability_engine=engine_name,
                combined_stability_score=combined_score,
                perturbation_score=round(sum(item.perturbation_score for item in mutations), 4),
                rationale=(
                    f"{len(mutations)} spatially distributed substitutions maximize charge/volume/hydropathy change "
                    "while satisfying protected-site, sequon and combined-risk constraints."
                ),
                warnings=warnings,
                mutant_ecd_sequence=mutant_sequence,
                mutant_coding_sequence=mutant_dna,
                saturation_scan={
                    str(position): [
                        {
                            "position": choice.position,
                            "wt": choice.wt,
                            "mutant": choice.mutant,
                            "stability_score": choice.stability_score,
                            "perturbation_score": choice.perturbation_score,
                            "stability_percentile": choice.stability_percentile,
                            "foldx_ddg": choice.foldx_ddg,
                            "stability_class": choice.stability_class,
                            "allowed": choice.allowed,
                            "reasons": choice.reasons,
                            "notation": choice.notation,
                        }
                        for choice in values
                    ]
                    for position, values in scans.items()
                },
            )
        )
    if inventory.foldx_path and structure_pdb and foldx_repaired:
        scored_designs = [design for design in designs if design.mutations]
        groups = [
            [
                (
                    mutation.wt,
                    by_position[mutation.position].auth_seq_id,
                    by_position[mutation.position].chain_id,
                    mutation.mutant,
                )
                for mutation in design.mutations
            ]
            for design in scored_designs
        ]
        foldx_scores = score_many_with_foldx(
            inventory.foldx_path,
            foldx_repaired,
            groups,
            workdir,
            runs=5,
            already_repaired=True,
        )
        for design, foldx_score in zip(scored_designs, foldx_scores):
            if foldx_score is None:
                design.status = "FOLDX_REVIEW_REQUIRED"
                design.warnings.append("FoldX BuildModel output could not be parsed; no stability conclusion was assigned.")
                continue
            design.combined_stability_score = foldx_score
            design.stability_engine = "FoldX 5.1 BuildModel median delta-delta-G"
            design.warnings.append(f"FoldX BuildModel median delta-delta-G: {foldx_score:.3f} kcal/mol over 5 runs.")
            if foldx_score > 0.46:
                design.status = "FOLDX_REVIEW_REQUIRED"
                design.warnings.append("FoldX combination delta-delta-G exceeds the configured +0.46 kcal/mol gate.")
    elif inventory.foldx_path and structure_pdb and not foldx_repaired:
        for design in designs:
            if design.mutations:
                design.status = "FOLDX_REVIEW_REQUIRED"
                design.warnings.append("FoldX was detected but RepairPDB did not produce a reusable repaired structure.")
    return designs, {
        "regions_considered": len(regions),
        "recommended_designs": sum(item.status == "RECOMMENDED_LOW_PREDICTED_RISK" for item in designs),
        "in_silico_substitutions_evaluated": all_scans,
        "final_library_size": len(designs),
        "score_interpretation": (
            "FoldX 5.1 single-point and combination ΔΔG in kcal/mol when available; otherwise a dimensionless within-target proxy"
        ),
    }
