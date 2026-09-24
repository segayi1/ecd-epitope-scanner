from __future__ import annotations

from ..models import ProteinContext, ResidueRecord
from .models import MutationOption, ProbeSector
from .structure_audit import residue_id

_CHARGED = set("DEKR")
_AROMATIC = set("FWY")
_HYDROPHOBIC = set("AILM V".replace(" ", ""))
_POLAR = set("STNQ")


def _menu(record: ResidueRecord) -> list[tuple[str, str, float]]:
    aa = record.aa
    options: list[tuple[str, str, float]] = []
    if aa in _CHARGED:
        options.append(("charge_removal", "A", 0.72))
        options.append(("charge_neutralization", "Q" if aa in "DE" else "N", 0.48))
    elif aa in _AROMATIC:
        options.append(("aromatic_removal", "A", 0.76))
        options.append(("sidechain_shortening", "S", 0.52))
    elif aa in _HYDROPHOBIC:
        options.append(("hydrophobicity_reduction", "A", 0.55))
        options.append(("polar_substitution", "S", 0.46))
    elif aa in _POLAR:
        options.append(("sidechain_shortening", "A", 0.42))
        options.append(("polarity_change", "V", 0.45))
    else:
        options.append(("alanine_scan", "A", 0.40))
    return [(kind, mutant, strength) for kind, mutant, strength in options if mutant != aa and mutant != "C"]


def _creates_glycosylation_sequon(sequence: str, position: int, mutant: str) -> bool:
    chars = list(sequence)
    index = position - 1
    if not 0 <= index < len(chars):
        return False
    original = set()
    for start in range(max(0, index - 2), min(len(chars) - 2, index + 1)):
        if chars[start] == "N" and chars[start + 1] != "P" and chars[start + 2] in "ST":
            original.add(start)
    chars[index] = mutant
    for start in range(max(0, index - 2), min(len(chars) - 2, index + 1)):
        if chars[start] == "N" and chars[start + 1] != "P" and chars[start + 2] in "ST":
            if start not in original:
                return True
    return False


def build_mutation_options(
    sector: ProbeSector,
    records: dict[int, ResidueRecord],
    protein: ProteinContext,
    options_per_site: int = 3,
    allowed_positions: set[int] | None = None,
) -> tuple[list[MutationOption], dict[int, list[str]]]:
    options: list[MutationOption] = []
    rejected: dict[int, list[str]] = {}
    for position in sector.residue_ids:
        record = records.get(position)
        if record is None:
            continue
        if allowed_positions is not None and position not in allowed_positions:
            rejected[position] = ["surface-unit-not-designable"]
            continue
        reasons = list(record.hard_exclusions)
        if record.disulfide or record.aa == "C":
            reasons.append("disulfide-cysteine")
        if record.glycosylation_site or record.glycosylation_sequon:
            reasons.append("glycosylation-sequon")
        if record.flexible_terminal:
            reasons.append("flexible-terminal")
        if record.interface:
            reasons.append("native-interface")
        if record.domain_interface:
            reasons.append("domain-interface")
        if record.salt_bridge:
            reasons.append("salt-bridge-core")
        if record.aa in {"G", "P"}:
            reasons.append("backbone-conformation-sensitive")
        if reasons:
            rejected[position] = sorted(set(reasons))
            continue
        for perturbation, mutant, strength in _menu(record)[:options_per_site]:
            if _creates_glycosylation_sequon(protein.sequence, position, mutant):
                rejected.setdefault(position, []).append("introduces-glycosylation-sequon")
                continue
            risk = 0.10
            risk += 0.18 if record.salt_bridge else 0.0
            risk += 0.15 if record.interface else 0.0
            risk += 0.12 if record.domain_interface else 0.0
            risk += 0.12 if record.glycan_shielded else 0.0
            risk += 0.08 if record.engineered_neighborhood else 0.0
            risk += max(0.0, 0.35 - record.rsasa) * 0.8
            options.append(
                MutationOption(
                    residue=residue_id(record), position=position, wild_type=record.aa,
                    mutant_aa=mutant, perturbation_class=perturbation,
                    rule_risk=round(min(1.0, risk), 4),
                    chemical_perturbation_strength=strength,
                    stability_value=None, stability_unit=None,
                    evidence_level="PROXY_ONLY",
                )
            )
    options.sort(key=lambda item: (-item.chemical_perturbation_strength, item.rule_risk, item.position, item.mutant_aa))
    return options, rejected
