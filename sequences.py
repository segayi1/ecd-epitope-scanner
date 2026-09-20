from __future__ import annotations

import re
from collections.abc import Iterable

from Bio.Seq import Seq

from .models import MutationChoice, ProteinContext
from .util import PipelineError

HUMAN_PREFERRED_CODONS = {
    "A": "GCC",
    "R": "CGC",
    "N": "AAC",
    "D": "GAC",
    "C": "TGC",
    "Q": "CAG",
    "E": "GAG",
    "G": "GGC",
    "H": "CAC",
    "I": "ATC",
    "L": "CTG",
    "K": "AAG",
    "M": "ATG",
    "F": "TTC",
    "P": "CCC",
    "S": "AGC",
    "T": "ACC",
    "W": "TGG",
    "Y": "TAC",
    "V": "GTG",
}


def normalize_cds(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z]", "", value or "").upper().replace("U", "T")
    if cleaned and (set(cleaned) - set("ACGT")):
        raise PipelineError("CDS contains characters other than A/C/G/T")
    if cleaned and len(cleaned) % 3:
        raise PipelineError("CDS length must be divisible by three")
    return cleaned


def translate_cds(cds: str) -> str:
    return str(Seq(cds).translate()).rstrip("*")


def reverse_translate(sequence: str) -> str:
    try:
        return "".join(HUMAN_PREFERRED_CODONS[aa] for aa in sequence)
    except KeyError as exc:
        raise PipelineError(f"Cannot reverse-translate amino acid {exc.args[0]!r}") from exc


def optimize_human_coding_fragment(sequence: str) -> tuple[str, str]:
    initial = reverse_translate(sequence)
    try:
        from dnachisel import (
            AvoidPattern,
            CodonOptimize,
            DnaOptimizationProblem,
            EnforceGCContent,
            EnforceTranslation,
        )

        constraints = [
            EnforceTranslation(translation=sequence),
            EnforceGCContent(mini=0.40, maxi=0.65, window=min(100, len(initial))),
            AvoidPattern("AAAAAAA"),
            AvoidPattern("CCCCCCC"),
            AvoidPattern("GGGGGGG"),
            AvoidPattern("TTTTTTT"),
        ]
        problem = DnaOptimizationProblem(
            sequence=initial,
            constraints=constraints,
            objectives=[CodonOptimize(species="h_sapiens_9606", method="use_best_codon")],
            logger=None,
        )
        problem.resolve_constraints()
        problem.optimize()
        optimized = str(problem.sequence).upper()
        if translate_cds(optimized) == sequence:
            return optimized, "DNA Chisel HEK293/human codon optimization"
    except Exception:
        pass
    return initial, "deterministic HEK293/human preferred codons (DNA Chisel fallback)"


def apply_mutations(
    sequence: str,
    mutations: Iterable[MutationChoice],
    sequence_start: int = 1,
) -> str:
    output = list(sequence)
    for mutation in mutations:
        index = mutation.position - sequence_start
        if not (0 <= index < len(output)):
            raise PipelineError(f"Mutation {mutation.notation} falls outside the output sequence")
        if output[index] != mutation.wt:
            raise PipelineError(
                f"Mutation {mutation.notation} does not match sequence residue {output[index]} at {mutation.position}"
            )
        if mutation.mutant == mutation.wt:
            raise PipelineError(f"WT-to-WT pseudo-mutation is not allowed: {mutation.notation}")
        output[index] = mutation.mutant
    return "".join(output)


def coding_fragment(
    protein: ProteinContext,
    mutant_ecd_sequence: str,
    original_cds: str | None = None,
) -> tuple[str, str]:
    if len(mutant_ecd_sequence) != protein.ecd_length:
        raise PipelineError("Mutant ECD sequence length differs from the canonical ECD")
    source = "deterministic HEK293/human preferred codons"
    if original_cds:
        normalized = normalize_cds(original_cds)
        translated = translate_cds(normalized)
        if translated == protein.sequence:
            fragment = normalized[(protein.ecd_start - 1) * 3 : protein.ecd_end * 3]
        elif translated == protein.ecd_sequence:
            fragment = normalized
        else:
            raise PipelineError(
                "Uploaded CDS does not translate to either the canonical full protein or the selected ECD"
            )
        codons = [fragment[index : index + 3] for index in range(0, len(fragment), 3)]
        for index, (wt, mutant) in enumerate(zip(protein.ecd_sequence, mutant_ecd_sequence)):
            if wt != mutant:
                codons[index] = HUMAN_PREFERRED_CODONS[mutant]
        fragment = "".join(codons)
        source = "uploaded CDS with only mutated codons replaced"
    else:
        fragment, source = optimize_human_coding_fragment(mutant_ecd_sequence)

    if translate_cds(fragment) != mutant_ecd_sequence:
        raise PipelineError("Internal reverse-translation validation failed")
    return fragment, source
