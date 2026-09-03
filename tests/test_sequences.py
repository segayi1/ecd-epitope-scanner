import pytest

from ecd_epitope_scanner.models import MutationChoice, ProteinContext
from ecd_epitope_scanner.sequences import apply_mutations, coding_fragment, translate_cds
from ecd_epitope_scanner.util import PipelineError


def protein():
    return ProteinContext(
        accession="PTEST1",
        protein_name="test",
        organism="human",
        sequence="MACDEFGHIKLMNPQRSTVWY",
        ecd_start=2,
        ecd_end=21,
        membrane_anchor_position=21,
        topology="test",
    )


def choice(position=2, wt="A", mutant="R"):
    return MutationChoice(position, wt, mutant, 0.2, 1.0, 0.1, True, [])


def test_ecd_mutation_and_reverse_translation_roundtrip():
    target = protein()
    mutant = apply_mutations(target.ecd_sequence, [choice()], target.ecd_start)
    assert mutant[0] == "R"
    dna, source = coding_fragment(target, mutant)
    assert translate_cds(dna) == mutant
    assert len(dna) == len(mutant) * 3
    assert "start" not in source.lower()


def test_wt_to_wt_and_reference_mismatch_are_rejected():
    target = protein()
    with pytest.raises(PipelineError, match="pseudo-mutation"):
        apply_mutations(target.ecd_sequence, [choice(mutant="A")], target.ecd_start)
    with pytest.raises(PipelineError, match="does not match"):
        apply_mutations(target.ecd_sequence, [choice(wt="V")], target.ecd_start)
