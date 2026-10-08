from ecd_epitope_scanner.scripts.map_cross_protein_validation import align_polymer


def test_semiglobal_mapping_preserves_absolute_reference_positions():
    mapping, audit = align_polymer("XXACDEFYY", "ACDEF")
    assert [mapping[index] for index in range(5)] == [3, 4, 5, 6, 7]
    assert audit["sequence_identity"] == 1.0


def test_internal_substitution_does_not_shift_downstream_residue_map():
    mapping, audit = align_polymer("ACDEFG", "ACDFFG")
    assert [mapping[index] for index in range(6)] == [1, 2, 3, 4, 5, 6]
    assert audit["sequence_identity"] == 5 / 6
