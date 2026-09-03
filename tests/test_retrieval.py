from ecd_epitope_scanner.models import FeatureRange, StructureCandidate
from ecd_epitope_scanner.retrieval import find_public_af3_candidates, parse_uniprot_entry
from ecd_epitope_scanner.structure_analysis import _feature_positions


def minimal_uniprot():
    return {
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "primaryAccession": "PTEST1",
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Test receptor"}}},
        "organism": {"scientificName": "Homo sapiens"},
        "sequence": {"value": "M" * 10 + "ACDEFGHIKLMNPQRSTVWY" * 5 + "L" * 24 + "K" * 10},
        "features": [
            {
                "type": "Topological domain",
                "description": "Extracellular",
                "location": {"start": {"value": 11}, "end": {"value": 110}},
            },
            {
                "type": "Transmembrane",
                "location": {"start": {"value": 111}, "end": {"value": 134}},
            },
            {
                "type": "Domain",
                "description": "Ig-like domain",
                "location": {"start": {"value": 20}, "end": {"value": 90}},
            },
        ],
        "uniProtKBCrossReferences": [{"database": "PDB", "id": "1ABC"}],
    }


def test_uniprot_ecd_and_topology_are_derived_from_annotations():
    protein = parse_uniprot_entry(minimal_uniprot())
    assert protein.ecd_start == 11
    assert protein.ecd_end == 110
    assert protein.membrane_anchor_position == 110
    assert protein.topology.startswith("type-I")
    assert protein.reviewed
    assert protein.pdb_ids == ["1abc"]


def test_non_single_pass_entry_is_rejected():
    entry = minimal_uniprot()
    entry["features"] = [item for item in entry["features"] if item["type"] != "Transmembrane"]
    try:
        parse_uniprot_entry(entry)
    except Exception as exc:
        assert "single-pass" in str(exc)
    else:
        raise AssertionError("Expected a single-pass validation error")


class FakeClient:
    def text(self, url):
        return "Model generated with AlphaFold 3 for publication"


def test_strict_af3_does_not_accept_alphafold_db_as_af3():
    afdb = StructureCandidate(
        identifier="AF-PTEST1-F1",
        category="AB-INITIO",
        provider="AlphaFold DB",
        model_url="https://example/afdb.cif",
    )
    archive = StructureCandidate(
        identifier="ma-test",
        category="AB-INITIO",
        provider="ModelArchive",
        model_url="https://example/model.cif",
        page_url="https://example/model",
    )
    accepted = find_public_af3_candidates([afdb, archive], FakeClient())
    assert [item.identifier for item in accepted] == ["ma-test"]
    assert archive.explicit_af3


def test_disulfide_feature_protects_only_the_two_bond_endpoints():
    protein = parse_uniprot_entry(minimal_uniprot())
    protein.features.append(FeatureRange("Disulfide bond", 30, 80, ""))
    assert _feature_positions(protein, [r"disulfide bond"]) == {30, 80}
