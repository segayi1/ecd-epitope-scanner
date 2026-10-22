from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, PDBIO, Residue, Structure

from ecd_epitope_scanner.assay_profiles import DEFAULT_ASSAY_PROFILE_ID
from ecd_epitope_scanner.cli import build_parser
from ecd_epitope_scanner.joint_design_v4.engine import _adapter
from ecd_epitope_scanner.pipeline import PipelineRequest, _pdb_only_context, _run_pipeline
from ecd_epitope_scanner.structure_analysis import best_chain_mapping, load_structure
from ecd_epitope_scanner.util import PipelineError


def _coordinates(tmp_path, chains):
    structure = Structure.Structure("binder")
    model = Model.Model(0)
    structure.add(model)
    for offset, (chain_id, start) in enumerate(chains):
        chain = Chain.Chain(chain_id)
        model.add(chain)
        for index in range(12):
            residue = Residue.Residue((" ", start + index, " "), "ALA", " ")
            residue.add(Atom.Atom("CA", np.array([index * 3.8, offset * 10., 0.]),
                                 0., 1., " ", " CA ", index + 1, element="C"))
            chain.add(residue)
    path = tmp_path / "binder.pdb"
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(path))
    return path


def test_ordinary_af3_archives_structure_and_confidence_before_normalization(tmp_path):
    import json
    from ecd_epitope_scanner.pipeline import _preserve_uploaded_inputs
    from ecd_epitope_scanner.util import sha256_file
    source = tmp_path / 'model.cif'
    source.write_bytes(b'data_original\n')
    confidence = tmp_path / 'full_data.json'
    confidence.write_text('{"atom_plddts":[90]}', encoding='utf-8')
    output = tmp_path / 'output'
    output.mkdir()
    request = PipelineRequest(input_type='af3', structure_path=str(source), confidence_path=str(confidence))
    identity = _preserve_uploaded_inputs(request, output)
    assert identity['files']['structure']['sha256'] == sha256_file(source)
    assert identity['files']['confidence']['sha256'] == sha256_file(confidence)
    assert (output / 'native_input.cif').read_bytes() == source.read_bytes()
    assert (output / 'source_confidence.json').read_bytes() == confidence.read_bytes()
    source.write_bytes(b'data_normalized\n')
    assert (output / 'native_input.cif').read_bytes() == b'data_original\n'
    assert json.loads((output / 'source_input_identity.json').read_text()) == identity


def test_input_archive_handles_no_local_file_and_already_archived_inputs(tmp_path):
    from ecd_epitope_scanner.pipeline import _preserve_uploaded_inputs
    assert _preserve_uploaded_inputs(PipelineRequest(input_type='uniprot'), tmp_path)['files'] == {}
    source = tmp_path / 'native_input.pdb'
    source.write_bytes(b'HEADER original\n')
    identity = _preserve_uploaded_inputs(PipelineRequest(input_type='structure', structure_path=str(source)), tmp_path)
    assert identity['files']['structure']['artifact'] == 'native_input.pdb'


def test_binder_context_preserves_explicit_non_one_numbered_chain(tmp_path):
    parsed = load_structure(_coordinates(tmp_path, [("A", 1), ("B", 78)]))
    request = PipelineRequest(input_type="structure", analysis_target="binder", chain="B",
                              assay_profile=DEFAULT_ASSAY_PROFILE_ID)
    protein, chain = _pdb_only_context(parsed, request, None)
    mapping = best_chain_mapping(parsed, protein, [chain], allowed_chains=[chain])
    assert chain == "B"
    assert mapping.uniprot_to_residue[1].id[1] == 78
    assert mapping.uniprot_to_residue[12].id[1] == 89
    assert protein.protein_name == "Coordinate-derived binder"
    assert protein.raw["canonical_numbering_available"] is False
    assert request.assay_profile is None
    assert request.assay_format == "unspecified"
    assert "K562" not in protein.topology


def test_binder_monomer_selects_its_only_chain(tmp_path):
    parsed = load_structure(_coordinates(tmp_path, [("Z", 94)]))
    protein, chain = _pdb_only_context(parsed, PipelineRequest(
        input_type="structure", analysis_target="binder"), None)
    assert chain == "Z"
    assert protein.ecd_sequence == "A" * 12


@pytest.mark.parametrize("sequence", [None, "A" * 12])
def test_binder_complex_requires_chain_even_with_supplied_sequence(tmp_path, sequence):
    path = _coordinates(tmp_path, [("A", 1), ("B", 78)])
    request = PipelineRequest(input_type="structure", analysis_target="binder",
        structure_path=str(path), sequence=sequence, ecd_start=1, ecd_end=12,
        output_dir=str(tmp_path / "run"))
    with pytest.raises(PipelineError, match="explicit binder chain"):
        _run_pipeline(request)


def test_binder_cli_and_disabled_epitope_adapter_override_environment(tmp_path, monkeypatch):
    args = build_parser().parse_args(["analyze", "--structure", "binder.pdb",
        "--analysis-target", "binder", "--chain", "A", "--output", "out"])
    assert args.analysis_target == "binder"
    monkeypatch.setenv("ECD_V4_DISCOTOPE_COMMAND", '["unwanted-worker"]')
    adapter = _adapter("discotope", {"enabled": False}, Path(tmp_path))
    assert adapter.command == []
    assert adapter.inventory().available is False
