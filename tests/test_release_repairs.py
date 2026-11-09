import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_public_assay_profiles_contain_no_private_construct_provenance():
    from ecd_epitope_scanner.assay_profiles import ASSAY_PROFILES, get_assay_profile

    for profile_id in ASSAY_PROFILES:
        profile = get_assay_profile(profile_id)
        for section in ('antigen_display', 'fab_display'):
            assert not any('plasmid' in key or 'sha256' in key or 'snapgene' in key
                           for key in profile[section])
        assert profile['antigen_display']['anchor_side'] == 'C_terminal'
        assert profile['design_policy']['membrane_proximal_penalty'] is False


def test_foldx_worker_uses_its_own_package(tmp_path):
    other = tmp_path / 'ecd_epitope_scanner'
    other.mkdir()
    (other / '__init__.py').write_text('', encoding='utf-8')
    worker = ROOT / 'integrations' / 'workers' / 'foldx_worker.py'
    script = (
        "import io, runpy, sys; "
        f"w=runpy.run_path({str(worker)!r}); "
        "sys.stdin=io.StringIO('{}'); w['main'](); "
        "import ecd_epitope_scanner.integrations.paired as p; print(p.__file__)"
    )
    env = dict(os.environ, PYTHONPATH=str(tmp_path))
    result = subprocess.run([sys.executable, '-c', script], cwd=tmp_path, env=env,
                            capture_output=True, text=True, check=True)
    assert str(ROOT / 'integrations' / 'paired.py') in result.stdout
    assert json.loads(result.stdout.splitlines()[0])['errors'][0]['message'] == (
        'FoldX received no mutation combinations'
    )


def test_missing_procaliper_is_unavailable(tmp_path):
    worker = ROOT / 'integrations' / 'workers' / 'procaliper_worker.py'
    # Isolate the specific missing optional package even if installed locally.
    script = (
        "import builtins,io,runpy,sys; original=builtins.__import__; "
        "exec(\"def blocked(name, *args, **kwargs):\\n"
        " if name.startswith('procaliper'): raise ModuleNotFoundError(name='procaliper')\\n"
        " return original(name, *args, **kwargs)\"); "
        "builtins.__import__=blocked; sys.stdin=io.StringIO('{}'); "
        f"runpy.run_path({str(worker)!r}, run_name='__main__')"
    )
    result = subprocess.run([sys.executable, '-c', script], cwd=tmp_path,
                            capture_output=True, text=True, check=True)
    payload = json.loads(result.stdout)
    assert payload['status'] == 'failed'
    assert payload['errors'][0]['code'] == 'TOOL_NOT_CONFIGURED'


def test_missing_confidence_does_not_use_bfactor(tmp_path, monkeypatch):
    from ecd_epitope_scanner import structure_analysis as analysis
    from ecd_epitope_scanner.models import ProteinContext

    parsed = analysis.load_structure(ROOT / 'examples' / 'pdl1_run' / 'structure.pdb')
    chain = next(parsed.model.get_chains())
    residue = next(chain.get_residues())
    for atom in residue:
        atom.bfactor = 99.0
    aa = analysis._one_letter(residue.resname)
    protein = ProteinContext('test', 'test', 'test', aa, 1, 1, 1, 'test')
    mapping = analysis.ChainMapping(chain.id, [residue], {0: 1}, {1: residue}, 1.0, 1.0, 0, 0)
    monkeypatch.setattr(analysis, 'calculate_sasa', lambda *args: ({}, 'test'))
    monkeypatch.setattr(analysis, 'calculate_sidechain_exposure', lambda *args: {})
    records, _ = analysis.analyze_residues(parsed, mapping, protein, tmp_path,
                                            predicted=True, confidence=None)
    assert records[0].plddt is None
    assert 'low-pLDDT' in records[0].review_reasons


def test_unmapped_confidence_cannot_be_assigned_by_atom_count():
    from ecd_epitope_scanner import structure_analysis as analysis

    parsed = analysis.load_structure(ROOT / 'examples' / 'pdl1_run' / 'structure.pdb')
    values = [99.0] * len(list(parsed.model.get_atoms()))
    assert analysis._af3_atom_plddt(parsed, {'atom_plddts': values}) == {}


def test_confidence_cif_validates_atom_mapping(tmp_path):
    import pytest
    from Bio.PDB import MMCIFIO
    from ecd_epitope_scanner import structure_analysis as analysis
    from ecd_epitope_scanner.util import PipelineError

    source = analysis.load_structure(ROOT / 'examples' / 'pdl1_run' / 'structure.pdb')
    cif = tmp_path / 'model.cif'
    writer = MMCIFIO()
    writer.set_structure(source.structure)
    writer.save(str(cif))
    parsed = analysis.load_structure(cif)
    values = [float(atom.bfactor) for atom in parsed.model.get_atoms()]
    mapped = analysis._af3_atom_plddt(parsed, {'atom_plddts': values})
    assert mapped
    values[0] = 101.0
    with pytest.raises(PipelineError):
        analysis._af3_atom_plddt(parsed, {'atom_plddts': values})


def test_preserved_input_identity_uses_artifact_and_digest(tmp_path):
    from ecd_epitope_scanner.pipeline import PipelineRequest, _preserve_uploaded_inputs
    from ecd_epitope_scanner.util import sha256_file

    source = tmp_path / 'private-owner' / 'model.pdb'
    source.parent.mkdir()
    source.write_text('ATOM\n', encoding='utf-8')
    output = tmp_path / 'output'
    output.mkdir()
    identity = _preserve_uploaded_inputs(PipelineRequest(input_type='structure',
                                         structure_path=str(source)), output)
    assert identity['files']['structure'] == {
        'artifact': 'native_input.pdb', 'sha256': sha256_file(source)
    }
    assert (output / 'native_input.pdb').read_bytes() == source.read_bytes()
    assert 'private-owner' not in (output / 'source_input_identity.json').read_text()


def test_recorded_worker_commands_omit_local_paths():
    from ecd_epitope_scanner.joint_design_v4.reporting import _portable_tool_commands

    commands = {'spurs': [('C:' + '\\Users\\private-owner\\python.exe'),
                          ('C:' + '\\Users\\private-owner\\wsl_bridge.py'), '--python',
                          ('/' + 'home/' + 'private-linux/env/bin/python'),
                          ('--worker=/' + 'home/' + 'private-linux/tool.py'), '--device', 'cuda']}
    safe = _portable_tool_commands(commands)['spurs']
    assert safe == ['<local>/python.exe', '<local>/wsl_bridge.py', '--python',
                    '<local>/python', '--worker=<local-path>', '--device', 'cuda']
    assert 'private' not in str(safe)
    assert commands['spurs'][0].startswith('C:')


def test_foldx_runs_original_command_but_records_no_private_path(tmp_path, monkeypatch):
    import subprocess
    from ecd_epitope_scanner import stability

    command = [('C:' + '\\Users\\private-owner\\foldx.exe'), '--command=BuildModel']
    seen = []
    def fake_run(args, **kwargs):
        seen.append(args)
        return subprocess.CompletedProcess(args, 0, stdout='', stderr='')
    monkeypatch.setattr(stability.subprocess, 'run', fake_run)
    assert stability._run_engine(command, tmp_path, 10).returncode == 0
    assert seen == [command]
    receipt = json.loads((tmp_path / 'engine_execution.json').read_text())
    assert receipt['command_paths_redacted'] is True
    assert 'private-owner' not in str(receipt)
    assert receipt['command'][1] == '--command=BuildModel'


def test_portable_projection_redacts_drive_home_and_mnt_paths(tmp_path):
    from ecd_epitope_scanner.util import portable_projection, sanitize_text_file

    value = {
        'windows': ('C:' + '\\Users\\private-owner\\run\\analysis.json'),
        'home': ('/' + 'home/' + 'private-linux/work'),
        'mnt': ('/' + 'mnt/e/private-owner/data'),
        'url': 'https://example.org/not-a-local-path',
    }
    portable = portable_projection(value)
    text = json.dumps(portable)
    assert 'private-owner' not in text and 'private-linux' not in text
    assert portable['url'] == value['url']
    artifact = tmp_path / 'analysis.json'
    artifact.write_text(json.dumps(value), encoding='utf-8')
    assert sanitize_text_file(artifact)
    assert 'private-owner' not in artifact.read_text(encoding='utf-8')
    escaped = json.dumps(value)
    from ecd_epitope_scanner.util import portable_text
    assert 'private-owner' not in portable_text(escaped)


def test_request_export_projection_keeps_runtime_path_out_of_serialized_identity(tmp_path):
    from ecd_epitope_scanner.pipeline import PipelineRequest

    private = tmp_path / 'private-owner' / 'run'
    request = PipelineRequest(input_type='structure', output_dir=str(private),
                              structure_path=str(private / 'model.pdb'),
                              metadata={'joint_v4_ml': {'custom': {'command': [
                                  ('C:' + '\\Users\\private-owner\\python.exe'),
                                  str(private / 'foldx_worker.py'),
                              ]}}})
    assert request.as_dict()['output_dir'] == str(private)
    exported = request.export_dict()
    assert 'private-owner' not in json.dumps(exported)
    assert exported['output_dir'] == '<local-path>'
