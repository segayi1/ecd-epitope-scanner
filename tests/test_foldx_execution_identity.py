import json
from pathlib import Path

import pytest

from ecd_epitope_scanner import stability


@pytest.mark.parametrize('workers', [1, 2])
def test_each_candidate_has_a_canonical_private_invocation(monkeypatch, tmp_path, workers):
    calls = []
    a = [('K', 29, 'A', 'F'), ('E', 23, 'A', 'Y')]
    b = [('R', 15, 'A', 'A')]

    def native(binary, pdb, groups, directory, runs, repaired, seeded):
        assert len(groups) == 1  # A neighbour must never consume this candidate's sampler.
        assert groups[0] == sorted(groups[0], key=lambda x: (x[2], x[1], x[0], x[3]))
        calls.append(groups)
        return [float(sum(item[1] for item in groups[0]))]

    monkeypatch.setattr(stability, '_run_foldx_batch', native)
    assert stability.score_many_with_foldx('foldx', tmp_path/'input.pdb', [a,b], tmp_path,
        max_workers=workers) == [52.,15.]
    assert stability.score_many_with_foldx('foldx', tmp_path/'input.pdb', [b,list(reversed(a))], tmp_path,
        max_workers=workers) == [15.,52.]
    assert len(calls) == 4


def test_execution_receipt_freezes_inputs_command_and_outputs(monkeypatch, tmp_path):
    binary = tmp_path/'foldx.exe'; binary.write_bytes(b'fixture binary')
    source = tmp_path/'source.pdb'; source.write_bytes(b'fixed repaired input')
    (tmp_path/'rotabase.txt').write_bytes(b'fixed rotabase')

    def native(command, directory, timeout):
        assert '--timeSeedRotabase=false' in command
        request = json.loads((directory/'buildmodel_request.json').read_text())
        assert request['ordered_mutations'] == ['EA23Y,KA29F;']
        assert request['execution_layout'] == 'candidate_isolated_v1'
        assert (directory/'rotabase.txt').read_bytes() == b'fixed rotabase'
        (directory/'Dif_input.fxout').write_text('Pdb\ttotal energy\ninput_1_0.pdb\t2.1\n')
        return type('Result', (), {'returncode':0})()

    monkeypatch.setattr(stability, '_run_engine', native)
    assert stability.score_many_with_foldx(str(binary), source,
        [[('K',29,'A','F'),('E',23,'A','Y')]],tmp_path/'work',runs=1,
        already_repaired=True,max_workers=1) == [2.1]
    receipt = json.loads(next((tmp_path/'work').glob('engine_runs/*/buildmodel_request.json')).read_text())
    assert receipt['rotabase_sha256']
    assert receipt['output_sha256']['Dif_input.fxout']
    assert receipt['replay_mode'] == 'deterministic_native_reexecution'


def test_isolated_receipt_validation_survives_portable_export_and_detects_damage(monkeypatch, tmp_path):
    import numpy as np
    from zipfile import ZipFile
    from ecd_epitope_scanner.surface_perturbation.electrostatics import write_atoms
    from ecd_epitope_scanner.surface_perturbation.paired import scan_pairs, validate_pair
    from ecd_epitope_scanner.surface_perturbation.export import evidence_archive

    binary=tmp_path/'foldx.exe';binary.write_bytes(b'fixture binary')
    source=tmp_path/'source.pdb'
    wt={('A',1,'','CA'):(np.array([0.,0.,0.]),'ALA','C')}
    write_atoms(source,wt)
    (tmp_path/'rotabase.txt').write_bytes(b'fixed rotabase')

    def native(command,directory,timeout):
        write_atoms(directory/'WT_input_1.pdb',wt)
        write_atoms(directory/'input_1.pdb',{('A',1,'','CA'):(np.array([0.,0.,0.]),'VAL','C')})
        (directory/'PdbList_input.fxout').write_text('input_1.pdb\nWT_input_1.pdb\n')
        (directory/'Dif_input.fxout').write_text('Pdb\ttotal energy\ninput_1.pdb\t.35\n')
        return type('Result',(),{'returncode':0})()

    monkeypatch.setattr(stability,'_run_engine',native)
    output=tmp_path/'run';output.mkdir()
    stability.score_many_with_foldx(str(binary),source,[[('A',1,'A','V')]],output,runs=1,already_repaired=True)
    pairs,errors=scan_pairs(output);assert not errors and len(pairs)==1
    validate_pair(pairs[0])
    archive=evidence_archive(output,{'candidate':pairs},output/'null')
    review=tmp_path/'portable'
    with ZipFile(archive) as z:z.extractall(review)
    row=json.loads((review/'portable_pairs.json').read_text())['candidate'][0]
    for key in ('wt_path','mut_path','pdb_list_path','ddg_file_path','batch_path','input_path','mutant_file_path','buildmodel_request_path'):
        row[key]=str(review/row[key])
    validate_pair(row)
    rotabase=Path(row['batch_path'])/'rotabase.txt';rotabase.write_bytes(b'changed')
    with pytest.raises(ValueError,match='rotabase changed'):validate_pair(row)
    rotabase.write_bytes(b'fixed rotabase')
    Path(row['buildmodel_request_path']).write_text('{}')
    with pytest.raises(ValueError,match='execution receipt changed'):validate_pair(row)
