"""Exercise worker, paired structures and adapter without executing a model binary."""
import io
import json
import shutil
from types import SimpleNamespace

import numpy as np

from ecd_epitope_scanner import stability
from ecd_epitope_scanner.integrations.foldx_adapter import FoldXAdapter
from ecd_epitope_scanner.integrations.workers import foldx_worker
from ecd_epitope_scanner.joint_design_v4.models import PredictionRequest, PredictionResult
from ecd_epitope_scanner.surface_perturbation.electrostatics import write_atoms
from ecd_epitope_scanner.surface_perturbation.paired import validate_pair


def test_worker_produces_valid_isolated_pairs_and_reuses_wt_repair(tmp_path, monkeypatch, capsys):
    binary = tmp_path / 'foldx.exe'
    binary.write_bytes(b'test executable; never run')
    source = tmp_path / 'input.pdb'
    wt_atoms = {('A', 1, '', 'CA'): (np.array([0., 0., 0.]), 'ALA', 'C'),
                ('A', 2, '', 'CA'): (np.array([4., 0., 0.]), 'ALA', 'C')}
    write_atoms(source, wt_atoms)
    repair_calls = []

    def fake_engine(command, cwd, timeout):
        if '--command=RepairPDB' in command:
            repair_calls.append(str(cwd))
            shutil.copy2(cwd / 'structure.pdb', cwd / 'structure_Repair.pdb')
        else:
            notation = (cwd / 'individual_list.txt').read_text().strip().rstrip(';')
            position, amino_acid = int(notation[2:-1]), notation[-1]
            mutant = dict(wt_atoms)
            old = wt_atoms[('A', position, '', 'CA')]
            mutant[('A', position, '', 'CA')] = (old[0], {'V': 'VAL', 'D': 'ASP'}[amino_acid], old[2])
            runs = int(next(value.split('=')[1] for value in command if value.startswith('--numberOfRuns=')))
            pairs, energies = [], []
            for index in range(runs):
                wild, change = f'WT_input_1_{index}.pdb', f'input_1_{index}.pdb'
                write_atoms(cwd / wild, wt_atoms)
                write_atoms(cwd / change, mutant)
                pairs.extend([wild, change])
                energies.append(f'{change}\t0.25')
            (cwd / 'PdbList_input.fxout').write_text('\n'.join(pairs) + '\n')
            (cwd / 'Dif_input.fxout').write_text('Pdb\ttotal energy\n' + '\n'.join(energies) + '\n')
            (cwd / 'sampling_execution.json').write_text(json.dumps({'protocol': 'time_rotabase_v3'}))
        return SimpleNamespace(returncode=0, stdout='', stderr='')

    monkeypatch.setattr(stability, '_run_engine', fake_engine)
    base_parameters = {'binary': str(binary), 'sampling_protocol': 'time_rotabase_v3',
        'sampling_cohort': 'ordinary_site_first_v2', 'sampling_role': 'production',
        'repair_work_root': str(tmp_path / 'shared' / 'work'),
        'engine_work_root': str(tmp_path / 'short'), 'emit_paired_models': True, 'max_workers': 2}
    for index, (slot_start, runs) in enumerate([(0, 3), (3, 7)]):
        parameters = {**base_parameters, 'sampling_slot_start': slot_start, 'number_of_runs': runs}
        request = PredictionRequest(schema_version='1.0', request_id=str(index), structure_path=str(source),
            chain_id='A', sequence='AA', output_dir=str(tmp_path / f'batch{index}'),
            mutations=[['A1V'], ['A2D']], parameters=parameters)
        monkeypatch.setattr('sys.stdin', io.StringIO(json.dumps(request.as_dict())))
        assert foldx_worker.main() == 0
        payload = json.loads(capsys.readouterr().out)
        result = FoldXAdapter().checked_result(request, PredictionResult(**payload))
        assert result.status == 'success', result.errors
        assert len(result.records) == 2
        for row in result.records:
            assert [pair['sampling_slot'] for pair in row['paired_models']] == list(range(slot_start, slot_start + runs))
            for pair in row['paired_models']:
                validate_pair(pair)
    assert len(repair_calls) == 1
