import hashlib
from pathlib import Path

import pytest

from ecd_epitope_scanner import stability
from ecd_epitope_scanner.integrations.foldx_adapter import FoldXAdapter
from ecd_epitope_scanner.integrations.workers import foldx_worker as worker
from ecd_epitope_scanner.joint_design_v4.models import PredictionRequest, PredictionResult


def test_shared_repair_runs_once_and_rejects_changed_cache(tmp_path, monkeypatch):
    source = tmp_path / 'input.pdb'
    source.write_bytes(b'WT')
    binary_hash = 'b' * 64
    workroot = tmp_path / 'shared' / 'work'
    calls = []

    def repair(binary, source, root, *, binary_sha256):
        calls.append(root)
        target = root.parent / 'foldx_repair_cache' / hashlib.sha256(source.read_bytes()).hexdigest()[:12] / binary_sha256[:12] / 'structure_Repair.pdb'
        target.parent.mkdir(parents=True)
        target.write_bytes(b'repaired WT')
        return target

    monkeypatch.setattr(stability, 'repair_with_foldx', repair)
    repaired, fresh = worker._site_first_repair('binary', source, workroot, binary_hash)
    assert fresh
    assert worker._site_first_repair('binary', source, workroot, binary_hash) == (repaired, False)
    assert len(calls) == 1
    repaired.write_bytes(b'corrupted')
    with pytest.raises(ValueError, match='identity changed'):
        worker._site_first_repair('binary', source, workroot, binary_hash)


@pytest.mark.parametrize('workers', [1, 2])
def test_native_candidate_exception_preserves_healthy_peer(tmp_path, monkeypatch, workers):
    def native(binary, source, groups, directory, runs, repaired, seeded):
        if groups[0][0][1] == 2:
            raise ValueError('one candidate crashed')
        return [.25]
    monkeypatch.setattr(stability, '_run_foldx_batch', native)
    groups = [[('A', 1, 'A', 'V')], [('A', 2, 'A', 'D')]]
    errors = []
    scores = stability.score_many_with_foldx('binary', tmp_path/'wt.pdb', groups, tmp_path,
        max_workers=workers, candidate_errors=errors)
    assert scores == [.25, None]
    assert len(errors) == 1 and errors[0]['candidate_index'] == 1


def test_native_empty_score_cannot_shift_peer_candidate_identity(tmp_path, monkeypatch):
    monkeypatch.setattr(stability, '_run_foldx_batch', lambda binary, source, groups, *args:
                        [] if groups[0][0][1] == 1 else [.25])
    groups = [[('A', 1, 'A', 'V')], [('A', 2, 'A', 'D')]]
    errors = []
    scores = stability.score_many_with_foldx('binary', tmp_path/'wt.pdb', groups, tmp_path,
        max_workers=2, candidate_errors=errors)
    assert scores == [None, .25]
    assert len(errors) == 1 and errors[0]['candidate_index'] == 0


@pytest.mark.parametrize('damaged', [None, 'run_identity', 'cross_batch', 'peer_scan'])
def test_worker_batch_retains_raw_pairs_for_successful_candidate(tmp_path, monkeypatch, damaged):
    from ecd_epitope_scanner.surface_perturbation import paired
    batch = tmp_path / 'native'
    batch.mkdir()
    (batch/'sampling_execution.json').write_text('{}')
    pairs = [{'mutations': ['A1V'], 'batch_path': str(batch), 'group': 0, 'run': index,
              'ddg_kcal_mol': .1} for index in range(3)]
    if damaged == 'run_identity':
        pairs[-1]['run'] = 4
    elif damaged == 'cross_batch':
        pairs[-1]['batch_path'] = str(tmp_path/'other_batch')
    monkeypatch.setattr(stability, 'score_many_with_foldx', lambda *_, **kwargs: [.1, None])
    monkeypatch.setattr(paired, 'scan_pairs', lambda _: (pairs, []))
    if damaged == 'peer_scan':
        def scan(root, *, listing_paths):
            if listing_paths[0].parent.name == 'bad':
                raise ValueError('Malformed peer PdbList')
            return pairs, []
        monkeypatch.setattr(paired, 'scan_pairs', scan)
        # The mocked engine returns immediately; list files are enough to exercise isolation.
        monkeypatch.setattr(worker, 'uuid4', lambda: type('Id', (), {'hex': 'fixture'})())
        for name in ['bad', 'good']:
            path = tmp_path/'fresh_fixture'/'engine_runs'/name/'PdbList_input.fxout'
            path.parent.mkdir(parents=True)
            path.write_text('fixture')
    monkeypatch.setattr(paired, 'validate_pair', lambda _: None)
    request = {'mutations': [['A1V'], ['A2D']], 'sequence': 'AA', 'chain_id': 'A',
        'parameters': {'number_of_runs': 3, 'sampling_slot_start': 0, 'sampling_protocol': 'time_rotabase_v3'}}
    identity = {'binary_sha256': 'b' * 64, 'structure_sha256': 'a' * 64}
    result = worker._site_first_batch(request, 'binary', tmp_path/'WT.pdb', tmp_path,
        tmp_path, tmp_path, identity, .2, True)
    assert result['status'] == 'partial'
    good, bad = result['records']
    if damaged in {'run_identity', 'cross_batch'}:
        assert good['status'] == 'failed'
    else:
        assert good['ddg'] == .1 and [p['sampling_slot'] for p in good['paired_models']] == [0, 1, 2]
    assert bad['status'] == 'failed'
    assert result['provenance']['attempted_slots'] == 6
    assert result['provenance']['repair_executed']
    assert (tmp_path/'paired_model_manifest.json').is_file()


@pytest.mark.parametrize(('aggregate', 'expected_status'), [
    (0.566, 'success'),     # Four-decimal aggregate of the saved raw FoldX median.
    (0.7, 'failed'),        # A genuinely mismatched worker aggregate stays rejected.
])
def test_worker_checks_raw_foldx_median_without_rounding(tmp_path, monkeypatch, aggregate, expected_status):
    from ecd_epitope_scanner.surface_perturbation import paired
    batch = tmp_path / 'native'
    batch.mkdir()
    (batch / 'sampling_execution.json').write_text('{}')
    energies = (0.565976, 0.615856, 0.287897)
    pairs = [{'mutations': ['K2A', 'E3W'], 'batch_path': str(batch), 'group': 1,
              'run': index, 'ddg_kcal_mol': value}
             for index, value in enumerate(energies)]
    monkeypatch.setattr(stability, 'score_many_with_foldx', lambda *_, **kwargs: [aggregate])
    monkeypatch.setattr(paired, 'scan_pairs', lambda *args, **kwargs: (pairs.copy(), []))
    monkeypatch.setattr(paired, 'validate_pair', lambda _: None)
    request = {'mutations': [['K2A', 'E3W']], 'sequence': 'AKE', 'chain_id': 'A',
        'parameters': {'number_of_runs': 3, 'sampling_slot_start': 0,
                       'sampling_protocol': 'time_rotabase_v3'}}
    result = worker._site_first_batch(request, 'binary', tmp_path/'WT.pdb', tmp_path,
        tmp_path, tmp_path, {'binary_sha256': 'b' * 64, 'structure_sha256': 'a' * 64}, .2, True)
    record = result['records'][0]
    assert record['status'] == expected_status, record
    if expected_status == 'success':
        assert record['ddg'] == 0.565976
    else:
        assert record['errors'][0]['message'] == 'FoldX aggregate differs from raw paired energy'


def test_adapter_validates_each_batch_candidate_without_erasing_peers(tmp_path):
    request = PredictionRequest(schema_version='1.0', request_id='test', structure_path=str(tmp_path/'WT.pdb'), chain_id='A', sequence='AA',
        mutations=[['A1V'], ['A2D']], parameters={'number_of_runs': 3, 'emit_paired_models': True,
            'sampling_protocol': 'time_rotabase_v3', 'sampling_cohort': 'ordinary_site_first_v2'})
    good = {'mutations': ['A1V'], 'ddg': .1, 'paired_models': [
        {'batch_path': 'batch', 'group': 0, 'run': i, 'sampling_slot': i, 'ddg_kcal_mol': .1}
        for i in range(3)]}
    bad = {'mutations': ['A2D'], 'status': 'failed', 'errors': [{'code': 'FAILED'}]}
    adapter = FoldXAdapter()
    result = adapter.checked_result(request, PredictionResult('partial', records=[bad, good]))
    assert result.status == 'partial' and len(result.records) == 2
    assert result.records[1] == good and result.errors
    good['ddg'] = 9
    rejected = adapter.checked_result(request, PredictionResult('success', records=[good]))
    assert all(row.get('status') == 'failed' for row in rejected.records)
