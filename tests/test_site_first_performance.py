from copy import deepcopy

import pytest

from ecd_epitope_scanner.joint_design_v4.site_first_performance import compare_performance, peak_process_bytes


def test_runtime_gate_requires_matched_inputs_cache_and_half_time():
    baseline = {'protocol': {'name': 'legacy'}, 'input': {'structure_path': 'hash'}, 'hardware': {'cpu': 'host'},
        'cache_condition': 'warm', 'score_cache_initially_empty': True, 'models_disabled': False, 'total_seconds': 100}
    optimized = deepcopy(baseline)
    optimized['protocol'] = {'name': 'ordinary_site_first_v2'}
    optimized['total_seconds'] = 50
    assert compare_performance(baseline, optimized)['runtime_gate'] == 'PROVISIONAL_PASS'
    optimized['total_seconds'] = 51
    assert compare_performance(baseline, optimized)['runtime_gate'] == 'NOT_MET'
    for field, value in [('cache_condition', 'cold'), ('input', {}), ('hardware', {}),
                         ('score_cache_initially_empty', False), ('models_disabled', True),
                         ('total_seconds', float('nan'))]:
        changed = {**optimized, field: value}
        with pytest.raises(ValueError):
            compare_performance(baseline, changed)
    assert peak_process_bytes() > 0


def test_new_protocol_keeps_three_or_ten_card_contract():
    from ecd_epitope_scanner.joint_design_v4.impact_completeness import _repeat_protocol
    assert _repeat_protocol({'execution_protocol': {'name': 'ordinary_site_first_v2'}}) == (
        'ordinary_site_first_v2', [list(range(3)), list(range(10))])


def test_post_export_receipt_does_not_change_scientific_archive(tmp_path):
    from types import SimpleNamespace
    from time import perf_counter
    from ecd_epitope_scanner.joint_design_v4.site_first_performance import write_performance
    archive = tmp_path/'ecd_epitope_design.zip'
    archive.write_bytes(b'frozen scientific archive')
    source = tmp_path/'source.pdb'
    source.write_bytes(b'WT')
    request = SimpleNamespace(structure_path=str(source), confidence_path=None, chain='A',
        input_type='structure', analysis_target='binder', metadata={'benchmark_cache_condition': 'warm'})
    receipt = write_performance(tmp_path, request, {'execution_protocol': {'name': 'ordinary_site_first_v2'}},
        perf_counter()-2, perf_counter()-1, True)
    assert receipt['total_seconds'] >= 2 and receipt['export_seconds'] >= 1
    assert archive.read_bytes() == b'frozen scientific archive'
    assert receipt['scientific_acceptance'] == 'NOT_ASSESSED_BY_RUNTIME_RECEIPT'
