import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from build_site_first_assessment import (audit_formal_display, audit_panel_signatures,
    compare_lanes, evaluate, scientific_status, selected, table_escape,
    unique_formal_combinations)


def test_performance_target_requires_matching_models_and_actual_cold_caches():
    gpu = {'model_files_sha256': {'weight': 'sha'}, 'actual_device': 'cuda', 'gpu_name': 'gpu'}
    cpu = {'version': '1', 'biopython_version': '1', 'algorithm': 'SASA', 'probe_radius_a': 1.4, 'sphere_points': 100}
    identity = {'worker_model_identities': {'spurs_single': gpu, 'spurs_multi': gpu, 'procaliper': cpu}}
    new = {'status': 'PASS', 'input_sha256': 'input', 'elapsed_seconds': 49,
           'cache_evidence': {'before_run': {'cold_verified': True}},
           'source_digest_checkpoint': 's', 'source_digest_run': 's', 'source_digest_end': 's'}
    old = {'status': 'COMPLETE', 'input_sha256': 'input', 'elapsed_pipeline_seconds': 100,
           'normalized_structure': {'sha256': 'structure'}, 'model_identity': {'foldx_binary': {'sha256': 'fx'}},
           'cold_cache_evidence': {'score_cache_entries_at_start': 0, 'all_possible_repair_cache_paths_absent_before_run': True}}
    ledger = {'spurs': {'batches': [{'provenance': gpu}]}, 'foldx_results': [{'paired_models': [{'binary_sha256': 'fx'}]}]}
    manifest = {'parameters': {'joint_v4_ml': {'discotope': {'enabled': False}}}}
    args = (new, old, identity, identity, ledger, manifest, manifest, {'status': 'PASS'}, 'structure')
    assert evaluate(*args)['performance_target'] == 'NOT_ASSESSED'
    assert evaluate(*args, os_page_cache_matched=True)['performance_target'] == 'PASS'
    mismatched = {**identity, 'worker_model_identities': {**identity['worker_model_identities'], 'spurs_single': {**gpu, 'gpu_name': 'other'}}}
    assert evaluate(new, old, mismatched, identity, ledger, manifest, manifest, {'status': 'PASS'}, 'structure')['performance_target'] == 'NOT_ASSESSED'
    assert evaluate({**new, 'status': 'NEEDS_REVIEW'}, old, identity, identity, ledger, manifest, manifest, {'status': 'PASS'}, 'structure')['performance_target'] == 'NOT_ASSESSED'


def test_scientific_acceptance_requires_post_validation_display_signatures_and_identity():
    checks = {
        'new_acceptance': True, 'same_input_sha256': True, 'same_normalized_structure': True,
        'spurs_single_weights_and_gpu': True, 'spurs_multi_weights_and_gpu': True,
        'procaliper_identity': True, 'discotope_disabled_in_both_binder_runs': True,
        'foldx_binary_identity': True, 'legacy_raw_foldx_binary_identity': True,
        'original_input_unchanged': True, 'new_source_stable': True,
    }
    new = {'validation': {'status': 'PASS', 'errors': []}}
    signatures = {'status': 'PASS', 'mismatches': [], 'unknown': []}
    status, failed, unknown = scientific_status(new, dict(checks), [], signatures, [])
    assert (status, failed, unknown) == ('PASS', [], [])

    new_review = {'validation': {'status': 'NEEDS_REVIEW', 'errors': ['failed gate']}}
    status, failed, _ = scientific_status(new_review, dict(checks), ['failed gate'], signatures, [])
    assert status == 'NEEDS_REVIEW'
    assert {'post_validation_pass', 'post_validation_errors_absent'} <= set(failed)

    status, failed, _ = scientific_status(new, {**checks, 'same_input_sha256': False}, [], signatures, [])
    assert status == 'NEEDS_REVIEW' and 'same_input_sha256' in failed
    status, failed, _ = scientific_status(new, dict(checks), [], signatures, ['display failure'])
    assert status == 'NEEDS_REVIEW' and 'formal_display_eligibility' in failed


def test_missing_signature_or_validation_evidence_is_unknown():
    audit = audit_panel_signatures({}, [])
    assert audit['status'] == 'UNKNOWN'
    assert audit['unknown']
    checks = {
        'new_acceptance': True, 'same_input_sha256': True, 'same_normalized_structure': True,
        'spurs_single_weights_and_gpu': True, 'spurs_multi_weights_and_gpu': True,
        'procaliper_identity': True, 'discotope_disabled_in_both_binder_runs': True,
        'foldx_binary_identity': True, 'legacy_raw_foldx_binary_identity': True,
        'original_input_unchanged': True, 'new_source_stable': True,
    }
    status, failed, unknown = scientific_status({}, checks, None, audit, None)
    assert status == 'NOT_ASSESSED'
    assert not failed
    assert {'post_validation_pass', 'post_validation_errors_absent',
            'formal_display_eligibility', 'panel_signature_consistency'} <= set(unknown)


def test_signature_audit_does_not_replace_missing_values_with_zero():
    panel = {'selected_candidate_ids': ['C1'], 'hypothesis_signatures': {'H1': [0.0, 1.0]}}
    candidates = [{'candidate_id': 'C1', 'hypothesis_perturbation': {'H1': 0.0}}]
    audit = audit_panel_signatures(panel, candidates)
    assert audit['status'] == 'UNKNOWN'
    assert audit['compared_values'] == 0
    assert audit['unknown']

    panel['hypothesis_signatures'] = {'H1': [0.0]}
    candidates[0]['hypothesis_response'] = {'H1': 0.5}
    audit = audit_panel_signatures(panel, candidates)
    assert audit['status'] == 'NEEDS_REVIEW'
    assert len(audit['mismatches']) == 1


def test_lane_comparison_reports_gaps_and_compares_shared_lanes_only():
    old = {('VS001', 'mechanism_pareto'): {'status': 'FORMAL_ELIGIBLE', 'mutations': [['A1V']]}}
    new = {
        ('VS001', 'mechanism_pareto'): {'status': 'FORMAL_ELIGIBLE', 'mutations': [['A1V']]},
        ('VS002', 'sector_coverage'): {'status': 'FORMAL_ELIGIBLE', 'mutations': [['B2L']]},
    }
    rows = compare_lanes(old, new)
    assert len(rows) == 2
    shared = next(row for row in rows if row['region'] == 'VS001')
    gap = next(row for row in rows if row['region'] == 'VS002')
    assert shared['comparable'] and shared['identical_mutations'] is True
    assert gap['comparable'] is False and gap['comparison_status'] == 'NOT_REPORTED_ON_ONE_SIDE'
    assert gap['old_status'] == 'NOT_REPORTED'


def test_lane_comparison_all_pass_and_historical_display_field_unknown(tmp_path):
    lanes = {
        ('VS001', 'mechanism_pareto'): {'status': 'FORMAL_ELIGIBLE', 'mutations': [['A1V']]},
        ('VS002', 'sector_coverage'): {'status': 'FORMAL_ELIGIBLE', 'mutations': [['B2L']]},
    }
    rows = compare_lanes(lanes, lanes)
    assert len(rows) == 2 and all(row['comparable'] for row in rows)
    assert all(row['identical_mutations'] and row['identical_sites'] for row in rows)

    def write(name, value):
        (tmp_path / name).write_text(__import__('json').dumps(value), encoding='utf-8')
    write('selected_panel_v4.json', {'selected_candidate_ids': ['C1']})
    write('panel_candidates.json', [{'candidate_id': 'C1', 'formal_eligible': True}])
    write('sector_recommendations_v4.json', [{'candidate_id': 'C1', 'formal_eligible': True,
                                               'principles': ['mechanism_pareto']}])
    audit = audit_formal_display(tmp_path, ['site-first display qualification unknown: C1'])
    assert audit['status'] == 'UNKNOWN'
    assert audit['unknown'] == ['C1: display_eligible field missing']
    checks = {
        'new_acceptance': True, 'same_input_sha256': True, 'same_normalized_structure': True,
        'spurs_single_weights_and_gpu': True, 'spurs_multi_weights_and_gpu': True,
        'procaliper_identity': True, 'discotope_disabled_in_both_binder_runs': True,
        'foldx_binary_identity': True, 'legacy_raw_foldx_binary_identity': True,
        'original_input_unchanged': True, 'new_source_stable': True,
    }
    status, failed, unknown = scientific_status({'validation': {'status': 'PASS', 'errors': []}},
        checks, [], {'status': 'PASS'}, audit)
    assert status == 'NOT_ASSESSED' and not failed
    assert unknown == ['formal_display_eligibility']

    candidate = {'candidate_id': 'C1', 'formal_eligible': True, 'display_eligible': True}
    write('panel_candidates.json', [candidate])
    audit = audit_formal_display(tmp_path, [])
    assert audit['status'] == 'PASS'
    assert not audit['errors'] and not audit['unknown']

    candidate['display_eligible'] = False
    write('panel_candidates.json', [candidate])
    audit = audit_formal_display(tmp_path, [])
    assert audit['status'] == 'NEEDS_REVIEW'


def test_lane_comparison_marks_explicit_not_applicable():
    rows = compare_lanes({('VS001', 'epitope_support'): {'status': 'NOT_REPORTED', 'mutations': []}},
                         {('VS001', 'epitope_support'): {'status': 'NOT_APPLICABLE', 'mutations': []}})
    assert rows[0]['comparison_status'] == 'NOT_APPLICABLE'


def test_actual_recommendation_shape_yields_notation_combinations_and_escaped_table_cells(tmp_path):
    recommendations = [{
        'sector_id': 'VS001', 'candidate_id': 'VS001|C1', 'formal_eligible': True,
        'principles': ['mechanism_pareto'], 'mutation_notations': ['A12V', 'B33K'],
        'alternatives': [{
            'sector_id': 'VS001', 'candidate_id': 'VS001|C2', 'formal_eligible': True,
            'principles': ['sector_coverage'], 'mutation_notations': ['C8Y'],
        }],
        'unavailable_principles': [{'principle': 'epitope_support', 'state': 'not_applicable'}],
    }]
    (tmp_path / 'sector_recommendations_v4.json').write_text(
        __import__('json').dumps(recommendations), encoding='utf-8')
    lanes = selected(tmp_path)
    assert lanes[('VS001', 'mechanism_pareto')]['mutations'] == [['A12V', 'B33K']]
    assert lanes[('VS001', 'sector_coverage')]['mutations'] == [['C8Y']]
    assert lanes[('VS001', 'epitope_support')]['status'] == 'NOT_APPLICABLE'
    assert unique_formal_combinations(lanes) == {('A12V', 'B33K'), ('C8Y',)}
    assert table_escape('VS001|C1') == r'VS001\|C1'
