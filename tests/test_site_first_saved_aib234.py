"""Replay the saved AIB234 structure through the pipeline with fake model scores."""
import hashlib
import json
from pathlib import Path

import pytest

from ecd_epitope_scanner.joint_design_v4.models import PredictionResult


@pytest.mark.parametrize("ten_slot_confirmation", [False, True])
def test_saved_aib234_pipeline_replay_without_model_calls(tmp_path, monkeypatch, ten_slot_confirmation):
    from ecd_epitope_scanner import pipeline
    from ecd_epitope_scanner.joint_design_v4 import engine

    root = Path(__file__).resolve().parents[1]
    saved = root / 'workspace_data/2b522bedd0f5'
    structure = saved / 'uploads/structure.pdb'
    request_path = saved / 'request.json'
    if not structure.is_file() or not request_path.is_file():
        pytest.skip('Saved AIB234 source input is unavailable')
    assert hashlib.sha256(structure.read_bytes()).hexdigest() == '8ff23af24604861cec56365dc835068f97fa714222b8313a71291da2dc7851e4'
    calls = []

    def predict(adapter, structure, chain, sequence, output, structure_type, parameters, supplied,
                key=None, mutations=None):
        key = key or adapter.name
        calls.append(key)
        if key == 'spurs':
            return PredictionResult('success', records=[{'position': index, 'wt': aa, 'scores': {
                mutant: .1 for mutant in 'ACDEFGHIKLMNPQRSTVWY' if mutant != aa}}
                for index, aa in enumerate(sequence, 1)])
        if key in {'spurs_multi', 'foldx'}:
            rows = [{'mutations': group, 'ddg': .1} for group in mutations]
            if key == 'foldx':
                for index, row in enumerate(rows):
                    start = parameters['sampling_slot_start']
                    row['paired_models'] = [{'sampling_slot': slot,
                        'ddg_kcal_mol': 2.5 if ten_slot_confirmation and slot == 0 else .1,
                        'batch_path': str(output / f'candidate{index}'), 'group': 1, 'run': slot}
                        for slot in range(start, start + parameters['number_of_runs'])]
            return PredictionResult('success', records=rows)
        if key == 'procaliper':
            return PredictionResult('success', records=[{'position': index, 'sasa': 100., 'rsasa': .5}
                for index in range(1, len(sequence) + 1)])
        if key == 'discotope':
            return PredictionResult('success', records=[{'position': index, 'score': index / len(sequence)}
                for index in range(1, len(sequence) + 1)])
        raise AssertionError(f'Unexpected model: {key}')

    monkeypatch.setattr(engine, '_prediction', predict)
    data = json.loads(request_path.read_text(encoding='utf-8'))
    data.update(structure_path=str(structure), output_dir=str(tmp_path),
        metadata={'joint_v4_ml': {'execution_protocol': {'name': 'ordinary_site_first_v2'}}})
    result = pipeline.run_pipeline(pipeline.PipelineRequest(**data))
    assert result['execution_protocol']['name'] == 'ordinary_site_first_v2'
    assert result['execution_protocol']['display_prefilter'] == 'known_review_v1'
    assert result['execution_protocol']['repeat_confirmation_policy'] == 'ten_slot_median_v1'
    assert result['source']['structure_sha256']
    assert (tmp_path / 'ecd_epitope_design.zip').is_file()
    assert not (tmp_path / 'ordinary_candidate_index.sqlite3').exists()
    receipt = json.loads((tmp_path / 'card_evidence_completeness.json').read_text(encoding='utf-8'))
    assert receipt['status'] == 'PASS' and receipt['card_count'] > 0
    assert {'spurs', 'spurs_multi', 'foldx'} <= set(calls)
    panel = json.loads((tmp_path / 'selected_panel_v4.json').read_text(encoding='utf-8'))
    candidates = {row['candidate_id']: row for row in json.loads((tmp_path / 'panel_candidates.json').read_text(encoding='utf-8'))}
    if ten_slot_confirmation:
        from ecd_epitope_scanner.joint_design_v4.site_first_display import confirmed_repeat_pass
        for candidate_id in panel['selected_candidate_ids']:
            audit = candidates[candidate_id]['impact_evidence']['stability_repeat_audit']
            assert confirmed_repeat_pass(audit)
            assert audit['sampling_decision_sensitive'] is True
        recommendations = json.loads((tmp_path / 'sector_recommendations_v4.json').read_text(encoding='utf-8'))
        assert recommendations and all(confirmed_repeat_pass(card['stability_repeat_diagnostic'])
            for main in recommendations for card in [main, *main.get('alternatives', [])])
    assert panel['hypothesis_signatures'] and panel['selected_candidate_ids']
    for hypothesis, signature in panel['hypothesis_signatures'].items():
        for index, candidate_id in enumerate(panel['selected_candidate_ids']):
            candidate = candidates[candidate_id]
            responses = candidate.get('hypothesis_perturbation') or candidate.get('hypothesis_response') or {}
            assert hypothesis in responses
            assert signature[index] == pytest.approx(responses[hypothesis], abs=1e-6)
    import importlib.util
    spec = importlib.util.spec_from_file_location('ordinary_rerun', root / 'scripts/run_ordinary_adaptive_rerun.py')
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    assert runner._site_first_display_errors(tmp_path) == []
