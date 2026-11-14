"""Display size preferences remain separate from hard geometry and science gates."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from ecd_epitope_scanner.config import JOINT_V4_ML_PARAMETERS
from ecd_epitope_scanner.joint_design_v4 import sectors as module
from ecd_epitope_scanner.joint_design_v4.models import EpitopeEvidence, ResidueKey, SurfaceEvidence
from ecd_epitope_scanner.models import ResidueRecord


def saved_inputs(name):
    data = json.loads((Path(__file__).parent / 'fixtures' / name).read_text(encoding='utf-8'))
    records = {row['uniprot_pos']: ResidueRecord(**row) for row in data['records']}
    graph = {p: {} for p in records}
    for row in data['edges']:
        left, right, distance = (row['source'], row['target'], row['distance_a']) if isinstance(row, dict) else row
        graph[left][right] = graph[right][left] = distance
    surface = ([SurfaceEvidence(**(row | {'residue': ResidueKey(**row['residue'])})) for row in data['surface']]
               if 'surface' in data else [SurfaceEvidence(
                   ResidueKey('0', 'A', row['uniprot_pos'], uniprot_pos=row['uniprot_pos']),
                   mapped_surface_area_a2=row['surface_area_a2']) for row in data['nodes']])
    epitope = [EpitopeEvidence(**(row | {'residue': ResidueKey(**row['residue'])})) for row in data['epitope']]
    return surface, epitope, graph, records


@pytest.fixture
def pdl1():
    return saved_inputs('pdl1_balanced_display.json')


def assert_partition(result, surface, graph, records, geo=30, ca=22, induced=True):
    assigned = [p for row in result for p in row.residue_positions]
    assert len(assigned) == len(set(assigned)) == len(graph) and set(assigned) == set(graph)
    areas = {row.residue.uniprot_pos: row.mapped_surface_area_a2 for row in surface}
    for row in result:
        group = set(row.residue_positions)
        assert len(group) >= 4 and sum(areas[p] for p in group) >= 50
        assert set(module._distance(graph, min(group), group)) == group
        for p in group:
            distances = module._distance(graph, p, group if induced else set(graph))
            for q in group:
                assert distances[q] <= geo
                assert np.linalg.norm(np.asarray(records[p].ca) - records[q].ca) <= ca


def line_inputs(count):
    records = {p: ResidueRecord(p, 'A', 'A', p, ca=(float(p), 0, 0)) for p in range(1, count + 1)}
    graph = {p: {} for p in records}
    for p in range(1, count):
        graph[p][p + 1] = graph[p + 1][p] = 1.0
    surface = [SurfaceEvidence(ResidueKey('0', 'A', p, uniprot_pos=p), mapped_surface_area_a2=80)
               for p in graph]
    return surface, [], graph, records


def test_new_defaults_only_change_display_scale_and_preferences():
    params = JOINT_V4_ML_PARAMETERS
    assert params['sectoring'] == {
        'policy': 'balanced_connected_v4', 'target_sector_residues': 12,
        'preferred_minimum_sector_residues': 10, 'preferred_maximum_sector_residues': 15,
        'maximum_geodesic_diameter_a': 30.0, 'maximum_euclidean_diameter_a': 22.0,
        'minimum_sector_area_a2': 50.0, 'minimum_sector_residues': 4}
    assert params['candidate_search']['maximum_geodesic_span_a'] == 24.0


def test_saved_pdl1_meets_induced_diameter(pdl1):
    surface, epitope, graph, records = pdl1
    audit = {}
    result = module.build_epitope_sectors(*pdl1, diagnostics=audit)
    assert len(graph) == 103
    assert_partition(result, surface, graph, records)
    assert audit['version'] == 'balanced_connected_v4'
    assert audit['optimization_status'] == 'LOCAL_SEARCH_COMPLETE'
    assert audit['search_states'] <= 100_000
    assert len(audit['initializations']) == 4
    assert audit['diameter_basis'] == 'induced_eligible_residue_graph_shortest_paths'
    assert len(audit['region_metrics']) == len(result)
    for row, metrics in zip(result, audit['region_metrics']):
        assert row.sector_id == metrics['sector_id']
        assert metrics['full_graph_geodesic_diameter_a'] <= 30
        assert metrics['induced_graph_geodesic_diameter_a'] >= metrics['full_graph_geodesic_diameter_a']
        assert metrics['induced_graph_geodesic_diameter_a'] <= 30
        assert row.geodesic_diameter_a == round(metrics['induced_graph_geodesic_diameter_a'], 3)
        assert metrics['ca_diameter_a'] <= 22


def test_saved_pdl1_order_is_irrelevant_including_audit(pdl1):
    surface, epitope, graph, records = pdl1
    first, second = {}, {}
    result = module.build_epitope_sectors(*pdl1, diagnostics=first)
    reversed_graph = {p: dict(reversed(list(neighbors.items()))) for p, neighbors in reversed(list(graph.items()))}
    other = module.build_epitope_sectors(list(reversed(surface)), list(reversed(epitope)), reversed_graph,
              dict(reversed(list(records.items()))), diagnostics=second)
    assert other == result
    assert second == first


def test_ctla4_balanced_defaults_keep_original_mapping_and_hard_rules():
    inputs = saved_inputs('ctla4_small_partition.json')
    snapshot = {p: (row.chain_id, row.auth_seq_id, tuple(row.ca), list(row.hard_exclusions))
                for p, row in inputs[3].items()}
    audit = {}
    result = module.build_epitope_sectors(*inputs, diagnostics=audit)
    assert len(inputs[2]) == 71
    assert_partition(result, inputs[0], inputs[2], inputs[3])
    assert audit['version'] == 'balanced_connected_v4'
    assert snapshot == {p: (row.chain_id, row.auth_seq_id, tuple(row.ca), list(row.hard_exclusions))
                        for p, row in inputs[3].items()}


def test_explicit_old_scale_is_hard_and_does_not_claim_full_size_success(pdl1):
    audit = {}
    result = module.build_epitope_sectors(*pdl1, parameters={
        'maximum_geodesic_diameter_a': 26, 'maximum_euclidean_diameter_a': 18}, diagnostics=audit)
    assert_partition(result, pdl1[0], pdl1[2], pdl1[3], 26, 18)
    assert len(result) >= 13
    assert audit['version'] == 'balanced_connected_v4'
    assert audit['size_exceptions']
    assert all(row['reason'] in {'LOCAL_BALANCE_LIMIT', 'SEARCH_BUDGET_EXHAUSTED'}
               for row in audit['size_exceptions'])


def test_legacy_policy_retains_saved_pdl1_and_distinguishes_path_metrics(pdl1):
    audit = {}
    result = module.build_epitope_sectors(*pdl1, max_sectors=7, parameters={'policy': 'connected_bounded_v3'}, diagnostics=audit)
    assert [len(row.residue_positions) for row in result] == [7, 8, 8, 4, 13, 8, 4, 4, 7, 13, 6, 9, 8, 4]
    assert_partition(result, pdl1[0], pdl1[2], pdl1[3], 26, 18, induced=False)
    assert audit['version'] == 'connected_bounded_v3'
    metrics = next(row for row in audit['region_metrics'] if row['sector_id'] == 'VS010')
    assert metrics['full_graph_geodesic_diameter_a'] < 26 < metrics['induced_graph_geodesic_diameter_a']


def test_saved_ctla4_diagnostic_replay_obeys_region_paths_and_is_deterministic():
    name = 'ctla4_1d812aa2f9be_partition.json'
    inputs = saved_inputs(name)
    params = json.loads((Path(__file__).parent / 'fixtures' / name).read_text(encoding='utf-8'))['parameters']
    audit, repeat = {}, {}
    result = module.build_epitope_sectors(*inputs, max_sectors=12, parameters=params, diagnostics=audit)
    other = module.build_epitope_sectors(*inputs, max_sectors=12, parameters=params, diagnostics=repeat)
    assert len(inputs[2]) == 74
    assert_partition(result, inputs[0], inputs[2], inputs[3])
    assert result == other and audit == repeat
    assert all(row.geodesic_diameter_a == round(metric['induced_graph_geodesic_diameter_a'], 3)
               for row, metric in zip(result, audit['region_metrics']))


def test_region_path_cannot_use_another_regions_shortcut():
    surface, _, graph, records = line_inputs(5)
    for left in range(1, 4):
        graph[left][left + 1] = graph[left + 1][left] = 10.0
    graph[4].pop(5)
    graph[5].pop(4)
    graph[1][5] = graph[5][1] = 1.0
    graph[4][5] = graph[5][4] = 1.0
    distances = {p: module._distance(graph, p, set(graph)) for p in graph}
    by_position = {row.residue.uniprot_pos: row for row in surface}
    args = ([{1, 2, 3, 4}], by_position, {}, graph, records, distances, 1, 4, 50, 25, 22)
    assert module._sector_outputs(*args)[0].geodesic_diameter_a == 12
    with pytest.raises(ValueError, match='diameter constraints'):
        module._sector_outputs(*args, strict_diameter=True)


@pytest.mark.parametrize('ca', [None, (float('nan'), 0, 0), (float('inf'), 0, 0), (1, 2),
                                   (1, 2, 3, 4), 'invalid'])
@pytest.mark.parametrize('policy', ['balanced_connected_v4', 'connected_bounded_v3'])
def test_invalid_ca_cannot_silently_skip_hard_diameter(ca, policy):
    inputs = line_inputs(12)
    inputs[3][1].ca = ca
    audit = {}
    with pytest.raises(module.PartitionConstraintError, match='INVALID_PARTITION_INPUT'):
        module.build_epitope_sectors(*inputs, parameters={'policy': policy}, diagnostics=audit)
    assert audit['status'] == 'FAILED'
    assert 'INVALID_CA_COORDINATE:1' in audit['input_errors']


@pytest.mark.parametrize('weight', [-1, 0, float('nan'), float('inf')])
def test_invalid_graph_weight_fails_before_shortest_path_search(weight):
    inputs = line_inputs(12)
    inputs[2][1][2] = inputs[2][2][1] = weight
    with pytest.raises(module.PartitionConstraintError, match='INVALID_PARTITION_INPUT'):
        module.build_epitope_sectors(*inputs)


def test_missing_area_is_not_fabricated():
    surface, epitope, graph, records = line_inputs(12)
    with pytest.raises(module.PartitionConstraintError, match='INVALID_PARTITION_INPUT'):
        module.build_epitope_sectors(surface[1:], epitope, graph, records)


def test_small_component_remains_complete_with_explained_soft_exception():
    inputs = line_inputs(6)
    audit = {}
    result = module.build_epitope_sectors(*inputs, diagnostics=audit)
    assert_partition(result, inputs[0], inputs[2], inputs[3])
    assert len(result) == 1 and len(result[0].residue_positions) == 6
    assert audit['size_exceptions'] == [{'sector_id': 'VS001', 'residue_count': 6,
        'reason': 'COMPONENT_SIZE_BELOW_PREFERRED_RANGE',
        'description': '该独立表面连通分量小于偏好大小，不能跨越不连通表面归并。'}]


def test_budget_exhaustion_keeps_valid_incumbent_and_reports_incomplete_balance(pdl1, monkeypatch):
    monkeypatch.setattr(module, 'PARTITION_SEARCH_LIMIT', 80)
    audit = {}
    result = module.build_epitope_sectors(*pdl1, diagnostics=audit)
    assert_partition(result, pdl1[0], pdl1[2], pdl1[3])
    assert audit['status'] == 'COMPLETE'
    assert audit['optimization_status'] == 'SEARCH_BUDGET_EXHAUSTED'
    assert audit['search_states'] == 80
    assert audit['size_exceptions']
    assert all(row['reason'] == 'SEARCH_BUDGET_EXHAUSTED' for row in audit['size_exceptions'])


def test_budget_exhaustion_before_any_feasible_partition_is_explicit(pdl1, monkeypatch):
    monkeypatch.setattr(module, 'PARTITION_SEARCH_LIMIT', 2)
    audit = {}
    with pytest.raises(module.PartitionConstraintError, match='SEARCH_BUDGET_EXHAUSTED'):
        module.build_epitope_sectors(*pdl1, diagnostics=audit)
    assert audit['reason'] == 'SEARCH_BUDGET_EXHAUSTED'
    assert audit['search_states'] == 2


def test_area_and_size_hard_minimum_cannot_be_hidden():
    inputs = line_inputs(3)
    audit = {}
    with pytest.raises(module.PartitionConstraintError, match='UNSATISFIED_PARTITION_CONSTRAINTS'):
        module.build_epitope_sectors(*inputs, diagnostics=audit)
    assert audit['draft_regions'][0]['residue_positions'] == [1, 2, 3]
    assert 'MINIMUM_RESIDUES_NOT_MET' in audit['violations']
