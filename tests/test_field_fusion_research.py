import json
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from types import SimpleNamespace
from ecd_epitope_scanner.scripts.run_ctla4_field_fusion import (
    LABELS, PANEL_IDS, _aggregate_numerical_bounds, _grid_contains,
)
from ecd_epitope_scanner.scripts import run_ctla4_field_fusion as research_runner
from ecd_epitope_scanner.surface_perturbation import engine

from ecd_epitope_scanner.surface_perturbation.field_fusion_research import (
    fuse_scope, geometry_summary, pair_null_blocks, validate_null_roles,
)


def evidence(n=2):
    train = np.tile(np.linspace(-.2, .2, 5)[:, None], (1, n))
    held = np.zeros((5, n))
    prod = np.zeros((5, n))
    error = np.full(n, .01)
    return prod, train, held, error


def test_sixth_candidate_label_does_not_redefine_the_frozen_five_candidate_panel():
    assert len(PANEL_IDS) == 5
    assert "CV2_5baddd4518a94929" not in PANEL_IDS
    assert LABELS["CV2_5baddd4518a94929"] == "S55F/I102R (coverage extension)"


def run_fusion(prod_a, prod_b, held_a=None):
    null_only, train, held_values, error = evidence()
    held = {'shape_sdf_a': held_values.copy(), 'nonpolar_support': held_values.copy()}
    if held_a is not None:
        held['shape_sdf_a'] = held_a
    return fuse_scope(
        {'shape_sdf_a': prod_a, 'nonpolar_support': prod_b},
        {'shape_sdf_a': train, 'nonpolar_support': null_only + train},
        held,
        {'shape_sdf_a': error, 'nonpolar_support': error},
        np.ones(2), scalar_names=('shape_sdf_a', 'nonpolar_support'),
        vector_groups={}, protocol='test',
    )


def test_pair_null_blocks_uses_disjoint_fixed_slot_pairs():
    assert np.array_equal(pair_null_blocks(np.arange(10)[:, None]).ravel(), [-1, -1, -1, -1, -1])


def test_apbs_grid_containment_rejects_nonfinite_and_out_of_bounds_atoms():
    grid = {'origin': [0., 0., 0.], 'lengths': [10., 10., 10.]}
    assert _grid_contains(grid, np.array([[0., 5., 10.]]))
    assert not _grid_contains(grid, np.array([[10.01, 5., 5.]]))
    assert not _grid_contains(grid, np.array([[np.nan, 5., 5.]]))


def test_global_apbs_grid_allows_only_append_only_extension_with_identical_geometry(tmp_path, monkeypatch):
    def atom_map(x):
        return {('A', 1, '', 'CA'): (np.array([x, 1., 1.]), 'ALA', 'C')}

    rows = {}
    for cid, x in (('c1', 2.), ('c2', 3.), ('c3', 12.)):
        rows[cid] = [{'wt_sha256': f'{cid}-wt-{i}', 'mut_sha256': f'{cid}-mut-{i}',
                      'wt_atoms': atom_map(x), 'mut_atoms': atom_map(x + .1)} for i in range(5)]

    monkeypatch.setattr(research_runner, '_role_rows', lambda manifest, cid, role: rows[cid] if role == 'production' else [])
    monkeypatch.setattr(research_runner, '_apbs_root', lambda study, manifest: tmp_path / 'apbs')
    monkeypatch.setattr(research_runner, 'reference_id', lambda _reference: 'frozen-reference')
    monkeypatch.setattr(research_runner, 'aligned_atoms',
                        lambda row, native, records: (row['wt_atoms'], row['mut_atoms'], {'core_rmsd_a': 0.}))

    def fake_grid(coords, spacing, margin):
        extent = 10. if float(np.max(coords[:, 0])) < 10. else 20.
        return {'counts': [int(round(extent / spacing))] * 3,
                'lengths': [extent] * 3, 'center': [extent / 2] * 3,
                'spacing': [spacing] * 3, 'origin': [0., 0., 0.]}

    monkeypatch.setattr(research_runner, 'make_grid', fake_grid)
    monkeypatch.setattr(research_runner, '_grid_contains', lambda _grid, _coords: True)
    native = atom_map(1.)
    bank = {'candidates': {cid: {'null_train': [], 'null_holdout': []} for cid in rows}}
    base = {'candidate_ids': ['c1'], '_independent_null_bank': bank}
    reference = {'area': np.ones(1)}
    original = research_runner._fixed_apbs_grids(tmp_path, base, native, [], reference)

    extended = {'candidate_ids': ['c1', 'c2'], '_independent_null_bank': bank}
    assert research_runner._fixed_apbs_grids(tmp_path, extended, native, [], reference) == original
    contract_path = tmp_path / 'apbs' / 'apbs_grid_contract.json'
    saved = json.loads(contract_path.read_text())
    assert saved['contract']['candidate_ids'] == ['c1', 'c2']
    assert saved['contract_history'][0]['contract']['candidate_ids'] == ['c1']

    changed = {'candidate_ids': ['c1', 'c2', 'c3'], '_independent_null_bank': bank}
    with pytest.raises(ValueError, match='panel extension requires identical grids'):
        research_runner._fixed_apbs_grids(tmp_path, changed, native, [], reference)
    assert json.loads(contract_path.read_text())['contract']['candidate_ids'] == ['c1', 'c2']


def test_null_role_contract_counts_distinct_runs_and_reports_repeated_coordinates():
    train = [{'wt_sha256': f't{i}', 'coordinate_sha256': 'same',
              'batch_path': f'train{i}', 'group': 1, 'run': i} for i in range(10)]
    held = [{'wt_sha256': f'h{i}', 'coordinate_sha256': 'same',
             'batch_path': f'held{i}', 'group': 1, 'run': i} for i in range(10)]
    production = [{'wt_sha256': 'p', 'coordinate_sha256': 'same',
                   'batch_path': 'prod', 'group': 1, 'run': i} for i in range(5)]
    audit = validate_null_roles(train, held, production)
    assert audit['null_train']['unique_run_identities'] == 10
    assert audit['null_train']['unique_coordinate_models'] == 1
    with pytest.raises(ValueError, match='reuses a FoldX run identity'):
        validate_null_roles([*train[:-1], train[0]], held, production)
    with pytest.raises(ValueError, match='overlap'):
        validate_null_roles(train, [*held[:-1], {'batch_path': 'train0', 'group': 1, 'run': 0,
            'wt_sha256': 'h9', 'coordinate_sha256': 'same'}], production)


def test_fusion_requires_same_channel_four_of_five_and_preserves_unknown():
    a = np.zeros((5, 2)); b = np.zeros((5, 2))
    a[:2, 0] = 4.; b[2:4, 0] = 4.
    a[:, 1] = 4.
    result = run_fusion(a, b)
    assert not result['positive'][0]  # two plus two cannot be pooled
    assert result['positive'][1]
    assert not result['unknown'][1]


def test_holdout_same_channel_quorum_marks_that_position_unknown():
    a = np.zeros((5, 2)); b = np.zeros((5, 2)); a[:, 0] = 4.; a[:, 1] = 4.
    held = np.zeros((5, 2)); held[:4, 1] = 4.
    result = run_fusion(a, b, held)
    assert result['positive'][0]
    assert result['unknown'][1]
    assert not result['positive'][1]
    assert result['channel_unknown_reasons']['shape_sdf_a']['holdout_quorum_failure'][1]


def test_unknown_channel_does_not_clear_independently_closed_positive():
    geometry = np.zeros((5, 2)); chemistry = np.zeros((5, 2))
    geometry[:4, 0] = 4.
    chemistry[:, 0] = np.nan
    result = run_fusion(geometry, chemistry)
    assert result['positive'][0]
    assert result['channel_unknown']['nonpolar_support'][0]
    assert not result['unknown'][0]


def test_vector_channels_use_joint_magnitude_and_direction_shapes():
    n = 2
    production = {key: np.zeros((5, n)) for key in ('nx', 'ny', 'nz')}
    for key in production:
        production[key][:4, 0] = 1.
    estimate = {key: pair_null_blocks(np.zeros((10, n))) for key in production}
    heldout = {key: pair_null_blocks(np.zeros((10, n))) for key in production}
    error = {key: np.full(n, .01) for key in production}
    result = fuse_scope(production, estimate, heldout, error, np.ones(n), scalar_names=(),
                        vector_groups={'surface_normal_delta': ('nx', 'ny', 'nz')}, protocol='test')
    assert result['channel_positive']['surface_normal_delta'][0]
    assert result['channel_negative']['surface_normal_delta'][1]


def test_geometry_summary_keeps_irregular_mask_and_reports_weighted_geodesic_quantiles():
    graph = csr_matrix(np.array([[0., 1., 0.], [1., 0., 2.], [0., 2., 0.]]))
    result = geometry_summary(np.array([True, True, False]), np.array([1., 3., 2.]), graph,
                             np.array([True, False, False]))
    assert result['positive_area_a2'] == 4.
    assert result['equivalent_area_radius_a'] > 0
    assert result['geodesic_distance_area_quantiles_a']['0.9'] == 1.


def test_geodesic_summary_does_not_bridge_disconnected_surface_sides():
    graph = csr_matrix(np.array([[0., 2., 0., 0.], [2., 0., 0., 0.],
                                 [0., 0., 0., 1.], [0., 0., 1., 0.]]))
    result = geometry_summary(np.array([True, False, True, False]), np.ones(4), graph,
                             np.array([True, False, False, False]))
    assert result['unreachable_positive_samples'] == 1


def test_nonclosed_numerical_mesh_is_cached_as_unknown_and_charged_once(monkeypatch, tmp_path):
    pair = {'wt_sha256': 'wt', 'mut_sha256': 'mut'}
    reference = {'area': np.ones(3)}
    tools = SimpleNamespace(config={}, nano='nano', identity=lambda _path: 'nano-hash')
    monkeypatch.setattr(engine, 'protocol_identity', lambda *_args, **_kwargs: 'protocol')
    calls = []

    def nonclosed(*_args):
        calls.append(1)
        raise RuntimeError('SES mesh is not closed')

    monkeypatch.setattr(engine, 'numerical_bounds', nonclosed)
    first, audit = engine.cached_numerical_bounds(pair, {}, [], reference, tools, {}, tmp_path,
                                                   return_audit=True)
    second, cached_audit = engine.cached_numerical_bounds(pair, {}, [], reference, tools, {}, tmp_path,
                                                           return_audit=True)
    assert audit == cached_audit == {'status': 'UNRESOLVED', 'reason': 'SES mesh is not closed'}
    assert set(first) == set(second) == set(engine.FIELD_CHANNELS)
    assert all(np.isnan(first[name]).all() and np.isnan(second[name]).all() for name in engine.FIELD_CHANNELS)
    assert calls == [1]


def test_four_valid_numerical_slots_survive_one_unknown_and_three_remain_unknown():
    rows = np.array([[.1, .1, .1], [.2, .2, .2], [.3, .3, .3],
                     [.4, np.nan, .4], [np.nan, np.nan, .5]])
    error_runs = [{channel: row.copy() for channel in engine.FIELD_CHANNELS} for row in rows]
    result = _aggregate_numerical_bounds(error_runs)
    assert all(np.isclose(value[0], .4) and np.isnan(value[1]) and np.isclose(value[2], .5)
               for value in result.values())
