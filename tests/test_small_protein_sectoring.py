"""Length-aware display preferences keep geometry and ownership hard gates."""
from __future__ import annotations

import pytest

from ecd_epitope_scanner.config import JOINT_V4_ML_PARAMETERS, display_sectoring_parameters
from ecd_epitope_scanner.joint_design_v4 import sectors
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, SurfaceEvidence
from ecd_epitope_scanner.models import ResidueRecord


def graph_inputs(component_sizes):
    records = {}
    graph = {}
    surface = []
    position = 1
    for component_index, size in enumerate(component_sizes):
        component = list(range(position, position + size))
        x_offset = component_index * 100.0
        for offset, p in enumerate(component):
            records[p] = ResidueRecord(p, 'A', 'A', p, ca=(x_offset + float(offset), 0.0, 0.0))
            graph[p] = {}
            surface.append(SurfaceEvidence(
                ResidueKey('0', 'A', p, uniprot_pos=p), mapped_surface_area_a2=80.0))
            if offset:
                previous = component[offset - 1]
                graph[p][previous] = graph[previous][p] = 1.0
        position += size + 1
    return surface, [], graph, records


def assert_complete_hard_constraints(result, surface, graph, records):
    assigned = [p for sector in result for p in sector.residue_positions]
    assert len(assigned) == len(set(assigned))
    assert set(assigned) == set(graph)
    area = {row.residue.uniprot_pos: row.mapped_surface_area_a2 for row in surface}
    for sector in result:
        group = set(sector.residue_positions)
        assert len(group) >= 4
        assert sum(area[p] for p in group) >= 50.0
        assert set(sectors._distance(graph, min(group), group)) == group
        full_distances = {p: sectors._distance(graph, p, set(graph)) for p in group}
        for p in group:
            for q in group:
                assert full_distances[p][q] <= 30.0
                assert sum((records[p].ca[i] - records[q].ca[i]) ** 2 for i in range(3)) ** 0.5 <= 22.0


def test_length_defaults_and_overrides_are_isolated():
    original = JOINT_V4_ML_PARAMETERS['sectoring'].copy()
    small = display_sectoring_parameters(100)
    large = display_sectoring_parameters(151)
    explicit = display_sectoring_parameters(100, {'target_sector_residues': 9})

    assert (small['target_sector_residues'], small['preferred_minimum_sector_residues'],
            small['preferred_maximum_sector_residues']) == (7, 6, 8)
    assert (small['preferred_minimum_sector_count'], small['preferred_maximum_sector_count']) == (4, 12)
    assert (small['minimum_sector_residues'], small['minimum_sector_area_a2'],
            small['maximum_geodesic_diameter_a'], small['maximum_euclidean_diameter_a']) == (4, 50.0, 30.0, 22.0)
    assert large == {**original, 'preferred_minimum_sector_count': 4, 'preferred_maximum_sector_count': 12}
    assert explicit['target_sector_residues'] == 9
    assert explicit['preferred_minimum_sector_residues'] <= 9 <= explicit['preferred_maximum_sector_residues']
    assert JOINT_V4_ML_PARAMETERS['sectoring'] == original


@pytest.mark.parametrize('length,target', [(64, 7), (119, 7), (120, 8), (150, 8), (151, 12), (614, 12)])
def test_automatic_size_uses_analyzed_length_only(length, target):
    params = display_sectoring_parameters(length)
    assert params['target_sector_residues'] == target
    assert (params['preferred_minimum_sector_count'], params['preferred_maximum_sector_count']) == (4, 12)


def test_small_chain_prefers_eight_residues_per_connected_sector():
    inputs = graph_inputs([80])
    params = display_sectoring_parameters(80)
    audit = {}
    result = sectors.build_epitope_sectors(*inputs, parameters=params, diagnostics=audit)

    assert_complete_hard_constraints(result, inputs[0], inputs[2], inputs[3])
    assert 10 <= len(result) <= 13
    assert 6 <= sum(len(row.residue_positions) for row in result) / len(result) <= 8
    assert audit['preferred_minimum_sector_count'] == 4
    assert audit['preferred_maximum_sector_count'] == 12
    assert audit['count_preference_met'] is True


def test_count_range_is_soft_when_components_cannot_be_joined():
    inputs = graph_inputs([4] * 13)
    params = display_sectoring_parameters(52)
    audit = {}
    result = sectors.build_epitope_sectors(*inputs, parameters=params, diagnostics=audit)

    assert_complete_hard_constraints(result, inputs[0], inputs[2], inputs[3])
    assert len(result) == 13
    assert audit['count_range_is_soft'] is True
    assert audit['count_preference_met'] is False
    assert audit['count_range_exception']['reason'] == 'SOFT_COUNT_RANGE_EXCEEDED'


def test_explicit_ten_residue_target_keeps_existing_size_preferences():
    params = display_sectoring_parameters(100, {'target_sector_residues': 10})
    assert params['target_sector_residues'] == 10
    assert params['preferred_minimum_sector_residues'] == 10
    assert params['preferred_maximum_sector_residues'] == 15


@pytest.mark.parametrize('invalid_bound', [True, False, 3.5, 4.0, '4'])
def test_sector_count_bounds_reject_non_integer_types(invalid_bound):
    inputs = graph_inputs([8])
    with pytest.raises(ValueError, match='Preferred sector-count bounds must be positive integers'):
        sectors.build_epitope_sectors(*inputs, parameters={
            'preferred_minimum_sector_count': invalid_bound,
            'preferred_maximum_sector_count': 12,
        })

