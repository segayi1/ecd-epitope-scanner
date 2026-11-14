import pytest

from ecd_epitope_scanner.joint_design_v4 import sectors
from ecd_epitope_scanner.joint_design_v4.models import ResidueKey, SurfaceEvidence
from ecd_epitope_scanner.models import ResidueRecord


def test_region_diameter_cannot_use_another_regions_shortcut():
    records = {p: ResidueRecord(p, 'A', 'A', p, ca=(float(p), 0, 0)) for p in range(1, 6)}
    graph = {p: {} for p in records}
    for left in range(1, 4):
        graph[left][left + 1] = graph[left + 1][left] = 10.0
    graph[1][5] = graph[5][1] = 1.0
    graph[4][5] = graph[5][4] = 1.0
    surface = [SurfaceEvidence(ResidueKey('0', 'A', p, uniprot_pos=p), mapped_surface_area_a2=80) for p in graph]
    distances = {p: sectors._distance(graph, p, set(graph)) for p in graph}
    by_position = {row.residue.uniprot_pos: row for row in surface}
    args = ([{1, 2, 3, 4}], by_position, {}, graph, records, distances, 1, 4, 50, 25, 22)
    assert sectors._sector_outputs(*args)[0].geodesic_diameter_a == 12
    with pytest.raises(ValueError, match='diameter constraints'):
        sectors._sector_outputs(*args, strict_diameter=True)
