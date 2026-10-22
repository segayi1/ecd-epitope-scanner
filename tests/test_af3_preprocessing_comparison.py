import pytest
from ecd_epitope_scanner.compare_af3_preprocessing import (
    _is_partition,
    adjusted_rand_index,
    jaccard,
    kabsch_rmsd,
)


def test_rigid_motion_has_zero_kabsch_rmsd():
    points = [[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3]]
    moved = [[-y + 4, x - 2, z + 9] for x, y, z in points]
    assert kabsch_rmsd(points, moved) < 1e-10


def test_partition_metrics_handle_identical_and_disjoint_sets():
    assert jaccard({1, 2}, {2, 3}) == 1 / 3
    assert jaccard(set(), set()) == 1.0
    assert adjusted_rand_index([[1, 2], [3]], [[1, 2], [3]]) == 1.0
    assert adjusted_rand_index([[1, 2], [3]], [[1], [2, 3]]) == pytest.approx(-0.5)
    assert _is_partition([[1, 2], [3]])
    assert not _is_partition([[1, 2], [2, 3]])
