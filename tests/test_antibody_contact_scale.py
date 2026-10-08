import numpy as np
import hashlib
import http.client
import io
import json
import urllib.error

from ecd_epitope_scanner.scripts.build_antibody_contact_scale_reference import (
    _candidate_mask_overlaps, _digest_matches, _download_cif, _sample_contact_vertices, _surface_geodesics, build,
    _replay_frozen_cohort, _resolved_geodesic,
    _write_incomplete_contact_scale_diagnostic,
    _weighted_quantile,
)
from ecd_epitope_scanner.surface_perturbation.geometry import Mesh


def test_contact_vertex_sampling_is_area_proportional_and_reproducible():
    contact = np.arange(256)
    vertex_area = np.ones(256)
    vertex_area[:8] = 0.0
    atom_map = {(i,): (np.array([float(i), 0., 0.]), "ALA", "C") for i in range(256)}
    first, complete = _sample_contact_vertices(contact, vertex_area, atom_map)
    second, _ = _sample_contact_vertices(contact, vertex_area, atom_map)
    assert not complete
    assert len(first) == 128 and len(set(first)) == len(first)
    assert np.array_equal(first, second)
    assert np.all(vertex_area[first] > 0)


def _tetrahedron(offset=0.0):
    vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    vertices[:, 0] += offset
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    return vertices, faces


def test_surface_geodesics_builds_nonzero_edges_and_positive_distances():
    vertices, faces = _tetrahedron()
    mesh = Mesh(vertices, faces)

    class Tools:
        @staticmethod
        def ses(_atom_map, scale):
            assert scale == 2.0
            return mesh

    result = _surface_geodesics(Tools(), {}, np.array([[0., 0., 0.]]))
    assert result["status"] == "RESOLVED"
    assert result["connected_pair_fraction"] == 1.0
    assert result["interface_surface_components"] == 1
    assert all(value > 0 for value in result["pairwise_surface_distance_quantiles_a"].values())


def test_disconnected_contact_components_report_within_component_distances():
    vertices_a, faces_a = _tetrahedron()
    vertices_b, faces_b = _tetrahedron(100.0)
    mesh = Mesh(np.concatenate((vertices_a, vertices_b)), np.concatenate((faces_a, faces_b + 4)))

    class Tools:
        @staticmethod
        def ses(_atom_map, scale):
            assert scale == 2.0
            return mesh

    result = _surface_geodesics(Tools(), {}, np.array([[0., 0., 0.], [100., 0., 0.]]))
    assert result["status"] == "RESOLVED_COMPONENTWISE"
    assert result["connected_pair_fraction"] < 1.0
    assert result["distance_scope"].startswith("within connected antigen SES components")
    assert all(value > 0 for value in result["pairwise_surface_distance_quantiles_a"].values())
    assert _resolved_geodesic({"geodesic": result})


def test_frozen_cohort_replay_rejects_a_different_input_manifest(tmp_path):
    selection = tmp_path / "selection.json"
    selection.write_text(json.dumps({"profile": "contact_scale_reference_v1",
                                     "source_manifest_sha256": "a" * 64,
                                     "selected_structure_count": 0, "samples": []}), encoding="utf-8")
    try:
        _replay_frozen_cohort(selection, [], tmp_path, "A", "b" * 64, None)
    except ValueError as exc:
        assert "does not match this contact input manifest" in str(exc)
    else:
        raise AssertionError("frozen selection from a different input identity must be rejected")


def test_weighted_contact_distance_quantile_uses_surface_area_weights():
    assert _weighted_quantile(np.array([0., 10.]), np.array([9., 1.]), .85) == 0.
    assert _weighted_quantile(np.array([0., 10.]), np.array([9., 1.]), .95) == 10.


def test_saved_fab_contact_masks_can_be_reused_for_an_extended_panel(tmp_path):
    geometry = tmp_path / "geometry" / "candidate"
    apbs = tmp_path / "apbs" / "candidate"
    geometry.mkdir(parents=True); apbs.mkdir(parents=True)
    np.savez(geometry / "geometry_fusion_masks.npz",
             geometry_chemistry__positive=np.array([True, False, False]),
             geometry_chemistry__unknown=np.array([False, True, False]))
    np.savez(apbs / "field_fusion_evidence.npz",
             three_mechanism_positive=np.array([False, True, False]),
             three_mechanism_unknown=np.array([False, False, True]))

    overlaps = _candidate_mask_overlaps(["candidate"], {10, 11}, np.array([10, 11, 12]),
                                        np.array([1., 2., 3.]), tmp_path / "geometry", tmp_path / "apbs")
    assert overlaps["candidate"]["geometry_chemistry"]["overlap_fraction"] == 1 / 3
    assert overlaps["candidate"]["geometry_chemistry_unknown_opportunity"]["overlap_fraction"] == 2 / 3
    assert overlaps["candidate"]["three_mechanism"]["overlap_fraction"] == 2 / 3
    assert overlaps["candidate"]["three_mechanism_unknown_opportunity"]["overlap_fraction"] == 0.


def test_frozen_sha_text_is_case_insensitive(tmp_path):
    path = tmp_path / "source.bin"; path.write_bytes(b"frozen bytes")
    assert _digest_matches(path, hashlib.sha256(path.read_bytes()).hexdigest().upper())


def test_structure_download_retries_transient_eof_and_commits_atomically(tmp_path, monkeypatch):
    calls = []

    def urlopen(_request, timeout):
        calls.append(timeout)
        if len(calls) == 1:
            raise urllib.error.URLError("transient TLS EOF")
        if len(calls) == 2:
            raise http.client.IncompleteRead(b"partial")
        return io.BytesIO(b"x" * 1000)

    monkeypatch.setattr("ecd_epitope_scanner.scripts.build_antibody_contact_scale_reference.urllib.request.urlopen", urlopen)
    result = _download_cif({"pdb_code": "test"}, tmp_path)
    assert result.read_bytes() == b"x" * 1000
    assert calls == [90, 90, 90]
    assert not (tmp_path / "test.cif.part").exists()


def test_open_contact_surface_is_local_unknown_but_other_solver_errors_propagate():
    class BrokenSurface:
        def __init__(self, error):
            self.error = error

        def ses(self, _atom_map, scale):
            assert scale == 2.0
            raise self.error

    result = _surface_geodesics(BrokenSurface(RuntimeError("SES mesh is not closed")), {}, np.zeros((1, 3)))
    assert result == {"status": "UNKNOWN_OPEN_SURFACE", "reason": "SES mesh is not closed"}

    try:
        _surface_geodesics(BrokenSurface(RuntimeError("solver identity mismatch")), {}, np.zeros((1, 3)))
    except RuntimeError as exc:
        assert str(exc) == "solver identity mismatch"
    else:
        raise AssertionError("non-topology solver errors must remain fatal")


def test_below_minimum_contact_set_persists_diagnostics_without_reference(tmp_path):
    class Tools:
        nano = "NanoShaper"
        calls = {"nanoshaper": 1}
        activity = [{"runtime_seconds": 0.5}]

        @staticmethod
        def identity(_executable):
            return "a" * 64

    path = _write_incomplete_contact_scale_diagnostic(
        tmp_path,
        target_clusters=30,
        minimum_clusters=15,
        candidates=4274,
        accepted=[{"cluster_id": "ANTIGEN_CLUSTER_01"}],
        excluded=[{"reason": "HOMOLOGOUS_TO_CLUSTER_01"}],
        source_manifest_sha256="b" * 64,
        tools=Tools(),
        started=0.0,
    )
    result = json.loads(path.read_text(encoding="utf-8"))
    assert result["status"] == "INSUFFICIENT_SEQUENCE_CLUSTERS"
    assert result["selected_structure_count"] == 1
    assert result["minimum_required"] == 15
    assert result["minimum_policy"]["default_minimum"] == 20
    assert result["minimum_policy"]["override_used"] is True
    assert result["samples"][0]["cluster_id"] == "ANTIGEN_CLUSTER_01"
    assert result["excluded_attempts"][0]["reason"] == "HOMOLOGOUS_TO_CLUSTER_01"
    assert not (tmp_path / "contact_scale_reference.json").exists()


def test_contact_cluster_bounds_are_checked_before_reading_study_inputs(tmp_path):
    try:
        build(tmp_path, target_clusters=14, minimum_clusters=15)
    except ValueError as exc:
        assert "1 <= minimum_clusters <= target_clusters" in str(exc)
    else:
        raise AssertionError("minimum cohort larger than target must be rejected")
