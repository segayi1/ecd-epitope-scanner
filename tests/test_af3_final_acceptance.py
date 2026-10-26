import importlib.util
import json
from pathlib import Path

import pytest


def verifier():
    path = Path(__file__).resolve().parents[1] / "workspace_data/af3_binder_batch_20261017/verify_af3_batch.py"
    spec = importlib.util.spec_from_file_location("af3_final_acceptance", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def geometry_fixture(root):
    sectors = [{"sector_id": "r1", "residue_positions": [1, 2, 3, 4]}]
    data = {"probe_sectors_v4.json": sectors,
            "analysis.json": {"joint_v4_ml": {"sectors": sectors},
                              "residues": [{"uniprot_pos": i, "ca": [i * 3, 0, 0]} for i in range(1, 5)]},
            "surface_graph_v41.json": {"nodes": [{"uniprot_pos": i, "surface_area_a2": 20} for i in range(1, 5)],
                                       "edges": [{"source": i, "target": i + 1, "distance_a": 3} for i in range(1, 4)]}}
    for name, value in data.items():
        (root / name).write_text(json.dumps(value), encoding="utf-8")
    return data


def test_final_acceptance_independently_checks_connected_unique_geometry(tmp_path):
    module = verifier()
    data = geometry_fixture(tmp_path)
    assert module._geometry_evidence(tmp_path)["surface_residues"] == 4
    data["surface_graph_v41.json"]["edges"].pop(1)
    (tmp_path / "surface_graph_v41.json").write_text(json.dumps(data["surface_graph_v41.json"]))
    with pytest.raises(ValueError, match="independent geometry"):
        module._geometry_evidence(tmp_path)


def test_final_acceptance_rejects_duplicate_coverage(tmp_path):
    module = verifier()
    data = geometry_fixture(tmp_path)
    data["probe_sectors_v4.json"].append(data["probe_sectors_v4.json"][0])
    (tmp_path / "probe_sectors_v4.json").write_text(json.dumps(data["probe_sectors_v4.json"]))
    with pytest.raises(ValueError, match="independent geometry"):
        module._geometry_evidence(tmp_path)


def test_table_cache_preserves_parser_results_and_rejects_changed_bytes(tmp_path, monkeypatch):
    from ecd_epitope_scanner import stability
    module = verifier()
    table = tmp_path / "Dif.fxout"
    table.write_text("original")
    calls = []
    def parser(path):
        calls.append(path)
        return [("model.pdb", 0.25)]
    monkeypatch.setattr(stability, "_parse_dif_rows", parser)
    with module.cached_foldx_tables():
        assert stability._parse_dif_rows(table) == stability._parse_dif_rows(table) == [("model.pdb", 0.25)]
    assert len(calls) == 1 and stability._parse_dif_rows is parser
    with pytest.raises(ValueError, match="changed during acceptance"):
        with module.cached_foldx_tables():
            stability._parse_dif_rows(table)
            table.write_text("changed")
    assert stability._parse_dif_rows is parser
