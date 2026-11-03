"""Replay archived FoldX pairs through the repaired site-first worker without invoking FoldX."""
import copy
import json
import statistics
from pathlib import Path

import pytest


def test_saved_foldx_pairs_pass_worker_with_legacy_rounded_aggregates(tmp_path, monkeypatch):
    from ecd_epitope_scanner import stability
    from ecd_epitope_scanner.integrations.workers import foldx_worker
    from ecd_epitope_scanner.surface_perturbation import paired

    root = Path(__file__).resolve().parents[1]
    ledger = root / "workspace_data/d14576822c6f/ordinary_candidate_ledger.json"
    if not ledger.is_file():
        pytest.skip("Saved real FoldX batch is unavailable")
    batches = json.loads(ledger.read_text(encoding="utf-8"))["foldx"]["batches"]
    records = [batch["stages"][0]["records"][0] for batch in batches]
    pairs = [copy.deepcopy(pair) for row in records for pair in row["paired_models"]]
    assert len(records) == 8 and len(pairs) == 24
    medians = [statistics.median(pair["ddg_kcal_mol"] for pair in row["paired_models"]) for row in records]
    monkeypatch.setattr(foldx_worker, "_parse_mutations", lambda mutations, chain, sequence: mutations)
    monkeypatch.setattr(stability, "score_many_with_foldx", lambda *args, **kwargs: [round(x, 4) for x in medians])
    monkeypatch.setattr(paired, "scan_pairs", lambda *args, **kwargs: (pairs, []))
    out, invocation = tmp_path / "output", tmp_path / "invocation"
    out.mkdir()
    invocation.mkdir()
    result = foldx_worker._site_first_batch(
        {"parameters": {"number_of_runs": 3, "sampling_protocol": "time_rotabase_v3"},
         "mutations": [row["mutations"] for row in records], "chain_id": "A", "sequence": ""},
        tmp_path / "foldx.exe", tmp_path / "repaired.pdb", tmp_path, out, invocation,
        {"binary_sha256": pairs[0]["binary_sha256"], "structure_sha256": pairs[0]["normalized_input_sha256"]}, 0.0, False)
    assert result["status"] == "success"
    assert [row["ddg"] for row in result["records"]] == medians
    assert all(len(row["paired_models"]) == 3 for row in result["records"])
