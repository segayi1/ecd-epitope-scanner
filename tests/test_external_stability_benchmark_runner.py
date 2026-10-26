import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "run_external_stability_benchmark", ROOT / "scripts/run_external_stability_benchmark.py"
)
benchmark = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(benchmark)
SUMMARY_SPEC = importlib.util.spec_from_file_location(
    "summarize_external_stability_benchmark", ROOT / "scripts/summarize_external_stability_benchmark.py"
)
summarizer = importlib.util.module_from_spec(SUMMARY_SPEC)
SUMMARY_SPEC.loader.exec_module(summarizer)


def test_runner_obeys_approved_pair_cap_and_counts_failed_start_reservations(tmp_path, monkeypatch):
    runner_manifest = tmp_path / "runner_execution_manifest.json"
    runner_manifest.write_text(json.dumps({"approved_pair_slots": 10, "max_pair_slots": 300}), encoding="utf-8")
    ledger = tmp_path / "execution_ledger.jsonl"
    ledger.write_text(json.dumps({"event": "foldx_start", "candidate_id": "started_without_result",
                                  "reserved_pair_slots": 5}) + "\n", encoding="utf-8")
    monkeypatch.setattr(benchmark, "RUNNER_MANIFEST", runner_manifest)
    monkeypatch.setattr(benchmark, "LEDGER", ledger)
    monkeypatch.setattr(benchmark, "OUT", tmp_path / "predictions")
    monkeypatch.setattr(benchmark, "run_worker", lambda *args: pytest.fail("budget check must precede worker"))

    assert benchmark.approved_pair_slot_cap() == 10
    with pytest.raises(RuntimeError, match=r"approved pair-slot budget exceeded: 5\+10>10"):
        benchmark.run_foldx([{"source_variant_id": "started_without_result"},
                             {"source_variant_id": "new_a"}, {"source_variant_id": "new_b"}], 3)

    assert len(ledger.read_text(encoding="utf-8").splitlines()) == 1


def test_foldx_result_pair_must_map_to_the_started_candidate():
    candidate = {"chain": "A", "mutations_normalized": ["D84H", "M129N"],
                 "normalized_pdb_sha256": "a" * 64}
    pair = {"chain": "A", "mutations": ["D84H", "M129N"],
            "normalized_input_sha256": "a" * 64}
    assert summarizer.candidate_pair_identity_matches(candidate, pair)
    assert not summarizer.candidate_pair_identity_matches(candidate, pair | {"mutations": ["D84H", "M130N"]})
    assert not summarizer.candidate_pair_identity_matches(candidate, pair | {"chain": "B"})
    assert not summarizer.candidate_pair_identity_matches(candidate, pair | {"normalized_input_sha256": "b" * 64})


def test_normalized_benchmark_requests_and_pair_audit_use_chain_a(tmp_path, monkeypatch):
    ledger = tmp_path / "execution_ledger.jsonl"
    runner_manifest = tmp_path / "runner_execution_manifest.json"
    runner_manifest.write_text(json.dumps({"approved_pair_slots": 5, "max_pair_slots": 300}), encoding="utf-8")
    monkeypatch.setattr(benchmark, "LEDGER", ledger)
    monkeypatch.setattr(benchmark, "RUNNER_MANIFEST", runner_manifest)
    monkeypatch.setattr(benchmark, "OUT", tmp_path / "predictions")
    requests = []

    def fake_worker(command, request, timeout, log_path):
        requests.append(request)
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        return {"status": "success", "records": []}, 0.01, "{}"

    monkeypatch.setattr(benchmark, "run_worker", fake_worker)
    row = {
        "source_variant_id": "example:B2V",
        "normalized_pdb": str(tmp_path / "normalized.pdb"),
        "normalized_pdb_sha256": "a" * 64,
        "source_pdb_sha256": "b" * 64,
        "chain": "B",
        "sequence": "ACDEFG",
        "mutations_normalized": ["A1V", "C2D"],
        "order": 2,
    }

    assert benchmark.run_spurs([row]) == 1
    assert len(benchmark.run_foldx([row], 1)) == 1
    assert [request["chain_id"] for request in requests] == ["A", "A"]
    assert ":" not in Path(requests[1]["output_dir"]).name
    assert summarizer.candidate_pair_identity_matches(
        row,
        {"chain": "A", "mutations": ["A1V", "C2D"], "normalized_input_sha256": "a" * 64},
    )
    assert not summarizer.candidate_pair_identity_matches(
        row,
        {"chain": "B", "mutations": ["A1V", "C2D"], "normalized_input_sha256": "a" * 64},
    )


def test_remaining_resume_matches_numeric_source_ids_to_string_ledger_ids(tmp_path, monkeypatch):
    ledger = tmp_path / "execution_ledger.jsonl"
    ledger.write_text(json.dumps({"event": "foldx_start", "candidate_id": "37",
                                  "reserved_pair_slots": 5}) + "\n", encoding="utf-8")
    runner_manifest = tmp_path / "runner_execution_manifest.json"
    runner_manifest.write_text(json.dumps({"approved_pair_slots": 10, "max_pair_slots": 300}), encoding="utf-8")
    monkeypatch.setattr(benchmark, "LEDGER", ledger)
    monkeypatch.setattr(benchmark, "RUNNER_MANIFEST", runner_manifest)
    rows = [
        {"source_variant_id": 37, "order": 1, "sequence_position": 1},
        {"source_variant_id": "new_variant", "order": 1, "sequence_position": 2},
    ]
    monkeypatch.setattr(benchmark, "read_candidates", lambda: rows)
    monkeypatch.setattr(benchmark, "run_spurs", lambda unused_rows: 0)
    selected = []

    def capture_foldx(selected_rows, limit):
        selected.extend(row["source_variant_id"] for row in selected_rows)
        return selected_rows

    monkeypatch.setattr(benchmark, "run_foldx", capture_foldx)
    benchmark.run(type("Args", (), {"pilot": False, "remaining": True, "skip_spurs": False})())

    assert selected == ["new_variant"]


def test_incomplete_benchmark_status_tracks_terminal_runner_receipt():
    assert summarizer.benchmark_status(False, False, "RUNNING") == "PILOT_OR_RUNNING"
    assert summarizer.benchmark_status(False, False, "INTERRUPTED_NO_TERMINAL_RECEIPT") == "INTERRUPTED_INCOMPLETE"
    assert summarizer.benchmark_status(True, False, "INTERRUPTED_NO_TERMINAL_RECEIPT") == "COMPLETE_WITH_MODEL_FAILURES"
    assert summarizer.benchmark_status(True, True, "FINISHED") == "COMPLETE"


def test_runner_lock_blocks_a_second_process(tmp_path):
    lock_path = tmp_path / "runner.lock"
    script = "\n".join([
        "import sys",
        "from pathlib import Path",
        f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
        "from run_external_stability_benchmark import exclusive_runner_lock",
        f"with exclusive_runner_lock(Path({str(lock_path)!r})): pass",
    ])
    with benchmark.exclusive_runner_lock(lock_path):
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode != 0
    assert "External benchmark runner is already active" in result.stderr


def test_benchmark_summary_keeps_worker_failures_and_separates_preflight_status(tmp_path, monkeypatch):
    row = {
        "source_variant_id": "candidate-1",
        "source": "test",
        "order": 2,
        "sequence_position": 1,
        "sequence": "ACDEFG",
        "mutations_normalized": ["A1V", "C2D"],
        "chain": "B",
        "normalized_pdb_sha256": "a" * 64,
        "experimental_ddg_kcal_mol": 0.4,
        "experimental_ddg_source_value": 0.4,
    }
    (tmp_path / "candidate_manifest.json").write_text(json.dumps({
        "summary": {"prediction_status": "NOT_RUN_COORDINATED_WITH_PRODUCTION_MODEL_WORK"},
        "candidates": [row],
    }), encoding="utf-8")
    (tmp_path / "candidates.csv").write_text(
        "source_variant_id,train_cluster_overlap,train_single_component_overlap\n"
        "candidate-1,{},{}\n", encoding="utf-8")
    events = [
        {"event": "spurs_start", "candidate_ids": ["candidate-1"]},
        {"event": "spurs_result", "candidate_ids": ["candidate-1"], "status": "failed",
         "key": "spurs:test", "elapsed_seconds": 1.0,
         "payload": {"errors": [{"code": "TOOL_CRASHED", "message": "chain not found"}]}},
    ]
    (tmp_path / "execution_ledger.jsonl").write_text(
        "\n".join(json.dumps(event) for event in events) + "\n", encoding="utf-8")
    (tmp_path / "runner_execution_manifest.json").write_text(json.dumps({
        "status": "RUNNING", "max_pair_slots": 300, "approved_pair_slots": 5,
        "prior_attempt_recovery": {"snapshot_file": "failed_attempt.json"},
    }), encoding="utf-8")
    monkeypatch.setattr(summarizer, "BASE", tmp_path)

    summarizer.main()

    summary = json.loads((tmp_path / "benchmark_results.json").read_text(encoding="utf-8"))
    assert summary["status"] == summary["prediction_status"] == "PILOT_OR_RUNNING"
    assert summary["preflight_prediction_status"] == "NOT_RUN_COORDINATED_WITH_PRODUCTION_MODEL_WORK"
    assert summary["spurs"]["failed_candidate_ids"] == ["candidate-1"]
    assert summary["spurs"]["failed_events"][0]["errors"][0]["message"] == "chain not found"
    assert summary["prior_attempt_recovery"]["snapshot_file"] == "failed_attempt.json"
