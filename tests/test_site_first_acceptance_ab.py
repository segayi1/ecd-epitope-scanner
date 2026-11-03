import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_site_first_acceptance_ab import compare
import run_site_first_acceptance_ab as controller


def test_comparison_requires_matched_evidence_and_scientific_acceptance():
    identity = {"worker_model_identities": {name: {"model_files_sha256": {"w": "x"}, "device": "cuda", "gpu_name": "GPU"}
        for name in ("spurs_single", "spurs_multi", "discotope")}}
    cache = {"score_cache_entries_at_start": 0, "all_possible_repair_cache_paths_absent_before_run": True}
    digest = "8ff23af24604861cec56365dc835068f97fa714222b8313a71291da2dc7851e4"
    new = {"status": "PASS", "input_sha256": digest, "cache_evidence": {"before_run": {"cold_verified": True}}, "elapsed_seconds": 50}
    old = {"status": "COMPLETE", "input_sha256": digest, "cold_cache_evidence": cache,
           "elapsed_pipeline_seconds": 100, "normalized_structure": {"sha256": "s"}, "model_identity": {"x": 1}}
    result = compare(new, old, {"x": 1}, identity, identity, "s")
    assert result["performance_target"] == "PASS" and result["reduction_percent"] == 50
    assert compare(new, old, {"x": 2}, identity, identity, "s")["performance_target"] == "NOT_ASSESSED"
    assert compare({**new, "status": "NEEDS_REVIEW"}, old, {"x": 1}, identity, identity, "s")["performance_target"] == "NOT_ASSESSED"


def test_explicit_new_attempt_keeps_old_results_and_stops_before_baseline(tmp_path, monkeypatch):
    old = tmp_path / "site_first_acceptance_previous"
    old.mkdir()
    (old / "report.json").write_text("original")
    new_dir = tmp_path / "site_first_acceptance_repair1"
    monkeypatch.setattr(controller.current, "DATA", tmp_path)
    monkeypatch.setattr(controller.current, "SOURCE_STRUCTURE", old / "input.pdb")
    (old / "input.pdb").write_text("input")
    monkeypatch.setattr(controller.current, "preflight", lambda **kwargs: {"execution_already_started": False})
    monkeypatch.setattr(controller.legacy, "preflight", lambda **kwargs: {"model_identity_paths": {}})
    monkeypatch.setattr(controller.legacy, "_fingerprints", lambda paths: {})
    def run_new(*, acceptance_dir):
        assert acceptance_dir == new_dir
        acceptance_dir.mkdir()
        return {"status": "NEEDS_REVIEW", "errors": ["failure"]}
    monkeypatch.setattr(controller.current, "execute", run_new)
    monkeypatch.setattr(controller.legacy, "execute", lambda **kwargs: (_ for _ in ()).throw(AssertionError("must not run")))
    monkeypatch.setattr(sys, "argv", ["controller", "--execute", "--output-dir", str(new_dir)])
    controller.main()
    assert (old / "report.json").read_text() == "original"
    assert (new_dir / "comparison_report.json").is_file()
