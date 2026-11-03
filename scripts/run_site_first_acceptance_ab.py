"""Run the authorized single AIB234 case and, if it passes, its cold legacy comparison."""
import argparse
import json
import math
from pathlib import Path

import run_site_first_legacy_baseline as legacy
import run_site_first_single_acceptance as current


def compare(new, old, initial_models, new_identity, old_identity, normalized_sha):
    checks = {
        "new_scientific_acceptance": new.get("status") == "PASS",
        "legacy_completed": old.get("status") == "COMPLETE",
        "same_input": new.get("input_sha256") == old.get("input_sha256") == current.EXPECTED_STRUCTURE_SHA256,
        "same_normalized_structure": bool(normalized_sha) and normalized_sha == (old.get("normalized_structure") or {}).get("sha256"),
        "model_files_unchanged": bool(initial_models) and initial_models == old.get("model_identity"),
        "cold_new_cache": ((new.get("cache_evidence") or {}).get("before_run") or {}).get("cold_verified") is True,
        "cold_legacy_cache": (old.get("cold_cache_evidence") or {}).get("score_cache_entries_at_start") == 0 and (old.get("cold_cache_evidence") or {}).get("all_possible_repair_cache_paths_absent_before_run") is True,
    }
    a = new_identity.get("worker_model_identities") or {}
    b = old_identity.get("worker_model_identities") or {}
    for name in ("spurs_single", "spurs_multi", "discotope"):
        x, y = a.get(name) or {}, b.get(name) or {}
        files = x.get("model_files_sha256")
        checks[name + "_weights_match"] = bool(files) and files == y.get("model_files_sha256")
        checks[name + "_device_match"] = bool(x.get("gpu_name")) and x.get("gpu_name") == y.get("gpu_name") and (x.get("actual_device") or x.get("device")) == (y.get("actual_device") or y.get("device")) == "cuda"
    new_time, old_time = new.get("elapsed_seconds"), old.get("elapsed_pipeline_seconds")
    valid_time = all(isinstance(x, (int, float)) and math.isfinite(x) and x > 0 for x in (new_time, old_time))
    ratio = new_time / old_time if valid_time else None
    checks["timing_available"] = valid_time
    matched = all(checks.values())
    return {"checks": checks, "new_seconds": new_time, "legacy_seconds": old_time,
            "ratio": ratio, "reduction_percent": 100 * (1 - ratio) if ratio is not None else None,
            "performance_target": ("PASS" if ratio <= .5 else "FAIL") if matched else "NOT_ASSESSED"}


def read(path):
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", required=True)
    parser.add_argument("--output-dir", type=Path, default=current.ACCEPTANCE_DIR)
    args = parser.parse_args()
    directory = args.output_dir.resolve()
    if directory.parent != current.DATA.resolve() or not directory.name.startswith("site_first_acceptance_"):
        raise SystemExit("Acceptance output must be a new site_first_acceptance_ directory within workspace_data")
    plan = current.preflight(acceptance_dir=directory)
    if plan["execution_already_started"]:
        raise SystemExit("Existing acceptance state: refusing to start again")
    legacy_plan = legacy.preflight(acceptance_dir=directory)
    models = legacy._fingerprints(legacy_plan["model_identity_paths"])
    new = current.execute(acceptance_dir=directory)
    if new.get("status") == "PASS":
        try:
            current.assert_no_active_jobs(current.DATA / "jobs.sqlite3")
            old = legacy.execute(acceptance_dir=directory)
        except Exception as exc:
            old = {"status": "NOT_STARTED", "errors": [f"Legacy dispatch failed: {type(exc).__name__}: {exc}"]}
    else:
        old = {"status": "NOT_STARTED", "errors": ["New workflow did not pass acceptance"]}
    workdir = Path(new.get("workdir") or ".")
    normalized = workdir / "structure_for_energy.pdb"
    old_output = directory / "baseline" / "legacy_output"
    result = compare(new, old, models, read(workdir / "run_source_identity_v1.json"),
                     read(old_output / "run_source_identity_v1.json"), current.sha256(normalized) if normalized.is_file() else None)
    result["checks"]["original_input_unchanged"] = current.sha256(current.SOURCE_STRUCTURE) == current.EXPECTED_STRUCTURE_SHA256
    if not result["checks"]["original_input_unchanged"]:
        result["performance_target"] = "NOT_ASSESSED"
    result.update(new_status=new.get("status"), legacy_status=old.get("status"), job_id=new.get("job_id"),
                  new_report=str(directory / "report.json"), legacy_report=str(directory / "baseline_report.json"),
                  new_region_principles=new.get("region_principles", []),
                  legacy_region_principles=current._region_principles(old_output, read(old_output / "analysis.json")),
                  errors={"new": new.get("errors", []), "legacy": old.get("errors", [])})
    current.atomic_json(directory / "comparison_report.json", result)
    text = ["# AIB234 单例验收与耗时对照", "", f"新流程验收：{result['new_status']}；旧流程：{result['legacy_status']}。",
            f"新流程墙钟：{result['new_seconds']} 秒；旧流程管线墙钟：{result['legacy_seconds']} 秒；耗时比例：{result['ratio']}。",
            f"50% 提速验收：**{result['performance_target']}**。", "",
            "旧流程为冻结的 pre-adaptive 版本，科学协议与新流程不同；单例不代表其他输入。代表替换限制了搜索空间。", "",
            "## 可比条件核验", "", "```json", json.dumps(result['checks'], indent=2, ensure_ascii=False), "```", "",
            "## 新流程逐区原则", "", "```json", json.dumps(new.get('region_principles', []), indent=2, ensure_ascii=False), "```", "",
            "## 旧流程逐区原则", "", "```json", json.dumps(result['legacy_region_principles'], indent=2, ensure_ascii=False), "```", "",
            "## 错误与不足", "", "```json", json.dumps(result['errors'], indent=2, ensure_ascii=False), "```", ""]
    (directory / "comparison_report.md").write_text("\n".join(text), encoding="utf-8")


if __name__ == "__main__":
    main()
