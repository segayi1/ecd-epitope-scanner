from __future__ import annotations

"""Run the frozen score-blind external benchmark through project workers.

The append-only ledger reserves a candidate before launching either predictor,
so an interrupted/failed invocation is never silently retried or double-counted.
"""
import argparse
import contextlib
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "workspace_data" / "p0_p3_upgrade_20261016" / "benchmark"
LEDGER = BASE / "execution_ledger.jsonl"
OUT = BASE / "predictions"
MANIFEST = BASE / "candidate_manifest.json"
RUNNER_MANIFEST = BASE / "runner_execution_manifest.json"
HARD_PAIR_SLOT_CAP = 300
RUNNER_LOCK = BASE / "runner_execution.lock"
FOLDX = ROOT / "integrations" / "workers" / "foldx_worker.py"
SPURS = ROOT / "integrations" / "workers" / "spurs_gpu_worker.py"
BRIDGE = ROOT / "integrations" / "wsl_bridge.py"
NORMALIZED_CHAIN_ID = "A"


def read_candidates():
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    rows = data.get("candidates", data if isinstance(data, list) else [])
    return rows


def append(row):
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with LEDGER.open("a", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(row, sort_keys=True, ensure_ascii=False) + "\n")
        f.flush()
        os.fsync(f.fileno())


@contextlib.contextmanager
def exclusive_runner_lock(path=RUNNER_LOCK):
    """Prevent concurrent resumes from reserving the same score-blind budget."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+b") as handle:
        handle.seek(0, os.SEEK_END)
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        if os.name == "nt":
            import msvcrt
            try:
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError as error:
                raise RuntimeError("External benchmark runner is already active") from error
        else:
            import fcntl
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as error:
                raise RuntimeError("External benchmark runner is already active") from error
        try:
            yield
        finally:
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def events():
    out = []
    if LEDGER.exists():
        for line in LEDGER.read_text(encoding="utf-8").splitlines():
            if line.strip():
                out.append(json.loads(line))
    return out


def run_worker(command, request, timeout, log_path):
    start = time.perf_counter()
    proc = subprocess.run(command, input=json.dumps(request), text=True,
                          capture_output=True, timeout=timeout, cwd=ROOT)
    elapsed = time.perf_counter() - start
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(proc.stderr or "", encoding="utf-8")
    try:
        payload = json.loads(proc.stdout)
    except Exception:
        payload = {"status": "failed", "errors": [{"code": "INVALID_WORKER_JSON",
                    "stdout_tail": proc.stdout[-3000:]}], "records": []}
    if proc.returncode:
        payload = {"status": "failed", "errors": [{"code": "WORKER_EXIT",
                    "returncode": proc.returncode, "stderr_tail": proc.stderr[-3000:]}],
                   "records": []}
    return payload, elapsed, proc.stdout


def mutation_list(row):
    value = row["mutations_normalized"]
    return json.loads(value) if isinstance(value, str) else list(value)


def approved_pair_slot_cap():
    data = json.loads(RUNNER_MANIFEST.read_text(encoding="utf-8"))
    approved = int(data["approved_pair_slots"])
    declared_max = int(data.get("max_pair_slots", HARD_PAIR_SLOT_CAP))
    if approved <= 0 or declared_max <= 0:
        raise ValueError("Benchmark FoldX slot limits must be positive")
    return min(approved, declared_max, HARD_PAIR_SLOT_CAP)


def run_spurs(rows):
    """Batch by exact normalized structure and model order; append every result."""
    existing = events()
    attempted = {str(cid) for e in existing if e.get("event") == "spurs_start"
                 for cid in e.get("candidate_ids", [])}
    rows = [r for r in rows if str(r["source_variant_id"]) not in attempted]
    grouped = defaultdict(list)
    for row in rows:
        # Single and multi SPURS weights are distinct, but the multi worker
        # batches doubles and triples together and reuses one WT encoding.
        grouped[(row["normalized_pdb_sha256"], "single" if int(row["order"]) == 1 else "multi")].append(row)
    distro = "Ubuntu-22.04"
    wsl_python = "/home/segayi1/miniconda3/envs/spurs_v4/bin/python"
    count = 0
    for (struct_sha, mode), items in sorted(grouped.items()):
        order = 1 if mode == "single" else 2
        ids = [str(r["source_variant_id"]) for r in items]
        identity_hash = hashlib.sha256(json.dumps(ids, sort_keys=True).encode()).hexdigest()[:10]
        key = f"spurs:{struct_sha}:{mode}:{identity_hash}"
        first = items[0]
        outdir = OUT / "spurs" / key.replace(":", "_")
        mutations = [mutation_list(r) for r in items]
        request = {"schema_version": "1.0", "request_id": hashlib.sha256(key.encode()).hexdigest()[:24],
                   "structure_path": str((ROOT / first["normalized_pdb"]).resolve()),
                   "chain_id": NORMALIZED_CHAIN_ID, "structure_type": "experimental_pdb_normalized",
                   "sequence": first["sequence"], "mutations": mutations if mode == "multi" else [],
                   "output_dir": str(outdir.resolve()), "device": "cuda",
                   "parameters": {"cuda_oom_retries": 1, "allow_cpu_fallback_on_oom": True,
                                  "reuse_wt_encoding": True}}
        append({"event": "spurs_start", "key": key, "candidate_ids": ids,
                "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "n_candidates": len(items), "orders": sorted({int(r["order"]) for r in items}),
                "structure_sha256": struct_sha, "weight_mode": mode})
        cmd = [sys.executable, str(BRIDGE), "--distro", distro, "--python", wsl_python,
               "--worker", str(SPURS)]
        payload, elapsed, stdout = run_worker(cmd, request, 7200, outdir / "worker.stderr.log")
        (outdir / "worker.stdout.json").write_text(stdout, encoding="utf-8")
        append({"event": "spurs_result", "key": key, "status": payload.get("status"),
                "candidate_ids": ids, "elapsed_seconds": elapsed,
                "payload": payload, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
        count += len(items)
    return count


def run_foldx(rows, limit):
    history = events()
    attempted = {str(e.get("candidate_id")) for e in history if e.get("event") == "foldx_start"}
    remaining = [r for r in rows if str(r["source_variant_id"]) not in attempted]
    remaining = remaining[:limit]
    used = sum(int(e.get("reserved_pair_slots", 5)) for e in history
               if e.get("event") == "foldx_start")
    cap = approved_pair_slot_cap()
    if used + 5 * len(remaining) > cap:
        raise RuntimeError(f"FoldX approved pair-slot budget exceeded: {used}+{5 * len(remaining)}>{cap}")
    for slot, row in enumerate(remaining):
        cid = str(row["source_variant_id"])
        outdir = OUT / "foldx" / f"candidate_{hashlib.sha256(cid.encode()).hexdigest()[:16]}"
        mutations = mutation_list(row)
        append({"event": "foldx_start", "candidate_id": cid, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "reserved_pair_slots": 5, "attempt_index": len(attempted), "order": int(row["order"]),
                "source_structure_sha256": row["source_pdb_sha256"],
                "normalized_structure_sha256": row["normalized_pdb_sha256"]})
        request = {"structure_path": str((ROOT / row["normalized_pdb"]).resolve()),
                   "chain_id": NORMALIZED_CHAIN_ID, "sequence": row["sequence"],
                   "mutations": [mutations], "output_dir": str(outdir.resolve()),
                   "parameters": {"number_of_runs": 5, "max_workers": 1,
                                  "emit_paired_models": True,
                                  "sampling_protocol": "time_rotabase_v3",
                                  "sampling_role": "production",
                                  "sampling_cohort": "external_stability_benchmark_20261016",
                                  "sampling_slot_start": 0}}
        payload, elapsed, stdout = run_worker([sys.executable, str(FOLDX)], request, 21600,
                                               outdir / "worker.stderr.log")
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "worker.stdout.json").write_text(stdout, encoding="utf-8")
        rec = payload.get("records", [{}])[0] if payload.get("records") else {}
        append({"event": "foldx_result", "candidate_id": cid, "status": payload.get("status"),
                "elapsed_seconds": elapsed, "reserved_pair_slots": 5,
                "ddg_kcal_mol": rec.get("ddg"), "ddg_runs_kcal_mol": rec.get("ddg_runs_kcal_mol"),
                "paired_models": rec.get("paired_models"), "provenance": payload.get("provenance"),
                "errors": payload.get("errors", []), "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
        attempted.add(cid)
    return remaining


def run(args):
    rows = read_candidates()
    rows.sort(key=lambda r: (not bool(r.get("timing_pilot")), int(r["order"]), int(r.get("sequence_position") or 0), str(r["source_variant_id"])))
    runner_manifest = json.loads(RUNNER_MANIFEST.read_text(encoding="utf-8"))
    runner_manifest.update({
        "status": "RUNNING",
        "runner_pid": os.getpid(),
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "candidate_count": len(rows),
        "approved_pair_slots": approved_pair_slot_cap(),
        "runner_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    })
    runner_manifest.pop("finished_utc", None)
    runner_manifest.pop("error", None)
    runner_manifest.pop("terminal_receipt_missing", None)
    runner_manifest.pop("last_checked_at_utc", None)
    runner_manifest.pop("runner_pid_observed_absent", None)
    RUNNER_MANIFEST.write_text(json.dumps(runner_manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    limit = 10 if args.pilot else (len(rows) if args.remaining else 0)
    if args.pilot:
        foldx_rows = rows[:10]
    elif args.remaining:
        foldx_done = {str(e.get("candidate_id")) for e in events() if e.get("event") == "foldx_start"}
        foldx_rows = [r for r in rows if str(r["source_variant_id"]) not in foldx_done]
    else:
        foldx_rows = []
    try:
        if not args.skip_spurs:
            spurs_rows = foldx_rows if (args.pilot or args.remaining) else rows
            spurs_count = run_spurs(spurs_rows)
        else:
            spurs_count = 0
        ran = run_foldx(foldx_rows, limit)
        history = events()
        starts = [e for e in history if e.get("event") == "foldx_start"]
        foldx_result_ids = {str(e.get("candidate_id")) for e in history if e.get("event") == "foldx_result"}
        foldx_started_ids = {str(e.get("candidate_id")) for e in starts}
        spurs_started_ids = {str(cid) for e in history if e.get("event") == "spurs_start"
                             for cid in e.get("candidate_ids", [])}
        spurs_scored_ids = {str(cid) for e in history if e.get("event") == "spurs_result"
                            and e.get("status") == "success" for cid in e.get("candidate_ids", [])}
        spurs_failed_ids = {str(cid) for e in history if e.get("event") == "spurs_result"
                            and e.get("status") != "success" for cid in e.get("candidate_ids", [])}
        runner_manifest.update({
            "status": "FINISHED",
            "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "foldx_attempted_candidates": len(foldx_started_ids),
            "foldx_pair_slots_reserved": sum(int(e.get("reserved_pair_slots", 5)) for e in starts),
            "foldx_result_events": sum(e.get("event") == "foldx_result" for e in history),
            "foldx_candidates_without_result": sorted(foldx_started_ids - foldx_result_ids),
            "foldx_failed_result_candidates": sorted({str(e.get("candidate_id")) for e in history
                                                        if e.get("event") == "foldx_result" and e.get("status") != "success"}),
            "spurs_attempted_candidates": len(spurs_started_ids),
            "spurs_scored_candidates": len(spurs_scored_ids),
            "spurs_failed_candidate_ids": sorted(spurs_failed_ids),
        })
        RUNNER_MANIFEST.write_text(json.dumps(runner_manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    except Exception as exc:
        runner_manifest.update({
            "status": "FAILED",
            "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "error": f"{type(exc).__name__}: {exc}",
        })
        RUNNER_MANIFEST.write_text(json.dumps(runner_manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        raise
    print(json.dumps({"spurs_candidate_count_in_new_calls": spurs_count,
                      "foldx_candidates_started": len(ran), "foldx_pair_slots_reserved": 5*len(ran),
                      "foldx_max_candidates": 60, "foldx_max_pair_slots": 300,
                      "approved_pair_slot_cap": approved_pair_slot_cap()}, ensure_ascii=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true", help="run frozen first 10 score-blind records")
    ap.add_argument("--remaining", action="store_true", help="run all not-yet-attempted records")
    ap.add_argument("--skip-spurs", action="store_true")
    args = ap.parse_args()
    with exclusive_runner_lock():
        run(args)


if __name__ == "__main__":
    main()
