from __future__ import annotations

import contextlib
import csv
import json
import os
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

from gpu_runtime import is_cuda_oom, recovery_devices, serialized_gpu, sha256_model_files


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def main() -> int:
    started = time.perf_counter()
    try:
        request = json.load(sys.stdin)
        import torch

        requested = str(request.get("device") or "auto").lower()
        cuda_available = torch.cuda.is_available()
        if requested in {"cuda", "gpu"} and not cuda_available:
            raise RuntimeError("CUDA was requested but is unavailable in the DiscoTope worker")
        if requested in {"cuda", "gpu"} and len(request["sequence"]) >= 1000:
            raise RuntimeError(
                "DiscoTope GPU embedding supports sequences shorter than 1000 residues"
            )
        use_cuda = cuda_available and requested != "cpu" and len(request["sequence"]) < 1000
        parameters = request.get("parameters", {})
        retries = max(0, int(parameters.get("cuda_oom_retries", 1)))
        allow_cpu_fallback = bool(parameters.get("allow_cpu_fallback_on_oom", True))
        root = _project_root()
        vendor = root / "vendor" / "discotope3"
        models = vendor / "models"
        structure_type = "alphafold" if request.get("structure_type") == "alphafold" else "solved"
        devices = recovery_devices(
            "cuda" if use_cuda else "cpu", retries, allow_cpu_fallback
        )
        warnings = []
        oom_retries = 0
        fallback_reason = None
        actual_device = "cuda" if use_cuda else "cpu"
        lock_context = serialized_gpu() if use_cuda else contextlib.nullcontext(0.0)
        with lock_context as gpu_lock_wait_seconds:
            for attempt, attempt_device in enumerate(devices):
                with tempfile.TemporaryDirectory(
                    prefix="discotope_gpu_", dir=request.get("output_dir") or None
                ) as temp_dir:
                    command = [
                        sys.executable,
                        "-m",
                        "discotope3.main",
                        "--pdb_or_zip_file",
                        str(request["structure_path"]),
                        "--struc_type",
                        structure_type,
                        "--out_dir",
                        temp_dir,
                        "--models_dir",
                        str(models),
                        "--verbose",
                        "0",
                    ]
                    if attempt_device == "cpu":
                        command.append("--cpu_only")
                    completed = subprocess.run(
                        command,
                        cwd=vendor,
                        env={**dict(os.environ), "PYTHONPATH": str(vendor)},
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    if completed.stderr:
                        print(completed.stderr, file=sys.stderr, end="")
                    if completed.returncode:
                        failure = completed.stderr + "\n" + completed.stdout[-2000:]
                        if attempt_device == "cuda" and is_cuda_oom(failure):
                            if attempt + 1 < len(devices):
                                if devices[attempt + 1] == "cuda":
                                    oom_retries += 1
                                    warnings.append(f"CUDA_OOM_RETRY_{oom_retries}")
                                    time.sleep(min(5.0, 1.5 * oom_retries))
                                else:
                                    fallback_reason = "CUDA_OUT_OF_MEMORY"
                                    warnings.append("CUDA_OOM_CPU_FALLBACK")
                                continue
                        raise RuntimeError(
                            f"DiscoTope exited with code {completed.returncode}: "
                            f"{completed.stdout[-2000:]}"
                        )
                    csv_files = sorted(Path(temp_dir).rglob("*_discotope3.csv"))
                    if len(csv_files) != 1:
                        raise RuntimeError(f"Expected one DiscoTope CSV, found {len(csv_files)}")
                    records = []
                    with csv_files[0].open(newline="", encoding="utf-8") as handle:
                        for position, row in enumerate(csv.DictReader(handle), 1):
                            records.append(
                                {
                                    "position": position,
                                    "score": float(row["DiscoTope-3.0_score"]),
                                    "calibrated_score": float(row["calibrated_score"]),
                                    "rsasa": float(row["rsa"]),
                                    "epitope": str(row.get("epitope", "")).lower() == "true",
                                }
                            )
                    actual_device = attempt_device
                    break
        if len(records) != len(request["sequence"]):
            raise RuntimeError(
                f"DiscoTope returned {len(records)} residues for sequence length "
                f"{len(request['sequence'])}"
            )
        device = actual_device
        provenance = {
            "worker": "discotope_gpu_worker",
            "device": device,
            "torch_version": torch.__version__,
            "cuda_version": torch.version.cuda,
            "gpu_max_sequence_length": 999,
            "requested_device": requested,
            "actual_device": device,
            "cuda_oom_retry_count": oom_retries,
            "cpu_fallback_reason": fallback_reason,
            "gpu_lock_wait_seconds": gpu_lock_wait_seconds,
        }
        model_paths = list(models.rglob("*")) + [
            root / "models" / "discotope" / "esm_if1_gvp4_t16_142M_UR50.pt"
        ]
        provenance["weights_sha256"], provenance["model_files_sha256"] = (
            sha256_model_files(model_paths, root)
        )
        if device == "cuda":
            provenance["gpu_name"] = torch.cuda.get_device_name(0)
        print(
            json.dumps(
                {
                    "status": "success",
                    "records": records,
                    "warnings": warnings,
                    "errors": [],
                    "provenance": provenance,
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
        return 0
    except Exception as exc:
        traceback.print_exc(file=sys.stderr)
        print(
            json.dumps(
                {
                    "status": "failed",
                    "records": [],
                    "warnings": [],
                    "errors": [{"code": "TOOL_CRASHED", "message": str(exc)}],
                    "provenance": {"worker": "discotope_gpu_worker"},
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
