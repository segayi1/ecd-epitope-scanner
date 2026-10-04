"""Isolated Paraplume 2.0.0 full-model worker with CUDA OOM CPU recovery."""

from __future__ import annotations

import hashlib
import json
import sys
import time
import traceback
from contextlib import redirect_stdout
from importlib.metadata import version
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
EMBEDDING_STAGE_ORDER = ("igbert", "igt5", "esm", "antiberty", "prot-t5")
EMBEDDING_STAGE_LABELS = {
    "igbert": "IgBert",
    "igt5": "IgT5",
    "esm": "ESM2",
    "antiberty": "AntiBERTy",
    "prot-t5": "ProtT5",
}


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def validate_request(request: dict[str, Any]) -> list[dict[str, str]]:
    if request.get("schema_version") != "1.0" or request.get("model_size") != "large":
        raise ValueError("Only schema 1.0 with the full large model is supported")
    if len(str(request.get("model_fingerprint") or "")) != 64:
        raise ValueError("A 64-character model fingerprint is required")
    constructs = request.get("constructs")
    if not isinstance(constructs, list) or not constructs or len(constructs) > 20:
        raise ValueError("constructs must contain between 1 and 20 records")
    result = []
    seen = set()
    for item in constructs:
        construct_id = str(item.get("construct_id") or "")
        vh = str(item.get("vh_sequence") or "").strip().upper()
        vl = str(item.get("vl_sequence") or "").strip().upper()
        if not construct_id or construct_id in seen:
            raise ValueError("Construct IDs must be non-empty and unique")
        if not vh or not vl or set(vh + vl) - AMINO_ACIDS:
            raise ValueError(
                f"Construct {construct_id} VH and VL must use standard amino acids"
            )
        seen.add(construct_id)
        result.append({"construct_id": construct_id, "vh_sequence": vh, "vl_sequence": vl})
    return result


def _weight_provenance() -> dict[str, Any]:
    package_version = version("paraplume")
    if package_version != "2.0.0":
        raise RuntimeError(f"Expected Paraplume 2.0.0, found {package_version}")
    model_dir = Path(sys.prefix) / "Lib" / "site-packages" / "paraplume" / "trained_models" / "large_1_0"
    checkpoint = model_dir / "checkpoint.pt"
    summary_path = model_dir / "summary_dict.json"
    if not checkpoint.is_file() or not summary_path.is_file():
        raise RuntimeError("The installed full Paraplume checkpoint is unavailable")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("embedding_models") not in {"all", "igbert,igT5,esm,antiberty,prot-t5"}:
        raise RuntimeError("The installed Paraplume checkpoint does not declare all five embeddings")
    return {
        "model_version": package_version,
        "model_size": "large",
        "model_checkpoint_sha256": _digest(checkpoint),
        "model_summary_sha256": _digest(summary_path),
        "embedding_models": ["igbert", "igT5", "esm", "antiberty", "prot-t5"],
    }


def _is_cuda_oom(exc: BaseException) -> bool:
    return "out of memory" in str(exc).lower() and "cuda" in str(exc).lower()


def run_inference(
    constructs: list[dict[str, str]], progress_callback=None
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    import pandas as pd
    import paraplume.infer as paraplume_infer
    import torch
    original_embedding_functions = dict(paraplume_infer.llm_to_func)

    def report(progress: float, message: str) -> None:
        if progress_callback is not None:
            progress_callback(progress, message)

    cuda_available = bool(torch.cuda.is_available())
    requested_device = "cuda:0" if cuda_available else "cpu"
    fallback_reason = None if cuda_available else "CUDA_UNAVAILABLE"
    attempts = [0, -1] if cuda_available else [-1]
    records = None
    actual_gpu = -1
    for gpu in attempts:
        try:
            report(0.14, "Loading the full Paraplume model and preparing five embeddings")
            frame = pd.DataFrame(
                [
                    {
                        "construct_id": item["construct_id"],
                        "sequence_heavy": item["vh_sequence"],
                        "sequence_light": item["vl_sequence"],
                    }
                    for item in constructs
                ]
            )
            embedding_functions = original_embedding_functions
            function_names = {
                str(name).casefold(): name for name in embedding_functions
            }
            if len(function_names) != len(EMBEDDING_STAGE_ORDER) or set(function_names) != set(
                EMBEDDING_STAGE_ORDER
            ):
                raise RuntimeError("Paraplume embedding functions do not match the verified five-stage model")
            embedding_count = len(EMBEDDING_STAGE_ORDER)
            progress_start = 0.16
            progress_span = 0.70
            instrumented_functions = {}
            for index, normalized_name in enumerate(EMBEDDING_STAGE_ORDER):
                name = function_names[normalized_name]
                function = embedding_functions[name]
                label = EMBEDDING_STAGE_LABELS[normalized_name]
                start = progress_start + progress_span * index / embedding_count
                end = progress_start + progress_span * (index + 1) / embedding_count

                def with_progress(*args, _function=function, _label=label,
                                  _start=start, _end=end, **kwargs):
                    report(_start, f"Computing {_label} residue embeddings")
                    value = _function(*args, **kwargs)
                    report(_end, f"Completed {_label} residue embeddings")
                    return value

                instrumented_functions[name] = with_progress
            paraplume_infer.llm_to_func = instrumented_functions
            predicted = paraplume_infer.predict_paratope(
                frame,
                gpu=gpu,
                emb_proc_size=1,
                compute_sequence_embeddings=False,
                single_chain=False,
                large=True,
                compute_shap=False,
            )
            records = []
            for item in constructs:
                row = predicted.loc[predicted["construct_id"] == item["construct_id"]]
                if len(row) != 1:
                    raise ValueError("Paraplume returned a missing or duplicate construct row")
                heavy = [float(value) for value in row.iloc[0]["model_prediction_heavy"]]
                light = [float(value) for value in row.iloc[0]["model_prediction_light"]]
                if len(heavy) != len(item["vh_sequence"]) or len(light) != len(item["vl_sequence"]):
                    raise ValueError("Paraplume residue counts do not match the input sequences")
                if any(not 0 <= value <= 1 for value in heavy + light):
                    raise ValueError("Paraplume scores must be finite probabilities in [0, 1]")
                records.append(
                    {
                        "construct_id": item["construct_id"],
                        "vh_length": len(heavy),
                        "vl_length": len(light),
                        "vh_probabilities": heavy,
                        "vl_probabilities": light,
                    }
                )
            actual_gpu = gpu
            break
        except Exception as exc:
            if gpu != 0 or not _is_cuda_oom(exc):
                raise
            fallback_reason = "CUDA_OUT_OF_MEMORY"
            report(0.14, "GPU memory was insufficient; clearing allocations and retrying the full model on CPU")
            print("CUDA out of memory; releasing device allocations and retrying full model on CPU", file=sys.stderr)
            import gc

            gc.collect()
            torch.cuda.empty_cache()
            torch.cuda.ipc_collect()
    if records is None:
        raise RuntimeError("No Paraplume inference attempt completed")
    report(0.90, "Validating residue probabilities and model provenance")
    provenance = _weight_provenance()
    actual_device = "cuda:0" if actual_gpu == 0 else "cpu"
    provenance.update(
        {
            "worker": "paraplume_worker",
            "requested_device": requested_device,
            "actual_device": actual_device,
            "cuda_available": cuda_available,
            "cpu_fallback_reason": fallback_reason,
            "torch_version": torch.__version__,
            "cuda_runtime": torch.version.cuda,
        }
    )
    if actual_device.startswith("cuda"):
        provenance["gpu_name"] = torch.cuda.get_device_name(0)
    return records, provenance


def main() -> int:
    started = time.perf_counter()
    try:
        request = json.load(sys.stdin)
        constructs = validate_request(request)

        def report_progress(progress: float, message: str) -> None:
            stream = sys.__stdout__ or sys.stdout
            stream.write(json.dumps({
                "event": "progress", "progress": progress, "message": message,
            }, separators=(",", ":")) + "\n")
            stream.flush()

        with redirect_stdout(sys.stderr):
            records, provenance = run_inference(constructs, progress_callback=report_progress)
        provenance["model_fingerprint"] = str(request.get("model_fingerprint") or "")
        print(
            json.dumps(
                {
                    "status": "success",
                    "records": records,
                    "warnings": [],
                    "errors": [],
                    "provenance": provenance,
                    "runtime_seconds": time.perf_counter() - started,
                },
                separators=(",", ":"),
                sort_keys=True,
            )
        )
    except Exception as exc:
        traceback.print_exc(file=sys.stderr)
        print(
            json.dumps(
                {
                    "status": "failed",
                    "records": [],
                    "warnings": [],
                    "errors": [{"code": "PARAPLUME_INFERENCE_FAILED", "message": str(exc)}],
                    "provenance": {"worker": "paraplume_worker"},
                    "runtime_seconds": time.perf_counter() - started,
                },
                separators=(",", ":"),
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
