from __future__ import annotations

import contextlib
import copy
import json
import sys
import time
import traceback
from pathlib import Path

from gpu_runtime import (
    cleanup_cuda,
    is_cuda_oom,
    recovery_devices,
    serialized_gpu,
    sha256_model_files,
    software_identity,
)

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _device(torch, request: dict) -> str:
    requested = str(request.get("device") or "auto").lower()
    if requested in {"cuda", "gpu"}:
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested but is unavailable in the SPURS worker")
        return "cuda"
    if requested == "cpu":
        return "cpu"
    return "cuda" if torch.cuda.is_available() else "cpu"


def _load_multi(model_dir: Path, device: str):
    from omegaconf import OmegaConf
    from spurs.inference import _load_model_state_dict, _seed_everything
    from spurs.models.stability.spurs_multi import SPURSMulti

    cfg = OmegaConf.load(model_dir / ".hydra" / "config.yaml")
    model_cfg = copy.deepcopy(cfg["model"])
    del model_cfg["_target_"]
    _seed_everything(int(cfg["train"]["seed"]))
    model = SPURSMulti(model_cfg).to(device)
    model = _load_model_state_dict(model, model_dir / "checkpoints" / "best.ckpt")
    return model, cfg


def _provenance(torch, device: str, mode: str) -> dict:
    payload = {
        "worker": "spurs_gpu_worker",
        "mode": mode,
        "device": device,
        "torch_version": torch.__version__,
        "cuda_version": torch.version.cuda,
    }
    if device == "cuda":
        payload.update(
            {
                "gpu_name": torch.cuda.get_device_name(0),
                "gpu_memory_allocated_bytes": int(torch.cuda.memory_allocated(0)),
                "gpu_memory_reserved_bytes": int(torch.cuda.memory_reserved(0)),
            }
        )
    return payload


def _predict_once(request: dict, root: Path, torch, device: str) -> tuple[list[dict], str]:
    from spurs.inference import get_SPURS, parse_pdb, parse_pdb_for_mutation

    structure = Path(request["structure_path"])
    chain = str(request["chain_id"])
    mutations = request.get("mutations") or []
    mode = "multi" if mutations else "single"
    model_dir = root / "models" / "spurs" / ("spurs_multi" if mutations else "spurs")
    if mutations:
        model, cfg = _load_multi(model_dir, device)
    else:
        model, cfg = get_SPURS(str(model_dir), device=device)
    pdb = parse_pdb(str(structure), structure.stem, chain, cfg, device=device)
    model.eval()
    with torch.inference_mode():
        if mutations:
            groups: dict[int, list[tuple[int, list[str]]]] = {}
            for index, mutation_list in enumerate(mutations):
                groups.setdefault(len(mutation_list), []).append((index, mutation_list))
            ordered_values: list[float | None] = [None] * len(mutations)
            for mutation_order in sorted(groups):
                indexed_group = groups[mutation_order]
                group_mutations = [item for _, item in indexed_group]
                mut_ids, append_tensors = parse_pdb_for_mutation(group_mutations)
                pdb["mut_ids"] = mut_ids.to(device)
                pdb["append_tensors"] = append_tensors.to(device)
                values = model(pdb).detach().float().cpu().reshape(-1)
                if len(values) != len(indexed_group):
                    raise RuntimeError(
                        f"SPURS multi returned {len(values)} scores for "
                        f"{len(indexed_group)} order-{mutation_order} requests"
                    )
                for (original_index, _), value in zip(
                    indexed_group, values.tolist(), strict=True
                ):
                    ordered_values[original_index] = float(value)
            if any(value is None for value in ordered_values):
                raise RuntimeError("SPURS multi did not score every requested combination")
            records = [
                {"mutations": mutation_list, "ddg": ordered_values[index]}
                for index, mutation_list in enumerate(mutations)
            ]
        else:
            values = model(pdb, return_logist=True).detach().float().cpu().squeeze()
            if values.ndim != 2 or values.shape != (
                len(request["sequence"]),
                len(ALPHABET),
            ):
                raise RuntimeError(
                    f"SPURS single returned shape {tuple(values.shape)}, expected "
                    f"({len(request['sequence'])}, {len(ALPHABET)})"
                )
            records = []
            for index, wt in enumerate(request["sequence"], 1):
                records.append(
                    {
                        "position": index,
                        "wt": wt,
                        "scores": {
                            aa: float(values[index - 1, aa_index])
                            for aa_index, aa in enumerate(ALPHABET)
                            if aa != wt
                        },
                    }
                )
    return records, mode


def main() -> int:
    started = time.perf_counter()
    try:
        request = json.load(sys.stdin)
        root = Path(request.get("parameters", {}).get("model_assets_root") or _project_root())
        sys.path.insert(0, str(root / "vendor" / "spurs"))
        # Some SPURS dependencies print informational messages to stdout during
        # import/inference. Keep stdout reserved for the worker JSON contract.
        with contextlib.redirect_stdout(sys.stderr):
            import torch
            requested_device = _device(torch, request)
            parameters = request.get("parameters", {})
            retries = max(0, int(parameters.get("cuda_oom_retries", 1)))
            allow_cpu_fallback = bool(parameters.get("allow_cpu_fallback_on_oom", True))
            warnings = []
            oom_retries = 0
            fallback_reason = None
            devices = recovery_devices(requested_device, retries, allow_cpu_fallback)
            records = []
            mode = "multi" if request.get("mutations") else "single"
            final_device = requested_device
            lock_context = serialized_gpu() if requested_device == "cuda" else contextlib.nullcontext(0.0)
            with lock_context as gpu_lock_wait_seconds:
                for attempt, device in enumerate(devices):
                    try:
                        records, mode = _predict_once(request, root, torch, device)
                        final_device = device
                        break
                    except Exception as exc:
                        if device != "cuda" or not is_cuda_oom(exc):
                            raise
                        cleanup_cuda(torch)
                        if attempt + 1 >= len(devices):
                            raise
                        if devices[attempt + 1] == "cuda":
                            oom_retries += 1
                            warnings.append(f"CUDA_OOM_RETRY_{oom_retries}")
                            time.sleep(min(5.0, 1.5 * oom_retries))
                        else:
                            fallback_reason = "CUDA_OUT_OF_MEMORY"
                            warnings.append("CUDA_OOM_CPU_FALLBACK")
        checkpoint = root / "models" / "spurs" / ("spurs_multi" if mode == "multi" else "spurs") / "checkpoints" / "best.ckpt"
        esm_weights = root / "models" / "spurs" / "esm2_t33_650M_UR50D.pt"
        weights_sha256, model_files = sha256_model_files([checkpoint, esm_weights], root)
        print(
            json.dumps(
                {
                    "status": "success",
                    "records": records,
                    "warnings": warnings,
                    "errors": [],
                    "provenance": {
                        **_provenance(torch, final_device, mode),
                        **software_identity("spurs", root / "vendor" / "spurs"),
                        "requested_device": requested_device,
                        "actual_device": final_device,
                        "cuda_oom_retry_count": oom_retries,
                        "cpu_fallback_reason": fallback_reason,
                        "gpu_lock_wait_seconds": gpu_lock_wait_seconds,
                        "weights_sha256": weights_sha256,
                        "model_files_sha256": model_files,
                    },
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
                    "provenance": {"worker": "spurs_gpu_worker"},
                    "runtime_seconds": time.perf_counter() - started,
                },
                sort_keys=True,
            )
        )
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
