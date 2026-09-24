from __future__ import annotations

import contextlib
import hashlib
import time
from collections.abc import Iterator
from pathlib import Path


def is_cuda_oom(exc: BaseException | str) -> bool:
    message = str(exc).lower()
    return "cuda out of memory" in message or "cuda error: out of memory" in message


def recovery_devices(
    requested_device: str,
    cuda_oom_retries: int,
    allow_cpu_fallback: bool,
) -> list[str]:
    if requested_device != "cuda":
        return [requested_device]
    devices = ["cuda"] * (max(0, int(cuda_oom_retries)) + 1)
    if allow_cpu_fallback:
        devices.append("cpu")
    return devices


@contextlib.contextmanager
def serialized_gpu(lock_path: str = "/tmp/ecd_v4_gpu.lock") -> Iterator[float]:
    """Serialize WSL GPU workers across tools and Python processes."""
    import fcntl

    path = Path(lock_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with path.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        waited = time.perf_counter() - started
        try:
            yield waited
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def cleanup_cuda(torch) -> None:
    import gc

    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        try:
            torch.cuda.ipc_collect()
        except RuntimeError:
            pass


def sha256_model_files(paths: list[Path], root: Path) -> tuple[str | None, dict[str, str]]:
    """Hash model files and their names so a changed checkpoint changes the run identity."""
    files = sorted({path for path in paths if path.is_file()})
    if not files:
        return None, {}
    combined = hashlib.sha256()
    individual = {}
    for path in files:
        name = path.relative_to(root).as_posix()
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            while chunk := handle.read(8 * 1024 * 1024):
                digest.update(chunk)
        individual[name] = digest.hexdigest()
        combined.update(name.encode("utf-8"))
        combined.update(bytes.fromhex(individual[name]))
    return combined.hexdigest(), individual
