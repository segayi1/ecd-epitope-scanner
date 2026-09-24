from ecd_epitope_scanner.integrations.workers.gpu_runtime import (
    is_cuda_oom,
    recovery_devices,
)


def test_cuda_oom_detection_accepts_pytorch_messages() -> None:
    assert is_cuda_oom("CUDA out of memory. Tried to allocate 20.00 MiB")
    assert is_cuda_oom(RuntimeError("CUDA error: out of memory"))
    assert not is_cuda_oom("CUDA was requested but is unavailable")


def test_gpu_recovery_plan_retries_then_falls_back_explicitly() -> None:
    assert recovery_devices("cuda", 2, True) == ["cuda", "cuda", "cuda", "cpu"]
    assert recovery_devices("cuda", 1, False) == ["cuda", "cuda"]
    assert recovery_devices("cpu", 9, True) == ["cpu"]
