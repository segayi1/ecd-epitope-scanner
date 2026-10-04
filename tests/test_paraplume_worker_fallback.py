from __future__ import annotations

import sys
import types
from typing import Any

from ecd_epitope_scanner.integrations.workers import paraplume_worker


class _Column(list):
    def __eq__(self, other: object) -> list[bool]:
        return [value == other for value in self]


class _Frame:
    def __init__(self, rows: list[dict[str, Any]]):
        self.rows = rows

    @property
    def loc(self):
        return self

    @property
    def iloc(self):
        return self

    def __getitem__(self, key):
        if isinstance(key, str):
            return _Column(row.get(key) for row in self.rows)
        if isinstance(key, list):
            return _Frame([row for row, selected in zip(self.rows, key, strict=True) if selected])
        if isinstance(key, int):
            return self.rows[key]
        raise TypeError(key)

    def __len__(self):
        return len(self.rows)


def _mock_model_runtime(monkeypatch, *, cuda_available, predict_paratope):
    pandas = types.ModuleType("pandas")
    pandas.DataFrame = _Frame

    torch = types.ModuleType("torch")
    torch.__version__ = "test-torch"
    torch.version = types.SimpleNamespace(cuda="test-cuda")
    cache_calls = []
    torch.cuda = types.SimpleNamespace(
        is_available=lambda: cuda_available,
        empty_cache=lambda: cache_calls.append("empty_cache"),
        ipc_collect=lambda: cache_calls.append("ipc_collect"),
    )

    paraplume = types.ModuleType("paraplume")
    paraplume.__path__ = []
    infer = types.ModuleType("paraplume.infer")
    infer.predict_paratope = predict_paratope
    infer.llm_to_func = {
        name: (lambda *args, **kwargs: None)
        # The installed 2.0.0 mapping stores IgT5 before IgBert, while
        # predict_paratope invokes IgBert first. Progress must follow calls.
        for name in ("igT5", "igbert", "esm", "antiberty", "prot-t5")
    }
    monkeypatch.setitem(sys.modules, "pandas", pandas)
    monkeypatch.setitem(sys.modules, "torch", torch)
    monkeypatch.setitem(sys.modules, "paraplume", paraplume)
    monkeypatch.setitem(sys.modules, "paraplume.infer", infer)
    monkeypatch.setattr(
        paraplume_worker,
        "_weight_provenance",
        lambda: {
            "model_version": "2.0.0",
            "model_size": "large",
            "model_checkpoint_sha256": "a" * 64,
        },
    )
    return cache_calls


def _predictions(constructs):
    return _Frame([{
        "construct_id": item["construct_id"],
        "model_prediction_heavy": [0.1] * len(item["sequence_heavy"]),
        "model_prediction_light": [0.2] * len(item["sequence_light"]),
    } for item in constructs.rows])


def test_cuda_oom_retries_full_model_on_cpu_and_records_provenance(monkeypatch):
    calls = []

    def predict(frame, **kwargs):
        calls.append(kwargs)
        assert kwargs["large"] is True
        if kwargs["gpu"] == 0:
            raise RuntimeError("CUDA out of memory. Tried to allocate 1 GiB")
        return _predictions(frame)

    cache_calls = _mock_model_runtime(monkeypatch, cuda_available=True, predict_paratope=predict)
    constructs = [{"construct_id": "c1", "vh_sequence": "ACD", "vl_sequence": "EFG"}]
    records, provenance = paraplume_worker.run_inference(constructs)

    assert [call["gpu"] for call in calls] == [0, -1]
    assert all(call["large"] is True for call in calls)
    assert cache_calls == ["empty_cache", "ipc_collect"]
    assert provenance["requested_device"] == "cuda:0"
    assert provenance["actual_device"] == "cpu"
    assert provenance["cpu_fallback_reason"] == "CUDA_OUT_OF_MEMORY"
    assert records[0]["vh_probabilities"] == [0.1, 0.1, 0.1]
    assert records[0]["vl_probabilities"] == [0.2, 0.2, 0.2]


def test_cuda_unavailable_runs_full_model_on_cpu(monkeypatch):
    calls = []

    def predict(frame, **kwargs):
        calls.append(kwargs)
        assert kwargs["large"] is True
        return _predictions(frame)

    _mock_model_runtime(monkeypatch, cuda_available=False, predict_paratope=predict)
    records, provenance = paraplume_worker.run_inference([
        {"construct_id": "c1", "vh_sequence": "ACD", "vl_sequence": "EFG"}
    ])

    assert [call["gpu"] for call in calls] == [-1]
    assert provenance["requested_device"] == "cpu"
    assert provenance["actual_device"] == "cpu"
    assert provenance["cpu_fallback_reason"] == "CUDA_UNAVAILABLE"
    assert len(records) == 1


def test_full_model_reports_embedding_stages_in_order(monkeypatch):
    calls = []

    def predict(frame, **kwargs):
        infer = sys.modules["paraplume.infer"]
        heavy = [row["sequence_heavy"] for row in frame.rows]
        light = [row["sequence_light"] for row in frame.rows]
        for name in ("igbert", "igT5", "esm", "antiberty", "prot-t5"):
            infer.llm_to_func[name](heavy, light, gpu=kwargs["gpu"], single_chain=False)
        return _predictions(frame)

    _mock_model_runtime(monkeypatch, cuda_available=False, predict_paratope=predict)
    records, _ = paraplume_worker.run_inference(
        [{"construct_id": "c1", "vh_sequence": "ACD", "vl_sequence": "EFG"}],
        progress_callback=lambda progress, message: calls.append((progress, message)),
    )

    completed = [message for _, message in calls if message.startswith("Completed ")]
    assert completed == [
        "Completed IgBert residue embeddings",
        "Completed IgT5 residue embeddings",
        "Completed ESM2 residue embeddings",
        "Completed AntiBERTy residue embeddings",
        "Completed ProtT5 residue embeddings",
    ]
    stage_progress = [
        progress for progress, message in calls
        if message.startswith(("Computing ", "Completed "))
    ]
    assert stage_progress == sorted(stage_progress)
    completed_progress = [
        progress for progress, message in calls if message.startswith("Completed ")
    ]
    assert all(after > before for before, after in zip(completed_progress, completed_progress[1:]))
    progress = [value for value, _ in calls]
    assert progress == sorted(progress)
    assert progress[-1] == 0.90
    assert len(records) == 1
