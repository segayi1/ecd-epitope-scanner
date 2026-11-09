from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from typing import Protocol

from ..joint_design_v4.models import (
    PredictionRequest,
    PredictionResult,
    ToolInventory,
)


class ToolError(RuntimeError):
    def __init__(self, code: str, message: str, *, details: dict | None = None):
        super().__init__(message)
        self.code = code
        self.details = details or {}


class ExternalPredictor(Protocol):
    name: str

    def inventory(self) -> ToolInventory: ...
    def validate_input(self, request: PredictionRequest) -> list[dict]: ...
    def predict(self, request: PredictionRequest) -> PredictionResult: ...


def cache_key(
    structure_path: Path,
    *,
    chain_id: str,
    sequence: str,
    structure_type: str,
    tool_version: str | None,
    model_revision: str | None,
    weights_sha256: str | None = None,
    relevant_parameters: dict,
) -> str:
    """Return a content-addressed key; filenames never participate."""
    digest = hashlib.sha256()
    digest.update(structure_path.read_bytes())
    canonical = {
        "chain_id": chain_id,
        "sequence": sequence,
        "structure_type": structure_type,
        "tool_version": tool_version,
        "model_revision": model_revision,
        "weights_sha256": weights_sha256,
        "relevant_parameters": relevant_parameters,
    }
    digest.update(json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode())
    return digest.hexdigest()


def load_cached_result(cache_dir: Path, key: str) -> PredictionResult | None:
    path = cache_dir / f"{key}.json"
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        return PredictionResult(**payload)
    except (OSError, TypeError, ValueError):
        return None


def save_cached_result(cache_dir: Path, key: str, result: PredictionResult) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{key}.json"
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(asdict(result), ensure_ascii=False, indent=2), encoding="utf-8")
    tmp.replace(path)
    return path
