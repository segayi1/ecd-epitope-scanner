from __future__ import annotations

from pathlib import Path
from typing import Any
from uuid import uuid4

from ..joint_design_v4.models import PredictionRequest, PredictionResult, ToolInventory
from .base import ToolError, cache_key, load_cached_result, save_cached_result
from .process_runner import JsonWorkerRunner


class JsonPredictorAdapter:
    name = "external"

    def __init__(
        self,
        *,
        command: list[str] | None = None,
        timeout_seconds: float = 300,
        cache_dir: Path | None = None,
        inventory: ToolInventory | None = None,
    ):
        self.command = list(command or [])
        self.timeout_seconds = timeout_seconds
        self.cache_dir = cache_dir
        self._inventory = inventory or ToolInventory(
            name=self.name,
            available=bool(self.command),
            executable=self.command[0] if self.command else None,
            failure_reason=None if self.command else "not configured",
        )

    def inventory(self) -> ToolInventory:
        return self._inventory

    def validate_input(self, request: PredictionRequest) -> list[dict]:
        errors = []
        if not Path(request.structure_path).is_file():
            errors.append(
                {"code": "UNSUPPORTED_STRUCTURE", "message": "Structure path does not exist"}
            )
        if not request.chain_id:
            errors.append({"code": "CHAIN_MAPPING_AMBIGUOUS", "message": "chain_id is required"})
        if not request.sequence:
            errors.append({"code": "OUTPUT_SCHEMA_INVALID", "message": "sequence is required"})
        return errors

    def validate_result(
        self, request: PredictionRequest, result: PredictionResult
    ) -> list[dict[str, Any]]:
        if not result.ok:
            return []
        if not isinstance(result.records, list) or any(
            not isinstance(row, dict) for row in result.records
        ):
            return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "records must be objects"}]
        return []

    def checked_result(
        self, request: PredictionRequest, result: PredictionResult
    ) -> PredictionResult:
        errors = self.validate_result(request, result)
        if errors:
            return PredictionResult(
                status="failed",
                records=[],
                warnings=list(result.warnings),
                errors=list(result.errors) + errors,
                provenance=dict(result.provenance),
                runtime_seconds=result.runtime_seconds,
            )
        return result

    def predict(self, request: PredictionRequest) -> PredictionResult:
        validation = self.validate_input(request)
        if validation:
            return PredictionResult(status="failed", errors=validation)
        key = None
        # Model digests are returned after execution. A config label or declared
        # digest cannot prove the on-disk files are still the same before lookup.
        cache_identity_known = self.name not in {"spurs", "discotope", "foldx"}
        if self.cache_dir and cache_identity_known:
            key = cache_key(
                Path(request.structure_path),
                chain_id=request.chain_id,
                sequence=request.sequence,
                structure_type=request.structure_type,
                tool_version=self._inventory.version,
                model_revision=self._inventory.model_revision,
                weights_sha256=self._inventory.weights_sha256,
                relevant_parameters={
                    **request.parameters,
                    "_request_device": request.device,
                    "_request_mutations": request.mutations,
                },
            )
            cached = load_cached_result(self.cache_dir, key)
            if cached is not None:
                cached.provenance = dict(cached.provenance) | {"cache_hit": True, "cache_key": key}
                return self.checked_result(request, cached)
        try:
            result = JsonWorkerRunner(self.command, self.timeout_seconds).run(
                request,
                log_path=Path(request.output_dir) / f"{self.name}.log"
                if request.output_dir
                else None,
            )
        except ToolError as exc:
            return PredictionResult(
                status="failed",
                errors=[{"code": exc.code, "message": str(exc), "details": exc.details}],
                provenance={"cache_hit": False},
            )
        result.provenance = dict(result.provenance) | {
            "cache_hit": False,
            **({"cache_disabled_reason": "model_digest_unavailable"}
               if self.cache_dir and not cache_identity_known else {}),
            **({"cache_key": key} if key else {}),
        }
        result = self.checked_result(request, result)
        if key and self.cache_dir and result.ok:
            save_cached_result(self.cache_dir, key, result)
        return result


def request_for(
    adapter: JsonPredictorAdapter,
    structure_path: Path,
    chain_id: str,
    sequence: str,
    output_dir: Path,
    *,
    structure_type: str = "unknown",
    mutations: list[list[str]] | None = None,
    parameters: dict | None = None,
) -> PredictionRequest:
    return PredictionRequest(
        schema_version="1.0",
        request_id=uuid4().hex,
        structure_path=str(structure_path.resolve()),
        chain_id=chain_id,
        structure_type=structure_type,
        sequence=sequence,
        mutations=mutations or [],
        output_dir=str(output_dir.resolve()),
        device=str((parameters or {}).get("device", "auto")),
        parameters=parameters or {},
    )
