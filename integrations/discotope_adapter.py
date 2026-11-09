from __future__ import annotations

import math

from ..joint_design_v4.models import PredictionRequest, PredictionResult
from ._adapter import JsonPredictorAdapter


class DiscoTopeAdapter(JsonPredictorAdapter):
    name = "discotope"

    def validate_result(self, request: PredictionRequest, result: PredictionResult) -> list[dict]:
        errors = super().validate_result(request, result)
        if errors or not result.ok:
            return errors
        positions = set()
        for row in result.records:
            try:
                position = int(row.get("position", row.get("residue_index")))
                value = float(
                    row.get(
                        "calibrated_score",
                        row.get("score", row.get("raw_score", row.get("propensity"))),
                    )
                )
            except (TypeError, ValueError):
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid DiscoTope record"}]
            if not 1 <= position <= len(request.sequence) or not math.isfinite(value):
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid DiscoTope values"}]
            if position in positions:
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Duplicate DiscoTope position"}]
            positions.add(position)
        if positions != set(range(1, len(request.sequence) + 1)):
            return [{"code": "RESIDUE_COUNT_MISMATCH", "message": "DiscoTope output is incomplete"}]
        return []
