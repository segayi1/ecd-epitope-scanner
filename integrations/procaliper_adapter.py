from __future__ import annotations

import math

from ..joint_design_v4.models import PredictionRequest, PredictionResult
from ._adapter import JsonPredictorAdapter


class ProCaliperAdapter(JsonPredictorAdapter):
    name = "procaliper"

    def validate_result(self, request: PredictionRequest, result: PredictionResult) -> list[dict]:
        errors = super().validate_result(request, result)
        if errors or not result.ok:
            return errors
        positions = set()
        for row in result.records:
            try:
                position = int(row.get("position", row.get("residue_index")))
                sasa = float(row.get("sasa", row.get("residue_sasa", row.get("sasa_a2"))))
                rsasa = float(row.get("rsasa", row.get("rSASA", row.get("relative_sasa"))))
            except (TypeError, ValueError):
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid ProCaliper record"}]
            if not (1 <= position <= len(request.sequence)) or not all(
                math.isfinite(value) for value in (sasa, rsasa)
            ):
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid ProCaliper values"}]
            if position in positions:
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Duplicate ProCaliper position"}]
            positions.add(position)
        if positions != set(range(1, len(request.sequence) + 1)):
            return [
                {"code": "RESIDUE_COUNT_MISMATCH", "message": "ProCaliper output is incomplete"}
            ]
        return []
