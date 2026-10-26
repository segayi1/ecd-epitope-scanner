from __future__ import annotations

import math
import re

from ..joint_design_v4.models import PredictionRequest, PredictionResult
from ._adapter import JsonPredictorAdapter


class FoldXAdapter(JsonPredictorAdapter):
    name = "foldx"

    def validate_result(self, request: PredictionRequest, result: PredictionResult) -> list[dict]:
        errors = super().validate_result(request, result)
        if errors or not result.ok:
            return errors
        if not result.records:
            return [{"code": "OUTPUT_MISSING", "message": "FoldX returned no scores"}]
        expected = {
            tuple(sorted(str(item).upper() for item in mutations))
            for mutations in request.mutations
        }
        expected_runs = max(1, int(request.parameters.get("number_of_runs", 5)))
        sampling = request.parameters.get("sampling_protocol") is not None
        slot_start = int(request.parameters.get("sampling_slot_start", 0))
        observed: set[tuple[str, ...]] = set()
        for row in result.records:
            mutations = row.get("mutations", row.get("mutation_list"))
            if isinstance(mutations, str):
                mutations = [part.strip() for part in mutations.split(",")]
            if (
                not isinstance(mutations, list)
                or not mutations
                or any(
                    not re.fullmatch(r"[A-Z][0-9]+[A-Z]", str(item).upper()) for item in mutations
                )
            ):
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid FoldX mutations"}]
            normalized = tuple(sorted(str(item).upper() for item in mutations))
            if normalized not in expected:
                return [
                    {
                        "code": "OUTPUT_SCHEMA_INVALID",
                        "message": "FoldX returned an unrequested combination",
                    }
                ]
            if normalized in observed:
                return [
                    {"code": "OUTPUT_SCHEMA_INVALID", "message": "Duplicate FoldX combination"}
                ]
            observed.add(normalized)
            for notation in normalized:
                match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", notation)
                assert match is not None
                wt, position, _ = match.groups()
                index = int(position)
                if not 1 <= index <= len(request.sequence):
                    return [{"code": "RESIDUE_COUNT_MISMATCH", "message": notation}]
                if request.sequence[index - 1] != wt:
                    return [{"code": "WT_RESIDUE_MISMATCH", "message": notation}]
            try:
                value = float(row.get("ddg", row.get("foldx_ddg_kcal_mol")))
            except (TypeError, ValueError):
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid FoldX ddG"}]
            if not math.isfinite(value):
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid FoldX ddG"}]
            if request.parameters.get("emit_paired_models"):
                pairs = row.get("paired_models")
                if not isinstance(pairs, list) or len(pairs) != expected_runs:
                    return [{"code": "OUTPUT_MISSING", "message": "FoldX raw paired-repeat evidence is incomplete"}]
                pair_values = []
                pair_ids = set()
                slots = []
                for pair in pairs:
                    if not isinstance(pair, dict):
                        return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid FoldX paired-repeat row"}]
                    try:
                        pair_value = float(pair.get("ddg_kcal_mol"))
                    except (TypeError, ValueError):
                        return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid FoldX paired-repeat ddG"}]
                    identity = (pair.get("batch_path"), pair.get("group"), pair.get("run"))
                    if not math.isfinite(pair_value) or None in identity or identity in pair_ids:
                        return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid or duplicate FoldX paired-repeat identity"}]
                    pair_ids.add(identity)
                    pair_values.append(pair_value)
                    if sampling:
                        try:
                            slots.append(int(pair.get("sampling_slot")))
                        except (TypeError, ValueError):
                            return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Missing FoldX logical repeat slot"}]
                if sampling and sorted(slots) != list(range(slot_start, slot_start + expected_runs)):
                    return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "FoldX logical repeat slots do not match the request"}]
                if not math.isclose(float(row.get("ddg")), float(__import__("statistics").median(pair_values)), rel_tol=0, abs_tol=1e-6):
                    return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "FoldX aggregate does not match its raw repeat evidence"}]
        if observed != expected:
            return [{"code": "OUTPUT_MISSING", "message": "FoldX combinations are incomplete"}]
        return []
