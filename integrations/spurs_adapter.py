from __future__ import annotations

import math
import re

from ..joint_design_v4.models import PredictionRequest, PredictionResult
from ._adapter import JsonPredictorAdapter


class SPURSAdapter(JsonPredictorAdapter):
    name = "spurs"

    def validate_result(self, request: PredictionRequest, result: PredictionResult) -> list[dict]:
        errors = super().validate_result(request, result)
        if errors or not result.ok:
            return errors
        if not result.records:
            return [{"code": "OUTPUT_MISSING", "message": "SPURS returned no scores"}]
        expected = {
            tuple(sorted(str(item).upper() for item in mutations))
            for mutations in request.mutations
        }
        observed: set[tuple[str, ...]] = set()
        observed_singles: set[tuple[int, str]] = set()
        for row in result.records:
            if request.mutations:
                mutations = row.get("mutations", row.get("mutation_list"))
                if isinstance(mutations, str):
                    mutations = [part.strip() for part in mutations.split(",")]
                if not isinstance(mutations, list) or not mutations:
                    return [
                        {
                            "code": "OUTPUT_SCHEMA_INVALID",
                            "message": "SPURS multi mutations missing",
                        }
                    ]
                if any(
                    not re.fullmatch(r"[A-Z][0-9]+[A-Z]", str(item).upper()) for item in mutations
                ):
                    return [
                        {"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid mutation notation"}
                    ]
                normalized = tuple(sorted(str(item).upper() for item in mutations))
                if normalized not in expected:
                    return [
                        {
                            "code": "OUTPUT_SCHEMA_INVALID",
                            "message": "SPURS returned an unrequested combination",
                        }
                    ]
                if normalized in observed:
                    return [
                        {
                            "code": "OUTPUT_SCHEMA_INVALID",
                            "message": "Duplicate SPURS multi combination",
                        }
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
                value = row.get("ddg", row.get("spurs_multi_ddg_kcal_mol"))
            else:
                try:
                    position = int(row.get("position", row.get("residue_index")))
                except (TypeError, ValueError):
                    return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "SPURS position missing"}]
                if not 1 <= position <= len(request.sequence):
                    return [
                        {"code": "RESIDUE_COUNT_MISMATCH", "message": "SPURS position out of range"}
                    ]
                wt = row.get("wt", row.get("wild_type"))
                if wt is not None and str(wt).upper() != request.sequence[position - 1]:
                    return [
                        {
                            "code": "WT_RESIDUE_MISMATCH",
                            "message": f"SPURS WT mismatch at {position}",
                        }
                    ]
                scores = row.get("scores")
                if isinstance(scores, dict) and scores:
                    expected_mutants = set("ACDEFGHIKLMNPQRSTVWY") - {request.sequence[position - 1]}
                    if set(scores) != expected_mutants:
                        return [{"code": "RESIDUE_COUNT_MISMATCH", "message": "SPURS single matrix is incomplete"}]
                    if any(key[0] == position for key in observed_singles):
                        return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Duplicate SPURS position"}]
                    try:
                        if not all(math.isfinite(float(score)) for score in scores.values()):
                            raise ValueError
                    except (TypeError, ValueError):
                        return [
                            {
                                "code": "OUTPUT_SCHEMA_INVALID",
                                "message": "Invalid SPURS score matrix",
                            }
                        ]
                    observed_singles.update((position, aa) for aa in expected_mutants)
                    continue
                mutant = str(row.get("mutant_aa", row.get("mutant", row.get("aa", "")))).upper()
                if len(mutant) != 1 or mutant == request.sequence[position - 1]:
                    return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "SPURS mutant missing"}]
                value = row.get("ddg", row.get("spurs_ddg_kcal_mol", row.get("score")))
                key = (position, mutant)
                if key in observed_singles:
                    return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Duplicate SPURS mutation"}]
                observed_singles.add(key)
            try:
                if value is None or not math.isfinite(float(value)):
                    raise ValueError
            except (TypeError, ValueError):
                return [{"code": "OUTPUT_SCHEMA_INVALID", "message": "Invalid SPURS ddG"}]
        if request.mutations and observed != expected:
            return [{"code": "OUTPUT_MISSING", "message": "SPURS multi combinations are incomplete"}]
        if not request.mutations:
            complete = {
                (position, aa)
                for position, wt in enumerate(request.sequence, 1)
                for aa in "ACDEFGHIKLMNPQRSTVWY"
                if aa != wt
            }
            if observed_singles != complete:
                return [{"code": "RESIDUE_COUNT_MISMATCH", "message": "SPURS single matrix is incomplete"}]
        return []
