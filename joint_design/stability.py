from __future__ import annotations

import datetime as dt
from pathlib import Path
from typing import Any

from .. import config
from ..stability import EngineInventory, engine_metadata, repair_with_foldx, score_many_with_foldx
from .models import ProbeCandidate


def _mutation_group(candidate: ProbeCandidate) -> list[tuple[str, int, str, str]]:
    return [
        (
            mutation.wild_type,
            int(mutation.residue.auth_seq_id),
            mutation.residue.chain_id,
            mutation.mutant_aa,
        )
        for mutation in candidate.mutations
    ]


def _proxy_annotation(candidates: list[ProbeCandidate], reason: str) -> dict[str, Any]:
    for candidate in candidates:
        candidate.evidence_level = "RULE_ONLY"
        candidate.stability_status = "NOT_RUN"
        candidate.stability_unit = None
        candidate.stability_value = None
        candidate.stability_provenance = {"status": "not_run", "reason": reason}
        candidate.warnings = [
            "No external combination-level energy was run; no kcal/mol value is reported.",
            "Experimental surface expression and conformation-sensitive controls are required.",
        ]
    return {
        "enabled": False,
        "engine": None,
        "evidence_level": "RULE_ONLY",
        "status": "NOT_RUN",
        "candidate_count": len(candidates),
        "submitted_count": 0,
        "scored_count": 0,
        "failed_count": 0,
        "deferred_count": 0,
        "priority_submitted_count": 0,
        "all_failed": False,
        "batch_constructs": 0,
        "reason": reason,
    }


def annotate_stability_evidence(
    candidates: list[ProbeCandidate],
    *,
    inventory: EngineInventory | None = None,
    structure_pdb: Path | None = None,
    workdir: Path | None = None,
    parameters: dict[str, Any] | None = None,
    priority_candidate_ids: list[str] | None = None,
) -> dict[str, Any]:
    """Run combination-level stability for V3's shortlisted candidates."""
    params = parameters or config.JOINT_V3_PARAMETERS.get("stability", {})
    inventory = inventory or EngineInventory(True, None, None)
    if params.get("enabled", True) is False:
        return _proxy_annotation(candidates, "V3_external_stability_disabled_by_configuration")
    if not candidates:
        return _proxy_annotation(candidates, "no_candidates")
    if structure_pdb is None or workdir is None:
        return _proxy_annotation(candidates, "structure_or_workdir_not_provided")
    if not inventory.foldx_path:
        return _proxy_annotation(candidates, "FoldX_not_available")

    binary = inventory.foldx_path
    runs = int(params.get("foldx_runs", 5))
    max_workers = params.get("foldx_max_workers")
    max_candidates = max(1, int(params.get("foldx_max_candidates", len(candidates))))
    priority_id_set = set(priority_candidate_ids or [])
    by_id = {candidate.candidate_id: candidate for candidate in candidates}
    scored_candidates: list[ProbeCandidate] = []
    selected_ids: set[str] = set()

    # Always reserve capacity for the proxy preselection panel (and its
    # alternatives) before filling the remaining budget.  The old first-N
    # truncation could spend the whole FoldX budget on early sector IDs and
    # leave every finally selected probe without an energy value.
    for candidate_id in priority_candidate_ids or []:
        candidate = by_id.get(candidate_id)
        if candidate is None or candidate_id in selected_ids:
            continue
        scored_candidates.append(candidate)
        selected_ids.add(candidate_id)
        if len(scored_candidates) >= max_candidates:
            break

    # Fill remaining slots round-robin by sector so the capped energy batch
    # remains spatially diverse and does not depend on lexical sector order.
    sector_buckets: dict[str, list[ProbeCandidate]] = {}
    for candidate in candidates:
        if candidate.candidate_id not in selected_ids:
            sector_buckets.setdefault(candidate.sector_id, []).append(candidate)
    while len(scored_candidates) < max_candidates:
        changed = False
        for sector_id in sorted(sector_buckets):
            bucket = sector_buckets[sector_id]
            if not bucket or len(scored_candidates) >= max_candidates:
                continue
            candidate = bucket.pop(0)
            scored_candidates.append(candidate)
            selected_ids.add(candidate.candidate_id)
            changed = True
        if not changed:
            break
    deferred_candidates = [
        candidate for candidate in candidates
        if candidate.candidate_id not in selected_ids
    ]
    for candidate in deferred_candidates:
        candidate.evidence_level = "PROXY_ONLY"
        candidate.stability_status = "NOT_RUN_CAP"
        candidate.stability_provenance = {
            "status": "not_run",
            "reason": "foldx_candidate_cap",
            "foldx_max_candidates": max_candidates,
        }
        candidate.warnings = [
            f"FoldX combination scoring deferred after configured cap of {max_candidates} candidates; proxy risk retained."
        ]
    repaired = repair_with_foldx(binary, structure_pdb, workdir)
    common_provenance = {
        "engine": "FoldX",
        "engine_path": str(binary),
        "engine_metadata": engine_metadata(inventory),
        "input_structure": str(structure_pdb),
        "repaired_structure": str(repaired) if repaired else None,
        "command": "RepairPDB + BuildModel",
        "buildmodel_runs": runs,
        "combination_level": True,
        "single_site_sum_used": False,
        "candidate_selection_policy": "proxy_panel_priority_then_sector_round_robin",
        "priority_candidate_ids": list(priority_candidate_ids or []),
        "started_at": dt.datetime.now(dt.UTC).isoformat(),
    }
    if repaired is None:
        for candidate in candidates:
            candidate.evidence_level = "PROXY_ONLY"
            candidate.stability_status = "REPAIR_FAILED"
            candidate.stability_provenance = common_provenance | {"status": "repair_failed"}
            candidate.warnings = [
                "FoldX was detected but RepairPDB did not produce a reusable structure; proxy risk retained.",
                "No kcal/mol value is reported for this candidate.",
            ]
        return {
            "enabled": True,
            "engine": "FoldX",
            "evidence_level": "PROXY_ONLY",
            "status": "REPAIR_FAILED",
            "candidate_count": len(candidates),
            "submitted_count": len(scored_candidates),
            "scored_count": 0,
            "failed_count": len(candidates),
            "deferred_count": len(deferred_candidates),
            "priority_submitted_count": sum(
                candidate.candidate_id in priority_id_set
                for candidate in scored_candidates
            ),
            "all_failed": True,
            "batch_constructs": len(candidates),
            "provenance": common_provenance,
        }

    groups = [_mutation_group(candidate) for candidate in scored_candidates]
    values = score_many_with_foldx(
        binary,
        repaired,
        groups,
        workdir,
        runs=runs,
        already_repaired=True,
        max_workers=int(max_workers) if max_workers is not None else None,
    )
    # The engine adapter normally returns one value per construct.  Pad a
    # truncated response explicitly so an incomplete batch is auditable and
    # cannot silently leave candidates in the default NOT_RUN state.
    values = list(values[: len(scored_candidates)])
    if len(values) < len(scored_candidates):
        values.extend([None] * (len(scored_candidates) - len(values)))
    gate = float(params.get("combination_gate_kcal_mol", config.FOLDX_COMBINATION_GATE))
    scored = failed = 0
    for candidate, value in zip(scored_candidates, values):
        provenance = common_provenance | {
            "status": "scored" if value is not None else "score_failed",
            "combination_gate_kcal_mol": gate,
        }
        candidate.stability_provenance = provenance
        if value is None:
            failed += 1
            candidate.evidence_level = "PROXY_ONLY"
            candidate.stability_status = "SCORE_FAILED"
            candidate.stability_value = None
            candidate.stability_unit = None
            candidate.warnings = [
                "FoldX BuildModel failed or returned incomplete runs; proxy risk retained.",
                "No kcal/mol value is reported for this candidate.",
            ]
            continue
        scored += 1
        candidate.stability_value = float(value)
        candidate.stability_unit = "kcal/mol"
        candidate.stability_status = "ABOVE_GATE" if value > gate else "PASS"
        candidate.evidence_level = "EXTERNAL_COMBINATION"
        candidate.warnings = [] if value <= gate else [
            f"FoldX combination ΔΔG {value:.3f} kcal/mol exceeds configured +{gate:.2f} gate; retained for transparent panel comparison."
        ]
    return {
        "enabled": True,
        "engine": "FoldX",
        "evidence_level": "EXTERNAL_COMBINATION" if scored else "PROXY_ONLY",
        "status": (
            "PARTIAL" if failed else
            "COMPLETE_CAPPED" if deferred_candidates else
            "COMPLETE"
        ),
        "candidate_count": len(candidates),
        "submitted_count": len(scored_candidates),
        "scored_count": scored,
        "failed_count": failed,
        "deferred_count": len(deferred_candidates),
        "priority_submitted_count": sum(
            candidate.candidate_id in priority_id_set
            for candidate in scored_candidates
        ),
        "all_failed": bool(scored_candidates) and failed == len(scored_candidates),
        "batch_constructs": len(scored_candidates),
        "combination_gate_kcal_mol": gate,
        "provenance": common_provenance,
    }
