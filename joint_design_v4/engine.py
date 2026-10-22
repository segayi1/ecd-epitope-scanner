from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any

from .. import config
from ..stability import AUTO_EXCLUDED_MUTANTS
from ..integrations._adapter import request_for
from ..integrations.discotope_adapter import DiscoTopeAdapter
from ..integrations.foldx_adapter import FoldXAdapter
from ..integrations.procaliper_adapter import ProCaliperAdapter
from ..integrations.spurs_adapter import SPURSAdapter
from ..joint_design.surface_units import build_surface_units
from ..models import (
    DesignRecord,
    MutationChoice,
    ProteinContext,
    RegionRecord,
    ResidueRecord,
    to_dict,
)
from ..sequences import apply_mutations, coding_fragment
from ..util import sha256_file, write_json
from .candidate_search import generate_candidates
from .diagnostics import recommendation_diagnostics
from .contexts import (
    context_manifest,
    normalized_prediction_surface_rows,
    residue_keys,
    validate_mapping,
    write_mapping_sidecar,
    write_normalized_prediction_chain,
)
from .evidence import (
    _spurs_rows,
    attach_epitope_islands,
    build_single_evidence,
    compressed_options,
    epitope_evidence,
    surface_evidence,
)
from .hypotheses import attach_hypothesis_responses, hypotheses_from_surface_v2
from .hierarchy import build_hierarchy
from .qualification import qualification_audit
from .models import SurfaceComplexitySummary
from .geometry_evidence import annotate_hypothesis_geometry
from .models import (
    CombinationStabilityEvidence,
    JointV4Result,
    PredictionRequest,
    PredictionResult,
    ToolInventory,
    V4Candidate,
    V4Panel,
    V4Sector,
)
from .optimizer import select_panel
from .reporting import write_v4_outputs
from .risk import annotate_candidate_risks
from .sectors import PartitionConstraintError, build_epitope_sectors
from .sector_recommendations import build_sector_recommendations
from .surface_snapshot import write_surface_graph_snapshot


def _complete_success(result: PredictionResult) -> bool:
    return result.status == "success" and not result.errors


def _validate_procaliper_input_scope(
    result: PredictionResult,
    structure_path: Path | None,
    chain_id: str,
    sequence: str,
) -> PredictionResult:
    """Fail closed if ProCaliper did not analyze the frozen normalized input."""
    if structure_path is None or not _complete_success(result):
        return result
    provenance = result.provenance
    expected = {
        "input_structure_sha256": sha256_file(structure_path),
        "chain_id": chain_id,
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "residue_count": len(sequence),
    }
    mismatches = [
        key for key, value in expected.items()
        if provenance.get(key) != value
    ]
    if not mismatches:
        return result
    message = "ProCaliper result is not bound to the normalized prediction input: " + ", ".join(mismatches)
    return PredictionResult(
        status="failed",
        records=[],
        warnings=list(result.warnings),
        errors=[*result.errors, {"code": "PROCALIPER_INPUT_SCOPE_UNVERIFIED", "message": message}],
        provenance={**provenance, "scope_validation": "FAILED", "scope_validation_fields": mismatches},
        runtime_seconds=result.runtime_seconds,
    )


def _prediction(
    adapter: Any,
    structure_pdb: Path | None,
    chain_id: str,
    sequence: str,
    output_dir: Path,
    structure_type: str,
    parameters: dict,
    supplied: dict[str, PredictionResult] | None,
    supplied_key: str | None = None,
    mutations: list[list[str]] | None = None,
) -> PredictionResult:
    from ..surface_perturbation.runtime import ACTIVE_BUDGET
    active_budget = ACTIVE_BUDGET.get()
    if active_budget is not None:
        adapter.timeout_seconds = min(adapter.timeout_seconds, active_budget.remaining())
    if mutations:
        # Worker results have one score per physical combination. Repeated
        # display memberships must never duplicate inference requests.
        mutations = [list(signature) for signature in dict.fromkeys(
            tuple(sorted(str(item).strip().upper() for item in group))
            for group in mutations
        )]
    if supplied:
        supplied_keys = (supplied_key,) if supplied_key else (adapter.name,)
        for key in supplied_keys:
            if key and key in supplied:
                request = PredictionRequest(
                    schema_version="1.0",
                    request_id=f"supplied-{key}",
                    structure_path=str(structure_pdb or ""),
                    chain_id=chain_id,
                    structure_type=structure_type,
                    sequence=sequence,
                    mutations=mutations or [],
                    output_dir=str(output_dir),
                    device=str(parameters.get("device", "auto")),
                    parameters=parameters,
                )
                return adapter.checked_result(request, supplied[key])
    if structure_pdb is None:
        return PredictionResult(
            status="failed",
            errors=[{"code": "UNSUPPORTED_STRUCTURE", "message": "No structure path"}],
        )
    request = request_for(
        adapter,
        structure_pdb,
        chain_id,
        sequence,
        output_dir,
        structure_type=structure_type,
        mutations=mutations,
        parameters=parameters,
    )
    if active_budget is not None and adapter.name == "spurs":
        from ..surface_perturbation.upstream_cache import predict_cached
        return predict_cached(adapter,request)
    return adapter.predict(request)


def _map_external_rows(
    rows: list[dict[str, Any]], records: list[ResidueRecord]
) -> list[dict[str, Any]]:
    """Map normalized worker positions back to canonical UniProt positions."""
    normalized = {index: item.uniprot_pos for index, item in enumerate(records, 1)}
    mapped = []
    for row in rows:
        copy_row = dict(row)
        key_name = next(
            (
                name
                for name in ("position", "uniprot_pos", "residue_index", "pos")
                if name in copy_row
            ),
            "position",
        )
        value = copy_row.get(key_name)
        try:
            position = int(value)
        except (TypeError, ValueError):
            mapped.append(copy_row)
            continue
        target = normalized.get(position, position)
        copy_row[key_name] = target
        mapped.append(copy_row)
    return mapped


def _map_mutation_notation(value: str, position_map: dict[int, int]) -> str:
    match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z])", value.strip().upper())
    if not match:
        return value
    wt, position, mutant = match.groups()
    return f"{wt}{position_map.get(int(position), int(position))}{mutant}"


def _physical_confirmation_priority(candidates, maximum, areas, initially_covered=()):
    """Order costly confirmations by proxy diversity; never assign real coverage."""
    remaining=list(candidates);covered=set(initially_covered);output=[]
    while remaining and len(output)<maximum:
        def priority(candidate):
            support=set(candidate.covered_positions_by_scenario.get('geodesic_12',[]))
            gain=sum(areas.get(position,0.) for position in support-covered)
            return (-gain,candidate.structural_risk,-candidate.epitope_score,candidate.candidate_id)
        winner=min(remaining,key=priority)
        output.append(winner);remaining.remove(winner)
        covered.update(winner.covered_positions_by_scenario.get('geodesic_12',[]))
    return output


def _foldx_shortlist(
    candidates: list[V4Candidate],
    maximum: int,
    epitope_priority_per_sector: int = 3,
    stability_fallback_per_sector: int = 3,
    hypotheses=None,
    sectors=None,
    surface_areas=None,
    graph=None,
    diagnostics=None,
) -> list[V4Candidate]:
    """Select a deterministic, sector-balanced FoldX confirmation set."""
    from .spatial_support import spatial_shortlist, site_signature
    maximum = max(1, int(maximum))
    if sectors is not None:
        from dataclasses import replace
        from .spatial_support import candidate_region_membership, regional_support, spatial_priority
        ranked = [c for c in candidates if c.order in {2, 3}
                  and c.search_role != "EXPLORATORY" and c.qualification != "REVIEW_REQUIRED"]
        memberships = {c.candidate_id: candidate_region_membership(c, sectors) for c in ranked}
        local = {}
        for sector in sorted(sectors, key=lambda s: s.sector_id):
            local[sector.sector_id] = {
                c.candidate_id: regional_support(c, sector, surface_areas or {}, graph=graph)
                for c in ranked if sector.sector_id in memberships[c.candidate_id]["region_ids"]}
        output, seen, sites = [], set(), set()

        def add(candidate):
            output.append(candidate)
            seen.add(candidate.candidate_id)
            sites.add(site_signature(candidate))

        def local_key(candidate, sector_id):
            return spatial_priority(replace(candidate, geometry_support=local[sector_id][candidate.candidate_id]))

        # One opportunity per region before any second opportunity.
        for offset in range(2):
            for sector in sorted(sectors, key=lambda s: s.sector_id):
                pool = [c for c in ranked if memberships[c.candidate_id]["exclusive_region_id"] == sector.sector_id
                        and c.candidate_id not in seen and site_signature(c) not in sites]
                if offset == 0:
                    pool = [c for c in pool if c.order == 2] or pool
                if pool and len(output) < maximum:
                    add(min(pool, key=lambda c: local_key(c, sector.sector_id)))
        for offset in range(2):
            for sector in sorted(sectors, key=lambda s: s.sector_id):
                served = {site_signature(c) for c in output if sector.sector_id in memberships[c.candidate_id]["region_ids"]}
                pool = [c for c in ranked if sector.sector_id in memberships[c.candidate_id]["region_ids"]
                        and c.candidate_id not in seen and site_signature(c) not in sites]
                if len(served) < offset + 1 and pool and len(output) < maximum:
                    add(min(pool, key=lambda c: local_key(c, sector.sector_id)))
        for hypothesis in sorted(hypotheses or [], key=lambda h: h.hypothesis_id):
            pool = [c for c in ranked if hypothesis.hypothesis_id in c.targeted_hypotheses
                    and c.candidate_id not in seen and site_signature(c) not in sites]
            if pool and len(output) < maximum:
                add(spatial_shortlist(pool, 1)[0])
        remaining = spatial_shortlist([c for c in ranked if c.candidate_id not in seen], len(ranked))
        remaining.sort(key=lambda c: site_signature(c) in sites)
        for candidate in remaining:
            if len(output) >= maximum:
                break
            add(candidate)
        if diagnostics is not None:
            diagnostics["foldx_region_resources"] = [{
                "sector_id": sector.sector_id, "requested_exclusive_groups": 2,
                "confirmable_exclusive_groups": len({site_signature(c) for c in ranked
                    if memberships[c.candidate_id]["exclusive_region_id"] == sector.sector_id}),
                "scheduled_exclusive_groups": len({site_signature(c) for c in output
                    if memberships[c.candidate_id]["exclusive_region_id"] == sector.sector_id}),
                "scheduled_associated_groups": len({site_signature(c) for c in output
                    if sector.sector_id in memberships[c.candidate_id]["region_ids"]}),
                "scheduled_candidate_ids": [c.candidate_id for c in output
                    if sector.sector_id in memberships[c.candidate_id]["region_ids"]],
            } for sector in sorted(sectors, key=lambda s: s.sector_id)]
        return output
    if hypotheses:
        ranked = sorted((c for c in candidates if c.order < 4 and c.search_role != "EXPLORATORY" and c.qualification != "REVIEW_REQUIRED"),
                        key=lambda c: (c.order > 2, c.structural_risk, -c.epitope_score, c.candidate_id))
        ranked = spatial_shortlist(ranked, len(ranked)) if any(c.geometry_support for c in ranked) else ranked
        output, seen, seen_sites = [], set(), set()
        # Confirmation capacity is a hard worker budget, unlike the expanded
        # generation shortlist. Reserve across hypotheses before global fill.
        for offset in range(2):
            for hypothesis in sorted(hypotheses, key=lambda h: h.hypothesis_id):
                pool = [c for c in ranked if hypothesis.hypothesis_id in c.targeted_hypotheses and c.candidate_id not in seen]
                if any(c.geometry_support for c in pool):
                    pool = [c for c in pool if site_signature(c) not in seen_sites] or pool
                if pool and len(output) < maximum:
                    winner = pool[0]
                    output.append(winner)
                    seen.add(winner.candidate_id)
                    seen_sites.add(site_signature(winner))
        remaining = [c for c in ranked if c.candidate_id not in seen]
        if any(c.geometry_support for c in remaining):
            remaining.sort(key=lambda c: site_signature(c) in seen_sites)
        output.extend(remaining)
        return output[:maximum]
    by_sector: dict[str, list[V4Candidate]] = {}
    for candidate in candidates:
        if candidate.order >= 4 or candidate.search_role == "EXPLORATORY":
            continue
        by_sector.setdefault(candidate.sector_id, []).append(candidate)
    for sector_id, pool in by_sector.items():
        ranked = sorted(
            pool,
            key=lambda item: (
                item.order > 2,
                item.structural_risk,
                -item.epitope_score,
                item.order,
                item.candidate_id,
            )
        )
        if any(c.geometry_support for c in ranked):
            ranked = spatial_shortlist(ranked, len(ranked))
        prioritized = [
            item for item in ranked if item.recommendation_class == "EPITOPE_PRIORITIZED"
        ]
        fallback = [
            item for item in ranked if item.recommendation_class == "STABILITY_FALLBACK"
        ]
        leading = (
            (spatial_shortlist(prioritized, max(0, int(epitope_priority_per_sector)))
             + spatial_shortlist(fallback, max(0, int(stability_fallback_per_sector))))
            if any(c.geometry_support for c in ranked) else
            prioritized[: max(0, int(epitope_priority_per_sector))]
            + fallback[: max(0, int(stability_fallback_per_sector))]
        )
        leading_ids = {item.candidate_id for item in leading}
        by_sector[sector_id] = leading + [
            item for item in ranked if item.candidate_id not in leading_ids
        ]
        if any(c.geometry_support for c in ranked):
            first, variants, sites = [], [], set()
            for candidate in by_sector[sector_id]:
                signature = site_signature(candidate)
                (variants if signature in sites else first).append(candidate)
                sites.add(signature)
            by_sector[sector_id] = first + variants
    output: list[V4Candidate] = []
    offset = 0
    while len(output) < maximum:
        added = False
        for sector_id in sorted(by_sector):
            pool = by_sector[sector_id]
            if offset < len(pool):
                output.append(pool[offset])
                added = True
                if len(output) >= maximum:
                    break
        if not added:
            break
        offset += 1
    return output


def _adapter(
    name: str,
    settings: dict[str, Any],
    cache_dir: Path,
    *,
    tool_name: str | None = None,
    disabled: bool = False,
) -> Any:
    tool_name = tool_name or name
    disabled = disabled or settings.get("enabled") is False
    if disabled:
        command: list[str] = []
    else:
        command = settings.get("command") or []
        env_names = [f"ECD_V4_{name.upper()}_COMMAND"]
        if name.startswith("spurs_"):
            env_names.append("ECD_V4_SPURS_COMMAND")
        env_value = next((os.environ[key] for key in env_names if os.environ.get(key)), None)
        if not command and env_value:
            try:
                value = json.loads(env_value)
                command = [str(item) for item in value] if isinstance(value, list) else []
            except json.JSONDecodeError:
                command = []
        if not command and settings.get("backend") == "wsl_gpu" and os.name == "nt":
            worker_name = {
                "spurs_single": "spurs_gpu_worker.py",
                "spurs_multi": "spurs_gpu_worker.py",
                "discotope": "discotope_gpu_worker.py",
            }.get(name)
            if worker_name:
                root = Path(__file__).resolve().parents[1]
                command = [
                    sys.executable,
                    str(root / "integrations" / "wsl_bridge.py"),
                    "--distro",
                    str(settings.get("wsl_distro", "Ubuntu-22.04")),
                    "--python",
                    str(settings["wsl_python"]),
                    "--worker",
                    str(root / "integrations" / "workers" / worker_name),
                ]
        if not command and settings.get("backend") == "local_worker":
            worker_name = {
                "procaliper": "procaliper_worker.py",
                "foldx": "foldx_worker.py",
            }.get(name)
            if worker_name:
                root = Path(__file__).resolve().parents[1]
                command = [
                    str(settings.get("python") or sys.executable),
                    str(root / "integrations" / "workers" / worker_name),
                ]
    inventory = ToolInventory(
        name=name,
        available=bool(command),
        executable=command[0] if command else None,
        environment_type=(
            "wsl_gpu_worker"
            if settings.get("backend") == "wsl_gpu"
            else "local_json_worker"
            if settings.get("backend") == "local_worker"
            else "local_process"
        ),
        version=settings.get("version"),
        device=settings.get("device"),
        model_revision=settings.get("model_revision"),
        weights_sha256=settings.get("weights_sha256"),
        license_id=settings.get("license_id"),
        failure_reason=(
            "disabled by request" if disabled else None if command else "not configured"
        ),
    )
    klass = {
        "procaliper": ProCaliperAdapter,
        "spurs_single": SPURSAdapter,
        "spurs_multi": SPURSAdapter,
        "discotope": DiscoTopeAdapter,
        "foldx": FoldXAdapter,
    }[tool_name]
    return klass(
        command=command,
        timeout_seconds=float(settings.get("timeout_seconds", 300)),
        cache_dir=cache_dir,
        inventory=inventory,
    )


def _annotate_stability(
    candidates: list[V4Candidate],
    single_rows: list[dict[str, Any]],
    multi: PredictionResult,
    foldx: PredictionResult,
    params: dict[str, Any],
    foldx_params: dict[str, Any],
    position_map: dict[int, int] | None = None,
) -> None:
    multi_map: dict[tuple[str, ...], float] = {}
    for row in multi.records if _complete_success(multi) else []:
        mutations = row.get("mutations") or row.get("mutation_list") or []
        if isinstance(mutations, str):
            mutations = [part.strip() for part in mutations.split(",")]
        if position_map:
            mutations = [_map_mutation_notation(str(item), position_map) for item in mutations]
        value = row.get("ddg", row.get("spurs_multi_ddg_kcal_mol"))
        try:
            if mutations and value is not None and math.isfinite(float(value)):
                multi_map[tuple(sorted(str(item) for item in mutations))] = float(value)
        except (TypeError, ValueError):
            pass
    foldx_map: dict[tuple[str, ...], float] = {}
    for row in foldx.records if _complete_success(foldx) else []:
        mutations = row.get("mutations") or row.get("mutation_list") or []
        if isinstance(mutations, str):
            mutations = [part.strip() for part in mutations.split(",")]
        if position_map:
            mutations = [_map_mutation_notation(str(item), position_map) for item in mutations]
        value = row.get("ddg", row.get("foldx_ddg_kcal_mol"))
        try:
            if mutations and value is not None and math.isfinite(float(value)):
                foldx_map[tuple(sorted(str(item) for item in mutations))] = float(value)
        except (TypeError, ValueError):
            pass
    single_map = _spurs_rows(single_rows)
    for candidate in candidates:
        notation = tuple(sorted(option.notation for option in candidate.mutations))
        spurs = multi_map.get(notation)
        foldx_value = foldx_map.get(notation)
        singles = [
            single_map.get((option.residue.uniprot_pos or -1, option.mutant_aa))
            for option in candidate.mutations
        ]
        single_sum = (
            sum(value for value in singles if value is not None)
            if all(value is not None and math.isfinite(value) for value in singles)
            else None
        )
        disagreement = (
            abs(spurs - foldx_value) if spurs is not None and foldx_value is not None else None
        )
        spurs_pass = spurs is not None and spurs <= float(params.get("preferred_max_ddg", 0.5))
        spurs_review = spurs is not None and spurs <= float(params.get("review_max_ddg", 1.0))
        foldx_pass = foldx_value is not None and foldx_value <= float(
            foldx_params.get("combination_gate_kcal_mol", 2.0)
        )
        disagreement_gate = float(params.get("disagreement_gate_kcal_mol", 1.5))
        if spurs is not None and foldx_value is not None:
            if disagreement is not None and disagreement > disagreement_gate:
                status = "MODEL_DISAGREEMENT"
            elif spurs_pass and foldx_pass:
                status = "CONSENSUS_PASS"
            elif spurs_review and foldx_pass:
                status = "CONSENSUS_REVIEW"
            else:
                status = "COMBINATION_REJECT"
        elif spurs is not None:
            status = "SPURS_ONLY_REVIEW_REQUIRED" if spurs_pass else "COMBINATION_REJECT"
        elif foldx_value is not None:
            status = "FOLDX_ONLY_LEGACY"
        else:
            status = "NOT_SCORED"
        if params.get("stability_contract") == "spurs_primary_foldx_sanity_v3":
            from .stability_v3 import decision
            status = decision(spurs, foldx_value, candidate.order)
        candidate.stability = CombinationStabilityEvidence(
            mutations=notation,
            spurs_multi_ddg_kcal_mol=spurs,
            foldx_ddg_kcal_mol=foldx_value,
            order=candidate.order,
            epistasis_estimate=(
                spurs - single_sum if spurs is not None and single_sum is not None else None
            ),
            extrapolation_level=(
                "WITHIN_PRIMARY_VALIDATION_SCOPE"
                if candidate.order in {1, 2}
                else "HIGHER_ORDER_EXTRAPOLATION_MEDIUM"
                if candidate.order == 3
                else "HIGHER_ORDER_EXTRAPOLATION_HIGH"
            ),
            consensus_status=status,
            disagreement=disagreement,
            provenance={"single_sum_used": single_sum is not None, "stability_contract": params.get("stability_contract", "legacy"),
                        "epistasis_interpretation": "multi minus independently trained single sum; diagnostic only",
                        "models": {"spurs_multi": dict(multi.provenance), "foldx": dict(foldx.provenance)}},
        )


def _formal_candidate_reasons(
    candidate: V4Candidate,
    *,
    allow_spurs_only: bool,
    minimum_coverage_gain: float,
    lower_order_coverage: float,
    higher_order_information_gain: float | None = None,
    regional_coverage_gain: float | None = None,
    regional_gain_required: bool = False,
    physical_coverage_pending: bool = False,
    sasa_conflict_positions: set[int] | None = None,
    sasa_context_sensitive_positions: set[int] | None = None,
    sasa_method_unverified_positions: set[int] | None = None,
    sector_quality_flags: set[str] | None = None,
) -> list[str]:
    reasons = []
    # Recheck the authoritative policy even for options loaded from old caches.
    if any(option.mutant_aa in AUTO_EXCLUDED_MUTANTS for option in candidate.mutations):
        reasons.append("AUTO_EXCLUDED_MUTANT")
    if any(option.residue.one_letter in "CGP" for option in candidate.mutations):
        reasons.append("AUTO_PROTECTED_NATIVE_RESIDUE")
    if candidate.qualification == "REVIEW_REQUIRED" or any(option.review_required for option in candidate.mutations):
        reasons.append("STRUCTURAL_CONFIDENCE_REVIEW_REQUIRED")
    if not candidate.connected:
        reasons.append("DISCONNECTED_CANDIDATE")
    if candidate.order >= 4:
        reasons.append("EXPLORATORY_FOUR_SITE_NOT_PRIMARY")
    # External ProCaliper disagreement remains diagnostic. The authoritative
    # surface gate is the native normalized-context classification.
    if any(
        option.residue.uniprot_pos in (sasa_context_sensitive_positions or set())
        for option in candidate.mutations
    ):
        reasons.append("SASA_STRUCTURE_CONTEXT_CLASSIFICATION_CHANGED")
    # ProCaliper settings/provenance are a diagnostic cross-check. Native
    # context and explicit surface geometry remain the eligibility source.
    if (sector_quality_flags or set()) & {
        "DISCONNECTED_SECTOR",
        "GEODESIC_DISTANCE_UNAVAILABLE",
        "GEODESIC_DIAMETER_EXCEEDED",
        "EUCLIDEAN_DIAMETER_EXCEEDED",
    }:
        reasons.append("SECTOR_GEOMETRY_QC_FAILED")
    if any(
        option.spurs_ddg_kcal_mol is None or option.hard_exclusions
        for option in candidate.mutations
    ):
        reasons.append("SINGLE_EVIDENCE_MISSING_OR_EXCLUDED")
    stability = candidate.stability
    if stability is None:
        reasons.append("COMBINATION_NOT_SCORED")
    elif stability.consensus_status == "CONSENSUS_PASS":
        pass
    else:
        reasons.append(stability.consensus_status)
    # A cached status or the former SPURS-only option cannot replace either model.
    for name, value in (
        ("SPURS", stability.spurs_multi_ddg_kcal_mol if stability else None),
        ("FOLDX", stability.foldx_ddg_kcal_mol if stability else None),
    ):
        if value is None or not math.isfinite(value):
            reasons.append(f"COMBINATION_{name}_MISSING")
    if candidate.order > 2:
        candidate_gain = regional_coverage_gain if regional_gain_required else (
            higher_order_information_gain
            if higher_order_information_gain is not None else
            max(candidate.coverage_by_scenario.values(), default=0.0) - lower_order_coverage
        )
        if not physical_coverage_pending:
            if candidate_gain is None or not math.isfinite(candidate_gain):
                reasons.append("HIGHER_ORDER_COVERAGE_BASELINE_MISSING")
            elif candidate_gain < minimum_coverage_gain:
                reasons.append("INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN")
        if stability is None or stability.foldx_ddg_kcal_mol is None:
            reasons.append("HIGHER_ORDER_REQUIRES_FOLDX")
    return sorted(set(reasons))


def _extra_region_color(index: int) -> str:
    import colorsys

    rgb = colorsys.hsv_to_rgb(((index - 1) * 0.61803398875) % 1.0, 0.65, 0.72)
    return "#" + "".join(f"{round(channel * 255):02X}" for channel in rgb)


def _legacy_regions(sectors: list[V4Sector], records: list[ResidueRecord]) -> list[RegionRecord]:
    by_position = {item.uniprot_pos: item for item in records}
    for record in records:
        record.region_id = None
    output = []
    for index, sector in enumerate(sectors, 1):
        present = [by_position[p] for p in sector.residue_positions if p in by_position]
        for record in present:
            record.region_id = index
        output.append(
            RegionRecord(
                id=index,
                name=f"{present[0].domain if present else 'ECD'} / V4 sector {index}",
                domain=present[0].domain if present else "ECD",
                membrane_position="unknown",
                structural_character="mixed surface",
                residues=list(sector.residue_positions),
                surface_area=sector.total_area_a2,
                geodesic_diameter=sector.geodesic_diameter_a or 0.0,
                interface_fraction=sum(item.interface for item in present) / max(1, len(present)),
                pae_median=None,
                auto_design_allowed=bool(sector.designable_area_a2),
                rationale=(
                    "Connected, diameter-bounded display region; overlapping hypotheses are separate. "
                    f"Geodesic diameter {sector.geodesic_diameter_a or 0.0:.1f} A; "
                    f"CA diameter {sector.euclidean_diameter_a or 0.0:.1f} A."
                ),
                color=(
                    ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#6F4E7C"][index - 1]
                    if index <= 7 else _extra_region_color(index)
                ),
                domain_id=present[0].domain if present else "ECD",
                epitope_prior=sector.epitope_score,
                prior_components={"epitope_score": sector.epitope_score},
                tier="HIGH" if sector.island_ids else "LOW",
                assay_accessible=True,
                predicted_accessibility="unknown",
                required_probe_coverage=1,
                partition_mode="joint_v4_ml",
                partition_basis="explicit_surface_graph",
                envelope_residues=list(sector.residue_positions),
                mutable_residues=[item.uniprot_pos for item in present if not item.hard_exclusions],
                footprint_residues=list(sector.residue_positions),
                surface_connectivity="connected",
                quality_flags=list(sector.quality_flags),
                legacy_projection=True,
                source_engine="joint_v4_ml",
                source_schema_version="4.0",
                semantic_projection=True,
            )
        )
    return output


def _legacy_designs(
    selected: list[V4Candidate],
    sectors: list[V4Sector],
    protein: ProteinContext,
    original_cds: str | None,
    recommendation_status: str,
) -> list[DesignRecord]:
    region_id = {sector.sector_id: index for index, sector in enumerate(sectors, 1)}
    output = []
    for index, candidate in enumerate(selected, 1):
        positions = {option.residue.uniprot_pos for option in candidate.mutations}
        owners = [sector for sector in sectors if positions.intersection(sector.residue_positions)]
        owner = min(owners, key=lambda sector: (-len(positions.intersection(sector.residue_positions)), sector.sector_id)) if owners else None
        patch_id = region_id.get(owner.sector_id) if owner else None
        choices = [
            MutationChoice(
                position=option.residue.uniprot_pos or 0,
                wt=option.residue.one_letter,
                mutant=option.mutant_aa,
                stability_score=round(1.0 - max(0.0, option.spurs_ddg_kcal_mol or 0.0), 4),
                perturbation_score=round(option.chemical_disruption_score, 4),
                stability_percentile=None,
                allowed=True,
                reasons=list(option.warnings),
                stability_class="spurs" if option.spurs_ddg_kcal_mol is not None else "proxy",
                menu_class=option.chemistry_class,
                disruption_score=option.chemical_disruption_score,
                accessibility_score=option.accessibility_score,
                chemical_disruption_score=option.chemical_disruption_score,
                epitope_prior=option.epitope_prior,
                structural_risk_flags=list(option.structural_risk_flags),
            )
            for option in candidate.mutations
        ]
        sequence = apply_mutations(protein.ecd_sequence, choices, protein.ecd_start)
        cds, source = coding_fragment(protein, sequence, original_cds)
        energy = candidate.stability.foldx_ddg_kcal_mol if candidate.stability else None
        saturation_scan = {
            str(option.residue.uniprot_pos): [
                {
                    "position": option.residue.uniprot_pos,
                    "wt": option.residue.one_letter,
                    "mutant": amino_acid,
                    "notation": f"{option.residue.one_letter}{option.residue.uniprot_pos}{amino_acid}",
                    "menu_class": "unscored_saturation_option",
                    "status": "NOT_SCORED",
                }
                for amino_acid in "ACDEFGHIKLMNPQRSTVWY"
                if amino_acid != option.residue.one_letter
            ]
            for option in candidate.mutations
        }
        if recommendation_status == "FORMAL_RECOMMENDATION":
            status = "RECOMMENDED_LOW_PREDICTED_RISK"
        elif recommendation_status == "PREVIEW_NOT_FOR_FORMAL_RECOMMENDATION":
            status = "PREVIEW_CANDIDATE"
        else:
            status = "BLOCKED_EXTERNAL_EVIDENCE"
        output.append(
            DesignRecord(
                id=f"V4_{index:02d}",
                name=f"Patch {patch_id} / {candidate.recommendation_class.lower()}",
                region_id=patch_id,
                mutations=choices,
                status=status,
                stability_engine=(
                    "SPURS + FoldX" if energy is not None else "SPURS/structure evidence"
                ),
                combined_stability_score=energy,
                perturbation_score=round(sum(item.perturbation_score for item in choices), 4),
                rationale=(
                    "DiscoTope-prioritized, stability-gated V4 candidate."
                    if candidate.recommendation_class == "EPITOPE_PRIORITIZED"
                    else "Stability-gated fallback candidate; no DiscoTope-supported site is required."
                ),
                warnings=list(candidate.warnings)
                + [candidate.recommendation_class]
                + [f"Coding fragment source: {source}."]
                + ([candidate.stability.consensus_status] if candidate.stability else [])
                + ["Saturation substitutions are enumerated options; they were not scored by V4."],
                mutant_ecd_sequence=sequence,
                mutant_coding_sequence=cds,
                saturation_scan=saturation_scan,
                score_origin="foldx" if energy is not None else "proxy",
                foldx_diagnostics={
                    **(to_dict(candidate.stability) if candidate.stability else {}),
                    "candidate_id": candidate.candidate_id,
                    "formal_eligible": candidate.formal_eligible,
                    "formal_exclusion_reasons": list(candidate.formal_exclusion_reasons),
                    "recommendation_class": candidate.recommendation_class,
                    "epitope_supported_sites": candidate.epitope_supported_sites,
                    "epitope_prior_score": candidate.epitope_score,
                },
                design_role="STAGE1_PROBE",
                disruption_score=round(sum(item.perturbation_score for item in choices), 4),
                anchor_notation=None,
                source_engine="joint_v4_ml",
                source_schema_version="4.0",
                semantic_projection=True,
            )
        )
    return output


def _sector_recommendation_designs(panel, candidates, primary_designs, sectors, protein, original_cds):
    """Reuse physical constructs across regions; review suggestions remain outside Stage 1."""
    by_candidate = {design.foldx_diagnostics["candidate_id"]: design for design in primary_designs}
    candidates_by_id = {candidate.candidate_id: candidate for candidate in candidates}
    region_ids = {sector.sector_id: index for index, sector in enumerate(sectors, 1)}
    output = []
    rows = [item for primary in panel.sector_recommendations
            for item in [primary, *primary.get("alternatives", [])]]
    for row in rows:
        candidate_id = row["candidate_id"]
        row["design_id"] = None
        if candidate_id is None:
            continue
        if candidate_id not in by_candidate:
            candidate = candidates_by_id[candidate_id]
            design = _legacy_designs([candidate], sectors, protein, original_cds,
                                     "FORMAL_RECOMMENDATION" if candidate.formal_eligible else "BLOCKED_EXTERNAL_EVIDENCE")[0]
            design.id = f"REGION_{candidate_id}"
            design.region_id = region_ids[row["sector_id"]]
            design.name = "Region combination " + "+".join(row["mutation_notations"])
            design.design_role = "SECTOR_RECOMMENDATION" if candidate.formal_eligible else "SECTOR_REVIEW_CANDIDATE"
            design.status = "PASSED_SECTOR_RECOMMENDATION" if candidate.formal_eligible else "NEAR_THRESHOLD_NOT_QUALIFIED"
            design.foldx_diagnostics.update({"failed_gates": row["failed_gates"],
                                            "gate_details": row["gate_details"],
                                            "gate_distance": row["gate_distance"]})
            if not candidate.formal_eligible:
                design.rationale = "Closest available gate comparison for this region; failed gates require review."
                design.warnings.extend(row["failed_gates"])
                for mutation in design.mutations:
                    mutation.allowed = False
                    mutation.reasons.extend(row["failed_gates"])
            by_candidate[candidate_id] = design
            output.append(design)
        row["design_id"] = by_candidate[candidate_id].id
    return output


def _deconvolution_presets(
    primary_designs: list[DesignRecord],
    records: list[ResidueRecord],
    protein: ProteinContext,
    original_cds: str | None,
) -> list[DesignRecord]:
    """Create conditional single and pair controls for selected multi-site probes."""
    by_position = {item.uniprot_pos: item for item in records}
    output: list[DesignRecord] = []
    for parent in primary_designs:
        if parent.design_role != "STAGE1_PROBE" or len(parent.mutations) < 2:
            continue
        choices = sorted(parent.mutations, key=lambda item: item.position)
        for choice in choices:
            sequence = apply_mutations(protein.ecd_sequence, [choice], protein.ecd_start)
            coding, coding_source = coding_fragment(protein, sequence, original_cds)
            output.append(DesignRecord(
                id=f"SING_{parent.id}_{choice.position}",
                name=f"{parent.name} single {choice.notation}",
                region_id=None,
                mutations=[choice],
                status="STAGE2_SINGLE",
                stability_engine=parent.stability_engine,
                combined_stability_score=None,
                perturbation_score=choice.disruption_score,
                rationale=(
                    "Conditional Stage-2 construct to test whether one parent-probe substitution "
                    "reproduces the functional phenotype; it does not establish a direct contact."
                ),
                warnings=[
                    "Test only after the linked Stage-1 probe shows the phenotype.",
                    f"Coding fragment source: {coding_source}.",
                ],
                mutant_ecd_sequence=sequence,
                mutant_coding_sequence=coding,
                score_origin=parent.score_origin,
                design_role="STAGE2_SINGLE",
                parent_design_id=parent.id,
                disruption_score=choice.disruption_score,
                anchor_notation=choice.notation,
                source_engine=parent.source_engine,
                source_schema_version=parent.source_schema_version,
                semantic_projection=parent.semantic_projection,
            ))

        pair_options = []
        for left, right in itertools.combinations(choices, 2):
            left_record = by_position.get(left.position)
            right_record = by_position.get(right.position)
            if left_record is None or right_record is None or left_record.ca is None or right_record.ca is None:
                continue
            distance = math.dist(left_record.ca, right_record.ca)
            if distance >= config.MIN_COMBO_SPACING_A:
                pair_options.append((left, right))
        pair_options.sort(
            key=lambda pair: -(pair[0].disruption_score + pair[1].disruption_score)
        )
        for left, right in pair_options[: config.MAX_PAIRWISE_CONFIRMS]:
            sequence = apply_mutations(
                protein.ecd_sequence, [left, right], protein.ecd_start
            )
            coding, coding_source = coding_fragment(protein, sequence, original_cds)
            output.append(DesignRecord(
                id=f"PAIR_{parent.id}_{left.position}_{right.position}",
                name=f"{parent.name} pair {left.notation}+{right.notation}",
                region_id=None,
                mutations=[left, right],
                status="STAGE2_PAIRWISE",
                stability_engine=parent.stability_engine,
                combined_stability_score=None,
                perturbation_score=left.disruption_score + right.disruption_score,
                rationale=(
                    "Conditional confirmatory subset for combination-only functional loss; "
                    "the result does not establish direct residue-residue contacts."
                ),
                warnings=[
                    "Test only when linked Stage-2 singles do not reproduce the phenotype.",
                    f"Coding fragment source: {coding_source}.",
                ],
                mutant_ecd_sequence=sequence,
                mutant_coding_sequence=coding,
                score_origin=parent.score_origin,
                design_role="STAGE2_PAIRWISE",
                parent_design_id=parent.id,
                disruption_score=left.disruption_score + right.disruption_score,
                source_engine=parent.source_engine,
                source_schema_version=parent.source_schema_version,
                semantic_projection=parent.semantic_projection,
            ))
    return output


def run_joint_v4_ml(
    records: list[ResidueRecord],
    protein: ProteinContext,
    output_dir: Path,
    original_cds: str | None = None,
    chain_id: str | None = None,
    structure_pdb: Path | None = None,
    structure_type: str = "unknown",
    parameters: dict[str, Any] | None = None,
    source_provenance: dict[str, Any] | None = None,
    prediction_results: dict[str, PredictionResult] | None = None,
    surface_geometry: Any | None = None,
    surface_v2_result: Any | None = None,
    membrane_anchor: tuple[float, float, float] | None = None,
    membrane_normal: tuple[float, float, float] | None = None,
    assay_orientation_geometry: dict[str, Any] | None = None,
    confidence: dict[str, Any] | None = None,
    assay_format: str | None = None,
    disable_external_engines: bool = False,
) -> JointV4Result:
    params = {**config.JOINT_V4_ML_PARAMETERS, **(parameters or {})}
    # Nested configuration is merged explicitly to keep caller overrides small.
    for key in config.JOINT_V4_ML_PARAMETERS:
        if isinstance(config.JOINT_V4_ML_PARAMETERS[key], dict) and isinstance(
            (parameters or {}).get(key), dict
        ):
            params[key] = {**config.JOINT_V4_ML_PARAMETERS[key], **(parameters or {})[key]}
    supplied_predictions = dict(prediction_results or {})
    from ..surface_perturbation import STRATEGY_ID
    physical_preview = params.get("surface_strategy") == STRATEGY_ID
    from .coverage_profile import PROFILE_V3, coverage_enabled, sequence_reasons
    coverage_preview = coverage_enabled(params.get("optimization_profile"))
    if coverage_preview and not physical_preview:
        raise ValueError("Coverage profile requires physical surface strategy")
    if params.get("surface_strategy") not in {None, "legacy", STRATEGY_ID}:
        raise ValueError("Unknown surface_strategy")
    if physical_preview:
        # Physical preview retains its frozen display policy independently of V4 defaults.
        params["sectoring"] = {"policy": "connected_bounded_v3",
            "maximum_geodesic_diameter_a": 26.0, "maximum_euclidean_diameter_a": 18.0,
            "minimum_sector_residues": 4, "minimum_sector_area_a2": 50.0,
            **(parameters or {}).get("sectoring", {})}
        from ..surface_perturbation.combo_contract import CONTRACT
        params['surface_perturbation']={**params.get('surface_perturbation',{}),'panel_contract':CONTRACT}
        params['surface_perturbation'].setdefault('candidate_route','unpartitioned_local')
        params["candidate_search"] = {**params["candidate_search"], "orders": [2, 3],
            "primary_orders":[2,3],"exploratory_orders":[],"global_span_paths":True}
        params["foldx_confirmation"] = {**params["foldx_confirmation"], "number_of_runs": 1, "emit_paired_models": True}
        params["foldx_confirmation"]["paired_cache_roots"] = params.get("surface_perturbation", {}).get("paired_cache_roots", [])
    if coverage_preview:
        params["candidate_search"] = {**params["candidate_search"], "coverage_scan_v2": True,
            "canonical_sequence": protein.sequence, "orders": [2], "exploratory_orders": [],
            "mandatory_candidate_ids": params.get('candidate_ids') or []}
    params['surface_perturbation'] = {**params.get('surface_perturbation', {}),
        'export_mode': params.get('export_mode', 'deliverable'),
        'optimization_profile': params.get('optimization_profile', 'legacy')}
    if params.get('calibration_repeat_schedule') is not None:
        params['surface_perturbation']['calibration_repeat_schedule'] = params['calibration_repeat_schedule']
    if params.get("optimization_profile") == PROFILE_V3:
        params["spurs_multi"]["stability_contract"] = "spurs_primary_foldx_sanity_v3"
        params["foldx_confirmation"].update(sampling_protocol="time_rotabase_v3", sampling_role="production", sampling_cohort="production_v3")
    if disable_external_engines:
        for name in ("procaliper", "spurs", "spurs_multi", "discotope", "foldx"):
            supplied_predictions.setdefault(name, PredictionResult(
                status="failed",
                errors=[{
                    "code": "ENGINE_DISABLED",
                    "message": "External engine execution was disabled by the request.",
                }],
                provenance={"disabled_by_request": True},
            ))
        prediction_results = supplied_predictions
    output_dir.mkdir(parents=True, exist_ok=True)
    if coverage_preview:
        params['foldx_confirmation'].setdefault('engine_work_root',str(output_dir/'fx'))
        params['foldx_confirmation'].setdefault('short_binary_alias',True)
        assets=Path(params.get('spurs_single',{}).get('model_assets_root') or Path(__file__).resolve().parents[1])
        binaries=sorted((assets/'vendor/foldx').glob('foldx_*.exe'))
        if binaries:params['foldx_confirmation'].setdefault('binary',str(binaries[-1]))
    mode = str(params.get("mode", "strict")).lower()
    if mode not in {"strict", "preview"}:
        raise ValueError("mode must be 'strict' or 'preview'")
    chain = chain_id or (records[0].chain_id if records else "")
    position_map = {index: item.uniprot_pos for index, item in enumerate(records, 1)}
    mapping_errors = validate_mapping(records, protein)
    if mapping_errors:
        raise ValueError("WT_RESIDUE_MISMATCH")
    keys = residue_keys(records)
    write_mapping_sidecar(output_dir / "residue_mapping.csv", records, keys)
    normalized_path = write_normalized_prediction_chain(
        structure_pdb, output_dir / "prediction_context.pdb", records, chain
    )
    prediction_structure = normalized_path or structure_pdb
    prediction_sequence = "".join(item.aa for item in records)
    procaliper_parameters = dict(params.get("procaliper", {}))
    procaliper_parameters["_normalized_context_contract_version"] = 2
    cache_dir = Path(params.get("runtime", {}).get("cache_dir") or (output_dir / "cache"))
    adapters = {
        "procaliper": _adapter(
            "procaliper", params.get("procaliper", {}), cache_dir,
            disabled=disable_external_engines,
        ),
        "spurs_single": _adapter(
            "spurs_single", params.get("spurs_single", {}), cache_dir,
            tool_name="spurs_single", disabled=disable_external_engines,
        ),
        "spurs_multi": _adapter(
            "spurs_multi", params.get("spurs_multi", {}), cache_dir,
            tool_name="spurs_multi", disabled=disable_external_engines,
        ),
        "discotope": _adapter(
            "discotope", params.get("discotope", {}), cache_dir,
            disabled=disable_external_engines,
        ),
        "foldx": _adapter(
            "foldx", params.get("foldx_confirmation", {}), cache_dir,
            disabled=disable_external_engines,
        ),
    }
    inventories = {name: adapter.inventory() for name, adapter in adapters.items()}
    graph_provenance: dict[str, Any] = {}
    native_units, graph, areas = build_surface_units(
        records,
        rsasa_threshold=float(params["surface"].get("rsasa_threshold", config.RSASA_THRESHOLD)),
        neighbor_cutoff_a=10.0,
        max_neighbors=8,
        surface_geometry=surface_geometry,
        graph_provenance=graph_provenance,
    )
    components = {unit.position: unit.component_id for unit in native_units}
    procaliper = _prediction(
        adapters["procaliper"],
        prediction_structure,
        chain,
        prediction_sequence,
        output_dir,
        structure_type,
        procaliper_parameters,
        prediction_results,
    )
    procaliper = _validate_procaliper_input_scope(
        procaliper, normalized_path, chain, prediction_sequence
    )
    spurs = _prediction(
        adapters["spurs_single"],
        prediction_structure,
        chain,
        prediction_sequence,
        output_dir,
        structure_type,
        params.get("spurs_single", {}),
        prediction_results,
        "spurs",
    )
    discotope = _prediction(
        adapters["discotope"],
        prediction_structure,
        chain,
        prediction_sequence,
        output_dir,
        structure_type,
        params.get("discotope", {}),
        prediction_results,
    )
    procaliper.records = _map_external_rows(procaliper.records, records)
    spurs.records = _map_external_rows(spurs.records, records)
    discotope.records = _map_external_rows(discotope.records, records)
    normalized_context_rows = None
    normalized_context_method = None
    if normalized_path is not None:
        normalized_context_rows, normalized_context_method = normalized_prediction_surface_rows(
            normalized_path, chain, records
        )
    native_rows = [
        {"position": item.position, "sasa": item.mapped_area_a2, "rsasa": item.rasa}
        for item in native_units
    ]
    surface = surface_evidence(
        records,
        keys,
        procaliper.records if _complete_success(procaliper) else [],
        areas,
        components,
        normalized_context_rows=normalized_context_rows,
        normalized_context_method=normalized_context_method,
        procaliper_provenance=procaliper.provenance,
        rsasa_threshold=float(params["surface"].get("rsasa_threshold", config.RSASA_THRESHOLD)),
    )
    epitope = epitope_evidence(
        records, keys, discotope.records if _complete_success(discotope) else [], graph
    )
    attach_epitope_islands(
        epitope,
        graph,
        float(params["discotope"].get("seed_percentile", 0.80)),
        float(params["discotope"].get("local_peak_percentile", 0.65)),
    )
    single = build_single_evidence(
        records,
        keys,
        spurs.records if _complete_success(spurs) or mode == "preview" else [],
        params["spurs_single"],
        allow_unscored=mode == "preview",
        surface_by_position={item.residue.uniprot_pos: item for item in surface},
        sequence=protein.sequence,
        coverage_profile=coverage_preview,
    )
    options = compressed_options(
        records,
        single,
        int(params["spurs_single"].get("max_options_per_site", 3)),
        epitope,
        float(params["candidate_search"].get("epitope_priority_weight", 1.0)),
        float(params["discotope"].get("seed_percentile", 0.80)),
        coverage_profile=coverage_preview,
    )
    if coverage_preview:
        if params.get('full_replacement_search',False):
            existing={(o.residue.uniprot_pos,o.mutant_aa) for o in options}
            for row in single:
                key=(row.residue.uniprot_pos,row.mutant_aa)
                if row.eligible and not row.review_required and key not in existing:
                    previous=(row.menu_selected,row.pareto_layer,row.selection_reason)
                    options.extend(compressed_options(records,[row],1,coverage_profile=True))
                    row.menu_selected,row.pareto_layer,row.selection_reason=previous
                    existing.add(key)
        for supplement in params.get('menu_supplements', []):
            matches=[row for row in single if row.residue.uniprot_pos==supplement['position']
                and row.mutant_aa==supplement['mutant'] and row.eligible and not row.review_required]
            if len(matches)!=1 or not supplement.get('reason'):
                details=[{'eligible':row.eligible,'review_required':row.review_required,
                    'exclusion_reasons':row.hard_exclusions} for row in single
                    if row.residue.uniprot_pos==supplement['position'] and row.mutant_aa==supplement['mutant']]
                raise ValueError(f'Ineligible or unexplained menu supplement {supplement}: '
                    f'{details}; SPURS status={spurs.status}, errors={spurs.errors}')
            row=matches[0]
            if not any(o.residue.uniprot_pos==supplement['position'] and o.mutant_aa==supplement['mutant'] for o in options):
                original_layer=row.pareto_layer
                options.extend(compressed_options(records,[row],1,coverage_profile=True))
                row.pareto_layer=original_layer
                row.selection_reason='SUPPLEMENT: '+supplement['reason']
    record_map = {item.uniprot_pos: item for item in records}
    partition_diagnostics: dict[str, Any] = {}
    graph_provenance["display_partition_policy"] = partition_diagnostics
    try:
        sectors = build_epitope_sectors(
            surface,
            epitope,
            graph,
            record_map,
            max(1, int(params.get("max_sectors", 7))),
            params.get("sectoring", {}),
            diagnostics=partition_diagnostics,
        )
    except PartitionConstraintError as error:
        write_json(output_dir / "partition_diagnostics_v4.json", error.diagnostics)
        # Save the original eligible graph, without certifying draft ownership.
        write_surface_graph_snapshot(
            output_dir / "surface_graph_v4.json", graph, areas, None,
            mapping_path=output_dir / "residue_mapping.csv",
            structure_path=structure_pdb, provenance=graph_provenance,
            include_display_ownership=False,
        )
        raise
    write_json(output_dir / "partition_diagnostics_v4.json", partition_diagnostics)
    graph_snapshot_path = output_dir / "surface_graph_v4.json"
    graph_snapshot_sha256 = write_surface_graph_snapshot(
        graph_snapshot_path,
        graph,
        areas,
        sectors,
        mapping_path=output_dir / "residue_mapping.csv",
        structure_path=structure_pdb,
        provenance=graph_provenance,
    )
    graph_provenance["snapshot_file"] = graph_snapshot_path.name
    graph_provenance["snapshot_sha256"] = graph_snapshot_sha256
    if surface_v2_result is not None:
        write_surface_graph_snapshot(
            output_dir / "surface_graph_v41.json", graph, areas, sectors,
            mapping_path=output_dir / "residue_mapping.csv",
            structure_path=structure_pdb,
            provenance=graph_provenance,
            include_display_ownership=False,
        )
        write_json(
            output_dir / "surface_display_regions_v41.json",
            {
                "schema_version": "4.1",
                "interpretation": "Display ownership only; scientific footprint memberships may overlap.",
                "regions": to_dict(sectors),
                "partition_policy": graph_provenance["display_partition_policy"],
                "surface_v2_regions": to_dict(surface_v2_result.surface_regions),
            },
        )
    raw_hypotheses, families = [], []
    hypotheses = (
        hypotheses_from_surface_v2(surface_v2_result, records, epitope)
        if bool(params.get("use_surface_hypotheses", True)) else []
    )
    if surface_v2_result is not None and params.get("use_hypothesis_hierarchy", False):
        raw_hypotheses, families, hypotheses = build_hierarchy(
            surface_v2_result, records, epitope, graph, params.get("hypothesis_hierarchy", {}),
            confidence=confidence, predicted=structure_type in {"alphafold", "predicted"},
        )
    candidate_pool_diagnostics = {}
    search_hypotheses=hypotheses or None
    if physical_preview:
        from ..surface_perturbation.routes import candidate_hypotheses
        search_hypotheses=candidate_hypotheses(params.get('surface_perturbation',{}).get(
            'candidate_route','overlapping_windows'),options,graph,params['candidate_search'],hypotheses)
    candidates = generate_candidates(
        sectors, [option for option in options if not option.review_required], graph,
        {**params["candidate_search"], "sector_resource_guarantees": not physical_preview},
        hypotheses=search_hypotheses,
        diagnostics=candidate_pool_diagnostics,
        surface_areas=areas if not physical_preview else None,
    )
    review_diagnostics = {}
    if any(option.review_required for option in options):
        candidates.extend(generate_candidates(
            sectors, options, graph, {**params["candidate_search"], "require_review": True,
                                    "sector_resource_guarantees": False},
            hypotheses=hypotheses or None, diagnostics=review_diagnostics,
            surface_areas=areas if not physical_preview else None,
        ))
    priority_ids=[c.candidate_id for c in candidates]
    if coverage_preview:
        backup_ids={row['candidate_id'] for row in params['candidate_search'].get('conditional_pair_backups',[])}
        for candidate in candidates:
            if candidate.candidate_id in backup_ids:candidate.search_role='STABILITY_FAILURE_BACKUP'
    if physical_preview and (not coverage_preview or any(c.order==3 for c in candidates)):
        from .candidate_search import add_pair_controls
        candidates=add_pair_controls(candidates,graph)
    if coverage_preview:
        candidates = [c for c in candidates if not sequence_reasons(protein.sequence,
            [(o.residue.uniprot_pos, o.residue.one_letter, o.mutant_aa) for o in c.mutations])]
        # Controls and triples share the same hard scoring-pool capacity.
        required=set(params.get('candidate_ids') or [])
        for triple in (c for c in candidates if c.order==3):
            mutations={o.notation for o in triple.mutations}
            required.update(c.candidate_id for c in candidates
                if c is triple or (c.order==2 and {o.notation for o in c.mutations}<mutations))
        cap=int(params['candidate_search'].get('global_shortlist',32))
        if len(required)>cap:raise ValueError('Required nodes exceed hard scoring pool capacity')
        ordinal={cid:i for i,cid in enumerate(priority_ids)}
        candidates=sorted(candidates,key=lambda c:(c.candidate_id not in required,
            ordinal.get(c.candidate_id,len(ordinal)),c.candidate_id))[:cap]
        from ..surface_perturbation.paired import atoms
        atom_rows=atoms(structure_pdb);heavy_by_residue={}
        for atom_key,value in atom_rows.items():
            if value[2].upper()!='H':heavy_by_residue.setdefault(atom_key[:3],[]).append(value[0])
        contact_cache={}
        for candidate in candidates:
            separations=[];contacts=[]
            for left,right in itertools.combinations(candidate.mutations,2):
                a,b=record_map[left.residue.uniprot_pos].ca,record_map[right.residue.uniprot_pos].ca
                if a is not None and b is not None:
                    separations.append(math.dist(a,b))
                keys=tuple(sorted((left.residue.stable_id,right.residue.stable_id)))
                if keys not in contact_cache:
                    sites=[heavy_by_residue.get((o.residue.chain_id,o.residue.auth_seq_id,
                        (o.residue.insertion_code or '').strip()),[]) for o in (left,right)]
                    contact_cache[keys]=min((math.dist(a,b) for a in sites[0] for b in sites[1]),default=None)
                contacts.append(contact_cache[keys])
            candidate.risk_components['minimum_ca_distance_a']=min(separations) if separations else None
            candidate.risk_components['maximum_ca_distance_a']=max(separations) if separations else None
            candidate.risk_components['minimum_native_heavy_atom_distance_a']=min(
                (d for d in contacts if d is not None),default=None)
            candidate.risk_components['native_heavy_contact_pair_count_4p5_a']=sum(
                d<=4.5 for d in contacts if d is not None) if all(d is not None for d in contacts) else None
            # C-alpha separation is geometry context, never a contact/stability gate.
        write_json(output_dir/'candidate_pool_v2.json', {'schema_version':'coverage_pool_2',
            'candidates':to_dict(candidates),'diagnostics':candidate_pool_diagnostics})
    if params.get('candidate_stage_limit') is not None and params.get('candidate_ids') is None:
        limit = int(params['candidate_stage_limit'])
        if limit < 0:
            raise ValueError('Negative candidate stage limit')
        candidates = candidates[:limit]
    if params.get("candidate_ids") is not None:
        whitelist = set(params["candidate_ids"])
        missing = whitelist - {c.candidate_id for c in candidates}
        if missing:
            raise ValueError("Unknown candidate_ids: " + ", ".join(sorted(missing)))
        candidates = [c for c in candidates if c.candidate_id in whitelist]
    if params.get('optimization_profile') == PROFILE_V3:
        for candidate in candidates:
            if candidate.order==3:
                candidate.triple_role=params['candidate_search'].get('triple_roles_by_candidate_id',{}).get(candidate.candidate_id,candidate.triple_role)
                candidate.parent_candidate_ids=params['candidate_search'].get('triple_parents_by_candidate_id',{}).get(candidate.candidate_id,candidate.parent_candidate_ids)
    annotate_hypothesis_geometry(
        hypotheses + raw_hypotheses, records,
        membrane_anchor=membrane_anchor,
        membrane_normal=membrane_normal,
        assay_orientation_geometry=assay_orientation_geometry,
        confidence=confidence,
        predicted=structure_type in {"alphafold", "predicted"},
    )
    response_radius = float(params["panel_optimizer"].get("hypothesis_response_radius_a", 15.0))
    response_decay = float(params["panel_optimizer"].get("hypothesis_response_decay_a", 6.0))
    attach_hypothesis_responses(
        candidates, hypotheses, graph,
        radius_a=response_radius, decay_length_a=response_decay,
    )
    if not physical_preview:
        from .spatial_support import annotate_higher_order_regional_gain
        annotate_higher_order_regional_gain(candidates, sectors, areas)
    normalized_positions = {item.uniprot_pos: index for index, item in enumerate(records, 1)}

    def candidate_notation(candidate: V4Candidate) -> list[str]:
        return [
            f"{option.residue.one_letter}{normalized_positions[option.residue.uniprot_pos]}"
            f"{option.mutant_aa}"
            for option in candidate.mutations
        ]

    multi = (
        _prediction(
            adapters["spurs_multi"],
            prediction_structure,
            chain,
            prediction_sequence,
            output_dir,
            structure_type,
            params.get("spurs_multi", {}),
            prediction_results,
            "spurs_multi",
            [candidate_notation(candidate) for candidate in candidates if not physical_preview or candidate.order > 1],
        )
        if candidates
        else PredictionResult(status="failed", errors=[{
            "code": "NO_CANDIDATES_TO_SCORE",
            "message": "SPURS multi was not run because candidate generation returned no combinations.",
        }])
    )
    if physical_preview and _complete_success(spurs):
        # Singles reuse the already computed stability landscape, with no multi forward pass.
        singles = _spurs_rows(spurs.records)
        reused = [{"mutations": candidate_notation(c), "ddg": singles[(c.mutations[0].residue.uniprot_pos,c.mutations[0].mutant_aa)],
                   "source": "reused_spurs_single"} for c in candidates if c.order == 1
                  and (c.mutations[0].residue.uniprot_pos,c.mutations[0].mutant_aa) in singles]
        if reused:
            multi = PredictionResult(status="success" if _complete_success(multi) or all(c.order == 1 for c in candidates) else multi.status,
                                     records=multi.records+reused, errors=multi.errors, warnings=multi.warnings, provenance=multi.provenance)
    _annotate_stability(
        candidates, spurs.records if _complete_success(spurs) else [],
        multi, PredictionResult(status="failed"),
        params["spurs_multi"], params["foldx_confirmation"], position_map,
    )
    # Confirm combinations that can still enter the formal panel. An unrelated
    # rejected/exploratory pool must not consume the limited FoldX budget.
    if (mode == "strict" and not physical_preview and _complete_success(multi)
            and prediction_results is None and params.get("candidate_ids") is None
            and params.get("candidate_stage_limit") is None):
        from .candidate_search import supplemental_exclusive_doubles
        supplemental, supplemental_audit = supplemental_exclusive_doubles(
            candidates, sectors, [o for o in options if not o.review_required],
            graph, params["candidate_search"], areas)
        candidate_pool_diagnostics["spurs_feedback_supplement"] = supplemental_audit
        if supplemental:
            extra = _prediction(adapters["spurs_multi"], prediction_structure, chain,
                prediction_sequence, output_dir, structure_type, params.get("spurs_multi", {}),
                None, "spurs_multi", [candidate_notation(c) for c in supplemental])
            candidate_pool_diagnostics["spurs_feedback_worker"] = to_dict(extra)
            # A failed optional round cannot erase successful original model evidence.
            multi = PredictionResult(status=multi.status,
                records=multi.records + (extra.records if _complete_success(extra) else []),
                warnings=multi.warnings + extra.warnings, errors=multi.errors,
                provenance={**multi.provenance, "supplemental_run": extra.provenance},
                runtime_seconds=multi.runtime_seconds + extra.runtime_seconds)
            attach_hypothesis_responses(supplemental, hypotheses, graph,
                radius_a=response_radius, decay_length_a=response_decay)
            candidates.extend(supplemental)
            annotate_higher_order_regional_gain(candidates, sectors, areas)
            _annotate_stability(candidates, spurs.records if _complete_success(spurs) else [],
                multi, PredictionResult(status="failed"), params["spurs_multi"],
                params["foldx_confirmation"], position_map)
    confirmable = [
        item for item in candidates
        if item.order < 4 and item.qualification != "REVIEW_REQUIRED" and item.stability is not None
        and item.stability.consensus_status == "SPURS_ONLY_REVIEW_REQUIRED"
    ] if mode == "strict" else [item for item in candidates if item.order < 4]
    if mode == "strict" and not physical_preview:
        threshold = float(params["panel_optimizer"].get("minimum_coverage_gain", .05))
        rejected = [item for item in confirmable if item.order > 2 and (
            item.risk_components.get("regional_higher_order_gain") is None
            or not math.isfinite(item.risk_components["regional_higher_order_gain"])
            or item.risk_components["regional_higher_order_gain"] < threshold)]
        candidate_pool_diagnostics["foldx_preconfirmation_exclusions"] = [{
            "candidate_id": item.candidate_id,
            "reason": "HIGHER_ORDER_COVERAGE_BASELINE_MISSING" if item.risk_components.get("regional_higher_order_gain") is None
                      else "INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN",
            "observed": item.risk_components.get("regional_higher_order_gain"), "threshold": threshold,
            "measurement_basis": "mean_regional_area_coverage_9_12_15",
        } for item in sorted(rejected, key=lambda c: c.candidate_id)]
        rejected_ids = {item.candidate_id for item in rejected}
        confirmable = [item for item in confirmable if item.candidate_id not in rejected_ids]
    foldx_candidates = _foldx_shortlist(
        confirmable,
        int(params["foldx_confirmation"].get("max_candidates", 60)),
        int(params["panel_optimizer"].get("epitope_priority_per_sector", 3)),
        int(params["panel_optimizer"].get("stability_fallback_per_sector", 3)),
        hypotheses=hypotheses,
        sectors=sectors if not physical_preview else None,
        surface_areas=areas, graph=graph, diagnostics=candidate_pool_diagnostics,
    )
    if physical_preview and not coverage_preview:
        maximum=int(params['foldx_confirmation'].get('max_candidates',60))
        # Reserve most first-pass capacity for doubles. A few spatially distinct
        # triples retain all exact pair controls without consuming the entire queue.
        triples=_physical_confirmation_priority([c for c in confirmable if c.order==3],max(1,maximum//12),areas)
        required={tuple(sorted(o.notation for o in pair)) for triple in triples for pair in itertools.combinations(triple.mutations,2)}
        controls=[c for c in confirmable if c.order==2 and tuple(sorted(o.notation for o in c.mutations)) in required]
        selected_controls=controls[:max(0,maximum-len(triples))]
        used={c.candidate_id for c in triples+selected_controls}
        initial_support={p for c in selected_controls+triples for p in c.covered_positions_by_scenario.get('geodesic_12',[])}
        rest=_physical_confirmation_priority([c for c in confirmable if c.order==2 and c.candidate_id not in used],maximum,areas,initial_support)
        foldx_candidates=(selected_controls+triples+rest)[:maximum]
    if coverage_preview:
        foldx_candidates = _physical_confirmation_priority(
            [c for c in confirmable if c.order == 2 or params.get("candidate_ids") is not None],
            int(params['foldx_confirmation'].get('max_candidates',60)), areas)
    foldx = (
        _prediction(
            adapters["foldx"],
            prediction_structure,
            chain,
            prediction_sequence,
            output_dir,
            structure_type,
            params.get("foldx_confirmation", {}),
            prediction_results,
            "foldx",
            [candidate_notation(candidate) for candidate in foldx_candidates],
        )
        if foldx_candidates
        else PredictionResult(status="failed", errors=[{
            "code": "NO_CONFIRMABLE_CANDIDATES",
            "message": "FoldX was not run because no primary candidate passed the SPURS multi confirmation gate.",
        }])
    )
    _annotate_stability(
        candidates,
        spurs.records if _complete_success(spurs) else [],
        multi,
        foldx,
        params["spurs_multi"],
        params["foldx_confirmation"],
        position_map,
    )
    annotate_candidate_risks(candidates, records, hypotheses)
    allow_spurs_only = False
    params["allow_spurs_only_final"] = False
    minimum_gain = float(params["panel_optimizer"].get("minimum_coverage_gain", 0.05))
    by_sector: dict[str, list[V4Candidate]] = {}
    sector_flags = {sector.sector_id: set(sector.quality_flags) for sector in sectors}
    sasa_conflicts = {
        item.residue.uniprot_pos
        for item in surface
        if item.procaliper_classification_changed is True
        and item.sasa_agreement_status == "SASA_DISAGREEMENT"
    }
    sasa_context_sensitive = {
        item.residue.uniprot_pos
        for item in surface
        if item.native_context_classification_changed is True
    }
    sasa_method_unverified = {
        item.residue.uniprot_pos
        for item in surface
        if normalized_context_rows is not None
        and item.procaliper_rsasa is not None
        and item.sasa_agreement_status == "METHOD_SETTINGS_UNVERIFIED"
    }
    for candidate in candidates:
        by_sector.setdefault(candidate.sector_id, []).append(candidate)
    for candidate in candidates:
        # These proxy comparisons are used only by higher-order legacy gates.
        # Physical preview requires actual combination-field gain downstream.
        lower_coverage = 0.0 if candidate.order <= 2 or physical_preview else max(
            (
                max(item.coverage_by_scenario.values(), default=0.0)
                for item in (candidates if hypotheses else by_sector[candidate.sector_id])
                if item.order < candidate.order
            ),
            default=0.0,
        )
        information_gain = None
        if hypotheses and candidate.order > 2 and not physical_preview:
            def separation(item: V4Candidate) -> float:
                values = list(item.hypothesis_response.values())
                count = len(values)
                if count < 2:
                    return 0.0
                return sum(abs(left - right) for left, right in itertools.combinations(values, 2)) / (count * (count - 1) / 2)

            candidate_sites = {item.residue.uniprot_pos for item in candidate.mutations}
            comparable = [
                item for item in candidates
                if item.order < candidate.order
                and {mutation.residue.uniprot_pos for mutation in item.mutations} <= candidate_sites
            ]
            if comparable:
                candidate_values = list(candidate.hypothesis_response.values())
                candidate_mean = sum(candidate_values) / max(1, len(candidate_values))
                lower_mean = max(
                    sum(item.hypothesis_response.values()) / max(1, len(item.hypothesis_response))
                    for item in comparable
                )
                information_gain = max(
                    candidate_mean - lower_mean,
                    separation(candidate) - max(separation(item) for item in comparable),
                )
        candidate.formal_exclusion_reasons = _formal_candidate_reasons(
            candidate,
            allow_spurs_only=allow_spurs_only,
            minimum_coverage_gain=minimum_gain,
            lower_order_coverage=lower_coverage,
            higher_order_information_gain=information_gain,
            regional_coverage_gain=candidate.risk_components.get("regional_higher_order_gain"),
            regional_gain_required=not physical_preview,
            physical_coverage_pending=physical_preview,
            sasa_conflict_positions=sasa_conflicts,
            sasa_context_sensitive_positions=sasa_context_sensitive,
            sasa_method_unverified_positions=sasa_method_unverified,
            sector_quality_flags=sector_flags.get(candidate.sector_id),
        )
        if candidate.order > 2 and not physical_preview:
            candidate.risk_components["higher_order_information_gain"] = information_gain
            candidate.risk_components["higher_order_coverage_gain"] = max(candidate.coverage_by_scenario.values(), default=0.) - lower_coverage
        candidate.formal_eligible = not candidate.formal_exclusion_reasons
        if physical_preview and candidate.order==3 and candidate.formal_eligible:
            candidate.warnings.append('REAL_HIGHER_ORDER_GAIN_PENDING')
    selectable = (
        candidates if mode == "preview" else [item for item in candidates if item.formal_eligible]
    )
    if physical_preview:
        from ..surface_perturbation.combo_contract import waiting_confirmation
        selectable=[c for c in candidates if c.formal_eligible or waiting_confirmation(c)]
        if params.get('optimization_profile') == PROFILE_V3:
            selectable += [c for c in candidates if c not in selectable and c.stability is not None
                and c.stability.consensus_status == 'STRUCTURAL_ENERGY_CONFLICT'
                and set(c.formal_exclusion_reasons) <= {'STRUCTURAL_ENERGY_CONFLICT'}]
    physical_summary = None
    if physical_preview:
        from ..surface_perturbation.engine import run_physical_preview
        supplement_index = 0
        seen_pair_ids = [(p['batch_path'],p['group'],p['run']) for row in foldx.records for p in row.get('paired_models',[])]
        def register_pairs(pairs):
            existing=set(seen_pair_ids)
            for pair in pairs:
                identity=(pair['batch_path'],pair['group'],pair['run'])
                if identity not in existing:
                    seen_pair_ids.append(identity);existing.add(identity)
        def supplement(candidate, runs, *, role="production", cohort="production_v3", slot_start=0):
            nonlocal supplement_index
            supplement_index += 1
            result = _prediction(adapters["foldx"], prediction_structure, chain, prediction_sequence,
                                 output_dir/"surface_perturbation"/f"supplement_{supplement_index:04d}", structure_type,
                                 {**params["foldx_confirmation"], "number_of_runs": runs,
                                  "paired_cache_exclude_ids": seen_pair_ids, **({"sampling_role": role, "sampling_cohort": cohort, "sampling_slot_start": slot_start} if params.get("optimization_profile") == PROFILE_V3 else {})}, None, "foldx", [candidate_notation(candidate)])
            if not _complete_success(result):
                raise RuntimeError("FoldX supplemental repeats failed: " + str(result.errors))
            pairs = result.records[0].get("paired_models", [])
            seen_pair_ids.extend((p['batch_path'],p['group'],p['run']) for p in pairs)
            return pairs
        def recheck(candidate, pairs):
            import numpy as np
            scores = [p.get("ddg_kcal_mol") for p in pairs]
            if not scores or any(v is None or not np.isfinite(v) for v in scores): return False
            combined = PredictionResult(status="success", records=[{"mutations":candidate_notation(candidate),"ddg":float(np.median(scores))}])
            _annotate_stability([candidate],spurs.records,multi,combined,params["spurs_multi"],params["foldx_confirmation"],position_map)
            passed=candidate.stability.consensus_status == "CONSENSUS_PASS"
            if passed:
                cleared={'SPURS_ONLY_REVIEW_REQUIRED','HIGHER_ORDER_REQUIRES_FOLDX','COMBINATION_FOLDX_MISSING'}
                if params.get('optimization_profile') == PROFILE_V3:cleared.add('STRUCTURAL_ENERGY_CONFLICT')
                candidate.formal_exclusion_reasons=[r for r in candidate.formal_exclusion_reasons if r not in cleared]
                candidate.formal_eligible=not candidate.formal_exclusion_reasons
            elif params.get('optimization_profile') == PROFILE_V3:
                candidate.formal_eligible=False
                reason=candidate.stability.consensus_status
                if reason not in candidate.formal_exclusion_reasons:candidate.formal_exclusion_reasons.append(reason)
            return passed and candidate.formal_eligible
        if params.get('surface_perturbation',{}).get('audit_candidate_routes',False):
            from ..surface_perturbation.routes import search_controls
            write_json(output_dir/'surface_perturbation_candidate_routes.json',search_controls(
                sectors,options,graph,params['candidate_search'],hypotheses,selectable))
        selected, physical_summary = run_physical_preview(structure=structure_pdb, records=records, options=options,
                candidates=selectable, paired_records=foldx.records if _complete_success(foldx) else [],
                notation=candidate_notation, supplement=supplement, recheck=recheck, output=output_dir,register_pairs=register_pairs,
                settings={**params.get("surface_perturbation", {}),'ecd_length':len(protein.ecd_sequence),"primary_budget":params.get("panel_budget") or {"maximum_constructs":params.get("max_primary_constructs",7),"adaptive_budget":False}})
        selected_ids = [c.candidate_id for c in selected]
        for candidate in selected: candidate.role="PRIMARY"
        curve = physical_summary.get("coverage_curve", [])
        coverage = curve[min(len(selected),12)-1] if selected and curve else {}
        panel = V4Panel(selected_candidate_ids=selected_ids,candidates_by_role={"PRIMARY":selected_ids},
                        achieved_coverage={},uncovered_components=[],solver_status=physical_summary["status"],
                        coverage_denominator_area_a2=physical_summary.get('reference_area_all_a2'),achieved_area_coverage={k:v for k,v in coverage.items() if k.startswith("C_") and v is not None})
    else:
      panel, selected = select_panel(
        selectable,
        max(1, int(params.get("max_primary_constructs", 7))),
        int(params["panel_optimizer"].get("higher_order_max_primary", 2)),
        int(params["panel_optimizer"].get("recommendations_per_sector", 6)),
        int(params["panel_optimizer"].get("epitope_priority_per_sector", 3)),
        int(params["panel_optimizer"].get("stability_fallback_per_sector", 3)),
        sectors=sectors,
        surface_areas={item.residue.uniprot_pos: item.mapped_surface_area_a2 for item in surface},
        hypotheses=hypotheses,
        risk_penalty_weights=params["panel_optimizer"].get("risk_penalty_weights", {}),
        panel_budget=params.get("panel_budget"),
        stop_parameters=params["panel_optimizer"],
    )
    if raw_hypotheses:
        from .identifiability import panel_metrics
        from dataclasses import replace
        attach_hypothesis_responses(selected, raw_hypotheses, graph, radius_a=response_radius, decay_length_a=response_decay)
        raw_metrics = panel_metrics(selected, raw_hypotheses)
        raw_by_id = {h.hypothesis_id: h for h in raw_hypotheses}
        family_hypotheses = [replace(raw_by_id[f.representative_footprint_id],
            hypothesis_id=f.family_id, residue_weights=f.residue_weights) for f in families]
        attach_hypothesis_responses(selected, family_hypotheses, graph, radius_a=response_radius, decay_length_a=response_decay)
        panel.hierarchy_resolution = {"raw": raw_metrics, "family": panel_metrics(selected, family_hypotheses),
            "interpretation": "Unresolved raw/family hypotheses remain audited; connectivity groups are not equivalence classes."}
        attach_hypothesis_responses(selected, hypotheses, graph, radius_a=response_radius, decay_length_a=response_decay)
    recommendation_source = {**(source_provenance or {}), "models": {
        name: {"provenance": worker.provenance, "status": worker.status}
        for name, worker in (("spurs_single", spurs), ("spurs_multi", multi), ("discotope", discotope), ("foldx", foldx), ("procaliper", procaliper))}}
    recommendation_source["resource_audits"] = {
        row["sector_id"]: row for row in candidate_pool_diagnostics.get("foldx_region_resources", [])}
    panel.sector_recommendations = build_sector_recommendations(candidates, sectors,
        {item.residue.uniprot_pos: item.mapped_surface_area_a2 for item in surface}, params, selected,
        graph=graph, records=records, source_provenance=recommendation_source, epitope_evidence=epitope)
    from .spatial_support import candidate_region_membership, site_signature
    foldx_ids = {c.candidate_id for c in foldx_candidates}
    for row in panel.sector_recommendations:
        associated = [c for c in candidates if row["sector_id"] in candidate_region_membership(c, sectors)["region_ids"]]
        exclusive = [c for c in associated if candidate_region_membership(c, sectors)["exclusive_region_id"] == row["sector_id"]]
        generated = next((r for r in candidate_pool_diagnostics.get("sector_resources", [])
                          if r["sector_id"] == row["sector_id"]), {})
        confirmation = next((r for r in candidate_pool_diagnostics.get("foldx_region_resources", [])
                             if r["sector_id"] == row["sector_id"]), {})
        row["resource_audit"] = {**generated, **confirmation,
            "pool_candidate_count": len(associated), "pool_exclusive_count": len(exclusive),
            "pool_exclusive_groups": len({site_signature(c) for c in exclusive}),
            "spurs_scored_count": sum(c.stability is not None and c.stability.spurs_multi_ddg_kcal_mol is not None for c in associated),
            "foldx_scheduled_count": sum(c.candidate_id in foldx_ids for c in associated),
            "foldx_success_count": sum(c.stability is not None and c.stability.foldx_ddg_kcal_mol is not None for c in associated),
            "formal_eligible_count": sum(c.formal_eligible for c in associated),
            "formal_exclusive_count": sum(c.formal_eligible for c in exclusive)}
        audit = row["resource_audit"]
        audit["reasons"] = []
        if any(audit.get("quota_deficits", {}).values()):
            audit["reasons"].append("现有搜索未找到足够的不同专属位点组；保留实际配额缺口，不放宽位点与几何规则。")
        confirmable_count = audit.get("confirmable_exclusive_groups", 0)
        scheduled_count = audit.get("scheduled_exclusive_groups", 0)
        if confirmable_count < 2:
            audit["reasons"].append("符合组合预审条件的专属位点组不足两个，无法安排完整的专属FoldX保障。")
        elif scheduled_count < 2:
            audit["reasons"].append("FoldX固定预算不足以安排两个专属位点组；剩余缺口如实报告。")
        if audit["foldx_success_count"] < audit["foldx_scheduled_count"]:
            audit["reasons"].append("部分已安排的FoldX组合没有返回有效评分，证据仍缺失。")
        if audit["formal_exclusive_count"] == 0:
            audit["reasons"].append("本区目前没有满足全部准入条件的专属组合；待审查建议不属于合格推荐。")
    candidate_pool_diagnostics["sector_resources_final"] = [row["resource_audit"] for row in panel.sector_recommendations]
    selected_positions = {option.residue.uniprot_pos for candidate in selected for option in candidate.mutations}
    panel.uncovered_components = [
        {"sector_id": row["sector_id"], "reason": (
            "NO_SELECTED_FEASIBLE_CANDIDATE" if row["formal_eligible"] else "NO_FEASIBLE_CANDIDATE")}
        for row in panel.sector_recommendations
        if not selected_positions.intersection(next(s.residue_positions for s in sectors if s.sector_id == row["sector_id"]))
    ]
    required = {"spurs_single": spurs, "spurs_multi": multi}
    missing = [name for name, result in required.items() if not _complete_success(result)]
    foldx_run_success = _complete_success(foldx)
    foldx_candidate_coverage = sum(
        item.stability is not None and item.stability.foldx_ddg_kcal_mol is not None
        for item in candidates
    )
    selected_panel_foldx_complete = bool(selected) and all(
        item.stability and item.stability.foldx_ddg_kcal_mol is not None
        and math.isfinite(item.stability.foldx_ddg_kcal_mol) for item in selected
    )
    multi_ok = bool(selected) and all(
        item.stability and item.stability.spurs_multi_ddg_kcal_mol is not None
        and math.isfinite(item.stability.spurs_multi_ddg_kcal_mol) for item in selected
    )
    if physical_preview or mode == "preview":
        recommendation_status = "PREVIEW_NOT_FOR_FORMAL_RECOMMENDATION"
    elif (
        missing
        or not selected
        or not multi_ok
        or (not selected_panel_foldx_complete and not allow_spurs_only)
        or (not allow_spurs_only and not foldx_run_success)
    ):
        recommendation_status = "BLOCKED_EXTERNAL_EVIDENCE"
    else:
        recommendation_status = "FORMAL_RECOMMENDATION"
    warnings = ["Joint V4-ML thresholds are uncalibrated engineering heuristics."]
    if len(hypotheses) > 30:
        warnings.append("More than 30 inference hypotheses remain after conservative merging; count is diagnostic only and no forced compression is applied.")
    if assay_format == "cell_surface":
        warnings.append(
            "Cell-surface assay accessibility is unassessed: membrane-relative geometry is a structural proxy, while intermembrane distance and receptor/Fab-specific reach remain unknown."
        )
    if sasa_conflicts:
        warnings.append(
            f"ProCaliper and native Shrake-Rupley disagree on the normalized prediction input at {len(sasa_conflicts)} positions; ProCaliper is diagnostic, and native surface eligibility remains authoritative."
        )
    if sasa_context_sensitive:
        warnings.append(
            f"Original-chain and normalized ECD SASA cross the surface threshold at {len(sasa_context_sensitive)} positions; affected candidates are excluded from strict recommendations because structure context changes classification."
        )
    if sasa_method_unverified:
        warnings.append(
            f"ProCaliper and native SASA settings differ at {len(sasa_method_unverified)} positions; ProCaliper values are diagnostic only there, while surface eligibility uses the native isolated-protein calculation."
        )
    if not _complete_success(discotope):
        warnings.append(
            "DiscoTope evidence is unavailable; candidates are stability-gated fallback "
            "recommendations without epitope-prior support."
        )
    if not _complete_success(procaliper):
        warnings.append(
            "ProCaliper SASA cross-check is unavailable; native SASA and explicit surface geometry remain authoritative."
        )
    else:
        warnings.append(
            "ProCaliper contributes an independently audited residue SASA cross-check only; this worker does not compute finite-probe steric accessibility."
        )
    if missing:
        warnings.append("Missing external evidence: " + ", ".join(missing) + ".")
    if not foldx_run_success:
        warnings.append("FoldX run did not complete successfully; no kcal/mol value was fabricated.")
    elif not selected:
        warnings.append(
            f"FoldX scored {foldx_candidate_coverage}/{len(candidates)} candidates, but no candidate passed panel selection."
        )
    elif not selected_panel_foldx_complete:
        warnings.append("FoldX run succeeded, but selected panel coverage is incomplete.")
    if mode == "strict" and not selected:
        warnings.append("No candidate satisfied all strict candidate-level evidence gates.")
    if mode == "preview":
        warnings.append("Preview output must not be treated as a formal recommendation.")
    if any(sector.sector_type == "COVERAGE_ONLY" for sector in sectors):
        warnings.append("At least one sector has LOW_EPITOPE_PRIOR.")
    if len(sectors) > max(1, int(params.get("max_sectors", 7))):
        warnings.append("Display region count exceeds the soft target to balance surface residue counts while preserving connectivity and diameter constraints; construct budget is unchanged.")
    legacy_regions = _legacy_regions(sectors, records)
    legacy_designs = _legacy_designs(
        selected, sectors, protein, original_cds, recommendation_status
    )
    legacy_designs.extend(
        _deconvolution_presets(legacy_designs, records, protein, original_cds)
    )
    legacy_designs.extend(_sector_recommendation_designs(panel, candidates,
        [design for design in legacy_designs if design.design_role == "STAGE1_PROBE"],
        sectors, protein, original_cds))
    result = JointV4Result(
        schema_version="4.0",
        recommendation_status=recommendation_status,
        tool_inventory=inventories,
        surface_evidence=surface,
        epitope_evidence=epitope,
        mutation_options=single,
        sectors=sectors,
        candidates=candidates,
        panel=panel,
        raw_hypotheses=raw_hypotheses,
        footprint_families=families,
        complexity=SurfaceComplexitySummary(
            eligible_residues=len(graph), eligible_area_a2=sum(areas.values()),
            display_regions=len(sectors), raw_footprints=len(getattr(surface_v2_result, "raw_footprints", []) or raw_hypotheses),
            footprint_families=len(families), inference_hypotheses=len(hypotheses),
            surface_components=len(set(components.values())),
            structural_domains=len({record_map[p].domain for p in graph}),
        ),
        surface_hypotheses=hypotheses,
        legacy_regions=legacy_regions,
        legacy_designs=legacy_designs,
        warnings=warnings,
        errors=[
            error
            for result in (procaliper, spurs, discotope, multi, foldx)
            for error in result.errors
        ],
        provenance={
            "sector_recommendation_policy": {
                "version": "coverage_first_regions_v2_spatial_complementarity",
                "membership": "actual_mutation_positions_in_final_display_region",
                "coverage_basis": "geodesic_0_eligible_surface_area",
                "regional_priority": "mean_9_12_15_sector_coverage_then_marginal_area_then_overlap_then_direct_area_then_risk",
                "geometric_support_interpretation": "neighborhood proxy only; no biological perturbation or minimum distance gate",
                "near_threshold": "independent_failed_gates_then_missing_evidence_then_known_normalized_excess",
                "budget_scope": "selected_primary_constructs_only; regional suggestions are alternatives",
            },
            "surface_perturbation": physical_summary,
            **dict(source_provenance or {}),
            "surface_graph": graph_provenance,
            "qualification_audit": qualification_audit(records, float(params["surface"]["rsasa_threshold"])),
            "hierarchy_raw_footprint_ids": [f.footprint_id for f in getattr(surface_v2_result, "raw_footprints", [])],
            "global_candidate_pool": candidate_pool_diagnostics,
            "review_candidate_pool": review_diagnostics,
            "recommendation_diagnostics": recommendation_diagnostics(
                candidates, selected,
                {"procaliper": procaliper, "spurs_single": spurs,
                 "spurs_multi": multi, "discotope": discotope, "foldx": foldx},
                mutation_options=single,
            ),
            "foldx_evidence_status": {
                "run_success": foldx_run_success,
                "scored_candidates": foldx_candidate_coverage,
                "total_candidates": len(candidates),
                "selected_panel_complete": selected_panel_foldx_complete,
            },
            "hypothesis_response_model": {
                "kind": "surface_geodesic_chemical_noisy_or_proxy",
                "radius_a": response_radius,
                "decay_length_a": response_decay,
                "source": "surface_v2_deduplicated_footprints",
                "hypothesis_count": len(hypotheses),
            },
            "assay_format": assay_format,
            "assay_accessibility_status": "UNASSESSED" if assay_format == "cell_surface" else "NOT_APPLICABLE",
        },
    )
    contexts = context_manifest(protein, records, chain)
    write_json(output_dir / "raw_surface_footprints_v411.json", to_dict(getattr(surface_v2_result, "raw_footprints", [])))
    contexts["prediction_context"]["normalized_path"] = (
        normalized_path.name if normalized_path else None
    )
    workers = {
        "procaliper": procaliper,
        "spurs_single": spurs,
        "spurs_multi": multi,
        "discotope": discotope,
        "foldx": foldx,
    }
    for name, worker in workers.items():
        provenance = worker.provenance
        if provenance.get("model_revision"):
            inventories[name].model_revision = str(provenance["model_revision"])
        if provenance.get("weights_sha256"):
            inventories[name].weights_sha256 = str(provenance["weights_sha256"])
        if provenance.get("version"):
            inventories[name].version = str(provenance["version"])
        if provenance.get("git_revision"):
            inventories[name].git_revision = str(provenance["git_revision"])
        if provenance.get("actual_device"):
            inventories[name].device = str(provenance["actual_device"])
    write_v4_outputs(
        output_dir,
        result,
        parameters=params,
        contexts=contexts,
        source_provenance=result.provenance,
        worker_results=workers,
        tool_commands={name: adapter.command for name, adapter in adapters.items()},
        normalized_structure_path=normalized_path,
    )
    write_json(output_dir / "native_surface_rows.json", native_rows)
    return result
