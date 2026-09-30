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
from .geometry_evidence import annotate_hypothesis_geometry
from .models import (
    CombinationStabilityEvidence,
    JointV4Result,
    PredictionRequest,
    PredictionResult,
    ToolInventory,
    V4Candidate,
    V4Sector,
)
from .optimizer import select_panel
from .reporting import write_v4_outputs
from .risk import annotate_candidate_risks
from .sectors import build_epitope_sectors
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


def _foldx_shortlist(
    candidates: list[V4Candidate],
    maximum: int,
    epitope_priority_per_sector: int = 3,
    stability_fallback_per_sector: int = 3,
) -> list[V4Candidate]:
    """Select a deterministic, sector-balanced FoldX confirmation set."""
    maximum = max(1, int(maximum))
    by_sector: dict[str, list[V4Candidate]] = {}
    for candidate in candidates:
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
        prioritized = [
            item for item in ranked if item.recommendation_class == "EPITOPE_PRIORITIZED"
        ]
        fallback = [
            item for item in ranked if item.recommendation_class == "STABILITY_FALLBACK"
        ]
        leading = (
            prioritized[: max(0, int(epitope_priority_per_sector))]
            + fallback[: max(0, int(stability_fallback_per_sector))]
        )
        leading_ids = {item.candidate_id for item in leading}
        by_sector[sector_id] = leading + [
            item for item in ranked if item.candidate_id not in leading_ids
        ]
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
            if mutations and value is not None:
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
            if mutations and value is not None:
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
            if any(value is not None for value in singles)
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
                if candidate.order == 2
                else "HIGHER_ORDER_EXTRAPOLATION_MEDIUM"
                if candidate.order == 3
                else "HIGHER_ORDER_EXTRAPOLATION_HIGH"
            ),
            consensus_status=status,
            disagreement=disagreement,
            provenance={"single_sum_used": single_sum is not None},
        )


def _formal_candidate_reasons(
    candidate: V4Candidate,
    *,
    allow_spurs_only: bool,
    minimum_coverage_gain: float,
    lower_order_coverage: float,
    higher_order_information_gain: float | None = None,
    sasa_conflict_positions: set[int] | None = None,
    sasa_context_sensitive_positions: set[int] | None = None,
    sasa_method_unverified_positions: set[int] | None = None,
    sector_quality_flags: set[str] | None = None,
) -> list[str]:
    reasons = []
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
    elif stability.consensus_status == "SPURS_ONLY_REVIEW_REQUIRED" and allow_spurs_only:
        pass
    else:
        reasons.append(stability.consensus_status)
    if candidate.order > 2:
        candidate_gain = (
            higher_order_information_gain
            if higher_order_information_gain is not None else
            max(candidate.coverage_by_scenario.values(), default=0.0) - lower_order_coverage
        )
        if candidate_gain < minimum_coverage_gain:
            reasons.append("INSUFFICIENT_HIGHER_ORDER_COVERAGE_GAIN")
        if stability is None or stability.foldx_ddg_kcal_mol is None:
            reasons.append("HIGHER_ORDER_REQUIRES_FOLDX")
    return sorted(set(reasons))


def _legacy_regions(sectors: list[V4Sector], records: list[ResidueRecord]) -> list[RegionRecord]:
    by_position = {item.uniprot_pos: item for item in records}
    output = []
    for index, sector in enumerate(sectors, 1):
        present = [by_position[p] for p in sector.residue_positions if p in by_position]
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
                rationale="Joint V4 connected sector; epitope tendency is a soft ranking target.",
                color=["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#6F4E7C"][
                    (index - 1) % 7
                ],
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
                mutable_residues=list(sector.residue_positions),
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
        patch_id = region_id.get(candidate.sector_id)
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
    if normalized_path is not None and _complete_success(procaliper):
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
    )
    options = compressed_options(
        records,
        single,
        int(params["spurs_single"].get("max_options_per_site", 3)),
        epitope,
        float(params["candidate_search"].get("epitope_priority_weight", 1.0)),
        float(params["discotope"].get("seed_percentile", 0.80)),
    )
    record_map = {item.uniprot_pos: item for item in records}
    sectors = build_epitope_sectors(
        surface,
        epitope,
        graph,
        record_map,
        min(7, max(1, int(params.get("max_sectors", 7)))),
        params.get("sectoring", {}),
    )
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
                "regions": to_dict(surface_v2_result.surface_regions),
            },
        )
    hypotheses = (
        hypotheses_from_surface_v2(surface_v2_result, records, epitope)
        if bool(params.get("use_surface_hypotheses", True)) else []
    )
    candidates = generate_candidates(
        sectors, options, graph, params["candidate_search"],
        hypotheses=hypotheses or None,
    )
    annotate_hypothesis_geometry(
        hypotheses, records,
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
            [candidate_notation(candidate) for candidate in candidates],
        )
        if candidates
        else PredictionResult(status="failed")
    )
    foldx_candidates = _foldx_shortlist(
        candidates,
        int(params["foldx_confirmation"].get("max_candidates", 42)),
        int(params["panel_optimizer"].get("epitope_priority_per_sector", 3)),
        int(params["panel_optimizer"].get("stability_fallback_per_sector", 3)),
    )
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
        else PredictionResult(status="failed")
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
    allow_spurs_only = bool(params.get("allow_spurs_only_final", False))
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
        lower_coverage = max(
            (
                max(item.coverage_by_scenario.values(), default=0.0)
                for item in by_sector[candidate.sector_id]
                if item.order < candidate.order
            ),
            default=0.0,
        )
        information_gain = None
        if hypotheses and candidate.order > 2:
            def separation(item: V4Candidate) -> float:
                values = list(item.hypothesis_response.values())
                count = len(values)
                if count < 2:
                    return 0.0
                return sum(abs(left - right) for left, right in itertools.combinations(values, 2)) / (count * (count - 1) / 2)

            candidate_sites = {item.residue.uniprot_pos for item in candidate.mutations}
            comparable = [
                item for item in by_sector[candidate.sector_id]
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
            sasa_conflict_positions=sasa_conflicts,
            sasa_context_sensitive_positions=sasa_context_sensitive,
            sasa_method_unverified_positions=sasa_method_unverified,
            sector_quality_flags=sector_flags.get(candidate.sector_id),
        )
        candidate.formal_eligible = not candidate.formal_exclusion_reasons
    selectable = (
        candidates if mode == "preview" else [item for item in candidates if item.formal_eligible]
    )
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
    )
    selected_sector_ids = {item.sector_id for item in selected}
    reported_uncovered = {item["sector_id"] for item in panel.uncovered_components}
    for sector in sectors:
        if sector.sector_id not in selected_sector_ids | reported_uncovered:
            panel.uncovered_components.append(
                {"sector_id": sector.sector_id, "reason": "NO_FEASIBLE_CANDIDATE"}
            )
    required = {"spurs_single": spurs, "spurs_multi": multi}
    missing = [name for name, result in required.items() if not _complete_success(result)]
    foldx_run_success = _complete_success(foldx)
    foldx_candidate_coverage = sum(
        item.stability is not None and item.stability.foldx_ddg_kcal_mol is not None
        for item in candidates
    )
    selected_panel_foldx_complete = bool(selected) and all(
        item.stability and item.stability.foldx_ddg_kcal_mol is not None for item in selected
    )
    multi_ok = bool(selected) and all(
        item.stability and item.stability.spurs_multi_ddg_kcal_mol is not None for item in selected
    )
    if mode == "preview":
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
    legacy_regions = _legacy_regions(sectors, records)
    legacy_designs = _legacy_designs(
        selected, sectors, protein, original_cds, recommendation_status
    )
    legacy_designs.extend(
        _deconvolution_presets(legacy_designs, records, protein, original_cds)
    )
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
            **dict(source_provenance or {}),
            "surface_graph": graph_provenance,
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
