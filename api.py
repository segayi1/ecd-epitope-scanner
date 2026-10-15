from __future__ import annotations

import asyncio
import csv
import json
import os
import subprocess
import uuid
from dataclasses import asdict
from pathlib import Path
from typing import Any, Literal

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field, model_validator

from .assay_profiles import DEFAULT_ASSAY_PROFILE_ID
from .experiment_projects import ExperimentStore, IntegrityError, RevisionConflict, SnapshotError
from .experimental_constructs import ConstructError
from .experimental_effects import (
    ProbeBatchEffect,
    WellCount,
    WellRate,
    calculate_antibody_variant_batch_effect,
    calculate_batch_effect,
    calculate_factorial_interaction,
    project_probe_effects,
    summarize_antibody_variant_effects,
    summarize_factorial_interactions,
)
from .experimental_gates import GateError, parse_compensation_matrix
from .flowjo_result_tables import FlowJoTableError, parse_flowjo_percentage_table
from .integrations.paraplume_service import (
    ParaplumeConfigurationError,
    ParaplumeRunner,
    model_fingerprint,
)
from .jobs import JobRunner, JobStore
from .models import (
    DesignRecord,
    FeatureRange,
    MutationChoice,
    ProteinContext,
    RegionRecord,
    ResidueRecord,
    to_dict,
)
from .pipeline import PipelineRequest
from .sequences import apply_mutations, coding_fragment
from .stability import evaluate_substitution
from .util import read_json, safe_name, write_json

ROOT = Path(__file__).resolve().parent
WEB_ROOT = ROOT / "web"
DATA_ROOT = ROOT / "workspace_data"
STORE = JobStore(DATA_ROOT)
RUNNER = JobRunner(STORE)
EXPERIMENTS = ExperimentStore(DATA_ROOT, STORE)
PARAPLUME_RUNNER = ParaplumeRunner(EXPERIMENTS)


class JobCreate(BaseModel):
    input_type: Literal["uniprot", "pdb_id"] = "uniprot"
    identifier: str = Field(min_length=4, max_length=32)
    accession: str | None = None
    chain: str | None = None
    ecd_start: int | None = None
    ecd_end: int | None = None
    original_cds: str | None = None
    strict_af3: bool = True
    allow_partial: bool = False
    structure_candidate_id: str | None = None
    assay_format: Literal["cell_surface", "purified"] = "cell_surface"
    design_engine: Literal["joint_v4_ml"] = "joint_v4_ml"
    ml_mode: Literal["strict", "preview"] = "strict"
    surface_strategy: Literal["legacy", "surface_perturbation_v1_preview"] = "legacy"
    optimization_profile: Literal["legacy", "coverage_scan_v2_preview", "coverage_scan_v3_preview"] = "legacy"
    candidate_ids: list[str] | None = None
    calibration_repeat_schedule: list[int] | None = None
    export_mode: Literal["checkpoint", "deliverable"] = "deliverable"
    model_assets_root: str | None = None
    case_time_budget_seconds: float | None = Field(default=None,gt=0,allow_inf_nan=False)
    case_byte_budget: int | None = Field(default=None,gt=0)
    max_sectors: int = 7
    target_sector_residues: int | None = Field(default=None, ge=10, le=15)
    maximum_sector_geodesic_diameter_a: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    maximum_sector_ca_diameter_a: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    max_primary_constructs: int | None = Field(default=None, ge=1, le=12)
    minimum_constructs: int | None = Field(default=None, ge=1, le=12)
    target_constructs: int | None = Field(default=None, ge=1, le=12)
    maximum_constructs: int | None = Field(default=None, ge=1, le=12)
    adaptive_budget: bool | None = None
    assay_profile: str | None = DEFAULT_ASSAY_PROFILE_ID
    allow_spurs_only_final: bool = False
    biological_context_mode: str = "provided_structure"

    @model_validator(mode="after")
    def validate_panel_budget(self):
        if self.optimization_profile in {"coverage_scan_v2_preview", "coverage_scan_v3_preview"}:
            self.surface_strategy = "surface_perturbation_v1_preview"
            if 'assay_profile' not in self.model_fields_set:
                self.assay_profile = None
            if self.maximum_constructs is None and self.max_primary_constructs is None:
                self.maximum_constructs = 12
        from .joint_design_v4.budget import resolve_budget, optimizer_budget
        from .assay_profiles import get_assay_profile
        profile = get_assay_profile(self.assay_profile)
        budget = resolve_budget(minimum=self.minimum_constructs, target=self.target_constructs,
                       maximum=self.maximum_constructs, adaptive=self.adaptive_budget,
                       legacy_max=self.max_primary_constructs,
                       policy=profile["design_policy"] if profile else None,
                       maximum_limit=12 if self.surface_strategy == "surface_perturbation_v1_preview" else 8)
        optimizer_budget(profile, budget)
        PipelineRequest(**self.model_dump())
        return self


class MutationInput(BaseModel):
    position: int
    mutant: str = Field(pattern="^[ACDEFGHIKLMNPQRSTVWY]$")


class CustomDesignInput(BaseModel):
    name: str = Field(default="Custom construct", min_length=1, max_length=80)
    mutations: list[MutationInput] = Field(min_length=1, max_length=12)
    force_high_risk: bool = False


class ExperimentProjectCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    source_job_id: str = Field(min_length=1, max_length=64)
    name: str = Field(min_length=1, max_length=120)


class StructureEvidenceFreezeRequest(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    expected_revision: int = Field(ge=1)


class AntigenMutationInput(BaseModel):
    position: int = Field(ge=1)
    wt: str = Field(pattern="^[ACDEFGHIKLMNPQRSTVWY]$")
    mut: str = Field(pattern="^[ACDEFGHIKLMNPQRSTVWY]$")


class AntigenConstructCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    name: str = Field(min_length=1, max_length=120)
    source: Literal["CUSTOM", "RECOMMENDED"] = "CUSTOM"
    recommended_construct_id: str | None = None
    mutations: list[AntigenMutationInput] | None = None
    full_sequence: str | None = None
    expected_revision: int = Field(ge=1)


class AntigenConstructPreview(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    mutations: list[AntigenMutationInput] | None = None
    full_sequence: str | None = None


class AntibodyConstructCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    name: str = Field(min_length=1, max_length=120)
    role: Literal["PARENT", "VARIANT"]
    parent_construct_id: str | None = None
    vh_sequence: str
    vl_sequence: str
    expected_revision: int = Field(ge=1)


class ParaplumeTaskCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    antibody_construct_ids: list[str] = Field(min_length=1, max_length=20)
    model_version: Literal["2.0.0"] = "2.0.0"
    model_size: Literal["large"] = "large"
    expected_revision: int = Field(ge=1)


class ConstraintEvidenceInput(BaseModel):
    constraint: str = Field(min_length=1, max_length=200)
    evidence_id: str = Field(min_length=1, max_length=64)
    rationale: str = Field(min_length=1, max_length=1000)
    details: dict[str, Any] = Field(default_factory=dict)


class ConstraintExportCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    antibody_construct_id: str = Field(min_length=1, max_length=64)
    experimental_ir: list[ConstraintEvidenceInput] = Field(default_factory=list, max_length=500)
    predicted_ir: list[ConstraintEvidenceInput] = Field(default_factory=list, max_length=500)
    manual_rpr: list[ConstraintEvidenceInput] = Field(default_factory=list, max_length=200)
    merge_sources: list[Literal["experimental_ir", "predicted_ir", "manual_rpr"]] = Field(default_factory=list)
    reviewer: str = Field(min_length=1, max_length=120)
    expected_revision: int = Field(ge=1)


class AssayBatchCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    name: str = Field(min_length=1, max_length=120)
    date: str = Field(pattern=r"^\d{4}-\d{2}-\d{2}$")
    assay: Literal["ACTIVATION", "ANTIGEN_MYC", "FAB_MYC", "OTHER"]
    cell_line: str = Field(min_length=1, max_length=120)
    dose: float | None = Field(default=None, ge=0)
    dose_unit: str | None = None
    positive_gate_id: str | None = None
    independent_batch: bool = True
    settings: dict[str, Any] = Field(default_factory=dict)
    expected_revision: int = Field(ge=1)


class SummaryMeasurementCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    assay: Literal["ANTIGEN_MYC", "FAB_MYC"]
    construct_id: str = Field(min_length=1, max_length=64)
    measured_on: str = Field(pattern=r"^\d{4}-\d{2}-\d{2}$")
    metric: Literal["POSITIVE_FRACTION", "MEDIAN_FLUORESCENCE", "GEOMETRIC_MEAN", "OTHER"]
    value: float = Field(ge=0)
    unit: str = Field(default="", max_length=80)
    source_label: str = Field(min_length=1, max_length=240)
    notes: str = Field(default="", max_length=1000)
    reviewer: str = Field(min_length=1, max_length=120)
    expected_revision: int = Field(ge=1)


class TemporarySampleIdentityInput(BaseModel):
    status: Literal["TEMPORARY_LABEL"] = "TEMPORARY_LABEL"
    source: Literal["USER_CONFIRMED_LABELS"] = "USER_CONFIRMED_LABELS"
    antigen_label: str | None = Field(default=None, min_length=1, max_length=120)
    antigen_role: Literal["REFERENCE", "VARIANT", "NONE"]
    panel_label: str = Field(min_length=1, max_length=120)
    panel_role: Literal["ANTIBODY", "POSITIVE_LIGAND_CONTROL"]


class SampleMatchInput(BaseModel):
    upload_id: str
    batch_id: str | None = None
    role: Literal["ACTIVATION", "NC", "ANTIGEN_MYC", "FAB_MYC", "OTHER_CONTROL"]
    antigen_construct_id: str | None = None
    antibody_construct_id: str | None = None
    temporary_identity: TemporarySampleIdentityInput | None = None
    replicate_type: Literal["SINGLE", "TECHNICAL", "BIOLOGICAL"]
    replicate_index: int = Field(ge=1)
    condition: dict[str, Any] = Field(default_factory=dict)


class SampleConfirmInput(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    upload_session_id: str
    matches: list[SampleMatchInput]
    expected_revision: int = Field(ge=1)


class SamplePlotRequest(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    x_channel: str = Field(min_length=1, max_length=120)
    y_channel: str = Field(min_length=1, max_length=120)
    maximum_points: int = Field(default=3000, ge=100, le=10_000)


class FlowJoGateRebuildInput(BaseModel):
    flowjo_import_id: str = Field(min_length=1, max_length=64)
    reviewer: str = Field(min_length=1, max_length=120)
    reason: str = Field(min_length=1, max_length=1000)


class GateRevisionCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    name: str = Field(min_length=1, max_length=120)
    parent_gate_revision_id: str | None = None
    gate_tree: list[dict[str, Any]] = Field(min_length=1, max_length=256)
    transforms: list[dict[str, Any]] = Field(default_factory=list, max_length=64)
    compensation_decision: str
    compensation_artifact_id: str | None = None
    flowjo_rebuild: FlowJoGateRebuildInput | None = None
    expected_revision: int = Field(ge=1)


class GateQCPolicyInput(BaseModel):
    minimum_gate_events: int = Field(default=500, ge=0, le=10_000_000)
    time_window_count: int = Field(default=20, ge=3, le=100)
    time_mad_multiplier: float = Field(default=5.0, gt=0, le=100)
    boundary_fraction_threshold: float = Field(default=0.01, gt=0, le=1)


class EventWindowRangeInput(BaseModel):
    event_start: int = Field(ge=0)
    event_end: int = Field(gt=0)


class SampleEventWindowExclusionInput(BaseModel):
    sample_id: str = Field(min_length=1, max_length=64)
    windows: list[EventWindowRangeInput] = Field(min_length=1, max_length=100)


class GateAnalysisCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    gate_revision_id: str = Field(min_length=1, max_length=64)
    fcs_artifact_ids: list[str] = Field(min_length=1, max_length=128)
    qc_policy: GateQCPolicyInput = Field(default_factory=GateQCPolicyInput)
    source_gate_analysis_id: str | None = Field(default=None, max_length=64)
    window_exclusions: list[SampleEventWindowExclusionInput] = Field(default_factory=list, max_length=128)
    exclusion_reviewer: str | None = Field(default=None, max_length=120)
    exclusion_reason: str = Field(default="", max_length=1000)
    expected_revision: int = Field(ge=1)


class EffectSampleInput(BaseModel):
    sample_id: str = Field(min_length=1, max_length=64)
    probe_id: str = Field(min_length=1, max_length=64)
    condition: Literal["MUTANT", "REFERENCE", "MUTANT_NC", "REFERENCE_NC"]


class EffectAnalysisCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    gate_analysis_id: str = Field(min_length=1, max_length=64)
    positive_gate_id: str = Field(min_length=1, max_length=120)
    samples: list[EffectSampleInput] = Field(min_length=4, max_length=128)
    comparison_mode: Literal["ANTIGEN_EFFECT", "ANTIBODY_EFFECT", "DOUBLE_MUTANT_INTERACTION"] = "ANTIGEN_EFFECT"
    min_reference_window: float = Field(default=0.1, ge=0, le=1)
    bootstrap_seed: int = Field(default=0, ge=0, le=2_147_483_647)
    bootstrap_iterations: int = Field(default=2000, ge=100, le=100_000)
    expected_revision: int = Field(ge=1)


class FlowJoPercentageEffectCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    min_reference_window: float = Field(default=0.10, ge=0, le=1)
    expected_revision: int = Field(ge=1)


class QCReviewCreate(BaseModel):
    schema_version: Literal["1.0"] = "1.0"
    gate_analysis_id: str | None = Field(default=None, max_length=64)
    flowjo_import_id: str | None = Field(default=None, max_length=64)
    flowjo_percentage_import_id: str | None = Field(default=None, max_length=64)
    sample_id: str = Field(min_length=1, max_length=64)
    status: Literal["PASS", "PASS_WITH_OVERRIDE", "HOLD", "EXCLUDE"]
    reason: str = Field(default="", max_length=1000)
    reviewer: str = Field(min_length=1, max_length=120)
    expected_revision: int = Field(ge=1)


app = FastAPI(
    title="ECD Epitope Workbench",
    version="1.0.0",
    description="Auditable ECD surface-region mutation design for single-pass membrane proteins",
)


def _job_or_404(job_id: str) -> dict[str, Any]:
    try:
        return STORE.get(job_id)
    except KeyError as exc:
        raise HTTPException(404, "Job not found") from exc


def _public_job(job: dict[str, Any], *, view: str = "full") -> dict[str, Any]:
    result = job.get("result")
    if result and job["status"] == "COMPLETE":
        result = dict(result)
        result["structure_url"] = f"/api/jobs/{job['id']}/files/ecd_structure.cif"
        result["download_url"] = f"/api/jobs/{job['id']}/files/ecd_epitope_design.zip"
        result["job_id"] = job["id"]
        if view == "display":
            from .joint_design_v4.reporting import project_display_analysis
            result = project_display_analysis(result)
    output = dict(job)
    output["result"] = result
    request = dict(output.get("request") or {})
    if request.get("original_cds"):
        request["original_cds"] = "<provided>"
    output["request"] = request
    return output


def _job_summary(job: dict[str, Any]) -> dict[str, Any]:
    """Small list payload; full residue/design records stay on the detail API."""
    result = job.get("result") or {}
    protein = result.get("protein") or {}
    summary = result.get("summary") or {}
    return {
        key: job[key]
        for key in (
            "id",
            "status",
            "stage",
            "progress",
            "message",
            "error",
            "created_at",
            "updated_at",
        )
    } | {
        "accession": protein.get("accession"),
        "regions": summary.get("regions"),
        "recommended_constructs": summary.get("recommended_constructs"),
    }


@app.get("/api/health")
def health() -> dict[str, Any]:
    return {"status": "ok", "service": "ECD Epitope Workbench", "version": "1.0.0"}


@app.get("/favicon.ico", include_in_schema=False)
def favicon() -> Response:
    return Response(status_code=204)


@app.get("/api/jobs")
def list_jobs() -> list[dict[str, Any]]:
    return [_job_summary(item) for item in STORE.list()]


@app.post("/api/jobs", status_code=202)
def create_job(payload: JobCreate) -> dict[str, Any]:
    request = PipelineRequest(**payload.model_dump())
    job = STORE.create(request)
    RUNNER.submit(job["id"])
    return _public_job(job)


async def _save_upload(upload: UploadFile, path: Path, maximum: int = 80 * 1024 * 1024) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with path.open("wb") as handle:
        while chunk := await upload.read(1024 * 1024):
            total += len(chunk)
            if total > maximum:
                handle.close()
                path.unlink(missing_ok=True)
                raise HTTPException(413, "Uploaded file is too large")
            handle.write(chunk)


@app.post("/api/jobs/upload", status_code=202)
async def create_upload_job(
    structure: UploadFile = File(...),
    accession: str | None = Form(None),
    chain: str | None = Form(None),
    sequence: str | None = Form(None),
    ecd_start: int | None = Form(None),
    ecd_end: int | None = Form(None),
    original_cds: str | None = Form(None),
    allow_partial: bool = Form(False),
    design_engine: Literal["joint_v4_ml"] = Form("joint_v4_ml"),
    ml_mode: Literal["strict", "preview"] = Form("strict"),
    surface_strategy: Literal["legacy", "surface_perturbation_v1_preview"] = Form("legacy"),
    optimization_profile: Literal["legacy", "coverage_scan_v2_preview", "coverage_scan_v3_preview"] = Form("legacy"),
    candidate_ids: list[str] | None = Form(None),
    calibration_repeat_schedule: list[int] | None = Form(None),
    export_mode: Literal["checkpoint", "deliverable"] = Form("deliverable"),
    model_assets_root: str | None = Form(None),
    case_time_budget_seconds: float | None = Form(None,gt=0,allow_inf_nan=False),
    case_byte_budget: int | None = Form(None,gt=0),
    max_sectors: int = Form(7),
    target_sector_residues: int | None = Form(None, ge=10, le=15),
    maximum_sector_geodesic_diameter_a: float | None = Form(None, gt=0, allow_inf_nan=False),
    maximum_sector_ca_diameter_a: float | None = Form(None, gt=0, allow_inf_nan=False),
    max_primary_constructs: int | None = Form(None),
    minimum_constructs: int | None = Form(None),
    target_constructs: int | None = Form(None),
    maximum_constructs: int | None = Form(None),
    adaptive_budget: bool | None = Form(None),
    allow_spurs_only_final: bool = Form(False),
    biological_context_mode: str = Form("provided_structure"),
) -> dict[str, Any]:
    suffix = Path(structure.filename or "structure.cif").suffix.lower()
    if optimization_profile in {"coverage_scan_v2_preview", "coverage_scan_v3_preview"}:
        surface_strategy = "surface_perturbation_v1_preview"
        if maximum_constructs is None and max_primary_constructs is None:
            maximum_constructs = 12
    try:
        from .joint_design_v4.budget import resolve_budget
        resolve_budget(minimum=minimum_constructs, target=target_constructs,
                       maximum=maximum_constructs, adaptive=adaptive_budget, legacy_max=max_primary_constructs,
                       maximum_limit=12 if surface_strategy == "surface_perturbation_v1_preview" else 8)
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    if suffix not in {".cif", ".mmcif", ".pdb", ".ent"}:
        raise HTTPException(422, "Structure must be PDB or mmCIF")
    # Validate stage semantics before allocating a task or saving an upload.
    try:
        JobCreate(identifier=accession or 'USER_STRUCTURE',optimization_profile=optimization_profile,
            candidate_ids=candidate_ids,calibration_repeat_schedule=calibration_repeat_schedule,
            export_mode=export_mode,model_assets_root=model_assets_root,
            case_time_budget_seconds=case_time_budget_seconds,case_byte_budget=case_byte_budget,
            surface_strategy=surface_strategy,minimum_constructs=minimum_constructs,
            target_constructs=target_constructs,maximum_constructs=maximum_constructs,
            max_primary_constructs=max_primary_constructs,adaptive_budget=adaptive_budget,
            target_sector_residues=target_sector_residues,
            maximum_sector_geodesic_diameter_a=maximum_sector_geodesic_diameter_a,
            maximum_sector_ca_diameter_a=maximum_sector_ca_diameter_a)
    except ValueError as exc:
        raise HTTPException(422,str(exc)) from exc
    placeholder = PipelineRequest(
        optimization_profile=optimization_profile,
        candidate_ids=candidate_ids, calibration_repeat_schedule=calibration_repeat_schedule,
        export_mode=export_mode, model_assets_root=model_assets_root,
        case_time_budget_seconds=case_time_budget_seconds, case_byte_budget=case_byte_budget,
        surface_strategy=surface_strategy,
        input_type="structure",
        identifier=structure.filename or "structure",
        accession=accession or None,
        chain=chain or None,
        sequence=sequence or None,
        ecd_start=ecd_start,
        ecd_end=ecd_end,
        original_cds=original_cds or None,
        allow_partial=allow_partial,
        assay_profile=None if optimization_profile in {"coverage_scan_v2_preview", "coverage_scan_v3_preview"} else DEFAULT_ASSAY_PROFILE_ID,
        minimum_regions=1,
        maximum_regions=6,
        design_engine=design_engine,
        ml_mode=ml_mode,
        max_sectors=max_sectors,
        target_sector_residues=target_sector_residues,
        maximum_sector_geodesic_diameter_a=maximum_sector_geodesic_diameter_a,
        maximum_sector_ca_diameter_a=maximum_sector_ca_diameter_a,
        max_primary_constructs=max_primary_constructs,
        minimum_constructs=minimum_constructs, target_constructs=target_constructs,
        maximum_constructs=maximum_constructs, adaptive_budget=adaptive_budget,
        allow_spurs_only_final=allow_spurs_only_final,
        biological_context_mode=biological_context_mode,
    )
    job = STORE.create(placeholder)
    target = Path(job["workdir"]) / "uploads" / f"structure{suffix}"
    await _save_upload(structure, target)
    request = STORE.request_for(job["id"])
    request.structure_path = str(target)
    STORE.update(job["id"], request_json=request.as_dict() | {"original_cds": request.original_cds})
    RUNNER.submit(job["id"])
    return _public_job(STORE.get(job["id"]))


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str, view: Literal["full", "display"] = "full") -> dict[str, Any]:
    return _public_job(_job_or_404(job_id), view=view)


@app.get("/api/jobs/{job_id}/events")
async def job_events(job_id: str) -> StreamingResponse:
    _job_or_404(job_id)

    async def stream():
        previous = ""
        while True:
            job = _public_job(_job_or_404(job_id))
            payload = json.dumps(
                {
                    "id": job["id"],
                    "status": job["status"],
                    "stage": job["stage"],
                    "progress": job["progress"],
                    "message": job["message"],
                    "error": job["error"],
                },
                ensure_ascii=False,
            )
            if payload != previous:
                yield f"data: {payload}\n\n"
                previous = payload
            if job["terminal"]:
                break
            await asyncio.sleep(0.6)

    return StreamingResponse(
        stream(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-store", "X-Accel-Buffering": "no"},
    )


@app.post("/api/jobs/{job_id}/af3", status_code=202)
async def upload_af3(
    job_id: str,
    structure: UploadFile = File(...),
    confidence: UploadFile = File(...),
    official_source_confirmed: bool = Form(...),
) -> dict[str, Any]:
    job = _job_or_404(job_id)
    if job["status"] != "WAITING_FOR_AF3_UPLOAD":
        raise HTTPException(409, "Job is not waiting for AF3 files")
    if not official_source_confirmed:
        raise HTTPException(422, "Official AF3 source confirmation is required")
    structure_name = (structure.filename or "model.cif").lower()
    confidence_name = (confidence.filename or "confidences.json").lower()
    if not structure_name.endswith((".cif", ".mmcif")):
        raise HTTPException(422, "AF3 structure must be mmCIF")
    if not confidence_name.endswith(".json"):
        raise HTTPException(422, "AF3 confidence file must be JSON")
    upload_dir = Path(job["workdir"]) / "af3_upload"
    structure_path = upload_dir / "official_af3_model.cif"
    confidence_path = upload_dir / "official_af3_confidences.json"
    await _save_upload(structure, structure_path)
    await _save_upload(confidence, confidence_path, maximum=200 * 1024 * 1024)
    try:
        confidence_data = read_json(confidence_path)
    except Exception as exc:
        raise HTTPException(422, f"Invalid confidence JSON: {exc}") from exc
    if not isinstance(confidence_data, dict) or not any(
        key in confidence_data for key in ("pae", "atom_plddts", "predicted_aligned_error")
    ):
        raise HTTPException(422, "Confidence JSON lacks AF3 PAE/pLDDT arrays")
    RUNNER.resume_with_af3(job_id, structure_path, confidence_path)
    return _public_job(STORE.get(job_id))


def _hydrate_for_custom(
    job: dict[str, Any],
) -> tuple[
    dict[str, Any], ProteinContext, list[ResidueRecord], list[RegionRecord], list[DesignRecord]
]:
    analysis = read_json(Path(job["workdir"]) / "analysis.json")
    protein_data = analysis["protein"]
    uniprot_path = Path(job["workdir"]) / "uniprot.json"
    if uniprot_path.exists():
        sequence = read_json(uniprot_path)["sequence"]["value"]
    else:
        sequence = (job.get("request") or {}).get("sequence") or protein_data["ecd_sequence"]
    features = [FeatureRange(**item) for item in protein_data.get("features", [])]
    protein = ProteinContext(
        accession=protein_data["accession"],
        protein_name=protein_data["protein_name"],
        organism=protein_data.get("organism", ""),
        sequence=sequence,
        ecd_start=protein_data["ecd_start"],
        ecd_end=protein_data["ecd_end"],
        membrane_anchor_position=protein_data["membrane_anchor_position"],
        topology=protein_data["topology"],
        reviewed=protein_data.get("reviewed", False),
        features=features,
    )
    records = [ResidueRecord(**item) for item in analysis["residues"]]
    regions = [RegionRecord(**item) for item in analysis["regions"]]
    designs: list[DesignRecord] = []
    for item in analysis["designs"]:
        data = dict(item)
        data["mutations"] = [MutationChoice(**mutation) for mutation in data.get("mutations", [])]
        designs.append(DesignRecord(**data))
    return analysis, protein, records, regions, designs


@app.post("/api/jobs/{job_id}/designs")
def add_custom_design(job_id: str, payload: CustomDesignInput) -> dict[str, Any]:
    job = _job_or_404(job_id)
    if job["status"] != "COMPLETE":
        raise HTTPException(409, "Custom designs require a completed analysis")
    analysis, protein, records, regions, designs = _hydrate_for_custom(job)
    by_position = {record.uniprot_pos: record for record in records}
    choices: list[MutationChoice] = []
    errors: list[str] = []
    for item in payload.mutations:
        record = by_position.get(item.position)
        if record is None:
            errors.append(f"Residue {item.position} is outside the mapped ECD")
            continue
        if record.hard_exclusions:
            errors.append(
                f"{record.aa}{item.position}: protected by {', '.join(record.hard_exclusions)}"
            )
            continue
        choice = evaluate_substitution(record, item.mutant, protein)
        if not choice.allowed and not payload.force_high_risk:
            errors.append(f"{choice.notation}: {', '.join(choice.reasons)}")
            continue
        choices.append(choice)
    if errors:
        raise HTTPException(
            422, detail={"message": "Custom design failed validation", "issues": errors}
        )
    if len({choice.position for choice in choices}) != len(choices):
        raise HTTPException(422, "Each residue may appear only once in a construct")
    mutant_sequence = apply_mutations(protein.ecd_sequence, choices, protein.ecd_start)
    # Preserve the original construct context: users may add a custom probe
    # after an analysis specifically to order a full-length cell-surface clone.
    original_cds = STORE.request_for(job_id).original_cds
    mutant_dna, dna_source = coding_fragment(protein, mutant_sequence, original_cds)
    design = DesignRecord(
        id=f"CUSTOM_{len([item for item in designs if item.custom]) + 1:02d}",
        name=payload.name,
        region_id=None,
        mutations=choices,
        status="USER_CONFIRMED_HIGH_RISK" if payload.force_high_risk else "USER_CONFIRMED",
        stability_engine="builtin structural-risk proxy",
        combined_stability_score=round(sum(item.stability_score for item in choices), 4),
        perturbation_score=round(sum(item.perturbation_score for item in choices), 4),
        rationale="User-defined cross-region construct validated against mapping and protected-site rules.",
        warnings=(
            ["High-risk substitutions were explicitly forced by the user."]
            if payload.force_high_risk
            else []
        )
        + [f"Coding fragment source: {dna_source}."],
        mutant_ecd_sequence=mutant_sequence,
        mutant_coding_sequence=mutant_dna,
        custom=True,
    )
    designs.append(design)
    analysis["designs"] = [to_dict(item) for item in designs]
    analysis["design_summary"]["final_library_size"] = len(designs)
    analysis["summary"]["custom_constructs"] = sum(item.custom for item in designs)
    output_dir = Path(job["workdir"])
    from .exporting import write_scientific_outputs

    files = write_scientific_outputs(
        output_dir, analysis, protein, records, regions, designs, original_cds=original_cds
    )
    analysis["files"] = files
    write_json(output_dir / "analysis.json", analysis)
    STORE.update(job_id, result_json=analysis, message=f"Custom design {design.id} saved")
    return {"design": to_dict(design), "job": _public_job(STORE.get(job_id))}


@app.get("/api/jobs/{job_id}/export/offline")
def export_offline(job_id: str, lang: str = "zh") -> Response:
    """Self-contained offline report: structure + data + 3Dmol in one HTML."""
    job = _job_or_404(job_id)
    if job.get("status") != "COMPLETE":
        raise HTTPException(409, "Job is not complete; nothing to export")
    workdir = Path(job["workdir"])
    analysis = read_json(workdir / "analysis.json")
    from .exporting import build_offline_html

    html = build_offline_html(analysis, workdir=workdir, lang=lang)
    return Response(
        html,
        media_type="text/html; charset=utf-8",
        headers={
            "Content-Disposition": (
                f'attachment; filename="ecd_epitope_offline_{safe_name(job_id)}.html"'
            )
        },
    )


@app.get("/api/jobs/{job_id}/files/{filename}")
def job_file(job_id: str, filename: str) -> FileResponse:
    job = _job_or_404(job_id)
    safe = safe_name(filename)
    if safe != filename:
        raise HTTPException(400, "Invalid filename")
    path = Path(job["workdir"]) / safe
    try:
        path.resolve().relative_to(Path(job["workdir"]).resolve())
    except ValueError as exc:
        raise HTTPException(400, "Invalid result path") from exc
    if not path.is_file():
        raise HTTPException(404, "Result file not found")
    media = "application/zip" if path.suffix == ".zip" else None
    return FileResponse(path, media_type=media, filename=path.name)


def _experiment_response(data: Any) -> dict[str, Any]:
    return {"schema_version": "1.0", "data": data, "warnings": [], "errors": []}


@app.post("/api/v4/experiment-projects", status_code=201)
def create_experiment_project(payload: ExperimentProjectCreate) -> dict[str, Any]:
    try:
        project = EXPERIMENTS.create(payload.source_job_id, payload.name)
    except KeyError as exc:
        raise HTTPException(
            404, {"code": "SOURCE_JOB_NOT_FOUND", "message": "Source job not found"}
        ) from exc
    except SnapshotError as exc:
        raise HTTPException(422, {"code": exc.code, "message": str(exc)}) from exc
    return _experiment_response(project)


@app.get("/api/v4/experiment-projects")
def list_experiment_projects() -> dict[str, Any]:
    return _experiment_response(
        [
            {
                key: project[key]
                for key in (
                    "project_id",
                    "name",
                    "source_job_id",
                    "revision",
                    "status",
                    "created_at",
                    "frozen_hashes",
                    "integrity_errors",
                )
            }
            for project in EXPERIMENTS.list()
        ]
    )


@app.get("/api/v4/experiment-projects/{project_id}")
def get_experiment_project(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.get(project_id))
    except KeyError as exc:
        raise HTTPException(
            404, {"code": "PROJECT_NOT_FOUND", "message": "Project not found"}
        ) from exc


@app.get("/api/v4/experiment-projects/{project_id}/support-bundle")
def experiment_support_bundle(project_id: str, include_originals: bool = False) -> Response:
    try:
        content = EXPERIMENTS.build_support_bundle(
            project_id, include_originals=include_originals
        )
    except KeyError as exc:
        raise HTTPException(
            404, {"code": "PROJECT_NOT_FOUND", "message": "Project not found"}
        ) from exc
    except IntegrityError as exc:
        raise HTTPException(409, {"code": exc.code, "message": str(exc)}) from exc
    return Response(
        content=content,
        media_type="application/zip",
        headers={
            "Content-Disposition": f'attachment; filename="experiment_{project_id}_support.zip"'
        },
    )


@app.get("/api/v4/experiment-projects/{project_id}/structure-context")
def experiment_structure_context(project_id: str) -> dict[str, Any]:
    try:
        project = EXPERIMENTS.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        artifacts = {item["kind"]: item for item in project["artifacts"]}
        required = {"structure", "mapping", "surface_graph"}
        if not required.issubset(artifacts):
            raise IntegrityError("Frozen structure context is incomplete")
        paths = {}
        for kind in required:
            paths[kind], _ = EXPERIMENTS.artifact_path(project_id, artifacts[kind]["id"])
        with paths["mapping"].open(newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        mapping = []
        for row in rows:
            try:
                position = int(row["uniprot_pos"])
                auth_seq_id = int(row["auth_seq_id"])
            except (KeyError, TypeError, ValueError):
                continue
            chain_id = str(row.get("chain_id") or "").strip()
            if not chain_id:
                continue
            mapping.append({
                "uniprot_pos": position,
                "normalized_position": int(row["normalized_position"])
                if str(row.get("normalized_position", "")).isdigit() else None,
                "chain_id": chain_id,
                "auth_seq_id": auth_seq_id,
                "insertion_code": str(row.get("insertion_code") or "").strip(),
                "aa": str(row.get("one_letter") or "").strip(),
            })
        return _experiment_response({
            "structure": {
                "artifact_id": artifacts["structure"]["id"],
                "media_type": artifacts["structure"]["media_type"],
                "format": paths["structure"].suffix.lower().lstrip("."),
                "url": f"/api/v4/experiment-projects/{project_id}/files/{artifacts['structure']['id']}",
            },
            "mapping": mapping,
            "surface_graph": read_json(paths["surface_graph"]),
            "structure_evidence": EXPERIMENTS.read_structure_evidence(project_id),
            "frozen_hashes": project["frozen_hashes"],
        })
    except KeyError as exc:
        raise HTTPException(404, {"code": "PROJECT_NOT_FOUND", "message": "Project not found"}) from exc
    except IntegrityError as exc:
        raise HTTPException(409, {"code": exc.code, "message": str(exc)}) from exc


@app.post("/api/v4/experiment-projects/{project_id}/structure-evidence")
def freeze_project_structure_evidence(
    project_id: str, payload: StructureEvidenceFreezeRequest
) -> dict[str, Any]:
    try:
        result = EXPERIMENTS.freeze_source_structure_evidence(
            project_id, expected_revision=payload.expected_revision
        )
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except SnapshotError as exc:
        raise HTTPException(
            422, {"code": exc.code, "message": str(exc)}
        ) from exc
    return _experiment_response(result)


@app.get("/api/v4/experiment-projects/{project_id}/files/{artifact_id}")
def experiment_artifact(project_id: str, artifact_id: str) -> FileResponse:
    try:
        path, artifact = EXPERIMENTS.artifact_path(project_id, artifact_id)
    except KeyError as exc:
        raise HTTPException(
            404, {"code": "ARTIFACT_NOT_FOUND", "message": "Artifact not found"}
        ) from exc
    except IntegrityError as exc:
        raise HTTPException(409, {"code": exc.code, "message": str(exc)}) from exc
    return FileResponse(path, media_type=artifact["media_type"], filename=path.name)


def _experiment_write_error(exc: Exception) -> HTTPException:
    if isinstance(exc, RevisionConflict):
        return HTTPException(409, {"code": exc.code, "message": str(exc)})
    if isinstance(exc, IntegrityError):
        return HTTPException(409, {"code": exc.code, "message": str(exc)})
    if isinstance(exc, KeyError):
        return HTTPException(
            404, {"code": "RESOURCE_NOT_FOUND", "message": "Project or parent construct not found"}
        )
    if isinstance(exc, GateError):
        code = (
            "UNKNOWN_COMPENSATION_SOURCE"
            if "UNKNOWN_COMPENSATION_SOURCE" in str(exc)
            else "GATE_INVALID"
        )
        return HTTPException(422, {"code": code, "message": str(exc)})
    return HTTPException(422, {"code": "CONSTRUCT_INVALID", "message": str(exc)})


@app.post("/api/v4/experiment-projects/{project_id}/antigen-constructs", status_code=201)
def add_antigen_construct(project_id: str, payload: AntigenConstructCreate) -> dict[str, Any]:
    try:
        if payload.source == "RECOMMENDED":
            if not payload.recommended_construct_id:
                raise ConstructError("recommended_construct_id is required for RECOMMENDED source")
            if payload.mutations is not None or payload.full_sequence is not None:
                raise ConstructError("Recommended constructs must use the frozen v4 sequence")
            result = EXPERIMENTS.add_recommended_antigen_construct(
                project_id,
                recommended_construct_id=payload.recommended_construct_id,
                name=payload.name,
                expected_revision=payload.expected_revision,
            )
        else:
            if payload.recommended_construct_id is not None:
                raise ConstructError("recommended_construct_id requires RECOMMENDED source")
            result = EXPERIMENTS.add_antigen_construct(
                project_id,
                name=payload.name,
                mutations=[item.model_dump() for item in payload.mutations]
                if payload.mutations
                else None,
                full_sequence=payload.full_sequence,
                expected_revision=payload.expected_revision,
            )
    except (KeyError, IntegrityError, RevisionConflict, ConstructError) as exc:
        raise _experiment_write_error(exc) from exc
    return _experiment_response(result)


@app.post("/api/v4/experiment-projects/{project_id}/antigen-constructs/preview")
def preview_antigen_construct(
    project_id: str, payload: AntigenConstructPreview
) -> dict[str, Any]:
    try:
        result = EXPERIMENTS.preview_antigen_construct(
            project_id,
            mutations=[item.model_dump() for item in payload.mutations]
            if payload.mutations
            else None,
            full_sequence=payload.full_sequence,
        )
    except (KeyError, IntegrityError, ConstructError) as exc:
        raise _experiment_write_error(exc) from exc
    return _experiment_response(result)


@app.get("/api/v4/experiment-projects/{project_id}/antigen-constructs")
def list_antigen_constructs(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_entities(project_id, "antigen_constructs"))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc


@app.get("/api/v4/experiment-projects/{project_id}/recommended-antigens")
def list_recommended_antigens(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_recommended_antigens(project_id))
    except (KeyError, IntegrityError) as exc:
        raise _experiment_write_error(exc) from exc


@app.post("/api/v4/experiment-projects/{project_id}/antibody-constructs", status_code=201)
def add_antibody_construct(project_id: str, payload: AntibodyConstructCreate) -> dict[str, Any]:
    try:
        result = EXPERIMENTS.add_antibody_construct(
            project_id,
            name=payload.name,
            role=payload.role,
            parent_construct_id=payload.parent_construct_id,
            vh_sequence=payload.vh_sequence,
            vl_sequence=payload.vl_sequence,
            expected_revision=payload.expected_revision,
        )
    except (KeyError, IntegrityError, RevisionConflict, ConstructError) as exc:
        raise _experiment_write_error(exc) from exc
    return _experiment_response(result)


@app.get("/api/v4/experiment-projects/{project_id}/antibody-constructs")
def list_antibody_constructs(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_entities(project_id, "antibody_constructs"))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc


@app.post("/api/v4/experiment-projects/{project_id}/paraplume", status_code=202)
def create_paraplume_task(
    project_id: str, payload: ParaplumeTaskCreate
) -> dict[str, Any]:
    try:
        fingerprint = model_fingerprint()
    except ParaplumeConfigurationError as exc:
        raise HTTPException(503, {"code": exc.code, "message": str(exc)}) from exc
    except (OSError, ValueError, KeyError) as exc:
        raise HTTPException(
            503, {"code": "PARAPLUME_MODEL_UNAVAILABLE", "message": str(exc)}
        ) from exc
    try:
        task = EXPERIMENTS.create_paraplume_task(
            project_id,
            antibody_construct_ids=payload.antibody_construct_ids,
            model_fingerprint=fingerprint,
            expected_revision=payload.expected_revision,
        )
    except KeyError as exc:
        raise HTTPException(
            404, {"code": "ANTIBODY_CONSTRUCT_NOT_FOUND", "message": str(exc)}
        ) from exc
    except (IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        raise HTTPException(
            422, {"code": "PARAPLUME_REQUEST_INVALID", "message": str(exc)}
        ) from exc
    PARAPLUME_RUNNER.submit(task["id"])
    return _experiment_response(
        {
            "task_id": task["id"],
            "status": task["status"],
            "progress": task["progress"],
            "events_url": f"/api/v4/experiment-tasks/{task['id']}/events",
            "result_scope": "antibody-side paratope probabilities; not antigen epitope evidence",
        }
    )


@app.get("/api/v4/experiment-projects/{project_id}/paraplume-runs")
def list_paraplume_runs(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_paraplume_runs(project_id))
    except (KeyError, IntegrityError) as exc:
        raise _experiment_write_error(exc) from exc


@app.post("/api/v4/experiment-projects/{project_id}/constraint-exports", status_code=201)
def create_constraint_export(
    project_id: str, payload: ConstraintExportCreate
) -> dict[str, Any]:
    try:
        return _experiment_response(
            EXPERIMENTS.create_constraint_export(
                project_id,
                antibody_construct_id=payload.antibody_construct_id,
                experimental_ir=[item.model_dump() for item in payload.experimental_ir],
                predicted_ir=[item.model_dump() for item in payload.predicted_ir],
                manual_rpr=[item.model_dump() for item in payload.manual_rpr],
                merge_sources=payload.merge_sources,
                reviewer=payload.reviewer,
                expected_revision=payload.expected_revision,
            )
        )
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        message = str(exc)
        prefix = message.split(":", 1)[0]
        code = getattr(exc, "code", None)
        if code == "FASTA_CONSTRAINT_MISMATCH" and prefix.isupper() and "_" in prefix:
            code = prefix
        code = code or (prefix if prefix.isupper() and "_" in prefix else "CONSTRAINT_EXPORT_INVALID")
        raise HTTPException(422, {"code": code, "message": str(exc)}) from exc


@app.get("/api/v4/experiment-projects/{project_id}/constraint-exports")
def list_constraint_exports(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_constraint_exports(project_id))
    except (KeyError, IntegrityError) as exc:
        raise _experiment_write_error(exc) from exc


@app.get("/api/v4/experiment-tasks/{task_id}")
def get_experiment_task(task_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.get_experiment_task(task_id))
    except KeyError as exc:
        raise HTTPException(
            404, {"code": "TASK_NOT_FOUND", "message": "Experiment task not found"}
        ) from exc


@app.get("/api/v4/experiment-tasks/{task_id}/events")
async def experiment_task_events(task_id: str) -> StreamingResponse:
    try:
        EXPERIMENTS.get_experiment_task(task_id)
    except KeyError as exc:
        raise HTTPException(
            404, {"code": "TASK_NOT_FOUND", "message": "Experiment task not found"}
        ) from exc

    async def stream():
        previous = None
        while True:
            task = EXPERIMENTS.get_experiment_task(task_id)
            output = task.get("output") or {}
            message = output.get("message")
            payload = json.dumps(
                {
                    "task_id": task["id"],
                    "status": task["status"],
                    "progress": task["progress"],
                    "message": message,
                    "error_code": task["error_code"],
                    "output": output if task["status"] in {"COMPLETE", "FAILED", "INTERRUPTED"} else None,
                },
                ensure_ascii=False,
            )
            if payload != previous:
                yield f"data: {payload}\n\n"
                previous = payload
            if task["status"] in {"COMPLETE", "FAILED", "INTERRUPTED"}:
                break
            await asyncio.sleep(1.0)

    return StreamingResponse(
        stream(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-store", "X-Accel-Buffering": "no"},
    )


@app.post("/api/v4/experiment-projects/{project_id}/batches", status_code=201)
def add_assay_batch(project_id: str, payload: AssayBatchCreate) -> dict[str, Any]:
    try:
        result = EXPERIMENTS.add_batch(project_id, **payload.model_dump(exclude={"schema_version"}))
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    return _experiment_response(result)


@app.get("/api/v4/experiment-projects/{project_id}/batches")
def list_assay_batches(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_entities(project_id, "assay_batches"))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc


@app.post(
    "/api/v4/experiment-projects/{project_id}/summary-measurements", status_code=201
)
def add_summary_measurement(
    project_id: str, payload: SummaryMeasurementCreate
) -> dict[str, Any]:
    try:
        result = EXPERIMENTS.add_summary_measurement(
            project_id, **payload.model_dump(exclude={"schema_version"})
        )
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        raise HTTPException(
            422, {"code": "SUMMARY_MEASUREMENT_INVALID", "message": str(exc)}
        ) from exc
    return _experiment_response(result)


@app.get("/api/v4/experiment-projects/{project_id}/summary-measurements")
def list_summary_measurements(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_summary_measurements(project_id))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc


def _inspect_uploaded_fcs(path: Path) -> dict[str, Any]:
    python = ROOT / ".tool_envs" / "flowkit" / "Scripts" / "python.exe"
    worker = ROOT / "integrations" / "workers" / "flowkit_worker.py"
    if not python.is_file() or not worker.is_file():
        raise HTTPException(
            503, {"code": "FLOWKIT_WORKER_UNAVAILABLE", "message": "FlowKit worker is unavailable"}
        )
    request = {
        "schema_version": "1.0",
        "request_id": uuid.uuid4().hex,
        "operation": "inspect_fcs",
        "path": str(path),
    }
    try:
        completed = subprocess.run(
            [str(python), str(worker)],
            input=json.dumps(request),
            capture_output=True,
            text=True,
            encoding="utf-8",
            timeout=90,
            check=False,
            env=os.environ | {"PYTHONIOENCODING": "utf-8"},
        )
        result = json.loads(completed.stdout)
    except (subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
        raise HTTPException(422, {"code": "FCS_PARSE_FAILED", "message": str(exc)}) from exc
    if completed.returncode or result.get("status") != "COMPLETE":
        errors = result.get("errors") or [{"message": completed.stderr[-1000:]}]
        raise HTTPException(
            422,
            {"code": "FCS_PARSE_FAILED", "message": errors[0].get("message", "FCS parse failed")},
        )
    return result["records"]


def _import_uploaded_wsp(path: Path, fcs_paths: list[Path]) -> dict[str, Any]:
    python = ROOT / ".tool_envs" / "flowkit" / "Scripts" / "python.exe"
    worker = ROOT / "integrations" / "workers" / "flowkit_worker.py"
    if not python.is_file() or not worker.is_file():
        raise HTTPException(
            503, {"code": "FLOWKIT_WORKER_UNAVAILABLE", "message": "FlowKit worker is unavailable"}
        )
    request = {
        "schema_version": "1.0",
        "request_id": uuid.uuid4().hex,
        "operation": "apply_wsp_gates",
        "path": str(path),
        "fcs_paths": [str(item) for item in fcs_paths],
    }
    try:
        completed = subprocess.run(
            [str(python), str(worker)],
            input=json.dumps(request),
            capture_output=True,
            text=True,
            encoding="utf-8",
            timeout=180,
            check=False,
            env=os.environ | {"PYTHONIOENCODING": "utf-8"},
        )
        result = json.loads(completed.stdout)
    except (subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
        raise HTTPException(422, {"code": "FLOWJO_IMPORT_FAILED", "message": str(exc)}) from exc
    if completed.returncode or result.get("status") in {"FAILED", None}:
        errors = result.get("errors") or [{"message": completed.stderr[-1000:]}]
        raise HTTPException(
            422,
            {
                "code": "FLOWJO_IMPORT_FAILED",
                "message": errors[0].get("message", "FlowJo import failed"),
            },
        )
    return result


@app.post("/api/v4/experiment-projects/{project_id}/samples/import", status_code=201)
async def import_experiment_samples(
    project_id: str,
    files: list[UploadFile] = File(...),
    batch_id: str | None = Form(None),
    expected_revision: int = Form(..., ge=1),
) -> dict[str, Any]:
    if not files or len(files) > 128:
        raise HTTPException(
            422,
            {"code": "FCS_UPLOAD_COUNT_INVALID", "message": "Upload between 1 and 128 FCS files"},
        )
    if EXPERIMENTS.get(project_id)["status"] != "COMPLETE":
        raise HTTPException(
            409,
            {
                "code": "INTEGRITY_ERROR",
                "message": "Project snapshot failed integrity verification",
            },
        )
    staging = EXPERIMENTS.projects_root / project_id / ".uploads"
    staging.mkdir(parents=True, exist_ok=True)
    imported = []
    temporaries = []
    try:
        for upload in files:
            original = (upload.filename or "sample.fcs").replace("\\", "/").split("/")[-1]
            if Path(original).suffix.lower() != ".fcs":
                raise HTTPException(
                    415,
                    {
                        "code": "UNSUPPORTED_SAMPLE_FILE",
                        "message": f"Only FCS files are accepted here: {original}",
                    },
                )
            target = staging / f"{uuid.uuid4().hex}.fcs"
            temporaries.append(target)
            await _save_upload(upload, target)
            metadata = await asyncio.to_thread(_inspect_uploaded_fcs, target)
            imported.append(
                {"path": str(target), "original_filename": original, "metadata": metadata}
            )
        result = EXPERIMENTS.import_fcs_files(
            project_id, imported, batch_id=batch_id, expected_revision=expected_revision
        )
    except (KeyError, IntegrityError, RevisionConflict, ValueError) as exc:
        raise _experiment_write_error(exc) from exc
    finally:
        for path in temporaries:
            path.unlink(missing_ok=True)
        for upload in files:
            await upload.close()
    return _experiment_response(result)


@app.post("/api/v4/experiment-projects/{project_id}/samples/import-flowjo", status_code=201)
async def import_flowjo_workspace(
    project_id: str,
    workspace: UploadFile = File(...),
    fcs_artifact_ids: str = Form(...),
    expected_revision: int = Form(..., ge=1),
) -> dict[str, Any]:
    original = (workspace.filename or "workspace.wsp").replace("\\", "/").split("/")[-1]
    if Path(original).suffix.lower() != ".wsp":
        raise HTTPException(
            415, {"code": "UNSUPPORTED_WORKSPACE_FILE", "message": "FlowJo workspace must use .wsp"}
        )
    try:
        artifact_ids = json.loads(fcs_artifact_ids)
    except json.JSONDecodeError as exc:
        raise HTTPException(
            422, {"code": "FCS_ARTIFACTS_INVALID", "message": "fcs_artifact_ids must be JSON array"}
        ) from exc
    if not isinstance(artifact_ids, list) or not artifact_ids or len(artifact_ids) > 128:
        raise HTTPException(
            422,
            {
                "code": "FCS_ARTIFACTS_INVALID",
                "message": "Select 1 to 128 registered FCS artifacts",
            },
        )
    paths = []
    try:
        for artifact_id in artifact_ids:
            path, artifact = EXPERIMENTS.artifact_path(project_id, str(artifact_id))
            if artifact["kind"] != "raw_fcs":
                raise ValueError("Only registered raw FCS artifacts can be matched to a workspace")
            paths.append(path)
    except (KeyError, IntegrityError, ValueError) as exc:
        raise _experiment_write_error(exc) from exc
    staging = EXPERIMENTS.projects_root / project_id / ".uploads"
    staging.mkdir(parents=True, exist_ok=True)
    temp = staging / f"{uuid.uuid4().hex}.wsp"
    try:
        await _save_upload(workspace, temp, maximum=80 * 1024 * 1024)
        preview = await asyncio.to_thread(_import_uploaded_wsp, temp, paths)
        preview["flowjo_input_artifact_ids"] = artifact_ids
        result = EXPERIMENTS.register_flowjo_import(
            project_id, temp, original, preview, expected_revision
        )
    except (KeyError, IntegrityError, RevisionConflict, ValueError) as exc:
        raise _experiment_write_error(exc) from exc
    finally:
        temp.unlink(missing_ok=True)
        await workspace.close()
    return _experiment_response(result)


@app.get("/api/v4/experiment-projects/{project_id}/flowjo-percentage-imports")
def list_flowjo_percentage_imports(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_flowjo_percentage_imports(project_id))
    except (KeyError, IntegrityError) as exc:
        raise _experiment_write_error(exc) from exc


@app.get(
    "/api/v4/experiment-projects/{project_id}/flowjo-percentage-imports/{import_id}/qc-worksheet"
)
def download_flowjo_percentage_qc_worksheet(project_id: str, import_id: str) -> Response:
    try:
        payload = EXPERIMENTS.build_flowjo_percentage_qc_worksheet(project_id, import_id)
    except (KeyError, IntegrityError, ValueError) as exc:
        raise _experiment_write_error(exc) from exc
    return Response(
        content=payload,
        media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={
            "Content-Disposition": f'attachment; filename="flowjo-qc-review-{import_id[:12]}.xlsx"'
        },
    )


@app.post(
    "/api/v4/experiment-projects/{project_id}/flowjo-percentage-imports/{import_id}/qc-worksheet",
    status_code=201,
)
async def import_flowjo_percentage_qc_worksheet(
    project_id: str,
    import_id: str,
    worksheet: UploadFile = File(...),
) -> dict[str, Any]:
    original = (worksheet.filename or "qc-review.xlsx").replace("\\", "/").split("/")[-1]
    if Path(original).suffix.lower() != ".xlsx":
        await worksheet.close()
        raise HTTPException(415, {
            "code": "QC_WORKSHEET_FORMAT_UNSUPPORTED",
            "message": "Upload the downloaded .xlsx QC review worksheet",
        })
    try:
        payload = await worksheet.read(2 * 1024 * 1024 + 1)
        result = EXPERIMENTS.import_flowjo_percentage_qc_worksheet(
            project_id, import_id, payload
        )
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        raise HTTPException(422, {
            "code": "QC_WORKSHEET_INVALID", "message": str(exc),
        }) from exc
    finally:
        await worksheet.close()
    return _experiment_response(result)


@app.get(
    "/api/v4/experiment-projects/{project_id}/flowjo-percentage-imports/{import_id}/effect-analyses"
)
def list_flowjo_percentage_effect_analyses(
    project_id: str, import_id: str
) -> dict[str, Any]:
    try:
        return _experiment_response(
            EXPERIMENTS.list_flowjo_percentage_effect_analyses(project_id, import_id)
        )
    except (KeyError, IntegrityError) as exc:
        raise _experiment_write_error(exc) from exc


@app.post("/api/v4/experiment-projects/{project_id}/flowjo-percentage-imports", status_code=201)
async def import_flowjo_percentage_table(
    project_id: str,
    table: UploadFile = File(...),
    batch_id: str = Form(...),
    numeric_scale: str = Form("AUTO"),
    decimal_separator: str = Form("."),
    compensation_status: str = Form("UNKNOWN"),
    compensation_note: str = Form(""),
    expected_revision: int = Form(..., ge=1),
) -> dict[str, Any]:
    original = (table.filename or "flowjo-results.csv").replace("\\", "/").split("/")[-1]
    suffix = Path(original).suffix.lower()
    if suffix not in {".xls", ".xlsx", ".csv", ".tsv", ".txt"}:
        raise HTTPException(
            415,
            {"code": "UNSUPPORTED_FLOWJO_RESULT_TABLE", "message": "上传已圈门的 Excel/CSV/TSV 结果表"},
        )
    if numeric_scale not in {"AUTO", "PERCENT", "FRACTION"}:
        raise HTTPException(422, {
            "code": "FLOWJO_RESULT_SCALE_INVALID",
            "message": "numeric_scale must be AUTO, PERCENT or FRACTION",
        })
    if decimal_separator not in {".", ","}:
        raise HTTPException(422, {
            "code": "FLOWJO_DECIMAL_SEPARATOR_INVALID",
            "message": "decimal_separator must be '.' or ','",
        })
    if compensation_status not in {"NOT_APPLIED", "APPLIED", "UNKNOWN"}:
        raise HTTPException(422, {
            "code": "FLOWJO_COMPENSATION_STATUS_INVALID",
            "message": "compensation_status must be NOT_APPLIED, APPLIED or UNKNOWN",
        })
    if len(compensation_note) > 500:
        raise HTTPException(422, {
            "code": "FLOWJO_COMPENSATION_NOTE_TOO_LONG",
            "message": "compensation_note must be 500 characters or fewer",
        })
    staging = EXPERIMENTS.projects_root / project_id / ".uploads"
    staging.mkdir(parents=True, exist_ok=True)
    temp = staging / f"{uuid.uuid4().hex}{suffix}"
    try:
        await _save_upload(table, temp, maximum=8 * 1024 * 1024)
        parsed = await asyncio.to_thread(
            parse_flowjo_percentage_table,
            temp.read_bytes(),
            original,
            numeric_scale=numeric_scale,
            decimal_separator=decimal_separator,
        )
        result = EXPERIMENTS.register_flowjo_percentage_import(
            project_id, batch_id=batch_id, source_path=temp, parsed_table=parsed,
            compensation_status=compensation_status,
            compensation_note=compensation_note,
            expected_revision=expected_revision,
        )
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except (FlowJoTableError, ValueError) as exc:
        raise HTTPException(422, {
            "code": "FLOWJO_RESULT_TABLE_INVALID", "message": str(exc),
        }) from exc
    finally:
        temp.unlink(missing_ok=True)
        await table.close()
    return _experiment_response(result)


@app.post("/api/v4/experiment-projects/{project_id}/samples/confirm", status_code=201)
def confirm_experiment_samples(project_id: str, payload: SampleConfirmInput) -> dict[str, Any]:
    try:
        result = EXPERIMENTS.confirm_fcs_matches(
            project_id,
            payload.upload_session_id,
            [item.model_dump() for item in payload.matches],
            payload.expected_revision,
        )
    except (KeyError, IntegrityError, RevisionConflict, ValueError) as exc:
        raise _experiment_write_error(exc) from exc
    return _experiment_response(result)


@app.get("/api/v4/experiment-projects/{project_id}/samples/import-sessions/{session_id}")
def get_sample_import_session(project_id: str, session_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.get_import_session(project_id, session_id))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc


@app.get("/api/v4/experiment-projects/{project_id}/samples")
def list_experiment_samples(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_entities(project_id, "samples"))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc


@app.post("/api/v4/experiment-projects/{project_id}/samples/{sample_id}/plot-data")
def sample_plot_data(
    project_id: str, sample_id: str, payload: SamplePlotRequest
) -> dict[str, Any]:
    try:
        sample = next(
            item
            for item in EXPERIMENTS.list_entities(project_id, "samples")
            if item["id"] == sample_id
        )
        path, artifact = EXPERIMENTS.artifact_path(project_id, sample["artifact_id"])
        if artifact["kind"] != "raw_fcs":
            raise ValueError("Sample is not linked to a registered raw FCS artifact")
        labels = {
            str(value)
            for channel in sample.get("fcs_metadata", {}).get("channels", [])
            for value in (channel.get("pnn"), channel.get("pns"))
            if value
        }
        if payload.x_channel not in labels or payload.y_channel not in labels:
            raise ValueError("Selected scatter channel is not present in this FCS sample")
    except StopIteration as exc:
        raise HTTPException(
            404, {"code": "SAMPLE_NOT_FOUND", "message": "Sample not found"}
        ) from exc
    except (KeyError, IntegrityError) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        raise HTTPException(
            422, {"code": "PLOT_DATA_INVALID", "message": str(exc)}
        ) from exc

    python = ROOT / ".tool_envs" / "flowkit" / "Scripts" / "python.exe"
    worker = ROOT / "integrations" / "workers" / "flowkit_worker.py"
    if not python.is_file() or not worker.is_file():
        raise HTTPException(
            503, {"code": "FLOWKIT_WORKER_UNAVAILABLE", "message": "FlowKit worker is unavailable"}
        )
    request = {
        "schema_version": "1.0",
        "request_id": uuid.uuid4().hex,
        "operation": "sample_events",
        "path": str(path),
        "x_channel": payload.x_channel,
        "y_channel": payload.y_channel,
        "maximum_points": payload.maximum_points,
    }
    try:
        completed = subprocess.run(
            [str(python), str(worker)], input=json.dumps(request), capture_output=True,
            text=True, encoding="utf-8", timeout=60, check=False,
            env=os.environ | {"PYTHONIOENCODING": "utf-8"},
        )
        worker_response = json.loads(completed.stdout)
    except OSError as exc:
        raise HTTPException(
            503, {"code": "FLOWKIT_WORKER_UNAVAILABLE", "message": str(exc)}
        ) from exc
    except (subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
        raise HTTPException(
            422, {"code": "PLOT_DATA_FAILED", "message": str(exc)}
        ) from exc
    if completed.returncode or worker_response.get("status") != "COMPLETE":
        errors = worker_response.get("errors") or [{"message": completed.stderr[-1000:]}]
        raise HTTPException(422, {
            "code": "PLOT_DATA_FAILED",
            "message": errors[0].get("message", "FCS display sampling failed"),
        })
    records = worker_response.get("records")
    if not isinstance(records, dict) or records.get("display_only") is not True:
        raise HTTPException(422, {
            "code": "PLOT_DATA_FAILED",
            "message": "FlowKit returned an invalid display-only plot payload",
        })
    return _experiment_response(records | {
        "sample": {
            "sample_id": sample_id,
            "label": sample.get("condition", {}).get("sample_label") or sample_id,
            "role": sample["role"],
        },
        "provenance": worker_response.get("provenance", {}),
    })


@app.get("/api/v4/experiment-projects/{project_id}/gates")
def list_gate_revisions(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_gate_revisions(project_id))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc


@app.post(
    "/api/v4/experiment-projects/{project_id}/compensation-matrices", status_code=201
)
async def upload_compensation_matrix(
    project_id: str,
    matrix: UploadFile = File(...),
    expected_revision: int = Form(..., ge=1),
) -> dict[str, Any]:
    original = (matrix.filename or "").replace("\\", "/").split("/")[-1]
    try:
        content = await matrix.read(1024 * 1024 + 1)
        matrix_document = parse_compensation_matrix(content, original)
        result = EXPERIMENTS.register_compensation_matrix(
            project_id, matrix_document, expected_revision
        )
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except (GateError, ValueError) as exc:
        raise HTTPException(
            422, {"code": "COMPENSATION_MATRIX_INVALID", "message": str(exc)}
        ) from exc
    finally:
        await matrix.close()
    return _experiment_response(result)


@app.post("/api/v4/experiment-projects/{project_id}/gates", status_code=201)
def create_gate_revision(project_id: str, payload: GateRevisionCreate) -> dict[str, Any]:
    try:
        result = EXPERIMENTS.save_gate_revision(
            project_id,
            name=payload.name,
            gates=payload.gate_tree,
            transforms=payload.transforms,
            compensation_decision=payload.compensation_decision,
            compensation_artifact_id=payload.compensation_artifact_id,
            parent_revision_id=payload.parent_gate_revision_id,
            flowjo_rebuild=payload.flowjo_rebuild.model_dump() if payload.flowjo_rebuild else None,
            expected_revision=payload.expected_revision,
        )
    except (KeyError, IntegrityError, RevisionConflict, GateError, ValueError) as exc:
        raise _experiment_write_error(exc) from exc
    return _experiment_response(result)


@app.get("/api/v4/experiment-projects/{project_id}/gate-analyses")
def list_gate_analyses(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_gate_analyses(project_id))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc


@app.get("/api/v4/experiment-projects/{project_id}/effect-analyses")
def list_effect_analyses(project_id: str) -> dict[str, Any]:
    try:
        return _experiment_response(EXPERIMENTS.list_effect_analyses(project_id))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc


@app.post("/api/v4/experiment-projects/{project_id}/effect-analyses", status_code=201)
def create_effect_analysis(project_id: str, payload: EffectAnalysisCreate) -> dict[str, Any]:
    try:
        project = EXPERIMENTS.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        gate_analysis = next((item for item in EXPERIMENTS.list_gate_analyses(project_id)
                              if item["id"] == payload.gate_analysis_id), None)
        if gate_analysis is None or gate_analysis["status"] != "COMPLETE":
            raise KeyError(payload.gate_analysis_id)
        result_path, _ = EXPERIMENTS.artifact_path(project_id, gate_analysis["result_artifact_id"])
        source_document = read_json(result_path)
        sample_result_by_id = {
            str(item.get("sample", {}).get("sample_id")): item
            for item in source_document.get("samples", [])
        }
        sample_records = EXPERIMENTS.list_entities(project_id, "samples")
        sample_by_id = {item["id"]: item for item in sample_records}
        batches = {item["id"]: item for item in EXPERIMENTS.list_entities(project_id, "assay_batches")}
        constructs = {item["id"]: item for item in EXPERIMENTS.list_entities(project_id, "antigen_constructs")}
        antibodies = {item["id"]: item for item in EXPERIMENTS.list_entities(project_id, "antibody_constructs")}
        review_by_sample: dict[str, dict[str, Any]] = {}
        for review in EXPERIMENTS.list_qc_reviews(project_id, "GATE_ANALYSIS", payload.gate_analysis_id):
            previous = review_by_sample.get(review["sample_id"])
            if previous is None or review["version"] > previous["version"]:
                review_by_sample[review["sample_id"]] = review

        if len({item.sample_id for item in payload.samples}) != len(payload.samples):
            raise ValueError("Each sample may be mapped only once")
        group_rows: dict[tuple[str, str, str], list[WellCount]] = {}
        well_records = []
        probe_constructs: dict[str, dict[str, Any]] = {}
        shared_base_stratum = None
        shared_antibody_id = None
        selected_antibody_ids: set[str] = set()
        for item in payload.samples:
            sample = sample_by_id.get(item.sample_id)
            result_sample = sample_result_by_id.get(item.sample_id)
            if sample is None or result_sample is None:
                raise ValueError(f"Sample {item.sample_id} is not part of the selected gate analysis")
            if item.probe_id not in constructs:
                raise ValueError(f"Unknown antigen probe construct {item.probe_id}")
            probe = constructs[item.probe_id]
            probe_constructs[item.probe_id] = probe
            if payload.comparison_mode == "ANTIBODY_EFFECT":
                if item.condition not in {"REFERENCE", "REFERENCE_NC"}:
                    raise ValueError("ANTIBODY_EFFECT accepts only WT antigen activation and each antibody's own NC")
                if probe["sequence"] != project["manifest"]["reference_sequence"]:
                    raise ValueError("ANTIBODY_EFFECT requires one fixed WT antigen construct")
                if sample.get("antigen_construct_id") != item.probe_id:
                    raise ValueError("ANTIBODY_EFFECT samples must use the selected fixed WT antigen construct")
            if item.condition in {"MUTANT", "MUTANT_NC"} and sample.get("antigen_construct_id") != item.probe_id:
                raise ValueError("Mutant antigen and NC wells must reference their probe construct")
            if item.condition in {"REFERENCE", "REFERENCE_NC"}:
                reference_construct = constructs.get(sample.get("antigen_construct_id"))
                if reference_construct is None or reference_construct["sequence"] != project["manifest"]["reference_sequence"]:
                    raise ValueError("Reference antigen and NC wells must use a registered WT construct")
            is_nc = item.condition.endswith("_NC")
            if (sample["role"] == "NC") != is_nc:
                raise ValueError("NC conditions must use NC samples; activation conditions must use ACTIVATION samples")
            if not is_nc and sample["role"] != "ACTIVATION":
                raise ValueError("Only ACTIVATION and NC roles can enter activation effect analysis")
            batch = batches.get(sample["batch_id"])
            if batch is None or batch["assay"] != "ACTIVATION":
                raise ValueError("Effect analysis requires samples from ACTIVATION assay batches")
            if not batch["independent_batch"]:
                raise ValueError("Effect analysis requires batches marked as independent")
            antibody_id = sample.get("antibody_construct_id")
            if antibody_id not in antibodies:
                raise ValueError("Every activation comparison requires a registered antibody construct")
            selected_antibody_ids.add(antibody_id)
            base_stratum = {
                "assay": batch["assay"],
                "cell_line": batch["cell_line"], "dose": batch["dose"],
                "dose_unit": batch["dose_unit"],
                "gate_revision_id": source_document["gate_revision_id"],
                "positive_gate_id": payload.positive_gate_id,
            }
            serialized_base_stratum = json.dumps(base_stratum, sort_keys=True, separators=(",", ":"))
            if shared_base_stratum is None:
                shared_base_stratum = serialized_base_stratum
            elif serialized_base_stratum != shared_base_stratum:
                raise ValueError("Activation wells differ in assay, dose, cell line or gate threshold")
            if payload.comparison_mode == "ANTIGEN_EFFECT":
                if shared_antibody_id is None:
                    shared_antibody_id = antibody_id
                elif antibody_id != shared_antibody_id:
                    raise ValueError("ANTIGEN_EFFECT requires one fixed antibody construct")
            stratum = json.dumps(
                base_stratum | {"antibody_construct_id": antibody_id},
                sort_keys=True, separators=(",", ":"),
            )
            gate_count = next((row for row in result_sample.get("gate_counts", [])
                               if row.get("gate_id") == payload.positive_gate_id), None)
            if gate_count is None:
                raise ValueError(f"Positive gate {payload.positive_gate_id} is missing for sample {item.sample_id}")
            review = review_by_sample.get(item.sample_id)
            if review is None:
                raise ValueError(f"QC_REVIEW_REQUIRED: sample {item.sample_id} has not been manually reviewed")
            positive, total = gate_count.get("count"), gate_count.get("parent_count")
            if not isinstance(positive, int) or not isinstance(total, int) or total < 1:
                raise ValueError(f"Positive gate has invalid counts for sample {item.sample_id}")
            well = WellCount(
                batch_id=sample["batch_id"], stratum_id=stratum,
                condition=item.condition, positive=positive, total=total,
                qc_status=review["status"],
            )
            key = (item.probe_id, sample["batch_id"], stratum)
            group_rows.setdefault(key, []).append(well)
            well_records.append({
                "sample_id": item.sample_id, "probe_id": item.probe_id,
                "batch_id": sample["batch_id"], "condition": item.condition,
                "positive_gate_id": payload.positive_gate_id,
                "positive_count": positive, "parent_count": total,
                "fraction": well.fraction, "qc_status": review["status"],
                "qc_review_id": review["id"], "qc_reason": review["reason"],
                "replicate_type": sample["replicate_type"],
                "replicate_index": sample["replicate_index"],
                "antigen_construct_id": sample.get("antigen_construct_id"),
                "antibody_construct_id": antibody_id,
            })

        factorial_parent_id = None
        factorial_variant_id = None
        if payload.comparison_mode == "DOUBLE_MUTANT_INTERACTION":
            if len(selected_antibody_ids) != 2:
                raise ValueError("DOUBLE_MUTANT_INTERACTION requires one parent and one variant antibody")
            variants = [
                antibodies[construct_id]
                for construct_id in selected_antibody_ids
                if antibodies[construct_id]["role"] == "VARIANT"
            ]
            if len(variants) != 1:
                raise ValueError("DOUBLE_MUTANT_INTERACTION requires exactly one antibody variant")
            factorial_variant_id = variants[0]["id"]
            factorial_parent_id = variants[0]["parent_id"]
            if (
                factorial_parent_id not in selected_antibody_ids
                or antibodies.get(factorial_parent_id, {}).get("role") != "PARENT"
            ):
                raise ValueError("The selected variant must be paired with its registered parent antibody")
            if any(
                not any(mutation.get("kind") == "SUBSTITUTION" for mutation in probe["mutations"])
                for probe in probe_constructs.values()
            ):
                raise ValueError("DOUBLE_MUTANT_INTERACTION requires antigen variant probes")

        antibody_effect_parent_id = None
        antibody_effect_variant_id = None
        if payload.comparison_mode == "ANTIBODY_EFFECT":
            if len(probe_constructs) != 1:
                raise ValueError("ANTIBODY_EFFECT requires exactly one fixed WT antigen construct")
            if len(selected_antibody_ids) != 2:
                raise ValueError("ANTIBODY_EFFECT requires one parent and one variant antibody")
            parents = [
                antibodies[construct_id] for construct_id in selected_antibody_ids
                if antibodies[construct_id]["role"] == "PARENT"
            ]
            variants = [
                antibodies[construct_id] for construct_id in selected_antibody_ids
                if antibodies[construct_id]["role"] == "VARIANT"
            ]
            if len(parents) != 1 or len(variants) != 1:
                raise ValueError("ANTIBODY_EFFECT requires one parent and one variant antibody")
            antibody_effect_parent_id = parents[0]["id"]
            antibody_effect_variant_id = variants[0]["id"]
            if variants[0]["parent_id"] != antibody_effect_parent_id:
                raise ValueError("The selected variant must be paired with its registered parent antibody")

        batch_effects = []
        projection_rows = []
        batch_effect_by_factor: dict[tuple[str, str, str], Any] = {}
        antibody_effect_wells: dict[tuple[str, str, str], list[WellCount]] = {}
        for (probe_id, batch_id, stratum), wells in sorted(group_rows.items()):
            if payload.comparison_mode == "ANTIBODY_EFFECT":
                antibody_id = json.loads(stratum)["antibody_construct_id"]
                antibody_effect_wells[(probe_id, batch_id, antibody_id)] = wells
                continue
            condition_names = {well.condition for well in wells}
            mutant_nc = "MUTANT_NC" if "MUTANT_NC" in condition_names else None
            if mutant_nc is None:
                raise ValueError(f"Probe {probe_id}, batch {batch_id} is missing MUTANT_NC")
            effect = calculate_batch_effect(
                wells, mutant="MUTANT", reference="REFERENCE", mutant_nc=mutant_nc,
                reference_nc="REFERENCE_NC" if "REFERENCE_NC" in condition_names else None,
                min_reference_window=payload.min_reference_window,
            )
            positions = tuple(sorted({
                int(mutation["position"])
                for mutation in probe_constructs[probe_id]["mutations"]
                if mutation.get("kind") == "SUBSTITUTION"
            }))
            antibody_id = json.loads(stratum)["antibody_construct_id"]
            batch_effect_by_factor[(probe_id, batch_id, antibody_id)] = effect
            batch_effects.append({
                "probe_id": probe_id, "mutation_positions": positions,
                "antibody_construct_id": antibody_id,
                "batch_effect": asdict(effect),
            })
            if positions and payload.comparison_mode == "ANTIGEN_EFFECT":
                projection_rows.append(ProbeBatchEffect(probe_id, positions, effect))

        antibody_effects = []
        antibody_summaries = []
        if payload.comparison_mode == "ANTIBODY_EFFECT":
            fixed_probe_id = next(iter(probe_constructs))
            effect_batch_ids = sorted({batch_id for _, batch_id, _ in antibody_effect_wells})
            if not effect_batch_ids:
                raise ValueError("ANTIBODY_EFFECT requires complete parent and variant groups")
            for batch_id in effect_batch_ids:
                parent_wells = antibody_effect_wells.get(
                    (fixed_probe_id, batch_id, antibody_effect_parent_id)
                )
                variant_wells = antibody_effect_wells.get(
                    (fixed_probe_id, batch_id, antibody_effect_variant_id)
                )
                if parent_wells is None or variant_wells is None:
                    raise ValueError(
                        f"Batch {batch_id} must contain parent and variant WT activation/NC groups"
                    )
                antibody_effects.append({
                    "probe_id": fixed_probe_id,
                    "parent_antibody_construct_id": antibody_effect_parent_id,
                    "variant_antibody_construct_id": antibody_effect_variant_id,
                    **calculate_antibody_variant_batch_effect(
                        parent_wells, variant_wells,
                        min_reference_window=payload.min_reference_window,
                    ),
                })
            antibody_summaries = [
                {
                    "probe_id": fixed_probe_id,
                    "parent_antibody_construct_id": antibody_effect_parent_id,
                    "variant_antibody_construct_id": antibody_effect_variant_id,
                    **summarize_antibody_variant_effects(
                        antibody_effects, seed=payload.bootstrap_seed,
                        bootstrap_iterations=payload.bootstrap_iterations,
                    ),
                }
            ]

        factorial_interactions = []
        factorial_summaries = []
        if payload.comparison_mode == "DOUBLE_MUTANT_INTERACTION":
            by_probe: dict[str, list[dict[str, Any]]] = {}
            probe_batch_ids = sorted({(probe_id, batch_id) for probe_id, batch_id, _ in batch_effect_by_factor})
            for probe_id, batch_id in probe_batch_ids:
                parent_effect = batch_effect_by_factor.get((probe_id, batch_id, factorial_parent_id))
                variant_effect = batch_effect_by_factor.get((probe_id, batch_id, factorial_variant_id))
                if parent_effect is None or variant_effect is None:
                    raise ValueError(
                        f"Probe {probe_id}, batch {batch_id} lacks one of the four antibody-by-antigen groups"
                    )
                interaction = calculate_factorial_interaction(parent_effect, variant_effect)
                positions = tuple(sorted({
                    int(mutation["position"])
                    for mutation in probe_constructs[probe_id]["mutations"]
                    if mutation.get("kind") == "SUBSTITUTION"
                }))
                interaction_record = {
                    "probe_id": probe_id,
                    "mutation_positions": positions,
                    "parent_antibody_construct_id": factorial_parent_id,
                    "variant_antibody_construct_id": factorial_variant_id,
                    **interaction,
                }
                factorial_interactions.append(interaction_record)
                by_probe.setdefault(probe_id, []).append(interaction)
            factorial_summaries = [
                {
                    "probe_id": probe_id,
                    "mutation_positions": sorted({
                        int(mutation["position"])
                        for mutation in probe_constructs[probe_id]["mutations"]
                        if mutation.get("kind") == "SUBSTITUTION"
                    }),
                    **summarize_factorial_interactions(
                        rows, seed=payload.bootstrap_seed,
                        bootstrap_iterations=payload.bootstrap_iterations,
                    ),
                }
                for probe_id, rows in sorted(by_probe.items())
            ]

        graph_path = EXPERIMENTS.projects_root / project_id / "reference" / "surface_graph_v4.json"
        graph = read_json(graph_path)
        sectors: dict[str, set[int]] = {}
        graph_positions = set()
        for node in graph.get("nodes", []):
            graph_positions.add(int(node["uniprot_pos"]))
            if node.get("sector_id") is not None:
                sectors.setdefault(str(node["sector_id"]), set()).add(int(node["uniprot_pos"]))
        projection = (
            project_probe_effects(
                projection_rows, sectors, seed=payload.bootstrap_seed,
                bootstrap_iterations=payload.bootstrap_iterations,
            )
            if projection_rows else None
        )
        if payload.comparison_mode == "ANTIGEN_EFFECT":
            unmapped_probe_ids = {
                item.probe_id for item in projection.probes
                if item.attribution == "UNMAPPED_OR_NOT_SURFACE"
            } if projection else set()
            unmapped_probe_ids.update(set(probe_constructs) - {row.probe_id for row in projection_rows})
        elif payload.comparison_mode == "DOUBLE_MUTANT_INTERACTION":
            unmapped_probe_ids = {
                probe_id for probe_id, probe in probe_constructs.items()
                if any(
                    mutation.get("kind") == "SUBSTITUTION"
                    and int(mutation["position"]) not in graph_positions
                    for mutation in probe["mutations"]
                )
            }
        else:
            unmapped_probe_ids = set()
        document = {
            "comparison_mode": payload.comparison_mode,
            "gate_analysis_sha256": next(
                artifact["sha256"] for artifact in project["artifacts"]
                if artifact["id"] == gate_analysis["result_artifact_id"]
            ),
            "gate_revision_id": source_document["gate_revision_id"],
            "gate_revision_version": source_document["gate_revision_version"],
            "frozen_hashes": project["frozen_hashes"],
            "parameters": {
                "min_reference_window": payload.min_reference_window,
                "bootstrap_seed": payload.bootstrap_seed,
                "bootstrap_iterations": payload.bootstrap_iterations,
            },
            "well_results": well_records,
            "batch_effects": batch_effects,
            "probe_mutations": {
                probe_id: probe["mutations"] for probe_id, probe in sorted(probe_constructs.items())
            },
            "surface_projection": asdict(projection) if projection else None,
            "factorial_interactions": factorial_interactions,
            "factorial_summaries": factorial_summaries,
            "antibody_effects": antibody_effects,
            "antibody_summaries": antibody_summaries,
            "unmapped_probe_ids": sorted(unmapped_probe_ids),
            "limitations": [
                "Positive rates are condition-specific functional readouts, not affinity or direct-contact probabilities.",
                *(["DOUBLE_MUTANT_INTERACTION is a descriptive coupling clue and is not assigned to individual antigen residues."]
                  if payload.comparison_mode == "DOUBLE_MUTANT_INTERACTION" else []),
                *( ["ANTIBODY_EFFECT is a functional WT-antigen readout and does not establish affinity, folding, surface expression or spatial antigen-residue attribution."]
                  if payload.comparison_mode == "ANTIBODY_EFFECT" else []),
            ],
        }
        result = EXPERIMENTS.save_effect_analysis(
            project_id, gate_analysis_id=payload.gate_analysis_id,
            positive_gate_id=payload.positive_gate_id, document=document,
            expected_revision=payload.expected_revision,
        )
        return _experiment_response(result | {
            "comparison_mode": payload.comparison_mode,
            "well_results": well_records, "batch_effects": batch_effects,
            "surface_projection": document["surface_projection"],
            "factorial_interactions": factorial_interactions,
            "factorial_summaries": factorial_summaries,
            "antibody_effects": antibody_effects,
            "antibody_summaries": antibody_summaries,
            "limitations": document["limitations"],
        })
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        code = "QC_HOLD" if "QC_REVIEW_REQUIRED" in str(exc) else "EFFECT_ANALYSIS_INVALID"
        raise HTTPException(422, {"code": code, "message": str(exc)}) from exc


@app.get("/api/v4/experiment-projects/{project_id}/qc-reviews")
def list_qc_reviews(
    project_id: str, source_type: str | None = None, source_id: str | None = None
) -> dict[str, Any]:
    try:
        if bool(source_type) != bool(source_id):
            raise ValueError("source_type and source_id must be supplied together")
        return _experiment_response(EXPERIMENTS.list_qc_reviews(project_id, source_type, source_id))
    except KeyError as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        raise HTTPException(422, {"code": "QC_REVIEW_SOURCE_INVALID", "message": str(exc)}) from exc


@app.post(
    "/api/v4/experiment-projects/{project_id}/flowjo-percentage-imports/{import_id}/effect-analyses",
    status_code=201,
)
def create_flowjo_percentage_effect_analysis(
    project_id: str,
    import_id: str,
    payload: FlowJoPercentageEffectCreate,
) -> dict[str, Any]:
    try:
        project = EXPERIMENTS.get(project_id)
        if project["status"] != "COMPLETE":
            raise IntegrityError("Project snapshot failed integrity verification")
        imported = next(
            (item for item in EXPERIMENTS.list_flowjo_percentage_imports(project_id)
             if item["id"] == import_id),
            None,
        )
        if imported is None:
            raise KeyError(import_id)
        source = imported["result"]
        flowjo_import_id = source.get("flowjo_import_id")
        missing_sample_ids = source.get("missing_sample_ids") or []
        if missing_sample_ids:
            example_ids = ", ".join(str(item) for item in missing_sample_ids[:3])
            raise ValueError(
                f"QC_REVIEW_REQUIRED: result table omits {len(missing_sample_ids)} registered samples; "
                f"re-import a complete table before effect analysis (examples: {example_ids})"
            )
        qc_source_type = "FLOWJO_PERCENTAGE_IMPORT"
        qc_source_id = import_id
        reviews = EXPERIMENTS.list_qc_reviews(project_id, qc_source_type, qc_source_id)
        # Read older reviews created before percentage tables had their own QC source.
        # New reviews always attach to this exact table import, whether or not a WSP exists.
        if not reviews and flowjo_import_id:
            legacy_reviews = EXPERIMENTS.list_qc_reviews(
                project_id, "FLOWJO_IMPORT", flowjo_import_id
            )
            if legacy_reviews:
                qc_source_type = "FLOWJO_IMPORT"
                qc_source_id = flowjo_import_id
                reviews = legacy_reviews
        latest_review: dict[str, dict[str, Any]] = {}
        for review in reviews:
            old = latest_review.get(review["sample_id"])
            if old is None or review["version"] > old["version"]:
                latest_review[review["sample_id"]] = review

        observations = source.get("observations") or []
        sample_ids = {item.get("sample_id") for item in observations}
        missing_reviews = sorted(sample_id for sample_id in sample_ids if sample_id not in latest_review)
        if missing_reviews:
            raise ValueError(
                f"QC_REVIEW_REQUIRED: {len(missing_reviews)} 个导入样本尚未人工审核，"
                f"例如 {missing_reviews[0]}"
            )

        by_panel: dict[str, list[dict[str, Any]]] = {}
        for item in observations:
            if item.get("panel_role") == "ANTIBODY" and item.get("panel_label"):
                by_panel.setdefault(str(item["panel_label"]), []).append(item)
        if not by_panel:
            raise ValueError("结果表没有可用于固定抗体比较的已登记临时面板标签")

        effect_rows = []
        qc_statuses = {"PASS", "PASS_WITH_OVERRIDE", "HOLD", "EXCLUDE"}
        for panel_label, panel_rows in sorted(by_panel.items()):
            references = [
                item for item in panel_rows
                if item.get("sample_role") == "ACTIVATION" and item.get("antigen_role") == "REFERENCE"
            ]
            controls = [item for item in panel_rows if item.get("sample_role") == "NC"]
            variants: dict[str, list[dict[str, Any]]] = {}
            for item in panel_rows:
                if item.get("sample_role") == "ACTIVATION" and item.get("antigen_role") == "VARIANT":
                    label = str(item.get("antigen_label") or item.get("condition_label") or "")
                    if label:
                        variants.setdefault(label, []).append(item)
            if not references or not controls:
                raise ValueError(f"面板 {panel_label} 缺少同期 WT 或 NC 孔")

            for variant_label, mutant_rows in sorted(variants.items()):
                selected = references + controls + mutant_rows
                malformed = [item for item in selected if item["sample_id"] not in latest_review]
                if malformed:
                    raise ValueError("QC_REVIEW_REQUIRED: 效应输入有样本未审核")
                replicate_types = {item.get("replicate_type") for item in selected}
                if len(replicate_types) != 1 or not replicate_types <= {"SINGLE", "TECHNICAL", "BIOLOGICAL"}:
                    raise ValueError(
                        f"面板 {panel_label} 的 WT、NC 和变体必须使用同一种有效重复类型"
                    )
                replicate_type = next(iter(replicate_types))
                if any(latest_review[item["sample_id"]]["status"] not in qc_statuses for item in selected):
                    raise ValueError("效应输入包含未知 QC 状态")
                stratum_id = f"{source['batch_id']}|{panel_label}|{source['gate_label']}"
                wells = []
                for item in mutant_rows:
                    wells.append(WellRate(
                        batch_id=source["batch_id"], stratum_id=stratum_id,
                        condition="MUTANT", fraction=float(item["fraction"]),
                        qc_status=latest_review[item["sample_id"]]["status"],
                    ))
                for item in references:
                    wells.append(WellRate(
                        batch_id=source["batch_id"], stratum_id=stratum_id,
                        condition="REFERENCE", fraction=float(item["fraction"]),
                        qc_status=latest_review[item["sample_id"]]["status"],
                    ))
                for item in controls:
                    wells.append(WellRate(
                        batch_id=source["batch_id"], stratum_id=stratum_id,
                        condition="MUTANT_NC", fraction=float(item["fraction"]),
                        qc_status=latest_review[item["sample_id"]]["status"],
                    ))
                try:
                    calculated = calculate_batch_effect(
                        wells, mutant="MUTANT", reference="REFERENCE", mutant_nc="MUTANT_NC",
                        min_reference_window=payload.min_reference_window,
                    )
                    effect_data = asdict(calculated)
                    analysis_status = "CALCULATED"
                except ValueError as exc:
                    effect_data = None
                    analysis_status = "NO_QC_APPROVED_DATA"
                    exclusion_note = str(exc)
                included = [
                    item for item in selected
                    if latest_review[item["sample_id"]]["status"] in {"PASS", "PASS_WITH_OVERRIDE"}
                ]
                comparison_group_by_sample = {
                    item["sample_id"]: group
                    for group, group_items in (
                        ("MUTANT", mutant_rows),
                        ("REFERENCE", references),
                        ("MUTANT_NC", controls),
                    )
                    for item in group_items
                }
                replicate_rows = [{
                    "sample_id": item["sample_id"], "well_id": item["well_id"],
                    "comparison_group": comparison_group_by_sample[item["sample_id"]],
                    "condition_label": item.get("condition_label"),
                    "replicate_type": item["replicate_type"],
                    "replicate_index": item["replicate_index"],
                    "fraction": item["fraction"], "percentage": item["percentage"],
                    "displayed_value": item.get("displayed_value"),
                    "qc_status": latest_review[item["sample_id"]]["status"],
                    "accepted_for_summary": latest_review[item["sample_id"]]["status"]
                    in {"PASS", "PASS_WITH_OVERRIDE"},
                    "qc_review_id": latest_review[item["sample_id"]]["id"],
                } for item in selected]
                row = {
                    "batch_id": source["batch_id"],
                    "panel_label": panel_label,
                    "panel_role": "ANTIBODY",
                    "probe_label": variant_label,
                    "reference_label": references[0].get("antigen_label") or references[0].get("condition_label"),
                    "replicate_type": replicate_type,
                    "effect_status": analysis_status,
                    "batch_effect": effect_data,
                    "measurement_count": len(selected),
                    "qc_accepted_count": len(included),
                    "replicate_count_by_group": {
                        "mutant": sum(row.condition == "MUTANT" for row in wells
                                      if row.qc_status in {"PASS", "PASS_WITH_OVERRIDE"}),
                        "reference": sum(row.condition == "REFERENCE" for row in wells
                                         if row.qc_status in {"PASS", "PASS_WITH_OVERRIDE"}),
                        "nc": sum(row.condition == "MUTANT_NC" for row in wells
                                  if row.qc_status in {"PASS", "PASS_WITH_OVERRIDE"}),
                    },
                    "replicate_observations": replicate_rows,
                    "exclusion_note": exclusion_note if analysis_status != "CALCULATED" else None,
                }
                effect_rows.append(row)

        if not effect_rows:
            raise ValueError("结果表没有可与同期 WT、NC 比较的抗原变体")
        replicate_types_used = sorted({item["replicate_type"] for item in effect_rows})
        compensation = source.get("compensation")
        if not isinstance(compensation, dict):
            compensation = {"status": "UNKNOWN", "source": "LEGACY_STATUS_UNRECORDED"}
        compensation_status = compensation.get("status")
        compensation_limitations = {
            "NOT_APPLIED": (
                "结果表按用户报告的未补偿输入生成；若通道串色影响门内分类，阳性率会受影响。"
            ),
            "APPLIED": (
                "结果表按用户报告的已补偿输入生成；仅凭百分比表无法独立验证所用补偿矩阵或设置。"
            ),
            "UNKNOWN": (
                "结果表所依据的补偿状态未知；效应仅描述导出百分比，不能假定不同补偿条件可直接比较。"
            ),
        }
        compensation_limitation = compensation_limitations.get(
            compensation_status, compensation_limitations["UNKNOWN"]
        )
        document = {
            "interpretation": "DESCRIPTIVE_FUNCTIONAL_READOUT",
            "comparison_basis": "within-batch, within-panel, equal-weight well medians with matched replicate type",
            "effects": effect_rows,
            "provenance": {
                "source_import_id": import_id,
                "flowjo_import_id": flowjo_import_id,
                "qc_source_type": qc_source_type,
                "qc_source_id": qc_source_id,
                "source_sha256": source["source_sha256"],
                "gate_label": source["gate_label"],
                "metric": source["metric"],
                "measurement_basis": "PERCENTAGE_ONLY",
                "denominator_available": False,
                "compensation": compensation,
                "replicate_types_used": replicate_types_used,
                "excluded_control_panels": [
                    {"panel_label": label, "panel_role": "POSITIVE_LIGAND_CONTROL",
                     "reason": "Positive ligand/control panel is not an antibody panel"}
                    for label in sorted({
                        str(item.get("panel_label"))
                        for item in observations
                        if item.get("panel_role") == "POSITIVE_LIGAND_CONTROL"
                        and item.get("panel_label")
                    })
                ],
                "construct_identities": "TEMPORARY_LABELS",
                "structural_attribution": "NOT_AVAILABLE",
                "independent_batch_count_per_panel": 1,
                "technical_replicates_are_not_biological_replicates": (
                    all(item["replicate_type"] == "TECHNICAL" for item in effect_rows)
                ),
                "confidence_intervals_or_significance": "NOT_CALCULATED",
            },
            "limitations": [
                "输入是 FlowJo 导出的百分比；事件数分子和分母未知，不推导二项分布误差或事件级置信区间。",
                "WT、NC 和变体必须使用一致的登记重复类型；中位数按孔等权计算，不跨不同抗体面板合并。",
                compensation_limitation,
                "E1–E7 和抗体身份仍是临时标签；结果只表示该门控条件下的功能阳性率变化。",
                "结果不代表结合亲和力、直接接触或结构位点效应。",
                "所有使用样本均需先完成人工 QC；HOLD/EXCLUDE 不参与条件中位数。",
            ],
        }
        saved = EXPERIMENTS.save_flowjo_percentage_effect_analysis(
            project_id, import_id=import_id, document=document,
            expected_revision=payload.expected_revision,
        )
        return _experiment_response(saved | {"effects": effect_rows})
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        code = "QC_HOLD" if "QC_REVIEW_REQUIRED" in str(exc) else "FLOWJO_PERCENTAGE_EFFECT_INVALID"
        raise HTTPException(422, {"code": code, "message": str(exc)}) from exc


@app.post("/api/v4/experiment-projects/{project_id}/qc-reviews", status_code=201)
def create_qc_review(project_id: str, payload: QCReviewCreate) -> dict[str, Any]:
    supplied_sources = [
        ("GATE_ANALYSIS", payload.gate_analysis_id),
        ("FLOWJO_IMPORT", payload.flowjo_import_id),
        ("FLOWJO_PERCENTAGE_IMPORT", payload.flowjo_percentage_import_id),
    ]
    selected_sources = [(kind, source_id) for kind, source_id in supplied_sources if source_id]
    if len(selected_sources) != 1:
        raise HTTPException(422, {
            "code": "QC_REVIEW_SOURCE_INVALID",
            "message": "Supply exactly one QC source ID",
        })
    source_type, source_id = selected_sources[0]
    try:
        result = EXPERIMENTS.create_qc_review(
            project_id, source_type=source_type, source_id=source_id,
            sample_id=payload.sample_id, status=payload.status, reason=payload.reason,
            reviewer=payload.reviewer, expected_revision=payload.expected_revision,
        )
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        raise HTTPException(422, {"code": "QC_REVIEW_INVALID", "message": str(exc)}) from exc
    return _experiment_response(result)


@app.post("/api/v4/experiment-projects/{project_id}/gate-analyses", status_code=201)
def apply_gate_analysis(project_id: str, payload: GateAnalysisCreate) -> dict[str, Any]:
    event_exclusion_ranges_by_artifact: dict[str, list[dict[str, int]]] = {}
    event_exclusion_document: dict[str, Any] | None = None
    try:
        gate_revision = next(
            (item for item in EXPERIMENTS.list_gate_revisions(project_id)
             if item["id"] == payload.gate_revision_id),
            None,
        )
        if gate_revision is None:
            raise KeyError(payload.gate_revision_id)
        external_matrix = None
        if gate_revision["compensation_decision"] == "EXTERNAL_MATRIX":
            if not gate_revision.get("compensation_artifact_id"):
                raise ValueError("External compensation matrix artifact is missing from the gate revision")
            matrix_path, matrix_artifact = EXPERIMENTS.artifact_path(
                project_id, gate_revision["compensation_artifact_id"]
            )
            if matrix_artifact["kind"] != "compensation_matrix":
                raise ValueError("Gate revision does not reference a compensation matrix artifact")
            external_matrix = read_json(matrix_path)
        EXPERIMENTS.validate_gate_analysis_inputs(
            project_id, payload.gate_revision_id, payload.fcs_artifact_ids
        )
        if payload.window_exclusions:
            if not payload.source_gate_analysis_id:
                raise ValueError("Manual time-window exclusions must cite the source gate analysis")
            if not payload.exclusion_reviewer or not payload.exclusion_reviewer.strip():
                raise ValueError("A reviewer is required for manual time-window exclusions")
            if not payload.exclusion_reason.strip():
                raise ValueError("A reason is required for manual time-window exclusions")
            source_analysis = next((item for item in EXPERIMENTS.list_gate_analyses(project_id)
                                    if item["id"] == payload.source_gate_analysis_id), None)
            if source_analysis is None or source_analysis["status"] != "COMPLETE":
                raise ValueError("The source gate analysis is missing or incomplete")
            if (source_analysis["gate_revision_id"] != payload.gate_revision_id
                    or source_analysis["input_artifact_ids"] != payload.fcs_artifact_ids):
                raise ValueError("Time-window exclusions must use the exact source gate revision and FCS artifact set")
            source_path, source_artifact = EXPERIMENTS.artifact_path(
                project_id, source_analysis["result_artifact_id"]
            )
            source_document = read_json(source_path)
            prior_exclusions = source_document.get("event_exclusions") or {}
            prior_window_rows = prior_exclusions.get("windows", [])
            prior_rows = {item["sample_id"]: item for item in prior_window_rows}
            if len(prior_rows) != len(prior_window_rows):
                raise ValueError("Source analysis contains duplicate sample exclusion records")
            source_samples = {
                str(item.get("sample", {}).get("sample_id")): item
                for item in source_document.get("samples", [])
            }
            registered_samples = {
                item["id"]: item for item in EXPERIMENTS.list_entities(project_id, "samples")
            }
            selected_sample_ids = [item.sample_id for item in payload.window_exclusions]
            if len(selected_sample_ids) != len(set(selected_sample_ids)):
                raise ValueError("A sample may appear only once in the time-window exclusion selection")
            exclusion_rows = {}
            for sample_id, prior_row in prior_rows.items():
                source_sample = source_samples.get(sample_id)
                registered_sample = registered_samples.get(sample_id)
                if source_sample is None or registered_sample is None:
                    raise ValueError("Source analysis contains an exclusion for an unknown sample")
                artifact_id = registered_sample["artifact_id"]
                if artifact_id not in payload.fcs_artifact_ids or prior_row.get("artifact_id") != artifact_id:
                    raise ValueError("Source analysis exclusion is not bound to the same FCS artifact set")
                original_count = source_sample.get("source_event_count", source_sample.get("event_count"))
                visible_count = source_sample.get("event_count")
                if (not isinstance(original_count, int) or isinstance(original_count, bool)
                        or original_count < 1 or not isinstance(visible_count, int)
                        or isinstance(visible_count, bool) or visible_count < 1):
                    raise ValueError("Source analysis has invalid original or remaining event counts")
                previous = sorted(prior_row.get("ranges", []),
                                  key=lambda item: (item.get("event_start", -1), item.get("event_end", -1)))
                removed = 0
                previous_end = 0
                for item in previous:
                    start, end = item.get("event_start"), item.get("event_end")
                    if (not isinstance(start, int) or isinstance(start, bool)
                            or not isinstance(end, int) or isinstance(end, bool)
                            or start < previous_end or end <= start or end > original_count):
                        raise ValueError("Source analysis contains invalid previous event exclusions")
                    removed += end - start
                    previous_end = end
                if (original_count - removed != visible_count
                        or prior_row.get("events_removed") != removed
                        or prior_row.get("events_after_exclusion") != visible_count):
                    raise ValueError("Source analysis exclusion counts do not match its remaining event count")
                event_exclusion_ranges_by_artifact[artifact_id] = previous
                exclusion_rows[sample_id] = {
                    **prior_row, "ranges": previous, "events_removed": removed,
                    "events_after_exclusion": visible_count,
                }
            for selection in payload.window_exclusions:
                source_sample = source_samples.get(selection.sample_id)
                registered_sample = registered_samples.get(selection.sample_id)
                if source_sample is None or registered_sample is None:
                    raise ValueError(f"Sample {selection.sample_id} is not part of the source analysis")
                artifact_id = registered_sample["artifact_id"]
                if artifact_id not in payload.fcs_artifact_ids:
                    raise ValueError("Excluded sample FCS artifact is not in the source analysis")
                available = {
                    (window.get("event_start"), window.get("event_end"))
                    for window in (source_sample.get("qc") or {}).get("time_windows", [])
                }
                ranges = []
                for item in selection.windows:
                    start, end = item.event_start, item.event_end
                    if end <= start:
                        raise ValueError("A time window must have event_end greater than event_start")
                    if (start, end) not in available:
                        raise ValueError("Each excluded range must match a complete time window in the source analysis")
                    ranges.append({"event_start": start, "event_end": end})
                ranges.sort(key=lambda item: (item["event_start"], item["event_end"]))
                if any(left["event_end"] > right["event_start"]
                       for left, right in zip(ranges, ranges[1:])):
                    raise ValueError("Selected time windows may not overlap")
                visible_event_count = source_sample.get("event_count")
                source_event_count = source_sample.get("source_event_count", visible_event_count)
                if not isinstance(source_event_count, int) or source_event_count < 1:
                    raise ValueError("Source analysis has an invalid original event count")
                if not isinstance(visible_event_count, int) or visible_event_count < 1:
                    raise ValueError("Source analysis has an invalid remaining event count")
                selected_count = sum(item["event_end"] - item["event_start"] for item in ranges)
                if selected_count >= visible_event_count:
                    raise ValueError("Manual exclusions must leave at least one event in every sample")
                previous = prior_rows.get(selection.sample_id, {}).get("ranges", [])
                previous = sorted(previous, key=lambda item: (item["event_start"], item["event_end"]))
                previous_removed = 0
                previous_end = 0
                for item in previous:
                    start, end = item.get("event_start"), item.get("event_end")
                    if (not isinstance(start, int) or not isinstance(end, int)
                            or start < previous_end or end <= start or end > source_event_count):
                        raise ValueError("Source analysis contains invalid previous event exclusions")
                    previous_removed += end - start
                    previous_end = end
                if source_event_count - previous_removed != visible_event_count:
                    raise ValueError("Source analysis exclusion counts do not match its remaining event count")
                original_ranges = []
                original_cursor = 0
                visible_cursor = 0
                for item in [*previous, {"event_start": source_event_count, "event_end": source_event_count}]:
                    segment_length = item["event_start"] - original_cursor
                    visible_segment_end = visible_cursor + segment_length
                    for selected in ranges:
                        start = max(selected["event_start"], visible_cursor)
                        end = min(selected["event_end"], visible_segment_end)
                        if start < end:
                            original_ranges.append({
                                "event_start": original_cursor + (start - visible_cursor),
                                "event_end": original_cursor + (end - visible_cursor),
                            })
                    original_cursor = item["event_end"]
                    visible_cursor = visible_segment_end
                cumulative_ranges = sorted(
                    [*previous, *original_ranges],
                    key=lambda item: (item["event_start"], item["event_end"]),
                )
                cumulative_removed = sum(item["event_end"] - item["event_start"]
                                         for item in cumulative_ranges)
                if cumulative_removed >= source_event_count:
                    raise ValueError("Manual exclusions must leave at least one source event")
                event_exclusion_ranges_by_artifact[artifact_id] = cumulative_ranges
                exclusion_rows[selection.sample_id] = {
                    "sample_id": selection.sample_id, "artifact_id": artifact_id,
                    "source_index_ranges": ranges, "newly_excluded_ranges": original_ranges,
                    "ranges": cumulative_ranges, "source_event_count": source_event_count,
                    "events_removed": cumulative_removed,
                    "events_after_exclusion": source_event_count - cumulative_removed,
                }
            event_exclusion_document = {
                "source_gate_analysis_id": payload.source_gate_analysis_id,
                "source_result_artifact_id": source_artifact["id"],
                "source_result_sha256": source_artifact["sha256"],
                "reviewer": payload.exclusion_reviewer.strip(),
                "reason": payload.exclusion_reason.strip(),
                "windows": [exclusion_rows[key] for key in sorted(exclusion_rows)],
            }
        elif (payload.source_gate_analysis_id or payload.exclusion_reviewer
              or payload.exclusion_reason.strip()):
            raise ValueError("Source analysis, reviewer, and reason are only valid with selected windows")
        fcs_paths = []
        for artifact_id in payload.fcs_artifact_ids:
            path, artifact = EXPERIMENTS.artifact_path(project_id, artifact_id)
            if artifact["kind"] != "raw_fcs":
                raise ValueError("Only registered raw FCS artifacts can be analyzed")
            fcs_paths.append(str(path))
    except (KeyError, IntegrityError) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        code = "SAMPLE_MAPPING_REQUIRED" if "SAMPLE_MAPPING_REQUIRED" in str(exc) else "GATE_ANALYSIS_INPUT_INVALID"
        raise HTTPException(422, {"code": code, "message": str(exc)}) from exc

    worker = ROOT / "integrations" / "workers" / "flowkit_worker.py"
    python = ROOT / ".tool_envs" / "flowkit" / "Scripts" / "python.exe"
    if not python.is_file() or not worker.is_file():
        raise HTTPException(503, {
            "code": "FLOWKIT_WORKER_UNAVAILABLE",
            "message": "FlowKit worker is unavailable",
        })
    request = {
        "schema_version": "1.0", "request_id": uuid.uuid4().hex,
        "operation": "apply_gate_revision", "paths": fcs_paths,
        "gate_revision": gate_revision,
        "compensation_matrix": external_matrix,
        "qc_policy": payload.qc_policy.model_dump(),
        "event_exclusion_ranges_by_path": [
            event_exclusion_ranges_by_artifact.get(artifact_id, [])
            for artifact_id in payload.fcs_artifact_ids
        ],
    }
    try:
        completed = subprocess.run(
            [str(python), str(worker)], input=json.dumps(request), capture_output=True,
            text=True, encoding="utf-8", timeout=180, check=False,
            env=os.environ | {"PYTHONIOENCODING": "utf-8"},
        )
        worker_response = json.loads(completed.stdout)
    except (subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
        raise HTTPException(422, {
            "code": "GATE_APPLICATION_FAILED", "message": str(exc),
        }) from exc
    if completed.returncode or worker_response.get("status") != "COMPLETE":
        errors = worker_response.get("errors") or [{"message": completed.stderr[-1000:]}]
        raise HTTPException(422, {
            "code": "GATE_APPLICATION_FAILED",
            "message": errors[0].get("message", "Gate application failed"),
        })
    worker_samples = worker_response.get("records", {}).get("samples", [])
    if event_exclusion_document is not None:
        source_rows = {item["artifact_id"]: item for item in event_exclusion_document["windows"]}
        if len(worker_samples) != len(payload.fcs_artifact_ids):
            raise HTTPException(422, {
                "code": "GATE_APPLICATION_FAILED", "message": "FlowKit result count does not match excluded FCS files",
            })
        for artifact_id, sample_result in zip(payload.fcs_artifact_ids, worker_samples, strict=True):
            expected = source_rows.get(artifact_id)
            if expected is None:
                continue
            if (sample_result.get("source_event_count") != expected["source_event_count"]
                    or sample_result.get("events_removed") != expected["events_removed"]
                    or sample_result.get("event_count") != expected["events_after_exclusion"]
                    or sample_result.get("event_exclusion_ranges") != expected["ranges"]):
                raise HTTPException(422, {
                    "code": "GATE_EXCLUSION_COUNT_MISMATCH",
                    "message": "FlowKit did not apply the reviewed event-window exclusions exactly",
                })
    try:
        result = EXPERIMENTS.save_gate_analysis(
            project_id,
            gate_revision_id=payload.gate_revision_id,
            input_artifact_ids=payload.fcs_artifact_ids,
            worker_result={"samples": worker_response["records"]["samples"],
                           "provenance": worker_response.get("provenance", {}),
                           "qc_policy": payload.qc_policy.model_dump(),
                           "event_exclusions": event_exclusion_document},
            expected_revision=payload.expected_revision,
        )
    except (KeyError, IntegrityError, RevisionConflict) as exc:
        raise _experiment_write_error(exc) from exc
    except ValueError as exc:
        raise HTTPException(422, {"code": "GATE_ANALYSIS_INVALID", "message": str(exc)}) from exc
    return _experiment_response(result)


@app.get("/", response_class=HTMLResponse)
def index() -> HTMLResponse:
    return HTMLResponse((WEB_ROOT / "index.html").read_text(encoding="utf-8"))


@app.get("/experiments", response_class=HTMLResponse)
def experiment_index() -> HTMLResponse:
    return HTMLResponse((WEB_ROOT / "experiments.html").read_text(encoding="utf-8"))


@app.get("/pd1-antibodies", response_class=HTMLResponse)
def pd1_antibodies_index() -> HTMLResponse:
    return HTMLResponse((WEB_ROOT / "pd1_antibodies.html").read_text(encoding="utf-8"))


app.mount("/assets", StaticFiles(directory=WEB_ROOT), name="assets")
