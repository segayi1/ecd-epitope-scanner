from __future__ import annotations

import asyncio
import json
from pathlib import Path
from typing import Any, Literal

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .jobs import JobRunner, JobStore
from .models import DesignRecord, FeatureRange, MutationChoice, ProteinContext, RegionRecord, ResidueRecord, to_dict
from .pipeline import PipelineRequest
from .sequences import apply_mutations, coding_fragment
from .stability import evaluate_substitution
from .util import read_json, safe_name, write_json


ROOT = Path(__file__).resolve().parent
WEB_ROOT = ROOT / "web"
DATA_ROOT = ROOT / "workspace_data"
STORE = JobStore(DATA_ROOT)
RUNNER = JobRunner(STORE)


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


class MutationInput(BaseModel):
    position: int
    mutant: str = Field(pattern="^[ACDEFGHIKLMNPQRSTVWY]$")


class CustomDesignInput(BaseModel):
    name: str = Field(default="Custom construct", min_length=1, max_length=80)
    mutations: list[MutationInput] = Field(min_length=1, max_length=12)
    force_high_risk: bool = False


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


def _public_job(job: dict[str, Any]) -> dict[str, Any]:
    result = job.get("result")
    if result and job["status"] == "COMPLETE":
        result = dict(result)
        result["structure_url"] = f"/api/jobs/{job['id']}/files/ecd_structure.cif"
        result["download_url"] = f"/api/jobs/{job['id']}/files/ecd_epitope_design.zip"
        result["job_id"] = job["id"]
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
        key: job[key] for key in ("id", "status", "stage", "progress", "message", "error", "created_at", "updated_at")
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
) -> dict[str, Any]:
    suffix = Path(structure.filename or "structure.cif").suffix.lower()
    if suffix not in {".cif", ".mmcif", ".pdb", ".ent"}:
        raise HTTPException(422, "Structure must be PDB or mmCIF")
    placeholder = PipelineRequest(
        input_type="structure",
        identifier=structure.filename or "structure",
        accession=accession or None,
        chain=chain or None,
        sequence=sequence or None,
        ecd_start=ecd_start,
        ecd_end=ecd_end,
        original_cds=original_cds or None,
        allow_partial=allow_partial,
    )
    job = STORE.create(placeholder)
    target = Path(job["workdir"]) / "uploads" / f"structure{suffix}"
    await _save_upload(structure, target)
    request = STORE.request_for(job["id"])
    request.structure_path = str(target)
    STORE.update(
        job["id"], request_json=request.as_dict() | {"original_cds": request.original_cds}
    )
    RUNNER.submit(job["id"])
    return _public_job(STORE.get(job["id"]))


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str) -> dict[str, Any]:
    return _public_job(_job_or_404(job_id))


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


def _hydrate_for_custom(job: dict[str, Any]) -> tuple[
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
            errors.append(f"{record.aa}{item.position}: protected by {', '.join(record.hard_exclusions)}")
            continue
        choice = evaluate_substitution(record, item.mutant, protein)
        if not choice.allowed and not payload.force_high_risk:
            errors.append(f"{choice.notation}: {', '.join(choice.reasons)}")
            continue
        choices.append(choice)
    if errors:
        raise HTTPException(422, detail={"message": "Custom design failed validation", "issues": errors})
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
            ["High-risk substitutions were explicitly forced by the user."] if payload.force_high_risk else []
        ) + [f"Coding fragment source: {dna_source}."],
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


@app.get("/", response_class=HTMLResponse)
def index() -> HTMLResponse:
    return HTMLResponse((WEB_ROOT / "index.html").read_text(encoding="utf-8"))


app.mount("/assets", StaticFiles(directory=WEB_ROOT), name="assets")
