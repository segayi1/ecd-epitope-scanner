from __future__ import annotations

import datetime as dt
import json
import platform
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import Bio
import networkx
import numpy as np
from Bio.PDB import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from . import config
from . import models
from .design import design_regions
from .exporting import write_scientific_outputs
from .models import ProteinContext, StructureCandidate, to_dict
from .partitioning import partition_surface
from .retrieval import (
    candidate_manifest,
    direct_pdb_candidate,
    fetch_protein_context,
    fetch_sifts_mapping,
    fetch_structure_candidates,
    find_public_af3_candidates,
    infer_uniprot_from_sifts,
    prefilter_experimental_candidates,
)
from .stability import engine_metadata, find_engine_inventory
from .sequences import coding_fragment
from .structure_analysis import (
    ChainMapping,
    ParsedStructure,
    analyze_residues,
    best_chain_mapping,
    load_structure,
    write_mapped_ecd_cif,
    write_mapped_ecd_pdb,
    write_chain_pdb,
)
from .util import PipelineError, RemoteClient, WaitingForAF3, read_json, safe_name, sha256_file, write_json


ProgressCallback = Callable[[str, float, str], None]


@dataclass(slots=True)
class PipelineRequest:
    input_type: str
    identifier: str = ""
    accession: str | None = None
    chain: str | None = None
    structure_path: str | None = None
    confidence_path: str | None = None
    sequence: str | None = None
    ecd_start: int | None = None
    ecd_end: int | None = None
    original_cds: str | None = None
    strict_af3: bool = True
    allow_partial: bool = False
    structure_candidate_id: str | None = None
    minimum_regions: int = 5
    maximum_regions: int = 8
    output_dir: str | None = None
    assay_format: str = "cell_surface"   # cell_surface | purified
    metadata: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        data = to_dict(self)
        if data.get("original_cds"):
            data["original_cds"] = "<provided>"
        return data


def _noop(stage: str, progress: float, message: str) -> None:
    del stage, progress, message


def _synthetic_context(request: PipelineRequest) -> ProteinContext:
    sequence = (request.sequence or "").strip().upper()
    if not sequence or request.ecd_start is None or request.ecd_end is None:
        raise PipelineError(
            "A structure without a UniProt accession requires sequence, ecd_start and ecd_end"
        )
    if request.ecd_start < 1 or request.ecd_end > len(sequence) or request.ecd_start > request.ecd_end:
        raise PipelineError("Invalid ECD range for the supplied sequence")
    return ProteinContext(
        accession=request.accession or "USER_STRUCTURE",
        protein_name="User-provided structure",
        organism="unspecified",
        sequence=sequence,
        ecd_start=request.ecd_start,
        ecd_end=request.ecd_end,
        membrane_anchor_position=request.ecd_end,
        topology="user-defined ECD boundary",
        source_url="user input",
    )


def _transmembrane_geometry(
    mapping: ChainMapping, protein: ProteinContext
) -> tuple[tuple[float, float, float] | None, tuple[float, float, float] | None]:
    """Derive membrane orientation only from an observed TM segment.

    A low-confidence ECD stalk is not a membrane landmark.  Missing TM
    coordinates therefore produce an explicit unknown orientation.
    """
    tm_residues = [
        residue for position, residue in mapping.uniprot_to_residue.items()
        if protein.membrane_anchor_position < position <= protein.membrane_anchor_position + 35
        and residue.has_id("CA")
    ]
    ecd_residues = [
        residue for position, residue in mapping.uniprot_to_residue.items()
        if protein.ecd_start <= position <= protein.ecd_end and residue.has_id("CA")
    ]
    if len(tm_residues) < 8 or len(ecd_residues) < 3:
        return None, None
    tm_coordinates = np.asarray([residue["CA"].coord for residue in tm_residues], dtype=float)
    tm_center = np.mean(tm_coordinates, axis=0)
    ecd_center = np.mean([residue["CA"].coord for residue in ecd_residues], axis=0)
    # The principal axis of a sufficiently long TM helix is the membrane
    # normal.  Orient its sign towards the ECD; do not use the ECD centroid as
    # the axis itself, because a tilted/asymmetric ectodomain would bias it.
    _, _, axes = np.linalg.svd(tm_coordinates - tm_center, full_matrices=False)
    normal = axes[0]
    if float(np.dot(normal, ecd_center - tm_center)) < 0:
        normal = -normal
    if float(np.linalg.norm(normal)) < 1e-6:
        return None, None
    return tuple(float(value) for value in tm_center), tuple(float(value) for value in normal)


def _domain_coverage(protein: ProteinContext, mapping: ChainMapping) -> tuple[bool, dict[str, float]]:
    observed = set(mapping.uniprot_to_residue)
    results: dict[str, float] = {}
    major_domains = [
        feature
        for feature in protein.features
        if feature.type.lower() in {"domain", "repeat"}
        and max(feature.start, protein.ecd_start) <= min(feature.end, protein.ecd_end)
        and min(feature.end, protein.ecd_end) - max(feature.start, protein.ecd_start) + 1 >= 15
    ]
    for feature in major_domains:
        start = max(feature.start, protein.ecd_start)
        end = min(feature.end, protein.ecd_end)
        coverage = sum(position in observed for position in range(start, end + 1)) / (end - start + 1)
        results[feature.description or f"{feature.type} {start}-{end}"] = round(coverage, 4)
    return all(value >= config.MIN_MAJOR_DOMAIN_COVERAGE for value in results.values()), results


def _extract_citation(path: Path) -> dict[str, Any]:
    if path.suffix.lower() not in {".cif", ".mmcif"}:
        return {}
    try:
        data = MMCIF2Dict(str(path))
    except Exception:
        return {}

    def first(key: str) -> str:
        value = data.get(key, "")
        if isinstance(value, list):
            return str(value[0]) if value else ""
        return str(value or "")

    title = first("_citation.title").strip("'\"")
    pmid = first("_citation.pdbx_database_id_PubMed")
    doi = first("_citation.pdbx_database_id_DOI")
    year = first("_citation.year")
    if not any((title, pmid, doi, year)):
        return {}
    return {"title": title, "pmid": pmid, "doi": doi, "year": year}


def _candidate_suffix(candidate: StructureCandidate) -> str:
    return ".pdb" if candidate.model_format.upper() == "PDB" or candidate.model_url.lower().endswith(".pdb") else ".cif"


def _evaluate_candidate(
    candidate: StructureCandidate,
    protein: ProteinContext,
    client: RemoteClient,
    candidate_dir: Path,
    preferred_chain: str | None = None,
) -> tuple[ParsedStructure, ChainMapping]:
    suffix = _candidate_suffix(candidate)
    local = candidate_dir / f"{safe_name(candidate.identifier)}{suffix}"
    if not local.exists():
        client.download(candidate.model_url, local)
    candidate.local_path = str(local.resolve())
    parsed = load_structure(local)
    preferred = [preferred_chain] if preferred_chain else candidate.chain_ids
    mapping = best_chain_mapping(parsed, protein, preferred)
    candidate.selected_chain = mapping.chain_id
    candidate.observed_ecd_coverage = round(mapping.observed_ecd_coverage, 4)
    candidate.max_internal_gap = mapping.max_internal_gap
    candidate.engineered_mismatches = mapping.engineered_mismatches
    domain_complete, domain_coverages = _domain_coverage(protein, mapping)
    candidate.metadata["aligned_sequence_identity"] = round(mapping.identity, 4)
    candidate.metadata["domain_coverages"] = domain_coverages
    candidate.metadata["citation"] = _extract_citation(local)
    candidate.rejection_reasons = []
    if mapping.observed_ecd_coverage < config.MIN_ECD_COVERAGE:
        candidate.rejection_reasons.append("observed ECD coverage < 95%")
    if mapping.max_internal_gap > config.MAX_INTERNAL_GAP_AA:
        candidate.rejection_reasons.append("internal coordinate gap > 8 aa")
    if mapping.identity < config.MIN_SEQUENCE_IDENTITY:
        candidate.rejection_reasons.append("aligned sequence identity < 95%")
    if not domain_complete:
        candidate.rejection_reasons.append("one or more major ECD domains have < 90% observed coverage")
    candidate.eligible = not candidate.rejection_reasons
    return parsed, mapping


def _select_from_uniprot(
    request: PipelineRequest,
    protein: ProteinContext,
    client: RemoteClient,
    output_dir: Path,
    progress: ProgressCallback,
) -> tuple[StructureCandidate, ParsedStructure, ChainMapping, list[StructureCandidate], bool]:
    progress("SELECTING_STRUCTURE", 0.12, "Querying 3D-Beacons for experimental and public model structures")
    candidates = fetch_structure_candidates(protein, client, output_dir)
    candidate_dir = output_dir / "candidates"
    candidate_dir.mkdir(parents=True, exist_ok=True)
    experimental = prefilter_experimental_candidates(protein, candidates, limit=10)
    if request.structure_candidate_id:
        experimental = [
            item for item in experimental if item.identifier.lower() == request.structure_candidate_id.lower()
        ]
        if not experimental:
            raise PipelineError(f"Requested structure candidate was not found: {request.structure_candidate_id}")

    evaluated: list[tuple[StructureCandidate, ParsedStructure, ChainMapping]] = []
    for index, candidate in enumerate(experimental):
        progress(
            "SELECTING_STRUCTURE",
            0.14 + 0.16 * (index / max(1, len(experimental))),
            f"Auditing experimental structure {candidate.identifier}",
        )
        try:
            parsed, mapping = _evaluate_candidate(
                candidate, protein, client, candidate_dir, request.chain
            )
            evaluated.append((candidate, parsed, mapping))
        except Exception as exc:
            candidate.eligible = False
            candidate.rejection_reasons.append(f"download/parse/mapping failed: {exc}")

    eligible = [item for item in evaluated if item[0].eligible]
    eligible.sort(
        key=lambda item: (
            -float(item[0].observed_ecd_coverage or 0),
            -float(item[0].metadata.get("aligned_sequence_identity", 0)),
            item[0].engineered_mismatches or 0,
            0 if item[0].metadata.get("citation") else 1,
            item[0].resolution if item[0].resolution is not None else 99.0,
            item[0].identifier,
        )
    )
    for rank, (candidate, _, _) in enumerate(eligible, start=1):
        candidate.rank = rank
    if eligible:
        selected, parsed, mapping = eligible[0]
        write_json(output_dir / "structure_candidates.json", candidate_manifest(candidates))
        return selected, parsed, mapping, candidates, False

    progress("SELECTING_STRUCTURE", 0.31, "No experimental structure passed the complete-ECD gate; checking explicit public AF3 models")
    public_af3 = find_public_af3_candidates(candidates, client)
    for candidate in public_af3:
        try:
            parsed, mapping = _evaluate_candidate(candidate, protein, client, candidate_dir, request.chain)
        except Exception as exc:
            candidate.rejection_reasons.append(f"AF3 model could not be audited: {exc}")
            continue
        if (
            mapping.observed_ecd_coverage >= config.MIN_ECD_COVERAGE
            and mapping.max_internal_gap <= config.MAX_INTERNAL_GAP_AA
            and mapping.identity >= config.MIN_SEQUENCE_IDENTITY
        ):
            candidate.eligible = True
            candidate.rank = 1
            write_json(output_dir / "structure_candidates.json", candidate_manifest(candidates))
            return candidate, parsed, mapping, candidates, True

    if not request.strict_af3:
        afdb = [
            item
            for item in candidates
            if item.provider.lower() == "alphafold db" and item.reported_coverage == 1.0
        ]
        if afdb:
            candidate = afdb[0]
            parsed, mapping = _evaluate_candidate(candidate, protein, client, candidate_dir, request.chain)
            candidate.eligible = True
            candidate.rank = 1
            candidate.metadata["scientific_warning"] = "AFDB fallback was explicitly enabled; this model is not asserted to be AF3."
            write_json(output_dir / "structure_candidates.json", candidate_manifest(candidates))
            return candidate, parsed, mapping, candidates, True

    write_json(output_dir / "structure_candidates.json", candidate_manifest(candidates))
    waiting = {
        "status": "WAITING_FOR_AF3_UPLOAD",
        "message": (
            "No experimental structure passed the >=95% ECD coverage/all-domain/max-gap gate, "
            "and no public model carried explicit AlphaFold 3 provenance. Upload the official AF3 model CIF and confidence JSON."
        ),
        "protein": {
            "accession": protein.accession,
            "ecd_start": protein.ecd_start,
            "ecd_end": protein.ecd_end,
        },
        "candidates": candidate_manifest(candidates),
    }
    write_json(output_dir / "waiting_for_af3.json", waiting)
    raise WaitingForAF3(waiting["message"], waiting["candidates"])


def _standardize_cif(parsed: ParsedStructure, output_path: Path) -> Path:
    if parsed.format == "mmcif":
        if parsed.path.resolve() != output_path.resolve():
            shutil.copy2(parsed.path, output_path)
    else:
        io = MMCIFIO()
        io.set_structure(parsed.structure)
        io.save(str(output_path))
    return output_path


def _software_versions() -> dict[str, str]:
    return {
        "ecd_epitope_workbench": models.__name__.split(".")[0] + " 1.0.0",
        "python": platform.python_version(),
        "biopython": Bio.__version__,
        "numpy": np.__version__,
        "networkx": networkx.__version__,
        "platform": platform.platform(),
    }


def run_pipeline(
    request: PipelineRequest,
    progress_callback: ProgressCallback | None = None,
) -> dict[str, Any]:
    progress = progress_callback or _noop
    output_dir = Path(request.output_dir or Path.cwd() / "results" / safe_name(request.identifier or "run"))
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "request.json", request.as_dict())
    client = RemoteClient()
    ecd_override = (
        (request.ecd_start, request.ecd_end)
        if request.ecd_start is not None and request.ecd_end is not None
        else None
    )
    progress("FETCHING", 0.03, "Resolving sequence, topology and structure provenance")

    parsed: ParsedStructure
    mapping: ChainMapping
    candidates: list[StructureCandidate] = []
    predicted = False
    source: dict[str, Any] = {}
    input_type = request.input_type.lower()
    if input_type == "uniprot":
        accession = request.identifier or request.accession or ""
        protein = fetch_protein_context(accession, client, output_dir, ecd_override)
        selected, parsed, mapping, candidates, predicted = _select_from_uniprot(
            request, protein, client, output_dir, progress
        )
        source = {
            "kind": "public database",
            "input_identifier": accession,
            "normalized_accession": protein.accession,
            "provider": selected.provider,
            "identifier": selected.identifier,
            "structure_url": selected.model_url,
            "page_url": selected.page_url,
            "category": selected.category,
            "explicit_af3": selected.explicit_af3,
            "downloaded_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
    elif input_type == "pdb_id":
        pdb_id = request.identifier
        sifts = fetch_sifts_mapping(pdb_id, client, output_dir)
        accession, sifts_chains = infer_uniprot_from_sifts(sifts, pdb_id)
        protein = fetch_protein_context(
            request.accession or accession, client, output_dir, ecd_override
        )
        selected = direct_pdb_candidate(pdb_id)
        selected.chain_ids = [request.chain] if request.chain else sifts_chains
        parsed, mapping = _evaluate_candidate(
            selected, protein, client, output_dir / "candidates", request.chain
        )
        candidates = [selected]
        if not selected.eligible and not request.allow_partial:
            raise PipelineError(
                f"PDB {pdb_id} failed the complete-ECD gate: {', '.join(selected.rejection_reasons)}"
            )
        source = {
            "kind": "PDB identifier",
            "provider": "PDBe",
            "identifier": pdb_id.lower(),
            "structure_url": selected.model_url,
            "page_url": selected.page_url,
            "downloaded_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
    elif input_type in {"structure", "af3_upload"}:
        structure_path = Path(request.structure_path or request.identifier)
        if request.accession:
            protein = fetch_protein_context(request.accession, client, output_dir, ecd_override)
        else:
            protein = _synthetic_context(request)
        parsed = load_structure(structure_path)
        mapping = best_chain_mapping(
            parsed, protein, [request.chain] if request.chain else None
        )
        predicted = input_type == "af3_upload" or bool(request.confidence_path)
        selected = StructureCandidate(
            identifier=structure_path.stem,
            category="AB-INITIO" if predicted else "USER-PROVIDED",
            provider="user upload",
            model_url=str(structure_path.resolve()),
            model_format=parsed.format.upper(),
            selected_chain=mapping.chain_id,
            observed_ecd_coverage=round(mapping.observed_ecd_coverage, 4),
            max_internal_gap=mapping.max_internal_gap,
            engineered_mismatches=mapping.engineered_mismatches,
            explicit_af3=input_type == "af3_upload",
            eligible=(
                mapping.observed_ecd_coverage >= config.MIN_ECD_COVERAGE
                and mapping.max_internal_gap <= config.MAX_INTERNAL_GAP_AA
                and mapping.identity >= config.MIN_SEQUENCE_IDENTITY
            ),
        )
        if not selected.eligible and not request.allow_partial:
            raise PipelineError(
                "Uploaded structure failed the complete-ECD mapping gate; enable allow_partial only for deliberate review runs"
            )
        candidates = [selected]
        source = {
            "kind": "AF3 upload" if predicted else "downloaded structure upload",
            "provider": "user",
            "identifier": structure_path.name,
            "original_path": str(structure_path.resolve()),
            "explicit_af3": input_type == "af3_upload",
        }
    else:
        raise PipelineError(f"Unsupported input_type: {request.input_type}")

    progress("ANALYZING", 0.42, "Computing residue-level SASA, interfaces and folding-risk features")
    structure_cif = _standardize_cif(parsed, output_dir / "structure.cif")
    write_mapped_ecd_cif(parsed, mapping, protein, output_dir / "ecd_structure.cif")
    # Folding ΔΔG is calculated for the mapped ECD heavy atoms alone. Other chains
    # remain in structure.cif for interface annotation and visual inspection.
    structure_pdb = write_mapped_ecd_pdb(parsed, mapping, protein, output_dir / "structure_for_energy.pdb")
    confidence = read_json(Path(request.confidence_path)) if request.confidence_path else None
    records, methods = analyze_residues(
        parsed,
        mapping,
        protein,
        output_dir,
        predicted=predicted,
        confidence=confidence,
    )
    mismatch_positions = list(getattr(mapping, "mismatch_positions", []) or [])
    mismatch_records = [
        record for record in records if record.uniprot_pos in set(mismatch_positions)
    ]
    for record in records:
        if any(
            abs(record.uniprot_pos - mismatch.uniprot_pos) <= 3
            or (
                record.ca is not None
                and mismatch.ca is not None
                and float(np.linalg.norm(np.asarray(record.ca) - np.asarray(mismatch.ca)))
                <= config.ENGINEERED_NEIGHBORHOOD_RADIUS_A
            )
            for mismatch in mismatch_records
        ):
            record.engineered_neighborhood = True
            if "engineered-neighborhood" not in record.risk_flags:
                record.risk_flags.append("engineered-neighborhood")
    methods["engineered_mismatch_positions"] = mismatch_positions
    methods["engineered_neighborhood_residues"] = sum(
        item.engineered_neighborhood for item in records
    )

    conservation_warnings: list[str] = []
    if input_type != "uniprot" and not request.accession:
        conservation_warnings.append(
            "No UniProt accession; per-residue conservation was not computed."
        )
    else:
        try:
            from .conservation import fetch_conservation

            progress("ANALYZING", 0.48, "Fetching reviewed orthologs for conservation scoring")
            scores = fetch_conservation(protein, client, output_dir)
            if scores is None:
                conservation_status = read_json(output_dir / "conservation.json") if (output_dir / "conservation.json").exists() else {}
                state = conservation_status.get("status", "unavailable")
                conservation_warnings.append(
                    f"Conservation was not computed ({state}); epitope priors use a neutral variability term."
                )
            else:
                matched = 0
                for record in records:
                    value = scores.get(record.uniprot_pos)
                    if value is not None:
                        record.conservation_score = value
                        matched += 1
                if not matched:
                    conservation_warnings.append(
                        "Ortholog alignment covered none of the mapped ECD positions; conservation was not applied."
                    )
        except Exception as exc:
            conservation_warnings.append(
                f"Conservation computation failed ({type(exc).__name__}); continuing without it."
            )

    progress("PARTITIONING", 0.62, "Building Fab-footprint surface patches with epitope priors")
    # Panel budget adapts to ECD size (config.adaptive_maximum_regions) and is
    # additionally capped by any explicit user request.
    adaptive_cap = config.adaptive_maximum_regions(protein.ecd_length)
    effective_max = min(request.maximum_regions, adaptive_cap)
    effective_min = min(request.minimum_regions, effective_max)
    membrane_anchor, membrane_normal = _transmembrane_geometry(mapping, protein)
    regions, partition_summary = partition_surface(
        records,
        protein,
        confidence=confidence,
        minimum_regions=effective_min,
        maximum_regions=effective_max,
        predicted=predicted,
        assay_format=request.assay_format,
        membrane_anchor=membrane_anchor,
        membrane_normal=membrane_normal,
    )
    inventory = find_engine_inventory()
    if request.metadata.get("disable_external_engines"):
        from .stability import EngineInventory

        inventory = EngineInventory(True, None, None)
    progress("DESIGNING", 0.76, "Evaluating in-silico saturation substitutions and regional combinations")
    designs, design_summary = design_regions(
        regions,
        records,
        protein,
        inventory,
        structure_pdb,
        output_dir,
        request.original_cds,
    )
    surface_count = sum(item.surface_accessible for item in records)
    if (
        protein.ecd_length <= config.SMALL_ECD_LENGTH_AA
        or surface_count <= config.SMALL_ECD_SURFACE_RESIDUES
    ):
        fine_windows: list[dict[str, Any]] = []
        for region in regions:
            positions = sorted(region.residues)
            for index in range(0, len(positions), config.SMALL_ECD_FINE_MAP_WINDOW_RESIDUES):
                window = positions[index : index + config.SMALL_ECD_FINE_MAP_WINDOW_RESIDUES]
                if window:
                    fine_windows.append(
                        {
                            "patch_id": region.id,
                            "positions": window,
                            "purpose": "Stage-2 overlapping fine-mapping window after a Stage-1 binding hit",
                        }
                    )
        design_summary.setdefault("epitope_mapping_plan", {})["resolution_mode"] = (
            "small-ECD fine mapping: Stage-2 single mutants and overlapping residue windows are prioritised"
        )
        design_summary["epitope_mapping_plan"]["fine_mapping_windows"] = fine_windows
        partition_summary["resolution_level"] = "enhanced_stage2"
    else:
        partition_summary["resolution_level"] = "standard"
    foldx_diag = design_summary.get("foldx_diagnostics") or {}
    if foldx_diag.get("all_failed"):
        # P1-3: a completely silent FoldX batch is an engine failure, not a
        # scientific conclusion.  Fail the job loudly instead of stamping every
        # region FOLDX_REVIEW_REQUIRED around proxy numbers.
        raise PipelineError(
            "ERROR_FOLDX_ENGINE: FoldX BuildModel returned no parsable scores for "
            f"{foldx_diag.get('batch_constructs', 0)} constructs "
            "(single-point screen included). See run_manifest.json foldx block."
        )

    selected_candidate = candidates[0] if input_type != "uniprot" else next(
        (item for item in candidates if item.identifier == source.get("identifier")), candidates[0]
    )
    selection = {
        "selected_identifier": source.get("identifier", selected_candidate.identifier),
        "selected_chain": mapping.chain_id,
        "observed_ecd_coverage": mapping.observed_ecd_coverage,
        "max_internal_gap": mapping.max_internal_gap,
        "aligned_sequence_identity": mapping.identity,
        "engineered_mismatches": mapping.engineered_mismatches,
        "engineered_mismatch_positions": mismatch_positions,
        "eligibility_gate": {
            "minimum_observed_ecd_coverage": config.MIN_ECD_COVERAGE,
            "maximum_internal_gap_aa": config.MAX_INTERNAL_GAP_AA,
            "minimum_aligned_identity": config.MIN_SEQUENCE_IDENTITY,
            "minimum_major_domain_coverage": config.MIN_MAJOR_DOMAIN_COVERAGE,
        },
        "candidate_count": len(candidates),
        "candidates_file": "structure_candidates.json",
    }
    warnings = list(partition_summary.get("warnings", [])) + list(conservation_warnings)
    warnings.append("Production SASA uses Bio.PDB Shrake-Rupley with pinned Tien/DSSP radii.")
    if not inventory.evoef2_path and not inventory.foldx_path:
        warnings.append(
            "Neither FoldX nor EvoEF2 is configured; combinations are ranked with the built-in dimensionless structural-risk proxy."
        )
    if predicted:
        warnings.append("Predicted structure confidence gates do not establish biological correctness.")
    if request.assay_format == "cell_surface":
        warnings.append(
            "cell_surface assay: membrane-proximal patches were down-weighted because IgG cannot reach them on intact cells."
        )

    wt_coding_sequence, wt_coding_source = coding_fragment(
        protein, protein.ecd_sequence, request.original_cds
    )

    analysis = {
        "schema_version": "1.0",
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "software": _software_versions(),
        "request": request.as_dict(),
            "source": source,
        "protein": {
            "accession": protein.accession,
            "protein_name": protein.protein_name,
            "organism": protein.organism,
            "reviewed": protein.reviewed,
            "sequence_length": len(protein.sequence),
            "ecd_start": protein.ecd_start,
            "ecd_end": protein.ecd_end,
            "ecd_length": protein.ecd_length,
            "ecd_sequence": protein.ecd_sequence,
            "ecd_coding_sequence": wt_coding_sequence,
            "ecd_coding_source": wt_coding_source,
            "topology": protein.topology,
            "membrane_anchor_position": protein.membrane_anchor_position,
            "features": [to_dict(item) for item in protein.features],
        },
        "structure_selection": selection,
        "parameters": {
            "scientific_thresholds": {
                "minimum_ecd_coverage": config.MIN_ECD_COVERAGE,
                "maximum_internal_gap_aa": config.MAX_INTERNAL_GAP_AA,
                "minimum_sequence_identity": config.MIN_SEQUENCE_IDENTITY,
                "minimum_major_domain_coverage": config.MIN_MAJOR_DOMAIN_COVERAGE,
                "rsasa_threshold": config.RSASA_THRESHOLD,
                "sidechain_sasa_threshold_A2": config.SIDECHAIN_SASA_THRESHOLD_A2,
                "minimum_orthologs": config.CONSERVATION_MIN_ORTHOLOGS,
                "engineered_neighborhood_radius_A": config.ENGINEERED_NEIGHBORHOOD_RADIUS_A,
                "membrane_accessibility_scale_A": config.MEMBRANE_ACCESSIBILITY_SCALE_A,
            },
            "strict_af3_fallback": request.strict_af3,
            "assay_format": request.assay_format,
            "minimum_regions": effective_min,
            "maximum_regions": effective_max,
            "adaptive_region_budget": {
                "requested_maximum": request.maximum_regions,
                "adaptive_cap_from_ecd_length": adaptive_cap,
                "policy": "~1 patch per 22 ECD residues, capped at 7 (config.ADAPTIVE_REGIONS_*)",
            },
            "patch_scale": "Fab footprint: geodesic radius 13 A, 15-25 residues (docs/改进方案.md)",
            "compact_single_domain_surface_face": {
                "enabled_only_after_legacy_quality_failure": True,
                "ecd_max_length_aa": config.COMPACT_ECD_MAX_LENGTH_AA,
                "dominant_domain_min_fraction": config.COMPACT_DOMINANT_DOMAIN_MIN_FRACTION,
                "max_surface_regions": config.COMPACT_MAX_SURFACE_REGIONS,
                "footprint_geodesic_radius_A": config.COMPACT_FOOTPRINT_GEODESIC_RADIUS_A,
                "footprint_euclidean_cap_A": config.COMPACT_FOOTPRINT_EUCLIDEAN_CAP_A,
                "combo_max_spread_A": config.COMPACT_COMBO_MAX_SPREAD_A,
            },
            "regional_mutations": "directed alphabet + anchor mutation; exhaustive enumeration of cleared menus",
            "functional_interfaces": "included with persistent risk flag",
            "output_sequence": "ECD coding fragment only; no start/stop/signal peptide/tag/TM/cytoplasmic region",
            "random_seed": None,
        },
        "methods": methods,
        "stability": engine_metadata(inventory),
        "partitioning": partition_summary,
        "design_summary": design_summary,
        "epitope_mapping_plan": design_summary.get("epitope_mapping_plan", {}),
        "summary": {
            "mapped_ecd_residues": len(records),
            "surface_accessible": sum(item.surface_accessible for item in records),
            "folding_screen_pass": sum(item.folding_screen_pass for item in records),
            "flexible_terminal_excluded": partition_summary.get(
                "excluded_flexible_terminal_residues", 0
            ),
            "engineered_neighborhood_residues": sum(
                item.engineered_neighborhood for item in records
            ),
            "functional_interface": sum(item.interface for item in records),
            "regions": len(regions),
            "recommended_constructs": sum(
                item.status == "RECOMMENDED_LOW_PREDICTED_RISK" for item in designs
            ),
        },
        "residues": [to_dict(item) for item in records],
        "regions": [to_dict(item) for item in regions],
        "designs": [to_dict(item) for item in designs],
        "warnings": warnings,
        "validation_gates": [
            "total and surface expression relative to WT",
            "non-competing conformational antibody or domain folding control",
            "monodispersity/thermal stability for purified ECD",
            "native ligand/function when interface mutations are interpreted",
        ],
        "files": {},
    }
    progress("EXPORTING", 0.91, "Writing auditable JSON, CSV, FASTA, GenBank and methods outputs")
    # structure.cif was written before residue analysis; hashing it here lets
    # the export run exactly once instead of twice (the old second pass existed
    # only to refresh this hash into the manifest).
    analysis["source"]["structure_sha256"] = sha256_file(structure_cif)
    files = write_scientific_outputs(
        output_dir,
        analysis,
        protein,
        records,
        regions,
        designs,
        request.original_cds,
    )
    analysis["files"] = files
    # Flexible terminal segments: excluded from patches, kept as linear-epitope
    # candidates with a ready-to-order alanine-scan block list.
    tail_ranges = methods.get("flexible_terminal_ranges") or []
    if tail_ranges:
        analysis["epitope_mapping_plan"]["linear_epitope_candidates"] = [
            {
                "range": f"{item['start']}-{item['end']}",
                "reason": (
                    "flexible terminal segment (no observed secondary structure): excluded from "
                    "patch probes because tail conformation adapts and negative results would be uninterpretable"
                ),
                "alanine_scan_blocks": [
                    f"{start}-{min(start + 2, item['end'])}"
                    for start in range(item["start"], item["end"] + 1, 3)
                ],
            }
            for item in tail_ranges
        ]
    progress("COMPLETE", 1.0, "Analysis and epitope-mapping panel design completed")
    return analysis
