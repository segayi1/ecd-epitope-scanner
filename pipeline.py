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
    return all(value >= 0.9 for value in results.values()), results


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
    if mapping.observed_ecd_coverage < 0.95:
        candidate.rejection_reasons.append("observed ECD coverage < 95%")
    if mapping.max_internal_gap > 8:
        candidate.rejection_reasons.append("internal coordinate gap > 8 aa")
    if mapping.identity < 0.95:
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
        if mapping.observed_ecd_coverage >= 0.95 and mapping.max_internal_gap <= 8 and mapping.identity >= 0.95:
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
                mapping.observed_ecd_coverage >= 0.95
                and mapping.max_internal_gap <= 8
                and mapping.identity >= 0.95
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
    progress("PARTITIONING", 0.62, "Building domain- and membrane-constrained geodesic surface regions")
    regions, partition_summary = partition_surface(
        records,
        protein,
        confidence=confidence,
        minimum_regions=request.minimum_regions,
        maximum_regions=request.maximum_regions,
        predicted=predicted,
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
        "eligibility_gate": {
            "minimum_observed_ecd_coverage": 0.95,
            "maximum_internal_gap_aa": 8,
            "minimum_aligned_identity": 0.95,
            "minimum_major_domain_coverage": 0.90,
        },
        "candidate_count": len(candidates),
        "candidates_file": "structure_candidates.json",
    }
    warnings = list(partition_summary.get("warnings", []))
    if methods["sasa_method"].startswith("Bio.PDB"):
        warnings.append(
            "FreeSASA was unavailable; the auditable Bio.PDB Shrake-Rupley fallback was used with the same probe/point settings."
        )
    if not inventory.evoef2_path and not inventory.foldx_path:
        warnings.append(
            "Neither FoldX nor EvoEF2 is configured; combinations are ranked with the built-in dimensionless structural-risk proxy."
        )
    if predicted:
        warnings.append("Predicted structure confidence gates do not establish biological correctness.")

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
            "strict_af3_fallback": request.strict_af3,
            "minimum_regions": request.minimum_regions,
            "maximum_regions": request.maximum_regions,
            "regional_mutations": "adaptive 3-5; may reduce to 2 if no larger combination passes",
            "functional_interfaces": "included with persistent risk flag",
            "output_sequence": "ECD coding fragment only; no start/stop/signal peptide/tag/TM/cytoplasmic region",
            "random_seed": None,
        },
        "methods": methods,
        "stability": engine_metadata(inventory),
        "partitioning": partition_summary,
        "design_summary": design_summary,
        "summary": {
            "mapped_ecd_residues": len(records),
            "surface_accessible": sum(item.surface_accessible for item in records),
            "folding_screen_pass": sum(item.folding_screen_pass for item in records),
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
    analysis["source"]["structure_sha256"] = sha256_file(structure_cif)
    write_json(output_dir / "analysis.json", analysis)
    # Refresh manifest after the final analysis file is written.
    write_scientific_outputs(
        output_dir,
        analysis,
        protein,
        records,
        regions,
        designs,
        request.original_cds,
    )
    progress("COMPLETE", 1.0, "Analysis and regional combination design completed")
    return analysis
