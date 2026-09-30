from __future__ import annotations

import datetime as dt
import platform
import shutil
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import Bio
import networkx
import numpy as np
from Bio.PDB import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from . import config, models
from .assay_profiles import get_assay_profile
from .design import design_regions
from .exporting import write_scientific_outputs
from .joint_design import run_joint_v3
from .joint_design_v4 import run_joint_v4_ml
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
from .sequences import coding_fragment
from .stability import engine_metadata, find_engine_inventory
from .structure_analysis import (
    ChainMapping,
    ParsedStructure,
    analyze_residues,
    best_chain_mapping,
    chain_sequence,
    infer_antigen_chain,
    load_structure,
    standard_residues,
    write_mapped_ecd_cif,
    write_mapped_ecd_pdb,
)
from .structure_reliability import classify_reliability, export_reliability, import_iupred
from .surface_patch import run_surface_patch_v2
from .util import (
    PipelineError,
    RemoteClient,
    WaitingForAF3,
    read_json,
    safe_name,
    sha256_file,
    write_json,
)

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
    assay_profile: str | None = None
    design_engine: str = "joint_v4_ml"
    ml_mode: str = "strict"              # strict | preview (joint_v4_ml only)
    max_sectors: int = 7
    max_primary_constructs: int = 7
    allow_spurs_only_final: bool = False
    biological_context_mode: str = "provided_structure"
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


def _pdb_only_context(
    parsed: ParsedStructure,
    request: PipelineRequest,
    assay_profile: dict[str, Any],
) -> tuple[ProteinContext, str]:
    """Build an ECD-only context from observed coordinates.

    No canonical numbering or unobserved residues are invented.  When no
    chain is specified, a complex-aware antigen-chain heuristic is used and
    its evidence and uncertainty are recorded in the output.
    """
    chains = {
        str(chain.id): chain_sequence(chain)
        for chain in parsed.model.get_chains()
        if len(standard_residues(chain)) >= 8
    }
    if not chains:
        raise PipelineError("Uploaded structure contains no protein chain with at least 8 residues")
    chain_selection = infer_antigen_chain(parsed, request.chain)
    chain_id = str(chain_selection["chain_id"])
    sequence = chains[chain_id]
    return ProteinContext(
        accession=f"PDB_ONLY_CHAIN_{chain_id}",
        protein_name="Coordinate-derived antigen ECD",
        organism="unspecified",
        sequence=sequence,
        ecd_start=1,
        ecd_end=len(sequence),
        membrane_anchor_position=len(sequence),
        topology="K562 display: N-terminal extracellular antigen ECD, C-terminal PDGFRB membrane anchor",
        source_url="user-provided coordinate file",
        raw={
            "sequence_source": "observed coordinate residues",
            "canonical_numbering_available": False,
            "auto_selected_chain": request.chain is None,
            "chain_selection": chain_selection,
            "available_chain_lengths": {key: len(value) for key, value in sorted(chains.items())},
            "assay_profile": assay_profile["id"],
        },
    ), chain_id


def _profile_orientation_geometry(
    mapping: ChainMapping,
    protein: ProteinContext,
) -> dict[str, Any]:
    """Construct a relative C-terminal anchor axis, never a fake membrane plane."""
    ordered = [
        (position, residue)
        for position, residue in sorted(mapping.uniprot_to_residue.items())
        if protein.ecd_start <= position <= protein.ecd_end and residue.has_id("CA")
    ]
    if len(ordered) < 4:
        return {
            "status": "unavailable",
            "axis_source": "C-terminal anchor surrogate",
            "reason": "fewer than four mapped ECD CA atoms",
        }
    tail = ordered[-3:]
    anchor = np.mean([residue["CA"].coord for _, residue in tail], axis=0)
    centroid = np.mean([residue["CA"].coord for _, residue in ordered], axis=0)
    axis = centroid - anchor
    norm = float(np.linalg.norm(axis))
    if norm < 1e-6:
        return {
            "status": "unavailable",
            "axis_source": "C-terminal anchor surrogate",
            "reason": "ECD centroid and C-terminal anchor coincide",
        }
    axis /= norm
    return {
        "status": "relative_axis_only",
        "axis_source": "C-terminal anchor surrogate",
        "anchor_residues": [position for position, _ in tail],
        "anchor_xyz": [round(float(value), 4) for value in anchor],
        "extracellular_axis": [round(float(value), 6) for value in axis],
        "physical_membrane_plane_inferred": False,
        "caveat": (
            "The axis supports relative membrane-proximal/distal labels only; it is not a "
            "membrane normal and does not estimate intercellular Fab reach."
        ),
    }


def _annotate_profile_orientation(
    regions: list[Any],
    records: list[Any],
    geometry: dict[str, Any],
) -> None:
    if geometry.get("status") != "relative_axis_only":
        for region in regions:
            region.membrane_position = "relative orientation unavailable"
            region.predicted_accessibility = "unknown"
            region.membrane_accessibility_score = None
            region.assay_accessible = True
        return
    anchor = np.asarray(geometry["anchor_xyz"], dtype=float)
    axis = np.asarray(geometry["extracellular_axis"], dtype=float)
    by_position = {item.uniprot_pos: item for item in records}
    heights: dict[int, float] = {}
    for region in regions:
        positions = region.envelope_residues or region.residues
        coords = [by_position[pos].ca for pos in positions if pos in by_position and by_position[pos].ca]
        if coords:
            heights[region.id] = float(np.dot(np.mean(np.asarray(coords), axis=0) - anchor, axis))
    if not heights:
        return
    lower, upper = np.quantile(list(heights.values()), [1 / 3, 2 / 3])
    for region in regions:
        value = heights.get(region.id)
        if value is None:
            label = "relative orientation unavailable"
        elif value <= lower:
            label = "C-anchor-proximal"
        elif value >= upper:
            label = "C-anchor-distal"
        else:
            label = "C-anchor-middle"
        region.membrane_position = label
        region.predicted_accessibility = "unknown"
        region.membrane_accessibility_score = None
        region.assay_accessible = True
        region.rationale += (
            f"; {label} on a C-terminal-anchor relative axis; no membrane-reach penalty applied"
        )


def _select_geometry_spread_panel(
    regions: list[Any], records: list[Any], maximum: int
) -> list[Any]:
    """Greedy maximin coverage without epitope-likelihood ranking."""
    if len(regions) <= maximum:
        return list(regions)
    by_position = {item.uniprot_pos: item for item in records}

    def center(region: Any) -> np.ndarray | None:
        positions = region.envelope_residues or region.residues
        coords = [by_position[pos].ca for pos in positions if pos in by_position and by_position[pos].ca]
        return np.mean(np.asarray(coords), axis=0) if coords else None

    centers = {region.id: center(region) for region in regions}
    pending = list(regions)
    selected = [max(pending, key=lambda item: (item.auto_design_allowed, item.surface_area, -item.id))]
    pending.remove(selected[0])
    while pending and len(selected) < maximum:
        def coverage_key(candidate: Any) -> tuple[Any, ...]:
            point = centers[candidate.id]
            distances = [
                float(np.linalg.norm(point - centers[item.id]))
                for item in selected
                if point is not None and centers[item.id] is not None
            ]
            return (
                candidate.auto_design_allowed,
                min(distances) if distances else -1.0,
                candidate.surface_area,
                -candidate.id,
            )
        chosen = max(pending, key=coverage_key)
        selected.append(chosen)
        pending.remove(chosen)
    return sorted(selected, key=lambda item: item.id)


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
    preferred = [preferred_chain] if preferred_chain else None
    if preferred_chain:
        mapping = best_chain_mapping(
            parsed, protein, preferred, allowed_chains=[preferred_chain]
        )
    else:
        try:
            mapping = best_chain_mapping(
                parsed,
                protein,
                allowed_chains=candidate.chain_ids or None,
            )
        except PipelineError:
            # Some providers report label-chain identifiers while the parsed model
            # exposes author-chain identifiers.  Keep the mapping-based fallback,
            # but make the identifier mismatch visible in candidate provenance.
            mapping = best_chain_mapping(parsed, protein)
            if candidate.chain_ids:
                candidate.metadata["chain_selection_warning"] = (
                    "SIFTS/provider chain IDs did not match coordinate author-chain IDs; "
                    "selected chain by sequence alignment instead."
                )
    _apply_admission_gate(candidate, protein, mapping)
    candidate.metadata["citation"] = _extract_citation(local)
    return parsed, mapping


def _apply_admission_gate(
    candidate: StructureCandidate,
    protein: ProteinContext,
    mapping: ChainMapping,
) -> None:
    domain_complete, domain_coverages = _domain_coverage(protein, mapping)
    candidate.rejection_reasons = []
    candidate.selected_chain = mapping.chain_id
    candidate.observed_ecd_coverage = round(mapping.observed_ecd_coverage, 4)
    candidate.max_internal_gap = mapping.max_internal_gap
    candidate.engineered_mismatches = mapping.engineered_mismatches
    candidate.metadata["aligned_sequence_identity"] = round(mapping.identity, 4)
    candidate.metadata["domain_coverages"] = domain_coverages
    candidate.metadata["engineered_mismatch_positions"] = list(mapping.mismatch_positions)
    if mapping.observed_ecd_coverage < config.MIN_ECD_COVERAGE:
        candidate.rejection_reasons.append("observed ECD coverage < 95%")
    if mapping.max_internal_gap > config.MAX_INTERNAL_GAP_AA:
        candidate.rejection_reasons.append("internal coordinate gap > 8 aa")
    if mapping.identity < config.MIN_SEQUENCE_IDENTITY:
        candidate.rejection_reasons.append("aligned sequence identity < 95%")
    if not domain_complete:
        candidate.rejection_reasons.append("one or more major ECD domains have < 90% observed coverage")
    candidate.eligible = not candidate.rejection_reasons


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
    try:
        assay_profile = get_assay_profile(request.assay_profile)
    except ValueError as exc:
        raise PipelineError(str(exc)) from exc
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
    chain_selection: dict[str, Any] = {}
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
            "downloaded_at": dt.datetime.now(dt.UTC).isoformat(),
        }
    elif input_type == "pdb_id":
        pdb_id = request.identifier
        sifts = fetch_sifts_mapping(pdb_id, client, output_dir)
        try:
            accession, sifts_chains = infer_uniprot_from_sifts(
                sifts, pdb_id, request.accession
            )
        except PipelineError:
            accession, sifts_chains = request.accession or "", []
        selected = direct_pdb_candidate(pdb_id)
        candidate_dir = output_dir / "candidates"
        candidate_dir.mkdir(parents=True, exist_ok=True)
        local = candidate_dir / f"{safe_name(selected.identifier)}{_candidate_suffix(selected)}"
        if not local.exists():
            client.download(selected.model_url, local)
        selected.local_path = str(local.resolve())
        parsed = load_structure(local)
        if not request.accession and sifts:
            provisional = infer_antigen_chain(parsed)
            if provisional["confidence"] in {"high", "medium"}:
                try:
                    accession, sifts_chains = infer_uniprot_from_sifts(
                        sifts, pdb_id, preferred_chain=provisional["chain_id"]
                    )
                except PipelineError:
                    pass
        if accession:
            protein = fetch_protein_context(accession, client, output_dir, ecd_override)
            selected.chain_ids = [request.chain] if request.chain else sifts_chains
            parsed, mapping = _evaluate_candidate(
                selected, protein, client, candidate_dir, request.chain
            )
            chain_selection = {
                "chain_id": mapping.chain_id,
                "method": (
                    "user specified" if request.chain else
                    "SIFTS UniProt-to-chain mapping" if sifts_chains else
                    "best sequence alignment across protein chains"
                ),
                "confidence": "explicit" if request.chain else "high" if sifts_chains else "medium",
                "sifts_candidate_chains": sifts_chains,
            }
        elif request.sequence:
            protein = _synthetic_context(request)
            selected.chain_ids = [request.chain] if request.chain else []
            mapping = best_chain_mapping(
                parsed,
                protein,
                [request.chain] if request.chain else None,
                allowed_chains=[request.chain] if request.chain else None,
            )
            _apply_admission_gate(selected, protein, mapping)
            chain_selection = {
                "chain_id": mapping.chain_id,
                "method": "user specified" if request.chain else "best sequence alignment across protein chains",
                "confidence": "explicit" if request.chain else "medium",
                "aligned_sequence_identity": round(mapping.identity, 4),
                "observed_ecd_coverage": round(mapping.observed_ecd_coverage, 4),
            }
        elif assay_profile:
            protein, preferred_chain = _pdb_only_context(parsed, request, assay_profile)
            selected.chain_ids = [preferred_chain]
            mapping = best_chain_mapping(
                parsed, protein, [preferred_chain], allowed_chains=[preferred_chain]
            )
            _apply_admission_gate(selected, protein, mapping)
            chain_selection = dict(protein.raw.get("chain_selection", {}))
            chain_selection["sifts_candidate_chains"] = []
        else:
            raise PipelineError(
                f"PDB {pdb_id} has no UniProt/SIFTS target mapping; provide an accession, "
                "sequence with ECD boundaries, or an assay profile for coordinate-only analysis"
            )
        selected.metadata["chain_selection"] = chain_selection
        if protein.raw.get("canonical_numbering_available") is False:
            selected.metadata["coordinate_only_context"] = True
        candidates = [selected]
        source = {
            "kind": "PDB identifier",
            "provider": "PDBe",
            "identifier": pdb_id.lower(),
            "structure_url": selected.model_url,
            "page_url": selected.page_url,
            "downloaded_at": dt.datetime.now(dt.UTC).isoformat(),
            "chain_selection": chain_selection,
        }
    elif input_type in {"structure", "af3_upload"}:
        structure_path = Path(request.structure_path or request.identifier)
        parsed = load_structure(structure_path)
        preferred_chain = request.chain
        if request.accession:
            protein = fetch_protein_context(request.accession, client, output_dir, ecd_override)
        elif request.sequence:
            protein = _synthetic_context(request)
        elif assay_profile:
            protein, preferred_chain = _pdb_only_context(parsed, request, assay_profile)
        else:
            protein = _synthetic_context(request)
        mapping = best_chain_mapping(
            parsed,
            protein,
            [preferred_chain] if preferred_chain else None,
            allowed_chains=[preferred_chain] if preferred_chain else None,
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
        )
        _apply_admission_gate(selected, protein, mapping)
        chain_selection = dict(protein.raw.get("chain_selection", {}))
        if not chain_selection:
            chain_selection = {
                "chain_id": mapping.chain_id,
                "method": "user specified" if request.chain else "best sequence alignment across protein chains",
                "confidence": "explicit" if request.chain else "medium",
                "aligned_sequence_identity": round(mapping.identity, 4),
                "observed_ecd_coverage": round(mapping.observed_ecd_coverage, 4),
            }
        selected.metadata["chain_selection"] = chain_selection
        candidates = [selected]
        source = {
            "kind": "AF3 upload" if predicted else "downloaded structure upload",
            "provider": "user",
            "identifier": structure_path.name,
            "original_path": str(structure_path.resolve()),
            "explicit_af3": input_type == "af3_upload",
            "sequence_source": protein.raw.get("sequence_source", "provided or canonical sequence"),
            "chain_selection": chain_selection,
        }
    else:
        raise PipelineError(f"Unsupported input_type: {request.input_type}")

    if candidates:
        write_json(output_dir / "structure_candidates.json", candidate_manifest(candidates))

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
        ignored_positions=set(mapping.mismatch_positions),
    )
    mismatch_positions = list(getattr(mapping, "mismatch_positions", []) or [])
    mismatch_coordinates = [
        (position, residue["CA"].coord)
        for position in mismatch_positions
        if (residue := mapping.uniprot_to_residue.get(position)) is not None
        and residue.has_id("CA")
    ]
    for record in records:
        if any(
            abs(record.uniprot_pos - mismatch_position) <= 3
            or (
                record.ca is not None
                and float(np.linalg.norm(np.asarray(record.ca) - np.asarray(mismatch_ca)))
                <= config.ENGINEERED_NEIGHBORHOOD_RADIUS_A
            )
            for mismatch_position, mismatch_ca in mismatch_coordinates
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
    # Region count is a Stage-1 experimental budget, not a geometry limit.
    # Surface V2 enumerates the complete surface first and the design stage
    # chooses a bounded panel afterwards.
    if assay_profile:
        policy = assay_profile["design_policy"]
        effective_max = int(policy["stage1_maximum_constructs"])
        stage1_target = min(int(policy["stage1_target_constructs"]), effective_max)
        # A construct count is an experimental ceiling/target, not permission
        # to split one coherent surface into artificial extra regions.
        effective_min = 1
    else:
        effective_max = min(7, max(1, request.maximum_regions))
        effective_min = min(max(1, request.minimum_regions), effective_max)
        stage1_target = effective_min
    membrane_anchor, membrane_normal = _transmembrane_geometry(mapping, protein)
    assay_geometry = _profile_orientation_geometry(mapping, protein) if assay_profile else None
    patch_engine = str(
        request.metadata.get("patch_engine", request.design_engine or config.PATCH_ENGINE_VERSION)
    ).lower()
    if patch_engine != "joint_v4_ml":
        raise PipelineError(
            f"Unsupported design engine: {patch_engine}; new analyses use joint_v4_ml"
        )
    surface_v2_payload: dict[str, Any] | None = None
    joint_v3_payload: dict[str, Any] | None = None
    joint_v4_payload: dict[str, Any] | None = None
    reliability = None
    if patch_engine in {"surface_v2", "joint_v3", "joint_v4_ml"}:
        sequence_scores, sequence_provenance = import_iupred(
            request.metadata.get('iupred_output'), protein.sequence
        )
        reliability = classify_reliability(
            records, protein, predicted=predicted, confidence=confidence,
            chain_length=len(mapping.residues), single_chain=len(list(parsed.model.get_chains())) == 1,
            parameters=request.metadata.get('structure_reliability'),
            sequence_scores=sequence_scores, sequence_provenance=sequence_provenance,
        )
        export_reliability(reliability, output_dir)
        ecd_residues = {
            position: residue
            for position, residue in mapping.uniprot_to_residue.items()
        if protein.ecd_start <= position <= protein.ecd_end
        and position not in set(mismatch_positions)
        }
    inventory = find_engine_inventory()
    if request.metadata.get("disable_external_engines"):
        from .stability import EngineInventory

        inventory = EngineInventory(True, None, None)
    if patch_engine == "surface_v2":
        surface_result = run_surface_patch_v2(
            parsed.model,
            mapping.chain_id,
            ecd_residues,
            records,
            protein,
            output_dir=output_dir,
            confidence=confidence,
            predicted=predicted,
            parameters=config.SURFACE_V2_PARAMETERS,
            reliability=reliability,
        )
        # Natural communities remain visible regardless of the panel budget.
        # The legacy design adapter selects its input later in this function.
        regions = surface_result.legacy_regions
        records_by_position = {item.uniprot_pos: item for item in records}
        for record in records:
            record.region_id = None
        for region in regions:
            for position in region.envelope_residues:
                if position in records_by_position and records_by_position[position].region_id is None:
                    records_by_position[position].region_id = region.id
        partition_summary = surface_result.summary
        partition_summary["surface_region_count"] = len(surface_result.surface_regions)
        partition_summary["actual_regions"] = len(regions)
        footprint_payload = []
        region_footprint_ids = {fid for region in surface_result.surface_regions for fid in region.footprint_ids}
        for footprint in surface_result.raw_footprints:
            if footprint.footprint_id not in region_footprint_ids:
                continue
            item = to_dict(footprint)
            item.pop("vertex_ids", None)
            footprint_payload.append(item)
        surface_v2_payload = {
            "footprints": footprint_payload,
            "regions": [to_dict(item) for item in surface_result.surface_regions],
            "adjacency": [to_dict(item) for item in surface_result.footprint_edges],
            "qc": {
                key: value
                for key, value in partition_summary.items()
                if key in {
                    "surface_vertices", "surface_edges", "surface_components",
                    "center_count", "raw_footprint_count", "deduplicated_footprint_count",
                    "surface_region_count", "coverage",
                }
            },
        }
    elif patch_engine == "joint_v3":
        progress("PARTITIONING", 0.62, "Building Joint V3 surface units, multi-scale windows and probe sectors")
        v3_surface = run_surface_patch_v2(
            parsed.model,
            mapping.chain_id,
            ecd_residues,
            records,
            protein,
            output_dir=None,
            confidence=confidence,
            predicted=predicted,
            parameters=config.SURFACE_V2_PARAMETERS,
            reliability=reliability,
        ).surface
        joint_result = run_joint_v3(
            records,
            protein,
            output_dir,
            original_cds=request.original_cds,
            chain_id=mapping.chain_id,
            tm_anchor_xyz=membrane_anchor,
            tm_axis_xyz=membrane_normal,
            parameters=config.JOINT_V3_PARAMETERS,
            surface_geometry=v3_surface,
            structure=parsed.structure,
            source_provenance={
                "identifier": source.get("identifier"),
                "provider": source.get("provider"),
                "structure_sha256": sha256_file(parsed.path),
                "selected_chain": mapping.chain_id,
            },
            structure_pdb=structure_pdb,
            engine_inventory=inventory,
        )
        regions = joint_result.legacy_regions
        designs = joint_result.designs
        partition_summary = joint_result.partition_summary
        design_summary = joint_result.design_summary
        joint_v3_payload = {
            "schema_version": "3.0",
            "structure_audit": to_dict(joint_result.audit),
            "topology_context": to_dict(joint_result.topology),
            "microfootprints": to_dict(joint_result.microfootprints),
            "probe_sectors": to_dict(joint_result.sectors),
            "candidates": to_dict(joint_result.candidates),
            "selected_panel": to_dict(joint_result.panel),
        }
    elif patch_engine == "joint_v4_ml":
        progress("PROCALIPER_ANALYSIS", 0.62, "Running Joint V4-ML evidence adapters and surface audit")
        v4_surface_result = run_surface_patch_v2(
            parsed.model, mapping.chain_id, ecd_residues, records, protein,
            output_dir=None, confidence=confidence, predicted=predicted,
            parameters=config.SURFACE_V2_PARAMETERS, reliability=reliability,
        )
        v4_params = dict(request.metadata.get("joint_v4_ml", {}))
        v4_params.update({
            "mode": request.ml_mode,
            "max_sectors": request.max_sectors,
            "max_primary_constructs": min(request.max_primary_constructs, effective_max),
            "allow_spurs_only_final": request.allow_spurs_only_final,
            "use_surface_hypotheses": (
                assay_profile is None or int(assay_profile.get("version", 1)) >= 2
            ),
        })
        joint_result = run_joint_v4_ml(
            records, protein, output_dir, original_cds=request.original_cds,
            chain_id=mapping.chain_id, structure_pdb=structure_pdb,
            structure_type="alphafold" if predicted else "solved",
            parameters=v4_params, surface_geometry=v4_surface_result.surface,
            surface_v2_result=v4_surface_result,
            membrane_anchor=membrane_anchor,
            membrane_normal=membrane_normal,
            assay_orientation_geometry=assay_geometry,
            confidence=confidence,
            assay_format=request.assay_format,
            source_provenance={
                "identifier": source.get("identifier"), "provider": source.get("provider"),
                "structure_sha256": sha256_file(parsed.path), "selected_chain": mapping.chain_id,
                "biological_context_mode": request.biological_context_mode,
            },
            disable_external_engines=bool(request.metadata.get("disable_external_engines")),
        )
        regions = joint_result.legacy_regions
        designs = joint_result.legacy_designs
        partition_summary = {
            "engine": "joint_v4_ml", "schema_version": "4.0",
            "schema_extension_version": "4.1",
            "method": "Joint V4-ML 4.1: Surface V2 overlapping footprint hypotheses with stability-gated probe selection",
            "sector_count": len(joint_result.sectors),
            "primary_construct_count": len(joint_result.panel.selected_candidate_ids),
            "recommendation_status": joint_result.recommendation_status,
            "warnings": list(joint_result.warnings),
            "selected_panel": to_dict(joint_result.panel),
            "surface_evidence_count": len(joint_result.surface_evidence),
            "candidate_count": len(joint_result.candidates),
        }
        design_summary = {
            "engine": "joint_v4_ml", "schema_version": "4.0",
            "schema_extension_version": "4.1",
            "recommendation_status": joint_result.recommendation_status,
            "tool_inventory": to_dict(joint_result.tool_inventory),
            "candidate_count": len(joint_result.candidates),
            "selected_candidate_ids": joint_result.panel.selected_candidate_ids,
            "recommended_by_class": {
                key: len(value)
                for key, value in joint_result.panel.candidates_by_role.items()
                if key != "PRIMARY"
            },
            "higher_order_candidate_count": sum(item.order > 2 for item in joint_result.candidates),
            "epitope_mapping_plan": {
                "sectors": to_dict(joint_result.sectors),
                "overlapping_hypothesis_count": len(joint_result.surface_hypotheses),
                "hypothesis_source": "surface_v2_footprints",
                "hypotheses": to_dict(joint_result.surface_hypotheses),
                "panel_identifiability": {
                    "uniform_hypothesis_response": joint_result.panel.mean_hypothesis_response,
                    "prior_weighted_hypothesis_response": joint_result.panel.prior_weighted_hypothesis_response,
                    "marginal_information_gain": joint_result.panel.marginal_information_gain,
                    "mean_pairwise_separation": joint_result.panel.mean_pairwise_separation,
                    "minimum_signature_distance": joint_result.panel.minimum_signature_distance,
                    "prior_weighted_mean_pairwise_separation": joint_result.panel.prior_weighted_mean_pairwise_separation,
                    "unresolved_hypothesis_groups": joint_result.panel.unresolved_hypothesis_groups,
                    "uncovered_hypotheses": joint_result.panel.uncovered_hypotheses,
                },
                "soft_target": True,
            },
        }
        stage2_presets = [
            item for item in joint_result.legacy_designs
            if item.design_role in {"STAGE2_SINGLE", "STAGE2_PAIRWISE"}
        ]
        design_summary["stage2_preset_count"] = len(stage2_presets)
        stage1_designs = [
            item for item in joint_result.legacy_designs
            if item.design_role == "STAGE1_PROBE"
        ]
        unscored_saturation_count = sum(
            len(options)
            for design in stage1_designs
            for options in design.saturation_scan.values()
        )
        design_summary["saturation_scan_status"] = "ENUMERATED_NOT_SCORED"
        design_summary["unscored_saturation_substitution_count"] = unscored_saturation_count
        selected_by_id = {item.candidate_id: item for item in joint_result.candidates}
        stage1_probes = []
        coverage_matrix = []
        for index, design in enumerate(stage1_designs):
            candidate_id = joint_result.panel.selected_candidate_ids[index]
            candidate = selected_by_id[candidate_id]
            mutation_notations = [choice.notation for choice in design.mutations]
            stage1_probes.append({
                "probe_id": design.id,
                "patch_id": design.region_id,
                "candidate_id": candidate_id,
                "mutations": mutation_notations,
                "targeted_hypotheses": candidate.targeted_hypotheses,
                "hypothesis_response_proxy": candidate.hypothesis_response,
                "status": design.status,
            })
            scenario_coverage = candidate.coverage_by_scenario or {}
            for scenario, coverage in sorted(scenario_coverage.items()):
                coverage_matrix.append({
                    "probe_id": design.id,
                    "patch_id": design.region_id,
                    "scenario": scenario,
                    "predicted_coverage": coverage,
                    "evidence_type": "GEOMETRIC_MODEL_ONLY",
                })
            if not scenario_coverage:
                coverage_matrix.append({
                    "probe_id": design.id,
                    "patch_id": design.region_id,
                    "scenario": None,
                    "predicted_coverage": None,
                    "evidence_type": "NO_GEOMETRIC_COVERAGE_AVAILABLE",
                })
        design_summary["epitope_mapping_plan"].update({
            "stage1_probes": stage1_probes,
            "coverage_matrix": coverage_matrix,
            "blind_spots": to_dict(joint_result.panel.uncovered_components),
            "assay_note": (
                "Stage 1 is a predicted probe panel. Geodesic neighborhood coverage is a model output, "
                "not a measured epitope or binding result."
            ),
        })
        design_summary["epitope_mapping_plan"]["stage2_presets"] = [{
            "design_id": item.id,
            "parent_design_id": item.parent_design_id,
            "role": item.design_role,
            "mutations": [choice.notation for choice in item.mutations],
            "conditional_on_parent_probe_hit": True,
        } for item in stage2_presets]
        joint_v4_payload = {
            "schema_version": "4.0", "recommendation_status": joint_result.recommendation_status,
            "schema_extension_version": "4.1",
            "tool_inventory": to_dict(joint_result.tool_inventory),
            "surface_evidence": to_dict(joint_result.surface_evidence),
            "epitope_evidence": to_dict(joint_result.epitope_evidence),
            "sectors": to_dict(joint_result.sectors), "candidates": to_dict(joint_result.candidates),
            "selected_panel": to_dict(joint_result.panel), "warnings": joint_result.warnings,
            "surface_hypotheses": to_dict(joint_result.surface_hypotheses),
        }
    else:
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
    if assay_profile:
        _annotate_profile_orientation(regions, records, assay_geometry or {})
    progress("DESIGNING", 0.76, "Selecting a deterministic probe panel from the configured design engine")
    if patch_engine == "joint_v3":
        # Joint V3 has already enforced the panel budget in its optimizer.
        stage1_regions = regions
        partition_summary["stage1_selected_region_ids"] = [item.id for item in regions]
        partition_summary["stage1_min_constructs"] = min(effective_min, len(designs))
        partition_summary["stage1_target_constructs"] = min(stage1_target, effective_max)
        partition_summary["stage1_max_constructs"] = effective_max
        partition_summary["stage1_selection_basis"] = "joint_v3_robust_surface_coverage"
        partition_summary["epitope_prior_used_for_panel_selection"] = False
    elif patch_engine == "joint_v4_ml":
        stage1_regions = regions
        partition_summary["stage1_selected_region_ids"] = [item.id for item in regions]
        partition_summary["stage1_min_constructs"] = min(len(regions), len(designs))
        partition_summary["stage1_target_constructs"] = len(designs)
        partition_summary["stage1_max_constructs"] = min(
            request.max_primary_constructs, effective_max
        )
        partition_summary["stage1_selection_basis"] = (
            joint_result.panel.solver_status.lower()
        )
        partition_summary["epitope_prior_used_for_candidate_pool"] = bool(
            any(item.recommendation_class == "EPITOPE_PRIORITIZED" for item in joint_result.candidates)
        )
        partition_summary["epitope_prior_used_for_panel_selection"] = False
    else:
        stage1_regions = regions
        if patch_engine == "surface_v2" and len(regions) > effective_max:
            if assay_profile:
                stage1_regions = _select_geometry_spread_panel(regions, records, effective_max)
            else:
                pending = list(regions)
                stage1_regions = []
                covered_domains: set[str] = set()
                while pending and len(stage1_regions) < effective_max:
                    selected = max(pending, key=lambda r: (
                        r.auto_design_allowed, r.domain not in covered_domains,
                        len(r.envelope_residues), -r.id,
                    ))
                    stage1_regions.append(selected)
                    covered_domains.add(selected.domain)
                    pending.remove(selected)
                stage1_regions.sort(key=lambda r: r.id)
        if patch_engine == "surface_v2":
            selected_ids = {r.id for r in stage1_regions}
            partition_summary["coverage_plan"] = [r for r in partition_summary["coverage_plan"] if r["patch_id"] in selected_ids]
            partition_summary["unselected_region_ids"] = [r.id for r in regions if r.id not in selected_ids]
        partition_summary["stage1_selected_region_ids"] = [item.id for item in stage1_regions]
        partition_summary["stage1_min_constructs"] = min(effective_min, len(stage1_regions))
        partition_summary["stage1_target_constructs"] = stage1_target
        partition_summary["stage1_max_constructs"] = effective_max
        if assay_profile:
            partition_summary["stage1_selection_basis"] = assay_profile["design_policy"]["selection_basis"]
            partition_summary["epitope_prior_used_for_panel_selection"] = False
        designs, design_summary = design_regions(
            stage1_regions,
            records,
            protein,
            inventory,
            structure_pdb,
            output_dir,
            request.original_cds,
            minimum_combo_sites=(
                int(assay_profile["design_policy"]["combination_minimum_sites"])
                if assay_profile else None
            ),
            maximum_combo_sites=(
                int(assay_profile["design_policy"]["combination_maximum_sites"])
                if assay_profile else None
            ),
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
        "chain_selection": chain_selection,
        "observed_ecd_coverage": mapping.observed_ecd_coverage,
        "max_internal_gap": mapping.max_internal_gap,
        "aligned_sequence_identity": mapping.identity,
        "engineered_mismatches": mapping.engineered_mismatches,
        "engineered_mismatch_positions": mismatch_positions,
        "ignored_engineered_mismatch_positions": mismatch_positions,
        "eligibility_gate": {
            "minimum_observed_ecd_coverage": config.MIN_ECD_COVERAGE,
            "maximum_internal_gap_aa": config.MAX_INTERNAL_GAP_AA,
            "minimum_aligned_identity": config.MIN_SEQUENCE_IDENTITY,
            "minimum_major_domain_coverage": config.MIN_MAJOR_DOMAIN_COVERAGE,
        },
        "candidate_count": len(candidates),
        "candidates_file": "structure_candidates.json",
        "admission": {
            "eligible": selected_candidate.eligible,
            "status": "PASS" if selected_candidate.eligible else "CONTINUED_WITH_WARNINGS",
            "warnings": list(selected_candidate.rejection_reasons),
        },
    }
    warnings = list(partition_summary.get("warnings", [])) + list(conservation_warnings)
    for reason in selected_candidate.rejection_reasons:
        warnings.append(
            f"PDB admission warning: {reason}. Analysis continued; review mapping coverage and chain provenance."
        )
    warnings.extend(chain_selection.get("warnings", []))
    if selected_candidate.metadata.get("chain_selection_warning"):
        warnings.append(str(selected_candidate.metadata["chain_selection_warning"]))
    if mismatch_positions:
        positions = ", ".join(str(position) for position in mismatch_positions)
        warnings.append(
            f"The selected PDB chain differs from the UniProt standard sequence at "
            f"{len(mismatch_positions)} position(s): {positions}. These may be engineered "
            "substitutions or sequence variants. The original structure is retained for context; "
            "those positions are omitted from "
            "per-residue scoring, surface partitioning, and mutation design."
        )
    if protein.raw.get("canonical_numbering_available") is False:
        warnings.append(
            "No canonical UniProt sequence or numbering was available; this analysis uses only observed "
            "coordinates from the selected chain."
        )
    warnings.append("Production SASA uses Bio.PDB Shrake-Rupley with pinned Tien/DSSP radii.")
    if patch_engine == "joint_v3":
        warnings.extend([
            "Joint V3 sectors are coarse localization areas; they are not exact antibody footprints.",
            "Joint V3 does not output binding probabilities or kcal/mol values without a combination-level external engine.",
        ])
    elif patch_engine == "joint_v4_ml":
        warnings.extend(joint_result.warnings)
    if not inventory.evoef2_path and not inventory.foldx_path:
        warnings.append(
            "Neither FoldX nor EvoEF2 is configured; combinations are ranked with the built-in dimensionless structural-risk proxy."
        )
    if predicted:
        warnings.append("Predicted structure confidence gates do not establish biological correctness.")
    if assay_profile:
        warnings.append(
            "Fixed dual-cell assay profile: membrane-proximal regions were not down-weighted. "
            "The ECD-only structure supplies a relative C-terminal anchor axis, not a physical "
            "membrane plane or a quantitative Fab-reach model."
        )
        warnings.append(
            "N-terminal Myc accessibility controls K562 antigen surface display but does not "
            "prove native antigen folding; use an independent conformation-sensitive control."
        )
    elif request.assay_format == "cell_surface" and patch_engine != "joint_v4_ml":
        warnings.append(
            "cell_surface assay: membrane-proximal patches were down-weighted because IgG cannot reach them on intact cells."
        )

    wt_coding_sequence, wt_coding_source = coding_fragment(
        protein, protein.ecd_sequence, request.original_cds
    )

    analysis = {
        "schema_version": "1.0",
        "created_at": dt.datetime.now(dt.UTC).isoformat(),
        "software": _software_versions(),
        "request": request.as_dict(),
        "source": source,
        "assay_profile": assay_profile,
        "assay_geometry": assay_geometry,
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
            "structure_reliability": reliability['qc']['parameters'] if reliability else None,
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
            "stage1_construct_budget": {
                "requested_maximum": request.maximum_regions,
                "minimum": effective_min,
                "target": stage1_target,
                "maximum": effective_max,
                "fixed_by_assay_profile": bool(assay_profile),
                "policy": (
                    "fixed six-construct coarse-mapping ceiling; surface regions are enumerated "
                    "first and selected by geometric spread"
                    if assay_profile else
                    "3-7 is applied only when selecting the Stage-1 experimental panel; surface footprints are not capped"
                ),
            },
            "patch_engine": patch_engine,
            "joint_v3": config.JOINT_V3_PARAMETERS if patch_engine == "joint_v3" else None,
            "joint_v4_ml": config.JOINT_V4_ML_PARAMETERS if patch_engine == "joint_v4_ml" else None,
            "surface_v2": config.SURFACE_V2_PARAMETERS,
            "patch_scale": "Overlapping surface-geodesic footprints at fixed 9/12/15 A and adaptive 600 A2 scales",
            "compact_single_domain_surface_face": {
                "enabled_only_after_legacy_quality_failure": True,
                "ecd_max_length_aa": config.COMPACT_ECD_MAX_LENGTH_AA,
                "dominant_domain_min_fraction": config.COMPACT_DOMINANT_DOMAIN_MIN_FRACTION,
                "max_surface_regions": config.COMPACT_MAX_SURFACE_REGIONS,
                "footprint_geodesic_radius_A": config.COMPACT_FOOTPRINT_GEODESIC_RADIUS_A,
                "footprint_euclidean_cap_A": config.COMPACT_FOOTPRINT_EUCLIDEAN_CAP_A,
                "combo_max_spread_A": config.COMPACT_COMBO_MAX_SPREAD_A,
            },
            "regional_mutations": (
                "Joint V3 deterministic beam search with 2-4 sites; anchor mutation is not forced"
                if patch_engine == "joint_v3" else
                "Joint V4-ML order-specific 2-4 site search and footprint hypothesis panel; external scores are versioned evidence"
                if patch_engine == "joint_v4_ml" else
                "2-4 sites under fixed assay profile; directed alphabet + anchor mutation; "
                "exhaustive enumeration of cleared menus"
                if assay_profile else
                "directed alphabet + anchor mutation; exhaustive enumeration of cleared menus"
            ),
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
                item.status in {
                    "RECOMMENDED_LOW_PREDICTED_RISK",
                    "RECOMMENDED_MULTI_SITE_PROBE",
                    "RECOMMENDED_SINGLE_SITE_PROBE",
                }
                for item in designs
            ),
        },
        "residues": [to_dict(item) for item in records],
        "regions": [to_dict(item) for item in regions],
        "surface_v2": surface_v2_payload,
        "joint_v3": joint_v3_payload,
        "joint_v4_ml": joint_v4_payload,
        "structure_reliability": reliability,
        "designs": [to_dict(item) for item in designs],
        "warnings": warnings,
        "validation_gates": [
            "K562 antigen surface-display normalization to WT using N-terminal Myc",
            "independent non-competing conformation-sensitive antibody or ligand-binding control",
            "Jurkat Fab surface-expression and H/L assembly control",
            "WT, no-antigen, irrelevant-antigen and signaling-positive controls in the same batch",
            "native ligand/function control when functional-interface mutations are interpreted",
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
    if patch_engine == "joint_v3":
        for name in (
            "manifest.json", "structure_audit.json", "topology_context.json",
            "surface_residues.csv", "microfootprints.json", "probe_sectors.json",
            "candidate_probes.csv", "coverage_by_scenario.csv", "selected_panel.json", "exclusions.csv",
            "scientific_report.md",
        ):
            if (output_dir / name).exists():
                analysis["files"][name] = name
    if patch_engine == "joint_v4_ml":
        for name in (
            "tool_inventory.json", "structure_contexts.json", "residue_mapping.csv", "surface_graph_v4.json",
            "prediction_context.pdb",
            "surface_evidence.csv", "discotope_residues.csv", "single_mutation_options.csv",
            "candidate_combinations.csv.gz", "probe_sectors_v4.json", "panel_candidates.json",
            "selected_panel_v4.json", "surface_hypotheses_v4.json",
            "hypothesis_probe_matrix_v4.json", "v4_manifest.json", "v4_manifest_v41.json", "v4_scientific_report.md",
            "surface_graph_v41.json", "surface_display_regions_v41.json",
            "epitope_hypotheses_v41.json", "hypothesis_membership_v41.csv",
            "probe_hypothesis_matrix_v41.csv", "panel_identifiability_v41.json",
        ):
            if (output_dir / name).exists():
                analysis["files"][name] = name
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
    # The exporter creates the initial analysis file before it knows the
    # complete file map.  Persist the final analysis once more so V3 file
    # provenance and the legacy view agree on disk and in the API response.
    write_json(output_dir / "analysis.json", analysis)
    progress("COMPLETE", 1.0, "Analysis and epitope-mapping panel design completed")
    return analysis
