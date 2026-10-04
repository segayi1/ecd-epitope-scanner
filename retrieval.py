from __future__ import annotations

import re
from pathlib import Path
from typing import Any
from urllib.parse import quote

from .models import FeatureRange, ProteinContext, StructureCandidate, to_dict
from .util import PipelineError, RemoteClient, safe_name, write_json

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_SEARCH_URL = (
    "https://rest.uniprot.org/uniprotkb/search?query={query}"
    "&fields=accession,gene_primary,protein_name,organism_name,reviewed&format=json&size=6"
)
BEACONS_URL = "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/{accession}.json"
SIFTS_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
PDB_CIF_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif"


def _location_value(location: dict[str, Any], side: str) -> int | None:
    value = location.get(side, {})
    raw = value.get("value")
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def parse_uniprot_entry(
    entry: dict[str, Any], ecd_override: tuple[int, int] | None = None
) -> ProteinContext:
    sequence = entry.get("sequence", {}).get("value", "")
    accession = entry.get("primaryAccession", "")
    if not accession or not sequence:
        raise PipelineError("UniProt entry is missing a primary accession or sequence")

    features: list[FeatureRange] = []
    for raw in entry.get("features", []):
        start = _location_value(raw.get("location", {}), "start")
        end = _location_value(raw.get("location", {}), "end")
        if start is None or end is None:
            continue
        evidence = ",".join(
            item.get("evidenceCode", "") for item in raw.get("evidences", []) if item
        )
        features.append(
            FeatureRange(
                type=raw.get("type", "Feature"),
                start=start,
                end=end,
                description=raw.get("description", "") or "",
                evidence=evidence,
            )
        )

    transmembranes = [item for item in features if item.type.lower() == "transmembrane"]
    if len(transmembranes) != 1:
        raise PipelineError(
            f"Expected a single-pass membrane protein, found {len(transmembranes)} transmembrane annotations"
        )

    extracellular = [
        item
        for item in features
        if item.type.lower() == "topological domain"
        and re.search(r"extracellular|lumenal|luminal", item.description, re.I)
    ]
    if ecd_override:
        ecd_start, ecd_end = ecd_override
    elif extracellular:
        chosen = max(extracellular, key=lambda item: item.end - item.start)
        ecd_start, ecd_end = chosen.start, chosen.end
    else:
        tm = transmembranes[0]
        signal = [item for item in features if item.type.lower() == "signal peptide"]
        if signal and tm.start > signal[0].end:
            ecd_start, ecd_end = signal[0].end + 1, tm.start - 1
        else:
            ecd_start, ecd_end = tm.end + 1, len(sequence)

    if not (1 <= ecd_start <= ecd_end <= len(sequence)):
        raise PipelineError(f"Invalid ECD range {ecd_start}:{ecd_end} for sequence length {len(sequence)}")

    tm = transmembranes[0]
    if tm.start > ecd_end:
        topology = "type-I / N-terminal extracellular"
        membrane_anchor = ecd_end
    else:
        topology = "type-II / C-terminal extracellular"
        membrane_anchor = ecd_start

    description = entry.get("proteinDescription", {})
    protein_name = (
        description.get("recommendedName", {}).get("fullName", {}).get("value")
        or description.get("submissionNames", [{}])[0].get("fullName", {}).get("value")
        or accession
    )
    pdb_ids = [
        item.get("id", "").lower()
        for item in entry.get("uniProtKBCrossReferences", [])
        if item.get("database") == "PDB" and item.get("id")
    ]
    organism = entry.get("organism", {}).get("scientificName", "")
    reviewed = "reviewed" in entry.get("entryType", "").lower()
    return ProteinContext(
        accession=accession,
        protein_name=protein_name,
        organism=organism,
        sequence=sequence,
        ecd_start=ecd_start,
        ecd_end=ecd_end,
        membrane_anchor_position=membrane_anchor,
        topology=topology,
        reviewed=reviewed,
        features=features,
        pdb_ids=pdb_ids,
        source_url=UNIPROT_URL.format(accession=accession),
        raw=entry,
    )


def _gene_candidates(identifier: str, client: RemoteClient, organism_id: int = 9606) -> list[dict[str, Any]]:
    query = quote(f'(gene_exact:"{identifier}") AND (organism_id:{organism_id})')
    payload = client.json(UNIPROT_SEARCH_URL.format(query=query))
    return [
        item for item in (payload.get("results") or [])
        if item.get("primaryAccession") and item.get("organism", {}).get("taxonId") == organism_id
    ]


def resolve_uniprot_identifier(
    identifier: str, client: RemoteClient, organism_id: int = 9606
) -> tuple[str, str]:
    """Resolve an accession or an unambiguous gene symbol.

    The returned tuple is (canonical_accession, input_kind).  Gene symbols are
    intentionally restricted to human by default: silently selecting a
    same-named gene from another species would be worse than a useful error.
    """
    raw = identifier.strip()
    if not raw:
        raise PipelineError("A UniProt accession or gene symbol is required")
    accession = raw.upper()
    try:
        entry = client.json(UNIPROT_URL.format(accession=accession))
    except PipelineError:
        entry = None
    if isinstance(entry, dict) and entry.get("primaryAccession"):
        return str(entry["primaryAccession"]).upper(), "accession"
    candidates = _gene_candidates(raw, client, organism_id)
    if len(candidates) == 1:
        return str(candidates[0]["primaryAccession"]).upper(), "gene_symbol"
    if not candidates:
        raise PipelineError(
            f"No human UniProt accession or unambiguous gene symbol match was found for {raw!r}."
        )
    choices = ", ".join(
        f"{item.get('primaryAccession')} ({item.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'candidate')})"
        for item in candidates
    )
    raise PipelineError(
        f"Gene symbol {raw!r} is ambiguous for Homo sapiens; choose one accession explicitly: {choices}"
    )


def fetch_protein_context(
    accession: str,
    client: RemoteClient,
    output_dir: Path,
    ecd_override: tuple[int, int] | None = None,
) -> ProteinContext:
    canonical, input_kind = resolve_uniprot_identifier(accession, client)
    entry = client.json(UNIPROT_URL.format(accession=canonical))
    if input_kind != "accession":
        entry["_workbench_input"] = {"original_identifier": accession, "kind": input_kind}
    write_json(output_dir / "uniprot.json", entry)
    return parse_uniprot_entry(entry, ecd_override)


def _candidate_from_beacon(item: dict[str, Any]) -> StructureCandidate:
    summary = item.get("summary", item)
    entities = summary.get("entities") or []
    chain_ids: list[str] = []
    for entity in entities:
        if entity.get("entity_type") == "POLYMER":
            chain_ids.extend(entity.get("chain_ids") or [])
    text = " ".join(
        str(summary.get(key, ""))
        for key in ("model_identifier", "model_type", "confidence_version", "provider")
    )
    explicit_af3 = bool(re.search(r"\b(?:AF3|AlphaFold\s*3)\b", text, re.I))
    resolution = summary.get("resolution")
    try:
        resolution = float(resolution) if resolution is not None else None
    except (TypeError, ValueError):
        resolution = None
    return StructureCandidate(
        identifier=str(summary.get("model_identifier", "unknown")),
        category=str(summary.get("model_category", "")),
        provider=str(summary.get("provider", "")),
        model_url=str(summary.get("model_url", "")),
        page_url=str(summary.get("model_page_url", "") or ""),
        model_format=str(summary.get("model_format", "MMCIF") or "MMCIF"),
        method=str(summary.get("experimental_method", "") or summary.get("model_type", "") or ""),
        resolution=resolution,
        sequence_identity=summary.get("sequence_identity"),
        uniprot_start=summary.get("uniprot_start"),
        uniprot_end=summary.get("uniprot_end"),
        reported_coverage=summary.get("coverage"),
        chain_ids=sorted(set(chain_ids)),
        created=str(summary.get("created", "") or ""),
        explicit_af3=explicit_af3,
        metadata=summary,
    )


def fetch_structure_candidates(
    protein: ProteinContext, client: RemoteClient, output_dir: Path
) -> list[StructureCandidate]:
    payload = client.json(BEACONS_URL.format(accession=protein.accession))
    write_json(output_dir / "3d_beacons.json", payload)
    candidates = [_candidate_from_beacon(item) for item in payload.get("structures", [])]
    return candidates


def prefilter_experimental_candidates(
    protein: ProteinContext, candidates: list[StructureCandidate], limit: int = 10
) -> list[StructureCandidate]:
    scored: list[tuple[tuple[Any, ...], StructureCandidate]] = []
    for candidate in candidates:
        if candidate.category.upper() != "EXPERIMENTALLY DETERMINED":
            continue
        start = candidate.uniprot_start or 1
        end = candidate.uniprot_end or 0
        overlap = max(0, min(end, protein.ecd_end) - max(start, protein.ecd_start) + 1)
        coarse_coverage = overlap / protein.ecd_length
        identity = float(candidate.sequence_identity or 0)
        if identity > 1.5:
            identity /= 100.0
        resolution = candidate.resolution if candidate.resolution is not None else 99.0
        candidate.metadata["coarse_ecd_coverage"] = round(coarse_coverage, 4)
        key = (-coarse_coverage, -identity, resolution, candidate.identifier)
        scored.append((key, candidate))
    return [item for _, item in sorted(scored, key=lambda pair: pair[0])[:limit]]


def find_public_af3_candidates(
    candidates: list[StructureCandidate], client: RemoteClient
) -> list[StructureCandidate]:
    accepted: list[StructureCandidate] = []
    for candidate in candidates:
        if candidate.category.upper() == "EXPERIMENTALLY DETERMINED":
            continue
        if candidate.provider.lower() == "alphafold db":
            continue
        if not candidate.explicit_af3 and candidate.provider.lower() == "modelarchive" and candidate.page_url:
            try:
                page = client.text(candidate.page_url)
                candidate.explicit_af3 = bool(re.search(r"\b(?:AF3|AlphaFold\s*3)\b", page, re.I))
            except PipelineError:
                candidate.explicit_af3 = False
        if candidate.explicit_af3:
            accepted.append(candidate)
    return accepted


def fetch_sifts_mapping(pdb_id: str, client: RemoteClient, output_dir: Path) -> dict[str, Any]:
    try:
        payload = client.json(SIFTS_URL.format(pdb_id=pdb_id.lower()))
    except PipelineError:
        payload = {}
    write_json(output_dir / f"sifts_{safe_name(pdb_id.lower())}.json", payload)
    return payload


def infer_uniprot_from_sifts(
    payload: dict[str, Any],
    pdb_id: str,
    preferred_accession: str | None = None,
    preferred_chain: str | None = None,
) -> tuple[str, list[str]]:
    root = payload.get(pdb_id.lower(), payload.get(pdb_id.upper(), {}))
    entries = root.get("UniProt", {})
    if preferred_accession:
        data = entries.get(preferred_accession)
        if not data:
            return preferred_accession, []
        chains = sorted({
            item.get("chain_id")
            for item in data.get("mappings", [])
            if item.get("chain_id")
        })
        return preferred_accession, chains
    if not entries:
        raise PipelineError(f"No UniProt mapping was found for PDB {pdb_id}")
    if preferred_chain:
        matching = [
            (accession, data)
            for accession, data in entries.items()
            if any(
                item.get("chain_id") == preferred_chain
                for item in data.get("mappings", [])
            )
        ]
        if matching:
            accession, data = max(
                matching,
                key=lambda pair: sum(
                    int(item.get("unp_end", 0)) - int(item.get("unp_start", 0)) + 1
                    for item in pair[1].get("mappings", [])
                ),
            )
            chains = sorted({
                item.get("chain_id")
                for item in data.get("mappings", [])
                if item.get("chain_id")
            })
            return accession, chains
    ranked = sorted(
        entries.items(),
        key=lambda pair: -sum(
            int(item.get("unp_end", 0)) - int(item.get("unp_start", 0)) + 1
            for item in pair[1].get("mappings", [])
        ),
    )
    accession, data = ranked[0]
    chains = sorted({item.get("chain_id") for item in data.get("mappings", []) if item.get("chain_id")})
    return accession, chains


def direct_pdb_candidate(pdb_id: str) -> StructureCandidate:
    pdb_id = pdb_id.strip().lower()
    if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", pdb_id):
        raise PipelineError(f"Invalid PDB identifier: {pdb_id}")
    return StructureCandidate(
        identifier=pdb_id,
        category="EXPERIMENTALLY DETERMINED",
        provider="PDBe",
        model_url=PDB_CIF_URL.format(pdb_id=pdb_id),
        page_url=f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}",
        model_format="MMCIF",
    )


def candidate_manifest(candidates: list[StructureCandidate]) -> list[dict[str, Any]]:
    return [to_dict(item) for item in candidates]
