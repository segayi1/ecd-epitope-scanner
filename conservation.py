"""Per-residue evolutionary conservation from reviewed UniProt orthologs.

Epitope surfaces tend to be variable (immune evasion) and variable sites are
cheaper to mutate functionally, so conservation feeds both the epitope prior
and the site ranking (docs/改进方案.md).  The score is the fraction of aligned
reviewed same-gene orthologs carrying the identical residue; >= 3 orthologs
are required before any score is reported.  Everything degrades to None
offline or when the gene is unknown -- the pipeline states that plainly
instead of inventing a number.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any
from urllib.parse import quote

from Bio.Align import PairwiseAligner, substitution_matrices

from .models import ProteinContext
from .util import RemoteClient, write_json

UNIPROT_SEARCH_URL = (
    "https://rest.uniprot.org/uniprotkb/search?query={query}"
    "&fields=accession,gene_primary,sequence,organism_name,reviewed"
    "&format=json&size={size}"
)


def _primary_gene(protein: ProteinContext) -> str | None:
    genes = (protein.raw or {}).get("genes") or []
    for gene in genes:
        name = (gene.get("geneName") or {}).get("value")
        if name:
            return name
    return None


def _aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def fetch_conservation(
    protein: ProteinContext,
    client: RemoteClient,
    output_dir: Path,
) -> dict[int, float] | None:
    """Return {uniprot_pos: conservation in 0..1}, or None when unavailable."""
    cache = output_dir / "conservation.json"
    if cache.exists():
        import json

        try:
            stored = json.loads(cache.read_text(encoding="utf-8"))
        except Exception:
            stored = None
        if isinstance(stored, dict) and stored.get("status") == "computed" and stored.get("scores") is not None:
            return {int(key): float(value) for key, value in stored["scores"].items()}

    gene = _primary_gene(protein)
    if not gene:
        write_json(cache, {"status": "gene_unknown", "scores": None})
        return None
    try:
        query_text = f'(gene_exact:"{gene}") AND reviewed:true'
        query = quote(query_text)
        payload = client.json(UNIPROT_SEARCH_URL.format(query=query, size=config_size()))
    except Exception as exc:
        write_json(
            cache,
            {"status": "query_failed", "gene": gene, "error": type(exc).__name__, "scores": None},
        )
        return None
    if not isinstance(payload, dict):
        write_json(cache, {"status": "invalid_payload", "gene": gene, "scores": None})
        return None
    results = payload.get("results")
    if not isinstance(results, list):
        write_json(cache, {"status": "invalid_payload", "gene": gene, "scores": None})
        return None
    orthologs: list[dict[str, Any]] = []
    for item in results:
        if not isinstance(item, dict) or item.get("primaryAccession") == protein.accession:
            continue
        sequence = item.get("sequence")
        if isinstance(sequence, dict) and isinstance(sequence.get("value"), str) and sequence["value"]:
            orthologs.append(item)
    from . import config

    if len(orthologs) < config.CONSERVATION_MIN_ORTHOLOGS:
        write_json(
            cache,
            {
                "status": "insufficient_orthologs",
                "gene": gene,
                "query": query_text,
                "orthologs_found": len(orthologs),
                "scores": None,
            },
        )
        return None

    aligner = _aligner()
    counts: dict[int, int] = {}
    totals: dict[int, int] = {}
    used = 0
    for ortholog in orthologs[: config.CONSERVATION_MAX_ORTHOLOGS]:
        other = str((ortholog.get("sequence") or {}).get("value") or "").replace("\n", "")
        if not other:
            continue
        try:
            alignment = next(iter(aligner.align(protein.sequence, other)))
            indices = alignment.indices
        except Exception:
            continue
        used += 1
        for target_index, source_index in zip(indices[0], indices[1]):
            if target_index < 0 or source_index < 0:
                continue
            position = target_index + 1
            totals[position] = totals.get(position, 0) + 1
            if protein.sequence[target_index] == other[source_index]:
                counts[position] = counts.get(position, 0) + 1
    if used < config.CONSERVATION_MIN_ORTHOLOGS:
        write_json(
            cache,
            {
                "status": "insufficient_aligned_orthologs",
                "gene": gene,
                "query": query_text,
                "orthologs_found": len(orthologs),
                "orthologs_used": used,
                "scores": None,
            },
        )
        return None
    scores = {
        position: round(counts.get(position, 0) / totals[position], 4)
        for position in totals
        if totals[position] >= max(config.CONSERVATION_MIN_ORTHOLOGS, used // 2)
    }
    if not scores:
        write_json(
            cache,
            {
                "status": "no_usable_alignment_columns",
                "gene": gene,
                "query": query_text,
                "orthologs_used": used,
                "scores": None,
            },
        )
        return None
    write_json(
        cache,
        {
            "status": "computed",
            "gene": gene,
            "query": query_text,
            "orthologs_used": used,
            "alignment": {
                "algorithm": getattr(aligner, "algorithm", "PairwiseAligner global"),
                "matrix": "BLOSUM62",
                "open_gap": -11,
                "extend_gap": -1,
            },
            "scores": {str(key): value for key, value in scores.items()},
        },
    )
    return scores


def config_size() -> int:
    from . import config

    return config.CONSERVATION_MAX_ORTHOLOGS + 4
