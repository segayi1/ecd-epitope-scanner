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
        if isinstance(stored, dict) and stored.get("scores") is not None:
            return {int(key): float(value) for key, value in stored["scores"].items()}

    gene = _primary_gene(protein)
    if not gene:
        return None
    try:
        query = f'(gene_exact:"{gene}") AND reviewed:true'
        payload = client.json(UNIPROT_SEARCH_URL.format(query=query, size=config_size()))
    except Exception:
        return None
    results = payload.get("results", [])
    orthologs = [
        item
        for item in results
        if item.get("primaryAccession") != protein.accession and item.get("sequence", {}).get("value")
    ]
    if len(orthologs) < 3:
        return None

    aligner = _aligner()
    counts: dict[int, int] = {}
    totals: dict[int, int] = {}
    used = 0
    for ortholog in orthologs[:8]:
        other = ortholog["sequence"]["value"].replace("\n", "")
        try:
            alignment = next(iter(aligner.align(protein.sequence, other)))
        except Exception:
            continue
        used += 1
        for target_index, source_index in zip(alignment.indices[0], alignment.indices[1]):
            if target_index < 0 or source_index < 0:
                continue
            position = target_index + 1
            totals[position] = totals.get(position, 0) + 1
            if protein.sequence[target_index] == other[source_index]:
                counts[position] = counts.get(position, 0) + 1
    if used < 3:
        return None
    scores = {
        position: round(counts[position] / totals[position], 4)
        for position in totals
        if totals[position] >= max(3, used // 2)
    }
    if not scores:
        return None
    write_json(
        cache,
        {
            "gene": gene,
            "orthologs_used": used,
            "scores": {str(key): value for key, value in scores.items()},
        },
    )
    return scores


def config_size() -> int:
    from . import config

    return config.CONSERVATION_MAX_ORTHOLOGS + 4
