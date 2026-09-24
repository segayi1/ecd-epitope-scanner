from __future__ import annotations

import heapq
from itertools import combinations

from ..models import ProteinContext, ResidueRecord
from .models import MutationOption, ProbeCandidate, ProbeSector
from .mutation_options import build_mutation_options


def _graph_distances(
    graph: dict[int, dict[int, float]], start: int, allowed: set[int]
) -> dict[int, float]:
    distances = {start: 0.0}
    queue = [(0.0, start)]
    while queue:
        distance, current = heapq.heappop(queue)
        if distance != distances.get(current):
            continue
        for neighbor, step in graph.get(current, {}).items():
            if neighbor not in allowed:
                continue
            proposed = distance + step
            if proposed < distances.get(neighbor, float("inf")):
                distances[neighbor] = proposed
                heapq.heappush(queue, (proposed, neighbor))
    return distances


def _combo_geometry(
    combo: tuple[MutationOption, ...],
    sector_positions: set[int],
    graph: dict[int, dict[int, float]],
) -> tuple[bool, float, float]:
    positions = sorted(item.position for item in combo)
    pair_distances: list[float] = []
    for index, left in enumerate(positions):
        distances = _graph_distances(graph, left, sector_positions)
        for right in positions[index + 1 :]:
            value = distances.get(right, float("inf"))
            if value == float("inf"):
                return False, value, value
            pair_distances.append(value)
    return True, min(pair_distances, default=0.0), max(pair_distances, default=0.0)


def _candidate_score(
    combo: tuple[MutationOption, ...],
    sector_positions: set[int],
    graph: dict[int, dict[int, float]],
) -> tuple[float, ...]:
    connected, minimum_distance, maximum_distance = _combo_geometry(combo, sector_positions, graph)
    if not connected:
        return (-float("inf"),)
    strength = sum(item.chemical_perturbation_strength for item in combo)
    risk = sum(item.rule_risk for item in combo)
    novelty = min(minimum_distance, 8.0) / 40.0
    compactness_penalty = max(0.0, maximum_distance - 12.0) / 30.0
    return (
        round(strength + novelty - compactness_penalty - risk, 6),
        -round(risk, 6),
        -round(maximum_distance, 6),
        -len(combo),
    )


def _signature(combo: tuple[MutationOption, ...]) -> tuple[tuple[int, str], ...]:
    return tuple(sorted((item.position, item.mutant_aa) for item in combo))


def _canonical(combo: tuple[MutationOption, ...]) -> tuple[MutationOption, ...]:
    return tuple(sorted(combo, key=lambda item: (item.position, item.mutant_aa)))


def _creates_new_sequon(combo: tuple[MutationOption, ...], protein: ProteinContext) -> bool:
    original = protein.sequence
    mutated = list(original)
    for item in combo:
        index = item.position - 1
        if 0 <= index < len(mutated):
            mutated[index] = item.mutant_aa

    def motifs(sequence: str) -> set[int]:
        return {
            start
            for start in range(len(sequence) - 2)
            if sequence[start] == "N" and sequence[start + 1] != "P" and sequence[start + 2] in "ST"
        }

    return bool(motifs("".join(mutated)) - motifs(original))


def _rank_unique(
    combos: list[tuple[MutationOption, ...]],
    sector_positions: set[int],
    graph: dict[int, dict[int, float]],
    maximum_geodesic_span_a: float,
    limit: int,
) -> list[tuple[MutationOption, ...]]:
    unique: dict[tuple[tuple[int, str], ...], tuple[MutationOption, ...]] = {}
    for raw in combos:
        combo = _canonical(raw)
        signature = _signature(combo)
        if signature in unique:
            continue
        connected, _, span = _combo_geometry(combo, sector_positions, graph)
        if connected and span <= maximum_geodesic_span_a:
            unique[signature] = combo
    ranked = sorted(
        unique.values(),
        key=lambda combo: (_candidate_score(combo, sector_positions, graph), _signature(combo)),
        reverse=True,
    )
    return ranked[:limit]


def _balanced_shortlist(
    layers: dict[int, list[tuple[MutationOption, ...]]], limit: int
) -> list[tuple[MutationOption, ...]]:
    output: list[tuple[MutationOption, ...]] = []
    offsets = {size: 0 for size in layers}
    while len(output) < limit:
        changed = False
        for size in sorted(layers):
            index = offsets[size]
            if index < len(layers[size]) and len(output) < limit:
                output.append(layers[size][index])
                offsets[size] += 1
                changed = True
        if not changed:
            break
    return output


def generate_probe_candidates(
    sectors: list[ProbeSector],
    records: list[ResidueRecord],
    protein: ProteinContext,
    graph: dict[int, dict[int, float]],
    designable_positions: set[int],
    site_pool_cap: int = 12,
    options_per_site: int = 3,
    beam_width: int = 200,
    shortlist_per_sector: int = 20,
    min_mutations: int = 2,
    max_mutations: int = 4,
    maximum_geodesic_span_a: float = 24.0,
) -> tuple[list[ProbeCandidate], dict[str, list[str]]]:
    by_position = {item.uniprot_pos: item for item in records}
    all_candidates: list[ProbeCandidate] = []
    exclusions: dict[str, list[str]] = {}
    for sector in sectors:
        sector_positions = set(sector.residue_ids)
        options, rejected = build_mutation_options(
            sector,
            by_position,
            protein,
            options_per_site,
            allowed_positions=designable_positions,
        )
        ranked_sites = sorted(
            {item.position for item in options},
            key=lambda position: (
                -max(item.chemical_perturbation_strength for item in options if item.position == position),
                min(item.rule_risk for item in options if item.position == position),
                position,
            ),
        )[:site_pool_cap]
        ranked_site_set = set(ranked_sites)
        option_pool = [item for item in options if item.position in ranked_site_set]
        if rejected:
            exclusions[sector.sector_id] = [
                f"{position}:{','.join(sorted(set(reasons)))}"
                for position, reasons in sorted(rejected.items())
            ]
        if len(ranked_sites) < min_mutations:
            exclusions.setdefault(sector.sector_id, []).append("insufficient-designable-sites")
            continue

        initial = [
            tuple(combo)
            for combo in combinations(option_pool, min_mutations)
            if len({item.position for item in combo}) == min_mutations
            and not _creates_new_sequon(tuple(combo), protein)
        ]
        beam = _rank_unique(
            initial, sector_positions, graph, maximum_geodesic_span_a, beam_width
        )
        layers: dict[int, list[tuple[MutationOption, ...]]] = {}
        if beam:
            layers[min_mutations] = list(beam)
        for size in range(min_mutations + 1, max_mutations + 1):
            expanded = [
                combo + (option,)
                for combo in beam
                for option in option_pool
                if option.position not in {item.position for item in combo}
                and not _creates_new_sequon(combo + (option,), protein)
            ]
            beam = _rank_unique(
                expanded, sector_positions, graph, maximum_geodesic_span_a, beam_width
            )
            if not beam:
                break
            layers[size] = list(beam)

        shortlisted = _balanced_shortlist(layers, shortlist_per_sector)
        for rank, combo in enumerate(shortlisted, start=1):
            connected, _, _ = _combo_geometry(combo, sector_positions, graph)
            all_candidates.append(
                ProbeCandidate(
                    candidate_id=f"PC_{sector.sector_id}_{rank:03d}",
                    sector_id=sector.sector_id,
                    mutations=list(combo),
                    connected=connected,
                    structural_risk=round(sum(item.rule_risk for item in combo) / len(combo), 4),
                    evidence_level="PROXY_ONLY",
                    warnings=[
                        "No external combination-level stability engine was run; risk is a dimensionless rule proxy."
                    ],
                )
            )
    all_candidates.sort(
        key=lambda item: (
            item.sector_id,
            tuple((mutation.position, mutation.mutant_aa) for mutation in item.mutations),
        )
    )
    globally_unique: dict[tuple[tuple[int, str], ...], ProbeCandidate] = {}
    for candidate in all_candidates:
        signature = tuple((item.position, item.mutant_aa) for item in candidate.mutations)
        globally_unique.setdefault(signature, candidate)
    return list(globally_unique.values()), exclusions
