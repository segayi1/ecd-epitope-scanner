"""Fab-footprint-scale surface patch partitioning for epitope mapping.

Replaces the old domain-quota geodesic partitioning.  Rationale
(docs/改进方案.md): conformational epitopes span 15-25 residues / ~600-900
A2 (Kringelum 2012) and routinely cross domain boundaries, so patches are
enumerated at antibody-footprint scale directly on the surface graph instead
of being forced into domain quotas.  Each patch carries an epitope prior
(exposure, protrusion, loop content, variability, glycan/membrane penalties)
that drives the Stage-1 probe panel: HIGH-prior patches get two probes,
others one, and anything untestable is declared a blind spot rather than
silently skipped.
"""
from __future__ import annotations

import math
from collections import Counter
from statistics import mean
from typing import Any

import networkx as nx
import numpy as np

from . import config
from .models import ProteinContext, RegionRecord, ResidueRecord
from .structure_analysis import median_region_pae


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))


def _build_graph(records: list[ResidueRecord]) -> nx.Graph:
    graph = nx.Graph()
    for record in records:
        graph.add_node(record.uniprot_pos)
    by_pos = {record.uniprot_pos: record for record in records}
    for index, left in enumerate(records):
        if left.ca is None:
            continue
        for right in records[index + 1 :]:
            if right.ca is None:
                continue
            distance = _distance(left.ca, right.ca)
            sequence_neighbor = abs(left.uniprot_pos - right.uniprot_pos) <= 2
            if distance <= 12.5 or (sequence_neighbor and distance <= 16.0):
                penalty = 4.0 if left.domain != right.domain else 0.0
                secondary_penalty = 1.5 if left.secondary_structure != right.secondary_structure else 0.0
                graph.add_edge(
                    left.uniprot_pos,
                    right.uniprot_pos,
                    weight=distance + penalty + secondary_penalty,
                    euclidean=distance,
                )

    # Surface islands can arise around flexible loops.  Connect each island to
    # its nearest neighbor so every eligible residue is auditable and assigned
    # exactly once.
    components = list(nx.connected_components(graph))
    while len(components) > 1:
        first = components[0]
        rest = set().union(*components[1:])
        best: tuple[float, int, int] | None = None
        for left_pos in first:
            for right_pos in rest:
                distance = _distance(by_pos[left_pos].ca, by_pos[right_pos].ca)  # type: ignore[arg-type]
                candidate = (distance, left_pos, right_pos)
                if best is None or candidate < best:
                    best = candidate
        if best is None:
            break
        distance, left_pos, right_pos = best
        domain_penalty = 4.0 if by_pos[left_pos].domain != by_pos[right_pos].domain else 0.0
        graph.add_edge(left_pos, right_pos, weight=distance + domain_penalty + 6.0, euclidean=distance)
        components = list(nx.connected_components(graph))
    return graph


def _exposure_priority(record: ResidueRecord) -> float:
    score = record.rsasa * 2.0 + min(1.2, record.sidechain_sasa / 100.0)
    return score + 0.6 / (1 + record.contacts / 6.0)


def _enumerate_patch_seeds(
    graph: nx.Graph, records: list[ResidueRecord], radius: float
) -> list[int]:
    """Farthest-point seeding until every residue is within reach of a seed.

    "Reach" = geodesic distance <= `radius` AND euclidean distance <=
    PATCH_EUCLIDEAN_CAP_A.  The euclidean condition stops patches from
    wrapping around thin sheets: a residue 12 A along the surface but 40 A
    through space belongs to a different face, not this patch.
    """
    by_pos = {record.uniprot_pos: record for record in records}
    nodes = sorted(by_pos)
    first = max(
        nodes,
        key=lambda position: (_exposure_priority(by_pos[position]), -position),
    )
    seeds = [first]
    distance_to_seeds = nx.single_source_dijkstra_path_length(graph, first, weight="weight")

    def effective_distances() -> dict[int, float]:
        output: dict[int, float] = {}
        for node in nodes:
            if node in seeds:
                continue
            reachable = [
                distance_to_seeds[node]
                for seed in seeds
                if _distance(by_pos[node].ca, by_pos[seed].ca) <= config.PATCH_EUCLIDEAN_CAP_A  # type: ignore[arg-type]
                and distance_to_seeds.get(node, math.inf) <= radius
            ]
            output[node] = min(reachable) if reachable else math.inf
        return output

    effective = effective_distances()
    while True:
        candidates = {node: value for node, value in effective.items() if node not in seeds}
        if not candidates:
            break
        next_seed, worst = max(candidates.items(), key=lambda item: (item[1], -item[0]))
        if worst <= radius:
            break
        seeds.append(next_seed)
        fresh = nx.single_source_dijkstra_path_length(graph, next_seed, weight="weight")
        distance_to_seeds = {
            node: min(distance_to_seeds.get(node, math.inf), fresh.get(node, math.inf))
            for node in nodes
        }
        effective = effective_distances()
    return seeds


def _assign_to_patches(
    graph: nx.Graph, records: list[ResidueRecord], seeds: list[int]
) -> dict[int, list[int]]:
    by_pos = {record.uniprot_pos: record for record in records}
    seed_distances = {
        seed: nx.single_source_dijkstra_path_length(graph, seed, weight="weight")
        for seed in seeds
    }
    assignments: dict[int, list[int]] = {seed: [] for seed in seeds}
    for record in records:
        local_seeds = [
            seed for seed in seeds
            if _distance(record.ca, by_pos[seed].ca) <= config.PATCH_EUCLIDEAN_CAP_A  # type: ignore[arg-type]
        ]
        if local_seeds:
            best_distance, best_seed = min(
                ((seed_distances[seed].get(record.uniprot_pos, math.inf), seed) for seed in local_seeds),
                key=lambda item: (item[0], item[1]),
            )
        else:
            # Stray residue beyond the euclidean cap of every seed: fall back
            # to the spatially nearest seed so coverage stays complete while
            # the patch stays locally compact.
            best_seed = min(
                seeds,
                key=lambda seed: (_distance(record.ca, by_pos[seed].ca), seed),  # type: ignore[arg-type]
            )
            best_distance = 0.0
        if math.isinf(best_distance):
            best_seed = min(
                seeds,
                key=lambda seed: (_distance(record.ca, by_pos[seed].ca), seed),  # type: ignore[arg-type]
            )
        assignments[best_seed].append(record.uniprot_pos)
    return assignments


def _merge_small_patches(
    assignments: dict[int, list[int]], by_pos: dict[int, ResidueRecord]
) -> dict[int, list[int]]:
    for seed in list(assignments):
        members = assignments.get(seed, [])
        if len(members) >= config.PATCH_MIN_RESIDUES or len(assignments) == 1:
            continue
        alternatives = [item for item in assignments if item != seed]
        if not alternatives:
            continue
        replacement = min(
            alternatives,
            key=lambda other: (_distance(by_pos[seed].ca, by_pos[other].ca), other),  # type: ignore[arg-type]
        )
        assignments[replacement].extend(members)
        del assignments[seed]
    return assignments


def _merge_down_to_limit(
    assignments: dict[int, list[int]],
    priors: dict[int, float],
    by_pos: dict[int, ResidueRecord],
    limit: int,
) -> tuple[dict[int, list[int]], list[int]]:
    """Fold the weakest-prior patch into its spatial nearest neighbor.

    Returns the new assignment map plus the seeds that absorbed a neighbour
    (so rationales can flag patches whose extent grew beyond footprint scale).
    """
    merged_targets: list[int] = []
    while len(assignments) > limit:
        weakest = min(assignments, key=lambda seed: (priors.get(seed, 0.0), seed))
        alternatives = [item for item in assignments if item != weakest]
        replacement = min(
            alternatives,
            key=lambda other: (_distance(by_pos[weakest].ca, by_pos[other].ca), other),  # type: ignore[arg-type]
        )
        priors[replacement] = min(priors.get(replacement, 0.0), priors.get(weakest, 0.0))
        assignments[replacement].extend(assignments[weakest])
        del assignments[weakest]
        merged_targets.append(replacement)
    return assignments, merged_targets


def _patch_prior(items: list[ResidueRecord]) -> tuple[float, dict[str, float]]:
    exposure = mean(min(1.0, item.rsasa / 0.5) for item in items)
    protrusion = mean(1.0 / (1.0 + item.contacts / 6.0) for item in items)
    loop = sum(item.secondary_structure == "C" for item in items) / len(items)
    conservation_values = [item.conservation_score for item in items if item.conservation_score is not None]
    if conservation_values:
        variability = 1.0 - mean(conservation_values)
    else:
        variability = 0.4  # neutral placeholder when conservation is unavailable
    glycan = sum(item.glycan_shielded for item in items) / len(items)
    interface = sum(item.interface for item in items) / len(items)
    prior = (
        config.PRIOR_W_EXPOSURE * exposure
        + config.PRIOR_W_PROTRUSION * protrusion
        + config.PRIOR_W_LOOP * loop
        + config.PRIOR_W_VARIABILITY * variability
        + config.PRIOR_W_INTERFACE * min(1.0, interface * 2.0)
        - config.PRIOR_PENALTY_GLYCAN * glycan
    )
    components = {
        "exposure": round(exposure, 3),
        "protrusion": round(protrusion, 3),
        "loop_fraction": round(loop, 3),
        "variability": round(variability, 3),
        "glycan_fraction": round(glycan, 3),
        "interface_fraction": round(interface, 3),
        "conservation_available": 1.0 if conservation_values else 0.0,
    }
    return round(prior, 4), components


def _assign_tiers(priors: dict[int, float]) -> dict[int, str]:
    ordered = sorted(priors, key=lambda seed: (-priors[seed], seed))
    total = len(ordered)
    tiers: dict[int, str] = {}
    for rank, seed in enumerate(ordered):
        if rank < max(1, round(total * config.TIER_HIGH_FRACTION)):
            tiers[seed] = "HIGH"
        elif rank < max(1, round(total * (config.TIER_HIGH_FRACTION + config.TIER_MEDIUM_FRACTION))):
            tiers[seed] = "MEDIUM"
        else:
            tiers[seed] = "LOW"
    return tiers


def _geodesic_diameter(graph: nx.Graph, members: list[int]) -> float:
    subgraph = graph.subgraph(members)
    maximum = 0.0
    for _, lengths in nx.all_pairs_dijkstra_path_length(subgraph, weight="weight"):
        finite = [float(value) for value in lengths.values() if math.isfinite(value)]
        if finite:
            maximum = max(maximum, max(finite))
    return round(maximum, 2)


def _structural_character(items: list[ResidueRecord]) -> str:
    counts = Counter(item.secondary_structure for item in items)
    total = max(1, len(items))
    if counts["C"] / total >= 0.55:
        return "loop-rich surface"
    if counts["E"] / total >= 0.45:
        return "beta-sheet face"
    if counts["H"] / total >= 0.45:
        return "helical face"
    return "mixed secondary-structure face"


def partition_surface(
    records: list[ResidueRecord],
    protein: ProteinContext,
    confidence: dict[str, Any] | None = None,
    minimum_regions: int = 5,
    maximum_regions: int = 8,
    predicted: bool = False,
    assay_format: str = "cell_surface",
) -> tuple[list[RegionRecord], dict[str, Any]]:
    eligible = [
        record for record in records
        if record.folding_screen_pass and record.ca is not None
        and not (config.EXCLUDE_FLEXIBLE_TERMINAL and record.flexible_terminal)
    ]
    excluded_tail = sum(
        1 for record in records
        if record.folding_screen_pass and record.ca is not None and record.flexible_terminal
    )
    if not eligible:
        return [], {
            "eligible_residues": 0,
            "assigned_residues": 0,
            "coverage": 0.0,
            "excluded_flexible_terminal_residues": excluded_tail,
            "warnings": ["No residues passed the configured surface and folding-risk screen."],
            "method": "Fab-footprint patch partitioning with epitope priors",
        }

    graph = _build_graph(eligible)
    by_pos = {record.uniprot_pos: record for record in eligible}
    seeds = _enumerate_patch_seeds(graph, eligible, config.PATCH_GEODESIC_RADIUS_A)
    assignments = _assign_to_patches(graph, eligible, seeds)
    assignments = _merge_small_patches(assignments, by_pos)

    priors: dict[int, float] = {}
    for seed, members in assignments.items():
        items = [by_pos[position] for position in members]
        prior, _ = _patch_prior(items)
        priors[seed] = prior

    patch_limit = max(1, maximum_regions)
    assignments, merged_targets = _merge_down_to_limit(assignments, priors, by_pos, patch_limit)

    # Membrane orientation, unchanged from the previous partitioning.  Computed
    # BEFORE tier assignment so the cell-surface accessibility penalty can
    # actually move a membrane-proximal patch out of the HIGH tier.
    anchor_record = min(
        eligible,
        key=lambda item: abs(item.uniprot_pos - protein.membrane_anchor_position),
    )
    anchor = anchor_record.ca
    membrane_distances = {
        record.uniprot_pos: _distance(record.ca, anchor)  # type: ignore[arg-type]
        for record in eligible
    }
    distance_values = np.asarray(list(membrane_distances.values()), dtype=float)
    lower = float(np.quantile(distance_values, 0.33))
    upper = float(np.quantile(distance_values, 0.67))

    membrane_positions: dict[int, str] = {}
    for seed, members in assignments.items():
        patch_mean = mean(membrane_distances[position] for position in members)
        if patch_mean <= lower:
            membrane_positions[seed] = "membrane-proximal"
        elif patch_mean >= upper:
            membrane_positions[seed] = "membrane-distal"
        else:
            membrane_positions[seed] = "middle ECD"

    if assay_format == "cell_surface":
        for seed in assignments:
            if membrane_positions[seed] == "membrane-proximal":
                priors[seed] = round(priors[seed] * config.MEMBRANE_PROXIMAL_PRIOR_SCALE, 4)

    tiers = _assign_tiers({seed: priors[seed] for seed in assignments})

    ordered = sorted(
        assignments.items(),
        key=lambda pair: (mean(membrane_distances[position] for position in pair[1]), min(pair[1])),
    )
    regions: list[RegionRecord] = []
    coverage_plan: list[dict[str, Any]] = []
    for index, (seed, members) in enumerate(ordered, start=1):
        members = sorted(set(members))
        items = [by_pos[position] for position in members]
        domain = Counter(item.domain for item in items).most_common(1)[0][0]
        mean_distance = mean(membrane_distances[position] for position in members)
        membrane_position = membrane_positions[seed]

        prior, components = _patch_prior(items)
        assay_accessible = not (
            assay_format == "cell_surface" and membrane_position == "membrane-proximal"
        )
        if not assay_accessible:
            prior = round(prior * config.MEMBRANE_PROXIMAL_PRIOR_SCALE, 4)
        tier = tiers.get(seed, "MEDIUM")
        required_probes = 2 if tier == "HIGH" else 1

        interface_fraction = sum(item.interface for item in items) / len(items)
        pae_median = median_region_pae(members, records, confidence) if predicted else None
        auto_allowed = not predicted or pae_median is None or pae_median <= 10.0

        rationale_parts = [
            f"patch scale: antibody-footprint (geodesic radius <= {config.PATCH_GEODESIC_RADIUS_A:g} A)",
            f"domain anchor: {domain}",
            f"orientation: {membrane_position}",
            f"surface character: {_structural_character(items)}",
            f"epitope prior={prior:.3f}",
        ]
        if seed in merged_targets:
            rationale_parts.append(
                "absorbed a weaker neighbour to respect the Stage-1 panel budget"
            )
        if not assay_accessible:
            rationale_parts.append(
                "membrane-proximal for a cell-surface assay; IgG accessibility is limited"
            )
        if predicted and pae_median is not None:
            rationale_parts.append(f"median intra-region PAE={pae_median:.2f} A")

        region = RegionRecord(
            id=index,
            name=f"{domain} / {membrane_position} / patch {index}",
            domain=domain,
            membrane_position=membrane_position,
            structural_character=_structural_character(items),
            residues=members,
            surface_area=round(sum(item.sidechain_sasa for item in items), 2),
            geodesic_diameter=_geodesic_diameter(graph, members),
            interface_fraction=round(interface_fraction, 3),
            pae_median=pae_median,
            auto_design_allowed=auto_allowed,
            rationale="; ".join(rationale_parts),
            color=config.LEGACY_REGION_COLORS[(index - 1) % len(config.LEGACY_REGION_COLORS)],
            domain_id=domain,
            epitope_prior=prior,
            prior_components=components,
            tier=tier,
            assay_accessible=assay_accessible,
            required_probe_coverage=required_probes,
        )
        regions.append(region)
        for position in members:
            by_pos[position].region_id = index
        coverage_plan.append(
            {
                "patch_id": index,
                "tier": tier,
                "epitope_prior": prior,
                "required_probes": required_probes,
                "residues": len(members),
                "assay_accessible": assay_accessible,
            }
        )

    assigned = sum(len(region.residues) for region in regions)
    coverage = assigned / len(eligible)
    warnings: list[str] = []
    if len(regions) < minimum_regions:
        warnings.append(
            f"Only {len(regions)} structurally defensible patches were produced; the ECD surface did not support the requested minimum of {minimum_regions}."
        )
    if coverage < 0.9:
        warnings.append(f"Eligible-surface assignment coverage is only {coverage:.1%}.")
    if not any(item.conservation_score is not None for item in eligible):
        warnings.append(
            "Sequence conservation was not computed; epitope priors used a neutral variability term."
        )
    if excluded_tail:
        warnings.append(
            f"{excluded_tail} flexible terminal residues (no observed secondary structure) were excluded from patch design; they remain candidates for a linear alanine scan."
        )
    return regions, {
        "eligible_residues": len(eligible),
        "assigned_residues": assigned,
        "coverage": round(coverage, 4),
        "patch_radius_a": config.PATCH_GEODESIC_RADIUS_A,
        "actual_regions": len(regions),
        "assay_format": assay_format,
        "excluded_flexible_terminal_residues": excluded_tail,
        "coverage_plan": coverage_plan,
        "method": "Fab-footprint patch partitioning with epitope priors (see docs/改进方案.md)",
        "warnings": warnings,
    }
