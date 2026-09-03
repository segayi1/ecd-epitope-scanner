from __future__ import annotations

import math
from collections import Counter, defaultdict
from statistics import mean
from typing import Any

import networkx as nx
import numpy as np

from .models import ProteinContext, RegionRecord, ResidueRecord
from .structure_analysis import median_region_pae


REGION_COLORS = [
    "#087E8B",
    "#E4572E",
    "#6A4C93",
    "#2E7D32",
    "#C47F00",
    "#1565C0",
    "#AD1457",
    "#5D4037",
]


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))


def _target_region_count(count: int, domains: int, minimum: int, maximum: int) -> int:
    if count < 8:
        return 1
    if count < 18:
        return min(3, max(2, count // 5))
    if count < 30:
        return min(4, max(3, count // 7))
    target = max(minimum, round(count / 18))
    target = max(target, min(domains, maximum))
    return min(maximum, target)


def _build_graph(records: list[ResidueRecord]) -> nx.Graph:
    graph = nx.Graph()
    for record in records:
        graph.add_node(record.uniprot_pos)
    for index, left in enumerate(records):
        assert left.ca is not None
        for right in records[index + 1 :]:
            assert right.ca is not None
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

    # Surface islands can arise around flexible loops. Connect each island to its nearest
    # neighbor so every eligible residue is auditable and assigned exactly once.
    components = list(nx.connected_components(graph))
    while len(components) > 1:
        best: tuple[float, int, int] | None = None
        first = components[0]
        rest = set().union(*components[1:])
        by_pos = {record.uniprot_pos: record for record in records}
        for left_pos in first:
            for right_pos in rest:
                distance = _distance(by_pos[left_pos].ca, by_pos[right_pos].ca)  # type: ignore[arg-type]
                candidate = (distance, left_pos, right_pos)
                if best is None or candidate < best:
                    best = candidate
        assert best is not None
        distance, left_pos, right_pos = best
        domain_penalty = 4.0 if by_pos[left_pos].domain != by_pos[right_pos].domain else 0.0
        graph.add_edge(left_pos, right_pos, weight=distance + domain_penalty + 6.0, euclidean=distance)
        components = list(nx.connected_components(graph))
    return graph


def _allocate_quotas(records: list[ResidueRecord], target: int) -> dict[str, int]:
    counts = Counter(record.domain for record in records)
    domains = sorted(counts, key=lambda name: (-counts[name], name))
    quotas = {domain: 1 for domain in domains[:target]}
    remaining = target - len(quotas)
    if remaining <= 0:
        return quotas
    total = sum(counts[domain] for domain in quotas)
    fractional: list[tuple[float, str]] = []
    for domain in quotas:
        exact = remaining * counts[domain] / max(1, total)
        whole = int(math.floor(exact))
        quotas[domain] += whole
        fractional.append((exact - whole, domain))
    assigned = sum(quotas.values())
    for _, domain in sorted(fractional, key=lambda item: (-item[0], item[1])):
        if assigned >= target:
            break
        quotas[domain] += 1
        assigned += 1
    return quotas


def _select_seeds(
    graph: nx.Graph, records: list[ResidueRecord], quotas: dict[str, int]
) -> list[int]:
    by_pos = {record.uniprot_pos: record for record in records}
    seeds: list[int] = []
    for domain in sorted(quotas):
        nodes = sorted(record.uniprot_pos for record in records if record.domain == domain)
        if not nodes:
            continue
        local = graph.subgraph(nodes)
        first = max(
            nodes,
            key=lambda position: (
                by_pos[position].rsasa
                + min(1.0, by_pos[position].sidechain_sasa / 120.0)
                + 1.0 / (1 + by_pos[position].contacts),
                -position,
            ),
        )
        local_seeds = [first]
        while len(local_seeds) < min(quotas[domain], len(nodes)):
            distances: dict[int, float] = {}
            for node in nodes:
                if node in local_seeds:
                    continue
                values = []
                for seed in local_seeds:
                    try:
                        values.append(nx.shortest_path_length(local, node, seed, weight="weight"))
                    except nx.NetworkXNoPath:
                        values.append(float("inf"))
                distances[node] = min(values)
            if not distances:
                break
            next_seed = max(distances, key=lambda position: (distances[position], -position))
            local_seeds.append(next_seed)
        seeds.extend(local_seeds)
    return seeds


def _assign_to_seeds(
    graph: nx.Graph, records: list[ResidueRecord], seeds: list[int]
) -> dict[int, list[int]]:
    by_pos = {record.uniprot_pos: record for record in records}
    assignments: dict[int, list[int]] = {seed: [] for seed in seeds}
    for record in records:
        same_domain = [seed for seed in seeds if by_pos[seed].domain == record.domain]
        candidates = same_domain or seeds
        ranked: list[tuple[float, int]] = []
        for seed in candidates:
            try:
                distance = float(nx.shortest_path_length(graph, record.uniprot_pos, seed, weight="weight"))
            except nx.NetworkXNoPath:
                distance = _distance(record.ca, by_pos[seed].ca) + 100.0  # type: ignore[arg-type]
            ranked.append((distance, seed))
        _, selected = min(ranked, key=lambda item: (item[0], item[1]))
        assignments[selected].append(record.uniprot_pos)

    # Merge regions too small to support a meaningful multi-site perturbation.
    for seed in list(assignments):
        members = assignments.get(seed, [])
        if len(members) >= 3 or len(assignments) == 1:
            continue
        alternatives = [item for item in assignments if item != seed]
        if not alternatives:
            continue
        replacement = min(
            alternatives,
            key=lambda other: (
                _distance(by_pos[seed].ca, by_pos[other].ca),  # type: ignore[arg-type]
                other,
            ),
        )
        assignments[replacement].extend(members)
        del assignments[seed]
    return assignments


def _geodesic_diameter(graph: nx.Graph, members: list[int]) -> float:
    subgraph = graph.subgraph(members)
    maximum = 0.0
    for source, lengths in nx.all_pairs_dijkstra_path_length(subgraph, weight="weight"):
        del source
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
) -> tuple[list[RegionRecord], dict[str, Any]]:
    eligible = [record for record in records if record.folding_screen_pass and record.ca is not None]
    if not eligible:
        return [], {
            "eligible_residues": 0,
            "assigned_residues": 0,
            "coverage": 0.0,
            "warnings": ["No residues passed the configured surface and folding-risk screen."],
            "method": "domain-constrained geodesic surface partitioning",
        }
    graph = _build_graph(eligible)
    domains = len({record.domain for record in eligible})
    target = _target_region_count(len(eligible), domains, minimum_regions, maximum_regions)
    quotas = _allocate_quotas(eligible, target)
    seeds = _select_seeds(graph, eligible, quotas)
    assignments = _assign_to_seeds(graph, eligible, seeds)
    by_pos = {record.uniprot_pos: record for record in eligible}

    anchor_record = min(
        eligible,
        key=lambda item: abs(item.uniprot_pos - protein.membrane_anchor_position),
    )
    anchor = anchor_record.ca
    assert anchor is not None
    membrane_distances = {
        record.uniprot_pos: _distance(record.ca, anchor)  # type: ignore[arg-type]
        for record in eligible
    }
    distance_values = np.asarray(list(membrane_distances.values()), dtype=float)
    lower = float(np.quantile(distance_values, 0.33))
    upper = float(np.quantile(distance_values, 0.67))

    ordered = sorted(
        assignments.items(),
        key=lambda pair: (
            mean(membrane_distances[position] for position in pair[1]),
            min(pair[1]),
        ),
    )
    regions: list[RegionRecord] = []
    for index, (_, members) in enumerate(ordered, start=1):
        members = sorted(set(members))
        items = [by_pos[position] for position in members]
        domain = Counter(item.domain for item in items).most_common(1)[0][0]
        mean_distance = mean(membrane_distances[position] for position in members)
        if mean_distance <= lower:
            membrane_position = "membrane-proximal"
        elif mean_distance >= upper:
            membrane_position = "membrane-distal"
        else:
            membrane_position = "middle ECD"
        character = _structural_character(items)
        interface_fraction = sum(item.interface for item in items) / len(items)
        pae_median = median_region_pae(members, records, confidence) if predicted else None
        auto_allowed = not predicted or pae_median is None or pae_median <= 10.0
        name = f"{domain} / {membrane_position} / surface {index}"
        rationale_parts = [
            f"domain anchor: {domain}",
            f"orientation: {membrane_position}",
            f"surface character: {character}",
            "boundary follows geodesic distance, low-accessibility bottlenecks and domain/secondary-structure penalties",
        ]
        if interface_fraction:
            rationale_parts.append(f"native-interface fraction={interface_fraction:.2f}")
        if predicted and pae_median is not None:
            rationale_parts.append(f"median intra-region PAE={pae_median:.2f} A")
        region = RegionRecord(
            id=index,
            name=name,
            domain=domain,
            membrane_position=membrane_position,
            structural_character=character,
            residues=members,
            surface_area=round(sum(item.sidechain_sasa for item in items), 2),
            geodesic_diameter=_geodesic_diameter(graph, members),
            interface_fraction=round(interface_fraction, 3),
            pae_median=pae_median,
            auto_design_allowed=auto_allowed,
            rationale="; ".join(rationale_parts),
            color=REGION_COLORS[(index - 1) % len(REGION_COLORS)],
        )
        regions.append(region)
        for position in members:
            by_pos[position].region_id = index

    assigned = sum(len(region.residues) for region in regions)
    coverage = assigned / len(eligible)
    warnings: list[str] = []
    if len(regions) < minimum_regions:
        warnings.append(
            f"Only {len(regions)} structurally defensible regions were produced; the ECD surface did not support the requested minimum of {minimum_regions}."
        )
    if coverage < 0.9:
        warnings.append(f"Eligible-surface assignment coverage is only {coverage:.1%}.")
    return regions, {
        "eligible_residues": len(eligible),
        "assigned_residues": assigned,
        "coverage": round(coverage, 4),
        "target_regions": target,
        "actual_regions": len(regions),
        "method": "domain- and membrane-orientation-constrained geodesic surface partitioning",
        "warnings": warnings,
    }
