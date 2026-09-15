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
            # A C-alpha contact alone joins the two sheets of an IgV sandwich.
            # Preserve local backbone continuity, but make non-local spatial
            # edges surface-facing: opposing side-chain vectors cannot define
            # one antibody-accessible face.
            same_face = False
            if left.sidechain_direction is not None and right.sidechain_direction is not None:
                same_face = float(np.dot(left.sidechain_direction, right.sidechain_direction)) > config.PATCH_OPPOSING_FACE_COSINE_MAX
            if (sequence_neighbor and distance <= 16.0) or (distance <= 12.5 and same_face):
                penalty = 4.0 if left.domain != right.domain else 0.0
                secondary_penalty = 1.5 if left.secondary_structure != right.secondary_structure else 0.0
                graph.add_edge(
                    left.uniprot_pos,
                    right.uniprot_pos,
                    weight=distance + penalty + secondary_penalty,
                    euclidean=distance,
                )

    # Deliberately retain disconnected surface islands.  A forced nearest-Cα
    # bridge can tunnel through a beta-sandwich core and merge opposite faces.
    # Seed enumeration handles disconnected components deterministically.
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
    seed_distances = {first: nx.single_source_dijkstra_path_length(graph, first, weight="weight")}

    def effective_distances() -> dict[int, float]:
        output: dict[int, float] = {}
        for node in nodes:
            if node in seeds:
                continue
            reachable = [
                seed_distances[seed].get(node, math.inf)
                for seed in seeds
                if _distance(by_pos[node].ca, by_pos[seed].ca) <= config.PATCH_EUCLIDEAN_CAP_A  # type: ignore[arg-type]
                and seed_distances[seed].get(node, math.inf) <= radius
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
        seed_distances[next_seed] = nx.single_source_dijkstra_path_length(graph, next_seed, weight="weight")
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
            # Seed enumeration must have created a local seed.  Keeping an
            # unassigned residue is safer than silently stretching a patch
            # across an inaccessible surface route.
            assignments.setdefault(record.uniprot_pos, []).append(record.uniprot_pos)
            continue
        if math.isinf(best_distance):
            assignments.setdefault(record.uniprot_pos, []).append(record.uniprot_pos)
            continue
        assignments[best_seed].append(record.uniprot_pos)
    return assignments


def _patch_is_compact(
    graph: nx.Graph, members: list[int], by_pos: dict[int, ResidueRecord]
) -> bool:
    del graph
    for index, left in enumerate(members):
        for right in members[index + 1 :]:
            if _distance(by_pos[left].ca, by_pos[right].ca) > config.PATCH_MAX_EUCLIDEAN_DIAMETER_A:  # type: ignore[arg-type]
                return False
    return True


def _split_noncompact_patches(
    assignments: dict[int, list[int]], graph: nx.Graph, by_pos: dict[int, ResidueRecord]
) -> dict[int, list[int]]:
    """Split, never stretch, a patch that violates a footprint hard limit."""
    output: dict[int, list[int]] = {}
    for seed, members in assignments.items():
        pending = sorted(set(members))
        while pending:
            local_seed = pending.pop(0)
            cluster = [local_seed]
            for position in list(pending):
                if _patch_is_compact(graph, cluster + [position], by_pos):
                    cluster.append(position)
                    pending.remove(position)
            output[local_seed] = cluster
    return output


def _merge_small_patches(
    assignments: dict[int, list[int]], graph: nx.Graph, by_pos: dict[int, ResidueRecord]
) -> dict[int, list[int]]:
    for seed in list(assignments):
        members = assignments.get(seed, [])
        if len(members) >= config.PATCH_MIN_RESIDUES or len(assignments) == 1:
            continue
        alternatives = [item for item in assignments if item != seed]
        if not alternatives:
            continue
        legal = [other for other in alternatives if _patch_is_compact(graph, assignments[other] + members, by_pos)]
        if not legal:
            continue
        replacement = min(legal, key=lambda other: (_distance(by_pos[seed].ca, by_pos[other].ca), other))  # type: ignore[arg-type]
        assignments[replacement].extend(members)
        del assignments[seed]
    return assignments


def _merge_down_to_limit(
    assignments: dict[int, list[int]],
    priors: dict[int, float],
    graph: nx.Graph,
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
        legal = [other for other in alternatives if _patch_is_compact(graph, assignments[other] + assignments[weakest], by_pos)]
        if not legal:
            break
        replacement = min(legal, key=lambda other: (_distance(by_pos[weakest].ca, by_pos[other].ca), other))  # type: ignore[arg-type]
        priors[replacement] = min(priors.get(replacement, 0.0), priors.get(weakest, 0.0))
        assignments[replacement].extend(assignments[weakest])
        del assignments[weakest]
        merged_targets.append(replacement)
    return assignments, merged_targets


def _patch_prior(items: list[ResidueRecord]) -> tuple[float, dict[str, float]]:
    exposure = mean(min(1.0, item.rsasa / 0.5) for item in items)
    protrusion = mean(1.0 / (1.0 + item.contacts / 6.0) for item in items)
    sse_available = any(item.secondary_structure in {"H", "E", "C"} for item in items)
    # AF3 mmCIF commonly contains no HELIX/SHEET assignment.  Unknown must
    # not masquerade as a loop-rich epitope or receive a loop bonus.
    loop = (
        sum(item.secondary_structure == "C" for item in items) / len(items)
        if sse_available else 1.0 / 3.0
    )
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
        "secondary_structure_available": 1.0 if sse_available else 0.0,
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
    if not any(item.secondary_structure in {"H", "E", "C"} for item in items):
        return "secondary structure unavailable"
    counts = Counter(item.secondary_structure for item in items)
    total = max(1, len(items))
    if counts["C"] / total >= 0.55:
        return "loop-rich surface"
    if counts["E"] / total >= 0.45:
        return "beta-sheet face"
    if counts["H"] / total >= 0.45:
        return "helical face"
    return "mixed secondary-structure face"


def _membrane_accessibility(
    records: list[ResidueRecord],
    anchor: tuple[float, float, float] | None,
    normal: tuple[float, float, float] | None,
) -> tuple[dict[int, float], dict[int, float], str]:
    """Estimate continuous cell-surface reachability from an ECD-facing plane.

    Experimental structures frequently omit the transmembrane helix.  We use
    the vector from the ECD anchor to the observed ECD centroid as a
    reproducible surrogate for the membrane normal, and retain an explicit
    fallback label when this geometry cannot be estimated.  This is a prior,
    not a claim that a membrane-proximal epitope is inaccessible.
    """
    if anchor is None or normal is None:
        return ({item.uniprot_pos: 1.0 for item in records}, {}, "unavailable")
    anchor_point = np.asarray(anchor, dtype=float)
    normal_vector = np.asarray(normal, dtype=float)
    norm = float(np.linalg.norm(normal_vector))
    if norm < 1e-6:
        return ({item.uniprot_pos: 1.0 for item in records}, {}, "unavailable")
    normal_vector /= norm
    plane_distances: dict[int, float] = {}
    scores: dict[int, float] = {}
    for item in records:
        if item.ca is None:
            continue
        distance = abs(float(np.dot(np.asarray(item.ca, dtype=float) - anchor_point, normal_vector)))
        plane_distances[item.uniprot_pos] = distance
        scaled = min(1.0, distance / config.MEMBRANE_ACCESSIBILITY_SCALE_A)
        score = config.MEMBRANE_ACCESSIBILITY_MIN_FACTOR + (
            1.0 - config.MEMBRANE_ACCESSIBILITY_MIN_FACTOR
        ) * scaled
        if item.glycan_shielded:
            score *= 1.0 - config.MEMBRANE_ACCESSIBILITY_GLYCAN_PENALTY
        scores[item.uniprot_pos] = round(float(np.clip(score, 0.0, 1.0)), 4)
    return scores, plane_distances, "transmembrane-helix-derived membrane-plane"


def _partition_surface_footprint(
    records: list[ResidueRecord],
    protein: ProteinContext,
    confidence: dict[str, Any] | None = None,
    minimum_regions: int = 5,
    maximum_regions: int = 8,
    predicted: bool = False,
    assay_format: str = "cell_surface",
    membrane_anchor: tuple[float, float, float] | None = None,
    membrane_normal: tuple[float, float, float] | None = None,
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
    assignments = _split_noncompact_patches(assignments, graph, by_pos)
    assignments = _merge_small_patches(assignments, graph, by_pos)

    priors: dict[int, float] = {}
    for seed, members in assignments.items():
        items = [by_pos[position] for position in members]
        prior, _ = _patch_prior(items)
        priors[seed] = prior

    patch_limit = max(1, maximum_regions)
    assignments, merged_targets = _merge_down_to_limit(assignments, priors, graph, by_pos, patch_limit)

    # Estimate a continuous cell-surface accessibility prior before tiering.
    # It is based on a membrane-plane surrogate rather than a binary
    # membrane-proximal exclusion.
    membrane_scores, membrane_distances, membrane_geometry = _membrane_accessibility(
        eligible, membrane_anchor, membrane_normal
    )
    distance_values = np.asarray(list(membrane_distances.values()), dtype=float)
    lower = float(np.quantile(distance_values, 0.33)) if len(distance_values) else 0.0
    upper = float(np.quantile(distance_values, 0.67)) if len(distance_values) else 0.0

    membrane_positions: dict[int, str] = {}
    for seed, members in assignments.items():
        patch_mean = mean(membrane_distances[position] for position in members) if membrane_distances else 0.0
        if not membrane_distances:
            membrane_positions[seed] = "orientation unavailable"
        elif patch_mean <= lower:
            membrane_positions[seed] = "membrane-proximal"
        elif patch_mean >= upper:
            membrane_positions[seed] = "membrane-distal"
        else:
            membrane_positions[seed] = "middle ECD"

    if assay_format == "cell_surface":
        for seed, members in assignments.items():
            accessibility = mean(membrane_scores[position] for position in members)
            priors[seed] = round(priors[seed] * accessibility, 4)

    tiers = _assign_tiers({seed: priors[seed] for seed in assignments})

    ordered = sorted(
        assignments.items(),
        key=lambda pair: (
            mean(membrane_distances[position] for position in pair[1])
            if membrane_distances else 0.0,
            min(pair[1]),
        ),
    )
    regions: list[RegionRecord] = []
    coverage_plan: list[dict[str, Any]] = []
    for index, (seed, members) in enumerate(ordered, start=1):
        members = sorted(set(members))
        items = [by_pos[position] for position in members]
        domain = Counter(item.domain for item in items).most_common(1)[0][0]
        membrane_position = membrane_positions[seed]

        prior, components = _patch_prior(items)
        membrane_accessibility = mean(membrane_scores[position] for position in members)
        predicted_accessibility = (
            "unknown" if not membrane_distances else
            "reduced" if assay_format == "cell_surface" and membrane_accessibility < config.MEMBRANE_ACCESSIBILITY_REDUCED_THRESHOLD
            else "normal"
        )
        assay_accessible = predicted_accessibility != "reduced"
        if assay_format == "cell_surface":
            prior = round(prior * membrane_accessibility, 4)
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
        if predicted_accessibility == "reduced":
            rationale_parts.append(
                f"predicted cell-surface accessibility reduced ({membrane_accessibility:.2f}); retain as a lower-priority coverage probe"
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
            predicted_accessibility=predicted_accessibility,
            membrane_accessibility_score=round(membrane_accessibility, 4),
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
                "predicted_accessibility": predicted_accessibility,
                "membrane_accessibility_score": round(membrane_accessibility, 4),
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
        "membrane_accessibility_model": membrane_geometry,
        "membrane_accessibility_scale_A": config.MEMBRANE_ACCESSIBILITY_SCALE_A,
        "excluded_flexible_terminal_residues": excluded_tail,
        "coverage_plan": coverage_plan,
        "method": "Fab-footprint patch partitioning with epitope priors (see docs/改进方案.md)",
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Compact single-domain rescue mode
# ---------------------------------------------------------------------------

def _compact_envelope(records: list[ResidueRecord]) -> list[ResidueRecord]:
    """Return a *geometric* surface envelope, not a mutation menu.

    Protected residues are intentionally retained: deleting a glycan, a
    disulfide or a high-risk site before geometry is built turns one physical
    face into artificial islands.  Low-confidence predicted coordinates and
    flexible termini remain outside the automatic 3-D envelope.
    """
    return [
        item for item in records
        if item.rsasa >= config.COMPACT_ENVELOPE_RSASA and item.ca is not None
        and item.surface_normal is not None
        and not (config.EXCLUDE_FLEXIBLE_TERMINAL and item.flexible_terminal)
        and (item.plddt is None or item.plddt >= config.COMPACT_MIN_CORE_CONFIDENCE)
    ]


def _outward_normal(record: ResidueRecord) -> np.ndarray:
    """Read the solvent-derived local normal; V2 has no radial fallback."""
    normal = np.asarray(record.surface_normal, dtype=float)
    return normal / max(float(np.linalg.norm(normal)), 1e-6)


def _compact_graph(
    envelope: list[ResidueRecord], normals: dict[int, np.ndarray]
) -> nx.Graph:
    graph = nx.Graph()
    for item in envelope:
        graph.add_node(item.uniprot_pos)
    for index, left in enumerate(envelope):
        for right in envelope[index + 1 :]:
            distance = _distance(left.ca, right.ca)  # type: ignore[arg-type]
            alignment = float(np.dot(normals[left.uniprot_pos], normals[right.uniprot_pos]))
            if distance <= config.COMPACT_GRAPH_CONTACT_A and alignment >= config.COMPACT_MIN_NORMAL_COSINE:
                graph.add_edge(
                    left.uniprot_pos,
                    right.uniprot_pos,
                    weight=distance * (1.0 + 0.35 * (1.0 - alignment)),
                )
    return graph


def _dominant_domain_fraction(envelope: list[ResidueRecord]) -> tuple[str, float]:
    labels = [item.domain for item in envelope if item.domain and item.domain != "ECD"]
    if not labels:
        return "ECD", 0.0
    domain, count = Counter(labels).most_common(1)[0]
    return domain, count / len(envelope)


def _beta_sheet_assignments(
    records: list[ResidueRecord], protein: ProteinContext, envelope: list[ResidueRecord]
) -> dict[int, list[int]] | None:
    """Group annotated beta strands by sheet-like backbone contacts.

    Reviewed UniProt beta-strand annotations are used only as segment bounds;
    sheet membership is derived from the supplied 3-D coordinates.  Adjacent
    sequence segments are excluded from the contact graph because a tight
    beta turn joins different sheets spatially without making them one sheet.
    """
    by_pos = {item.uniprot_pos: item for item in records if item.ca is not None}
    segments: list[list[int]] = []
    for feature in protein.features:
        if feature.type.lower() != "beta strand":
            continue
        positions = [position for position in range(feature.start, feature.end + 1) if position in by_pos]
        if len(positions) >= 2:
            segments.append(positions)
    if len(segments) < 2 * config.COMPACT_BETA_SHEET_MIN_STRANDS:
        return None
    strand_graph = nx.Graph()
    strand_graph.add_nodes_from(range(len(segments)))
    for left_index, left in enumerate(segments):
        for right_index, right in enumerate(segments[left_index + 1 :], start=left_index + 1):
            sequence_gap = min(abs(a - b) for a in left for b in right)
            if sequence_gap <= 5:
                continue
            contacts = sum(
                _distance(by_pos[a].ca, by_pos[b].ca) <= config.COMPACT_BETA_STRAND_CONTACT_A  # type: ignore[arg-type]
                for a in left for b in right
            )
            if contacts >= config.COMPACT_BETA_STRAND_MIN_CONTACTS:
                strand_graph.add_edge(left_index, right_index)
    components = [
        component for component in nx.connected_components(strand_graph)
        if len(component) >= config.COMPACT_BETA_SHEET_MIN_STRANDS
    ]
    if not 1 <= len(components) <= config.COMPACT_MAX_SURFACE_REGIONS:
        return None
    annotated_positions = {position for segment in segments for position in segment}
    covered_positions = {
        position for component in components for index in component for position in segments[index]
    }
    if len(covered_positions) / max(1, len(annotated_positions)) < config.COMPACT_BETA_SHEET_MIN_COVERAGE:
        return None
    scaffolds = [
        sorted({position for index in component for position in segments[index]})
        for component in components
    ]
    assignments: dict[int, list[int]] = {min(scaffold): [] for scaffold in scaffolds}
    for item in envelope:
        seed, _ = min(
            (
                (min(scaffold), min(_distance(item.ca, by_pos[position].ca) for position in scaffold))  # type: ignore[arg-type]
                for scaffold in scaffolds
            ),
            key=lambda pair: (pair[1], pair[0]),
        )
        assignments[seed].append(item.uniprot_pos)
    return {seed: sorted(members) for seed, members in assignments.items() if members}


def _add_within_face_geometry_edges(
    graph: nx.Graph, assignments: dict[int, list[int]], by_pos: dict[int, ResidueRecord]
) -> nx.Graph:
    """Fill local surface gaps only after opposite beta sheets are separated."""
    output = graph.copy()
    for members in assignments.values():
        for index, left in enumerate(members):
            for right in members[index + 1 :]:
                distance = _distance(by_pos[left].ca, by_pos[right].ca)  # type: ignore[arg-type]
                if distance <= config.COMPACT_GRAPH_CONTACT_A:
                    output.add_edge(left, right, weight=distance)
    return output


def _compact_trigger(
    records: list[ResidueRecord], protein: ProteinContext, legacy: list[RegionRecord]
) -> tuple[bool, list[str], list[ResidueRecord]]:
    envelope = _compact_envelope(records)
    domain, fraction = _dominant_domain_fraction(envelope)
    if (
        protein.ecd_length > config.COMPACT_ECD_MAX_LENGTH_AA
        or len(envelope) < config.COMPACT_MIN_ENVELOPE_RESIDUES
        or fraction < config.COMPACT_DOMINANT_DOMAIN_MIN_FRACTION
    ):
        return False, [], envelope
    mutable = sum(item.folding_screen_pass for item in envelope)
    ratio = mutable / len(envelope)
    reasons: list[str] = []
    if len(legacy) > config.COMPACT_MAX_SURFACE_REGIONS:
        reasons.append(f"legacy produced {len(legacy)} regions for compact domain {domain}")
    if any(region.geodesic_diameter > 2.0 * config.PATCH_GEODESIC_RADIUS_A for region in legacy):
        reasons.append("legacy region exceeds two footprint radii along the surface")
    if ratio < config.COMPACT_MUTABLE_ENVELOPE_RATIO_TRIGGER:
        reasons.append(f"mutable/envelope ratio {ratio:.2f} indicates sparse-site fragmentation")
    # A compact protein is re-routed only when its legacy result has a clear
    # pathology, never merely because it is short.
    return bool(reasons), reasons, envelope


def _compact_seed_positions(
    items: list[ResidueRecord], normals: dict[int, np.ndarray], count: int
) -> list[int]:
    first = max(items, key=lambda item: (_exposure_priority(item), -item.uniprot_pos))
    seeds = [first]
    while len(seeds) < count:
        candidate = max(
            (item for item in items if item.uniprot_pos not in {seed.uniprot_pos for seed in seeds}),
            key=lambda item: min(
                1.0 - float(np.dot(normals[item.uniprot_pos], normals[seed.uniprot_pos]))
                for seed in seeds
            ),
        )
        seeds.append(candidate)
    return [item.uniprot_pos for item in seeds]


def _compact_assignments(
    graph: nx.Graph, envelope: list[ResidueRecord], normals: dict[int, np.ndarray], seed_count: int
) -> dict[int, list[int]] | None:
    components = list(nx.connected_components(graph))
    if len(components) > config.COMPACT_MAX_SURFACE_REGIONS:
        return None
    seed_count = max(seed_count, len(components))
    seed_count = min(seed_count, config.COMPACT_MAX_SURFACE_REGIONS, len(envelope))
    # Guarantee that every pre-existing surface component gets a seed before
    # angular farthest-point seeding adds optional subdivisions.
    by_pos = {item.uniprot_pos: item for item in envelope}
    seeds = [
        max(component, key=lambda position: (_exposure_priority(by_pos[position]), -position))
        for component in sorted(components, key=lambda group: min(group))
    ]
    if len(seeds) < seed_count:
        extra_items = [item for item in envelope if item.uniprot_pos not in seeds]
        while len(seeds) < seed_count and extra_items:
            chosen = max(
                extra_items,
                key=lambda item: min(
                    1.0 - float(np.dot(normals[item.uniprot_pos], normals[seed]))
                    for seed in seeds
                ),
            )
            seeds.append(chosen.uniprot_pos)
            extra_items.remove(chosen)
    lengths = {seed: nx.single_source_dijkstra_path_length(graph, seed, weight="weight") for seed in seeds}
    assignments: dict[int, list[int]] = {seed: [] for seed in seeds}
    for item in envelope:
        choices = []
        for seed in seeds:
            route = lengths[seed].get(item.uniprot_pos, math.inf)
            angular = 1.0 - float(np.dot(normals[item.uniprot_pos], normals[seed]))
            # Surface direction defines the face; geodesic distance only
            # regularises its boundary.  This prevents left/middle/right
            # spatial slices through one flat beta sheet.
            choices.append((12.0 * angular + 0.15 * route, seed))
        _, selected = min(choices, key=lambda pair: (pair[0], pair[1]))
        assignments[selected].append(item.uniprot_pos)
    assignments = {seed: sorted(members) for seed, members in assignments.items() if members}
    # A surface face must be graph-connected.  A disconnected result is not
    # "repaired" by a nearest-C-alpha bridge; the caller safely falls back.
    if any(len(list(nx.connected_components(graph.subgraph(members)))) != 1 for members in assignments.values()):
        return None
    for members in assignments.values():
        mean_normal = np.mean([normals[position] for position in members], axis=0)
        mean_normal /= max(float(np.linalg.norm(mean_normal)), 1e-6)
        coherence = mean(float(np.dot(normals[position], mean_normal)) for position in members)
        if coherence < config.COMPACT_SURFACE_NORMAL_MIN_COHERENCE:
            return None
    return assignments


def _compact_footprints(
    members: list[int], graph: nx.Graph, by_pos: dict[int, ResidueRecord]
) -> list[tuple[int, list[int]]]:
    mutable = [position for position in members if by_pos[position].folding_screen_pass]
    if not mutable:
        return []
    ranked = sorted(mutable, key=lambda position: (-_exposure_priority(by_pos[position]), position))
    footprints: list[tuple[int, list[int]]] = []
    for center in ranked:
        paths = nx.single_source_dijkstra_path_length(graph.subgraph(members), center, weight="weight")
        footprint = sorted(
            position for position in mutable
            if paths.get(position, math.inf) <= config.COMPACT_FOOTPRINT_GEODESIC_RADIUS_A
            and _distance(by_pos[position].ca, by_pos[center].ca) <= config.COMPACT_FOOTPRINT_EUCLIDEAN_CAP_A  # type: ignore[arg-type]
        )
        if not footprint:
            continue
        candidate = set(footprint)
        if any(len(candidate & set(existing)) / len(candidate | set(existing)) >= 0.70 for _, existing in footprints):
            continue
        footprints.append((center, footprint))
        if len(footprints) >= config.COMPACT_MAX_FOOTPRINTS_PER_FACE:
            break
    return footprints


def _partition_compact_surface_face(
    records: list[ResidueRecord], protein: ProteinContext, confidence: dict[str, Any] | None,
    predicted: bool, assay_format: str, membrane_anchor: tuple[float, float, float] | None,
    membrane_normal: tuple[float, float, float] | None, envelope: list[ResidueRecord], trigger_reasons: list[str],
) -> tuple[list[RegionRecord], dict[str, Any]] | None:
    normals = {item.uniprot_pos: _outward_normal(item) for item in envelope}
    graph = _compact_graph(envelope, normals)
    assignments = _beta_sheet_assignments(records, protein, envelope)
    partition_basis = "beta_sheet_backbone_topology" if assignments is not None else "solvent_surface_normals"
    selected_face_count = None
    if assignments is None:
        minimum_alignment = min(
            (float(np.dot(normals[left], normals[right])) for i, left in enumerate(normals) for right in list(normals)[i + 1 :]),
            default=1.0,
        )
        first_face_count = 2 if minimum_alignment < -0.35 else 1
        for requested_faces in range(first_face_count, config.COMPACT_MAX_SURFACE_REGIONS + 1):
            assignments = _compact_assignments(graph, envelope, normals, requested_faces)
            if assignments is not None:
                selected_face_count = len(assignments)
                break
    else:
        selected_face_count = len(assignments)
    if assignments is None:
        return None
    by_pos = {item.uniprot_pos: item for item in envelope}
    graph = _add_within_face_geometry_edges(graph, assignments, by_pos)
    priors = {seed: _patch_prior([by_pos[position] for position in members])[0] for seed, members in assignments.items()}
    membrane_scores, membrane_distances, membrane_geometry = _membrane_accessibility(envelope, membrane_anchor, membrane_normal)
    if assay_format == "cell_surface":
        for seed, members in assignments.items():
            priors[seed] *= mean(membrane_scores[position] for position in members)
    tiers = _assign_tiers(priors)
    ordered = sorted(assignments.items(), key=lambda pair: (-priors[pair[0]], min(pair[1])))
    regions: list[RegionRecord] = []
    coverage_plan: list[dict[str, Any]] = []
    domain, _ = _dominant_domain_fraction(envelope)
    for index, (seed, members) in enumerate(ordered, start=1):
        items = [by_pos[position] for position in members]
        footprints = _compact_footprints(members, graph, by_pos)
        center, footprint = footprints[0] if footprints else (None, [])
        protected = [position for position in members if not by_pos[position].folding_screen_pass]
        mutable = [position for position in members if by_pos[position].folding_screen_pass]
        prior, components = _patch_prior(items)
        accessibility = mean(membrane_scores[position] for position in members)
        if assay_format == "cell_surface":
            prior *= accessibility
        mean_normal = np.mean([normals[position] for position in members], axis=0)
        mean_normal /= max(float(np.linalg.norm(mean_normal)), 1e-6)
        measured_normal_coherence = mean(float(np.dot(normals[position], mean_normal)) for position in members)
        normal_coherence = (
            None if partition_basis == "beta_sheet_backbone_topology"
            else round(measured_normal_coherence, 4)
        )
        membrane_position = "orientation unavailable"
        if membrane_distances:
            value = mean(membrane_distances[position] for position in members)
            lower, upper = np.quantile(list(membrane_distances.values()), [0.33, 0.67])
            membrane_position = "membrane-proximal" if value <= lower else "membrane-distal" if value >= upper else "middle ECD"
        pae = median_region_pae(footprint, records, confidence) if predicted and footprint else None
        flags = [] if len(footprint) >= 2 else ["insufficient_local_mutable_sites"]
        region = RegionRecord(
            id=index,
            name=f"{domain} / {'beta-sheet zone' if partition_basis == 'beta_sheet_backbone_topology' else 'surface face'} {index}",
            domain=domain,
            membrane_position=membrane_position, structural_character=_structural_character(items),
            residues=footprint, surface_area=round(sum(item.sidechain_sasa for item in items), 2),
            geodesic_diameter=_geodesic_diameter(graph, members),
            interface_fraction=round(sum(item.interface for item in items) / len(items), 3),
            pae_median=pae, auto_design_allowed=(not predicted or pae is None or pae <= 10.0) and len(footprint) >= 2,
            rationale=(f"compact single-domain continuous surface face; envelope={len(members)} residues; "
                       f"mutable={len(mutable)}; local footprint={len(footprint)}; "
                       + (
                           "beta-sheet membership derived from backbone contact topology; "
                           if partition_basis == "beta_sheet_backbone_topology"
                           else f"surface-normal coherence={measured_normal_coherence:.2f}; "
                       )
                       + "; ".join(trigger_reasons)),
            color=config.LEGACY_REGION_COLORS[(index - 1) % len(config.LEGACY_REGION_COLORS)],
            domain_id=domain, epitope_prior=round(prior, 4), prior_components=components,
            tier=tiers[seed], assay_accessible=accessibility >= config.MEMBRANE_ACCESSIBILITY_REDUCED_THRESHOLD,
            predicted_accessibility=("unknown" if not membrane_distances else "reduced" if assay_format == "cell_surface" and accessibility < config.MEMBRANE_ACCESSIBILITY_REDUCED_THRESHOLD else "normal"),
            membrane_accessibility_score=round(accessibility, 4), required_probe_coverage=1,
            partition_mode="compact_surface_face", partition_basis=partition_basis,
            envelope_residues=members,
            mutable_residues=mutable, protected_surface_residues=protected,
            footprint_residues=footprint, footprint_center=center,
            candidate_footprints=[positions for _, positions in footprints],
            surface_connectivity=(
                "beta_sheet_scaffold_connected"
                if partition_basis == "beta_sheet_backbone_topology" else "connected"
            ),
            normal_coherence=normal_coherence, quality_flags=flags,
        )
        regions.append(region)
        for position in members:
            by_pos[position].region_id = index
        coverage_plan.append({"patch_id": index, "tier": tiers[seed], "epitope_prior": round(prior, 4),
                              "required_probes": 1, "envelope_residues": len(members),
                              "mutable_residues": len(mutable), "footprint_residues": len(footprint),
                              "surface_connectivity": "connected"})
    assigned = sum(len(region.envelope_residues) for region in regions)
    return regions, {
        "eligible_residues": sum(item.folding_screen_pass for item in envelope),
        "assigned_residues": assigned, "coverage": round(assigned / len(envelope), 4),
        "patch_radius_a": config.COMPACT_FOOTPRINT_GEODESIC_RADIUS_A,
        "actual_regions": len(regions), "assay_format": assay_format,
        "membrane_accessibility_model": membrane_geometry,
        "membrane_accessibility_scale_A": config.MEMBRANE_ACCESSIBILITY_SCALE_A,
        "excluded_flexible_terminal_residues": sum(item.flexible_terminal for item in records),
        "coverage_plan": coverage_plan, "partition_mode": "compact_surface_face",
        "surface_geometry_version": "solvent_atomic_points_v2",
        "partition_basis": partition_basis,
        "selected_face_count": selected_face_count,
        "trigger_reasons": trigger_reasons,
        "method": "Compact single-domain continuous surface-face partitioning with local Fab footprints",
        "warnings": ["Protected surface residues remain in geometry but are excluded from automatic mutation menus."],
    }


def partition_surface(
    records: list[ResidueRecord], protein: ProteinContext, confidence: dict[str, Any] | None = None,
    minimum_regions: int = 5, maximum_regions: int = 8, predicted: bool = False,
    assay_format: str = "cell_surface", membrane_anchor: tuple[float, float, float] | None = None,
    membrane_normal: tuple[float, float, float] | None = None,
) -> tuple[list[RegionRecord], dict[str, Any]]:
    """Route only pathological compact single-domain outputs to the rescue mode."""
    legacy, legacy_summary = _partition_surface_footprint(
        records, protein, confidence, minimum_regions, maximum_regions, predicted,
        assay_format, membrane_anchor, membrane_normal,
    )
    triggered, reasons, envelope = _compact_trigger(records, protein, legacy)
    if not triggered:
        legacy_summary["partition_mode"] = "footprint"
        legacy_summary["trigger_reasons"] = []
        return legacy, legacy_summary
    # The legacy trial assigned sparse mutable residues.  Clear it before
    # painting compact envelope membership so no stale legacy colour survives
    # on a low-confidence/protected residue in the viewer or offline report.
    for record in records:
        record.region_id = None
    compact = _partition_compact_surface_face(
        records, protein, confidence, predicted, assay_format, membrane_anchor,
        membrane_normal, envelope, reasons,
    )
    if compact is None:
        # Rebuild the legacy state because its residue IDs were cleared above.
        legacy, legacy_summary = _partition_surface_footprint(
            records, protein, confidence, minimum_regions, maximum_regions, predicted,
            assay_format, membrane_anchor, membrane_normal,
        )
        legacy_summary["partition_mode"] = "footprint"
        legacy_summary["trigger_reasons"] = reasons
        legacy_summary.setdefault("warnings", []).append(
            "Compact single-domain rescue was considered but failed connectivity quality gates; retained legacy footprint result."
        )
        return legacy, legacy_summary
    return compact
