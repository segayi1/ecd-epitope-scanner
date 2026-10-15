from __future__ import annotations

from itertools import combinations, product
import hashlib
import math
from typing import Any

from .models import V4Candidate, V4MutationOption, V4Sector, V4SurfaceHypothesis


def _distances(
    graph: dict[int, dict[int, float]], start: int, allowed: set[int]
) -> dict[int, float]:
    import heapq

    result = {start: 0.0}
    queue = [(0.0, start)]
    while queue:
        distance, current = heapq.heappop(queue)
        if distance != result.get(current):
            continue
        for neighbor, step in graph.get(current, {}).items():
            if neighbor not in allowed:
                continue
            proposed = distance + step
            if proposed < result.get(neighbor, float("inf")):
                result[neighbor] = proposed
                heapq.heappush(queue, (proposed, neighbor))
    return result


def _connected(positions: set[int], graph: dict[int, dict[int, float]], allowed: set[int]) -> bool:
    if not positions:
        return False
    seen = {min(positions)}
    changed = True
    while changed:
        changed = False
        for left in tuple(seen):
            reachable = _distances(graph, left, allowed)
            newly_seen = positions.intersection(reachable).difference(seen)
            if newly_seen:
                seen.update(newly_seen)
                changed = True
    return seen == positions


def generate_candidates(
    sectors: list[V4Sector],
    options: list[V4MutationOption],
    graph: dict[int, dict[int, float]],
    parameters: dict[str, Any],
    hypotheses: list[V4SurfaceHypothesis] | None = None,
    diagnostics: dict[str, Any] | None = None,
    surface_areas: dict[int, float] | None = None,
) -> list[V4Candidate]:
    if parameters.get("coverage_scan_v2"):
        return coverage_candidates(options, graph, parameters, diagnostics)
    from .spatial_support import (annotate_spatial_support, diverse_sites, spatial_shortlist,
                                  candidate_region_membership, regional_support, site_signature)
    guarantee = surface_areas is not None and parameters.get("sector_resource_guarantees", True) and not parameters.get("require_review")
    base_quotas = {int(order): int(quota) for order, quota in parameters.get("exclusive_groups_by_order", {"2": 4, "3": 2}).items()}
    coverage_trigger = float(parameters.get("sector_coverage_search_trigger", 0.9))
    by_position: dict[int, list[V4MutationOption]] = {}
    for option in options:
        by_position.setdefault(option.residue.uniprot_pos or -1, []).append(option)
    orders = list(dict.fromkeys(
        [int(item) for item in parameters.get("orders", [2, 3])]
        + [int(item) for item in parameters.get("exploratory_orders", [])]
    ))
    maximum_span = float(parameters.get("maximum_geodesic_span_a", 24.0))
    site_cap = int(parameters.get("site_pool_per_sector", 12))
    if hypotheses:
        site_cap = min(site_cap, int(parameters.get("site_pool_per_hypothesis", 8)))
    global_distances = {position: _distances(graph, position, set(graph))
                        for position in {option.residue.uniprot_pos for option in options}} if surface_areas is not None else {}
    per_sector = int(parameters.get("shortlist_per_sector", 30))
    epitope_site_fraction = min(
        1.0, max(0.0, float(parameters.get("epitope_site_pool_fraction", 0.5)))
    )
    fallback_candidate_fraction = min(
        1.0, max(0.0, float(parameters.get("fallback_candidate_fraction", 0.5)))
    )
    display_owner = {
        position: sector.sector_id
        for sector in sectors for position in sector.residue_positions
    }
    if hypotheses:
        windows = []
        for hypothesis in hypotheses:
            covered = [position for position in hypothesis.residue_weights if position in by_position and position in graph]
            all_positions = sorted(position for position in hypothesis.residue_weights if position in graph)
            if not covered or not all_positions:
                continue
            owner_counts: dict[str, int] = {}
            for position in all_positions:
                owner = display_owner.get(position)
                if owner:
                    owner_counts[owner] = owner_counts.get(owner, 0) + 1
            owner = min(owner_counts, key=lambda item: (-owner_counts[item], item)) if owner_counts else "UNASSIGNED"
            window = V4Sector(
                sector_id=owner, residue_positions=all_positions, component_id=None,
                designable_area_a2=sum(hypothesis.residue_weights[p] for p in covered),
                total_area_a2=sum(hypothesis.residue_weights[p] for p in all_positions),
                epitope_score=hypothesis.discotope_prior or 0.0,
            )
            windows.append((window, hypothesis.hypothesis_id))
    else:
        windows = [(sector, None) for sector in sectors]
    output = []
    generation_audit = []
    for sector, hypothesis_id in windows:
        scope_areas = {p: surface_areas[p] for p in sector.residue_positions if p in surface_areas} if surface_areas is not None else None
        positions = [position for position in sector.residue_positions if position in by_position]
        positions.sort(key=lambda p: (
            not any(option.review_required for option in by_position[p]) if parameters.get("require_review") else False,
            -max(option.accessibility_score for option in by_position[p]),
            -max(option.chemical_disruption_score for option in by_position[p]),
            min(option.spurs_ddg_kcal_mol if option.spurs_ddg_kcal_mol is not None else float("inf") for option in by_position[p]),
            p,
        ))
        supported_positions = [
            position
            for position in positions
            if any(option.epitope_supported for option in by_position[position])
        ]
        fallback_positions = [position for position in positions if position not in supported_positions]
        if supported_positions and fallback_positions:
            supported_cap = min(
                len(supported_positions), max(1, round(site_cap * epitope_site_fraction))
            )
            fallback_cap = min(len(fallback_positions), max(1, site_cap - supported_cap))
            supported = diverse_sites(supported_positions, supported_cap, global_distances, scope_areas)
            fallback = diverse_sites(fallback_positions, fallback_cap, global_distances, scope_areas, supported)
            positions = supported + fallback
            remaining = [
                position
                for position in supported_positions + fallback_positions
                if position not in positions
            ]
            positions.extend(diverse_sites(remaining, max(0, site_cap - len(positions)), global_distances, scope_areas, positions))
        else:
            positions = diverse_sites(positions, site_cap, global_distances, scope_areas)
        candidates: list[V4Candidate] = []
        order_shortlists: list[V4Candidate] = []
        if hypotheses:
            positions = positions[:min(site_cap, int(parameters.get("site_pool_per_hypothesis", 8)))]
        allowed_positions = set(sector.residue_positions)
        distance_by_site = {
            position: _distances(graph, position, set(graph) if parameters.get('global_span_paths') else allowed_positions)
            for position in positions
        }
        for order in orders:
            candidates = []
            if len(positions) < order:
                continue
            for site_combo in combinations(positions, order):
                if not set(site_combo) <= distance_by_site[site_combo[0]].keys():
                    continue
                distances = [
                    distance_by_site[left] for left in site_combo
                ]
                span = max(
                    (
                        distances[i].get(right, float("inf"))
                        for i, left in enumerate(site_combo)
                        for right in site_combo[i + 1 :]
                    ),
                    default=0.0,
                )
                if span > maximum_span:
                    continue
                # Coverage depends on sites and geometry, not amino-acid
                # choices. Reusing it avoids repeating Dijkstra for every
                # member of the Cartesian mutation menu (up to 81 per set).
                coverage, covered_positions = _coverage_from_distances(distances, allowed_positions)
                pools = [
                    sorted(by_position[position], key=lambda item: (
                        -item.chemical_disruption_score,
                        item.spurs_ddg_kcal_mol if item.spurs_ddg_kcal_mol is not None else float("inf"),
                        -item.accessibility_score,
                        item.notation,
                    ))[:3]
                    for position in site_combo
                ]
                for selected in product(*pools):
                    if parameters.get("require_review") and not any(item.review_required for item in selected):
                        continue
                    if len({item.mutant_aa for item in selected}) == 0:
                        continue
                    risk_values = [
                        item.spurs_ddg_kcal_mol
                        for item in selected
                        if item.spurs_ddg_kcal_mol is not None
                    ]
                    candidate = V4Candidate(
                        candidate_id=(
                            f"V4_{sector.sector_id}_{hypothesis_id}_{order}_{len(candidates) + 1:04d}"
                            if hypothesis_id else
                            f"V4_{sector.sector_id}_{order}_{len(candidates) + 1:04d}"
                        ),
                        sector_id=sector.sector_id,
                        mutations=list(selected),
                        order=order,
                        epitope_score=round(sum(item.epitope_prior for item in selected), 6),
                        structural_risk=round(sum(max(0.0, value) for value in risk_values), 6),
                        connected=True,
                        warnings=(
                            []
                            if order in {1, 2}
                            else [
                                f"HIGHER_ORDER_EXTRAPOLATION_{'MEDIUM' if order == 3 else 'HIGH'}"
                            ]
                        ),
                        epitope_supported_sites=sum(
                            1 for item in selected if item.epitope_supported
                        ),
                        recommendation_class=(
                            "EPITOPE_PRIORITIZED"
                            if any(item.epitope_supported for item in selected)
                            else "STABILITY_FALLBACK"
                        ),
                        anchor_region_id=sector.sector_id,
                        search_role="EXPLORATORY" if order >= 4 else "PRIMARY_POOL",
                        targeted_hypotheses=[hypothesis_id] if hypothesis_id else [],
                        qualification="REVIEW_REQUIRED" if any(item.review_required for item in selected) else "FORMAL_POOL",
                    )
                    candidate.coverage_by_scenario = dict(coverage)
                    candidate.covered_positions_by_scenario = {
                        scenario: list(hit) for scenario, hit in covered_positions.items()
                    }
                    candidates.append(candidate)
            annotate_spatial_support(candidates, graph, surface_areas, global_distances, allowed_positions,
                                     hypothesis_id or sector.sector_id)
            beam_by_order = parameters.get("beam_width_by_order", {})
            beam_cap = beam_by_order.get(str(order), beam_by_order.get(order, parameters.get("beam_width", 300)))
            candidates = spatial_shortlist(candidates, max(1, int(beam_cap)))
            prioritized = [
                item for item in candidates if item.recommendation_class == "EPITOPE_PRIORITIZED"
            ]
            fallback = [
                item for item in candidates if item.recommendation_class == "STABILITY_FALLBACK"
            ]
            shortlist_by_order = parameters.get("shortlist_by_order", {})
            default_order_cap = parameters.get(
                "shortlist_per_order", max(1, per_sector // max(1, len(orders)))
            )
            order_cap = max(
                1,
                int(shortlist_by_order.get(str(order), shortlist_by_order.get(order, default_order_cap))),
            )
            if prioritized and fallback:
                fallback_cap = min(
                    len(fallback), max(1, round(order_cap * fallback_candidate_fraction))
                )
                fallback_cap = min(fallback_cap, max(1, order_cap - 1))
                priority_cap = min(len(prioritized), max(1, order_cap - fallback_cap))
                shortlisted = spatial_shortlist(prioritized, priority_cap) + spatial_shortlist(fallback, fallback_cap)
                remaining = [
                    item
                    for item in prioritized + fallback
                    if item not in shortlisted
                ]
                shortlisted.extend(spatial_shortlist(remaining, max(0, order_cap - len(shortlisted))))
            else:
                shortlisted = spatial_shortlist(candidates, order_cap)
            order_shortlists.extend(shortlisted[:order_cap])
        output.extend(order_shortlists)
    if guarantee:
        for sector in sorted(sectors, key=lambda item: item.sector_id):
            exclusive = [c for c in output if candidate_region_membership(c, sectors)["exclusive_region_id"] == sector.sector_id]
            counts = {order: len({site_signature(c) for c in exclusive if c.order == order}) for order in (2, 3)}
            best = max((regional_support(c, sector, surface_areas, global_distances)["target_mean_coverage_fraction"] or 0.0
                        for c in exclusive if c.order in (2, 3)), default=0.0)
            supplemental_orders = [order for order in orders if order in (2, 3)]
            supplement = bool(supplemental_orders) and (any(counts[order] < base_quotas.get(order, 0) for order in supplemental_orders) or best < coverage_trigger)
            if supplement:
                supplemental_parameters = {**parameters, "sector_resource_guarantees": False,
                    "orders": supplemental_orders, "exploratory_orders": [],
                    "shortlist_by_order": {str(order): max(base_quotas.get(order, 0), int(parameters.get("shortlist_per_order", 10)))
                                           for order in supplemental_orders}}
                output.extend(generate_candidates([sector], options, graph, supplemental_parameters,
                                                  surface_areas=surface_areas))
            generation_audit.append({"sector_id": sector.sector_id, "supplemental_window": supplement,
                                     "initial_exclusive_site_sets": {str(order): count for order, count in counts.items()},
                                     "initial_best_mean_coverage": best})
    if hypotheses or guarantee:
        # A construct can arise from overlapping footprints with different
        # display anchors. Its physical identity is global, not sector-local.
        unique: dict[tuple[tuple[int, str], ...], V4Candidate] = {}
        for candidate in sorted(output, key=lambda item: (
            item.structural_risk, -item.epitope_score, item.candidate_id
        )):
            signature = tuple(sorted(
                (option.residue.uniprot_pos or -1, option.mutant_aa)
                for option in candidate.mutations
            ))
            existing = unique.get(signature)
            if existing is None:
                unique[signature] = candidate
            else:
                existing.targeted_hypotheses = sorted(set(
                    existing.targeted_hypotheses + candidate.targeted_hypotheses
                ))
        for signature, candidate in unique.items():
            candidate.candidate_id = "V41_" + hashlib.sha256(repr(signature).encode()).hexdigest()[:16]
        # Global ranking compares one surface denominator, not the arbitrary
        # first overlapping window retained by physical-identity deduplication.
        annotate_spatial_support(list(unique.values()), graph, surface_areas, global_distances,
                                 scope_id="full_eligible_surface")
        if hypotheses:
            from .hypotheses import attach_hypothesis_responses
            attach_hypothesis_responses(list(unique.values()), hypotheses, graph)
        generated = list(unique.values())
        output = hypothesis_balanced_shortlist(generated, hypotheses or [], parameters, diagnostics,
                    sectors=sectors if guarantee else None, surface_areas=surface_areas, distances=global_distances)
        if diagnostics is not None and guarantee:
            for row in generation_audit:
                for label, pool in (("generated", generated), ("pool", output)):
                    exclusive = [c for c in pool if candidate_region_membership(c, sectors)["exclusive_region_id"] == row["sector_id"]]
                    row[label + "_exclusive_candidates"] = len(exclusive)
                    row[label + "_exclusive_site_sets"] = {str(order): len({site_signature(c) for c in exclusive if c.order == order}) for order in (2, 3)}
                row["quota_deficits"] = {str(order): max(0, base_quotas.get(order, 0) - row["pool_exclusive_site_sets"][str(order)]) for order in orders if order in (2, 3)}
                row["generation_count_basis"] = "deduplicated_local_shortlists_after_beam_before_global_pool"
            diagnostics["sector_resources"] = generation_audit
    if diagnostics is not None and surface_areas is not None:
        diagnostics.update(spatial_selection_policy="pareto_target_coverage_marginal_area_single_risk_then_site_set_round_robin_v2",
                           geometric_support_basis="full_eligible_surface_graph_9_12_15_equal_weight_mean",
                           geometric_support_interpretation="proxy only; not biological perturbation",
                           retained_unique_site_sets=len({tuple(sorted(o.residue.uniprot_pos for o in c.mutations)) for c in output}))
    return sorted(output, key=lambda item: (item.order, item.candidate_id))


def supplemental_exclusive_doubles(candidates, sectors, options, graph, parameters, areas):
    """One bounded search round for regions with no SPURS-acceptable exclusive pair."""
    from .spatial_support import candidate_region_membership, spatial_shortlist
    known = {tuple(sorted(o.notation for o in c.mutations)) for c in candidates}
    output, audit = [], []
    for sector in sorted(sectors, key=lambda s: s.sector_id):
        exclusive = [c for c in candidates if c.order == 2 and
            candidate_region_membership(c, sectors)["exclusive_region_id"] == sector.sector_id]
        if any(c.stability is not None and c.stability.consensus_status in
               {"CONSENSUS_PASS", "SPURS_ONLY_REVIEW_REQUIRED"} for c in exclusive):
            continue
        pool = generate_candidates([sector], options, graph, {
            **parameters, "orders": [2], "exploratory_orders": [],
            "shortlist_by_order": {"2": 300}, "global_shortlist_by_order": {"2": 300},
            "sector_resource_guarantees": True,
        }, surface_areas=areas)
        available = [c for c in pool if tuple(sorted(o.notation for o in c.mutations)) not in known]
        chosen = spatial_shortlist(available, 12)
        output.extend(chosen)
        known.update(tuple(sorted(o.notation for o in c.mutations)) for c in chosen)
        audit.append({"sector_id": sector.sector_id, "trigger": "NO_SPURS_ACCEPTABLE_EXCLUSIVE_DOUBLE",
            "previous_exclusive_doubles": len(exclusive), "new_candidates": len(chosen),
            "maximum_new_candidates": 12, "rounds": 1})
    return output, audit


def add_pair_controls(candidates, graph):
    """Actual double-mutant controls for every retained triple, not unions."""
    from dataclasses import replace
    existing={tuple(sorted(o.notation for o in c.mutations)):c for c in candidates}
    for parent in list(candidates):
        if parent.order!=3 or parent.qualification=='REVIEW_REQUIRED':continue
        for options in combinations(parent.mutations,2):
            key=tuple(sorted(o.notation for o in options))
            if key in existing:
                control=existing[key]
                if parent.candidate_id not in control.control_for_candidate_ids:
                    control.control_for_candidate_ids.append(parent.candidate_id)
                continue
            signature=tuple(sorted((o.residue.uniprot_pos,o.mutant_aa) for o in options))
            child=replace(parent,candidate_id='V41_'+hashlib.sha256(repr(signature).encode()).hexdigest()[:16],
                mutations=list(options),order=2,stability=None,formal_eligible=False,formal_exclusion_reasons=[],
                warnings=['REAL_TRIPLE_LOWER_CONTROL'],structural_risk=sum(max(0.,o.spurs_ddg_kcal_mol or 0.) for o in options),
                epitope_score=sum(o.epitope_prior for o in options),triple_role=None,parent_candidate_ids=[],
                control_for_candidate_ids=[parent.candidate_id],search_role='EXACT_TRIPLE_PAIR_CONTROL',risk_components={})
            distances=[_distances(graph,o.residue.uniprot_pos,set(graph)) for o in options]
            child.risk_components['geodesic_span_a']=distances[0].get(options[1].residue.uniprot_pos,math.inf)
            child.coverage_by_scenario,child.covered_positions_by_scenario=_coverage_from_distances(distances,set(graph))
            existing[key]=child
    return sorted(existing.values(),key=lambda c:(c.order,c.candidate_id))


def hypothesis_balanced_shortlist(candidates, hypotheses, parameters, diagnostics=None,
                                 sectors=None, surface_areas=None, distances=None):
    """Preserve hypothesis capacity and add distinct exclusive-region reserves."""
    from .spatial_support import (spatial_shortlist, candidate_region_membership,
                                  site_signature, regional_priority, regional_support)
    output = []
    ids = sorted(h.hypothesis_id for h in hypotheses)
    base_quotas = {int(order): int(quota) for order, quota in parameters.get("exclusive_groups_by_order", {"2": 4, "3": 2}).items()}
    ownership = {c.candidate_id: candidate_region_membership(c, sectors)["exclusive_region_id"] for c in candidates} if sectors else {}
    local_priority = {}
    local_coverage = {}
    if sectors:
        for sector in sectors:
            for candidate in candidates:
                if ownership[candidate.candidate_id] == sector.sector_id:
                    local_priority[candidate.candidate_id] = regional_priority(candidate, sector, surface_areas, distances)
                    local_coverage[candidate.candidate_id] = -local_priority[candidate.candidate_id][0]
    audit = []
    groups = {(c.qualification, c.order) for c in candidates}
    qualification = "REVIEW_REQUIRED" if parameters.get("require_review", False) else "FORMAL_POOL"
    groups.update((qualification, int(order)) for order in parameters.get("orders", [2, 3]))
    for qualification, order in sorted(groups):
        ranked = sorted((c for c in candidates if c.order == order and c.qualification == qualification),
                        key=lambda c: (c.structural_risk, -c.epitope_score, c.candidate_id))
        quota = 2 if order == 2 else 1 if order in {1, 3} else 0
        chosen = {}
        missing = {}
        for hypothesis_id in ids:
            available = [c for c in ranked if hypothesis_id in c.targeted_hypotheses]
            for candidate in spatial_shortlist(available, quota):
                chosen[candidate.candidate_id] = candidate
            if len(available) < quota:
                missing[hypothesis_id] = quota - len(available)
        overrides = parameters.get("global_shortlist_by_order", {})
        requested = max(1, int(overrides.get(str(order), overrides.get(order, parameters.get("global_shortlist_per_order", 60)))))
        capacity = max(requested, len(chosen))
        reserved_count = len(chosen)
        base_missing = 0
        if sectors and order in (2, 3):
            for sector in sorted(sectors, key=lambda s: s.sector_id):
                exclusive = [c for c in ranked if ownership[c.candidate_id] == sector.sector_id]
                representatives = {}
                for candidate in sorted(exclusive, key=lambda c: local_priority[c.candidate_id]):
                    representatives.setdefault(site_signature(candidate), candidate)
                for candidate in list(representatives.values())[:base_quotas.get(order, 0)]:
                    if candidate.candidate_id not in chosen:
                        base_missing += 1
                    chosen[candidate.candidate_id] = candidate
            capacity += base_missing
        # Fill by perturbation diversity within a stability-qualified pool.
        vectors = {c.candidate_id: [c.hypothesis_response.get(hid, 0.0) for hid in ids] for c in ranked}
        remaining = [c for c in ranked if c.candidate_id not in chosen]
        if sectors and order in (2, 3):
            # Recompute deficits after each distinct regional opportunity.
            def region_need(sector):
                exclusive = [c for c in chosen.values() if ownership[c.candidate_id] == sector.sector_id]
                quota_gap = max(0, base_quotas.get(order, 0) - len({site_signature(c) for c in exclusive}))
                best = max((local_coverage[c.candidate_id] for c in exclusive), default=0.0)
                area = sum(surface_areas.get(p, 0.0) for p in sector.residue_positions)
                return (-quota_gap, -(1 - best), -(1 - best) * area, sector.sector_id)
            while len(chosen) < capacity:
                existing_sets = {site_signature(c) for c in chosen.values()}
                available_by_region = [(sector, [c for c in remaining if site_signature(c) not in existing_sets and
                    ownership[c.candidate_id] == sector.sector_id])
                    for sector in sectors]
                available_by_region = [(sector, available) for sector, available in available_by_region
                                       if available and region_need(sector)[1] < 0]
                if not available_by_region:
                    break
                sector, available = min(available_by_region, key=lambda pair: region_need(pair[0]))
                winner = min(available, key=lambda c: local_priority[c.candidate_id])
                chosen[winner.candidate_id] = winner
                remaining.remove(winner)
        def vector_distance(left, right):
            return math.sqrt(sum((a-b)**2 for a,b in zip(vectors[left], vectors[right])) / max(1, len(ids)))
        nearest = {c.candidate_id: min((vector_distance(c.candidate_id, other) for other in chosen), default=math.inf) for c in remaining}
        if any(c.geometry_support for c in remaining):
            for candidate in spatial_shortlist(remaining, max(0, capacity - len(chosen))):
                chosen[candidate.candidate_id] = candidate
            remaining = []
        while remaining and len(chosen) < capacity:
            def key(candidate):
                distance = nearest[candidate.candidate_id] if chosen else 0.0
                return (-distance, candidate.structural_risk, -candidate.epitope_score, candidate.candidate_id)
            winner = min(remaining, key=key)
            chosen[winner.candidate_id] = winner
            remaining.remove(winner)
            for candidate in remaining:
                nearest[candidate.candidate_id] = min(nearest[candidate.candidate_id], vector_distance(candidate.candidate_id, winner.candidate_id))
        output.extend(chosen.values())
        audit.append({"qualification": qualification, "order": order, "requested_capacity": requested,
                      "effective_capacity": capacity, "reserved_unique": reserved_count,
                      "exclusive_base_expansion": base_missing,
                      "selected": len(chosen), "quota_deficits": missing})
    if diagnostics is not None:
        diagnostics.update(global_order_pools=audit, physical_candidates_before_global_shortlist=len(candidates))
    return output


def _coverage(
    candidate: V4Candidate, sector: V4Sector, graph: dict[int, dict[int, float]]
) -> tuple[dict[str, float], dict[str, list[int]]]:
    positions = set(sector.residue_positions)
    distances = [
        _distances(graph, mutation.residue.uniprot_pos or -1, positions)
        for mutation in candidate.mutations
    ]
    return _coverage_from_distances(distances, positions)


def _coverage_from_distances(
    distances: list[dict[int, float]], positions: set[int],
) -> tuple[dict[str, float], dict[str, list[int]]]:
    covered: dict[str, float] = {}
    covered_positions: dict[str, list[int]] = {}
    for radius in (0.0, 9.0, 12.0, 15.0):
        hit = set()
        for distance in distances:
            hit |= {position for position, value in distance.items() if value <= radius}
        scenario = f"geodesic_{int(radius)}"
        covered[scenario] = round(len(hit.intersection(positions)) / max(1, len(positions)), 6)
        covered_positions[scenario] = sorted(hit)
    return covered, covered_positions


def coverage_candidates(options, graph, parameters, diagnostics=None):
    """Enumerate full eligible site graph, then select a bounded proxy pool."""
    from .coverage_profile import sequence_reasons, mechanisms, descriptors
    by_site = {}
    for option in options:
        if (option.eligible and not option.review_required and option.residue.uniprot_pos in graph
                and option.spurs_ddg_kcal_mol is not None and math.isfinite(option.spurs_ddg_kcal_mol)):
            by_site.setdefault(option.residue.uniprot_pos, []).append(option)
    allowed = set(graph)
    distances = {p: _distances(graph, p, allowed) for p in by_site}
    pool = []
    rejected = 0
    span_limit = float(parameters.get("maximum_geodesic_span_a", 24))
    sequence = parameters["canonical_sequence"]
    for sites in combinations(sorted(by_site), 2):
        span = distances[sites[0]].get(sites[1], math.inf)
        if span > span_limit:
            continue
        coverage, support = _coverage_from_distances([distances[p] for p in sites], allowed)
        supports=[{p for p,d in distances[site].items() if d<=12} for site in sites]
        overlap=len(supports[0]&supports[1])/max(1,len(supports[0]|supports[1]))
        for selected in product(*(by_site[p] for p in sites)):
            if sequence_reasons(sequence, [(o.residue.uniprot_pos, o.residue.one_letter, o.mutant_aa) for o in selected]):
                rejected += 1
                continue
            signature = "|".join(sorted(o.notation for o in selected))
            identifier = "CV2_" + hashlib.sha256(signature.encode()).hexdigest()[:16]
            candidate = V4Candidate(identifier, "GLOBAL", list(selected), order=2,
                structural_risk=sum(max(0, o.spurs_ddg_kcal_mol) for o in selected),
                coverage_by_scenario=coverage, covered_positions_by_scenario=support,
                risk_components={"geodesic_span_a": span,"proxy_support_overlap_fraction":overlap}, search_role="DOUBLE_SCAN")
            pool.append(candidate)
    enumerated_count=len(pool)
    full_position_pair_count=len({tuple(sorted(o.residue.uniprot_pos for o in c.mutations)) for c in pool})
    all_site_pairs=sorted({tuple(sorted(o.residue.uniprot_pos for o in c.mutations)) for c in pool})
    excluded_ids=set(parameters.get('exclude_candidate_ids',[]))
    pool=[c for c in pool if c.candidate_id not in excluded_ids]
    excluded_count=enumerated_count-len(pool)
    position_pair_count=len({tuple(sorted(o.residue.uniprot_pos for o in c.mutations)) for c in pool})
    if parameters.get('site_pair_priority',False):
        from .coverage_profile import vector
        excluded_pairs={tuple(sorted(pair)) for pair in parameters.get('exclude_site_pairs',[])}
        grouped={}
        def representative_rank(candidate):
            descriptions=[descriptors(o.residue.one_letter,o.mutant_aa) for o in candidate.mutations]
            return (-len(set().union(*(mechanisms(d) for d in descriptions))),
                -sum(sum(abs(v) for v in vector(d)) for d in descriptions),
                candidate.structural_risk,candidate.candidate_id)
        for candidate in pool:
            pair=tuple(sorted(o.residue.uniprot_pos for o in candidate.mutations))
            if pair in excluded_pairs:continue
            if pair not in grouped or representative_rank(candidate)<representative_rank(grouped[pair]):
                grouped[pair]=candidate
        pool=list(grouped.values())
    cap = int(parameters.get("global_shortlist", 32))
    if not 1 <= cap <= 256:
        raise ValueError("Coverage pool capacity must be in [1, 256]")
    selected = []
    covered = set()
    mandatory={tuple(sorted(v)) for v in parameters.get('mandatory_candidate_signatures',[])}
    requested_ids=set(parameters.get('mandatory_candidate_ids',[]))
    mandatory.update(tuple(sorted(o.notation for o in c.mutations)) for c in pool if c.candidate_id in requested_ids)
    for signature in sorted(mandatory):
        matches=[c for c in pool if tuple(sorted(o.notation for o in c.mutations))==signature]
        if len(matches)!=1:raise ValueError('Mandatory candidate failed full eligibility/geometry audit')
        matches[0].search_role='STAGE_WHITELIST' if matches[0].candidate_id in requested_ids else 'SCIENTIFIC_VALIDATION_CONTROL'
        selected.append(matches[0]);pool.remove(matches[0])
        covered.update(matches[0].covered_positions_by_scenario['geodesic_12'])
    if len(selected)>cap:raise ValueError('Mandatory controls exceed hard pool capacity')
    chemistry_by_id = {c.candidate_id: len(set().union(*(mechanisms(descriptors(o.residue.one_letter, o.mutant_aa)) for o in c.mutations))) for c in pool}
    # Coverage diversity must precede chemistry variants at already sampled sites.
    while pool and len(selected) < cap:
        def rank(c):
            hit = set(c.covered_positions_by_scenario["geodesic_12"])
            return (-len(hit - covered), -chemistry_by_id[c.candidate_id], c.structural_risk, c.candidate_id)
        winner = min(pool, key=rank)
        pool.remove(winner)
        selected.append(winner)
        covered.update(winner.covered_positions_by_scenario["geodesic_12"])
    if diagnostics is not None:
        diagnostics.update(route="coverage_scan_v2", eligible_sites=len(by_site),
                           sequence_rejected=rejected, full_enumerated_count=enumerated_count,
                           excluded_screened_count=excluded_count,
                           position_pair_count_before_exclusion=position_pair_count,
                           full_position_pair_count=full_position_pair_count,
                           all_eligible_site_pairs=[list(pair) for pair in all_site_pairs],
                           remaining_position_pair_count=len({tuple(sorted(o.residue.uniprot_pos for o in c.mutations)) for c in pool+selected}),
                           site_pair_priority=bool(parameters.get('site_pair_priority',False)),
                           hard_pool_capacity=cap, selected_count=len(selected),
                           proxy_denominator="full_surface_graph_residues")
    if parameters.get("triple_parent_ids"):
        from ..surface_perturbation.compression import merge_parents
        parent_ids = set(parameters["triple_parent_ids"])
        parents = [c for c in selected if c.candidate_id in parent_ids]
        if parent_ids - {c.candidate_id for c in parents}:
            raise ValueError("Triple parents missing from frozen pool")
        triples = []
        for left, right in combinations(parents, 2):
            mutations = merge_parents(left, right)
            if mutations is None or sequence_reasons(sequence, [(o.residue.uniprot_pos,o.residue.one_letter,o.mutant_aa) for o in mutations]):
                continue
            positions = [o.residue.uniprot_pos for o in mutations]
            span = max(distances[p].get(q, math.inf) for p,q in combinations(positions,2))
            remote=span>span_limit
            if remote and not parameters.get('allow_remote_triples',False):
                continue
            signature = "|".join(sorted(o.notation for o in mutations))
            coverage, support = _coverage_from_distances([distances[p] for p in positions], allowed)
            triples.append(V4Candidate("CV2_"+hashlib.sha256(signature.encode()).hexdigest()[:16],
                "GLOBAL", mutations, order=3, coverage_by_scenario=coverage,
                covered_positions_by_scenario=support,
                triple_role=parameters.get("triple_role", "TRIPLE_COMPRESSION"),
                parent_candidate_ids=[left.candidate_id,right.candidate_id],
                risk_components={'geodesic_span_a':span},
                warnings=['DISTANT_LAYOUT_SPLIT_LOCALIZATION_REQUIRED'] if remote else [],
                search_role='REMOTE_TRIPLE_BACKUP' if remote else 'LOCAL_TRIPLE'))
        unique = {c.candidate_id:c for c in triples}
        selected.extend(list(unique.values())[:min(2,int(parameters.get("triple_pool_cap", 2)))])
    return selected
