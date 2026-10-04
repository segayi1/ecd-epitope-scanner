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
) -> list[V4Candidate]:
    if parameters.get("coverage_scan_v2"):
        return coverage_candidates(options, graph, parameters, diagnostics)
    by_position: dict[int, list[V4MutationOption]] = {}
    for option in options:
        by_position.setdefault(option.residue.uniprot_pos or -1, []).append(option)
    orders = list(dict.fromkeys(
        [int(item) for item in parameters.get("orders", [2, 3])]
        + [int(item) for item in parameters.get("exploratory_orders", [])]
    ))
    maximum_span = float(parameters.get("maximum_geodesic_span_a", 24.0))
    site_cap = int(parameters.get("site_pool_per_sector", 12))
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
    for sector, hypothesis_id in windows:
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
            positions = supported_positions[:supported_cap] + fallback_positions[:fallback_cap]
            remaining = [
                position
                for position in supported_positions[supported_cap:] + fallback_positions[fallback_cap:]
                if position not in positions
            ]
            positions.extend(remaining[: max(0, site_cap - len(positions))])
        else:
            positions = positions[:site_cap]
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
            candidates.sort(
                key=lambda item: (
                    item.structural_risk,
                    -item.epitope_score,
                    item.order,
                    item.candidate_id,
                )
            )
            beam_by_order = parameters.get("beam_width_by_order", {})
            beam_cap = beam_by_order.get(str(order), beam_by_order.get(order, parameters.get("beam_width", 300)))
            candidates = candidates[: max(1, int(beam_cap))]
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
                shortlisted = prioritized[:priority_cap] + fallback[:fallback_cap]
                remaining = [
                    item
                    for item in prioritized[priority_cap:] + fallback[fallback_cap:]
                    if item not in shortlisted
                ]
                shortlisted.extend(remaining[: max(0, order_cap - len(shortlisted))])
            else:
                shortlisted = candidates[:order_cap]
            order_shortlists.extend(shortlisted[:order_cap])
        output.extend(order_shortlists)
    if hypotheses:
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
        from .hypotheses import attach_hypothesis_responses
        attach_hypothesis_responses(list(unique.values()), hypotheses, graph)
        output = hypothesis_balanced_shortlist(list(unique.values()), hypotheses, parameters, diagnostics)
    return sorted(output, key=lambda item: (item.order, item.candidate_id))


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


def hypothesis_balanced_shortlist(candidates, hypotheses, parameters, diagnostics=None):
    """Quota is scientific-hypothesis based; display ownership is irrelevant."""
    output = []
    ids = [h.hypothesis_id for h in hypotheses]
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
            for candidate in available[:quota]:
                chosen[candidate.candidate_id] = candidate
            if len(available) < quota:
                missing[hypothesis_id] = quota - len(available)
        overrides = parameters.get("global_shortlist_by_order", {})
        requested = max(1, int(overrides.get(str(order), overrides.get(order, parameters.get("global_shortlist_per_order", 60)))))
        capacity = max(requested, len(chosen))
        reserved_count = len(chosen)
        # Fill by perturbation diversity within a stability-qualified pool.
        vectors = {c.candidate_id: [c.hypothesis_response.get(hid, 0.0) for hid in ids] for c in ranked}
        remaining = [c for c in ranked if c.candidate_id not in chosen]
        def vector_distance(left, right):
            return math.sqrt(sum((a-b)**2 for a,b in zip(vectors[left], vectors[right])) / max(1, len(ids)))
        nearest = {c.candidate_id: min((vector_distance(c.candidate_id, other) for other in chosen), default=math.inf) for c in remaining}
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
