"""Area-weighted geometric support; never measured mutation perturbation."""
from itertools import combinations
import math

from .candidate_search import _distances


def site_signature(candidate):
    return tuple(sorted(option.residue.uniprot_pos for option in candidate.mutations))


def candidate_region_membership(candidate, sectors):
    """Actual mutated positions determine ownership, never a window anchor."""
    sites = set(site_signature(candidate))
    owners = {p: sector.sector_id for sector in sectors for p in sector.residue_positions}
    region_ids = sorted({owners[p] for p in sites if p in owners})
    exclusive = region_ids[0] if len(region_ids) == 1 and sites <= owners.keys() else None
    return {"region_ids": region_ids, "exclusive_region_id": exclusive,
            "assignment": "EXCLUSIVE" if exclusive else "CROSS_REGION" if len(region_ids) > 1 else "UNASSIGNED"}


def regional_support(candidate, sector, areas, distances=None, graph=None):
    if areas is None:
        return {}
    if distances is None:
        distances = {p: _distances(graph, p, set(graph)) for p in site_signature(candidate)}
    result = support_metrics(site_signature(candidate), distances, areas,
                             sector.residue_positions, sector.sector_id)
    local = support_metrics(site_signature(candidate), distances,
                            {p: areas[p] for p in sector.residue_positions if p in areas},
                            sector.residue_positions, sector.sector_id)
    result["full_surface_mean_marginal_area_a2"] = result["mean_marginal_area_a2"]
    result["mean_marginal_area_a2"] = local["mean_marginal_area_a2"]
    result["mean_overlap_fraction"] = local["mean_overlap_fraction"]
    for row in result["scenarios"].values():
        row["target_uncovered_positions"] = sorted(set(sector.residue_positions) - set(row["covered_positions"]))
    for name, row in result["scenarios"].items():
        row["local_marginal_area_a2"] = local["scenarios"][name]["marginal_area_a2"]
        row["local_mean_pair_overlap_fraction"] = local["scenarios"][name]["mean_pair_overlap_fraction"]
    return result


def regional_priority(candidate, sector, areas, distances=None, graph=None):
    support = regional_support(candidate, sector, areas, distances, graph)
    return (-(support.get("target_mean_coverage_fraction") or 0.0),
            -support.get("mean_marginal_area_a2", 0.0), support.get("mean_overlap_fraction", 0.0),
            -support.get("direct_area_a2", 0.0), _initial_risk(candidate),
            -candidate.epitope_score, candidate.candidate_id)


def support_metrics(positions, distances, areas, target_positions=None, scope_id=None):
    positions = tuple(sorted(set(positions)))
    scenarios = {}
    target = set(areas) if target_positions is None else set(target_positions) & set(areas)
    def area(nodes):
        return sum(areas.get(position, 0.0) for position in sorted(nodes))
    target_area = area(target)
    for radius in (9, 12, 15):
        supports = [{p for p, distance in distances[position].items() if distance <= radius and p in areas}
                    for position in positions]
        union = set().union(*supports)
        lower_areas = [area(set().union(*(support for j, support in enumerate(supports) if j != i)))
                       for i in range(len(supports))]
        lower = max(lower_areas, default=0.0)
        overlaps = [area(left & right) / area(left | right) if area(left | right) else 0.0
                    for left, right in combinations(supports, 2)]
        scenarios[f"geodesic_{radius}"] = {
            "union_area_a2": area(union), "best_lower_order_area_a2": lower,
            "covered_positions": sorted(union),
            "target_covered_area_a2": area(union & target),
            "target_coverage_fraction": area(union & target) / target_area if target_area else None,
            "marginal_area_a2": max(0.0, area(union) - lower),
            "leave_one_out_gain_a2": {str(position): max(0.0, area(union) - lower_areas[i])
                                      for i, position in enumerate(positions)},
            "mean_pair_overlap_fraction": sum(overlaps) / max(1, len(overlaps)),
            "max_pair_overlap_fraction": max(overlaps, default=0.0),
        }
    return {
        "basis": "surface_geodesic_neighborhood_proxy", "area_basis": "eligible_surface_area_a2",
        "target_positions": sorted(target),
        "scenarios": scenarios,
        "mean_marginal_area_a2": sum(row["marginal_area_a2"] for row in scenarios.values()) / 3,
        "mean_overlap_fraction": sum(row["mean_pair_overlap_fraction"] for row in scenarios.values()) / 3,
        "direct_area_a2": area(positions),
        "target_scope_id": scope_id,
        "target_scope_area_a2": target_area,
        "target_mean_coverage_fraction": sum(row["target_coverage_fraction"] for row in scenarios.values()) / 3 if target_area else None,
    }


def annotate_spatial_support(candidates, graph, areas, distances=None, target_positions=None, scope_id=None):
    if areas is None:
        return
    if any(not math.isfinite(value) or value < 0 for value in areas.values()):
        raise ValueError("Spatial support requires finite nonnegative surface areas")
    if distances is None:
        distances = {position: _distances(graph, position, set(graph))
                     for position in {p for c in candidates for p in site_signature(c)}}
    cache = {}
    for candidate in candidates:
        signature = site_signature(candidate)
        if signature not in cache:
            cache[signature] = support_metrics(signature, distances, areas, target_positions, scope_id)
        candidate.geometry_support = cache[signature]


def pointwise_perturbation_field(candidate, areas, epitope_priors=None):
    """Build an auditable per-surface-point perturbation proxy.

    The physical channel is mutation utility carried over the candidate's
    geodesic-12 support. It is a surface support proxy, not a measured field.
    The epitope channel uses the target surface residue's score; missing
    scores retain the neutral value 0.5.
    """
    from .regional_perturbation import PERTURBATION_AXIS_SCALES, candidate_impact

    areas = {int(position): float(area) for position, area in (areas or {}).items()}
    if any(not math.isfinite(area) or area < 0 for area in areas.values()):
        raise ValueError("Perturbation fields require finite nonnegative surface areas")
    geometry = getattr(candidate, "geometry_support", None) or {}
    scenario = (geometry.get("scenarios") or {}).get("geodesic_12", {})
    support = set(scenario.get("covered_positions") or
                  (getattr(candidate, "covered_positions_by_scenario", None) or {}).get("geodesic_12", ()))
    target_positions = geometry.get("target_positions")
    if target_positions is not None:
        support &= set(target_positions)
    impact = candidate_impact(candidate)
    utility = None
    if impact.get("complete"):
        strength = sum(float(impact["mean_absolute_delta"][axis]) / scale
                       for axis, scale in PERTURBATION_AXIS_SCALES.items())
        exposure = float(impact["exposure_mean"])
        utility = (strength * max(1, int(impact["order"])) /
                   len(PERTURBATION_AXIS_SCALES) * (.5 + .5 * min(1., max(0., exposure))))
    priors = {int(position): float(value) for position, value in (epitope_priors or {}).items()
              if value is not None and math.isfinite(float(value))}
    supported = support & set(areas)
    scored = supported & set(priors)
    epitope_status = ("SCORED" if scored and scored == supported else
                      "PARTIAL" if scored else "NEUTRAL_MISSING")
    physical = {str(position): (float(utility) if utility is not None and position in support else 0.0)
                for position in sorted(areas)}
    weighted = {position: value * (.5 + priors.get(int(position), .5))
                for position, value in physical.items()}
    return {
        "schema_version": "surface_perturbation_field_v2",
        "basis": "geodesic_12_surface_support_proxy",
        "status": "SCORED" if utility is not None else "INCOMPLETE",
        "epitope_status": epitope_status,
        "epitope_mapping_basis": "target_surface_residue_prior",
        "physical_field": physical,
        "epitope_weighted_field": weighted,
        "surface_area_a2": {str(position): areas[position] for position in sorted(areas)},
        "supported_positions": sorted(support & set(areas)),
    }


def annotate_pointwise_perturbation_fields(candidates, areas, epitope_priors=None):
    """Attach point fields without changing eligibility or scientific gates."""
    for candidate in candidates:
        candidate.impact_evidence = dict(getattr(candidate, "impact_evidence", None) or {})
        candidate.impact_evidence["surface_perturbation_field"] = pointwise_perturbation_field(candidate, areas, epitope_priors)
    return candidates


def annotate_higher_order_regional_gain(candidates, sectors, areas):
    """Compare each higher-order candidate with the best lower-order geometry in its regions."""
    owners = {}
    sector_positions = {}
    for sector in sectors:
        positions = set(sector.residue_positions)
        sector_positions[sector.sector_id] = positions
        for position in positions:
            owners.setdefault(position, set()).add(sector.sector_id)

    def area_for_scope(scope):
        if areas is None:
            return None
        values = []
        for position in sorted(scope):
            if position not in areas or areas[position] is None:
                return None
            try:
                value = float(areas[position])
            except (TypeError, ValueError):
                return None
            if not math.isfinite(value) or value < 0:
                return None
            values.append(value)
        total = sum(values)
        return total if total > 0 else None

    def coverage(candidate, scope, scope_area):
        scenarios = (getattr(candidate, "geometry_support", None) or {}).get("scenarios", {})
        fractions = []
        for radius in (9, 12, 15):
            row = scenarios.get(f"geodesic_{radius}")
            covered = row.get("covered_positions") if isinstance(row, dict) else None
            if not isinstance(covered, (list, tuple, set)):
                return None
            covered = set(covered)
            fractions.append(sum(float(areas[position]) for position in scope & covered) / scope_area)
        return sum(fractions) / len(fractions)

    for candidate in candidates:
        if candidate.order <= 2:
            continue
        components = dict(getattr(candidate, "risk_components", None) or {})
        positions = {option.residue.uniprot_pos for option in candidate.mutations}
        region_ids = sorted({region_id for position in positions for region_id in owners.get(position, ())})
        scope = set().union(*(sector_positions[region_id] for region_id in region_ids)) if region_ids else set()
        components.update({
            "regional_higher_order_gain": None,
            "regional_lower_order_coverage": None,
            "regional_candidate_coverage": None,
            "regional_gain_baseline_candidate_id": None,
            "regional_gain_region_ids": region_ids,
            "regional_gain_basis": "mean_regional_area_coverage_9_12_15",
            "regional_gain_status": None,
            "regional_gain_interpretation": "counterfactual_geometry_only_not_candidate_qualification",
        })
        candidate.risk_components = components
        if not positions or any(not owners.get(position) for position in positions):
            components["regional_gain_status"] = "MISSING_REGION_OWNERSHIP"
            continue
        scope_area = area_for_scope(scope)
        if scope_area is None:
            components["regional_gain_status"] = "MISSING_AREA"
            continue
        candidate_coverage = coverage(candidate, scope, scope_area)
        if candidate_coverage is None:
            components["regional_gain_status"] = "MISSING_CANDIDATE_GEOMETRY"
            continue
        components["regional_candidate_coverage"] = candidate_coverage

        lower_order = [item for item in candidates
                       if item is not candidate and item.order in {2, 3}
                       and item.order < candidate.order and item.connected is True
                       and item.mutations
                       and {option.residue.uniprot_pos for option in item.mutations} <= scope]
        if not lower_order:
            components["regional_gain_status"] = "NO_LOWER_ORDER_BASELINE"
            continue
        ranked = []
        for item in lower_order:
            value = coverage(item, scope, scope_area)
            if value is None:
                ranked = None
                break
            ranked.append(( -value, str(item.candidate_id), value, item.candidate_id))
        if ranked is None:
            components["regional_gain_status"] = "MISSING_LOWER_ORDER_GEOMETRY"
            continue
        _, _, lower_coverage, baseline_id = min(ranked)
        components.update({
            "regional_higher_order_gain": candidate_coverage - lower_coverage,
            "regional_lower_order_coverage": lower_coverage,
            "regional_gain_baseline_candidate_id": baseline_id,
            "regional_gain_status": "COMPUTED",
        })


def spatial_priority(candidate):
    support = candidate.geometry_support
    return (-(support.get("target_mean_coverage_fraction") or 0.0),
            -support.get("mean_marginal_area_a2", 0.0), support.get("mean_overlap_fraction", 0.0),
            -support.get("direct_area_a2", 0.0), _initial_risk(candidate),
            -candidate.epitope_score, candidate.order, candidate.candidate_id)


def _initial_risk(candidate):
    if any(option.spurs_ddg_kcal_mol is None or not math.isfinite(option.spurs_ddg_kcal_mol)
           for option in candidate.mutations):
        return math.inf
    return candidate.structural_risk


def spatial_shortlist(candidates, maximum):
    """Keep coverage/risk Pareto fronts and one representative per site set first."""
    if maximum <= 0 or not candidates:
        return []
    if not any(candidate.geometry_support for candidate in candidates):
        return sorted(candidates, key=lambda c: (c.structural_risk, -c.epitope_score, c.order, c.candidate_id))[:maximum]
    groups = {}
    for candidate in sorted(candidates, key=spatial_priority):
        groups.setdefault(site_signature(candidate), []).append(candidate)
    remaining = [group[0] for group in groups.values()]
    representatives = []
    while remaining:
        frontier = [c for c in remaining if not any(
            (other.geometry_support.get("target_mean_coverage_fraction") or 0.0) >= (c.geometry_support.get("target_mean_coverage_fraction") or 0.0)
            and
            other.geometry_support["mean_marginal_area_a2"] >= c.geometry_support["mean_marginal_area_a2"]
            and _initial_risk(other) <= _initial_risk(c)
            and ((other.geometry_support.get("target_mean_coverage_fraction") or 0.0) > (c.geometry_support.get("target_mean_coverage_fraction") or 0.0)
                 or other.geometry_support["mean_marginal_area_a2"] > c.geometry_support["mean_marginal_area_a2"]
                 or _initial_risk(other) < _initial_risk(c))
            for other in remaining)]
        representatives.extend(sorted(frontier, key=spatial_priority))
        ids = {candidate.candidate_id for candidate in frontier}
        remaining = [candidate for candidate in remaining if candidate.candidate_id not in ids]
        if len(representatives) >= maximum:
            break
    output = []
    for offset in range(max((len(groups[site_signature(c)]) for c in representatives), default=0)):
        for candidate in representatives:
            group = groups[site_signature(candidate)]
            if offset < len(group):
                output.append(group[offset])
            if len(output) == maximum:
                return output
    return output


def diverse_sites(positions, maximum, distances, areas, initially_selected=()):
    if areas is None:
        return positions[:maximum]
    selected = list(initially_selected)
    remaining = list(positions)
    output = []
    covered = {radius: {p for site in selected for p, d in distances[site].items() if d <= radius}
               for radius in (9, 12, 15)}
    rank = {position: index for index, position in enumerate(positions)}
    while remaining and len(output) < maximum:
        def key(position):
            gain = sum(areas.get(p, 0.0) for radius in covered
                       for p, d in distances[position].items() if d <= radius and p not in covered[radius])
            return (-gain, rank[position])
        winner = min(remaining, key=key)
        output.append(winner)
        remaining.remove(winner)
        for radius in covered:
            covered[radius].update(p for p, d in distances[winner].items() if d <= radius)
    return output
