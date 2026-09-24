from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .models import ProbeCandidate, ProbePanel, SurfaceResidueUnit


@dataclass
class _Selection:
    candidate: ProbeCandidate
    mean_gain: float
    robust_gain: float


def _candidate_risk(candidate: ProbeCandidate) -> float:
    if candidate.stability_value is not None:
        return max(0.0, float(candidate.stability_value))
    return float(candidate.structural_risk)


def _area_map(units: list[SurfaceResidueUnit]) -> dict[int, float]:
    return {unit.position: unit.mapped_area_a2 for unit in units if unit.designable}


def _selection_scenarios(candidates: list[ProbeCandidate], requested: list[str] | None) -> tuple[str, ...]:
    available = set.intersection(
        *(set(item.coverage_units_by_scenario) for item in candidates)
    ) if candidates else set()
    ordered = [item for item in (requested or []) if item in available and item != "direct_site"]
    if not ordered:
        ordered = [
            item for item in (
                "geodesic_tight", "geodesic_medium", "geodesic_wide",
                "footprint_hit_1", "footprint_hit_2",
            ) if item in available
        ]
    return tuple(ordered or (["direct_site"] if "direct_site" in available else []))


def _robust_gain(
    candidate: ProbeCandidate,
    covered: dict[str, set[int]],
    areas: dict[int, float],
    scenarios: tuple[str, ...],
) -> tuple[float, float]:
    gains = []
    for scenario in scenarios:
        positions = set(candidate.coverage_units_by_scenario.get(scenario, [])) - covered[scenario]
        gains.append(sum(areas.get(position, 0.0) for position in positions))
    robust = min(gains) if gains else 0.0
    mean = sum(gains) / max(1, len(gains))
    return robust, mean - _candidate_risk(candidate) * max(1.0, mean) * 0.05


def _milp_select(
    candidates: list[ProbeCandidate],
    units: list[SurfaceResidueUnit],
    budget: int,
    timeout_seconds: float,
    mip_relative_gap: float,
    scenarios: tuple[str, ...],
) -> tuple[list[ProbeCandidate], str, float | None] | None:
    """Two-pass lexicographic robust/mean coverage MILP."""
    try:
        from scipy.optimize import Bounds, LinearConstraint, milp
        from scipy.sparse import lil_matrix, vstack
    except Exception:
        return None
    designable_units = [unit for unit in units if unit.designable]
    if not candidates or not designable_units or not scenarios:
        return None
    n_candidates = len(candidates)
    n_units = len(designable_units)
    n_scenarios = len(scenarios)
    y_offset = n_candidates
    t_index = n_candidates + n_scenarios * n_units
    variable_count = t_index + 1
    unit_index = {unit.position: index for index, unit in enumerate(designable_units)}

    rows = 1 + n_scenarios * n_units + n_scenarios
    matrix = lil_matrix((rows, variable_count), dtype=float)
    lower = np.full(rows, -np.inf, dtype=float)
    upper = np.zeros(rows, dtype=float)
    matrix[0, :n_candidates] = 1.0
    upper[0] = min(budget, n_candidates)
    row = 1
    for scenario_index, scenario in enumerate(scenarios):
        for unit in designable_units:
            y_index = y_offset + scenario_index * n_units + unit_index[unit.position]
            matrix[row, y_index] = 1.0
            for candidate_index, candidate in enumerate(candidates):
                if unit.position in candidate.coverage_units_by_scenario.get(scenario, []):
                    matrix[row, candidate_index] = -1.0
            row += 1
        matrix[row, t_index] = 1.0
        for unit in designable_units:
            y_index = y_offset + scenario_index * n_units + unit_index[unit.position]
            matrix[row, y_index] = -unit.mapped_area_a2
        row += 1

    # Do not select multiple candidates that mutate exactly the same sites.
    signatures: dict[tuple[int, ...], list[int]] = {}
    for index, candidate in enumerate(candidates):
        signatures.setdefault(tuple(sorted(item.position for item in candidate.mutations)), []).append(index)
    exclusion_groups = [values for values in signatures.values() if len(values) > 1]
    site_sets = [{item.position for item in candidate.mutations} for candidate in candidates]
    for left in range(n_candidates):
        for right in range(left + 1, n_candidates):
            overlap = len(site_sets[left] & site_sets[right]) / max(1, min(len(site_sets[left]), len(site_sets[right])))
            if overlap >= 0.75 and [left, right] not in exclusion_groups:
                exclusion_groups.append([left, right])
    if exclusion_groups:
        extra = lil_matrix((len(exclusion_groups), variable_count), dtype=float)
        for group_index, values in enumerate(exclusion_groups):
            for candidate_index in values:
                extra[group_index, candidate_index] = 1.0
        matrix = vstack((matrix.tocsr(), extra.tocsr()))
        lower = np.r_[lower, np.full(len(exclusion_groups), -np.inf)]
        upper = np.r_[upper, np.ones(len(exclusion_groups))]
    else:
        matrix = matrix.tocsr()

    integrality = np.r_[np.ones(n_candidates + n_scenarios * n_units), 0.0]
    lb = np.zeros(variable_count)
    ub = np.ones(variable_count)
    ub[t_index] = sum(unit.mapped_area_a2 for unit in designable_units)
    bounds = Bounds(lb, ub)
    options = {"time_limit": float(timeout_seconds), "mip_rel_gap": float(mip_relative_gap)}

    first_objective = np.zeros(variable_count)
    first_objective[t_index] = -1.0
    try:
        first = milp(
            c=first_objective,
            integrality=integrality,
            bounds=bounds,
            constraints=LinearConstraint(matrix, lower, upper),
            options=options,
        )
    except Exception:
        return None
    if first.x is None or not np.isfinite(first.x).all():
        return None

    robust_optimum = float(first.x[t_index])
    robust_row = lil_matrix((1, variable_count), dtype=float)
    robust_row[0, t_index] = 1.0
    second_matrix = vstack((matrix, robust_row.tocsr()))
    second_lower = np.r_[lower, robust_optimum - 1e-6]
    second_upper = np.r_[upper, np.inf]
    second_objective = np.zeros(variable_count)
    for index, candidate in enumerate(candidates):
        second_objective[index] = _candidate_risk(candidate) * 1e-3
    for scenario_index in range(n_scenarios):
        for unit in designable_units:
            y_index = y_offset + scenario_index * n_units + unit_index[unit.position]
            second_objective[y_index] = -unit.mapped_area_a2 / n_scenarios
    try:
        second = milp(
            c=second_objective,
            integrality=integrality,
            bounds=bounds,
            constraints=LinearConstraint(second_matrix, second_lower, second_upper),
            options=options,
        )
    except Exception:
        second = first
    solution = second if second.x is not None and np.isfinite(second.x).all() else first
    selected = [candidate for index, candidate in enumerate(candidates) if solution.x[index] >= 0.5]
    selected.sort(key=lambda item: item.candidate_id)
    if not selected:
        return None
    status = "MILP_OPTIMAL" if bool(solution.success) else "MILP_FEASIBLE"
    gap = getattr(solution, "mip_gap", None)
    return selected, status, float(gap) if gap is not None else None


def _panel_from_selected(
    selected: list[ProbeCandidate],
    units: list[SurfaceResidueUnit],
    status: str,
    optimality_gap: float | None,
    selection_scenarios: tuple[str, ...],
) -> ProbePanel:
    areas = _area_map(units)
    report_scenarios = tuple(dict.fromkeys(
        ["direct_site", *selection_scenarios]
        + [name for candidate in selected for name in candidate.coverage_units_by_scenario]
    ))
    covered = {scenario: set() for scenario in report_scenarios}
    for candidate in selected:
        for scenario in report_scenarios:
            covered[scenario].update(candidate.coverage_units_by_scenario.get(scenario, []))
    total = max(1.0, sum(areas.values()))
    achieved = {
        scenario: round(sum(areas.get(position, 0.0) for position in covered[scenario]) / total, 4)
        for scenario in report_scenarios
    }
    robust = min((achieved[item] for item in selection_scenarios), default=0.0)
    objective_vector = [
        round(robust * total, 4),
        round(sum(achieved[item] for item in selection_scenarios) / max(1, len(selection_scenarios)), 4),
        round(sum(_candidate_risk(item) for item in selected), 4),
        float(sum(len(item.mutations) for item in selected)),
    ]
    uncovered = [
        {"component_id": unit.component_id, "area_a2": round(unit.mapped_area_a2, 3), "position": unit.position, "reason": "not_covered_in_all_selection_scenarios"}
        for unit in units
        if unit.designable and any(unit.position not in covered[scenario] for scenario in selection_scenarios)
    ]
    uncovered.extend(
        {
            "component_id": unit.component_id,
            "area_a2": round(unit.mapped_area_a2, 3),
            "position": unit.position,
            "reason": "protected_surface",
            "exclusion_reasons": unit.exclusion_reasons,
        }
        for unit in units if not unit.designable
    )
    return ProbePanel(
        panel_id="PANEL_V3_001",
        selected_candidate_ids=[item.candidate_id for item in selected],
        objective_vector=objective_vector,
        achieved_coverage=achieved,
        uncovered_components=uncovered,
        solver_status=status,
        optimality_gap=optimality_gap,
    )


def _greedy_select(
    candidates: list[ProbeCandidate],
    areas: dict[int, float],
    budget: int,
    scenarios: tuple[str, ...],
    excluded_ids: set[str] | None = None,
) -> list[ProbeCandidate]:
    selected: list[ProbeCandidate] = []
    covered = {scenario: set() for scenario in scenarios}
    remaining = [item for item in candidates if item.candidate_id not in (excluded_ids or set())]
    used_site_sets: list[set[int]] = []
    while remaining and len(selected) < budget:
        ranked: list[_Selection] = []
        for candidate in remaining:
            site_set = {item.position for item in candidate.mutations}
            if any(
                len(site_set & used) / max(1, min(len(site_set), len(used))) >= 0.75
                for used in used_site_sets
            ):
                continue
            robust, mean = _robust_gain(candidate, covered, areas, scenarios)
            ranked.append(_Selection(candidate, mean, robust))
        if not ranked:
            break
        ranked.sort(
            key=lambda item: (
                item.robust_gain,
                item.mean_gain,
                -_candidate_risk(item.candidate),
                item.candidate.candidate_id,
            ),
            reverse=True,
        )
        best = ranked[0]
        if best.mean_gain <= 0 and selected:
            break
        selected.append(best.candidate)
        used_site_sets.append({item.position for item in best.candidate.mutations})
        remaining = [item for item in remaining if item.candidate_id != best.candidate.candidate_id]
        for scenario in covered:
            covered[scenario].update(best.candidate.coverage_units_by_scenario.get(scenario, []))
    return selected


def select_probe_panel(
    candidates: list[ProbeCandidate],
    units: list[SurfaceResidueUnit],
    budget: int,
    timeout_seconds: float = 60.0,
    mip_relative_gap: float = 0.01,
    coverage_scenarios: list[str] | None = None,
) -> tuple[ProbePanel, list[ProbeCandidate]]:
    areas = _area_map(units)
    budget = max(1, int(budget))
    scenarios = _selection_scenarios(candidates, coverage_scenarios)
    milp_result = _milp_select(
        candidates, units, budget, timeout_seconds, mip_relative_gap, scenarios
    )
    if milp_result is not None:
        selected, status, gap = milp_result
    else:
        selected = _greedy_select(candidates, areas, budget, scenarios)
        status, gap = "HEURISTIC_SOLUTION", None
    panel = _panel_from_selected(selected, units, status, gap, scenarios)

    alternatives: list[dict] = []
    seen = {tuple(panel.selected_candidate_ids)}
    for excluded in panel.selected_candidate_ids[:3]:
        alternative = _greedy_select(candidates, areas, budget, scenarios, {excluded})
        ids = tuple(item.candidate_id for item in alternative)
        if not ids or ids in seen:
            continue
        seen.add(ids)
        alt_panel = _panel_from_selected(alternative, units, "HEURISTIC_ALTERNATIVE", None, scenarios)
        alt_id = f"PANEL_V3_ALT_{len(alternatives) + 1:02d}"
        alternatives.append({
            "panel_id": alt_id,
            "candidate_ids": list(ids),
            "achieved_coverage": alt_panel.achieved_coverage,
            "objective_vector": alt_panel.objective_vector,
        })
    panel.alternative_panels = alternatives
    panel.alternative_panel_ids = [item["panel_id"] for item in alternatives]
    return panel, selected
