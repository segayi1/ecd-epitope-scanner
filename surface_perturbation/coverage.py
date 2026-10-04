"""Area-conserving center coverage with explicit unknowns and nested panels."""
from __future__ import annotations

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix, hstack, vstack, eye
from scipy.sparse.csgraph import connected_components, dijkstra


def scan_support(graph, changed, owner, *, radius=12., minimum_sites=1):
    """Unweighted contact-support scenario, never a binding probability."""
    n = len(changed)
    if n != len(owner) or graph.shape != (n, n):
        raise ValueError("canonical reference dimension mismatch")
    counts = np.zeros(n, dtype=np.int32)
    # One residue may own multiple samples but must count only once.
    for residue in np.unique(owner[changed]):
        sources = np.flatnonzero(changed & (owner == residue))
        distance = dijkstra(graph, directed=False, indices=sources, limit=radius, min_only=True)
        # Accumulate distinct direct-change residue support.
        counts += distance <= radius
    return counts >= minimum_sites


def largest_blind(graph, covered, area):
    ids = np.flatnonzero(~covered & (area > 0))
    if not len(ids):
        return 0.
    _, labels = connected_components(graph[ids][:, ids], directed=False)
    return float(np.bincount(labels, weights=area[ids]).max(initial=0))


def nested_panel(hard, area, risk, graph, ids, maximum=12):
    hard = np.asarray(hard, dtype=bool)
    covered = np.zeros(len(area), bool)
    chosen, remaining = [], set(range(len(ids)))
    while remaining and len(chosen) < maximum:
        def rank(i):
            proposed = covered | hard[i]
            gain = float(area[hard[i] & ~covered].sum())
            # Distinguish panel signatures after coverage/connected blind area.
            # This is an area-based engineering tie-break, not a prediction of
            # antibody response or an arbitrary weighted objective.
            distinction = min((float(area[hard[i] ^ hard[j]].sum()) for j in chosen),default=0.)
            return (-gain, largest_blind(graph, proposed, area), -distinction, float(risk[i]), ids[i])
        winner = min(remaining, key=rank)
        if float(area[hard[winner] & ~covered].sum()) <= 0:
            break
        chosen.append(winner); remaining.remove(winner); covered |= hard[winner]
    return chosen


def budget_curve(hard, extended, order, area, designable, graph):
    hard = np.asarray(hard, bool)
    extended = np.asarray(extended, float)  # nan is unknown, not false.
    area, designable = np.asarray(area), np.asarray(designable, bool)
    denominator = float(area.sum())
    design_denominator = float(area[designable].sum())
    rows = []
    h = np.zeros(len(area), bool); e = np.zeros(len(area), bool); known = np.ones(len(area), bool)
    counts=np.zeros(len(area),int)
    for k in range(1, 13):
        if k <= len(order):
            i = order[k-1]
            h |= hard[i]
            counts += hard[i]
            e |= np.nan_to_num(extended[i], nan=0).astype(bool)
            known &= np.isfinite(extended[i])
        def ratio(mask, design=False):
            denom = design_denominator if design else denominator
            value = float(area[mask & designable].sum()) if design else float(area[mask].sum())
            return value / denom if denom else None
        # Existing positive evidence is known even if another candidate's field is absent.
        resolved_extended = known | e | h
        ext = e | h
        rows.append({"K": k, "actual_constructs": min(k, len(order)),
                     "C_hard_all": ratio(h), "C_hard_designable": ratio(h, True),
                     "C_extended_all": ratio(ext) if resolved_extended.all() else None,
                     "C_extended_designable": ratio(ext, True) if resolved_extended[designable].all() else None,
                     "extended_lower_bound_all": ratio(ext),
                     "extended_unknown_area_a2": float(area[~resolved_extended].sum()),
                     "largest_blind_patch_a2": largest_blind(graph, h, area),
                     "area_weighted_redundancy":float(np.sum(area*np.maximum(counts-1,0))/denominator) if denominator else None,
                     "curve_kind": "nested_verified_panel_not_independent_K_optimum"})
    return rows


def independent_bounds(hard, area, k, seconds=10):
    """Maximum coverage MILP bounds on THIS measured matrix only."""
    hard = csr_matrix(np.asarray(hard, dtype=float).T)
    centers, candidates = hard.shape
    if candidates == 0 or centers == 0:
        return {"K": k, "lower_area_a2": 0., "upper_area_a2": 0., "gap": 0., "status": "EMPTY"}
    objective = np.r_[np.zeros(candidates), -np.asarray(area)]
    constraints = vstack([hstack([-hard, eye(centers)]),
                          csr_matrix(np.r_[np.ones(candidates), np.zeros(centers)][None])]).tocsr()
    result = milp(objective, integrality=np.ones(candidates+centers),
                  bounds=Bounds(0, 1), constraints=LinearConstraint(constraints, -np.inf, np.r_[np.zeros(centers), k]),
                  options={"time_limit": seconds, "mip_rel_gap": .01})
    dual = getattr(result, "mip_dual_bound", None)
    return {"K": k, "lower_area_a2": -float(result.fun) if result.fun is not None else 0.,
            "upper_area_a2": min(float(np.asarray(area).sum()), -float(dual)) if dual is not None and np.isfinite(dual) else None,
            "gap": float(result.mip_gap) if getattr(result, "mip_gap", None) is not None else None,
            "status": int(result.status), "selected_indices": np.flatnonzero(result.x[:candidates] > .5).tolist() if result.x is not None else [],
            "scope": "measured_candidate_pool_only"}
