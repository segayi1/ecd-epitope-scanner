"""Panel replacement audit, distinct from higher-order perturbation gain."""
from __future__ import annotations

import numpy as np

from .coverage import largest_blind


def compression_audit(panel_ids, parent_ids, triple_id, hard, area, graph, baseline_mask):
    area = np.asarray(area, float)
    if not np.isfinite(area).all() or area.sum() <= 0:
        raise ValueError("Invalid frozen reliable surface denominator")
    reasons = []
    if len(set(parent_ids)) != 2 or not set(parent_ids) <= set(panel_ids):
        reasons.append("PARENTS_NOT_IN_PANEL")
    after_ids = [cid for cid in panel_ids if cid not in parent_ids] + [triple_id]
    if len(after_ids) != len(set(after_ids)) or not 4 <= len(after_ids) <= 12:
        reasons.append("CONSTRUCT_COUNT_INVALID")
    before = np.any([hard[cid] for cid in panel_ids], axis=0)
    after = np.any([hard[cid] for cid in after_ids], axis=0)
    denominator = float(area.sum())
    fraction = float(area[after].sum() / denominator)
    loss = float((area[np.asarray(baseline_mask, bool)].sum() - area[after].sum()) / denominator)
    old_blind, new_blind = largest_blind(graph, before, area), largest_blind(graph, after, area)
    if fraction < .8:
        reasons.append("COVERAGE_BELOW_80")
    if loss > .01 + 1e-12:
        reasons.append("CUMULATIVE_LOSS_EXCEEDS_ONE_POINT")
    if new_blind > old_blind + 1e-12:
        reasons.append("LARGEST_BLIND_AREA_INCREASED")
    return {"accepted": not reasons, "reasons": reasons, "before_ids": list(panel_ids),
        "after_ids": after_ids, "hard_coverage_after": fraction,
        "cumulative_coverage_loss": loss, "added_area_a2": float(area[after & ~before].sum()),
        "lost_area_a2": float(area[before & ~after].sum()),
        "largest_blind_before_a2": old_blind, "largest_blind_after_a2": new_blind,
        "construct_count_change": len(after_ids) - len(panel_ids),
        "qualification": "PANEL_REPLACEMENT_ONLY_REQUIRES_INDEPENDENT_TRIPLE_EVIDENCE"}


def merge_parents(left, right):
    """Only AB+AC with exactly the same concrete A substitution yields ABC."""
    if left.order != 2 or right.order != 2:
        return None
    a = {o.residue.uniprot_pos: o for o in left.mutations}
    b = {o.residue.uniprot_pos: o for o in right.mutations}
    shared = a.keys() & b.keys()
    if len(shared) != 1 or len(a.keys() | b.keys()) != 3:
        return None
    p = next(iter(shared))
    if a[p].notation != b[p].notation:
        return None
    return [o for _, o in sorted({**a, **b}.items())]


def expansion_audit(triple_id, double_panel_ids, hard, area):
    """Panel blind-spot gain is additional to the unchanged lower-order gate."""
    area=np.asarray(area,float)
    baseline=np.any([hard[cid] for cid in double_panel_ids],axis=0) if double_panel_ids else np.zeros(len(area),bool)
    added=float(area[np.asarray(hard[triple_id],bool)&~baseline].sum())
    return {'candidate_id':triple_id,'baseline_double_panel_ids':list(double_panel_ids),
        'actual_panel_added_area_a2':added,'accepted':added>0,
        'reasons':[] if added>0 else ['NO_ACTUAL_PANEL_BLINDSPOT_GAIN'],
        'qualification':'ADDITIONAL_PANEL_GAIN_GATE_NOT_SUPERADDITIVITY'}
