"""Versioned panel budget resolution shared by requests and the optimizer."""
from __future__ import annotations


def resolve_budget(*, minimum=None, target=None, maximum=None, adaptive=None, legacy_max=None, policy=None, maximum_limit=8):
    policy = policy or {}
    for name, value in (("minimum", minimum), ("target", target), ("maximum", maximum), ("legacy_max", legacy_max)):
        if value is not None and (isinstance(value, bool) or not isinstance(value, int)):
            raise ValueError(f"{name} construct count must be an integer")
    explicit = any(value is not None for value in (minimum, target, maximum, adaptive))
    if legacy_max is not None and not explicit:
        cap = min(int(legacy_max), int(policy.get("stage1_maximum_constructs", legacy_max)))
        result = dict(minimum_constructs=min(4, cap), target_constructs=cap,
                      maximum_constructs=cap, adaptive_budget=False)
    else:
        cap = int(maximum if maximum is not None else policy.get("stage1_maximum_constructs", 8))
        aim = int(target if target is not None else min(7, policy.get("stage1_target_constructs", 7), cap))
        lower = int(minimum if minimum is not None else min(4, policy.get("stage1_minimum_constructs", 4), aim))
        result = dict(minimum_constructs=lower, target_constructs=aim,
                      maximum_constructs=cap,
                      adaptive_budget=bool(adaptive if adaptive is not None else policy.get("adaptive_budget", not bool(policy))))
    if maximum_limit not in {8,12}: raise ValueError("Unknown budget ceiling")
    if not 1 <= result["minimum_constructs"] <= result["target_constructs"] <= result["maximum_constructs"] <= maximum_limit:
        raise ValueError(f"Panel budget must satisfy 1 <= minimum <= target <= maximum <= {maximum_limit}")
    return result


def optimizer_budget(profile, budget):
    """V1/V2 retain their original selection semantics, V3 opts into stopping."""
    if profile and int(profile.get("version", 1)) < 3:
        if budget["adaptive_budget"]:
            raise ValueError("Adaptive selection requires assay profile V3; V1/V2 retain fixed-budget semantics")
        return None
    return budget
