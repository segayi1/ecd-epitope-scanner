"""Opt-in coverage policy; chemical descriptors are scheduling proxies only."""
from __future__ import annotations

import math
from typing import Any

from ..stability import AA_PROPERTIES, sequons

PROFILE = "coverage_scan_v2_preview"
PROFILE_V3 = "coverage_scan_v3_preview"
PROFILES = {PROFILE, PROFILE_V3}


def coverage_enabled(profile):
    return profile in PROFILES
DESCRIPTOR_VERSION = "sidechain_descriptors_v1"
DONORS = {"R": 3, "K": 1, "N": 1, "Q": 1, "S": 1, "T": 1, "W": 1, "Y": 1}
ACCEPTORS = {"D": 2, "E": 2, "N": 1, "Q": 1, "S": 1, "T": 1, "Y": 1, "M": 1}


def resolve_threshold(profile: str, metadata: dict[str, Any]) -> float:
    if profile not in {"legacy", *PROFILES}:
        raise ValueError("Unknown optimization_profile")
    from ..config import RSASA_THRESHOLD
    nested = metadata.get("joint_v4_ml", {}).get("surface", {}).get("rsasa_threshold")
    upstream = metadata.get("surface_rsasa_threshold")
    supplied = [float(v) for v in (upstream, nested) if v is not None]
    if any(not math.isfinite(v) or not 0 <= v <= 1 for v in supplied):
        raise ValueError("rSASA threshold must be finite and in [0, 1]")
    if len(set(supplied)) > 1:
        raise ValueError("Conflicting explicit rSASA thresholds")
    return supplied[0] if supplied else (0.20 if coverage_enabled(profile) else RSASA_THRESHOLD)


def sequence_reasons(sequence: str, mutations: list[tuple[int, str, str]]) -> list[str]:
    """Check the entire canonical sequence, including synergistic sequon creation."""
    reasons: set[str] = set()
    changed = list(sequence)
    seen: set[int] = set()
    for position, wt, mutant in mutations:
        if position in seen:
            reasons.add("duplicate-position")
        seen.add(position)
        if not 1 <= position <= len(sequence):
            reasons.add("numbering-out-of-range")
            continue
        if sequence[position - 1] != wt:
            reasons.add("WT_RESIDUE_MISMATCH")
        if wt in "CGP":
            reasons.add("native-protected-residue")
        if mutant in "CGP" or mutant not in AA_PROPERTIES or mutant == wt:
            reasons.add("protected-or-invalid-substitution")
        changed[position - 1] = mutant
    before, after = sequons(sequence), sequons("".join(changed))
    if after - before:
        reasons.add("glycosylation-sequon-created")
    if before - after:
        reasons.add("glycosylation-sequon-destroyed")
    return sorted(reasons)


def descriptors(wt: str, mutant: str) -> dict[str, Any]:
    a, b = AA_PROPERTIES[wt], AA_PROPERTIES[mutant]
    unknown = "H" in (wt, mutant)
    return {
        "version": DESCRIPTOR_VERSION,
        "reference": "stability.AA_PROPERTIES; nominal sidechain functional-group counts",
        "charge_delta": None if unknown else b[2] - a[2],
        "volume_delta_a3": b[1] - a[1],
        "hydropathy_delta": b[0] - a[0],
        "donor_delta": None if unknown else DONORS.get(mutant, 0) - DONORS.get(wt, 0),
        "acceptor_delta": None if unknown else ACCEPTORS.get(mutant, 0) - ACCEPTORS.get(wt, 0),
        "aromatic_delta": b[3] - a[3],
        "protonation_status": "UNKNOWN_HISTIDINE" if unknown else "NOMINAL_REFERENCE_STATE",
    }


def vector(d: dict[str, Any]) -> tuple[float, ...]:
    return tuple((d[name] or 0) / scale for name, scale in (
        ("charge_delta", 2), ("volume_delta_a3", 115), ("hydropathy_delta", 9),
        ("donor_delta", 3), ("acceptor_delta", 2), ("aromatic_delta", 1)))


def mechanisms(d: dict[str, Any]) -> set[str]:
    return {name for name in ("charge_delta", "volume_delta_a3", "hydropathy_delta",
                             "donor_delta", "acceptor_delta", "aromatic_delta")
            if d[name] is not None and d[name] != 0}


def select_menu(options: list[Any], maximum: int) -> list[Any]:
    if maximum <= 0 or not options:
        return []
    def energy(item):
        return item.spurs_ddg_kcal_mol if item.spurs_ddg_kcal_mol is not None else math.inf
    def axes(item):
        return (energy(item),) + tuple(-abs(x) for x in vector(item.chemical_descriptors))
    remaining = list(options)
    layer = 0
    while remaining:
        front = [x for x in remaining if not any(
            all(a <= b for a, b in zip(axes(y), axes(x))) and
            any(a < b for a, b in zip(axes(y), axes(x))) for y in remaining if y is not x)]
        for x in front:
            x.pareto_layer = layer
        remaining = [x for x in remaining if all(x is not y for y in front)]
        layer += 1
    selected = [min(options, key=lambda x: (energy(x), x.mutant_aa))]
    selected[0].selection_reason = "STABILITY_REPRESENTATIVE"
    while len(selected) < min(maximum, len(options)):
        used = set().union(*(mechanisms(x.chemical_descriptors) for x in selected))
        remaining = [x for x in options if all(x is not y for y in selected)]
        best_layer = min(x.pareto_layer for x in remaining)
        remaining = [x for x in remaining if x.pareto_layer == best_layer]
        def rank(x):
            v = vector(x.chemical_descriptors)
            distance = min(sum((a-b)**2 for a,b in zip(v, vector(y.chemical_descriptors))) for y in selected)
            return (-len(mechanisms(x.chemical_descriptors) - used), -distance, energy(x), x.mutant_aa)
        winner = min(remaining, key=rank)
        winner.selection_reason = "PARETO_MECHANISM_COMPLEMENT"
        selected.append(winner)
    for x in selected:
        x.menu_selected = True
    return selected
