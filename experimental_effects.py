"""Pure, auditable calculations for activation readouts and V4 surface attribution.

FCS parsing, gating and QC approval happen upstream. Nothing here treats cell
events as independent experimental replicates or splits combination effects.
"""

from __future__ import annotations

import random
from collections import defaultdict
from dataclasses import dataclass
from heapq import heappop, heappush
from math import exp, isfinite
from statistics import median
from typing import Protocol

ACCEPTED_QC = frozenset({"PASS", "PASS_WITH_OVERRIDE"})


@dataclass(frozen=True)
class WellCount:
    batch_id: str
    stratum_id: str  # fixed antibody, assay, concentration and cell type
    condition: str
    positive: int
    total: int
    qc_status: str = "PASS"

    @property
    def fraction(self) -> float:
        if self.total <= 0 or not 0 <= self.positive <= self.total:
            raise ValueError("Invalid positive/total event counts")
        return self.positive / self.total


@dataclass(frozen=True)
class WellRate:
    """A well-level percentage when the event numerator/denominator is unavailable."""

    batch_id: str
    stratum_id: str
    condition: str
    fraction: float
    qc_status: str = "PASS"

    def __post_init__(self) -> None:
        if not isfinite(self.fraction) or not 0 <= self.fraction <= 1:
            raise ValueError("Well percentage fraction must be finite and between 0 and 1")


class _WellMeasurement(Protocol):
    batch_id: str
    stratum_id: str
    condition: str
    qc_status: str

    @property
    def fraction(self) -> float: ...


@dataclass(frozen=True)
class BatchEffect:
    batch_id: str
    stratum_id: str
    mutant_fraction: float
    reference_fraction: float
    mutant_nc_fraction: float
    reference_nc_fraction: float
    difference_pp: float
    mutant_control_window_pp: float
    reference_control_window_pp: float
    retention: float | None
    decline_pct: float | None
    qc_flags: tuple[str, ...]


def calculate_batch_effect(
    wells: list[_WellMeasurement],
    *,
    mutant: str,
    reference: str,
    mutant_nc: str,
    reference_nc: str | None = None,
    min_reference_window: float = 0.10,
) -> BatchEffect:
    """Use equal-weight well medians within one comparable independent batch."""
    if not 0 <= min_reference_window <= 1:
        raise ValueError("min_reference_window must be between 0 and 1")
    if not wells:
        raise ValueError("No wells supplied")
    keys = {(well.batch_id, well.stratum_id) for well in wells}
    if len(keys) != 1:
        raise ValueError("Wells must share one batch and comparable assay stratum")
    grouped: dict[str, list[WellCount]] = defaultdict(list)
    for well in wells:
        _ = well.fraction  # validate excluded wells as well as accepted wells
        if well.qc_status in ACCEPTED_QC:
            grouped[well.condition].append(well)
    reference_nc = reference_nc or mutant_nc

    def p(condition: str) -> float:
        if not grouped[condition]:
            raise ValueError(f"Missing QC-approved condition: {condition}")
        return float(median(well.fraction for well in grouped[condition]))

    pm, pr, pnm, pnr = p(mutant), p(reference), p(mutant_nc), p(reference_nc)
    window = pr - pnr
    flags = []
    if any(well.qc_status == "PASS_WITH_OVERRIDE" for group in grouped.values() for well in group):
        flags.append("PASS_WITH_OVERRIDE")
    if window <= 0 or window < min_reference_window:
        flags.append("INSUFFICIENT_REFERENCE_WINDOW")
        retention = decline = None
    else:
        retention = (pm - pnm) / window
        decline = 100 * (1 - retention)
        if retention < 0 or retention > 1:
            flags.append("OUT_OF_RANGE_INTERPRETATION")
    batch_id, stratum_id = next(iter(keys))
    return BatchEffect(
        batch_id, stratum_id, pm, pr, pnm, pnr,
        100 * (pm - pr), 100 * (pm - pnm), 100 * (pr - pnr),
        retention, decline, tuple(flags),
    )


def calculate_factorial_interaction(
    parent_antibody: BatchEffect, variant_antibody: BatchEffect
) -> dict:
    """Calculate raw percentage-point difference-in-differences for one batch."""
    if parent_antibody.batch_id != variant_antibody.batch_id:
        raise ValueError("Factorial antibody contrasts must come from the same batch")
    flags = set(parent_antibody.qc_flags) | set(variant_antibody.qc_flags)
    if "INSUFFICIENT_REFERENCE_WINDOW" in parent_antibody.qc_flags:
        flags.add("PARENT_REFERENCE_WINDOW_INVALID")
    if "INSUFFICIENT_REFERENCE_WINDOW" in variant_antibody.qc_flags:
        flags.add("VARIANT_REFERENCE_WINDOW_INVALID")
    parent_effect = parent_antibody.difference_pp
    variant_effect = variant_antibody.difference_pp
    return {
        "batch_id": parent_antibody.batch_id,
        "parent_antigen_effect_pp": parent_effect,
        "variant_antigen_effect_pp": variant_effect,
        "interaction_pp": variant_effect - parent_effect,
        "four_cells": {
            "parent_antibody": {
                "mutant_antigen": parent_antibody.mutant_fraction,
                "wt_antigen": parent_antibody.reference_fraction,
                "mutant_antigen_nc": parent_antibody.mutant_nc_fraction,
                "wt_antigen_nc": parent_antibody.reference_nc_fraction,
            },
            "variant_antibody": {
                "mutant_antigen": variant_antibody.mutant_fraction,
                "wt_antigen": variant_antibody.reference_fraction,
                "mutant_antigen_nc": variant_antibody.mutant_nc_fraction,
                "wt_antigen_nc": variant_antibody.reference_nc_fraction,
            },
        },
        "qc_flags": tuple(sorted(flags)),
    }


def summarize_factorial_interactions(
    rows: list[dict], *, seed: int = 0, bootstrap_iterations: int = 2000
) -> dict:
    """Summarize independent-batch interaction estimates without pooling events."""
    if not rows:
        raise ValueError("No factorial interaction batches supplied")
    batch_ids = [row["batch_id"] for row in rows]
    if len(set(batch_ids)) != len(batch_ids):
        raise ValueError("Factorial interaction rows must have unique independent batches")
    values = [float(row["interaction_pp"]) for row in rows]
    if any(not isfinite(value) for value in values):
        raise ValueError("Factorial interaction values must be finite")
    center = float(median(values))
    return {
        "batch_count": len(rows),
        "median_interaction_pp": center,
        "interaction_interval_pp": _batch_bootstrap_interval(
            values, seed=seed, iterations=bootstrap_iterations
        ),
        "positive_batch_count": sum(value > 0 for value in values),
        "negative_batch_count": sum(value < 0 for value in values),
        "zero_batch_count": sum(value == 0 for value in values),
        "qc_flags": tuple(sorted({flag for row in rows for flag in row.get("qc_flags", ())})),
        "interpretation": "DESCRIPTIVE_COUPLING_CLUE",
    }


def _condition_median(wells: list[_WellMeasurement], condition: str) -> float:
    if not wells:
        raise ValueError(f"Missing {condition} wells")
    if len({well.batch_id for well in wells}) != 1 or len({well.stratum_id for well in wells}) != 1:
        raise ValueError("Condition wells must share one batch and antibody stratum")
    accepted = []
    for well in wells:
        _ = well.fraction
        if well.condition != condition:
            raise ValueError(f"Expected only {condition} wells")
        if well.qc_status in ACCEPTED_QC:
            accepted.append(well.fraction)
    if not accepted:
        raise ValueError(f"Missing QC-approved {condition} wells")
    return float(median(accepted))


def calculate_antibody_variant_batch_effect(
    parent_wells: list[_WellMeasurement],
    variant_wells: list[_WellMeasurement],
    *,
    min_reference_window: float = 0.10,
) -> dict:
    """Compare variant-vs-parent WT activation windows within one batch."""
    if not 0 <= min_reference_window <= 1:
        raise ValueError("min_reference_window must be between 0 and 1")
    if not parent_wells or not variant_wells:
        raise ValueError("Both parent and variant antibody wells are required")
    parent_batch_ids = {well.batch_id for well in parent_wells}
    variant_batch_ids = {well.batch_id for well in variant_wells}
    if len(parent_batch_ids) != 1 or parent_batch_ids != variant_batch_ids:
        raise ValueError("Parent and variant antibody wells must share one independent batch")
    parent_wt = _condition_median(
        [well for well in parent_wells if well.condition == "REFERENCE"], "REFERENCE"
    )
    parent_nc = _condition_median(
        [well for well in parent_wells if well.condition == "REFERENCE_NC"], "REFERENCE_NC"
    )
    variant_wt = _condition_median(
        [well for well in variant_wells if well.condition == "REFERENCE"], "REFERENCE"
    )
    variant_nc = _condition_median(
        [well for well in variant_wells if well.condition == "REFERENCE_NC"], "REFERENCE_NC"
    )
    parent_window = parent_wt - parent_nc
    variant_window = variant_wt - variant_nc
    flags: set[str] = set()
    if any(well.qc_status == "PASS_WITH_OVERRIDE" for well in parent_wells + variant_wells):
        flags.add("PASS_WITH_OVERRIDE")
    parent_valid = parent_window > 0 and parent_window >= min_reference_window
    variant_valid = variant_window > 0 and variant_window >= min_reference_window
    if not parent_valid:
        flags.add("PARENT_REFERENCE_WINDOW_INVALID")
    if not variant_valid:
        flags.add("VARIANT_REFERENCE_WINDOW_INVALID")
    return {
        "batch_id": next(iter(parent_batch_ids)),
        "parent_wt_fraction": parent_wt,
        "parent_wt_nc_fraction": parent_nc,
        "parent_control_window_pp": 100 * parent_window,
        "variant_wt_fraction": variant_wt,
        "variant_wt_nc_fraction": variant_nc,
        "variant_control_window_pp": 100 * variant_window,
        "variant_minus_parent_pp": 100 * (variant_window - parent_window),
        "valid": parent_valid and variant_valid,
        "qc_flags": tuple(sorted(flags)),
        "interpretation": "DESCRIPTIVE_ANTIBODY_EFFECT",
    }


def summarize_antibody_variant_effects(
    rows: list[dict], *, seed: int = 0, bootstrap_iterations: int = 2000
) -> dict:
    """Summarize valid independent-batch WT-window differences."""
    if not rows:
        raise ValueError("No antibody-effect batches supplied")
    batch_ids = [row["batch_id"] for row in rows]
    if len(set(batch_ids)) != len(batch_ids):
        raise ValueError("Antibody-effect rows must have unique independent batches")
    values = [float(row["variant_minus_parent_pp"]) for row in rows if row.get("valid")]
    if any(not isfinite(value) for value in values):
        raise ValueError("Antibody-effect values must be finite")
    center = float(median(values)) if values else None
    return {
        "batch_count": len(rows),
        "valid_batch_count": len(values),
        "median_variant_minus_parent_pp": center,
        "interval_pp": _batch_bootstrap_interval(
            values, seed=seed, iterations=bootstrap_iterations
        ) if values else None,
        "positive_batch_count": sum(value > 0 for value in values),
        "negative_batch_count": sum(value < 0 for value in values),
        "zero_batch_count": sum(value == 0 for value in values),
        "qc_flags": tuple(sorted({flag for row in rows for flag in row.get("qc_flags", ())})),
        "interpretation": "DESCRIPTIVE_ANTIBODY_EFFECT",
    }


@dataclass(frozen=True)
class ProbeBatchEffect:
    probe_id: str
    mutated_positions: tuple[int, ...]
    effect: BatchEffect


@dataclass(frozen=True)
class ProbeSummary:
    probe_id: str
    mutated_positions: tuple[int, ...]
    batch_count: int
    valid_decline_batch_count: int
    difference_pp: float
    decline_pct: float | None
    decline_interval_pct: tuple[float, float] | None
    positive_decline_batch_count: int
    negative_decline_batch_count: int
    zero_decline_batch_count: int
    attribution: str
    sector_id: str | None
    qc_flags: tuple[str, ...]


@dataclass(frozen=True)
class SurfaceProjection:
    probes: tuple[ProbeSummary, ...]
    sector_decline_pct: dict[str, float]
    sector_probe_ids: dict[str, tuple[str, ...]]
    single_site_probe_ids: dict[int, tuple[str, ...]]


def _percentile(values: list[float], q: float) -> float:
    values.sort()
    offset = (len(values) - 1) * q
    lower = int(offset)
    fraction = offset - lower
    return values[lower] * (1 - fraction) + values[min(lower + 1, len(values) - 1)] * fraction


def _batch_bootstrap_interval(
    values: list[float], *, seed: int, iterations: int
) -> tuple[float, float] | None:
    if len(values) < 3:
        return None
    if iterations < 100:
        raise ValueError("Bootstrap requires at least 100 iterations")
    rng = random.Random(seed)
    estimates = [
        float(median(rng.choices(values, k=len(values)))) for _ in range(iterations)
    ]
    return _percentile(estimates.copy(), 0.025), _percentile(estimates, 0.975)


def project_probe_effects(
    rows: list[ProbeBatchEffect],
    sectors: dict[str, set[int]],
    *,
    seed: int = 0,
    bootstrap_iterations: int = 2000,
) -> SurfaceProjection:
    """Assign whole-probe effects to sectors; only true singles get site entries."""
    owner: dict[int, str] = {}
    for sector_id, positions in sectors.items():
        for position in positions:
            if position in owner:
                raise ValueError(f"Overlapping V4 sectors at position {position}")
            owner[position] = sector_id
    grouped: dict[str, list[ProbeBatchEffect]] = defaultdict(list)
    for row in rows:
        grouped[row.probe_id].append(row)
    probes = []
    sector_values: dict[str, list[tuple[str, float]]] = defaultdict(list)
    site_ids: dict[int, list[str]] = defaultdict(list)
    for probe_id, probe_rows in sorted(grouped.items()):
        positions = tuple(sorted(probe_rows[0].mutated_positions))
        if not positions or len(set(positions)) != len(positions):
            raise ValueError(f"Probe {probe_id} has invalid mutation positions")
        if any(tuple(sorted(row.mutated_positions)) != positions for row in probe_rows):
            raise ValueError(f"Probe {probe_id} has inconsistent mutation positions")
        if len({row.effect.batch_id for row in probe_rows}) != len(probe_rows):
            raise ValueError(f"Probe {probe_id} repeats an independent batch")
        if len({row.effect.stratum_id for row in probe_rows}) != 1:
            raise ValueError(f"Probe {probe_id} mixes assay strata")
        declines = [
            row.effect.decline_pct for row in probe_rows if row.effect.decline_pct is not None
        ]
        sectors_touched = {owner.get(position) for position in positions}
        sector_id = None
        if None in sectors_touched:
            attribution = "UNMAPPED_OR_NOT_SURFACE"
        elif len(sectors_touched) > 1:
            attribution = "CROSS_SECTOR_UNASSIGNED"
        else:
            sector_id = next(iter(sectors_touched))
            attribution = "SINGLE_SITE" if len(positions) == 1 else "SAME_SECTOR_COMBINATION"
        decline = float(median(declines)) if declines else None
        if sector_id is not None and decline is not None:
            sector_values[sector_id].append((probe_id, decline))
        if attribution == "SINGLE_SITE":
            site_ids[positions[0]].append(probe_id)
        probes.append(
            ProbeSummary(
                probe_id=probe_id,
                mutated_positions=positions,
                batch_count=len(probe_rows),
                valid_decline_batch_count=len(declines),
                difference_pp=float(median(row.effect.difference_pp for row in probe_rows)),
                decline_pct=decline,
                decline_interval_pct=_batch_bootstrap_interval(
                    declines, seed=seed, iterations=bootstrap_iterations
                ),
                positive_decline_batch_count=sum(value > 0 for value in declines),
                negative_decline_batch_count=sum(value < 0 for value in declines),
                zero_decline_batch_count=sum(value == 0 for value in declines),
                attribution=attribution,
                sector_id=sector_id,
                qc_flags=tuple(sorted({flag for row in probe_rows for flag in row.effect.qc_flags})),
            )
        )
    return SurfaceProjection(
        probes=tuple(probes),
        sector_decline_pct={
            sector_id: float(median(value for _, value in values))
            for sector_id, values in sorted(sector_values.items())
        },
        sector_probe_ids={
            sector_id: tuple(probe_id for probe_id, _ in values)
            for sector_id, values in sorted(sector_values.items())
        },
        single_site_probe_ids={position: tuple(ids) for position, ids in sorted(site_ids.items())},
    )


def visual_fade(
    graph: dict[int, dict[int, float]],
    site: int,
    *,
    bandwidth_a: float = 5.0,
    cutoff_a: float = 12.0,
) -> dict[int, float]:
    """Gaussian display weights along a frozen graph; never use as evidence."""
    if site not in graph:
        raise ValueError("Site is not a mapped surface node")
    if bandwidth_a <= 0 or cutoff_a <= 0:
        raise ValueError("Fade distances must be positive")
    distances = {site: 0.0}
    queue = [(0.0, site)]
    while queue:
        distance, current = heappop(queue)
        if distance != distances[current]:
            continue
        for neighbor, step in graph[current].items():
            if not isfinite(step) or step <= 0:
                raise ValueError("Surface graph needs finite positive edge weights")
            proposed = distance + step
            if proposed > cutoff_a or proposed >= distances.get(neighbor, float("inf")):
                continue
            distances[neighbor] = proposed
            heappush(queue, (proposed, neighbor))
    return {
        position: exp(-distance * distance / (2 * bandwidth_a * bandwidth_a))
        for position, distance in sorted(distances.items())
    }
