"""Radius sensitivity of the saved scan kernel; read-only, zero new solves.

The 12 A scan kernel is a scenario assumption, not a calibrated antibody
footprint radius. This audit recomputes the achieved hard coverage of an
already delivered panel under alternative radii so the reported target status
can be read against its own scenario assumption.

It never writes to the case folder, never changes the frozen 12 A contract and
never re-runs FoldX, SPURS or APBS. Only the delivered coverage matrix, its
canonical reference graph and the panel summary are read.
"""
import argparse
from pathlib import Path
import json

import numpy as np
from scipy.sparse import load_npz

from ecd_epitope_scanner.surface_perturbation.coverage import scan_support
from ecd_epitope_scanner.util import write_json

# Radii are reported side by side; 12.0 is the frozen contract value and is
# never written back into any result. 6.0 approximates the direct footprint
# scale, 9.0 and 15.0 bracket the V4.1 core and legacy proxy cuts.
RADI_A = (6.0, 9.0, 12.0, 15.0, 18.0)
CONTRACT_RADIUS_A = 12.0


def load_reference(folder):
    """The reference graph is a scipy sparse file; the arrays a plain npz.

    Mirrors how the engine wrote them: save_npz for the graph and
    savez_compressed for every ndarray field of the canonical reference.
    """
    with np.load(folder / "canonical_reference.npz", allow_pickle=False) as data:
        area, owner = data["area"], data["owner"]
    return load_npz(folder / "reference_graph.npz"), area, owner


def seeds_for(output, candidate_ids):
    """Re-read the actual verified per-site seeds from the delivered evidence."""
    rows = []
    for candidate_id in candidate_ids:
        path = output / "surface_perturbation" / f"{candidate_id}_direct_site_evidence.npz"
        if not path.is_file():
            raise FileNotFoundError(f"missing direct site evidence for {candidate_id}")
        with np.load(path, allow_pickle=False) as data:
            rows.append(data["site_changes"])
    if not rows:
        raise ValueError("no delivered seed masks")
    return np.concatenate(rows, axis=0)


def coverage_at(seeds, owner, graph, area, radius):
    """Fraction of the fixed reliable WT denominator reached at this radius."""
    support = scan_support(graph, seeds.any(axis=0), owner, radius=radius)
    return float(area[support].sum() / area.sum()), support


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="delivered results folder")
    parser.add_argument("--report", type=Path, help="audit output path")
    args = parser.parse_args()
    output = args.output.resolve()

    summary = json.loads((output / "surface_perturbation_preview.json").read_text(encoding="utf-8"))
    candidate_ids = summary["selected_candidate_ids"]
    graph, area, owner = load_reference(output / "surface_perturbation")
    seeds = seeds_for(output, candidate_ids)

    curve = []
    for radius in RADI_A:
        fraction, support = coverage_at(seeds, owner, graph, area, radius)
        curve.append({
            "radius_a": radius,
            "is_frozen_contract": radius == CONTRACT_RADIUS_A,
            "C_screen_lower_bound": fraction,
            "covered_area_a2": float(area[support].sum()),
            "blind_area_a2": float(area[~support].sum()),
        })

    denominator = float(area.sum())
    report = {
        "audit": "scan_kernel_radius_sensitivity",
        "constructs": len(candidate_ids),
        "selected_candidate_ids": candidate_ids,
        "denominator_area_a2": denominator,
        "frozen_contract_radius_a": CONTRACT_RADIUS_A,
        "frozen_contract_coverage": next(r["C_screen_lower_bound"] for r in curve
                                          if r["is_frozen_contract"]),
        "curve": curve,
        "interpretation": (
            "Coverage is recomputed from the delivered seeds under alternative "
            "scan radii. This measures how much the achieved target status "
            "depends on the scenario assumption; it is not a calibrated "
            "footprint radius and not an experimental detection rate."),
        "no_new_scientific_calls": True,
    }
    destination = args.report or output.parent / "scan_kernel_radius_sensitivity.json"
    write_json(destination, report)
    print(destination)


if __name__ == "__main__":
    main()
