"""Replay saved APBS evidence with its preparation, medium and numerical QC.

No solve or preparation waiver. Raw conditional fields are exported separately
from usable evidence; missing prerequisites remain unknown with explicit reasons.
"""
import argparse
from pathlib import Path
import json

import numpy as np
from scipy.sparse import load_npz
from scipy.sparse.csgraph import dijkstra

from ecd_epitope_scanner.surface_perturbation.calibration import calibrate_channel, resolved_changes
from ecd_epitope_scanner.util import write_json, sha256_file


def owner_index(owner, wanted):
    return np.flatnonzero(np.isin(owner, sorted(wanted)))


def surface_distances(graph, owner, wanted):
    if graph.shape != (len(owner), len(owner)):
        raise ValueError('surface graph/reference shape mismatch')
    sites = owner_index(owner, wanted)
    distance = (np.asarray(dijkstra(graph, directed=False, indices=sites, min_only=True)).ravel()
                if len(sites) else np.full(len(owner), np.inf))
    if distance.shape != (len(owner),):
        raise ValueError('per-point distance shape mismatch')
    return distance


def replay_saved(data, report, protocol):
    raw = np.asarray(data['raw_preparation_conditional_offset_deltas'], float)
    if raw.ndim != 3 or raw.shape[:2] != (5, 3):
        raise ValueError('five production fields and three offsets required')
    deltas = raw[:, 1, :]
    numerical = np.asarray(data['numerical_error'], float)
    solvent = np.asarray(data['solvent_valid'], bool)
    null = np.asarray(data['null_blocks'], float)
    if numerical.shape != deltas.shape[1:] or solvent.shape != deltas.shape:
        raise ValueError('saved numerical/medium reference shape mismatch')
    if null.ndim != 2 or null.shape[1:] != deltas.shape[1:]:
        raise ValueError('saved null reference shape mismatch')
    preparation = (report.get('evidence_states', {}).get('APBS_PREPARATION_VALID') is True
                   and report.get('preparation_gate_rejected') is False)
    medium = bool(np.asarray(data['medium_evaluated']).item())
    usable_runs = solvent & np.isfinite(deltas) & preparation & medium
    verified = np.where(usable_runs, deltas, np.nan)
    calibration = calibrate_channel(null, numerical, protocol, 'apbs')
    changed, unknown = resolved_changes(verified, calibration, protocol)
    reasons = {
        'PREPARATION_NOT_VALIDATED': len(numerical) if not preparation else 0,
        'MEDIUM_NOT_EVALUATED': len(numerical) if not medium else 0,
        'MISSING_OR_INVALID_NUMERICAL_ERROR': int((~np.isfinite(numerical) | (numerical < 0)).sum()),
        'INSUFFICIENT_COMMON_SOLVENT_REPEATS': int((~usable_runs.all(axis=0)).sum()),
        'NULL_CALIBRATION_UNRESOLVED': int((~calibration.get('valid', np.zeros(len(numerical), bool))).sum()),
    }
    return deltas, verified, changed, unknown, reasons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="delivered results folder")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError('Use a new report directory; preserve previous reanalysis')
    args.report.parent.mkdir(parents=True, exist_ok=True)
    output = args.output.resolve()
    physical = output / "surface_perturbation"

    summary = json.loads((output / "surface_perturbation_preview.json").read_text(encoding="utf-8"))
    for cid in summary['selected_candidate_ids']:
        if args.report.with_name(f'{cid}_electrostatic_decay_input.npz').exists():
            raise FileExistsError('preserve existing decay inputs; use a new directory')
    candidates = {c["candidate_id"]: c for c in json.loads(
        (output / "panel_candidates.json").read_text(encoding="utf-8"))}
    with np.load(physical / "canonical_reference.npz", allow_pickle=False) as data:
        owner, area = data["owner"], data["area"]
    graph = load_npz(physical / "reference_graph.npz")

    rows = []
    for candidate_id in summary["selected_candidate_ids"]:
        path = physical / f"{candidate_id}_electrostatics.npz"
        report_path = physical / f"{candidate_id}_electrostatics.json"
        original_report = json.loads(report_path.read_text(encoding='utf-8'))
        with np.load(path, allow_pickle=False) as data:
            deltas, verified, changed, unknown, reasons = replay_saved(
                data, original_report, candidate_id + ':apbs:paired-grid-v1')
        if deltas.shape[1] != len(owner):
            raise ValueError('APBS/canonical reference length mismatch')
        finite = np.isfinite(deltas).all(axis=0)
        usable = ~unknown
        correlation = (float(np.corrcoef(deltas[0][finite], deltas[1][finite])[0, 1])
                       if finite.sum() >= 2 and np.all(np.std(deltas[:2, finite], axis=1) > 0) else None)

        positions = sorted({m["residue"]["uniprot_pos"] for m in candidates[candidate_id]["mutations"]})
        distance = surface_distances(graph, owner, positions)

        rows.append({
            "candidate_id": candidate_id,
            "mutated_uniprot_positions": positions,
            "resolved_points": int(changed.sum()),
            "unknown_points": int(unknown.sum()),
            "usable_fraction": float(usable.mean()),
            "scientific_status": 'UNRESOLVED' if unknown.all() else ('PARTIAL' if unknown.any() else 'RESOLVED'),
            "unknown_reason_counts": reasons,
            "raw_statistics_status": 'CONDITIONAL_NOT_VERIFIED_EVIDENCE',
            "input_sha256": sha256_file(path),
            "preparation_report_sha256": sha256_file(report_path),
            "repeat_correlation_run0_run1": correlation,
            "max_abs_delta_phi_kT_e": float(np.max(np.abs(deltas[np.isfinite(deltas)]))) if np.isfinite(deltas).any() else None,
            "median_abs_delta_phi_kT_e": float(np.median(np.abs(deltas[np.isfinite(deltas)]))) if np.isfinite(deltas).any() else None,
        })
        np.savez_compressed(
            args.report.with_name(f"{candidate_id}_electrostatic_decay_input.npz"),
            delta_phi=deltas.astype(np.float32),
            verified_delta_phi=verified.astype(np.float32),
            geodesic_distance_a=distance.astype(np.float32),
            owner=owner, area=area, resolved=changed, unknown=unknown)

    write_json(args.report, {
        "audit": "electrostatic_reanalysis_after_preparation_fix",
        "constructs": len(rows),
        "rows": rows,
        "new_apbs_calls": 0,
        "new_pdb2pqr_calls": 0,
        "source_sha256": sha256_file(Path(__file__)),
        "reference_sha256": sha256_file(physical/'canonical_reference.npz'),
        "graph_sha256": sha256_file(physical/'reference_graph.npz'),
        "interpretation": (
            "Original preparation verdict, common solvent masks and measured numerical "
            "error retained. Raw finite fields do not imply scientific validity. "
            "Reason counts overlap. Legacy APBS null is an empirical audit, not the "
            "independent production/train/holdout formal contract. No hard coverage change."),
    })
    print(args.report)


if __name__ == "__main__":
    main()
