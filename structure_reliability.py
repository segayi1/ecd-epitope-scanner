"""Conservative, auditable reliability routing; evidence scores are not probabilities.

Experimental B factors are never interpreted as pLDDT. PAE is accepted only
with a validated token mapping. Missing evidence cannot trigger exclusion.
"""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np

from .util import sha256_file, write_json

DEFAULTS = {
    "enabled": True,
    "structured_plddt": 70.0,
    "low_plddt": 50.0,
    "core_pae_threshold": 15.0,
    "contact_radius": 8.0,
    "low_contacts": 2,
    "minimum_idr_length": 15,
    "embedded_loop_max_length": 14,
    "smoothing_window": 7,
    "transition_length": 3,
    "sequence_threshold": 0.5,
    "minimum_low_fraction": 0.8,
    "boundary_residues": 2,
    "surface_mode": "strict",
    "include_conditional": False,
    "idr_window": 12,
    "idr_stride": 5,
    "weights": {
        "plddt": 0.30,
        "pae": 0.25,
        "sequence": 0.25,
        "packing": 0.15,
        "terminal": 0.05,
    },
}


def import_iupred(path, sequence):
    """Read complete standard IUPred output (index, AA, score[, anchor])."""
    if path is None:
        return {}, {"status": "unavailable", "method": "iupred2a", "version": None}
    path = Path(path)
    scores = {}
    headers = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            headers.append(line)
            continue
        columns = line.split()
        if len(columns) < 3:
            raise ValueError("Invalid IUPred row: expected position, amino acid, score")
        position, aa, value = int(columns[0]), columns[1], float(columns[2])
        if (
            position in scores
            or not 1 <= position <= len(sequence)
            or sequence[position - 1] != aa
        ):
            raise ValueError(
                "IUPred positions must uniquely match the full canonical sequence"
            )
        if not np.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("IUPred score must be finite and between zero and one")
        scores[position] = value
    if set(scores) != set(range(1, len(sequence) + 1)):
        raise ValueError("IUPred output must cover the full canonical sequence")
    return scores, {
        "status": "imported",
        "method": "iupred2a",
        "version": None,
        "headers": headers,
        "sha256": sha256_file(path),
        "note": "Version is only established by the imported header; no predictor was run.",
    }


def _pae_mapping(confidence, records, chain_length, single_chain):
    if not isinstance(confidence, dict):
        return None, {}, "unavailable"
    raw = confidence.get("pae", confidence.get("predicted_aligned_error"))
    if raw is None:
        return None, {}, "unavailable"
    try:
        matrix = np.asarray(raw, dtype=float)
    except (ValueError, TypeError):
        return None, {}, "invalid_matrix"
    if (
        matrix.ndim != 2
        or matrix.shape[0] != matrix.shape[1]
        or not np.isfinite(matrix).all()
        or (matrix < 0).any()
    ):
        return None, {}, "invalid_matrix"
    chains, resids = confidence.get("token_chain_ids"), confidence.get("token_res_ids")
    mapping = {}
    if (
        isinstance(chains, list)
        and isinstance(resids, list)
        and len(chains) == len(resids) == len(matrix)
    ):
        lookup = {}
        for i, key in enumerate(zip(map(str, chains), map(str, resids))):
            lookup.setdefault(key, []).append(i)
        for r in records:
            keys = [(r.chain_id, str(r.auth_seq_id))]
            if r.label_chain_id and r.label_seq_id is not None:
                keys = [(r.label_chain_id, str(r.label_seq_id))]
            hits = lookup.get(keys[0], [])
            if len(hits) == 1 and not r.insertion_code:
                mapping[r.uniprot_pos] = hits[0]
        if len(set(mapping.values())) != len(mapping):
            return None, {}, "ambiguous_token_mapping"
        return matrix, mapping, "explicit_tokens"
    if single_chain and chain_length == len(matrix):
        mapping = {
            r.uniprot_pos: int(r.token_index)
            for r in records
            if r.token_index is not None and 0 <= r.token_index < len(matrix)
        }
        return matrix, mapping, "validated_single_chain"
    return None, {}, "unmapped_pae_ignored"


def _runs(rows, key):
    output = []
    for row in rows:
        if (
            output
            and row["residue_id"] == output[-1][-1]["residue_id"] + 1
            and key(row) == key(output[-1][-1])
        ):
            output[-1].append(row)
        else:
            output.append([row])
    return output


def classify_reliability(
    records,
    protein,
    *,
    predicted=False,
    confidence=None,
    chain_length=None,
    single_chain=False,
    parameters=None,
    sequence_scores=None,
    sequence_provenance=None,
):
    cfg = DEFAULTS | (parameters or {})
    if cfg["surface_mode"] not in {"strict", "occlusion_sensitivity"}:
        raise ValueError("Unknown reliability surface_mode")
    if cfg["smoothing_window"] < 1 or cfg["smoothing_window"] % 2 == 0:
        raise ValueError("smoothing_window must be positive and odd")
    if not 1 <= cfg["idr_stride"] <= cfg["idr_window"] or cfg["minimum_idr_length"] < 1:
        raise ValueError("Invalid reliability segment/window configuration")
    by_pos = {r.uniprot_pos: r for r in records}
    scores = sequence_scores or {}
    matrix, token, pae_status = _pae_mapping(
        confidence if predicted else None, records, chain_length, single_chain
    )
    core = [
        r.uniprot_pos
        for r in records
        if predicted and r.plddt is not None and r.plddt >= cfg["structured_plddt"]
    ]
    rows = []
    for pos in range(protein.ecd_start, protein.ecd_end + 1):
        r = by_pos.get(pos)
        valid = r is not None and r.ca is not None
        plddt = r.plddt if valid and predicted else None
        contacts = (
            sum(
                abs(other.uniprot_pos - pos) > 2
                and np.linalg.norm(np.asarray(r.ca) - other.ca) <= cfg["contact_radius"]
                for other in records
                if other.ca is not None
            )
            if valid
            else 0
        )

        def pae_to(positions):
            if matrix is None or pos not in token:
                return None
            ids = [token[p] for p in positions if p in token and p != pos]
            return (
                float(np.mean((matrix[token[pos], ids] + matrix[ids, token[pos]]) / 2))
                if ids
                else None
            )

        # Prefer own-domain reference, to avoid calling whole mobile domains IDRs.
        domain_core = (
            [p for p in core if p != pos and by_pos[p].domain == r.domain]
            if valid
            else []
        )
        ref = (
            domain_core
            if len(domain_core) >= 5
            else [p for p in core if abs(p - pos) > 7]
        )
        local_pae = pae_to(
            range(max(protein.ecd_start, pos - 3), min(protein.ecd_end, pos + 3) + 1)
        )
        core_pae = pae_to(ref)
        terminal = (
            min(pos - protein.ecd_start, protein.ecd_end - pos)
            < cfg["minimum_idr_length"]
        )
        evidence = {
            "plddt": (
                (
                    1.0
                    if plddt < cfg["low_plddt"]
                    else 0.5 if plddt < cfg["structured_plddt"] else 0.0
                )
                if plddt is not None
                else None
            ),
            "pae": (
                float(core_pae >= cfg["core_pae_threshold"])
                if core_pae is not None
                else None
            ),
            "sequence": scores.get(pos),
            "packing": float(contacts <= cfg["low_contacts"]) if valid else None,
            "terminal": float(terminal),
        }
        state = (
            "UNRESOLVED"
            if not valid
            else (
                "STRUCTURED"
                if not predicted
                or (plddt is not None and plddt >= cfg["structured_plddt"])
                else "UNCERTAIN"
            )
        )
        rows.append(
            dict(
                residue_id=pos,
                plddt=plddt,
                mean_local_pae=local_pae,
                mean_core_pae=core_pae,
                secondary_structure=r.secondary_structure if valid else None,
                contact_count=int(contacts),
                disorder_sequence_score=scores.get(pos),
                disorder_evidence_score=round(
                    sum(
                        cfg["weights"][k] * v
                        for k, v in evidence.items()
                        if v is not None
                    ),
                    4,
                ),
                available_evidence_weight=round(
                    sum(
                        cfg["weights"][k] for k, v in evidence.items() if v is not None
                    ),
                    4,
                ),
                state=state,
                terminal=terminal,
                quality_flags=[],
                evidence=evidence,
            )
        )
    for run in _runs(rows, lambda x: x["state"] != "UNRESOLVED"):
        half = cfg["smoothing_window"] // 2
        for i, row in enumerate(run):
            row["smoothed_evidence_score"] = (
                float(
                    np.median(
                        [
                            x["disorder_evidence_score"]
                            for x in run[max(0, i - half) : i + half + 1]
                        ]
                    )
                )
                if row["state"] != "UNRESOLVED"
                else None
            )
    if cfg["enabled"] and predicted:
        # Smooth only contiguous low-confidence candidates; never bridge missing residues.
        for run in _runs(rows, lambda x: x["state"] == "UNCERTAIN"):
            if run[0]["state"] != "UNCERTAIN":
                continue
            start, end = run[0]["residue_id"], run[-1]["residue_id"]
            left = by_pos.get(start - 1)
            right = by_pos.get(end + 1)
            embedded = (
                len(run) <= cfg["embedded_loop_max_length"]
                and left
                and right
                and left.domain == right.domain
                and all(
                    x.plddt is not None and x.plddt >= cfg["structured_plddt"]
                    for x in (left, right)
                )
            )
            if embedded:
                for row in run:
                    row["state"] = "FLEXIBLE_STRUCTURED"
                    row["quality_flags"].append("EMBEDDED_LOW_CONFIDENCE_LOOP")
                continue
            low = [
                x["plddt"] is not None and x["plddt"] < cfg["low_plddt"] for x in run
            ]
            packing = np.mean([x["contact_count"] <= cfg["low_contacts"] for x in run])
            sequence_support = (
                np.mean(
                    [
                        x["disorder_sequence_score"] is not None
                        and x["disorder_sequence_score"] >= cfg["sequence_threshold"]
                        for x in run
                    ]
                )
                >= 0.8
            )
            pae_support = (
                np.mean(
                    [
                        x["mean_core_pae"] is not None
                        and x["mean_core_pae"] >= cfg["core_pae_threshold"]
                        for x in run
                    ]
                )
                >= 0.8
            )
            sse_conflict = any(x["secondary_structure"] in {"H", "E"} for x in run)
            supported = (
                len(run) >= cfg["minimum_idr_length"]
                and np.mean(low) >= cfg["minimum_low_fraction"]
                and packing >= 0.8
                and (sequence_support or pae_support)
                and not sse_conflict
            )
            if supported:
                for row in run:
                    row["state"] = "DISORDERED"
                    row["quality_flags"].append(
                        "MULTI_EVIDENCE_LOW_RELIABILITY_SEGMENT"
                    )
                # Retain a flagged transition at a structured attachment.
                if (
                    cfg["transition_length"] > 0
                    and left
                    and left.plddt is not None
                    and left.plddt >= cfg["structured_plddt"]
                ):
                    for row in run[: cfg["transition_length"]]:
                        row["state"] = "UNCERTAIN"
                        row["quality_flags"].append("TRANSITION")
                if (
                    cfg["transition_length"] > 0
                    and right
                    and right.plddt is not None
                    and right.plddt >= cfg["structured_plddt"]
                ):
                    for row in run[-cfg["transition_length"] :]:
                        row["state"] = "UNCERTAIN"
                        row["quality_flags"].append("TRANSITION")
        # Conditional folding requires independent sequence disorder support;
        # high inter-domain PAE alone is insufficient.
        for row in rows:
            if (
                row["state"] == "STRUCTURED"
                and row["secondary_structure"] == "H"
                and row["disorder_sequence_score"] is not None
                and row["disorder_sequence_score"] >= cfg["sequence_threshold"]
                and row["mean_core_pae"] is not None
                and row["mean_core_pae"] >= cfg["core_pae_threshold"]
                and row["contact_count"] <= cfg["low_contacts"]
            ):
                row["state"] = "CONDITIONALLY_FOLDED"
                row["quality_flags"].append("CONDITIONAL_FOLD_HYPOTHESIS")
    # Coordinate eligibility is separate from a biological disorder diagnosis.
    # AF3 can place a poorly anchored terminal chain near a folded domain even
    # when there is not enough evidence to call that chain an IDR. Preserve its
    # state and sequence screening, but do not use that placement as a footprint.
    uncertain_geometry = set()
    if cfg["enabled"] and predicted:
        for run in _runs(rows, lambda r: r["state"] == "UNCERTAIN"):
            if run[0]["state"] != "UNCERTAIN":
                continue
            terminal_run = (
                run[0]["residue_id"] == protein.ecd_start
                or run[-1]["residue_id"] == protein.ecd_end
            )
            supported = (
                terminal_run
                and all(r["plddt"] is not None and r["plddt"] < cfg["structured_plddt"] for r in run)
                and np.mean([r["contact_count"] <= cfg["low_contacts"] for r in run]) >= 0.8
                and any(r["mean_core_pae"] is not None and r["mean_core_pae"] >= cfg["core_pae_threshold"] for r in run)
            )
            if supported:
                for row in run:
                    uncertain_geometry.add(row["residue_id"])
                    row["quality_flags"].append("TERMINAL_COORDINATES_UNRELIABLE_FOR_3D")
    excluded = set()
    for row in rows:
        row["use_for_surface_mesh"] = row["state"] not in {
            "DISORDERED",
            "UNRESOLVED",
        } and (row["state"] != "CONDITIONALLY_FOLDED" or cfg["include_conditional"])
        row["use_for_3d_footprint"] = row["use_for_surface_mesh"]
        row["use_for_epitope_screening"] = True
        if not row["use_for_surface_mesh"]:
            excluded.add(row["residue_id"])
    geometry_excluded = excluded | uncertain_geometry
    boundaries = {
        row["residue_id"]
        for row in rows
        if row["use_for_surface_mesh"]
        and row["residue_id"] not in geometry_excluded
        and any(
            abs(row["residue_id"] - p) <= cfg["boundary_residues"] for p in geometry_excluded
        )
    }
    for row in rows:
        row["boundary_cap_flag"] = row["residue_id"] in boundaries
        row["use_for_3d_footprint"] = (
            row["use_for_surface_mesh"] and not row["boundary_cap_flag"]
            and row["residue_id"] not in uncertain_geometry
        )
        if row["boundary_cap_flag"]:
            row["quality_flags"].append("ARTIFICIAL_BOUNDARY_POSSIBLE")
    segments = []
    candidates = []
    for number, run in enumerate(_runs(rows, lambda r: r["state"]), 1):
        start, end = run[0]["residue_id"], run[-1]["residue_id"]
        state = run[0]["state"]

        def avg(key):
            values = [r[key] for r in run if r[key] is not None]
            return float(np.mean(values)) if values else None

        segment = dict(
            segment_id=f"SEG{number:03d}",
            start=start,
            end=end,
            length=len(run),
            state=state,
            mean_plddt=avg("plddt"),
            mean_core_pae=avg("mean_core_pae"),
            mean_disorder_score=avg("disorder_sequence_score"),
            use_for_surface_mesh=run[0]["use_for_surface_mesh"],
            use_for_idr_screening=state in {"DISORDERED", "CONDITIONALLY_FOLDED"},
            evidence=sorted({flag for r in run for flag in r["quality_flags"]}),
        )
        segments.append(segment)
        if segment["use_for_idr_screening"]:
            width = min(cfg["idr_window"], len(run))
            starts = sorted(
                set(range(start, end - width + 2, cfg["idr_stride"]))
                | {end - width + 1}
            )
            for s in starts:
                e = s + width - 1
                positions = list(range(s, e + 1))
                sequon_positions = {
                    q
                    for i in range(len(protein.sequence) - 2)
                    if protein.sequence[i] == "N"
                    and protein.sequence[i + 1] != "P"
                    and protein.sequence[i + 2] in "ST"
                    for q in range(i + 1, i + 4)
                }
                protected = [
                    p
                    for p in positions
                    if protein.sequence[p - 1] in "CGP"
                    or p in sequon_positions
                    or any(
                        (
                            p in {f.start, f.end}
                            if f.type.lower() in {"disulfide bond", "cross-link"}
                            else f.contains(p)
                        )
                        and f.type.lower()
                        in {
                            "glycosylation",
                            "disulfide bond",
                            "modified residue",
                            "cross-link",
                        }
                        for f in protein.features
                    )
                    or (
                        p in by_pos
                        and (
                            by_pos[p].glycosylation_sequon
                            or by_pos[p].ptm
                            or by_pos[p].disulfide
                        )
                    )
                ]
                candidates.append(
                    dict(
                        candidate_id=f'{segment["segment_id"]}_{s}_{e}',
                        candidate_type=(
                            "CONDITIONAL_FOLD"
                            if state == "CONDITIONALLY_FOLDED"
                            else "IDR_LINEAR"
                        ),
                        start=s,
                        end=e,
                        sequence=protein.sequence[s - 1 : e],
                        sasa=None,
                        protected_positions=protected,
                        review_positions=[p for p in positions if p not in protected],
                        run_foldx=False,
                        status="REVIEW_REQUIRED",
                        note="Sequence windows only; not validated epitopes or approved mutations.",
                    )
                )
    qc = dict(
        classification_version="reliability_rules_v2",
        parameters=cfg,
        enabled=cfg["enabled"],
        predicted=predicted,
        pae_status=pae_status,
        sequence_predictor=sequence_provenance or {"status": "unavailable"},
        state_counts=dict(Counter(r["state"] for r in rows)),
        excluded_positions=sorted(excluded),
        footprint_excluded_positions=sorted(geometry_excluded),
        uncertain_terminal_geometry_positions=sorted(uncertain_geometry),
        boundary_positions=sorted(boundaries),
        idr_candidate_count=len(candidates),
        smoothing="centered median evidence diagnostic; classification uses segment consensus, not score thresholds",
        score_interpretation="uncalibrated evidence sum; missing evidence contributes no support",
    )
    return dict(residues=rows, segments=segments, idr_candidates=candidates, qc=qc)


def export_reliability(result, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for name, rows in [
        ("residue_reliability.csv", result["residues"]),
        ("structural_segments.csv", result["segments"]),
        (
            "surface_mask.csv",
            [
                {
                    k: r[k]
                    for k in (
                        "residue_id",
                        "state",
                        "use_for_surface_mesh",
                        "use_for_3d_footprint",
                        "boundary_cap_flag",
                    )
                }
                for r in result["residues"]
            ],
        ),
        ("idr_candidates.csv", result["idr_candidates"]),
    ]:
        fields = (
            list(rows[0])
            if rows
            else [
                "candidate_id",
                "candidate_type",
                "start",
                "end",
                "sequence",
                "status",
            ]
        )
        with (output / name).open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {
                        k: json.dumps(v) if isinstance(v, (list, dict)) else v
                        for k, v in row.items()
                    }
                )
    write_json(output / "reliability_qc.json", result["qc"])
