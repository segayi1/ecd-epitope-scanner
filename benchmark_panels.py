"""Evaluate Stage-1 panels against held-out structural contact labels.

This is a benchmark of direct mutation-site hits, not an assay of binding loss.
The manifest format and interpretation limits are documented in
docs/joint_v4_validation_protocol.md.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Any

ENGINES = ("surface_v2", "joint_v3", "joint_v4_ml")


def _load_run(path: Path, engine: str, budget: int) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    analysis = json.loads(path.read_text(encoding="utf-8"))
    observed_engine = analysis.get("parameters", {}).get("patch_engine")
    if observed_engine != engine:
        raise ValueError(f"{path}: expected {engine}, found {observed_engine}")
    if engine == "joint_v4_ml" and analysis.get("joint_v4_ml", {}).get(
        "recommendation_status"
    ) != "FORMAL_RECOMMENDATION":
        return analysis, []
    selected = [
        design
        for design in analysis.get("designs", [])
        if design.get("design_role") == "STAGE1_PROBE"
        and str(design.get("status", "")).startswith("RECOMMENDED")
    ]
    return analysis, selected[:budget]


def _run_identity(analysis: dict[str, Any]) -> tuple[Any, ...]:
    protein = analysis.get("protein", {})
    identity = (
        protein.get("sequence_length"),
        protein.get("ecd_start"),
        protein.get("ecd_end"),
        protein.get("ecd_sequence"),
        analysis.get("source", {}).get("structure_sha256"),
        analysis.get("request", {}).get("assay_format"),
    )
    if any(value is None or value == "" for value in identity):
        raise ValueError("benchmark run lacks sequence, ECD range, structure hash or assay")
    return identity


def _quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    index = round(probability * (len(ordered) - 1))
    return ordered[index]


def _cluster_interval(
    rows: list[dict[str, Any]], metric: str, draws: int, seed: int
) -> dict[str, float | int]:
    by_antigen: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_antigen[row["antigen_id"]].append(row)
    antigens = sorted(by_antigen)
    observed = mean(float(row[metric]) for row in rows) if rows else 0.0
    if len(antigens) < 2 or draws < 1:
        return {"mean": observed, "ci_low": observed, "ci_high": observed, "antigen_count": len(antigens)}
    rng = random.Random(seed)
    sampled = []
    for _ in range(draws):
        draw_rows = [
            row for antigen in (rng.choice(antigens) for _ in antigens) for row in by_antigen[antigen]
        ]
        sampled.append(mean(float(row[metric]) for row in draw_rows))
    return {
        "mean": observed,
        "ci_low": _quantile(sampled, 0.025),
        "ci_high": _quantile(sampled, 0.975),
        "antigen_count": len(antigens),
    }


def evaluate_manifest(
    manifest_path: Path, *, split: str = "frozen_test", bootstrap_draws: int = 1000, seed: int = 0
) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    budget = int(manifest.get("construct_budget", 7))
    if budget < 1:
        raise ValueError("construct_budget must be positive")
    cases = manifest.get("cases", [])
    if not isinstance(cases, list) or not cases:
        raise ValueError("manifest requires at least one case")
    family_split: dict[str, str] = {}
    seen_cases: set[tuple[str, str]] = set()
    rows: list[dict[str, Any]] = []
    for case in cases:
        antigen = str(case["antigen_id"])
        antibody = str(case["antibody_id"])
        if (antigen, antibody) in seen_cases:
            raise ValueError(f"duplicate benchmark case: {antigen}/{antibody}")
        seen_cases.add((antigen, antibody))
        family = str(case.get("antigen_family") or antigen)
        case_split = str(case["split"])
        if family in family_split and family_split[family] != case_split:
            raise ValueError(f"antigen family {family} occurs in multiple splits")
        family_split[family] = case_split
        if case_split != split:
            continue
        contacts = {int(position) for position in case["true_contact_positions"]}
        if not contacts or min(contacts) < 1:
            raise ValueError(f"{antigen}/{antibody}: true_contact_positions must be positive")
        runs = case["runs"]
        if set(runs) != set(ENGINES):
            raise ValueError(f"{antigen}/{antibody}: all three engine runs are required")
        identities = []
        for engine in ENGINES:
            run_path = (manifest_path.parent / runs[engine]).resolve()
            analysis, selected = _load_run(run_path, engine, budget)
            identities.append(_run_identity(analysis))
            positions = {
                int(mutation["position"])
                for design in selected
                for mutation in design.get("mutations", [])
            }
            contact_hits = positions & contacts
            rows.append({
                "antigen_id": antigen,
                "antibody_id": antibody,
                "antigen_family": family,
                "engine": engine,
                "construct_count": len(selected),
                "budget_shortfall": max(0, budget - len(selected)),
                "mutated_position_count": len(positions),
                "direct_contact_hit_count": len(contact_hits),
                "hit_at_least_one": int(len(contact_hits) >= 1),
                "hit_at_least_two": int(len(contact_hits) >= 2),
                "run_path": str(run_path),
            })
        if len(set(identities)) != 1:
            raise ValueError(f"{antigen}/{antibody}: runs do not share input, structure and assay")
        ecd_start, ecd_end = identities[0][1:3]
        if min(contacts) < ecd_start or max(contacts) > ecd_end:
            raise ValueError(f"{antigen}/{antibody}: contacts fall outside the ECD range")
    if not rows:
        raise ValueError(f"no cases in split {split!r}")
    summary = {
        engine: {
            metric: _cluster_interval(
                [row for row in rows if row["engine"] == engine], metric, bootstrap_draws, seed
            )
            for metric in ("hit_at_least_one", "hit_at_least_two", "construct_count", "budget_shortfall")
        }
        for engine in ENGINES
    }
    paired_differences = {}
    keyed = {(row["antigen_id"], row["antibody_id"], row["engine"]): row for row in rows}
    for baseline in ENGINES[:-1]:
        differences = []
        for case in cases:
            if case["split"] != split:
                continue
            antigen, antibody = str(case["antigen_id"]), str(case["antibody_id"])
            current = keyed[(antigen, antibody, "joint_v4_ml")]
            control = keyed[(antigen, antibody, baseline)]
            differences.append({
                "antigen_id": antigen,
                "hit_at_least_one": current["hit_at_least_one"] - control["hit_at_least_one"],
                "hit_at_least_two": current["hit_at_least_two"] - control["hit_at_least_two"],
            })
        paired_differences[baseline] = {
            metric: _cluster_interval(differences, metric, bootstrap_draws, seed)
            for metric in ("hit_at_least_one", "hit_at_least_two")
        }
    return {
        "benchmark_version": "1.0",
        "split": split,
        "construct_budget": budget,
        "case_count": len(rows) // len(ENGINES),
        "endpoint": "direct mutated-site overlap with known structural antibody contacts",
        "interpretation": "A contact hit is not evidence of binding loss or preserved expression.",
        "case_rows": rows,
        "summary": summary,
        "paired_v4_minus_baseline": paired_differences,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--split", default="frozen_test")
    parser.add_argument("--bootstrap-draws", type=int, default=1000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = evaluate_manifest(
        args.manifest, split=args.split, bootstrap_draws=args.bootstrap_draws
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Evaluated {result['case_count']} cases at budget {result['construct_budget']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
