#!/usr/bin/env python
"""Create read-only, explicitly status-aware charts for the P0-P3 six-case report.

No model calls are made. Incomplete cases are emitted as PENDING rows and are
never coerced to zero. The input report is the single current audit snapshot.
"""
from __future__ import annotations

import argparse
import csv
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _number(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _sizes_json(value: Any) -> list[int] | None:
    if isinstance(value, list):
        try:
            return [int(x) for x in value]
        except (TypeError, ValueError):
            return None
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
            return _sizes_json(parsed)
        except (json.JSONDecodeError, TypeError):
            return None
    if isinstance(value, dict) and isinstance(value.get("sizes"), list):
        return _sizes_json(value["sizes"])
    return None


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _audit_counts(task: dict[str, Any]) -> dict[str, Any]:
    audit = task.get("audit") or {}
    counts = audit.get("counts") or {}
    return counts


def build_sources(report: dict[str, Any], old_rows: list[dict[str, str]]) -> dict[str, list[dict[str, Any]]]:
    old_by_protein = {r.get("binder", ""): r for r in old_rows}
    regions: list[dict[str, Any]] = []
    budgets: list[dict[str, Any]] = []
    funnels: list[dict[str, Any]] = []
    biases: list[dict[str, Any]] = []
    for task in report.get("tasks", []):
        protein = task.get("protein") or task.get("input_file") or task.get("new_job_id", "unknown")
        status = str(task.get("status", "PENDING"))
        old = old_by_protein.get(protein, {})
        old_part = _sizes_json(old.get("delivery_partition"))
        new_part = _sizes_json(task.get("new_region_sizes_surface_nodes"))
        regions.extend([
            {"protein": protein, "series": "原始流程/旧分区", "status": old.get("job_status", "UNKNOWN"),
             "region_index": i + 1, "surface_residues": n}
            for i, n in enumerate(old_part or [])
        ])
        if new_part is None:
            regions.append({"protein": protein, "series": "新流程/新分区", "status": "PENDING", "region_index": "", "surface_residues": ""})
        else:
            regions.extend({"protein": protein, "series": "新流程/新分区", "status": status,
                            "region_index": i + 1, "surface_residues": n}
                           for i, n in enumerate(new_part))

        seq = task.get("new_sequential_foldx_budget") or {}
        allocation = task.get("new_foldx_allocation") or {}
        budgets.append({
            "protein": protein, "status": status,
            "old_unique_combinations": _number(old.get("foldx_candidates")),
            "old_attempt_slots_60x5": (_number(old.get("foldx_candidates")) * 5) if _number(old.get("foldx_candidates")) is not None else None,
            "new_unique_combinations": _number(seq.get("unique_combinations")),
            "new_initial_slots": _number(seq.get("initial_slots")),
            "new_supplemental_slots_attempted": _number(seq.get("supplemental_slots")),
            "new_total_slots_attempted": _number(seq.get("total_slots")),
            "new_max_unique": 80, "new_max_pair_slots": 300,
            "requested_local": _number(allocation.get("requested_local", allocation.get("local_requested"))),
            "actual_local": _number(allocation.get("actual_local")),
            "actual_global": _number(allocation.get("actual_global")),
            "local_exploration": _number(allocation.get("local_exploration")),
        })

        counts = _audit_counts(task)
        # Each number has its own denominator: initial schedule is unique
        # combinations, candidate audit stages use audited candidates, and the
        # main panel is a subset of formal candidates.
        seq_budget = task.get("new_sequential_foldx_budget") or {}
        initial_slots = _number(seq_budget.get("initial_slots"))
        unique_combos = _number(seq_budget.get("unique_combinations"))
        pool_den = _number(counts.get("candidate_denominator"))
        vals = [
            ("eligible candidate pool", counts.get("candidate_denominator"), pool_den, "candidates"),
            ("scheduled initial 3 slots", initial_slots / 3 if initial_slots is not None else None, unique_combos, "unique combinations"),
            ("confirmed 5 slots", counts.get("foldx_complete_five"), pool_den, "audited candidate pool"),
            ("formal eligible", counts.get("formal_eligible"), pool_den, "audited candidate pool"),
            ("main panel", counts.get("selected_panel_count"), _number(counts.get("formal_eligible")), "formal eligible"),
        ]
        for stage, val, denominator, denominator_type in vals:
            funnels.append({"protein": protein, "status": status, "stage": stage,
                            "count": _number(val), "denominator": _number(denominator),
                            "denominator_type": denominator_type})

        diversity = task.get("audit", {}).get("selection_diversity_and_chemistry") or {}
        for cohort in ("candidate_pool", "foldx_scheduled", "formal_eligible", "main_panel"):
            item = diversity.get(cohort) or {}
            biases.append({"protein": protein, "status": status, "cohort": cohort,
                           "chemical_identities": _number(item.get("unique_chemical_identities")),
                           "unique_site_sets": _number(item.get("unique_site_sets")),
                           "mean_hydropathy_shift": _number((item.get("hydropathy_shift") or {}).get("mean_sum_per_candidate")),
                           "mean_charge_shift": _number((item.get("charge_shift") or {}).get("mean_sum_per_candidate")),
                           "mean_aromatic_shift": _number((item.get("aromatic_shift") or {}).get("mean_sum_per_candidate")),
                           "native_interface_fraction": _number(item.get("native_interface_mutation_fraction")),
                           "interface_enrichment": _number(item.get("interface_enrichment_ratio_vs_accessible_surface"))})
    return {"regions": regions, "budgets": budgets, "funnels": funnels, "biases": biases}


def pending_tasks(report):
    return [t for t in report.get("tasks", []) if t.get("status") != "COMPLETE_ARTIFACTS_PRESENT"
            or (t.get("audit") or {}).get("status") not in {"PASS", "AUDIT_COMPLETE"}]


def render(sources: dict[str, list[dict[str, Any]]], outdir: Path, report: dict[str, Any]) -> dict[str, Any]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import font_manager
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is unavailable in the selected Python environment; no installation attempted") from exc

    outdir.mkdir(parents=True, exist_ok=True)
    pending = pending_tasks(report)
    label = "PARTIAL — 有任务未完成" if pending else "FINAL — 六例均完成"
    # Readable Chinese defaults; the output remains useful if a CJK font is not installed.
    plt.rcParams["axes.unicode_minus"] = False
    cjk_font = next((p for p in (Path("C:/Windows/Fonts/msyh.ttc"), Path("C:/Windows/Fonts/simhei.ttf"),
                      Path("C:/Windows/Fonts/simsun.ttc")) if p.is_file()), None)
    if cjk_font:
        font_manager.fontManager.addfont(str(cjk_font))
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=str(cjk_font)).get_name()
    else:
        # Keep generation usable on headless systems without CJK fonts; suppress
        # repetitive glyph warnings and retain Latin labels/data CSVs.
        warnings.filterwarnings("ignore", message="Glyph .* missing from current font")

    # Figure 1: old vs new per-protein partition sizes, unavailable stays missing.
    fig, ax = plt.subplots(figsize=(11, 5.8))
    proteins = list(dict.fromkeys(r["protein"] for r in sources["regions"]))
    for y, protein in enumerate(proteins):
        rows = [r for r in sources["regions"] if r["protein"] == protein]
        series_styles = {"原始流程/旧分区": ("#5B8FF9", -0.12), "新流程/新分区": ("#E8684A", 0.12)}
        for series, (color, dy) in series_styles.items():
            points = [r for r in rows if r["series"] == series and _number(r.get("surface_residues")) is not None]
            if points:
                xs = [float(r["surface_residues"]) for r in points]
                ys = [y + dy + ((i % 5) - 2) * .018 for i in range(len(points))]
                ax.scatter(xs, ys, label=series if y == 0 else None, color=color, s=42, alpha=.85)
            elif series == "新流程/新分区":
                ax.text(0.02, y + dy, "待完成", transform=ax.get_yaxis_transform(), color=color, va="center", fontsize=8)
    ax.axvspan(6, 8, color="#5DAA68", alpha=.12, label="新小蛋白偏好 6–8")
    ax.set_yticks(range(len(proteins)), proteins)
    ax.set_xlabel("每区表面残基数（aa）")
    ax.set_title(f"每例分区大小对照｜{label}")
    ax.grid(axis="x", alpha=.2); ax.legend(loc="best", fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / "01_region_sizes_comparison.png", dpi=170); plt.close(fig)

    # Figure 2: old fixed 60x5 and new sequential protocol, not conflating slots with unique candidates.
    fig, ax = plt.subplots(figsize=(11, 5.8))
    x = range(len(sources["budgets"]))
    names = [r["protein"] for r in sources["budgets"]]
    oldslots = [_number(r["old_attempt_slots_60x5"]) for r in sources["budgets"]]
    initial = [_number(r["new_initial_slots"]) for r in sources["budgets"]]
    supplemental = [_number(r["new_supplemental_slots_attempted"]) for r in sources["budgets"]]
    ax.bar([i - .18 for i in x], [v if v is not None else 0 for v in oldslots], width=.36, label="旧：60 unique × 5 slots", color="#5B8FF9")
    # Unknown bars remain visually absent and are labeled pending.
    ax.bar([i + .18 for i in x], [v if v is not None else 0 for v in initial], width=.36, label="新：已尝试初始槽（3/组）", color="#61DDAA")
    ax.bar([i + .18 for i in x], [v if v is not None else 0 for v in supplemental], width=.36,
           bottom=[v if v is not None else 0 for v in initial], label="新：已尝试补算槽", color="#F6BD16")
    for i, r in enumerate(sources["budgets"]):
        if r["new_total_slots_attempted"] is None:
            ax.text(i + .18, 6, "待完成", ha="center", va="bottom", fontsize=8, rotation=90)
    ax.axhline(300, color="#E8684A", ls="--", lw=1, label="新协议最多300配对槽")
    ax.set_xticks(list(x), names); ax.set_ylabel("FoldX配对运行槽（不是模型调用数）")
    ax.set_title(f"FoldX预算：unique组合与重复槽分开｜{label}")
    ax.grid(axis="y", alpha=.2); ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout(); fig.savefig(outdir / "02_foldx_budget_comparison.png", dpi=170); plt.close(fig)

    # Figure 3: explicitly unequal funnel counts and descriptive-only chemistry/interface trends.
    fig, (axf, axb) = plt.subplots(2, 1, figsize=(12, 9), gridspec_kw={"height_ratios": [1.2, 1]})
    stages = ["eligible candidate pool", "scheduled initial 3 slots", "confirmed 5 slots", "formal eligible", "main panel"]
    for protein in proteins:
        rows = [r for r in sources["funnels"] if r["protein"] == protein]
        bystage = {r["stage"]: r for r in rows}
        yy = [i for i, st in enumerate(stages)]
        xx = [_number(bystage.get(st, {}).get("count")) for st in stages]
        present = [(xxi, yyi) for xxi, yyi in zip(xx, yy) if xxi is not None]
        if present:
            axf.plot([p[0] for p in present], [p[1] for p in present], marker="o", label=protein, alpha=.8)
            for xv, yv in present:
                stage_row = bystage[stages[yv]]
                den = _number(stage_row.get("denominator"))
                suffix = f"/{int(den)}" if den is not None else ""
                axf.annotate(f"{int(xv)}{suffix}", (xv, yv), xytext=(4, 3), textcoords="offset points", fontsize=7)
        if any(v is None for v in xx):
            axf.text(.995, .03 + proteins.index(protein) * .04, f"{protein}: 待完成", transform=axf.transAxes, ha="right", fontsize=7)
    axf.set_yticks(range(len(stages)), stages); axf.invert_yaxis(); axf.set_xlabel("候选数（各阶段按任务各自候选全集计数）")
    axf.set_title("候选漏斗（点标注=数量/该阶段分母；分母定义不同，不作阶段转化率）")
    axf.grid(axis="x", alpha=.2)
    if proteins:
        axf.legend(ncol=3, fontsize=7, loc="best")

    # Plot available pool-to-panel interface fraction and chemistry shift, never infer missing as zero.
    cohort_order = ["candidate_pool", "foldx_scheduled", "formal_eligible", "main_panel"]
    axchem = axb.twinx()
    for protein in proteins:
        rows = [r for r in sources["biases"] if r["protein"] == protein]
        bycohort = {r["cohort"]: r for r in rows}
        available = [(i, _number(bycohort.get(c, {}).get("native_interface_fraction"))) for i, c in enumerate(cohort_order)]
        avail = [(i, v) for i, v in available if v is not None]
        if avail:
            axb.plot([i for i, _ in avail], [v for _, v in avail], marker="o", label=protein, alpha=.8)
        for field, axis_label, linestyle in (("mean_hydropathy_shift", "Δhydropathy", "--"),
                                             ("mean_charge_shift", "Δcharge", ":"),
                                             ("mean_aromatic_shift", "Δaromatic", "-.")):
            vals = [(i, _number(bycohort.get(cohort, {}).get(field))) for i, cohort in enumerate(cohort_order)]
            vals = [(i, v) for i, v in vals if v is not None]
            if vals:
                axchem.plot([i for i, _ in vals], [v for _, v in vals], linestyle=linestyle,
                            marker=".", alpha=.42, linewidth=.8,
                            label=f"{protein} {axis_label}")
    axb.set_xticks(range(4), ["pool", "FoldX scheduled", "formal", "main panel"])
    axb.set_ylabel("native-interface mutation fraction")
    axchem.set_ylabel("每候选平均化学变化（描述值）")
    axb.set_title("描述性界面/化学诊断：界面比例与平均疏水/电荷/芳香变化（缺失不补零；不代表亲和力/因果）")
    axb.set_ylim(bottom=0); axb.grid(alpha=.2)
    h1, l1 = axb.get_legend_handles_labels(); h2, l2 = axchem.get_legend_handles_labels()
    axb.legend(h1 + h2, l1 + l2, ncol=3, fontsize=6, loc="best")
    fig.suptitle(f"候选路径与偏倚诊断｜{label}", y=.995)
    fig.tight_layout(); fig.savefig(outdir / "03_candidate_funnel_and_bias.png", dpi=170); plt.close(fig)

    return {"status": "PARTIAL_PENDING" if pending else "FINAL", "pending_jobs": [t.get("new_job_id") for t in pending],
            "task_count": len(report.get("tasks", [])), "matplotlib_version": matplotlib.__version__,
            "figures": ["01_region_sizes_comparison.png", "02_foldx_budget_comparison.png", "03_candidate_funnel_and_bias.png"]}


def generate(report_path: Path, baseline_csv: Path, outdir: Path) -> dict[str, Any]:
    report = _read_json(report_path)
    sources = build_sources(report, _read_csv(baseline_csv))
    fields = {
        "regions": ["protein", "series", "status", "region_index", "surface_residues"],
        "budgets": ["protein", "status", "old_unique_combinations", "old_attempt_slots_60x5", "new_unique_combinations", "new_initial_slots", "new_supplemental_slots_attempted", "new_total_slots_attempted", "new_max_unique", "new_max_pair_slots", "requested_local", "actual_local", "actual_global", "local_exploration"],
        "funnels": ["protein", "status", "stage", "count", "denominator", "denominator_type"],
        "biases": ["protein", "status", "cohort", "chemical_identities", "unique_site_sets", "mean_hydropathy_shift", "mean_charge_shift", "mean_aromatic_shift", "native_interface_fraction", "interface_enrichment"],
    }
    outdir.mkdir(parents=True, exist_ok=True)
    for name, rows in sources.items():
        _write_csv(outdir / f"figure_{name}_source.csv", rows, fields[name])
    receipt = render(sources, outdir, report)
    receipt.update({"generated_at_utc": datetime.now(timezone.utc).isoformat(),
                    "input_report": str(report_path), "baseline_csv": str(baseline_csv),
        "scope_note": "Read-only summaries; pending values are blank/annotated, not zero. Interface/chemistry are descriptive diagnostics.",
        "supplemental_references": {
            "external_benchmark_summary": "workspace_data/p0_p3_upgrade_20261016/benchmark/benchmark_data_review_zh.md",
            "external_benchmark_acceptance": "workspace_data/p0_p3_upgrade_20261016/benchmark/benchmark_acceptance.json",
            "environment_lock": "workspace_data/p0_p3_upgrade_20261016/environment_lock_20261017T013955Z.json",
            "environment_lock_scope": "Supplemental lock captured after the first model call; does not replace a pre-start environment capture."
        }})
    (outdir / "figure_generation_receipt.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--baseline-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(generate(args.report, args.baseline_csv, args.outdir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
