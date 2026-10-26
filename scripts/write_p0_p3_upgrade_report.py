"""Build an evidence-linked Chinese report for the six fresh P0-P3 V4 reruns.

The generator is safe to rerun while jobs are pending: incomplete tasks remain
PENDING and are never treated as failures or silently filled from old results.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import sys
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))


CHECKLIST = [
    ("P0", "字段与状态口径", "主代表/正式候选/替代数分开；真实表面门来自运行metadata；请求区数、有效预算、候选配额、FoldX组合与重复槽分开；响应近邻与精确重复分开；原始搜索与派生重分类资源分开。", "代码与展示回归通过；逐任务资格、预算、ZIP和离线状态在下方独立列出。", "不同历史schema继续显示原语义；在线和离线投影应保持一致。"),
    ("P1", "候选证据与预算公平", "通用/机制候选共用当前区归属、序列、描述符、暴露、几何及敏感性证据；未知化学轴不补零；三点按本区增益重算；探索按区轮转，其余名额优先不同位点集的合格双/三点，允许短缺回流但不降门；保存调度用途/未调度原因；分别统计化学身份与位点集多样性。", "调度策略和字段已接入；逐任务审计报告候选漏斗、未调度与多样性证据。", "本批诊断只描述选择结果，不因面板差异拟合校正。"),
    ("P1", "模型差与重复诊断", "按蛋白/阶数/SPURS窗口/重复层级报告有符号及绝对差、重复值/中位数/范围/SD、3到5次门状态变化、缺失/未调度/高风险/证据不完整的独立分母。", "独立审计器 scripts/audit_p0_p3_batch.py 已实现并通过合成正反例；真实分层结果逐任务纳入下方审计。", "重复SD是计算采样描述，不是实验置信区间。"),
    ("P1", "外部实测小基准", "FireProtDB单点、PTMUL多点；映射/条件/单位/符号/重复/同源/训练重叠审查；最多60例，每阶最多20；不按预测值预筛；单变体FoldX五次，独立记账。", "保留score-blind候选清单、作者编号与结构hash审查、runner账本及独立模型配对验证；实际分层统计从benchmark_results.json读取。", "训练集独立性或实验条件无法核实时只能作描述性小样本；不据此自动调阈值或宣称binder通用校准。"),
    ("P2", "序贯FoldX与重复证据", "80个唯一组合请求配额为局部60/全局20，局部强扰动探索最多16；池短缺可在不放宽资格的前提下回流，80总组仍硬上限。每组先3次，最多30组追加2次，总计最多300配对槽；主面板和区主代表正式输出须五次。", "新协议已写入任务manifest/配对slot contract；独立审计逐任务验证回流、slot完整性、80/300/30及正式主面板/区代表五槽。", "候选供给不足按实际值报告，不放宽门槛凑数。"),
    ("P2", "严格缓存与阶段剖析", "身份应绑定模型/源码/结构/链/序列/变体/环境；读取前后验证身份和内容；FoldX复用要有配对和输入证据；记录锁等待、模型加载、解析/编码/解码和端到端。", "缓存与FoldX复用身份缺口会拒绝命中；任务级execution trace、worker runtime和cache receipt进入逐例审计。", "时长是批处理累计，不是CPU核时或实验精度；身份不完整结果不计缓存命中。"),
    ("P2", "SPURS多点编码复用", "同一multi模型/WT只编码一次，双/三点复用解码；不跨单点模型；候选值、顺序、mask、混批与显存一致性回归。", "真实GPU编码复用数值回归已验证最大差值0；新任务记录实际调用/缓存诊断。", "收益受设备与缓存状态影响；不据单个benchmark推定常驻GPU服务收益。"),
    ("P2", "全表面与界面对照", "正式默认仍全表面；界面对照固定候选池/资格/预算，只变目标，不触发模型；比较候选交并、界面比例、直接面积、9/12/15Å代理覆盖、分区覆盖和缺口；描述化学富集不作为惩罚。", "已有同池只读对照实现；逐任务审计报告候选池→FoldX→formal→主面板的界面突变富集、疏水/电荷/芳香变化。", "界面面积和覆盖是描述性代理，不代表结合能或亲和力。"),
    ("P3", "源码、时间与几何复现", "逐任务源文件快照与hash、配置/环境/commit状态；追加式阶段日志；逐残基几何owner采样、贡献、入图和拒绝原因。", "新任务保存source snapshot、execution trace和geometry owner diagnostics；逐任务身份/归档/导出哈希纳入独立审计。", "AIB193运行模块源码哈希依赖同进程导入链推断，未直接记录runtime bytecode digest；见逐例限制。"),
    ("P3", "完整离线网站", "每任务ZIP带根index.html、独立报告页和资源；共享在线布局/结果/结构/下载；本地突变草稿不显示计算分；禁用模型提交。", "逐任务ZIP成员/文件hash与主站、独立报告资源检查纳入报告；在线浏览器证明与file://策略阻塞逐例区分。", "file://交互验收被浏览器安全策略阻止的任务仍标未验证；静态ZIP与资源扫描不替代运行时离线UX验收。"),
]


def read_json(path: Path, default=None):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default


def load_csv(path: Path):
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def _task_status(task_dir: Path, trace, job_store_row=None):
    artifacts_present = (task_dir / "run_manifest.json").is_file() and (task_dir / "ecd_epitope_design.zip").is_file()
    if job_store_row:
        status = str(job_store_row["status"])
        stage = str(job_store_row["stage"] or status)
        if status == "COMPLETE":
            return "COMPLETE_ARTIFACTS_PRESENT" if artifacts_present else "COMPLETE_ARTIFACTS_MISSING"
        if status == "FAILED":
            return "FAILED_REQUIRES_REVIEW"
        if status in {"QUEUED", "RUNNING"}:
            return f"{status}:{stage}"
        return status
    if not trace:
        return "ARTIFACTS_PRESENT_TERMINAL_STATUS_UNKNOWN" if artifacts_present else "PENDING_NO_TRACE"
    last = trace[-1]
    event = str(last.get("event", ""))
    stage = str(last.get("stage", "UNKNOWN"))
    if event in {"PIPELINE_FAILED", "JOB_FAILED", "FAILED"}:
        return "FAILED_REQUIRES_REVIEW"
    if event in {"PIPELINE_COMPLETED", "JOB_COMPLETED", "COMPLETE"}:
        return "COMPLETE_ARTIFACTS_PRESENT" if artifacts_present else "TERMINAL_EXPORT_PENDING"
    return f"RUNNING:{stage}"


def _job_store_rows(root: Path, job_ids):
    if not job_ids:
        return {}
    path = root / "workspace_data" / "jobs.sqlite3"
    if not path.is_file():
        return {}
    try:
        with sqlite3.connect(path.as_uri() + "?mode=ro", uri=True) as conn:
            conn.row_factory = sqlite3.Row
            rows = conn.execute(
                f"SELECT id,status,stage FROM jobs WHERE id IN ({','.join('?' for _ in job_ids)})",
                tuple(job_ids),
            ).fetchall()
    except sqlite3.Error:
        return {}
    return {row["id"]: dict(row) for row in rows}


def _trace(task_dir: Path):
    rows = []
    path = task_dir / "execution_trace.jsonl"
    if path.is_file():
        for line in path.read_text(encoding="utf-8").splitlines():
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return rows


def _runtime_cache_rows(value, prefix="analysis", output=None):
    output = output if output is not None else []
    if isinstance(value, dict):
        keys = set(value)
        if keys & {"runtime_seconds", "cache_hit", "cache_disabled_reason", "gpu_lock_wait_seconds"}:
            output.append({"path": prefix, **{key: value.get(key) for key in (
                "tool", "tool_name", "runtime_seconds", "device", "gpu_lock_wait_seconds",
                "cache_hit", "cache_disabled_reason", "cache_key", "cache_identity_status") if key in value}})
        for key, child in value.items():
            _runtime_cache_rows(child, f"{prefix}.{key}", output)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _runtime_cache_rows(child, f"{prefix}[{index}]", output)
    return output


def _foldx_execution_timing(task_dir: Path, v4_manifest, request):
    receipts = []
    for path in sorted(task_dir.rglob("engine_execution.json")):
        value = read_json(path)
        if value is not None:
            command = " ".join(str(part) for part in value.get("command", []))
            receipts.append({"path": str(path), "command": command,
                "status": value.get("status"), "runtime_seconds": value.get("runtime_seconds")})
    build = [row for row in receipts if "BuildModel" in row["command"]]
    repair = [row for row in receipts if "RepairPDB" in row["command"]]
    worker = (((v4_manifest.get("worker_runs") or {}).get("foldx") or {}).get("provenance") or {})
    repair_path = worker.get("repair_cache_path")
    repair_receipt = Path(repair_path).parent / "engine_execution.json" if repair_path else None
    repair_execution = read_json(repair_receipt) if repair_receipt else None
    repair_file_exists = bool(repair_path and Path(repair_path).is_file())
    queued_at = (request.get("metadata", {}).get("execution_queue", {}).get("created_at")
                 if isinstance(request, dict) else None)
    repair_state = "UNKNOWN"
    if repair_execution and repair_path:
        try:
            receipt_time = repair_receipt.stat().st_mtime
            queued_ts = datetime.fromisoformat(queued_at.replace("Z", "+00:00")).timestamp() if queued_at else None
            repair_state = ("REUSED_PREEXISTING_REPAIR_CACHE" if queued_ts is not None and receipt_time < queued_ts
                            else "REPAIR_PDB_EXECUTED_DURING_TASK_OR_AFTER_QUEUE")
        except (OSError, ValueError):
            repair_state = "REPAIR_RECEIPT_PRESENT_TIME_SCOPE_UNKNOWN"
    elif repair_file_exists:
        repair_state = "REPAIR_CACHE_PRESENT_EXECUTION_TIME_UNKNOWN"
    elif repair_path:
        repair_state = "REPAIR_CACHE_PATH_RECORDED_BUT_FILE_NOT_FOUND"
    return {
        "buildmodel_receipt_count_in_task": len(build),
        "buildmodel_status_counts": {status: sum(row["status"] == status for row in build)
                                     for status in sorted({str(row["status"]) for row in build})},
        "buildmodel_worker_seconds_sum_across_batches": sum(float(row["runtime_seconds"] or 0) for row in build),
        "buildmodel_max_batch_seconds": max((float(row["runtime_seconds"] or 0) for row in build), default=None),
        "buildmodel_receipts": build,
        "repair_cache_path": repair_path,
        "repair_cache_state": repair_state,
        "repair_receipt_seconds_not_added_to_buildmodel_sum": repair_execution.get("runtime_seconds") if repair_execution else None,
        "repair_receipt_status": repair_execution.get("status") if repair_execution else None,
        "repair_runtime_scope_note": "RepairPDB cache receipt is shown separately; prior cache runtime is not added to current task BuildModel time.",
        "foldx_worker_runtime_seconds": ((v4_manifest.get("worker_runs") or {}).get("foldx") or {}).get("runtime_seconds"),
        "foldx_worker_reused_model_counts": worker.get("reused_model_counts"),
        "foldx_worker_new_model_counts": worker.get("new_model_counts"),
        "sum_is_cpu_core_time": False,
    }


def _load_audit(task_dir: Path):
    required = ["panel_candidates.json", "v4_manifest_v41.json", "selected_panel_v4.json",
                "sector_recommendations_v4.json", "residue_mapping.csv"]
    if not all((task_dir / name).is_file() for name in required):
        return {"status": "PENDING_REQUIRED_OUTPUTS", "checks": None}
    try:
        from audit_p0_p3_batch import audit_task
        value = audit_task(task_dir)
        return {"status": "AUDIT_COMPLETE", "checks": value["checks"],
                "counts": value["counts"], "budgets": value["budgets"],
                "candidate_strata": value["candidate_strata"],
                "selection_diversity_and_chemistry": value["selection_diversity_and_chemistry"],
                "panel_repeat_sensitivity": value.get("panel_repeat_sensitivity"),
                "paired_integrity_errors": value["paired_integrity_errors"],
                "source_hash_audit": value["source_hash_audit"],
                "export_hash_audit": value["export_hash_audit"]}
    except Exception as exc:  # Keep the report usable while preserving the failure.
        return {"status": "AUDIT_ERROR_REQUIRES_REVIEW", "error": f"{type(exc).__name__}: {exc}"}


def _external_benchmark(root: Path):
    value = read_json(root / "workspace_data" / "p0_p3_upgrade_20261016" / "benchmark" / "benchmark_results.json")
    if not value:
        return None
    keys = ("status", "score_blind_candidate_count", "chain_length_aa", "source_order_counts",
            "actual_spurs_training_overlap", "PTMUL_reference_train_clusters", "foldx_development_overlap",
            "foldx", "spurs", "metrics_by_order", "metrics_by_length_group", "identity_or_pair_failures",
            "interpretation", "evidence_sha256")
    return {key: value.get(key) for key in keys if key in value}


def _offline_site_static_scan(task_dir: Path):
    from html.parser import HTMLParser

    offline_pages = ("index.html", "offline_report.html")
    offline_assets = (*offline_pages, "assets/app.bundle.js", "assets/workbench.css", "assets/app.css",
                      "assets/renderer.js", "assets/vendor/3Dmol-min.js",
                      "assets/vendor/pdbe-molstar-plugin-3.10.1.js",
                      "assets/vendor/pdbe-molstar-light-3.10.1.css",
                      "assets/vendor/PDBe-Molstar-LICENSE")
    offline_present = all((task_dir / name).is_file() for name in offline_assets)
    remote_refs = []

    class ResourceScanner(HTMLParser):
        remote_pattern = re.compile(
            r"^(?:(?:https?:)?//(?:(?:[a-z0-9-]+\.)+[a-z0-9-]+|localhost|127\.0\.0\.1)(?::\d+)?"
            r"|https?%3A%2F%2F(?:(?:[a-z0-9-]+\.)+[a-z0-9-]+|localhost|127\.0\.0\.1))",
            re.IGNORECASE,
        )
        css_pattern = re.compile(r"(?:url\(\s*['\"]?|@import\s+['\"]?)([^\s)'\"]+)", re.IGNORECASE)

        def __init__(self):
            super().__init__()
            self.references = []
            self.in_style = False

        @classmethod
        def is_remote(cls, value):
            return bool(value and cls.remote_pattern.match(value.strip()))

        def handle_starttag(self, tag, attrs):
            values = dict(attrs)
            if tag == "link" and self.is_remote(values.get("href")):
                self.references.append(values["href"])
            if tag in {"script", "img", "iframe", "source", "video", "audio", "embed"}:
                for attr in ("src", "poster"):
                    if self.is_remote(values.get(attr)):
                        self.references.append(values[attr])
            if tag == "object" and self.is_remote(values.get("data")):
                self.references.append(values["data"])
            if tag == "form" and self.is_remote(values.get("action")):
                self.references.append(values["action"])
            if tag == "meta" and values.get("http-equiv", "").lower() == "refresh":
                target = re.search(r"url\s*=\s*([^;]+)", values.get("content", ""), re.IGNORECASE)
                if target and self.is_remote(target.group(1).strip(" '\"")):
                    self.references.append(target.group(1))
            if self.is_remote(values.get("srcset")):
                self.references.extend(part.strip().split()[0] for part in values["srcset"].split(",")
                                       if part.strip() and self.is_remote(part.strip().split()[0]))
            for match in self.css_pattern.findall(values.get("style", "")):
                if self.is_remote(match):
                    self.references.append(match)
            if tag == "style":
                self.in_style = True

        def handle_endtag(self, tag):
            if tag == "style":
                self.in_style = False

        def handle_data(self, data):
            if self.in_style:
                for match in self.css_pattern.findall(data):
                    if self.is_remote(match):
                        self.references.append(match)

    for name in offline_pages:
        page = task_dir / name
        if page.is_file():
            scanner = ResourceScanner()
            scanner.feed(page.read_text(encoding="utf-8", errors="replace"))
            remote_refs.extend(scanner.references)
    return offline_present, sorted(set(remote_refs))


def _summarize_job(item, baseline, derived, root: Path, job_store_row=None):
    task_dir = Path(item.get("workdir") or root / "workspace_data" / item["job_id"])
    trace = _trace(task_dir)
    status = _task_status(task_dir, trace, job_store_row)
    analysis = read_json(task_dir / "analysis.json", {}) or {}
    manifest = read_json(task_dir / "v4_manifest_v41.json", {}) or {}
    regions = analysis.get("regions", []) or []
    if not regions:
        regions = read_json(task_dir / "surface_display_regions_v41.json", []) or []
    if isinstance(regions, dict):
        regions = list(regions.values())
    if not isinstance(regions, list):
        regions = []
    sizes = []
    for region in regions:
        if not isinstance(region, dict):
            continue
        positions = region.get("residue_positions") or region.get("positions") or region.get("residues") or []
        if positions:
            sizes.append(len(positions))
    pool = manifest.get("provenance", {}).get("global_candidate_pool", {})
    repeat = pool.get("foldx_repeat_budget") or {}
    allocation = pool.get("foldx_allocation") or {}
    stage = trace[-1].get("stage") if trace else None
    elapsed = trace[-1].get("elapsed_monotonic_seconds") if trace else None
    model_sequence = analysis.get("protein", {}).get("sequence")
    residues = analysis.get("residues", []) or []
    chain_length = len(model_sequence) if model_sequence else len({r.get("uniprot_pos") for r in residues if r.get("uniprot_pos") is not None}) or None
    runtime_rows = []
    runtime_cache_diagnostics = _runtime_cache_rows(analysis)[:100] if analysis else []
    for filename in ("model_runtime_audit.json", "runtime_audit.json", "run_audit.json"):
        value = read_json(task_dir / filename)
        if value is not None:
            runtime_rows.append({"file": filename, "data": value})
    offline_present, remote_refs = _offline_site_static_scan(task_dir)
    offline_network_free = offline_present and not remote_refs
    browser_receipt = read_json(root / "workspace_data" / "p0_p3_upgrade_20261016" /
                                "offline_browser_audit" / "acceptance.json", {}) or {}
    browser_result = next((row for row in browser_receipt.get("results", []) if row.get("job_id") == item["job_id"]), None)
    validation_dir = root / "workspace_data" / "p0_p3_upgrade_20261016"
    run_acceptance = read_json(validation_dir / f"{item['job_id']}_acceptance.json", {}) or {}
    bundle_acceptance = read_json(validation_dir / f"{item['job_id']}_bundle_audit.json", {}) or {}
    browser_policy_block = read_json(validation_dir / "offline_browser_audit" / "file_protocol_policy_block.json", {}) or {}
    run_validation = {
        "status": run_acceptance.get("status"),
        "errors": run_acceptance.get("errors", []),
        "failed_checks": [name for name, passed in (run_acceptance.get("checks") or {}).items() if not passed],
        "pair_count": run_acceptance.get("pair_count"),
        "unique_combinations": run_acceptance.get("unique_combinations"),
        "source_provenance": run_acceptance.get("source_provenance"),
    } if run_acceptance else None
    bundle_validation = {key: bundle_acceptance.get(key) for key in ("status", "members", "case_status", "scope")}
    if not bundle_acceptance:
        bundle_validation = None
    offline_browser_audit = None
    if browser_result:
        offline_browser_audit = {
            "status": browser_result.get("status"), "opened_from": browser_result.get("opened_from"),
            "network_request_count": len(browser_result.get("network_requests", [])),
            "browser_error_count": len(browser_result.get("browser_errors", [])),
            "region_count": browser_result.get("checks", {}).get("regionCount"),
            "design_count": browser_result.get("checks", {}).get("designCount"),
            "offline_draft_status": browser_result.get("offline_draft", {}).get("status"),
            "offline_draft_download_count": len(browser_result.get("offline_draft", {}).get("downloads", [])),
            "screenshot": browser_result.get("screenshot"),
        }
    elif (validation_dir / "offline_browser_audit" / f"{item['job_id']}.png").is_file():
        offline_browser_audit = {
            "status": "PRIOR_BROWSER_ARTIFACTS_ONLY",
            "screenshot": f"offline_browser_audit/{item['job_id']}.png",
            "draft_artifacts_present": all((validation_dir / "offline_browser_audit" / name).is_file()
                for name in (f"{item['job_id']}_local_draft.json", f"{item['job_id']}_local_fasta.fasta",
                             f"{item['job_id']}_local_cds.fasta")),
            "limitation": "No structured browser acceptance receipt is retained for this earlier run.",
        }
    elif browser_policy_block:
        offline_browser_audit = {
            "status": "POLICY_BLOCKED_FILE_URL",
            "reason": browser_policy_block.get("reason"),
            "scope": browser_policy_block.get("scope"),
        }
    return {
        "protein": item["file"].split("_")[-1].replace(".pdb", ""),
        "input_file": item["file"], "source_job_id": item["source_job_id"],
        "new_job_id": item["job_id"], "url": item.get("url"),
        "source_sha256": item.get("source_sha256"), "task_directory": str(task_dir),
        "status": status, "job_store_status": job_store_row,
        "latest_stage": stage, "pipeline_elapsed_monotonic_seconds": elapsed,
        "chain_length_aa": chain_length,
        "historical_source_result": baseline,
        "historical_saved_pool_repartition": derived,
        "new_region_count": len(sizes) if sizes else None,
        "new_region_sizes_surface_nodes": sizes if sizes else None,
        "new_parameters": manifest.get("parameters") or None,
        "new_sequential_foldx_budget": repeat or None,
        "new_foldx_allocation": allocation or None,
        "runtime_receipts": runtime_rows,
        "runtime_cache_diagnostics": runtime_cache_diagnostics,
        "foldx_execution_timing": _foldx_execution_timing(task_dir, manifest, read_json(task_dir / "request.json", {})),
        "offline_site_files_present": offline_present if status == "COMPLETE_ARTIFACTS_PRESENT" else None,
        "offline_site_static_network_scan": ("PASS" if offline_network_free else "FAIL" if offline_present else "PENDING")
            if status == "COMPLETE_ARTIFACTS_PRESENT" else None,
        "offline_site_network_reference_examples": sorted(set(remote_refs))[:20],
        "offline_site_browser_audit": offline_browser_audit,
        "run_acceptance": run_validation,
        "bundle_acceptance": bundle_validation,
        "source_archive_present": (task_dir / "executed_source_snapshot.zip").is_file(),
        "export_manifest_present": (task_dir / "run_manifest.json").is_file(),
        "audit": _load_audit(task_dir),
    }


def build_report(manifest_path: Path, baseline_dir: Path, output_dir: Path, root: Path = ROOT):
    fresh_manifest = read_json(manifest_path)
    if not fresh_manifest:
        raise FileNotFoundError(f"Fresh-run manifest unavailable: {manifest_path}")
    baseline_rows = {row["source_job_id"]: row for row in load_csv(baseline_dir / "case_summary.csv")}
    historical_manifest = read_json(root / "workspace_data" / "binder_input_preflight_20261016" / "article_batch_manifest.json", {}) or {}
    derived_rows = {row["source_job_id"]: row for row in historical_manifest.get("small_protein_deliveries", [])}
    job_rows = _job_store_rows(root, [item["job_id"] for item in fresh_manifest.get("jobs", [])])
    jobs = [_summarize_job(item, baseline_rows.get(item["source_job_id"]),
                           derived_rows.get(item["source_job_id"]), root, job_rows.get(item["job_id"]))
            for item in fresh_manifest.get("jobs", [])]
    ready = [job for job in jobs if job["audit"]["status"] == "AUDIT_COMPLETE"]
    independent_audits_complete = len(ready) == len(jobs)
    terminal = len([job for job in jobs if job["status"] == "COMPLETE_ARTIFACTS_PRESENT"])
    core_validation_complete = terminal == len(jobs) and all(
        job.get("run_acceptance", {}).get("status") == "PASS" and
        job.get("bundle_acceptance", {}).get("status") == "PASS"
        for job in jobs)
    browser_validation_complete = all(
        job.get("offline_site_browser_audit", {}).get("status") == "PASS_FILE_URL_NO_NETWORK"
        for job in jobs)
    browser_policy_limited = any(
        job.get("offline_site_browser_audit", {}).get("status") == "POLICY_BLOCKED_FILE_URL"
        for job in jobs)
    benchmark = _external_benchmark(root)
    benchmark_status = (benchmark or {}).get("status")
    benchmark_terminal = benchmark_status in {"COMPLETE", "COMPLETE_WITH_MODEL_FAILURES"}
    benchmark_limited = benchmark_status == "COMPLETE_WITH_MODEL_FAILURES" or bool(
        (benchmark or {}).get("identity_or_pair_failures"))
    browser_policy_block = read_json(root / "workspace_data" / "p0_p3_upgrade_20261016" /
                                     "offline_browser_audit" / "file_protocol_policy_block.json", {}) or {}
    if any(job["status"] in {"FAILED_REQUIRES_REVIEW", "COMPLETE_ARTIFACTS_MISSING"} for job in jobs):
        report_status = "SIX_TASK_FAILURE_REQUIRES_REVIEW"
    elif terminal < len(jobs):
        report_status = "PENDING_SIX_FRESH_TASKS"
    elif not core_validation_complete or not independent_audits_complete:
        report_status = "SIX_TASK_ARTIFACTS_PRESENT_AUDIT_REVIEW_REQUIRED"
    elif not benchmark_terminal:
        report_status = "PENDING_EXTERNAL_BENCHMARK"
    elif browser_policy_limited and benchmark_limited:
        report_status = "COMPLETE_ENGINEERING_AUDIT_WITH_FILE_BROWSER_AND_BENCHMARK_LIMITATIONS"
    elif browser_policy_limited:
        report_status = "COMPLETE_ENGINEERING_AUDIT_WITH_FILE_BROWSER_LIMITATION"
    elif benchmark_limited:
        report_status = "COMPLETE_ENGINEERING_AUDIT_WITH_BENCHMARK_LIMITATION"
    elif browser_validation_complete:
        report_status = "COMPLETE_ENGINEERING_AUDIT_WITH_LIMITATIONS"
    else:
        report_status = "OFFLINE_BROWSER_AUDIT_REVIEW_REQUIRED"
    result = {
        "schema_version": "p0_p3_upgrade_evidence_report_v1",
        "generated_at_local": datetime.now().astimezone().isoformat(),
        "report_status": report_status,
        "fresh_manifest": str(manifest_path.resolve()), "fresh_manifest_protocol": fresh_manifest.get("protocol"),
        "fresh_manifest_status": fresh_manifest.get("status"), "fresh_manifest_sha256": fresh_manifest.get("source_manifest_sha256"),
        "fresh_task_count": len(jobs), "fresh_terminal_artifact_count": terminal,
        "independent_audits_complete": independent_audits_complete,
        "external_benchmark_terminal": benchmark_terminal,
        "external_benchmark_limitations": benchmark_limited,
        "old_baseline_report_dir": str(baseline_dir.resolve()), "tasks": jobs,
        "external_benchmark": benchmark,
        "offline_browser_policy_block": browser_policy_block or None,
        "p0_p3_checklist": [dict(zip(("phase", "area", "implementation_scope", "evidence", "remaining"), row)) for row in CHECKLIST],
        "scientific_boundary": [
            "文章来源或模型名称不能证明结构是实验解析；坐标结构/编号与任务来源按文件证据报告。",
            "FoldX/SPURS、几何和表面面积均为计算代理；工程验收不代表折叠、结合或实验效果。",
            "模型分歧、重复离散度、界面富集和化学描述均为运行内诊断，不是实验置信区间或因果解释。",
            "对仍运行任务，候选、预算、模型/导出核验结果必须保留PENDING，不用旧任务或派生重分类数据补齐。",
        ],
        "evidence_notes": {
            "original_baseline": "从binder_run_budget_review_20261016/case_summary.csv读取既有完整六例，不覆盖。",
            "old_derived": "从article_batch_manifest.json.small_protein_deliveries读取保存池重分类身份；这是派生展示任务，不是重新调用模型。",
            "new_tasks": "六个新job ID来自rerun_manifest.json；每个任务独立核对参数、源快照、FoldX配对、导出和离线站点。",
            "runtime": "阶段elapsed来自execution_trace；worker runtime来自工具回执。两者不能混成CPU核时或端到端墙钟。",
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "p0_p3_upgrade_report.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_dir / "p0_p3_upgrade_report_zh.md").write_text(render_markdown(result), encoding="utf-8")
    return result


def _fmt(value):
    if value is None:
        return "待产物生成"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def render_markdown(result):
    status = result["report_status"]
    lines = [
        "# P0–P3优化实施与六例新运行审查（滚动报告）", "",
        f"生成时间：{result['generated_at_local']}  ",
        f"当前状态：**{status}**；六例终态产物 {result['fresh_terminal_artifact_count']}/{result['fresh_task_count']}。", "",
        "## 本轮范围与解释", "",
        "本报告把已确认的P0–P3工程优化与六个新任务的实际运行证据分开呈现。历史完整任务和历史保存分区派生均保留为对照；新任务是另一次完整模型流程。五个AF3单链不在本轮。任务仍运行时，结果和验收保留待完成状态。", "",
        f"新任务唯一清单：`{result['fresh_manifest']}`（协议 `{result['fresh_manifest_protocol']}`）。", "",
        "本报告不把工程验收、稳定性代理、界面面积或本批模型重复离散度称为实验折叠/结合保证或实验校准。", "",
        "## 六例原始、历史派生与新任务映射", "",
        "| Binder | 原完整任务 | 历史保存池派生 | 新完整任务 | 长度 aa | 新分区数/大小 | 状态 | 最新trace阶段偏移 s |", "|---|---|---|---|---:|---|---|---:|",
    ]
    for task in result["tasks"]:
        old = task.get("historical_source_result") or {}
        derived = task.get("historical_saved_pool_repartition") or {}
        old_len = old.get("chain_length_aa")
        new_len = task.get("chain_length_aa")
        sizes = task.get("new_region_sizes_surface_nodes")
        lines.append(f"| {task['protein']} | `{task['source_job_id']}` | `{derived.get('delivery_job_id','—')}` | `{task['new_job_id']}` | {_fmt(new_len or old_len)} | {len(sizes) if sizes is not None else '待定'} / {_fmt(sizes)} | {task['status']} | {_fmt(task['pipeline_elapsed_monotonic_seconds'])} |")
    lines += ["", "长度取新任务分析结构中的规范化残基数；尚未生成分析时仅展示历史坐标链长，并明确不等同于天然完整构建长度。分区大小按表面图显示区成员计，不是连续序列片段。运行中的最新trace偏移不是任务墙钟或worker累计耗时；以完成后的独立时间收据分开统计。", "",
        "## 新协议与预算口径", "",
        "- <120 aa：分区目标7个表面残基，偏好6–8；总分区4–12为软范围。120–150 aa仍为目标8/偏好6–10；>150 aa仍为目标12/偏好10–15。连通/完整覆盖以及每区至少4残基、50 Å²、测地直径30 Å、Cα直径22 Å为硬门。最终值以每任务保存的实际parameters与分区诊断为准。",
        "- 新常规任务只保留双点/三点评分，四点不进入新模型请求。三点按本区增益和完整组合证据评估，不因阶数增加而放宽稳定性线。",
        "- SPURS完整组合≤1、FoldX≤2为正式稳定性代理门；SPURS 1–2或FoldX 2–3保留待审、不进正式推荐。取消模型差值和必须同号作为淘汰门；差值只记录诊断。实际正式状态仍由每个任务的完整资格审计决定。",
        "- 每蛋白最多80个唯一FoldX组合，请求60局部/20全局、局部强扰动探索至多16；本地/全局供给不足允许回流，不能把60/20误读成实际选中硬上限。每组先3次、最多30组补至5次，最多300配对槽。主面板与各区正式主代表必须五次；其他候选按序贯策略补算，缺失仍待审。",
        "- 主面板构建预算独立记录，不折算为分区数或FoldX数。下方逐任务列配置请求和实际选择。", "",
        "## 每任务实际过程、证据与运行状态", ""]
    for task in result["tasks"]:
        old = task.get("historical_source_result") or {}
        derived = task.get("historical_saved_pool_repartition") or {}
        audit = task.get("audit") or {}
        lines += [f"### {task['protein']}：{task['new_job_id']}", "",
            f"原完整任务 `{task['source_job_id']}`；历史派生任务 `{derived.get('delivery_job_id','未登记')}`；新任务 `{task['new_job_id']}`。输入 `{task['input_file']}`，来源SHA-256 `{task.get('source_sha256')}`。当前 **{task['status']}**，最近阶段 `{task.get('latest_stage') or '暂无trace'}`。", "",
            f"新链长：{_fmt(task.get('chain_length_aa'))} aa；新分区成员数：{_fmt(task.get('new_region_sizes_surface_nodes'))}；参数：`{json.dumps(task.get('new_parameters'), ensure_ascii=False, sort_keys=True) if task.get('new_parameters') else '尚未生成v4 manifest'}`。", "",
            f"FoldX预算：`{json.dumps(task.get('new_sequential_foldx_budget'), ensure_ascii=False, sort_keys=True) if task.get('new_sequential_foldx_budget') else '尚未生成'}`；局部/全局回流：`{json.dumps(task.get('new_foldx_allocation'), ensure_ascii=False, sort_keys=True) if task.get('new_foldx_allocation') else '尚未生成'}`。", "",
            f"历史原完整任务：链长 {old.get('chain_length_aa','未知')} aa，{old.get('candidate_count','未知')} 个保存候选、{old.get('foldx_candidates','未知')} 个FoldX组合、{old.get('panel_size','未知')} 个面板构建，状态 {old.get('job_status','未知')}。历史保存池派生分区：{derived.get('region_sizes','未知')}；该操作只重分类旧候选池，不是本轮新计算。", "",
            f"历史worker阶段耗时合计（非端到端墙钟）：{old.get('worker_elapsed_seconds_sum_not_pipeline_wall','未记录')} 秒；历史FoldX在已记录worker耗时中占比：{old.get('foldx_fraction_of_worker_elapsed_sum','未记录')}。新任务最新trace阶段偏移：{_fmt(task.get('pipeline_elapsed_monotonic_seconds'))} 秒（运行中并非当前总耗时）；新工具回执：`{json.dumps(task.get('runtime_cache_diagnostics'), ensure_ascii=False, sort_keys=True) if task.get('runtime_cache_diagnostics') else '分析导出前待生成'}`。", "",
            f"BuildModel执行回执（并行批次秒数求和，非墙钟/CPU核时）：`{json.dumps(task.get('foldx_execution_timing'), ensure_ascii=False, sort_keys=True)}`。RepairPDB缓存准备单列；如果复用历史缓存则不把历史RepairPDB耗时算入本任务。", "",
            f"本任务独立预算/导出审计：**{audit.get('status','PENDING')}**；身份/配对/源码验收：`{json.dumps(task.get('run_acceptance'), ensure_ascii=False, sort_keys=True) if task.get('run_acceptance') else '待完成'}`；ZIP验收：`{json.dumps(task.get('bundle_acceptance'), ensure_ascii=False, sort_keys=True) if task.get('bundle_acceptance') else '待完成'}`。离线站点文件已检查：{task.get('offline_site_files_present') if task.get('offline_site_files_present') is not None else '待导出'}；解压ZIP浏览器验收：`{json.dumps(task.get('offline_site_browser_audit'), ensure_ascii=False, sort_keys=True) if task.get('offline_site_browser_audit') else '待完成'}`；运行源码归档存在：{task.get('source_archive_present')}；导出manifest存在：{task.get('export_manifest_present')}。", ""]
        if audit.get("status") == "AUDIT_COMPLETE":
            lines.append(f"独立核对计数：`{json.dumps(audit.get('counts'), ensure_ascii=False, sort_keys=True)}`。")
            lines.append(f"预算硬检查：`{json.dumps(audit.get('checks'), ensure_ascii=False, sort_keys=True)}`。")
            lines.append(f"分层gap/重复统计：`{json.dumps(audit.get('candidate_strata'), ensure_ascii=False, sort_keys=True)}`。")
            lines.append(f"池→FoldX→formal→main-panel化学/位点/界面对照：`{json.dumps(audit.get('selection_diversity_and_chemistry'), ensure_ascii=False, sort_keys=True)}`。")
            lines.append(f"固定五次池的3vs5面板复选（与候选可用性分离）及全3次调度池对照（含可用性差异）：`{json.dumps(audit.get('panel_repeat_sensitivity'), ensure_ascii=False, sort_keys=True)}`。此项是本地读取已保存FoldX槽后重跑选择器，不触发模型。")
        elif audit.get("error"):
            lines.append(f"审计器错误需复核：`{audit['error']}`。")
        else:
            lines.append("候选池、SPURS、FoldX配对、正式候选、主面板、在线/离线导出的一致性尚未齐备，全部维持PENDING；不从历史任务填补。")
        lines.append("")
    lines += ["## P0–P3实施清单：已实现、证据与剩余优化", "",
        "| 类别 | 本轮实施内容 | 当前证据 | 剩余工作/限制 |", "|---|---|---|---|"]
    for phase, area, scope, evidence, remaining in CHECKLIST:
        lines.append(f"| {phase} · {area} | {scope} | {evidence} | {remaining} |")
    lines += ["", "## 偏倚和预算应如何解释", "",
        "候选漏斗是在多重结构、序列、表面、SPURS、几何和预算条件下筛出的条件分布。面板中化学身份、位点集或界面残基的富集只能说明选择结果与候选池构成不同，不能证明优化器造成了真实生物学提升。新批完成后，报告分别计算化学身份和位点集多样性，并比较候选池、实际FoldX调度池、正式合格池和主面板；供给和门槛造成的差别会单独呈现。", "",
        "FoldX重复槽只描述当前运行协议下的模型采样。重复内容相同、模型身份不完整、未调度、失败和有限重复都会单独保留，不把重复行当成独立实验样本。SPURS与FoldX差值不再是淘汰条件；仍需报告其符号、绝对值和重复稳定性，避免只看平均方向。", "",
        "本轮另外的数据基准工作与六例模型预算分开记账。FireProtDB/PTMUL原始记录下载不是验证完成；须先完成标签单位/符号/条件、作者编号、结构映射、去重、同源与训练重叠审查，再按最多60例协议运行。样本不足时只报告未充分校准。", "",
        "## 本阶段可以确认与不能确认", "",
        "可以确认：六个新job身份、输入哈希和任务路径来自唯一rerun_manifest；历史原完整任务和历史保存池派生身份可追溯；对已写出的新任务文件可读取实际运行参数、当前阶段和源快照。", "",
        "仍不能确认：只要新任务未写完全部证据，就不能确认80组配额实际使用、每组3/5次、正式门状态、偏倚统计、离线包完整或在线/离线一致。文件在目录出现也不等于验证通过。", "",
        "## 稳定性门槛的文献背景", "",
        "[SPURS（Nature Communications, 2025）](https://www.nature.com/articles/s41467-025-67609-4)报告了该模型在其指定的热稳定性与其他数据集上的结果；这些数据集、结构来源和评估切分并不等同于本批binder小蛋白。该文没有把本项目的SPURS≤1 kcal/mol规定为通用安全线。",
        "[FoldX力场更新研究（Bioinformatics, 2025）](https://academic.oup.com/bioinformatics/article/41/2/btaf064/8003679)重新整理了蛋白稳定性突变数据并提出力场更新，说明软件版本和评估数据会影响解释；本项目FoldX≤2 kcal/mol是操作门槛，不是预测误差带或置信区间。",
        "[CAGI5 frataxin盲测](https://pmc.ncbi.nlm.nih.gov/articles/PMC6744338/)比较多个方法对单一蛋白有限突变集合的稳定性预测，支持用预先冻结的外部数据评估模型，但不足以推导跨蛋白通用门槛。故本轮按用户确认的SPURS≤1、FoldX≤2执行，SPURS 1–2或FoldX 2–3列待审，外部小基准只作描述性检查，不据此自动调门。", "",
    ]
    benchmark = result.get("external_benchmark")
    lines += ["## 外部实测基准进度与分层结果", ""]
    if not benchmark:
        lines.append("尚无外部基准预测汇总；候选数据预审和运行账本单独保留。")
    else:
        foldx = benchmark.get("foldx") or {}
        spurs = benchmark.get("spurs") or {}
        lines += [f"状态 **{benchmark.get('status')}**；score-blind候选 {benchmark.get('score_blind_candidate_count')} 条。SPURS 已尝试 {spurs.get('attempted_candidate_count')} 条/成功评分 {spurs.get('scored_candidate_count')} 条；FoldX 已启动 {foldx.get('attempted_candidates')} 条/有结果 {foldx.get('successful_candidates')} 条，预占 {foldx.get('pair_slots_reserved')} / 批准 {foldx.get('approved_pair_slots')} 个配对槽；配对证明 {foldx.get('pair_proofs_validated')} 个，失败 {foldx.get('pair_validation_failures')} 个。", ""]
        recovery = benchmark.get("prior_attempt_recovery") or {}
        if recovery:
            lines += [f"续跑前有一次FoldX预算预检在启动新模型前中止，原收据保存在 `{recovery.get('snapshot_file')}`；根因：{recovery.get('root_cause')}。该次新增FoldX配对槽为0。", ""]
        for failure in spurs.get("failed_events", []):
            errors = failure.get("errors") or []
            reason = errors[0].get("message") if errors and isinstance(errors[0], dict) else "worker未返回具体错误"
            lines.append(f"SPURS未得分（保留失败、不重试）：`{json.dumps(failure.get('candidate_ids', []), ensure_ascii=False)}`；{reason}。")
        if spurs.get("failed_events"):
            lines.append("")
        lines += ["| 突变阶数 | SPURS n / MAE / RMSE / signed bias / Spearman | FoldX n / MAE / RMSE / signed bias / Spearman |", "|---:|---|---|"]
        by_order = benchmark.get("metrics_by_order")
        for order in ("1", "2", "3"):
            metrics = by_order.get(order) if isinstance(by_order, dict) else None
            metrics = metrics if isinstance(metrics, dict) else {}
            columns = []
            for model in ("SPURS", "FoldX"):
                row = metrics.get(model) or {}
                columns.append(" / ".join(_fmt(row.get(key)) for key in (
                    "n", "mae_kcal_mol", "rmse_kcal_mol", "mean_signed_error_kcal_mol", "spearman_rho")))
            lines.append(f"| {order} | {columns[0]} | {columns[1]} |")
        lines += ["", f"FoldX已启动但缺终态结果的候选：`{json.dumps(foldx.get('missing_result_candidates', []), ensure_ascii=False)}`；这些槽保持消耗且不重试。SPURS训练集关系：{benchmark.get('actual_spurs_training_overlap')}；PTMUL参考簇：{benchmark.get('PTMUL_reference_train_clusters')}。", "",
            "上述误差和排序只描述当前可用的小样本；训练集独立性未充分证明，PTMUL实验条件有缺失，单点与双/三点覆盖也不均衡。因此本轮不据此自动调整门槛，不称作独立通用校准。", ""]
    if result.get("offline_browser_policy_block"):
        lines += ["离线文件页的实时交互验收受浏览器安全策略阻止：策略仅开放HTTP/HTTPS，拒绝`file://`，并禁止改走其他浏览器入口或CDP。结构、ZIP、静态资源/网络引用与源代码回归仍可核验；受阻任务的七步切换、Mol*渲染和本地草稿下载保持未验证。", ""]
    lines += ["", "科学边界：文章来源不等于实验解析结构。SPURS、FoldX、结构几何、界面面积/接触与表面覆盖均为计算证据或代理；不能据此宣称折叠实测、亲和力、结合概率、功能或实验效果。", "",
        "## 证据文件", "",
        f"- 新任务清单：`{result['fresh_manifest']}`",
        f"- 历史六例基线：`{result['old_baseline_report_dir']}`",
        f"- 本报告结构化JSON：`{Path(result['fresh_manifest']).parent / 'p0_p3_upgrade_report.json'}`",
        "- 独立批次审计器：`scripts/audit_p0_p3_batch.py`（完整任务齐备后自动纳入本报告）。",
        f"- 外部基准分层结果：`{Path(result['fresh_manifest']).parent / 'benchmark' / 'benchmark_results.json'}`。", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=ROOT / "workspace_data/p0_p3_upgrade_20261016/rerun_manifest.json")
    parser.add_argument("--baseline-dir", type=Path, default=ROOT / "workspace_data/binder_run_budget_review_20261016")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "workspace_data/p0_p3_upgrade_20261016")
    args = parser.parse_args()
    result = build_report(args.manifest, args.baseline_dir, args.output_dir)
    print(json.dumps({"status": result["report_status"], "terminal": result["fresh_terminal_artifact_count"],
                      "tasks": result["fresh_task_count"], "report": str(args.output_dir / "p0_p3_upgrade_report_zh.md")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
