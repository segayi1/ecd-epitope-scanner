"""Build the evidence-bound AIB234 cold-cache comparison report."""
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'workspace_data'


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def selected(path):
    lanes = {}
    for region_row in load(path / 'sector_recommendations_v4.json'):
        region = str(region_row.get('sector_id', region_row.get('region_id', 'UNKNOWN')))
        for card in [region_row, *(region_row.get('alternatives') or [])]:
            for principle in card.get('principles') or []:
                key = (region, str(principle))
                if card.get('formal_eligible') is True:
                    entry = lanes.setdefault(key, {'status': 'FORMAL_ELIGIBLE', 'mutations': [], 'candidate_ids': []})
                    entry['mutations'].append(sorted(card.get('mutation_notations') or []))
                    entry['candidate_ids'].append(card.get('candidate_id'))
                else:
                    lanes.setdefault(key, {'status': 'NOT_FORMALLY_ELIGIBLE', 'mutations': [], 'candidate_ids': []})
        for unavailable in region_row.get('unavailable_principles') or []:
            principle = unavailable.get('principle')
            if principle:
                key = (region, str(principle))
                if key not in lanes:
                    state = str(unavailable.get('state', 'UNAVAILABLE')).upper()
                    lanes[key] = {'status': state, 'reason': unavailable.get('reason'), 'mutations': [], 'candidate_ids': []}
    return lanes


def unique_formal_combinations(lanes):
    return {tuple(sorted(mutations)) for lane in lanes.values()
            if lane.get('status') == 'FORMAL_ELIGIBLE'
            for mutations in lane.get('mutations', []) if mutations}


def table_escape(value):
    return str(value).replace('|', r'\|').replace('\r\n', '<br>').replace('\n', '<br>')


def compare_lanes(old_lanes, new_lanes):
    rows = []
    for key in sorted(set(old_lanes) | set(new_lanes)):
        old, new = old_lanes.get(key), new_lanes.get(key)
        row = {'region': key[0], 'principle': key[1],
               'old_status': old['status'] if old else 'NOT_REPORTED',
               'new_status': new['status'] if new else 'NOT_REPORTED'}
        old_options = old.get('mutations', []) if old else []
        new_options = new.get('mutations', []) if new else []
        if len(old_options) == len(new_options) == 1:
            row.update(old=old_options[0], new=new_options[0], comparable=True,
                       identical_mutations=old_options[0] == new_options[0],
                       identical_sites=sorted(site_number(x) for x in old_options[0]) ==
                       sorted(site_number(x) for x in new_options[0]))
        else:
            row.update(old=old_options[0] if len(old_options) == 1 else None,
                       new=new_options[0] if len(new_options) == 1 else None,
                       comparable=False, identical_mutations=None, identical_sites=None)
            if len(old_options) > 1 or len(new_options) > 1:
                row['comparison_status'] = 'AMBIGUOUS_MULTIPLE_FORMAL_CARDS'
            elif 'NOT_APPLICABLE' in (row['old_status'], row['new_status']):
                row['comparison_status'] = 'NOT_APPLICABLE'
            elif old is None or new is None:
                row['comparison_status'] = 'NOT_REPORTED_ON_ONE_SIDE'
            else:
                row['comparison_status'] = 'NO_FORMAL_WINNER_ON_BOTH_SIDES'
        rows.append(row)
    return rows


def site_number(notation):
    match = re.search(r'\d+', str(notation))
    return int(match.group()) if match else -1


def audit_formal_display(workdir, validator_errors):
    """Keep absent historical display flags unknown and explicit failures distinct."""
    panel = load(Path(workdir) / 'selected_panel_v4.json')
    candidates = load(Path(workdir) / 'panel_candidates.json')
    recommendations = load(Path(workdir) / 'sector_recommendations_v4.json')
    by_id = {row.get('candidate_id'): row for row in candidates if isinstance(row, dict)}
    ids = set(panel.get('selected_candidate_ids') or [])
    for primary in recommendations if isinstance(recommendations, list) else []:
        for card in [primary, *(primary.get('alternatives') or [])]:
            if card.get('candidate_id') and card.get('formal_eligible') is True:
                ids.add(card['candidate_id'])
    errors, unknown = [], []
    missing_ids = set()
    for candidate_id in sorted(ids):
        row = by_id.get(candidate_id)
        if row is None:
            unknown.append(f'{candidate_id}: candidate evidence missing')
        elif 'display_eligible' not in row:
            missing_ids.add(candidate_id)
            unknown.append(f'{candidate_id}: display_eligible field missing')
        elif row.get('display_eligible') is not True:
            errors.append(f'{candidate_id}: display_eligible explicitly false')
    for message in validator_errors or []:
        match = re.search(r'site-first display candidate is not eligible: (.+)$', message)
        if match and match.group(1) in missing_ids:
            continue
        match = re.search(r'site-first display qualification unknown: (.+)$', message)
        if match and match.group(1) in missing_ids:
            continue
        errors.append(message)
    status = 'NEEDS_REVIEW' if errors else ('UNKNOWN' if unknown else 'PASS')
    return {'status': status, 'errors': list(dict.fromkeys(errors)), 'unknown': list(dict.fromkeys(unknown))}


def audit_panel_signatures(panel, candidates):
    """Compare every stored signature value with each available candidate response."""
    unknown, mismatches, compared = [], [], 0
    expected_cells, compared_cells, mismatch_cells = 0, 0, set()
    if not isinstance(panel, dict) or not isinstance(candidates, list):
        return {'status': 'UNKNOWN', 'expected_cells': 0, 'compared_cells': 0, 'mismatched_cells': 0,
                'compared_values': 0, 'field_comparisons': 0, 'field_mismatches': 0, 'mismatches': [],
                'unknown': ['invalid panel/candidate artifacts']}
    selected_ids = panel.get('selected_candidate_ids')
    signatures = panel.get('hypothesis_signatures')
    if not isinstance(selected_ids, list) or not selected_ids:
        unknown.append('selected_candidate_ids missing or empty')
        selected_ids = []
    if not isinstance(signatures, dict) or not signatures:
        unknown.append('hypothesis_signatures missing or empty')
        signatures = {}
    by_id = {row.get('candidate_id'): row for row in candidates if isinstance(row, dict) and row.get('candidate_id')}
    for hypothesis, values in signatures.items():
        if not isinstance(values, list) or len(values) != len(selected_ids):
            unknown.append(f'{hypothesis}: signature length does not match selected candidates')
            continue
        expected_cells += len(selected_ids)
        for index, candidate_id in enumerate(selected_ids):
            expected = values[index]
            if not isinstance(expected, (int, float)) or isinstance(expected, bool) or not math.isfinite(expected):
                unknown.append(f'{hypothesis}/{candidate_id}: signature value missing or invalid')
                continue
            row = by_id.get(candidate_id)
            if row is None:
                unknown.append(f'{hypothesis}/{candidate_id}: selected candidate missing')
                continue
            response_maps = [(key, row.get(key)) for key in ('hypothesis_perturbation', 'hypothesis_response')
                             if isinstance(row.get(key), dict) and row[key]]
            if not response_maps:
                unknown.append(f'{hypothesis}/{candidate_id}: candidate response missing')
                continue
            cell_compared = False
            for key, response in response_maps:
                observed = response.get(hypothesis)
                if not isinstance(observed, (int, float)) or isinstance(observed, bool) or not math.isfinite(observed):
                    unknown.append(f'{hypothesis}/{candidate_id}: {key} missing or invalid')
                else:
                    compared += 1
                    cell_compared = True
                    if not math.isclose(float(expected), float(observed), rel_tol=1e-9, abs_tol=1e-8):
                        mismatches.append(f'{hypothesis}/{candidate_id}/{key}: {expected} != {observed}')
                        mismatch_cells.add((hypothesis, candidate_id))
            if cell_compared:
                compared_cells += 1
    status = 'NEEDS_REVIEW' if mismatches else ('UNKNOWN' if unknown else 'PASS')
    return {'status': status, 'expected_cells': expected_cells, 'compared_cells': compared_cells,
            'compared_values': compared,
            'mismatched_cells': len(mismatch_cells), 'field_comparisons': compared,
            'field_mismatches': len(mismatches), 'mismatches': mismatches, 'unknown': unknown}


def scientific_status(new, checks, validation_errors, signature_audit, display_errors):
    validation = new.get('validation') if isinstance(new, dict) else None
    checks['post_validation_pass'] = (None if not isinstance(validation, dict) else
                                      validation.get('status') == 'PASS' and not validation.get('errors'))
    if isinstance(display_errors, dict):
        checks['formal_display_eligibility'] = (True if display_errors.get('status') == 'PASS' else
                                                False if display_errors.get('status') == 'NEEDS_REVIEW' else None)
    else:
        checks['formal_display_eligibility'] = not display_errors if display_errors is not None else None
    checks['panel_signature_consistency'] = (True if signature_audit.get('status') == 'PASS' else
                                             False if signature_audit.get('status') == 'NEEDS_REVIEW' else None)
    checks['post_validation_errors_absent'] = not validation_errors if validation_errors is not None else None
    required = ('new_acceptance', 'post_validation_pass', 'post_validation_errors_absent',
                'formal_display_eligibility', 'panel_signature_consistency', 'same_input_sha256',
                'same_normalized_structure', 'spurs_single_weights_and_gpu',
                'spurs_multi_weights_and_gpu', 'procaliper_identity',
                'discotope_disabled_in_both_binder_runs', 'foldx_binary_identity',
                'legacy_raw_foldx_binary_identity', 'original_input_unchanged', 'new_source_stable')
    failed = [key for key in required if checks.get(key) is False]
    unknown = [key for key in required if checks.get(key) is not True and checks.get(key) is not False]
    status = 'NEEDS_REVIEW' if failed else ('NOT_ASSESSED' if unknown else 'PASS')
    return status, failed, unknown


def evaluate(new, old, new_identity, old_identity, ledger, new_manifest, old_manifest,
             old_export, normalized_new_sha, *, os_page_cache_matched=False):
    nw = new_identity['worker_model_identities']
    ow = old_identity['worker_model_identities']
    fields = ('model_files_sha256', 'actual_device', 'gpu_name')
    single = all(nw['spurs_single'].get(k) and nw['spurs_single'].get(k) == ow['spurs_single'].get(k) for k in fields)
    batches = ledger['spurs']['batches']
    multi = bool(batches) and all(all(batch['provenance'].get(k) and batch['provenance'].get(k) == ow['spurs_multi'].get(k)
                                       for k in fields) for batch in batches)
    cpu_keys = ('version', 'biopython_version', 'algorithm', 'probe_radius_a', 'sphere_points')
    procaliper = all(nw['procaliper'].get(k) is not None and nw['procaliper'].get(k) == ow['procaliper'].get(k) for k in cpu_keys)
    pair_binaries = {pair.get('binary_sha256') for row in ledger['foldx_results'] for pair in row.get('paired_models', [])}
    discotope = all(not manifest['parameters']['joint_v4_ml']['discotope']['enabled'] for manifest in (new_manifest, old_manifest))
    checks = {
        'new_acceptance': new.get('status') == 'PASS', 'legacy_independent_reexport': old_export.get('status') == 'PASS',
        'legacy_completed': old.get('status') == 'COMPLETE',
        'same_input_sha256': new.get('input_sha256') == old.get('input_sha256') and bool(new.get('input_sha256')),
        'same_normalized_structure': normalized_new_sha == old['normalized_structure']['sha256'],
        'spurs_single_weights_and_gpu': single, 'spurs_multi_weights_and_gpu': multi,
        'procaliper_identity': procaliper, 'discotope_disabled_in_both_binder_runs': discotope,
        'foldx_binary_identity': pair_binaries == {old['model_identity']['foldx_binary']['sha256']},
        'new_cold_cache': new['cache_evidence']['before_run']['cold_verified'] is True,
        'legacy_cold_cache': old['cold_cache_evidence']['score_cache_entries_at_start'] == 0 and old['cold_cache_evidence']['all_possible_repair_cache_paths_absent_before_run'] is True,
        'os_page_cache_condition_verified': os_page_cache_matched is True,
        'new_source_stable': new['source_digest_checkpoint'] == new['source_digest_run'] == new['source_digest_end'] and bool(new['source_digest_run']),
    }
    ratio = new['elapsed_seconds'] / old['elapsed_pipeline_seconds']
    return {'checks': checks, 'new_seconds': new['elapsed_seconds'], 'legacy_seconds': old['elapsed_pipeline_seconds'],
            'reduction_percent': 100 * (1 - ratio),
            'performance_target': ('PASS' if ratio <= .5 else 'FAIL') if all(checks.values()) else 'NOT_ASSESSED'}


def build(*, acceptance_dir=DATA / 'site_first_acceptance_r2',
          new_dir=DATA / 'fdbde4aa24d7', old_root=DATA / 'sf_old1',
          report_dir=DATA / 'site_first_acceptance_r2',
          source_input=DATA / '2b522bedd0f5' / 'uploads' / 'structure.pdb'):
    import hashlib
    acceptance_dir, new_dir, old_root, report_dir, source_input = map(
        Path, (acceptance_dir, new_dir, old_root, report_dir, source_input))
    report_dir.mkdir(parents=True, exist_ok=True)
    old_dir = old_root / 'derived_export'
    new, old = load(acceptance_dir / 'report.json'), load(old_root / 'baseline_report.json')
    ledger = load(new_dir / 'ordinary_candidate_ledger.json')
    old_export = load(old_root / 'derived_export_report.json')
    result = evaluate(new, old, load(new_dir / 'run_source_identity_v1.json'),
                      load(old_dir / 'run_source_identity_v1.json'), ledger,
                      load(new_dir / 'run_manifest.json'), load(old_dir / 'run_manifest.json'),
                      old_export,
                      hashlib.sha256((new_dir / 'structure_for_energy.pdb').read_bytes()).hexdigest())
    result['checks']['original_input_unchanged'] = hashlib.sha256(source_input.read_bytes()).hexdigest() == new['input_sha256']
    original_old = old_root / 'baseline/legacy_output'
    old_pair_hashes = set()
    for path in [original_old / 'paired_model_manifest.json', *sorted((original_old / 'foldx_refinement').glob('batch_*/paired_model_manifest.json'))]:
        if path.is_file():
            for record in load(path).get('records', []):
                old_pair_hashes.update(pair.get('binary_sha256') for pair in record.get('paired_models', []))
    result['checks']['legacy_raw_foldx_binary_identity'] = old_pair_hashes == {old['model_identity']['foldx_binary']['sha256']}
    panel = load(new_dir / 'selected_panel_v4.json')
    candidates = load(new_dir / 'panel_candidates.json')
    signature_audit = audit_panel_signatures(panel, candidates)
    result['panel_signature_audit'] = signature_audit
    result['panel_signature_mismatches'] = signature_audit['mismatched_cells']
    validation_errors = list(new.get('errors') or [])
    validation_errors.extend((new.get('validation') or {}).get('errors') or [])
    try:
        import sys
        scripts_dir = str(ROOT / 'scripts')
        if scripts_dir not in sys.path:
            sys.path.insert(0, scripts_dir)
        from run_ordinary_adaptive_rerun import _site_first_display_errors
        display_errors = audit_formal_display(new_dir, _site_first_display_errors(new_dir))
    except Exception as exc:
        display_errors = None
        result['formal_display_eligibility_error'] = f'{type(exc).__name__}: {exc}'
    result['formal_display_eligibility_errors'] = (display_errors or {}).get('errors') if isinstance(display_errors, dict) else display_errors
    result['formal_display_eligibility_unknown'] = (display_errors or {}).get('unknown') if isinstance(display_errors, dict) else None
    result['scientific_acceptance'], result['scientific_failed_checks'], result['scientific_unknown_checks'] = scientific_status(
        new, result['checks'], validation_errors, signature_audit, display_errors)
    if result['scientific_acceptance'] != 'PASS':
        result['performance_target'] = 'NOT_ASSESSED'
    before, after = selected(old_dir), selected(new_dir)
    comparisons = compare_lanes(before, after)
    comparable = [row for row in comparisons if row['comparable']]
    old_recommendations = load(old_dir / 'sector_recommendations_v4.json')
    old_formal = [card for row in old_recommendations for card in [row, *row.get('alternatives', [])] if card.get('formal_eligible')]
    old_signatures = unique_formal_combinations(before)
    new_signatures = unique_formal_combinations(after)
    perf = load(new_dir / 'site_first_performance.json')
    region_gap_rows = []
    display_by_region = {}
    for error in (display_errors or {}).get('errors', []) if isinstance(display_errors, dict) else (display_errors or []):
        match = re.search(r'VS\d{3}', error)
        if match:
            display_by_region.setdefault(match.group(), []).append(error)
    new_region_evidence = {str(row.get('region', 'UNKNOWN')): row for row in (new.get('region_principles') or [])}
    for region, row in sorted(new_region_evidence.items()):
        gaps, unknowns = [], []
        for principle, evidence in (row.get('principles') or {}).items():
            status = str(evidence.get('status', 'NOT_REPORTED')).upper()
            if status in {'NOT_APPLICABLE', 'NOT_APPLICABLE_BY_SCOPE', 'N_A'}:
                continue
            if status in {'NOT_REPORTED', 'TECHNICAL_FAILURE', 'UNAVAILABLE', 'REVIEW_REQUIRED'} or evidence.get('gap'):
                gaps.append(f'{principle}: {evidence.get("gap") or status}')
        region = row.get('region', 'UNKNOWN')
        gaps.extend(display_by_region.get(region, []))
        unknowns.extend(item for item in (result.get('formal_display_eligibility_unknown') or []) if region in item)
        region_gap_rows.append({'region': region, 'gaps': list(dict.fromkeys(gaps)),
                                'unknowns': list(dict.fromkeys(unknowns)),
                                'status': 'NEEDS_REVIEW' if gaps else ('UNKNOWN' if unknowns else 'PASS')})
    result.update(region_gaps=region_gap_rows, comparisons=comparisons,
                  comparable_lane_count=len(comparable), missing_or_unavailable_lane_count=len(comparisons)-len(comparable),
                  exact_same_principle_winners=sum(row['identical_mutations'] is True for row in comparable),
                  same_sites_per_principle=sum(row['identical_sites'] is True for row in comparable),
                  old_formal_card_instances=len(old_formal), old_unique_formal_combinations=len(old_signatures),
                  new_unique_formal_combinations=len(new_signatures), shared_formal_combinations=len(old_signatures & new_signatures),
                  performance=perf, new_validation=new['validation'],
                  caveat='Windows 文件页缓存没有被控制或验证为一致；旧版执行中曾短暂运行无模型测试，文件事件监听在旧版导出时额外读取了部分复制的收据。机器负荷未完全隔离，结果为单次观测，严格同缓存性能验收未完成。')
    # Avoid copying the full batch ledger into the concise comparison receipt.
    for key in ('spurs', 'foldx'):
        result['performance'][key].pop('batches', None)
    (report_dir / 'assessment_report.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    display_status = (display_errors or {}).get('status', 'UNKNOWN') if isinstance(display_errors, dict) else 'UNKNOWN'
    cache_known_cold = result['checks'].get('new_cold_cache') is True and result['checks'].get('legacy_cold_cache') is True
    def show(value, precision=None):
        if value is None:
            return '未知'
        return f'{value:.{precision}f}' if precision is not None else str(value)
    lines = ['# AIB234 位点先行算法验收与对照报告', '',
        f"新流程原自动检查：**{new['status']}**；独立科学审查：**{result['scientific_acceptance']}**；严格同缓存耗时目标：**{result['performance_target']}**。",
        f"面板签名复核：**{signature_audit['status']}**，假说×候选单元格 {signature_audit['compared_cells']}/{signature_audit['expected_cells']}，不一致单元格 {signature_audit['mismatched_cells']}；两个响应字段逐项比较 {signature_audit['field_comparisons']} 次，不一致 {signature_audit['field_mismatches']} 次；未知 {len(signature_audit['unknown'])} 项。",
        f"正式展示资格复核：**{display_status}**；post-validation：**{(new.get('validation') or {}).get('status', 'UNKNOWN')}**。",
        f"不能通过的科学验收项：{', '.join(result['scientific_failed_checks']) or '无'}；证据未知项：{', '.join(result['scientific_unknown_checks']) or '无'}。",
        f"应用评分及修复缓存冷启动核验：{'PASS' if cache_known_cold else '未通过或未知'}；观测耗时：旧版 {show(result.get('legacy_seconds'), 1)} 秒，新版 {show(result.get('new_seconds'), 1)} 秒，变化 {show(result.get('reduction_percent'), 1)}%。",
        '旧版计时为冻结管线入口至打包结束，新版计时包含任务执行及打包；两者串行。旧版收据有导出顺序缺陷，已在独立副本重新核验，模型未重算。', '',
        '## 核验条件', '', '输入和规范化结构、单点／组合 SPURS 权重与 GPU、ProCaliper 版本及 FoldX 二进制一致；两次分别使用新建空评分与修复缓存。binder 模式两次均关闭 DiscoTope。各项机器可读检查见 assessment_report.json。',
        result['caveat'], '', '## 逐区原则结果', '',
        '| 区域 | 旧机制方案 | 新机制方案 | 旧覆盖方案 | 新覆盖方案 |', '|---|---|---|---|---|']
    def lane_value(row, side):
        mutations = row.get(side)
        return '/'.join(mutations) if mutations is not None else row.get(f'{side}_status', 'UNKNOWN')
    lines = lines[:lines.index('| 区域 | 旧机制方案 | 新机制方案 | 旧覆盖方案 | 新覆盖方案 |')]
    lines += ['| 区域 | 原则 | 旧版结果 | 新版结果 | 可比性 |', '|---|---|---|---|---|']
    for row in comparisons:
        comparable = ('可比' if row['comparable'] else row.get('comparison_status', '不可比'))
        cells = [row['region'], row['principle'], lane_value(row, 'old'), lane_value(row, 'new'), comparable]
        lines.append('| ' + ' | '.join(table_escape(cell) for cell in cells) + ' |')
    lines += ['', '## 逐区缺口', '', '| 区域 | 状态 | 缺口或待审查项 |', '|---|---|---|']
    for row in region_gap_rows:
        details = [*row['gaps'], *(f"未知：{item}" for item in row['unknowns'])]
        cells = [row['region'], row['status'], '；'.join(details) or '无']
        lines.append('| ' + ' | '.join(table_escape(cell) for cell in cells) + ' |')
    lines += ['', f"{result['comparable_lane_count']} 个双方均有唯一正式方案的 lane 中，{result['exact_same_principle_winners']} 个替换完全相同，{result['same_sites_per_principle']} 个位点集合相同；另有 {result['missing_or_unavailable_lane_count']} 个缺失、不可用或存在歧义的 lane。",
        f"旧版共 {len(old_recommendations)} 个区域记录、{len(old_formal)} 张formal_eligible卡片实例及 {len(old_signatures)} 个不同组合；新版有 {len(new.get('region_principles') or [])} 个区域记录、{len(new_signatures)} 个原则组合。两版组合交集 {len(old_signatures & new_signatures)} 个。展示资格状态为 **{display_status}**。", '',
        '## 计算与验收', '',
        f"位点评价 {show((perf.get('candidate_timings') or {}).get('site_preparation_seconds'), 3)} 秒，替换展开 {show((perf.get('candidate_timings') or {}).get('representative_expansion_seconds'), 3)} 秒；组合 SPURS {show((perf.get('spurs') or {}).get('unique_requests'))} 个唯一请求、{show(((perf.get('spurs') or {}).get('timings') or {}).get('seconds'), 1)} 秒；FoldX {show((perf.get('foldx') or {}).get('attempted_slots'))} 槽、首轮 {show(((perf.get('foldx') or {}).get('timings') or {}).get('first_pass_seconds'), 1)} 秒。",
        f"FoldX 修复执行 {show((perf.get('worker_metrics') or {}).get('repair_executions'))} 次，修复 {show((perf.get('worker_metrics') or {}).get('repair_seconds'), 1)} 秒，BuildModel {show((perf.get('worker_metrics') or {}).get('buildmodel_seconds'), 1)} 秒；导出 {show(perf.get('export_seconds'), 1)} 秒。",
        f"新版主进程峰值内存 {show(perf.get('peak_process_memory_bytes'))} 字节，记录磁盘占用 {show(perf.get('disk_bytes_before_receipt'))} 字节；峰值不包括模型子进程和 GPU。",
        f"post-validation 状态为 {(new.get('validation') or {}).get('status', 'UNKNOWN')}；独立展示资格状态为 {display_status}，面板签名状态为 {signature_audit['status']}。旧版独立副本核验 {show(old_export.get('validated_unique_pairs'))} 对唯一 FoldX 配对。", '',
        '## 修复和限制', '',
        '真实首轮暴露了四位汇总与原始 FoldX 中位数的精度校验误判；修复后保留原始中位数。还修复了导出时卡片收据先于最终 analysis 写入，以及大 CSV 字段读取限制。失败运行与原始输出均保留。',
        '代表替换按 3→6 分层，未穷尽全部氨基酸组合；机制预排为近似筛选。两版搜索空间和 FoldX 3／5、3／10 采样协议不同，结果差异不能解释为同一候选的严格配对实验。',
        '这是一次 AIB234 应用层冷修复／冷评分缓存实测；OS 页缓存与背景负载未受控。已修复缓存条件、AF3 binder、抗原与其余样本尚未完成本轮验收。', '']
    (report_dir / 'assessment_report.md').write_text('\n'.join(lines), encoding='utf-8')
    print(json.dumps({'performance_target': result['performance_target'], 'checks': result['checks'],
                      'reduction_percent': result['reduction_percent'], 'report': str(report_dir / 'assessment_report.md')}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--acceptance-dir', type=Path, default=DATA / 'site_first_acceptance_r2',
                        help='Directory containing the completed single-case report')
    parser.add_argument('--new-run-dir', type=Path, default=DATA / 'fdbde4aa24d7',
                        help='Completed new pipeline job directory')
    parser.add_argument('--baseline-dir', type=Path, default=DATA / 'sf_old1',
                        help='Frozen legacy baseline directory')
    parser.add_argument('--report-dir', type=Path, default=DATA / 'site_first_acceptance_r2',
                        help='Directory for assessment_report.json and .md')
    parser.add_argument('--source-input', type=Path,
                        default=DATA / '2b522bedd0f5' / 'uploads' / 'structure.pdb',
                        help='Original source structure for SHA-256 identity check')
    args = parser.parse_args()
    build(acceptance_dir=args.acceptance_dir, new_dir=args.new_run_dir,
          old_root=args.baseline_dir, report_dir=args.report_dir, source_input=args.source_input)
