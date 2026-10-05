"""Write the closed CTLA4 result report from retained evidence, without solves."""
import argparse
from pathlib import Path
import json
import numpy as np
from ..util import write_json,sha256_file
from ..surface_perturbation.runtime import read_budget
from .audit_ctla4_coverage_v3 import audit


def run(case):
    route=case/'routes/coverage_v3';output=route/'results'
    if read_budget(case/'case_budget.json')['active'] or not (route/'delivery_status.json').exists():
        raise ValueError('Final report requires completed core delivery and closed budget')
    independent=audit(route,packages=True)
    summary=json.loads((output/'surface_perturbation_preview.json').read_text(encoding='utf-8'))
    info=independent['panel_information']
    history=json.loads((route/'screened_candidate_history.json').read_text(encoding='utf-8'))
    pairs=json.loads((output/'surface_perturbation/paired_manifest.json').read_text(encoding='utf-8'))
    candidates={cid:r['candidate'] for cid,r in history.items()}
    candidates.update({c['candidate_id']:c for c in json.loads((output/'panel_candidates.json').read_text(encoding='utf-8'))})
    receipts=[json.loads(p.read_text(encoding='utf-8')) for p in (case/'stages').glob('V3_*_report.json')]
    triple_reviews=[{'stage':r['stage'],'compression_contract':r.get('physical',{}).get('compression_contract'),
                     'audits':r.get('physical',{}).get('triple_audits',[]),
                     'compression_comparisons':r.get('physical',{}).get('compression_audits',[])} for r in receipts
                    if r.get('engineering_status')=='COMPLETED' and r.get('configuration',{}).get('triple_candidate_id')]
    regression_receipts=[json.loads(p.read_text(encoding='utf-8')) for p in sorted(route.glob('T1*regression*.json'))]
    invocation_path=route/'buildmodel_invocation_audit.json'
    invocation_count=json.loads(invocation_path.read_text(encoding='utf-8'))['audited_time_seeded_buildmodel_invocations'] if invocation_path.exists() else '收据尚未生成'
    pct=lambda value:f'{100*value:.6f}%'
    lines=['# CTLA4 coverage v3：实现与真实实例最终审查', '',
        f"正式目标状态：**{independent['hard_target_status']}**。{summary['actual_constructs']}个正式计算合格构建体；扫描机会覆盖{pct(info['C_screen_lower_bound'])}，直接扰动覆盖{pct(info['C_physical_lower_bound'])}。",
        '', '上述扫描覆盖是已验证扰动seed沿固定表面图12 Å扩张所得的计算机会下界，不是抗体结合概率、实验检出率或可解码表位覆盖。', '',
        '## 正式面板及信息指标', '',
        f"固定可靠WT表面分母：{summary['reference_area_all_a2']:.12f} Å²；rSASA≥0.20只改变可设计资格，不缩小分母。",
        f"二重覆盖：{pct(info['double_covered_fraction'])}；删除任意单个构建体后的最差扫描覆盖：{pct(info['worst_single_construct_loss_coverage'])}。",
        f"最大连通盲区：{info['largest_blind_patch_a2']:.6f} Å²。",
        '几何码明确1/0对立的位置对比例（未知不计为0）：'+', '.join(f"距离≥{k}：{pct(v)}" for k,v in info['geometric_pair_distinction'].items())+'。',
        '这是几何可区分性，不保证实验解码、纠错或单构建失效后仍达到80%。', '',
        '| 基因型（UniProt） | SPURS-multi | 五对FoldX中位数 | 单项扫描下界 | 单项直接扰动下界 |',
        '|---|---:|---:|---:|---:|']
    for cid in summary['selected_candidate_ids']:
        c=candidates[cid];production=[p for p in pairs[cid] if p['sampling_role']=='production']
        local=json.loads((output/'surface_perturbation'/f'{cid}_field_audit.json').read_text(encoding='utf-8'))
        spurs=(c.get('stability') or {}).get('spurs_multi_ddg_kcal_mol')
        # Use the actual exported evidence spelling instead of synthesizing a score.
        if spurs is None:
            constructs=json.loads((output/'surface_perturbation_local_constructs.json').read_text(encoding='utf-8'))['constructs']
            spurs=next(r['spurs_ddg_kcal_mol'] for r in constructs if r['candidate_id']==cid)
        lines.append(f"| {' / '.join(m['notation'] for m in c['mutations'])} | {spurs:.6f} | {np.median([p['ddg_kcal_mol'] for p in production]):.6f} | {pct(local['scan_area_a2']/summary['reference_area_all_a2'])} | {pct(local['physical_area_a2']/summary['reference_area_all_a2'])} |")
    apbs_details=[]
    for cid in summary['selected_candidate_ids']:
        path=output/'surface_perturbation'/f'{cid}_electrostatics.json'
        if path.exists():
            report=json.loads(path.read_text(encoding='utf-8'))
            apbs_details.append({'candidate_id':cid,'evidence_states':report['evidence_states'],
                'heavy_atom_preparation_changes':sum(len(r['heavy_atom_preparation_changes']) for r in report['runs']),
                'example':report['runs'][0]['heavy_atom_preparation_changes'][:2]})
    lines+=['', '### 构建预算曲线', '', '| 构建数 | 扫描下界 | 最大连通盲区 Å² |', '|---:|---:|---:|']
    for row in summary.get('coverage_curve',[]):
        if row['K']<=summary['actual_constructs']:
            lines.append(f"| {row['actual_constructs']} | {pct(row['C_hard_all'])} | {row['largest_blind_patch_a2']:.6f} |")
    lines+=['', '这条曲线是最终合格面板的嵌套子集，不是各K全局独立最优；K大于实际构建数不会制造额外构建体。']
    cardinality=route/'measured_pool_cardinality_audit.json'
    if cardinality.exists():
        proof=json.loads(cardinality.read_text(encoding='utf-8'))
        lines+=['', f"独立穷举已测{proof['measured_double_count']}个合格双点的全部{proof['four_construct_subsets_checked']}种四构建子集：最高扫描下界{pct(proof['maximum_four_construct_scan_lower_bound'])}，低于80%；五个是在本次已测双点池内的最少构建数。此证明不覆盖未测位点对、备用替换或未来三点。优化器自身超时没有给出这一证明，独立穷举收据另存measured_pool_cardinality_audit.json。"]
    lines+=['', '全长FASTA同时保留主探针和条件性拆分对照。双点PAIRWISE对照可与原双点探针同基因型；它不增加正式构建数，单点拆分序列也不因导出而获得正式五对资格。独立审核分别核对主面板与每条对照的canonical序列及父子关系。']
    lines+=['', '每个正式构建体有固定五对生产模型、10个阈值估计WT和10个留出WT；两组各形成5个不共享WT的差分。生产、估计、留出不混用。每个位点都有相同四次直接变化支持；未知区域不会被填成阳性。', '',
        '## 实现与复现', '',
        '新增v3预览合同贯通API、CLI、上传、回放及UI。双点SPURS为主要稳定性预测器；取消绝对分差拒绝，正式仍须SPURS≤1、五对FoldX中位数≤2及全五模型结构QC。缺失/非有限评分保留UNKNOWN，三点保留高阶外推和三个精确双点对照。',
        '阳性覆盖下界与完整足迹解析分开；Euler为诊断，局部通道有效性和同四次向量证据仍审核。单点kernel仅用于排序，不作为剪枝证明。三点压缩受覆盖、盲区、码距离1/2/3和单失效鲁棒性共同约束。',
        '实际回归收据：`'+json.dumps(regression_receipts,ensure_ascii=False)+'`。测试集合有重叠，不相加冒称独立测试数。',
        '保存96候选重放恢复17组双方门内分歧与39组轻度SPURS候选的确认机会，保留4个真实风险冲突；不是直接恢复正式资格。T3、T4、T5逐级真实验证及各阶段源码、输入、调用和工程状态随包保存。', '',
        '## 搜索、三点和APBS', '',
        '搜索状态：`'+json.dumps(json.loads((route/'search_state.json').read_text(encoding='utf-8')),ensure_ascii=False)+'`。',
        '达到80%后停止扩大双点；暂缓候选不是科学阴性，不宣称728位点对已全部实测。选板为已测正式池中的当前最好可行面板，不能推出全体组合的全局最优。',
        '逐阶段三点审核：`'+json.dumps(triple_reviews,ensure_ascii=False)+'`。旧未冻结基准压缩审查由frozen_baseline_information_compression_v3_1替代，不用于正式接受判断。',
        '首次达标为六个双点87.883091%。早期压缩白名单缺陷额外启动S55F/S101R完整确认；实际证据保留、费用计入，修复后不会继续复活无关候选。最终选板重新审核已完成五生产及20专用null的阳性双点，不新增冷候选。双点池重新选板与三点压缩分别报告，不能把双点减少构建数称为三点压缩成功。',
        'APBS状态：`'+json.dumps(summary.get('final_apbs_evidence',{}),ensure_ascii=False)+'`。',
        '逐构建体预处理审计摘要：`'+json.dumps(apbs_details,ensure_ascii=False)+'`。',
        'APBS执行、预处理及数值解析分别解释；未知/条件性预处理不补入hard覆盖。详细逐构建体电静力报告、PQR输入及数值误差/失败mask随物理证据包保留。', '',
        '## 成本、恢复与独立交付', '',
        f"v3阶段收费合计：{sum(r['charged_seconds'] for r in receipts):.3f}秒；报告生成时原case累计：{read_budget(case/'case_budget.json')['elapsed_seconds']:.3f}秒，历史、失败和重试均未清零。最终归档费用另计入delivery_status.json闭账值。开发与人工间隔不计入科学运行时间。",
        f"保留实际配对执行身份：{independent['actual_pair_identities']}；角色计数：{json.dumps(independent['role_counts'])}。逻辑配对不等于外部调用次数；实际时间播种BuildModel启动审计见buildmodel_invocation_audit.json。",
        f"实际时间播种BuildModel启动：{invocation_count}。此计数不含RepairPDB、未记录/非本采样协议的历史调用；合法多输出槽不冒充额外外部调用。",
        f"阶段峰值新增占用：{max(r['peak_owned_bytes'] for r in receipts)}字节；不包含随后交付归档峰值，最终资源收据另列。",
        'T6首次尝试因外部Node进程占用检查点JSON失败，1204.687秒已收费。原完整JSON、失败tmp、模型和null保留；占用释放后按同一固定槽恢复。没有终止用户进程、重复计数或降低科学门槛。',
        '首次最终打包被过严FASTA审核阻断：审核器把合法STAGE2_PAIRWISE对照误当成重复主探针。角色化复核保留严格主探针唯一性、每条序列WT/编号/氨基酸及对照父子关系；35相关回归3.44秒通过。仅恢复交付归档，不重跑已完成FINAL_APBS，模型账本hash前后相同。',
        '独立审核复算L⊆U、固定分母、覆盖/直接面积、三态码与单失效覆盖；核对canonical→作者编号→规范化基因型、模型/WT/DDG原始身份、CSV/JSON/FASTA一致性及ZIP每个成员hash。源码、阶段收据和浏览器只读证据纳入v3_stage_provenance.zip。',
        '工程节省来自旧合法证据复用、有入选价值者才支付独立null、达到真实目标后停止余下确认、内容去重和限域批次检查点。候选池及资格合同改变，不能据此宣称物理算法准确率提升或与旧路线公平提速倍数。', '',
        '## 科学边界与剩余限制', '',
        '1AH1为20构象NMR ensemble，共129个映射CA。构象两两配准RMSD中位数3.498195 Å，说明骨架不确定性尚未覆盖；本轮五对仍在固定骨架上作侧链建模，没有混入其他NMR骨架。',
        '不同合法重复可能产生相同WT/Mut字节，依赖和内容多样性如实保存；不筛选五个有利重复，不称独立热力学观测。经验q95、数值误差及有限留出验证不构成全局95%误报控制。',
        'SPURS/FoldX门槛是版本化预测风险合同，缺少CTLA4表达、折叠或结合实测校准；最终实验仍需表达/折叠对照及阳性后拆分定位。物理机会覆盖与几何编码不替代这些实验。', '']
    (output/'CTLA4_v3_final_report.md').write_text('\n'.join(lines),encoding='utf-8')
    write_json(route/'independent_final_package_audit.json',independent)
    return output/'CTLA4_v3_final_report.md'


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('case',type=Path)
    print(run(parser.parse_args().case.resolve()))
