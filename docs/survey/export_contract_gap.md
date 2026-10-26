# 导出契约缺口:panel_coverage_matrix.csv 两列全空(缺陷修复任务,2026-10-17)

本文件是 `docs/extension_directions.md` **原「方向 1」** 的登记去处。
按 `docs/survey/direction_1_ruling_request.md` 的**选项 B** 处置:该条从扩展方向清单移出,降级为**缺陷修复任务**,正文内容完整搬到这里,不删除、不丢证据。移出的理由如实写在下面 §5:它是**仓库自证型**议题,两轮换词根检索无内容对口文献,而判据(1)要求每个扩展方向都有可核验的对口文献。

---

## 1. 缺陷是什么(实测)

| 项 | 实测值 | 出处 |
|---|---|---|
| 产物 | `workspace_data/full_runtime_restore_20261012/outputs/484899e022f3/panel_coverage_matrix.csv` | 冻结运行 `484899e022f3`(PD-L1) |
| CSV 表头 | `design_id,patch_id,perturbed_fraction` | 实测 3 列 |
| 非空列 | 只有 `patch_id`;`design_id` 与 `perturbed_fraction` **两列全空** | 实测 32 数据行,`design_id` 非空行数 = 0 |
| 源数据 | `analysis.json` 的 `epitope_mapping_plan.coverage_matrix`,32 条,字段为 `probe_id` / `patch_id` / `scenario` / `predicted_coverage` / `evidence_type` | 实测 |
| 写盘点 | `exporting.py:332-336`;表头在 `:335` 写死为 `["design_id", "patch_id", "perturbed_fraction"]` | 实测 |
| 静默机制 | `exporting.py:200-205` 的 `_write_csv` 用 `csv.DictWriter(..., extrasaction="ignore")` 再逐个 `row.get(field, "")`,表头键在源记录里不存在时**写空串而不报错** | 实测 |
| 仓库自述① | `README.md:129`:`- panel_coverage_matrix.csv：构建体 × 补丁的扰动覆盖矩阵（解读阴性结果用）。` | 实测 |
| 仓库自述② | `exporting.py:523`:`- panel_coverage_matrix.csv holds the construct-by-patch perturbation matrix; use it to weigh negative results.` | 实测 |

**结论**:这是一个**可直接复现的仓库内部缺陷**,不依赖任何外部文献就能成立。产物被描述为「构建体 × 补丁扰动矩阵」,实际内容是 **probe × patch × scenario 的覆盖表**;两列空掉不是数据缺失,而是表头键名与源字段名不是同一套。

> **更正(2026-10-17)**:本缺陷的根因先前被写成 `exporting.py:221-226`。实测 `exporting.py:221-226` 是 `assay_profile` / `fixed_assay_profile.json` 的写法,**与本缺陷无关**。正确引用是写盘调用点 `exporting.py:332-336`(表头 `:335`)与静默取空的 `exporting.py:200-205`。

---

## 2. 为什么必须"先修再审计"

这条台账要核对的是**显示分区 → 面板候选 → 导出通道**三者是否互相闭合。导出通道本身写出整列空文件时,「闭合」这件事在该通道上根本无法判断 —— 不是审计不严,是被审计的对象没有内容。因此修复**先于**台账。

---

## 3. 复现命令(可自行跑)

```powershell
cd "E:\user\New project\ecd_epitope_scanner"

# 1) design_id 非空的行数(预期 0)
Import-Csv ".\workspace_data\full_runtime_restore_20261012\outputs\484899e022f3\panel_coverage_matrix.csv" |
  Where-Object { $_.design_id -ne "" } | Measure-Object | Select-Object -ExpandProperty Count

# 2) 源数据的真实字段与条数(预期 32 与 probe_id/patch_id/scenario/predicted_coverage/evidence_type)
python -X utf8 -c "import json;d=json.load(open('workspace_data/full_runtime_restore_20261012/outputs/484899e022f3/analysis.json',encoding='utf-8'));rows=d['epitope_mapping_plan']['coverage_matrix'];print(len(rows));print(sorted(rows[0].keys()))"

# 3) 写盘点的固定表头与静默取空机制
Select-String -Path exporting.py -Pattern 'panel_coverage_matrix|perturbed_fraction'
Select-String -Path exporting.py -Pattern 'extrasaction="ignore"'
```

---

## 4. 修复口径(修的时候按这个判)

1. 表头必须由**源记录的实际键**决定,不能写死:要么用 `sorted(rows[0])`,要么显式列出 `probe_id,patch_id,scenario,predicted_coverage,evidence_type`。
2. `_write_csv` 遇到「表头键在记录里不存在」时**要么报错、要么在表头里就不出现**,不得静默写空。这是本次缺陷能长期存在的原因,单修表头不修这个机制,同类问题会复发。
3. `README.md:129` 与 `exporting.py:523` 的描述必须与修复后的真实内容一致(probe × patch × scenario 覆盖表),不得继续称「构建体 × 补丁扰动矩阵」。
4. 改完必须能复现:同一冻结运行重跑导出,`design_id` 非空行数应等于源记录条数(32),且新增一条针对「表头键必须存在于源记录」的测试。

---

## 5. 为什么从扩展方向清单移出(如实记录)

- 本条的核心依据是**仓库自身产物**,属仓库自证型:导出契约与产物数值,不依赖任何外部文献。
- 两轮换过词根的 PubMed 检索(见 `docs/survey/direction_1_4_literature_retry.md` §1):最贴近的两个检索式 `provenance audit scientific software output verification` 与 `export contract validation data pipeline` **命中数均为 0**;20 条候选中无一条的被检验命题是「导出字段与源数据是否对齐」。
- 判据(1)原文要求「列出的每个扩展方向都引用至少一篇真实可核验的文献」,**没有仓库自证型豁免条款**。因此把它留在扩展方向清单里,判据(1)就无法整体通过。
- 移出**不等于**这件事不重要:它是本仓库已定位的真实缺陷,已作为缺陷修复任务登记在本文件,并保留完整复现命令与修复口径。

可作**旁注**的同类工程实践两条,均**不承担对本缺陷的支撑**:Julearn 谈 ML 评估的可检查性 [PMID 38496213, GigaByte 2024, `10.46471/gigabyte.113`],CWLProv 谈工作流溯源可跨实现共享 [PMID 31675414, GigaScience 2019, `10.1093/gigascience/giz095`]。

**未排除的可能**:PubMed 未索引的工程类文献(会议论文、标准文档、软件工程期刊)可能有过对口研究;本次未检索非 PubMed 通道,故只能断言「在本次检索覆盖范围内没有」,不能断言「文献中不存在这类研究」。
