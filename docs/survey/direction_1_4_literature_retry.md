# 方向 1 / 方向 4 的文献补检(第二轮,2026-10-16)

目的:独立评估者指出判据(1)"每个方向至少引一篇可核验文献"对方向 1 与方向 4 未满足,而**判据原文没有"仓库自证型可豁免"的例外**。本文档记录本轮的补救检索:换词根重检、对被引条目全部做 esummary 二次拉取并落盘、并逐条判断"直接对口 / 仅机制相近"。

本文件只记录检索与核对事实。是否把这两个方向判定为"不适用",仍属需要人裁决的事项(见 §5)。

---

## 1. 本轮检索式与命中数

全部经 PubMed E-utilities `esearch.fcgi`(`retmode=json`),`retmax=5~6`。原始返回存于 `docs/survey/_tmp_retry_esearch.json`。

| # | 检索式(原文) | 命中数 | 用于 |
|---|---|---|---|
| R1 | `reproducibility workflow identical results` | 143 | 方向 4 |
| R2 | `bitwise reproducibility scientific software` | **0** | 方向 4 |
| R3 | `provenance audit scientific software output verification` | **0** | 方向 1 |
| R4 | `export contract validation data pipeline` | **0** | 方向 1 |
| R5 | `computational reproducibility analysis pipelines research` | 1670 | 方向 4(过宽,只取前 6 条) |
| R6 | `same software same data results differ` | 5775 | 方向 4(过宽,只取前 6 条) |
| R7 | `regression test detection of output change` | 65 | 方向 1 |
| R8 | `reproducible research computational results cannot be reproduced` | 38 | 方向 4 |
| R9 | `software reproducibility crisis computational` | 34 | 方向 4 |
| R10 | `verifiable reproducible bioinformatics framework` | 74 | 方向 4 |
| R11 | `reproducibility of published computational results` | 4893 | 方向 4(过宽) |

**观察**:R3 与 R4 这两个最贴近方向 1(导出字段与源数据是否对齐、产物契约是否被验证)的检索式**返回 0**;方向 1 找不到任何直接对口的文献,与本轮之前一轮的结论一致。

---

## 2. 二次拉取(全部落盘)

| 文件 | 内容 | 请求 |
|---|---|---|
| `docs/survey/_tmp_retry_esummary.xml` | 10 条 DocSum(第一轮候选) | `esummary.fcgi?db=pubmed&retmode=xml&id=42861721,42801184,42800986,42796516,42787826,42850345,42656234,42440088,42430271,42322893` |
| `docs/survey/_tmp_retry_esummary2.xml` | 10 条 DocSum(第二轮候选) | `esummary.fcgi?db=pubmed&retmode=xml&id=38514204,27401684,32286507,35483271,33205122,30971826,41107899,35277501,29141660,31240330` |
| `docs/survey/_tmp_retry_efetch.xml` | 5 条摘要原文 | `efetch.fcgi?db=pubmed&retmode=xml&id=27401684,33205122,35277501,35483271,42800986` |

三份文件均为服务端原始返回,未做删改。

**第一轮候选的处置**:全部淘汰。`42861721`(CAD-score-LT)是结构比较方法、`42800986`(FairFlow)标题对口但**未取摘要**(当时只取到 esummary 标题层,全文与摘要都没有拉取;它不在「R1 前 6 之外」这个说法所指的范围内 —— 那份排序属第二轮候选表,本句此前把两轮写混了)、其余 8 条为医学影像、脊柱外科、视网膜影像、呼吸阈值、心肾、放疗质控等无关领域。

> 更正(2026-10-17):本句原写「排在 R1 前 6 之外未取摘要」,与本文件 `:35` 记「`_tmp_retry_esummary.xml` 含 42800986」自相矛盾。
> 实测:该文件确有 10 条 DocSum,42800986 是第 3 条,标题为 "FairFlow: A transparency-first framework for verif…",即**元数据已取到、摘要未取**。
> 两处现在说法一致:42800986 在第一轮候选里、被淘汰的原因是没有内容对口证据(只有标题层),不是「没取到」。它的元数据后来由第二轮 efetch(快照 D)与单独补拉的 esummary(快照 E)双通道补齐。

**第二轮候选的元数据(来自 `_tmp_retry_esummary2.xml`)**:

| PMID | PubDate | 期刊 | 标题(截断) |
|---|---|---|---|
| 27401684 | 2016 Jul 11 | GigaScience | Tools and techniques for computational reproducibility. |
| 33205122 | 2020 Aug 14 | Patterns (N Y) | dtoolAI: Reproducibility for Deep Learning. |
| 35277501 | 2022 Mar 11 | Scientific data | FAIRly big: A framework for computationally reproducible processing of large-scale data. |
| 35483271 | 2022 Jun | Comput Methods Programs Biomed | ClinicaDL: An open-source deep learning software for reproducible neuroimaging processing. |
| 42800986 | 2026 Sep | Array (N Y) | FairFlow: A transparency-first framework for verifiable and reproducible bioinformatics. |

---

## 3. 逐条判"直接对口 / 仅机制相近"

判断标准(本轮事先写下):该文献陈述的**被检验命题**必须与方向的对象至少同构才算直接对口。
- 方向 1 的对象:**一份计算产物的导出字段与其源数据是否对齐,内部不变量是否闭合**。
- 方向 4 的对象:**同一份冻结输入、同一份代码,两次独立执行的结果集合是否逐条相等**。

| PMID | 摘要中与本方向相关的原文(逐字) | 判定 | 理由 |
|---|---|---|---|
| 27401684 | "The deterministic nature of most computer programs means that the same analysis tasks, applied to the same data, will often produce the same outputs. However, in practice, computational findings often cannot be reproduced because of complexities in how software is packaged, installed, and executed-and because of limitations associated with how scientists document and ..." | **直接对口(方向 4)** | 它把"同代码 + 同数据 → 同输出"当作可检验的默认预期,并指出实践中失败的来源。这正是方向 4 的命题(集合相等)及其替代说明,属同一可检验命题。 |
| 35277501 | "It affords the capture of machine-actionable computational provenance records that can be used to retrace and **verify the origins of research outcomes**, as well as be **re-executed independent of the original computing infrastructure**." | **直接对口(方向 4,且部分覆盖方向 1)** | "verify the origins of research outcomes"与方向 1 的"导出字段能否追溯到源数据"同构;"re-executed independent of the original infrastructure"与方向 4 同构。 |
| 42800986 | "Because tools are updated to fix bugs and add features, **identical code can produce different results depending on which versions are installed**... we present FairFlow, a framework that makes rigorous, reproducible analysis accessible..." | **直接对口(方向 4)** | 它正面陈述"同代码 + 版本差异 → 结果不同",与方向 4 的检验对象同构;且是标题就自述为 "verifiable and reproducible bioinformatics" 的框架论文。 |
| 35483271 | "it checks the absence of data leakage when inferring the results of new data with trained networks, and saves all necessary information to guarantee the reproducibility of results" | 仅机制相近 | 对象是深度学习推理的泄漏检查与结果复现,**不是**导出字段与源数据的对齐审计。 |
| 33205122 | "Model outputs are critically dependent on the data and processing approach used to initially generate the model, but this provenance information is usually lost during model training." | 仅机制相近 | 关注模型训练阶段的溯源丢失,不是产物导出契约。 |

---

## 4. 结论

- **方向 4 现有一篇内容直接对口、且经二次拉取落盘的文献:27401684(GigaScience 2016)**。辅助:`35277501`、`42800986`。三者都来自本轮真实拉取,原始返回见 §2。
- **方向 1 仍然没有直接对口文献**:R3、R4 两个最贴近的检索式返回 0,本轮 20 条候选中无一条的对象是"导出字段与源数据对齐"。如实记录,不以机制相近充数。

**重要限定**:上表"摘要中与本方向相关的原文"一栏全部取自 `_tmp_retry_efetch.xml` 的摘要**原文**,不是底账转述;但**本文档未通读全文**,故这些文献对该方向的支持强度只到"摘要层陈述的命题与本方向同构"。

---

## 5. 仍需人裁决的事项

判据(1)对方向 1 的适用性。三种可选处置:

1. **接受现状**:方向 1 保留在交付物中,并**明确标注**它不满足判据(1)的文献要求。这意味着判据(1)整体上不能判为 pass。
2. **由人书面裁定方向 1 不适用该判据**(理由:其对象是纯工程闭环,要求外部文献等于要求它改变自身性质),裁定记录写入 `docs/extension_directions.md`,则判据(1)可对剩余 8 个方向判。
3. **把方向 1 从"扩展方向"降级为"已知缺陷的修复任务"**,从判据(1)的适用范围中移出。这样交付物仍有 8 个方向,每个都有直接对口文献。

我不自行选定其中任何一项。
