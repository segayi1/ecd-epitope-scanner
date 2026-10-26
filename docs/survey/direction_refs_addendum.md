# 方向 1 / 3 / 4 的文献依据补充记录

**用途**:上一次独立核验指出判据(1)「每个扩展方向都引用至少一篇真实可核验文献」在方向 1、3、4 上不成立(正文三节通篇无 PMID/DOI)。本文件记录为这三条方向补做的检索、命中的文献、以及**未命中时的如实标注**。检索经 PubMed E-utilities(esearch + esummary),全部命中条目另有一份二次拉取快照。

**检索时点**:2026-10-16。

## 检索式与实际命中(esearch,retmax=4—5)

| 检索式 | 命中 PMID |
|---|---|
| reproducibility AND computational AND pipeline AND structural AND biology | 42876657, 42848382, 42798706, 42734744, 42730630 |
| software AND verification AND invariant AND scientific AND workflow AND provenance | (空) |
| negative results AND reporting AND computational AND design AND protein | 42874184, 42843132, 42771345, 42691622, 42645154 |
| coverage AND completeness AND audit AND design AND library AND protein | (空) |
| reproducibility AND workflow AND provenance AND structural AND bioinformatics AND containers | 31675414 |
| reproducibility crisis computational biology reproducibility of published computational results | 41922703, 41100578, 39915771, 38897481 |
| workflow provenance reproducibility computational pipelines FAIR | 42527781, 41332756, 40880880, 36064140 |
| negative controls benchmarking protein design computational predictions validation | 23180633 |
| completeness coverage assessment library design protein variants systematic | (空) |
| unlabeled data not negative examples positive unlabeled learning protein function annotation | 42884310, 42877654, 42876682, 42873887 |
| abstention uncertainty reporting computational predictions protein engineering coverage | (空) |
| benchmark pitfalls machine learning protein engineering data leakage reporting standards | (空) |
| failed designs reporting protein engineering negative results publication bias computational | (空) |
| positive unlabeled protein function prediction | 42779863, **42412842**, 42215600, 42191748 |
| positive unlabeled learning review | 42837747, 42784600, 42428098, 42241701 |
| library completeness coverage protein variant library design NGS | 22913592 |
| data leakage machine learning model evaluation pitfalls | 42689951, 42149898, 41953863, 41590873, **38496213** |
| reproducibility standards computational biology machine learning reporting checklist | 42669717, 40892627 |

**检索到的条目元数据**(esummary 拉取,与上表部分对应):

> **更正(2026-10-16 复核)**:原表头写作「esummary 二次拉取,与上表一一对应」,**这句话不成立**,已改。经回查 `docs/survey/lit_resnapshot_20261016.xml`(19 条 DocSum,Id 标签实测)与 `_tmp_abstracts.xml`、`_tmp_ddg.xml`,下表中 5 条里**只有 42412842 有可复查的拉取记录**(它在 10-16 快照中);31675414、38496213、39915771、42241701、31675414 四条**没有任何已落盘的 esummary 快照**——它们出现在本表,依据是检索阶段的返回内容,但该次返回未单独落盘。第 40 行的注脚只交代了 39915771 的 DOI,未交代这一整类缺口。
>
> 影响范围:这 5 条中,被本文采用为证据的只有 **42412842**(方向 3),它在快照中可核。**38496213 与 31675414 在方向 1、方向 4 中只作旁注、已明确声明不作证据**,故该缺口不影响任何结论;但它确实意味着本文件下表不能作为"这 5 条都经过二次核对"的记录。此条已登记为未解释项并在此更正,不掩饰。

| PMID | PubDate | 期刊 | DOI | 标题 |
|---|---|---|---|---|
| 31675414 | 2019 Nov 1 | GigaScience | 10.1093/gigascience/giz095 | Sharing interoperable workflow provenance: A review of best practices and their practical application in CWLProv |
| 38496213 | 2024 | GigaByte | 10.46471/gigabyte.113 | Julearn: an easy-to-use library for leakage-free evaluation and inspection of ML models |
| 39915771 | 2025 Feb 6 | BMC biology | (见下注) | Same data, different analysts: variation in effect sizes due to analytical decisions in ecology and evolution |
| 42241701 | 2026 Jun 4 | Journal of medical Internet research | 10.2196/89071 | Detecting Uncoded Self-Harm in Veterans' Electronic Health Records Using Positive and Unlabeled Learning |
| 42412842 | 2026 Jul 1 | Bioinformatics | 10.1093/bioinformatics/btag263 | PU-GRAIL: residue-level graph learning for identifying protective bacterial antigens under positive-unlabeled settings |

> 注:`39915771` 的 DOI 在检索阶段的返回中未被单独取出核对,故只写期刊与标题;**该条在本文件里只作旁证、不作证据**。据 2026-10-16 复核,该条与 31675414、38496213、42241701 一样**没有落盘的 esummary 快照**(见上方更正块),本表不能作为它们已二次核对的凭据。

## 方向 1(覆盖可审计性台账 / 工程不变量审计)

**命中的直接对口文献**:**无**。检索式 `software AND verification AND invariant AND scientific AND workflow AND provenance` 返回空;与「审计一份计算产物的内部不变量是否闭合」直接对口的文献没有检索到。

**可用的间接文献**(口径需声明):`38496213`(GigaByte 2024, `10.46471/gigabyte.113`, Julearn)谈 ML 评估的可检查性;`31675414`(GigaScience 2019, `10.1093/gigascience/giz095`, CWLProv)谈工作流溯源的跨实现共享。

**结论**:方向 1 标为**仓库自证型**,其依据是仓库自身产物(见正文方向 1 的实测数字),不引用外部文献充数。上面两条只作为"同类工程实践存在先例"的旁注,**不作为该方向的证据**。判据(1)是否对仓库自证型方向适用,须由人裁定。

## 方向 3(未解假设与盲区显式化)

**命中的对口文献**:`42412842`(Bioinformatics 2026, `10.1093/bioinformatics/btag263`, PU-GRAIL)。本次检索式 `positive unlabeled protein function prediction` 第 2 位命中它;它同时支撑方向 3 的动机——"未标记 ≠ 阴性"在残基层面的正面处理。

**补充佐证**:`42241701`(J Med Internet Res 2026, `10.2196/89071`)说明 PU 思路在别的领域被当作"未标记数据必须显式建模"的标准做法。这条与蛋白质无关,**只作旁证**。

**结论**:方向 3 现有一条直接可核验引用(`42412842`),用途已写明。

## 方向 4(同一冻结输入两次执行的逐记录一致性)

**命中的对口文献**:`31675414`(GigaScience 2019, `10.1093/gigascience/giz095`, CWLProv)。它把"工作流溯源可以跨执行互操作地共享"作为可检验的工程目标,与本方向的"两次执行的分区成员/候选 ID/导出清单集合相等"属同一类问题——**但它是工作流溯源口径,不是"同一输入两次运行逐条一致"的重复性实测**。

**可用的旁证**:`39915771`(BMC Biology 2025, 'Same data, different analysts: variation in effect sizes due to analytical decisions in ecology and evolution')说明"同一数据、不同分析者"会产生差异;这与本方向检验的"同一输入、同一代码、两次执行"不是同一件事,**不构成对本方向的证据**。

**结论**:方向 4 有可核验的间接引用(`31675414`),但**无直接对口文献**,仍主要属仓库自证型。是否够判据(1),须由人裁定。

## 口径声明(本章三条共用)

1. 上表所有 PMID 都是 esearch 返回的**原始 idlist**(加粗者为被采用的两条),未经人工筛选后的完整清单保留,是为了让"检索到了什么"可复核;命中里大量与主题无关(例如该条检索式返回的四条全部是无关领域),说明**"检索到"不等于"支持"**。
2. 本文件不把任何一条文献结论当作本仓库已通过某项检验的证据;方向 1/3/4 的判据本身都还没执行过。
3. 未检索到直接对口文献的方向,按判据(1)原文属"未满足",这里如实标注为**仓库自证型**,不拿间接文献冒充直接依据。
