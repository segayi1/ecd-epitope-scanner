# 引用核对记录(reference reconciliation,2026-10-16)

**做的是什么**:把 `docs/extension_directions.md` 里引用的每一个 PubMed 编号,与两份落盘的 esummary 原始返回逐条对照,列出全部不一致。这份文件是那次对照的记录;机器可核版本在 `docs/survey/reference_reconciliation.json`。

## 1. 对照的两侧

| 侧 | 内容 | 规模 |
|---|---|---|
| 权威侧 | `docs/survey/lit_resnapshot_20261016.xml`(19 条 DocSum)+ `docs/survey/lit_resnapshot_gap_20261016.xml`(本表此前缺快照的 4 条) | **23 条 PMID**,每条都有原始 XML 返回 |
| 被检侧 | `docs/extension_directions.md` 全文出现的编号(正则 `\b\d{7,8}\b` 抽取) | **21 条 PMID** |

**集合差**:`doc_pmids_without_snapshot = []` —— 正文引用的每一个编号,现在都能落到可复查的拉取产物上。补拉前为 17/21(缺 `31675414`、`38496213`、`39915771`、`42241701`)。

## 2. 逐条不符清单(6 处)

判定口径:拿快照里的 `Title` / `PubDate` / `FullJournalName` / `DOI` 与正文写法逐字对比。**"改写"不算不符,"改变主张"或"改变文献类型"才算。**

| # | PMID | 正文原写法 | 快照实际 | 性质 | 处置 |
|---|---|---|---|---|---|
| 1 | 39872490 | 「Pythia 零样本 ΔΔG 与 105× 加速」 | "Structure-based self-supervised learning enables ultrafast protein stability prediction upon mutation." | **数值与关键词均不在标题中**。"105×" 只在底账摘要层转述里 | §5 表改为"标题无零样本/105×;加速比系底账转述,未回原文核",打 †;方向 5 已据 §6 限制 7 降级 |
| 2 | 42620036 | 「偏好有序区」 | "… a Bias for **Orderness**." | 同义改写,不改变主张 | 正文补注原文用词,表内标注**预印本** |
| 3 | 39455825 | §5 表「GeoDDG 的标注数据不足自述」 | "Improving the prediction of protein stability changes upon mutations by geometric learning and a **pre-training strategy**." | 标题主张是**方法**,"数据不足"是摘要层作者自陈 | 表内改为方法描述 + 注明自陈为摘要层,打 † |
| 4 | 42252708 | §5 表缩写作「SNAC-DB」 | "SNAC-DB: An ML-ready database for antibody and **NANOBODY® VHH**-antigen complexes…" | 丢掉**范围限定**(原文明文含纳米抗体) | 表内与方向 6 正文均补上范围 |
| 5 | 25048123 | 方向 5 写「AgAbDb 的条目也记着…」 | 快照确认为 *Methods Mol Biol* 的 AgAbDb 介绍章节 | **文献类型不符**:数据库介绍章节被当作独立评测研究 | 方向 5 补「引用强度说明」,该条降为范围说明、不承担方向 5 的立论 |
| 6 | 19055730 | §5 表「ElliPro 不连续表位基准集」 | "ElliPro: a new structure-based tool for the prediction of antibody epitopes." | **文献类型不符**:工具论文被当作基准集论文 | 表内改为"工具论文;其数据集常被用作基准" |

## 3. 同一编号跨文件的不一致(3 处)

不只正文与快照之间有差异,派生文件之间也互不一致:

| PMID | 各文件写法 |
|---|---|
| 42620036 | `lit_resnapshot_20261016_summary.json:127` 作 "Bias for **Ordered** Regions";快照 txt 与 addendum 作 "Bias for **Orderness**";`method_landscape.json` 为更短的截断写法 |
| 42252708 | 快照为全称含 NANOBODY® VHH;底账与正文均缩略 |
| 38197311 / 31243847 / 38177393 / 19055730 | 正文表记年份为 2024 / 2020 / 2024 / 2009,快照 `PubDate` 为 2023 Nov 22 / 2019 Sep / 2022 Dec / 2008 Dec 2 —— 前者用卷期正式年、后者用 PubMed 在线先发年,**两种记法都不假但混用使机械核对失效** |

## 4. 对结论的影响

- **第 1 处最相关**:方向 5 用 Pythia 的加速比论证"零样本 ΔΔG 可作快速筛选",该数值未回原文核,故方向 5 的**性价比论据降级为未核实**;方向 5 的归类(需外部数据)不因此改变,因为它靠的是输入契约与标签现状两条仓库自身事实。
- 第 2、3、4 处只改措辞强度,不改变任何方向的归属。
- 第 5、6 处改变的是**引用类型**:这两条不再被用来支撑"某研究实测了 X"这类主张。第 5 处影响方向 5 的立论负担(已转由 `docs/project_memory.md` 的标签现状承担)。
- 跨文件不一致(第 3 节)不改变本文任何结论,但意味着**任何以 `_summary.json` 为源的二次核对都会继承 "Ordered"/"Orderness" 这个错**,后续若自动化核对须以快照 XML 为准。

## 5. 本次核对没做到什么

1. 只核到 **esummary 层**(标题/年份/期刊/DOI),**没有**取全文,也没有重取摘要。"105×""标注数据不足"这类摘要层或正文层的量,本文一律标注为"底账转述,未回原文核",不得当作已验证数值。
2. 底账里 `cross_check.mismatches_found = 0` 的字段**没有被推翻**——它核的是各底账自身的范围,不覆盖本文件做的"正文 vs 快照"对照;两者不是同一件事,不能互相代替。
3. 编号之外的内容(例如方向 1 的 `exporting.py:221-226` 键名不匹配)属代码事实,不在本记录范围,见 `docs/extension_directions.md` 方向 1 与 `docs/survey/self_check.json`。
