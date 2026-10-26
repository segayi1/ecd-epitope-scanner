# 引用底账(2026-10-16)

本文记录 `docs/extension_directions.md` 的**可机械核对层**:每一个被引用的 PMID,它在哪个落盘快照里、快照原文给出的年份/期刊/DOI、以及文档里的写法是否与之逐字一致。

**范围**:只覆盖有落盘快照的条目。无快照的条目一律标注 `NO_SNAPSHOT`,不在此列编造任何"已核"。

**年份口径(2026-10-16 第三轮定稿)**:出版年**一律取 PubMed 记录的 `PubDate`**(无 PubDate 时退回 `JournalIssue/Year`)。`DateCompleted` 是 MEDLINE 的加工日,**任何情况下都不作出版年**——本轮 R1 的错误正是把 `DateCompleted` 读成了出版年。唯一例外是 `19055730`,正文取期刊卷期年并按行内注明。
逐条实测值见 `docs/survey/date_field_audit.json`(由 `docs/survey/_audit_dates.py` 直接从四份落盘 XML 解析,不经人工转录;复核命令:`python docs/survey/_audit_dates.py`)。

## 0. 快照清单(实测)

| 代号 | 文件 | 请求 | docsum 数 | 实测 PMID 数 |
|---|---|---|---|---|
| A | `docs/survey/lit_resnapshot_20261016.xml` | `esummary.fcgi?db=pubmed&retmode=xml&id=<19 PMIDs>` | 19 `<DocSum>` | 19 |
| B | `docs/survey/lit_resnapshot_gap_20261016.xml` | `id=31675414,38496213,39915771,42241701` | 4 `<DocSum>` | 4 |
| C | `docs/survey/_tmp_retry_esummary2.xml` | 第二轮补检 esummary | 10 `<DocSum>` | 10 |
| D | `docs/survey/_tmp_retry_efetch.xml` | 第二轮补检 efetch(含摘要与 DOI) | — | 5 条 `PubmedArticle`(PMID `27401684`、`33205122`、`35277501`、`35483271`、`42800986`) |
| E | `docs/survey/_fetch_esummary_42800986.xml` | 第三轮单条补拉 `esummary.fcgi?db=pubmed&retmode=xml&id=42800986`(2026-10-17,HTTP 200) | 1 `<DocSum>` | 1 |

**计数口径(2026-10-17 定稿)**:上表「数」一律是 **XML 元素计数**(`<DocSum>` / `<PubmedArticle>` 子节点数),由 `docs/survey/_count_snapshot_articles.py` 用 `xml.etree.ElementTree` 实测,不经物理行数。
D 文件 105,868 字符但**只有 4 个物理行**(PubMed 不在记录之间换行),按行 grep 会把它数成 1 条 —— 这是 D 曾被误报「只有 1 个 PubmedArticle」的原因。复核命令:

```powershell
python -X utf8 docs/survey/_count_snapshot_articles.py
```

A 与 B **不相交**(`snapshot_accounting.json` 的 `snapshot_overlap = []`,已核)。四份快照共 34 个不同 PMID。

D 的 5 条中 4 条(`27401684`、`33205122`、`35277501`、`35483271`)在 C 的 esummary 里也有;`42800986` **不在 C 里**,但已于 2026-10-17 由快照 E 单独补拉一次 esummary 并与 D 互校 —— 标题、年、DOI、卷四项逐字一致。D 是完整的 `PubmedArticle` 记录(含年/期刊/DOI/卷期页),E 是同一条目的 esummary 视图,两者来自**不同的 E-utilities 接口**,故现在可以这么说:`42800986` 经两处独立拉取核对。比对脚本与结果:`docs/survey/_check_42800986_two_channels.py` → `docs/survey/_esummary_42800986_check.json`(`verdict = consistent`)。

> **解析口径提醒(2026-10-17)**:从 efetch 取 DOI **必须**限定在 `PubmedData/ArticleIdList` 下,不能用 `.//ArticleId` —— 后者会连 `ReferenceList` 里参考文献的 ArticleId 一起匹配,取到最后一条参考文献的 DOI。本项第一次运行就是这样把 `42800986` 的 DOI 读成 `10.1038/nature21406`(Nature 一篇参考文献),`verdict` 误报 `inconsistent`;改成定向取后两通道一致。

**只存在于 C/D 的条目**:`27401684`、`35277501`、`42800986` 见于 D 的 efetch;**`42800986` 不在 C 的 esummary 返回里**(实测 C 的 10 个 `Id` 为 38514204、27401684、32286507、35483271、33205122、30971826、41107899、35277501、29141660、31240330),它的标题/期刊/DOI 由 D 与 E 两处支撑。

## 1. 逐条对照表

**六列口径(2026-10-17 第四轮改定)**:`快照原文` 照抄快照;`文档现值` 是 2026-10-17 从 `docs/extension_directions.md` §5 表**实测抄出**的当前值,不是记忆中的值;`更正前写法(已过期)` 保留历史写法以便回溯,这一列**任何一条都不得被当作现值引用**;`判定` 只允许 OK / FIXED / NO_SNAPSHOT / UNSUPPORTED。

**本表与文档 §5 表的关系(2026-10-17 实测,第五次改定)**:文档 §5 表 **24 行,24 个唯一 PMID**;本表 §1 **24 行,24 个唯一 PMID**,两侧**完全一致**。两条「有快照但未被文档引用」的记录(`21846404`、`1537828`)改为**非表行记法**写在 §1 末尾,不计入行数 —— 先前 `21846404` 被写成表行,使台账比文档多一行,脚本读出的 `in_ledger_but_not_cited` 非空,与「本表只记被引条目」的口径自相矛盾。
文档中每一个被引 PMID 都在本表中有对应行(`cited_but_not_in_ledger = []`),24 条的**年份**(`year_mismatches = []`)、**期刊**(`journal_mismatches = []`,按 NLM 缩写规则匹配)与**标题片段**(`title_fragment_missing = []`)均逐条机械核对通过。
核对命令(数表行前先剔除 fenced 代码块,否则脚本文本里的 `| NNNNNNNN |` 模式会被当成表行 —— 2026-10-17 早先的命令就是这样把本表数成 0 行/25 行的):

```powershell
python -X utf8 docs/survey/_check_ledger_vs_doc.py
```

结果落在 `docs/survey/ledger_vs_doc_check.json`(`pass = true`)。该脚本除比 PMID 集合外,还逐条比对文档 §5 的年份与本表 `文档现值` 列的年份,差异写入 `year_mismatches`;`19055730` 是唯一登记在 `year_exceptions` 里的口径例外。`1537828` 是 1992 年 JBC 的 RNA helicase 论文,属检索噪声,**在本表中不用表行格式书写,不计入行数**。

| PMID | 快照 | 快照原文(年 / 期刊 / DOI) | 文档现值(2026-10-17 实测自 `docs/extension_directions.md`) | 更正前写法(已过期) | 判定 |
|---|---|---|---|---|---|
| 27401684 | D | **2016** Jul 11(PubDate;ArticleDate 同为 2016-07-11)/ GigaScience / 10.1186/s13742-016-0135-4(卷 5(1):30) | 2016 / GigaScience / 10.1186/s13742-016-0135-4 | 2017 / —(取自 `DateCompleted`) | **OK**(R1 曾误取 `DateCompleted` 2017-10-30,现已改回 2016,见 §2 R1) |
| 35277501 | D | 2022 / Scientific data / 10.1038/s41597-022-01163-2 | 2022 / Scientific Data / —(未列 DOI) | 同现值 | OK |
| 42800986 | D + E | 2026 / Array (New York, N.Y.) / 10.1016/j.array.2026.101150(卷 31) | 2026 / Array / —(未列 DOI) | 同现值 | **OK(D 与 E 两通道一致)**;比对见 `docs/survey/_esummary_42800986_check.json` |
| 42241701 | B | 2026 / Journal of medical Internet research / 10.2196/89071 | 2026 / J Med Internet Res / —(未列 DOI) | 摘要层措辞「已被当作标准做法」 | OK(标题已在 §5 逐字重写、该措辞已删) |
| 39915771 | B | 2025 / BMC biology / 10.1186/s12915-024-02101-x | 2025 / BMC Biol / —(未列 DOI) | 同现值 | OK |
| 31675414 | B | 2019 / GigaScience / 10.1093/gigascience/giz095 | 2019 / GigaScience / `10.1093/gigascience/giz095` | 同现值 | OK |
| 38496213 | B | 2024 / GigaByte (Hong Kong, China) / 10.46471/gigabyte.113 | 2024 / GigaByte / `10.46471/gigabyte.113` | 同现值 | OK |
| 42412842 | A | 2026 / Bioinformatics (Oxford, England) / 10.1093/bioinformatics/btag263 | 2026 / Bioinformatics / `10.1093/bioinformatics/btag263` | 同现值 | OK |
| 42620036 | A | 2026 Jul 31 / bioRxiv / 10.64898/2026.07.30.741792 | 2026 / bioRxiv(预印本)/ —(未列 DOI) | 同现值 | OK(**年份粗细差异已登记**:正文取年,快照为 2026 Jul 31) |
| 42252708 | A | 2026 / Protein science / 10.1002/pro.70655 | 2026 / Protein Sci / —(未列 DOI) | 同现值 | OK |
| 42719814 | A | 2026 / Bioinformatics advances / 10.1093/bioadv/vbag223 | 2026 / Bioinform Adv / —(未列 DOI) | 同现值 | OK |
| 38197311 | A | **2023** Nov 22 / Briefings in bioinformatics / 10.1093/bib/bbad491 | **2023** / Brief Bioinform / —(未列 DOI) | 2024 / Brief Bioinform | **OK(曾写 2024,已改 2023,见 §2 R2)** |
| 31243847 | A | **2019** Sep / Human mutation / 10.1002/humu.23852 | **2019** / Hum Mutat / —(未列 DOI) | 2020 / Hum Mutat | **OK(曾写 2020,已改 2019,见 §2 R3)** |
| 32761082 | A | 2021 Apr 20 / Bioinformatics / 10.1093/bioinformatics/btaa708 | 2021 / Bioinformatics / —(未列 DOI) | 同现值 | OK |
| 39455825 | A | 2024 Nov / Nature computational science / 10.1038/s43588-024-00716-2 | 2024 / Nat Comput Sci / —(未列 DOI) | 同现值 | OK |
| 39872490 | A | 2025 Jan 6 / Innovation (Cambridge (Mass.)) / 10.1016/j.xinn.2024.100750 | 2025 / Innovation (Camb) / —(未列 DOI) | 同现值 | OK |
| 38227011 | A | 2025 Mar / Human genetics / 10.1007/s00439-023-02623-4 | 2025 / Hum Genet / —(未列 DOI) | 同现值 | OK |
| 38844044 | A | 2024 Aug 15 / Journal of molecular biology / 10.1016/j.jmb.2024.168640 | 2024 / J Mol Biol / —(未列 DOI) | 同现值 | OK |
| 36611255 | A | 2023 Jan 19 / Briefings in bioinformatics / 10.1093/bib/bbac567 | 2023 / Brief Bioinform / —(未列 DOI) | 同现值 | OK |
| 25048123 | A | 2014 / Methods in molecular biology / 10.1007/978-1-4939-1115-8_8 | 2014 / Methods Mol Biol / —(未列 DOI) | 同现值 | OK |
| 38177393 | A | **2022** Dec / Nature computational science / 10.1038/s43588-022-00372-4 | **2022** / Nat Comput Sci / —(未列 DOI) | 2024 / Nat Comput Sci | **OK(曾写 2024,已改 2022,见 §2 R4)** |
| 42824475 | A | 2026 / Frontiers in immunology / 10.3389/fimmu.2026.1873135 | 2026 / Front Immunol / —(未列 DOI) | 同现值 | OK |
| 19055730 | A | 2008 Dec 2 / BMC bioinformatics / 10.1186/1471-2105-9-514 | **2009** / BMC Bioinformatics / —(未列 DOI) | 同现值 | **OK(唯一口径例外,见 §2 R5)** |
| 27216816 | A | 2016 Jul 8(PubDate;EPubDate 2016 May 23)/ Nucleic acids research / 10.1093/nar/gkw458 | 2016 / Nucleic Acids Res / —(未列 DOI) | 不在 §5 表中(2026-10-17 补入) | **OK**(此前被误记为本表「多拉的条目」;实为正文方向 5 已引用、§5 表却漏列,现已补齐) |

**以下两条有快照但文档从未引用,故不写成表行、不计入本表行数**(写成 `| NNNNNNNN |` 表行会让 `_check_ledger_vs_doc.py` 的 `ROW` 正则把它们算成台账行,从而与「本表只记被引条目」的口径自相矛盾):

- `21846404` — 快照 A;2011 Aug 17 / BMC bioinformatics / `10.1186/1471-2105-12-341`;标题「Prediction of conformational B-cell epitopes from 3D structures by random forests with a distance-based feature.」。文档未引用(主题相关但未被任何方向采用)。
- `1537828` — 快照 A;1992 Mar 5 / The Journal of biological chemistry / 无 DOI;标题「A new RNA helicase isolated from HeLa cells that catalytically translocates in the 3' to 5' direction.」。文档未引用(明显无关命中)。


## 2. 本轮修正记录

| 编号 | 对象 | 原写法 | 快照原文 | 处置 |
|---|---|---|---|---|
| R1 | 27401684 年份(**纠错,非修正**) | 我曾写 2017 | efetch `PubDate` 2016 Jul 11、`ArticleDate` 2016-07-11;2017 来自 `<DateCompleted>` 2017-10-30 | 年份**保持 2016**;2017 是 MEDLINE 加工日,不作出版年;DOI 补为 `10.1186/s13742-016-0135-4` |
| R2 | 38197311 年份 | 2024 | PubDate **2023** Nov 22 | 年份改 2023 |
| R3 | 31243847 年份 | 2020 | PubDate **2019** Sep | 年份改 2019 |
| R4 | 38177393 年份 | 2024 | PubDate **2022** Dec(EPubDate 2022 Dec 19) | 年份改 2022,并注明"在线先发于 2022-12" |
| R5 | 19055730 年份 | 2009 | PubDate **2008** Dec 2(= EPubDate,正式卷期 2009) | 保留 2009,行内标明这是本表唯一的口径例外 |

> R1–R5 的 `PubDate` 一律取快照**原文字符串**,不做任何语义重排。**R1 是纠错**:原写法 2017 不是快照给的出版年,而是加工日字段读错,故本条判定为 OK 而非 FIXED。R5 属"在线先发 vs 正式卷期"的已知差异,`benchmark_landscape.json` 早有记载,这里只加注不改数。

## 3. 仍然存在的限制

1. `NO_SNAPSHOT` 的条目(第二轮第一次检索命中、但只落了 esummary 而未落 efetch 的那批)**不在此表判定范围内**;文档正文凡引到它们的,已标注"未经二次拉取核对"。
2. 本表只核**书目元数据**,不核**内容命题**。内容层面的对口性判定在 `docs/survey/direction_1_4_literature_retry.md`,结论是:方向 4 有直接对口(27401684),方向 1 无。
3. R1 **不是**「esummary 与 efetch 两个接口不一致」:两处快照(`_tmp_retry_esummary2.xml` 的 PubDate「2016 Jul 11」与 `_tmp_retry_efetch.xml` 的 `<PubDate>` 2016 Jul 11)对该条都写 2016,并不矛盾。是我在取年份时读到了 efetch 的 `<DateCompleted>` 2017-10-30(Medline 加工日)并把它当成了出版年。该条已改回 2016;`DateCompleted` 在本台账所有取年份的地方一律不作出版年。逐条实测值见 `docs/survey/date_field_audit.json`。
