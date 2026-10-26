# 台账修正记录(2026-10-16 第三轮)

本文件针对 `docs/survey/citation_ledger.md` 与 `docs/survey/extension_directions_reconciliation.json`
的三处错误写法,逐条给出**改前原文 / 改后原文 / 依据**。每条都附一条可直接运行的
PowerShell 复核命令;命令里的行号以本文件写作时的文件状态为准,若行号漂移,用其中的
关键词重新定位。

依据的统一来源是 `docs/survey/date_field_audit.json`(由 `docs/survey/_audit_dates.py`
直接从四份落盘 XML 解析生成,不经人工转录)。可随时重跑:

```powershell
python docs/survey/_audit_dates.py
```

---

## C1 — 27401684 的年份:2017 是错的,应为 2016

**改了哪里**:`docs/survey/citation_ledger.md:57`(§2 修正记录的 R1 行)与
`docs/survey/citation_ledger.md:26`(§1 对照表 27401684 行)。

**改前**(R1 行,原文):

```
| R1 | 27401684 年份 | 2016 | efetch `2017` | 改为 2017;DOI 补为 `10.1186/s13742-016-0135-4` |
```

**改后**:

```
| R1 | 27401684 年份(**本条为纠错,非修正**) | efetch `PubDate` 与 `ArticleDate` 均为 2016 | 我曾误取 `DateCompleted` 2017-10-30 当作出版年 | 年份**保持 2016**;2017 是 MEDLINE 加工日,不作出版年;DOI 补为 `10.1186/s13742-016-0135-4` |
```

**改前**(§1 对照表 27401684 行,原文):

```
| 27401684 | D | 2017 / GigaScience / 10.1186/s13742-016-0135-4 | 2016 / GigaScience / — | **FIXED**(见 §2 R1) |
```

**改后**:

```
| 27401684 | D | **2016** Jul 11 / GigaScience / 10.1186/s13742-016-0135-4(卷 5(1):30) | 2016 / GigaScience / 10.1186/s13742-016-0135-4 | **OK**(见 §2 R1 纠错说明) |
```

**依据**:`docs/survey/_tmp_retry_efetch.xml` 中该条 `<PubDate><Year>2016</Year><Month>Jul</Month><Day>11</Day></PubDate>`,
`<ArticleDate>` 同为 2016-07-11,而 `<DateCompleted>` 为 2017-10-30;
`docs/survey/_tmp_retry_esummary2.xml` 的 PubDate 同样是「2016 Jul 11」。
即两个接口一致说 2016,只有**加工日字段**是 2017。

**复核命令**:

```powershell
python -c "import json;a=json.load(open('docs/survey/date_field_audit.json',encoding='utf-8'))['by_pmid']['27401684']['efetch'];print(a['pub_date_raw'],'|',a['article_date_raw'],'|',a['date_completed_raw'])"
```

预期输出形如 `<Year>2016</Year><Month>Jul</Month><Day>11</Day> | <Year>2016</Year><Month>07</Month><Day>11</Day> | <Year>2017</Year><Month>10</Month><Day>30</Day>`。

---

## C2 — §3.3 的归因是错的:差异不来自「两个接口不一致」

**改了哪里**:`docs/survey/citation_ledger.md:69`(§3 待补 / 残留限制 第 3 条)。

**改前**(原文):

```
3. 三处 `FIXED` 中 R1 的年份差异(2016 vs 2017)是 esummary 与 efetch 两个接口对同一 DOI 的不同 `PubDate` 表达所导致,**不是任一方"编造"**;本文以 efetch 的 `<PubmedData>` DOI 为准并保留两处原文。
```

**改后**:

```
3. R1 不是「两个接口不一致」:两处快照(`_tmp_retry_esummary2.xml` 与 `_tmp_retry_efetch.xml`)对该条 `PubDate` 都写 **2016 Jul 11**,并不矛盾。是我在取年份时读到了 efetch 的 `<DateCompleted>` 2017-10-30(Medline 加工日)并把它当成了出版年。本条已改回 2016;`DateCompleted` 在本书所有取年份的地方一律不作出版年。
```

**依据**:`docs/survey/date_field_audit.json` 的 `by_pmid["27401684"]` 下,`esummary` 数组三条记录的 `pubdate_esummary`
与 `efetch.pub_date_raw` 均指向 2016;只有 `date_completed_raw` 是 2017。

**复核命令**:

```powershell
python -c "import json;a=json.load(open('docs/survey/date_field_audit.json',encoding='utf-8'))['by_pmid']['27401684'];print([e['pubdate_esummary'] for e in a['esummary']]);print(a['efetch']['pub_date_raw']);print(a['efetch']['date_completed_raw'])"
```

---

## C3 — D 文件的记录数与「5 条」口径

**现状**:`docs/survey/citation_ledger.md:14` 记 D 为「5 条 `PubmedArticle`」。实测
`docs/survey/_tmp_retry_efetch.xml` **确实是 5 条** `<PubmedArticle>`,PMID 为
`27401684`、`33205122`、`35277501`、`35483271`、`42800986`。这一行**数字正确,不改**,
但需要把「哪几条只靠 D 支撑」写清,以免被读成「5 条都有 esummary 互校」。

**改后**(在 §0 该表下方已有说明之外补一句):

```
D 的 5 条中,4 条(`27401684`、`33205122`、`35277501`、`35483271`)在 C 的 esummary 里也有;
**`42800986` 只在 D 里**,没有任何 esummary 互校。D 是完整 `PubmedArticle` 记录(含年/期刊/DOI/卷期页),证据强度不低于 esummary,但不可声称「经两处独立核对」。
```

**依据**:`docs/survey/date_field_audit.json` 的 `summary.efetch_only` 为 `["42800986"]`;
`per_file_counts` 给出 D=5、C=10 与各自的 PMID 列表。

**复核命令**:

```powershell
python -c "import json;r=json.load(open('docs/survey/date_field_audit.json',encoding='utf-8'));print(r['per_file_counts']['efetch_D_tmp_retry_efetch.xml']);print('efetch_only',r['summary']['efetch_only'])"
```

---

## C4 — 核对记录 JSON 里 27401684 的年份

**改了哪里**:`docs/survey/extension_directions_reconciliation.json:33`。

**改前**:该条 `found` / 备注里写「GigaScience 2016」而正文当时写 2017,两者不一致
(现正文已改回 2016,本条只需把 JSON 判据与正文对齐)。

**改后**:`27401684` 条目的 `found` 保留「GigaScience 2016」(与快照一致),并在同条目的
`status` 旁注明「正文已改回 2016;原写 2017 系误取 DateCompleted」。

**依据**:同 C1。

**复核命令**:

```powershell
python -c "import json;d=json.load(open('docs/survey/extension_directions_reconciliation.json',encoding='utf-8'));print(json.dumps(d,ensure_ascii=False)[:400])"
```

---

## 未变的既有事实

- §0 的 A=19、C=10、并集 21 条被正文引用(23 行表 = 21 被引 + `[A2]` 两条 35277501/42800986 中已在并集外的……实际以
  `docs/survey/date_field_audit.json` 的 `summary.distinct_pmids` = 34 与 §5 表 23 行为准,两者口径不同:34 是**全部落盘快照**的并集,23 是**本文实际引用**的行数)。
- `docs/survey/direction_1_4_literature_retry.md` 里「请求 5 个 PMID、返回 5 条」的记录与 D 文件实测一致,无需改。
