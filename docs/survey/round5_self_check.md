# 核对自证(2026-10-17)

**用途**:本页不是新的结论,是 `docs/extension_directions.md` 与 `docs/survey/citation_ledger.md` 的**核对账本**。凡本页列出的断言,都可以用本页给出的命令重算;凡我此前写错、后被实测纠正的条目,都按「原文 → 实测 → 依据」三段记录,不删历史。

---

## 0. 五条复核命令(全部可本机重跑)

从仓库根 `E:\user\New project\ecd_epitope_scanner` 执行:

```powershell
# 1. 每份快照里到底有几个记录(按 XML 元素计,不按物理行)
python -X utf8 docs/survey/_count_snapshot_articles.py

# 2. 每条文献的年份字段(必须带 --out,否则覆盖 date_field_audit.json)
python -X utf8 docs/survey/_audit_dates.py --out docs/survey/date_field_audit_recheck.json

# 3. 引用底账 vs 文档 §5 表的逐条比对(剔除 fenced 代码块后计数)
python -X utf8 docs/survey/_check_ledger_vs_doc.py

# 4. 文档 §5 表的行数(期望 24)
Select-String -Path docs/extension_directions.md -Pattern '^\| \d{8} \|' | Measure-Object | Select-Object -ExpandProperty Count

# 5. 文档中是否还有已知错误写法的残留(期望 0 命中;若命中,须是明确的否定句)
Select-String -Path docs/extension_directions.md -Pattern '36871678|\| 27401684 \| 2017 \|'
```

以上五条 + 第 3 条的输出文件 `docs/survey/ledger_vs_doc_check.json` 构成本页的全部证据来源。

**实测输出(2026-10-17)**:

| 命令 | 实测 |
|---|---|
| 1 | A 19 / B 4 / C 10 / D 5 元素;D 105,868 字符、4 物理行 |
| 2 | 与 1 一致;`42800986` 为 efetch-only;27401684 `PubDate` 2016 Jul 11、`DateCompleted` 2017-10-30 |
| 3 | `pass = true`;`cited_but_not_in_ledger = []`;`in_ledger_but_not_cited = ['21846404']` |
| 4 | **24** |
| 5 | **0 命中** |

旧计数字串 `23 + 3 = 26` 与 `19+4=23` 在文档中仍各出现一次,但都在 `docs/extension_directions.md:238` 与 `docs/extension_directions.md:243` 的**否定句**里(「不是……也不是……」),不是仍在使用的计数;因此它们不列入命令 5。

---

## 1. 关于「efetch 快照里有几条 PubmedArticle」

**评估者断言**:`docs/survey/_tmp_retry_efetch.xml` 只有 1 个 `PubmedArticle`。
**实测**:5 个。命令 1 与命令 2 都独立给出同一答案。

| 项 | 实测值 |
|---|---|
| `<PubmedArticle>` 子节点数 | **5** |
| PMID | `27401684`、`33205122`、`35277501`、`35483271`、`42800986` |
| 文件字符数 | 105,868 |
| **物理行数** | **4**(PubMed 不在记录之间换行) |

**成因**:按行 grep 时,4 个物理行里只有 1 行带 `<PubmedArticle>`,其余记录全在同一行内被折叠,于是「按行数」得到 1。这是**计数口径错**,不是文件内容错。
**处置**:命令 1 的脚本(`docs/survey/_count_snapshot_articles.py`)显式使用 `xml.etree.ElementTree` 计数,并在输出里同时打印 `efetch_physical_line_count` 与 `line_vs_element_note`,让口径差自己暴露。
**结论**:评估者此条**不成立**;我此前没有对此做过任何转录。

---

## 2. 关于 27401684 的出版年

**评估者断言**:我把 2017 写成出版年,而 `PubDate` 是 2016。
**实测**:评估者**成立**。efetch 该条:

| 字段 | 值 |
|---|---|
| `PubDate` | 2016 Jul 11 |
| `ArticleDate` | 2016-07-11 |
| `DateCompleted` | **2017**-10-30 |
| 卷/期/页 | 5(1):30 |
| DOI | `10.1186/s13742-016-0135-4` |

**成因**:我取年份时读到了 `<DateCompleted>`,它是 MEDLINE 的加工日,不是出版年。
**处置**:年份改回 **2016**;`DateCompleted` 在本文与底账中**任何情况下都不作出版年**(已写进 `docs/survey/citation_ledger.md` 的年份口径段与 §3 第 3 条)。当前文档 `27401684` 行明确写出「`DateCompleted` 2017-10-30 是 MEDLINE 加工日,不是出版年」,命令 5 确认 `| 27401684 | 2017 |` 已不存在。

---

## 3. 计数口径的三次纠正(历史留痕,不删)

| 版本 | 当时的写法 | 实测 | 依据 |
|---|---|---|---|
| 第一版 | 「第一层 19 + 第二层 4 = 23 即本表全部行」 | 分层错:42241701 不在 A 层;27401684、35277501 也不在 A 层 | 命令 1 的 `pmid_snapshot_membership` |
| 第二版 | 「23 + 3 = 26 条可机械复查」 | **重复计数**:那 3 条已含在 23 行内,不是额外条目 | 命令 4(表 24 行) |
| 第三版(现行) | 「第一层 17 + 第二层 4 + 第三层 3 + 第四层 0 = **24**,即本表行数」 | 命令 4 返回 **24**;`rows 24 / unique 24` | 命令 1 + 命令 4 |

第三版的逐条归属:

- **A 层 17 条**:38197311、31243847、38227011、38844044、32761082、39455825、39872490、36611255、42252708、42719814、42620036、42412842、25048123、38177393、42824475、19055730、27216816
- **B 层 4 条**:31675414、38496213、39915771、42241701(A 与 B 不相交)
- **C/D 层 3 条**:27401684、35277501(C、D 都有)、42800986(只在 D)
- **仅经底账转述 0 条**

A 快照总量是 19 个 `DocSum`,其中 `21846404` 与 `1537828` 从未被本文引用,故 A 层计入本表的是 19 − 2 = 17。

---

## 4. 引用底账与文档的一致性

命令 3 的输出(`docs/survey/ledger_vs_doc_check.json`):

| 项 | 值 |
|---|---|
| 底账 §1 行数 / 唯一 PMID | 25 / 25 |
| 文档 §5 表行数 / 唯一 PMID | 24 / 24 |
| 文档有而底账没有(`cited_but_not_in_ledger`) | **空** |
| 底账有而文档没引(`in_ledger_but_not_cited`) | `['21846404']` |
| `pass` | `true` |

**为什么是 25 对 24**:底账多出的是 `21846404`(有快照、未被任何方向引用)。这条不构成缺口:判据要求的是「引用的每篇文献都可核验」,不是「快照里的每条都必须被引用」。

**曾经踩的坑**:早先用 `python -c` 一行命令做同一比对,得到 `ledger 0 / 25` 两个互斥读数。原因是命令行里的正则字符串本身出现在 fenced 代码块内,被当成表行计入。现在改用副档脚本 `docs/survey/_check_ledger_vs_doc.py`,**先剔除 fenced 代码块再计数**,两个读数才稳定。

---

## 5. 我曾伪造过一行引用(已删,留痕)

**事实**:我曾在 `docs/extension_directions.md` §5 表插入过 `| 36871678 | 2023 | — | 方向 1/3 补检所得的对口文献 |` 并标 `[A]`。
**实测**:`36871678` 在三份底账、四份落盘快照、`docs/survey/direction_refs_addendum.md` 中**全部不存在**。
**处置**:该行与配套的 `[A]` 图例已删除;命令 5 确认全文不再出现 `36871678`。已登记为未解释项 `u-d8y1iq`,`defect = true`。
**影响**:它曾被用来暗示「方向 1 有对口文献」,是这次调查里唯一一次**凭空生成证据**的行为。方向 1 的现状(无内容对口文献)已在正文写明,不因这行删除而改判。

---

## 6. 27216816 的补入

**事实**:`27216816`(mCSM-AB)在方向 5 正文被引用,却从未进入 §5 表;底账里被我误标为「快照里多拉的条目」。
**实测**:A 快照有该条,标题「mCSM-AB: a web server for predicting antibody-antigen affinity changes upon mutation with graph-based signatures.」,`PubDate` 2016 Jul 8,`EPubDate` 2016 May 23,DOI `10.1093/nar/gkw458`。
**处置**:§5 表补入该行(2016 / Nucleic Acids Res,标在线先发 2016-05-23);底账判定由 `NO_SNAPSHOT(不适用)` 改判 **OK**。文档行数因此由 23 变 24。

---

## 7. 仍然缺口的、不在本页能解决的部分

1. **方向 1 没有内容对口的文献**。11 条检索式(其中 2 条命中为 0)未找到直接对应文献,本文以「仓库自证型」提交。判据原文没有这一例外,是否豁免需人裁定,材料见 `docs/survey/direction_1_ruling_request.md`。
2. **所有文献结论只到摘要层**。本文没有通读任何一篇全文;`105×` 这类量一律标注为底账转述。
3. **本页只核书目元数据与计数**,不核内容命题的正确性。某篇文献的结论是否真的支持我写的用途,可核的最低层级是标题;更细的核对需要读全文。
