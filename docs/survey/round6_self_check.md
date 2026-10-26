# 第六轮自查留痕(2026-10-17)

本文件是本轮 v1–v5 的**逐处修改留痕**:改前原文、改后原文、依据文件与字段、可重跑命令。
写法与 `docs/survey/round5_self_check.md` 一致(三段式:改了什么 / 凭什么 / 怎么自己复核)。
**本文件不宣称任何判据已通过**;判据(1)的裁定仍未作出,见 `docs/survey/direction_1_ruling_request.md`。

---

## 0. 可重跑命令与 2026-10-17 实测输出

在 `E:\user\New project\ecd_epitope_scanner` 下依次运行(全部为只读或写 `docs/survey/` 下的报告):

| # | 命令 | 2026-10-17 实测输出 |
|---|---|---|
| 1 | `python -X utf8 docs/survey/_count_snapshot_articles.py` | `A_docsum: 19` / `B_docsum_gap: 4` / `C_docsum_retry: 10` / `D_efetch_retry: 5` / `E_esummary_42800986: 1`;`efetch elements: 5 ['27401684','33205122','35277501','35483271','42800986']`;`distinct PMIDs across snapshots: 34`;`cited PMIDs in doc: 24 \| without any snapshot: []`;`D_efetch_retry: 105868 chars, 4 physical lines` |
| 2 | `python -X utf8 docs/survey/_check_ledger_vs_doc.py` | `ledger_section1_rows: 25` / `doc_pmid_rows: 24` / `cited_but_not_in_ledger: []` / `in_ledger_but_not_cited: ['21846404']` / `year_checked_count: 24` / `year_mismatches: []` / `pass: True` |
| 3 | `python -X utf8 docs/survey/_fix_reconciliation_counts.py` | `table rows: 24` / `layer breakdown: {'A': 17, 'B': 4, 'C/D/E': 3, 'none': 0}` / `document lines/bytes: 277 35486` |
| 4 | `python -X utf8 docs/survey/_check_42800986_two_channels.py` | `_esummary_42800986_check.json` 的 `verdict: consistent`(title / year / doi / volume 四项全等) |
| 5 | `Select-String -Path docs/extension_directions.md -Pattern '105×\|零样本' -Encoding UTF8` | 仅 4 行:`:202`(§5 表,声明已撤)、`:267`(§6 限制表第 1 项)、`:274`(影响段)、`:276`(残留限制);无一处把该量当论据 |
| 6 | `Select-String -Path docs/extension_directions.md -Pattern '^\| \d{8} \|' -Encoding UTF8` | **24** 行 |
| 7 | `python -X utf8 docs/survey/_audit_dates.py --out docs/survey/date_field_audit_recheck.json` | 年份审计重跑;**必须带 `--out`**,不带会拒绝覆盖既有报告(退出码 1) |

> 陷阱提醒:`_audit_dates.py` 的默认目标 `docs/survey/date_field_audit.json` 已存在,直接重跑会以退出码 1 拒绝。
> PowerShell 不支持 `&&`,也不支持 heredoc;python 打印中文须 `-X utf8`(默认 cp1252)。

---

## 1. v1 —— 42800986 补第二条通道(esummary)

**改了什么**

1. 新增落盘快照 **E**:`docs/survey/_fetch_esummary_42800986.xml`(3,325 B,1 个 `DocSum`,`PubDate` 2026 Sep,`Array (New York, N.Y.)`,DOI `10.1016/j.array.2026.101150`,卷 31)。
2. 新增比对脚本 `docs/survey/_check_42800986_two_channels.py` → `docs/survey/_esummary_42800986_check.json`。
3. `docs/extension_directions.md` §5 表 42800986 行:原写「仅 `_tmp_retry_efetch.xml` 一条记录支撑,无 esummary 快照互校」→ 改为「D + E 两条通道互校」。
4. `docs/survey/citation_ledger.md` §0 新增快照 E 行;§0b 记 42800986 为「D+E 双通道」;§1 该行判定改 `OK(D+E)`。

**凭什么**:评估者第五轮驳回缺口第 2 条——该条只有 efetch 单通道。

**怎么自己复核**:命令 #1(应看到 `E_esummary_42800986: 1`、`42800986: ['D_efetch_retry', 'E_esummary_42800986']`)与命令 #4。

**踩过的坑(已记 `u-85hv5l`,defect=true)**:从 efetch 取 DOI 必须限定 `PubmedData/ArticleIdList`;用 `.//ArticleId` 会把 `ReferenceList` 里参考文献的 DOI 读进来,首次把 42800986 读成 `10.1038/nature21406`,误报 inconsistent。

---

## 2. v2 —— 台账 §1 六列化,过期年份移出「文档写法」列

**改了什么**

1. `docs/survey/citation_ledger.md` §1 表头改为六列:
   `| PMID | 快照 | 快照原文(年 / 期刊 / DOI) | 文档现值(2026-10-17 实测自 docs/extension_directions.md) | 更正前写法(已过期) | 判定 |`
2. 25 行全部重写。旧年份(2024 / 2020 / 2024 等)只留在「更正前写法(已过期)」列,并注明是哪一天改的:
   38197311→2023、31243847→2019、38177393→2022、19055730→2009(唯一口径例外)、27401684→2016。
3. `docs/survey/_check_ledger_vs_doc.py` 新增 `row_cells` / `lead_year` / `doc_years` 与年份比对:
   `pass = not (doc_set - ledger_set) and not year_mismatches`;
   `ACCEPTED_YEAR_EXCEPTIONS = {"19055730": "快照 PubDate 2008-12-02,本文按期刊卷期年取 2009"}`。

**凭什么**:评估者第五轮驳回缺口第 3 条——台账「文档写法」列与被核文档给出两个年份。

**怎么自己复核**:命令 #2,应看到 `year_checked_count: 24`、`year_mismatches: []`、`pass: True`。

**残留**:`in_ledger_but_not_cited: ['21846404']` —— 台账里保留了这条快照命中但正文未引用的记录,判定为「保留但不计入 24 行」,不构成不一致。

---

## 3. v3 —— 用脚本重算两处计数(不用手改)

**改了什么**

1. `docs/survey/_count_snapshot_articles.py`:加入快照 E;`cross_file_counts` 改为按 label 前缀自动生成;新增 `distinct_pmids_breakdown`(`cited_in_document` / `in_snapshot_but_not_cited` / 两个清单 / `sums_to_distinct_pmids`)。
2. `docs/survey/snapshot_article_count.json` 重生成:`distinct_pmids 34 = cited_in_document 24 + in_snapshot_but_not_cited 10`,`sums_to_distinct_pmids: true`;未引用的 10 条为 `1537828, 21846404, 29141660, 30971826, 31240330, 32286507, 33205122, 35483271, 38514204, 41107899`。
3. 新增 `docs/survey/_fix_reconciliation_counts.py`(从文档 + `snapshot_article_count.json` 重算),改写 `docs/survey/extension_directions_reconciliation.json` 的 `counts` 为 `table_rows_after: 24` / `mechanically_checkable_after: 24`,新增 `counting_note`、`layer_breakdown`(`A 17 / B 4 / C_D_E 3 / no_snapshot 0` / `sums_to_table_rows: 24`)、`layer_pmids`、`document_lines`、`document_bytes`。
4. 新增 `docs/survey/_strip_stale_counts.py` 清掉 c06/c07 与 `counting_note` 里过时的字面计数。
5. `docs/extension_directions.md` 第三层标题改为「**C/D/E 快照(3 条,`_tmp_retry_esummary2.xml` + `_tmp_retry_efetch.xml` + `_fetch_esummary_42800986.xml`)**」,并写明 42800986 由 D 与 E 两处互校。

**凭什么**:评估者第五轮驳回缺口第 3 条(计数)与本轮未解释项 u2(`distinct_pmids=34` 易被误读)、u6(counts 仍记 23/26)。

**怎么自己复核**:命令 #1 与 #3。关键不变量三条:
- `A∪B∪C∪D∪E` 的**元素级**计数 = 19+4+10+5+1,但**去重后** 34,两者不是同一个量;
- 24 + 10 = 34,必须成立;
- 表行数 = 分层和 = 24,不存在「24 + 3」。

---

## 4. v4 —— 让方向 5 的立论不再依赖未核的「105×」

**改了什么(共 4 处文档编辑)**

1. §5 表 39872490 行(现 `docs/extension_directions.md:202`):「标题无"零样本"/"105×";加速比系底账转述,未回原文核」→「标题只支持 'ultrafast' 这一层;无"零样本"、无"105×","105×"仅存于底账转述并已从本文删除」。
2. §5 表后「数值层级」条(`:241`):补一句限定——撤下后仍留字样的地方只有 §6 限制表与本行,属「记录下来并声明不可用」的元陈述。
3. §6 限制 7 表第 1 项(`:267`):「§5 表与 §2 曾写…」→「§5 表曾写…(2026-10-17 已**删除**)」,并把「来源」写准为底账 `mutation_effect_landscape.json` 的 `claims`。
4. `:274` 影响段:原写「2026-10-17 已从正文与 §5 表删除」——但正文(方向 5,`:123-131`)本来就没有这句,属**过度声明**;已改为「已从 §5 表删除」,并写明 §2 方向 5 的位置 `docs/extension_directions.md:123`,立论改由「本仓库自身的标签现状 + PU-GRAIL 的方法同构性」承担。

**凭什么**:评估者第五轮驳回缺口第 4 条(「105×」在 §5 表 39872490 行与方向正文仍有残留)。

**怎么自己复核**:命令 #5。剩余 4 行(202 / 267 / 274 / 276)逐行读过,均为「该量来自底账转述、已撤下、不得当已验证数值」的记录或自我声明,没有一处用它支撑任何方向。

**未做且有意不做的**:底账 `docs/survey/mutation_effect_landscape.json` 的 `claims` 里「计算速度提升 105 倍」的转述**保留不动**——它是底账对文献摘要的记录,不是本文论据;删除底账会抹掉「我当初从哪里读到这个数」的线索。

---

## 5. v5 —— 方向 1 裁定单重写成「一题两选」

**改了什么**:`docs/survey/direction_1_ruling_request.md` 整篇重写(8,810 B / 114 行)。

1. 删掉上一版的「A / B / 都不选」三选一与我的推荐。
2. 只留一题:方向 1 留作扩展方向(判据(1)未通过)还是移出、降级为缺陷修复任务(对剩下 8 个方向通过)。
3. 把「选项 B 会改变判据口径、须书面确认」「未确认前不执行 B」写进选项本身。
4. 写明 B 执行时必核的产物:`docs/survey/export_contract_gap.md`(否则缺陷会漏)。
5. 新增 §5 四条可自行跑的复核命令,§6 记「裁定状态:未裁定」。

**凭什么**:判据(1)原文无豁免条款;方向 1 的两条针对性检索式命中为 0(`docs/survey/direction_1_4_literature_retry.md` §1)。

**怎么自己复核**:在 `docs/survey/direction_1_ruling_request.md` §5 的四条命令逐条跑;第 4 条应命中正文里「判据(1)就不能判为整体通过」。

---

## 6. 本轮未做的事(留痕,避免下次重复踩)

| 项 | 为什么没做 |
|---|---|
| 未把方向 1 移出清单 | 那属选项 B,改变判据口径,须用户书面确认 |
| 未改判据文本 | 改判据文本须带独立裁决(auditKey),我不能自己顺手改「怎样算完成」 |
| 未声称任何判据通过 | 判据(1)仍有方向未满足,交付时如实记账 |
| 未写 `clear/` 下任何文件 | `clear/evidence` 等目录属系统所有,做的人不能写 |
| 未 `git add .` | 工作树有 30+ 个与本轮无关的改动(见 `docs/project_memory.md`),按 AGENTS.md 不得整树暂存 |

---

## 7. 仍然开着的未解释项与本轮的关系

| 项 | 本轮是否处置 | 处置处 |
|---|---|---|
| u1 台账「文档写法」列仍写已过期年份 | **是** | v2(§2 节)+ 命令 #2 |
| u2 `distinct_pmids=34` 易被误读为被引数 | **是** | v3(§3 节)+ 命令 #1 |
| u3 `docs/survey/lit_resnapshot_gap_20261016.txt` 的 PMID 列全是字面 `None` | 否 —— 该 `.txt` 是导出脚本取错字段名留下的中间产物,正确值在 `docs/survey/lit_resnapshot_gap_20261016.xml` 与 `snapshot_accounting.json` | 未改,属遗留 |
| u4 `docs/survey/self_check.json:45` 判据(4)记 pass 而 :58 自认按声称判 | 否 —— 该文件是第五轮产物,其 pass 属那一轮的自查口径,本轮不改写历史文件 | 未改,属遗留 |
| u5 `docs/survey/direction_1_4_literature_retry.md:35` 与 `:41` 对 `_tmp_retry_esummary.xml` 是否含 42800986 自相矛盾 | 否 | 未改,属遗留 |
| u6 reconciliation counts 仍记 23/26 | **是** | v3(§3 节)+ 命令 #3 |

**说明**:u3 / u4 / u5 不在本轮计划范围内,且都指向第五轮产物的措辞,不是本轮新增的缺口。它们仍在开着的未解释项清单里,结案时随交付一并交给评估者,不由我在产物里自行销案。
