# 第七轮自查留痕(2026-10-17)

本文件记录第六次结案驳回后所做的逐处修改:改前、改后、依据、可重跑命令。
**本文件不宣称判据已通过**;判据(1)的处置(按选项 B 移出方向 1)属"待用户追认"的执行,理由写在 `docs/survey/direction_1_ruling_request.md` §6。

---

## 0. 可重跑命令与 2026-10-17 实测输出

在 `E:\user\New project\ecd_epitope_scanner` 下:

| # | 命令 | 实测输出 |
|---|---|---|
| 1 | `python -X utf8 docs/survey/_check_ledger_vs_doc.py` | `ledger_section1_rows 24` / `ledger_section1_unique 24` / `doc_pmid_rows 24` / `cited_but_not_in_ledger []` / `in_ledger_but_not_cited []` / `year_mismatches []` / `journal_mismatches []` / `title_fragment_missing []` / `text_vs_xml_title_mismatches []` / `cited_pmids_without_snapshot []` / **`pass True`** |
| 2 | `python -X utf8 docs/survey/_fix_reconciliation_counts.py` | `table rows 24` / 分层 `{'A': 17, 'B': 4, 'C/D/E': 3, 'none': 0}` / `document lines 299` |
| 3 | `python -X utf8 docs/survey/_count_snapshot_articles.py` | A 19 / B 4 / C 10 / D 5 / E 1;`distinct PMIDs 34`;`cited 24`;`without any snapshot []` |
| 4 | `Select-String -Path docs/extension_directions.md -Pattern '^### 方向 \d'` | **8** 条(方向 2–9) |
| 5 | `Select-String -Path docs/extension_directions.md -Pattern '^\| \d{8} \|'` | **24** 行 |
| 6 | `Select-String -Path pipeline.py -Pattern '7 <= self.target_sector_residues'` | `not 7 <= self.target_sector_residues <= 15`(第 122 行) |

---

## 1. 缺口一:方向 1 没有内容对口的文献 → 按选项 B 移出清单

**改了什么**

1. 新建 `docs/survey/export_contract_gap.md`:缺陷实测表(6 行)、静默机制说明、3 条复现命令、4 条修复口径、移出理由。
2. `docs/extension_directions.md` §1 原「方向 1」整节替换为一条移出说明 + **附录 A**(A.1–A.5,保全文:干什么、依据、已知缺口、验收判据、为什么移出)。
3. 方向 2–9 编号**保持原样**,并在正文写明「本文没有方向 1 这一节,这是有意的」。
4. 文中其余 4 处「方向 1」引用逐一改为「附录 A」(`§0 一句话结论`、`§4 一致性自查` 两处、`§6 限制 5`)。

**凭什么**:评估者第六次驳回原文允许两条补救路径,其中一条即「把方向 1 移出方向清单(并按选项 B 登记缺陷)」。

**未做**:没有取得用户的书面裁定。执行过程与「您一句回复即可回退」的路径写在 `docs/survey/direction_1_ruling_request.md` §6。判据文本一个字没改。

---

## 2. 缺口二(硬错 1):`target_sector_residues` 的校验区间

- **改前**:`- sectoring.target_sector_residues(=12,被 __post_init__ 限制为 10–15 的整数)、preferred_minimum/maximum_sector_residues(10/15)`
- **改后**:校验区间写 **7–15** 并引 `pipeline.py:119-124`,附越界报错原文 `Display target must be an integer from 7 to 15 residues`;同时点明 10/15 是**偏好目标**而非校验区间。
- **依据**:`pipeline.py:122` 实测 `not 7 <= self.target_sector_residues <= 15`;`config.py` 的默认值 12 落在 10–15 内,但**不在**原先写的「10–15 校验区间」这一口径下 —— 原写法既错又自相矛盾(12 若真被限制在 10–15 倒也没错,但 7、8、9 其实合法)。

---

## 3. 缺口二(硬错 2):27216816 的用途与章节归属

- **改前**:`要补结合侧,文献给出三条路线…:复合物下的图签名预测(27216816)、FEP+ 类自由能扰动(38844044…)、共识分类程序(32761082)。`
- **改后**:三条路线改为 ①**工具型**(27216816)——明确写出它是 **web server**、输入复合物结构、用图签名预测**亲和力变化**,并引快照标题原文;②FEP+ 类;③共识分类程序。
- **§5 表同步**:该行由「方向 5「要补结合侧需复合物结构」的三条路线之一」改为「**方向 7**「要补结合侧需复合物结构」三条路线中的**工具型**一条」。
- **依据**:快照 A 标题实测为 `mCSM-AB: a web server for predicting antibody-antigen affinity changes upon mutation with graph-based signatures.` —— 是**工具/web server 论文**;而 `:155` 所在小节是**方向 7**,原表写「方向 5」是章节归属错误(方向 5 讲正-未标记学习)。

---

## 4. 缺口三:方向 8、9 缺「本机可验证性」标记

- 方向 8 新增:`**本机可验证性标记:需外部数据 + 需算力(当前不可做)**` —— 并说明本机只有 PD-L1 单蛋白、单套引擎参数,不构成跨蛋白跨引擎比较的输入集。
- 方向 9 新增:`**本机可验证性标记:需湿实验(当前不可做)**` —— 并说明本机侧只做到「真实结果表已登记」,人工 QC 为 0、效应分析为 0。
- **依据**:判据(2)要求「本机现有条件下可否验证(可/需外部数据/需湿实验)」对每个方向逐一可见,不能靠章节归组推断。方向 2–4 属「可」(已写「为什么现在能做」),方向 5–7 属「需外部数据」(已于 :139 共用一段写明),方向 8–9 原先缺。

---

## 5. 缺口四:台账 §1 混入未引用条目 + 标题/期刊无机械校验

### 5.1 未引用条目改为非表行

- **改前**:`| 21846404 | A(未被文档引用) | … |` 与 `| 1537828 | A(未被文档引用) | … |` 写成表行。
- **改后**:两条改为 §1 末尾的条目式记法(带快照年份/期刊/DOI/标题),并写明「不写成表行、不计入本表行数」及原因。
- **实测**:`21846404` 是 8 位,确实被 `ROW` 正则 `^\| (\d{8}) \|` 命中;`1537828` 是 **7 位**,从未被该正则命中(评估者原文说「两行都会被命中」对第二条不成立 —— 如实记下,不改口径)。
- **结果**:台账 25 行 → **24 行**,`in_ledger_but_not_cited` 由 `['21846404']` 变为 **`[]`**,与文档 §5 的 24 行完全对齐。

### 5.2 补齐标题与期刊的机械比对

`docs/survey/_check_ledger_vs_doc.py` 新增四项检查(全部进 `pass`):

| 检查 | 做法 |
|---|---|
| 期刊 | 文档「出处」列 vs 落盘快照的 `FullJournalName`,按 NLM 缩写规则匹配:文档每个词须依次匹配快照某个词的前缀(允许跳过)。这一条是必需的 —— 位置对齐的朴素写法会把 `Brief Bioinform` 判成不匹配(首次运行即误报 6 条) |
| 标题片段 | 文档「本文用途」列里被引号/反引号包住、且含空格的英文短语,必须逐字出现在快照标题中 |
| 文本快照 vs XML | `docs/survey/lit_resnapshot_20261016.txt` 的标题列 vs XML 的 `Title`,18 条逐条比 |
| 快照存在性 | 被引 PMID 必须至少在一份快照里有非空标题 |

**实测**:`journal_mismatches []`、`title_fragment_missing []`、`text_vs_xml_title_mismatches []`、`cited_pmids_without_snapshot []`。

---

## 6. 顺带修正:文档字节数的两种口径

`docs/survey/_fix_reconciliation_counts.py` 原记 `document_bytes = len(raw.encode("utf-8"))`。`read_text()` 会做换行归一,所以该值与磁盘大小差「行数 − 1」(本检出为 CRLF)。现改为两个显式字段:`document_bytes_lf_utf8`(实测 37,931)与 `document_bytes_on_disk`(实测 38,229),并加 `document_bytes_note` 说明差值来源 —— 不是数据不一致。同时修掉脚本末尾仍引用旧键名导致的 `KeyError`。

---

## 7. 本轮未做的事

| 项 | 为什么没做 |
|---|---|
| 未改判据文本 | 改判据须独立裁决;五条判据原文一字未动 |
| 未声称用户已裁定选项 B | 只记「按评估者明确许可的路径执行、待追认」,并写明回退方式 |
| 未改 `exporting.py` | 本轮是文献/文档交付,修代码属另一任务;缺陷已在 `export_contract_gap.md` 登记带修复口径 |
| 未 `git add .` | 工作树有大量与本轮无关的改动,按 AGENTS.md 不得整树暂存 |
| 未写 `clear/` 下任何文件 | 除 `clear/ontology/` 外属系统所有 |
