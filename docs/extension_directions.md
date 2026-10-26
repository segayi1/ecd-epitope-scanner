# ECD Epitope Workbench:扩展方向与后续研究思路

**文档性质**:文献调研 + 扩展方向建议。本文只提方向、依据与可检验的判据,**不代表任何方向已被实施或已获新增执行授权**。文中所有数值均标注来源文件,未实测的一律写明"未实测"。

**调研时点**:2026-10-13,并于 2026-10-16 补做一轮二次拉取核对(检索通道:NCBI E-utilities;`web_search` 工具不可用,缺口见 `docs/survey/web_access_channels.md`)

**检索底账**(本文引用编号 [L]-[M] 均指这三份):
- `docs/survey/method_landscape.json` — 方法与评估格局,32 条,4 条检索式
- `docs/survey/benchmark_landscape.json` — 公开基准集、指标定义、独立评测结论,5 个数据集 + 4 条争议
- `docs/survey/mutation_effect_landscape.json` — 突变效应与稳定性预测,7 条精度证据

---

## 0. 一句话结论

这个工作台当前最强的能力是**几何与覆盖的合规审计**,最缺的是**外部定位**。方向可以分成三档:能用本机已有冻结数据立刻做的(方向 2–4;原方向 1 已移入附录 A)、需要外部数据或算力才能做的(方向 5–7)、属于长期科学验证、短期不可能结清的(方向 8–9)。

---

## 1. 现在就能做的干实验验证(不需要新湿实验数据)

这**三条**改动量都在本机可承受范围内,判据都能写成布尔值或数值阈值。它们**检验的是工程与几何不变量,不是生物学有效性**——这一点必须在结论里保留。

### (原方向 1 已移出扩展方向清单,改为缺陷修复任务)

原「方向 1:用现有未跟踪的验收收据做一份"覆盖可审计性"台账 —— 但先修一个已知缺口」已按 `docs/survey/direction_1_ruling_request.md` 的**选项 B** 移出本清单,降级为**缺陷修复任务**,完整内容(实测表、复现命令、修复口径、移出理由)见 **`docs/survey/export_contract_gap.md`** 与本文 **附录 A**。

**为什么移出**:该条是仓库自证型议题(导出字段与源数据是否对齐),两轮换词根 PubMed 检索无内容对口文献,而判据(1)要求每个扩展方向都有可核验的对口文献。移出后 §1 剩下 **方向 2–4** 三个当前可执行方向,每条都有内容对口的可核验文献。

**编号不重排**:方向 2–9 的编号**故意保持原样**,以免打断 §5 表、§6 限制表与 `docs/survey/` 下其他报告里对「方向 5」「方向 7」等编号的引用。本文因此没有「方向 1」这一节,这是有意的,不是漏写。

---

### 方向 2:科学门与软目标的单因素敏感性消融(覆盖几何链)

**干什么**:把配置里那些"决定结果形状"的数逐个单独扰动,看结论是否跨过阈值或换序。要扰动的是(全部取自 `config.py` 的 `JOINT_V4_ML_PARAMETERS` 与 `pipeline.py:109-124` 的 `__post_init__` 校验):
- `sectoring.target_sector_residues`(=12,被 `__post_init__` 限制为 **7–15** 的整数,见 `pipeline.py:119-124`,越界报错原文 `Display target must be an integer from 7 to 15 residues`)、`preferred_minimum/maximum_sector_residues`(10/15,这是**偏好目标**不是校验区间,12 落在 10–15 内)
- `maximum_geodesic_diameter_a`(30.0)、`maximum_euclidean_diameter_a`(22.0)、`minimum_sector_area_a2`(50.0)、`minimum_sector_residues`(4)
- Surface V2 侧:`minimum_global_surface_coverage`(0.95)、`min_surface_area_A2`(8.0)、`min_surface_fraction`(0.20)、去重加权 Jaccard(0.75) [`surface_patch/engine.py:104-109,162-173`]
- 可靠性阈值:`structured_plddt=70.0`、`low_plddt=50.0`、`core_pae_threshold=15.0` [`structure_reliability.py:18-43`]

**依据**:
- 这个工作台的**敏感度证据目前只有零散痕迹**:PD-L1 那一例记录了最小配对签名距离 `8.2e-05`,而同一次运行的平均配对分离度是 `0.072872`(`panel_identifiability_v41.json`)——两者差三个数量级;`partition_diagnostics_v4.json` 自己写着 `size_exceptions` 四条(9/9/9/8 区)原因为 `LOCAL_BALANCE_LIMIT`,并注明"不代表全局无解"。
- 关于"参数扰动会不会改变结论"的一般性理由,文献里对同类计算给出的答案是会:`docs/survey/mutation_effect_landscape.json` 记录 SKEMPI 2.0 基准上的性能"seem overly optimistic",八种预测器泛化能力有限 [M, PMID 38197311];`benchmark_landscape.json` 记录七个结构预测模型上基于置信度的排序无法挑出最佳预测 [PMID 42252708]。这些是**别的系统**的证据,不能直接搬到本仓库,但足以说明"参数与排序的稳定性"是一个值得单独测量的量。
- 仓库侧的证据是:软目标与硬约束混在同一次运行里报告,`count_is_soft_target=true`、`size_is_soft_target=true` 只以布尔字段存在,`search_scope` 自陈是局部搜索。

**可重跑命令**(现有脚本,几何重放):
```powershell
cd "E:\user\New project"; .\ecd_epitope_scanner\.venv\Scripts\python.exe ecd_epitope_scanner/scripts/replay_joint_v41_run.py workspace_data/fa85119485cd <new_dir> --geometry-only
```
每次扰动必须用**新目录**,不得覆盖历史运行目录。

**通过判据**:对每个被扰动的参数,给出"结论不变 / 结论改变"的判定,并写明改变发生在哪个量上(分区数、分区大小、直径、面板成员)。若某参数在 ±10% 内即改变面板成员,记为该链条的敏感点。

**不通过/反面结果的价值**:如果全部参数在合理范围内都不改变结论,说明结论由数据主导而非参数主导——这本身是对该链条稳健性的证据,应当同样记录(它不构成对生物学有效性的证据)。

---

### 方向 3:把未解假设与覆盖盲区做成显式的"下一步实验清单",并验证它与面板的对应关系

**干什么**:用已算好的假设与盲区数据,核对"面板覆盖了哪些假设、漏了哪些",再把漏掉的转成下一阶段候选,而不是悄悄丢掉。

**依据(已实测)**:
- `outputs/484899e022f3/panel_identifiability_v41.json`:179 个几何/化学扰动签名,**1 个未解假设组包含全部 179 个**,`uncovered_hypotheses` **54 个**;`interpretation` 字段自陈"separation is not experimental binding evidence"。
- `epitope_hypotheses_v41.json`:179 个假设,`generation_mode` 全部为 `fixed_radius`,`quality_flags` 中 `WARNING_FEW_RESIDUES` 出现 **164 次**(即 164/179 个假设都带这条警告)。
- `analysis.json` 的 `epitope_mapping_plan.blind_spots` 只有一条:`{"sector_id":"VS007","reason":"NO_SELECTED_FEASIBLE_CANDIDATE"}` —— 有一个分区没有可用候选。
- `probe_hypothesis_matrix_v41.csv`:476 候选 × 179 假设 = 85,204 行,其中非零 **30,874**,零 **54,330**(稀疏度 64%)。
- `epitope_mapping_plan.stage2_presets`:28 条去卷积预设,全部 `conditional_on_parent_probe_hit=true`。
- **文献侧的直接对口依据**:PU-GRAIL 在残基层面用正-未标记学习处理"未标记 ≠ 阴性" [PMID 42412842, Bioinformatics 2026, `10.1093/bioinformatics/btag263`]。**旁证及其边界**:另一条 PU 学习研究报告的是"Veterans 电子病历中未编码自伤的检出"(快照标题原文:*Detecting Uncoded Self-Harm in Veterans' Electronic Health Records Using Positive and Unlabeled Learning: Retrospective Cohort Study*, PMID 42241701)——它只能说明该学习方法被用在另一个领域的未标记问题上,**生物学语境不可迁移,不构成本方向的证据**。本文件**不声称**该思路在别的领域"已被当作标准做法":那句表述是摘要层的归纳,快照只给出标题,按 §6 限制 7 的口径属"底账转述,未回原文核",已删除。

**可重跑命令**:
```powershell
cd "E:\user\New project"; .\ecd_epitope_scanner\.venv\Scripts\python.exe -c "import json;d=json.load(open('ecd_epitope_scanner/workspace_data/full_runtime_restore_20261012/outputs/484899e022f3/panel_identifiability_v41.json',encoding='utf-8'));print(len(d['uncovered_hypotheses']),len(d['hypothesis_signatures']))"
```

**通过判据**:54 个未覆盖假设与 VS007 盲区,每一个都要在产物里有明确去处(被某候选覆盖、被判定为不可设计、或被显式登记为待外部证据)。不允许"既没覆盖也没登记"。

**注意保真**:`WARNING_FEW_RESIDUES` 出现在 164/179 个假设上,说明这些假设的残基很少——在把它当成"待验证表位"之前必须先看这条警告,否则会把采样产物当生物学实体。

---

### 方向 4:同一冻结输入在两次独立执行下的逐记录一致性(重跑一致性)

**干什么**:用已有的重放能力,对同一冻结输入跑两次,比较产物中应当确定性的部分是否逐条一致。

**依据**:
- 本机已有重放脚本族,例如 `scripts/replay_joint_v41_run.py`、`scripts/replay_v4_selection.py`,以及 `docs/project_memory.md` 记录的用法。
- 文献侧对"可重复性"的一般性支持见 `docs/survey/method_landscape.json` 的检索式 B/C 记录;注意本文不把任何文献结论当作本仓库已通过该项测试的证据。
- **内容直接对口的文献(第二轮补检所得,2026-10-16)**:`27401684`(GigaScience 2016, *Tools and techniques for computational reproducibility*)——摘要原文把"同一分析任务、同一数据"应产生同一输出当作可检验的默认预期,并指出实践中失败源于软件打包/安装/执行与记录不足;这正是本方向要检验的集合相等命题及其替代说明。辅助:`35277501`(Scientific Data 2022, FAIRly big,"re-executed independent of the original computing infrastructure")、`42800986`(Array 2026, FairFlow,"identical code can produce different results depending on which versions are installed")。三者的摘要原文与逐条对口判定见 `docs/survey/direction_1_4_literature_retry.md` §3,原始返回已落盘为 `docs/survey/_tmp_retry_efetch.xml` 与 `_tmp_retry_esummary2.xml`。
- **仅作旁注、不承担支撑**:CWLProv 的工作流溯源互操作 [PMID 31675414, GigaScience 2019, `10.1093/gigascience/giz095`] 与"同一数据、不同分析者"的差异研究 [PMID 39915771, BMC Biology 2025, `10.1186/s12915-024-02101-x`]——前者是溯源口径而非重复性实测,后者对象是**不同分析者**而非同一代码的两次执行。

  **附带一条自我约束(2026-10-16 起):不拿"机制相近"充作"问题相同"。** 判断标准是:该文献陈述的被检验命题是否与本方向的对象同构。按此标准,上述三篇属于同构(直接对口),31675414 与 39915771 不属于。

**通过判据**:两次运行在**分区成员集合、面板候选 ID 集合、导出清单文件集合**三处逐条一致(集合相等,不是数值近似)。只要有一处不同,就按"未通过"记录并定位差异来源,不改判据。

**边界**:这项检验只证明计算可重复,不证明计算正确,更不证明推荐在实验室里成立。

---

## 2. 需要外部数据或额外算力才能做的

### 方向 5:把负样本问题正面处理掉(正-未标记学习)

**现状边界**:`docs/project_memory.md` 明写"实验标签稀疏且未知表位不是负样本";`benchmark_landscape.json` 就 `AgAbDb` 记着"以已知复合物结构为准,天然偏向已被解析的抗体-抗原对,对未知表位不提供负样本" [PMID 25048123]。**引用强度说明**:`25048123` 是一篇 *Methods Mol Biol* 的**数据库介绍章节**,不是独立评测研究 —— 它陈述的是该数据库自身的收录口径,可作为"某数据库如何界定正样本"的**范围说明**,但不能当作"负样本缺失问题的实测证据";方向 3 的立论主要由本仓库自身的标签现状(`project_memory.md`)承担。

**文献给的路线**:PU-GRAIL 用正-未标记学习在残基层面处理同构问题 [PMID 42412842, `method_landscape.json`]。这与本仓库"未知 ≠ 阴性"的立场是同一问题的两种表述。

**为什么现在做不到**:需要一套带部分标签的外部数据集,本机没有;而且要把残基级 PU 方法与扇区级覆盖目标对接,中间还缺一个定义映射。

**若要启动需要什么**:外部标注数据 + 一套明确的"未标记"采样的定义。

---

### 方向 6:给"几何可达 ≠ 细胞可达"这条边界接上复合物层面的证据

**现状边界**:`分析产物 warnings` 里逐字写着"Cell-surface assay accessibility is unassessed: membrane-relative geometry is a structural proxy, while intermembrane distance and receptor/Fab-specific reach remain unknown";`deliverables/project_introduction_20261013/项目介绍.md` 第 5 节也把 ProCaliper 限定为"不构成独立物理方法验证;当前 worker 没有计算有限尺寸抗体的空间可达性"。

**文献给的路线**:复合物层面的几何精度本身已成瓶颈 —— AlphaFold3 在 50 个非冗余人源化抗体-抗原复合物上只有 46% 能给出中等或高质量结合界面,且氢键与精细接触强度的复现一致性更差 [PMID 42719814, `benchmark_landscape.json`];SNAC-DB(收录抗体**与纳米抗体 VHH**-抗原复合物)对七个 SOTA 模型的评测显示成功率很少超过 25%,且基于置信度的排序无法挑出最佳预测 [PMID 42252708]。另有 AF3 在抗体复合物上的预印本结论,其**快照标题原文**为 "A Structural Antibody Benchmark of AlphaFold3 reveals **Hallucinated Epitopes** and **a Bias for Orderness**" [PMID 42620036]。**说明**:本文只读到该条的 esummary 标题层,未取摘要或全文,因此这里引用的是**标题自身的措辞**,不引任何数值,也不把它当作已验证的实测结论(与 §5 表该行、§6 限制 7 的记载一致)。它与本仓库"AF3 必须带置信度文件、不得当实验结构"的要求方向一致。

  > **该条为预印本**:`42620036` 的出版方确认为 `bioRxiv : the preprint server for biology`(DOI 前缀 `10.64898` 为 bioRxiv 专用),**未经同行评审**;它经 PubMed 检索式 A「epitope prediction benchmark antibody」入库,底账 `method_landscape.json` 的 `journal` 字段本就写着 "bioRxiv (preprint)"。它对本条归因的**削弱程度为低**:方向 6 归入"需外部数据"靠的是本仓库自身的输入契约(无抗体/复合物字段)与 `analysis.json` warnings 原文,该预印本只是"结论方向一致"的补充,**其数值在本文中未被引作数值证据**。详见 `docs/survey/preprint_status_note.md`。

**为什么现在做不到**:本仓库的输入契约不含抗体序列或复合物结构,这一条是设计边界而非实现缺口。

**若要做**:需要引入复合物结构或对接姿态,并且必须先定义"有限尺寸探针可达性"的算法口径——existing ProCaliper worker 明确不做这件事。

---

### 方向 7:把 FoldX 的折叠稳定性量与结合量分开,并补上结合侧证据

**现状边界**:`docs/project_memory.md` 与 `tasks/lessons.md` 都明写"FoldX DDG 不是 binding ΔΔG"。

**文献给的依据与路线** [`mutation_effect_landscape.json`]:
- 两个量在定义式、测量方式与所回答的问题上不同:折叠稳定性 ΔΔG = ΔG_mut − ΔG_WT,测自化学/热变性曲线;结合 ΔΔG_bind = ΔG_bind(mut) − ΔG_bind(WT),测自 SPR/BLI/ITC。**互不可换算**。
- 要补结合侧,文献给出三条路线,全部需要抗体-抗原复合物结构:①**工具型**(27216816)—— mCSM-AB 是一个 **web server**,输入抗体-抗原复合物结构、用图签名预测突变引起的**亲和力变化**(快照标题即 `mCSM-AB: a web server for predicting antibody-antigen affinity changes upon mutation with graph-based signatures.`),属「拿来就能算」的实现,不是方法学新路线;②FEP+ 类自由能扰动(38844044,成本高且需离群校正);③共识分类程序(32761082)。三者的定位差异见 §5 表该三行。
- 即便走 FEP+,CAGI5 盲测显示没有一种方法在所有指标上严格占优 [PMID 31243847],而 CAGI6 中 AUC 在部分挑战上达 1 [PMID 38227011] —— 小样本挑战的高指标不能当普适精度。

**结论**:这一条属"需外部数据",不是"现在可做"。在拿到复合物结构前,产物的命名纪律(禁止裸 ΔΔG、FoldX 输出只能写 ΔΔG_fold)必须保持。

**口径提醒(本章三条共用)**:方向 5–7 之所以归入"需外部数据",不是因为实现难,而是因为本仓库的输入契约(UniProt/PDB/AF3 结构 + ECD 区间)**不含抗体序列或复合物结构**。这一限制写在 `pipeline.py:66` 起的 `PipelineRequest` 字段表里(无任何抗体/复合物字段),也写在 `analysis.json` 的 `warnings` 中("membrane-relative geometry is a structural proxy, while intermembrane distance and receptor/Fab-specific reach remain unknown")。任何绕过该契约的做法都必须先改输入契约,属于设计决策,不是脚本改动。

---

## 3. 属于长期科学验证的(短期不可能结清)

### 方向 8:固定预算下的方法比较(跨蛋白、跨引擎)

**现状**:`deliverables/project_introduction_20261013/项目介绍.md` 第 12 节把它列为方向第 2 条;`docs/project_memory.md` 记录跨蛋白验证尚未完成,MBP/T4_lysozyme 冻结未完成。`benchmark_landscape.json` 的 `contrast_with_this_repo` 明确说明:本仓库评的是扇区级的几何与覆盖合规,公开基准评的是残基级排序质量或复合物结构精度,**口径不同,不可换算**——把覆盖率高写成 AUC 高属于跨口径换算。

**为什么慢**:需要冻结输入、构建预算、模型资源与判定合同,并且必须区分"有结构/候选准备/真实评分/物理解析/实验确认"五个层次,避免把准备完成写成验证完成。

**本机可验证性标记:需外部数据 + 需算力(当前不可做)**。本机虽有 PD-L1 冻结运行,但只有一个蛋白、一套引擎参数,不构成"跨蛋白、跨引擎"比较所需的输入集;冻结新的蛋白输入与重跑评分都需要模型资源。**不得**把它写成"现在可立即执行"。

### 方向 9:前瞻性的实验校准

**现状**:`项目介绍.md` 第 9.3 节记录,CTLA4 activation 项目 `e0f2c26650f04681` 登记状态 COMPLETE、revision 16、162 个样本、3 份 FlowJo 百分比导入,而**人工 QC 为 0,常规及百分比专用效应分析均为 0**;第 9.4 节的分层评价把"实验闭环"标为"真实结果表已登记,QC 与效应路径可用 / 尚需人工 QC、身份确认、真实效应与复测"。

**为什么慢**:这是唯一能把"几何代理"变成"真实结合后果"的路径,但它需要湿实验,且登记状态 COMPLETE 只说明登记事务状态,不能称实验分析完成。

**本机可验证性标记:需湿实验(当前不可做)**。本机侧只做到"真实结果表已登记"(162 个样本、3 份 FlowJo 百分比导入),人工 QC 为 0、效应分析为 0;把几何代理校准成真实结合后果必须新做湿实验,本机任何脚本改动都不构成这条方向的验证。

---

## 4. 与仓库现有表述的一致性自查

写本文时逐条对照了 `deliverables/project_introduction_20261013/项目介绍.md`,未发现冲突:
- 该文第 5 节的分工表(原生 SASA / ProCaliper / DiscoTope / SPURS / FoldX / 表面图 / Paraplume)与本文各方向的前提一致。
- 该文第 9.1 节的 PD-L1 数字(215 规范化残基、476 组合、60×5、105 节点 10 区、10 个专属首选、8 个正式面板)与本文**附录 A**引用的验收收据一致。
- 该文第 11 节的边界(SASA ≠ 抗体可达性、邻域覆盖 ≠ 实验结合覆盖)与本文方向 6 的前提一致。
- 该文第 12 节的五条发展方向与本文方向 5–8 部分重叠,本文补充的是**文献依据与可执行判据**,不改变其"未实施"的性质。

**本文新增但该文未含的一项**:附录 A 记录的导出契约缺口(`panel_coverage_matrix.csv` 两列为空)。这是本次侦察的实测结果,已登记为缺陷修复任务并定位到 `exporting.py:332-336`(表头 `:335`)与静默取空的 `exporting.py:200-205`,完整登记见 `docs/survey/export_contract_gap.md`。(此前把根因写成 `exporting.py:221-226` 是错的,该处是 `assay_profile` 的写法。)

---

## 5. 文献清单(全部可在 `docs/survey/` 底账中查到)

| PMID | 年份 | 出处 | 本文用途 |
|---|---|---|---|
| 38197311 | 2023 | Brief Bioinform | 突变效应预测基准"过于乐观"、泛化有限 [M] |
| 31243847 | 2019 | Hum Mutat | CAGI5 盲测无方法全指标占优 [M] |
| 38227011 | 2025 | Hum Genet | CAGI6 相关性 0.74、小样本 AUC=1 的解读边界 [M] |
| 38844044 | 2024 | J Mol Biol | FEP+ 需离群校正,结合侧成本与限制 [M] |
| 32761082 | 2021 | Bioinformatics | 结合亲和力变化数据稀缺 [M] |
| 39455825 | 2024 | Nat Comput Sci | 几何学习 + 预训练策略改进稳定性变化预测;摘要层自陈标注数据不足 [M]† |
| 39872490 | 2025 | Innovation (Camb) | 自监督学习做超快折叠稳定性预测(标题只支持 "ultrafast" 这一层;无"零样本"、无"105×","105×"仅存于底账转述并已从本文删除)[M]† |
| 36611255 | 2023 | Brief Bioinform | 九种 SOTA 表位预测器性能低、部分不优于随机 [B] |
| 42252708 | 2026 | Protein Sci | SNAC-DB:抗体**与纳米抗体 VHH**-抗原复合物的 ML-ready 数据库;七个模型成功率 <25%、置信度排序失效 [B] |
| 42719814 | 2026 | Bioinform Adv | AF3 复合物界面精度瓶颈(46% 中等/高质量) [B] |
| 42620036 | 2026 | bioRxiv(**预印本,未经同行评审**) | AF3 抗体基准:快照标题原文含 "Hallucinated Epitopes" 与 "a Bias for Orderness"(本文未通读全文,不据此引任何数值)[L]† |
| 42412842 | 2026 | Bioinformatics | PU-GRAIL 正-未标记学习,负样本问题的正面处理 [L] |
| 25048123 | 2014 | Methods Mol Biol | AgAbDb 数据库条目综述(**类型为数据库介绍章节,非独立评测研究**)[B]† |
| 38177393 | 2022 | Nat Comput Sci | Absolut! 合成真值范式(方法学排序可迁移性;在线先发于 2022-12)[B] |
| 42824475 | 2026 | Front Immunol | 148 个 AI 表位模型、评测实践碎片化 [B] |
| 19055730 | 2009 | BMC Bioinformatics | ElliPro 基于结构的抗体表位预测**工具**(标题为工具论文;其数据集常被用作不连续表位基准)。**年份口径说明**:快照 `PubDate` 为 **2008 Dec 2**(即 `EPubDate` 2008 Dec 2,等同于在线先发日),本表取**期刊卷期所属年 2009**;口径不同于本表其余各行(见下方"年份口径"一段)[B]† |
| 31675414 | 2019 | GigaScience | 工作流溯源可跨执行互操作共享(CWLProv);方向 4 的**间接**文献 [A] |
| 38496213 | 2024 | GigaByte | Julearn:面向无泄漏评估与检查的 ML 库;仅作旁注 [A] |
| 39915771 | 2025 | BMC Biol | 同一数据、不同分析者的效应量差异;仅作旁注 [A] |
| 42241701 | 2026 | J Med Internet Res | 临床电子病历未编码自伤检出中的 PU 学习;**生物学语境不可迁移,不构成本文任何方向的证据** [A] |
| 27401684 | 2016 | GigaScience | 计算可重复性的工具与技术:**同一分析、同一数据应得同一输出**及其失败原因;方向 4 的**直接对口**文献 [A2](DOI `10.1186/s13742-016-0135-4`;efetch 的 `PubDate` 为 **2016 Jul 11**,卷 5(1):30,`ArticleDate` 同为 2016-07-11。`DateCompleted` 2017-10-30 是 MEDLINE 加工日,**不是出版年**) |
| 35277501 | 2022 | Scientific Data | FAIRly big:脱离原计算基础设施独立重跑并核验结果来源;方向 4 辅助 [A2] |
| 42800986 | 2026 | Array | FairFlow:代码相同、依赖版本不同即可产生不同结果;方向 4 辅助 [A2](**efetch 与 esummary 双通道核对一致**:2026-10-17 补拉 `esummary.fcgi?db=pubmed&retmode=xml&id=42800986`,原始返回存 `docs/survey/_fetch_esummary_42800986.xml`;标题/年 2026/DOI `10.1016/j.array.2026.101150`/卷 31 四项与 `_tmp_retry_efetch.xml` 逐字相同,比对见 `docs/survey/_esummary_42800986_check.json`,`verdict = consistent`) |
| 27216816 | 2016 | Nucleic Acids Res | mCSM-AB:用图签名在**抗体-抗原复合物**上预测突变引起的亲和力变化;**工具/web server 论文**,不是方法学新路线;方向 7「要补结合侧需复合物结构」三条路线中的**工具型**一条 [A](在线先发于 2016-05-23) |

来源标记:[L] = `docs/survey/method_landscape.json`;[B] = `docs/survey/benchmark_landscape.json`;[M] = `docs/survey/mutation_effect_landscape.json`;[A] = `docs/survey/direction_refs_addendum.md`(该文件记录了方向 1/3/4 的第一次补检过程;**其表内条目并非全部有落盘的 esummary 快照,该文件顶部已就此更正**);[A2] = `docs/survey/direction_1_4_literature_retry.md`(2026-10-16 第二轮补检;**三条的来源不同**:`27401684`、`35277501` 在 `_tmp_retry_esummary2.xml` 中;`42800986` **不在该文件内**,只由 `_tmp_retry_efetch.xml` 的完整 PubmedArticle 记录支撑)。
**† 标记的含义**:该行的"本文用途"含**摘要层转述或性质归类**,而标题层不足以坐实 —— 具体逐条见 §6 限制 7。打 † 的行不得被单独引用为已验证结论。

**年份口径(2026-10-16 二次核验后定稿)**:本表年份**一律取 PubMed 记录的 `PubDate`**(无 PubDate 时退回 `JournalIssue/Year`);`DateCompleted` 是 MEDLINE 加工日,**任何情况下都不作出版年**。唯一例外是 `19055730`,取期刊卷期年 2009 而快照 `PubDate` 为 2008-12-02,已在行内标明。逐条实测值见 `docs/survey/date_field_audit.json`(由 `docs/survey/_audit_dates.py` 直接从四份落盘 XML 解析生成,不经人工转录)。复核命令:

```powershell
python -X utf8 docs/survey/_audit_dates.py --out docs/survey/date_field_audit_recheck.json
Select-String -Path docs/extension_directions.md -Pattern '^\| \d{8} \|' | Measure-Object | Select-Object -ExpandProperty Count
```

> 第一条命令**必须带 `--out`**:不带它会把 `docs/survey/date_field_audit.json` 覆盖掉(脚本会在运行前打印 `NOTE: writing to the default path REPLACES the existing report.`)。第二条实测返回 `24`。

**核对方式(2026-10-17 按实际快照归属重写)**:本条**不再声称"上述条目均经 esummary 核对"**。分层**按每个 PMID 实际落在哪份快照**,不按快照总数硬套;**每行只归一层,合计等于本表行数 24**。逐条归属由 `docs/survey/_count_snapshot_articles.py` 用 XML 解析实测(复核命令见本节末):

- **第一层 · A 快照 esummary(17 条,`lit_resnapshot_20261016.xml`)**:`38197311`、`31243847`、`38227011`、`38844044`、`32761082`、`39455825`、`39872490`、`36611255`、`42252708`、`42719814`、`42620036`、`42412842`、`25048123`、`38177393`、`42824475`、`19055730`、`27216816`。标题/年份/期刊可逐字对照;DOI 需另经 efetch 或出版物页核对。
  - **本层已有更正**:`42241701` 曾被写成 A 层,实测它**不在 A 快照内**(只出现在 B),已按实测移到第二层;`39915771` 的期刊为 **BMC Biol**;`31675414`、`38496213` 两行的期刊与标题已据快照逐字复核。
  - **`27216816` 是本轮补入的**:它在方向 5 正文被引用却未进本表,属清单覆盖不全;快照标题为「mCSM-AB: a web server for predicting antibody-antigen affinity changes upon mutation with graph-based signatures.」,与正文用途一致。
- **第二层 · B 快照 esummary(4 条,`lit_resnapshot_gap_20261016.xml`)**:`31675414`、`38496213`、`39915771`、`42241701`。A 与 B **不相交**,故这 4 条不在第一层内,不能与第一层相加成"19+4=23 即 A 的 19 条"。
- **第三层 · C/D/E 快照(3 条,`_tmp_retry_esummary2.xml` + `_tmp_retry_efetch.xml` + `_fetch_esummary_42800986.xml`)**:`27401684`、`35277501`(C 与 D **两处都有**)、`42800986`(**不在 C**,由 D 的 efetch 记录与 E 的 esummary 单条补拉两处互校,标题/年/DOI/卷四项一致)。这一层**不在 A 也不在 B**,此前被误并入第一层 19 条。
- **第四层 · 仅经底账转述(本表 0 条)**:`docs/survey/direction_refs_addendum.md` 表内除上述条目外的检索命中均无落盘快照;本表**没有**引用这一层的条目。
- **数值层级**:本表所列"论文自述"均来自**摘要层或底账转述**,未通读全文;凡属摘要层/正文层的量(该 Pythia 加速比即为一例),本表一律只按"底账转述、未回原文核"记录,且已从**表里**撤下、不再作为任何方向论据(见 §6 限制 7)。**撤下后仍留着字样的地方只有 §6 的限制表与本行**,都属"记录下来并声明不可用"的元陈述,不构成对该量的引用。

**可机械复查的条目数 = 24**(即本表行数)= 第一层 17 + 第二层 4 + 第三层 3 + 第四层 0。不是"23+3=26",也不是"19+4=23"。

**归属复核**(两条命令,前者列出每个 PMID 实际在哪份快照、后者数表行):
```powershell
python -X utf8 docs/survey/_count_snapshot_articles.py
python -X utf8 -c "import json;d=json.load(open('docs/survey/snapshot_article_count.json',encoding='utf-8'));print(d['efetch_pubmedarticle_count']);print(d['pmid_snapshot_membership'])"
```

各底账的 `cross_check.mismatches_found` 均为 0,但该字段只覆盖各底账**自身**的核对范围,不能代替本表分层。唯一标注的年份口径差是 `19055730`(快照 `PubDate` 2008-12-02 / 本表取卷期年 2009),已在行内标明,不算矛盾。

---

## 6. 本次调研的限制

1. 检索仅经 PubMed E-utilities,未使用任何非 PubMed 的预印本服务器(bioRxiv/arXiv)直检通道,也未做非英文文献检索。PubMed 索引中的预印本记录**会被检索到**(本次即命中 1 条,见 `42620036`,已在底账中标注为 "bioRxiv (preprint)"、未经同行评审),但**未被 PubMed 索引的预印本不在本次覆盖范围内**。部分检索式的 `sort=date` 未被 PubMed 接受(返回 `Unknown sort schema 'date' ignored`),时间偏倚不可控。
2. 关键词召回含明显噪声(例如把电化学传感器论文拉进来),已人工筛除。
3. 所有文献结论均来自摘要层级,未通读全文;标注为"论文自述"的条目均以摘要可见部分为准。
4. 本文**未在本地复算任何公开基准指标**,也未下载任何基准数据集;所有基准数字都是文献报告值。
5. 方向 2–4 的判据是**本次新写的**,尚未执行过;执行它们需要另立任务并获得授权。
6. `web_search` 工具在本机不可用(缺 `DEEPSEEK_API_KEY`),检索只走 PubMed E-utilities;补上该凭据属用户决定,未代为填写(见 `docs/survey/web_access_channels.md`)。
7. **结论 7(2026-10-16 补充,内容级限制,不是检索级限制)**:上述限制 1–6 说的都是"检索通道覆盖到哪",但调查后另发现一类不同性质的问题 —— **本文件对已取得文献的内容转述与实际标题/年份不符**,共 6 处。这不是检索没覆盖到,而是转录时引入了原文没有的含义:

   | # | PMID | 本文件写法 | 快照实际 | 性质 |
   |---|---|---|---|---|
   | 1 | 39872490 | §5 表曾写「Pythia 零样本 ΔΔG 与 105× 加速」(2026-10-17 已**删除**;本文件正文与 §5 表现均无 "105×"、无 "零样本" 对该文的定性) | 标题为「Structure-based self-supervised learning enables ultrafast protein stability prediction upon mutation.」,无 "105×"、无 "zero-shot" | **数值不在标题中**,"105×" 只存在于底账 `mutation_effect_landscape.json` 的 `claims` 里(摘要层转述),本文未回原文核 → 已按"不得当作已验证的数值"删除 |
   | 2 | 42620036 | 「偏好有序区」(方向 6)与「幻觉表位、偏好有序区」(§5 表) | "a Bias for **Orderness**" | 同义改写;**该文为预印本** |
   | 3 | 39455825 | §5 表写「GeoDDG 的标注数据不足自述」 | "Improving the prediction of protein stability changes upon mutations by geometric learning and a **pre-training strategy**." | 标题未提数据不足;"标注数据不足"是摘要层的作者自陈,本文未回原文核 |
   | 4 | 42252708 | §5 表写「SNAC-DB」 | "SNAC-DB: An ML-ready database for antibody and **NANOBODY® VHH**-antigen complexes..." | 原标题明确含纳米抗体 VHH;本文件缩略写法丢了这层范围限定 |
   | 5 | 25048123 | 方向 8 写「AgAbDb 的条目也记着…」 | 快照标题确认为 AgAbDb 综述条目 | **类型不符**:该条是 *Methods Mol Biol* 的 **数据库介绍章节**,不是独立评测研究;本文件把它当作"某研究的实测记录"引用过重 |
   | 6 | 19055730 | §5 表写「ElliPro 不连续表位基准集」 | "ElliPro: a new structure-based tool for the prediction of antibody epitopes." | 标题是**工具论文**;"基准集"是工具附带数据集的性质,不是该文标题的主张 |

   **这 6 处对结论的影响**:第 1 处最相关 —— 曾用 Pythia 的加速比支撑"零样本 ΔΔG 可作快速筛选"的性价比论据,该量与"零样本"定性**都不在本文可核的快照里**;2026-10-17 已从 §5 表删除,方向 5(§2 方向 5,`docs/extension_directions.md:123`)的立论因此**不依赖任何未核数值**(它只依赖本仓库自身的标签现状与 PU-GRAIL 的方法同构性)。其余 5 处不改变方向归属,只改变措辞强度(第 5 处已把 AgAbDb 的引用强度下调,第 3、4 处补上原标题范围)。
   **处置**:上述 6 处已逐条改写在正文各方向处与 §5 表中;改写只做"退回到快照可核内容",不新增任何未经二次拉取的断言。
   **残留限制**:本次二次拉取只取到 esummary 层的**标题/年份/期刊/DOI**,未取全文;"105×""标注数据不足"这类**摘要层或正文层的量**,在本文中一律标注为"底账转述,未回原文核",不得当作本文已验证的数值。

---

## 附录 A:已移出扩展方向清单的条目 —— 原「方向 1」(缺陷修复任务)

**状态**:按 `docs/survey/direction_1_ruling_request.md` 的选项 B 处置 —— 从扩展方向清单移出,降级为**缺陷修复任务**。完整登记(含复现命令与修复口径)见 `docs/survey/export_contract_gap.md`。本附录保留原有的方向定义与验收判据,便于追溯。

### A.1 干什么

把已经存在但散落在 `workspace_data/` 各次运行里的验收 JSON 收成一份可复查的台账,核对每一例的**显示分区 → 面板候选 → 导出通道**三者是否互相闭合。

### A.2 依据(已实测的具体数字,全部来自 PD-L1 冻结运行 `484899e022f3`)

- `workspace_data/full_runtime_restore_20261012/full_run_medium/scientific_acceptance.json`:显示分区 10 个,大小 `[10, 9, 9, 10, 9, 8, 11, 15, 13, 11]`,共 105 个表面节点;`passed_exclusive_regions=10`,`formal_main_count=8`;worker 计数 `{procaliper:215, spurs_single:215, spurs_multi:476, discotope:215, foldx:60}`。
- `workspace_data/full_runtime_restore_20261012/outputs/484899e022f3/partition_diagnostics_v4.json`:最大测地直径 **29.727 Å**(上限 30.0),最大 Cα 直径 **20.999 Å**(上限 22.0)。
- `full_run_medium/export_acceptance.json`:65 个清单文件、1022 条显示引用、`region_membership_complete_unique=true`、`recommendation_json_matches_api_disk_zip=true`。
- `full_run_medium/foldx_pair_acceptance.json`:60 候选 × 5 次 = 300 个 fresh 配对,复用 0,`all_pair_hash_contracts_valid=true`。

### A.3 已知缺口(必须先修,详见 `docs/survey/export_contract_gap.md`)

`outputs/484899e022f3/panel_coverage_matrix.csv` 只有 3 列 `design_id,patch_id,perturbed_fraction`,32 行中 `design_id` 与 `perturbed_fraction` **两列全为空**,只有 `patch_id` 有值。根因在 `exporting.py:332-336`(表头 `:335` 写死成三个与源字段全不匹配的键名)与 `exporting.py:200-205`(`_write_csv` 用 `csv.DictWriter(..., extrasaction="ignore")` + `row.get(field, "")`,表头键不存在时静默写空串)。源数据 `analysis.json` 的 `epitope_mapping_plan.coverage_matrix` 有 32 条,真实字段是 `probe_id` / `patch_id` / `scenario` / `predicted_coverage` / `evidence_type`。`README.md:129` 与 `exporting.py:523` 都把它描述为「构建体 × 补丁扰动矩阵」,与实际内容(probe × patch × scenario 覆盖表)不符。

> **根因行号更正(2026-10-17)**:先前写作 `exporting.py:221-226`,实测该处是 `assay_profile` / `fixed_assay_profile.json` 的写法,与本缺口无关。正确引用为 `exporting.py:332-336` 与 `exporting.py:200-205`。

### A.4 验收判据(台账的通过条件)

台账里每一例都同时满足 ——①分区成员对合格表面残基图完整唯一;②每个分区的最大测地直径 ≤ 配置上限、最大 Cα 直径 ≤ 配置上限;③导出通道之间的候选身份与数值一致。任一条为否则记为不通过,不修改上限来凑通过。

**不通过意味着什么**:说明某次运行的「已通过」只覆盖了它自己声明的字段,不代表导出链整体闭合。

**为什么现在能做**:数据已在本机,不需要任何新计算。

### A.5 为什么移出(文献层面的理由)

本条属**仓库自证型**——核心依据是仓库自身产物(导出契约与产物数值)。第二轮为它专门换词根再检一轮(见 `docs/survey/direction_1_4_literature_retry.md` §1):最贴近的两个检索式 `provenance audit scientific software output verification` 与 `export contract validation data pipeline` **命中数均为 0**,20 条候选中无一条的被检验命题是「导出字段与源数据是否对齐」。

判据(1)原文是「列出的每个扩展方向都引用至少一篇真实可核验的文献」,**没有「仓库自证型可豁免」的例外**。把它留在扩展方向清单里,判据(1)就无法整体通过;移出后 §1 的每个方向都有内容对口的可核验文献,而这条缺陷并未丢失 —— 它完整登记在 `docs/survey/export_contract_gap.md`,带复现命令与修复口径。

**未排除的可能**:PubMed 未索引的工程类文献(会议论文、标准文档、软件工程期刊)可能有过对口研究;本次未检索非 PubMed 通道,故只能断言「在本次检索覆盖范围内没有」。
