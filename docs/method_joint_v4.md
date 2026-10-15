# Method

## 2026-10-09 正式面板双模型与证据展示

正式主构建体必须同时具有有限的SPURS组合和FoldX组合评分，并通过现行准入门槛；旧SPURS-only配置和缓存PASS不构成例外。规则检查与数值模型检查分开展示，前者提供实际规则事实，数值门槛不适用。直接上传PDB/CIF不自动证明实验来源：从原始头部读取已知实验方法，预测来源另行识别。NMR结构不适用pLDDT或晶体学分辨率，局部实验精度未评估；坐标与序列映射可用不等于已证明精度高。完整追溯标识与序列仍保存于折叠详情及导出。详见[本轮修复记录](recommendation_detail_polish_20261009.md)。

## 2026-10-09 专属分区搜索与推荐解释

常规V4按实际突变位置区分专属与跨区组合，hypothesis锚定只代表生成来源。局部短名单形成后，如某区独立专属双点不足4组、三点不足2组，或最佳专属组合9/12/15Å平均本区代理覆盖低于90%，补充本区窗口；90%只触发搜索，不是资格门槛。保持12位点/3替换、300 beam、24Å跨度及科学门。全局池保留hypothesis软容量及专属基础保障，额外扩张仅用于补足基础组数；跨窗口物理去重、同位点不同替换不重复计独立组数。

FoldX默认最多60个物理组合（2026-10-12由42提升），先逐区一个专属组，再补第二组；专属预审候选不足时保留缺口，跨区候选只作补充，实际模型成本不重复。审计区分生成短名单、入池、SPURS评分、可确认、安排/成功FoldX、最终资格。专属首选优先合格者，否则最近门槛安全候选；跨区不顶替专属首选。更低预测风险、更强表位先验和跨区覆盖三类建议独立资格，同一物理组合合并原则。区域指标独立计算，不按遍历顺序扣减。

正式与离线页面共用摘要/展开详情组件，中文说明代替内部状态码，三尺度覆盖表格与完整门槛、位点、互补性、证据来源和序列同一数据；候选哈希仅作追踪。几何代理、模型分歧、缺失证据与真实实验扰动分开解释。旧任务只兼容显示，不伪称重新实施专属保障；物理预览沿用原专用策略。

## 1. 目的与适用范围

基于实验室之前的流程，使用一系列目前已有的生物信息分析，筛除明显缺乏结构映射或预测稳定性证据的候选，给出一系列用于实验粗定位的表面突变构建体，在给定蛋白结构、ECD 序列和外部模型证据的条件下，生成可用于实验粗定位的表面突变探针面板。

V4 组合使用互补工具：**ProCaliper** 是 2025 年发表的开源蛋白结构分析与可视化库，在本项目中用于计算 SASA/rSASA；其 SASA 实现与本地 **Bio.PDB Shrake–Rupley** 使用同一算法，因此两者用于一致性核对，而非独立方法验证（[ProCaliper 论文](https://doi.org/10.1093/bioadv/vbaf275)、[代码库](https://github.com/PNNL-Predictive-Phenomics/ProCaliper)、[Bio.PDB 文档](https://biopython.org/docs/latest/api/Bio.PDB.SASA.html)）。**DiscoTope-3.0** 是较新的结构型 B 细胞表位倾向预测器，利用逆折叠表示改进残基评分；在 V4 中只作软排序先验，不代表特定抗体的结合概率（[原始论文](https://doi.org/10.3389/fimmu.2024.1322712)、[官方代码](https://github.com/Magnushhoie/DiscoTope-3.0)）。**SPURS** 是近期发表的稳定性深度学习模型，融合 ESM2 与 ProteinMPNN 表征以预测单点和组合突变的稳定性变化，代表较新的数据驱动方法，但应用经验积累仍少于传统工具（[原始论文](https://doi.org/10.1038/s41467-025-67609-4)、[官方代码](https://github.com/luo-group/SPURS)）。**FoldX** 则是发展多年的经验力场工具，用于突变稳定性评估；V4 将它作为 SPURS 之外的组合复核证据，两者一致可增强计算筛选的可信度，但仍不能替代实验验证（[FoldX 论文](https://doi.org/10.1093/nar/gki387)）。

V4 的输出不是精确抗体表位预测，也不是结合概率，只是在当前结构上下文和配置的模型证据下，哪些表面区域可以用少量、可审计的突变组合进行实验区分，并且哪些组合没有通过稳定性或映射安全门控。

支持两种运行模式：

- `preview`：用于调试、UI 预览和候选审查；允许缺失单点分数，但不产生正式推荐。
- `strict`：用于正式计算筛选；必需证据和候选级硬门控必须全部满足，否则返回 `BLOCKED_EXTERNAL_EVIDENCE`。

外部模型权重、二进制和许可证不随仓库分发。仓库提供本地隔离 worker 协议、结果校验、缓存、审计和安全门控。

## 2. 输入与数据边界

V4 接收：

1. ECD 蛋白序列和 ECD 起止坐标；
2. 选择的 PDB/mmCIF 目标链；
3. 原始结构类型，例如实验结构或预测结构；
4. 可选原始全长 CDS；
5. 表面几何和结构可靠性信息；
6. ProCaliper、SPURS、DiscoTope-3.0 和 FoldX 的本地 worker 配置；
7. V4 参数快照，包括阈值、sector 数、候选阶数、面板预算和超时。

结构解析层同时保留三类位置：UniProt 位置、author residue number/插入码和规范化预测链位置。任何位置映射歧义、WT 氨基酸不一致、链不存在或结构残基缺失，均作为错误处理，不按数组顺序静默补位。

所有科学阈值由 `config.JOINT_V4_ML_PARAMETERS` 管理，并写入 V4 manifest。V4 的阈值是工程筛选阈值，不是经过实验校准的概率或放行标准。

## 3. 处理流程

V4 的确定性处理流程如下：

```text
原始结构/序列
      |
      v
链选择与 WT 映射校验
      |
      v
规范化连续编号预测链 + residue_mapping.csv
      |
      +--> ProCaliper surface/SASA
      +--> SPURS single
      +--> DiscoTope-3.0
      |
      v
表面证据融合与表位岛识别
      |
      v
连通表面 Voronoi sector 分区
      |
      v
单点替换菜单与 2/3 点候选搜索
      |
      +--> SPURS multi
      +--> FoldX confirmation
      |
      v
候选级稳定性/映射/外推门控
      |
      v
确定性 panel selection
      |
      v
V4 原生输出 + legacy UI 兼容投影
```

## 4. 规范化预测上下文与残基映射

### 4.1 规范化目标链

V4 将目标链写为单链、连续编号的 `prediction_context.pdb`。该文件只作为外部预测工具的 prediction context；原始结构仍作为 biological context 保存，不被静默替换。

预测序列为输入记录中按规范化顺序排列的 WT 氨基酸序列。规范化位置从 1 开始，结构内部缺失残基不会被虚构插入。

### 4.2 映射表

`residue_mapping.csv` 是所有外部结果坐标回溯的权威表，至少包含：

- `normalized_position`
- `uniprot_pos`
- `chain_id`
- `auth_seq_id`
- `insertion_code`
- `label_seq_id`
- `one_letter`

规范化结构写出时逐个检查结构残基是否存在以及三字母残基名是否与输入 WT 一致。错误使用以下稳定错误码：

- `CHAIN_MAPPING_AMBIGUOUS`
- `RESIDUE_COUNT_MISMATCH`
- `WT_RESIDUE_MISMATCH`
- `UNSUPPORTED_STRUCTURE`

外部 worker 只能使用规范化位置。返回结果通过位置范围和 WT 校验后，才映射回 UniProt 位置进入 V4 证据对象。

## 5. 四条外部证据线

### 5.1 ProCaliper：表面可及性

ProCaliper 用于提供第二次表面/SASA 计算。它和本地 native SASA 均采用 Shrake–Rupley 采样算法，因此二者的一致性是实现与参数敏感性检查，不能当作独立物理方法的交叉验证。V4 同时保留本地 native SASA、ProCaliper SASA、相对 SASA 和差异状态：

- `AGREE`
- `SASA_DISAGREEMENT`
- `NOT_COMPARED`

ProCaliper 成功结果必须覆盖规范化目标链全部位置。缺少位置、位置超界、数值非有限或输出不完整时，该证据线失败。

候选位点若处于 `SASA_DISAGREEMENT`，其表面分类缺少一致证据，不能进入 strict 正式面板。

### 5.2 DiscoTope-3.0：表位倾向先验

DiscoTope 分数只作为软表位先验，不代表真实抗体接触。V4 分别保存原始分数、校准分数、模型二元标签、在当前表面图内计算的分位数、图邻域局部峰值和 provenance；这些字段不能互相替代。成功结果必须覆盖规范化目标链全部位置。

高分位置和局部峰值用于建立 epitope islands。DiscoTope 只提供候选优先级，不是突变位点或正式推荐的硬门槛。缺少 DiscoTope 不会伪造分数；候选会全部标记为 `STABILITY_FALLBACK`，只要稳定性必需证据和组合门控完整，strict 模式仍可形成正式推荐，同时写入缺失先验的审计 warning。

### 5.3 SPURS single：单点稳定性筛选

每个结构映射位点枚举 19 个非 WT 替换。首先应用硬排除：

- 二硫键或 Cys 风险；
- Pro/Gly 主链构象风险；
- 糖基化 sequon 或已知糖基化位点；
- 其他输入中明确标记的保护位点。

SPURS single 分数按配置阈值分带：

| 分带           | 条件                                     | strict 资格      |
| ------------ | -------------------------------------- | -------------- |
| `PREFERRED`  | ddG ≤ `preferred_max_ddg`              | 可以进入候选         |
| `ACCEPTABLE` | preferred < ddG ≤ `acceptable_max_ddg` | 可以进入候选         |
| `REVIEW`     | acceptable < ddG ≤ `review_max_ddg`    | 不进入 strict 主候选 |
| `REJECT`     | ddG > review                           | 淘汰             |
| `NOT_SCORED` | 没有分数                                   | 仅 preview 可保留  |

单点矩阵和所有通过/未通过信息都必须导出。不能因为没有分数而把替换当作低风险。

### 5.4 SPURS multi 与 FoldX：组合复核

候选组合发送给 SPURS multi 和 FoldX。两条工具线使用独立配置和独立缓存键；SPURS single 的命令不会隐式替代 SPURS multi。

组合记录计算：

<div align="center"><i>epistasis</i> = ΔΔG<sub>SPURS,multi</sub> − Σ ΔΔG<sub>SPURS,single</sub></div>

组合输出必须满足：突变格式可解析、位置在规范化序列范围内、WT 与输入序列一致、组合属于请求集合、每个请求组合都返回且 ddG 为有限数值。worker 超时、崩溃、空输出、stdout 污染和 schema 错误都作为失败处理。即使失败输出携带部分评分，也不得作为 strict 放行证据。

## 6. 表面图、epitope islands 与 sector

常规显示策略为 `balanced_connected_v4`，以每区12个合格表面残基为目标，10–15为大小偏好。最少4个残基、50 Å²、诱导图连通、完整唯一覆盖及配置直径仍是硬约束；大小偏好及请求分区数是软目标。物理覆盖预览保留 `connected_bounded_v3` 的26/18 Å尺度，24 Å候选组合跨度与科学准入门槛保持原规则。

常规初始化尝试总表面残基数/目标大小附近的中心数及原请求数，复用原有Voronoi和小区修复，再执行合法边界迁移、交换、合并与局部重建。按偏好范围偏离、区数、目标大小偏离依次比较，使用稳定残基编号去重和决策。所有初始化与平衡共享100,000状态预算；有硬合格方案但优化预算耗尽时保留方案并明确标注，尚无硬合格方案时停止后续设计。局部搜索结束不是全局最优证明。大小例外记录独立小连通分量、局部搜索未改善或预算耗尽；不能为了凑大小跨不连通表面归并。

硬直径使用完整合格残基图最短路，诱导分区图直径另行报告供用户了解区内路径；两者不混用。正式结果默认保留完整API；网页显式请求轻量展示投影，完整候选和大规模成对诊断通过原始文件引用获取。投影保留推荐、门槛、位点、覆盖、资源审计和来源身份，离线导出复用相同投影。

V4 从显式溶剂表面点图构建残基级测地图。高暴露残基是 sector 节点；未达到突变暴露阈值、但在点图上具有真实溶剂表面的边界残基可作为测地路径中继，避免 rSASA 阈值在连续表面上打孔。若稀疏点采样遗漏肽键处的 owner transition，则仅在两个残基连续、都具有表面点且 Cα 距离不超过 4.5 Å 时加入主链桥。中继路径最长 15 Å，不允许形成跨结构核心的远程捷径。

DiscoTope 高分位置形成 epitope islands，但 island 是分区先验而不是独立 sector 的硬边界。每个物理表面连通分量先取得一个高先验中心，剩余中心使用确定性的最远点采样分配；随后执行最多 8 轮 minimax medoid 更新和测地 Voronoi 分配。Voronoi 分配后先拆分不连通或超径区域，再归并不满足最小大小的区域；k-center 本身不能保证最小残基数。

1. 每个物理表面连通分量至少获得一个中心；
2. 请求区域数（默认 7）是软目标；物理连通性和大小约束优先，构建体预算独立；
3. 按残基数、面积、最小残基编号处理小区域。相邻候选按两区残基之间的最短全局表面测地距离排序，同距按最小残基编号排序；优先直接并入最近的可行区域；
4. 整区合并超径时，保留包含小区域的合格连通部分，将冲突边界残基通过确定性回溯分配给周边区域。每次分配检查连通性和直径，完整方案通过后一次提交；失败则尝试下一候选，成功后重新识别小区域；
5. 最小残基数 `minimum_sector_residues`（默认 4）、最小表面积（默认 50 Å²）、区域连通性、全局表面测地直径（常规默认 ≤30 Å）和 Cα 欧氏直径（常规默认 ≤22 Å）都是最终硬约束。所有高暴露表面节点必须完整、唯一归属，不能删除低先验位置改变覆盖分母。

边界回溯最多访问 100,000 个状态；预算耗尽报告搜索未解决，不能解释为不存在可行分区。无法满足最终约束时保存分区诊断并停止后续设计，不返回带 1–2 残基小区的正常结果。导出记录分区策略版本、实际限制和残基迁移；所有导出使用同一最终分区。以上是显示分区约束，不改变 Surface V2 科学假设、突变评分或科学证据门槛。

若所有完整归并候选均无法安置冲突节点，则按最近邻顺序尝试重建受体边界：小区与受体的一部分形成合格区域，受体剩余部分也必须连通且满足最小大小和直径。此回退保持区域数，完整归并则减少区域数；两者共享同一个搜索预算和最终验收门。

每个 `V4Sector` 保存：

- sector ID；
- UniProt 残基位置；
- 表面连通分量；
- 总表面积与可设计表面积；
- epitope score；
- 中心残基、表面测地直径和欧氏直径；
- island ID；
- `EPITOPE_PRIOR` 或 `COVERAGE_ONLY` 类型；
- `LOW_EPITOPE_PRIOR`、`LOW_SURFACE_AREA`、`DISCONNECTED_SECTOR`、直径超限等质量标记。

sector 是实验粗定位单元，不是抗体真实 footprint。

## 7. 候选搜索

每个 sector 从有效单点选项中生成 2、3 点候选。候选必须满足：

- 所有位点来自同一生成窗口；启用hypothesis时窗口可重叠并跨显示分区；
- 位点在有效表面图中连通；
- 组合测地跨度不超过配置上限；
- 不包含单点硬排除；
- 保存逐位点化学类别和单点风险；
- 生成直接位点、9 Å、12 Å 和 15 Å 测地覆盖情景。

候选位点池和shortlist保留DiscoTope支持/回退配额。常规路径按9/12/15Å等权平均新增窗口邻域面积选位点；beam按目标范围覆盖比例、组合内部低阶边际面积和单点风险代理建立Pareto前沿，每个位点组先一个代表，再轮换替换。前沿按覆盖比例、边际面积、重叠、直接面积、风险/先验/ID截断。全局物理身份去重后重算全合格表面分母，避免先保留的重叠窗口改变全局排序；科学假设配额及FoldX队列使用同一位点多样性策略。单点风险不是组合稳定性，未知单点风险不会按0参与Pareto比较。

几何支持以全表面图测地邻域和残基面积计算。双点边际为并集减最佳单点；三点减最佳双点，另存每个位点leave-one-out增益。保存各尺度覆盖位置、并集面积、边际面积、pairwise重叠及目标范围覆盖。全部是几何代理，不新增最低距离或最低边际面积科学门。无面积的旧API沿用原排序；物理coverage_scan专用生成/验证合同不变。

每个候选归入两个互斥类别之一：

- `EPITOPE_PRIORITIZED`：至少一个组成位点达到 DiscoTope seed 分位数或是局部峰值；
- `STABILITY_FALLBACK`：不要求 DiscoTope 支持，但仍必须通过完全相同的表面、硬排除和稳定性门控。

组合阶数解释为：

- 2 点：`WITHIN_PRIMARY_VALIDATION_SCOPE`；
- 3 点：`HIGHER_ORDER_EXTRAPOLATION_MEDIUM`；
- 4 点：`HIGHER_ORDER_EXTRAPOLATION_HIGH`。

高阶候选不能仅因为扰动分数更高而自动成为主推荐。

## 8. 候选级稳定性门控

组合共识状态包括：

- `CONSENSUS_PASS`
- `CONSENSUS_REVIEW`
- `MODEL_DISAGREEMENT`
- `SPURS_ONLY_REVIEW_REQUIRED`
- `FOLDX_ONLY_LEGACY`
- `COMBINATION_REJECT`
- `NOT_SCORED`

默认判定逻辑：

1. SPURS multi 和 FoldX 都存在，且两者差值超过 `disagreement_gate_kcal_mol`：`MODEL_DISAGREEMENT`；
2. SPURS multi 和 FoldX 同时通过各自配置门槛：`CONSENSUS_PASS`；
3. 两者都有结果但未同时通过：`CONSENSUS_REVIEW` 或 `COMBINATION_REJECT`；
4. 只有 SPURS multi：`SPURS_ONLY_REVIEW_REQUIRED`；
5. 只有 FoldX：`FOLDX_ONLY_LEGACY`；
6. 两者都没有：`NOT_SCORED`。

正式候选还必须满足：

- 每个组成单点有 SPURS single 分数，且分带为 `PREFERRED` 或 `ACCEPTABLE`；
- 没有硬排除；
- 组合状态为 `CONSENSUS_PASS`；
- 正式面板全部组合必须同时有有限的 SPURS multi 与 FoldX 实测模型输出；旧 `allow_spurs_only_final` 不能替代缺失证据；
- 常规3/4点组合比较实际突变所属显示分区范围内的9/12/15 Å面积覆盖比例均值，相对该范围最佳低阶几何对照至少增加 `minimum_coverage_gain`（默认5个百分点）；跨区组合使用所属分区并集。对照仅代表反事实几何，不声明对照已通过稳定性审查；缺失对照或几何时不能准入。全局hypothesis响应差仅保存为诊断，不与该面积门槛混用；
- 常规FoldX确认前检查上述已知几何门，第一轮专属机会优先双点，避免无新增覆盖三点消耗预算；物理覆盖预览仍沿用其独立合同；
- 候选必须来自连通 sector。
- native源链与规范化ECD统一仅用重原子（去H/D）、960点/1.4 Å探针及Tien/DSSP半径计算SASA；原始结构保留。真实源链/ECD分类变化仍审查，不能因为统一口径就豁免。ProCaliper参数不一致时仅作为交叉核对。所在sector无直径或连通性质量失败。

每个候选保存 `formal_eligible` 和 `formal_exclusion_reasons`，使“为什么没有进入正式面板”可以直接审查。

## 9. 面板优化

V4 主面板遵守独立的构建预算、候选资格、连通性和高阶数量限制。FoldX 复核预算独立，默认最多42个候选；显示分区数量是软目标，逐区建议不增加主面板预算。

正常流水线按新增直接表面覆盖面积 → 几何互补 → 信息增益 → 风险及确定性ID排序，solver为 `COVERAGE_FIRST_GREEDY`。直接面积按合格表面分母上 `geodesic_0` 位点并集计算，不作邻域膨胀。正覆盖增益不会被原信息/风险净目标否决；满足最小数量后没有新增覆盖则停止，原因 `NO_POSITIVE_COVERAGE_GAIN`。最小数量只能使用当前模式允许的候选，不能放宽科学门槛。缺少表面面积的兼容调用采用明确标注的假设响应代理。

每个显示分区按实际突变位置匹配建议，不使用生成窗口anchor断言归属。从通过者中先选9/12/15Å等权平均本区面积覆盖比例最高者，再比较内部边际、重叠、直接面积和风险。覆盖区域必须与实际分区取交集，邻区面积不能代替本区覆盖。没有通过者时先按独立未通过门数、缺失证据数、已知归一化超限量选最近者，再比较覆盖。保存每尺度本区覆盖比例、未覆盖残基及模型门槛。零面积/缺失支持不虚构覆盖比例。跨区组合可被多个分区共用，只投影一份物理构建体。

`PASSED`逐区建议不自动占用主面板名额；`NEAR_THRESHOLD`为待审查组合，不能获得正式资格或生成正式Stage2控制。安全候选全部不可构建时输出 `NO_COMBINATION`和原因。网页、JSON、CSV、离线报告和ZIP使用同一最终建议列表。旧无几何/假设的API保留系列选择兼容路径。

直接面积与邻域互补分别报告，二者均不能证明生物学扰动或抗体结合覆盖。实现、批判性取舍及验证见 [分区内优化验收](intra_sector_design_optimization_20261008.md)。FoldX最大60，主面板预算独立控制；全局假设配额的保留数仍可超过请求的软容量，保留更多不同组合可能增加SPURS multi工作量，需查看有效候选数。


## 10. 运行模式与状态

### 10.1 strict

只有当 ProCaliper、SPURS single、SPURS multi 的整个请求成功且输出完整、FoldX 按配置完成全部短名单（显式开启 `allow_spurs_only_final` 的例外见候选门控），并且最终选中的每个候选都通过候选级门控时，顶层状态才可能是：

```text
FORMAL_RECOMMENDATION
```

以下任一情况会导致：

```text
BLOCKED_EXTERNAL_EVIDENCE
```

- 稳定性必需工具未配置、超时、崩溃或输出非法；
- 必需证据缺失；
- 没有严格合格的单点或组合候选；
- 选中候选存在模型分歧或缺少 required FoldX；
- 结构/WT/残基坐标映射失败。

### 10.2 preview

preview 顶层状态固定为：

```text
PREVIEW_NOT_FOR_FORMAL_RECOMMENDATION
```

preview 可以保留 `NOT_SCORED` 单点并生成调试候选，但旧版兼容设计状态为 `PREVIEW_CANDIDATE`，不能被 UI 或汇总计为正式推荐。preview 不填充伪造 kcal/mol。

## 11. Worker 接口

worker 通过 stdin 接收一个 JSON 请求，通过 stdout 返回一个 JSON 对象，日志只能写 stderr。stdout 不能包含进度文本、调试信息或多个 JSON 文档。

请求的核心字段包括：

```json
{
  "schema_version": "1.0",
  "request_id": "...",
  "structure_path": "...",
  "chain_id": "A",
  "structure_type": "solved",
  "sequence": "...",
  "mutations": [["D12A", "W18A"]],
  "output_dir": "...",
  "device": "auto",
  "parameters": {}
}
```

返回对象至少包含：

```json
{
  "status": "success",
  "records": [],
  "warnings": [],
  "errors": [],
  "provenance": {},
  "runtime_seconds": 0.0
}
```

`status` 只能是 `success`、`partial` 或 `failed`。V4 会对不同工具执行专门校验：表面工具检查全链覆盖，SPURS/FoldX 检查突变组合、位置和 WT，所有数值必须是有限值。

工具命令可以写入配置，也可以使用：

```text
ECD_V4_PROCALIPER_COMMAND
ECD_V4_SPURS_SINGLE_COMMAND
ECD_V4_SPURS_MULTI_COMMAND
ECD_V4_DISCOTOPE_COMMAND
ECD_V4_FOLDX_COMMAND
```

环境变量值必须是 JSON argv 数组。`ECD_V4_SPURS_COMMAND` 可作为 single/multi 的兼容后备。

### 11.1 当前工作站的 GPU 执行路径

Windows 主进程通过 `integrations/wsl_bridge.py` 把规范化结构路径和输出目录转换为 WSL 路径，再调用隔离的 Linux Python worker。默认配置为：

| 证据线              | 执行环境                                                        | 设备                 |
| ---------------- | ----------------------------------------------------------- | ------------------ |
| SPURS single     | `/home/segayi1/miniconda3/envs/spurs_v4/bin/python`         | CUDA               |
| SPURS multi      | `/home/segayi1/miniconda3/envs/spurs_v4/bin/python`         | CUDA               |
| DiscoTope-3.0    | `/home/segayi1/miniconda3/envs/discotope_v4/bin/python`     | CUDA（预测链 <1000 aa） |
| ProCaliper 0.2.1 | `G:\anaconda\envs\ecd_procaliper\python.exe`，本地 JSON worker | CPU                |
| FoldX 5.1        | `vendor/foldx/foldx_1_20270131.exe`，本地 JSON worker          | CPU，最多 8 路并行       |

GPU worker 在请求指定 `device=cuda` 但 CUDA 本身不可用时仍返回失败。对于运行期间的 CUDA OOM，SPURS single/multi 和 DiscoTope 共用 WSL 文件锁，确保同一时刻只有一个 worker 占用 GPU；worker 会释放 Python/CUDA 缓存并重试 1 次，第二次仍 OOM 时才按 `allow_cpu_fallback_on_oom` 显式回退到 CPU。该回退不是静默降级：结果会加入 `CUDA_OOM_RETRY_1`、`CUDA_OOM_CPU_FALLBACK` warning，并在 provenance 中记录请求设备、实际设备、重试次数、回退原因和 GPU 锁等待时间。默认超时为 SPURS single 900 秒、DiscoTope 900 秒、SPURS multi 1800 秒，使 CPU 回退有足够时间完成，同时仍保持有界执行。

DiscoTope 遵循其 IF1 实现的长度边界：预测链达到 1000 aa 时，V4 的 CUDA 请求会显式失败；如需处理该长度，应把该工具请求的设备改为 `cpu`。SPURS 依赖产生的普通提示和 warning 被重定向到 stderr，stdout 始终只保留一个 JSON 文档。

每个成功结果的 provenance 保存 worker 名称、PyTorch 版本、CUDA 版本、请求设备、实际设备、GPU 名称、OOM 重试次数、CPU 回退原因和 GPU 锁等待时间；SPURS 还记录推理结束时的 CUDA allocated/reserved bytes。模型 worker 计算实际权重或二进制摘要并写入 provenance。当前无法在查缓存前独立验证模型文件摘要，因此 SPURS、DiscoTope 与 FoldX 禁用结果缓存，以免更换本地文件后误用旧评分；ProCaliper 的缓存键仍包含结构内容、目标链、序列、结构类型、工具参数和请求设备。

ProCaliper worker 通过包内的 Bio.PDB ShrakeRupley 实现计算全链 SASA（探针半径 1.4 Å、100 个球面点），再以 Tien 等人的最大 ASA 表换算 rSASA；它会核对规范化结构链的序列和返回长度。FoldX worker 先运行 `RepairPDB`，并按输入结构与 FoldX 二进制哈希复用 `workspace_data/foldx_repair_cache/` 中的修复结构；随后对候选运行 `BuildModel --numberOfRuns=5` 并取中位数。为控制严格模式成本，SPURS multi 仍评估完整候选集，而 FoldX 只复核最多 60 个确定性短名单候选；短名单先在每个 sector 内优先选择二点、高表位倾向、低结构风险候选，再按 sector 轮询抽取，避免单个大 sector 占满预算。FoldX 必须返回短名单全部组合，缺失任一结果即视为该证据线失败。

本机以 290 aa PD-L1 结构完成的冒烟测试结果为：SPURS single 约 58 秒、SPURS multi 单个双突变组合约 60 秒、DiscoTope-3.0 约 39 秒，设备均为 NVIDIA RTX 3060 Laptop GPU、PyTorch 2.5.1+cu124。另一次 215 aa 严格模式完整运行中，ProCaliper 处理 215 个残基约 2.28 秒，FoldX 对 42 个候选各运行 5 次约 1077.55 秒；5 条证据线均可用，顶层状态为 `FORMAL_RECOMMENDATION`，最终选择 4 个 PRIMARY 候选。OOM 恢复修复后，CTLA-4 的 125 个结构残基也完成了 SPURS single CUDA 冒烟测试，约 71.6 秒并返回全部 125 条记录。该数字只用于运行量级估算；冷启动、序列长度、结构大小、候选组合数量、显存/CPU 占用、GPU 锁等待和缓存命中都会改变总耗时。

网页服务在启动时加载配置。修改 `config.py`、bridge 或 worker 后必须停止并重新启动服务，已经运行的旧进程不会自动获得新配置。每次正式运行应检查 `tool_inventory.json` 以及 worker provenance，不能只凭系统存在 NVIDIA GPU 推断模型实际使用了 GPU。

## 12. 输出与审计

V4 主要输出如下：

| 文件                                   | 内容                                                                      |
| ------------------------------------ | ----------------------------------------------------------------------- |
| `v4_manifest.json`                   | 经 JSON Schema 校验的模式、工具、参数、上下文、代码/结构哈希、实际 worker 状态与 provenance、命令、警告和错误 |
| `tool_inventory.json`                | 工具可用性、命令、版本、模型 revision、权重摘要和失败原因                                       |
| `structure_contexts.json`            | biological context 与 prediction context                                 |
| `prediction_context.pdb`             | 连续编号的预测目标链                                                              |
| `residue_mapping.csv`                | normalized/UniProt/author/WT 权威映射                                       |
| `surface_evidence.csv`               | native 与 ProCaliper 表面证据比较                                              |
| `sasa_comparison.csv`                | rSASA 差异与一致性状态                                                          |
| `spurs_single_matrix.csv.gz`         | 单点完整评分矩阵                                                                |
| `single_mutation_options.csv`        | 单点分带、化学类别、资格和警告                                                         |
| `discotope_residues.csv`             | DiscoTope 分数、分位数和岛屿                                                     |
| `epitope_islands.json`               | 表位先验岛屿                                                                  |
| `probe_sectors_v4.json`              | V4 sector 及质量标记                                                         |
| `panel_candidates.json`              | 所有候选、覆盖、稳定性和排除原因                                                        |
| `sector_recommendations_v4.json/csv` | 每个显示分区的组合建议、资格、门槛及失败原因                                                  |
| `candidate_combinations.csv.gz`      | 候选组合审计表                                                                 |
| `combination_stability.csv`          | SPURS multi/FoldX/ddG/共识状态                                              |
| `selected_panel_v4.json`             | 最终面板与未覆盖 sector                                                         |
| `v4_scientific_report.md`            | 面向人工审查的科学报告                                                             |

旧版 `regions.json` 和 `designs.json` 是兼容投影。V4 解释应以 V4 原生输出为准，不能仅凭旧版字段名把 sector 解释为真实抗体 footprint。

## 13. 可重复性与验证

V4 运行应保存：

- 输入结构和序列摘要；
- 结构链和位置映射；
- 全部 V4 参数；
- 工具命令、版本、模型 revision、权重摘要和设备；
- cache key；
- 每个 worker 的 warnings/errors/runtime；
- 候选排除原因和最终面板状态。

当前软件回归测试覆盖：worker 超时、worker 崩溃、stdout 污染、single/multi 配置隔离、预览状态隔离、strict 候选级门控、模型分歧阻断、每补丁优先/回退候选系列、DiscoTope 缺失时的审计降级、规范化位置 sidecar 和完整表面分区。测试验证的是软件契约与确定性逻辑，不替代 ProCaliper、SPURS、DiscoTope 或 FoldX 自身的外部 benchmark。固定预算、抗原家族隔离的比较方案见 `docs/joint_v4_validation_protocol.md`；目前没有由该方案产生的实测性能结论。

## 14. 科学解释边界

- **不具抗体特异性：**&#x5F53;前输入不包含目标抗体序列、抗原–抗体复合物或对接姿势。DiscoTope 先验不能回答指定抗体是否结合某 patch。
- **几何可及性不等于细胞可达：**`cell_surface` 模式目前不模拟膜平面上的 Fab/IgG 体积排斥、糖萼或全长蛋白邻近结构；孤立 ECD 表面的残基不能直接称为活细胞可达表位。
- **计算稳定性不是表达/折叠实测：**&#x53;PURS 与 FoldX 的共识只代表指定结构和模型中的预测一致，不能证明突变体表达、膜定位、构象或功能正常。
- **实验标签稀疏与未知表位：**&#x5DF2;知抗体复合物只标出被观察到的表位，不代表其他未观察表面绝非表位。Cia 等对 250 余个抗体-抗原结构评估九种上一代构象表位工具，报告总体性能偏低且部分接近随机表面 patch；该研究早于 DiscoTope-3.0，不能直接作为 V4 或 DiscoTope-3.0 的性能估计，但说明独立基准和盲测不可省略。评估必须按抗原家族隔离训练/测试并报告覆盖与失败率。
- **尚无 v4 性能优势的实验结论：**&#x5F53;前测试覆盖软件契约和确定性逻辑。结构接触命中只能用于开发期代理评估；优越性需固定预算外部对照和前瞻实验（表达/膜展示、构象对照、目标抗体结合）支持。

组会汇报建议把结论表述为：**v4 把残基级表位倾向、表面几何和突变稳定性约束整合成一个可审计的实验候选面板生成流程；它的贡献在流程整合、候选可解释性和失败门控，当前尚未证明比成熟方法更准确地定位真实抗体表位。**

DiscoTope 高分是表位倾向先验；SPURS/FoldX ddG 是模型预测；`CONSENSUS_PASS` 只是配置意义上的计算共识。正式实验仍需 WT、表达/膜表面展示、构象敏感性、折叠或稳定性以及必要的功能对照。preview 和 `BLOCKED_EXTERNAL_EVIDENCE` 结果不得作为正式实验放行依据。

对于 `cell_surface` 测定，V4 当前没有显式 Fab/IgG 与膜平面的可达性建模；结果会提示这一限制，不能将孤立 ECD 的表面暴露直接解释成活细胞抗体可达。

## 相关文献

1. Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms. *J Mol Biol*. 1973. <https://doi.org/10.1016/0022-2836(73)90011-9>
2. Tien MZ, et al. Maximum Allowed Solvent Accessibilites of Residues in Proteins. *PLOS ONE*. 2013. <https://doi.org/10.1371/journal.pone.0080635>
3. Høie MH, et al. DiscoTope-3.0: improved B-cell epitope prediction using inverse folding latent representations. *Front Immunol*. 2024. <https://doi.org/10.3389/fimmu.2024.1322712>
4. Hsu C, et al. Learning inverse folding from millions of predicted structures. *ICML/PMLR*. 2022. <https://proceedings.mlr.press/v162/hsu22a.html>
5. Li Z, Luo Y. Generalizable and scalable protein stability prediction with rewired protein generative models. *Nat Commun*. 17:891 (2026; online 2025). <https://doi.org/10.1038/s41467-025-67609-4>
6. Lin Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*. 2023. <https://doi.org/10.1126/science.ade2574>
7. Dauparas J, et al. Robust deep learning–based protein sequence design using ProteinMPNN. *Science*. 2022. <https://doi.org/10.1126/science.add2187>
8. Guerois R, Nielsen JE, Serrano L. Predicting changes in the stability of proteins and protein complexes: a study of more than 1000 mutations. *J Mol Biol*. 2002. <https://doi.org/10.1016/S0022-2836(02)00442-4>
9. Schymkowitz J, et al. The FoldX web server: an online force field. *Nucleic Acids Res*. 2005. <https://doi.org/10.1093/nar/gki387>
10. Cia G, Pucci F, Rooman M. Critical review of conformational B-cell epitope prediction methods. 2023. [<https://doi.org/10.1093/bib/bbac567\\\](>
