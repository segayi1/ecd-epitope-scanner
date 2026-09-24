# ECD Epitope Scanner Joint V4-ML 详细实施规格

> 状态：Draft for decisions and subsequent code generation  
> 目标：在现有 `joint_v3` 基础上接入 ProCaliper、SPURS 和 DiscoTope-3.0，生成不超过七个空间连通、稳定性风险受控、对表面表位热点具有互补覆盖的组合突变探针。  
> 非目标：从抗原结构直接断言真实抗体表位，或通过计算“确保”突变体正确折叠。  
> 推荐引擎名：`joint_v4_ml`，在完成基准和实验校准前保持 opt-in。

## 1. 决策摘要

V4-ML 不重写现有几何与优化框架，而是在 `joint_v3` 上增加三个独立证据适配器和一套版本化证据融合规则：

1. ProCaliper提供残基 SASA、距离网络及结构特征的独立计算和审计。
2. 现有显式溶剂表面点图继续作为连通性、测地距离和表面覆盖的唯一权威几何来源。
3. SPURS single 对全部映射位点进行一次饱和突变稳定性预测，生成残基-替换级安全菜单。
4. DiscoTope-3.0 对完整目标抗原链计算逐残基表位倾向，只作为排序和分区软目标。
5. 图约束分区在 `K=1..7` 中选取最合理的连通区域数，不强迫每个结构生成七区。
6. 每区使用确定性 beam search 生成 2–4 点组合，并计算多尺度几何覆盖。
7. SPURS multi 对短名单组合评分；最终少量候选使用 FoldX作正交复核。
8. 全局优化器在构建预算内选择覆盖互补、稳定性风险低、跨区干扰少的主面板。
9. 所有外部模型失败、缺失、外推和降级必须在顶层结果中可见。

## 2. 待确认决策

以下选项不会阻止计划编写。未确认前，代码生成按标注的推荐默认值设计成可配置行为。

### D1：新引擎的发布方式

- A（推荐）：新增 `joint_v4_ml`，保留 `surface_v2` 默认和 `joint_v3` 回归入口；通过基准后再切换默认。
- B：直接将 `joint_v3` 升级为新语义，不新增引擎名。
- C：完全替换旧引擎，只保留 V4-ML。

影响：A 的兼容性和可比较性最好；B 会导致同一引擎名在不同版本产生不同科学语义；C 风险最高。

### D2：七个分区与七个实验构建是否为同一预算

- A（推荐）：最多七个 `ProbeSector`，每区最多一个主构建；备选和单点去卷积不计入主面板。
- B：最多七个分区，但主构建总数可少于或多于七个。
- C：分区数不受七限制，仅最终主构建不超过七个。

推荐 A，因为实验解释最直接：一个主探针对应一个粗定位区域。

### D3：ProCaliper在几何层的权限

- A（推荐）：ProCaliper负责 SASA/距离特征和交叉审计；现有显式表面图负责连通与测地距离。
- B：ProCaliper残基距离图也参与候选邻接，但必须经过表面法向量、净空和同表面分量过滤。
- C：完全用ProCaliper距离图替换现有表面图。

C 不建议实施，因为原始欧氏距离图可能跨越蛋白内部或连接相反表面。

### D4：外部模型缺失时的运行政策

- A（推荐）：提供 `strict` 和 `preview` 两种模式。`strict` 缺少 SPURS 或 DiscoTope 即停止；`preview` 允许降级，但不输出正式推荐。
- B：始终允许降级，并使用旧代理分数生成推荐。
- C：任何外部模型缺失均终止整个分析，包括结构和表面结果。

### D5：组合稳定性共识

- A（推荐）：SPURS multi为主评分；最终候选必须通过FoldX，或显式标记 `SPURS_ONLY_REVIEW_REQUIRED`。
- B：只使用SPURS，不再运行FoldX。
- C：SPURS与FoldX都作为硬门槛，任何一方失败都淘汰。

A 能利用SPURS效率，同时保留与训练范式不同的正交证据。

### D6：三点和四点组合政策

- A（推荐）：二至四点均允许；三/四点必须带 `HIGHER_ORDER_EXTRAPOLATION`，且只有覆盖明显改善并通过FoldX才可成为主推荐。
- B：正式推荐最多二点；三/四点仅作为探索性备选。
- C：二至四点同等处理。

SPURS multi公开训练和主要验证集中在双突变，因此不建议 C。

### D7：DiscoTope高倾向的分区约束

- A（推荐）：高倾向是软目标；无法满足时允许低倾向区或少于七区，并报告原因。
- B：每区必须包含一个超过固定阈值的高倾向残基，否则不生成该区。
- C：只围绕高倾向残基分区，不追求完整可设计表面覆盖。

### D8：突变字母表

- A（推荐）：先由SPURS生成全20种替换，再由化学扰动、保护规则和稳定性选择每个位点最多3种。
- B：维持现有定向字母表，只用SPURS为其评分。
- C：组合阶段保留全部19种非WT替换。

A 既保留发现能力，又控制组合数量。

### D9：运行环境

- A（推荐）：主程序、SPURS worker、DiscoTope worker使用三个隔离环境，通过JSON/CSV和子进程协议通讯。
- B：尝试把所有依赖安装进现有Python 3.11环境。
- C：调用远程Web/API服务。

A 最适合当前Windows本地应用，也便于固定模型版本。C 会引入数据上传、服务可用性和复现问题。

### D10：DiscoTope许可与产品用途

DiscoTope-3.0代码采用 CC BY-NC 4.0。若项目可能用于商业或收费场景，需要在实施前确定授权路径。

- A（推荐，非商业研究）：本地独立安装，记录引用和许可，不把其源码或权重打包进本项目。
- B（商业可能）：先取得额外授权；适配器先实现，正式模型调用保持禁用。
- C：不用DiscoTope代码，改用许可兼容的替代表位倾向模型。

## 3. 科学声明边界

### 3.1 系统可以声明

- 残基满足所配置的表面可及性标准。
- 某个具体替换在SPURS/FoldX下属于较低或较高预测稳定性风险。
- 某个残基具有较高或较低DiscoTope相对倾向。
- 某个组合在指定表面图和覆盖情景下覆盖多少可设计表面积。
- 某个面板相对其他候选具有更好的预测覆盖、稳定性和互补性。

### 3.2 系统不得声明

- 高DiscoTope分数残基就是真实表位残基。
- 低分区域不是抗体表位。
- SPURS/FoldX通过即可确保正确折叠、表达或细胞表面展示。
- 一个突变能物理影响整个分区。
- 三点或四点SPURS预测具有与双点相同的验证水平。
- 预测的低结合信号一定来自直接表位破坏，而不是变构、糖基化或表达变化。

## 4. 总体架构

```text
PDB/mmCIF + sequence/context
  |
  +-- StructureAudit / residue mapping / confidence / protection
  |
  +-- StructureContextBuilder
  |     +-- isolated_target_chain.pdb
  |     `-- biological_context.pdb
  |
  +-- ProCaliperAdapter
  |     +-- SASA/rSASA
  |     `-- residue distance features
  |
  +-- Existing Surface V2 geometry
  |     +-- explicit solvent surface point graph
  |     +-- residue surface graph
  |     `-- geodesic distances and areas
  |
  +-- SPURSSingleAdapter
  |     `-- L x 20 mutation stability matrix
  |
  +-- DiscoTopeAdapter
  |     `-- per-residue propensity and calibrated/rank scores
  |
  +-- EvidenceFusion
  |     `-- allowed MutationOptions + exclusion evidence
  |
  +-- EpitopeAwareSectorBuilder (1..7 connected sectors)
  |
  +-- CandidateSearch (2..4 mutations, deterministic beam)
  |
  +-- CoverageEngine (direct/geodesic/footprint scenarios)
  |
  +-- SPURSMultiAdapter -> FoldX confirmation
  |
  +-- PanelOptimizer
  |
  `-- reports / API / UI / audit manifest
```

## 5. 目录与模块变更

推荐新增：

```text
ecd_epitope_scanner/
  joint_design_v4/
    __init__.py
    models.py
    engine.py
    contexts.py
    evidence_fusion.py
    epitope_sectors.py
    candidate_ranking.py
    panel_optimizer.py
    reporting.py
  integrations/
    __init__.py
    base.py
    process_runner.py
    procaliper_adapter.py
    spurs_adapter.py
    discotope_adapter.py
    foldx_adapter.py
  schemas/
    joint_v4_manifest.schema.json
```

推荐复用而不是复制：

- `surface_patch/*`：显式表面点云、表面图、测地距离和面积。
- `joint_design/structure_audit.py`：结构审计。
- `joint_design/topology.py`：连接端和膜方向不确定性。
- `joint_design/microfootprints.py`：多尺度窗口。
- `joint_design/candidate_search.py`：确定性beam search核心。
- `joint_design/coverage.py`：覆盖矩阵。
- `joint_design/optimizer.py`：MILP和greedy fallback。
- `stability.py`：FoldX执行与解析，可在统一adapter后继续调用。

不建议直接在 `joint_design` 内大量条件分支堆叠V4字段。V3结果需要继续可读和可回归，V4使用新schema，再通过legacy adapter投影到当前UI所需字段。

## 6. 统一外部工具协议

### 6.1 Adapter接口

```python
class ExternalPredictor(Protocol):
    name: str

    def inventory(self) -> ToolInventory: ...
    def validate_input(self, request: PredictionRequest) -> ValidationResult: ...
    def predict(self, request: PredictionRequest) -> PredictionResult: ...
```

`ToolInventory`至少包含：

- `available`
- `executable`
- `environment_type`
- `version`
- `git_revision`
- `model_revision`
- `weights_sha256`
- `device`
- `license_id`
- `failure_reason`

### 6.2 子进程JSON合同

主程序通过stdin或临时请求文件传入：

```json
{
  "schema_version": "1.0",
  "request_id": "...",
  "structure_path": "absolute path",
  "chain_id": "A",
  "structure_type": "solved|alphafold|unknown",
  "mutations": [["W31A", "K55E"]],
  "output_dir": "absolute path",
  "device": "auto"
}
```

worker必须返回：

```json
{
  "status": "success|partial|failed",
  "records": [],
  "warnings": [],
  "errors": [],
  "provenance": {},
  "runtime_seconds": 0.0
}
```

禁止从stdout混入进度日志；日志写stderr或独立文件。所有超时、非零退出码、缺失输出和残基数不一致都必须转换为稳定错误码。

### 6.3 缓存键

缓存键必须至少包含：

```text
sha256(
  normalized_structure_bytes
  + selected_chain
  + sequence
  + structure_type
  + tool_version
  + model_revision
  + relevant_parameters
)
```

单点SPURS和DiscoTope结果可跨候选组合复用；不得仅以文件名作为缓存键。

## 7. 数据模型

### 7.1 残基标识

内部主键必须保留完整结构身份：

```python
@dataclass(frozen=True, order=True)
class ResidueKey:
    model_id: str
    chain_id: str
    auth_seq_id: int
    insertion_code: str
    label_seq_id: int | None
    uniprot_pos: int | None
    one_letter: str
```

所有外部工具的连续位置编号必须通过显式映射表转换为 `ResidueKey`。映射出现缺失、重复或WT氨基酸不一致时，该记录不得进入自动设计。

### 7.2 新增记录

```python
@dataclass
class SurfaceEvidence:
    residue: ResidueKey
    procaliper_sasa_a2: float | None
    procaliper_rsasa: float | None
    native_sasa_a2: float | None
    native_rsasa: float | None
    sidechain_sasa_a2: float | None
    mapped_surface_area_a2: float
    surface_component_id: int | None
    sasa_agreement_status: str
    provenance: dict[str, Any]

@dataclass
class EpitopeEvidence:
    residue: ResidueKey
    raw_score: float | None
    calibrated_score: float | None
    surface_percentile: float | None
    local_peak: bool
    neighborhood_mean_by_radius: dict[str, float]
    status: str
    provenance: dict[str, Any]

@dataclass
class SingleMutationEvidence:
    residue: ResidueKey
    mutant_aa: str
    spurs_ddg_kcal_mol: float | None
    stability_band: str
    chemistry_class: str
    hard_exclusions: list[str]
    soft_warnings: list[str]
    eligible: bool
    provenance: dict[str, Any]

@dataclass
class CombinationStabilityEvidence:
    mutations: tuple[str, ...]
    spurs_multi_ddg_kcal_mol: float | None
    foldx_ddg_kcal_mol: float | None
    order: int
    extrapolation_level: str
    consensus_status: str
    disagreement: float | None
    provenance: dict[str, Any]
```

### 7.3 分数语义

每个分数都必须同时存储：

- 原始值；
- 单位；
- 方向定义；
- 模型和权重版本；
- 是否经过校准；
- 是否属于训练/验证范围内使用；
- 缺失原因。

不得用 `0` 表示未运行或失败。

## 8. 阶段A：输入、结构上下文与审计

### 8.1 输入标准化

1. 解析PDB/mmCIF全部model、chain、altloc、插入码和非标准残基。
2. 默认选第一个model，但写入警告；NMR ensemble暂不自动平均。
3. 建立结构编号、连续链编号、输入序列和UniProt编号的双向映射。
4. 为外部工具导出规范化单链PDB，连续编号从1开始，并保存sidecar映射CSV。
5. 另保留生物学装配体或用户上传复合体，用于遮挡、SASA和界面保护。

### 8.2 两种结构上下文

- `prediction_context`：目标抗原单链，用于DiscoTope和SPURS，匹配其常见训练/推理方式。
- `biological_context`：目标链加相关亚基、配体或糖链，用于计算真实遮挡和保护界面。

若无法确定生物学装配体，状态记为 `BIOLOGICAL_ASSEMBLY_UNKNOWN`，不得把晶体接触自动当作生理界面。

### 8.3 硬保护规则

默认不可突变：

- 二硫键Cys和配对残基；
- 创建或破坏N-X-S/T糖基化sequon的替换；
- 已知糖基化、切割、催化、配体结合和关键PTM位点；
- 主链不完整、非标准或残基映射不唯一；
- 低置信度/高度柔性且无法可靠解释的区域；
- 天然寡聚界面的关键埋藏残基；
- 用户显式保护位点；
- 结构数据不足以支持SPURS输入的位点。

保守但非绝对的规则放入软惩罚，不应全部硬删除，例如高进化保守性、膜近端、功能界面邻域。

## 9. 阶段B：ProCaliper与表面图

### 9.1 ProCaliper输出

需要保存：

- 原子和残基SASA；
- rSASA及所用归一化尺度；
- 可选电荷；
- 全残基距离矩阵或阈值后的稀疏边；
- 运行参数和版本。

不要求把ProCaliper对象pickle作为长期交换格式；长期格式使用CSV/JSON/NPZ。

### 9.2 SASA一致性检查

对当前Shrake–Rupley与ProCaliper结果进行逐残基比较：

- WT残基和编号必须完全对齐；
- 记录绝对差、相对差和Pearson/Spearman相关；
- 初始审计阈值建议为 `abs(rsasa_delta) <= 0.05` 或 `relative_sasa_delta <= 20%`；
- 超阈值不立即淘汰残基，而是标记 `SASA_DISAGREEMENT`；
- 若全链系统性偏差，停止自动设计并报告半径、采样点、归一化尺度或结构上下文差异。

阈值在基准阶段冻结前属于待校准启发式。

### 9.3 表面候选条件

推荐初始条件：

```text
native_rSASA >= 0.25
AND sidechain_SASA >= 20 A2
AND mapped_surface_area >= 8 A2
AND explicit_surface_component exists
```

`0.20 <= rSASA < 0.25` 可进入 `BORDERLINE_SURFACE`，仅在严格候选不足时参与备选。

### 9.4 图权威关系

- 结点：通过表面条件或作为连接桥保留的可靠残基。
- 边：必须由显式溶剂表面点图支持，并满足现有CA长边限制。
- 权重：表面测地距离。
- ProCaliper欧氏距离：仅作为紧凑性、跨区干扰和近邻特征，不作为连通证明。

## 10. 阶段C：SPURS单点模式

### 10.1 执行

1. 对规范化目标单链运行一次single模型。
2. 得到 `L x 20` 预测矩阵。
3. 验证WT列已归零或符合模型输出合同。
4. 映射回 `ResidueKey`。
5. 保存完整矩阵，而不是只保存通过的突变。

### 10.2 初始稳定性分带

假设正值表示去稳定化：

```text
PREFERRED: ddG <= 0.0
ACCEPTABLE: 0.0 < ddG <= +0.5 kcal/mol
REVIEW: +0.5 < ddG <= +1.0 kcal/mol
REJECT: ddG > +1.0 kcal/mol
```

规则：

- 正负方向必须在adapter自检和固定fixture中验证，不能只依赖变量名。
- `strict`模式默认只让 `PREFERRED` 和 `ACCEPTABLE` 进入主候选池。
- `REVIEW`只在某区没有足够候选时进入备选池。
- 阈值必须配置化，并在实验校准后更新版本。

### 10.3 替换菜单压缩

每个位点从19个非WT替换中选择最多3个：

1. 删除硬保护和SPURS `REJECT`。
2. 按化学扰动类型分组：体积减小、体积增大、电荷改变、芳香性改变、极性改变。
3. 优先保留不同扰动类型，避免三个高度相似替换。
4. 在同类型中按SPURS风险、扰动强度和可合成性排序。
5. 保存所有被删除替换的reason code。

## 11. 阶段D：DiscoTope-3.0

### 11.1 执行与映射

1. 对完整规范化目标链运行，而不是只提交安全表面残基。
2. 根据结构来源设置 `solved` 或 `alphafold`；未知来源不得伪装。
3. 验证输出残基数、WT氨基酸和连续编号。
4. 将原始分数、官方校准分数和二元标签全部保存。
5. 设计逻辑不使用官方二元标签作硬门槛。

### 11.2 派生特征

仅在可靠表面残基集合上计算：

- `surface_percentile`；
- top 10%、20%、30%标记；
- 是否为表面图局部极大值；
- 9/12/15 Å测地邻域加权均值；
- 高倾向连通岛ID；
- 高倾向岛总表面积。

### 11.3 高倾向种子

推荐把以下任一条件作为seed：

```text
surface_percentile >= 0.80
OR local_peak AND surface_percentile >= 0.65
```

seed只是分区目标，不是候选资格。阈值须在公开复合物benchmark中做敏感性分析。

## 12. 阶段E：表位倾向感知的连通分区

### 12.1 K的选择

对 `K=1..min(7, connected_designable_components)` 生成候选分区方案，使用以下字典序目标选择：

1. 所有sector严格连通，违反则方案不可行；
2. 最大化被至少一个sector纳入的可设计表面积；
3. 最大化被sector捕获的高倾向岛权重；
4. 最小化高倾向must-link边的切割权重；
5. 最小化最大未覆盖连通表面积；
6. 平衡sector可设计面积，避免极小区；
7. 最小化跨sector潜在扰动重叠；
8. 在近似等价时选择更少sector。

不建议仅通过最大化模块度自动决定K，因为模块度存在分辨率限制，且未直接表达实验预算。

### 12.2 Must-link与cannot-link

Must-link候选：

- 位于同一高倾向连通岛；
- 表面测地距离不超过初始12 Å；
- 之间存在可靠表面路径；
- 合并后不超过sector直径上限。

Cannot-link：

- 不同表面连通分量；
- 合并后测地直径 >35 Å；
- 合并后CA欧氏跨度 >30 Å；
- 相反表面或显式穿模检查失败；
- 被结构域间柔性断点分隔且没有可靠表面连接。

### 12.3 无高倾向区域

若某个重要可设计表面分量没有高倾向seed：

- 仍可建立 `COVERAGE_ONLY` sector；
- 其优先级低于包含高倾向岛的sector；
- UI和报告必须显示 `LOW_EPITOPE_PRIOR`；
- 若七个预算已满，优先保留高倾向和大未覆盖面积，而不是强制完整分区。

## 13. 阶段F：候选组合生成

### 13.1 位点评分

单个位点不生成不可解释的总概率，保存组成项：

```text
epitope_rank_component
local_peak_component
direct_surface_area_component
geodesic_centrality_component
chemical_perturbation_component
single_stability_penalty
functional_risk_penalty
uncertainty_penalty
```

用于beam排序的初始启发式分数可以是标准化线性组合，但所有权重必须配置化并经消融实验。

### 13.2 几何限制

- 每个组合2–4个不同位点；
- 位点间建议最小CA距离5.5 Å；
- 最大表面测地跨度24 Å；
- 最大欧氏跨度22–24 Å；
- 所有位点必须属于同一sector；
- 组合在sector诱导图上的局部影响区域必须连通；
- 不允许创建新的糖基化sequon或破坏硬保护位点。

### 13.3 Beam search配额

推荐初始值：

```text
site_pool_per_sector = 12
mutation_options_per_site = 3
beam_width_per_order = 300
pre_spurs_multi_shortlist_per_sector = 150
post_spurs_multi_shortlist_per_sector = 30
foldx_confirmation_per_sector = 6
```

每个突变阶数使用独立配额，避免四点组合因为启发式总分较高而吞掉二点候选。

### 13.4 覆盖情景

每个组合至少计算：

- `direct_site`；
- `geodesic_9a`；
- `geodesic_12a`；
- `geodesic_15a`；
- `microfootprint_hit_at_least_1`；
- `microfootprint_hit_at_least_2`；
- `high_propensity_island_coverage`。

这些是实验设计代理，不解释为实际影响半径。

## 14. 阶段G：SPURS多点与FoldX共识

### 14.1 SPURS multi批处理

- 同一蛋白的基础表示只计算一次。
- 候选按突变数和sector批处理。
- 保存组合预测及其 constituent single sum。
- 计算 `epistasis_estimate = multi_ddg - sum(single_ddg)`。
- 极端epistasis值进入 `MODEL_OUTLIER_REVIEW`，不静默接受。

### 14.2 组合阶数不确定性

```text
order=2 -> WITHIN_PRIMARY_VALIDATION_SCOPE
order=3 -> HIGHER_ORDER_EXTRAPOLATION_MEDIUM
order=4 -> HIGHER_ORDER_EXTRAPOLATION_HIGH
```

三/四点组合只有在下列条件全部满足时进入主推荐：

1. 比最佳低阶候选产生预先规定的最小覆盖增益；
2. SPURS multi通过配置门槛；
3. 每个组成单点均不高于单点门槛；
4. FoldX复核通过；
5. 不存在功能或糖基化风险升级。

初始最小覆盖增益建议为5个百分点，必须在benchmark后校准。

### 14.3 共识状态

```text
CONSENSUS_PASS
CONSENSUS_REVIEW
MODEL_DISAGREEMENT
SPURS_ONLY_REVIEW_REQUIRED
FOLDX_ONLY_LEGACY
COMBINATION_REJECT
NOT_SCORED
```

初始建议门槛：

- SPURS multi主推荐：≤ +0.5 kcal/mol；+0.5至+1.0进入review；
- FoldX主推荐：沿用项目现有组合门槛 +2.0 kcal/mol，但报告其较宽松性质；
- 两模型差值超过1.5 kcal/mol时标记模型分歧。

这些门槛均属于待校准配置，不能硬编码在adapter中。

## 15. 阶段H：全局面板优化

### 15.1 硬约束

- 主构建数不超过配置预算，推荐默认7；
- 每个sector最多选择一个PRIMARY；
- 候选必须满足运行模式所要求的稳定性证据；
- 不能选择几乎相同的位点集合；
- 三/四点外推候选受到独立数量上限，推荐主面板最多2个；
- 选择的候选必须来自有效连通sector。

### 15.2 字典序目标

1. 最大化包含高倾向岛的sector覆盖数；
2. 最大化多个覆盖情景中的最差可设计表面积覆盖；
3. 最小化最大未覆盖连通块；
4. 最大化高倾向面积覆盖；
5. 最小化探针间覆盖重叠和潜在交叉干扰；
6. 最小化组合稳定性风险；
7. 最小化突变总数和高阶外推数量。

应使用多阶段MILP或等价的严格优先级方法，避免一个加权总分因量纲变化改变科学优先级。

### 15.3 输出层级

每个sector输出：

- `PRIMARY`：主面板候选；
- `CONSERVATIVE`：更少位点或更低稳定性风险；
- `AGGRESSIVE`：覆盖或扰动更强；
- `DECONVOLUTION_SINGLES`：用于阳性探针后续拆解的单点构建；
- `NO_FEASIBLE_CANDIDATE`：无安全组合时明确报告。

## 16. 运行模式与失败降级

### 16.1 strict模式

正式推荐必须满足：

- ProCaliper完成或SASA一致性已有批准的本地替代；
- SPURS single完成；
- DiscoTope完成；
- SPURS multi完成；
- 最终主候选完成FoldX，或用户显式允许SPURS-only。

否则结构分析可继续，但 `recommendation_status = BLOCKED_EXTERNAL_EVIDENCE`。

### 16.2 preview模式

允许使用现有代理和旧表位先验生成候选，用于UI预览和调试，但：

- 不得标为正式推荐；
- 不得填充伪造的kcal/mol；
- 顶层显示缺失工具；
- 导出文件名和schema仍与strict一致。

### 16.3 稳定错误码

至少定义：

```text
TOOL_NOT_CONFIGURED
TOOL_ENVIRONMENT_INVALID
MODEL_WEIGHTS_MISSING
MODEL_REVISION_MISMATCH
TOOL_TIMEOUT
TOOL_CRASHED
OUTPUT_MISSING
OUTPUT_SCHEMA_INVALID
RESIDUE_COUNT_MISMATCH
WT_RESIDUE_MISMATCH
CHAIN_MAPPING_AMBIGUOUS
UNSUPPORTED_STRUCTURE
LICENSE_RESTRICTED
NO_SURFACE_CANDIDATES
NO_STABLE_SINGLE_OPTIONS
NO_FEASIBLE_COMBINATIONS
NO_CONSENSUS_STABLE_COMBINATION
```

## 17. 配置结构

新增 `JOINT_V4_ML_PARAMETERS`，禁止散落magic number：

```python
JOINT_V4_ML_PARAMETERS = {
    "schema_version": "4.0",
    "mode": "strict",
    "max_sectors": 7,
    "max_primary_constructs": 7,
    "surface": {...},
    "procaliper": {...},
    "spurs_single": {...},
    "discotope": {...},
    "sectoring": {...},
    "candidate_search": {...},
    "spurs_multi": {...},
    "foldx_confirmation": {...},
    "panel_optimizer": {...},
    "runtime": {...},
}
```

每次运行将解析后的完整配置、稳定排序后的JSON和SHA-256写入manifest。

## 18. 输出合同

### 18.1 新增文件

```text
tool_inventory.json
structure_contexts.json
residue_mapping.csv
surface_evidence.csv
sasa_comparison.csv
spurs_single_matrix.csv.gz
single_mutation_options.csv
discotope_residues.csv
epitope_islands.json
probe_sectors_v4.json
candidate_combinations.csv.gz
combination_stability.csv
panel_candidates.json
selected_panel_v4.json
v4_manifest.json
v4_scientific_report.md
external_tools/*.log
```

### 18.2 顶层结果必须显示

- `recommendation_status`；
- 实际sector数和主构建数；
- 各工具运行状态；
- 结构与模型版本；
- 未覆盖表面和原因；
- 模型外推候选数量；
- SPURS/FoldX分歧；
- 是否存在低表位倾向覆盖区；
- 所有待实验验证声明。

### 18.3 兼容投影

继续生成现有UI可消费的 `regions.json` 和 `designs.json`，但附加：

```text
source_engine = joint_v4_ml
source_schema_version = 4.0
semantic_projection = true
```

不得让旧UI把 `ProbeSector` 展示成已预测的真实抗体footprint。

## 19. API、CLI与UI变更

### 19.1 CLI

```powershell
ecd-epitope-scanner analyze `
  --structure input.pdb `
  --sequence sequence.fasta `
  --design-engine joint_v4_ml `
  --ml-mode strict `
  --max-sectors 7 `
  --max-primary-constructs 7 `
  --output results/run
```

外部环境路径通过专用参数或环境变量配置，不写死到源码。

### 19.2 API

请求增加：

- `design_engine = joint_v4_ml`
- `ml_mode`
- `max_sectors`
- `max_primary_constructs`
- `allow_spurs_only_final`
- `biological_context_mode`

作业状态增加明确阶段：

```text
PROCALIPER_ANALYSIS
SPURS_SINGLE
DISCOTOPE
EPITOPE_AWARE_PARTITIONING
COMBINATION_SEARCH
SPURS_MULTI
FOLDX_CONFIRMATION
PANEL_OPTIMIZATION
```

### 19.3 UI

UI至少能切换显示：

- SASA；
- DiscoTope倾向连续色带；
- 高倾向连通岛；
- sector边界；
- 主/备选突变；
- 未覆盖表面；
- 模型分歧和外推警告。

界面中“安全”统一改为“低预测稳定性风险”；“表位”在预测层统一写为“表位倾向”。

## 20. 测试策略

### 20.1 Adapter合同测试

- 工具未安装；
- 错误环境路径；
- 超时和崩溃；
- stdout污染；
- 输出字段缺失；
- WT残基不一致；
- 链长度不一致；
- 插入码和不连续编号；
- 缓存命中与模型版本变化失效。

测试使用固定小型fixture和伪worker，不要求CI下载大型模型。

### 20.2 科学逻辑单元测试

- 高DiscoTope分数不能绕过硬保护；
- 低分残基仍可作为连接桥或覆盖候选；
- ProCaliper欧氏近邻不能创建无表面路径的边；
- 相邻高倾向seed在满足直径限制时不被拆开；
- 无高倾向seed时可产生 `COVERAGE_ONLY` sector；
- 每个sector始终连通；
- 每种突变阶数有独立候选配额；
- 三/四点候选必须带外推等级；
- SPURS和FoldX分歧被正确标记；
- 外部工具失败不填0值。

### 20.3 优化器测试

- 小图与穷举oracle一致；
- 构建数严格不超过预算；
- 增加预算不降低同目标下最优覆盖；
- 同一sector不选择多个PRIMARY；
- 高阶组合数量上限有效；
- 候选输入顺序不影响结果；
- MILP超时后fallback确定性一致。

### 20.4 集成测试

- 仅单链PDB；
- AlphaFold结构；
- 多链生物学装配体；
- 强糖基化ECD；
- β-sandwich相反表面；
- 低置信度长末端；
- 无可行二点组合；
- strict和preview模式；
- 完整API到离线报告链路。

## 21. 科学基准与消融实验

### 21.1 数据分层

- 开发集：用于参数和阈值调整。
- 冻结测试集：仅用于最终比较。
- 前瞻实验集：本项目真实抗体和抗原数据。

公开复合物benchmark只向算法提供抗原结构，真实抗体接触残基只用于事后评估。

### 21.2 对照方法

- 当前 `surface_v2`；
- 当前 `joint_v3`；
- 随机安全表面位点；
- 均匀测地采样；
- V4不含DiscoTope；
- V4不含SPURS；
- V4仅SPURS、不做FoldX复核；
- 完整V4-ML。

### 21.3 主要指标

- 在最多七个主构建下，真实抗体footprint被至少一个突变命中的比例；
- 被至少两个位点命中的比例；
- 高倾向岛覆盖；
- 可设计表面积覆盖；
- 最大未覆盖连通块；
- 主构建间覆盖冗余；
- 候选的SPURS/FoldX风险分布；
- 运行时间和峰值内存；
- 失败/降级率。

### 21.4 前瞻实验验收

至少记录：

- Myc或等价表面表达归一化；
- 一个以上构象/折叠控制；
- WT和无关抗原；
- 主组合；
- 阳性组合的单点去卷积；
- 必要时天然配体或功能读出。

算法成功标准应是“在相同构建预算下提高表位粗定位信息量，同时不增加普遍表达/折叠失败”，而不是单纯提高预测分数。

## 22. 性能预算

Phase 0先在目标Windows硬件上实测，再冻结硬门槛。建议初始目标：

- 几何、ProCaliper和分区：典型500 aa目标CPU下少于2分钟；
- 单点SPURS：GPU单目标少于1分钟，模型加载可缓存；
- DiscoTope：GPU单目标少于1分钟，模型加载可缓存；
- 组合SPURS：短名单少于2分钟；
- FoldX复核：按每区最多6个候选限制总成本；
- 总流程：不含首次模型下载时，典型目标在10分钟级完成。

以上为工程预算，不是已验证性能承诺。CPU-only模式必须报告预计显著变慢。

## 23. 安全、许可与复现

- 不自动将未公开结构上传远程服务。
- 模型权重首次下载必须显式记录来源、revision和SHA-256。
- 不把CC BY-NC 4.0的DiscoTope代码或权重直接并入本项目发行包，除非用途和再分发条件已经确认。
- SPURS、ProCaliper、DiscoTope和FoldX分别记录许可证与使用限制。
- 任何模型自动更新都禁止；升级必须改变manifest中的版本和缓存键。
- 所有外部命令使用参数列表启动，禁止拼接shell字符串。
- 临时目录按单作业隔离，输出文件名不能由未清洗的用户输入直接构造。

## 24. 分阶段代码生成计划

### Phase 0：冻结V3和性能基线

工作：

- 保存V3代表性输入、输出哈希和运行时间；
- 建立大型模型之外的快速CI测试集；
- 建立V2/V3 shadow comparison基线。

完成条件：当前工作树相关测试通过；同一输入产生稳定的非时间字段哈希。

### Phase 1：V4 schema和通用adapter框架

工作：

- 新建 `joint_design_v4/models.py`；
- 实现tool inventory、JSON子进程合同、超时和错误码；
- 增加伪worker合同测试。

完成条件：不安装任何外部模型也可完整测试adapter生命周期。

### Phase 2：结构上下文和映射

工作：

- 生成规范化单链PDB和映射sidecar；
- 区分prediction/biological context；
- 实现严格WT和编号检查。

完成条件：插入码、缺失残基、多链和AlphaFold fixture全部通过。

### Phase 3：ProCaliper接入

工作：

- 实现inventory和worker；
- 导入SASA/距离特征；
- 实现双SASA一致性报告；
- 保持显式表面图权威关系。

完成条件：ProCaliper不能改变无表面路径残基的连通性。

### Phase 4：SPURS single接入

工作：

- 实现全矩阵推理与缓存；
- 验证符号方向；
- 建立稳定性分带和替换菜单压缩；
- 输出全部拒绝原因。

完成条件：固定fixture上的矩阵、映射、阈值边界和缓存测试通过。

### Phase 5：DiscoTope接入

工作：

- 实现完整链推理；
- 导入原始/校准分数；
- 计算表面百分位、局部峰和高倾向岛；
- 实现许可状态检查。

完成条件：DiscoTope分数只影响排序和软目标，不改变硬保护结果。

### Phase 6：V4证据融合与分区

工作：

- 合并表面、保护、SPURS和DiscoTope证据；
- 实现must-link/cannot-link；
- 枚举K=1..7并使用字典序目标选区；
- 输出低倾向区和未覆盖组件。

完成条件：所有sector连通；相反表面不误并；相邻热点在可行时不拆分。

### Phase 7：候选与SPURS multi

工作：

- 复用V3 beam search并加入阶数配额；
- 批量运行SPURS multi；
- 计算epistasis估计和外推等级；
- 保留每区多样化短名单。

完成条件：候选量受控；2/3/4点均有公平机会；高阶外推不可隐藏。

### Phase 8：FoldX共识与面板优化

工作：

- 对每区短名单运行FoldX；
- 生成共识状态；
- 扩展MILP加入sector、外推和干扰约束；
- 输出PRIMARY/CONSERVATIVE/AGGRESSIVE。

完成条件：最终面板不超过预算；strict模式不选择缺少所需证据的候选。

### Phase 9：pipeline/API/CLI/report/UI集成

工作：

- 新增引擎选择和运行阶段；
- 输出V4原生文件与旧格式投影；
- 在UI显示倾向、sector、稳定性证据和警告；
- 离线导出保持可用。

完成条件：旧任务不受影响；V4任务在无模型时给出清晰配置指导。

### Phase 10：benchmark和默认策略

工作：

- 运行公开复合物benchmark和消融；
- 冻结阈值与权重；
- 运行前瞻小规模实验；
- 决定是否替换默认引擎。

完成条件：预注册主要指标优于简单基线，且稳定性/表达失败率没有不可接受增加。

## 25. 每个代码生成任务的固定提示合同

后续每次生成代码只处理一个Phase内的一个可独立验证任务，并必须包含：

1. 允许修改的文件列表；
2. 不得改变的公共接口；
3. 输入/输出类型；
4. 错误码与失败模式；
5. 确定性和缓存要求；
6. 新增配置项；
7. 单元测试和集成测试；
8. 验收命令；
9. 科学声明边界；
10. 与V2/V3兼容要求。

禁止一次性提示模型“实现整个V4”，否则外部依赖、schema和科学语义会难以审阅。

## 26. 建议的首批代码任务

用户确认D1–D10后，按以下顺序开始：

1. `V4-001`：新增V4数据模型、配置骨架和schema序列化测试。
2. `V4-002`：新增外部工具inventory与通用JSON worker runner。
3. `V4-003`：新增规范化单链PDB与残基映射sidecar。
4. `V4-004`：ProCaliper adapter及SASA一致性报告。
5. `V4-005`：SPURS single adapter、缓存与fixture解析。
6. `V4-006`：单点证据融合和每位点最多3种替换菜单。
7. `V4-007`：DiscoTope adapter及高倾向岛派生。
8. `V4-008`：表位倾向感知sector builder。
9. `V4-009`：带阶数配额的V4 candidate search。
10. `V4-010`：SPURS multi和组合不确定性。
11. `V4-011`：FoldX共识与V4 panel optimizer。
12. `V4-012`：pipeline/API/CLI/export接入。
13. `V4-013`：UI与离线报告。
14. `V4-014`：benchmark harness与消融报告。

## 27. 开始编码前的Definition of Ready

至少应确认：

- D1、D2、D4、D5、D6、D9、D10；
- 最终实验主构建预算；
- 是否允许非商业许可依赖；
- 目标运行机器是否有NVIDIA GPU及可用显存；
- SPURS/DiscoTope模型是否允许首次联网下载，还是必须离线安装；
- 最终推荐是否必须有FoldX结果；
- 现有未提交V3代码是否作为V4基线保留。

若这些问题尚未确认，可以先完成 `V4-001` 至 `V4-003`，但不应冻结模型门槛、部署脚本或最终面板政策。

## 28. 依据与已知限制

- ProCaliper使用Shrake–Rupley计算SASA，并提供残基欧氏距离网络；这适合作为结构特征与审计层，不足以替代表面测地连通图。
- DiscoTope-3.0输出逐残基表位倾向并包含长度/表面校准，但其外部测试AUC-PR提示它更适合排序而非硬筛选。
- SPURS single可一次预测全部单点替换；multi模型显式建模非加和效应，但训练和主要公开验证集中在双突变。
- DiscoTope-3.0源码许可证为CC BY-NC 4.0，商业用途必须单独评估。
- 所有初始阈值均须在开发集和前瞻实验中校准；在此之前只能作为版本化启发式配置。

参考：

1. ProCaliper: https://pmc.ncbi.nlm.nih.gov/articles/PMC12607263/
2. SPURS: https://www.nature.com/articles/s41467-025-67609-4
3. SPURS repository: https://github.com/luo-group/SPURS
4. DiscoTope-3.0: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1322712/full
5. DiscoTope-3.0 repository and license: https://github.com/Magnushhoie/DiscoTope-3.0

## 29. A方案本机环境配置记录（RTX 3060 Laptop）

本节记录当前机器上已经完成的依赖与工具安装，作为后续代码生成和复现实验的运行基线。模型权重和缓存文件放在项目的 `models/` 目录，已加入 `.gitignore`，不得提交到版本库。

### 29.1 机器与运行层

- GPU：NVIDIA GeForce RTX 3060 Laptop GPU，6 GB显存；Windows驱动可用，WSL2 Ubuntu-22.04中 `nvidia-smi` 可见同一GPU。
- 推荐运行方式：Windows负责项目文件、源码和下载；WSL2 Ubuntu负责SPURS/DiscoTope的Linux推理环境。ProCaliper使用Windows Conda环境即可。
- GPU推理必须设置显存上限和批次；当输入结构或模型峰值超过6 GB时自动回退CPU，并在结果中记录 `device_fallback`。

### 29.2 已安装环境

| 用途 | 环境/路径 | 版本与状态 |
|---|---|---|
| 项目现有运行环境 | `.venv\\Scripts\\python.exe` | Python 3.12；保留现有V2/V3依赖，不混装深度学习依赖 |
| ProCaliper | `G:\\anaconda\\envs\\ecd_procaliper` | Python 3.11；`procaliper 0.2.1`，含 BioPython 1.88、OpenBabel、SciPy 等依赖 |
| SPURS single/multi | WSL `/home/segayi1/miniconda3/envs/spurs_v4` | Python 3.11.16；`torch 2.5.1+cu124`；CUDA可用；源码为 `vendor/spurs` 的 `fix/py311-inference-beta` 分支 |
| DiscoTope-3.0 | WSL `/home/segayi1/miniconda3/envs/discotope_v4` | Python 3.11；`torch 2.5.1+cu124`、PyG 2.7、BioPython 1.86、BioTite 1.6、XGBoost 3.1.2、pyGAM 0.12、scikit-learn 1.7；7c4s样例已通过 |

源码位置：

- SPURS：`vendor/spurs/`
- DiscoTope-3.0：`vendor/discotope3/`
- SPURS模型目录：`models/spurs/`（包含 `spurs/` 和 `spurs_multi/`）

### 29.3 模型权重与离线部署

- SPURS single：`models/spurs/spurs/checkpoints/best.ckpt`，官方仓库文件大小约2.68 GB。
- SPURS multi：`models/spurs/spurs_multi/checkpoints/best.ckpt`，官方仓库文件大小约2.70 GB。
- SPURS single SHA-256：`9A1FD0E0B123DC74CC394591034288348594AA9B2F4D152AD808A24D4CD5C642`。
- SPURS multi SHA-256：`018081B900A77958F67355F0D323DFDC3E18D63F61C8764F0B268E7C1E32080C`。
- SPURS运行时依赖的ESM-2：`models/spurs/esm2_t33_650M_UR50D.pt`，SHA-256：`EA9D0522B335A8778DEA6535A65301F10208DECE28CD5865482B0B1FC446168C`。
- 两个模型均通过SPURS提供的本地接口加载：`spurs.inference.get_SPURS(<模型目录>)` 和 `get_SPURS_multi_from_hub` 对应的本地化适配器；实现时不得在每次任务运行时联网下载。
- DiscoTope模型文件来自源码仓库的 `models.zip`（解压到 `vendor/discotope3/models/`）；ESM-IF1权重如首次运行需要下载，应预先放入离线缓存并记录校验值。
- DiscoTope ESM-IF1：`models/discotope/esm_if1_gvp4_t16_142M_UR50.pt`，SHA-256：`BE4BA36EDEC22A9BFAA4946FF6B2815F1F19D8A3D7E0EADA8B796D5A0EAE9FD4`。
- 权重下载完成后执行 `Get-FileHash -Algorithm SHA256`，将结果填入运行清单；权重本身不提交Git。

### 29.4 安装/验证命令

Windows PowerShell：

```powershell
conda run -n ecd_procaliper python -c "import procaliper, Bio; print(procaliper.__version__)"
Get-FileHash .\models\spurs\spurs\checkpoints\best.ckpt -Algorithm SHA256
Get-FileHash .\models\spurs\spurs_multi\checkpoints\best.ckpt -Algorithm SHA256
```

WSL2 Ubuntu：

```bash
source /home/segayi1/miniconda3/etc/profile.d/conda.sh
conda activate spurs_v4
python -c "import torch, spurs; print(torch.__version__, torch.cuda.is_available()); print(torch.cuda.get_device_name(0))"
conda activate discotope_v4
python -c "import torch, torch_geometric, Bio, biotite, xgboost; print(torch.__version__, torch.cuda.is_available())"
```

已完成验收：环境导入、GPU可见性、SPURS single本地checkpoint加载、SPURS multi本地checkpoint加载、SPURS single在RTX 3060上加载（约2.53 GiB显存）、DiscoTope 7c4s六链样例PDB推理均通过。项目adapter尚未生成，故在adapter完成前工具状态记为 `installed_verified_external`，而不是项目级 `available`。

### 29.5 兼容性和许可备注

- SPURS和DiscoTope使用隔离环境，避免PyTorch、BioPython和NumPy版本互相污染现有项目环境。
- DiscoTope-3.0源码为CC BY-NC 4.0；本项目默认仅用于非商业研究。若项目用途变化，必须在发布前重新审核许可。
- 当前DiscoTope官方依赖声明为PyTorch 2.10系列，而本机RTX 3060验证环境采用可用的2.5.1+cu124组合；因此在代码中记录实际版本，并把端到端样例结果作为兼容性门槛，而不是仅依据版本字符串放行。
