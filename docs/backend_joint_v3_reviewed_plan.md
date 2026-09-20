# ECD Epitope Scanner 后端联合设计 V3：科学审阅与实施规格

> 文档状态：Reviewed implementation plan
> 目标版本：`joint_v3`
> 适用实验：Jurkat 表达膜锚定 Fab；K562 表达 `N-Myc–antigen ECD–PDGFRB TM` 野生型及突变体；细胞混合后以结合/报告信号进行抗体表位粗定位
> 核心输入：抗原结构文件（PDB/mmCIF）；序列、链、膜锚点等信息尽量从结构与构建设定推导，无法可靠推导时必须显式报告不确定性

## 1. 审阅结论

项目的实验目标是可实现的，但算法不应声称“预测精确表位”。它应解决一个更窄、更可验证的问题：在最多约 6 个抗原组合突变探针的预算下，尽量均匀地扰动整个可设计 ECD 表面，使 10 条抗体的实验响应能够把结合区域粗分到若干空间区域。

当前“表面残基 → 折叠稳定性 → footprint → 分区 → 组合突变”的方向没有根本错误，问题在于若按单向硬筛选实现，会把多个未经验证的先验串联起来：前级误判会永久删除候选，最终得分看似精确，实际却无法解释。因此 V3 应改成“结构审计 + 多尺度表面表示 + 候选探针生成 + 联合面板优化 + 不确定性报告”，稳定性、可及性、空间覆盖和实验可解释性应联合决策，而不是逐级一票否决。

### 1.1 上一版方案中需要修正的地方

| 原设想 | 科学或工程问题 | V3 修正 |
|---|---|---|
| 用一个固定 footprint 大小代表抗体接触面 | 抗体接触面积是分布，不是常数；不同抗体和结构状态差异明显 | 使用多个可配置面积尺度做敏感性分析；不得把任一尺度称为真实表位 |
| 再创建一层任意 `SurfaceUnit` 图分区 | 与 footprint、region 的含义重叠，会引入新的分区先验 | 优化基本单元直接采用“可设计的表面残基”，以残基映射表面积作为权重 |
| 每个 surface unit 的最小覆盖最大化 | 只要存在一个不可突变或不可达单元，最小值就会退化为 0 | 改用分层目标：先最大化可设计表面积覆盖，再最小化最大未覆盖连通块，然后控制冗余和风险 |
| 指数距离衰减直接表示突变影响 | 这是启发式核函数，不是结合概率，也不是物理能量 | 同时计算位点覆盖、局部测地覆盖、footprint 拦截三类指标，并跨参数情景报告稳健性 |
| 枚举所有 2–4 位点组合 | 候选数呈组合爆炸，无法稳定运行 | 每区预筛位点与替换选项，再用确定性 beam search；高成本能量计算只用于短名单 |
| 无外部稳定性引擎就拒绝候选 | 会令默认流程不可用，也会掩盖“缺证据”和“高风险”的区别 | 允许规则代理结果，但标记 `PROXY_ONLY` / `REVIEW_REQUIRED`，绝不伪造 kcal/mol |
| 每个组合必须含一个“锚点突变” | 没有实验或结构依据，是人为先验 | 删除强制锚点规则；仅以空间分散、化学扰动与稳定性共同评分 |
| 从孤立 ECD PDB 推出唯一膜平面 | 质粒结构只给出连接拓扑，不能唯一决定游离 ECD 的膜相对姿态 | 默认仅输出 C 端连接点、相对轴和方向不确定性；只有存在实际 TM 坐标时才启用膜几何模式 |
| Myc 阳性等同正确折叠 | Myc 主要报告表面展示/表达，不能证明天然折叠 | 输出实验解释要求：Myc 为表达归一化，折叠需额外构象性对照或多抗体一致性证据 |
| 固定 95% 为科学合格线 | 未经基准集校准的阈值只是人为规定 | 总是输出实际覆盖率；放行阈值需在独立基准集上冻结后才可作为 QC |

## 2. 系统边界与明确非目标

### 2.1 V3 必须完成

1. 仅导入抗原 PDB/mmCIF 即可运行默认流程。
2. 识别结构链、缺失、非标准残基、重复链、结构置信度与 C 端锚定信息。
3. 基于真实三维表面生成多尺度假想抗体接触窗口。
4. 生成 2–4 个突变位点组成的局部组合探针候选。
5. 在探针数上限内联合选择面板，默认上限为 6，但作为配置而非科学常数。
6. 输出覆盖、不确定性、候选淘汰原因、稳定性证据来源与可复现配置。
7. 兼容现有 Surface V2 结果，并支持 shadow mode 对比。

### 2.2 V3 不应声称或默认执行

- 不从抗原单体结构预测某条具体抗体的精确表位。
- 不把一个粗定位探针区称为一个真实抗体 footprint。
- 不从质粒拓扑或孤立 ECD 结构计算唯一膜平面。
- 不把 FoldX/EvoEF 分数当作细胞表面正确折叠的证明。
- 不把 N 端 Myc 信号当作正确折叠证明。
- 不在 MVP 中模拟完整 linker、糖链、跨膜区或细胞间受力过程。
- 不在没有实验数据时输出“结合概率”。

## 3. 术语与数据语义

为防止代码和 UI 混淆，V3 固定以下术语：

| 术语 | 定义 | 是否重叠 |
|---|---|---|
| `SurfaceResidueUnit` | 一个满足表面条件且可映射的残基，是覆盖优化的最小单元；权重为该残基表面积或近似表面积 | 否 |
| `MicroFootprint` | 从表面中心向外按测地距离生长得到的抗体尺度假想窗口，用于敏感性分析 | 是 |
| `ProbeSector` | 一个候选探针意图扰动的连续表面区域，可由相邻多个 microfootprint 的并集构成 | 区内连续 |
| `MutationOption` | 某个残基的一种允许替换及其风险/证据 | 不适用 |
| `ProbeCandidate` | 2–4 个突变组成的具体构建体候选 | 不适用 |
| `ProbePanel` | 最终选中的不超过预算的多个构建体 | 不适用 |

关键限制：`MicroFootprint` 是用于模拟可能抗体接触范围的分析窗口；`ProbeSector` 是为节省实验数量而扩大的粗定位区域。两者不得在 API、报告或 UI 中共用同一个字段名。

## 4. 目标架构

```text
PDB/mmCIF
  │
  ├─ StructureAudit ── chain/model/residue mapping/quality warnings
  │
  ├─ TopologyContext ── C-terminal anchor + relative axis + uncertainty
  │
  ├─ SurfaceModel ── SASA/rSASA + surface graph + geodesic distances
  │       │
  │       ├─ SurfaceResidueUnits
  │       └─ Multi-scale MicroFootprints
  │
  ├─ ProbeSectorBuilder ── connected coarse regions
  │
  ├─ CandidateGenerator ── mutation options + deterministic beam search
  │       │
  │       ├─ cheap structural/rule proxy
  │       └─ optional combination-level FoldX/EvoEF shortlisting
  │
  ├─ PanelOptimizer ── MILP lexicographic objectives
  │       └─ deterministic greedy/local-swap fallback
  │
  └─ ScientificReport ── selected panel, coverage scenarios, provenance,
                          warnings, alternative panels and machine-readable files
```

## 5. 模块与文件修改方案

尽量复用当前模型与流程，避免一次性重写。

### 5.1 新增模块

```text
ecd_epitope_scanner/
  joint_design/
    __init__.py
    models.py                 # V3 数据结构与 schema version
    structure_audit.py        # 链/模型/缺失/编号/重复链检查
    topology.py               # C 端锚点、相对轴、方向不确定性
    surface_units.py          # 表面残基单元与权重
    microfootprints.py        # 多尺度测地 footprint
    sectors.py                # 连通粗定位区生成
    mutation_options.py       # 替换规则和局部风险特征
    candidate_search.py       # 确定性 beam search
    stability.py              # 代理与外部引擎适配器
    coverage.py               # 三类覆盖矩阵及敏感性情景
    optimizer.py              # MILP 和 fallback
    reporting.py              # JSON/CSV/Markdown 输出
```

### 5.2 修改现有模块

- `pipeline.py`
  - 为 `PipelineRequest` 增加 `design_engine: Literal["surface_v2", "joint_v3"]`。
  - `joint_v3` 走新入口；`surface_v2` 保留，作为回归与 shadow comparison。
  - 默认配置切换应在 V3 通过基准与验收后单独提交，不能与算法首次落地同时切换。
- `config.py`
  - 新增 `JOINT_V3_PARAMETERS`，所有尺度、阈值、预算、beam width 和 solver timeout 均配置化。
  - 配置写入结果 manifest，并计算稳定哈希。
- `models.py`
  - 保留旧 `RegionRecord`；增加 V3 adapter，避免旧 UI/导出立即失效。
- `exporting.py`
  - 增加 V3 JSON、CSV、Markdown 科学报告；旧格式继续输出兼容视图。
- `design.py`
  - 旧 `design_regions` 不直接承载 V3 联合优化逻辑。
  - 删除或绕过 V3 中的强制 anchor mutation 行为。
- `surface_patch/*`
  - 复用表面点云、表面图和测地距离实现。
  - 不在原 V2 类型中硬塞 `ProbeSector` 新语义。

## 6. 数据模型

以下字段是实施所需的最小合同，具体可用 `dataclass` 或 Pydantic 实现。

```python
@dataclass(frozen=True)
class ResidueId:
    model_id: str
    chain_id: str
    auth_seq_id: int
    insertion_code: str
    one_letter: str

@dataclass
class StructureAudit:
    selected_model: str
    selected_chains: list[str]
    chain_decision: str
    residue_mapping_status: str
    missing_segments: list[dict]
    nonstandard_residues: list[dict]
    duplicate_chain_groups: list[list[str]]
    quality_warnings: list[str]
    fatal_errors: list[str]

@dataclass
class TopologyContext:
    mode: Literal["topology_only", "tm_geometry"]
    anchor_residue: ResidueId | None
    anchor_point_xyz: tuple[float, float, float] | None
    relative_axis_xyz: tuple[float, float, float] | None
    membrane_plane: dict | None
    uncertainty: Literal["low", "medium", "high", "not_estimable"]
    assumptions: list[str]

@dataclass
class SurfaceResidueUnit:
    residue: ResidueId
    rasa: float
    sasa_a2: float | None
    designable: bool
    exclusion_reasons: list[str]
    component_id: int

@dataclass
class MicroFootprint:
    footprint_id: str
    center: ResidueId
    scale_id: str
    target_area_a2: float
    residue_ids: list[ResidueId]
    realized_area_a2: float | None

@dataclass
class MutationOption:
    residue: ResidueId
    mutant_aa: str
    perturbation_class: str
    rule_risk: float
    stability_value: float | None
    stability_unit: str | None
    evidence_level: Literal["EXTERNAL_ENGINE", "PROXY_ONLY"]
    rejection_reasons: list[str]

@dataclass
class ProbeCandidate:
    candidate_id: str
    sector_id: str
    mutations: list[MutationOption]
    connected: bool
    coverage_by_scenario: dict[str, dict[str, float]]
    structural_risk: float
    evidence_level: str
    warnings: list[str]

@dataclass
class ProbePanel:
    panel_id: str
    selected_candidate_ids: list[str]
    objective_vector: list[float]
    achieved_coverage: dict[str, float]
    uncovered_components: list[dict]
    solver_status: str
    optimality_gap: float | None
    alternative_panel_ids: list[str]
```

所有残基必须由完整结构标识映射，不能只用 `K123` 之类字符串作为内部主键。展示层可再格式化。

## 7. 算法详细规格

### 7.1 阶段 A：结构审计与链选择

#### 输入检查

1. 解析所有 model、polymer chain、标准与非标准残基、原子占有率、altloc、缺失坐标。
2. 默认使用第一个模型，但必须在报告中说明；多模型/NMR ensemble 可在后续版本扩展。
3. 生成 PDB/mmCIF 编号到连续序列编号的双向映射。
4. 记录缺失片段、未解析 C 端、断链和只有 Cα 的低信息区域。

#### 链决策

1. 首先识别高序列同一性且高结构相似的重复链，视为可能的同源拷贝。
2. 若多个链是等价重复拷贝，选择一个代表链运行，同时保留 oligomer context 用于遮挡/SASA 计算的可选情景。
3. 若存在多个非等价、均可能是目标抗原的链，不以“长度差 10%”之类任意规则静默选择；默认阻断自动设计并返回 `CHAIN_SELECTION_REQUIRED`。
4. 水、离子、常见结晶添加剂和明显短肽不作为抗原主链候选。

#### 结构置信度

- 实验结构：记录分辨率、缺失和局部 B factor（若存在），不把 B factor 单独解释为折叠置信度。
- AlphaFold 类结构：若 pLDDT 存放于 B-factor 字段，可作为局部置信度；必须先确认文件来源/字段语义。
- 未知来源：质量字段标记为 `UNKNOWN_PROVENANCE`，不能假设 B-factor 是 pLDDT。

### 7.2 阶段 B：膜拓扑与方向

默认 `topology_only`：

1. 由构建体知识固定抗原链方向为 `N-Myc–ECD–C-terminal linker/TM`。
2. 在结构中寻找最后一个有可靠坐标的 C 端 ECD 残基，记录为 anchor residue。
3. 相对轴可由 C 端局部主链方向与 ECD 质心构造，仅用于“靠近/远离连接端”的描述。
4. 未解析尾部长度、linker 柔性和结构来源均进入 uncertainty。
5. 不生成唯一膜平面，不以相对轴直接删除候选表面。

可选 `tm_geometry`：只有输入结构真实包含可识别的跨膜螺旋坐标时，才可调用 PPM/OPM 类方法估计膜位置。该模式的结果仍须带方法、版本和输入链 provenance。

### 7.3 阶段 C：表面残基与表面图

1. 计算单体情景和（若存在生物组装信息）组装体情景下的 SASA/rSASA。
2. 默认表面入选规则可延续当前 `rSASA >= 0.25` 且绝对 SASA 达到最低值，但它只是初始配置。
3. 对临界残基保留连续权重，避免阈值附近的硬跳变：

```text
surface_weight = mapped_area × accessibility_confidence
```

4. 构建表面残基图：节点为表面残基；边由表面点云连通性或可靠的近邻关系定义；边权为表面测地距离。
5. 分离互不连通的表面组件，记录小组件和缺失结构导致的断裂。

表面残基的可设计性与表面性分开：一个残基可以是表面残基但因二硫键、糖基化 sequon、配体接触、催化/功能位点、极端保守性或局部结构风险而不可设计。不可设计不等于不存在，仍参与“总表面”和“不可覆盖表面”的报告。

### 7.4 阶段 D：多尺度 MicroFootprint

以候选中心沿表面图按测地距离生长窗口，直到达到目标表面积。初始建议尺度为约 600、800、1000 Å²，仅作为文献支持的合理起点和敏感性区间，而非固定真值。

必须满足：

- 每个 footprint 记录目标面积、实际面积、残基数和边界截断原因。
- 对缺失/低置信区域相交的 footprint 降低置信度，不自动补齐虚构表面。
- 去重只能减少数值冗余，不能删除独特表面；去重阈值配置化并输出。
- footprint 中心采样使用确定性 farthest-point 或等价方法，随机算法必须固定 seed。

### 7.5 阶段 E：ProbeSector 构建

ProbeSector 的任务是允许少量探针进行更粗的定位。

1. 将空间相邻且高度重叠的 microfootprint 合并成连通 sector。
2. sector 必须是单个表面连通面，禁止跨越蛋白两侧或以欧氏距离“穿过蛋白内部”合并。
3. 合并目标兼顾：面积均衡、形状紧凑、可设计位点充足和全局覆盖。
4. sector 数量可多于最终探针预算；最终由 panel optimizer 选候选，而不是先强制切成正好 6 区。
5. 输出 sector 面积、直径、包含的 footprint、可设计面积比例和形状紧凑度。

不要把 sector 面积直接解释成抗体接触面积。UI 建议显示“粗定位覆盖区”。

### 7.6 阶段 F：突变选项生成

#### 先验保护与默认排除

默认不修改：

- 二硫键 Cys；
- 明确糖基化位点或可能破坏 N-X-S/T sequon 的位点；
- 金属/配体关键配位位点；
- 已知催化或必需功能位点（若有可靠注释）；
- 埋藏核心、主链氢键关键位置、明显盐桥核心；
- 结构缺失或坐标质量不足、无法进行风险计算的位置。

这些规则应允许高级用户覆盖，但覆盖操作必须进入 manifest。

#### 替换策略

不使用“所有位点统一丙氨酸”作为唯一规则。每个候选位点生成少量、具有不同扰动机制的替换选项，例如：

- 去除电荷；
- 电荷反转（风险通常更高）；
- 缩短侧链；
- 芳香/疏水侧链去除；
- 在可接受位置引入空间阻碍。

替换策略必须避免引入新的糖基化 sequon、游离 Cys、明显蛋白酶位点或不可解释的多重理化改变。具体替换表放入版本化配置并用基准集调整，而不是散落在代码条件分支中。

#### 评分分离

至少拆分以下字段，禁止先合成一个“科学总分”后丢失来源：

- `surface_accessibility`
- `local_structural_risk`
- `functional_annotation_risk`
- `chemical_perturbation_strength`
- `spatial_novelty`
- `stability_evidence`

### 7.7 阶段 G：候选组合搜索

不得完整枚举全局 2–4 位点组合。采用两级、确定性搜索：

1. 每个 sector 先按硬排除规则过滤。
2. 依据可及性、空间代表性、化学扰动互补和低结构风险选出有限位点池；池大小配置化。
3. 每个位点仅保留有限替换选项；数量配置化。
4. beam search 从 2 位点扩展到 3、4 位点。
5. 每一步先用便宜代理评分，保留固定 beam width。
6. 对候选使用稳定排序键，确保相同输入产生完全相同输出：

```text
(objective tuple, chain_id, residue_number, insertion_code, mutant_aa)
```

7. 对短名单执行可选外部稳定性计算。
8. 若使用 FoldX/EvoEF，应计算完整组合而不是简单相加单点 ΔΔG；单点相加只能作为明确标注的代理。

建议伪代码：

```python
for sector in sectors:
    site_pool = rank_and_cap_sites(sector, config.site_pool_cap)
    option_pool = build_options(site_pool, config.options_per_site_cap)
    beam = initialize_pairs(option_pool)
    for size in range(2, config.max_mutations_per_probe + 1):
        beam = expand_without_duplicate_sites(beam, option_pool, target_size=size)
        beam = reject_hard_failures(beam)
        beam = stable_top_k(score_proxy(beam), config.beam_width)
    shortlist.extend(diversify(beam, config.shortlist_per_sector))

shortlist = optional_combination_stability(shortlist)
```

候选至少满足：位点都位于同一连通表面 sector；不要求各位点在序列上连续；组合不能仅靠欧氏近距离跨越结构内部。

### 7.8 阶段 H：覆盖定义与敏感性分析

每个候选探针预计算三类覆盖向量：

#### A. 直接位点覆盖

若残基本身被突变则为 1，否则为 0。它最可解释，但低估突变对邻近界面的影响。

#### B. 局部测地覆盖

以每个突变位点为中心，在表面图上使用若干离散的距离情景。例如紧、中、宽三种核。可用阈值核或指数核：

```text
K(d; λ) = exp(-d / λ)
```

该函数只表示“空间影响权重”，绝不命名为结合概率。`λ` 必须来自配置，并在结果中分别报告，不能只输出某个默认 λ 的结论。

#### C. Footprint 拦截

对每个 microfootprint 计算：

- 是否包含至少 1 个突变位点；
- 是否包含至少 2 个突变位点；
- 被突变位点覆盖的 footprint 表面积比例。

面板优化使用预先阈值化的情景矩阵，结果报告同时展示所有情景。若不同情景导致面板显著变化，应降低推荐置信度并输出替代面板。

### 7.9 阶段 I：面板联合优化

#### 决策变量

- `x_q ∈ {0,1}`：是否选择候选探针 `q`。
- `y_us ∈ {0,1}`：在覆盖情景 `s` 下，表面单元 `u` 是否被至少一个已选探针覆盖。

对预计算的二元覆盖矩阵 `c[q,u,s]`：

```text
y_us <= Σ_q c[q,u,s] x_q
Σ_q x_q <= probe_budget
```

同一 sector 候选数量、互斥设计、突变总数和结构风险上限均可写成线性约束。连续覆盖必须先变成多个明确阈值情景，不能把非线性函数直接塞入线性求解器。

#### 分层目标

按以下顺序做多次求解，并把上一层最优值固定或允许极小容差后再优化下一层：

1. 最大化跨主要情景均能覆盖的“可设计表面积”。
2. 最小化最大未覆盖表面连通块的面积；若难以直接线性化，可先枚举/预计算组件约束，或在 fallback 阶段做局部交换优化。
3. 最大化 sector/表面方向多样性，降低多个探针覆盖同一位置的冗余。
4. 最小化结构风险。
5. 最小化突变总数和构建复杂度。

不再最大化所有残基中的最小覆盖，因为不可设计残基会让目标退化。

#### 求解器与 fallback

- 首选 `scipy.optimize.milp`/HiGHS；记录 SciPy 与 HiGHS 版本、状态、耗时和 optimality gap。
- 超时或不可用时使用确定性 greedy：每轮选择“新增稳健覆盖面积 / 风险成本”最高的候选，然后执行固定轮数的 1-swap/2-swap 局部改进。
- fallback 输出必须标记 `HEURISTIC_SOLUTION`，不能标为最优。
- 若预算不足以达到某覆盖水平，返回最佳可行解及未覆盖区域，而不是失败或暗中放宽规则。

### 7.10 阶段 J：稳定性证据

稳定性分为三级：

| 等级 | 含义 | 输出规则 |
|---|---|---|
| `RULE_ONLY` | 只有几何/理化规则代理 | 不输出 kcal/mol；标记需要实验表达与折叠验证 |
| `EXTERNAL_SINGLE` | 有单点外部能量结果 | 可用于预筛，但组合风险仍不确定 |
| `EXTERNAL_COMBINATION` | 对完整组合运行外部引擎 | 输出数值、单位、工具版本、参数、修复步骤与失败日志摘要 |

外部能量是候选优先级证据，不是稳定折叠的保证。若不同能量工具分歧，不做无依据平均；分别报告并降低置信度。

## 8. 默认参数原则

初始默认值只用于启动开发，不能在文档中被称为已验证阈值。

```yaml
joint_v3:
  max_probes: 6
  mutations_per_probe: [2, 4]
  footprint_target_areas_a2: [600, 800, 1000]
  coverage_scenarios:
    - direct_site
    - geodesic_tight
    - geodesic_medium
    - geodesic_wide
    - footprint_hit_1
    - footprint_hit_2
  candidate_search:
    site_pool_cap_per_sector: 12
    mutation_options_cap_per_site: 3
    beam_width: 200
    shortlist_per_sector: 20
  solver:
    backend: scipy_milp
    time_limit_seconds: 60
    mip_relative_gap: 0.01
    deterministic_seed: 0
```

这些数值均须：

1. 写入 manifest；
2. 可由 CLI/API 覆盖；
3. 在开发集上校准；
4. 在独立测试集冻结评估；
5. 进行至少一轮敏感性分析。

## 9. 输出合同

每次任务创建目录：

```text
results/<job_id>/
  manifest.json
  structure_audit.json
  topology_context.json
  surface_residues.csv
  microfootprints.json
  probe_sectors.json
  candidate_probes.csv
  selected_panel.json
  coverage_by_scenario.csv
  exclusions.csv
  scientific_report.md
  legacy_regions.json
```

### 9.1 `manifest.json` 必须包含

- `schema_version: "3.0"`
- 输入文件 SHA-256；
- 软件 commit、依赖版本、外部工具版本；
- 完整配置及配置哈希；
- 随机 seed；
- 运行开始/结束时间；
- 运行模式与 fallback；
- 质粒/实验 profile 名称与版本；
- 每项关键假设和用户覆盖项。

### 9.2 每个推荐探针必须解释

- 突变列表与原始残基编号；
- 所属 probe sector；
- 三类覆盖贡献；
- 对已选面板的新增覆盖；
- 结构风险与证据级别；
- 被选择原因；
- 最接近但未选择的替代候选及淘汰原因；
- 与锚定端的相对几何描述，但不伪造膜面可达概率。

### 9.3 UI 兼容

- 默认仍可将最终探针的 sector 以表面模型着色。
- surface 灰色；粗定位 sector 黄色；实际突变位点橙/红色。
- tooltip 同时显示：`Probe sector`、`Microfootprint sensitivity`、`Mutation sites`，避免三者混称 footprint。
- 旧 `RegionRecord` 由 adapter 生成，并显式带 `legacy_projection: true`。

## 10. 实验解释与数据分析约束

算法报告必须随结果给出以下解释框架：

1. 每个抗体对每个探针的结合/报告信号应先按抗原表面表达量归一化；N 端 Myc 可作为表达控制之一。
2. Myc 阳性只能支持“构建体到达细胞表面且标签可及”，不能独立证明 ECD 正确折叠。
3. 推荐设置：WT 抗原、空载体/无关抗原、每个突变探针、已知构象抗体或配体（若有）。
4. 若一个探针对多条原本差异很大的抗体均造成普遍信号丢失，同时 Myc 下降或构象对照下降，应优先判断为表达/折叠问题，而非共同表位。
5. 一个探针阴性只能说明该组合未产生可检测影响，不能证明对应区域绝不属于表位。
6. 一个探针阳性也可能来自变构、糖基化改变或局部折叠变化；因此结果是粗定位证据，不是接触残基证明。

实验数据回流、混合效应模型和自动表位归类属于 Phase 2，不阻塞 V3 设计引擎 MVP。但 V3 输出 ID 和 schema 必须为后续数据回流预留稳定主键。

## 11. 实施阶段

### Phase 0：基线冻结

- 保存当前 Surface V2 的代表性输入和结果快照。
- 固定现有测试集并加入 3–5 个结构类型不同的小型基准。
- 定义性能基线：运行时间、内存、候选数、覆盖率。

完成条件：相同代码与输入可重复生成相同 V2 结果哈希（排除时间戳字段）。

### Phase 1：结构审计与 V3 schema

- 实现 `StructureAudit`、`TopologyContext`、稳定 ResidueId。
- 完成重复链/非等价链决策与 fatal/warning 分级。
- 增加 manifest、schema version 和旧格式 adapter。

完成条件：链歧义、缺失 C 端、插入码、altloc、非标准残基均有测试。

### Phase 2：表面单元与多尺度 footprint

- 复用现有 surface point cloud/graph。
- 输出残基级面积权重。
- 实现多面积尺度、确定性中心采样与 footprint 去重。

完成条件：每个 footprint 连通；面积误差有界；相同输入完全可复现。

### Phase 3：sector 与候选生成

- 实现连通 sector。
- 建立版本化突变替换规则。
- 实现 site pool、beam search、硬排除与代理评分。

完成条件：不进行全组合枚举；大结构候选数受明确上限控制；所有淘汰均可追踪原因。

### Phase 4：稳定性适配器

- 建立统一接口与 `RULE_ONLY` fallback。
- 接入 FoldX/EvoEF 时保存完整 provenance。
- 对短名单运行组合级评估。

完成条件：外部工具缺失或失败时流程仍能完成，但报告降级且不产生伪能量值。

### Phase 5：覆盖矩阵与面板优化

- 预计算三类覆盖及多个情景。
- 实现分层 MILP。
- 实现 deterministic greedy/local-swap fallback。
- 输出未覆盖最大连通块和替代面板。

完成条件：预算 1–6 单调测试通过；MILP 与穷举小样本最优解一致；fallback 可复现。

### Phase 6：报告、API 与 UI 接入

- 导出全部机器可读文件和科学报告。
- API 增加 V3 参数与状态。
- UI 使用表面模型展示 sector、footprint 和突变点。
- 增加 shadow comparison 页面或报告段落。

完成条件：旧任务可打开；V3 任务不依赖旧字段语义；界面不把 sector 标成真实 footprint。

### Phase 7：基准与默认切换

- 使用有已知抗体–抗原复合物的公开结构构造离线 benchmark。
- 仅使用抗原侧结构生成探针，再评估其对真实抗体 footprint 的拦截率。
- 开发集调参，独立测试集只做最终评估。
- 比较 V2、V3、随机表面分区和简单均匀采样基线。

完成条件：V3 在预先定义的主要指标上优于简单基线，且没有不可接受的稳定性风险增加；之后才考虑成为默认引擎。

## 12. 测试计划

### 12.1 单元测试

- 残基编号：负编号、插入码、多模型、重复 residue number。
- 链选择：同源重复链、异源复合物、短肽配体、无明确主链。
- 表面图：连通性、缺失片段、跨内部近邻不得误连。
- footprint：目标面积、去重、边界截断、确定性。
- mutation rules：二硫键、糖基化 sequon、新增 Cys、非标准残基。
- beam search：组合位点唯一、2–4 位限制、候选上限、稳定排序。
- coverage：手工小图上的直接/测地/footprint 拦截结果。
- optimizer：预算约束、互斥约束、风险约束、不可覆盖单元不导致目标退化。

### 12.2 性质测试

- 增加 probe budget 不应降低同一目标定义下的最优可覆盖面积。
- 重命名等价 chain ID 不应改变几何结果，仅改变展示标识。
- 刚体旋转/平移结构不应改变表面分区、候选和覆盖。
- 候选输入顺序变化不应改变最终结果。
- 表面图距离不得小于明显不合理的穿模直线捷径。

### 12.3 集成测试

- 仅 PDB 输入完成全流程。
- 外部稳定性引擎不存在、超时、部分失败。
- MILP 不可用或超时后 fallback。
- 150 aa 单域 ECD、多个结构域、同源二聚体、强糖基化 ECD。
- API → 后端 → 导出 → UI 表面着色一致。

### 12.4 科学基准指标

至少报告：

- 已知真实 footprint 被任一探针至少 1 个突变位点拦截的比例；
- 被至少 2 个突变位点拦截的比例；
- 可设计表面积覆盖率；
- 最大未覆盖连通表面积；
- 探针间覆盖冗余；
- 每探针突变数与结构风险分布；
- 对 footprint 尺度、SASA 阈值、测地核和稳定性阈值的敏感性；
- 与随机、均匀分区、V2 的差值和置信区间。

不能只选择一个最有利的指标报告。

## 13. 验收标准

### 13.1 工程验收

- 同一输入、配置和版本产生同一排序与面板。
- 所有结果均可追溯到输入哈希和配置。
- 无 FoldX/EvoEF 时可完成运行且正确降级。
- 150 aa 典型 ECD 在约定硬件上的运行时间满足预设预算；具体数值在 Phase 0 测得后冻结。
- 最多 6 个探针是严格约束，不通过事后截断实现。
- 每个 probe 为 2–4 个突变，除非返回明确的“无可行设计”。
- 任何不可行、超时、fallback 都在顶层报告可见。

### 13.2 科学验收

- 输出明确区分测量、计算结果、启发式假设和未知项。
- 不输出不存在证据支持的膜平面、结合概率或折叠成功结论。
- 粗定位区始终保持表面连通。
- footprint 尺度与覆盖核至少有三个情景，并报告面板稳健性。
- 外部稳定性工具的版本、参数与组合级计算状态完整。
- 公开 benchmark 中与简单基线比较，而非只比较旧算法。

## 14. 禁止性实现规则

以下做法应在 code review 中直接拒绝：

1. 将 `sector_score > 常数` 写死并称为科学阈值。
2. 用孤立 ECD 的 PCA 最小轴作为确定膜法向量而不标不确定性。
3. 用欧氏球直接定义跨越蛋白内部的表面邻域。
4. 把所有单点 ΔΔG 相加冒充完整组合 ΔΔG，而不标代理。
5. 因 Myc 阳性而输出 `folded = true`。
6. 在没有校准的情况下输出“表位概率 0.xx”。
7. 用最终截断 `candidates[:6]` 代替预算约束优化。
8. 让不可设计残基从覆盖分母中静默消失；必须同时报告总表面覆盖与可设计表面覆盖。
9. 把 `ProbeSector` 命名为 `Footprint`。
10. 外部工具失败后填充 0 kcal/mol 或其他哨兵数值。

## 15. 推荐的首批代码任务

按依赖顺序生成代码：

1. `joint_design/models.py` 与 JSON schema。
2. `structure_audit.py` 和 residue mapping 测试。
3. `topology.py` 的 `topology_only` 实现。
4. 从 Surface V2 抽取可复用的表面图接口。
5. `surface_units.py` 与 `microfootprints.py`。
6. `sectors.py`。
7. `mutation_options.py` 和版本化替换配置。
8. `candidate_search.py` 的 deterministic beam search。
9. `coverage.py` 的三类覆盖矩阵。
10. `optimizer.py` 的小样本 MILP、穷举 oracle 和 fallback。
11. `stability.py` 适配器。
12. pipeline/API/export/UI 集成。
13. benchmark harness 和 shadow comparison。

每个任务的代码生成提示必须包含：输入/输出类型、失败模式、确定性要求、单元测试、不得改变的旧接口。不要一次提示模型重写整个项目。

## 16. 代码生成执行合同

每个实现 PR/提交应满足：

- 只完成一个可独立测试的里程碑；
- 新增公共字段必须带 schema 版本与迁移/adapter；
- 算法参数必须集中配置，禁止散落 magic numbers；
- 每个候选淘汰都有机器可读 reason code；
- 每个评分保留组成项，不只保留总分；
- 所有图算法使用表面测地关系；
- 所有随机行为固定 seed 并稳定排序；
- 优化器必须有小问题穷举 oracle 测试；
- 外部工具必须隔离在 adapter 中并可缺省运行；
- 新科学阈值必须附来源、校准记录或明确标记为待校准启发式；
- 任何 benchmark 调参都只能使用开发集，测试集结果不得反向改参数。

当以下任一条件出现时，代码生成应停止并返回明确错误，而不是猜测：

- 多个非等价链均可能是目标抗原；
- 结构无法建立稳定残基映射；
- 没有足够表面残基形成任何 2 位点候选；
- 用户要求的硬约束彼此冲突；
- 输出 schema 与旧任务兼容策略未定义。

## 17. 文献与工具依据

以下资料用于限定合理范围和工具能力，不应被解读为本项目参数已经完成验证：

1. Kringelum JV et al. *Structural analysis of B-cell epitopes in antibody:protein complexes*. Molecular Immunology, 2013. 该研究支持蛋白抗体表位通常由多个空间邻近残基构成、接触面积存在较宽分布，因此适合使用多尺度窗口而非固定 footprint。
   https://pmc.ncbi.nlm.nih.gov/articles/PMC3461403/
2. Lomize AL et al. *PPM 3.0: enhanced spatial positioning of proteins in membranes*. Nucleic Acids Research, 2022. 膜定位方法需要与膜相互作用相关的结构信息；孤立 ECD 与质粒连接拓扑本身不足以唯一确定膜平面。
   https://pmc.ncbi.nlm.nih.gov/articles/PMC8740824/
3. Delgado J et al. *FoldX 5.0: working with RNA, small molecules and a new graphical interface*. Bioinformatics, 2019. FoldX 是经验能量函数工具，适合作为稳定性证据来源之一，而非实验折叠结论。
   https://academic.oup.com/bioinformatics/article/35/20/4168/5381539
4. SciPy documentation: `scipy.optimize.milp`. 提供基于 HiGHS 的混合整数线性规划接口，可实现带预算和覆盖约束的确定性面板选择，并报告求解状态。
   https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.milp.html

## 18. 最终判断

V3 最重要的变化不是换一个聚类算法，而是改变问题定义：从“先预测 footprint，再在每个区挑突变”改为“在多个合理表面与影响模型下，联合选择少量、局部、低风险且覆盖互补的实验探针”。这更符合当前实验能力，也能把无法从 PDB 得到的知识明确留给实验，而不是用先验假装成计算结论。
