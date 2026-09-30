# ECD Epitope Workbench：论文级方法学（Methods）

## 1. 研究目的与总体设计

本研究建立一个面向膜蛋白胞外结构域（extracellular domain, ECD）构象表位定位的计算—实验设计工作台。系统的输入为目标蛋白序列、ECD 边界、可用的编码序列（CDS）及结构/置信度证据；输出为一组可直接用于结合实验的 ECD 或全长突变构建体。每个构建体被视为测量表位位置的实验仪器：其突变应在尽量不破坏整体折叠、表达和膜定位的前提下，对局部抗体足迹产生足够扰动。

分析流程由结构准入、残基映射、表面可及性计算、表面补丁划分、突变菜单生成、稳定性与功能风险过滤、组合优化、序列重建和可重复导出组成。系统当前以本机单用户模式运行，计算核心不依赖登录或云服务；数据模型、工件接口和作业归属边界为未来公开部署预留替换点。

## 2. 输入数据与质量控制

### 2.1 序列和边界

目标蛋白序列由用户输入或由 UniProt accession 解析获得。ECD 起止位置以 UniProt 绝对坐标表示，并在进入结构分析前验证：起点、终点均为正整数，起点不大于终点，边界位于目标序列内，ECD 长度与结构链覆盖范围可计算。原始 CDS 若存在，则保留为不可变输入快照；任何 ECD 优化序列均从该快照重新生成，避免仅保存 ECD 片段而无法恢复全长构建体。

### 2.2 结构候选和映射

结构候选优先来自 PDB/mmCIF 或 AF3 预测结构。结构链与目标序列通过序列比对和残基编号映射建立一一对应关系，并记录 accession、链标识、覆盖区间、序列一致性、内部缺口和结构域覆盖率。当前准入参数由运行配置统一提供：ECD 覆盖率下限为 0.95，允许的最大内部缺口为 8 个残基，结构—目标序列一致性下限为 0.95，主要结构域覆盖率下限为 0.90。未满足准入条件的结构不进入后续补丁和突变设计；所有拒绝原因写入诊断事件和 manifest。

AF3 confidence JSON 仅在其残基/原子数组能够与目标链完成映射时使用。多链、缺失数组、长度错配或链标识不一致均产生明确 warning；无法确认映射时回退到结构文件中的可靠 B-factor，若仍无法确认则终止该结构的置信度加权并要求人工复核，不静默地产生高置信度结论。

## 3. 结构预处理与 SASA 计算

结构解析采用 Bio.PDB。生产环境的溶剂可及表面积（solvent-accessible surface area, SASA）固定使用 Shrake–Rupley 算法，探针半径为 1.4 Å，采样点数为 960，原子半径采用 Tien/DSSP 参数集：C=1.80 Å、N=1.60 Å、O=1.40 Å、S=1.85 Å。FreeSASA 不参与生产筛选、分区或突变面板，仅可作为诊断性对照，并且不得改变最终结果。

对每个映射残基计算总 SASA 和侧链 SASA，并以残基类型的最大可及面积归一化：

<div align="center"><i>rSASA</i><sub>i</sub> = <i>SASA</i><sub>i</sub> / <i>ASA</i><sub>max</sub>(<i>aa</i><sub>i</sub>)</div>

残基被标记为表面暴露的条件为 native total `rSASA ≥ 0.25`。侧链 SASA 为描述证据，统一 `20 Å²` 条件已移除，避免 Gly 和小残基的面积偏差。坐标可靠性、保护位点与突变稳定性分别审核。SASA 方法、版本、半径、采样点数、归一化表和实际阈值均写入 manifest 与报告；0.15/0.20/0.25 敏感性比较不自动降低默认值。

同时提取以下结构特征：二级结构、环区/转角、凸出度、膜邻近性、结构域边界、糖基化 sequon、已知或预测糖基、二硫键、翻译后修饰、蛋白—蛋白界面标记及预测模型的 pLDDT。Cys、Pro、二硫键残基、糖基化位点和翻译后修饰位点默认进入突变保护集合。实验结构的 B-factor 不被解释为 pLDDT 或无序证据。

### 3.1 结构可靠性和无序区分流

在 SASA 计算完成后、Surface V2 几何建模开始前，系统以保守的多证据规则对 ECD 每个 UniProt 位置分类。该层只决定残基是否进入三维表面分区；它不改写输入结构、原始 SASA、残基暴露判定、突变菜单、FoldX/EvoEF2、Stage-2 或序列重建。其输出不是经校准的无序概率，也不能用于断言结合诱导折叠。

预测结构中的 pLDDT 仅作为局部结构置信度：`≥70` 为高置信结构候选，`<50` 为低置信候选。PAE 只有在 AF3 token 与结构链/残基可验证映射时才使用：优先使用显式链—残基 token；没有显式 token 时仅接受维度与已映射单链一致的矩阵。PAE 无法安全映射、多链含义不明确、矩阵无效或证据缺失时，PAE 记为不可用而非正向无序证据。核心 PAE 优先相对于同一结构域的高置信残基计算，以避免将整体可动但局部折叠的结构域误判为无序。

除 pLDDT 与核心 PAE 外，系统计算 8 Å 半径内、序列距离大于 2 的非局部 Cα 接触数；可选导入完整 canonical 序列的 IUPred 输出（位置、氨基酸、分数）。导入文件必须完整覆盖全序列、编号和氨基酸严格一致、分数在 `[0,1]`，其 SHA-256 与头信息写入 QC；系统不自动安装或运行 IUPred。缺失任一证据不会被填补为“无序”。

分类状态为：`STRUCTURED`（有坐标的实验结构，或预测 pLDDT≥70）、`FLEXIBLE_STRUCTURED`（由同域高置信两端支撑、长度≤14 aa 的低置信嵌入 loop）、`UNCERTAIN`（证据不足，保留）、`DISORDERED`、`CONDITIONALLY_FOLDED` 和 `UNRESOLVED`（缺失可用坐标）。`DISORDERED` 需是连续候选段，默认长度≥15 aa、至少 80% pLDDT<50、至少 80% 残基低非局部接触，并且至少有序列无序或高核心 PAE 支持，且不含 H/E 二级结构冲突。结构连接端最多保留 3 个 `UNCERTAIN` 过渡残基。`CONDITIONALLY_FOLDED` 需同时具有高置信螺旋、序列无序支持、高核心 PAE 和低接触，只作为条件折叠假设，默认不进入表面几何。7 aa 居中中位数证据分数仅供诊断；最终判定使用连续片段共识，不按该分数阈值切分。

默认 `strict` 几何模式下，`DISORDERED`、`UNRESOLVED` 与默认排除的 `CONDITIONALLY_FOLDED` 残基既不贡献采样顶点，也不作为 ECD 遮挡原子；`occlusion_sensitivity` 模式保留它们的遮挡作用但仍不采样顶点，用于敏感性比较。屏蔽或未解析位置相邻 2 个残基保留遮挡原子、但不贡献表面顶点，以降低人为切口成为采样中心的风险；此为整残基边界保护，而不是网格切口重建。`UNCERTAIN` 和 `FLEXIBLE_STRUCTURED` 保留在几何中，并在候选足迹 QC 中标记。原始 SASA 始终取自未屏蔽结构，不能解释为屏蔽后的可及面积。

## 4. Surface V2 表面补丁划分和表位先验

### 4.1 几何图和补丁枚举

Surface V2 使用独立于生产 SASA 的确定性 `sr_pointcloud` 表面几何，仅为分区提供图结构，绝不替换第 3 节的 SASA。对可靠性分流后保留的 ECD 重原子，以 1.4 Å 探针扩张原子半径，在每原子 48 个 Fibonacci 球面方向上采样；落入其他原子排除体积的点被剔除。顶点间在距离≤2.5 Å、法向余弦≥0.30 且原子间线段未穿越排除体积时连边，形成显式溶剂表面点图。顶点面积、所属残基/原子、法向、连通分量和置信度随结果导出。

仅来自 SASA 暴露残基、pLDDT≥70（若有预测置信度）且所在表面分量面积≥25 Å² 的顶点可作为中心候选。最远点优先采样使用 7 Å 中心间距与 7 Å 覆盖半径。每个中心生成 9、12、15 Å 测地窗口以及目标面积 600 Å² 的自适应窗口（允许半径 8–16 Å、面积 350–1000 Å²）；窗口可重叠。残基纳入 footprint 的默认条件是该残基拥有的表面面积≥8 Å²，且其表面面积中落入窗口的比例≥0.20。该定义以原子表面测地距离而非 Cα 邻接替代旧的“固定残基图、13 Å 单一补丁”规则。

每个 footprint 计算连通性、测地直径、法向离散度、低置信表面积比例、残基数和可突变位点等 QC。按加权顶点面积 Jaccard 相似度≥0.75 去除冗余 footprint；原始与去重后的窗口均导出，避免把去重误解为未枚举。Surface V2 的“所有自然区域”与 Stage-1 实验预算分离：前者完整保留，后者最多选择 7 个设计区域。

去重 footprint 以重叠、共享表面积、结构域一致性和表面图的测地关系构图，再以固定随机种子 0 的、主尺度受约束 Louvain 社区方法聚合。聚合限制包括中心间测地距离、社区中心跨度、最小区域面积和连通性修复；最终每个可显示表面残基仅归属一个区域颜色，重叠 footprint 仍作为单独对象可查看。区域聚类参数、种子、实际区域数、覆盖率、中心数、足迹数、连通分量面积和所有 QC 均写入 `surface_v2_qc.json`。

### 4.2 兼容的紧凑单结构域异常修复模式

通用足迹分区首先照常运行。只有在 ECD 不超过 180 aa、一个主要结构域覆盖可信表面包络至少 70%，且通用结果出现区域数超过 3、测地跨度超过两个足迹半径或可突变位点相对表面包络过度稀疏时，才进入紧凑单结构域模式。未触发时，通用分区及其突变输出保持不变。

兼容模式使用与生产 SASA 相同的 1.4 Å 探针和 Tien/DSSP 原子半径，在每个重原子扩张球面上取 96 个确定性 Fibonacci 采样点。被邻近原子排除体积遮挡的点被删除，剩余点的单位方向向量经残基内合成后得到局部溶剂表面法向。该法向不依赖蛋白整体质心，因而适用于扁平 β-sheet、凹面和非球形结构。常规新作业默认采用 4.1 的 Surface V2；本节仅定义旧的兼容回退路径。

对于具有可信 UniProt Beta strand 注释的紧凑 β-sandwich，分区优先使用三维 β-sheet 骨架拓扑。两条非序列相邻的 β-strand 若至少具有 3 对 Cα 距离不超过 6.5 Å 的接触，则建立 sheet 邻接；序列间隔不超过 5 个残基的紧邻片段不用于跨 sheet 连接，以避免 β-turn 将两张 sheet 错误合并。未被主 sheet 直接覆盖的近表面残基按最近 sheet 骨架归属。缺少可信 β-strand 拓扑时，才使用溶剂表面法向图，并要求区域连通且平均法向一致性不低于 0.65。

每个紧凑表面区域最多生成 3 个非冗余局部 footprint。footprint 测地半径不超过 13 Å，中心欧氏距离不超过 18 Å；保护残基参与表面几何但不进入突变菜单。多点探针必须完整属于同一个 footprint，突变间最大跨度不超过 15 Å。

### 4.3 先验评分

补丁先验由可解释特征线性组合得到：暴露度权重 0.30、凸出度 0.20、环区权重 0.15、序列变异性权重 0.15、界面权重 0.05；糖基化/糖链遮挡施加 0.10 惩罚。`cell_surface` 试验格式对膜近端补丁施加 0.50 的乘法惩罚，并在报告中注明 IgG 可及性限制；`purified` 格式不施加该惩罚。先验只用于候选排序，不作为单独的阳性证据。

## 5. 保守性、语言模型与功能风险

保守性分析优先使用 UniProt reviewed、同基因正交同源序列。序列按 PairwiseAligner 映射到目标序列；至少 3 条有效同源序列才计算位点保守性，否则返回缺失并记录原因。位点变异性可增加表位先验，但对高度保守或功能关键位点施加风险惩罚。

ESM-1v 为可选 zero-shot 证据，仅在显式提供 `ECD_ESM_MODEL` 且 torch/esm 可用时执行。对每个候选替换计算突变相对似然分数，并与结构稳定性证据共同用于排序。模型不可用、超时或输出异常时优雅降级为缺失，不阻断可审计的结构规则流程。

## 6. Surface V2/legacy 突变菜单与单点筛选

本节描述默认 `surface_v2`/legacy 工作流。Joint V3 使用独立的突变选项和组合搜索规则，见第 7.1 节；两套引擎的结果不能仅凭字段名称互相解释。

系统在允许的表面位点上枚举氨基酸替换。默认定向菜单为：体积较大或芳香/带正电残基（W/Y/F/R/K/M/I/L）优先替换为 Ala；小体积或极性残基（G/A/S/T/D/N）优先替换为 Trp；带电残基可进行电荷反转。每个位点最多保留 2 个候选。锚点突变定义为电荷反转，或 W/Y/F→A；每个主构建尽量包含至少一个锚点。

以下残基或情形默认排除：Cys、Pro、二硫键、糖基化 sequon/已知糖基化位点、已知翻译后修饰位点、埋藏位点和跨越结构域核心的替换。可靠性分流排除的残基不进入三维 footprint，也不进入自动突变菜单；`UNCERTAIN` 和 `FLEXIBLE_STRUCTURED` 残基若落入 footprint 只保留 QC 标记，不因低 pLDDT 单独剔除。单点候选需同时满足表面暴露、结构映射、置信度和功能风险门槛。稳定性门槛按暴露度联动：`rSASA≥0.35` 时允许预测 ΔΔG 至少 0.80 kcal/mol，`0.25≤rSASA<0.35` 时至少 0.60 kcal/mol；无法获得稳定性评分时标记为不确定而非假定安全。

## 7. Surface V2/legacy 的 FoldX/EvoEF2 评估和组合设计

FoldX 与 EvoEF2 属于可替换的稳定性评估引擎。生产流程对候选突变进行确定性解析、超时控制、空输出校验和缓存；失败仅保留诊断信息，不将空输出解释为低风险。跨补丁候选先汇总单点结果，再按统一批次执行，确保评分、阈值、顺序和审计记录一致。成功作业完成后自动清理 `engine_runs` 中间目录，失败作业保留诊断目录。

通用模式组合设计在 3–6 个突变位点范围内进行，并尽量包含至少一个锚点。候选之间的最小空间间隔为 5.5 Å，组合最大空间跨度为 22 Å；每个补丁最多保留 5 个关键位点。紧凑单域模式在每个区域最多扫描 12 个优先位置，多点组合必须属于同一候选 footprint，最大跨度为 15 Å。若所有局部多点组合均失败，但存在通过 FoldX 门槛的强扰动单点，则输出 `RECOMMENDED_SINGLE_SITE_PROBE`；该单点阴性不能排除整个区域，报告中必须保留此解释限制。组合稳定性门槛为预测总 ΔΔG 不高于 +2.0 kcal/mol（或按配置快照中的等价组合门槛执行）。候选在“结合扰动潜力—结构稳定性风险—补丁覆盖度”三维空间形成 Pareto 前沿，选取膝点作为主探针，保留一个更激进的备选探针。

### 7.1 Joint V3：联合表面表示、候选搜索与面板优化

Joint V3 是显式选择的可选后端，默认不替换 Surface V2。它将问题定义为：在不超过配置预算的少量构建体内，尽量覆盖可设计 ECD 表面，使实验响应能够把结合区域粗分到空间区域；它不预测某条抗体的精确表位。

#### 结构审计和拓扑上下文

V3 复用前述结构解析和残基映射，但额外生成 `StructureAudit` 和 `TopologyContext`。审计记录选定 model/chain、缺失片段、非标准残基、映射状态和质量警告。残基内部主键包含 model、chain、author residue number、insertion code 和氨基酸，展示层才格式化为简写编号。多链或非等价链的选择不能被静默解释为生物学唯一性。

对于孤立 ECD，默认拓扑模式只记录 C 端锚点、相对轴和高不确定性。相对轴可用于“靠近/远离连接端”的描述，但不是膜法向、膜平面或 Fab 可达概率。只有输入结构真实包含可识别 TM 坐标时，才允许输出 `tm_geometry` 上下文，并记录方法和输入 provenance。

#### 表面残基单元和多尺度窗口

V3 的覆盖优化基本单元是 `SurfaceResidueUnit`，而不是人为创建的固定大小 `SurfaceUnit` 分区。表面残基保留总表面积、rSASA、可设计状态、排除原因和连通组件。二硫键 Cys、糖基化 sequon、缺失/低质量位置、明显界面或结构风险可被保护，但保护残基仍保留在总表面和不可覆盖表面的报告中。

表面图使用表面候选残基间的保守三维近邻关系和测地距离；不会用穿过蛋白内部的欧氏球直接连接两侧表面。沿图按测地距离生长多个 `MicroFootprint`，默认敏感性面积为约 600、800、1000 Å²。每个窗口记录中心、目标面积、实际面积、残基集合和边界截断；这些面积是敏感性情景，不是所有抗体的固定真实接触面积。

相邻且重叠的窗口合并为连通 `ProbeSector`。合并必须满足同一残基图连通分量、测地直径不超过 35 Å 和 Cα 三维跨度不超过 30 Å；显式表面点图折叠为残基图时，超过 10 Å 的 Cα 长边不作为连接。这样可避免一个分区包含空间上相距很远的不同片段。Sector 用于少量构建体的粗定位，不能在 API、报告或网页中被称为真实 antibody footprint。现有网页仍接收兼容投影的 `RegionRecord`；该投影带有 `legacy_projection=true`，不改变 V3 的内部语义。

#### 突变选项和确定性 beam search

V3 默认保护二硫键、糖基化位点/sequon、明显功能位点和结构缺失区域，并检查单点及组合突变是否引入新的 N-X-S/T sequon 或游离 Cys。每个位点只生成少量有明确扰动机制的替换，例如去除电荷、降低疏水性、芳香侧链去除或有限的极性变化；不把统一 Ala 扫描作为唯一策略。

每个 sector 先按表面可及性、空间代表性、化学扰动和规则风险截取有限位点池，再从 2 位点扩展到 3、4 位点。候选不进行全局 2–4 位点组合枚举，不强制包含 anchor mutation。beam 的排序键包含评分、链、残基位置、插入码和替换氨基酸，确保同一输入和配置得到相同候选顺序。

没有外部组合级稳定性工具时，V3 只保留 `RULE_ONLY`/`PROXY_ONLY` 证据等级和无量纲风险，不输出 kcal/mol。当前 V3 FoldX 适配器先用代理风险生成预选面板，优先将预选及替代面板候选纳入组合级 `RepairPDB + BuildModel`，剩余配额按 sector 轮转分配；最终面板只从成功取得组合 ΔΔG 的候选中选择，组合值同步写入 `DesignRecord.combined_stability_score`。默认 `require_pass_for_final_panel=true`：`ABOVE_GATE` 结果只保留在审计和排序信息中，不得进入最终推荐面板；若 PASS 候选不足预算，输出短缺而不以不合格组合凑数。适配器保存工具路径/版本、参数、输入与修复结构、运行次数、候选状态、失败摘要和组合级结果；单点分数相加不参与 V3 能量证据。超出配置候选上限、修复失败、超时或空输出的候选保留明确降级状态。

#### 多情景覆盖与面板优化

每个候选预计算三类覆盖：直接突变位点覆盖、表面图上的 tight/medium/wide 测地覆盖，以及 `MicroFootprint` 的至少 1 位点或至少 2 位点拦截。测地指数核或距离阈值只表示空间影响权重，不能命名为结合概率。

面板优化的主要目标是跨多个覆盖情景的可设计表面积稳健覆盖，其次控制未覆盖区域、冗余、结构风险和突变数；不可设计残基不会从总表面分母中静默消失。首选 SciPy/HiGHS MILP，并记录 solver 状态和 optimality gap；求解器不可用、超时或无可行解时使用确定性 greedy/local fallback，并标记 `HEURISTIC_SOLUTION`。预算是优化约束，不通过事后截断候选列表实现。

若没有任何满足 2–4 位点和保护规则的候选，V3 返回明确的 `NO_FEASIBLE_PROBE` 状态，并在 design summary 中保存无可行设计原因、淘汰原因和未覆盖表面，不制造一个看似完整的探针。V3 输出的覆盖率是结构情景下的计算结果，不是结合概率或科学合格线。

### 7.2 Joint V3 稳定性与实验解释边界

V3 的稳定性证据分为 `RULE_ONLY`、`EXTERNAL_SINGLE` 和 `EXTERNAL_COMBINATION`。外部能量工具仅用于候选优先级和风险分层，不能证明细胞表面正确折叠。实验解释仍要求 WT、无抗原/无关抗原、突变探针、表达归一化和独立构象敏感性对照；Myc 阳性只能支持表面展示/标签可及性，不能单独证明天然折叠。

探针阴性只能说明该组合未产生可检测影响，不能证明对应区域不属于抗体表位。探针阳性也可能由变构、糖基化改变或局部折叠变化造成。因此 V3 的输出应称为“粗定位证据”和“实验面板设计”，而不是“表位预测”。

### 7.3 Joint V4-ML：多工具证据、规范化映射与正式推荐门控

V4 的完整独立方法学见 [`method_joint_v4.md`](method_joint_v4.md)；本节保留与通用 Methods 的接口和解释边界摘要。

Joint V4-ML 是独立于默认 Surface V2 和 Joint V3 的 opt-in 后端。它把 ProCaliper 的表面可及性、DiscoTope-3.0 的表位倾向、SPURS single/multi 的稳定性证据和 FoldX 的组合复核作为四条可审计证据线；这些分数只用于候选筛选和实验面板排序，不被解释为抗体结合概率或真实接触预测。V4 不捆绑外部模型权重、二进制或许可证。

#### 预测上下文与位置映射

对输入结构选定的目标链生成连续编号的 `prediction_context.pdb`。预测请求使用该规范化链及其连续序列；原始结构、author residue number、插入码和 UniProt 位置仍作为生物学上下文保存。`residue_mapping.csv` 至少记录 `normalized_position`、`uniprot_pos`、`chain_id`、`auth_seq_id`、`insertion_code`、`label_seq_id` 和 WT 氨基酸。若目标链不存在、残基数量不一致、结构残基 WT 与输入序列不一致或映射存在歧义，运行失败并返回稳定错误码，不静默按枚举顺序补齐。

规范化链上的外部结果先通过位置范围、WT 和输出字段校验，再映射回 UniProt 坐标。ProCaliper/DiscoTope 成功输出必须覆盖目标链的全部规范化位置；SPURS single 必须提供每个位置的 19 个非 WT 替换分数；SPURS multi 与 FoldX 必须完整返回每个请求过的组合，并且组合中的 WT、位置和突变格式必须一致。worker 进程的 stdout 只能包含一个 JSON 文档，日志写入 stderr；超时、崩溃、空输出、污染输出或 schema 错误均保留诊断并使该证据线失败。

#### 表面分区与候选生成

V4 在已计算的自然表面残基图上建立表位先验。DiscoTope 分数只用于高分种子和局部岛屿识别；表面图的每个连通分量随后按到表位岛种子的最短路径进行 Voronoi 分配，从而覆盖整个连通表面。没有表位种子的分量仍保留为 `COVERAGE_ONLY`，不能从覆盖分母中删除。超过配置 sector 数时只合并图上相邻的分区；不连通分量不被强行合并。

每个位点枚举 19 个非 WT 替换，先应用 Cys/Pro、二硫键、糖基化 sequon、功能保护和结构风险等硬排除，再按 SPURS single 分数分为：

- `PREFERRED`：ddG ≤ 配置的 preferred 阈值；
- `ACCEPTABLE`：preferred < ddG ≤ acceptable 阈值；
- `REVIEW`：acceptable < ddG ≤ review 阈值；
- `REJECT`：超过 review 阈值；
- `NOT_SCORED`：没有 SPURS single 结果。

strict 模式只有 `PREFERRED` 和 `ACCEPTABLE` 单点进入正式候选池；preview 可以保留 `NOT_SCORED`，但所有输出只能标记为 `PREVIEW_CANDIDATE`。候选按 2、3、4 点组合生成，组合位点必须在同一有效连通 sector 内并满足配置的测地跨度。每个候选保存组成单点风险、覆盖情景、组合阶数、外推等级和排除原因。

#### 组合稳定性与候选级放行

对候选计算：

<div align="center"><i>epistasis</i> = ΔΔG<sub>SPURS,multi</sub> − Σ ΔΔG<sub>SPURS,single</sub></div>

当 SPURS multi 与 FoldX 均有结果时，差值超过配置的 1.5 kcal/mol 默认标记为 `MODEL_DISAGREEMENT`；两者同时通过组合门槛才为 `CONSENSUS_PASS`。仅有 SPURS multi 时为 `SPURS_ONLY_REVIEW_REQUIRED`，只有在运行参数显式允许时才可作为正式候选。其余情况分别标记 `CONSENSUS_REVIEW`、`FOLDX_ONLY_LEGACY`、`COMBINATION_REJECT` 或 `NOT_SCORED`。

三/四点组合还必须比同一 sector 的最佳低阶候选提高至少配置的覆盖增益，并且必须有 FoldX 结果；否则不能进入正式系列。DiscoTope 仅作为软排序先验：达到 seed 分位数或局部峰值的组合标记为 `EPITOPE_PRIORITIZED`，其余通过同一稳定性门控的组合标记为 `STABILITY_FALLBACK`。面板默认每个 sector 最多选择 6 个构建体（优先 3+3 分层），总数最多 7 个；FoldX 最多复核 42 个候选，高阶组合受独立数量上限约束。选择优先组合稳定性和低阶证据，再考虑覆盖与 epitope prior；没有通过候选时保留 `NO_FEASIBLE_CANDIDATE`，不用不合格候选补齐预算。

#### 运行模式和推荐语义

`strict` 的顶层状态只有在稳定性必需证据线完成且最终推荐系列均通过候选级门控时才可能为 `FORMAL_RECOMMENDATION`。DiscoTope 缺失会产生审计 warning 并使候选归入 `STABILITY_FALLBACK`，但不单独阻断；缺少稳定性证据、没有可行候选或出现模型分歧时为 `BLOCKED_EXTERNAL_EVIDENCE`。`preview` 始终为 `PREVIEW_NOT_FOR_FORMAL_RECOMMENDATION`，不填充伪造 kcal/mol，也不把兼容投影标记为推荐。旧网页中的 `regions.json`/`designs.json` 只是语义兼容投影；V4 的正式解释应以 `v4_manifest.json`、`panel_candidates.json`、`selected_panel_v4.json`、组合稳定性表和工具 inventory 为准。

V4 的 worker 是本地隔离进程协议，不等同于已安装的模型。可配置 `ECD_V4_PROCALIPER_COMMAND`、`ECD_V4_SPURS_SINGLE_COMMAND`、`ECD_V4_SPURS_MULTI_COMMAND`、`ECD_V4_DISCOTOPE_COMMAND` 和 `ECD_V4_FOLDX_COMMAND`；`ECD_V4_SPURS_COMMAND` 可作为 single/multi 的兼容后备。没有这些工具或对应权重时，系统应明确阻断 strict，而不是用代理值填充外部模型结果。

## 8. Stage-1/Stage-2 实验面板

对于 Surface V2/legacy，Stage-1 为跨越主要表面补丁的主组合，优先覆盖不同结构域及域间边界；备选组合用于在主组合表达或折叠不佳时提供替代。对于 Joint V3，Stage-1 由面板优化器在最多配置预算内联合选择，所属对象是 `ProbeSector`，不强制 anchor mutation，也不将 sector 称为真实 footprint。Stage-2 在抗体结合下降后，对贡献最大的单点和相邻配对进行精细化验证，属于后续实验回流，不改变 V3 的结构覆盖定义。每个设计记录 `design_role`、`parent_design_id`、突变集合、所属区域、稳定性证据和风险标记；旧网页使用的 `region_id` 是兼容投影，不应覆盖 V3 的 sector/micro-footprint 语义。

## 9. ECD、全长蛋白和 CDS 重建

对每个设计，先在 ECD 氨基酸序列上应用绝对坐标突变并进行长度、字符和突变一致性校验。若存在可解析的原始全长 CDS，则根据原始 CDS 的密码子边界将对应氨基酸替换回全长编码序列，同时重新生成全长蛋白序列和 GenBank 注释；原始 CDS 永久保留，派生产物采用原子化写入。无可用 CDS 时，仅导出优化 ECD 片段，并明确标记其不能代表全长膜蛋白构建体；全长蛋白可作为合成公司的反向翻译输入，但不宣称保留原始密码子优化策略。

术语和复制语义严格区分：ECD 氨基酸按钮显示“复制 ECD”，内容为当前 WT 或构建体的 ECD 序列；“复制 CDS”仅复制对应 ECD 编码序列，缺少 CDS 时禁用并解释原因。未来若提供全长复制，使用独立的“复制全长蛋白”和“复制全长 CDS”入口。

## 10. 报告、离线 HTML 与预览状态

每个作业生成 manifest、methods report、JSON、CSV、FASTA、GenBank 和离线 HTML。manifest 包含输入哈希、结构 accession/链、配置快照、SASA 完整参数、引擎版本、缓存命中、警告和失败诊断；输出文件采用临时文件加原子替换，避免半写入文件被下载。结构可靠性层另导出 `residue_reliability.csv`、`structural_segments.csv`、`surface_mask.csv`、`idr_candidates.csv` 与 `reliability_qc.json`，并纳入 ZIP 和 manifest。IDR/条件折叠候选使用 12 aa 窗口、步长 5 并补齐末端；SASA 明确为空，候选状态为 `REVIEW_REQUIRED`，不自动生成突变或调用 FoldX。C/G/P、N-X-S/T、既有 PTM 和二硫键位置在这些序列窗口中仍保留保护标记。

启用 Joint V3 时，额外生成 `manifest.json`、`structure_audit.json`、`topology_context.json`、`surface_residues.csv`、`microfootprints.json`、`probe_sectors.json`、`candidate_probes.csv`、`coverage_by_scenario.csv`、`selected_panel.json`、`exclusions.csv` 和 `scientific_report.md`。其中 `manifest.json` 保存 V3 schema、输入结构校验和、运行环境、完整配置及配置哈希；`selected_panel.json` 保存候选 ID、覆盖情景、solver 状态、替代面板和未覆盖/保护表面；`candidate_probes.csv` 保存每个候选的突变、sector、风险、覆盖与稳定性证据级别。旧网页使用选中面板 sector 的兼容投影 `regions.json`/`designs.json`，不能据此把 ProbeSector 重新解释为真实 footprint。

启用 Joint V4-ML 时，额外生成 `tool_inventory.json`、`structure_contexts.json`、`prediction_context.pdb`、`residue_mapping.csv`、`surface_evidence.csv`、`sasa_comparison.csv`、`spurs_single_matrix.csv.gz`、`single_mutation_options.csv`、`discotope_residues.csv`、`epitope_islands.json`、`probe_sectors_v4.json`、`panel_candidates.json`、`candidate_combinations.csv.gz`、`combination_stability.csv`、`selected_panel_v4.json`、`v4_manifest.json` 和 `v4_scientific_report.md`。`tool_inventory.json` 保存工具可用性、命令、版本/模型 revision、设备、权重摘要和失败原因；manifest 保存模式、配置快照、上下文、缓存键、警告和错误。`residue_mapping.csv` 是外部结果坐标回溯的权威表。

实时工作台和离线 HTML 共享三种互斥预览状态：`WT`、`推荐构建(id)` 和 `自定义草稿`。推荐构建由其 ID 对应的 `mutations` 直接生成，不依赖此前点击或界面重绘副作用；点击突变 chip 只定位残基，不改变 active construct ID。自定义草稿独立保存，仅在用户主动进入“预览自定义草稿”时生效。构建卡片、下拉框、突变高亮、序列坞、复制内容、三维着色和标题均绑定同一 active construct ID，并显示当前来源。

序列区采用顶部工具栏—可收缩主工作区—固定序列区的 Grid/Flex 布局，序列区独立滚动；低高度或窄宽度窗口允许页面滚动，不裁切序列、复制按钮、三维视图或分析面板。

## 11. 本机运行、安全和数据治理

默认运行模式为 `local`，服务仅绑定 `127.0.0.1`、`localhost` 或 `::1`。内部 Principal 固定解析为唯一 `local` 用户；作业、事件、文件、导出和 custom design 操作均先通过作业归属接口。历史作业迁移时补齐本机归属字段。API 仅返回受控下载端点，不暴露物理绝对路径；文件解析执行路径纵深校验。

作业列表返回摘要，完整结果需访问单作业详情。成功作业清理 FoldX/EvoEF2 中间目录，失败作业保留诊断；repair cache 按容量或最近使用时间淘汰。历史结果和测试产物只提供预览式清理，实际删除必须显式确认。启动健康检查验证服务身份；停止操作通过 PID 文件和命令行校验，避免误杀其他进程。公开网站所需的认证、HTTPS、允许域名、下载大小限制、异步队列和多租户授权属于后续运行模式，不在本机模式中启用。

## 12. 可重复性与验证标准

配置文件是全部科学阈值的唯一来源；修改配置应在测试中引起可观察的筛选变化。生产结果不得因 FreeSASA 是否安装而改变。测试覆盖配置运行时行为、SASA 固定参数、Surface V2 表面图/去重/区域确定性、可靠性长尾屏蔽、短嵌入 loop 保留、实验结构 B-factor 不作为无序、PAE 映射失败降级、严格与遮挡敏感性模式、边界保护和空表面、AF3 多链/错配、CDS 继承、FoldX 超时/空输出/缓存/清理、构建连续点击、草稿隔离、ECD/CDS 复制语义和低高度布局；Joint V3 另覆盖确定性候选/面板、2–4 位候选保留与全局去重、保护掩码一致性、表面捷径防护、覆盖分母与 footprint 去重、多情景 MILP、拓扑不确定性、FoldX 配额优先级和组合能量投影；Joint V4 覆盖 worker 超时/崩溃/输出污染、规范化残基映射、single/multi 配置隔离、严格候选级证据门控、模型分歧阻断、预览/阻断兼容投影和完整表面分区。具体回归结果以当前测试命令输出为准。

所有结果必须随 manifest 保存运行配置、输入摘要、方法参数和诊断事件，从而允许第三方在相同输入、依赖版本和配置快照下重建结果。计算排名不等同于表位确证；最终结论需通过抗体结合实验，并结合表达、折叠、膜定位或纯化质量对照进行解释。

## 13. 方法学边界与局限

该流程定位的是“可被当前结构和突变面板检出的表面表位”，不能保证发现柔性、瞬态、糖链依赖或构象重排形成的全部表位。SASA、Surface V2 点云、FoldX/EvoEF2、ESM 和结构可靠性规则均为近似证据，不能替代实验结构生物学。严格屏蔽可消除不可靠片段对本轮几何分区的影响，但不会重建其真实构象；遮挡敏感性比较也不能代替多种真实构象或 AF 种子的验证。`CONDITIONALLY_FOLDED` 只是假设，IDR 窗口也不是确证表位或获批准突变。结构覆盖不足、AF3 置信度不可靠、CDS 缺失或关键位点被保护时，系统应输出限制说明和人工复核提示，而不是生成无条件的确定性结论。

Joint V3 进一步把不确定性显式放入多个覆盖情景和面板目标中，但这不会消除结构模型、糖基化状态、构象变化、细胞表面表达和抗体几何取向的未知项。V3 的“稳健覆盖率”只表示在配置的结构情景中被计算覆盖的可设计表面积比例；它不是抗体结合概率、表位置信度，也不是未经独立 benchmark 校准的 QC 合格线。只有在公开复合物 benchmark 和实验回流数据完成独立校准后，才能冻结更强的放行标准。

Joint V4 的正式推荐还受外部工具实际安装、权重、版本、许可证和 worker 输出质量限制。通过候选级门控只表示“在给定结构和工具证据下未发现配置的计算风险”，不表示正确折叠、表达、膜表面展示、抗体结合或功能保持。DiscoTope 的高分是表位倾向先验；SPURS/FoldX 的 ddG 是模型预测；模型共识也不能消除工具系统误差。`preview` 适合调试和 UI 审查，不得用于正式实验放行；`BLOCKED_EXTERNAL_EVIDENCE` 应由人工检查缺失工具、输出错误或没有安全候选的具体原因。

## 14. 参考方法来源

参数选择与解释参考 Kringelum 等关于构象表位几何的研究、Bogan 与 Thorn 关于界面热点的分析、Tien 等关于 DSSP 原子半径的工作，以及 Meier 等关于 ESM-1v zero-shot 突变效应预测的研究。pLDDT/PAE 的解释遵循 [AlphaFold Database FAQ](https://alphafold.ebi.ac.uk/faq)：局部 pLDDT 与相对位置 PAE 不能单独证明蛋白无序。具体版本、DOI/PMID 和软件许可证应在项目发布时与本 Methods 一并固化。
