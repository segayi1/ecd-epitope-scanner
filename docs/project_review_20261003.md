**ECD Epitope Scanner 全项目第一性原理审阅（2026-10-03）**

本次审阅的结论是：项目已形成较好的科学证据保留机制，但还不能把“约 8 个构建体扫描绝大多数可靠 ECD 表面并帮助实验定位”作为已验证能力。当前最需要的是补齐映射、环境交付、执行闭环和独立验证，再决定哪些更复杂的物理通道值得保留。

审阅对象为当前工作树，包括未提交的 SES 扰动预览、实验工作台、FlowJo/FlowKit、Paraplume 与 GRASP 功能。Git 基点为 `9f01e12f0a7924b5b6253e4bb50ebc956cfc56f5`，不将当前工作树称为 clean HEAD。本次只新增审阅文档与审阅证据，没有修改算法、阈值、profile 或启动新的科学生产计算。既有四例运行保持原状。

已阅读根 AGENTS.md、project_memory、PD1 AF3 修复审计及当前物理实现/四例运行审计；审查了输入检索与映射、结构可靠性、SAS/SES、候选搜索、稳定性门控、物理校准、覆盖优化、worker/缓存、API/任务、实验统计、序列/导出、网页与测试配置。覆盖方式为架构和关键路径审查、最小反例、构建与回归验证，不等于对每条代码路径的形式化证明。没有重新运行真实模型权重、四例完整物理计算或浏览器视觉验收。

**1. 应当从什么科学问题出发**

真正的任务是受预算约束的实验设计：在表达、折叠和实验可执行性足够可靠的前提下，用少量突变构建体获得对未知结合区域有用的观测。至少需要区分下面五层。

| 层次 | 可回答的问题 | 不能直接推出的结论 |
|---|---|---|
| 结构与坐标 | 结构表示的是哪条序列、哪种构象和上下文？ | 预测坐标正确、膜方位唯一 |
| SASA / SES | 水尺度暴露、分子边界和连通路径如何？ | Fab 或双细胞受体真实可达 |
| WT/Mut 配对场 | 给定建模协议下哪里发生可解析变化？ | 抗体结合必然下降 |
| 探针覆盖与区分 | 预设足迹中心有多少被触及，探针模式是否不同？ | 实验检出率或残基级定位精度 |
| 实验观测 | 特定体系的 reporter/结合/表达发生怎样的变化？ | 未控制混杂时的直接表位归因 |

建议将目标写成分层合同：固定可靠目标表面 S 和面积权重 a，面板 P 的计算扫描覆盖为 `C(P)=Σ[a(x)·I(存在合格探针覆盖x)]/Σa(x)`。未知量不能作为已解析阴性。优化优先级建议为：硬资格 → 主预算下扫描覆盖 → 最大连通盲区 → 噪声下可辨识性 → 成本和风险。每一层分别输出，不混成一个没有单位的“综合可信分数”。

定位与检出必须分别验证。例如所有位置都被至少一个探针触及，仍可能出现大片区域具有同一响应模式；增加面积覆盖并不保证可以区分这些区域。应报告相同/近似签名的最大连通面积及混淆组，而不只报告平均签名距离。

FoldX 的多次 BuildModel 是算法收敛/局部解检查，不能视为独立生物学重复。官方说明支持这一用途：[numberOfRuns](https://foldxsuite.crg.eu/parameter/numberOfRuns)。配对参考 WT 与突变体的能量差及结构产物也应按官方配对关系解释：[BuildModel](https://foldxsuite.crg.eu/command/BuildModel)。

**2. 应保留的已有设计**

- author chain、auth_seq_id、insertion code 的导出映射与 AF3 confidence provenance；显示分区连通和直径约束；区域数、科学假设数、构建体预算分离。
- native SASA、稳定性模型、化学扰动和实验读出的证据分层；缺失评分不填假数；无合格面板可以是正确结果。
- 配对 WT/Mut、共同刚体变换、原单位通道、unknown mask、null 与数值敏感性分开记录；禁止补洞制造闭合 SES。
- 保护位置保留在主要覆盖分母中；同四次模型同时满足直接支撑与阈值；真实组合重算；候选池内求解上下界明确限定范围。
- 实验工作台的不可变来源快照、版本/QC 审核、按孔与独立批次计算、不把细胞事件当独立重复；组合效应不直接平摊到残基。
- 保存原始模型、置信度、DDG 表和 ZIP hash 的方向是正确的。缺少完整验收不能通过放松这些约束解决。

**3. 已确认的问题与具体修改建议**

优先级 P1 表示影响核心正确性、复现或完整执行，P2 表示局部正确性/解释/维护问题。以下既有最小反例，也有可直接确认的调用链缺口；影响范围分别限定。

| 编号 | 优先级 | 问题 | 主要位置 |
|---|---|---|---|
| R1 | P1 | wheel 漏新增包、脚本和网页资源 | pyproject.toml:21、26；joint_design_v4/engine.py:804 |
| R2 | P1 | 多个最优序列比对直接选择第一个 | structure_analysis.py:280–309 |
| R3 | P1 | 新物理原子读取丢弃糖链等 HETATM，上下文语义不统一 | surface_perturbation/paired.py:21–31；surface_patch/backend/sr_pointcloud.py:41 |
| R4 | P1 | 普通 API 预览入口没有完整校准/生产协调流程 | api.py:65–99、419；surface_perturbation/engine.py:169–180；scripts/continue_surface_validation.py |
| R5 | P1 | 通用 worker 超时缺少后代进程生命周期保证 | integrations/process_runner.py:19–38；surface_perturbation/tools.py:39–47 |
| R6 | P1 | Paraplume 上游权重实际内容不进入完整验证 | integrations/paraplume_service.py:34–95 |
| R7 | P1 | 成功后的基线对照被限制到物理方法已经选中的面板 | surface_perturbation/engine.py:342、531；surface_perturbation/benchmark.py |
| R8 | P1 | 测试和工作区没有完全隔离，CI 当前也不通过 | api.py:54–60；tests/test_api.py:7、77；.github/workflows/quality.yml:17 |
| R9 | P2 | unreviewed 条目被标为 reviewed | retrieval.py:105 |
| R10 | P2 | 离线报告的标题/品牌文本没有 HTML 转义 | exporting.py:65–75 |
| R11 | P2 | API、CLI、profile 的预算/策略入口不一致 | cli.py:59–81；api.py:77–85；joint_design_v4/budget.py |

R1：在临时源码副本构建 wheel，产物共 111 个 ZIP 成员，没有 `surface_perturbation/`、`scripts/`、`web/`。移除 editable-install 的回退查找器后，导入新包实际得到 ModuleNotFoundError。默认 V4 也无条件导入 STRATEGY_ID，因此不只是预览用户受到影响；网页服务还依赖缺失的静态资源。建议改为明确覆盖全部运行时包的 discovery/清单，将运行时审计函数移出 scripts，并打包所需 web/schema 数据。验收：在源码目录之外安装 wheel，运行默认 V4 的禁用外部模型 fixture、导入物理模块、访问首页/静态资源及离线导出；不能只测 `pip install -e .`。

R2：`ACDEFGHIKLACDEFGHIKL` 对 `ACDEFGHIKL` 有两个等分最优映射，当前函数静默返回后半段 11–20，identity 为 1.0。后续 WT 字母核验无法发现重复序列中的等价位置错误。该反例证明存在歧义处理缺口，不证明现有四例已经错位。建议优先使用逐残基 SIFTS/mmCIF 对应关系，核验与坐标/序列一致；fallback alignment 必须检查最优或近最优映射对目标位点的影响，无法唯一定位则标记 `MAPPING_AMBIGUOUS`。验收包含重复域、内部缺口、插入码、工程替换及链顺序扰动，不能只验证字母匹配。

R3：标准蛋白原子与 NAG 原子各一个的最小 PDB，Bio.PDB 读取 2 个原子，当前 `atoms()` 只保留 1 个。此过滤用于 native 物理参考及配对场；旧 SAS 点云却允许结构上下文原子遮挡。因而有糖链的输入可能在候选生成与 SES 计量阶段使用不同空间环境。rSASA 本来明确在孤立选中蛋白链计算，这本身不是错误；错误风险来自没有将不同上下文作为一等合同一致传递。建议区分 target、occluding proteins、resolved glycans、waters/ions 等角色，按明确配置选择遮挡体，保留上下文 hash；糖链可以只提供几何遮挡，未知化学类型继续弃权，不能直接把蛋白供受体规则套到糖上。验收比较同一结构含/不含糖链，校验缓存身份、遮挡和分母的预期变化。不要将已有 glycan-shield 注释当作真实糖链几何已经纳入 SES 的证明。

R4：普通新作业默认 `phase=production`，校准目录默认为新任务自己的空目录；JobCreate/上传表单没有校准来源或协调器参数，JobRunner 只运行一次 pipeline。现有完整扩池闭环主要在实验脚本中。因而用户在网页选择物理预览，并没有获得与四例验证脚本相同的闭环服务；此处按调用链确认，没有为复现而启动昂贵生产作业。建议把 `PRECHECK → CALIBRATING → VERIFYING → EXPANDING → EXPORTING → COMPLETE/INCOMPLETE` 做成可恢复的应用任务流程，共享现有科学引擎和累计预算；GUI 提供实际缺口与恢复入口。校准仍不足时必须保留未完成状态。验收：空缓存 API 作业必须先获得匹配校准，或明确返回可操作的前提缺失；不得把没有面板的通用 COMPLETE 文案当作设计已成功。

R5：已用与当前代码相同的 `subprocess.run(timeout=...)` 模式做本机无害反例：父进程在 0.5 秒超时，子进程仍在 1.5 秒后完成写标记。它证明直接子进程超时不等于进程树终止，不证明当前正在运行的队列已经存在孤儿。独立 benchmark 的 phase guard 已做 owned-tree 管理，但常规 API/worker 路径并未统一采用它。建议抽出统一执行器，Windows 使用 Job Object 或受控进程树，Linux/WSL 使用明确进程组/会话及取消协议；失败也保存 stderr/请求/退出信息。验收检查 timeout、取消、服务关闭和监督器异常后均无本任务后代存活，且不触及其他用户进程。

R6：服务检查 Hugging Face snapshot 文件存在与 size，但返回 fingerprint 使用清单里声明的 `commit/sha256_or_blob`。临时假环境中把 `AAAA` 改为同大小的 `BBBB`，两次 fingerprint 完全相同。模型主 checkpoint 的实际 hash 已计算，上游编码器权重仍可能漂移，不能据此声称“全部模型内容已核验”。建议核对实际 revision，并计算实际文件 hash；worker 源码、预处理/运行配置及关键版本也应进入预测缓存身份。可用只读内容寻址和校验记录降低重复读取成本。验收：同大小改写、refs 指向另一 snapshot、worker 行为变更均应使缓存失效或拒绝运行。

R7：初始 `controls(measured, ...)` 在最终面板产生后被 `controls(final, ...)` 覆盖。这样几何/SASA 基线只能从物理方法已经挑出的少数构建体中选，不是同候选池独立选面板的公平比较。它可用于描述最终面板的代理覆盖，但不能证明复杂物理路线优于基线。建议保留两类不同结果：`selected_panel_descriptive_controls` 和冻结全池下各方法独立选择的 `same_pool_method_comparison`。候选不足或未计算的物理数据保持未测；独立路线另做等预算比较，不能给未知候选填 0。验收使用一个“基线最优候选不在物理最终面板”的反例，确保公平对照不会漏掉它。

R8：API 在 import 时构造 JobStore/ExperimentStore，路径固定在源码 `workspace_data`；部分测试只 mock pipeline、直接向该全局 STORE 创建任务。初次测试确实新增了模拟任务 `ebc56ba63ff7`，发现后停止本次测试进程，改用临时源码副本完成后续验证；该模拟任务保留，未清理既有历史。建议 `create_app(settings, store, runner)`、环境可配置数据目录、pytest fixture 强制临时 store/runner，禁止测试接触日常工作区。Ruff 本轮实际为 1044 项，主要为 E701/E702/I001；这是发布检查失败，不等于 1044 个算法缺陷。验收先让既定 CI 通过，再加入 wheel 安装测试和 Windows 生命周期/文件锁测试。

R9：`"reviewed" in "unreviewed"` 为真，真实函数反例已确认。建议解析完整 entryType，至少明确区分 Swiss-Prot/TrEMBL/未知值；加 unreviewed 与缺失字段回归。影响是来源审核标签失真；没有证据表明它目前直接改变稳定性评分。

R10：给蛋白名称 `<b id="review_marker">text</b>`，HTMLParser 在导出报告中读到两个真实 b 元素，而不是文字，证明原始元数据被解释为 HTML。建议标题/品牌用 html.escape，所有 script 内 JSON 统一用 `_json_inline`；目前 PDB 文本单独使用 json.dumps，也应统一。后者的实际外部可控范围受上游结构重写限制，不能据静态代码夸大为已验证攻击链。验收覆盖尖括号、引号、`</script>` 和占位符样式文本，不访问外网、不执行注入内容。

R11：CLI 未提供 surface strategy 和新的 min/target/max/adaptive 参数，仍默认 `--max-primary-constructs 7`；网页采用另一套默认合同。相同输入通过不同入口，可能得到不同预算/算法，影响可复现性。建议共用请求 schema 与解析器，并导出“用户值 → profile 值 → 最终生效值”；添加 API/CLI 等价请求测试。不要在修复时追溯改写历史 profile。

**4. 方法学修改建议：证据不足与程序错误分开处理**

| 模块 | 当前判断 | 建议及验收 |
|---|---|---|
| 默认设计目标 | legacy V4 greedy 仍将 coverage、平均区分、风险、冗余合为净收益并早停；与广覆盖目标不完全一致 | 保留旧版本解释，新增明确覆盖优先 policy；输出 K=1…8/12 的可行性、下界、候选池内上界和盲区；新默认必须经验证再切换 |
| 候选空间 | 有按 order 独立搜索、全局去重，但每位点替换菜单、位点池、beam/shortlist 仍可能提前排除有价值解 | 扩池分别诊断“无可用替换、未搜索、未评分、未校准、被门控”；测候选池收敛，不能用当前池的上界证明全局不可行 |
| SAS/SES 与图 | 新 SES 已实现，早期“只有 SAS”的记忆已过时；两种分母与半径集不同 | 分开命名和存储；在同 native 输入上做参考面积、图路径、稠密真实 hard 场收敛，不能用固定锚点的面积/图代理一致代替全部物理场一致 |
| 坐标可靠性 | pLDDT 与 PAE 的角色已部分解耦；单构象仍限制长程解释 | 高 PAE 的域间方向做场景敏感性；不给高局部 pLDDT 自动授予域间几何可靠性；低置信度未知面积单列 |
| 化学场 | donor/acceptor/芳香方向和核是可追溯描述符，并非接触自由能 | 单独验证质子化、His、方向几何、半径/带宽敏感性；复杂通道只有在等池消融有增益时才进入主路线 |
| null 与阈值 | Q95/Q99、4/5 和数值敏感性是经验计算规则；相邻 WT 差分相关 | 按可恢复执行错误、查询缺口、几何拓扑、质子化不确定、有限重复不收敛分类；不可解析问题及早标记，不能靠机械追加20次解决；不要赋予统计FPR含义 |
| FoldX/SPURS | 适合稳定性风险预筛，非表达/正确折叠保证；模型共识不是两份独立实验 | 记录适用域与分歧，三/四点外推单列。SPURS 原文 multi 模型以单/双突变数据训练，不能把更高阶实测可靠性视为已有：[原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12830971/) |
| 覆盖/定位 | 12 Å 扫描规则是场景假设；一个直接变化点可覆盖多个中心 | 同时报 one-site/two-site、半径/权重场景、签名混淆和定位分辨率；这些曲线仍不是未知抗体的保证 |
| APBS | 已拆分执行、准备、数值与 extended；这一点正确 | 明确“hard 合格但电势未知”；不要将同一 `VERIFIED_PANEL` 标签理解为全部电势已解析。若调整主面板与12项对照的成本合同，应前瞻版本化，不能追溯退款/变更旧验收 |
| 物理验证 | 当前四例未提供最终验收，不能证明复杂模型有实际收益 | 先小规模冻结对照：同候选池、同K、同资格、同输入，分别比较距离/SASA、原子占据、化学、SES、可选电势；再做独立候选路线比较 |
| 模型盲测 | 隐藏复合物接触标签有价值，但接触几何并非功能检出真值 | 按蛋白/结构域/抗体家族 holdout；避免同抗原上的多个抗体都算独立泛化证据；受体/配体案例独立汇报 |
| 实验统计 | 当前用孔中位数和独立批次 bootstrap，方向合理 | 明确技术孔、生物学重复、日期/批次层级；小批次数显示原始点与不确定性，不给精密显著性印象；预先规定主要终点与多重比较策略 |
| 实验归因 | reporter 下降包含结合、表达、膜接触/方位、信号阈值等路径 | 为“可导出定位证据”增加表达/构象控制与复测的结构化状态；阴性不能自动排除表位。抗体逃逸研究分别控制表达与功能，说明此类混杂需要独立处理：[原始研究](https://www.nature.com/articles/s41467-021-24435-8) |
| GRASP 导出 | 目前单点来源、映射与人工 rationale 被校验；一个批次且非空 decline 就可选择 experimental_ir | 将“实验关联线索”与“可作为直接接触约束”区分。记录效应方向、大小、不确定性、表达/构象控制；人工强制假设可保留，但不能与经过验证的接触证据同等级 |
| 序列/CDS | 有 WT/位置、翻译核验，保留；密码子优化可退化到确定性选码 | fallback 附带异常原因、已满足/未验证约束清单；每个设计保存完整构建背景、标签/锚定/信号肽策略，避免ECD序列与实验全构建体混淆 |
| UI/报告 | 现有 unknown 和预览警示应保留；静态字符串测试不能证明三维实际着色 | 顶部同时显示执行状态、hard证据、electrostatics和实验验证状态；无已验收面板时正式覆盖显示不可用。补真实浏览器 auth/label 选择与下载产物一致性验证 |

AlphaFold 的局部置信度和相对位置置信度必须分别理解，见 [EMBL-EBI 官方 AF3 解释](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/how-to-assess-the-quality-of-alphafold-3-predictions/)。上述 ensemble/场景建议是基于这些证据边界提出的设计建议，不声称当前已经实现。

**5. 架构与交付建议**

建议继续保持单机研究工作台架构，不急于拆微服务。先按职责划开：纯科学计算、证据状态与持久化、外部工具执行、Web/API 表示层。`api.py` 约 2.6k 行、`experiment_projects.py` 约 4.7k 行，包含大量事务、上传、QC、统计与导出责任；物理 engine 也高度嵌套。拆分应以稳定接口和测试约束为前提，而不是整体重写。

推荐统一四种合同：

1. `StructureContext`：canonical/author/label/normalized 映射、target/occluders、原结构和 confidence hash。
2. `EvidenceArtifact`：真实输入、工具/权重/预处理身份、原始值、单位、未知理由、运行身份与依赖图。
3. `Campaign/Attempt`：不可变 attempt_id、阶段/预算/取消/恢复、运行源代码身份，避免仅依赖文件mtime判断当前结果。
4. `PanelAssessment`：资格、hard/extended、分母、盲区、可辨识性和实验验证分别表达，由 UI/CSV/JSON/ZIP 共享。

源码 hash 不是科学有效性的证明，但对缓存身份至关重要。当前部分数值身份将源码字节串联，格式化也会整体失效；另外一些模型身份又未核实实际权重。应先保证依赖闭包完整，再考虑版本化语义合同/可审计 raw 缓存迁移；不能以减少重算为由复用旧阈值或旧资格。

保留当前累计预算合同，同时分别记录排队时间、active solver 时间、wall time、休眠/中断保守计费、缓存命中和失败成本。性能结论用这些分量与同输入对照，不能把长壁钟间隔直接当 solver 性能。先依据 profile 优化重复序列化、重复 scope 查询和无信息调度，再考虑GPU或更多并发；本机已有运行时尤其不应盲目增加并发。

**6. 建议实施顺序与完成标准**

| 阶段 | 工作 | 完成标准 |
|---|---|---|
| A：可信交付 | R1/R2/R5/R6/R8/R9/R10；保留用户改动并按依赖做小提交 | wheel 外目录烟测、所有必需回归、Ruff；映射/哈希/子进程反例全部变为受控行为 |
| B：统一结构与任务合同 | R3/R4/R11；统一输入上下文、应用校准闭环及状态展示 | 空缓存 API 全流程可恢复；原AF3/confidence与author映射守恒；表/图/JSON/ZIP一致；失败状态准确 |
| C：证明方法增益 | R7与等池对照、候选池/密集场收敛、四例合规验收 | 完整发布每例结果或未完成原因；覆盖和上下界可重放；复杂通道相对廉价基线的收益/成本明确 |
| D：实验闭环 | 单点/组合复测、表达/构象控制、实验定位终点与独立holdout | 有可追溯真实实验效果；检出和定位误差分别评价，才考虑校准概率或更强模型 |

不建议现在同时加黑箱 binding predictor、全蛋白 MD 或扩大组合阶数。它们未解决映射、流程完成率和独立验收的前提问题；是否加入应由阶段 C/D 的具体误差归因决定。也不建议直接降低稳定性/校准门槛来填满8项。

**7. 当前运行与本次验证证据**

本次只读状态快照见 [review_case_snapshot.json](project_review_20261003/review_case_snapshot.json)。快照时 CTLA4 校准不足，CTLA4/PD1/PDL1 生产均为 `CALIBRATION_EXPANSION_REQUIRED`；PD1/PDL1 校准正在运行，B7-H6 新批次尚未开始。四例 `accepted_panel=false`。这些是有时间戳的观察，既有队列后续可以继续变化；不是对最终性能的结论。

最小反例见 [review_repros.json](project_review_20261003/review_repros.json)，wheel 导入失败见 [wheel_import.txt](project_review_20261003/wheel_import.txt)。所有假权重、子进程与 HTML 反例在临时目录内完成，没有修改真实权重或运行历史。

全套测试结果为 **391 passed、8 failed、11 skipped、45 warnings，耗时 772.19 秒**，详见 [review_pytest.log](project_review_20261003/review_pytest.log)。测试从临时源码副本的父目录、使用原项目 venv 执行，副本含 271 个源码/必要资源/fixture 文件，未复制外部模型环境或历史 workspace。命令为 `python -m pytest -p no:cacheprovider ecd_epitope_scanner/tests -q --tb=short`。

8 项失败应分别解释，不能全部算作科学算法缺陷，也不能忽略后宣称测试全绿：

- 2 项 assay profile 测试仍断言最大构建数为 6，当前默认 profile 给出 8。需要按 profile 版本明确行为与回归，不能为通过旧断言直接把新默认改回 6。
- 4 项实验 API 集成测试返回 `FLOWKIT_WORKER_UNAVAILABLE`（503）；隔离副本未带 `.tool_envs`，因此这些路径未获得运行验证。建议区分有固定工具环境的集成测试与使用明确替身的 API 合同测试，并注册相应 marker。
- 2 项工作表相关断言涉及空文本重载后为 `None` / `inlineStr`，而测试要求空串 / `s`。当前 openpyxl 为 3.1.2；额外直接核验 `parse(build(rows)) == rows` 通过，见 [worksheet_probe.json](project_review_20261003/worksheet_probe.json)。这支持优先修正空单元格表示的合同与测试，而非据此断定数据丢失；完整 API 场景仍需修复断言后重跑。

Ruff 的 1044 项分布为 E702 516、E701 377、I001 118、E401 20、F401 7、E731 3、F841 2、UP017 1。现有 CI 先执行 Ruff，因此当前工作树未达到发布检查要求。此次没有批量格式化或自动修复，以免混入无关用户改动。

本审阅仍未覆盖的验证：真实 FlowKit 固定环境/FCS 全部样例、Paraplume 真实推理、全新浏览器视觉检查、新四例最终科学面板及实验检出/定位性能。后续应按验收表补齐，不能把此次代码检查或单元回归写成这些项目已通过。
