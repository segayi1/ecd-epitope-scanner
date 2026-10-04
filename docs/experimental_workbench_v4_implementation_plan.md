# v4 实验工作区实施计划

## 1. 目标、范围与实现基线

本计划把 FCS 原始数据、构建体、门控、质控、效应统计、v4 结构定位、抗体序列证据和 GRASP 输入导出接成一个可恢复、可审计的工作流。一个已完成的 v4 结构作业可建立一个冻结实验项目；项目内可累计多个实验批次和多个不可变分析版本。

实现遵循以下边界：

- 项目固定抗原参考序列、结构文件及 SHA-256、`residue_mapping.csv`、v4 sector 文件和 `surface_graph_v4.json`；用新结构/新分区重分析须显式创建项目新版本或克隆项目，不能改写已有证据坐标。
- 原始 FCS 只读保存；样本归属、门控、QC 决定、统计结果、人工审核和导出物分开版本化。
- 从 FCS 计算的统计单位是孔/样本，再到独立实验批次；绝不把事件数用作独立重复数。
- 缺少结构映射的构建体可以保留实验数据；只是不生成该位点的结构定位。
- 组合突变的功能效应不分摊成组内单点效应。结构渐淡只用于显示，不产生证据或 GRASP 约束。
- 单浓度阳性率下降描述给定门控、染色条件和浓度下的功能变化，不表示亲和力变化或直接接触概率。
- GRASP 本阶段是经过校验的文件接口；不安装、不运行 GRASP 模型。

计划编写时的代码基线（2026-09-25；现况以 1.4 节的当前开发交接为准）：v4 已提供结构上下文、权威 `residue_mapping.csv`、表面分区及后台 v4 作业；新增 `experimental_effects.py` 实现孔/批次效应、探针归属和可视化测地渐淡，新增 v4 表面图快照输出。当时尚无实验项目数据库、FCS 解析与门控界面、Paraplume 任务和 GRASP 导出界面。已实现的效应/映射边界见 [`experimental_surface_mapping_v4.md`](experimental_surface_mapping_v4.md)。

### 1.1 交接快照（2026-09-25）

本表记录 2026-09-25 初版计划形成时的基线，只作历史对照；当前实现状态以 1.4 节的当前开发交接为准。下文 API 表和 Phase 清单是目标设计，不能单凭接口表推断功能已经上线。

| 能力 | 当前状态 | 可核对的实现/证据 | 仍需开发 |
| --- | --- | --- | --- |
| v4 结构作业与分区 | 已有代码 | `joint_design_v4/engine.py`、`joint_design_v4/sectors.py` | 实验项目冻结、哈希核验、API 读取 |
| 表面图快照 | 已接入 v4 输出 | `joint_design_v4/surface_snapshot.py` 写 `surface_graph_v4.json`；`engine.py` 调用 | 项目复制、版本化读取、结构视图 |
| 阳性率效应与表面归属 | 独立计算内核已实现 | `experimental_effects.py`；7 个相关测试于 2026-09-25 通过 | FCS 计数输入、QC 审核、序列化、批次面板区间、产品 UI |
| FlowKit 1.3.2 | 独立环境已安装 | `.tool_envs/flowkit`、`tooling/requirements-flowkit.lock` | FCS/FlowJo worker、真实 fixture、门控编辑、补偿和 QC 流程 |
| Paraplume 2.0.0 `--large` | 环境、权重与单次 GPU 推理已验证 | `.tool_envs/paraplume`、`tooling/manifest.json`、`.tool_cache/paraplume_smoke/user_example_test_summary.json` | 受控 worker、任务 API、缓存、结果审核和页面；CPU 回退尚未实测 |
| GRASP-JAX 输入 | 仅格式计划 | 本文第 8 节及官方格式链接 | FASTA/IR/RPR 校验器、导出清单与 API；模型按范围要求不安装 |
| 完整实验工作区 | 未实现 | 当前 `api.py` 仍只有已有 `/api/jobs` 路由 | Phase 0–7 的数据库、接口、页面和重启验收 |

已验证的 Paraplume 结果是**抗体侧**逐残基 paratope 概率。示例序列与已发表的 FMC63 VH/VL 相同，CD19 的表位信息来自独立的 [FMC63–CD19 复合物结构](https://www.rcsb.org/structure/7URV)和[界面分析](https://pmc.ncbi.nlm.nih.gov/articles/PMC10228544/)；不能把这些文献中的抗原位点写成 Paraplume 的预测输出。产品中的抗原定位仍须走本文第 6 节的实验效应与结构映射流程。

### 1.2 原计划的开发起点（2026-09-25）

先完成 Phase 0 的最小闭环：从一个**真实已完成**的 v4 job 创建实验项目，复制并校验参考序列、结构、`residue_mapping.csv`、sector 和 `surface_graph_v4.json`，写入 SQLite 项目/修订/artifact 记录；提供项目创建、列表、详情和按 artifact ID 下载四类接口。重启服务后读取同一项目并复核所有 SHA-256，作为首个可交付的集成验收。此步不依赖 FCS fixture 或 Paraplume worker。

随后按 Phase 1→4 接入构建体、FCS、门控、QC 和结构图；Phase 5 把已验证的 Paraplume CLI 封装成后台任务；Phase 6 实现 GRASP 文件校验与导出；最后按 Phase 7 做全链路重启复现。每个 Phase 的退出条件在第 11 节，逐项需求的验收证据在第 12 节。

### 1.3 开发续记（2026-09-26）

以下是本轮工作树状态；它覆盖上方 2026-09-25 的历史交接快照，Phase 退出条件仍以第 11 节为准。

| Phase | 当前状态 | 已核验内容 | 仍待完成 |
| --- | --- | --- | --- |
| 0 冻结项目 | 已实现并对真实 v4 job 验收 | 项目快照、artifact 哈希/只读登记、SQLite 迁移、详情/列表/下载与重启读取；项目快照、结构证据、FCS/WSP/结果表原件、补偿矩阵、FlowKit 与结果表效应、Paraplume 结果及 GRASP ZIP 发布均有硬退出重启恢复测试；启动时从已提交 registry 修复派生 manifest，未提交记录不出现在项目查询中；Paraplume 任务中断标记为 `INTERRUPTED`，以新任务重试；矩阵按内容 SHA 寻址且重复上传幂等；跨进程与 API 层同 revision 批次写入均验证为一胜一冲突且只提交一个 batch | 已覆盖的硬崩溃会留下未登记文件供检查，不自动删除；其余 artifact 边界需随新增写入类型继续审计 |
| 1 构建体/样本 | 已实现 | 推荐/自定义抗原、参考 WT 序列、抗体 parent/variant、批次、FCS/WSP 导入/匹配 API 与页面；完整序列可先对冻结 WT 做无写入差异预览，成功后才可登记；批次表单包含剂量/单位、阳性门、独立批次和备注；逐 FCS 匹配可选抗原/抗体构建体、重复类型/编号，零选择不会伪装成确认；后端拒绝未关联双侧构建体的 ACTIVATION 样本；ACTIVATION 批次返回并显示缺 activation/NC/WT/构建体关联的设计警告；历史 antigen/Fab Myc 可作为独立 `HISTORICAL_SUMMARY_ONLY_NOT_CONCURRENT_CONTROL` 记录录入，强制构建体、日期、指标、来源和审核，且不进入 FCS/activation 分析；实际 H1–H3 标签和 H4=NC 已登记；Phase 1 退出条件均有 API/重启/浏览器测试 | — |
| 2 FCS/FlowJo/门控 | 部分实现 | FlowKit 1.3.2 worker；260903 三个独立批次各 54 个 FCS，逐批完成 162 样本/486 条门重算；通过 `$FIL/$WELLID/$PLATENAME/$SMNO` 唯一配对 WSP 样本，并只在临时副本中重命名重复样本键，避免 FlowKit 折叠 WSP 内同名孔；缓存计数差异作为可见诊断，不阻断同批门控方法一致的 FlowKit 计数；逐批检查门层级、原始门型、FlowKit 门型、所用通道变换、WSP 补偿以及门成员语义（补集和 Boolean 运算/引用）；FlowKit 支持的 WSP BooleanGate 可预览和执行；FlowKit 将 WSP EllipsoidGate 近似为 PolygonGate，预览会明确标注近似；未知、参数不完整或锁定版不支持的通道变换会在 workspace 构造前标为 `WAITING_FOR_USER`；FCS/WSP 原始哈希保持不变；跨批重复 `$FIL` 拒绝合并导入；不一致方法和未匹配样本仍阻断计数；手动门控重建入口保留；不可变 gate revision、版本化变换、层级执行、抽样、外部矩阵登记、时间窗排除和 provenance 有实现/测试 | 用户已确认 260903 采集后未设置补偿；WSP 无补偿矩阵且 FCS `$SPILLOVER` 为 10×10 identity，与未补偿状态一致。按用户选择，以已圈门结果表作为主路径，不处理或重算 FCS；FlowJo/WSP 更多门型与变换覆盖属于可选的方法复核，不是当前结果表流程的前置条件；同批方法检查不比较逐样本几何边界，但会比较补集/Boolean 成员逻辑 |
| 3 QC/效应 | 部分实现 | 采集率/边界/非有限值/低事件 QC；手工 QC 版本；effect-analysis API/UI 从 gate parent count 计算孔/批次效应并写入只读 artifact；单点、抗体效应和双侧交互模式均保留逐批数据、QC 与技术重复汇总；已圈门结果表效应保存逐孔原值、人工 QC 状态、比较组和是否纳入摘要；UI 提供变体/WT/NC 分层逐孔散点与明细表，固定 0–100% 坐标，不显示推断区间；用户已确认 260903 标签：CTLA4-WT、E1–E7（抗原构建体临时标签）、NC（仅抗体-SNIPR 阴性对照）、L32/L49/Y159/W202/W212（抗体）、CD80（阳性配体）；三个目录名定义三个独立批次，每批内两个试剂组各 9 条件×3 技术重复，condition map 可从首份 Excel 的 Sheet1 读取 | 尚无已审核的 260903 QC 决定或正式效应；每种抗体/配体只出现在其中一个独立批次，不能把三批当成同一抗体的重复或跨批推断；已确认数据未补偿，事件分母未知，结果表百分比只能按导出观测解释；首份 Excel 更新后 W212 首个技术重复已与 table-3 对应，不再有已知冲突；E1–E7 的具体突变清单未知，暂不做结构位点着色 |
| 4 结构效应图 | 部分实现 | 冻结结构上下文校验完整性并返回 structure/mapping/surface graph/hash；保留 direct/sector/single/fade 图层、组合效应不拆分和仅可视化渐淡；来源 run manifest 逐文件 SHA 核验后将 238 位点可靠性、SASA 和 QC 只读冻结到真实项目 revision 8；可靠性与 SASA 分类图层及残基详情、历史 Myc、人工 QC 豁免和 UNMAPPED 独立状态卡；VH/VL 序列突变轨道与抗原 3D 坐标分离；当前图层 PNG 与包含映射、graph、效应/构建体、Myc/QC 和限制说明的证据 JSON 导出；浏览器 smoke 与 API/重启用例通过 | 用具备真实 activation/QC 结果的数据核验效应图层着色和逐孔回溯；更丰富的复合物结构视图（有可用结构时） |
| 5 Paraplume | 部分实现 | 隔离 worker 强制 Paraplume 2.0.0 `large` 并禁止网络缺权重时降级；项目序列/模型指纹缓存、异步 task/SSE、失败码和重启中断状态；预测 artifact/manifest 哈希及 `PREDICTED_UNREVIEWED` 标签；真实用户 VH/VL 经 RTX 3060 GPU 和隔离 worker CPU 分别执行完整模型成功，均返回 120/107 个有限 `[0,1]` 概率；CPU 与 GPU 结果最大绝对差 `2.15e-6`；CPU worker provenance 标为 `CUDA_UNAVAILABLE`；CUDA OOM 注入→CPU 重试单测通过；新增 5 个完整模型 embedding 步骤的实际进度，经任务记录、SSE 和页面显示，候选仍标为 PREDICTED_UNREVIEWED 且导出要求审核人/理由 | 硬件真实 CUDA OOM 触发尚未实测；新增进度链路由 worker 单测和浏览器模拟 SSE 验证，本轮尚未通过真实模型任务端到端验证；Paraplume 许可条件需正式部署前确认 |
| 6 GRASP | 部分实现 | 按官方 `generate_restr.py` 固定 commit `db1715f51039e27fcc8ec4ca05e3217cb0c4ea67` 锁定 1=ECD/2=VH/3=VL 与 IR/RPR 语法；上游 parser 小型 fixture round-trip 通过并锁定 SHA；本项目 parser 校验 1 起始位置和 FASTA 氨基酸；固定时间戳 ZIP 与独立清单；导出 API/下载、revision/audit/artifact hash；预测 IR 关联同一抗体 Paraplume run 并保留 `PREDICTED_UNREVIEWED`；实验 IR 只能来自 effect artifact 中有效 `SINGLE_SITE` decline，并经冻结 CSV `uniprot_pos→normalized_position` 映射；页面可预览/生成两种 IR；无复合物结构时 RPR 稳定拒绝并输出 `IR_ONLY`；serializer、存储/API 下载、错误码及 ZIP/manifest 硬退出恢复 fixture 通过；真实浏览器打开现有实验项目，并以 mock 候选走通两种 IR 的预览/提交 | 当前没有复合物结构/RPR 证据 |
| 7 全链路恢复 | 部分实现 | 快照、合成四格 activation→QC→effect→GRASP ZIP、重启恢复和可复现下载已有验收；浏览器覆盖未解析 FlowJo 手动重建、逐孔/技术重复和时间窗排除；真实 260921 Myc 与 260903 activation FCS/WSP 只读回放和哈希核验；260903 162 个样本逐板导入、同批方法一致性、WSP 缓存差异显示及跨板重复 `$FIL` 拒绝有真实 fixture；JobStore/Paraplume task 和项目 artifact/manifest 崩溃恢复已覆盖多处硬退出点；FCS 原件、补偿矩阵、FlowKit/结果表效应、Paraplume 结果和 GRASP ZIP 提交边界现有硬退出恢复测试；批次 revision 的跨进程与 API 并发冲突均有回归 | 260903 尚未完成项目登记→逐样本映射→人工 QC→效应的真实 activation 全链路；Paraplume 进度真实模型端到端、真实硬件 CUDA OOM 与部署许可仍待验收 |

当前真实项目 `4823bdbae91440bb` 保持 COMPLETE、revision 8，已冻结 17 个 artifacts；四个 H 样本中 H4 是 NC，H1–H3 使用临时突变体标签。该项目属于 Myc 染色，不提供 activation/WT 比较，不能生成 activation effect。真实 FlowJo QC flags 仍需实验人员逐样本审核；开发过程没有替代操作者填写 PASS/HOLD。

2026-09-27 在原验收基础上，实验项目集成用例扩展为合成四格 activation 全链路：effect 后重启、由有效单点下降生成实验 IR、冻结映射生成 `1-2-C`、ZIP 内证据回溯 effect artifact SHA，再次重启后下载字节与 SHA 不变；Paraplume task 重启用例也改用确定性 PID 存活注入，整组 Paraplume task/fallback 8 passed。全仓最终 `pytest --collect-only` 收集 169 项；由于当前执行器约 60 秒会清理单一长进程，按相同文件集合拆为 54、44、33、8、30 五个分片执行，合计 169 passed。真实 260921 数据仍是 Myc 染色，H1–H3 保持临时突变体标签、H4=NC，不作为 activation/WT 比较；真实 activation 全链路仍需后续配套数据。

同日继续补齐 Phase 1 的实验安排边界：ACTIVATION 批次由 API 计算 `READY/INCOMPLETE` 和缺 activation、NC、WT 构建体、WT 参考样本、构建体关联等警告；未关联抗原和抗体构建体的 ACTIVATION 样本在后端拒绝。页面批次表单增加剂量、单位、阳性门、独立批次和备注，逐文件匹配增加抗原/抗体、单孔/技术/生物重复及编号，并明确保留未匹配文件不能分析。合成闭环测试覆盖警告从不完整转为 READY。

Phase 1 随后增加数据库迁移 v12、`summary-measurements` API 和页面表单，历史 antigen/Fab Myc 汇总强制关联对应构建体、测量日期、指标/单位、来源标签和审核人，固定返回 `HISTORICAL_SUMMARY_ONLY_NOT_CONCURRENT_CONTROL`，不会创建样本或进入 activation 效应。比例越界会拒绝，记录随项目 revision/audit 和支持包保存，服务重启后仍可读取。完整序列新增无写入 preview API，页面在显示冻结 WT 差异、映射状态、风险和序列 SHA 后才启用登记；实验项目现为 14 passed，浏览器 smoke 扩展到 23 项并全部通过且无 JavaScript 错误。

Phase 2 随后补入登记样本的事件绘图抽样：FlowKit 从原始未变换事件中按确定性均匀索引取最多 10,000 对通道值，返回完整/有限/抽样事件数及 1–99 百分位显示范围；API 只从项目样本解析已登记 FCS artifact，并要求结果带 `display_only` 标记。工作台提供散点与 64×64 密度显示，切换图形复用同一抽样。真实 H1 FCS 两次抽样得到相同点和范围（13,863 个完整事件中显示 257 点），前后 SHA-256 相同；API/实验项目测试 14 passed、真实 FlowKit fixture 5 passed；真实项目浏览器 smoke 31 项通过且无 JS 错误。未解析 FlowJo 门新增显式人工重建路径：用户需提供完整内部门树、审核人和理由；修订与重算结果固定保存预览哈希及同一组 FCS artifact IDs，换用其他 FCS 会被拒绝；原始 WSP 预览仍保留 `WAITING_FOR_USER`，不会伪称 FlowJo 原生门已解析。登记状态由内层状态和未解析节点共同核验；内外状态不一致 fixture 与人工重建闭环均已通过实验项目测试。该页面路径用模拟的 EllipsoidGate 节点验收，尚无真实未支持门型样本。真实项目仍是 Myc 染色，不产生 activation 结论。绘图不用于门控计数或 QC。

2026-09-27 继续完成 Phase 2 外部补偿矩阵接入：API 接受 UTF-8 带通道行列标签的 CSV 和 JSON，校验矩阵维度、唯一通道名、有限数值与数值稳定性；矩阵以只读 artifact 登记并写入项目 revision、audit 和 manifest。`EXTERNAL_MATRIX` gate revision 只能引用同项目 `compensation_matrix` artifact；FlowKit 1.3.2 worker 使用 `flowkit.Matrix` 对全事件执行门控；分析结果固定保存 artifact ID、SHA-256 和 detector 列表，页面支持上传、选择与版本回看。API/门控集成与 parser 测试 29 passed；真实本机 H1 FCS 的 worker 外部矩阵执行 fixture、其他真实 FlowKit fixture 共 6 passed；浏览器 smoke 32 项通过且无 JS 错误，Ruff 与 JavaScript 语法检查通过。全仓 178 项测试按 48、17、71、20、22 项分组执行，178 passed。外部矩阵的 H1 fixture 使用显式标注的合成矩阵验证软件路径，不代表实验矩阵或实验结论。检查用户提供的 FlowJo WSP 时，H1–H4 均未返回 workspace compensation matrix；因此真实外部矩阵下的四样本验收仍等待实验矩阵文件。项目仍是 Myc 染色，H1–H3 是临时突变体标签，H4=NC，不产生 activation/WT 比较。

同日继续补入 Phase 4 当前结构视图 PNG 导出：3Dmol 导出用户选择的冻结结构图层，文件名含项目 ID、图层和冻结结构 SHA 前缀；浏览器实际下载 PNG 并确认内容非空，新增回归后 smoke 共 33 项通过且无 JavaScript 错误。该 PNG 是当前视图的下载件，不改写项目或其分析证据。VH/VL 突变轨道、结构可靠性/SASA/Myc/人工豁免徽标及效应证据图导出仍待完成。

同日继续补入 Phase 3 逐探针×独立批次热图：格内只显示后端计算的每探针/每批下降比例；无效对照窗口保留为“无效”，不插补，也不把孔作为独立批次。浏览器 fixture 覆盖有效与无效两个批次格；最新结构图与热图 smoke 共 34 项通过且无 JavaScript 错误。

2026-09-27 继续完成 Phase 4 结构证据联接：可靠性、SASA 和 QC 文件必须同时出现在来源 `run_manifest.files` 并逐一匹配 SHA-256，所有 UniProt 位置必须与冻结 mapping 完全一致，分类计数与逐残基记录核对后才能接入。真实项目 `4823bdbae91440bb` 的 238 个位置可靠性（预测模型；STRUCTURED 214、FLEXIBLE_STRUCTURED 7、UNCERTAIN 3、DISORDERED 14）及 SASA 比较（AGREE 234、SASA_DISAGREEMENT 4）已作为只读 artifacts 冻结到 revision 8；structure-context 仅从项目已登记 artifacts 读取，重启后验证可读。页面把可靠性、SASA、summary-only Myc、人工 QC/override、UNMAPPED 分开显示；新增可靠性/SASA 类别图层、点击残基细节、独立 VH/VL 序列突变轨道、当前视图 PNG 和结构效应证据 JSON（含 hash、mapping、graph、构建体、效应、Myc/QC 与限制）。实验项目测试 18 passed；全仓按五组 36、34、42、31、38 项执行，合计 181 passed；浏览器 smoke 验证徽标、图层、显式冻结操作、VH/VL 轨道、PNG 和 JSON 下载，无 JavaScript 错误。真实样本仍为 Myc 染色；H1–H3 是临时突变体标签、H4=NC，不据此推断 activation 或 WT 效应。真实 activation/QC 数据图层验收仍待配套样本。

2026-09-27 继续补入 Phase 3 双侧交互分析：`DOUBLE_MUTANT_INTERACTION` 要求一个已登记变体 Fab 与其亲本、一个抗原变体探针，并由真实样本构建 parent/variant Fab × mutant/WT antigen 对照；每个 Fab 下还需突变抗原 NC 和 WT NC。每个独立批次分别计算两种 Fab 的突变抗原−WT 阳性率百分点差，再输出差中之差、四格 activation/NC 原始比例、QC 旗标及 probe 级批次中位数；至少 3 个批次才给描述性 bootstrap 区间。即使 WT−NC 窗口无效，原始 DID 仍保留并单独加无效窗口旗标；该模式不生成位点效应或 antigen 表面着色。页面要求逐样本指派、显示样本抗原/Fab 身份，并用独立 DID 表替代普通 antigen-effect 热图。纯计算及 8 样本 API fixture、固定抗体拒绝路径通过；本轮 `test_api + test_experiment_projects + test_experimental_gates` 37 passed，效应计算/项目测试 26 passed，浏览器 smoke 验证 DID、独立批次区间缺失提示、无位点归属及无 JS 错误。全仓当前收集 184 项；其余未改动分组沿用本轮前刚通过的分片结果。真实 activation 与 Fab×antigen 四格数据仍待用户提供。

2026-09-27 继续完成 Phase 3 逐孔图与技术重复汇总：逐孔阳性比例图按探针、独立批次和抗体构建体分面，逐点保留样本 ID、阳性/父门计数和 QC 审核状态；技术重复中位数/范围仅用 `PASS` 与 `PASS_WITH_OVERRIDE` 孔计算，`HOLD/EXCLUDE` 数量单独列出，不与生物重复混合。浏览器 fixture 特意放入 HOLD 技术孔，确认图上保留但不进入汇总；数值中位数和范围断言、DID 和结构证据路径均通过，浏览器 smoke 无 JavaScript 错误。当前全仓收集 184 项；本增量仅改前端显示、夹具和计划文档，Python 分片回归证据见上文。

2026-09-27 继续补齐 Phase 2/3 人工采集时间窗排除：操作人只能从来源分析的完整 QC 窗口中选择，必须提交审核人和理由；门控重跑保持原 gate revision 与 FCS artifact 集不变，按 FlowKit membership 对所有门的 count/parent_count 重新计数，重新计算 QC，并新建只读 gate-analysis artifact。manifest/audit 记录来源 analysis artifact SHA、选择的事件区间、累计原始索引区间、审核信息和排除前/后事件数；连续多次排除会将当前过滤序列索引映射回原始事件索引，原始 FCS 不改写。API 测试覆盖不完整窗口拒绝及连续三次排除的区间映射；真实 H1 worker 回放通过，原文件 SHA 保持不变；浏览器验证标记窗口、审核表单、提交和继承记录显示，无 JS 错误。全仓当前收集 186 项，按 38、9、25、85、29 五组执行，186 passed；工作台浏览器 smoke 与 Ruff/JavaScript 语法检查均通过。

2026-09-27 继续完善 Phase 3 单探针跨批次汇总：后端的原始百分点差中位数包含全部纳入批次；下降中位数、每批方向计数和 bootstrap 95% 区间只按有效 WT−NC 批次计算，且至少 3 个有效独立批次才提供区间。无效 WT−NC 窗口计入总批次数但不进入下降统计。工作台为 antigen-effect 展示该汇总表，与逐批效应表和热图并列；旧版 artifacts 缺方向计数时显示“旧版未记录”。效应纯函数测试及 effect→GRASP 集成测试共 8 passed；浏览器断言中位数区间/方向数量、双侧交互仍使用独立摘要表，整个工作区 smoke 无 JS 错误。此前 186 项全仓分片回归保持通过，本次增量的相关 Python/UI 测试已复跑。

同日继续补齐 Phase 3 每批对照窗口面板：后端从同一批已 QC 样本计算突变/WT 抗原率、各自 NC 率，以及 antigen−NC 百分点窗口；页面逐批展示四种原始阳性率、两个校正窗口和 QC flags，与效应/下降表并列。无效 WT−NC 批次仍显示其原始对照率与窗口，但保留无效旗标；不将孔或事件并为独立重复。新增计算断言、effect→GRASP 集成测试共 8 passed，浏览器确认面板行数、率和窗口显示，工作区 smoke 无 JS 错误。

2026-09-27 收尾回归：全仓 186 项按 38、9、25、85、29 五组均通过；真实 FlowKit H1 时间窗排除、浏览器全工作区 smoke、Ruff、两份 JavaScript 语法检查通过，`git diff --check` 无空白错误。真实项目 `4823bdbae91440bb` 仍为 COMPLETE revision 8；现有 FCS/WSP 是 Myc 染色，H1–H3 仅用临时突变体标签、H4=NC。真实 activation/WT 效应及真实实验外部补偿矩阵仍需相应配套数据，因此不能用这组 Myc 文件替代 Phase 2/3 的真实 activation 验收。

实现时复用现有 `api.py` 的 FastAPI 路由/文件返回模式、`jobs.py` 的 SQLite 作业与重启处理模式，以及 `joint_design_v4/surface_snapshot.py` 的冻结图格式。新增实验数据表和迁移独立于既有 job 表；文件写入先临时落盘、校验哈希、原子改名，再提交完成状态。首个集成测试须确认原有 `/api/jobs` 行为不变。详细接口契约见第 4 节。

### 1.4 当前开发交接（2026-09-29，覆盖前述历史状态）

**下一位开发者先读本节。** 本节结合当前源码、`workspace_data/experiments.sqlite3` 和既有验收记录核实；1.1–1.3 及文末续记保留历史证据。第 2–12 节描述目标设计，不代表所有退出条件已经满足。结果表主路径不处理原始 FCS；最近一次开发补齐了 V4 SASA 输入范围/方法比较的分层记录及 FCS/补偿矩阵 artifact 崩溃恢复，并完成 240 项全仓回归。没有代填人工 QC。

#### 当前完成度与真实数据位置

| 模块 | 当前进度 | 下一步 / 限制 |
| --- | --- | --- |
| Phase 0–1：项目、构建体、样本 | 冻结项目、版本、哈希、重启恢复、构建体差异预览、临时标签登记已实现；真实 activation 项目已有 3 批、162 样本 | 临时标签不是已核验序列；需补齐构建体身份后才能进入对应效应/结构定位路径 |
| Phase 2：FCS/FlowJo/门控 | 三批已各登记 54 个 FlowJo `% of Parent` 逐孔结果；原始 XLS、标准化百分比、来源 SHA-256 和可选的同批 workspace 关联均可追溯。已圈门结果表路径不读取或重算 FCS；文本小数分隔符显式选择并写入 parser provenance，避免逗号文本被静默改值；未提供的事件分子/分母保持未知 | 162 个样本仍需人工 QC；QC 决定绑定到具体结果表导入，不要求关联 FlowJo workspace；导入的表须与已登记批次样本孔位唯一匹配，缺孔表不能计算效应 |
| Phase 3–4：QC、效应与结构图 | 已有百分比专用数据模型、导入 API/UI 和按同批、同抗体面板计算的描述性效应入口；结果类型支持没有事件计数的 `WellRate`；结果表 QC 逐行保存、不重建整页，串行使用最新 revision、显示审核进度、只对空白行提供审核人填充快捷项，并识别 QC 改动导致的旧效应版本 | 当前 162 个样本尚无人工 QC，未生成真实效应；E1–E7 序列和抗体身份仍是待定临时标签，当前不可作结构残基归因；不报告事件级误差或生物学显著性 |
| Phase 5：Paraplume | 完整 large 模型 GPU/CPU 推理、缓存、后台任务、进度和恢复已有实现与既往验收；另以隔离临时项目和合成 VH/VL 完成 worker→任务进度→预测 artifact 持久化端到端运行，五种 embedding 阶段均产生进度，结果保持 `PREDICTED_UNREVIEWED`；用户确认仅本机非商业研究、不分发 | 硬件真实 CUDA OOM 尚未实测；用户确认目前没有与 MIT 项目组合的书面许可；上游“不与不兼容许可证软件组合”条款未限定为分发场景，故本机集成调用适用范围仍未确认；商业集成需事先书面授权 |
| Phase 6：GRASP | 固定 commit 的官方 `generate_restr.py` parser fixture 可解析并往返序列化导出的 IR/RPR；导出包格式/API、审核来源和重启下载已验证；模型未安装或运行 | 当前无可用配对距离证据，维持 IR-only；不安装或运行 GRASP 模型 |
| Phase 7：恢复与发布 | 结果表原文件与标准化 artifact 已写入 manifest；路径 schema 修复后重启载入通过；百分比表导入在原件发布和 manifest 发布两个硬崩溃点均可重启并重试；支持包重启/API 下载可复现审计与项目版本记录，默认省略原件、仅显式选择时附入；真实浏览器只读检查显示 162 个待审 QC 表单且无写请求；合成 API/浏览器闭环验证 QC 门槛、效应 provenance、逐孔审核连续保存、逐次 revision 递增及 QC 改动后的效应重算入口，所有浏览器写请求由 stub 截获且真实 ledger 未变化 | 真实 activation 尚未完成人工 QC→效应分析；其他 artifact 写入路径和真实人工 QC 决定后的浏览器验收仍需完成 |

本次数据库核对：Myc 项目 `4823bdbae91440bb` 为 COMPLETE / revision 8，1 批、4 样本；CTLA4 activation 项目 `e0f2c26650f04681` 为 COMPLETE / revision 16，3 批、162 样本。后者有 3 个 `flowjo_percentage_imports`，各含 54/54 个唯一孔并与各自 FlowJo workspace import 唯一关联。`gate_analyses`、`qc_reviews`、`effect_analyses` 及百分比专用效应分析仍无记录。这里 COMPLETE 是项目登记状态，**不表示实验分析全部完成**。原始样本分别位于 `D:/MFL lab/FACS/260921/` 和 `D:/MFL lab/FACS/260903/`。

用户已确认且后续须沿用：H1/H2/H3 是 Myc 染色的三个突变体，暂用临时标签，H4=NC；260903 中 CTLA4-WT 为野生型，E1–E7 为七个突变体，NC 为仅抗体-SNIPR 的自激活阴性对照，L32/L49/Y159/W202/W212 为抗体，CD80 为阳性配体。三目录是独立批次，按目录标签和工作表对应；每种试剂仅在一批出现，不能把三个目录合并成同一抗体的独立重复。用户确认采集未补偿；结果表中显示的 FlowJo `% of Parent` 作为百分比观测输入，不能推导事件分母。用户明确无需开发者处理原始 FCS：已经圈门的结果表可直接导入；WSP 仅作可选方法来源关联，不是结果表 QC 或效应分析的前提。E1–E7 序列和抗体身份当前待定，本轮不核验这两项。

三份结果表现已登记为只读输入：`table-1`→`260903PLATE1-l32-49`、`table-2`→`260903PLATE2-Y159-W202`、`table-3`→`260903PLATE3-W212-CD80`。三份文件表头都标注 `Live/GFP-A, SSC-A subset/mCherry-A, SSC-A subset | Freq. of Parent`；parser 只取每孔前两列，第三张表的汇总区不重复计入。当前 table-3 的 W212 `B5`（CTLA4-E2 技术重复 1）为 `21.0 %`，与该表的三孔汇总首值一致；前次“table-1 与 table-3 冲突”的提醒不适用于当前三份文件，因为它们分别属于三个独立批次。人工 QC 尚未登记。每份表现在将独立承载其逐孔人工 QC；同批 WSP 关联只作额外来源记录。

#### 已完成的结构分析规则调整

这些调整属于结构分析基础能力，后续实验项目冻结和结构定位须沿用：

- **直接 PDB 输入软准入。** 覆盖率、序列一致性等准入指标不达标时显示警告并继续；无法解析结构或建立有效映射等执行错误仍明确报告。抗原链优先根据 SIFTS/UniProt 映射选择，复合物结合实体注释和抗体接触证据识别，并记录选择依据及不确定性。
- **标准序列与工程突变。** 以 UniProt 标准序列为参考，保留原结构坐标，标注工程突变并从设计候选中排除；不把工程突变当作天然 WT。工程突变的坐标仍参与原生 SASA 的几何遮挡，这与“排除候选位点”是不同处理步骤。
- **只用选定蛋白链计算 SASA。** 排除结晶水、离子、非蛋白配体和其他链，使用固定探针/采样点/元素半径；原结构保留，蛋白伙伴和糖链的接触风险另作注释。氯离子不再误标为天然蛋白界面。
- **SASA 输入范围与方法比较。** 现分别记录：原始选定链的原生 SASA（保留工程突变及链侧翼）、规范化映射 ECD 上按相同参数重算的原生参考、以及 ProCaliper 对同一规范化输入的结果。ProCaliper 必须回报并匹配输入结构 SHA、链、序列哈希/长度；采样半径、采样点数和 rSASA 参考表须与原生参考一致。Biopython 版本单独记录；版本不同会显示 `VERSION_DIFFERS`，但只要算法参数相符仍直接比较实际结果。输入/算法设置无法核实时不把结果标成“一致”，严格候选也会阻断。原链与规范化 ECD 的分类变化作为 `CONTEXT_SENSITIVE` 单独保存并阻断受影响位点的严格推荐。既有冻结 artifact 没有这些列时兼容读取，但 UI 明确显示 `LEGACY_SCOPE_UNVERIFIED`。当前 V4 默认配置使用 `G:\anaconda\envs\ecd_procaliper\python.exe`；本轮通过 V4 adapter 和该环境完成规范化单链 worker 冒烟，ProCaliper 0.2.1 成功返回 109/109 残基并通过结构 SHA、链 ID、序列 SHA/长度范围校验；这验证 worker 接线，不代表完整 V4 作业或真实实验身份已端到端验收。

相关入口：`retrieval.py`、`pipeline.py`、`structure_analysis.py`、`joint_design_v4/engine.py`、`integrations/workers/procaliper_worker.py`。冻结项目的结构、映射、分区与表面图不得原地替换；采用更新的分析证据须走明确的新项目/新版本流程并重新核验哈希。

#### 新作业引擎策略

用户决定新作业默认使用 `joint_v4_ml`，不再兼容旧引擎入口。当前 API、CLI、网页和 pipeline 已将新作业限制为 V4；旧参数会被拒绝。V4 仍复用 `surface_patch` 表面几何和 `joint_design.surface_units` 共用计算模块，因此后续精简源码时须保留这些运行依赖或先迁移它们。历史方法文档和已生成结果暂留作审计记录，不代表旧引擎仍可运行。结果表入口、V4 SASA 范围比较、FCS 与补偿矩阵 artifact 崩溃恢复均已进入最新全仓回归；最近一次全仓测试 240 项通过，明细见文末。

#### 验证证据与下一步顺序

最近一次完整验证为 `pytest tests` **252 passed**（626.32 秒，6 条依赖/缓存警告）。新增 FCS 与补偿矩阵 artifact 硬退出用例已进入全仓回归；ProCaliper V4 worker 实跑 109 个规范化残基并通过完整输入指纹校验。Ruff、相关 Python 编译、浏览器脚本语法和 `git diff --check` 均通过。2026-09-29 完整工作台浏览器 smoke 再次通过结构证据、效应图层、导出、Paraplume 进度、FlowJo 状态和临时身份样本匹配等检查，无 JavaScript 异常；样本匹配后等待新项目视图中的文件输入框可见，修复了原 smoke 在旧 DOM 脱离与新视图渲染之间的等待竞态。真实 CTLA4 结果表路径的只读 smoke 显示 162 个空白 QC 表单、未补偿状态及分母未知，无网络写入。当前工作树仍有大量已修改/未跟踪文件，继续改动前先核对 `git status --short`，避免覆盖已有工作。

本轮新增路径验证：`tests/test_flowjo_result_tables.py` 与 `tests/test_experimental_effects.py` **21 passed**；相关 Python `py_compile`、`node --check web/experiments.js`、Ruff 检查及 `git diff --check` 通过。API 已将三张真实 FlowJo 表各登记为 54/54 个孔（项目 revision 16）；在未人工 QC 时效应 API 返回 `QC_HOLD`，数据库仍无人工 QC 或效应记录。没有调用 FlowKit worker 处理这些表对应的原始 FCS。

2026-09-28 ProCaliper worker 复核：V4 `_adapter()` 按默认配置选用 `G:\anaconda\envs\ecd_procaliper\python.exe`，该环境可导入 ProCaliper 0.2.1 和 Biopython 1.88。使用项目已有 PDB 的单链规范化临时输入，经 `ProCaliperAdapter`/worker 实算 109 个残基，输出完整且 `_validate_procaliper_input_scope()` 对结构 SHA-256、链、序列 SHA-256 和长度全部通过；provenance 记录 Shrake-Rupley、探针 1.4 Å、100 点。直接把多链原始复合物交给单链 worker 会因其返回全文件残基数而失败，符合 V4 先规范化所选 ECD 链的输入契约。该检查未运行完整 `run_joint_v4_ml`，也未触碰用户实验 FCS/QC。

2026-09-29 浏览器回归：`tests/browser_experiments.mjs` 连接现有 V4 服务 `127.0.0.1:8766` 后全项通过。样本匹配成功后，旧 smoke 仅等待旧区块 detached，早于整页项目工具重新渲染；等待条件现同时要求项目的 FCS 区域和可见文件输入已重建。临时身份样本由浏览器路由 stub 验证请求内容，没有写入真实 CTLA4 样本或 QC。

2026-09-29 增加 FCS 原件发布后的硬退出恢复：若只读、content-addressed FCS 文件已发布但数据库尚未登记，重试仅在目标为普通文件且大小/SHA-256 与再次上传内容一致时复用；符号链接或不匹配目标拒绝覆盖。临时子进程在只读发布后 `os._exit(80)`，重启后重试成功并只产生一个 raw-FCS artifact 和 import session。使用合成字节 fixture，没有导入或重处理真实 CTLA4 FCS。

2026-09-29 补齐补偿矩阵 artifact 提交恢复：矩阵 JSON 按规范化字节 SHA-256 命名；硬退出后重试验证并复用未登记的只读原件；相同内容再次上传返回既有 artifact，不增加 revision。分别在原件发布后及 manifest 写入后/SQLite 提交前执行 `os._exit`，重启后 manifest 修复、重试登记和项目完整性均通过。`tests/test_experiment_projects.py` 矩阵/FCS 恢复与矩阵 gate 集成聚焦用例 **4 passed**；最新全仓 `pytest tests` **240 passed**（489.52 秒）；Ruff、py_compile、`git diff --check` 通过。合成 spill matrix 只验证软件路径，不作为 260903 实验补偿数据；用户已确认当前采集未补偿。

另增加百分比结果表导入的硬崩溃恢复验收：子进程分别在 content-addressed 原表只读发布后、manifest 更新后且 SQLite 尚未提交时退出。重启会按已提交 artifact registry 修复派生 manifest；只有原件路径非符号链接且字节长度与 SHA-256 完全一致时才复用孤立原件。两种崩溃场景均成功重试并留下可读的 revision 5 项目，新增聚焦用例 **2 passed**；Ruff、相关 Python `py_compile` 与 `git diff --check` 通过。该测试使用临时合成项目，没有更改真实 CTLA4 项目；本轮未跑全仓回归。

后续完整实验项目模块回归 **35 passed**，覆盖上述崩溃恢复、跨进程同 revision 冲突、FlowJo 样本 QC 关联及百分比效应 API 闭环。并发批次写入在读取 revision 前使用 `BEGIN IMMEDIATE`，确保第二个进程得到明确 `REVISION_CONFLICT`；FlowJo QC 关联在缺少显式 `well_id` 时仅回退到可唯一解析的 `sample_label` 孔位。合成闭环验证未完成 QC 时返回 `QC_HOLD`，QC 后产生预期 `-35 pp` 对比，并将事件分母未知、未补偿状态和 CD80 `POSITIVE_LIGAND_CONTROL` 排除原因写入可回读 artifact。只读浏览器 smoke 检查真实 CTLA4 项目，显示 162 个空白 QC 表单、W212 B5 的 `21.0 %`、分母未知和未补偿说明；没有执行任何 POST 或填入 QC。全仓测试仍未运行。

2026-09-28 根据用户明确的结果表工作方式，修正原有隐藏依赖：已圈门百分比表现在以 `FLOWJO_PERCENTAGE_IMPORT` 自身作为人工 QC 来源，无需同批 WSP；每份决定仍须精确绑定该表内的样本，非表内样本会拒绝。WSP 关联仅作为可选方法来源，不会自动替代对结果表的审核。数据库 schema v14 保留既有 gate/FlowJo QC 记录并开放新的 QC 来源类型。若结果表缺少任一已登记批次孔位，效应 API 返回 `QC_HOLD`，页面显示缺失孔并要求重新导入完整结果，不把缺失值视为阴性或零。无 WSP 的合成结果表→QC hold→逐样本审核→效应完整 API 回归通过；包含实验项目、结果表 parser 和效应内核的测试组 **56 passed**，Ruff、Python 编译、前端语法与 `git diff --check` 通过。真实项目工作台浏览器 smoke 将三个导入的 WSP 关联置空，仍显示 **162** 个空白 QC 表单、无页面异常、无网络写入；项目保持 revision 16、3 张结果表、0 条 QC、0 个效应。合成样本登记仍沿用现有样本/孔位映射；该百分比分析路径不调用 FlowKit。随后全仓回归已于 2026-09-28 完成并通过；实际人工 QC 和效应仍等待用户审核决定。

下一阶段主线为 **结果表→人工 QC→描述性效应→（身份核验后）结构回溯及恢复**。沿用 Phase 0–1 已有项目/样本基础，按以下顺序推进：

1. **人工 QC。** 使用每份结果表导入记录和孔位匹配审核 162 个样本；不要求读取或重算原始 FCS，也不要求关联 WSP。开发者不得代填 `PASS/PASS_WITH_OVERRIDE/HOLD/EXCLUDE`。显示值保持原始百分比精度，事件分母未知。
2. **生成面板内描述性效应。** 所有纳入样本均需先有人工 QC；按批次及固定抗体面板分别比较临时抗原标签、WT 和同期 NC。三孔是技术重复；不跨不同抗体合并、不计算事件级置信区间或显著性，CD80 保持阳性配体控制角色。
3. **序列/身份后续核验。** E1–E7 突变清单与抗体构建体身份保持待定，当前不据临时标签生成结构位点归因；等用户补充并授权进入该部分后，再校验构建体与冻结结构映射。

2026-09-28 完成结果表路径与 V4 候选状态的全仓验收：根据用户更新后的工作簿，W212 首个技术重复已与 table-3 对应孔一致；fixture 现按所有面板逐孔比较工作簿与原始表，不再保留已失效的冲突断言。V4 离线配置现在会同时关闭环境变量和设置文件中的外部引擎命令，防止预览测试误启动 WSL/GPU worker，并遵守 assay profile 的探针上限。Stage-2 单点/配对构建体保持独立条件状态、关联 parent Stage-1，明确仅在表型出现后才作后续验证，不能证明直接接触；饱和替换只列为 `NOT_SCORED` 备选，不生成模型分数。更新相应历史断言后，完整 `pytest tests` **232 passed**（353.91 秒）；全仓 Ruff、Python compileall、`node --check web/experiments.js` 和 `git diff --check` 通过。真实 CTLA4 项目只读浏览器 smoke 在去除 WSP 关联后显示 162 个未预填 QC 表单、`21.0 %` W212 B5、事件分母未知与未补偿状态，无网络写入或页面错误。当前项目仍为 revision 16、3 张表、0 条人工 QC、0 个效应；该百分比表主路径未处理原始 FCS。真实样本 QC/效应、序列/抗体身份核验仍未完成。

2026-09-28 完成 SASA 比较范围拆分：ProCaliper worker 现在返回规范化结构 SHA-256、链 ID、序列 SHA-256/长度、采样参数、rSASA 参考表及 Biopython 版本；V4 验证输入指纹，并让旧缓存因上下文契约版本变化而失效。应用在与 worker 相同的规范化映射 ECD 输入上重算原生 Shrake-Rupley SASA，再单独比较原始选定链与规范化链的范围效应。方法不一致、输入范围造成的表面分类变化分开记录；严格模式会阻断受影响候选，旧冻结表仍能读取并显式标记范围未核。结构证据 API 与 UI 分别显示方法比较和结构范围影响。聚焦测试 62 passed；完整 `pytest tests` **236 passed**（267.90 秒）；Ruff、compileall、前端/浏览器脚本语法、git diff 检查以及全工作台浏览器 smoke 通过。当时这条记录将当前 shell 环境误当作 V4 worker 环境；后续已按 V4 默认 adapter 所配置的 Python 环境实测 ProCaliper 0.2.1，确认 V4 可用（见上文 §1.4）。CTLA4 真实项目保持 revision 16、3 张百分比表、0 条 QC 和 0 个效应，未处理其原始 FCS。人工 QC、效应、身份核验仍等待用户输入。
4. **补齐发布剩余项。** 完成 Paraplume 结果与 GRASP 导出的 artifact 中断恢复；以用户人工 QC 决定完成真实结果表→描述性效应→重启浏览器验收；补齐 Paraplume 任务进度真实推理端到端验证；有可靠配对距离证据后再扩展 RPR。ProCaliper 的 V4 调用和同范围 SASA 比较已完成，发布前继续跑全仓回归。

2026-09-28 用户明确已圈门结果表可直接作为本项目的输入，不需要开发者处理原始 FCS。为使这个工作方式真正闭环，新增 `FLOWJO_PERCENTAGE_IMPORT` 人工 QC 来源：每个决定绑定具体结果表及孔位，不能再因缺少同批 WSP 而阻断 QC/描述性效应。百分比表与既有登记批次样本唯一匹配；WSP 保留为可选方法来源。数据库迁移 v14 保留历史 QC 并扩展来源类型。端到端合成验收将移除 FlowJo workspace，仅用已登记孔位和用户导出的表格完成 QC 门槛→百分比效应，覆盖未审核阻断和非表内样本拒绝；真实 CTLA4 项目尚未写入任何 QC。

2026-09-29 增补 gate/effect artifact 写入恢复：对分析文件只读发布后、manifest 写入后/SQLite 提交前分别注入硬退出。重启后已提交 project revision 与分析 registry 保持一致，manifest 从已提交 registry 修复；未提交分析不会出现在项目查询中，重试可成功。stale revision 回归同时确认 `RevisionConflict` 不被只读文件清理错误遮蔽。硬崩溃遗留的未登记分析文件仍保留在磁盘，不会自动清理。相同安全清理函数已用于结构证据、FlowJo 百分比导入/效应和 GRASP ZIP 的异常回滚。最后一轮清理函数复用及新增 crash fixtures 前的全仓 `pytest tests` **244 passed**（469.31 秒，6 条依赖/缓存警告）。随后 Paraplume artifact、GRASP ZIP 和 FlowJo 百分比效应各在只读产物发布后及 manifest 发布后执行硬退出；重启均恢复已提交 registry，重试成功，孤儿文件保留。最新 `tests/test_experiment_projects.py` **42 passed**（118.37 秒，2 条依赖/缓存警告），`tests/test_paraplume_tasks.py` **11 passed**（33.11 秒，2 条依赖/缓存警告）；Ruff、相关 Python `py_compile`、changed-file trailing whitespace 和 `git diff --check` 通过。测试使用合成数据；只读复核真实 CTLA4 项目仍为 COMPLETE/revision 16，162 个样本、3 张百分比表、0 条人工 QC、0 个效应分析；没有重算原始 FCS。真实 activation 项目 QC→效应仍等待人工决定；Paraplume 进度真实模型端到端、真实硬件 CUDA OOM、部署许可和有复合物证据后的 RPR 仍待验收。

2026-09-29 为结果表 effect artifact 增加 `comparison_group` 和 `accepted_for_summary` 明细，工作台效应表新增可展开的逐孔散点图及原始值/QC/纳入状态表；散点按变体、WT、NC 分层，纵轴固定 0–100%，不绘制小样本密度、不重算统计量，也明确技术重复不等于生物学 n。API 回归 **1 passed**；只读真实项目浏览器 smoke 用 route stub 注入三点合成 effect，核对点数、原始百分比、PASS_WITH_OVERRIDE/HOLD 显示，无页面错误或网络写入；实际项目仍有 162 个未预填 QC 表单。Ruff、Python 编译、两个 JS 脚本语法检查通过；此视图只呈现用户已圈门百分比表，不读取 FCS。

## 2. 系统结构和持久化

### 2.1 组件划分

1. **实验项目服务**：项目/构建体/批次/样本及分析版本的 CRUD、输入校验和审计记录。
2. **项目文件存储**：原始 FCS、FlowJo workspace、门控修订导出、结果表、截图和约束文件。由 manifest 记录哈希与相对路径。
3. **实验分析服务**：FlowKit 读 FCS、补偿/变换/门控、事件计数、QC、效应及批次统计。
4. **结构映射服务**：根据项目冻结的 v4 残基映射和 surface graph 产生三层效应数据及可视化响应。
5. **抗体任务服务**：完整 Paraplume `--large` 推理、缓存、状态、硬件和模型 provenance。
6. **约束导出服务**：生成实验 IR、预测辅助 IR、人工 RPR 清单、合并 TXT 与校验报告。
7. **前端实验工作区**：按“项目 → 构建体 → 实验安排 → FCS 与门控 → 质控 → 效应分析 → 结构图 → 抗体证据 → GRASP 导出”导航。

后端仍保持当前仅回环访问的部署模型。文件 API 只允许项目根目录中的已登记相对路径；上传文件名不作为文件系统目标路径。FCS 和 ZIP 读取设文件大小/解压总量上限，导出用安全文件名。

### 2.2 数据库存储方案

新增实验数据库（建议 `workspace_data/experiments.sqlite3`），不要把大块事件数据塞入 SQLite；事件数据和原始文件保存在项目目录。SQLite 保存实体关系、版本号、任务状态、操作者/时间戳和相对路径。每项分析在事务中登记输入哈希和状态；生成文件先写临时路径，通过校验后原子改名，再提交为完成版本。崩溃留下的临时文件不得被当作完成结果。

建议初始表及主要唯一约束：

| 表 | 主要字段与约束 |
| --- | --- |
| `experiment_projects` | `id`, `source_job_id`, `root_path`, `reference_sequence`, `structure_sha256`, `mapping_sha256`, `surface_graph_sha256`, `sector_sha256`, `created_at`; 冻结结构上下文 |
| `project_revisions` | 项目参数、结构/映射/分区哈希和创建原因；只追加，不覆盖 |
| `antigen_constructs` | 项目、名称、参考/突变序列、突变 JSON、来源推荐/自定义、风险标记、映射状态、版本 |
| `antibody_constructs` | 亲本/变体、VH/VL 序列、链别、相对亲本的替换列表、单点/组合标记、版本 |
| `assay_batches` | 日期、批次标签、实验设置、负责人备注、是否独立批次 |
| `samples` | batch、抗原/抗体构建体、role、重复类别/编号、条件、FCS artifact、FCS 解析状态 |
| `artifacts` | 项目相对路径、类型、SHA-256、字节数、只读标志、创建时间；FCS 哈希唯一 |
| `gate_revisions` | parent gate、类型、通道、参数、补偿/变换引用、版本、操作者和备注 |
| `qc_reviews` | 样本/分析版本、PASS/PASS_WITH_OVERRIDE/HOLD/EXCLUDE、理由、审核人和时间 |
| `analysis_revisions` | 输入版本集合、参数 JSON、代码版本、统计输出路径/哈希、状态、创建者、父版本 |
| `paraplume_runs` | VH/VL 序列哈希、模型/权重版本、命令、设备、状态、输出路径和错误 |
| `constraint_sets` | 约束类型/条目、审核、输入 FASTA 哈希、校验报告和导出版本 |
| `audit_events` | 实体 ID、操作类型、操作者、时间、前后版本及备注 |

所有 JSON 设置使用带 `schema_version` 的 JSON Schema。构建体、门控、QC、分析、人工审核及约束修订都保留父版本，编辑操作创建新版本。重复 FCS 依据内容哈希只存一份，可由多个样本引用；不得因相同文件名判定重复。

### 2.3 项目目录约定

```text
workspace_data/experiments/<project_id>/
  project_manifest.json
  reference/{antigen_sequence.fasta, structure.cif, residue_mapping.csv,
             probe_sectors_v4.json, surface_graph_v4.json}
  constructs/{antigens.json, antibodies.json}
  raw/fcs/<sha256>.fcs
  raw/flowjo/<sha256>.wsp
  gates/revision_<n>.json
  analyses/<analysis_id>/{manifest.json, per_well.csv, batch_effects.csv,
                         probe_surface_effects.json, qc_report.json}
  paraplume/<run_id>/
  exports/<export_id>/
```

manifest 中用相对路径，避免工作区移动后绝对路径失效。启动时校验 manifest 与文件存在性/哈希；坏文件标记 `INTEGRITY_ERROR`，不自动覆盖或修复原始文件。

### 2.4 依赖和工具获取计划

本工作区现有 venv 锁定 NumPy 1.26.4。FlowKit 1.3.2 的依赖要求 NumPy >2，并包含用于 FCS/FlowJo 的 `flowio`、`flowutils` 等包，因此不得安装进现有 venv；使用 `.tool_envs/flowkit` 独立 Python 3.12 环境，通过 JSON stdin/stdout worker 接口调用。Paraplume 按用户指定版本 2.0.0 安装在 `.tool_envs/paraplume`，完整模型所需 PLM 权重已缓存在 `.tool_cache/hf`。PyPI 当前将 Paraplume 标为专有许可，正式运行前保存其发行包、LICENSE/metadata、权重来源和 SHA-256，并由管理员确认使用条件；本地工作区已按用户请求获取发行包不代表可再分发。

GRASP-JAX 不安装运行环境、不下载模型权重、不调用推理；为锁定文件接口，仅保留官方 `generate_restr.py` parser fixture 和 Apache-2.0 LICENSE，源码固定 commit/hash 见 `tooling/manifest.json`，并用其 parser 通过 serializer round-trip。FlowJo Desktop 是商业软件；计划不下载或复制 FlowJo 软件本身，只导入用户提供的 `.wsp`。FlowKit 对 FlowJo 10 workspace 的支持有限，故实现阶段应下载官方 `.wsp` 示例，明确通过 fixture 支持的门型与未支持项目，未知节点继续走用户处理流程；现已取得 260903 activation 配套 `.wsp`/FCS 并完成只读回放；正式项目门控/QC/效应链仍待完成（见 1.4 节）。

当前机器为 Python 3.12.4、RTX 3060 Laptop 6 GB。环境下载位置 `.tool_envs/`、wheel/cache 下载位置 `.tool_cache/` 均加入 `.gitignore`。依赖获取产物另写 `tooling/manifest.json`（工具、版本、下载 URL、wheel/source SHA-256、Python/CUDA、环境路径、许可证元数据、权重缓存状态）；不把环境目录或大权重提交到源码库。升级时需显式更新该 manifest 和锁文件，不能无版本安装“latest”。

**本工作区获取状态（2026-09-25）**：FlowKit 1.3.2 已装入 `.tool_envs/flowkit`，numpy 2.5.3、scipy 1.18.1、pandas 2.3.3；wheel SHA-256 已与 PyPI 发布信息核对。Paraplume 2.0.0 完整 `--large` 所需 IgBert、IgT5、ESM2、ProtT5 权重及 AntiBERTy 包内权重均已取得；四个 Hugging Face 仓库固定 commit、文件大小及摘要记录于 `tooling/manifest.json` 和 `.tool_cache/paraplume_weight_manifest.json`。PyTorch 2.6.0+cu124 官方 wheel SHA-256 已核验，RTX 3060 Laptop GPU CUDA 运算通过；以用户提供的 VH/VL 完成完整 `paraplume-infer file-to-paratope --large --gpu 0` 推理，输出分别为 120 和 107 个逐残基概率，均为有限的 [0,1] 数值。推理记录在 `.tool_cache/paraplume_smoke/user_example_gpu.log`，结果 pickle 在 `.tool_cache/paraplume_smoke/out-user/paratope_user_example.pkl`。首轮曾被 Transformers 对 torch<2.6 的安全版本门槛阻止；升级到官方 torch 2.6.0 后推理成功。GPU 路径无需 CPU 回退。Paraplume 输出为抗体 paratope 概率，不能单独作为抗原 epitope 定位；需要指定抗原并提供复合物结构、抗原突变扫描或其他界面证据。详细硬件、版本及所有权重摘要见 `tooling/manifest.json`。

`pip check` 仍报告 ESM 3.0.8 所声明的 `torchtext`、`torchvision` 未安装以及其 Biotite 精确版本要求与现装 1.7.1 不符；这是**依赖声明未完全满足**，不能记录为依赖检查通过。已测试的 Paraplume 五模型推理路径成功；若后续直接调用 ESM 的其他功能，应另建兼容环境并单独测试，不要为消除提示而直接改动已验证的 Paraplume 环境。`requirements-paraplume.lock` 中的 ESM/PyTorch 是本机 `file:///E:/...` wheel 引用；移机交接须按 `tooling/manifest.json` 的来源、版本和摘要重新获取，不能把该锁文件当成跨机器可直接安装的包。

## 3. 主要数据对象和业务规则

### 3.1 构建体

**抗原**可以从已完成 v4 结果的推荐构建导入，也可通过突变清单或完整序列录入。突变清单须校验位置范围、WT 氨基酸、重复位点、替换字母；完整序列须与项目参考序列全局对齐，展示替换/插入/缺失差异后由用户保存。系统记录设计风险而不阻止登记已经实际制作的高风险构建体。位点按权威映射投影；未映射/结构缺失位点保存在构建体与样本中，结果标记不参与定位。

**抗体**先保存亲本 VH 和 VL，再保存变体序列以及逐链相对亲本的变化。VH/VL 变化、单点/组合突变分别登记。序列修改后新建构建体版本，不覆盖旧序列。

**构建体级历史 Myc**是单独证据对象：记录检测日期、FCS/历史比例来源、构建体、WT/亲本比较对象、阳性率、MFI、门控/汇总来源和备注。它不得自动填到某次激活批次，也不得显示为同期表达；没有 FCS 的历史手填数标 `SUMMARY_ONLY`。

### 3.2 样本和批量匹配

样本 role 枚举至少包括 `ACTIVATION`、`NC`、`ANTIGEN_MYC`、`FAB_MYC`、`OTHER_CONTROL`。样本必需字段为批次、抗原构建体、抗体构建体（不适用时显式为 null）、role、技术/生物学重复编号、文件引用及条件。批量上传流程分为上传/哈希 → 显示文件名与元数据预览 → 文件-样本匹配表 → 用户确认 → 入库。未匹配文件进入 `UNASSIGNED` 暂存区，不进入分析。历史百分比可无 FCS 录入，但来源强制 `SUMMARY_ONLY`，无法参加基于事件计数的质控。

## 4. API 和任务协议

沿用现有 FastAPI、SQLite 作业状态/SSE 进度、项目文件下载和后台执行模式；增加实验路由，执行中的大任务返回 `202 + task_id`。建议接口：

| 路由 | 用途 |
| --- | --- |
| `POST/GET /api/v4/experiment-projects` | 从已完成 v4 job 建项目、列出项目 |
| `GET /api/v4/experiment-projects/{id}` | 读取 manifest 和版本摘要 |
| `POST /api/v4/experiment-projects/{id}/antigen-constructs` | 导入推荐或建立自定义抗原构建体 |
| `POST /api/v4/experiment-projects/{id}/antibody-constructs` | 建亲本/变体抗体构建体 |
| `POST /api/v4/experiment-projects/{id}/batches` | 建实验批次与条件设置 |
| `POST /api/v4/experiment-projects/{id}/samples/import` | 批量登记 FCS/FlowJo 文件，返回待确认匹配表 |
| `POST /api/v4/experiment-projects/{id}/samples/confirm` | 确认样本映射并入库 |
| `POST /api/v4/experiment-projects/{id}/analyses` | 建立新的门控/QC/效应分析修订并排入队列 |
| `GET /api/v4/experiment-tasks/{task_id}` / `events` | 查看后台状态、进度、日志和错误 |
| `GET/POST /api/v4/experiment-projects/{id}/gates` | 获取/保存门控修订 |
| `POST /api/v4/experiment-projects/{id}/qc-reviews` | 审核 QC 状态与豁免理由 |
| `GET /api/v4/experiment-projects/{id}/analyses/{analysis_id}` | 获取 QC、孔/批效应、图层与 provenance |
| `POST /api/v4/experiment-projects/{id}/paraplume` | 建 VH/VL 预测任务 |
| `POST /api/v4/experiment-projects/{id}/constraint-exports` | 生成指定审核记录的 GRASP 输入包 |
| `GET /api/v4/experiment-projects/{id}/files/{artifact_id}` | 安全下载登记文件 |

任务生命周期固定为 `QUEUED → VALIDATING → RUNNING → COMPLETE/PARTIAL/FAILED`，可有 `WAITING_FOR_USER`（未知 FlowJo 门或文件匹配）和 `CANCELLED`。任务每步写入进度、输入修订号、输出 artifact 哈希和标准错误码。重复提交相同输入哈希时可复用完整缓存结果，但要产生新的分析记录，标出 `cache_source_analysis_id`。重启时恢复数据库任务状态；运行中断的任务标 `INTERRUPTED`，由用户显式续跑/克隆，不将半成品显示为完成。

### 4.1 HTTP 请求/响应契约

JSON 响应外层统一为 {schema_version, data, warnings, errors}；异步创建接口额外返回 task_id、status、events_url。成功的 artifact 返回不可变 artifact_id、sha256、size_bytes、media_type，不暴露磁盘路径。稳定错误码至少包括 PROJECT_SNAPSHOT_INCOMPLETE、REVISION_CONFLICT、FCS_HASH_MISMATCH、SAMPLE_MAPPING_REQUIRED、UNKNOWN_COMPENSATION_SOURCE、UNRESOLVED_FLOWJO_NODE、GATE_INVALID、QC_HOLD、REFERENCE_WINDOW_INSUFFICIENT、RESIDUE_UNMAPPED、PARAPLUME_MODEL_UNAVAILABLE、FASTA_CONSTRAINT_MISMATCH、RPR_EVIDENCE_REQUIRED。

写接口使用 Pydantic 输入模型、JSON Schema schema_version 和 expected_revision 乐观锁；创建资源返回 201，排队任务返回 202，修订冲突返回 409。时间为 UTC ISO-8601；比例小数 0–1，百分点字段后缀 _pp，百分比字段后缀 _pct。patch 不允许改变原始 FCS artifact 或已完成分析版本。

主要请求和响应字段：

    POST /api/v4/experiment-projects
    request: {source_job_id, name}
    response: {project_id, revision, status, frozen_hashes:{structure,mapping,sectors,surface_graph}}

    POST /api/v4/experiment-projects/{id}/antigen-constructs
    request: {name, source, mutations:[{position,wt,mut}], full_sequence?, expected_revision}
    response: {construct_id, sequence_sha256, mutations, mapping_status, risk_flags}

    POST /api/v4/experiment-projects/{id}/batches
    request: {name, date, assay, cell_line, dose, dose_unit, positive_gate_id, expected_revision}
    response: {batch_id, revision, status}

    POST /api/v4/experiment-projects/{id}/samples/import
    multipart: files[] and optional batch_id
    response: {upload_session_id, files:[{upload_id,sha256,channels,event_count,metadata,match_status}]}

    POST /api/v4/experiment-projects/{id}/samples/confirm
    request: {upload_session_id, matches:[{upload_id,batch_id,role,antigen_construct_id,antibody_construct_id,replicate_type,replicate_index}], expected_revision}
    response: {sample_ids, unmatched_upload_ids}

    POST /api/v4/experiment-projects/{id}/analyses
    request: {batch_ids, gate_revision_id, qc_policy_revision_id, analysis_parameters:{minimum_reference_window,bootstrap_iterations,seed}, expected_revision}
    response: {task_id,status,events_url}

    POST /api/v4/experiment-projects/{id}/antibody-constructs
    request: {name,parent_construct_id?,vh_sequence,vl_sequence,role:"PARENT|VARIANT",expected_revision}
    response: {construct_id,parent_id,sequence_hashes,heavy_mutations,light_mutations,mutation_class:"SINGLE|COMBINATION"}

    POST /api/v4/experiment-projects/{id}/gates
    request: {name,parent_gate_id?,kind:"RANGE_1D|RECTANGLE_2D|POLYGON_2D",channels,parameters,transform_id,compensation_id,expected_revision}
    response: {gate_revision_id,gate_id,version,validated,gate_tree}

    POST /api/v4/experiment-projects/{id}/qc-reviews
    request: {analysis_id,sample_id,status:"PASS|PASS_WITH_OVERRIDE|HOLD|EXCLUDE",reason,reviewer,expected_revision}
    response: {qc_review_id,status,audit_event_id}

    POST /api/v4/experiment-projects/{id}/paraplume
    request: {antibody_construct_ids:["abc1"],model_version:"2.0.0",model_size:"large",expected_revision}
    response: {task_id,status,events_url}

    POST /api/v4/experiment-projects/{id}/constraint-exports
    request: {analysis_id,include_experimental_ir:true,include_predicted_ir:false,manual_rpr_ids:[],merge:["experimental_ir"],expected_revision}
    response: {export_id,status:"VALID|REJECTED",artifact_ids,validation_errors}

分析查询响应按 per_well → per_batch → per_probe → per_sector/per_site 层级组织。每个效应记录提供 n_positive/n_parent、原始比例、NC/WT 比例、difference_pp、retention、decline_pct、有效批次数、区间方法/种子、QC flags 和输入版本 hash；无效窗口的 decline_pct 为 null 并带原因。响应必须保留逐孔和逐批数据，不能只返回最后着色矩阵。

### 4.2 内部 Python 与隔离 worker 接口

实验领域层保持 API 无关的纯函数：

    create_experiment_project(source_job_id, name) -> ExperimentProject
    validate_antigen_construct(reference, mutations, full_sequence, mapping) -> ConstructValidation
    register_antibody_construct(parent_id, vh, vl, parent_vh, parent_vl) -> AntibodyConstruct
    inspect_fcs(path, compensation_decision) -> FcsInspection
    import_flowjo_wsp(path, fcs_by_sample_id) -> WorkspaceImportPreview
    apply_gate_revision(sample, gate_revision) -> GateCounts
    calculate_sample_qc(sample, gate_counts, controls, qc_policy) -> list[QcFinding]
    calculate_batch_effects(sample_counts, batch_design, effect_policy) -> AnalysisResult
    map_effects_to_surface(result, constructs, frozen_mapping, frozen_graph, frozen_sectors) -> SurfaceEvidenceView
    build_grasp_export(selection, fasta_chains, mapping) -> ConstraintExport

隔离进程沿用项目已有的 JSON-over-stdin/stdout 协议：stdout 只含一个 JSON object，日志写 stderr，父进程负责超时、退出码、schema、位置/WT 和 hash 校验。输入文件必须由项目内 artifact_id 解析，不接受 UI 传来的任意绝对路径。

FlowKit worker 请求字段为 schema_version、request_id、operation（inspect_fcs/import_wsp/apply_gates）、artifact_id、gate_revision、compensation_decision。响应为 status、records、warnings、errors、provenance、runtime_seconds。inspect_fcs 返回 FCS 版本、通道元数据、补偿来源和总事件数；apply_gates 返回每门的完整事件数、父门数和比例。绘图抽样只可出现在 plot_sample 字段，不能参与 counts。

Paraplume worker 接收 heavy_sequence、light_sequence、model_version=2.0.0、model_size=large、device_preference=cuda_then_cpu 和 run_id。worker 在临时目录写出包含 sequence_heavy/sequence_light 列的 CSV，调用官方 paraplume-infer file-to-paratope <csv> --large --gpu 0 --result-folder <task-dir>；CPU 调用时使用 --gpu -1，绝不传 --small。worker 只读取自己刚生成的 CLI pickle 输出，将 model_prediction_heavy/light 转换成 position/aa/probability JSON，并附模型/权重 hash、实际设备、版本和日志 provenance。不得把用户上传的 pickle 当可信输入；任务目录每次新建，避免 CLI 覆盖项目文件。

GRASP serializer 暴露 build_fasta(chains=[ECD,VH,VL])、serialize_ir(chain_index,residue_index,aa)、serialize_rpr(residue1,residue2,max_distance_a,fdr)、validate_against_fasta(fasta,lines)。GRASP residue 语法是 chain_index-residue_index-residue_type，均为 1 起始；RPR 行为 residue1, residue2, distance_cutoff[, fdr]，IR 行只有一个 residue。若固定三链 FASTA 含重复序列并触发上游链重排，必须检测并拒绝歧义导出，或实现经官方源码 fixture 验证的稳定 chain map。

## 5. FCS、门控、质控及实验统计

### 5.1 FCS 导入与完整事件处理

FlowKit 作为 FCS 读写、补偿、变换和门控执行依赖。固定兼容版本并记录实际版本；导入解析 FCS 版本、事件数、参数/通道名、时间、补偿矩阵、仪器元数据和文件完整性。样本声明补偿来源：FCS 内矩阵、用户指定外部矩阵、明确无补偿。未知来源阻塞统计并要求选择，防止重复补偿。变换和补偿步骤按版本保存在 gate/analysis revision。

外部补偿矩阵目前接受 UTF-8 带标签方阵 CSV（首行是 detector 列名、首列是相同行标签）或 JSON `{schema_version, detectors, matrix, fluorochromes?}`；通道名必须唯一、矩阵有限且数值稳定。矩阵作为项目内只读 artifact 登记，门控修订固定引用 artifact ID，结果记录矩阵 SHA-256。执行前 FlowKit 还会检查 detector 是对应 FCS 的荧光通道。

原始 FCS 复制到 content-addressed 只读目录并校验复制前后哈希。绘图允许固定随机种子的事件抽样或密度下采样；门后计数、百分比、事件率和 QC 必须基于完整事件数组。相同哈希已存在时复用原始文件并另建样本引用。

### 5.2 FlowJo 10 workspace 导入

第一版 `.wsp` 导入解析 workspace 的样本/组、补偿和 gate tree，映射通道名与 FlowKit 样本标识。导入结果先展示可识别与不支持的门。范围门、二维矩形/多边形及父子门映射成内部 gate schema；未知变换、补偿或门类型标 `UNRESOLVED_FLOWJO_NODE` 并使分析进入 `WAITING_FOR_USER`。不得静默略过门，也不在导入时覆盖已有人工门控版本。

**已圈门结果表是另一种受支持的输入路径。** 用户提供 FlowJo `.xls/.xlsx` 或 CSV/TSV 后，系统只解析逐孔孔位和百分比、保留原文件与门路径、按选定批次映射现有样本；不读取或重算 FCS。文本数值的小数分隔符需显式选择点或逗号；混合/歧义文本会拒绝，原显示值及所选解析尺度写入 parser provenance；Excel 数值单元格按自身 number format 解释。表中若为 `% of Parent`，事件分子/分母保持未知，不能反推事件计数或事件级误差。每份决定直接绑定该结果表和孔位，无需 WSP。每份表必须记录其门控所依据的补偿状态（未应用、已应用或未知），不得从批次、FCS 原始状态或其他表格静默继承；报告按该表状态说明解释限制。260903 现有结果按用户确认的未补偿状态保留；若通道串色影响门内分类，阳性率也可能受影响。

### 5.3 可编辑门控模型

内部 gate schema 支持：

- `RANGE_1D`：通道、下/上界和是否含端点；
- `RECTANGLE_2D`：两通道和各自界限；
- `POLYGON_2D`：通道、闭合顶点序列及自交检查；
- `PARENT_CHILD`：稳定 ID、parent ID、执行顺序与版本；
- 补偿与变换的显式引用，禁止一份样本重复补偿。

标准层级以 `ACQUISITION_VALID → CELLS → SINGLES → LIVE → JURKAT_OR_K562 → TARGET_POSITIVE` 表达。通用实验可调整层名/父子关系，但激活读数需明确采用选定 Jurkat 父门事件数作分母、其 report-positive 数作分子。细胞系选择显式记录为 Jurkat、K562 或其它自定义身份，不从通道名推断。新 gate revision 重算生成新 analysis revision，并同时展示变更前后各层计数。

### 5.4 QC 算法

QC 以“提出异常、不静默删除”为默认策略：

1. **文件/通道**：检查哈希、通道名映射、事件数、时间参数、补偿状态；缺失/矛盾记录为错误或待处理。
2. **采集时间**：按时间等事件数窗口计算事件速率及关键通道中位数，使用全样本中位数和 MAD 标记速率骤变/荧光漂移。建议窗口大小、MAD 倍数写入参数快照并可配置；窗口标记只供审核。手动删窗需保存窗口范围、事件前后计数、理由和审核人。
3. **饱和/边界**：统计每通道位于 ADC 上下界的比例和异常峰值堆积，产生 detector saturation/boundary flags。
4. **阳性阈值**：在批次设计中指定用于判阳的合适阴性对照及阈值方法。支持从该 NC 通道分布按用户设置的分位数定阈值，或在已变换的 2D NC 图上由用户审核多边形边界；保留 NC 分布、阈值坐标/算法参数和审核人。阈值一经批次确认，对该批所有可比样本共用；改阈值必须创建新 gate/analysis revision。跨批次不强制共用荧光坐标，但同协议批次须记录可比性。
5. **事件量**：按层记录完整事件数，配置最小门后细胞数。门后事件不足导致该样本 `HOLD` 建议，不自动丢弃。
6. **同批对照**：每批建立 NC、WT 对照角色，按同一阳性阈值分析。NC 与 WT 的阳性率窗口不足、WT/NC 靠近 0 或 100、对照重复差异过大时标记窗口/动态范围风险。
7. **重复和门控一致性**：报告技术重复离散度、生物重复间差异、gate revision 不同、同批条件缺项；不同设置不自动合并。
8. **决定留痕**：每个样本/分析用 `PASS`、`PASS_WITH_OVERRIDE`、`HOLD`、`EXCLUDE`。人工豁免必须有理由、操作者和时间；`HOLD/EXCLUDE` 不参与主结果，豁免及来源持续出现在下游图和导出文件。

Myc 仅作构建体级表面展示提示。历史 Myc 可在项目中挂到构建体，不能模拟为该激活批次的同期表达；Myc 也不能独自证明天然折叠。没有同期表达/构象控制可按人工流程放行，但必须 `PASS_WITH_OVERRIDE` 并记录所缺控制。

### 5.5 效应统计和结果

事件计数路径将门后阳性/总数及 QC 状态传入 `WellCount`；FlowJo 结果表路径使用 `WellRate` 只记录观测百分比，不合成事件数。两个路径都按等权孔中位数、百分点差和 NC 校正窗口计算；同一对比中的 WT、NC、变体必须登记为同一种重复类型。只有有足够独立生物批次且数据类型允许时才做批次 bootstrap。当前百分比结果表每个面板只有一个独立批次、三孔是技术重复，因此只报告面板内描述性效应，不给 CI、事件级标准误或显著性；未补偿状态和被排除的 CD80 阳性配体控制需保留在来源说明中。

实验比较按 role 固定：

- **抗原突变**：固定抗体，WT 与突变抗原使用同批对照；输出 `Δ百分点`、NC 校正保留率/下降幅度。
- **抗体突变**：固定 WT 抗原，突变 Fab 与亲本 Fab 比较；各自适用 NC，抗体位点展示在 VH/VL 轨道。
- **双侧组合**：要求亲本/突变抗体 × WT/突变抗原四组完整，输出百分点差中之差和原始四格；只标为耦合线索。

WT−NC 窗口非正或未达到配置阈值时，保留原始阳性率和百分点差，保留率/下降写 null。负保留率或超过 100% 不截断，附解释旗标。批次是独立重复单位；至少 3 个批次才提供描述性 batch-bootstrap CI，并同屏显示每批点及一致方向数。只有技术重复展示均值/中位数和离散度，不报告生物学显著性。热图支持抗体×抗原探针；跨抗体相关性只列复测线索。多突变下降探针自动建议同一区域单点替换/独立替换设计，不自动拆分已观测效应。

## 6. 结构定位、效应图层及审计

### 6.1 冻结结构上下文

项目创建时从源 v4 job 拷贝抗原参考序列、结构、`residue_mapping.csv`、`probe_sectors_v4.json`、`surface_graph_v4.json`，保存哈希和源 job/schema/code revision。创建失败时不得生成缺项的“冻结项目”。分析输入结构哈希与项目结构不匹配时拒绝按旧编号映射，并提示建立项目修订。

### 6.2 效应映射

每个抗原构建体通过参考序列对齐取得实际差异位点，校验 WT aa、UniProt position、规范化位置和 auth 编号。没有图节点的位点仍参加统计，标 `UNMAPPED/NOT_SURFACE_NODE`。

界面严格分三层：

1. **直接观测**：真实突变点高亮，点开显示完整构建体效应、批次、FCS/孔、QC、Myc 来源。组合每个点只标记 `COMBINATION_EFFECT`，禁止变成各自实测下降。
2. **分区效应**：只有所有突变同属一个 v4 sector、实验窗口合格且 QC 通过的探针才参加分区汇总。先按探针汇总跨批中位数，再对不同探针取中位数；同一探针多个批次不增加其权重。附探针数、批次覆盖、范围/可用时的批次重抽样区间。跨 sector 组合只连线/描边所及区域，状态 `CROSS_SECTOR_UNASSIGNED`。未测试与已测试零效应使用不同图例。
3. **位点敏感性**：只由真单点或明确成对去卷积实验产生，按氨基酸替换保留分项，显示效应/批次一致性/构象或展示证据。缺组合内单点数据时停止在区域级。

可选渐淡由冻结图的 Dijkstra 最短路径计算，高斯带宽默认 5 Å、截断 12 Å；仅 `VISUALIZATION_ONLY`，不进入区域统计、证据等级或 GRASP。抗原层和抗体 VH/VL 层独立。效应层叠加结构可信度、SASA 方法分歧、历史 Myc、人工豁免等图标，不将它们揉成单个未经验证的概率。

证据标签建议由可观察设计驱动：`UNTESTED`、`REGION_PERTURBATION`、`SITE_SENSITIVE_SINGLE`、`SITE_REPLICATED`，后者要求至少 3 个方向一致的独立批次及规定的展示/构象控制。技术阈值先定义为配置项；界面将工程阈值说明为复测排序规则，经本体系验证后再定稿。证据摘要导出结构截图、结构及映射版本、逐孔/逐批效应、区域/单点层、限制和 provenance。

## 7. Paraplume 抗体证据

新增独立 worker/adapter，锁定 Paraplume 2.0.0、完整 `--large` 模型及模型权重 SHA-256。启动先探测 GPU 和显存；尝试完整模型 GPU。若当前 6 GB 显存设备不足以加载完整模型，或运行 OOM，记录原因后清理并重试一次；仍失败时按同一完整模型 CPU 运行。禁止静默改小模型。对于硬件不可用或命令/权重缺失，任务明确失败，并允许 UI 显示尚未计算。

按 VH、VL 序列和模型权重散列缓存逐残基概率。每个 run 保存链、序列哈希、Paraplume 版本、权重版本/hash、命令、GPU/CPU 型号、实际设备、开始/结束时间、stdout/stderr 摘要、错误、输出文件和状态。接口要校验输出逐残基数量、输入序列氨基酸及 1 起始编号。

有突变实验时，实验展示/门控 QC/单点重复效应/人工复核是主证据，Paraplume 仅为独立排序证据并提示冲突。只有预测时按概率展示候选，标签必须是 `PREDICTED_UNREVIEWED`，审核确认后仍注明预测来源。实验阴性不自动否决预测位点；组合突变实验不把效应拆到多个抗体位点。

## 8. GRASP 输入接口与校验

从项目冻结序列导出 FASTA，固定 record 顺序与链序号：`1=ECD、2=VH、3=VL`。输出独立文件：

- `experimental_ir.tsv`：实验证据支持的 IR；
- `predicted_ir.tsv`：用户审核选择的预测辅助 IR；
- `manual_rpr.tsv`：用户明确录入的 RPR；
- 用户选择后的 `grasp_constraints.txt` 合并文件；
- `export_manifest.json`、校验报告、导出 FASTA 和相关证据记录。

IR 使用 `chain_index-residue_index-residue_type`（例 `1-10-G`）。RPR 使用 `1-10-G, 2-25-Y, 8.0, 0.05`，末尾 FDR 可省略且默认 0.05。导出前把项目位置映射到 ECD/VH/VL 各自链内位置，核对序列长度、chain 编号、氨基酸、重复约束、字段数量和 FASTA 哈希。歧义或不匹配直接拒绝导出并列出行号/理由。

RPR 只有人工审核的残基配对、距离上限和错误率齐全且验证通过时才输出；配对距离须由明确配对证据或人工指定来源支持，并按项目结构计算/校验。表面邻近、同区效应、抗原-抗体双突变交互及 Paraplume 概率都不能自动生成 RPR。证据不足时导出 `IR_ONLY` 并在 manifest 解释缺项。每行约束保留 evidence source、构建体、批次、QC/豁免状态、审核人、分析版本和映射版本。

GRASP-JAX 规范由官方 [GRASP-JAX `generate_restr.py`](https://github.com/xiergo/GRASP-JAX/blob/main/generate_restr.py) 与其 [README 格式说明](https://github.com/xiergo/GRASP-MS#restraints-file-format) 核对：FASTA 记录顺序决定从 1 开始的 chain index；IR 每行是一个 `chain-residue-aa`；RPR 每行是 `chain-residue-aa, chain-residue-aa, distance_cutoff[, fdr]`，默认 FDR 为 0.05。固定 commit 的 parser 源码与 Apache-2.0 许可证保存在 `tests/fixtures/grasp/`；`tests/test_grasp_official_format.py` 已将导出的 IR/RPR 通过该官方 parser 转为 restraint arrays 并往返序列化校验。此验收覆盖输入格式 parser，不代表 GRASP-JAX 推理环境或模型已安装/运行。缺 RPR 时导出包状态固定为 `IR_ONLY`，不把“文件校验通过”误报为有距离约束。

## 9. 网页流程

1. **项目列表/创建**：从完成的 v4 job 创建，显示结构/分区哈希，明确若源结果不完整不能建冻结项目。
2. **构建体页**：推荐导入、自定义突变/全序列比较、抗原风险、VH/VL parent/variant、单点/组合标记、历史 Myc。
3. **实验安排**：建批次、选择构建体、角色、浓度、细胞系、重复，预览预期比较组是否完整并提示 NC/WT 缺项。
4. **上传与匹配**：拖入 FCS/WSP；展示文件/hash/metadata、重复哈希、批量样本映射表；用户确认后写样本记录。
5. **门控**：FlowJo 导入状态和未解释节点；编辑 1D 范围门、2D 矩形/多边形、父子关系；散点/密度图可抽样，图层显示完整事件计数。
6. **质控**：按文件、时间、通道、门后事件数、阈值和重复展示 flags；支持覆盖/排除/手工窗口但要求原因，显示前后事件数。
7. **效应分析**：显示逐孔点、每批 NC/WT 窗口、异常和批次 bootstrap、抗体×探针热图、抗体突变 VH/VL 轨道。
8. **结构图**：切换直接观测/分区/单点/显示渐淡图层；点击结构节点回到证据与 QC 来源；Myc 及结构可信度作为独立徽标。
9. **抗体证据**：亲本/变体序列、Paraplume worker 状态和 VH/VL 逐残基概率，标记实验/预测/审核来源。
10. **GRASP 导出**：选择 IR 清单、预测辅助清单和 RPR；预览 FASTA 编号、每行证据、错误列表；成功后下载单项和完整审计包。

顶层显示当前项目、活动批次、活动 gate/analysis revision 和后台任务。任一步关闭浏览器后状态由服务端恢复；按钮在后台处理时可继续访问其它项目，但不能编辑正在被该任务读取的不可变输入版本。

## 10. 文件/分析复现和导出

单个完整分析从原始 FCS 和锁定 gate revision 重跑，写出：

- 每孔计数、每层 gate 计数及阳性率 CSV；
- 每批对照窗口和效应 CSV；
- 构建体、探针、sector 和位点证据 JSON；
- QC 标记/人工审核表；
- gate JSON、FlowJo 导入报告和图像；
- 结构映射、surface graph/sector/structure 哈希引用；
- 任务日志、参数、依赖版本、算法/代码 revision 和随机种子 manifest。

项目级 ZIP 包含原始 FCS（可选用户勾选并再次提示体积）、FlowJo 文件、全部版本表、manifest、图像、Paraplume 输出和 GRASP 文件。原始文件默认不加入精简分析导出，但必须可从 artifact ID 单独下载。CSV/JSON 和界面显示采用同一分析输出对象，避免各自重算或四舍五入差异。

## 11. 分期实施及依赖

### Phase 0 — schema、迁移、项目快照

- 新建实验数据库、迁移版本表和 artifact store；写 JSON Schema。
- 从 v4 job 创建冻结项目，拷贝序列、结构、映射、sector、surface graph，记录哈希。
- 建立审计事件、通用 task 状态与文件安全路由。
- 完成后端 API 基本 CRUD、项目页面骨架和服务重启恢复。

**退出条件：** 从真实完整 v4 job 建项目；所有快照 hash 可校验；服务重启后列表和文件关系完整；旧 pipeline job 行为不变。

### Phase 1 — 构建体、批次、样本和上传匹配

- 完成 antigen 推荐/自定义清单/完整序列解析及映射风险；antibody parent/variant/VH/VL/单点组合。
- 加入批次、角色、条件/剂量、重复元数据、历史 Myc 汇总。
- FCS/WSP 哈希去重、只读复制、文件解析元数据、批量匹配预览与确认。
- 做构建体及实验安排 UI。

**退出条件：** 推荐与自定义 antigen 共存；完整序列能核对差异；错 WT/重复位置能拒绝；无法映射的突变仍可保存；文件匹配未确认前不可分析；哈希重复不多复制。

### Phase 2 — FCS/FLOWJO 和门控编辑

- 锁定 FlowKit 版本，实现 FCS 元数据/通道/补偿/transform 读取和事件表缓存。
- 实现 gate schema、1D/2D gate editor、FlowJo 10 WSP 导入及待处理节点。
- 接受带标签 CSV/JSON 外部补偿矩阵，并将其不可变 artifact ID/hash 绑定到 gate revision 和分析结果。
- 执行门控树、分层事件计数和阳性统计；图像只抽样。
- 保存 gate revisions 和样本级门控差异。

**退出条件：** 测试 fixture 覆盖真实 FCS 与 WSP；未解释门会阻止分析并能用户处理；重复补偿被拒绝；gate 修改产生新版本和重算；统计计数等于完整事件结果。

### Phase 3 — QC、人工审核、效应统计

- 文件、采集稳定性、饱和边界、事件数、同批阈值、重复和 gate 一致性 QC。
- UI 标记异常，人工选择删除时间窗须有理由/前后计数。
- 连接 `experimental_effects.py`，加入 API 序列化和 analysis manifest；批次整体 bootstrap 对多探针区域按同一批次重抽样。
- 支持 PASS/PASS_WITH_OVERRIDE/HOLD/EXCLUDE、历史 Myc 和 summary-only 输入。
- 构建逐孔点图、重复离散度、batch control dashboard、条件效应表和热图。

**退出条件：** 覆盖原始要求中的抗原变异、抗体变异和四格交互；无效 WT−NC 时 ratio=null；只用批次作 n；豁免原因在结果序列化中可追溯；相同分析输入重跑输出一致。

### Phase 4 — 结构效应图和证据分级

- 从项目快照加载 mapping、sector 和 graph，严格核 hash 与 WT aa。
- 实现 direct/sector/single 三层图层、颜色规则、点击联动、VH/VL mutation tracks、渐淡最短路。
- 实现组合跨 sector 不归属、只组合证据停区域、探针/批次计数和数据导出。
- 叠加结构可靠性/SASA/Myc/override 分立状态及 `UNMAPPED` 数据。

**退出条件：** 所有初始要求的组合/单点边界通过 fixture；改变结构映射/图 hash 会阻止分析复用；sector 色值可以逐条回溯至孔/批次；渐淡从不变成位点证据/GRASP IR。

### Phase 5 — Paraplume 2.0.0 执行和缓存

- 实现预检、完整 `--large` GPU → 完整模型 CPU 策略、任务进度和错误。
- VH/VL 按序列和权重 hash 缓存；校验位置/氨基酸；前端区分预测、实验和人工审阅。
- 将预测辅助 IR 候选列入审核页面但不自动升级实验等级。

**退出条件：** 可在真实配置运行完整模型；GPU 不足会转 CPU 并留 provenance；不可能静默降小模型；失败可在页面诊断；仅预测不会显示实验验证。

### Phase 6 — GRASP constraint 输入和导出包

- 固定官方/目标 GRASP-JAX 语法 fixture，建立 FASTA chain index 与各链位置映射。
- 输出实验 IR、预测辅助 IR、人工 RPR 独立清单和选定合并文本。
- 逐行校验序列、1 起始位置、氨基酸、配对残基、距离/error 字段、重复和审核 provenance。
- 生成导出预览、错误明细、FASTA/constraints/manifest ZIP。

**退出条件：** 故意不匹配 FASTA 必定拒绝；无配对/距离证据时只有 IR-only；表面视觉扩散/功能交互/预测概率不自动成 RPR；每条约束可追到 QC 和分析版本。

### Phase 7 — 全链路恢复、发布验收

- 按最初用户操作脚本从项目创建跑到 GRASP 导出，再重启服务。
- 用同一原始 FCS、gate revision、参数和结构快照重建同一数据版本；比对事件计数、效应表、结构渲染输入和约束文件内容/hash。
- 测任务失败/中断、重复文件、部分映射、FlowJo 不支持门、PASS_WITH_OVERRIDE、历史 Myc、无效窗口和 CPU fallback。
- 增加用户级可见错误码、支持包和完整审计下载。

**退出条件：** 全部验收表通过；原始 FCS SHA 保持不变；分析和导出版本可复现；服务重启后工作可继续；UI 无“已完成”误报。

## 12. 验收矩阵（需求逐项映射）

| 需求 | 验收证据 |
| --- | --- |
| 项目冻结 antigen/structure/mapping/sector/版本并可加批次 | 项目 manifest 哈希、迁移/重启测试、追加批次不改既有 revision |
| 推荐/自定义抗原；替换表或完整序列 | 差异预览、WT/位置核验、风险旗标、未映射保留 |
| parent VH/VL 与突变抗体 | 链级差异和单点/组合标志可复核 |
| 样本 role/重复/文件匹配/summary-only | 上传确认测试、未匹配阻断、历史手填来源标签 |
| FCS 只读/哈希/门控/分析/QC 版本 | 文件 hash 前后相同；所有修订新增版本并有审计记录 |
| FlowKit、FlowJo、任意 1D/2D 和层级门 | 测试文件导入修改重算；不支持节点显式待用户处理 |
| 文件/采集/补偿/阈值/QC | 双补偿阻断；MAD 窗口提示；排除仅接受来源分析完整窗口且要求审核人/理由；新结果记录来源 hash、选择和前后计数，真实 FCS SHA 不变 |
| 完整事件统计和激活分母 | 抽样绘图与完整事件计数 fixture 数值一致；分母明确为 Jurkat 父门 |
| 同批阈值/重复/不足窗口 | 未通过窗口不产 ratio；重复不合并事件；异常可见且不静默剔除 |
| Myc 历史与同期语义 | 日期和构建体关联；历史值不出现同批标识；不单独证明折叠 |
| 效应和跨批统计 | 已列的抗原/抗体/四格公式；batch median/bootstrap；结果表路径保留并分层显示逐孔原值、QC 与摘要纳入状态；技术重复不做生物显著性 |
| 结构直观层、sector 层、单点层 | 效应归属单测覆盖组合不拆分、跨区不分摊、单点层限制；浏览器运行时返回冻结结构并生成 3Dmol canvas；真实项目无 activation 效应时只显示构建体位点 |
| 渐淡、QC/可信度、截图与导出 | 固定 graph shortest path 且标为 `VISUALIZATION_ONLY`；可靠性/SASA/Myc/QC/override/UNMAPPED 分立，VH/VL 突变轨道独立于抗原坐标；PNG 与结构效应证据 JSON 已从真实项目页面下载并核验；API/浏览器/重启用例通过；activation 关联着色及逐孔回溯仍需真实 activation/QC 数据 |
| Paraplume 2.0.0 全模型/设备 fallback | 同一权重/模型 GPU 和 CPU 路径 provenance；缓存和错误可追踪 |
| 实验与预测证据来源不同 | 预测只标预测；单点/审核主导实验结论；阴性不自动删除候选 |
| GRASP FASTA、IR/RPR 清单、合并和拒绝规则 | 固定 chain order；IR/RPR 分文件；序列不符拒绝；缺证据 IR-only |
| 关闭/重启恢复和重建 | 端到端重启前后结果 hash 对比；原始 FCS 保持只读 |

## 13. 实施约定和风险登记

- **先锁数据模型和文件版本，再做大页面。** 结构图层和约束都依赖稳定 position/analysis revision；过早在 UI 中另算会造成不可复现。
- **先做 fixture，再连接真实任务。** 本机已有用户提供的 260921 Myc 染色 FCS 与 FlowJo 10 WSP，并已只读回放；其 4 个文件不能代替 activation 配套样本或实验补偿矩阵。仍需覆盖更多 WSP gate/transform 类型的 fixture 和多色补偿样例；真实 activation 与实验矩阵作为外部验收输入登记。
- **依赖版本锁定。** FlowKit/FlowJo importer/Paraplume/GRASP serializer 均在 lock 和 manifest 中固定。GRASP 输出规范、FlowJo 未支持语法以及 Paraplume 权重许可/下载由 Phase 2/5/6 的进入条件检查。
- **计算和 UI 分离。** 统计算法为纯函数，接受已审核的 counts/QC/构建体差异/冻结映射；worker 不读写浏览器状态。前端呈现 API 返回对象，不重复实现公式。
- **不把尚未测量的参数硬编码成科学事实。** 最小参考窗口、效应排序阈值、采集 MAD 倍数和事件最低数均作为 schema 参数；初始值按配置显示为工程建议，并需用本实验室重复性/阳性阴性样本前瞻校准。
- **范围递增发布。** Phase 0–4 形成可用的流式实验与结构定位；Phase 5 和 6 是外部模型与文件接口，可作为明确状态的可选能力，不能阻塞或伪装成已有运行结果。

## 14. 本机复现与交接物

以下路径相对本仓库根目录。`.tool_envs/`、`.tool_cache/` 被 `.gitignore` 排除；交接到另一台机器时须重新取得环境和大型权重，不能假定克隆源码后即有模型。当前这些计划和内核文件在工作树中仍有未提交修改/未跟踪文件；后续开发先检查 `git status --short`，避免覆盖既有工作。

| 交接物 | 用途 |
| --- | --- |
| `docs/experimental_surface_mapping_v4.md` | 科学解释边界、孔/批次公式、表面归属与证据等级 |
| `experimental_effects.py`、`tests/test_experimental_effects.py` | 孔/批次纯计算内核，已连接不可变 FlowKit gate-analysis、QC、effect artifact/API 与实验工作台；相关计算、API、浏览器和合成 activation→GRASP 重启用例覆盖。真实 activation 数据已取得，正式门控/QC/效应链的科学验收仍待完成 |
| `joint_design_v4/surface_snapshot.py`、`tests/test_surface_snapshot_v4.py` | v4 表面图冻结格式及哈希测试 |
| `tooling/manifest.json`、`tooling/requirements-flowkit.lock`、`tooling/requirements-paraplume.lock` | 本机依赖版本、权重 commit/摘要、已验证硬件与安装来源 |
| `.tool_cache/paraplume_weight_manifest.json` | 四个 Hugging Face 仓库逐文件摘要及快照 commit；仅在当前机器存在 |
| `.tool_cache/paraplume_smoke/user_example.csv`、`user_example_test_summary.json`、`user_example_paratope_residues.csv`、`user_example_gpu.log`、`cpu-worker-response.json` | 示例 VH/VL 的 GPU/CPU 完整模型运行输入、响应和日志；仅在当前机器存在 |

在仓库根目录用 PowerShell 核验当前已有内容：

```powershell
.\.venv\Scripts\python.exe -m pytest tests\test_experimental_effects.py tests\test_surface_snapshot_v4.py -q
.\.tool_envs\paraplume\Scripts\python.exe -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available(), torch.cuda.get_device_name(0))"
.\.tool_envs\paraplume\Scripts\python.exe -m json.tool .tool_cache\paraplume_weight_manifest.json > $null
```

2026-09-25 的聚焦测试结果为 `7 passed`；pytest 因现有 `.pytest_cache` 写入权限产生一条缓存警告，不影响测试断言。完整 GPU 推理已针对示例运行一次，命令及输出哈希记录在 `tooling/manifest.json`。需要重跑时先配置离线缓存，并用新输出目录保存本次运行：

```powershell
$env:HF_HOME = (Resolve-Path '.tool_cache/hf').Path
$env:HF_HUB_OFFLINE = '1'
$env:HF_HUB_DISABLE_XET = '1'
$env:REQUESTS_CA_BUNDLE = (Resolve-Path '.tool_cache/windows_roots.pem').Path
& '.\.tool_envs\paraplume\Scripts\paraplume-infer.exe' file-to-paratope '.tool_cache/paraplume_smoke/user_example.csv' --large --gpu 0 --emb-proc-size 1 --result-folder '.tool_cache/paraplume_smoke/out-handoff'
```

该示例的预期输出是 VH 120 位、VL 107 位，且每位为有限的 `[0,1]` 概率。GPU 与实际 CPU worker 均已完整运行一次；CUDA OOM 重试路径由注入 OOM 的 worker 单测覆盖，未在硬件上强制制造真实 OOM。用户提供的 260921 FlowJo 10 `.wsp` 与 4 个 FCS 已完成只读回放验收；数据属于 Myc 染色，尚不能代替真实 activation/QC/effect 的科学验收。全仓回归本轮 196 passed（123.87 秒）；合成 activation 到 GRASP 流程已验证两次新 store 重启后的 artifact 回溯与 ZIP 字节一致。Phase 7 的真实 activation 跨重启验收及更多写入中断/崩溃注入仍待完成。


2026-09-27 继续补齐 Phase 3 独立抗体变体效应：新增 `ANTIBODY_EFFECT`，固定同一 WT 抗原并要求同一独立批次中的亲本/变体 Fab 各有 activation 与自身 NC；逐批保存两条 NC 校正窗口及差值，无效窗口保留原始对比但不计入有效汇总，至少 3 个有效批次才给 bootstrap 区间。API 验证固定构建体、亲本配对、批次齐全及 QC；结果不产生抗原结构投影。工作台提供 WT/NC 条件限定、固定 WT 探针选择、每批窗口表和跨批摘要，并注明读数不能说明亲和力、折叠、表面表达或空间位点。纯函数/API 相关测试通过（30 passed）；工作区浏览器 smoke 通过且无 JavaScript 错误；全仓 189 项测试通过，Ruff、两份 JavaScript 语法检查和 `git diff --check` 通过。真实 activation/WSP 与实际外部补偿矩阵仍待用户提供；已有 260921 Myc 数据不替代这些验收输入。

2026-09-27 继续补齐 Phase 5 推理进度：Paraplume worker 在模型加载、5 个完整模型 embedding、分数/来源校验和结果持久化阶段写出实际状态；service 将 worker JSONL 事件同步到任务进度，工作台通过 SSE 展示百分比和中文阶段说明，并在 SSE 不可用时回退轮询。进度不基于 ETA；模型步骤只在对应 embedding 完成后推进。worker 五阶段顺序/完成比例与 service 事件转发单测通过；浏览器使用模拟 SSE 核验大模型请求及 14%/86% 中文进度显示；浏览器工作区 smoke、Ruff、JavaScript 语法检查和 git diff --check 均通过，全仓 192 passed（130.41 秒）。本轮未通过完整 Paraplume 推理任务验证进度链路；既有 GPU/CPU 完整推理验收仍有效。真实硬件 CUDA OOM、部署许可、真实 activation 和真实外部补偿矩阵验收仍待完成。

2026-09-27 补强 Phase 0/7 项目快照失败恢复：若快照目录已从 staging 原子发布、随后 SQLite 登记失败，创建调用现在会移除未登记的 staging/final 目录；若数据库状态无法核实则保留已发布目录，避免误删可能已登记的数据。新增一次性数据库故障注入，确认失败后目录与项目列表均为空，并在新 ExperimentStore 启动后再次确认无残留。针对新失败场景和正常快照重启的聚焦用例 2 passed；全仓回归 192 passed（130.41 秒），Ruff 和 git diff --check 通过。硬进程终止时的恢复及更广泛写入中断仍待后续验收。

2026-09-27 Phase 0/7 硬崩溃恢复验收：测试子进程在快照目录原子发布后、SQLite 项目登记前调用 os._exit(73)。新 ExperimentStore 启动时只恢复 ID/manifest/schema 匹配、必需 artifact 完整且逐文件大小和 SHA-256 全部正确的快照；在 SQLite savepoint 内恢复项目、初始 revision、artifact registry 和 CREATE_RECOVERED 审计事件。无效/损坏目录保留供检查。崩溃恢复、数据库登记异常回滚、普通项目创建/重启聚焦用例 3 passed；全仓 193 passed（127.23 秒），Ruff、JavaScript 语法检查和 git diff --check 通过。其他 artifact 写入中断和真实 activation 全链路验收仍待完成。

2026-09-28 继续补齐结构证据写入恢复：对结构证据第一份文件原子发布后、以及 project_manifest.json 更新后/SQLite commit 前分别用 spawn 子进程 os._exit 注入硬崩溃。启动恢复会从 revision 1 冻结快照和已提交 artifact registry 重建派生 artifact 清单；只有不可变身份/冻结 hash 匹配且所有已登记文件仍通过 SHA-256/大小验证时才修复，并写 MANIFEST_RECOVERED 审计。若清单 identity 被篡改则不自动修复，项目报告 INTEGRITY_ERROR。重试时只复用与来源 SHA/大小一致的未登记结构证据文件，避免重写或覆盖可疑内容。硬退出两种场景均在重启后成功完成冻结；相关恢复、身份完整性和正常路径共 6 passed。全仓 196 passed（123.87 秒）。全仓 Ruff 首次检查发现并修复 14 条导入顺序、未使用导入和冗余 f-string 问题；随后全仓 Ruff、相关 Python 文件 py_compile、JavaScript 语法及 git diff --check 均通过。


2026-09-28 首次使用用户提供的 260903 CTLA4 activation fixture 做 FlowJo 集成验收：3 个板目录各含 54 个 FCS，另有 FlowJo 10.8.1 `.wsp` 和三张 54 孔百分比表。三组逐板回放均为 54/54 样本、每样本 3 层门、0 个未解析节点，累计 486 条 gate count；所有事件纳入 QC，原始 FCS 与 WSP 的 SHA-256 前后不变。发现三板重复使用相同 FlowJo `$FIL` 孔 ID，一次性加载 162 个文件会静默折叠成 54 个样本；已在 FlowKit worker 增加明确失败保护并提示逐板导入，加入真实 fixture 回归。三份 XLS 百分比表按板目录比较时，每孔结果尚未完全一致，最大偏差分别约 2.001、0.903、2.116 个百分点，已向用户核对是否为同一版门控/补偿设置。FlowJo workspace 未提供 compensation matrix，所有 FCS 嵌入的 10×10 `$SPILLOVER` 都是单位阵；没有据此宣称真实荧光补偿通过。三组重算的 162 个样本均为 QC `REVIEW_REQUIRED`，其中包含 293 条上限饱和 REVIEW 和 40 条低门事件 REVIEW；下限聚集 INFO 671 条。未代替实验人员执行 PASS/HOLD。真实构建体/对照映射和 activation effect 尚未登记或分析。聚焦 FlowKit fixture（含 260921）9 passed；全仓 pytest 199 passed（179.17 秒），全仓 Ruff、相关 Python py_compile、`node --check web/experiments.js`、`git diff --check` 和工作台浏览器 smoke（全部检查通过、无 JS 异常）均通过。

2026-09-28 按用户补充的验收标准更新 FlowJo 回放：自动重算值可与手动圈选不同；同批门层级、原始/FlowKit 门型、所用通道、门控通道变换及 WSP 补偿设置一致时，FlowKit 完整事件计数可用于 QC，FlowJo 缓存差异只作可见诊断；未知门、样本无法唯一匹配或同批方法不一致仍阻断计数。手动门控修订与 WSP 未解析门的人工重建功能保留。发现 FlowKit 还会折叠 WSP 内重复的 B3.fcs 等样本名；worker 现在按 FCS `$FIL/$WELLID/$PLATENAME/$SMNO` 唯一匹配记录，在临时 WSP 副本中重命名内部键后回放，原始文件只读。260903 每批 54 个样本均通过同批方法检查；260921 Myc 的缓存差异也按同一规则显示，不再误报为未知门。工作簿 Sheet1 确认三个目录对应的试剂组为 L32/L49、Y159/W202、W212/CD80；用户确认 CTLA4-WT、E1–E7、NC 与各试剂角色，以及三目录为独立批次。FCS/WSP 元数据里的 `$SMNO` 与后两个目录标签不一致，因此试剂/批次归属按用户提供的目录标签和首份表格解释，metadata 仅用于匹配 WSP 样本。每个面板包含 9 个条件×3 个技术重复；自动 QC 162 个样本均需人工 `PASS/HOLD`，目前没有登记人工决定或生成效应。WSP 无补偿矩阵，嵌入 spillover 为 identity，不能据此确认实验补偿；首份表格的 W212 第一个技术重复与 table-3 原始导出不一致，等用户指定权威值。真实原始文件 SHA-256 回放前后不变。最终回归：260903/260921 真实 FlowKit fixture 11 passed；FlowJo QC 导入 API 1 passed；全仓 pytest 201 passed（381.37 秒）；新增 API 同 revision 并发回归 1 passed（两个并发请求一胜一冲突，只提交一个 batch）；全仓 Ruff、相关 Python `py_compile`、`node --check`、`git diff --check` 和工作台浏览器 smoke 均通过，无 JavaScript 异常。

2026-09-28 补齐 FlowJo WSP 发布后的硬崩溃重试：子进程在 WSP artifact 已写入并设为只读、SQLite 登记前调用 `os._exit(77)`；新 store 重启后状态仍为 COMPLETE，重试会校验并复用相同 SHA-256 的只读 WSP，而不会尝试覆盖。不同内容或不安全目标仍拒绝复用。`tests/test_experiment_projects.py` 全模块 29 passed（50.33 秒）；全仓 Ruff、相关 Python `py_compile`、两份 JavaScript 语法检查和 `git diff --check` 均通过。全仓 pytest 的 201 项回归结果早于本次两项 project-store 增量，本轮没有重跑全仓；真实 activation 的项目登记、样本映射及人工 QC→效应链仍未完成。

2026-09-28 将用户确认的标签映射登记到真实 CTLA4 activation 项目 `e0f2c26650f04681`：逐批用 Sheet1 的试剂区块与 FCS 文件名孔位交叉核对 `$FIL`，每批 54/54 孔唯一匹配、无缺失或重复，共登记 162 个技术重复样本。L32、L49、Y159、W202、W212 作为抗体面板；CD80 面板保留为配体阳性对照 `OTHER_CONTROL`；NC 保留独立角色。CTLA4-WT/E1–E7、试剂均以用户确认的临时标签保存，未制造序列或构建体链接；批次分别出现 `TEMPORARY_IDENTITIES_UNVERIFIED`、`WT_IDENTITY_UNVERIFIED` 提示。样本导入页增加临时标签路径，并明确标签不能用于结构归属或效应分析；原有手动门控修订继续可用。临时标签存储/无序列冒充路径测试 2 passed；完整 `tests/test_experiment_projects.py` 32 passed（44.65 秒）。该回归首次触发 Windows `os.kill(pid, 0)` 对已退出服务 PID 的 `SystemError`；任务恢复现改用 `OpenProcess`/`WaitForSingleObject` 核对 Windows PID，并对无法确定的进程保守视为仍存活；新增测试通过。最终 Ruff、相关 Python 编译、`node --check web/experiments.js` 和 `git diff --check` 通过；实验工作台浏览器 smoke 全部通过且无 JavaScript 异常。未自动作出任何 QC 决定或计算效应。W212 table-1/table-3 数值冲突、真实补偿矩阵、序列身份核验及 162 孔人工 QC 仍待解决。

2026-09-28 补齐 FlowJo 百分比结果表导入后的两处硬崩溃恢复：分别在原始表 content-addressed 文件只读发布后、以及 manifest 发布后/SQLite 提交前通过 spawn 子进程 `os._exit` 注入退出。重启后项目保持 COMPLETE 和原 revision；manifest 发布中断会按已提交 artifact registry 修复并记 `MANIFEST_RECOVERED`。重试只在现存原件是普通文件且大小、SHA-256 与再次上传内容一致时复用；校验失败仍拒绝。两种场景聚焦测试 2 passed；Ruff、相关 py_compile 和 `git diff --check` 通过，未重跑全仓。真实项目数据没有被这些临时项目测试改动；activation 人工 QC→效应链仍等待用户审核。

2026-09-28 补充复核并修复跨进程写入边界与历史 FlowJo 样本孔位回退：批次写事务在 revision 读取前获取 SQLite immediate write lock；旧记录没有 `well_id` 时，仅对 `sample_label` 做末尾孔位解析并要求与 FlowJo workspace 唯一匹配。新增合成 FlowJo 百分比 API 闭环，同时暴露并修复百分比效应响应中 `WellRate` 孔数统计的对象访问错误，并把 CD80 排除来源与实际 `POSITIVE_LIGAND_CONTROL` 标签对齐。完整 `tests/test_experiment_projects.py` **35 passed**（77.48 秒）；FlowJo parser/效应单测 21 passed，Ruff、相关 Python `py_compile`、JavaScript 语法及 `git diff --check` 通过。真实项目只读浏览器检查确认 162 个 QC 表单为空、table-3 W212 B5 显示 `21.0 %`、未补偿/分母未知提示在场，网络中没有非 GET 请求且无页面错误。新版服务已在 `127.0.0.1:8766` 重新加载最新代码；旧 `8765` 服务保持不动。真实 CTLA4 项目仍为 revision 16、3 份结果表、0 条 QC 和 0 条效应。全仓测试及用户人工 QC→效应仍未完成。

2026-09-28 用户确认已更新工作簿中的 W212 首个技术重复，使其与 table-3 对应孔的原始导出相符；按修订后的 fixture 断言逐板/逐重复一致，移除先前失效的“仍有冲突”判断。数据入口继续以已圈门百分比表为主路径：未补偿状态如实标注，事件分母未知；FlowKit/WSP 重算仅为可选验证路径。

2026-09-29 修正 FlowJo 同批门控方法比较的成员语义遗漏：一致性签名现在包含矩形/多边形门的 `use_complement`、BooleanGate 的运算类型和引用路径；BooleanGate 无维度时不再触发通道遍历异常，并在 WSP 预览中按 BOOLEAN 标识。FlowKit 的 WSP EllipsoidGate 会被转换为近似 PolygonGate；预览新增明确说明，避免把近似重算说成精确椭圆重放。另加 WSP 预检：带 `transformation-ref` 的椭圆门在构造 FlowKit workspace 前成为 `WAITING_FOR_USER` 节点，包含样本、父路径、通道和原因，避免笼统 `FAILED`。新增只使用合成 workspace 对象和合成 XML 的回归，不读取或改写实验 FCS；该测试 1 passed，Ruff、Python 编译和 JavaScript 语法检查通过。真实项目、人工 QC 和效应记录不变。Phase 2 仍待覆盖更多 WSP 门型/变换并完成完整系统验收。

2026-09-29 扩展 FlowJo WSP 导入预检：锁定版 FlowKit 不支持的 transform 类型、`biex` 非 256 `length`、缺失/非有限参数及整个 Transformations 块缺失，均在建 workspace 前转成 `WAITING_FOR_USER` 诊断，包含样本、变换类型和通道；界面给出表格导入或修复 WSP 的指引。复核发现预检将属性键转为小写后按 `maxRange` 原大小写读取，会误报有效 `biex`；现保留 FlowJo XML 属性大小写，并在合成 XML 中覆盖非 256 长度、缺失/非有限参数、未知类型和受支持 logicle/fasinh。`tests/test_flowkit_worker_method_consistency.py` 1 passed；Ruff、py_compile、JS 语法及 `git diff --check` 通过。真实项目 smoke 只读显示 162 个空白 QC 表单，无写请求；没有读取或改写实验 FCS，真实项目/QC/effect 未更动。

2026-09-29 继续补齐 WSP 门型预检：根据锁定版 FlowKit 1.3.2 转换器，矩形、多边形、椭圆和 Boolean 门可进入预览；QuadrantGate 虽可被 XML 解析器识别，但实际转换路径抛出 NotImplemented，未知门型会在构造 workspace 时失败。现在这两类及缺失/多重 Gate 结构会在 workspace 构造前生成 `WAITING_FOR_USER` 节点，保留样本、门路径、门型和通道；页面提示优先导入用户已圈门结果表。合成 XML 覆盖嵌套 QuadrantGate 和未知 SectorGate，确认不产生门计数；FlowJo worker 定向回归 1 passed，Ruff、py_compile、JS 语法检查通过。未读取 FCS；真实项目、人工 QC 和效应未更动。

2026-09-29 为结果表 QC 提交补充合成浏览器验收：`tests/browser_flowjo_percentage_qc_stub.mjs` 验证 HOLD 缺少理由时浏览器阻止提交；补足理由后 POST 精确携带结果表 ID、样本 ID、审核人、决定、理由和当前 project revision。POST 由浏览器 route stub 截获，重载回显也使用合成审核记录；测试随后直接 GET 真实 QC ledger，确认记录数未变化。浏览器无页面异常；真实项目的 162 个样本仍无开发者代填的决定，未处理 FCS。

2026-09-29 修正已圈门结果表的补偿状态默认值：此前所有新导入表都被写成未补偿，可能把 260903 的用户确认外推到别的批次。现在每份新表可显式登记 `NOT_APPLIED`、`APPLIED` 或 `UNKNOWN`，默认 `UNKNOWN` 并保存用户说明；不从批次、FCS 或其他结果表继承补偿状态。效应 provenance 和限制文字按状态分别说明，表格 UI 使用中性的“导出百分比”而不把未知 metric 统称为 `% of Parent`。现有 260903 artifact 不变，人工 QC 决定不变。API 定向回归 1 passed，覆盖新导入默认 `UNKNOWN` 与用户声明 `APPLIED` 的 provenance/限制；只读真实工作台 smoke 显示 162 个 QC 表单未预填、页面无错误且没有非 GET 请求；Ruff、Python 编译与 JS 语法检查通过。

2026-09-29 按用户确认将当前 Paraplume 使用范围限定为本机非商业研究、不分发，并确认目前没有与 MIT 项目组合的书面许可；用户选择保留工作台集成，仅供本机运行。核对 `cli.py` 拒绝非回环 host，启动脚本绑定 `127.0.0.1`。复核 Paraplume 2.0.0 wheel 内的 `LICENCE`（SHA-256 `47fc1150…3107a0b9`）及上游原文：自定义许可写明 private/academic/institutional 用途、免费再分发的通知要求、商业集成需事先书面授权，以及不得与不兼容许可软件组合；该组合条款未限定为分发场景。仓库根许可证为 MIT，`.tool_envs/` 与 `.tool_cache/` 不随源码提交。工具 manifest、README 和抗体 UI 已记录许可条件与未决点。非商业/不分发范围没有解决本机集成调用是否属于被限制的“组合”；在此澄清前，不宣称集成路径已获许可。

2026-09-29 改进 260903 结果表逐孔人工 QC：每孔保存后不再重建整个项目页面或丢失审核位置；本次 UI 会即时显示审核计数，多个孔的提交串行绑定逐次返回的项目 revision，并阻止同一标签页并发写入。每张结果表新增只填空白行审核人字段的快捷操作，不触碰 QC 状态或发出写请求。若已有描述性效应引用的 QC review ID 与当前审核不同，明确标记旧分析并允许基于最新 QC 创建新的不可变效应版本。合成浏览器回归通过 54 孔连续审核、审核人快捷填充不写入、HOLD 理由必填、revision 逐项递增、旧效应过期提示与重算入口；所有审核 POST 均由 route stub 截获，真实 QC ledger 保持不变。`node --check`、manifest JSON 解析与 `git diff --check` 通过。真实 162 孔仍为空白，未提交或推断任何 QC 决定。

2026-09-29 完成 Paraplume 任务进度端到端验收：在 `.tool_cache` 下创建并自动清理隔离临时 V4 项目，使用 120 aa + 120 aa 的重复氨基酸字母合成输入（无真实抗体身份，不作生物学解释），经 `ParaplumeRunner` 调用完整 Paraplume 2.0.0 large worker。RTX 3060 Laptop GPU / CUDA 12.4 实际完成 IgBert、IgT5、ESM2、AntiBERTy、ProtT5 五个 embedding 阶段；任务记录捕获 15 次单调进度更新，随后完成来源校验和 artifact 持久化，task 为 COMPLETE、run 为 `PREDICTED_UNREVIEWED`，VH/VL 各返回 120 个有限 `[0,1]` 概率。仅使用临时数据库与项目目录，无序列、预测数组或 artifact 写入真实 CTLA4 项目。该运行闭合了实际 worker→service task/progress→持久化路径；真实硬件 CUDA OOM 及 Paraplume 与 MIT 项目本机组合许可范围仍未解决，真实抗体身份不在本轮测试范围。

复核上述 Paraplume 实际进度发现 `llm_to_func` 的字典顺序为 IgT5、IgBert，而模型调用顺序为 IgBert、IgT5，导致第二个 embedding 的百分比在单调保护下停滞。worker 现按锁定模型的运行顺序建立阶段百分比，并在模型 embedding 集不匹配时明确失败；定向回归将函数表顺序故意设为与调用顺序不同，3 项通过，Ruff 与 `py_compile` 通过。重新运行相同类型的临时合成任务后，五阶段结束进度为 30%、44%、58%、72%、86%，均严格递增；task COMPLETE、概率长度和值域检查通过，临时项目自动清理，真实项目未改动。

核对 Phase 6–7 现有验收：`tests/test_grasp_official_format.py` **1 passed**，固定 commit 的官方 parser 对 serializer 生成的 IR/RPR 执行解析、数组往返和 round-trip；不涉及真实抗体序列或 GRASP 推理模型。`tests/test_experiment_projects.py::test_support_bundle_audit_restarts_and_includes_originals_only_on_request` **1 passed**：合成项目的默认 ZIP 包含审计/修订记录并省略原始文件，显式勾选才带入原件；进程重启后的 API 下载字节与原 bundle 一致，项目 revision 不变。没有读写真实实验 FCS/QC；输出仅为临时 fixture。pytest 显示既有 `httpx` 弃用与 `.pytest_cache` 权限警告，断言通过。

2026-09-29 审计结果表主分析路径：`tests/test_flowjo_result_tables.py tests/test_experimental_effects.py` **21 passed**，覆盖百分比单位/精度、QC 门槛及以独立批次为 n 的描述性效应。以 SQLite `mode=ro` 检查真实 CTLA4 项目，并逐份核验三张原始结果表与标准化 JSON 的 SHA-256；再将三份已圈门表重新解析，逐孔比较 54 个孔位的比例、原始显示值、精度和来源行，均与已登记结果完全一致。每表样本 artifact 与重复类型/编号关联完整，范围为有限 0–100 百分比，未推断事件分子/分母，补偿 provenance 为 `NOT_APPLIED`。项目仍为 revision 16 / 162 samples / 0 QC / 0 effects；只读取用户已导出的结果表，未读取 FCS、提交写请求或更改数据库。

真实工作台 `tests/browser_flowjo_percentage_readonly.mjs` 只读 smoke 通过：API route stub 去除所有 WSP 关联并提供独立合成效应用于图表呈现；真实结果表页面仍加载 162 个空白 QC 表单，显示 W212 `21.0 %`、未补偿和分母未知提示，新导入默认补偿为 `UNKNOWN`。合成散点仅验证技术重复逐孔呈现和 `PASS_WITH_OVERRIDE/HOLD` 标记，不是 CTLA4 真实效应；浏览器检查确认零写请求、零页面异常，没有读取 FCS。

2026-09-29 修正 CSV/TSV 文本百分比的小数逗号歧义：旧 parser 无条件移除逗号，可能把 `2,1%` 静默改成 `21%`；还会移除数值内部空白，把格式错误的 `2 1%` 改成 `21%`。导入界面/API 现在要求点或逗号小数分隔符；未声明逗号、分隔符混用和数值内部空白均拒绝，只允许数值两侧及 `%` 前的空格。原显示文本、数字尺度与所选分隔符保存到 provenance，Excel 数值仍以单元格 number format 为准。新增 semicolon CSV 端到端 API 验收确认 `2,1%`→2.1%、单位/精度原样保留、QC→描述性效应数值一致及非法选项不改变 revision。解析器与导入 API 聚焦测试 **14 passed**；Ruff、py_compile、两份 JS 语法和实际结果表页面只读 smoke 通过，后者确认新下拉选项可选、162 个 QC 留空、零写请求。升级后再次只读重解析真实三张 XLS：每表 54 个孔位、原表和标准化结果 SHA-256 及所有原显示值均吻合，项目为 revision 16 / 162 samples / 0 QC / 0 effects。没有读取 FCS；真实 CTLA4 artifact/QC/effect 未更动。

2026-09-29 当前状态复核：全仓 `pytest tests -q` **252 passed**（626.32 秒，6 条依赖/缓存警告）。本机 `127.0.0.1:8766` API 对 CTLA4 activation 项目返回 HTTP 200；只读 SQLite 核对为 revision 16、162 samples、3 份百分比结果表、0 条人工 QC、0 条效应。六个原始结果表/标准化结果 artifacts 的文件均存在且 SHA-256 与 registry 相符。`tests/browser_flowjo_percentage_readonly.mjs` 对真实项目页面 smoke 通过：162 个 QC 表单均为空，W212 B5 显示 `21.0 %`，未补偿与事件分母未知提示可见；WSP 及 effect API 用 route stub 隔离，浏览器零写请求、零页面异常。测试与 smoke 均未写入真实项目 QC/effect，也未读取 FCS。真实 activation 的人工 QC→效应、E1–E7 序列/抗体身份、Paraplume 许可兼容范围及有复合物证据后的 RPR 仍待对应输入或授权；本次验收不据此标记目标完成。
