# 项目记忆：V4.1 实现、故障与维护约束

更新日期：2026-10-01。用于后续开发查阅；历史运行是证据，不是当前算法规范。
本文件记录已经核对的修复与已知边界，不能据此宣称所有输入都能产生正式推荐。

## 1. 当前架构与术语

```text
原始结构 + UniProt 映射 + AF3 confidence（预测输入）
    → native SASA / 结构可靠性
    → Surface V2：点表面、重叠足迹、Surface V2 区域
    → FootprintFamily → InferenceHypothesis + 全局独立双点/三点候选搜索
    → SPURS single / multi、FoldX 确认与 formal gates
    → hypothesis response matrix → identifiability panel
    → schema 4.0 兼容投影 + 4.1 审计文件
```

- **显示区域**：互斥颜色归属，必须空间连通且直径受限。不是已验证表位。
- **科学假设**：Surface V2 raw footprint 经 family 与 conservative inference hierarchy 解释；同一残基可参与多个假设。
- **扰动矩阵**：几何与化学规则的 proxy，不是结合下降概率。
- **面板**：greedy 信息增益与风险约束的结果，不宣称全局最优。
- 区域数、假设数、构建体预算相互独立。固定 assay profile V2 的 Stage 1
  上限为 6；未指定固定 profile 的旧请求可保留其独立预算。
- 膜方向、糖链、PAE 是结构证据；不能据此声称预测 reporter 或亲和力。

## 2. 本次 PD1 / AF3 分区故障

输入：`workspace_data/fa85119485cd`，Q15116，ECD 25–170。
原模型 SHA256：
`8e4c7e0b922377d19359d2f4ad1bd598dadb69d0cdcb3b991aa1c5ea07b7c496`。

### 原因与证据

0. **页面着色的直接原因是导出编号错位。** ECD CIF 的 label_seq_id 从 1
   开始，而 record 保留 AF3 原编号 25–170；旧选择器优先用 label，导致
   本输入颜色选择偏移 24 位。实际导出前两项 auth/label 为 25/1、26/2。
   输入映射本身正确并不意味着导出结构的 label 编号仍一致。
1. `sectors.py` 为接近 7 区而合并小区域，未在合并前检查直径；最终检查只
   添加警告，没有修复。原 VS004 测地 / Cα 直径为 42.410 / 35.499 Å；
   VS007 为 39.432 / 34.631 Å，超出当时 35 / 30 Å 上限。
2. 原 VS001 的 25–31 与 129–130 在序列上相隔远，但 Cα 最大距离为
   16.821 Å；不能把序列跨度直接当作三维错分。原 VS002 Cα 直径为
   25.029 Å，显示范围也明显大于页面宣称的 ≤18 Å。
3. AF3 末端 25–31、150–170 的低置信度坐标没有达到严格 IDR 分类条件，
   因而仍进入三维几何。这混淆了“不能诊断无序”和“坐标足够可靠”。
   原 reliability QC 的 excluded_positions 为空，88 个残基进入表面图。
4. UI 的固定文案称所有区域都是 13 Å / 18 Å Fab 补丁，实际显示的是 V4
   sector。旧颜色只有 7 种循环复用；区域数增加时会造成远处分区同色。
5. 重放脚本把所有结构设为 `structure`；用于 AF3 会丢失预测模型语义和
   confidence 的作用。必须保留 `af3_upload` 与原 confidence。

### 修复规则

- `web/app.js.residueQuery` 使用导出后保持一致的 auth chain、auth_seq_id
  与 insertion code；只有 author 标识缺失时才回退到 label。不能直接
  将原 AF3 label 用在经过 ECD 截取 / 写出的 CIF 上。
- 显示 regions 时，不叠加未归属位点的通用 SASA 蓝色或 folding 绿色；
  不可靠尾部为灰色。否则通用蓝色会被误认为区域 1。
- `web/renderer.js` 等待实际结构 hierarchy 与 assembly 可用才首次着色，
  不能仅等待 managers 初始化，防止第一次 paint 在加载完成前被跳过。
- `structure_reliability.py` 的 `reliability_rules_v2` 将坐标资格与生物学
  无序分类分开。ECD 端部连续 UNCERTAIN 段只有在以下证据同时满足时，
  才退出三维足迹：所有已测 pLDDT <70、至少 80% 位点非局部接触数 ≤2、
  至少一个位点有经过合法 token 映射的 core PAE ≥15 Å。
- 该段仍为 UNCERTAIN，保留原结构、序列和审核能力；不自动标成 IDR，
  不产生自动 IDR 突变。高置信度移动结构域和嵌在核心中的低置信度短环
  不因这条规则退出。缺少 PAE 不能触发这条退出规则。
- `surface_patch/engine.py` 真正消费 `use_for_3d_footprint`；不可靠末端
  不生成表面样本，仍作为当前构象的遮挡原子保留。既有明确 DISORDERED
  段继续遵守 strict / occlusion_sensitivity 策略。
- 相邻 2 个保留位点作为边界 guard；不可靠末端自身不重复算作边界。
- `sectors.py` 在输出前拆分不连通 / 超径组，合并只接受满足全部空间约束
  的并集。默认显示直径：测地 ≤26 Å、Cα ≤18 Å。这是明确可配置的显示
  紧凑性规则，不是所有抗体足迹的物理尺寸定律。
- 区域数是软目标。不能为了 7 色牺牲空间约束；小区域若无法安全合并，
  保留并标记。额外区域采用不同颜色。
- `surface_display_regions_v41.json.regions` 保存实际 V4 显示归属，
  `surface_v2_regions` 保存原 Surface V2 聚合区域，两者不能混称。
- 兼容投影的 `mutable_residues` 不再包含硬保护残基；空列表不能在 UI 中
  被 `||` 回退成“全部可突变”。
- V4 兼容投影同步 `ResidueRecord.region_id`，使珠子视图与 Molstar 的区域
  membership 着色一致；不属于有效区域的位点必须清空旧归属。

### 验证资产

- `tests/fixtures/pd1_af3_partition.json`：原坐标、PAE、原图和原分区；
  PAE 按明确 token 映射截取 ECD，不用数组下标猜测 UniProt 编号。
- `tests/test_pd1_af3_partition_regression.py`：原故障图、端部资格、缺失
  PAE、实验结构、短环保留、mask 到 footprint/graph 的传递、软计数、
  颜色唯一和 AF3 重放语义；实际执行 JS 选择器并验证 ECD 导出编号。
- `workspace_data/pd1_af3_partition_repair_20260930/verified_geometry_complete`：
  原模型与 confidence 的真实几何重算；外部评分显式禁用，无正式突变推荐。
- 新图保留 58 个符合表面与可靠性规则的节点，仍形成 7 个有效区域；
  最大测地直径 25.226 Å，最大 Cα 直径 17.885 Å。
- 28 个端部位置退出三维足迹；边界 guard 为 32、33、148、149；全部
  146 个 ECD 坐标和原输入结构保留。当前输入未被诊断为 IDR。
- 页面验证副本任务 `da27bc835241`，明确标为“仅重算几何，未运行突变评分”；
  原失败任务保留。172 个假设均无不合资格位点；47 个 bundle 文件一致。

详见 `docs/pd1_af3_partition_repair_20260930.md`。

## 3. 之前的 V4 → V4.1 实现问题

| 原问题 | 已实现的修复 / 维护要求 |
|---|---|
| Surface V2 已计算但只传 `.surface` | 保留完整结果，将重叠 footprint/region 传入 V4 adapter；不要另造科学分区引擎 |
| sector 被当作科学表位 | 显示与科学 membership 分离；schema 4.0 是兼容投影，4.1 是扩展审计 |
| 共享 beam 挤占双/三点预算 | 每 order 独立 beam 和 shortlist；四点为探索池 |
| option.score 混合所有证据 | 保留暴露、化学、prior、稳定性和风险分量；内部排序分数不是校准概率 |
| panel 按 sector 轮询 | 4.1 使用所有假设的 signature、覆盖、风险、冗余与信息增益 |
| 只报告平均区分 | 同时报告最小距离、未区分组、未覆盖假设；不得掩盖 minimum=0 |
| DiscoTope 缺失导致整体阻断 | generic epitope prior 是可选软先验；保留 stability fallback |
| 将 ProCaliper 当抗体可达性 | 当前 worker 是 residue SASA cross-check；不是 Fab 尺寸或进入凹槽的模拟 |
| 外部 SASA 方法参数缺失阻断所有候选 | native SASA + 显式几何为主；ProCaliper 方法差异为诊断，真实 native/normalized 分类冲突仍审查 |
| `_formal_candidate_reasons` 与测试参数漂移 | 兼容参数与契约测试同步；不要简单删掉失败断言 |
| ProCaliper atom radii provenance 缺失 | worker 与 context 保存物理设置 / 来源；不能编造 matched provenance |
| selected=[] 被误报为 FoldX 未运行 | 分开工具状态、候选覆盖、面板覆盖，成功评分 42 个不能报 missing |
| PAE / glycan / membrane 信息被丢弃 | 接入已有结构可靠性与 geometry evidence；实验 B factor 不当作 pLDDT |
| 几何均值访问 footprint 外残基导致 KeyError | 只对权重 membership 与有效 records 的交集累加 |
| 用声明的模型 revision 代替真实缓存身份 | 权重与模型文件 hash 必须由实际输入建立；不能无依据重用旧评分 |

## 4. CTLA4 无推荐故障与修复

详见 `docs/ctla4_v41_repair_20260930.md`；关键证据如下：

| 原问题 | 原因、修复与回归 |
|---|---|
| 新 CTLA4 无 single / DiscoTope | `.venv/python.exe` 为 Low Mandatory Level，WSL 返回 `Wsl/E_ACCESSDENIED`；仅本机解释器标签恢复 Medium |
| 桥接日志只有空 RuntimeError | WSL 错误在 UTF-16LE stdout；桥接解码 stdout + stderr 并返回结构化错误 |
| 含空格路径 / shell 参数问题 | WSL 使用 argv 与 `--exec`；不拼接不受控 shell 字符串 |
| ProCaliper 导入错误 | `Select` 从 `Bio.PDB` 导入，不能从不存在的 `Bio.PDB.Select` 导入 |
| 同组合来自不同显示区域，multi 整批拒绝 | 全局按位点与替换去重，保留多 hypothesis 归属；推理边界再防重 |
| 空候选被当工具失败而无原因 | 保存 NO_CANDIDATES_TO_SCORE / NO_CONFIRMABLE_CANDIDATES 与阶段诊断 |
| native normalized SASA 依赖可选 ProCaliper 成功 | 始终运行 native 对照；可选工具故障不跳过结构语境审计 |
| FoldX 预算被不可确认 / 四点候选占用 | 先做 SPURS multi，主确认预算只给通过该门的双/三点候选 |
| 同一位点组合反复跑图距离 | 复用窗口内单位置距离及位置集合覆盖，氨基酸菜单不改变空间距离 |
| 完成消息覆盖推荐阻断原因 | pipeline、diagnostics、jobs 保留真实原因；计算完成不等于正式推荐 |
| ZIP 漏掉 V4.1 文件 | 4.1 文件纳入 manifest 与 ZIP；最终 analysis 写入后重新封装 / 哈希 |

真实 CTLA4 修复重跑：`8c068ad701a6` →
`workspace_data/ctla4_v41_repair_20260930/verified_pipeline`。
210 个唯一组合、34 个正式通过、5 个正式入选；所有外部 worker 实际成功。
48 个导出文件的磁盘 / manifest / ZIP 一致。页面副本任务为 `05ecca7e4cc7`。
这证明特定实现故障修复，不证明其他输入必须产生组合。

本机完整性标签修复不随 Git 同步；程序不自动提升权限或重置整个目录 ACL。
长时间 multi 推理不要反复启动重复进程；先检查进度与模型实际日志。

## 5. 后续必须执行的检查

1. 先确认原始输入类型、chain/auth/label/UniProt 映射、confidence token
   映射和文件 hash。PAE 无合法映射即忽略，不能凭序列下标猜测。
2. 显示 membership 对 eligible graph **完整、唯一**；每区 induced graph
   连通，所有两两距离满足直径上限；不能只打 quality flag 后照常输出。
3. geometry-ineligible 的样本不能经 graph fallback、transit contraction、
   sequence bridge 或 hypothesis adapter 被重新引入。
4. 用真实故障 fixture 验证，再核对完整输入重算；历史运行不可原地改写。
5. 不因“必须推荐”而放宽稳定性、保护位点或外部证据门槛。有效空结果应
   明确区分真实不通过、工具缺失、几何资格不足与优化未获得信息增益。
6. 只更新文案是不够的；反之修后端也必须核对 UI 归属、颜色、位点计数、
   JSON、CSV、manifest、离线报告和下载 ZIP。
7. 当前工作树包含大量其他功能的未提交编辑，不能 `git add .` 混合提交。

### 测试与重放

从仓库父目录 `E:\user\New project` 运行，避免包导入错误：

```powershell
python -m pytest -p no:cacheprovider ecd_epitope_scanner/tests/test_joint_v4.py ecd_epitope_scanner/tests/test_joint_v4_sasa_context.py ecd_epitope_scanner/tests/test_surface_snapshot_v4.py ecd_epitope_scanner/tests/test_surface_patch_v2.py ecd_epitope_scanner/tests/test_structure_reliability.py ecd_epitope_scanner/tests/test_v4_hypothesis_panel.py ecd_epitope_scanner/tests/test_ctla4_v41_regression.py ecd_epitope_scanner/tests/test_pd1_af3_partition_regression.py -q
```

从仓库目录验证原 AF3 几何，不运行评分：

```powershell
.venv/Scripts/python.exe scripts/replay_joint_v41_run.py workspace_data/fa85119485cd workspace_data/new_pd1_geometry_validation --geometry-only
```

去掉 `--geometry-only` 才会运行真实评分、候选门槛和面板。
重放必须使用新目录。`replay_input_provenance.json` 记录实际计算范围与 hash。

## 6. 仍需谨慎的边界

- 缺少 PAE 的预测输入不会由这条 AF3 末端规则自动删除；应保留不确定性，
  不能称作已验证的表位几何。
- 高置信度但相对方向不确定的结构域仍需 PAE / 多构象审查；本次 PD1
  修复不等同于完整结构 ensemble 验证。
- 低置信度末端保留为遮挡原子，因此“移走尾部”仍可改变表面可见性；
  本次采用保守的给定构象，不声称对此类变化严格不变。
- 受保护 / 埋藏残基造成的颜色空隙不等同于空间断裂；连接判定基于已
  记录的显式表面图及有界 transit 路径，不是序列连续性。
- 70 项相关回归通过不等于整个仓库所有实验工作台 / 流式功能均已验证。
- V4.2 的 ensemble / retrospective benchmark、V5 的训练不是本次已完成
  项目；没有真实完成的数据不能写成已实现结果。


## 7. 2026-10-01：V4.1 下一阶段资格、层级和自适应预算

完整实施与证据见 [v411_upgrade_20260930.md](v411_upgrade_20260930.md)。

### 原因和政策变化

- 原 surface gate 将 total rSASA 与统一20 Å²侧链面积串联，误伤 Gly/小侧链；现在 native total rSASA 默认0.25决定暴露，侧链只描述。0.15/0.20比较不会自动改变默认。
- 低 pLDDT 表示坐标不确定，不能解释为突变必然破坏折叠。内部几何可用环区生成复核候选，严格面板拒绝复核候选；不可靠 AF3 末端仍退出图。具体替换继续受 Gly/Pro、二硫键、修饰、糖基化、SPURS/FoldX 政策约束。
- 显示 ownership、raw footprint、family、inference hypothesis 和构建体数独立。9/12/15工程权重1/0.7/0.3，adaptive独立保留；complete-link保守合并，不跨组件/域或不兼容可靠性。缺少合法预测 PAE 必须 UNAVAILABLE，不能以 None→0 伪装为LOW。
- 候选全局物理去重，来源假设合并；每order独立shortlist，60默认可为保留配额扩容，空pool也记录配额缺口。复核不占正式预算；四点仍探索。
- V3默认自适应4/7/8，固定√8距离；4只是目标，没有正净增益允许更早停止。4以后绝对信息增益0.02、相对0.05及分离改善条件必须满足。V1/V2原固定策略保留，新参数优先于旧max。最小距离为0不能否决其他有效分离改善；近距离连通组仅为待解析组，不代表全组两两不可区分。

### 回归与交付注意

- CTLA4 V3 strict 0.25：245唯一组合、26正式通过、3入选，worker无失败；0.15：258/28/3。均因 NO_POSITIVE_NET_INFORMATION_GAIN 提前停止，不能强补4个。
- PD1原AF3几何：62节点、7显示区，连通且26/18 Å受限；32个几何不合格位置不入图。该重放不运行评分，不能称作真实突变推荐验证。
- 必须核对最终analysis、manifest、ZIP哈希，避免先打包再改analysis；raw/family/inference映射、复核池、候选缺口、逐步预算曲线均需入包。CTLA4 bundle61文件一致，浏览器只读验证无JS异常。
- 除工作树测试，还需从Git独立检出测试。本次发现旧pipeline已有infer_antigen_chain引用但实现未提交，必要结构分析依赖补入单独提交。避免只在脏工作树可运行的交付。
- 使用项目venv运行包含API的测试，系统Anaconda缺FastAPI；从父目录执行。Pydantic显式要求2.x，与model_validator一致。
- 服务进程加载旧模块时网页新资产不等于后端已升级。本次自动审批拒绝服务重启，需用户正常重启；不要绕过拒绝或停止未知任务。

### 剩余限制

- surface centers约7 Å、合并上限6 Å，保守层级仍保留大量hypotheses（CTLA4 125、PD1 94）；不强制压缩，定位收益待真实benchmark。
- benchmark目前仅输入接口、隐藏真值隔离和指标规范，未完成真实抗体回顾性评估；ensemble/noise/training后置。
- 真实重放导出时采样源码hash，运行期间后续编辑不能冒充已加载代码；保留启动快照及独立检出证据。

- CTLA4三阈值比较已完整完成：0.15/0.20/0.25唯一候选258/250/245，正式通过28/28/26，均入选3，worker无失败，均NO_POSITIVE_NET_INFORMATION_GAIN。0.20新增9个表面位点，0.15新增12个；入选组合发生改变，默认保持0.25。0.20显示13区仍连通有界，61个bundle文件一致；区数非构建体预算。

- 最终独立Git检出回归 **97 passed / 33.07 s**，一个Starlette/httpx弃用提示；专用升级模块23项。Python编译、JavaScript语法检查及任务提交diff检查通过。
