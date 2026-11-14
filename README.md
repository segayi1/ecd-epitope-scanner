# ECD Epitope Workbench

本机单用户的单次跨膜蛋白胞外区（ECD）结构分区与组合突变设计工具。它生成可审计的抗原突变构建，不依据抗体序列预测表位，也不把计算筛选误写成实验表位结论。服务仅允许本机回环地址；公开部署需要未来的专用运行配置、认证与作业隔离，不能通过修改 `--host` 直接开启。

## 已实现的科学边界

新建分析统一使用 `joint_v4_ml`。历史版本的实现说明保留作研究记录，旧引擎不再作为新作业的可选入口。

- 输入 UniProt accession、PDB ID，或从相关网站下载的 PDB/mmCIF 文件。
- UniProt 输入自动下载 canonical 序列、拓扑、结构域、二硫键、糖基化和功能注释，并通过 3D-Beacons/PDBe 搜索实验结构。
- UniProt 搜索得到的结构候选仍按完整 ECD 质量门槛筛选：观察覆盖率 >=95%、最大内部缺口 <=8 aa、与 canonical 序列一致性 >=95%、主要结构域覆盖 >=90%。直接输入 PDB ID 或上传 PDB/mmCIF 时，这些门槛改为结果中的准入警告，不会单独终止分析；逐残基编号和警告原因会保存在结果中。
- PDB ID 优先使用 SIFTS 指向目标 UniProt 的链；上传结构在有 accession 时按序列映射选择目标链。无 accession 的抗原–抗体复合物会结合结构实体名称、免疫球蛋白序列特征和链间接触自动选择抗原链，并记录方法、证据与不确定性。支持显式指定链覆盖自动选择。检测到 PDB 与 UniProt 标准序列不一致的位置时，原结构保留为上下文，差异位点及警告写入结果，这些位点不参与逐残基评分、表面分区和突变设计。
- 没有合格实验结构时，只接受带明确 AF3 来源的公开模型；否则状态为 `WAITING_FOR_AF3_UPLOAD`，网页要求上传官方 AF3 mmCIF 与 confidences JSON。AFDB 不会被伪装为 AF3。
- 统一用 Bio.PDB Shrake-Rupley 与固定原子半径、Tien 最大 ASA 归一化计算 rSASA；默认表面标准为 `total rSASA >=0.25`，侧链 SASA 仅作描述，不因 Gly 零侧链面积排除表面资格。
- 表面、三维坐标资格和突变审核分别记录。二硫键、糖基化/PTM、缺失坐标与埋藏位点维持保护政策；低置信度内部环区可产生 `REVIEW_REQUIRED` 候选，不进入严格正式面板。不可靠 AF3 末端仍退出三维几何。
- 新任务默认 assay profile V3：多尺度 footprint family → inference hypothesis，全局假设均衡候选池，自适应 4–8 构建体，常见目标 7。最少数量不强制填充；预算、停止原因和未解析假设均进入审计。旧 V1/V2 profile 保留原预算。
- 以 Fab 足迹尺度（测地半径 13 Å，约 15-25 残基，Kringelum 2012）在表面测地图上枚举补丁并允许跨结构域；每个补丁带表位先验（暴露度/外凸度/环区占比/进化可变性/糖屏蔽/膜近端惩罚）并分级 HIGH/MEDIUM/LOW；面板预算随 ECD 大小自适应（~1 补丁/22 残基，上限 7），结构不支持时如实输出更少补丁与原因。
- 实验形式 `assay_format=cell_surface|purified` 是一等参数：细胞表面实验时膜近端补丁先验降权并标记 IgG 不可及；`purified` 时全部保留。
- 每个补丁按 WT 残基类型定向选择替换字母表（大侧链→Ala、小侧链→Trp 插入、带电→电荷反转），Stage-1 探针必须含至少一个锚点突变（电荷反转或 W/Y/F/R/K→A）；在 FoldX/proxy 清空的候选菜单上**全枚举** 3-6 点组合（两两间距 ≥5.5 Å、最大展开 ≤22 Å），按**帕累托前沿**（破坏势 vs 稳定性风险）选点：主探针取膝点（单位风险破坏力最大），备选探针取激进端，完整前沿写入 `epitope_mapping_plan.pareto_fronts` 审计；无满足约束的组合则输出 `NO_LOW_RISK_COMBINATION`。
- 每个探针自动生成 Stage-2 预设：单点去卷积构建体（`STAGE2_SINGLE`）与确认性双点组合（`STAGE2_PAIRWISE`），并输出面板覆盖矩阵、决策规则与显式盲区清单（`epitope_mapping_plan`）。
- FoldX 单点门限与暴露度联动（rSASA≥0.35 → +0.80；0.25-0.35 → +0.60 kcal/mol；文献 RMSE 0.85-1.1），组合门限维持 +2.0 kcal/mol。
- 可选证据线：UniProt 同源保守性（自动，离线优雅降级）与 ESM zero-shot LLR（设置 `ECD_ESM_MODEL` 且本机装有 torch/fair-esm 时启用，LLR<-3 共识否决）。
- 输出仅含 ECD 的氨基酸序列和编码片段，不自动添加信号肽、标签、起始/终止密码子、TM 或胞内区；`original_cds` 可解析为全长 CDS 时额外导出**全长**膜蛋白突变编码序列（细胞流式实验需要），无 CDS 时导出全长蛋白序列供基因合成反向翻译。无原始 CDS 时使用 DNA Chisel 的人源/HEK293 密码子优化。

### Joint V3 历史实现（不用于新作业）

旧版 `joint_v3` 使用“结构审计 + 残基级表面单元 + 多尺度 micro-footprint + 连通 ProbeSector + 联合面板优化”流程：

- `ProbeSector` 是用于少量探针粗定位的连通表面区域，不是真实抗体 footprint；`MicroFootprint` 只用于多个面积尺度的敏感性分析。
- 候选探针固定为 2–4 个突变位点，采用确定性 beam search，不完整枚举全局组合，也不强制加入 anchor mutation。
- 面板优先使用 SciPy/HiGHS MILP 最大化多个测地覆盖情景下的稳健可设计表面积；求解器不可用或超时则降级为确定性 greedy，并标记 `HEURISTIC_SOLUTION`。
- 二硫键 Cys、糖基化 sequon、结构缺失、明显结构风险等进入保护/排除审计；不可设计表面仍保留在覆盖分母中，不会被静默删除。
- 没有外部组合级稳定性引擎时只输出 `RULE_ONLY` 和无量纲风险代理，不填充 kcal/mol，也不输出结合概率。
- 孤立 ECD 只提供 C 端锚点和相对方向；没有真实 TM 坐标时不推断唯一膜平面。Myc/表达信号也不能单独证明 ECD 正确折叠。

V3 的详细实现和解释边界见 [`docs/joint_v3.md`](docs/joint_v3.md)。

Joint V4-ML 的独立方法学见 [`docs/method_joint_v4.md`](docs/method_joint_v4.md)。

实验阳性率下降到结构表面的归属规则与统计边界见 [`docs/experimental_surface_mapping_v4.md`](docs/experimental_surface_mapping_v4.md)。

从实验项目与 FCS 到结构证据、Paraplume 和 GRASP 输入的完整实施路线见 [`docs/experimental_workbench_v4_implementation_plan.md`](docs/experimental_workbench_v4_implementation_plan.md)。

实验工作区中的 Paraplume 2.0.0 安装在 `.tool_envs/paraplume`，不会随源码包分发。用户确认当前仅作本机非商业研究且不分发；尚无 Paraplume 与 MIT 项目组合的书面许可。上游自定义许可含“不与不兼容许可证软件组合”的限制，未说明该限制是否只针对分发，因此本机集成调用的许可适用范围也未确认。商业产品或服务集成需权利方事先书面授权；许可及已安装文本摘要见 [`tooling/manifest.json`](tooling/manifest.json)。

整个系统的入口、任务队列、结构分析、V4 证据线、分区、候选门控和导出关系见 [`docs/system_flowchart_v4.md`](docs/system_flowchart_v4.md)。

Joint V4-ML 4.1（SPURS single/multi、DiscoTope-3.0、FoldX 和可选 ProCaliper SASA cross-check）是新作业的默认后端。V4.1 把 Surface V2 的重叠 footprints 作为候选表位兼容表面假设，互斥 sector 只用于兼容显示与候选组织；探针—假设响应是测地距离与化学扰动构成的规则代理量，不是结合变化概率。DiscoTope 是软先验，缺失时不单独阻断推荐。ProCaliper 当前 worker 只返回 SASA/rSASA，用于与 native SASA 交叉审计；它没有计算有限尺寸探针的空间可达性。命令可写入配置，或通过环境变量 `ECD_V4_PROCALIPER_COMMAND`、`ECD_V4_SPURS_SINGLE_COMMAND`、`ECD_V4_SPURS_MULTI_COMMAND`、`ECD_V4_DISCOTOPE_COMMAND`、`ECD_V4_FOLDX_COMMAND` 设置，环境变量值必须是 JSON argv 数组。为兼容旧配置，`ECD_V4_SPURS_COMMAND` 可同时作为 single/multi 的后备值。

V4.1 保留确定性 V4 sector 作为显示归属，同时单独导出 Surface V2 footprint 假设成员关系；一个残基可同时属于多个假设。偏好的 12 Å 与 adaptive-area footprints 构成主假设池，同中心 9 Å 区域标为 core、9–15 Å 外缘标为 boundary。输出另含 geometry-only `surface_graph_v41.json`、显示区域、假设成员表、探针—假设矩阵与面板可辨识性摘要。旧 `surface_graph_v4.json` 继续保留 sector 字段以兼容实验项目。

V4 默认配置声明以下推理与确认环境，系统重装后需要分别检查可用性。SPURS single/multi 使用 WSL2 环境 `/home/segayi1/miniconda3/envs/spurs_v4`，DiscoTope-3.0 使用 `/home/segayi1/miniconda3/envs/discotope_v4`，通过 CUDA 12.4 在 NVIDIA RTX 3060 Laptop GPU 上推理；ProCaliper 0.2.1 使用 `G:\anaconda\envs\ecd_procaliper\python.exe` 在 CPU 上计算 SASA/rSASA cross-check；FoldX 5.1 使用 `vendor/foldx/foldx_1_20270131.exe`，执行 `RepairPDB + BuildModel`，修复结构跨作业缓存，每个候选默认运行 5 次并取中位数。FoldX 对确定性短名单最多 60 个候选启用 8 路 CPU 并行。所有 WSL GPU worker 使用跨进程锁串行占用显卡；CUDA OOM 时清理缓存后重试 1 次，仍失败则按配置显式回退到 CPU，并在 warning/provenance 中记录实际设备和回退原因。SPURS single 与 DiscoTope 的超时为 900 秒，SPURS multi 为 1800 秒。历史本机运行曾测得 290 aa PD-L1 的 SPURS single 约 58 秒、SPURS multi 单组合约 60 秒、DiscoTope 约 39 秒；另一次 215 aa strict V4 中 ProCaliper 约 2.28 秒，FoldX 的 42 候选 × 5 次重复约 17.96 分钟。这些数字是旧版本机量级，不是 V4.1 性能基准。DiscoTope 的 GPU 路径限短于 1000 aa 的预测链。修改配置或 worker 后需要重启网页服务；`tool_inventory.json` 和各 worker 返回的 provenance 是每次作业实际环境的审计依据。

V4 示例：

```powershell
.\.venv\Scripts\python.exe -m ecd_epitope_scanner analyze `
  --structure .\input\model.pdb --accession Q9NZQ7 --chain A `
  --design-engine joint_v4_ml --ml-mode preview `
  --max-sectors 12 --max-primary-constructs 7 --output .\results\joint_v4
```

`strict` 模式只选择单点证据合格且组合级共识通过的候选；缺少稳定性必需证据、模型分歧或没有候选通过门控时会输出 `BLOCKED_EXTERNAL_EVIDENCE`。DiscoTope 缺失会写入 warning 并使候选全部归入 `STABILITY_FALLBACK`，但不再单独阻断正式推荐。`preview` 允许调试和审查，但其设计状态为 `PREVIEW_CANDIDATE`，不会被旧 UI/汇总计为正式推荐。V4 原生结果写入 `v4_manifest.json`、`tool_inventory.json`、`surface_evidence.csv`、`discotope_residues.csv`、`candidate_combinations.csv.gz`、`selected_panel_v4.json` 和 `v4_scientific_report.md`，并同时生成旧 UI 可读取的 `regions.json`/`designs.json` 投影。

V4 worker 使用 stdin/stdout JSON 协议：程序接收一个 `PredictionRequest` JSON，stdout 必须只有一个 JSON 对象，至少包含 `status`（`success`、`partial` 或 `failed`）和 `records`；日志只能写入 stderr。请求中的预测链是连续编号的规范化目标链，`residue_mapping.csv` 保存 `normalized_position`、UniProt 位置、author 编号、插入码和 WT 氨基酸之间的映射。若 ProCaliper 或 DiscoTope 返回成功，其结果必须覆盖整条规范化链；ProCaliper 缺失或 provenance 不全只影响 SASA cross-check 状态。SPURS/FoldX 的组合结果必须只返回请求过的突变组合，并通过位置范围及 WT 校验。

严格模式下，正式入选组合须有合格的单点 SPURS 证据、有限的组合 SPURS 与 FoldX 评分，并通过稳定性共识及其他候选门；SPURS-only 结果不能替代 FoldX。ProCaliper 是可选 SASA cross-check；方法 provenance 不全、SASA 方法设置不同或 ProCaliper 不可用时只写诊断提示，native normalized-context SASA 与 Surface V2 几何负责 surface eligibility。DiscoTope 是软先验，缺失不单独阻断。四点组合保留为 exploratory candidates，不进入主面板；三点组合须满足所选协议的高阶增益门。

计算预筛只能说明“低预测风险”。进行表位归因时仍需同时验证表达/膜表面表达、构象对照抗体或折叠控制、稳定性，以及必要时的配体或功能读出。
V4.1 保留 ECD anchor-axis 膜相对几何代理，但没有模拟实际膜平面、细胞间距或特定 Fab 的有限尺寸 steric reach；`cell_surface` 结果会显式提示这些未知量。SASA 只描述溶剂暴露，ProCaliper 当前只提供 SASA/rSASA 交叉审计，不能据此宣称抗体可达。固定预算的 V2/V3/V4 对照与前瞻验证方案见 [`docs/joint_v4_validation_protocol.md`](docs/joint_v4_validation_protocol.md)，当前未提供经验性能比较结果。

## 安装与启动

在项目根目录执行：

```powershell
python -m venv .venv
$env:PIP_CACHE_DIR = "$PWD\.pip-cache"
.\.venv\Scripts\python.exe -m pip install -e . -r .\requirements-dev.txt
.\.venv\Scripts\python.exe -m ecd_epitope_scanner serve --host 127.0.0.1 --port 8765
```

浏览器打开 `http://127.0.0.1:8765`。首页会先显示冻结的 PD-L1 示例；提交真实作业后可下载完整结果和保存自定义组合。

### 双击启动（无需命令行）

整个项目自包含在本文件夹内（源码 + `.venv` 虚拟环境 + FoldX/EvoEF2 引擎 + 前端资源）。
双击本文件夹里的 **`run_ecd_epitope_workbench.bat`**：

- 在最小化窗口后台启动服务，并**轮询等待**端口就绪（最多 25 秒）后再打开浏览器，
  避免服务还没起来就打开页面；
- 已经在运行时不会重复启动，直接复用并打开浏览器；
- 找不到虚拟环境或启动失败时，窗口会停在错误信息上而不是一闪而过。

用完后双击 **`stop_ecd_epitope_workbench.bat`** 结束服务（释放 8765 端口）。

文件夹可整体移动或复制到本机其他位置后继续使用（`.venv` 依赖本机的 Anaconda 基础解释器）；
换电脑则需要重建虚拟环境：`python -m venv .venv` 后执行
`.\.venv\Scripts\python.exe -m pip install -e . -r requirements-dev.txt`。

## CLI 调用

```powershell
# 只输入 UniProt accession；严格 AF3 回退策略
.\.venv\Scripts\python.exe -m ecd_epitope_scanner analyze --uniprot Q9NZQ7 --output .\results\Q9NZQ7

# 使用指定 PDB，自动从 SIFTS 推断 UniProt 映射
.\.venv\Scripts\python.exe -m ecd_epitope_scanner analyze --pdb-id 8AOM --chain A --output .\results\8AOM

# 使用已下载结构。建议提供 accession；否则必须提供完整序列和 ECD 边界。
.\.venv\Scripts\python.exe -m ecd_epitope_scanner analyze --structure .\input\model.cif --accession Q9NZQ7 --chain A --output .\results\uploaded

# 上传的官方 AF3 结果可用 confidence JSON 启用 AF3 confidence/PAE 门槛
.\.venv\Scripts\python.exe -m ecd_epitope_scanner analyze --structure .\input\af3_model.cif --accession Q9NZQ7 --chain A --confidence-json .\input\af3_confidences.json --output .\results\af3

```

如果严格 AF3 策略找不到合格模型，CLI 返回退出码 `2` 和 `WAITING_FOR_AF3_UPLOAD`；网页可直接上传 AF3 文件后恢复该作业。

## 输出内容

每个作业目录包含：

- `structure.cif`：完整来源结构，用于复现和外部能量引擎。
- `ecd_structure.cif`：只含映射 ECD 的网页显示结构。
- `mapping.csv`、`residues.csv`：UniProt/PDB 编号与逐位点审计表（含保守性分数）。
- `regions.json`：Fab 足迹尺度补丁、表位先验、tier、覆盖方案与原因。
- `designs.json`、`designs.csv`：Stage-1 探针、备选探针与 Stage-2 去卷积/双点确认预设。
- `panel_coverage_matrix.csv`：构建体 × 补丁的扰动覆盖矩阵（解读阴性结果用）。
- `sequences.fasta`、`coding_sequences.fasta`、`sequences.gb`：ECD WT/突变蛋白与编码片段。
- `full_length_protein_sequences.fasta` / `full_length_coding_sequences.gb`：全长膜蛋白突变序列（细胞表面结合实验用）。
- `methods_report.md`：阈值、方法、补丁 tier、Stage-1/Stage-2 实验计划、盲区与实验确认门槛。
- `run_manifest.json`：来源 URL、参数、软件版本和 SHA-256。
- `structure_reliability.json`、`residue_reliability.csv`、`structural_segments.csv`、`surface_mask.csv`、`idr_candidates.csv`：结构可靠性、表面掩膜和无序区审计输出。
- `surface_v2_qc.json`、`surface_v2_footprints.json`：表面足迹枚举、去重、区域归属与 QC 输出。
- `fixed_assay_profile.json`：固定双细胞实验方位和解释边界的 assay profile 快照。
- 历史 Joint V3 输出：`manifest.json`、`structure_audit.json`、`topology_context.json`、`surface_residues.csv`、`microfootprints.json`、`probe_sectors.json`、`candidate_probes.csv`、`coverage_by_scenario.csv`、`selected_panel.json`、`exclusions.csv` 和 `scientific_report.md`。
- Joint V4-ML 额外输出：`tool_inventory.json`、`structure_contexts.json`、`prediction_context.pdb`、`residue_mapping.csv`、`surface_evidence.csv`、`sasa_comparison.csv`、`spurs_single_matrix.csv.gz`、`single_mutation_options.csv`、`discotope_residues.csv`、`epitope_islands.json`、`probe_sectors_v4.json`、`panel_candidates.json`、`candidate_combinations.csv.gz`、`combination_stability.csv`、`selected_panel_v4.json`、`v4_manifest.json` 和 `v4_scientific_report.md`。V4.1 另外输出 `surface_graph_v41.json`、`surface_display_regions_v41.json`、`epitope_hypotheses_v41.json`、`hypothesis_membership_v41.csv`、`probe_hypothesis_matrix_v41.csv` 和 `panel_identifiability_v41.json`。
- `ecd_epitope_design.zip`：完整导出包。

## 外部稳定性引擎

以下为历史 `surface_v2`/legacy 与 Joint V3 的 FoldX/EvoEF2 评估说明，不用于当前新作业。

本机开发目录可放置用户自行取得的 FoldX 5.1 Windows 程序；FoldX 二进制及其许可文件不会随 Git 仓库分发。默认先对每区关键位点的候选单点替换做 FoldX 单点筛选，只有通过 `ΔΔG <= +0.46 kcal/mol` 的替换才进入组合搜索。最终区域组合再统一用 `BuildModel --numberOfRuns=5` 直接计算并取中位数。

若要覆盖仓库内程序，可设置环境变量：

```powershell
$env:EVOEF2_PATH = 'C:\tools\EvoEF2\EvoEF2.exe'
$env:FOLDX_PATH = 'C:\tools\foldx\foldx.exe'
```

如果使用项目本机开发包中的 FoldX，请将可执行文件放在 `vendor/foldx/`，或显式设置 `FOLDX_PATH`。从 GitHub 克隆后需要用户自行下载并接受 FoldX 许可；项目不会把该二进制上传到公共仓库。

当 FoldX 不可用时，结果会明确降级为无量纲代理并标记为不可推荐；不会把代理分数写成 kcal/mol。FoldX 官方说明建议先 `RepairPDB`，再运行 `PositionScan`/`BuildModel`；程序输出能量项单位为 kcal/mol。

`vendor/foldx/foldx_1_20270131.exe`（FoldX 5.1）已实测可用：该构建**无需 rotabase.txt**（rotamer 库内嵌），PD-L1 ECD 的 `RepairPDB` 约 4 分钟（按结构哈希**跨作业缓存**于 `workspace_data/foldx_repair_cache/`，同一靶点重复运行直接复用）。执行性能设计：

- BuildModel 批次按近等分块**并行**跑（`config.FOLDX_MAX_WORKERS`，默认 8 进程），结果按输入顺序聚合；
- 单点筛选 runs=1（门限判断足够）；Stage-1 探针/备选组合复评 runs=5 取中位数；**Stage-2 预设只做 runs=1 筛查**（订购前可要求补 5-run 中位数，警告会写明）；
- 每位点候选 top-2 进入 FoldX（`MAX_CANDIDATES_PER_SITE`），关键位点数 4（`KEY_SITES_PER_PATCH`）；
- 批次超时按 60 s/构建-run 自动缩放（上限 4 小时），引擎挂起只降级不炸作业。

150 aa ECD 的完整面板（含真实 FoldX）实测约 10 分钟量级；同一靶点第二次运行因修复缓存再省数分钟。

## 关键实现模块

- `retrieval.py`：UniProt、3D-Beacons、PDBe/SIFTS 结构来源与候选审计。
- `structure_analysis.py`：mmCIF/PDB 解析、三套残基编号、SASA、界面、糖链、二硫键、二级结构和置信度。
- `conservation.py`：UniProt reviewed 直系同源比对，逐残基保守性打分（离线优雅降级）。
- `partitioning.py`：Fab 足迹尺度表面补丁枚举 + 表位先验 + 覆盖方案。
- `design.py`、`stability.py`：定向替换字母表、锚点突变、组合全枚举、Stage-2 去卷积设计、EvoEF2/FoldX 适配器。
- `joint_design/`：Joint V3 结构审计、拓扑上下文、表面单元、多尺度窗口、ProbeSector、候选 beam search、覆盖矩阵、MILP/greedy 面板优化和科学报告。
- `joint_design_v4/`、`integrations/`：Joint V4-ML 的规范化预测上下文、证据融合、表面分区、候选级稳定性门控、确定性面板优化、JSON worker 运行器和工具输出校验。
- `esm_scores.py`：可选 ESM zero-shot PLM 证据线（opt-in，见模块说明）。
- `config.py`：全部科学阈值集中管理（含文献依据，输出到 manifest/methods report）。
- `api.py`、`jobs.py`：FastAPI、SQLite、本机单队列、SSE、AF3 恢复和自定义构建。
- `web/`：固定版本 PDBe Mol* 结构显示、步骤筛选、位点审查、分区、AA/CDS 对照和导出。

方法学与实验设计依据见 `docs/改进方案.md`。

## 验证

```powershell
.\.venv\Scripts\python.exe -m pytest .\tests -q --basetemp .\test_output\pytest_tmp -p no:cacheprovider
```

当前完整回归测试为 106 项通过；V4 专项覆盖 worker 超时/崩溃/输出污染、预览与 strict 状态隔离、模型分歧阻断、候选级正式推荐门控、DiscoTope 优先/稳定性回退系列、缺失 DiscoTope 时的审计降级、规范化残基映射、GPU 请求/缓存隔离、CUDA OOM 识别与重试/CPU 回退计划、FoldX 短名单上限与 sector 均衡、表面中继路径收缩、k-center 分区、孤立表位先验及完整表面分区。测试只验证软件契约和确定性规则，不替代外部模型自身的 benchmark 或实验验证。

测试包含 SIFTS/拓扑解析、AF3 严格识别、AF3 恢复回归（UniProt 作业回填 accession）、序列反向翻译、冻结 PD-L1 全流水线（补丁分区/探针锚点/Stage-2 预设/实验计划）、Joint V3 确定性面板与保护位点测试、Joint V4 worker 契约/证据门控/规范化映射/模型分歧阻断/每补丁候选系列/DiscoTope 审计降级/GPU 缓存隔离/CUDA OOM 恢复计划/FoldX 短名单/表面中继路径/k-center 分区/孤立表位先验测试、补丁分区与定向菜单单测、区域构建约束与 FastAPI 静态资源检查。当前工作区回归测试为 106 项通过；前端 JavaScript 语法检查和 Ruff 检查通过。
