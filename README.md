# ECD Epitope Workbench

本机单用户的单次跨膜蛋白胞外区（ECD）结构分区与组合突变设计工具。它生成可审计的抗原突变构建，不依据抗体序列预测表位，也不把计算筛选误写成实验表位结论。服务仅允许本机回环地址；公开部署需要未来的专用运行配置、认证与作业隔离，不能通过修改 `--host` 直接开启。

## 已实现的科学边界

- 输入 UniProt accession、PDB ID，或从相关网站下载的 PDB/mmCIF 文件。
- UniProt 输入自动下载 canonical 序列、拓扑、结构域、二硫键、糖基化和功能注释，并通过 3D-Beacons/PDBe 搜索实验结构。
- PDB 候选需要满足：ECD 观察覆盖率 >=95%、最大内部缺口 <=8 aa、与 canonical 序列一致性 >=95%、主要结构域覆盖 >=90%。逐残基编号使用 UniProt、PDB author 和 mmCIF label 映射保存。
- 没有合格实验结构时，只接受带明确 AF3 来源的公开模型；否则状态为 `WAITING_FOR_AF3_UPLOAD`，网页要求上传官方 AF3 mmCIF 与 confidences JSON。AFDB 不会被伪装为 AF3。
- 统一用 Bio.PDB Shrake-Rupley 与固定 Tien/DSSP 原子半径计算 rSASA/侧链 SASA；默认候选标准为 `rSASA >=0.25` 且 `side-chain SASA >=20 A2`，不因可选依赖而改变生产结果。
- 自动排除二硫键 Cys、注释糖基化/PTM、低置信度 AF3 位点、坐标缺失和埋藏位点；功能界面保留在候选中但显示风险标签。
- 以 Fab 足迹尺度（测地半径 13 Å，约 15-25 残基，Kringelum 2012）在表面测地图上枚举补丁并允许跨结构域；每个补丁带表位先验（暴露度/外凸度/环区占比/进化可变性/糖屏蔽/膜近端惩罚）并分级 HIGH/MEDIUM/LOW；面板预算随 ECD 大小自适应（~1 补丁/22 残基，上限 7），结构不支持时如实输出更少补丁与原因。
- 实验形式 `assay_format=cell_surface|purified` 是一等参数：细胞表面实验时膜近端补丁先验降权并标记 IgG 不可及；`purified` 时全部保留。
- 每个补丁按 WT 残基类型定向选择替换字母表（大侧链→Ala、小侧链→Trp 插入、带电→电荷反转），Stage-1 探针必须含至少一个锚点突变（电荷反转或 W/Y/F/R/K→A）；在 FoldX/proxy 清空的候选菜单上**全枚举** 3-6 点组合（两两间距 ≥5.5 Å、最大展开 ≤22 Å），按**帕累托前沿**（破坏势 vs 稳定性风险）选点：主探针取膝点（单位风险破坏力最大），备选探针取激进端，完整前沿写入 `epitope_mapping_plan.pareto_fronts` 审计；无满足约束的组合则输出 `NO_LOW_RISK_COMBINATION`。
- 每个探针自动生成 Stage-2 预设：单点去卷积构建体（`STAGE2_SINGLE`）与确认性双点组合（`STAGE2_PAIRWISE`），并输出面板覆盖矩阵、决策规则与显式盲区清单（`epitope_mapping_plan`）。
- FoldX 单点门限与暴露度联动（rSASA≥0.35 → +0.80；0.25-0.35 → +0.60 kcal/mol；文献 RMSE 0.85-1.1），组合门限维持 +2.0 kcal/mol。
- 可选证据线：UniProt 同源保守性（自动，离线优雅降级）与 ESM zero-shot LLR（设置 `ECD_ESM_MODEL` 且本机装有 torch/fair-esm 时启用，LLR<-3 共识否决）。
- 输出仅含 ECD 的氨基酸序列和编码片段，不自动添加信号肽、标签、起始/终止密码子、TM 或胞内区；`original_cds` 可解析为全长 CDS 时额外导出**全长**膜蛋白突变编码序列（细胞流式实验需要），无 CDS 时导出全长蛋白序列供基因合成反向翻译。无原始 CDS 时使用 DNA Chisel 的人源/HEK293 密码子优化。

计算预筛只能说明“低预测风险”。进行表位归因时仍需同时验证表达/膜表面表达、构象对照抗体或折叠控制、稳定性，以及必要时的配体或功能读出。

## 安装与启动

在项目根目录执行：

```powershell
python -m venv ecd_epitope_scanner\.venv --system-site-packages
$env:PIP_CACHE_DIR = "$PWD\ecd_epitope_scanner\.pip-cache"
.\ecd_epitope_scanner\.venv\Scripts\python.exe -m pip install -r .\ecd_epitope_scanner\requirements-dev.txt
.\ecd_epitope_scanner\.venv\Scripts\python.exe -m ecd_epitope_scanner serve --host 127.0.0.1 --port 8765
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
换电脑则需要重建虚拟环境：`python -m venv .venv --system-site-packages` 后安装 requirements-dev.txt。

## CLI 调用

```powershell
# 只输入 UniProt accession；严格 AF3 回退策略
.\ecd_epitope_scanner\.venv\Scripts\python.exe -m ecd_epitope_scanner analyze --uniprot Q9NZQ7 --output .\results\Q9NZQ7

# 使用指定 PDB，自动从 SIFTS 推断 UniProt 映射
.\ecd_epitope_scanner\.venv\Scripts\python.exe -m ecd_epitope_scanner analyze --pdb-id 8AOM --chain A --output .\results\8AOM

# 使用已下载结构。建议提供 accession；否则必须提供完整序列和 ECD 边界。
.\ecd_epitope_scanner\.venv\Scripts\python.exe -m ecd_epitope_scanner analyze --structure .\input\model.cif --accession Q9NZQ7 --chain A --output .\results\uploaded

# 上传的官方 AF3 结果可用 confidence JSON 启用 AF3 confidence/PAE 门槛
.\ecd_epitope_scanner\.venv\Scripts\python.exe -m ecd_epitope_scanner analyze --structure .\input\af3_model.cif --accession Q9NZQ7 --chain A --confidence-json .\input\af3_confidences.json --output .\results\af3
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
- `ecd_epitope_design.zip`：完整导出包。

## 外部稳定性引擎

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
- `esm_scores.py`：可选 ESM zero-shot PLM 证据线（opt-in，见模块说明）。
- `config.py`：全部科学阈值集中管理（含文献依据，输出到 manifest/methods report）。
- `api.py`、`jobs.py`：FastAPI、SQLite、本机单队列、SSE、AF3 恢复和自定义构建。
- `web/`：固定版本 PDBe Mol* 结构显示、步骤筛选、位点审查、分区、AA/CDS 对照和导出。

方法学与实验设计依据见 `docs/改进方案.md`。

## 验证

```powershell
.\ecd_epitope_scanner\.venv\Scripts\python.exe -m pytest .\ecd_epitope_scanner\tests -q --basetemp .\ecd_epitope_scanner\test_output\pytest_tmp -p no:cacheprovider
```

测试包含 SIFTS/拓扑解析、AF3 严格识别、AF3 恢复回归（UniProt 作业回填 accession）、序列反向翻译、冻结 PD-L1 全流水线（补丁分区/探针锚点/Stage-2 预设/实验计划）、补丁分区与定向菜单单测、区域构建约束与 FastAPI 静态资源检查。
