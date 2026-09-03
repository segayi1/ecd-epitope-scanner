# ECD Epitope Workbench

本机单用户的单次跨膜蛋白胞外区（ECD）结构分区与组合突变设计工具。它生成可审计的抗原突变构建，不依据抗体序列预测表位，也不把计算筛选误写成实验表位结论。

## 已实现的科学边界

- 输入 UniProt accession、PDB ID，或从相关网站下载的 PDB/mmCIF 文件。
- UniProt 输入自动下载 canonical 序列、拓扑、结构域、二硫键、糖基化和功能注释，并通过 3D-Beacons/PDBe 搜索实验结构。
- PDB 候选需要满足：ECD 观察覆盖率 >=95%、最大内部缺口 <=8 aa、与 canonical 序列一致性 >=95%、主要结构域覆盖 >=90%。逐残基编号使用 UniProt、PDB author 和 mmCIF label 映射保存。
- 没有合格实验结构时，只接受带明确 AF3 来源的公开模型；否则状态为 `WAITING_FOR_AF3_UPLOAD`，网页要求上传官方 AF3 mmCIF 与 confidences JSON。AFDB 不会被伪装为 AF3。
- 用 FreeSASA（若可用）或同参数的 Bio.PDB Shrake-Rupley 回退计算 rSASA/侧链 SASA；默认候选标准为 `rSASA >=0.25` 且 `side-chain SASA >=20 A2`。
- 自动排除二硫键 Cys、注释糖基化/PTM、低置信度 AF3 位点、坐标缺失和埋藏位点；功能界面保留在候选中但显示风险标签。
- 以结构域、膜近端/远端、二级结构边界和表面测地距离分成目标 5-8 区；结构不支持时如实输出更少区域与原因。
- 每区内部评估最多 8 个关键位点的 19 种替换，最终只输出一条 2-5 点组合突变体；无满足约束的组合则输出 `NO_LOW_RISK_COMBINATION`。
- 输出仅含 ECD 的氨基酸序列和编码片段，不自动添加信号肽、标签、起始/终止密码子、TM 或胞内区。无原始 CDS 时使用 DNA Chisel 的人源/HEK293 密码子优化。

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
- `mapping.csv`、`residues.csv`：UniProt/PDB 编号与逐位点审计表。
- `regions.json`：结构语义分区、表面积、测地直径、界面比例与原因。
- `designs.json`、`designs.csv`：关键位点饱和扫描及每区推荐组合。
- `sequences.fasta`、`coding_sequences.fasta`、`sequences.gb`：ECD WT/突变蛋白与编码片段。
- `methods_report.md`：阈值、方法、实验确认门槛。
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

## 关键实现模块

- `retrieval.py`：UniProt、3D-Beacons、PDBe/SIFTS 结构来源与候选审计。
- `structure_analysis.py`：mmCIF/PDB 解析、三套残基编号、SASA、界面、糖链、二硫键、二级结构和置信度。
- `partitioning.py`：结构域/膜方向约束的表面测地分区。
- `design.py`、`stability.py`：19 种替换扫描、组合突变选择、EvoEF2/FoldX 适配器。
- `api.py`、`jobs.py`：FastAPI、SQLite、本机单队列、SSE、AF3 恢复和自定义构建。
- `web/`：固定版本 PDBe Mol* 结构显示、步骤筛选、位点审查、分区、AA/CDS 对照和导出。

## 验证

```powershell
.\ecd_epitope_scanner\.venv\Scripts\python.exe -m pytest .\ecd_epitope_scanner\tests -q --basetemp .\ecd_epitope_scanner\test_output\pytest_tmp -p no:cacheprovider
```

测试包含 SIFTS/拓扑解析、AF3 严格识别、序列反向翻译、冻结 PD-L1 全流水线、区域构建约束与 FastAPI 静态资源检查。
