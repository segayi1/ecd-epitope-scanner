# 结构可靠性分流（规则 v1）

## 采用范围与判断

此层只作用于 Surface V2 的几何输入与分区候选，不改原始结构文件、SASA 计算、残基筛选阈值、突变设计器、FoldX、Stage 2 或序列生成器。几何输入变化会合理地改变分区及下游候选数量，不意味着旧设计结果逐项保持一致。

低 pLDDT 不是无序的充分条件；高跨域 PAE 不是局部展开的充分条件；实验 B 因子不作为 pLDDT。参考 [AlphaFold 官方说明](https://alphafold.ebi.ac.uk/faq)。本实现是保守的工程启发式，不是经实验校准的无序预测器。

- STRUCTURED：有坐标的实验结构，或预测 pLDDT ≥70。
- FLEXIBLE_STRUCTURED：两端有同域高置信结构支撑的短低置信片段（≤14 aa），保留几何。
- UNCERTAIN：证据不足，保留几何；参与足迹时输出质量标记。
- DISORDERED：连续低置信候选段 ≥15 aa，至少 80% pLDDT <50、低非局部接触，并有序列无序或高核心 PAE 支持，且没有已观测 H/E 冲突。结构连接处保留 3 aa 过渡，因此最终屏蔽核心可能短于 15 aa。
- CONDITIONALLY_FOLDED：高置信螺旋同时有序列无序、高核心 PAE 和低接触证据。只是条件折叠假设，默认分流，不声称存在结合诱导折叠。
- UNRESOLVED：没有可用坐标，不等同于无序，也不强行生成三维表面。

PAE 优先校验显式链/残基 token；无显式映射时仅接受维度匹配的单链矩阵。默认优先同域核心，降低活动结构域被误判的风险。缺失证据不补为阳性。7 aa 中位数证据分数只用于诊断，分类采用片段共识而非任意加权分数阈值。

## 几何与边界

默认 `strict`：屏蔽残基既不生成表面点，也不作为遮挡原子。`occlusion_sensitivity`：仍保留其遮挡，但不生成表面点，供敏感性比较。

屏蔽/缺失位置邻近 2 aa 保留遮挡原子，不生成表面点，防止切断边界成为伪中心。这是保守的整残基边界保护，不是网格切口重建。`use_for_surface_mesh` 表示结构分支归属，`use_for_3d_footprint` 进一步扣除边界；以 `surface_mask.csv` 审计。

当前后端仍是 `sr_pointcloud`。原 SASA 保持原始结构计算值，不能解读为屏蔽后表面积；IDR 候选的 SASA 明确为空。

## 测试与使用

重启本地服务，新建任务。分区页展开“结构可靠性 / Reliability”，可切换独立着色及查看片段、IDR 窗口。类别颜色使用配色技能推荐的 Okabe–Ito 色系并附文字标签，不与原分区着色叠加。没有坐标的 UNRESOLVED 仅在列表出现。

命令行 `analyze` 新增：

- `--iupred-output <完整序列预测文件>`：导入 IUPred 的 position/AA/score 文本，严格校验全长 canonical 序列、编号和分数，保存输入哈希及头信息；不自动安装或运行 IUPred。
- `--reliability-surface-mode strict|occlusion_sensitivity`。
- `--disable-reliability`：几何分流旁路，用于与旧结果比较。

Python `PipelineRequest.metadata` 对应 `iupred_output` 和 `structure_reliability` 配置字典。网页当前使用默认规则，尚未提供 IUPred 上传控件。

导出新增 `residue_reliability.csv`、`structural_segments.csv`、`surface_mask.csv`、`idr_candidates.csv`、`reliability_qc.json`，纳入 ZIP/manifest。IDR 默认 12 aa 窗口、步长 5、补齐末端；仅待审核，不自动设计突变或调用 FoldX。C/G/P、N-X-S/T 与已有 PTM/二硫键保护保留。

验证命令（仓库目录 PowerShell）：

```powershell
$env:PYTHONPATH='E:\user\New project'
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
```

合成双构象验证只证明严格屏蔽后的表面采样不受被排除尾部坐标影响，不能替代多种真实 AF 种子/蛋白的生物学基准测试。条件折叠分流也需要人工审核；尚未实现 IDR 专用免疫原性或聚集风险预测。
