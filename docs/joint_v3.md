# Joint V3 后端设计引擎

Joint V3 是一个可选的后端设计引擎，现有 Surface V2 仍是默认引擎。启用方式：

```powershell
ecd-epitope-scanner analyze --structure <structure.pdb> --sequence <sequence.fasta> `
  --ecd 1:220 --design-engine joint_v3 --output results/joint_v3
```

V3 从现有残基级结构分析结果开始，依次生成结构审计、拓扑不确定性、残基级表面单元、多尺度 micro-footprint、连通 `ProbeSector`、突变候选和有限预算面板。生产流水线复用显式溶剂表面点图及其穿模检查；仅在直接调用 V3 且未提供表面模型时使用带法向量和空间净空检查的残基级降级图，并在 manifest 中标记。`ProbeSector` 合并现在同时要求同一残基图连通分量、测地直径不超过 35 Å、Cα 三维跨度不超过 30 Å；显式表面点边折叠到残基图时还拒绝超过 10 Å 的 Cα 长边，避免把空间上分离的片段合成一个分区。候选组合使用确定性、规范化去重的 beam search 并分别保留 2、3、4 位点候选。面板选择使用两阶段 MILP：先最大化多情景最差覆盖，再在保持最差覆盖的条件下优化平均覆盖和风险；MILP 不可用时使用确定性 greedy fallback。检测到 FoldX 时，V3 先用规则风险生成代理预选面板，优先保证该面板及替代面板进入最多 81 个候选的组合级 `RepairPDB + BuildModel` 批次，其余名额按 sector 轮转分配；最终面板仅从成功取得组合 ΔΔG 的候选中选择，并将真实 kcal/mol 值写入最终设计。外部引擎不可用时明确降级为 `RULE_ONLY`/`PROXY_ONLY`，不生成伪 kcal/mol 或结合概率。

网页可在 Surface V2 与 Joint V3 间显式选择，并继续使用选中面板 sector 的兼容投影 `regions.json` 和 `designs.json`；所有未入选 sector 仍保留在 V3 原生输出。V3 的机器可读结果另写入 `structure_audit.json`、`topology_context.json`、`surface_residues.csv`、`microfootprints.json`、`probe_sectors.json`、`candidate_probes.csv`、`coverage_by_scenario.csv`、`exclusions.csv`、`selected_panel.json` 和 `scientific_report.md`。

重要解释边界：ProbeSector 是粗定位区域，不是真实抗体 footprint；孤立 ECD 不用于推断唯一膜平面；Myc/表达信号不能单独证明 ECD 正确折叠；所有覆盖率都是结构几何情景下的代理指标，不能解释为结合概率。
