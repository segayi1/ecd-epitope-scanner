# CTLA4 分区与算法修复交接

> 本文件是开工前的交接记录。2026-11-12 完成状态与最终验证证据见 [`handoff_ctla4_completion_20261112.md`](handoff_ctla4_completion_20261112.md)。

更新时间：2026-11-12。工作目录：`E:\user\New project\ecd_epitope_scanner`（本机路径）。当前分支 `codex/joint-v4-ml`，HEAD `8ae5891`。交接用途：在同一工作目录的新对话继续用户已要求的修复和验收；本文件本身不新增授权。

## 用户目标与约束

用户要求检查并修复 CTLA4 诊断和新算法验收中发现的问题，尤其是为何配置最多 12 个分区却只得到 8 个，以及四个区域的区内最短路径直径超过 30 Å。最新明确方向是：**区域内路径的直径应满足配置限制**，同时保留连通性、唯一覆盖、最小残基数与面积等硬约束；可行时通过较小的区域得到更多分区。区域数量仍是软目标，构建体预算与分区数量分别管理；找不到可行解时要明确失败或搜索预算耗尽，不放宽几何门或伪造证据。

更广的算法目标：结构稳定筛选前提下最大化表面理化扰动及可设计表面覆盖；DiscoTope 用作表位优先级而非硬排除，缺失时中性，并保留非表位探索。每区尽量尝试三点，比较其三个精确氨基酸双点父本，允许覆盖或局部扰动的正边际；三点本身仍须组合级 SPURS/FoldX。保守性用于候选来源、风险复核及溯源，不作为固定惩罚或主目标。SNIPR 数据只用于总结规律，不进入生产评分。

遵守用户提供的 AGENTS.md：先读 `docs/project_memory.md` 与相关审计；保留用户已有改动和历史运行；AF3 必须按 AF3 及其置信度来源回放；几何验证与经过科学门的突变推荐分开；从仓库**父目录**运行测试；修复复发问题后更新项目记忆。当前工作树有大量既有及本轮未提交更改和未跟踪文件，不要 reset、整树暂存或覆盖历史任务。提交、部署、删除或新建另一对话均未由本次交接请求授权。

## 实际状态与资料

- 真实诊断任务 `1d812aa2f9be` 已 COMPLETE，目录 `workspace_data/1d812aa2f9be/`。输入是 UniProt P16410、实验 NMR 结构 1AH1 chain A，而非 AF3；125 个映射 ECD 残基、76 个表面残基，产生 8 个分区和 12 个三点主构建体。保存结果的 12 组合均有 SPURS multi/FoldX 且记录 `CONSENSUS_PASS`，仅是计算筛选，不是实验验证。该任务是修复前的开发行为，不能当作新协议的验收结果。
- `docs/ctla4_diagnostic_review_20261112.md` 保存诊断详情，但其中 FoldX runtime 一项写作未修复，已被 `docs/project_memory.md` 尾部 “CTLA4 partition and runtime follow-up” 更新：runtime 已做源码修复和测试，历史结果不变。历史条目按时间记录，不应将旧状态当作当前状态。
- 几何证据见任务内 `partition_diagnostics_v4.json`、`surface_graph_v4.json`、`surface_display_regions_v41.json`、`analysis.json` 和结构来源文件；查看 `run_manifest.json`、`execution_trace.jsonl` 以核对来源和流程。输出体积很大，避免无必要将完整 `analysis.json` 等打印到终端。旧的 `docs/ctla4_small_region_version_audit_20261008.md` 记录另一份 CTLA4 输入的较小分区见证，**不等于**对本次任务的严格区内直径回放。
- 现行 `joint_design_v4/sectors.py` 的 `balanced_connected_v4` 对区域大小偏好排序后倾向更少区域，12 是软目标；直径检查使用全 eligible surface graph 最短路径，路径可穿过其他区域。任务的 VS002、VS004、VS005、VS007 保存的 induced-region 直径分别约 37.177、36.600、33.549、39.349 Å。`region_metrics` 已同时输出 full-graph 与 induced-graph 直径，但当前硬门和 `V4Sector.geodesic_diameter_a` 按前者计算。必须核实完整调用链和图中 contracted surface transit 的含义，再改约束与报告；不要把不同图定义混在一起。
- `tests/test_balanced_display_partitions.py` 当前 `assert_partition` 按全图计算；测试还包括保存的 PD-L1 和历史 `connected_bounded_v3` 行为。更改时区分当前新要求与旧协议复现，增加真实保存 CTLA4 输入回放和合成反例，并检查兼容性；新的严格划分若实际不可行，报告原因与搜索范围，不预先断言可行或固定必须 12 区。

## 已做修复与验证

- `conservation.py` 对没有相同残基的比对列使用 0 而非触发 `KeyError`；`tests/test_mapping_design.py -k conservation` 5 passed。
- `joint_design_v4/identifiability.py` 对新选择器将面板边际扰动代理排在三点父本差值之前，保留新分区优先及几何正增益路径；`pipeline.py` 的 V4 表位先验标记根据选择器及有评分证据推导。相关定向测试 97 passed、2 skipped（有 Starlette 弃用警告）。
- 诊断后 `joint_design_v4/engine.py` 修复 adaptive/site-first 汇总 `PredictionResult` 遗失 FoldX worker runtime 的问题，记录 repair 与 BuildModel 合计及 worker metrics；`tests/test_ordinary_adaptive_engine.py` 7 passed。历史 CTLA4 报告仍显示 0，不重写历史。修改后尚未重新进行真实模型运行。
- `joint_design_v4/site_first.py` / `regional_perturbation.py` 已加强三点的三个精确父本完整性检查，相关测试属于上述定向回归；仍需沿真实候选生产路径、完整模型门和导出审查，不能据单项测试宣称整体验收。

## 待修复与接续顺序

1. **优先处理几何。** 从 `docs/project_memory.md`、`docs/ctla4_diagnostic_review_20261112.md` 和相应修复审计开始，核对当前代码与工作树；用任务 `1d812aa2f9be` 的实际表面图、结构、残基映射构造只读回放，先复现 8 区与四个超限。实现带缓存的区内最短路径硬检查，并相应处理初始划分、拆分、边界移动、最终验收与输出口径；保留合理的有限搜索和无法找到解时的诊断。回放必须独立写入新的测试临时目录，不能覆盖原任务。验收覆盖每点恰好一次、连通、最小尺寸/面积、区内测地直径 ≤30 Å、Cα 直径 ≤22 Å、重复运行确定性以及搜索失败行为。按实际输出报告分区数。
2. **新旧协议隔离。** 当前 `ordinary_site_first_v2` 已直接混入表位扰动 selector、三点增强门和预算变化，历史协议语义并未隔离。`site_first.py`、`engine.py`、`config.py`、`pipeline.py`、API、导出、前端及 FoldX cohort/adapter/worker 均需检查。建议用明确新协议名，仅新协议启用新行为，恢复/保留旧协议与旧文件语义，并保留 site-first FoldX 批次身份及完整稳定性门；新协议名称和切换策略是实现建议，须用测试证明后再陈述已隔离。先 shadow replay，未通过前不宣称新协议正式验收。
3. **输出体积。** `joint_design_v4/sector_recommendations.py` 每张卡将完整 `source_provenance` 重复内嵌，ledger 和模型元数据占主要体积；本次 `selected_panel_v4.json` 约 196 MB、`analysis.json` 约 190 MB。`web/recommendations.js` 消费卡片级模型和结构信息；`v4_manifest.json` 含运行级 provenance。去重时先确认 UI、API、离线导出和引用能找到原始证据，量化体积并测试一致性；不能简单丢弃科学来源或修改历史文件。
4. **科学性及面板验收。** `joint_design_v4/spatial_support.py` 当前是候选整体化学强度在 12 Å 邻域上均匀传播，并把突变位点表位先验传播到邻域；它是粗糙代理，不是已验证的逐点物理场或结合损失。逐点 DiscoTope 证据和局部响应若要改进，先明确可复现的模型、缺失行为和证据，不能通过任意权重宣称已校准；持续保留非表位探索、完整三点模型证据和候选不足时的缺口。检查搜索截断对父本/三点完整性的影响，并做相同输入、预算和硬门条件下的新旧 shadow replay、导出一致性及真实模型小范围复测；实验结合机制仍需湿实验。

## 验证入口

从 `E:\user\New project` 执行，例如：

```powershell
& '.\ecd_epitope_scanner\.venv\Scripts\python.exe' -m pytest ecd_epitope_scanner/tests/test_balanced_display_partitions.py ecd_epitope_scanner/tests/test_ordinary_adaptive_engine.py -q
```

按具体修改扩展到 site-first、FoldX、panel、API、reporting、UI 和保存输入回归；旧报告中所列通过数不替代修改后的测试。当前 checkout 混合 CRLF 与已有空白问题，`git diff --check` 会产生大量既有噪声；不要以此声称本轮 diff 清洁或擅自批量格式化。服务或任务状态在新对话中重新查询，避免误用旧进程或重复启动长时间任务。

建议接手第一步：核对工作目录、`git status` 和上述任务文件，只读复现实际 CTLA4 分区，再针对区内直径实现与回归；此前记录的是事实与边界，不代表尚未完成的工作已验收。
