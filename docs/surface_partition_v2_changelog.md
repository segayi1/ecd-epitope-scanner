# Surface partition V2 历史修改日志

> 当前新作业固定使用 Joint V4-ML。本文记录历史分区开发；旧引擎回退选项已从新作业入口移除。

## 2026-09-16：有限 SurfaceRegion 聚合与唯一显示颜色

- 分区输入改为每中心一个 primary footprint，优先合格 R12、回退 adaptive；敏感性窗口与现有去重结果完整保留。
- 使用固定 seed=0 的 NetworkX Louvain，加入真实共享 surface vertices、测地邻近和 domain 软惩罚；拆分不连通/跨度过大的社区，仅沿共享边合并小社区。
- 默认测地跨度上限 60 Å，定义为 primary centers 在原 surface graph 上的最大距离，不是全区域 vertex diameter；该初始阈值尚待跨目标 benchmark 校准。
- Region 面积按 vertex union 精确去重。显示归属按最近 primary center 的表面面积投票确定；footprint 科学 membership 仍允许重叠。
- 所有自然区域均保留，有限实验预算只用于旧设计器输入选择；不再用暴露度排序截取显示区域。
- 代表 footprint 选择优先可突变性、残基数和社区 medoid；adapter 向设计器传递可突变子集，完整成员保留在 Footprint 导出。
- 最终区域顺序使用唯一调色板；点击区域淡化其他区域，再次点击恢复；Mol* 基础色与 beads/离线报告共用 residue.region_id。
- 使用 e81db31a230d 已保存数据重放：216 → 6 个区域；保留 567 个去重 footprints；183 个表面残基具有唯一显示归属；三个重复运行结果完全一致。
- 重放入口：`python -m ecd_epitope_scanner.surface_patch.replay_regions RUN_DIR OUTPUT_DIR`。重放只生成分区报告，不能将旧设计结果冒充为新区域设计。
- 本轮未修改 structure_analysis.py、design.py、stability.py、序列处理或生产 SASA/footprint 生成计算；未引入 mesh backend 或新膜评分。

## Surface-geodesic V2（当前实现）

- 新增独立 `surface_patch` 包；底层对象为 `SurfaceGeometry`、`Footprint` 与 `SurfaceRegion`。
- 默认 `surface_v2` 使用 Shrake–Rupley 原子球未遮挡点构造显式表面点云；该几何采样不改写生产 SASA 数值。
- surface graph 边由局部点距离、法向一致性与排除体积线段检查共同确定，不再用 Cα 距离建立主要表面连通性。
- 使用 SciPy CSR 图和 Dijkstra 计算真实表面路径，并按连通组件执行完全确定性的 geodesic farthest-point sampling。
- 每个中心生成 9/12/15 Å 固定半径及 600 Å² 自适应面积 footprint；所有 footprint 允许重叠。
- residue membership 使用 footprint 内面积（≥8 Å²）或面积比例（≥20%）；保护残基保留在几何和 footprint 中，但不进入可突变列表。
- 已实现几何 QC、area-weighted Jaccard 去重、footprint adjacency、抗链式合并的 SurfaceRegion grouping，以及旧 `RegionRecord` adapter。
- `3–7` 仅用于 Stage-1 兼容面板选择；完整 footprints/SurfaceRegions 均保留并导出。
- 新增 surface vertices/edges、raw/deduplicated footprints、membership、adjacency、regions 和 QC 输出；在线页面可选择并聚焦任一 footprint。
- 历史版本曾通过 `ECD_PATCH_ENGINE=legacy_v1` 或请求 metadata 的 `patch_engine` 回退旧引擎；当前版本不再接受这些入口。
- MolecularMesh backend 已保留接口；当前默认实现为无需外部表面二进制的 `sr_pointcloud`。

## 回退基线

- 分支：`codex/surface-partition-v2`
- 提交：`6a53412 checkpoint: preserve pre-surface-partition-v2 workbench`
- 回退方式：`git switch codex/surface-partition-v2` 后使用 `git revert <阶段提交>`；若仅查看旧版本，可使用 `git show 6a53412:<文件>`。

## V2-P1：局部溶剂表面法向（进行中）

- 停止使用蛋白整体质心到残基的径向量作为生产分面依据。
- 使用与生产 SASA 相同的 1.4 Å 探针和固定原子半径，在原子球面上采样未遮挡点并汇总局部外向法向。
- `ResidueRecord` 新增 `surface_normal`；旧结果缺少该字段时不得启用 V2 自动分面。
- 新增紧凑模式独立配置，不修改原有通用 `PATCH_*` 阈值。

## V2-P2：β-sheet 骨架分区与多足迹设计（进行中）

- 对带有可信 Beta strand 注释的紧凑 ECD，根据实际 Cα 接触识别 β-sheet 骨架；紧邻序列的转角片段不得作为跨 sheet 连接。
- Q15116 重放结果从 3 个质心径向扇区改为 2 个 β-sheet zone。
- 每个 zone 最多记录 3 个非冗余局部 footprint。
- 紧凑模式将候选搜索从单个 footprint 的 5 个位置扩大到整个 zone 的最多 12 个优先位置，仍统一批量执行 FoldX。
- 多点组合失败但存在 FoldX 通过的强扰动单点时，输出 `RECOMMENDED_SINGLE_SITE_PROBE`，不再直接将整个区域声明为盲区。
- 多点组合必须完整落在同一个候选 footprint 内；仅仅满足欧氏距离上限不足以跨 footprint 组合。
- 实时工作台、离线 HTML 和 methods report 显示分区依据与候选 footprint 数量。
- 单点回退若来自扩展的全区候选搜索，会将其局部单点 footprint 写回区域数据，保证“推荐突变属于某个候选 footprint”的审计约束成立。

## Q15116 验证记录

### 中间验证（仅保留作缺陷审计）

- 目录：`results/q15116_surface_v2_validation`
- 得到 2 个 β-sheet zone，两个区域均产生构建。
- 发现区域 1 多点组合仅满足欧氏距离约束、未完整落入同一个候选 footprint；该结果不作为验收版本。

### 修正后验证

- 目录：`results/q15116_surface_v2_validation_r2`
- 结构输入：原任务 `54ce64e6227e` 保存的官方 AF3 CIF 与 confidence JSON。
- 分区：2 个区域，依据为 `beta_sheet_backbone_topology`；表面包络共 82 个可信近表面残基。
- 区域 1：FoldX 通过的三点构建 `L138D/T36K/S38K`，完整属于候选 footprint `[36,38,39,138]`，5 次 BuildModel 中位数 -1.64 kcal/mol。
- 区域 2：多点 footprint 无合格组合，输出 FoldX 通过的单点回退 `A40R`，5 次 BuildModel 中位数 -0.57 kcal/mol。
- 两个区域均有可实验构建；FoldX 批处理 8 个最终构建，全部返回可解析结果。

## 当前回退点

- `6a53412`：Surface partition V2 开始前的完整工作台。
- `aecfc34`：溶剂表面法向、β-sheet 分区、多 footprint 与初版单点回退。
- `fac2fac`：强制多点组合属于单一 footprint，补齐网页/离线报告与 Q15116 验证日志。

## 文档同步

- `docs/method.md` 已加入紧凑单域模式的触发条件、原子表面法向、β-sheet 接触拓扑、多 footprint 和单点回退方法。
- 同时修正表面阈值文字为 rSASA 与侧链 SASA 同时通过，将每位点候选数更新为当前配置的 2，并修正 ΔΔG 门槛方向。
