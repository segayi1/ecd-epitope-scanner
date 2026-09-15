# Surface partition V2 修改日志

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
