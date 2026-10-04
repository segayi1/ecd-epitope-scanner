# PD1 AF3 分区修复审计

## 输入与故障

- 原任务：`fa85119485cd`，AF3 upload，Q15116，ECD 25–170。
- 原结构与置信度文件保持不变；具体 SHA256 见重放 provenance。
- **直接造成远处分区着色的是导出 label 编号错位**：ECD CIF 的 auth 25
  对应新 label 1，旧页面却用 record 的原 AF3 label 25 选中 auth 49。
  本输入因此偏移 24 个位置。现改用保持一致的 author 标识。
- 原 VS004 / VS007 超径仍输出；分区合并未检查几何资格。
- 原 VS001 跨序列段并非编号映射错误，但包含低置信度末端；VS002 显示
  尺度达 25.029 Å，超过页面固定文案宣称的 18 Å。
- reliability 只将明确 IDR 从几何中排除，未区分端部三维位置不可靠。

## 代码变化

| 文件 | 变化 |
|---|---|
| `structure_reliability.py` | 独立的 AF3 末端三维资格；保留 UNCERTAIN 分类、原坐标和序列审核；显式记录 omitted 与 boundary |
| `surface_patch/engine.py` | 消费 `use_for_3d_footprint`，防止无资格位置进入 footprint 和 V4 图；保留其遮挡原子 |
| `joint_design_v4/sectors.py` | 拆分不连通 / 超径分区；合并前检查并集；区域数为软目标；无有效 merge 时保留小区域 |
| `config.py` | 显示紧凑性直径默认值：测地 26 Å、Cα 18 Å |
| `joint_design_v4/engine.py` | 显示政策 provenance；实际 display 文件与 Surface V2 区域分开；额外区域颜色；同步 record.region_id 和保护位点计数 |
| `pipeline.py` | 记录显示政策、eligible 图覆盖口径与几何省略；明确几何验证的实际计算范围 |
| `web/app.js` | 修正导出编号选择；未归属位点不叠加 SASA / folding 颜色；去掉不适用的 Fab 描述；空 mutable 列表保持 0 |
| `web/renderer.js` | 等待真实结构 hierarchy / assembly 后再首次着色 |
| `scripts/replay_joint_v41_run.py` | 保留 AF3 类型、原模型、confidence、strict 标志；增加显式 geometry-only 模式 |

26 / 18 Å 是可配置显示范围约束，不是表位尺寸定律；科学 hypothesis 仍
来自 Surface V2 的多尺度重叠足迹。序列连续不是三维空间连通的替代判据。

## 实际验证

`workspace_data/pd1_af3_partition_repair_20260930/verified_geometry_complete`
保存完整几何重算与导出。此次评分 worker 显式禁用，不是一次正式推荐重跑。

| 指标 | 结果 |
|---|---|
| 原始模型 SHA256 一致 | 是 |
| 原始 ECD 坐标保留 | 146 / 146 |
| 不可靠端部退出三维 | 25–31、150–170，共 28 个 |
| 保留位点边界 guard | 32、33、148、149 |
| 有效表面图节点 | 58 |
| 显示区域 | 7 |
| 最大测地直径 | 25.226 Å |
| 最大 Cα 直径 | 17.885 Å |
| 每区 induced graph 连通 | 全部通过 |
| graph 节点归属唯一且完整 | 通过 |
| hypothesis 含不合资格位点 | 0 / 172 |
| 相关回归测试 | 70 passed |
| manifest / 磁盘 / ZIP 一致 | 47 个文件；最终 analysis 也一致 |

测试 fixture 保存原始坐标、PAE 与故障表面图，新增 6 个回归场景；同时
运行原 V4、Surface V2、可靠性、hypothesis panel 与 CTLA4 故障测试。
Python 编译、JavaScript 语法检查及 diff 空白检查通过。
实际浏览器同时验证 7 区、58 个归属位点、结构 overpaint 和无 JavaScript
错误；保存 `browser_verification.json` 与 `pd1_repaired_partitions.png`。

原任务目录未改写。验证副本登记为页面任务 `da27bc835241`：
`http://127.0.0.1:8765/?job=da27bc835241`。
本地服务已重启，后续提交使用修复代码。

机器可读结果：`workspace_data/pd1_af3_partition_repair_20260930/validation_summary.json`。
持续维护约束与此前修复见 `docs/project_memory.md`；根 `AGENTS.md` 引导后续
开发先读取该文件。

## 边界

核心内低置信度短环保留；实验 B factor 不解释成 pLDDT。没有合法 PAE
映射不能由本条末端规则排除位置。高置信度多结构域的相对运动仍需要
PAE / ensemble 审核；本次修复不宣称所有 AF3 构象均可靠。

正式组合仍需要独立的真实评分与现有门槛；几何验证的空推荐是主动禁用
worker 的预期结果，不能作为新的模型失败或“蛋白没有可用突变”的证据。
