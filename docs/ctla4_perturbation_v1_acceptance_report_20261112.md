# CTLA4 表面扰动算法 v1：验收与算法报告

日期：2026-11-12

## 验收结论

- **工程验收：PASS。** `ordinary_site_first_perturbation_v1` 已在真实 CTLA4 产物上完成几何、证据、预算和展示回归。
- **科学有效性：未建立。** 表面扰动场、表位加权和扰动效用均为未校准的计算代理，不是结合自由能、结合损失或实验效应。
- **湿实验验证：无。** 本报告只确认算法执行和证据链，不确认候选抗体的真实结合、稳定性或功能。

本次真实运行输入为 UniProt P16410、结构链 A、`strict` 模式、目标抗原、最多 4 个构建体。协议参数为 100 次 SPURS 唯一请求、30 个 FoldX 槽位，显示区域上限为 30 Å 诱导图测地直径和 22 Å Cα 直径。

## 验收证据

| 检查项 | 结果 |
| --- | --- |
| 核心/改动/API/展示/导出回归 | 456 passed, 2 skipped，38.27 s |
| 完整预览管线 | 4 passed，219.43 s |
| Node 语法检查 | `web/app.js`、`web/recommendations.js`、`web/job_center.js` 全部通过 |
| Job-center 单元测试 | 2/2 passed |
| 表面扰动浏览器回归 | PASS，拦截 3 次提交，页面错误 0 |
| 真实产物几何/证据校验 | PASS，11 区、74 节点、4 个选中构建体、108 个可解析引用 |
| 真实产物预算校验 | PASS，8 个完成队列、25 个预算耗尽队列 |

真实产物位于 `workspace_data/ctla4_budget_gap_validation_20261112/`。两份独立只读校验均通过：最大诱导测地直径为 23.331 Å，最大 Cα 直径为 20.0209 Å，低于配置上限；区域连通且唯一覆盖 74 个节点；推荐状态为 `FORMAL_RECOMMENDATION`。

运行共形成 36 个候选，其中 7 个通过正式资格门，29 个因 30 槽 FoldX 预算未覆盖而保持未资格。25 个 `budget_limited` 队列均存在于 `unavailable_principles`，没有被错误表示为已完成结果。

## 新算法

### 1. 严格连通分区

`balanced_connected_v4` 在诱导后的可设计残基图上计算连通性和最短路径。区域必须整体连通，且诱导测地直径和 Cα 直径同时满足上限，否则直接拒绝，而不是放宽门限生成展示区域。

相关实现：`joint_design_v4/sectors.py:482`、`joint_design_v4/sectors.py:559`。

### 2. 与默认 v2 隔离的扰动协议

默认协议仍是 `ordinary_site_first_v2`。只有请求 `ordinary_site_first_perturbation_v1` 时，候选才携带 `epitope_weighted_surface_perturbation_v1` 选择器标记，因此新算法不会改变默认 v2 的行为。

相关实现：`joint_design_v4/site_first.py:20`、`joint_design_v4/site_first.py:191`、`joint_design_v4/site_first.py:458`、`joint_design_v4/site_first.py:674`。

### 3. 精确三突变证据与弹性队列

高阶候选使用可恢复的站点队列。三突变除自身的化学轴、暴露度和表位证据外，还要求三个精确的、保留氨基酸替换的双突变父体全部具备完整证据。三突变裕量为三突变扰动效用减去最强精确父体效用；缺证据的候选不能通过该门。

相关实现：`joint_design_v4/regional_perturbation.py:134`、`joint_design_v4/regional_perturbation.py:167`、`joint_design_v4/regional_perturbation.py:193`。

当 FoldX 槽位耗尽时，队列标记为 `budget_limited`，并记录未探索原因。预算不会伪造候选分数，也不会把未评分候选当作失败候选。

### 4. 物理表面场与表位加权表面场

对每个三突变构建体，算法在 `geodesic_12` 表面支持上生成逐表面点的两个通道：

- 物理通道：三突变扰动效用传播到支持残基，面积加权后得到 `physical_field_gain_area_a2`。
- 表位通道：每个目标表面残基的局部场乘以表位先验，得到 `epitope_weighted_field_gain_area_a2`；缺失先验使用中性值 0.5。

场状态明确记录为 `SCORED`、`PARTIAL` 或 `NEUTRAL_MISSING`。这些值是表面支持代理，不是分子动力学场或实验测得的变构效应。

相关实现：`joint_design_v4/spatial_support.py:107`、`joint_design_v4/spatial_support.py:141`、`joint_design_v4/spatial_support.py:158`。

### 5. 按证据层级选择构建体

选择器使用固定优先级，而不是把不同量纲的指标压成一个不可解释的分数：

1. 是否来自尚未覆盖的区域；
2. 表位加权表面场增益；
3. 物理表面场增益；
4. 三突变相对精确双突变父体的增益裕量；
5. 扰动效用；
6. 可设计面积覆盖；
7. 边际面积；
8. 重叠度；
9. 信息增益；
10. 目标函数增益；
11. 候选 ID。

达到最小构建体数量且没有新区域候选时，只有仍能增加面积覆盖、物理场或表位加权场的候选才继续进入选择。更新后的预算曲线保存每一步的候选、场增益、状态、扰动效用、区域首增标记和三突变裕量。

相关实现：`joint_design_v4/identifiability.py:89`、`joint_design_v4/identifiability.py:145`、`joint_design_v4/identifiability.py:153`、`joint_design_v4/identifiability.py:182`。

### 6. 证据导出与页面展示

当前导出代码已将覆盖面积、物理场增益、表位加权场增益、两个场状态、扰动效用、区域首增、三突变裕量和三突变角色写入 `panel_budget_curve_v411.csv`。当前页面在检测到扰动选择器时，显示选择依据和逐构建体场增益表。

相关实现：`joint_design_v4/reporting.py:238`、`web/app.js:1436`、`web/app.js:1465`。

## 本次真实运行的选择结果

| 顺序 | 候选 | 区域 | 突变 | 三突变裕量 | 扰动效用 | 物理场增益 (Å²) | 表位加权场增益 (Å²) | 新增覆盖 (Å²) | 区域首增 |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | `VS003\|V41_55371857cde180cd` | VS003 | S50F / R51F / D153W | 0.439433 | 1.872479 | 1201.353 | 1099.617 | 369.148 | 是 |
| 2 | `VS001\|V41_ecc2d8426f45f117` | VS001 | V40W / A41R / A44R | 0.075424 | 1.135717 | 916.940 | 685.031 | 170.170 | 是 |
| 3 | `VS002\|V41_d331294ca903f41a` | VS002 | I53Y / S55F / Q117F | 0.278320 | 1.085682 | 588.451 | 555.192 | 206.518 | 是 |
| 4 | `VS003\|V41_00979e1845c1ca82` | VS003 | S50F / I152K / D153I | 0.401367 | 1.554256 | 0.000 | 0.000 | 89.891 | 否 |

四个候选均为三突变，角色均为 `ENHANCEMENT`，推荐类别均为 `STABILITY_FALLBACK`。前三个分别覆盖 VS003、VS001 和 VS002 三个新区域；第四个没有新增区域，也没有新增物理场或表位加权场，只增加 89.891 Å² 可设计面积和相关信息，因此其科学优先级低于前三个。

该第四个构建体来自固定的目标数量 4：在三个区域都已覆盖后，选择器仍要求达到目标构建体数，且该候选具有正面积覆盖。它符合当前算法规则，但不能表述为新增表面扰动覆盖。其净目标函数变化为 -0.051658，已在预算曲线中保留。

## 验收边界

- 本次验收证明算法可运行、几何约束生效、证据链可解析、预算耗尽可追踪，不证明新算法相对默认 v2 具有更高的实验命中率。
- 29 个未资格候选不是“失败候选”，而是没有被 30 槽 FoldX 预算探索；需要更多预算才能判断。
- 表面场、表位先验、扰动效用的尺度和权重尚未用实测结合数据校准。当前排序只适合计算筛选，不能作为结合损失预测。
- 真实产物生成于最新 CSV 字段和页面明细表合入之前。其 `panel_budget_curve_v411.csv` 仍是旧 12 列，离线页面也不包含最新明细表；历史产物未被重写。当前代码的扩展导出和页面渲染已由回归测试覆盖，下次运行会写入新字段。
- 项目全量测试收集仍排除 5 个既有历史阻断文件：缺少 `scripts.audit_p0_p3_batch`、缺少可选 `psutil`，以及 3 个依赖已忽略 AF3 工作目录的测试。`test_ordinary_adaptive_rerun.py` 还依赖当前工作区不存在的历史夹具，单独运行时出现 4 个夹具缺失失败，不涉及本次算法断言。
- 本次没有进行湿实验、结合亲和力测量或独立实验复现。

## 复现入口

- 真实运行产物：`workspace_data/ctla4_budget_gap_validation_20261112/`
- 设计摘要：`workspace_data/ctla4_budget_gap_validation_20261112/sector_recommendations_v4.json`
- 选择面板：`workspace_data/ctla4_budget_gap_validation_20261112/selected_panel_v4.json`
- 算法实现：`joint_design_v4/site_first.py`、`joint_design_v4/regional_perturbation.py`、`joint_design_v4/spatial_support.py`、`joint_design_v4/identifiability.py`、`joint_design_v4/sectors.py`
- 导出与展示：`joint_design_v4/reporting.py`、`web/app.js`
