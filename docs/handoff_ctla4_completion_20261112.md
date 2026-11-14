# CTLA4 分区与算法修复完成交接

更新时间：2026-11-12。工作目录：`E:\user\New project\ecd_epitope_scanner`。分支：`codex/joint-v4-ml`，HEAD：`8ae5891`。当前修改尚未提交，工作树还包含此前已有的未跟踪文件；不要整树暂存，也不要改写历史任务输出。

## 结论

本轮指定的算法修复已完成，并有当前代码上的回归和真实模型产物支持：CTLA4 分区改用区内路径直径硬约束；新扰动协议与默认协议隔离；FoldX/SPURS 预算耗尽后的队列状态和导出一致；推荐卡片不再重复嵌入大账本；在线与离线 UI 能解析共享引用并解释预算缺口。

这里的“可用”指工程流程、几何门、证据门和导出链路通过验收。预测分数和空间扰动仍是计算代理量，没有湿实验验证或定量校准，不能解释为已验证的结合改变。

## 主要修改

- `joint_design_v4/sectors.py`：`balanced_connected_v4` 在构建、修复和最终验证中使用严格 induced eligible-region 路径；历史 `connected_bounded_v3` 保留全图语义。保存 CTLA4 回放从 8 区扩展为 12 个连通区，最大区内直径 28.2675 Å。
- `joint_design_v4/site_first.py`、`regional_perturbation.py`：默认协议保持 `ordinary_site_first_v2`；扰动选择和高阶增益门隔离到显式协议 `ordinary_site_first_perturbation_v1`。三点候选要求三个精确、序列合法的父本，逐代表通过门槛，并保留协议与证据来源。
- `joint_design_v4/site_first_validation.py`：FoldX 预算无法继续预留时，将剩余 ready/exhausted 队列标为 `budget_limited`；不再留下无候选却显示 ready 的队列。
- `joint_design_v4/sector_recommendations.py`、`reporting.py`、`exporting.py`：卡片引用共享 `ordinary_candidate_ledger.json` 和 `v4_manifest.json#/worker_runs/{name}`；离线导出包含引用文件并移除会轮询在线 API 的 `job_center.js`。
- `pipeline.py`、`api.py`、`cli.py`、`web/`：显式协议、12/12/12 构造预算和表位先验使用状态贯通；在线/离线推荐渲染共享预算缺口解释。
- `README.md`、`docs/method_joint_v4.md`：正式组合必须有有限的组合 SPURS 和 FoldX 证据；旧 SPURS-only 设置不能替代 FoldX。

## 真实模型验证

最新验收入口：`workspace_data/ctla4_budget_gap_validation_20261112/`。

- 输入：保存的 UniProt P16410、实验结构 1AH1 chain A；不是 AF3。
- 运行协议：`ordinary_site_first_perturbation_v1`；最多 4 个主构建体；100 个 SPURS 请求；30 个 FoldX 槽；分区 30/22 Å。
- 分区：11 个连通区，唯一覆盖 74 个 eligible nodes；最大 induced geodesic diameter 23.331 Å；最大 Cα diameter 20.0209 Å。
- 推荐：4 个选中构建体，108 个推荐 provenance 引用全部可解析并包含在归档中。
- 队列：8 个 `complete`，25 个 `budget_limited`，没有遗留 ready；25 个缺口全部出现在 ledger gaps 和推荐卡片 `unavailable_principles`。
- 科学门：36 个候选中 7 个正式合格；29 个因 30 槽 FoldX 预算未覆盖而保持非合格，没有补分或放宽门槛。

独立只读校验结果：

```json
{"status":"PASS","protocol":"ordinary_site_first_perturbation_v1","regions":11,"nodes":74,"max_induced_geodesic_a":23.331,"max_ca_a":20.0209,"selected":4,"refs":108,"recommendation_status":"FORMAL_RECOMMENDATION"}
{"status":"PASS","queue_states":{"complete":8,"budget_limited":25},"budget_gaps_verified":25}
```

产物大小：

- `analysis.json`：5,189,210 bytes
- `sector_recommendations_v4.json`：2,313,234 bytes
- `selected_panel_v4.json`：30,588,095 bytes
- `ecd_epitope_design.zip`：27,683,591 bytes
- 离线 `index.html`：2,783,935 bytes，不含在线 `job_center.js` 引用

此外，默认 `ordinary_site_first_v2` 的独立同输入/同预算真实运行也完成：11 个区、4 个双点构建体、7 个正式候选。它用于行为隔离对照，不是新旧协议的完全同证据等价证明。

## 回归证据

从仓库父目录 `E:\user\New project` 运行当前代码：

- 本轮修改文件加项目记忆规定的核心 V4/CTLA4/PD-L1 套件：`382 passed, 2 skipped`，37.83 秒。
- 导出变更后的预览管线 fixture：`1 passed`，291.22 秒。
- `node --check` 通过 `web/app.js`、`web/recommendations.js`、`web/job_center.js`。
- `node --test tests/job_center.test.cjs`：2/2 通过。
- `node tests/browser_surface_perturbation.mjs <playwright-core/index.mjs>`：`{"status":"PASS","intercepted":3,"errors":[]}`。

包级全量 pytest 仍被 5 个历史测试文件阻断，原因不是当前算法失败：

- `test_audit_p0_p3_prior_attempts.py` 缺少 `scripts.audit_p0_p3_batch`
- `test_case_coordinator_boundary.py` 缺少可选依赖 `psutil`
- `test_finalize_af3_completion.py`、`test_af3_batch_recovery.py`、`test_af3_final_acceptance.py` 引用已从活动树移除的 `workspace_data/af3_binder_batch_20261017` 文件

共享证据调度回放另有逐队列比较：`perturbation_v1` 的 33 个队列状态和完成候选身份与最新真实运行全部一致；`v2` 的 33 个完成候选身份与修复前真实运行一致，另有一条 `VS001/epitope_support` 期望从错误 ready 修复为 `budget_limited`。这验证限定预算下的调度一致性，不替代科学验证。

## 运行状态

最终检查时 18765 和 8767 均无监听进程，因此没有旧服务需要安全重载，也没有活动作业被中断。需要网页交互时，应按项目正常启动方式重新启动服务，再用新任务目录验证入口；不要复用已经停止的 PID。

历史任务 `workspace_data/1d812aa2f9be/` 和所有旧输出未重写。当前改动未提交、未推送。

## 已知边界

- 表面场、表位权重和扰动排序是未校准的几何/化学代理量，不是结合亲和力或实验效应预测。
- 29/36 候选未获得 FoldX 结果，保持非合格；没有插值、补分或放宽门。
- 全量历史测试因缺失的忽略目录和可选依赖无法收集；当前生产相关回归已通过。
- 真实模型复核证明当前输入、预算和导出链路可运行，不构成湿实验或跨靶点普适性验证。
