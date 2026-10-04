# 表面扰动运行问题根因核验（2026-10-02）

## 结论与边界

已确认并修复实际查询范围缺失、Windows短暂文件占用导致原子替换失败，以及后台计算未申请系统保持运行的问题。状态机此前的局部校准、候选池闭环与APBS状态修复继续有效。修复通过隔离检出回归与保存的真实配对模型重放核验，不能据此宣布五例科学验收完成。

原7200秒生产/14400秒校准累计预算、全部历史模型、原始场、旧ZIP及失败文件保留。新的工程重放单列账本，不用于补足已耗尽的正式预算。当前五例状态以 `workspace_data/surface_perturbation_20261001/frozen_benchmark/current_acceptance.json` 为准。

## 1. 数小时wall-time间隔：Modern Standby

System日志的Kernel-Power 506/507确认：2026-10-02本地02:13:31进入Modern Standby（SC_MONITORPOWER），03:19:25有内部退出/再进入，09:02:57因键盘输入退出。02:13:32的Kernel-General时间调整为0ms。此前只检索42/107而未检索506/507，未能识别这段待机；现在有事件证据，不能将数小时累计wall-time解释为APBS或FoldX计算耗时。

`surface_perturbation/power_guard.py`在活动任务期间通过PowerCreateRequest申请SystemRequired及ExecutionRequired，不要求屏幕常亮、不改电源计划。`pipeline.py`和外层`scripts/surface_phase_guard.py`在完成或异常后释放请求；原累计时间和超时监督仍保留。

原生Windows请求注册/释放烟测通过。请求不保证覆盖手动睡眠或电池Modern Standby策略，详见[Microsoft PowerSetRequest文档](https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-powersetrequest)。`powercfg /requests`读取权限不足，因此不声称已用其核验。原事件XML存于 `frozen_benchmark/posthoc_audits/power_events_20261002.json`。

## 2. WinError5：原子替换被不共享DELETE的读句柄阻挡

使用保存的PD1配对manifest副本，实际CreateFileW读句柄允许READ/WRITE共享但不允许DELETE共享，`os.replace`复现同样的WinError5/errno13。关闭句柄后写入成功，原PD1文件hash不变。历史运行未记录句柄所有者，不能指认具体UI、杀毒软件或进程。

`util.write_json`已有独立临时文件、关闭写句柄、有限短暂重试及保留旧完整JSON的处理；新增真实Windows句柄回归确认锁释放后重试成功。持续占用仍必须报错，不能静默丢失manifest或无限等待。复现报告为 `posthoc_audits/runtime_cause_replay/report.json`。

## 3. 校准局部NaN：遗漏查询与真实未知必须分开

全量检查保存的20次局部场：PD1最差实际WT查询面积比例99.2272%，PDL1旧cohort为95.8446%，B7-H6实际查询均完整。查询是否执行以总是定义的WT极性支撑分量判断；不能以化学数组每分量有限判断，否则会把His质子化或供受体方向的主动弃权误认成未查询。

`surface_perturbation/scoped_fields.py`在实际候选direct-support并集内补做缺失的配对WT/Mut几何和化学查询，同一骨架对齐用于两者。结果写独立scoped缓存，不改原始raw数组；已有拓扑拒绝继续保留NaN。缓存验证参考、两模型、实现、查询mask、原始数组和NanoShaper身份。

`engine.py`使用新query合同及独立null协议，不能沿用旧threshold。5→10→20的请求范围扩大时，数值误差查询mask进入缓存身份并重新审核；只减少范围外查询，不减少真实配对重复或更改SES网格。导出字段审计记录范围hash及每次补算点数。配置为 `frozen_settings_20261002_v3.json`，方法合同为 `paired_WT_native_SES_preview_1.2`，预算不变。

保存的PD1候选 `V41_fac41aa107967037` 20对真实模型重放：实际查询完整率99.2272%→100%，25.315秒，全部原场hash未变；无新增FoldX/APBS/PDB2PQR求解，NanoShaper网格合法缓存复用。这仅验证查询分配修复，不证明null收敛或最终面板。PDL1最差旧cohort不在当前paired manifest，严格拒绝猜测配对，未把它宣称为重放通过。

证据：`posthoc_audits/null_query_scope_20261002.json`、`posthoc_audits/PD1/query_scope_repair/report.json`。诊断及重放脚本分别为 `scripts/audit_surface_null_scope.py`、`scripts/replay_surface_scope_repair.py`。

## 4. 保留的科学未解析事项

- His质子化、部分氢键方向投影仍未知；有限数据的Q95不稳定也不能靠补查询宣称收敛。
- PD1一组20次参考均有拓扑改变；新增查询不得解除全局拓扑拒绝。保存失败OFF包含真实边界/非流形边，不属于仅翻面可修的问题，不填洞制造合法SDF。
- PD1旧8项APBS执行成功但准备审核均不通过，存在非突变重原子变化，如A107 ND1约2.245 Å。准备差异和未证明的数值误差保持unknown；不把代表项误差界转移给全panel。证据为 `posthoc_audits/apbs_unresolved_causes_20261002.json`。
- 旧PD1/CTLA4含违规Cys构建体，继续失效；五例均未正式验收，不能恢复旧覆盖数字。

## 5. 回归

最新版源码覆盖到独立检出（249文件）。从检出父目录使用项目venv执行pytest且禁用cacheprovider：127项表面场、候选资格、状态、API、AF3映射、原子JSON和电源请求相关回归全部通过（30.92秒）。测试含实际缺失查询、unknown保留、拓扑拒绝、查询范围扩展缓存失效、原数组不变以及真实Windows文件锁。

流程fixture、Web设计与SASA上下文另12项通过（300.96秒），本轮共139项通过。源码身份及最终结果见 `posthoc_audits/runtime_repair_regression_20261002.json`。测试通过表示实现与所测反例一致，不代表抗体功能验证或五例生产验收。
