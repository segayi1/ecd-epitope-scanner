# 普通 site-first 重复加载复盘与修复方案（2026-11-04）

## 范围与结论

检查普通 `ordinary_site_first_v2` 的当前工作区源码，以及用户停止的 PD1 任务 `workspace_data/pd1_ctla4_site_first_20261103/jobs/9b4d300f3c2b`。此处为只读分析；没有重启任务、改变模型、门槛或历史数据。当前工作区含未提交修改，运行时来源以任务的 `executed_source_snapshot.zip` 和 `run_source_snapshot_v1.json` 为准。启动收据记录源码摘要 `756cbc3c90fc57459c128aa55229e15702b3cd9eeff13e9b30e298d1cc7242f1`。

问题是调度、模型生命周期和配置三者共同造成：探索阶段可能只有少数队首可立即评分；64 是上限；而 PD1 请求未启用常驻 SPURS，每个未命中批次都启动一个 worker 并重新加载模型。我启动时沿用旧 PD1 的执行设置，未预先验证新版探索策略与小批次的交互，是本次执行配置上的失误。

## PD1 实际证据

- 中断收据 `workspace_data/pd1_ctla4_site_first_20261103/stopped_by_user.json`：PD1 未完成，CTLA4 未启动；正在处理的请求不应计作成功。
- `site_first_score_cache/verified_spurs` 中断时有 71 个成功缓存合同：33 个一组合批、34 个二组合批、2 个三组合批、1 个 19 组合批，另一个零组合请求是单点模型。缓存文件数不等于唯一组合数、也不是最终预算账本。
- 对 33 个一组合缓存，`stage_timings.model_load_seconds` 中位数 34.169 秒、合计 1149.232 秒；worker 运行中位数 61.001 秒；组合解码中位数 0.016 秒；WT 编码中位数 0.279 秒。缓存预检查 `fingerprint_seconds` 中位数 4.070 秒。上述分类依据已落盘记录，求和仅说明这些请求里的重复开销，不是完整案例的总运行时间或可直接兑现的提速。
- `execution_trace.jsonl` 记录后段连续 SPURS 阶段，间隔大致 65–80 秒；数据库进度长期 0.62 因整个联合分析共享阶段进度，不证明卡死。中断前曾有 FoldX 评分；未生成最终候选账本与验收报告。
- `spurs.log` 仅是最后一次 worker 的 stderr，当前实现覆盖写入，不是完整历史；缓存和执行轨迹是主要依据。
- 关键八个源码文件（调度、队列、引擎、SPURS 运行器、SPURS/FoldX worker、stability、缓存）与本例 `run_source_snapshot_v1.json` 的文件 SHA-256 均相符。
- FoldX 已落盘的 18 份批次 manifest 均为 success，含 31 个 success 候选及各 3 份配对记录；批次规模为 8、7 和 16 个单候选批。共享 RepairPDB 执行回执为一次成功、49.666 秒。其余队列和终态仍未知。

## 按执行链检查重复工作

| 环节 | 当前行为及条件 | 风险和建议 |
|---|---|---|
| 队列与批次 | `joint_design_v4/site_first_validation.py`：`SPURS_BATCH_SIZE=64` 是上限；有未探索位置时 `frontier_limit=1`，只合并各活动队列队首，去重、预筛、已完成队列和缓存都可能使实际批次很小。 | 维持新位点优先与公平性；可评估对多个独立位点组的有界前瞻，不能简单填满 64 个同一位点的替换而提前吃掉预算。预取属于计费评分，改变预算分配，需单独版本、对照和验证召回。 |
| SPURS 多点进程与模型 | `integrations/process_runner.py` 默认每次调用 `subprocess.run`；`joint_design_v4/engine.py` 仅在 `runtime_optimization.spurs_persistent=true` 时创建 case 级常驻作用域。未开启时每个缓存未命中批次重启桥与 worker、重新导入和加载模型。此次 PD1 未开启。 | 优先把已实现的受控常驻模式作为新版默认候选，覆盖所有普通 site-first 入口，启动时明确记录实际模式和 worker 实例 ID；先完成故障、资源与科学等价验收，再推广默认。不要误以为填写 64 就消除了冷启动。 |
| SPURS 单点模型 | `integrations/workers/spurs_gpu_worker.py` 的无突变请求走 `get_SPURS`，有突变走 `_load_multi`。普通流程通常一次单点扫描；现有常驻 worker 只缓存 `_load_multi`，不缓存 `get_SPURS`。 | 先统计单点调用次数，若一个 case 只有一次则优先级低；多次场景再考虑单点模型缓存，不要因少量固定开销引入额外显存占用。 |
| SPURS 结构和野生型编码 | 每次 `_predict_once` 都 `parse_pdb`，multi 每请求调用一次 `encode_wt`；现有复用仅在一个批次内部，不跨请求。 | 在常驻模型可靠后，按结构内容、链、配置、设备及模型指纹考虑 case 级复用；核查解码对输入字典的原地修改、设备和显存生命周期，证明数值一致。此次单组合 WT 编码中位数 0.279 秒，优先级低于模型冷启动。 |
| 环境身份与缓存预检 | `surface_perturbation/upstream_cache.py` 每次请求计算代码/配置等身份；非驻留模式每次启动 WSL Python 查询包和 torch，并启动 `nvidia-smi`；任务级模型文件摘要可复用。常驻模式在 worker 内缓存环境身份并在收尾验证。 | case 内复用身份与安全终验，记录预检、模型加载、推理和缓存命中原因。不得以文件名或时间戳单独取代内容/环境身份；终验失败应使证据不可正式采用。 |
| WSL 路径转换 | `integrations/wsl_bridge.py` 中普通及常驻桥在每请求对结构/输出路径调用 `wslpath` 子进程。 | 常驻 case 作用域内可缓存经验证的路径转换；变更路径或环境时失效。先计时，避免把小开销当主要瓶颈。 |
| FoldX worker 与修复 | 每个 site-first FoldX 批次经 JSON adapter 启动一次 Python worker；`integrations/workers/foldx_worker.py` 的 `_site_first_repair` 在同一个 case 按输入、二进制和修复结构摘要复用 RepairPDB，并校验收据；每批仍做二进制/结构身份检查和来源复制。 | 先依据逐批 `repair_executed` 和耗时核算，禁止假设每批都 RepairPDB；可减少重复身份检查/复制，但需保留溯源、短路径处理和并发原子性。新的隔离目录意味着跨 case 复用不能假设存在。 |
| FoldX BuildModel 与线程池 | `stability.py::score_many_with_foldx` 明确每个候选一个原生进程，避免邻近候选改变采样序列；3 槽及必要时追加 7 槽均是实际物理重复。`foldx_execution_scope` 在主进程设置 ContextVar，但常规 FoldX worker 由 `JsonWorkerRunner` 新建子进程；上下文不会跨进程传递，现有单元测试只证明同进程作用域复用。 | 不把候选合并到同一 FoldX 进程，也不能减少/伪造槽数。若要复用 Python 层池，需要在实际 worker 内实现有界生命周期并做跨请求集成验证；先测收益。当前 `foldx_pool` 还嵌套在 SPURS 常驻开关内，独立开关语义需核对。 |
| 其他输入模型 | `joint_design_v4/engine.py` 对 ProCaliper、SPURS single、DiscoTope 分别调用适配器；正常普通流程这些是 case 级输入证据，非每个突变体都跑一次。通用 JSON worker 缺省是每次请求一个进程；DiscoTope worker 内还启动其独立推理进程。 | 对重试、不同入口及缓存失效分别记录次数；按实际消耗排序。Discotope 和 SPURS 的通用缓存不能仅凭声明模型标签判断，现有适配器明确禁用不安全缓存；SPURS 另有内容验证缓存（包括 single），但默认 case 隔离。ProCaliper 通用缓存使用声明的工具身份，扩大共享前应核实已安装版本身份。 |
| 失败、显存、关闭 | 现有 `spurs_persistent_worker.py` 的常驻路径强制 `cuda_oom_retries=0`、`allow_cpu_fallback_on_oom=False`；与单次模式默认恢复语义不同。出错/部分结果会关闭 worker，作用域退出有身份终验。 | 推广前需明确保守失败策略或实现可审计恢复；不能悄悄回落 CPU 或更改原结果语义。异常、超时、取消、用户停止均应回收 Windows 和 WSL 子进程，并明确未完成请求和预算账本。 |

## 建议实施顺序与验收

1. **先补可观测性和启动合同。** 日志逐批记录队列数、探索开关、实际组合数、命中/请求数、worker 实例、冷/热加载、环境指纹、阶段耗时、异常及成本；任务级状态分离执行中、已完成和中断。启动时输出已解析的执行配置，而不是只显示版本名或 64 上限。把持久 worker stderr 改成可追踪记录，不能覆盖旧批日志。
2. **让 SPURS multi 在单个案例内稳定常驻。** 复用现有受管 worker，但先补上默认入口传播、后端支持检查、错误恢复/不支持时明确失败或显式告知、结束身份验证、取消清理。保持独立 case 与 GPU 资源边界；单点与多点各自明确生命周期。现有小型真实握手证据 `workspace_data/site_first_runtime_opt_20261029/smoke_03/receipt.json` 显示一次冷加载 65.836 秒、第二次热加载 0.008 秒，但不证明当前完整任务已验收或能达到同样总体加速。
3. **单独处理调度与剩余开销。** 以当前策略为对照，先保持候选选择顺序不变，测常驻效果；然后才比较有界、跨位点的预取策略。优化 WSL 路径、身份预检、结构/WT 编码与 FoldX Python 生命周期前分别测量；不为填满批次牺牲位置覆盖、队列公平或 1000 SPURS/300 FoldX 槽合同。
4. **构造可重复测试闭环。** 用模拟 worker 驱动真实调度：同一 case 多个小批次只启动一次、请求顺序及 ID 匹配、预算与决策不变；两 case 隔离；失败、超时、终止、模型变更、显存不足、缓存异常和子进程清理均可检出；检查真实 FoldX 子进程边界，验证池是否确实复用。必要时把请求回放与原始结构/链/AF3 confidence 保持绑定。
5. **做受控真实验证后再运行整批。** 冷/热启动各重复测量；固定输入、来源、预算及硬件，对照逐批请求/分数、区域与原则终态、科学/展示资格、3/10 槽配对和导出；报告队列覆盖、新位点触达、时间、GPU/CPU 内存和磁盘，并把新旧采样差异与预取改变单独列出。对 SPURS 数值给预先约定的容差和跨门复核；FoldX 时间种子不能要求逐槽完全相同，应核验槽身份、资格和采样波动。单批短测不能代替两个完整真实案例。

FoldX 的修复缓存还应处理结构文件存在但 `repair_identity.json` 缺失的中断情况：当前 site-first 会读缺失收据并失败，通用函数则只看文件是否存在。建议将修复产物及收据作为完整原子项发布，校验失败时重新计算，不把文件存在等同有效。site-first 分支在通用配对缓存检索前返回；如要支持同槽位中断恢复，需要另行实现完整身份核验及明确的槽位和预算账本。

## 仍未知

中断任务没有完整账本、最终面板或实验真值，无法估算完整 Pareto 前沿召回或最后总耗时。现有缓存按完整批次身份键控，拆分/合并批次或切换持久协议通常改变缓存合同；不能承诺中断的 71 份缓存都能被新模式直接复用。任何恢复必须验证结构、模型、运行环境、协议、顺序和已消费预算，未完成请求单独记账；不要覆盖原始中断证据。
