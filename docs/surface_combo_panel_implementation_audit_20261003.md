# 2–3点真实组合面板实现与四例验证审计

2026-10-03，用户授权独立组合批次；历史单点局部结果不作为本目标成功结果。

## 根因与实现

- 旧物理预览强制orders=[1,2]、单点占FoldX名额。`combo_2_3_hard80_v1`限定最终4–12项，每项2/3个不同位点、每个位点真实组合五重复同组4/5直接变化；未确认候选是待建模，不是已不稳定。
- 重叠窗口全局去重，跨度沿完整可靠图，24 Å上限保留。FoldX按新增面积代理分散排序、双点优先、少量三点及全部精确双点对照。代理仅排序，不授予覆盖资格。
- 三点必须有三个精确双点子集的真实五重复、匹配校准与局部解析，再比较真实场的现有5%额外面积。解析为零变化的双点仍是合法对照，不因此成为目标构建体；单点并集不作证明。
- site_changes按native author-chain/auth-seq/insertion身份绑定位点。重复位点、任一位点无真实贡献均拒绝。goal_status与search_status独立，80%代表计算设计目标；VERIFIED_PANEL仍要求原搜索闭合。
- null缓存导入原始样本/误差/失败mask，重新审核当前资格，导入hash清单进入双层导出。合法模型与原始场复用，旧资格标签不继承。正式DDG只用对应最终前五对，null扩展20对不替代正式重复。
- 回归发现计时`carried+clock-started`浮点抵消可低于已累计值；改为`carried+max(0,clock-started)`，历史时间不退还。
- APBS只服务已达到80%的最终嵌套集合，保持五对、统一求解、代表收敛与未知合同；中间不足目标面板不算电势。页面/API展示组合、逐位面积及目标/搜索状态，区分局部候选与选中面板。

## 批次与验证证据

新批次 `workspace_data/surface_combo_four_case_20261003/batch_20261003_155853`，保留四例原输入/映射/confidence，EGFR排除。每例校准14400秒、生产7200秒，同批重试累计，旧账本保留。A队列CTLA4→PD1，B队列PDL1→B7-H6，已有受限工具并行保留。

v1在运行校准不覆盖；退出阶段边界后源hash guard阻止混代。handoff_v2必须等待旧进程退出、账本关闭，且当前校准完整交付才能直接继续production；同预算不重置。

v2隔离检出 `E:/user/New project/surface_combo_verification_20261003_v2/ecd_epitope_scanner`：detached HEAD 9f01e12加261文件显式overlay，并非clean HEAD声明。项目venv从检出父目录169项科学/V4/API/FoldX/UI测试全通过（14.57秒）；已有Starlette warning。浏览器PASS，2次写请求均拦截，未创建测试任务。

回归覆盖不同位点、逐位贡献、真实低阶对照缺失、零变化对照、未建模和真实拒绝区分、跨显示边界路径、代理多样性、原始null再审核、真实非加和、目标与搜索状态独立。

另补首轮组合校准的同协议20次上限调度，已达上限的原始null仍未知、优先处理未尝试任务；首轮/扩池两种路径回归均通过。安装包遗漏surface_perturbation包已修复，wheel构建及42模块逐字节匹配PASS，证据位于本批package_audit/final_v2。

## 当前限制

记录时CTLA4/PDL1校准仍运行，已产生真实组合场；PD1/B7-H6随后执行。尚无四例完整目标面板，不宣称达到80%。最终清单、ECD、五对模型、逐位证据、K曲线、闭环和双ZIP必须依本批实际产物核验。预算耗尽或20次null未稳仍为未完成，不新增批次规避限时，不放松资格。


另16项AF3保存输入/全流程/SASA回归通过（180.67秒）。队列完成后只读finalize_campaign执行双层ZIP、目标审计、ECD/组合/K曲线导出及四例native浏览器；输出final_combo_acceptance.json/md。此处记录待执行，不预先宣布PASS。

## 实时监控（2026-10-03 17:34）

用户要求持续监控与记录问题。本批monitor_campaign.py在独立隐藏进程中每30秒读取四例阶段、累计预算、进程身份/CPU及产物更新；不改运行源码、阈值、缓存或预算。monitoring/status.md与latest.json更新现况，snapshots.jsonl保留时间线，issues.jsonl去重记录问题，evidence/保存原文与SHA256。记录历史继承失败与本批新失败，阶段handoff的FROZEN_SOURCE_CHANGED记INFO；无产物15分钟仅提示疑似慢任务/停滞，不自动杀进程或断言死锁。验收finalizer结束后监控自动停止。

已核对：CTLA4首轮12项null均达到20次，0项local_scope_converged，CALIBRATION_INSUFFICIENT（校准任务运行结束不等于校准成功）；生产已在v2继续。PDL1本批三点L27E/Y123W/R125L（V41_2e40c586a0687ea9）null任务记录SES mesh is not closed，不能推断突变本身不合法，需进一步网格拓扑诊断。两项原始证据已保存，监控进程存活、多轮采样、事件去重及证据hash核验通过。尚无四例目标验收。


## 2026-10-03 18:16：资源中断后同批并行恢复

lane_B_v2.log与launcher记录明确OSError(28)/Windows错误112磁盘空间不足，随后物理ZIP读取BadZipFile；不能将该包视为完整交付，亦不能据用户描述断言RAM OOM。恢复前C盘约31 GB、E盘24 GB可用，WSL根162 GB可用，主机可用内存约17.8 GB；未发现仍在运行的科学worker。旧pid19632已复用，未按旧pid杀进程。

用户明确要求并行，保留CTLA4→PD1与PDL1→B7-H6两条队列及现有工具限流。原v2冻结源码261项hash核验匹配，历史metadata/results/验收输出复制至同批resume_history/20261003_181600。四例仍沿batch_20261003_155853原累计账本，从校准恢复以重建损坏/不完整交付，合法缓存复用，20次同协议null上限不重复，阈值/重复/预算不变，不新建批次。没有额外删除文件。

本批monitor_campaign.py仅增加资源采样及低内存3 GiB/低磁盘5 GiB事件，每30秒保存monitoring/resources.json；语法检查通过。科学计算源码未改。finalizer等待两队列退出后重审当前双层ZIP及真实目标证据。两条launcher已经RUNNING并实际输出FETCHING/ANALYZING。恢复启动不是目标验收；null不足和SES不闭合仍由原科学合同报告未知。低资源事件仅记录，不承诺自动扩容或证明未来不会耗尽。

资源监控原生复现本机psutil.disk_usage的SystemError(impossible bad format char)，标准库shutil.disk_usage正常；监控仅替换磁盘读取接口，完整snapshot重放通过后重启。科学worker未中断、源码和预算未改。

2026-10-03 20:52用户集中CTLA4：核验command后停止旧两队列及其worker，关闭账本不退款；其他病例不排队。CTLA4校准已完成当前交付，直接恢复production，冻结科学源码/配置不变；现有工具并行保留、CTLA4 worker允许全部逻辑CPU/AboveNormal优先级。外部launcher、monitor、验收/浏览器仅CTLA4，附K=8与实际项数，旧报告/metadata保留。cal累计13472.943秒，prod停止时2415.829秒，原14400/7200上限不变；目标尚未验收。用户允许精简非必要要求：优先hard目标阶段性结果，不以搜索未闭合/APBS未解析否认已完整验证的hard证据，但不将阶段性结果冒称全流程验收，不放松逐位/稳定性/null/4of5资格。
