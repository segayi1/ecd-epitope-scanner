# 四例完整运行与目标核验（2026-10-02）

用户本轮要求审查全部目标，并完整运行CTLA4、PD1、PDL1、B7-H6；EGFR明确排除。
冻结计划及对话为验收依据。旧计划中的讨论状态不覆盖用户已授权的开发/验证请求。

## 起始目标审查

| 目标 | 起始核验 | 缺口或边界 |
| --- | --- | --- |
| paired WT/Mut、canonical SES、SDF/occupancy/化学场 | 已实现并有保存输入重放/回归 | 固定骨架、拓扑/方向未知限制保留 |
| SASA/SPURS复用、单/双点为主、共享保护政策 | 已接入 | 旧含Cys面板失效，不能使用旧正式覆盖率 |
| 每通道null、5→10→20、最终4/5一致 | 已实现，局部范围修复已验证 | 所有声明通道的完整校准未在四例达成；上界关闭仅证明hard盲区问题 |
| 真实组合复算与重新优化 | 已实现并测试 | 四例尚无本轮完整合格面板 |
| APBS仅最终集合、配对固定网格/准备审计 | 已实现并有真实求解 | 执行完成不等于准备/数值/extended解析完成 |
| 四项覆盖、嵌套K曲线、盲区、冗余、独立K bounds | 已实现 | 未完成病例不得报告正式数字；unknown仍为null |
| 约8项覆盖绝大多数表面 | 未证明 | 90%仅参考线，不放宽科学门槛制造达标 |
| 欧氏/测地/V4/SASA对照及物理消融 | 工具已实现 | 尚缺四例完整实际结果与参数判断 |
| 三条独立候选路线的同预算物理对照 | 未完成 | 已有同测量池过滤及候选生成对照，不能冒充独立物理搜索 |
| 隐藏接触标签sanity check | 标签已冻结，核验工具已实现 | 尚缺当前合格真实面板后的结果 |
| API/UI/报告/完整证据ZIP | 已接入，有隔离回归和浏览器记录 | 四例新结果仍需实际导出、解压核验 |
| 参数冻结、五例方法充分性 | 未完成 | 本轮排除EGFR；四例不能替代完整五例或抗体功能验证 |

四例原保存输入已独立核对模型hash、结构语义和confidence：CTLA4结构126 aa；PD1原AF3 146 aa；PDL1结构220 aa；B7-H6原AF3 238 aa，链均A。详情为 `workspace_data/surface_perturbation_four_case_20261002/initial_goal_review.json`。

## 预算与历史

旧四例生产累计预算均已耗尽，旧账本、失效面板、模型及失败产物不改。原计划明确重试不能清零，故本轮完整续算的预算选择已向用户提出；未取得该约束前不启动会改变旧限时合同的依赖计算。

为避免新批次重复生成已有证据，新增显式`raw_field_cache_roots`入口。复用必须通过相同参考、数值实现、NanoShaper及内容hash审核；只迁移原始物理数组，不迁移旧threshold、校准完成声明或面板资格。FoldX已有严格fingerprint和缺口追加机制保持。

各目标以当前源码、当前运行、累计账本及导出一致性联合验收，不能用回归通过替代实际完整运行。

## 本轮预算无关验证

- 新增跨批次合法raw缓存入口后的75项相关回归通过，另从249文件源码覆盖的隔离检出验证；浏览器提交、12项曲线/unknown呈现及legacy独立上限检查通过，无page error。浏览器请求被截获，不创建真实病例。
- 对旧交付目录复查：CTLA4旧物理ZIP漏带原DDG表，源配对本身通过；用当前导出器另目录重打包后1127项/90对全部通过，原manifest hash不变。该工程输出明确标为历史证据，不授予新面板资格。
- PD1旧外层ZIP与失败续跑后改写的启动输入清单不一致，是历史包和新尝试混用，不重新授予旧panel资格。PDL1旧包一致性通过；B7-H6缺完整生产交付包。结果存于 `historical_delivery_audit.json`，其中ZIP内历史status不能冒充当前病例验收。
- 以上不计作四例新生产验收、完整校准、独立路线物理对照或新的正式覆盖率。

## 用户确认后的独立重跑

用户随后明确确认修复并要求“重新跑四个样例”。据此在 `workspace_data/surface_perturbation_four_case_20261002/rerun_20261002_191524` 启动独立新批次；原账本保持不变，新批次仍使用7200秒生产/14400秒校准、同批重试累计计时。不是为旧失败病例清零或追加旧验收时间，也不改变科学资格。

启动前127项相关回归、2项冻结源码保护回归通过，250文件源码覆盖检出。两条独立队列为CTLA4→PD1和PDL1→B7-H6；`--frozen-source`在每阶段启动前逐文件核对源码及检出hash，不并发覆盖共享源码。GPU沿用既有锁，APBS按现有每例最多2和内存约束。

当前进度看该批次 `lane_A.log`、`lane_B.log`、分病例累计账本、`lane_A_results.json`/`lane_B_results.json`。启动不代表完成；最终状态以新启动manifest和当前summary联合审核，不用历史面板补缺。

### 重跑中发现的依赖路径问题

首轮CTLA4的FoldX失败：默认二进制绝对路径来自隔离检出的config，而该检出故意不复制vendor二进制，故未获得配对模型。原模型/工具没有缺失或科学失败；这是benchmark依赖解析缺陷。实际错误日志和已消耗账本保留，停止已核验的本任务进程树、关闭活动账本，不退还时间。

`benchmark_surface_perturbation.frozen_worker_parameters`现在在开始上游计算前将FoldX明确解析到保存输入的原工作区，配置缺失/路径错误/hash改变均快速拒绝，并将依赖路径及hash写入frozen_input_manifest。两项路径/hash回归和冻结源码、病例状态、phase guard共13项隔离通过。实际二进制hash为 `a5d18d39fef02fdf86ccd7eead8644fb9c2965a70d7d087c1823ca9adff0e844`，与合法历史缓存一致。

`settings_v5.json`单列修复身份，两条队列在原新批次账本内重启，日志为 `lane_A_v5.log`/`lane_B_v5.log`；不是另开预算。源码覆盖251文件。原失败配置和停止进程/账本记录见 `dependency_failure_stop.json` 和campaign_manifest.repair_history。

### 实测补充的重复追加路径回归

v5实测两例均进入42条真实配对模型后，校准补至20次再调用ensure已有重复，触发UnboundLocalError：repeated_field只在追加循环中局部导入，而相同重复数时循环不执行。engine改为模块级导入。新回归真实经过校准控制流，覆盖5次及5→10→20并再次使用相同重复数，明确无新建模型、缓存按坐标共享但保留真实run身份。

同时修正coordinator对旧外层ZIP的判断：新summary不意味着旧run_manifest/ZIP是新交付，两者各自mtime必须属于当前尝试；旧包输出UNAVAILABLE_STALE_OR_MISSING_BUNDLE。case_status接入当前driver失败状态和新启动后旧summary拒绝，防止物理VERIFIED但完整导出失败仍被接触审计视为合格。

87项相关隔离回归通过，252文件冻结。settings_v6与日志lane_A_v6/lane_B_v6在同一账本内恢复，数值七模块及既有null协议未改；合法当前模型/场/null保留。停止记录为repeated_field_failure_stop.json，故障耗时不退款。

### v6端到端工程验证

在252文件冻结检出的父目录、原项目venv运行V4/SASA-context/API/页面/完整流程fixture/PD1 AF3回归，53项通过（360.62秒）。唯一warning是已有Starlette httpx弃用提示。此结果核验代码路径及导出一致性，不代表四例科学面板已合格；两条真实队列继续使用v6及原新批次累计账本。

### v7原生工具路径与预算读竞争

CTLA4候选A30Y/K29W补算失败保存输入复现：FoldX在239字符cwd返回3221226505（0xc0000409），无分数；232和87字符cwd各4次成功（82.18/82.86秒）。该原生二进制有路径buffer限制，不能当作突变不稳定。worker现在为过长路径创建持久化短执行目录，模型manifest仍写原输出目录、记录真实短路径和工具hash；原文件不删除或移动。默认benchmark显式使用保存输入工作区的foldx_short，合法RepairPDB文件按hash复制并单列种子记录。

实际guard还遇到Windows预算读PermissionError，外层退出但CTLA4 child继续。新增read_budget有限共享冲突重试，永久监控失败停止自己启动的进程树、计费并输出INCOMPLETE_GUARD_FAILURE。80项相关隔离回归通过；含短路径、读拒绝重试及失败无孤儿进程。253文件冻结，settings_v7、lane_A_v7/lane_B_v7在同一账本恢复，预算和已消耗时间保留。原生重放数据见workspace_data/foldx_path_replay_20261002，停止记录foldx_path_failure_stop.json。

修复worker真实重放使用同一normalized input和相同A30Y/K29W、4次重复，在原长output成功，4对模型/PdbList/DDG完整，原output manifest保留；provenance记录execution_root短目录、binary hash及new_model_counts=4。证据worker_request.json/worker_response.json。该工程重放不计作完整面板或额外null收敛。

v7的PDL1 V41_108bf6b953d9bddb SES任务失败独立核验：147128顶点/294267三角面，3条边界边、0非流形边、0闭合流形边winding矛盾、0退化面，是实际缺面而非仅orientation问题。保存原OFF及hash，见PDL1_true_open_mesh_audit.json；按合同unknown并隔离任务，不补洞制造SDF资格。此数值限制保留，不计为已修复有效物理证据。

PDL1 v7首轮校准已实际完成：物理ZIP 2426项/235对模型/187303038未压缩字节审计PASS，完整流程ZIP 66项审计PASS；科学状态CALIBRATION_INSUFFICIENT，不能称已完成正式面板。队列自动进入production。当前paired/cache及新旧任务历史都已进入本次证据导出。

### v8扩池调度：不重做已到上限的匹配null

PDL1第二轮真实校准12项有10项已在首轮到20重复，仅2项新任务。CTLA4第二轮物理包仍236对，说明这次扩池未增加模型。原调度始终优先同一proxy面板，且coordinator把请求集合不变视为无进一步任务，可使未处理候选遗留。

新增calibration_schedule：仅在扩池阶段、经过原strict loader指纹审核的同required-scope完整null或19个相邻差分块（20模型）不再机械追加。部分校准保持unknown、继续阻止VERIFIED；范围扩大或10次证据仍需追加。未处理任务先于失败任务重试。summary明确selected/remaining/deferred及扩池scope是否完整，空追加不把未收敛的deferred变成COMPLETE。coordinator依据尚未处理且当前仍请求的任务继续，即使请求集合相同；全部已尝试且匹配20仍不足则停止无信息追加、保持CALIBRATION_EXPANSION_REQUIRED。

91项隔离回归通过，新增包含真实engine控制流的20次partial不重做也不宣布完整测试。保存实际数据重放：CTLA4 32请求/20 fresh/12 deferred；PDL1 35请求/23 fresh/1失败重试/11 deferred。v7已运行阶段按原255之前快照完成并导出，其后source guard主动阻止代际混用，无孤儿child；重新冻结255文件，settings_v8在原账本继续。数值七模块与阈值合同未改，模型和原始fields/null保留，详见各case_calibration_schedule_replay.json。

### 校准缓存对应源模型的恢复与导出保留

扩池不再重算max20任务后，仍须保留该null对应的源模型；不能用当前上游1次记录覆盖旧20次身份。新增resume_pairs按相同mutation、repaired input、normalized input与binary指纹逐对validate，恢复原repeat顺序及独立身份，覆盖前保存paired manifest历史；V4注册恢复身份至补算排除集合，仅补不足。正式稳定性复核使用当前参与场计算的模型前缀（最终5次），不以缓存20的median替代真实5次资格。37项模型/追加/调度/V4回归通过。保存实际数据重放恢复CTLA4 203额外模型、PDL1 198额外模型，原保存身份缺失0；源文件不变，native_model_resume_replay.json记录证明。原源模型manifest在v8阶段覆写前另存4份历史，不改正在运行冻结阶段；下一阶段切换后使用新恢复机制。该工程恢复不增加科学资格。

v9冻结于独立detached HEAD 9f01e12检出surface_perturbation_verification_20261002_v9/ecd_epitope_scanner，257文件原源码overlay；122项surface/状态/调度/恢复/V4隔离回归全部通过（12.87秒）。不是clean HEAD声明。独立检出使A继续时B仍按旧v8快照完成，不并发覆盖其源码；数值七模块未改，同案例账本保持累计。settings_v9和lane_A_v9记录恢复，B待其当前阶段完成后切换。

### 2026-10-03 同批次续跑核验

CTLA4当前v9校准和生产流程已正常结束：校准物理ZIP 5034项/526对、生产4276项/530对、两层完整包均PASS。科学状态仍为CALIBRATION_EXPANSION_REQUIRED，23个正式候选的局部证据不能关闭；正式构建数0，不能引用历史覆盖率或称四项覆盖完成。33个保存null协议的逐通道诊断记录于CTLA4_calibration_resolution_diagnostic.json；部分required-scope Q95到20次仍不稳定，少数数值界未知，按合同保留未知。

PD1已完成v9首批校准和生产导出（生产3887项/470对，外层67项PASS），进入剩余校准。核验无存活本批worker后，使用同一v9隔离检出恢复PD1及PDL1→B7-H6两队列，不重开账本或退款；PD1未关闭的活动账本由RunBudget按原规则计费。EGFR不进入队列。日志lane_A_v9_resume_20261003.log与lane_B_v9_20261003.log，当前状态current_acceptance.json；后续状态覆盖此起始快照。

known_defect_verification.json记录已知修复的源码hash、关联回归和未解析科学边界。122项隔离回归不能替代四病例正式验收；实际破洞、null不稳定及准备漂移不会被标作修复后有效物理结果。

### 2026-10-03 关机根因与全包回归

System Kernel-General事件13确认00:44:25关机，事件12确认08:32:44启动；不是算法数小时求解或卡死。session_interruption_20261003.json保存核验。PD1/PDL1活动校准账本按冻结合同保守计费到35782/41308秒并关闭，不退款或新开预算。两例生产复核已实际完成：PD1物理包3993项/470对、外层67项PASS；PDL1物理包4150项/487对、外层66项PASS。状态均CALIBRATION_EXPANSION_REQUIRED，正式panel未验收。B7-H6由独立Hidden后台队列接续，进程归属durable_lane_B_launcher.json；power_guard不能阻止用户手动关机。

完整包回归发现旧默认profile测试仍固定6。当前已版本化默认8，历史V1仍6；只修改原tracked test_assay_profile，加入历史6保持断言，不改profile或科学参数。v10隔离257文件全包393通过/6失败/11跳过（748.05秒）。其中4个FlowKit API在原工作区通过，隔离环境缺少.tool_envs/flowkit；只用junction复用原工具环境后4项隔离通过（7.55秒），API未改。余2项在原工作区也失败，属于用户已有未跟踪FlowJo XLSX空格None/inlineStr断言，原文件保留；不声称全包全通过。详见full_package_regression_v10.json。表面专项此前122通过，数值七模块未改。

### 正常退出计时与阶段边界续跑

原生Windows最小重放确认旧guard正常退出会漏记启动开销：初始账本7秒，实际执行约0.517秒，结束仍为7秒。`settle_completed_phase`现将子账本与整个启动区间的累计下界取较大者，不能把同一区间重复相加；结束超过限时也拒绝验收。修复重放约0.528秒、账本7.5115秒且关闭。真实证据为本批`guard_clock_replay/old_v10`、`fixed_v11`，只证明计时修复，不是科学病例。

协调器新增首病例`--start-phase production`，只允许校准账本已关闭、summary与manifest/ZIP均属于当前尝试时使用。部分null可以交给engine逐项复核，不能绕过局部范围、匹配协议、真实重复、稳定性或未知候选的资格门槛。此入口避免已完成校准交付因guard源码更新再跑一轮；后续确有扩池任务仍按原闭环追加。

v12隔离overlay共257文件，位于`surface_perturbation_verification_20261003_v11/ecd_epitope_scanner`；父目录原venv验证138项专项回归通过（13.02秒）。相对v10只改变guard、协调器及对应两份测试，数值模块不变。B7-H6当前v10校准仍按其冻结源码完成；source guard在下一阶段阻止混代，watcher只在账本关闭后切至v12生产，预算不重置。阶段切换记录为`guard_fix_v12_handoff.json`；截至本次记录仍在校准，不把预定交接当作已经成功。

只读`diagnose_calibration.py`对四例已保存协议按每份required_mask与canonical面积逐通道核验，输出`four_case_calibration_resolution_diagnostic.json`；数组内容hash不符时跳过，避免并发写入形成误诊。该诊断包括历史scope协议，不按独立生物样本统计、不用算术均值充当最终面板面积，也不改变任何阈值/资格/账本。当前失败不只由是否执行求解决定：部分Q95在20次仍不稳定，部分投影或数值界未知；真实破洞与非流形mesh保持unknown。B7-H6仍在追加，最终以关闭的阶段证据为准。

扩池新增失败mesh见`B7-H6_mesh_failure_audit_expansion.json`：仍是3条边界边或2条非流形边，闭合流形边的方向冲突均0。未完成的工具cache没有complete身份，重试会重新生成OFF；其中此前一份OFF内容hash已改变，不能把旧hash当作现文件的身份。`failed_mesh_snapshot_manifest.json`公开该差异，`failed_mesh_snapshots/<sha256>.off`另存本次读取的原始字节。此备份不修几何、不将失败cache标complete；已成功的物理cache仍按完整hash审核。

11:00阶段交接发现本次一次性后台launcher的日志包装器把`subprocess.Popen`替换成函数；Windows asyncio需要继承该类，导入即TypeError，生产子进程尚未启动。改为Popen子类，另防止watcher读取上次失败metadata的启动竞争；新解释器同顺序导入通过，唯一新lane B13已真实进入FETCHING/ANALYZING。失败metadata保留，本批cal7413.594秒/prod358.704秒在恢复前核对且不清零。修复只涉及campaign helper，不更改冻结物理源码；证据`launcher_type_repair_20261003.json`。

v12第一次B7-H6生产真实完成：3302项/398对物理包与67项外层包均PASS，状态仍CALIBRATION_EXPANSION_REQUIRED；隐藏接触审计明确UNAVAILABLE，无合格真实面板不运行标签比较。engine摘要累计712.782秒，外层正常退出账本772.1916秒，补计约59.41秒启动区间，没有把engine报告直接当作完整wall-time。随后进入下一批校准。

当前四例保存生产模型池逐组合替换政策核查：CTLA4 33候选/530对、PD1 36/470、PDL1 37/487、B7-H6 37/398均未发现新C/P/G或原C/P/G修改；manifest hash单列于`four_case_model_pool_policy_audit.json`。这是当前已建模型的政策证据，不推断未建候选或正式面板资格。B7-H6后续补模型时应按最终manifest再核验。

### 方向资格与局部构建体交付修复（v13）

四例原队列现均结束，完整面板仍为CALIBRATION_EXPANSION_REQUIRED；B7-H6最终生产物理包4454项/551对、外层67项均PASS。包一致不等于科学收敛。已保存严格候选子池包含CTLA4四项、PD1三项、PDL1两项、B7-H6一项；旧页面只展示正式面板零项，不能据此推断所有个体构建失败。

实际发现`calibration_scope.local_scope`同时要求向量三个坐标分量与联合模长收敛，造成坐标依赖：同一六块null和五次变化只做45度刚体旋转，旧资格从false变true，联合分类却不变。数值反例保存在同批`vector_scope_rigid_rotation_defect_before.json`。v13以标量通道和联合向量通道作为物理资格，坐标分量仍完整保留为诊断；要求每个声明物理通道在全部实际局部支持上有效，数值未知或联合模长未收敛仍阻断，不修改Q95/Q99或4/5方向一致合同。

新增`surface_perturbation/local_constructs.py`独立审核正式替换资格、当前政策、五个真实不同配对身份、逐对原始模型验证、稳定性与同五模型DDG一致、局部闭合及实际变化mask；输出VERIFIED_LOCAL_CONSTRUCT，与完整panel状态分开。JSON/CSV、物理ZIP、外层manifest及页面均接入；突变用UniProt编号，另列author编号和FoldX规范化编号。失败尝试不保留旧局部资格作为当前结果，历史导出另存attempt_snapshots。

隔离v13是detached HEAD加259文件明确overlay；父目录原venv专项149通过（16.93秒）、API/UI七项通过（3.07秒），浏览器新表格及原独立预算检查PASS，所有请求拦截不创建真实任务。七份物理数值源码与v12逐hash一致，证据`v13_regression_and_source_audit.json`。同批两队列正在生产重核CTLA4→PD1、PDL1→B7-H6，保留原输入/AF3/confidence与累计账本，不重新校准已经耗尽的病例、不重置预算。当前运行不是四例已经完成；以`local_v13_*_results.json`及最终独立交叉核验为准。

本轮v13四例生产现在均已结束，11个局部候选独立核验PASS：CTLA4 T82W/Q117R/S49T/A77Y；PD1 R112M/N/Q（同一位置的备选）；PDL1 L27V/E45D/K46N；B7-H6 Q132W+M31A。新增L27V来自修正方向资格，不是放宽阈值或新增虚构模型。每例物理ZIP和外层包PASS；完整原始模型池政策另由`four_case_model_pool_policy_audit_v13.json`复核。`verified_local_constructs_v13.csv/fasta`是独立局部交付，FASTA按保存ECD逐位WT核对，原未解析结构部分的可靠性限制不因此改变。

真实浏览器另外暴露旧模拟fixture没有覆盖的接入遗漏：pipeline/API将完整扰动对象存于`design_summary.surface_perturbation`，旧客户端只查询`joint_v4_ml.surface_perturbation`，真实页面不显示表格。v14只修改app.js读取真实字段优先及browser fixture匹配真实schema，隔离259文件中其余全部与v13一致；API/UI七项复测通过（3.65秒），四例原生analysis浏览器均PASS、无错误和API写入，证据`native_browser_v13.json`（物理数据v13，客户端v14）。无需重算已核验物理证据。

最终状态见同批`final_four_case_acceptance_v13.md/json`、`four_case_local_construct_verification_v13.json`及`current_acceptance.json`；原current另存before_v13，旧final报告不覆盖。四例生产预算均关闭且未超7200秒，全局仍CALIBRATION_EXPANSION_REQUIRED；不能把局部合格11项称为完整四例panel或约8项覆盖绝大多数表面。最终APBS、正式K曲线、独立全物理路线比较、参数冻结和功能验证仍未完成。
