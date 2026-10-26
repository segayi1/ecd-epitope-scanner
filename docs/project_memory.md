# 项目记忆：有效决策、故障证据与维护约束

## 2026-10-26 运行剩余时间与FoldX / AF3状态核查

新增job_timing.py与JobStore.time_estimate，详情API/SSE发布同owner/同输入协议的历史阶段耗时区间；当前execution trace识别FoldX首轮/追加等真实阶段，不用百分比倒推。至少3个不同执行参考，历史缓存60秒，排队/等待/结束不伪报倒计时，样本不足/超范围明确未知。恢复取最新execution_id，发布重放副本不重复计样本；既有数据库无迁移、科学结果不改。运行条中英/手机及新离线完成态已覆盖。

ETA/API/生命周期36项通过；FoldX诊断/执行身份与AF3轻量接入58项通过、1Windows权限跳过/1完整集成未选。Chrome1440/390px中英四场景模拟进度渲染PASS，零HTTP/pageerror；实际只读历史AF3 binder5例、structure binder6例可估，部分抗原组仅1例不给数字。证据test_output/job_eta_20261026/，详见[估时及状态审计](job_eta_and_integration_audit_20261026.md)。

核对0 RUNNING/QUEUED、WSL /bin/true预检、旧PID12876及完整启动时间后安全隐藏重载。当前服务PID18924，health、time_estimate响应字段与HTML新字段均PASS，service_reload.json终态PASS；前面统一工作台和逐区规则也由新进程加载。未启动新模型或新增自动化，未跑全仓。

用户询问的FoldX大波动：实际准入取五槽median，采样隔离与跨门/逐槽风险显示已实现；并未新增方差硬拒绝，median通过且部分槽超门的组合仍可能保留资格并标采样需复核。AF3多模型预处理已接入共用流程且既有完整等价收据PASS；同目录普通Server模型的sample_id仍以父目录fallback，遗留未修，不能称来源身份所有布局均正确。本轮没有改变上述科学门或重跑完整AF3模型链。

## 2026-10-25 当前主合同：4–12分区，各原则独立给出合格推荐

用户明确4–12指显示分区，各区分别保留不同设计原则下通过门槛的推荐，不以全局4–8构建体数量截断。常规长度分区函数现对所有链长应用4–12软偏好，每区目标仍按长度7/8/12；连通/直径/面积硬门优先。全局精简面板仅作附加选择，预算不再占据逐区结果顶部。物理预览与80组合/300槽评分合同保留。

最终原则备选不再回退未通过项；缺口保存最近候选ID及失败/缺失门，CSV/JSON同步。主代表的近门槛诊断和冻结旧备选放入折叠审查详情。评分前/追加确认的engine两处调用显式confirmation_planning，仍保留未确认原则候选优先级。展示描述补充先复制单点记录，避免改源科学单点元数据。

18文件179项回归通过，最后CSV/显示专项20项通过（重叠不相加）；12区×全局4/8容量各原则完整与CTLA4源分数/资格/hash不变通过。125残基抗原7区、77残基binder7区新离线站点真实Chrome断网PASS，未通过卡折叠、精确目标、伙伴/表示/七步/单viewer/390px与零HTTP/pageerror。证据test_output/regional_contract_20261025/，详见[本次审计](regional_recommendation_contract_20261025.md)。未重算模型/重启服务/跑全仓；历史包不回写，评分不足的原则仍明示缺口。以下既有全局面板预算说明仅是附加面板历史，不能替代本条主要产品结果。

## 2026-10-25 Binder / 抗原共用工作台与长度分区

用户要求沿用修改后的binder前端统一普通模式，并明确每区残基数根据实际分析链长度而非binder标签分配。V4原本共用`display_sectoring_parameters(len(records))`，<120/120–150/>150 aa分别目标7/8/12；本次保留算法，修正前端及旧API测试的8起步限制为7。DiscoTope仍是主要评分差别，binder关先验并保留assay/膜方向未知、目标链和参考映射语义。

app统一完整结构来源/重试和第四步伙伴开关；renderer用既有getBundle建立精确分析残基组件，避免普通抗原同链ECD外片段露出，表示重建/复位沿用精确集合，离开第四步隐藏伙伴。在线缓存版本更新，新离线站点使用相同源码，历史包不回写。

父目录venv八文件72项回归通过；真实保存125残基抗原与77残基binder新离线站点Chrome断网均PASS：目标集合准确、伙伴/跨步隐藏、表示切换、单viewer、七步、390px及零HTTP/pageerror。收据test_output/unified_target_20261025/regression.xml与browser_acceptance.json；详见[本次审计](unified_target_workbench_20261025.md)。未跑全仓/重验在线服务、未重启或重算模型。早期相关binder受体可见性修复仍须保留。

## 2026-10-22 AF3多模型预处理交付与五例效果评估（完成）

实现已接入目录CLI/API、公共流程、在线/离线显示、原始输入归档和重放；只选原始构象，AF3语义/置信度和既有位点规则/80组合300槽合同保留。最终比较完全重算两边Surface V2几何，不使用早期临时缓存，源码与25模型输入SHA核对通过。

最终仅AIB10从model0选为model1，其余AIB203/AIB26/AIB56/AIB60仍model0。AIB10核心RMSD0.252Å、全目标0.616Å、pTM与ranking不变、平均pLDDT+0.0166分；接触Jaccard0.967，合格暴露位点67→66，Surface V2区域3→4，成员Jaccard0.953但分区ARI0.204，反映边界敏感。早期“2/5换模、AIB26m002”为把并列对照标签误当selected_model_id，已更正，不能继续引用。

当前专项92通过/1Windows symlink权限跳过（test_output/af3_final_focused_20261022.xml）；最终完整等价1通过（af3_final_equivalence_v4_20261022.xml），映射/可靠性/表面/显示区域/AA候选与两个导出包审核通过。修复np.mean与原流程statistics.mean导致的0.01分差异；同义密码子随机优化仅在测试冻结。保存PD1映射100%，288个残基confidence和PAE已映射；最大单模型2048目标全矩阵容量通过，未实测100个最大模型组合。

权威交付：test_output/af3_final_delivery_20261022/evaluation_zh.md与acceptance.json；最终原始比较test_output/af3_preprocessing_comparison_20261022_32de0d8a/。不回写历史运行；没有新AF3/SPURS/FoldX或正式突变推荐评分，未重启服务。浏览器file://受URL策略阻止，静态与归档检查完成；真实多模型ECD和实验准确度基准未完成。以下2026-10-21条目是实施历史，已由本条替代。

## 2026-10-21 AF3多模型预处理（实施与验收中，未交付）

用户批准仅在突变筛选前增加文件夹多模型输入、原始CIF/confidence身份验证、目标映射、质量近似池和来源组等权核心medoid；0.05Å集合内用合法pTM等置信度选一个原始构象，保留原文件与AF3语义。新增诊断不回写位点资格；workers/普通80组合与300槽合同不改。实现已覆盖core、API/CLI/JobRunner、pipeline、在线离线显示、输入归档、重放；原工作树和历史结果保留。

真实源验证暴露并修复原CIF省略冗余auth列与同目录冻结复制自身的问题。五例25模型前处理和源SHA不变已通过，收据test_output/af3_preprocessing_20261021/real_validation_v2/real_binder_audit.json；该收据早于最终基线/元数据修订，最终需要刷新。2048tokens、各4194304项PAE/contact_probs容量测试峰值281501696 bytes，见capacity_receipt.json。旧Binder/可靠性/PD1 22项及核心/接口专项通过，分组有重叠不求和；Windows无权限创建symlink的用例跳过。

尚需：真实文件名job_request/data与seed/上下文归档、原始基线范围/并列、解码前限额和缺少成员的逐项诊断最终核验；目录与选中单对完整公共流程等价测试尚无最终PASS。无新AF3/SPURS/FoldX调用、无服务重启；在线浏览器因本地服务不可用未验证。接续唯一状态tasks/af3_preprocessing_state.md，不得将此条当完成验收。

## 2026-10-18 AF3交付清单遗漏修复（不重算模型）

追加29组/58槽已于2026-10-18T02:11:22Z全部成功，四例派生卡片均PASS，但终验真实发现派生输出遗漏活动foldx_refinement清单：保存评分重放没有复制子目录，各例60条原始配对因此未进入最终目录。源结构/原模型/新增评分均完整，不能将此报错称为模型失败，也不能放宽集合核验。finalize_af3_completion.copy_refinement_evidence现在明确复制五份阶段JSON，绑定analysis.files与run_manifest.files后重建归档；根补齐清单保持原样，并要求最终真实配对身份集合严格等于源加58补充。首次失败收据、派生旧ZIP与元数据备份保留在before_refinement_manifest_repair。

实际终验还暴露validate_pair逐行重复解析同批所有Dif表的平方级I/O。仅在验收器内使用生命周期内解析复用，每份表解析前记录hash、退出前重新校验所有表hash；原始配对的结构字节、坐标身份、突变、槽、DDG匹配检查保持。新增回归验证相同结果只解析一次、表内容变化拒绝通过、离开作用域恢复解析器；导出/独立几何/发布相关8项通过。终验与真实五站ZIP网页检查仍以final_acceptance_execution.json终态为准。

## 2026-10-17 AF3最终验收补齐（独立追加58槽，接续中）

五例已COMPLETE；复查确认4个FAIL都来自真实展示卡五次重复不足，并非仍有80/300溢出。AIB10/AIB26/AIB56原验收把区域主代表缺五槽混称预算失败，现区分源协议不变量与派生补齐后的真实重复证据。29个组合（13/1/9/6）各补槽3/4，共58槽；`supplement_displayed_recommendations.execute_job`允许显式独立预算，默认历史20/43保持，原1500活动槽不扩大改写。任务下`run_af3_completion.py`冻结计划并一次执行，`finalize_af3_completion.py`复用已有重导出器；原目录保留，新增派生COMPLETE记录由无模型发布器注册。

独立核验必须同时读取root首轮240槽与foldx_refinement60槽，排除prior_attempts/paired_manifest_history；只读root清单会错误断言缺60槽。最终验证器核对原协议80/300、追加账本、每条实际配对、完整候选池身份/SPURS不变、逐槽FoldX值和中位数，再审全卡片。五例几何独立重算全部PASS（35区/277节点）；29项相关回归通过及专项12项通过有重叠不求和。AIB60实际file网页/独立报告/七步/结构/草稿JSON与AA/CDS下载/零HTTP通过，其余四站需最终ZIP形成后实测。

`finish_af3_acceptance.py`已一次性后台等待追加进程43536退出（launcher48476），随后核验并发布四例派生记录，五ZIP解压到unpacked_offline_sites后复用浏览器脚本实测。所有原任务、失败产物及早期验收收据保留；没有sleep轮询/自动化。终态、报告、UI收据分别见final_acceptance_execution.json、final_acceptance_zh.md、final_offline_browser_audit/acceptance.json；只有终态PASS才是完整交付。运行服务核对PID/启动时间且任务全空闲后才重载，重载失败单列，不改科学结果。

## 2026-10-17 AF3 FoldX首轮预算溢出修复（等待原队列结束后补跑）

AF3 AIB10 `87f726c6fce8`与AIB26 `2788c19f05c8`在首轮FoldX后进入refinement时失败，分别已有81组/243配对、86组/258配对。AIB10按保存区划核对为61局部+20全局。根因是`select_v4_foldx_first_pass`给逐区优先候选分配局部名额后，探索份额仅扣除已入探索的优先候选，没有同时扣除所有已占局部容量。`candidate_search.py`现将探索限制为剩余局部容量和剩余探索份额的较小值；`engine.py`在启动首轮FoldX前另行检查80组上限，避免昂贵计算后才报错。45优先局部+16探索局部+20全局的回归覆盖本次81组复发。

普通协议/恢复筛选与历史/任务证据历史/推荐补算18项回归通过（5.37秒）。`workspace_data/af3_binder_batch_20261017/budget_repair_evidence.json`绑定失败实际消耗、源码hash与测试记录。原进程PID48180已加载旧模块，磁盘修复不热改其活跃任务；后续失败使用同任务ID在新进程加载修复代码，独立保存失败数据库记录及顶层产物至`prior_attempts`，严格复用身份及槽合同一致的原始配对。已超额的旧尝试单独留账，不能把新运行80/300口径当累计消耗。

一次性后台恢复脚本`recover_failed_af3.py`用OS进程退出等待与文件锁，原队列终态后只恢复FAILED，成功任务不重跑；拒绝重复恢复和存在活跃任务时启动。无sleep轮询、无新增自动化；终态调用独立五例验收器，计算成功与验收状态分别保留。接续已启动（launcher16256/实际43456），收据`recovery_receipt.json`为`WAITING_FOR_ORIGINAL_QUEUE`；最终19项回归4.55秒通过，新增活跃任务禁止服务重载的断言。在线服务仍加载旧Python模块；接续只在五例验收通过且任务全空闲后核对原服务PID27912及启动时间再重载，不能提前声称修复已在线生效。重载只记录新服务启动，健康检查待用户回来后单次核对。

## 2026-10-17 AF3 binder 五例全流程（后台运行）

用户在文章来源六例完成后授权继续处理此前暂缓的五个 AF3 binder。使用其新提供的 `20261016_raw data_binder` 原始AF3文件，不调用AF3推理：每例取 model_0 CIF、同索引 `full_data_0.json`（含 atom_plddts/PAE）及 summary_confidences_0；model_0 在五模型排名分数中均最高或并列最高。AIB10/AIB203/AIB26/AIB56/AIB60 五条A链序列分别为128/64/108/123/66 aa，均与既有预检序列完全一致；单链与 confidence/token/PAE 映射预检PASS。源哈希、五个种子hash和预算规则写入 `workspace_data/af3_binder_batch_20261017/batch_manifest.json`。

五任务以 `input_type=af3_upload`、`analysis_target=binder`、chain A、strict V4运行；目标区大小按长度用7（<120 aa）或8（120–150 aa），显示请求上限12、4–12为软范围。沿用普通V4双/三点、80唯一FoldX组合（60 local/20 global）、80×3首轮及最多30组补两次（总300配对槽）；DiscoTope对binder关闭。通过单独后台JobRunner串行提交，不重启在线服务；首个AIB10 `87f726c6fce8` RUNNING，其余四个QUEUED。完成后需按清单核对结果和每例离线包；旧六例与失败/历史目录保留。

## 2026-10-17 六例推荐补算、卡片完整性与结构总览（交付完成）

20个已展示但FoldX槽不足的唯一组合追加43槽后，六例在线结果以新COMPLETE任务发布，原运行目录、源哈希及已保存值均保留。新任务ID与来源映射记录在`workspace_data/recommendation_completion_20261017/delivery_manifest.json`。FoldX新增用途：19个组合补第4/5槽，AIB201单个无初始槽组合完成0–4五槽；总共43次WT/Mut配对槽，不增加其他候选计算。

推荐卡完整性硬门现在遍历primary、alternatives及全部正式本区展开候选，按每次显示实例记录父区/角色；每例52–63实例均PASS，三输入文件hash绑定。His的质子化在原始AA属性中未知；新增展示专用独立微态枚举，His按中性(电荷0/供体1/受体1)和质子化(1/2/0)分别计算场景、轴/逐位范围和近邻方向分类，不指定pH。原始排序描述符null及资格不改。`web/recommendations.js`将有效范围投影到对应计算轴/逐位单元格，避免范围已可追溯而旁边仍写“未计算”。

遇到两类导出复发点并最小修复：`reexport_recommendation_supplements.py`按ledger仓库相对路径读取时固定以`__file__`推导根目录；补算候选的display-only字段不得传入`V4Candidate`严格schema；发布候选身份校验按notation/位点/替换字母，允许新增逐位证据但逐单核对SPURS、FoldX槽及原值。首次失败/中间目录保留。

结构选择保留同一viewer实例及现有representation，聚焦相机600ms过渡（系统减少动态效果时为0）；多选摘要标注位点重合/替换冲突，可清空且按任务隔离。真实在线AIB239页检查His范围、两组合冲突、清空、跨任务隔离及单canvas定位；六个完整ZIP 71成员且`audit_run_bundle`均PASS。在线与离线共用显示代码并重建归档哈希。相关后端41项、前端14项、生命周期6项及本轮渲染2项回归通过。file://断网交互未在此轮重新实测；真实缺少来源记录的pLDDT/表位先验仍明确呈现，不填造。

## 2026-10-17 P0–P3终态复审（计算结束，验收缺口明确保留）

本轮用户询问运行是否结束并要求再审查。一次性终态核对：六例新任务均COMPLETE、源/配对/预算/导出验收PASS；42区6–8个表面残基、2,664组合SPURS、480 FoldX、180五槽、177正式、48主面板；每例80/300，local60/global20，强探索五例16、AIB239为13。旧283正式与新177不可直接归因为门槛，五槽确认池已由每例60变30。外部runner 2026-10-17T11:10:23Z FINISHED，48尝试、44 SPURS、47 FoldX、235有效配对、240已预占槽；4 SPURS失败及1启动无结果保留不重试。最终汇总已重建，数据预审脚本仍只负责源hash/映射/符号，不能把旧0槽预审收据当终态验收。

最终报告`workspace_data/p0_p3_review_20261017/final_review_zh.md`及review_acceptance.json（35项断言）核对指标/账本/身份hash；独立只读六例与基准审查一致。全候选high_risk是SASA_METHOD_SETTINGS_UNVERIFIED旗标传播，应补方法证据且分拆语义，不能称全部稳定性失败。外部双模型数值线内有1SAKA +8.4、1MGRA +2.7源实验标签，但实验条件缺失/结构上下文需核对、训练同源未知，不据此自动调阈值。仅AIB239有本轮结构化file页PASS；前两例保留早期浏览器产物，后三例file策略受阻，不绕行重试。全部ZIP静态一致性通过仍不代表5/6完整离线交互验收。

报告图表错误仅接受audit.status=PASS，实际报告生产者是AUDIT_COMPLETE，导致全终态图仍写PARTIAL；write_p0_p3_comparison_figures.py最小修正为联合终态与两种审计成功枚举，测试保留RUNNING+旧PASS为pending。本轮无模型调用、无服务重启、无自动化恢复；用户明确后台运行与成本约束继续有效。

## 2026-10-16 P0–P3新协议（代码实施及完整新批次运行中，未交付）

任务证据根目录`workspace_data/p0_p3_upgrade_20261016/`。用户完整批准：普通新V4仅双/三点；<120aa目标7/6–8，120–150保留8/6–10；SPURS单点/组合≤1（≤.5优先）、FoldX≤2，差/符号仅诊断。80 unique=60local/20global软资源池，可回流；最多16local强扰动窗口(.5,1]；3初始+最多30组合各追加2，≤300尝试槽，正式候选需5唯一槽。特殊物理合同显式恢复原阈值及60组合，历史不回写。

首轮部分子任务中断后主任务接续：真实SPURS GPU6个双/三点编码复用数值差0；严格cache实际权重/源码/环境/参数绑定、损坏保留拒绝复用。序贯调度、同证据本区候选、近邻/精确重复、同池界面对照、模型差/重复审计、源码ZIP/阶段事件/owner诊断、完整同源经典脚本离线网站均已接入。203相关回归PASS，之后专项回归有重叠不相加。真实禁模型AIB234冒烟7区，79文件/ZIP哈希及唯一成员通过。预览无评分不得贴机制policy已修复；P0–P3相关pipeline及专项回归继续验证。在线历史页结构与rSASA-only说明实际核验。当前两例最终任务ZIP先核验CRC/成员路径，再真实解压并由Chrome以`file://`断网打开；七步、7区、内嵌结构及`MolstarView.isReady()`均通过，本地未评分突变草稿JSON/氨基酸与CDS FASTA可下载，且无HTTP请求或浏览器异常。离线主站和独立报告共用分析投影的回归也通过；Mol* wrapper不公开先前假设的layout `isBusy`属性。其余四例待任务结束后同法验收。

新任务唯一清单`rerun_manifest.json`：f606506c8e4a、5ac94808fb61、55ac1dce3c04、f53da1c7dac0、d4613b14ee6b、9b232e35663e。六个源结构hash及A链身份沿原manifest；AF3五例暂缓。禁止重复提交、重启有活跃任务的服务或改门。恢复不移原始配对模型路径，先复制元数据及顶层产物到prior_attempts，再追加执行UUID与新原始目录。AIB234/AIB226/AIB239已通过独立身份/预算/配对/导出审核及离线网站验收；较早FAIL核验日志早于最终PASS收据，原文保存在`*_verification_initial.log`，当前日志是复核后的PASS。面板重复敏感性审计曾错误读取根层门槛；实际`run_manifest`将SPURS门槛放在`spurs_multi`，已修复读取和阈值路径，三例均完成3/5槽选择重放（30个五槽候选池、80个首轮候选，五槽面板精确复现保存面板）。三份80成员ZIP的完整成员/磁盘/manifest/CRC审计通过；通用审计器现在容许安全多层离线资产路径并拒绝路径逃逸。AIB239的根离线页和独立报告页均经断网浏览器核验，比较科学显示payload时排除在线/离线路径、文件hash及执行elapsed字段，其他结果字段一致；源码回归也直接验证两种导出共用同一投影。外部runner新增OS文件锁，防止并行恢复重复预留槽。AIB193实际FoldX子进程及输出活动、数据库阶段标签滞后；AIB201/AIB15仍排队。外部预审48变体（8单/20双/20三），训练参考簇重叠明显、单点条件有缺口，不作为通用校准。首次10条SPURS试点完成、5条FoldX尝试预占25/240槽，4条结果与20个配对证明有效，1条启动无结果且不重试；旧runner PID不存在且无终态收据，状态已标记`INTERRUPTED_NO_TERMINAL_RECEIPT`并保留原manifest快照。六例终态后续跑其余43条（215槽），总额不超240。


2026-10-17增量审计：安全的多层ZIP资产路径（如`assets/vendor/viewer.js`）此前被通用bundle检查器误拒，已支持安全相对路径并继续拒绝绝对路径、盘符、反斜线及`..`逃逸。AIB234/AIB226/AIB239三份80成员包经磁盘/manifest/ZIP SHA与CRC检查通过；runner锁、3/5槽敏感性读取及离线报告双页浏览器回归已补。

## 2026-10-16 Binder 六例预算与选择偏倚（只读报告完成）

完整报告与可复核CSV/图/脚本位于`workspace_data/binder_run_budget_review_20261016/`。本轮未改变生产代码、服务或历史运行，模型调用0；六例9,348单点→3,669合格替换→2,699组合SPURS→360 FoldX×5→283全局formal，48主面板；当前36区/306节点，平均8.5个表面残基。FoldX elapsed总12,807.42秒占记录模型阶段和94.1%，不能叫CPU核时；single/multi是批次计算，记录数不等于独立模型启动数。

修正区级缺口解释：7个NEAR_THRESHOLD机制primary区实际有14个全局formal且位点全落当前区的双点（2/4/3/2/1/1/1），10个已列PASSED alternatives、4个未展示。它们缺候选级regional policy而未参加机制primary竞争，并非全区无合格组合；当前历史primary仍待审，不直接改标签。可先在独立目录复核当前区ownership/几何/完整资格与impact/敏感性，无需新模型；三点额外重算当前区增益。此前13/9合计错误已由独立断言纠正为14/10/4。

预算层次：区域路径全菜单枚举后每区20双/10三（最多12双亲本），通用路径每位点最多3选项；合并物理身份去重。每蛋白FoldX60 unique组合、5重复、内部最多8并发，不是每区60；新36区专属131/360，跨区229/360。四点360组只作SPURS探索，不能FoldX正式确认。有效主面板min4/target7/max8，实际COVERAGE_FIRST_GREEDY，各例6双+2三；requested_maximum误存显示区请求，不能当有效构建上限。

偏倚信号：净疏水增加41/48主面板；界面突变槽保存池28.5%→主面板16.7%，混合覆盖/风险/菜单目标，未做单权重因果对照。77模型分歧中75各自门通过且FoldX更负；实际FoldX取5槽中位数，前3槽比较仅4/360联合门翻转、主面板0/48，存在尾部变异，不支持直接削减5次/放宽分歧门。`indistinguishable`是RMS<0.1近邻数，非向量完全相同；仅AIB201一对精确重复。部分UI仍称sidechain≥20联合表面门，但本批实际总rSASA-only、侧链描述用途；P0文案/字段修复尚未实施。菜单与图7个位点差异来自几何owner要求，缺原子级记录不能进一步归因。

`report_acceptance.json` PASS：74主输入+96原始FoldX文件重验hash不变，独立合计/14候选身份/本地链接/四图嵌入通过，两路代码/偏倚复核。四图已看，HTML浏览器file协议预览被策略阻断未绕过；静态HTML已查。文章来源仍不等于实验结构，训练集同源/外部稳定性校准/实验结合未验；报告是审查与优化提案，不是生产优化或实验完成。

## 2026-10-16 小蛋白分区尺度（实现及六例派生验收完成）

用户明确每区约8个表面残基、总区数4–12软范围；按分析链≤150aa用8/6–10偏好，长链12/10–15，显式覆盖优先。config.display_sectoring_parameters由engine传原始未合并overrides，避免默认12被误当显式；UI留空自动、API支持8。sectors保持4残基/50Å²/连通/唯一覆盖及30/22硬门，4–12仅软优化；物理预览冻结旧策略不改。

服务重载规则：不在有活跃任务时重启。六例全部终态后确认无QUEUED/RUNNING，按仓库launcher要求从包目录父目录用项目venv启动；/api/health、既有派生job读取和小链8/6–10/4–12、长链12/10–15默认检查通过，收据`service_reload_acceptance.json`。

文章来源批次最终6/6原始验收PASS，每例60 FoldX候选×5次与300对WT/Mut结构审计通过；六个小蛋白派生分区也通过结构、来源、评分保持及导出核验，共36区29 PASSED/7 NEAR_THRESHOLD，无新增模型调用。详表以`workspace_data/binder_input_preflight_20261016/article_batch_summary_20261016.md`、`batch_acceptance.json`、`small_protein_delivery_acceptance.json`为准。

## 2026-10-16 Binder 输入适配（文章来源六例批次完成）

上传/API/CLI新增analysis_target=binder（默认antigen不变）；单链自动、多链显式binder链，坐标派生序列保留作者链/编号与插入码。binder不套用固定抗原assay、膜方向/膜锚保持未知，不读取未知B-factor为pLDDT；DiscoTope抗体表位先验关闭，原SASA/ProCaliper/SPURS/FoldX与稳定性/空间门沿用。独立V4来源、序列及验证文案保留binder身份。真实浏览器发现旧app.js/i18n缓存不匹配新HTML，资源版本参数与单一动态说明修复后入口正常，窄屏表单无横向溢出。

用户确认本批6个文章来源双链复合物全流程，5个AF3单链暂缓；用户最终清单为`workspace_data/binder_input_preflight_20261016/article_batch_manifest.json`。六个原始运行与六个保存评分池分区派生全部完成，原始与派生验收细节见批次总结文件。11文件/17链映射、11个binder A链SASA/能量导出及两份禁模型预览通过，原始哈希保持。AIB203早期冒烟已停止且不计入本批。

## 2026-10-16 机制排序批判修订（实现与保存评分池验收完成）

逐区策略族保留regional_mechanism_pareto_v1，新排序显式mechanism_channels_v2：原始六轴+暴露Pareto、缺失掩码/冻结参考保留，电荷/亲疏水/氢键及体积/芳香两组内取最大效用，最小暴露单独成组，组间leximin；搜索多样性仍按六轴。增加去除芳香轴敏感性且重算前沿，单候选不标稳健。signed净变化/绝对和/代数反向比例、FWY数、WT Cα近邻与距离仅诊断，不能作净零惩罚、CB/盐桥或碰撞证明。SPURS单点REVIEW窗口审计保留其他排除及两种暴露定义，门槛不放宽，不新增8–18Å窗口；表位仍只在新逐区路线注释。

涉及regional_perturbation.py、engine.py、sector_recommendations.py、reporting.py、共享recommendations.js及显式新目录replay_v4_selection.py。重放重建单点资格与源CSV逐项相同，并保留描述符/侧链证据；这是当前算法重排旧评分池，不是复现旧排序或完整新候选模型链。CSV补齐resource_audit及实际sensitivity，旧报告无新ranking_version不套新解释。原始CTLA4/PD-L1芳香富集与近邻异向现象是审计信号，不能推断聚集或自中和。

58项相关回归；两组556/753固定候选新评分0，评分/资格/源哈希保持，各60文件与ZIP一致；在线冻结显示投影和实际离线桌面/手机8项通过。CTLA4分区2换A44R/S55F，PD-L1分区10换K178V/S176H仍EVIDENCE_LIMITED；原4个模型差异超限区仍保留。新完整成员窗口225/463中其余条件合格125/320，高侧链暴露33/120，不能与原mutable-menu窗口67/141混为同一分母。证据workspace_data/regional_mechanism_critique_20261016及[修订验收](regional_mechanism_critique_20261016.md)。服务空闲重启加载新代码，未重跑完整新模型链/全仓/实验或提交推送，历史产物不改。

## 2026-10-15 全站 UI 重设计（展示验收完成）

用户经 grill-me 确认七步向导、浅色工作台、笔记本优先、按步骤调整结构/右栏比例、序列默认紧凑摘要后批准实施。首页新建/上传/高级参数移入对话框，默认 WT 全结构；序列展开、选择/视角保留、键盘与分区/组合联动保持。共享 workbench.css 统一主站、实验、抗体和结构图册外围；新离线模板自包含，旧报告不重写。主代码 web/index.html、app.css、app.js、renderer.js、experiments.html/js、offline_template.html 和图册外围；分区报告生成器内联共享 CSS。

修复时不能把每行 QC 表单也套入通用折叠容器，须保留表格操作和预览输出的相邻关系；不能用 overflow-x:clip 隐藏手机溢出，应由表格内部滚动与正文换行解决。截图必须等待真实结构层级和 Molstar 布尔 isBusy=false，单纯存在 canvas 不足以证明已渲染。标题父级字号/字重不能传播到说明段落。

真实首页 be18bf6290fb 在 1366×768/1440×900/1920×1080/390×844 通过；1366 下结构画布约902×472、结果步骤779×472，切步骤相机目标/选择/实例保持，保存结果到冻结示例真实载入只销毁旧实例一次。四个附属在线页各三尺寸、实际 file:// 离线断网、完整实验页浏览器回归及18项相关pytest通过。实验页写操作用模拟响应；真实调用仅包含已核对的只读绘图/序列预览。结构图册内嵌脚本与 HEAD 哈希一致，冻结资产无diff。收据 workspace_data/ui_redesign_20261015/main_final_verified、pages_final_verified、offline、pytest.log、experiment_regression.log；详见 [UI 重设计验收](ui_redesign_20261015.md)。未重算科学模型、跑全仓测试或修改历史科学结果；不将 UI 验收称为实验验证。

## 2026-10-15 逐区机制帕累托推荐（实现与计算验收完成）

新常规默认regional_mechanism_pareto_v1：全量合格单点、原生侧链相对暴露、六项逐位带符号化学描述，组合按平均绝对变化、缺失掩码分层帕累托和冻结分布分组折中选工程代表；每区20双/10三、最多12双亲本扩展。DiscoTope仅在新逐区路线作注释，不参与排序/调度；固定分区不变性不等于整个旧全局流程无表位先验。共享物理身份只评分一次，原始完整序列/24Å跨度/双模型/高阶增益门及全局预算保留。主卡显示直接面积与本区12Å覆盖/未覆盖，详情9/12/15Å及敏感性，在线/离线/CSV/JSON共用证据。

核心代码joint_design_v4/regional_perturbation.py、engine.py、sector_recommendations.py、reporting.py及structure_analysis.py。Bio.PDB逐残基deepcopy会复制20模型NMR父级，改实体copy脱离parent，冻结1AH1五个完整侧链新旧数值一致。敏感性每种变体重算相应前沿；未评分preview单点不进入新机制组合；未知His化学轴保留缺失。历史方法、科学门、运行目录和物理预览合同保持。

基线6999657/文献d092c12先提交。真实最终CTLA4任务2fd8fa7fb244：7区6完整通过、1模型差异超限；PD-L1 be18bf6290fb：10区6完整通过、1稳定性通过但His轴未知、3模型差异超限。SPURS组合556/753；各FoldX60×5，CTLA4复用本轮真实中间任务300对，PD-L1新算300对。两组300配对均独立身份/映射/哈希验证PASS，各65清单文件与ZIP一致；最终科学analysis哈希未改。新机制27项、既有99项、最终展示20项及较早物理等187项通过，有重叠不相加。两任务在线真实API和file离线桌面/手机通过，冻结渲染模拟证据单独标注。

收据workspace_data/regional_mechanism_20261015；[实现与逐区面积验收](regional_mechanism_pareto_20261015.md)列明任务、组合、面积、测试、失败和模型复用来源。新实现保留工作区未推送；未跑全仓库测试/可选物理场预览/抗体结合实验。计算代表和几何邻域不能称亲和力下降最大或实测接触面积，未通过分区不能为凑齐输出而放行。

## 2026-10-12 完整科学评分链恢复（完成）

新完整任务首轮07bbf0dd36ca因venv Python继承Low完整性标签，Windows→WSL桥拒绝访问而失败；独立Linux CUDA推理证明SPURS/DiscoTope环境与权重仍有效。仅修复venv的python/pythonw为Medium，实际生产Python→WSL预检通过，新任务484899e022f3为Medium。旧Low网站进程已重启，启动器增加WSL实际执行预检。独立新任务COMPLETE：原生SASA/ProCaliper/SPURS single/DiscoTope各215，SPURS multi476，FoldX60×5=300fresh/0复用；10专属首选全通过、主面板8全有有限双模型。25相关回归、300配对身份/哈希、65manifest文件/1022展示引用/API/JSON/CSV/ZIP与正式网页单canvas/定位/详情全部通过。保守性仅1同源记录保留未知，IUPred无输入及可选物理预览未运行，四区大小9/9/9/8偏离软目标，不称实验完成；详见[完整恢复审计](full_runtime_restore_20261012.md)。

## 2026-10-12 网站恢复与FoldX默认60

系统重装后8765未监听，旧venv继承G:/anaconda而FastAPI缺失；当前系统TEMP不可写。旧venv保留，使用新Python3.12.10独立环境按lock装53依赖。启动器预检、后台隐藏进程、固定stdout/stderr、项目runtime_tmp及60秒等待已验收；正式历史PD-L1画布/10专属首选和健康API恢复，预算专项3项、相关站点58项回归通过。共享V4 FoldX默认/缺省回退为60、每组合5次/科学门不变，V3独立81不改，历史42次评分仍原样。当前工具环境WSL执行拒绝访问、ProCaliper导入25秒超时；未修改VHD权限、未重跑新科学评分，不能宣称推理链恢复。见docs/system_restore_budget60_20261012.md及workspace_data/system_restore_20261012。

## 2026-10-10至11 PD-L1画布及均衡分区（实施与真实验收完成）

用户批准实施。常规balanced_connected_v4目标12、偏好10–15、全图30/Cα22，硬4/50/连通/完整唯一及有限坐标保留；物理预览v3仍26/18，24组合跨度/42FoldX预算/科学门未改。固定PD-L1 103节点9区11/10/10/12/12/11/13/13/11，全合格；CTLA4/PD1允许软大小例外，局部搜索非全局最优。公共失败保留draft及具体违规，预算耗尽区分硬无方案与已合法但优化未完成。

画布单次初始化、旧实例一次销毁、过期回包令牌、单CIF请求、分阶段重试；tooltip失败保留结构。网页view=display，完整API默认保留，正式/离线共用精确文件路径引用投影；旧PD-L1首屏2,614,581bytes，完整科学审计/历史源保持。实际Chrome旧PD-L1/CTLA4桌面手机4项、离线2项及任务切换通过；手机隐式grid列4.5px溢出用minmax修复。83/91/44项分组回归（重叠不相加），完整工程PD-L1流水线4项及AF3模型/原confidence完整几何重放通过。

新真实任务a1b44daa4c7e SPURS multi458成功但FoldX RepairPDB Windows0xc0000409，FAILED保留；c244781d83a8运行配置键误放，在单点完成无子进程时停止并保留原因。正确配置30aa57c17c81修复成功后执行进程随会话中断而消失，未验收，保留FAILED。持久后台dcac41e103e8的BuildModel内部1800秒限制截断6组，36组180条完整评分；9对半成品PDB没有Dif评分，不能复用。无活动worker时停止失败后处理，保留FAILED。修正执行时限为至少3600/默认外层7200秒，记录原生超时/退出诊断；科学门、42组及5次不改，补充51项回归通过。

2026-10-11最终b4b51cda70f6 COMPLETE：ProCaliper215/SPURS单点215/DiscoTope215/SPURS组合458/FoldX42全部success；42组×5次210条完整配对（180身份验证复用+30新评分），原8个冻结源hash和映射/表面图不变。9个专属首选2通过7待审查，10个备选，8正式主构建体均双模型有限。基础4双/2三均达标；分区2/3各只有1可确认专属组，因此明确验证机会缺口。API/JSON/19行CSV/ZIP一致，65文件hash和485展示引用核对，首屏2,207,862bytes。真实正式/实际离线×1500/390最终4项PASS，无模拟路由/异常/溢出，正式单canvas/初始化/CIF、键盘展开及结构序列定位通过。root补分区3上下文暴露分类原因中文映射，2项渲染回归及重打包后导出/浏览器均复验，分析SHA不变。收据见workspace_data/pdl1_balanced_acceptance_20261010/recovery_final与workspace_data/pdl1_recovery_browser_20261011/final2。具体跨盘/文件名崩溃触发未独立证明；显示大小软目标及几何代理不代表实验验证。详见[pdl1\_canvas\_partition\_fix\_20261010.md](pdl1_canvas_partition_fix_20261010.md)。

## 2026-10-10 最新PD-L1画布与分区审查（仅计划）

最新真实b80f5cbd8784为Q9NZQ7/8AOM链A，103节点/14区，大小7/8/8/4/13/8/4/4/7/13/6/9/8/4；5个worker成功，SPURS组合480、FoldX42。旧17区/117节点是不同工程测试输入。当前4/50/26/18硬验收全部通过，但没有用户已确认的表面残基10–15大小目标。固定图冲突团13且合法13区见证证明现尺度至少13区，全区≥8已不可能；CA-only18也至少12区。30/22、9种子只读对照10区/min6，不能称全区10+或直接改生产默认。

真实Chrome冷加载桌面/390及表示/珠子/导航均正常，用户原载入失败未稳定复现；已实际确认每次并发两个Molstar实例（adoptResult/renderAll重复init、旧实例未dispose）与108MB接口/212MBanalysis/117MB离线。复用现有诊断外置助手可将内存紧凑analysis缩至29.77MB，但仍需轻量投影。下一步计划含单实例生命周期、错误阶段、诊断外置、目标大小与边界平衡、尺度选择及新归属/候选资源/正式科学验收；不是已经实施。源文件哈希保持、模型新增0、服务未重启。详见[pdl1\_canvas\_partition\_audit\_plan\_20261010.md](pdl1_canvas_partition_audit_plan_20261010.md)，证据在workspace_data/pdl1_canvas_audit_20261010和pdl1_partition_readonly_audit_20261010。

## 2026-10-09 模型软件版本识别

真实WSL安装包SPURS1.0.0，当前vendor提交70cbb99af5cdf1f6da392389aaf93fd245610efb，仅requirements.inference.txt有跟踪修改；DiscoTope模型系列3.0、Python包discotope3版本1.0，vendor没有Git信息。108个独立模型文件与6aef冻结哈希一致，原运行软件版本/代码提交未记录，当前补充身份不能替代历史源码冻结。worker新增实际包版本、源码指纹与vendor Git标识；禁止把Torch运行时或父工作台提交当模型工具版本。正式与离线展示模型系列/包版本分开并标注补充身份。新任务f2b213d9b0ad，评分未变、模型调用0、历史265条任务保留，API/JSON/CSV/ZIP一致。详见[版本审计](model_version_identification_20261009.md)。

## 2026-10-09 推荐详情与双模型准入补修

风险码按事实解释并去重；规则检查与模型数值分开，布尔规则不能显示为未计算。正式主面板强制SPURS multi与FoldX组合值均有限，缓存PASS/旧SPURS-only不得绕过；原8个主构建体其实本已双模型齐全。原始1AH1元数据为Solution NMR，pLDDT与晶体学分辨率不适用；局部实验精度未评估，坐标/映射核对单独披露。来源/序列简化为主表与折叠详情；暴露度准确标RSASA，先验未达高支持条件不等于缺失。基于6aef485cc7c1真实保存评分，新任务e824509801c2，9区7通过2待审查、8主构建体；本轮新模型调用0，历史264条任务保留。73项相关最终回归、94项物理/V4.1回归、完整pipeline4项及后续17项界面/风险补验通过（有重叠不相加）；API/JSON/CSV/ZIP哈希一致；正式与离线真实Chrome1500/390像素通过，无模拟接口、页面错误或溢出，完整风险解释/规则事实/具体审查理由已核对。详见[修复记录](recommendation_detail_polish_20261009.md)。

## 2026-10-09 专属资源与解释面板

按实际成员识别专属/跨区，4双点/2三点独立专属组及90%搜索触发、额外仅基础配额扩张，FoldX42先逐区两个专属机会；不足保留缺口。专属首选+风险/先验/跨区备选独立资格，完整通过/失败门槛、位点、互补、来源和序列，正式与离线共用渲染，未知不补零。jobs进度不提前发布COMPLETE，结果与终态原子提交。资源审计整数键导致返回与磁盘不一致已修复，168项覆盖回归与完整pipeline4项补验通过；末次45项细节补验通过。

真实评分任务6aef485cc7c1复用冻结1AH1坐标，输入CIF/规范PDB/映射相同、71节点211边及面积相同；SPURS346组合/FoldX42组合真实成功。展示最终任务8acec7ea445b，9专属首选7通过2待审查、8主构建体，保留真实源任务与原262条记录。分区7收益0.011079<0.05，分区8 SPURS0.802133>0.5且缺FoldX；不声称已证实不稳定。分区2/3/8基础资源不足及8无可进入FoldX专属组如实审计，不能为凑齐放行。详见[验收记录](exclusive_sector_acceptance_20261009.md)及workspace_data/ctla4_exclusive_acceptance_20261009；正式与离线真实Chrome在1500/390像素均通过，无模拟API、无页面异常或溢出；序列栏内部编号导致窄屏撑宽已移除正文编号并添加换行约束，重打包后API下载/ZIP/清单哈希及冻结评分源一致。浏览器收据见workspace_data/ctla4_exclusive_browser_20261009_final3与ctla4_exclusive_offline_browser_20261009_final。

更新日期：2026-10-08。用于后续开发查阅；历史运行是证据，不是当前算法规范。
  
2026-10-08 CTLA4六构建体结构查看页见 web/ctla4_six_mutants.html，冻结结构资产和清单见 web/data/ctla4_six_mutants/，由 tooling/export_ctla4_six_mutants_web.py 从原研究逐条配对清单重建。核验覆盖30个production WT/Mut槽位：PDB SHA-256匹配，链A序列差异恰为两个设计突变；折叠模型局部编号按完整序列比对映射至P16410 37–161。各候选WT唯一字节数为5/5、5/5、3/5、5/5、1/5、5/5；报告不得把逻辑槽位称作独立骨架。该页面是计算结构展示，不改变原正式面板、86.0262%旧合同或field-fusion unknown状态。
  
本日筛查条件专项审查见 [ctla4\_screening\_conditions\_audit\_20261004.md](ctla4_screening_conditions_audit_20261004.md)。固定96候选检查点实际44组合拒绝均无FoldX，其中39仅SPURS在0.5–1.0；21模型分歧中17双方各自稳定门均通过。30逐位审查中21未确认，其中18前五次≥2拓扑拒绝；fields.py全结构Euler不等将所有通道置NaN，保存11粗细网格对照10个Euler敏感。calibration.py/directions.py在四次阳性之外要求五次全有限，零求解反例[2,2,2,2,NaN]拒绝已复现。9逐位阳性仅1正式的主要合同瓶颈为全局部全通道完整或扫描上界闭合；shadow8组79.381610%不可升级。case/audits/screening_96保存4来源字节/hash、96行CSV、可重放脚本与PASS断言，新增科学调用0。此次仅报告/审查证据新增，活跃冻结源码、原科学门及正式17.660761%不改；局部拓扑、四次quorum、新下界资格、联合null与review调度为待验证建议，不称已修复或已科学达标。
  
本文件记录有效方向、已核对修复与边界；已放弃的实施建议从活动待办移除。历史输入、profile和审计保留，不能据此宣称所有输入都能产生正式推荐。

## 0. 当前用户目标与任务状态

覆盖路线后续修复：`engine.py`压缩审核通过后写独立compression标记，`local_constructs.py`按角色导出，不再错误要求补盲5%额外增益；对应与物理91项回归通过。ZIP以pending原子替换，失败保留完整旧包；独立解包核验新增真实残基/DDG、原始场/null内部指纹，校准/物理105项相关回归通过（数量有重叠）。计量扫描改为scandir加Windows os.stat，保留NTFS hardlink身份并测试原子消失；固定快照5次字节计量一致，中位耗时75→36 ms，仅文件计量提速。新五模型cohort只对审查白名单追加10/20；观测topology/WT查询/质子化/方向原因不改变资格。正式80%目标尚未达标，最新科学状态仍以完成生产检查点为准。

35组后续生产检查点：10组稳定通过、3组逐位真实改变成立，正式面板0；shadow联合下界51.7766%仅为对照。2375替换/125位点、作者到规范化模型及95个实际pair/DDG独立审计PASS；隔离273文件快照160项相关回归通过，655个内容对象只读hash审计通过。账本约96.6分钟、峰值1.44 GiB，后续合法局部三点E59T/A61W/M91W及精确缺失双点在同一总账继续验证，不能重置早期60分钟窗口或把shadow当作正式达标。

局部三点及缺失A61W/M91W各新增5对，共10对，均稳定通过但3/5及2/5拓扑拒绝，未获同四次逐位扰动资格。生产重放复现冷启动缓存恢复缺陷：外层预查空all_pairs使未被提案优先选中的缓存候选只有1对；重复数选择改到ensure恢复和身份验证之后，对新配置生产已有≥5对恢复5对，旧配置单对预筛保持。参数化及相关99项回归通过，实际零新增模型重放在同一总账继续验证；没有因shadow覆盖临时下降选择性沿用历史最高值。

上述实际重放完成，105个pair身份不变且新增0，三个逐位阳性/51.7766% shadow恢复，正式仍0。最新隔离162项相关回归通过。前五对≥2次topology rejection停止无效null追加，null-only后续失败不改变原五对资格；20模型cohort与确定性未知不占新扩展白名单。APBS分别报告空正式面板、执行/准备未请求、数值未评估，不误标完成。约113.4分钟/峰值1.59 GiB后继续八个新双点，未宣布科学或完整交付验收成功。

覆盖路线续验：隔离150项与24项编号/AF3/SASA/流程导出回归通过，真实检查点浏览器只读PASS。存储扫描复现账本原子替换临时文件消失，`case_resources.owned_bytes`只忽略FileNotFoundError、保留其他失败；竞态回归和实际恢复完成，971.11秒失败保留。最新27组筛选（另1组SASA跨阈值退池）、6组稳定通过、2组逐位真实改变成立，正式面板仍0；shadow联合下界26.85%仅作对照。当前账本约76.1分钟、峰值1.17 GiB。可靠WT分母实为6874.363277945836 Å²，14个阶段全部一致；早前报告误用旧路线6873.12，已在专项审计更正。保存网格11组scale2/3比较中10组Euler不同，仅为分辨率敏感诊断，不放宽拓扑资格；两组阳性shadow候选未进入相同资格/图/池容量的最低能量三项消融池，不能据此称物理算法提速或科学达标。详情仍见本轮专项审计。

2026-10-04覆盖优先CTLA4计划已明确获准实施，当前仍在逐级验证，详见 [ctla4\_coverage\_v2\_implementation\_audit\_20261004.md](ctla4_coverage_v2_implementation_audit_20261004.md)。独立任务 `E:\user\New project\ctla4_cov2_20261004` 使用P16410/1AH1作者链A、14400秒统一累计账本及20 GiB新增预算；不复制整套权重、不退款失败成本。本轮修复阈值传播/全序列保护/完整替换表、全图有界双点池、压缩与补盲分审、缓存与检查点对象所有权；实际复现FoldX长名启动崩溃、缩池丢白名单节点和失败阶段串用旧摘要，已分别以相同字节短名执行副本、先保留请求节点、当次mtime检查修复。隔离116项回归及真实失败批次恢复通过，实际原始场17组数值/NaN/mask等价审核通过。最近完整检查点20组双点、4组稳定通过，正式构建0项；存在null统计未解析和确定性拓扑/方向未知。后续合法三点、APBS、最终可移植包/真实浏览器及80%科学验收尚未完成，不能将工程回归称作整体成功。

- 用户目标：约8个突变构建体系统扫描绝大部分可靠ECD表面，不能只覆盖少数局部区域。
- **当前授权：用户随后明确要求实施冻结计划、继续五例实测并修复测试问题。** 同日较早的“先不执行”只描述当时计划阶段，已被后续实施请求替代。
- 活动算法计划：[surface\_perturbation\_algorithm\_next\_plan\_20261001.md](surface_perturbation_algorithm_next_plan_20261001.md)；[总体覆盖计划](surface_perturbation_coverage_plan_20261001.md)保留目标与计量合同，算法细节由新计划修订。90%足迹中心可扫描面积是待讨论工程目标，不是已确认政策或已验证实验检出率。
- 当前预览路线：native WT SES → 复用SASA/可靠性/SPURS → 重叠局部组合 → 配对WT/Mut和分通道校准 → 全局hard覆盖 → 真实组合补重复、资格复核与再优化 → 最终集合APBS审计 → 最多12项嵌套曲线与事后分区解释。详细合同与阶段结果见实现审计。
- 三个量独立：物理表面化学改变面积、footprint位置可扫描面积、真实实验检出与定位。直接接触与propagation-only分开；无抗体时不能从化学差异推出binding-loss probability。
- 主分母保留所有可靠暴露表面，包括保护位点；未知结构单列；不能删除困难位置、放大核、放宽稳定性或强选四点以制造达标。
- 8组是否足够逐蛋白评估；候选池上限、确认不足、solver下界与有效上界分别报告。greedy未达标不证明8组不可能。
- 现有V3净收益早停保留为历史/兼容实现；已另建`surface_perturbation_v1_preview`，不能将预览验证写成实验功能验收。

## 1. 当前已实现架构与术语

```text
原始结构 + UniProt 映射 + AF3 confidence（预测输入）
    → native SASA / 结构可靠性
    → Surface V2：点表面、重叠足迹、Surface V2 区域
    → FootprintFamily → InferenceHypothesis + 全局独立双点/三点候选搜索
    → SPURS single / multi、FoldX 确认与 formal gates
    → hypothesis response matrix → identifiability panel
    → schema 4.0 兼容投影 + 4.1 审计文件
```

- **显示区域**：互斥颜色归属，必须空间连通且直径受限。不是已验证表位。
- **科学假设**：Surface V2 raw footprint 经 family 与 conservative inference hierarchy 解释；同一残基可参与多个假设。
- **扰动矩阵**：几何与化学规则的 proxy，不是结合下降概率。
- **当前V3面板**：greedy信息增益与风险约束的结果，不宣称全局最优；CTLA4证明其广覆盖不足，下一步覆盖优先计划尚未实施。
- 区域数、假设数、构建体预算相互独立。固定 assay profile V2 的 Stage 1
    
  上限为 6；未指定固定 profile 的旧请求可保留其独立预算。
- 膜方向、糖链、PAE 是结构证据；不能据此声称预测 reporter 或亲和力。

## 2. 本次 PD1 / AF3 分区故障

输入：`workspace_data/fa85119485cd`，Q15116，ECD 25–170。
  
原模型 SHA256：
  
`8e4c7e0b922377d19359d2f4ad1bd598dadb69d0cdcb3b991aa1c5ea07b7c496`。

### 原因与证据

1. **页面着色的直接原因是导出编号错位。** ECD CIF 的 label_seq_id 从 1
     
   开始，而 record 保留 AF3 原编号 25–170；旧选择器优先用 label，导致
     
   本输入颜色选择偏移 24 位。实际导出前两项 auth/label 为 25/1、26/2。
     
   输入映射本身正确并不意味着导出结构的 label 编号仍一致。
2. `sectors.py` 为接近 7 区而合并小区域，未在合并前检查直径；最终检查只
     
   添加警告，没有修复。原 VS004 测地 / Cα 直径为 42.410 / 35.499 Å；
     
   VS007 为 39.432 / 34.631 Å，超出当时 35 / 30 Å 上限。
3. 原 VS001 的 25–31 与 129–130 在序列上相隔远，但 Cα 最大距离为
     
   16.821 Å；不能把序列跨度直接当作三维错分。原 VS002 Cα 直径为
     
   25.029 Å，显示范围也明显大于页面宣称的 ≤18 Å。
4. AF3 末端 25–31、150–170 的低置信度坐标没有达到严格 IDR 分类条件，
     
   因而仍进入三维几何。这混淆了“不能诊断无序”和“坐标足够可靠”。
     
   原 reliability QC 的 excluded_positions 为空，88 个残基进入表面图。
5. UI 的固定文案称所有区域都是 13 Å / 18 Å Fab 补丁，实际显示的是 V4
     
   sector。旧颜色只有 7 种循环复用；区域数增加时会造成远处分区同色。
6. 重放脚本把所有结构设为 `structure`；用于 AF3 会丢失预测模型语义和
     
   confidence 的作用。必须保留 `af3_upload` 与原 confidence。


### 修复规则

- `web/app.js.residueQuery` 使用导出后保持一致的 auth chain、auth_seq_id
    
  与 insertion code；只有 author 标识缺失时才回退到 label。不能直接
    
  将原 AF3 label 用在经过 ECD 截取 / 写出的 CIF 上。
- 显示 regions 时，不叠加未归属位点的通用 SASA 蓝色或 folding 绿色；
    
  不可靠尾部为灰色。否则通用蓝色会被误认为区域 1。
- `web/renderer.js` 等待实际结构 hierarchy 与 assembly 可用才首次着色，
    
  不能仅等待 managers 初始化，防止第一次 paint 在加载完成前被跳过。
- `structure_reliability.py` 的 `reliability_rules_v2` 将坐标资格与生物学
    
  无序分类分开。ECD 端部连续 UNCERTAIN 段只有在以下证据同时满足时，
    
  才退出三维足迹：所有已测 pLDDT <70、至少 80% 位点非局部接触数 ≤2、
    
  至少一个位点有经过合法 token 映射的 core PAE ≥15 Å。
- 该段仍为 UNCERTAIN，保留原结构、序列和审核能力；不自动标成 IDR，
    
  不产生自动 IDR 突变。高置信度移动结构域和嵌在核心中的低置信度短环
    
  不因这条规则退出。缺少 PAE 不能触发这条退出规则。
- `surface_patch/engine.py` 真正消费 `use_for_3d_footprint`；不可靠末端
    
  不生成表面样本，仍作为当前构象的遮挡原子保留。既有明确 DISORDERED
    
  段继续遵守 strict / occlusion_sensitivity 策略。
- 相邻 2 个保留位点作为边界 guard；不可靠末端自身不重复算作边界。
- `sectors.py` 在输出前拆分不连通 / 超径组，合并只接受满足全部空间约束
    
  的并集。默认显示直径：测地 ≤26 Å、Cα ≤18 Å。这是明确可配置的显示
    
  紧凑性规则，不是所有抗体足迹的物理尺寸定律。
- 区域数是软目标。不能为了 7 色牺牲空间约束；小区域若无法安全合并，
    
  保留并标记。额外区域采用不同颜色。
- `surface_display_regions_v41.json.regions` 保存实际 V4 显示归属，
    
  `surface_v2_regions` 保存原 Surface V2 聚合区域，两者不能混称。
- 兼容投影的 `mutable_residues` 不再包含硬保护残基；空列表不能在 UI 中
    
  被 `||` 回退成“全部可突变”。
- V4 兼容投影同步 `ResidueRecord.region_id`，使珠子视图与 Molstar 的区域
    
  membership 着色一致；不属于有效区域的位点必须清空旧归属。

### 验证资产

- `tests/fixtures/pd1_af3_partition.json`：原坐标、PAE、原图和原分区；
    
  PAE 按明确 token 映射截取 ECD，不用数组下标猜测 UniProt 编号。
- `tests/test_pd1_af3_partition_regression.py`：原故障图、端部资格、缺失
    
  PAE、实验结构、短环保留、mask 到 footprint/graph 的传递、软计数、
    
  颜色唯一和 AF3 重放语义；实际执行 JS 选择器并验证 ECD 导出编号。
- `workspace_data/pd1_af3_partition_repair_20260930/verified_geometry_complete`：
    
  原模型与 confidence 的真实几何重算；外部评分显式禁用，无正式突变推荐。
- 新图保留 58 个符合表面与可靠性规则的节点，仍形成 7 个有效区域；
    
  最大测地直径 25.226 Å，最大 Cα 直径 17.885 Å。
- 28 个端部位置退出三维足迹；边界 guard 为 32、33、148、149；全部
    
  146 个 ECD 坐标和原输入结构保留。当前输入未被诊断为 IDR。
- 页面验证副本任务 `da27bc835241`，明确标为“仅重算几何，未运行突变评分”；
    
  原失败任务保留。172 个假设均无不合资格位点；47 个 bundle 文件一致。

详见 `docs/pd1_af3_partition_repair_20260930.md`。

## 3. 之前的 V4 → V4.1 实现问题

| 原问题                                 | 已实现的修复 / 维护要求                                                        |
| ----------------------------------- | -------------------------------------------------------------------- |
| Surface V2 已计算但只传 `.surface`        | 保留完整结果，将重叠 footprint/region 传入 V4 adapter；不要另造科学分区引擎                 |
| sector 被当作科学表位                      | 显示与科学 membership 分离；schema 4.0 是兼容投影，4.1 是扩展审计                       |
| 共享 beam 挤占双/三点预算                    | 每 order 独立 beam 和 shortlist；四点为探索池                                   |
| option.score 混合所有证据                 | 保留暴露、化学、prior、稳定性和风险分量；内部排序分数不是校准概率                                  |
| panel 按 sector 轮询                   | 已由4.1 signature优化替换；广覆盖仍不足，活动方向见第0/7节                                |
| 只报告平均区分                             | 同时报告最小距离、未区分组、未覆盖假设；不得掩盖 minimum=0                                   |
| DiscoTope 缺失导致整体阻断                  | generic epitope prior 是可选软先验；保留 stability fallback                   |
| 将 ProCaliper 当抗体可达性                 | 当前 worker 是 residue SASA cross-check；不是 Fab 尺寸或进入凹槽的模拟               |
| 外部 SASA 方法参数缺失阻断所有候选                | native SASA + 显式几何为主；ProCaliper 方法差异为诊断，真实 native/normalized 分类冲突仍审查 |
| `_formal_candidate_reasons` 与测试参数漂移 | 兼容参数与契约测试同步；不要简单删掉失败断言                                               |
| ProCaliper atom radii provenance 缺失 | worker 与 context 保存物理设置 / 来源；不能编造 matched provenance                 |
| selected=[] 被误报为 FoldX 未运行          | 分开工具状态、候选覆盖、面板覆盖，成功评分 42 个不能报 missing                                |
| PAE / glycan / membrane 信息被丢弃       | 接入已有结构可靠性与 geometry evidence；实验 B factor 不当作 pLDDT                   |
| 几何均值访问 footprint 外残基导致 KeyError     | 只对权重 membership 与有效 records 的交集累加                                    |
| 用声明的模型 revision 代替真实缓存身份            | 权重与模型文件 hash 必须由实际输入建立；不能无依据重用旧评分                                    |

## 4. CTLA4 无推荐故障与修复

详见 `docs/ctla4_v41_repair_20260930.md`；关键证据如下：

| 原问题                                       | 原因、修复与回归                                                                               |
| ----------------------------------------- | -------------------------------------------------------------------------------------- |
| 新 CTLA4 无 single / DiscoTope              | `.venv/python.exe` 为 Low Mandatory Level，WSL 返回 `Wsl/E_ACCESSDENIED`；仅本机解释器标签恢复 Medium |
| 桥接日志只有空 RuntimeError                      | WSL 错误在 UTF-16LE stdout；桥接解码 stdout + stderr 并返回结构化错误                                  |
| 含空格路径 / shell 参数问题                        | WSL 使用 argv 与 `--exec`；不拼接不受控 shell 字符串                                                |
| ProCaliper 导入错误                           | `Select` 从 `Bio.PDB` 导入，不能从不存在的 `Bio.PDB.Select` 导入                                    |
| 同组合来自不同显示区域，multi 整批拒绝                    | 全局按位点与替换去重，保留多 hypothesis 归属；推理边界再防重                                                   |
| 空候选被当工具失败而无原因                             | 保存 NO_CANDIDATES_TO_SCORE / NO_CONFIRMABLE_CANDIDATES 与阶段诊断                            |
| native normalized SASA 依赖可选 ProCaliper 成功 | 始终运行 native 对照；可选工具故障不跳过结构语境审计                                                         |
| FoldX 预算被不可确认 / 四点候选占用                    | 先做 SPURS multi，主确认预算只给通过该门的双/三点候选                                                      |
| 同一位点组合反复跑图距离                              | 复用窗口内单位置距离及位置集合覆盖，氨基酸菜单不改变空间距离                                                         |
| 完成消息覆盖推荐阻断原因                              | pipeline、diagnostics、jobs 保留真实原因；计算完成不等于正式推荐                                           |
| ZIP 漏掉 V4.1 文件                            | 4.1 文件纳入 manifest 与 ZIP；最终 analysis 写入后重新封装 / 哈希                                       |

真实 CTLA4 修复重跑：`8c068ad701a6` →
  
`workspace_data/ctla4_v41_repair_20260930/verified_pipeline`。
  
210 个唯一组合、34 个正式通过、5 个正式入选；所有外部 worker 实际成功。
  
48 个导出文件的磁盘 / manifest / ZIP 一致。页面副本任务为 `05ecca7e4cc7`。
  
这证明特定实现故障修复，不证明其他输入必须产生组合。

本机完整性标签修复不随 Git 同步；程序不自动提升权限或重置整个目录 ACL。
  
长时间 multi 推理不要反复启动重复进程；先检查进度与模型实际日志。

## 5. 后续必须执行的检查

1. 先确认原始输入类型、chain/auth/label/UniProt 映射、confidence token
     
   映射和文件 hash。PAE 无合法映射即忽略，不能凭序列下标猜测。
2. 显示 membership 对 eligible graph **完整、唯一**；每区 induced graph
     
   连通，所有两两距离满足直径上限；不能只打 quality flag 后照常输出。
3. geometry-ineligible 的样本不能经 graph fallback、transit contraction、
     
   sequence bridge 或 hypothesis adapter 被重新引入。
4. 用真实故障 fixture 验证，再核对完整输入重算；历史运行不可原地改写。
5. 不因“必须推荐”而放宽稳定性、保护位点或外部证据门槛。有效空结果应
     
   明确区分真实不通过、工具缺失、几何资格不足与优化未获得信息增益。
6. 只更新文案是不够的；反之修后端也必须核对 UI 归属、颜色、位点计数、
     
   JSON、CSV、manifest、离线报告和下载 ZIP。
7. 当前工作树包含大量其他功能的未提交编辑，不能 `git add .` 混合提交。

### 测试与重放

从仓库父目录 `E:\user\New project` 运行，避免包导入错误：

```powershell
python -m pytest -p no:cacheprovider ecd_epitope_scanner/tests/test_joint_v4.py ecd_epitope_scanner/tests/test_joint_v4_sasa_context.py ecd_epitope_scanner/tests/test_surface_snapshot_v4.py ecd_epitope_scanner/tests/test_surface_patch_v2.py ecd_epitope_scanner/tests/test_structure_reliability.py ecd_epitope_scanner/tests/test_v4_hypothesis_panel.py ecd_epitope_scanner/tests/test_ctla4_v41_regression.py ecd_epitope_scanner/tests/test_pd1_af3_partition_regression.py -q
```

从仓库目录验证原 AF3 几何，不运行评分：

```powershell
.venv/Scripts/python.exe scripts/replay_joint_v41_run.py workspace_data/fa85119485cd workspace_data/new_pd1_geometry_validation --geometry-only
```

去掉 `--geometry-only` 才会运行真实评分、候选门槛和面板。
  
重放必须使用新目录。`replay_input_provenance.json` 记录实际计算范围与 hash。

## 6. 仍需谨慎的边界

- 缺少 PAE 的预测输入不会由这条 AF3 末端规则自动删除；应保留不确定性，
    
  不能称作已验证的表位几何。
- 高置信度但相对方向不确定的结构域仍需 PAE / 多构象审查；本次 PD1
    
  修复不等同于完整结构 ensemble 验证。
- 低置信度末端保留为遮挡原子，因此“移走尾部”仍可改变表面可见性；
    
  本次采用保守的给定构象，不声称对此类变化严格不变。
- 受保护 / 埋藏残基造成的颜色空隙不等同于空间断裂；连接判定基于已
    
  记录的显式表面图及有界 transit 路径，不是序列连续性。
- 70 项相关回归通过不等于整个仓库所有实验工作台 / 流式功能均已验证。
- 结构场景与真实benchmark尚未完成；没有真实完成的数据不能写成已实现结果。黑箱binding/report训练不列入当前活动计划。

## 7. 2026-10-01：资格修正结果与低覆盖审查

已完成实现与证据见 [v411\_upgrade\_20260930.md](v411_upgrade_20260930.md)；新的广覆盖方案是待讨论计划，不混写成已完成能力。

### 已实现且继续有效

- native total rSASA默认0.25；侧链面积只描述。Gly零侧链面积不排除暴露；GP具体替换及骨架风险政策仍保留。
- 内部低pLDDT但几何可用位点进入REVIEW_REQUIRED池，不进入严格面板；不可靠AF3末端仍退出三维图。
- raw/family/inference及显示ownership独立，family工程权重1/0.7/0.3；complete-link不跨组件/域或不兼容可靠性。不为降数量放宽合并。缺少合法预测PAE为UNAVAILABLE，不能None→0伪造LOW。
- 全局按位点+替换去重，来源假设合并；double/triple独立shortlist，60默认可为保留配额扩容，空order也报告缺口。复核不占正式名额；四点exploratory。
- 现有V3默认4/7/8、固定√8距离；V1/V2历史选择语义保留，新参数优先。下一步广覆盖策略不会原地改写这些版本。

### CTLA4覆盖不足的证据

- V3 strict 0.25：11显示区、125 inference、245唯一候选、26正式通过、3入选，worker无失败。12 Å旧proxy为15/71节点=21.1%，面积22.5%；84个假设低于当前响应阈值。
- 第四个候选最高净收益：信息+0.04610，风险−0.00369，冗余−0.05842，净值−0.01601。因此停止符合旧规则，但未实现用户的广覆盖目标。
- 共同零响应和近重复hypotheses可使互补探针被过罚；后续优化应先满足覆盖/盲区目标，再考虑区分度和合格候选间成本。不是简单删除风险或强补8个。
- 显示区anchor不是实际覆盖；无anchored正式候选不证明其他探针不能接触该区。应按实际包络/接触支撑报告缺口。
- 0.15/0.20/0.25真实strict比较：唯一候选258/250/245、正式通过28/28/26、均3入选，worker无失败。降低门槛改变组合但未解决广覆盖，默认不变。

### 新扰动体系必须遵守的科学边界

- 具体替换先用可追溯化学描述；范围/形状需局部结构与成对WT/Mut场支撑。没有通用A→N或S→N真实扰动半径表。
- EpiScope约12 Å是组合位点间距；耦联研究15–20 Å不是抗体结合效应范围。现有6 Å衰减/15 Å截断仅保留作旧工程对照。
- FoldX需保存对应、相同邻域处理的WT/Mut结构，DDG不是物理差异场或binding ΔΔG；WT-null只能校验数值噪声。
- 非圆形多点支撑不填凸包空白；组合非加和和长程结构变化需单独风险审查，不自动计为额外覆盖。
- 面积权重守恒，重复尺度不重复计面积；在独立密集中心集验证覆盖，不只在优化用的125个hypotheses自证90%。
- 几何截获不保证命中未知抗体关键能量接触；可能结合增强或降低。集中/分散贡献场景是先验，不是所有抗体的保证。
- 复合物真值只能验证几何；真实functional detection需实测突变、表达/构象数据及family holdout校准。当前benchmark仅输入契约，尚无真实性能证据。

### 同日算法科学复审：当前核对与活动政策

- 当前Surface V2是SAS点云，默认48几何点/原子；native SASA与之独立。`molecular_mesh.py`仅报未安装，SES尚未实现。新增SES必须并行验证面积/图收敛，不混用SAS/SES分母，不宣称水探针等于Fab可达性。
- `stability.py`的`engine_runs`持久保留原始文件；不能称FoldX模型已被临时清理。当前worker只返回DDG，缺少逐候选/run、PdbList、对应WT/Mut、hash/映射的结构产物合同。先审计已有文件，不能凭路径或排序假配对。
- fixed WT表面用于统一计量，实际形状/化学变化在共同局部三维壳层查询。仅最近点匹配会跨凹槽/漏新增侧链；占据差、合法法向追踪及未知mask优先，不逐下标减两张网格。
- 对应WT核心对齐到native参考后，同一变换用于其Mut；Repair/WT重排漂移另报，不计突变覆盖。工具的骨架自由度与重复采样范围要声明，固定骨架场不能证明没有骨架/远程响应。
- 暴露、折叠稳定性、接触变化分别输出，但不声称物理或统计完全独立。DDG、SASA、形状/电势单位不混加；有符号通道、run一致性、null与收敛独立保存。
- 单点场max并集只是筛选proxy，不是组合真实效应下界；需真实组合配对。直接变化、邻域重排、长程模型变化与未解析分层，长程与未知不自动计强覆盖。
- 方向性虚拟接触探针是待验证增强，不是整只抗体或binding DDG。MaSIF embedding距离、任意核尺度、独立重跑AF3和全蛋白MD不列入起步方案；复杂模型必须经同池公平消融/holdout证明收益。
- 本轮仅完成代码只读核对、原文核查与文档一致性检查；没有运行新物理场、新增算法测试或宣称提升覆盖。新计划及上文科学合同均待实施。

### 已验证资产与维护限制

- PD1资格解耦后原AF3几何为62节点/7区、32个末端/guard位置不入图；先前58节点为早期修复baseline，不能混为同一次。此验证禁用评分，不证明突变推荐。
- CTLA4三阈值显示区11/13/11均连通有界，61个bundle文件一致；UI只读fixture无JS错误。现有独立Git检出97项回归通过，是实现规则回归，不是高覆盖/实际检出率验证。
- surface centers约7 Å、合并上限6 Å，保守inference仍多（CTLA4 125、PD1 94）；不强制压缩。
- 长期worker源码hash在导出时采样，不能单凭hash声称加载了后续编辑；保存启动快照并做独立Git检出。
- 已遇到pipeline引用未提交链识别函数的依赖缺口；交付必须脱离脏工作树验证。仅提交授权内容与必要依赖。
- API测试用项目venv，从父目录执行；系统Anaconda缺FastAPI，Pydantic要求2.x。本地服务旧进程不自动加载新模块；重启曾被自动审批以blocked by policy拒绝，不得绕过。

## 8. 2026-10-01：WT参考扰动场独立预览实现

实现与当前验证范围见 [surface\_perturbation\_v1\_implementation\_audit\_20261001.md](surface_perturbation_v1_implementation_audit_20261001.md)。
  
上节“SES尚未实现/worker只返回DDG”是先前只读审查的历史状态；新策略已新增真实NanoShaper SES、FoldX配对manifest、分通道null和真实组合再优化，不原地改写旧策略。

- `surface_perturbation_v1_preview`独立接入API/UI/报告和导出，原正式门槛及预览身份保留。APBS只最终集合审计，缺失/噪声失配/介质切换均未知。
- FoldX单次BuildModel的文件名可能只有group而无run后缀；按实际PdbList、突变列表及DDG模型名联结，不按排序配对。
- NanoShaper可能含精确重合顶点与零面积面；只允许精确焊接/删除零面积面，禁止填洞使非闭合网格伪装合法。
- 重复证据身份是batch/group/run，不能把相同坐标文件去重后误称未完成重复；合法缓存必须核对结构、分数及完整运行配置。
- 五例独立检出重放在进行，曾发生WSL执行失败与连接重置；保留累计预算和失败目录，串行恢复，不伪造模型结果。
- 最新独立源码覆盖检出69项相关回归通过，工具烟测及CTLA4面积/图代理采样审计不证明覆盖达标。三条独立路线物理对照、接触标签事后验证、密集中心真实场收敛与盲区扩池尚未完成。
- 校准预算与生产7200秒限制分离；只更正误套用的校准上限并保留已计时间。中断未关闭间隔保守计入预算，不把保守wall-time当作worker性能。
- 三维方向通道按联合null幅值及同四次方向一致判断。介质切换按实际PQR半径SES检查；APBS平滑网格界面依然需数值误差审计。
- 直接原子支撑与变化超阈值必须共享同至少四次模型；被排除三点候选退出正式池。独立K对照与其矩阵采用全真实验证池。
- 跨校准/生产阶段合法原始场复用必须核对参考、数值实现、NanoShaper及内容hash，分类总是重新审核；替代路线比较共享每例限时账本。
- 五例隐藏接触标签已经冻结，但事后验证仍待正式面板。B7-H6/NKp30与EGFR/EGF是受体/配体几何，不能冒充抗体表位；EGFR原始7SYD声明D256N冲突，标签映射排除并记录。
- 2026-10-01加速修复：CTLA4原始author1A先于1；PQR准备必须按native记录顺序，不能排序author编号。非突变准备漂移保持未知；旧CTLA4耗尽生产预算，未完成面板，不追加时间制造验收。
- `surface_perturbation/execution.py`只对独立APBS WT/Mut做最多2任务并发，按host/WSL内存限制、WSL单solver线程1，相同PQR去重避免缓存竞争；`runtime.py`预算写入加锁。实际PQR SES溶剂查询缓存核验完整数值身份，保留未知mask，不能用缓存替代重复证据。
- `deadlines.py`在三角查询内检查累计预算，数值算法不变。最近隔离检出75项相关回归通过。单对真实APBS工程试验串行41.37秒/并行28.45秒（1.45倍），采样差0；SES查询缓存3.51秒/0.015秒且数组一致。非整个流程加速率，也不是CTLA4科学验收追加运行。
- CTLA4最近尝试中APBS记录约465秒，不能从总累计wall-time推断它为最大瓶颈；配对场、溶剂SDF、APBS任务对分别计时。进一步全流程性能验证与五例正式验收仍待完成。


## 9. 2026-10-01：继续实测期间的修复与证据

- **PD1旧生产结果已撤销当前资格：发现8项中R88C和S38C违反原有自动替换政策。** 旧累计账本6512.53秒、主5项90.47%和嵌套8项92.35%仅为保存历史，不能作为当前正式覆盖率。八项APBS审计均完成但准备合同未通过，extended未知。必须按修复后资格重建designable分母、重新优化并复核；不得删除两项后继续使用旧分母宣称验收。
- `case_status.py`组合物理summary与全流程预算验收。物理阶段完整但导出超时不得记作病例完成；接触标签事后检查遵守同一状态。独立检查PD1外层68文件manifest/ZIP和内层1233项/102配对模型均通过；一致性检查不证明科学充分性。
- CTLA4生产预算已耗尽。旧ZIP缺少部分旧配对行的DDG表引用；只在真实模型名与原始DDG唯一匹配时补入已有表。行政性修复ZIP1127项/90配对模型通过独立核验，保留原ZIP和未完成状态，没有恢复生产计算或重新授予正式面板。
- `classification.py`修复每个直接变化位点的四次支撑与超阈值集合可能仅交三次的问题。必须同一位点共享同四次；标量/方向回归及CTLA4、PD1各八项保存场重放通过。CTLA4单组合独立密集场重分类不变，中心覆盖差0.37个百分点；不推广为全五例/全板收敛。
- `query_cache.py`保存原始查询批次，并允许严格逐点状态缓存：坐标、批次旋转结果、最近面编号及mesh/数值库身份均一致才复用；部分未命中仍计算原始整批。保存CTLA4原生mesh三组查询距离、valid与法向逐元素相同，缓存查询41–59倍；只是局部重复查询的性能证据。非闭合mesh和未知mask保持未知。
- `runtime.py`采用独立临时checkpoint与Windows短暂PermissionError重试，旧完整账本保留；不能把并发读文件造成的重命名失败当成科学模型失败。最终清理即使保存失败也释放锁/恢复预算上下文。
- `buffers.py`仅为同候选、相同WT/Mut文件哈希和链共享数值数组，所有真实run身份仍保留。校准任务结束后释放大数组，原始场、null与模型文件不删除。
- 扩展测试暴露旧断言把显示区数绑到构建预算、要求UI写死7。改查实际连通/直径和独立预算；fixture已复验，浏览器拒绝legacy上限12并接受8，新策略接受12，不收紧显示数来使旧断言通过。
- B7-H6原请求为已停用`joint_v3`。benchmark显式选择本次V4预览方法，仍保留原AF3、confidence、chain、输入及原请求hash，并记录旧/新引擎；失败时间不退款。
- benchmark恢复时将本病例输出加入合法配对缓存。`foldx_worker.py`按实际缺口补模型、保留已有模型身份与manifest历史，新模型写独立目录；完整缓存、部分缓存和输入不匹配分别处理。真实工程回归请求5对、复用3对、新增2对，全五对身份/hash/DDG通过；不是CTLA4生产重跑或验收。
- 主机/WSL可用内存约16/15GB时，将余下实测安排为两个独立检出队列，GPU继续共用`/tmp/ecd_v4_gpu.lock`，APBS每例最多2任务且按内存限流。只终止已核验的本任务协调/worker进程，校准恢复保留累计时间与源码启动快照；数值七模块未改动。并行报告按病例分文件避免临时文件竞争，不宣称整体加速倍数。

原始证据位于`workspace_data/surface_perturbation_20261001`下的`frozen_benchmark`、`delivery_repairs`、`dense_field_audits`、`triangle_point_state_engineering_audit`、`foldx_append_engineering_audit`和`cache_resume_before/after.json`。PDL1/B7-H6/EGFR及跨路线独立物理对照尚未完成；最终状态以每例当前启动manifest、累计账本与验收记录联合判定。

- 资格漂移修复：V4单点菜单和正式缓存复核共享stability.AUTO_EXCLUDED_MUTANTS，禁止新引入C/P/G；原C/G/P保护继续保留。PD1两个旧组合违规，当前case_status拒绝旧VERIFIED_PANEL。配对产物及旧ZIP保留，重新生成设计资格和分母后才允许复核。
- 实际失败OFF仍在physical_cache：B7-H6有3条边界边或2条非流形边，PDL1有3条边界边加1条非流形边或单独1条非流形边，PD1有3条边界边；没有可凭翻面修复的纯方向问题。失败任务隔离，不填洞、不赋予未知SDF已验证覆盖。详见posthoc_audits/failed_mesh_topology.json。
- CTLA4旧hard草案亦含K29C（V41_6e6722e23d0da1b4），违反新Cys禁入政策；本来就因预算耗尽未正式验收，现再明确该旧草案不能被用作当前合格面板。密集场工程单例33ff本身未涉及该违规项。独立检出234文件、103项相关回归通过（surface62、V4 30、AF3 6、API 5）。
- portable evidence导出改为递归包含calibration/task_failures及恢复历史。否则仅打包旧null数组而漏掉失败屏蔽会使解压重放错误地重新启用阈值；对应ZIP成员内容回归通过。物理预览跳过正式资格根本不消费的高阶代理全池比较；真实组合额外收益门槛仍保留，旧V4模式仍计算原代理比较。

## 10. 2026-10-02：verification state machine修复

- 用户再次要求审查建议、修复并跑实例。校准valid.any不等于请求范围收敛：calibrate_channel增加required_mask和面积，默认全部点；在候选direct-support并集检查所有请求位置，分别报告有效/未知面积比例。每个声明hard通道均完成局部校准才能授予候选完整资格，部分解析数据仍保留为预览证据。
- 全池收敛不能由已校准子池身份不变证明：任何未完成局部校准的正式候选均阻止VERIFIED_PANEL，正/零proxy都不是排除上界。返回CALIBRATION_EXPANSION_REQUIRED并保存请求；协调器在原校准/生产预算内定向补校准再优化。最大重复后请求不变或预算耗尽则报告未解析，不能改名为完成。
- APBS状态拆成执行完成、准备有效、数值已解析、extended已解析；未证明代表项误差界可覆盖整个面板，因此不将其误差界转移给其他候选，仍保留未知。
- 修改calibration/null_store后数值七模块身份变化。显式保存修改前七源码与hash，用legacy_fields校验五个物理计算源码字节完全一致和numerical_bounds AST一致后，仅迁移内容hash相同的raw fields/数值误差数组；旧threshold和converged不迁移。历史文件不改；新contract保留迁移来源。配置单列frozen_settings_20261002_v2.json，原7200生产和14400校准累计上限不变。
- 对中断后已消失的本任务进程，关闭账本并保守计入未关闭间隔，不退还时间。新版PDL1、B7-H6、PD1→EGFR分别从独立检出继续。病例结果仍待实际运行，不能据新增回归宣称五例通过。
- 增加严格hard扫描上界：至少五对真实配对模型，按同一直接原子4 Å支撑、同四次重复条件和canonical图12 Å扫描规则计算保守种子上界（不同残基并集合并只会扩大）。仅当未知通道不能增加该候选已验证hard足迹，或该候选上界完全被当前面板覆盖，才允许关闭对应hard盲区问题；未知通道不赋0、不变成完成校准。校准范围仍严格100%请求位置，各通道及联合向量的未解析比例单列。上界不证明实验响应、不证明风险/可辨识性或各K最优。
- 对配对物理模型异常明确标formal exclusion及原因，不能从候选矩阵中静默消失仍显示正式合格。直接支撑也复核可靠核心RMSD≤0.5 Å。候选上界、局部校准和闭环请求进入独立版本身份与导出，历史ZIP及reference/manifest在重试前保留快照。

### 2026-10-02继续运行后的收尾

- PD1 production真实触发WinError 5：`paired_manifest.json.tmp`替换被拒绝，科学复核未完成。`util.write_json`改为每次独立临时文件、关闭句柄后替换、有限PermissionError重试；持续拒绝保留旧完整JSON与待写payload。并发40次写入、瞬时占用、持续占用回归通过。不重置PD1预算，恢复尝试仍未验收。
- 新增`scripts/summarize_surface_validation.py`联合当前launch时间、summary、账本和case_acceptance汇总；新尝试期间旧VERIFIED_PANEL不得提供覆盖率。`case_status`另从保存blindspot请求拒绝虚假VERIFIED_PANEL。主K覆盖取真实嵌套曲线中actual_constructs匹配行，不冒充K=12结果。
- 运行期间观察到长壁钟间隔，不能把它当作求解器耗时或加速证据；原因未从系统事件确认。按原合同保守计时，PD1、PDL1、B7-H6生产及EGFR校准已经超限。只停止核验过的本人任务进程树、关闭账本、保留检查点，不退还时间。外层`surface_phase_guard.py`进一步对阻塞子worker监督累计限时；已超限不启动，关闭账本不重复累加同一执行间隔，尚未打开账本的子进程启动等待也计费。
- EGFR校准未产生任何匹配null数组而校准预算耗尽；生产尝试约461秒后明确关闭为INCOMPLETE_CALIBRATION_BUDGET_EXHAUSTED，尚有生产余额并不补偿缺失校准。仅在没有null数组且校准预算耗尽的确切前提下省去无从完成正式验证的生产；存在数组仍由engine逐项审核身份，不能靠预检推断有效。
- 当前五例均未正式验收。CTLA4、PD1、PDL1、B7-H6为生产预算耗尽；EGFR为校准前提无法补齐。PDL1此前确实返回CALIBRATION_EXPANSION_REQUIRED，B7-H6为CALIBRATION_INSUFFICIENT，证明修复不会把局部子池或稀少有效点包装成完成。PD1旧Cys面板继续失效，所有未完成病例不报告正式覆盖率。事实汇总见`workspace_data/surface_perturbation_20261001/frozen_benchmark/current_acceptance.json`，关闭原因见posthoc_audits。
- 最新隔离源码覆盖检出119项相关回归与34项扩展回归通过，新增缺失校准预算预检所在4项回归亦通过（其中3项与119项重叠）。科学物理架构与数值七模块未在收尾阶段改变；回归通过不代表五例或抗体功能验证通过。局部null覆盖缺失/重采样误差未解析、部分网格不闭合、完整独立路线物理对照仍是方法冻结的限制。

## 11. 2026-10-02：运行根因复现与局部查询修复

详见 [surface\_perturbation\_runtime\_root\_cause\_audit\_20261002.md](surface_perturbation_runtime_root_cause_audit_20261002.md)。本节更新上节“数小时间隔原因未确认”的历史状态，不修改旧账本。

- System的506/507事件确认02:13至09:02为Modern Standby，时间调整0ms；不是数小时solver运行。新增power_guard在pipeline及外层phase guard期间申请System/Execution并于异常/退出释放，不改显示或全局电源计划。手动睡眠及电池策略仍可覆盖，原预算不退款。
- 保存PD1配对manifest副本使用真实不共享DELETE的Windows读句柄复现WinError5；释放句柄后有限重试成功，原文件hash不变。具体历史持锁进程未记录，不能指认。原子JSON回归新增原生句柄案例。
- WT化学部分NaN可为His/方向弃权，不等于未查询。按极性支撑判断实际查询：PD1最差99.2272%、PDL1旧cohort95.8446%、B7-H6均完整。scoped_fields按真实direct-support并集补配对查询，旧raw不可变、新null协议不复用旧threshold，查询scope进入缓存和误差身份。PD1真实20对重放补至100%，25.315秒、原hash不变且无新增FoldX/APBS/PDB2PQR求解；不代表null或面板完成。PDL1旧cohort缺当前配对manifest，拒绝推断。
- 新合同paired_WT_native_SES_preview_1.2及v3冻结配置保持7200/14400预算；数值误差仅查询所需支撑区，范围扩大即重算。His/投影、20次拓扑变化、真实mesh破洞及PQR重排仍unknown，不放宽资格。隔离检出127项相关与12项流程/导出回归通过，源码身份见运行根因audit。五例正式验收仍未完成。


## 12. 2026-10-02：四例完整验收请求与目标复核

用户排除EGFR，要求审查目标并完整运行CTLA4/PD1/PDL1/B7-H6。目标缺口见 [surface\_perturbation\_four\_case\_acceptance\_20261002.md](surface_perturbation_four_case_acceptance_20261002.md)，原始输入、AF3/confidence和hash已核对。旧四例生产账本全耗尽；按原计划重试不清零，完整新批次预算约束已向用户提问，未取得前不启动改变该合同的生产计算。

- engine与benchmark新增显式raw_field_cache_roots，只按完整原始物理身份复用，不迁移旧threshold/资格或重置预算。75项相关回归及浏览器检查通过。
- 原CTLA4 ZIP漏DDG表，但源配对合法；当前导出器另目录工程重打包1127项/90对独立通过，原manifest hash不变。旧包保留。PD1旧ZIP和失败续跑的新输入清单不可混用；PDL1旧包一致、B7-H6缺完整生产包。工程导出不代表正式新验收。
- 约8项覆盖绝大多数表面、四例新正式panel、独立路线物理比较、接触标签核验、电势完整解析及参数冻结尚未达成；不把接口/回归通过称作所有科学目标完成。
- 用户后续明确“确认之前测试中出现的问题已全部修改，然后重新跑四个样例”。据此独立新批次rerun_20261002_191524启动四例，原耗尽账本不改；新批次仍7200生产/14400校准、同批重试累计。启动前127+2回归通过、250文件冻结。continue_surface_validation增加逐文件验证的--frozen-source，两个队列共享只读检出，避免并发freeze写入worker源码。当前运行不等于四例验收，状态看新批次各phase而非旧ZIP或旧summary。
- 重跑实测发现FoldX默认二进制路径错误地来自隔离检出，那里不复制vendor工具。benchmark新增frozen_worker_parameters，启动前从保存输入工作区解析真实工具、核对配置hash并写入input manifest；缺工具快速拒绝，不再浪费上游计算。真实binary hash与原合法缓存一致。13项依赖/冻结源码/状态/guard回归通过；已消耗时间保留，核验后只停止本人队列进程，settings_v5在同一新批次账本重启，251文件冻结。详情见四例验收审计，四例仍待结果，不能称为全部已完成。
- v5两例真实配对后复现同重复数ensure的UnboundLocalError：追加循环内局部导入repeated_field在空循环路径未执行。改模块级导入，新回归覆盖5、5→10→20及同数复用；87项相关隔离回归通过。coordinator仅审计当次manifest+ZIP；case_status拒绝当前driver失败及旧summary，旧包仍保留。数值七模块、合法fields/null不改，252文件冻结，settings_v6在原新批次账本继续，不退款。详见四例验收审计。
- v6冻结检出父目录+原venv的53项V4/API/页面/完整流程/PD1 AF3回归通过（360.62秒）；唯一既有Starlette httpx warning。与87项专项回归分别记录，不把工程通过计为四例科学验收。
- v6真实发现FoldX长路径buffer崩溃：相同A30Y/K29W输入，cwd239返回0xc0000409且无产物，cwd232/87各4次成功。worker长目录改持久化短执行根，原output manifest不变、保留模型路径/hash。guard预算读取PermissionError可能留下孤儿child，新增读重试及异常owned-tree清理/计费。80项隔离回归通过、253文件冻结，settings_v7同账本恢复，旧失败和花费保留；详情见四例验收审计。
- v7实测扩池调度重复最大20 null：PDL1第二轮12项中10项重复，仅2项新增；CTLA4第二轮236模型数不增。v8新增同scope/correct-fingerprint最大20或完整null的调度defer（仍unknown、不授予资格），优先未处理任务，coordinator即使requests不变只要尚有未处理项则继续。91项隔离回归通过、255文件冻结；实际保存数据CTLA4 32请求/20 fresh/12 deferred，PDL1 35/23 fresh/1 retry/11 deferred。旧快照阶段完成后source guard阻止下一阶段混用，v8同预算恢复，科学七模块未改。详见四例验收审计。
- 扩池跳过最大20 null时，补充源模型恢复/manifest历史保留，避免当前1次记录覆盖原20模型用于null和导出。resume_pairs严格mutation/input/normalized/binary+逐对validate；V4注册其真实身份至supplement exclusion，稳定性只核查当前field模型前缀（最终5）。37项相关回归通过；实际恢复CTLA4 203、PDL1 198额外模型，保存身份缺失0，源不变。见native_model_resume_replay.json/四例验收审计；当前v8冻结阶段按原源码完成后再切换，不混代。
- v9在独立detached HEAD 9f01e12+257文件overlay检出验证122项全部通过（12.87秒）；A先继续、B旧v8阶段完成后切换，同账本、科学七模块不变。原隔离检出及已完成v8产物不覆盖。进度见四例验收审计和新批次lane_A_v9。

2026-10-03续跑：CTLA4 v9两层ZIP审计PASS，但23项局部证据仍未关闭，正式panel未验收；详见[surface\_perturbation\_four\_case\_acceptance\_20261002.md](surface_perturbation_four_case_acceptance_20261002.md)及同批CTLA4_calibration_resolution_diagnostic.json。PD1/PDL1/B7-H6由同批v9队列继续，预算不重置。已知缺陷源码/回归清单known_defect_verification.json，不能把工程PASS理解为全部科学目标完成。

2026-10-03补充：System事件确认00:44关机/08:32启动，是本批worker消失的原因。保守未关闭计时不退款，PD1/PDL1校准超限，随后生产真实导出均PASS但科学仍CALIBRATION_EXPANSION_REQUIRED。B7-H6 Hidden独立后台队列继续。默认profile8/历史V1 6的旧tracked测试已修正；全包393通过/6失败/11跳过，FlowKit隔离依赖alias后4项复测通过；余2为已有未跟踪FlowJo空单元格断言，未改其源码。见四例acceptance audit/full_package_regression_v10.json，不称所有目标已完成。

2026-10-03正常退出计时补修：原生guard重放初始7秒/执行0.517秒后仍7，现按启动区间累计下界与子账本取max、正常退出超时亦拒绝；修复账本7.5115秒关闭。`continue_surface_validation --start-phase production`只在当前校准完整交付且账本关闭时从阶段边界恢复，部分null仍由engine审核，不跳科学门槛。v12隔离257文件专项138通过（13.02秒），仅guard/协调器及两份测试相对v10改变，数值模块不变。B7-H6当前阶段完成后再安全交接，同账本无退款；真实交接/最终验收以campaign记录为准。详见四例acceptance audit与guard_clock_replay、v12_regression.json。

2026-10-03交接helper修复：后台日志包装Popen为函数会破坏Windows asyncio子类导入，尚未启动科学子进程即失败。改为Popen子类，watcher拒绝旧启动metadata；新解释器重放通过，lane B13真实生产已启动。冻结物理源码与cal7413.594/prod358.704累计不变，旧失败保留于launcher_type_repair_20261003.json及failed metadata。B7-H6仍待最终状态；失败OFF另存hash副本，未完成工具cache重试会覆写OFF，不能把旧hash当作现文件身份。

2026-10-03 v13：修复向量资格坐标依赖（45度旋转旧scope false→true、同联合分类不变）。物理资格统一为标量及联合向量，分量仍诊断，要求全部实际支持有效；七份数值源码与v12完全相同。新增严格VERIFIED_LOCAL_CONSTRUCT的JSON/CSV/UI/双ZIP交付，局部个体资格不冒充全局panel收敛，UniProt/auth/FoldX编号分别保存。259文件隔离overlay专项149通过、API/UI七项及浏览器PASS。同账本两条生产队列重核四例，结果与限制见[surface\_perturbation\_four\_case\_acceptance\_20261002.md](surface_perturbation_four_case_acceptance_20261002.md)；当前不能把启动视为完成。

2026-10-03本轮四例重核已结束：CTLA4四项、PD1三项、PDL1三项（新增L27V）、B7-H6一项，11个严格局部构建体独立通过编号/原序列/五对模型/DDG/mask/双ZIP核验，CSV及ECD FASTA已交付。PD1 R112M/N/Q是同一位点备选，不能当作三区覆盖。native浏览器发现旧fixture漏测真实字段位置：API已有`design_summary.surface_perturbation`，页面只查joint_v4_ml；v14客户端读取真实字段优先并保留旧位置兼容，fixture改真实schema。v14仅app.js与browser fixture不同于v13；隔离API/UI七项复测通过，四例真实页面全部PASS且无API写入。生产账本均关闭且未超限，全局面板仍CALIBRATION_EXPANSION_REQUIRED，未获得约8项多数表面覆盖或最终APBS解析。当前交付为同批`final_four_case_acceptance_v13.md/json`、`four_case_local_construct_verification_v13.json`及`verified_local_constructs_v13.csv/fasta`，历史报告保留；详见四例审计。

2026-10-03组合目标修复：用户明确授权2–3点、长度自适应4–12项、hard/all≥80%的独立四例批次，旧11项局部单点结果不当作达标。combo_contract/engine/candidate_search加入逐位同组4/5贡献、真实三点双点对照、完整图跨度、双点优先确认及goal/search独立状态；null从原始样本重审、DDG正式前五对，双ZIP保留导入hash。回归复现累计计时浮点抵消已修复。隔离v2 detached9f01e12+261文件overlay，父目录项目venv167项通过（14.61秒），浏览器逐位证据PASS。新批batch_20261003_155853两队列CTLA4→PD1/PDL1→B7-H6，校准14400/生产7200同批累计；v1校准阶段完成后安全交接v2不覆盖正在运行源码。记录时前两例校准仍运行，尚无80%目标验收；详见[surface\_combo\_panel\_implementation\_audit\_20261003.md](surface_combo_panel_implementation_audit_20261003.md)。

组合批次后续补修：首轮校准现同样跳过同协议/同scope完整null或20次上限（未知不授予资格），首轮 capped组合回归通过；真实三点低阶对照允许严格hard扫描上界已闭合的实测对照，仍须三个精确双点子集，不借单点并集。v2最终专项169项通过（14.57秒），另AF3/完整流程/SASA16项通过（180.67秒），native UI浏览器PASS；安装包补收surface_perturbation，42模块wheel源码字节核验PASS。两条有限队列及阶段交接保持同账本，finalize_campaign只读验收/导出和真实浏览器在队列结束后生成final_combo_acceptance，代码通过不代表四例目标完成。记录时CTLA4/PDL1仍校准，PD1/B7-H6待队列，详见组合审计。

2026-10-03 18:16资源恢复：本批日志明确磁盘空间不足OSError28/Win112及损坏ZIP，尚无RAM OOM证据。用户要求继续并行，两条v2病例队列原账本恢复，源码hash一致、旧metadata另存，无新增删除。monitor加入30秒内存/磁盘采样；finalizer只认本次完整交付。科学目标尚未验收，详见组合审计的资源恢复节。

资源监控原生复现本机psutil.disk_usage的SystemError(impossible bad format char)，标准库shutil.disk_usage正常；监控仅替换磁盘读取接口，完整snapshot重放通过后重启。科学worker未中断、源码和预算未改。

2026-10-03 20:52用户集中CTLA4：核验command后停止旧两队列及其worker，关闭账本不退款；其他病例不排队。CTLA4校准已完成当前交付，直接恢复production，冻结科学源码/配置不变；现有工具并行保留、CTLA4 worker允许全部逻辑CPU/AboveNormal优先级。外部launcher、monitor、验收/浏览器仅CTLA4，附K=8与实际项数，旧报告/metadata保留。cal累计13472.943秒，prod停止时2415.829秒，原14400/7200上限不变；目标尚未验收。用户允许精简非必要要求：优先hard目标阶段性结果，不以搜索未闭合/APBS未解析否认已完整验证的hard证据，但不将阶段性结果冒称全流程验收，不放松逐位/稳定性/null/4of5资格。

2026-10-04续验：45候选/145实际pair，首次正式局部双点D78A/M122W，下界17.6608%，目标仍未达；shadow五组62.8345%单独报告。2375替换、45完整序列及145实际映射/DDG独立审计PASS，隔离v10相关138项PASS。累计140.46分钟/峰值2.02 GiB。SPURS11份实际合同权重一致，补录独立资产pin但不改评分缓存身份；FoldX新断言误读旧T0缺失字段启动失败6.453秒保留计费，以44份已记录模型工具指纹补录独立T0_frozen_tools读取，保留历史清单。非空未达标面板补齐APBS三个未请求/未评估状态。后续53节点相同请求v11继续验证，浏览器关闭失败尚不能记PASS，正式80%与最终可移植包仍未通过。

WSL后续故障：wslpath超时导致SPURS缺评分/补充菜单拒绝；诊断必须区分环境故障与科学排除，不丢历史模型。两失败照计后账本8680.1406秒。定向终止Ubuntu成功但启动/bin/true仍超时，重启WslService被OS权限拒绝；不得继续盲重试消耗预算。Chrome真实83 MB JSON的CDP注入超过100 MB缓冲，改本地只读HTTP后45节点浏览器PASS（写0、异常0）。checkpoint重复proxy near_pairs/距离矩阵改完整保存artifact的hash/路径/count引用，仅新配置；独立恢复、拒绝错hash和科学字段不变回归必须通过。需要管理员恢复运行环境及继续本路线剩余预算，科学验收未达。

当前交付核验完成：隔离v13相关139项/10.15秒通过；实际checkpoint诊断投影41,716,435→2,030,335字节且科学摘要SHA不变。当前完整可移植包184,327,815字节，物理145pair/405field/26null、1838内部hash、另27保存/预筛模型及源码阶段包独立核验通过，897内容对象hash通过；数量可能共享历史身份，不累加。归档与审核计入后账本8714.6406秒，当前拥有2,518,288,511字节。包内artifact一致性PASS，但当前阶段/科学验收INCOMPLETE、53节点请求待WSL恢复，正式1项17.6608%，APBS未请求。验收矩阵见独立任务acceptance_audit.md。原checkout一轮56项回归受当前WSL无响应拖延，由精确已核验PID取消，不能记通过；隔离139项才是完成回归证据。仍需管理员恢复WslService后沿同一账本/证据库继续，未关闭科学目标。

WSL第二次目标续验20秒/bin/true仍超时，没有新增求解/重置账本。保存证据逐候选诊断candidate_failure_diagnosis.csv/json确认45组稳定通过17、模型分歧11、组合拒绝17；5组逐位改变成立、正式1。4组阳性未正式的共同缺口是局部校准未闭合且hard上界未闭合，其中V81W/D100F5对可审查统计追加，其余V81W/D99V、E59R/M91W、E59T/M91W20对已达cohort上限，不能机械扩展。诊断四个源hash固定、未知不填0；科学状态未变，不以shadow替代。当前管理员权限恢复WSL服务仍为下一科学批次前提。

后续自动续验WSL /bin/true实际返回0，环境已恢复；沿原8714.6406秒账本和v13隔离源码恢复同一53节点请求（T6_resume_45_4_calibration），没有重置4小时路线。替换表实测漏斗125映射位点/2375替换，63位点831替换合格，初始最多3项菜单实际184（四位点仅1–2项），加5项盲区补充共189；本批局部双点轻量枚举6514，池容量53，不全量FoldX。资产pin和原资格正常通过，实际求解阶段正在运行，最终科学结果以新生产检查点为准，不能沿用上一轮blocked状态或提前宣称达标。

恢复53节点校准664.031秒、生产291.171秒，账本9669.8266秒；53候选、165实际pair、稳定通过18、分歧13、组合拒绝22，正式仍1/17.6608%，shadow70.1543%只作对照。新轻量checkpoint在实际pipeline完成导出，APBS未请求三个状态明确。发现池被53历史白名单占满，后续会误早停；修复保存菜单/完整节点图纯搜索补池，排除全部历史已筛选身份但保留模型和轻量历史审查。下一阶段白名单保留所有实际物理pair候选再加新双点，不让无物理证据的既往稳定性拒绝项持续挤占活跃池；跨阶段已筛选集合取所有完成生产请求并集，防止重新调度旧失败组合。新纯逻辑回归验证有界补池、候选身份不重复、完整序列保护及真实全池耗尽与当前池耗尽的区别。初次相关回归55过/2失败是C盘可用15 GiB不足测试默认20 GiB预留，非科学缺陷；改使用任务所属E盘独立basetemp复测，不修改预算门槛或删除旧目录。

补池v16隔离58项回归通过8.14秒（测试临时目录E盘任务所有，C盘20 GiB预留失败未放宽）。真实保存菜单/完整图重放生成32个与53历史身份不重复的新双点，FoldX/SPURS新增调用0、165pair manifest SHA不变。兼容早期pilot缺candidate_ids，以实际candidate_site_evidence识别已筛选身份；跨生产阶段并集保留历史，不重新调度已拒绝组合。53节点完整序列/编号/DDG/菜单独立核验PASS，实际浏览器只读PASS。T6_resume_53_1校准沿同一总账继续，阶段白名单18个已测物理候选+8未测双点=26、硬池上限32，53条历史筛选审查单独保存；不删除旧模型/资格失败，重新审核所有保留物理候选。新checkpoint实际导出使用完整诊断文件hash路径引用，旧包是53阶段之前的历史交付，后续完整包尚待更新。

2026-10-04最终部分面板APBS调度修复：原engine仅在4项/80%目标成立后执行APBS，导致到预算终点的部分合格面板无法满足最终审计合同。增加显式final_panel_audit，仅改变最终阶段求解调度；不改变稳定性、五对/同四次、局部校准、hard分母或80%资格。普通中间阶段仍不提前求解。协调器最终审计复用所有实际pair身份、schedule=[5]且extension IDs为空，保留至少300秒归档审查；零新pair独立身份审查通过才可交付。APBS准备或求解失败明确INCOMPLETE，不能声称已审计；已成立hard证据保留并与全流程资格区分。原checkout相关100项回归通过5.71秒，含低覆盖即使APBS完成也仍构建不足、空面板不求解、最终阶段不新增模型/扩null。真实APBS尚待当前批次结束及隔离源码冻结验证。新增只读累计候选漏斗脚本区分历史预筛状态与当前≥5对证据，缺失证据不填0。

2026-10-04用户合同修订：取消4小时累计科学时间上限，不重置已计费时间与失败历史；停止条件改为正式hard≥80%或合格局部位点对穷尽。用户进一步明确不穷尽每位点三个替换的3×3，而是每对先选一个稳定门内机制/化学描述代表，失败或实测盲区才用备用替换；先前“全部合格替换”选择已被本次修订替代。三点仍按补盲/压缩规则生成，不宣称穷尽三点。允许当前磁盘可用空间，保留10 GiB空闲，另预留2 GiB归档峰值。新RunBudget仅显式None表示无限累计总额，JSON写null；旧有限预算保持min语义，修改只在锁定的关闭阶段边界激活并保留原账本快照。每阶段仍有7200秒工程watchdog和原工具超时，失败计费并保存节点，不结束无限总路线。StorageGuard各层及CAS写入统一读取已激活合同，不能仍被默认18/20 GiB拦截。新位点对去重与已计算位点对排除分别记录，未计算的备用替换不能冒充已筛选。160项相关原checkout回归通过10.54秒，隔离冻结及真实新路线随后验证。T6_resume_61_1校准978.969秒/生产278.375秒完成，累计64筛选、23物理组合/220 pair；正式1项17.6608%，shadow7组70.1594%，账本12434.2006秒。新合同未改变任何科学资格、分母、扫描核或旧历史合同。

2026-10-04穷尽协调器v21修复：无限路线若一次多批，旧history只在入口更新，会把刚测位点对的其他替换再排队。改为每个完成生产阶段立即记录累计候选及其原始科学状态，失败阶段不能标已处理；每次补池用已处理位点对排除，并冻结全体位点对清单，空池只有逐项差集为空才称EXHAUSTED。保存累计原始paired索引供恢复和归档，活跃评分池最多248已有节点+8新节点；保留全部正式/逐位阳性/可改善统计及三点精确对照，其他原始模型不删除且资格不继承。超过必要节点上限会明确停止，不能丢掉必要证据。三点使用活跃实测池而非整个历史失败白名单；远距最多2个不同确认身份，失败重试不清计数；所有阶段尝试检查request/log等半写历史，不能只有报告缺失便覆写。纯菜单消融沿用当前账本有效上限，防止无限账本被脚本悄悄改回14400秒。新增累计索引、同字节不同run、变更身份拒绝、批次间历史、冻结scope及半写重试回归；原checkout165项相关回归通过10.09秒。实际v20新合同批次T6_resume_64_1校准1074.953秒、生产331.094秒，累计72筛选（65个不同双点位点对）、27物理组/240 pair，正式仍1/17.6608%，shadow7/70.1594%。728位点对余663；累计13840.9506秒、原14400秒上限已取消，峰值3395983402字节。下一隔离v21批次将核验新增协调器及scope冻结，正式目标仍未完成。

2026-10-04实际v21完成：T6_resume_72_1校准856.046秒/生产387.907秒；累计80候选、73/728位点对、27物理组/255实际pair，新增15对仅为已有候选校准，8个新代表均上游组合筛选未通过。正式仍1/17.6608%，shadow7/70.1594%；账本15085.56秒已真实跨越原4小时且未停止。累计索引与逐批history更新在真实路径通过。独立72检查点2375菜单/125位点/31活跃候选/240pair/DDG/正式mask PASS，实际浏览器只读PASS；31活跃与72累计不可混用。v22按用户确认的失败备用规则，最多将一批8组中的1组用于已实际组合稳定性预筛失败且尚未试备用的位点对，从当前合法菜单选择稳定性预测优先的一个具体替换，其余优先新位点对；不穷尽3×3、不保证备用实测稳定。来源和父候选进入配置、history及候选search_role，资格全量重审。补充阶段工程状态审核：物理Fatal INCOMPLETE或无当前summary不能只因pipeline返回0就标COMPLETED或加入穷尽历史，校准不足等明确科学未知仍可工程完成。原checkout167项相关回归通过10.60秒，包含备用唯一选择/最多一次与Fatal不计完成；隔离v22真实验证随后。此轮未改变物理数值协议、NaN/valid mask、稳定门、五对/同四次、分母或扫描半径。

2026-10-04实际v22批次完成：T6_resume_80_1校准789.437秒、生产418.391秒；累计88具体候选、80/728合格局部位点对（余648）、28物理组/260实际pair。正式仍1项17.6608%，shadow7组仅作对照，不能计入80%验收。组合稳定性累计28通过、20分歧、40拒绝；新备用V40I/M134I因真实预筛拒绝保留失败身份，不能标合格。当前总账16294.074642秒，关闭阶段后沿同一无限总账启动连续resume-batches（T6_resume_88_1起），逐批更新历史与累计索引，不需每批人工确认。v22隔离167项相关回归9.86秒通过，当前新增5对模型身份、峰值任务所属3598857007字节。后续科学结果及完整最终包尚待运行，不以工程通过宣称目标完成。

2026-10-04连续运行期间补充隔离v22接口/序列回归：从source父目录使用项目venv、-p no:cacheprovider和独立E盘basetemp，test_api.py、test_sequences.py、test_mapping_design.py共40项通过5.69秒；API实际HTTP与AF3恢复使用测试替身，不新增FoldX/SPURS/APBS科学调用。日志T1_v22_api_sequences_regression.log；与167项回归不可直接相加，亦不证明当前CTLA4科学目标。当前T6_resume_88_1校准进程仍存活并持续更新局部场文件，未重启或重置账本。

2026-10-04连续运行性能观察：当前父协调器CPU约505秒、科学worker约890秒；owned_bytes三次真实完整遍历分别3.1648/3.8323/3.7336秒，CPU3.0/3.7969/3.6563秒，测量对象约3.68GB，见live_storage_scan_overhead.json。StorageGuard已有2秒缓存，launch虽每250ms轮询但不是每次遍历；完整遍历仍在监督进程反复占用约一个CPU核，不能将其算成物理算法开销或据此声称物理提速。该项为待阶段边界优化：保持写入前检查、每次实时disk_usage/安全余量及证据身份不变，可将外层owned inventory监督频率分开控制。当前v22科学进程与冻结源码不改，不因观察到存活/无新summary而重启。

2026-10-04存储监督优化（当前检出待下次关闭阶段冻结，运行中的v22不变）：case_resources.StorageGuard增加独立inventory_interval（默认2秒、最高30秒）与force_inventory；每次仍实时disk_usage核验，非零reserved写入前仍完整统计。无限路线的外层监督改30秒inventory，阶段结束强制重新统计；旧有限路线仍2秒。EvidenceStore每个检查点复用同一guard，immutable hardlink零新增字节可共享两秒缓存，mutable副本每次非零reserved仍重新统计。三项回归验证缓存时磁盘空间下降仍拒绝、强制/写入前统计、不可变inode/字节及可变副本独立inode/每次新统计；45项针对性通过3.24秒，相关170项通过14.01秒。真实目录只读三次admission比较：旧每次新guard11.3048秒CPU11.0313，新共享guard3.3959秒CPU3.3594；只支持存储准入开销下降，不能称物理算法提速。报告storage_admission_reuse_benchmark_v23.json，源码SHA ba5c434e4fa2db1dcbf488cbb293f556530b07678cb315e4f4fd61d8ac2396ce。benchmark首次因sha256_file传str而未保存报告，修为Path后复测，失败未改变科学结果。实际v22批次仍在计算，正式17.6608%尚未达标；根目录源码新修改不能冒充正在运行的冻结版本或已完成最终包。

2026-10-04已设置协调器安全切换监视器（实际切换尚未发生）：scripts/restart_coordinator_at_boundary.py仅核验明确PID/启动身份/当前case及--resume-batches命令；生产报告完成后短暂suspend协调器以关闭spawn竞争，复核没有科学子进程、最新阶段为production、case JSON inactive且原生账本锁可获取，再仅终止协调器。检查不成立立即resume并等待；不停止任何科学worker、不写/reset账本。6项独立回归1.82秒通过，覆盖新cal报告阻止旧prod冒充边界、活worker、active账本、JSONclosed但锁仍持有、真实关闭时仅终止协调器且账本字节不变、无关PID拒绝。监视器原early-continue路径补finally内sleep避免错过边界后频繁suspend；旧监视器在目标生产报告尚不存在的前提下仅终止其自身后重新启动。当前科学协调器PID21828、worker24396未停止，v22仍计算；边界结果将写coordinator_boundary_restart.json。下一步在实际PASS/原进程终态后保留v22源码、冻结v23并沿原账本恢复；不因watcher在等待便宣称已启用新版。

2026-10-04待冻结v23复核：已修订路线的StorageGuard默认inventory_interval=30秒，未修订legacy默认仍2秒，显式配置2秒仍有效。该缓存只作用于零reserved查询；非零新增字节准入仍立即完整统计，每次均查实时disk_usage，最终阶段force_inventory重新统计。新增worker默认传播回归后，菜单/物理证据/Joint V4/组合合同/局部构建/null schedule/边界控制综合177项通过16.95秒，日志T1_resource_worker_v23.log；当前运行的v22仍未切换。已对正在使用的冻结v22源码只读打包，276源码+1overlay manifest共277对象逐一核验，source_snapshot_v22.zip为3495645字节、SHA d148da8c28e49515dc210206d0098872c6779f9075c207d326acdac35181710f；source_snapshot_v22_audit.json PASS，科学新增调用0。新final provenance明确纳入case/source_snapshot\_*.zip及其audit JSON，防止切换时只保留旧源码hash而丢失实际字节。边界监视器与科学任务实际handle均仍存活，生产终态/闭账条件尚未满足；不会提前冻结或宣称正式80%通过。

2026-10-04实际v22闭合阶段与安全切换：T6_resume_88_1校准3382.453秒、生产627.344秒均工程COMPLETED；累计96具体候选、87/728位点对（余641）、32物理组/295合法pair身份。A41Y/A61R从5补至20；正式仍1项D78A/M122W、17.660761%，shadow8组79.381610%不计正式目标。WT可靠分母6874.363277945836平方埃保持不变。总账20304.558642秒、inactive、无限上限null，阶段峰值3989584365字节。监视器PASS复核PID21828创建身份、最新production报告、无科学子进程和原生账本锁后仅终止协调器；终止科学worker=0，账本SHA 8947eb540abf9c2cdaa53d580bc21291f914f320e60a8e211eff7c1992f8f77c 保持。外层进程exit1是此明确协调器切换，不是科学失败；实际两个阶段COMPLETED报告保留。关闭后用原v22函数补写累计pair索引/已完成history，防止切换窗口在报告之后、history之前而漏记或重排旧位点对。后续冻结v23仅5个已知文件修改及2个新增控制/回归模块；科学数值、候选资格及所有稳定/五对/同四次门槛未改。仍需正式80%或约定scope逐项穷尽及最终APBS/可移植包核验。

2026-10-04隔离v23门通过：278源码文件hash一致，冻结源综合177项通过16.73秒。独立96检查点语义审查PASS：2375替换/125位置、36活跃候选（不是96累计）、295实际pair身份、原作者到规范编号/WT字母/DDG及正式mask、固定可靠分母6874.363277945836平方埃均一致。累计诊断96候选稳定31通过/21分歧/44拒绝，32物理组保留原始身份但不等同当前稳定合格组；逐位同四次9组，shadow8组79.381610%，正式1项17.660761%。输出independent_case_semantic_checkpoint_96.json、T6_96_funnel_audit.log；APBS中间状态仍未请求/未评估。沿20304.558642秒原账本恢复v23一批真实确认（--max-new-batches 1），准备以实测验证工程监督改动，不将存储benchmark声称为物理算法提速。新版实际T6_resume_96阶段结果待计算，源码冻结后本节只追加根目录连续性文档，不改活跃冻结字节。

2026-10-04隔离v23真实批次通过：T6_resume_96_1校准2176.531秒、生产616.625秒均工程COMPLETED，累计104具体筛选候选、94/728合格局部位点对（余634）。关闭检查点累计收费23098.715642148974秒、阶段峰值4237065026字节；正式仍1项D78A/M122W、17.660761%，shadow8项79.284987%仅作对照，未达到正式80%。保持原可靠WT分母6874.363277945836平方埃、稳定门及五对同四次要求。已沿冻结v23及同一总账启动不设max-new-batches的连续位点对优先路线，当前T6_resume_104_1_calibration正在计算；未将当前池耗尽误称全scope穷尽，未更新历史最终包冒充当前交付。存储使用已激活当前可用空间合同，至少10 GiB磁盘安全余量及2 GiB交付预留；此次真实阶段确认工程监督改动可用，不能据单批耗时宣称物理算法提速。

2026-10-04位点对已完成scope独立核验：case/audits/audit_completed_pair_frontier.py逐项核验104累计候选身份、位点唯一性、last_completed_stage对应真实COMPLETED生产报告及其SHA、3个显式备用替换的同位点父候选/实际失败来源；与冻结728对逐项差集得到94已处理/634未处理，超scope双点=0。94对中85只测1个具体替换、9保留2个历史替换；不将历史额外替换冒充穷尽3×3。completed_pair_frontier_104.json工程PASS、搜索INCOMPLETE，科学新增调用0；审查脚本字节及SHA已单独收据保存，最终通用case JSON归档可携带。当前T6_resume_104_1_calibration仍在生成场证据，不能将未完成阶段计入已处理或科学目标。

2026-10-04连续v23下一批完成：T6_resume_104_1校准1415.797秒、生产643.14秒均工程COMPLETED；累计112具体候选、101/728位点对（余627）、35物理组/340合法pair身份。累计收费25158.637642148977秒，阶段峰值4401776268字节。正式仍1项/17.660761%，shadow8项79.284987%；下一T6_resume_112_1_calibration已自动启动，同一handle23942，不重置账本。独立完成范围审查completed_pair_frontier_112.json PASS：全部完成生产报告身份/位点唯一性/冻结差集，无超scope双点。

新增只读完成阶段调用收据：T6_resume_104_1_calibration_v4_manifest.json及production对应manifest原字节分别SHA b97f66e8a0379de590ddc626414f8dd1d12d7c2f34bb81a8d541d644253a5435、3424932c41dbde105a5e72378cd17ecb5f5cc91d0de9f78555be45f5c941586b；核验导出mtime处于该请求与result之间及完整白名单相等，不把后续阶段manifest错配。初始FoldX37候选中4个各新增1对，生产初始新增全0；后续仅2组进入五对实际物理库，其他初始预筛模型不冒充合格；初始counts不包括补充调用，不能据此宣称总FoldX subprocess次数。SPURS single缓存命中、multi未命中；cal/production物理solver调用分别NanoShaper21/1，APBS和PDB2PQR0，缓存hit26/4。较早阶段未保存相应worker运行收据的缺口不得追溯虚构。

独立语义快照核验independent_case_semantic_checkpoint_112.json PASS：从完成生产复制9类语义输入、前后源字节不变且summary与报告相等后审查，实际41活跃候选（不混同112累计）、2375替换/125位置、340模型身份/WT字母/作者到规范编号/DDG、正式mask及固定可靠6874.363277945836平方埃分母均通过。快照位于case/audits/semantic_checkpoint_112，不改运行中输出，科学新增0。

独立全替换scope核验full_replacement_pair_scope_audit.json PASS：从完整19项替换表中的831合格替换、63位点及完整保存图做24埃局部最短路径，再对每对寻找满足完整序列保护的实际替换存在见证；得到728对，与冻结scope逐项相同，漏对0/无见证多余对0。该证据证明初始3菜单没有缩小本case位点对穷尽范围；存在见证不等同实测稳定或扰动资格。审查脚本字节/SHA及报告以case JSON收据保存，可进入最终通用provenance；当前正式80%与最终APBS/可移植当前包仍未完成。

2026-10-04连续路线T6_resume_112_1校准工程COMPLETED：1747.641秒，累计收费26907.184642148975秒，新增222729718字节、阶段峰值4629721387字节；累计原始物理库38组/355合法pair身份，较前完成批次增加3组各5对。校准科学状态CALIBRATION_INSUFFICIENT，不能升级正式资格。已只读保存对应v4_manifest原字节及worker receipt，SHA bcace438d1e9a645a1f5ab9c4d44bdff21ec7723d0f7b8336dec068bb78ce412；初始FoldX3个候选各新1对（后续确认另计），NanoShaper31、APBS/PDB2PQR0、缓存hit32。4个独立审查脚本收据已检查保存source_utf8编码字节与文件/SHA完全一致；科学新增调用0。当前生产审查已自动启动，同一handle23942；累计完成仍112具体候选、101/728位点对，正式17.660761%/shadow79.284987%以已完成T6_resume_104_1_production为准，本批未提前计完成。不改运行中的冻结v23源码或原科学门。

2026-10-04连续v23完成T6_resume_112_1：cal1747.641秒、production541.203秒均工程COMPLETED；累计120具体候选、108/728位点对（余620），38实际物理组/355合法pair身份。累计收费27448.387642148977秒，阶段峰值4640933514字节。正式仍1项D78A/M122W、17.660761%；shadow9项84.828267%首次超过80%，不能升级原正式合同。生产初始FoldX新增model_counts全0，完整worker收据/v4_manifest保存，SHA f90a5b803d1f622df0ba921cd3639750b639521f2ec4974fbd2667357b90b888；物理NanoShaper1、APBS/PDB2PQR0、cache_hits4。

completed_pair_frontier_120.json独立完成报告/冻结scope审查PASS：108真实完成对/620剩余、无超scope对。independent_case_semantic_checkpoint_120.json从完成生产冻结9类语义输入后审查PASS：43活跃（不混同120累计）、2375替换/125位置、355模型身份/序列/WT/作者到规范化编号/DDG、正式mask及6874.363277945836平方埃固定分母均正确。shadow_curve_independent_checkpoint_120.json从该完成阶段的实际shadow_positive_matrix.npz逐K独立OR复算，5候选82.959455%、9候选84.828267%，与每行导出差≤1e-12；矩阵area与冻结canonical逐项相等、dtype bool、全部选入候选同四次逐位证据complete；原始压缩矩阵字节base64/SHA随独立JSON报告保存，便于可移植复核。shadow与正式状态始终分开，科学新增调用0。

本次新增shadow候选CV2_ca1bf7d4af3aaa56及CV2_e8cabda1d5453f12（并非净增2个面板位；重选排掉1个冗余/当前资格变化候选）均measurement_local_complete=true，但全通道local complete=false、hard_scan_bound_resolved=false，未闭合上界56.1254575/58.9022829平方埃。校准unknown记录均为QUANTILE_UNSTABLE、appendable=true，未观察到确定性拓扑、缺失查询、质子化或方向未知。已自动排入当前T6_resume_120_1_calibration的两个10/20扩展白名单；不承诺追加一定获得正式资格。当前handle23942真实存活，新批继续运行，正式80%/最终APBS/当前完整交付包仍未完成。

2026-10-04 v24收尾与审查修复：T6_resume_120_1校准1934.25秒/生产512.547秒完成，128具体候选、115/728位点对；T6_resume_128_1校准2252.406秒/生产578.297秒完成，136具体候选、122/728位点对（余606），44原始物理组/445合法pair身份。累计收费32727.247642148977秒、峰值5460749997字节。正式仍1项D78A/M122W、17.660761%，shadow9项84.834436%，可靠WT分母6874.363277945836平方埃不变。累计pair身份不等同新FoldX调用次数；两阶段初始worker manifest分别SHA45d1a01ab83eae1de322c2860b0f9971537aef237053ff8de102d402d24d9143/4c730616b56559fabb09ae928d89d55d45dbd20f41facba194df35dc68c029e6，production初始新模型数全0。completed_pair_frontier_136.json独立PASS、无超scope对。两个旧shadow新增候选达到20模型仍未闭合正式上界；下批没有重复扩展这两个cohort。

审查纠错：independent_case_semantic_checkpoint_128.json原审查读取下一校准阶段临时1对零值mask，却保留前一生产summary，输出official_hard_fraction=0仍误标工程PASS。旧audit只在ACHIEVED时检查mask与summary，因此前后字节不变检查不足以防跨阶段混读。本次保留原文件/源hash，semantic_checkpoint_128_invalidation.json明确使其正式mask一致性结论无效；此前112/120快照的实际mask与17.660761%摘要一致。audit_ctla4_case新增verify_hard_curve，在未达目标时也逐K检查mask覆盖和实际项数，反例已通过。关闭生产的当前mask1745个点对应17.660761%，科学结果没有被零值snapshot改写。shadow128独立矩阵审查未涉及该临时mask，仍仅为shadow证据。

源码修改仅控制与独立审核：位点对穷尽后新增按规则的局部三点补盲审查，然后最多2个远距备选；正式达标时压缩只从正式面板父双点提案。已完成production的三点按具体基因型与角色去重，旧合同可由实际父双点恢复身份；失败/仅calibration不计完成。每个三点完成生产后立即保存pair索引与历史，精确三个双点对照继续强制生成。压缩首次接受后停止跨阶段继续压缩，避免把原累计1个百分点损失合同重置到新基准；4项面板不执行无构建数收益的压缩。全局几何、资格、稳定性、五对/同四次和数值协议未修改。新相关回归83项通过3.31秒（先前测试末尾误置的NameError已修正复测，失败不作通过）。该收尾路线纯逻辑门已过，真实后期三点仍待位点对穷尽或正式达标后验证，不称全流程验收。

安全切换：coordinator_boundary_restart.json确认在T6_resume_128_1_production真实完成、账本inactive、无科学子进程、原生锁可得时，仅终止owned PID42600，终止科学worker=0；账本SHA4d623c44987926cc32fc231960e03edb3ec846601178ffd2adb555c33a958aa9保持。旧v23实际源码逐文件核验/归档后冻结v24，沿相同总账、证据库、728对合同继续，无退款/历史删除。最终APBS、当前完整包、只读浏览器及正式80%尚未验收。

2026-10-04完成检查点144：T6_resume_136_1_production及对应校准均工程COMPLETED，cal1433.922秒/production599.609秒，累计144具体候选、129/728局部位点对。累计收费34761.46764214898秒、峰值5645101499字节。独立范围、完整替换表/序列/编号/DDG/475实际pair身份及逐K正式mask一致性审查PASS；正式1项/17.660761%，固定WT分母6874.363277945836平方埃。shadow9项/84.834436%独立矩阵OR复算PASS，不能替代正式资格。production工具计数{'nanoshaper': 1, 'apbs': 0, 'pdb2pqr': 0, 'cache_hits': 4}，初始FoldX模型counts只代表初始worker，不冒充补充调用总数。证据分别为independent_case_semantic_checkpoint_144.json、completed_pair_frontier_144.json、shadow_curve_independent_checkpoint_144.json及对应阶段worker收据。保持当前冻结源码、原科学门、真实run身份与统一账本；正式目标和最终APBS/当前可移植包仍待完成。

2026-10-04完成检查点152：T6_resume_144_1_production及对应校准均工程COMPLETED，cal669.782秒/production795.047秒，累计152具体候选、136/728局部位点对。累计收费36226.99864214898秒、峰值5664972430字节。独立范围、完整替换表/序列/编号/DDG/475实际pair身份及逐K正式mask一致性审查PASS；正式1项/17.660761%，固定WT分母6874.363277945836平方埃。shadow9项/84.834436%独立矩阵OR复算PASS，不能替代正式资格。production工具计数{'nanoshaper': 1, 'apbs': 0, 'pdb2pqr': 0, 'cache_hits': 4}，初始FoldX模型counts只代表初始worker，不冒充补充调用总数。证据分别为independent_case_semantic_checkpoint_152.json、completed_pair_frontier_152.json、shadow_curve_independent_checkpoint_152.json及对应阶段worker收据。保持当前冻结源码、原科学门、真实run身份与统一账本；正式目标和最终APBS/当前可移植包仍待完成。

2026-10-04双点/三点稳定性及物理证据研究（仅审查，未改活跃合同）：新增docs/ctla4_stability_and_physical_evidence_research_20261004.md。固定96检查点的17个双方门内MODEL_DISAGREEMENT逐项核对SPURS multi缓存/权重/规范组合、FoldX模型/Dif/Raw/PdbList/突变列表hash；均SPURS近中性、FoldX更强稳定，差值中位数2.087601，17组全部引入Arg、14组兼有Trp，16组仅1对预筛、1组5对，不能直接升级正式。Raw与Dif一致、符号/组合未见错配；规范原始与Repair输入Cα完全相同，但931公共原子中282个位移>0.05埃。主要涉及joint_design_v4/engine.py分歧优先状态机及strict FoldX名单，建议primary SPURS双点风险分层与FoldX结构辅助审查，绝对能量损失须结合WT裕量和按阶数实验误差校准；三点输出仍外推，不混用独立single与multi之差作为精确epistasis-head。

第3–6项新证据：保存三个同原子scale2/3 closed mesh的384个固定查询均得到有限SDF，Euler不匹配而距离差中位数约0.01埃；支持局部逐通道QC，不证明全部局部场合格。两个20-model null cohort（T104W/Q111R、A41Y/L47R）各只有10种WT/Mut字节内容、首五pair只有4种；合法run身份保留，不能解释为20独立采样。固定当前20模型scope/error做前4/9/19块重放，A41Y/L47R shape校准有效面积97.7716%/100%/53.4778%；新噪声信息可以正确降低资格，不允许挑旧前缀通过。建议区间式同四次阳性下界、阳性资格与全范围解析分开、模型级依赖/多通道OR/12埃扫描核误覆盖校准；不把当前shadow直接升为正式。case/audits/scientific_methods_20261004/analyze_saved_evidence.py保存重放全部断言PASS，93份输入快照hash一致，新增外部求解0、固定mesh查询384。研究附件research_data.json/stability_17.csv/replay.log及sources保存；未新跑FoldX/SPURS、未建立独立null误报或CTLA4实验安全阈值，建议仍待新合同和逐级验证。完成152检查点仍正式17.660761%、shadow84.834436%，分别报告。

2026-10-04完成检查点160：T6_resume_152_1_production及对应校准均工程COMPLETED，cal2712.547秒/production717.109秒，累计160具体候选、143/728局部位点对。累计收费39657.65564214897秒、峰值5995817091字节。独立范围、完整替换表/序列/编号/DDG/495实际pair身份及逐K正式mask一致性审查PASS；正式1项/17.660761%，固定WT分母6874.363277945836平方埃。shadow10项/87.043152%独立矩阵OR复算PASS，不能替代正式资格。production工具计数{'nanoshaper': 1, 'apbs': 0, 'pdb2pqr': 0, 'cache_hits': 4}，初始FoldX模型counts只代表初始worker，不冒充补充调用总数。证据分别为independent_case_semantic_checkpoint_160.json、completed_pair_frontier_160.json、shadow_curve_independent_checkpoint_160.json及对应阶段worker收据。保持当前冻结源码、原科学门、真实run身份与统一账本；正式目标和最终APBS/当前可移植包仍待完成。

2026-10-04待冻结v25最终复用审查修复：final_panel_audit原身份集合仅包含WT/Mut hash、group、run，遗漏batch_path；不同真实批次输出相同字节时可能误报零新增模型。新增paired_execution_identities保留batch/group/run执行身份并拒绝缺批次与重复执行；最终APBS前先审核保存身份，后逐候选比较完整执行与模型字节身份。新增相同字节不同批次、重复UID、最终阶段新增批次拒绝/不变接受反例；旧测试fixture补充实际模型具有的batch_path。相关86项回归3.73秒通过；首次扩展回归85通过/1缺batch fixture失败已修复复测，不计作通过。当前科学计算仍使用v24，资格/数值/累计账本不变；v25须在关闭生产边界冻结后才生效，真实最终APBS与完整交付未验收。

2026-10-04完成检查点168：T6_resume_160_1_production及对应校准均工程COMPLETED，cal2080.328秒/production761.125秒，累计168具体候选、150/728局部位点对。累计收费42499.77964214898秒、峰值6256289166字节。独立范围、完整替换表/序列/编号/DDG/530实际pair身份及逐K正式mask一致性审查PASS；正式1项/17.660761%，固定WT分母6874.363277945836平方埃。shadow10项/87.043152%独立矩阵OR复算PASS，不能替代正式资格。production工具计数{'nanoshaper': 1, 'apbs': 0, 'pdb2pqr': 0, 'cache_hits': 4}，初始FoldX模型counts只代表初始worker，不冒充补充调用总数。证据分别为independent_case_semantic_checkpoint_168.json、completed_pair_frontier_168.json、shadow_curve_independent_checkpoint_168.json及对应阶段worker收据。保持当前冻结源码、原科学门、真实run身份与统一账本；正式目标和最终APBS/当前可移植包仍待完成。

2026-10-04真实保存模型补证：semantic_checkpoint_168固定paired_manifest SHA 9c2a9514c79f75fa4fe9c03059549c7f9559fe3ef819d0f4ef70135d4801de36，49物理候选/530合法batch-group-run身份全部通过新严格审核；原WT/Mut hash+group+run集合实际仅315身份，49候选发生身份合并。v25_saved_execution_identity_audit_168.json保存逐候选差异；这是独立重复执行被旧审查集合合并的实证，不代表模型缺失、重复计算可删除或正式覆盖改变。新增科学调用0。v25仍待安全关闭生产边界切换。

2026-10-04已授权coverage v3实现及真实T3：详见 [ctla4\_coverage\_v3\_implementation\_audit\_20261004.md](ctla4_coverage_v3_implementation_audit_20261004.md)。独立v3 profile、SPURS主要双点风险分层、时间播种固定生产槽、10训练/10留出专用WT、同四次阳性下界、局部Euler诊断、三态几何信息/失效鲁棒性和三点损失审核已实现；旧shadow不升级。176相关回归通过。保存96审查恢复17+39确认机会，4真实冲突保留。真实D78A/M122W五对生产+20 null工程COMPLETED，正式1项扫描17.6170307579%、直接2.8816261241%，固定分母6874.363277945836；未达80%。收费1154.594秒包含旧检查点额外开销，统一总账不重置。发现逐文件存储预留反复全盘扫描的平方级开销，改批次预留、限域检查点、独立hash核验，待第二冻结快照复测。20 NMR作者链A129公共CA核验通过；密图代理差异0.669646个百分点，只是敏感性报告。第一真实源码保留，后续T4/T5/full/APBS/包及浏览器仍待验收，不称全流程完成。

2026-10-04v3真实T4闭合：4候选中3正式下界合格，扫描49.6439461279%、直接5.0355311710%、二重18.5601742357%、单构建失效最差32.0269153700%，仍未达80%；80真实pair身份（production20/train30/held30）独立核验PASS。T67R/L141W五模型稳定性排除不计正式。收费1929.75秒，总账46421.66964214897秒保留；峰值6.97GB。218相关回归通过，真实25-pair便携探针151文件hash及浏览器只读DOM通过。源码新增风险注释保留跨度/重叠、五模型失败同步资格、全五模型配准门与同四次局部未知分开、压缩1/2/3码距保护、精确对照先行、终包显式纳入物理/source/audit ZIP及便携原始执行身份。T4关闭边界只终止owned coordinator，科学worker0；将由第三冻结源码继续T5/full，不重算合法缓存、不继承legacy资格。T3/部分null内容重复依赖限制详见v3审计，模型重复不称独立热力学观测。

2026-10-04v3 T5闭合并独立审核PASS：真实T67R/M134R精确对照及A61W/T67R/M134R三点各5 production+20专用null；三点自身与低阶upper比较合格进入Expansion，当前5构建扫描51.2083397385%、直接6.5141603282%、二重32.4715507925%、单失效最差33.5913089806%，未达80%。130实际pair（production30/train50/held50）不混角色；收费1799.875秒，统一累计48221.560642148965秒。安全关闭边界worker0终止，旧阶段保留。下一冻结修复全阶段三点role/父对照绑定、已测max吞吐缩批≤8、真实可行4–12双点80%中途停付余下昂贵确认、四项面板无收益压缩跳过及按实测位点未知备用。正式pool不得因三点自身合格就绕过角色审核。科学核心数值/null协议不变，最终full/APBS/包仍待验收。

2026-10-05v3 T6文件占用恢复：source_v3_4独立221回归通过20.29秒。第一T6尝试因检查点/summary原子替换WinError5而工程失败，已收费1204.687秒，累计49426.263642148966秒未退款。Windows Restart Manager定位外部node.exe36144持有两个JSON，非只读属性；不终止用户进程、不绕过原子性截断文件。占用解除后零持有者、同内容替换0.016秒通过；12项JSON/阶段保护回归3.26秒通过含真实Windows占用。已存Q117R/E155W及L47W/Q117R局部下界合格，暂存面板并集79.7089343475%仍未达目标；pending-stage按原固定槽恢复，原模型/物理/null/失败tmp和阶段源码不删除或升级。持久外部占用仍为工程发布限制；不能按科学排除清空证据。详见v3审计与T6_file_occupancy_recovery.json；full/APBS/最终包仍待验收。

2026-10-05v3首次达标：V3_T6_6_2闭合六个双点87.8830909626%扫描/13.9969181247%直接，下界包含/角色审核通过；二重40.6391445651%、单失效最差76.9530388452%、最大盲区300.6027476329 Å²。真实达标见证停止后续昂贵双点确认，未冒称728穷尽。六项APBS执行全部完成，PDB2PQR在作者A56 Asn交换OD1/ND2坐标约2.1–2.3 Å，准备和数值解析保留未知，不计电静力扩展。压缩时发现legacy bounded_retained_ids复活无关暂缓S55F/S101R，以及engine_v3未使用声明的固定compression_baseline_ids。已启动额外计算保留并收费；v3专用白名单及固定4–12合格双点/80%基准校验修复，原始基准不随新池优化重置，压缩角色审核同步。合同frozen_baseline_information_compression_v3_1必须以新快照、同实际三点/原基准重放替代旧未固定审查；34项针对性回归2.60秒通过。科学字段/null/采样/稳定性不变，闭合阶段恢复不得重置总账；最终复审、可移植包仍待验收。

2026-10-05冻结压缩复核闭合：source_v3_5独立239回归18.38秒通过，V3_T6_COMPRESSION_FROZEN_1收费964.891秒；原六双点87.883091%基准保存，实际三点因码距1/2/3及单失效鲁棒性下降拒绝。关闭边界仅终止owned coordinator28000，科学worker0，原账本53366.810642148965秒保持。最终选板需审核全部已付费且五production+10train+10held闭合的阳性双点；不能只保留当时selected而漏掉已测正式池，也不能恢复不完整冷候选。final_panel_audit不触发普通中途80%早停，避免漏审已付费池；source_v3_6独立240项20.76秒通过，新增完整/不完整/阴性/三点调度反例。已请求8个已测双点的FINAL_APBS，不新建基因型，实际最终包/浏览器尚待完成。来源与物理数值合同不变，各阶段源码冻结不改。

2026-10-05最终科学审查闭合：V3_FINAL_APBS_1复用8已测阳性双点，正式5项A61W/M134R、T67R/M134R、Q117R/E155W、S55F/S101R、D78A/M122W。扫描86.0262353051%、直接11.5744840011%、二重32.3636162624%、单失效最差72.3038190447%、盲区376.8339883881 Å²。独立核验235真实pair身份（production55/train90/held90，含未入选/稳定排除历史），正式每项5+10+10，固定分母6874.363277945836及L⊆U/三态/编号/DDG/序列PASS。最终阶段工具新增NanoShaper/APBS/PDB2PQR均0、cachehit107；模型身份保持，不称物理算法提速。独立枚举当前8双点全部70种4构建子集，最高79.6660147108%，证明5是该实测双点池内最少，不证明728对全局最少。目标ACHIEVED、搜索TARGET_REACHED且scope_exhausted=false；APBS执行COMPLETE但准备及数值UNRESOLVED，不补入电静力覆盖。

最终包审核复现一次工程缺陷：full_length_protein_sequences.fasta合法拆分对照PAIRWISE可同原双点主探针同基因型，独立审核器却全局拒绝重复。修复audit_ctla4_coverage_v3.verify_selected_exports，分别核查主探针唯一性、每条canonical WT/编号/完整序列、analysis设计身份和条件对照父子子集；35相关回归3.44秒通过，错误控制序列仍拒绝。package_delivery将已关闭科学阶段的归档恢复与科学求解分开，只恢复打包，paired_manifest hash前后不变，源6原字节和失败费用保留。真实headless Chrome最终DOM核验PASS，直接11.57%及测地扫描86.03%分开呈现。核心语义、便携模型/source及ZIP成员hash审核PASS；阶段/provenance补充另有收据。正式报告routes/coverage_v3/results/CTLA4_v3_final_report.md；NMR骨架、建模重复内容依赖、无全局95%控制及表达/折叠/结合未经实测仍是科学限制。

2026-10-06既往日志全审及支线修复：117日志/64阶段/四包hash审查，原正式5项86.026235%仍通过。当前未冻结APBS旋转分类二元/三元残基键导致IndexError，锚距原表达式又误比较位移；electrostatics统一三元身份/真实锚距，缺锚保留未知、旋转为条件诊断不放宽门。reanalysis_apbs_preparation_fix修复min_only距离又取min导致(1,)零数组，恢复逐10199点/断开inf；沿保存准备/介质/numerical_error审核，撤销无依据全NaN设置并区分原始与可用电势及unknown原因。新增费用对账脚本，历史1106.213642秒单列未归属不造分摊；新launch用保存闭账减carried计费，含准入/检查点，另列wall。audit输出跨角色WT内容重合，蛋白链由entity_poly识别。160相关回归5.99秒通过（首轮锚距失败修复复测），5×5保存PDB/PQR重放、5衰减附件形状及235pair/主ZIP再核验通过，新科学求解0。原包/冻结模型/总账不改；APBS仍准备/介质/数值未知，单骨架/样本依赖/实验未校准不是工程修补可宣称解决。详见[修复审计](prior_task_repairs_20261006.md)，真实补充证据在E:/user/New project/cross_protein_validation_20261006/repairs。

2026-10-07 field-fusion continuation: contact-scale cohort downloads retry transient TLS EOF and `IncompleteRead` into atomic `.part` files. A non-closed SES raises an exact topology error inside the shared mesh helper; the contact-scale builder records that complex as `UNKNOWN_OPEN_SURFACE` and continues, while unrelated runtime/identity errors remain fatal. Regression: contact-scale and cross-protein candidate-freeze scoped suite 11 passed. Frozen author-ID preflight passed for all 7 non-CTLA4 structures. 1AH1 first-model selection remains fixed; its 20 NMR models share 124 mapped Cα residues, with 190 pairwise aligned RMSDs (median 2.63 Å, p90 3.63 Å, max 4.91 Å); these are conformer sensitivity, not independent proteins or thermodynamic replicates. Contact-scale build remains active; last accepted-console checkpoint was 12/30, with 1,704 unique CIFs and 13 completed SES cache objects (one local topology failure). Builder SHA256 matches its frozen source snapshot; E: has about 36.0 GiB free. The 22,548-row source contains mapped candidates but no antigen sequences in the frozen chain annotations, so homolog filtering must parse the coordinate structures. Entity-sequence preclustering cannot be claimed equivalent because deposited sequences may include unresolved/truncated residues and each later homolog remains eligible if earlier geometric QC fails. Cross-protein SPURS candidate scoring waits for this shared WSL/NanoShaper job. No formal CTLA4 result or scientific gate changed.

## 8. 2026-10-07 CTLA4 field-fusion closeout

本轮 CTLA4 field-fusion 研究分支已完成，详细逐项数据与限制见 \[CTLA4完成报告\](E:/user/New project/ctla4_field_fusion_20261006/deliverables/CTLA4_field_fusion_completion_20261007.md)，可复核包为 `E:/user/New project/ctla4_field_fusion_20261006/deliverables/CTLA4_field_fusion_v1_portable_review.zip`。原正式五构建体与 86.026235% hard/12 Å 结果不改；第六构建体 S55F/I102R 在新研究清单中追加，旧hard mask增157.61 Å²/2.293百分点，新三机制raw阳性并集增194.49 Å²/2.829百分点。两种指标定义不同，不合并成同一coverage。

六候选均完成全表面 geometry/chemistry 和 APBS Δφ执行；每候选固定5 production、10 estimate WT、10 held-out WT槽位。PDB2PQR heavy-atom-preserving准备的生产与两组null身份匹配；APBS和融合仍保留大量unknown，尤其电静力，研究状态是“计算完成、科学场部分解析”。12 Å geodesic机会模型在panel5和panel6都饱和，不能把100%扩张机会解释为整个表面直接受突变改变或实验结合下降。15个冻结Fab–抗原结构样本15/15测地尺度解析，形成12.403/13.147/14.827 Å三档结构先验；非实验校准。另核验12个CTLA4 Fab复合物交叠，只有几何对应关系。

修复根因：接触图无向边排序沿错轴曾生成自边，且抽样允许重复；已改为逐边端点排序、无放回抽样及非正距离拒绝。APBS网格合同对候选ID变化过度严格；实算证明第六候选网格与原五项相同后才允许append-only更新，并保存旧合同历史。测试 `test_antibody_contact_scale.py`、`test_field_fusion_research.py`、`test_pdb2pqr_research_protocol.py` 共26项通过；完成包1564成员独立hash核验通过。累计solver ledger为NanoShaper/APBS/PDB2PQR=379/307/148。

原 `study_manifest.json` 保持不变，其 `pair_manifest` 已有六候选各25行，但冗余 `paired_models`摘要仅含原五项；交付收据依逐行pair记录补足六项计数，没有修改历史身份。冻结的八蛋白候选科学评分尚未开始；本分支完成不能被用来宣称跨蛋白普适性。

## 9. 2026-10-08：CTLA4 收尾与跨蛋白暂停

用户要求暂停后续跨蛋白工作，先完成 CTLA4 收尾。CTLA4 field-fusion 报告、便携包与完成清单已交付；六个设计候选为 A61W/M134R、T67R/M134R、Q117R/E155W、S55F/S101R、D78A/M122W、S55F/I102R。原正式五项及 86.026235% 结果未改；第六项仅作为研究扩展候选，旧 hard-mask 面积对照为 88.3190%，不等同正式资格或实验结合证据。

跨蛋白候选冻结进程已停止并保留部分结果。`freeze_manifest.json` 当前列有 CTLA4 固定集及 PD1、PDL1、B7H6、hen_lysozyme 各6项；MBP 与 T4_lysozyme 未完成冻结。该部分只属候选输入准备，没有跨蛋白科学评分、物理场确认或 APBS，也不能据此宣称普适性。后续若恢复，必须从保存清单检查点继续，不重复或覆盖既有输出。

## 10. 2026-10-08：两组 CTLA4 三点压缩复审

用户提出 1+2→A61W/T67R/M134R 与 4+6→S55F/S101R/I102R。完整审查及保存证据路径见 [三点压缩可行性报告](ctla4_triple_compression_feasibility_20261008.md)。前者已有五对生产、各10个估计/留出WT和三个精确双点对照；以当前六项反事实基准重算，扫描覆盖 88.3190%→88.9012%，但删任一构建体后的最差覆盖 75.9758%→57.4396%，几何区分1/2/3位均下降，按原压缩合同拒绝。后者缺三点自身与 S101R/I102R 精确对照；即使其 mask 理想，删除新三点后剩余 1/2/3/5 仅覆盖73.9350%，已低于基准最差75.9758%。两组同时压缩时相应上界74.3720%。这些是基于固定旧V3 hard-mask 的面板必要上界，不把第六项升级为原正式面板、不推断三点生物学表达或结合效果；新增科学求解0。

## 11. 2026-10-08：小分区历史代码审查（未修复）

用户指出1–2残基小区不符合分区预期。最新作业ee5b17b0f7ea有11区，其中VS007/008/009/010分别2/1/2/2残基。Git同图重放定位到01b52bb的几何拆分/有界整区合并：minimum_sector_residues=4仅是合并触发阈值，最终只验连通/直径并允许SMALL_DISPLAY_REGION；当前sectors.py仍与该提交同blob。旧35/30版本7区无小区，却违反当前26/18合同，不能回退。14次当前整区合并均不合法，但9个边界节点重分配可使11区全部≥4残基/≥50Å²、71节点完整唯一及连通有界，因此本例小区并非几何必然。只验几何就称正常的判断已更正；方法文档小区必合并/35与30默认值过时。审查证据、重放脚本、可行见证及限制见 [小分区历史版本审查](ctla4_small_region_version_audit_20261008.md)。生产算法/历史运行未改、科学新增0、7区可行性与普适修复尚未证明。

## 12. 2026-10-08：小分区归并修复已实现

用户批准最近可行邻区归并及边界重分配方案。`joint_design_v4/sectors.py`在现有有界整理后增加确定性局部回溯，公共返回门强制最少4残基/50Å²、连通、26/18直径及完整唯一覆盖。失败抛出带诊断异常；`engine.py`保存诊断和显式无归属原图，停止候选设计；策略版本connected_bounded_v3在各导出一致。CTLA4保存输入重放得到9区/71节点，大小按新编号12/5/4/13/9/9/9/4/6，最小面积133.224Å²、最大测地23.291Å/Cα17.991Å；完整PD1 AF3保存图仍7区/58节点。原结构、映射、analysis、图与confidence前后hash一致，科学新增0。

初轮相关86项回归通过；扩展PD-L1完整流程暴露合法受体边界重建遗漏，随后补充回退：小区与受体部分成为合格区域，剩余受体也保持连通及4/50/26/18，区域数不变。PD-L1现17区/117节点，CA-only合成PD1也恢复正向分区，不改原正向流程测试预期。最终114项回归398.48秒全部通过（1条已有弃用提示），包含完整pipeline设计序列/报告/ZIP；随后补齐exporting.py共享文件清单，完整pipeline与导出补验8项249.77秒通过，诊断/graph/display/sectors四类文件的ZIP、磁盘及最终manifest SHA一致。局部搜索失败不证明全局无解；总预算100000状态，7区仍软目标，真实突变推荐未重算、历史目录与页面未改、运行服务未重启。详细验收、重放证据与限制见 [小分区修复验收](small_partition_repair_20261008.md)。

## 13. 2026-10-08：覆盖优先与逐区组合建议

用户要求新增覆盖为最高优先、每区合格或最近门槛组合。identifiability/optimizer的公共主选择器使用合格表面直接面积增量优先，原信息净目标不再否决正覆盖；独立主预算及正式门保持。新增sector_recommendations按实际突变残基匹配区域，旧sector_id/anchor只是生成窗口。缺少通过者时按独立未过门数、缺失证据数、已知归一化超限量选择审查组合，展示模型值/阈值/原因；共享组合只投影一次。engine与pipeline保存身份及policy，网页、reporting、offline与ZIP新增相同JSON/CSV结果。近门槛者不加入正式预算或Stage2。

CTLA4作业5068b62450ec保存证据只读重放到 `workspace_data/ctla4_coverage_first_regions_20261008_final`：8个通过资格门的主组合、9条逐区建议（7通过、2待审查），直接面积coverage=0.256274。VS008 M134W/Y140W的SPURS0.648233>0.5且FoldX缺失；VS009 E155D/S50I模型差1.978935>1.5。不是已证实不稳定或实验结合效应。97项相关回归与真实Chrome/多格式一致性通过，新增科学调用0，来源9文件hash未变；完整pipeline补验另见tasks/todo.md。[验收记录](coverage_first_sector_recommendations_20261008.md)。

完整PD-L1流水线补验4项通过（183.32秒），逐区JSON及ZIP与最终结果一致；初次断言误用panel字段，修正为selected_panel后重跑通过。自动审批拒绝本地服务重启（blocked by policy），命令未执行，后端新源码待正常重启加载，现服务health正常。

用户随后引用相邻双点覆盖冗余讨论。只读审查确认candidate_search beam与engine FoldX队列仍先风险，同一位点组多替换可抢占预算；geodesic_0直接残基面积不能判断邻域重叠。该上游优化尚未实施，具体位点多样性、增量对照、真实mask与几何代理分离及固定预算验收见 [分区内设计优化](intra_sector_design_optimization_20261008.md)。引用旧会话只有预览可用，不声称已读完整旧回答。

## 2026-10-08：正式网页交付

用户随后指出分区6缺失：VS006推荐已存在且通过，V4_07（Q76Y/V81W）旧region_id=3表示生成窗口，Step5构建体卡片按此分组。web/app.js共享regionsForDesign改按实际突变所在当前分区，卡片组及单区定位同步；跨区构建体复用展示，预算保持8。真实Chrome无模拟路由回归新增“补丁6包含V4_07且可选择”断言通过，证据workspace_data/ctla4_region6_browser_20261008。历史JSON/ZIP仍保留原生成窗口编号，未修改冻结证据；本次修复为正式网页显示映射。

保存评分重选此前未登记到生产JobStore，报告/模拟接口验收不能证明正式网页可见。本次新增独立COMPLETE任务 `6e0451896ef4`，状态与结果同一事务提交，指向重选目录；260条已有任务逐字段保持一致、源运行9文件hash不变、科学新增调用0。正式入口 `http://127.0.0.1:8765/?job=6e0451896ef4`；真实API列表/任务/结构/ZIP以及无路由模拟Chrome验收通过，9区、8主组合、2待审查、覆盖与未覆盖残基/门槛原因显示正常，无pageerror。证据 `workspace_data/ctla4_live_delivery_20261008/registration.json`、`workspace_data/ctla4_live_browser_20261008/browser_verification.json`。结果为旧已评分池重选，不包含新生成未评分组合。

## 14. 2026-10-08：分区内互补与整区覆盖排序已实现

用户提供完整建议并要求尽量覆盖整个分区。新增spatial_support.py按全合格表面图/面积计算9/12/15Å并集、低阶边际、每位点leave-one-out及重叠；三点比较最佳双点。常规candidate_search位点池先扩大目标窗口覆盖，beam/Pareto先不同位点组再替换；全局去重后统一全表面分母避免显示归属/任意保留窗口影响科学池。engine FoldX在原预算与假设配额内选不同位点组，主面板direct面积同分时优先互补。逐区按真实成员面积交集的三尺度等权平均覆盖比例优先；近门槛仍先门距离。JSON/CSV/候选审计/报告/网页/离线同步coverage和uncovered_positions，unknown不伪造为0/100%。科学门及物理profile路径不改，不加最短距离或新几何硬门。

固定单点输入old HEAD与新生成对比：主候选184→259、独立位点组129→254、平均内部几何边际189.10→345.33Å²；12Å每区有100%代理候选（6→9区），尚未评分。全候选244→319，既有假设配额软容量随更少重复而扩张，因此实际SPURS multi工作量可能增加；FoldX42及主面板预算不变。9Å并非所有区域改善，VS001最佳代理87.7→86.2%。证据 `workspace_data/intra_sector_geometry_audit_20261008/generation_comparison_final.json` 及新未评分池，历史输入hash未变。

已评分池重选至 `workspace_data/ctla4_region_coverage_20261008`，仍8主组合及9逐区建议（7通过2待审查），VS003各尺度本区coverage51.9%，不放宽门补齐。16项最终空间/完整pipeline补验通过（含4项完整PD-L1、261.64秒），最终相关全套结果见tasks/todo.md；真实Chrome及JSON/CSV/ZIP一致，模型新增0。服务最初未监听且G盘环境缺FastAPI；无活跃任务/子任务时项目虚拟环境启动并最终重启，health正常，最新截图证据 `workspace_data/ctla4_region_coverage_browser_20261008_final2`。旧运行仍展示原结果，最终新算法用于后续任务。详细接受/拒绝建议与边界见 [分区内优化验收](intra_sector_design_optimization_20261008.md)。

# 2026-10-11 PD-L1待审查算法修复（验收完成）

审计与收据见 `docs/pdl1_review_gate_repair_20261011.md` 和 `workspace_data/pdl1_review_gate_audit_20261011`。旧b4b51cda70f6保留；局部面积门被全局hypothesis响应稀释、FoldX先确认无增益三点、含H源结构与无H评分输入SASA口径不一已确认修复。新增一次最多12个双点的无SPURS通过分区反馈补筛，保留真实组合模型门与42个FoldX上限。反馈任务119b97d86de9真实SPURS476/FoldX42成功，105节点10区专属首选全部通过、主面板8；交付18c0b8ee73fa仅重放说明与共享网页，评分不变、科学调用0。15备选12通过/3待审，对应2个未安排FoldX的跨区物理组合。60文件/1407引用/API/JSON/CSV/ZIP及正式/实际离线×1500/390浏览器验收通过；4区9/9/9/8偏离10–15软偏好，不称全局无解。全图测地30/Cα22是直径门，区内路径单独诊断，修正网页混用标签；原结构/映射及全部历史证据保留。


## AF3 binder 五例最终验收（完成）

五个原始任务及两例失败恢复全部COMPLETE；原始80组/300槽核对通过。29个展示组合追加58个缺失重复槽，独立交付目录和新COMPLETE记录保留原结果。五例全卡片/真实配对/几何/源hash/ZIP验收通过；五个ZIP实际解压后file网页七步、结构与本地序列导出、零HTTP请求实测通过。最终报告：workspace_data/af3_binder_batch_20261017/final_acceptance_zh.md；收据：final_acceptance_execution.json、af3_batch_acceptance.json及final_offline_browser_audit/acceptance.json。后台进程使用一次性OS退出等待，无循环监控或自动化。

## 2026-10-18：11例binder审查与普通V4复发防护

审查入口 `workspace_data/20261016_raw data_binder/full_review_zh.md`，旧交付快照独立证据为evidence_audit.json。11例77区/583表面节点、689卡片完整性与几何通过；正式327、主面板88。发现7个旧包漏420条refinement清单行，AIB60还漏顶层240条；原始模型未丢失。新增派生包及注册结果以article_repacked_manifest.json终态为准，不能用旧ZIP的hash通过代替源配对集合闭合。

公共修复：sequential_foldx.display_priority优先全部主代表再轮转备选；engine严格普通任务从complete_display_pool生成卡，保留raw池与未展示ID；pipeline导出后写全卡完整性收据，缺证据卡不发布COMPLETE，无卡明确诊断。普通输入处理归档native与confidence并保存SHA身份；exporting与replay携带；普通export_names包含top paired及active refinement清单；共享supplement重导出复制阶段清单、拒绝冲突。特殊物理评分合同与科学门不变。

78项相关回归通过（42定向、31含完整pipeline、5离线/标准ZIP；最后补top清单后5项再通过）；闲置服务经DB无活跃任务、进程身份与无子进程检查后加载，health通过。详见regression_results.json和service_activation.json。历史原件保留，审查新增模型0。原AF3五站file浏览器5/5；本轮7个重建站未重新逐站file交互，来源实验类型、SASA设置与外部训练独立性也未补造为通过。后续审查以review_acceptance.json声明缺口为准。

上述重包终态已完成：article_repacked_manifest.json七条均注册COMPLETE（该清单顶层阶段名REPACKED），2143条活动pair核验、420 refinement加240初始清单恢复、评分/序列/源hash不变；最终11例DB均COMPLETE；两项发布元数据差异已在11条当前交付记录同步，前值备份及核对PASS；本轮未宣称11站API/浏览器全验。交付入口full_review_zh.md、review_acceptance.json、review_package.zip和final_delivery_check.json，状态PASS_WITH_DECLARED_GAPS。七个新网站不以旧file浏览器证据代替新版本交互验收。

## Binder复合物受体显示（后台交付重建）
共享renderer新增独立binder/partner组件，默认隐藏整体组件和partner，复用同一viewer；自建组件保留polymer标签以沿用表示切换和overpaint。第四步显示/隐藏受体，切换保留视角并平滑聚焦。12项生命周期/离线导出回归通过，AIB234真实file页面切换/表示/着色及无外网请求通过。11例更新脚本workspace_data/20261016_raw data_binder/update_partner_delivery.py正在后台更新六例离线产物、检查11站和一站在线并生成v2分享ZIP；只有partner_visibility_delivery.json为PASS才代表新版离线交付完成。旧分享ZIP保留。旧结果目录更名后以article_repacked目录junction维持已注册在线路径兼容。

Binder受体开关最终交付PASS：partner_visibility_delivery.json绑定v2 ZIP及SHA，11例file网页（六复合物开关/五AF3无开关）和AIB234在线共12项浏览器检查通过。在线此前加载ecd_structure.cif仅含binder，现复合物显示使用structure.cif完整坐标并隐藏受体；不能假设在线与离线加载同一结构上下文。六例analysis/组合评分/序列/主面板共24文件与历史ZIP字节一致。模型调用0，旧包及before_partner_toggle历史保留。验收脚本修正在线模块私有变量访问为DOM/真实组件状态。

## 2026-10-18 受体仅限蛋白蛋白界面（修复）

用户指定v2解压站点真实复现AIB234离开第四步进入第五步仍showPartner=true、受体组件可见。根因setStep只切面板，没有清除共享viewer的受体状态。web/app.js在离开第四步隐藏受体并同步按钮；开关回调也检查当前步骤，异步失败不在其他步骤恢复显示。

14项生命周期/离线导出回归通过；六复合物真实file站逐一验证4→1/2/3/5/6/7、突变定位/多选/清空、表示切换、相机重置与同viewer，零HTTP/pageerror。仅更新指定解压目录的app.bundle.js及manifest对应hash，逐文件清单证明科学文件不变；原v2 ZIP和旧脚本保留，另存v3。根因、验证和交付收据见docs/binder_partner_step_repair_20261018.md及workspace_data/20261016_raw data_binder/partner_step_repair_20261018；最终归档以delivery.json为准。未重算模型、未跑全仓或重验在线/五个未改AF3站。

## 2026-10-18 定位结构后受体复发：表示层修复

上次父组件隐藏检查假阳性：AIB193定位S1D/R4A重建分子表面后，Polymer/partner父组件仍hidden=true，新表示与renderObject却均可见。renderer.setRepresentation现重施完整原组件与partner的子树可见性，color和开关共用同步方法，完整重载回退重建链组件。selectDesign删除先行重复paint，避免与删建表示争用state ref。

15项回归通过（4.52秒），新两项检查分别能使旧renderer与旧selectDesign失败。真实六例验收改核对representation、renderObject、可见原子链及捕获viewerError，含截图。证据partner_representation_repair_20261018与docs/binder_partner_step_repair_20261018.md；v4最终以delivery.json为准。原v2/v3保留，科学文件逐字节核对不变，新增模型0。此前step修复PASS不得当作实际绘制层的通过；未重验在线或五个未改AF3。

## 2026-10-18 用户授权清理本任务中间产物

已按明确清单将70个中间文件（1,543,570,015字节）移入回收站：旧网页备份、DB修复前result_json备份、一次性workspace脚本、旧显示验收与调试文件。用户授权覆盖历史中间备份保留默认约束；原始输入、科学模型运行目录和当前在线来源article_repacked保留。v4 ZIP/文件夹、报告/汇总、最终显示验收及两份回归夹具保留，关键文件SHA前后一致、12个网页入口在位。cleanup_manifest_20261018.json记录PASS、清理路径和保护文件SHA。旧文档/收据仍可能引用已回收的中间路径，以该清单说明其去向；回收站未清空。

## 2026-10-20 开题报告与实验进展确认

用户本次明确确认：部分PD-1、CTLA-4构建已完成测序及Myc表达检测。应视为已有实验进展；具体构建清单、原始测序/流式图及表达数值尚待归档，不扩展为天然折叠或表位确认。当前只读数据库activation e0f2c26650f04681仍revision16、162样本/3份导入，人工QC/effect均0，只代表工作台分析登记状态，不能否定用户已完成的线下实验。4样本Myc项目也不代表线下全部工作。

开题报告已按现有结构映射、连通分区、六轴化学机制、SPURS/FoldX、序贯预算与实际案例重写，新增完成状态和证据索引。交付E:/download/开题报告-chw/项目进展补充版（21页正文、15页答辩内容）；原稿保留，未运行新模型或修改数据库。报告区分CTLA-4旧评分、PD-1局部候选、PD-L1流程及11例binder迁移，不作为同批表位性能。


## 2026-10-22 CTLA4 R51A/E59Y/K65F 三点风险直算

同一P16410/1AH1连续编号输入上，FoldX 5.1/timeRotabase五槽直接计算三单、三双及R51A/E59Y/K65F三点；35条配对身份和原始Dif验证PASS。目标三点中位数+2.537 kcal/mol，范围+0.582至+2.924，3/5高于项目+2.0复核线；风险主要跟E59Y一致，故标高风险复核项，不推断确定失折叠/表达或结合。旧V4 E59Y/K65F +1.21为legacy记录且缺逐槽输出；本次同二进制新五槽双点为+2.951，差异无法由旧原始数据解析，未混合均值。完整报告与输入/模型hash见 `workspace_data/ctla4_triple_risk_51RA_59EY_65KF_20261022/`。调用经历作者/连续编号误配和路径长度失败两次未出分，均保留；纠正规则见 `tasks/lessons.md`。历史正式面板未改。


## 2026-10-23 CTLA4 E59Y独立第二批复算
E59Y与E59Y/K65F各新增5槽，旧配对复用0，运行106.99秒；原始输入、RepairPDB产物及FoldX二进制SHA均与第一批相同。10/10原始Dif/模型hash/完整序列/采样槽核验PASS；两组WT/Mut均5种内容且与第一批无内容重合。第二批中位数分别+2.954（4/5>2）与+2.517（5/5>2）；第一批分别+2.821与+2.951。维持FoldX稳定性风险复核，不合并批次为新正式合同，不改历史面板。旧legacy双点+1.21因缺原始槽仍无法精确归因。报告workspace_data/ctla4_e59y_recheck_20261022/recheck_report_zh.md。未重算三点/SPURS/其他骨架。


## 2026-10-23 CTLA4旧FoldX推荐与严格复现审计（更正）
本条替代前述CTLA4条目中“旧逐槽数据缺失、旧值无法复现”的结论。complete派生报告的selection_replay.source_run指向5068b62450ec，原数据完整；9来源hash相符。旧E59Y/K65F在foldx_batch_4afa2724第5组，五槽2.99488/0.276247/3.45390/1.21002/0.710114，中位1.21002→1.21。42组210配对及全部中位聚合独立通过。原batch同input/binary/rotabase/列表顺序默认调用严格复现25 DDG与50模型字节；单独调用该双点中位2.07290。批次上下文属于计算身份，新timeSeed采样不能当旧调用复现。旧推荐按中位<=2等合同正确放行，但标签未体现采样敏感性；不临时改为五槽全过门。新增受控实验30配对/210.24秒；两近期风险报告及清单已勘误，更正前副本保存。原源任务和正式资格/面板未回写。报告workspace_data/ctla4_legacy_repro_audit_20261023/audit_report_zh.md。

## 2026-10-23 FoldX重复门限敏感性诊断
原因：五槽median按现行合同放行并不说明五个采样槽都位于门内；旧E59Y/K65F即有2/5槽>2，而变更候选批次组成后同输入/工具的中位数也可跨过2。`joint_design_v4/stability_repeat_diagnostics.py`被engine、reporting、本区推荐和保存评分重导出复用；导出最小值、最大值、门限两侧计数、`sampling_decision_sensitive`与`recommendation_risk_label`。只在五个预期槽有效完整时输出敏感性布尔值，缺失/NaN/错误槽身份为未知；不改变median资格合同。E59Y/K65F五槽显示2/5超门且仍按median维持共识通过。父目录venv相关6个测试文件92项通过。未实现历史目录自动扫描或把已核验跨执行证据自动关联当前候选；报告字段保持`cross_execution_evidence_status=NOT_ASSESSED`，来源身份需结构化绑定后再显示候选级冲突。详情见`tasks/lessons.md`与`tasks/todo.md`。
UI续修：`web/recommendations.js`和`web/app.js`同序读取本区`stability_repeat_diagnostic`、稳定性provenance、impact evidence及`design.foldx_diagnostics`；采样敏感/重复不完整/旧审计缺诊断时卡片标签分别显示需复核、不完整或诊断未知，并覆盖兼容旧低风险enum；完整数据展示配置门限、范围、门两侧数量和编号槽值。引擎导出槽值另保留带缺失位置的对齐列表，避免UI把缺失槽后的值错标。全槽在门内只标“本批未见门限跨越”，不写跨批稳定。Node渲染、离线及生命周期相关6个文件32项通过；不做自动跨执行历史搜索。


## 2026-10-23 FoldX候选隔离、缓存与可重放证据修复（验收完成）
公共score_many现为candidate_isolated_v1，一组合一进程、组合内固定排序；模型次数/科学阈值不变。worker保存请求/源码/原始结构/命令/输入与rotabase及输出hash，旧布局或协议不兼容缓存不复用。默认流只取一个完整原生batch，不能重启补后缀或混外部三槽与自身完整五槽；时间播种保留固定槽增量。推荐/CSV/报告/网页增加跨门诊断与未知标签；新旧来源、可移植配对和执行依赖进入普通与重放包，旧CLI/source未知不补造。234项相关回归PASS，真实新增15对：默认双点2.07290及10模型匹配历史单独调用；换序并发2→1复用10对/新增0；timeSeed1→3→5新增1/2/2总5。旧210/新10对包hash与移包核验PASS。未重启在线主服务、未全仓或真浏览器验收。审计docs/foldx_reproducibility_repair_20261023.md和workspace_data/ctla4_isolation_acceptance_20261023/final_acceptance.json。


## 2026-10-25 PD1独立突变图统一seg1–8
此前新增seg总览但突变卡片/膜图/说明仍用旧R分区，现独立atlas全部使用统一8段定义及灰色范围外位点。保留mut1–5和原mut4/5/7，按所选位点方向设相机；历史距离审计另存。11画布/288位置颜色/逐位映射/CB标记/零外部资源及PNG核验通过，独立只读审查通过。无管线或评分变化；seg8为预测茎区。修复与验收记录：[seg_only_update_20261025.md](pd1_antibody_epitope_atlas/seg_only_update_20261025.md)。


## 2026-10-25 PD1文献8嵌合体与现有8设计共同结构对照
新增独立PD1_mutant_comparison.html（原seg图保留），报告逐表与ECD序列验证45替换，现有8设计23替换；58独立位点。采用统一AFQ15116v6 ECD25–170预测坐标覆盖缺失N端/茎区，中性表面、构建体多色定位、共享位点多标记，单canvas离线。主面板L1–8+C3/C4/C5=11，其他5保留补充；原mut4/5/7保留A132L背景并建议单独对照。文献MPER译名/Seg3坐标计数纠正在网站标明，不用错误SASA汇总筛选。桌面/390px、逐位/交互/零网络验收PASS，数据来源hash未变；用户报告目录已复制网站、PNG、建议与验收JSON。详见docs/pd1_antibody_epitope_atlas/PD1_mutant_comparison_acceptance.json及用户目录PD1_mutant_comparison_delivery.json。未调用模型评分或重算突变结构。
