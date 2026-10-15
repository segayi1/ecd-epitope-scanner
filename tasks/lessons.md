# Lessons

## 2026-10-13 离线结构展示与正式网站表示方式

- 错误模式：沿用配色，却默认使用卡通结构，而用户需要正式网站的表面表示。
- 预防：对齐参考页面的表示方式；单文件离线包内嵌渲染库和结构，表面计算不依赖外部worker资源。WT表面上的位点标记不能称为突变后的模型。

## 2026-10-12 完整模型链与Windows进程身份

- 错误模式：原生WSL命令成功，被误认为应用Python也可访问WSL；新venv继承Low完整性标签，启动后模型桥被拒绝。修改exe标签后继续复用旧Low服务，仍会失败。
- 预防：用实际生产Python执行WSL预检；只修复确认受影响的启动文件标签，从正常用户会话重启服务并核实进程身份。模型worker单独成功、网站打开、完整流水线成功分开验收。
- 触发条件：Windows重装、复制项目/重建venv、WSL_ACCESS_DENIED。不以递归修改WSL磁盘ACL绕过进程身份问题。
- 诊断：冷导入探测时间应与生产限制匹配；默认忽略models目录的rg结果不能证明权重丢失，应核实实际路径、大小与哈希。

## 2026-10-12 系统重装后的环境与启动验收

- 错误模式：venv文件仍存在，被误认为依赖仍有效；继承旧Anaconda全局包或无效/不可写TEMP；把端口打开和旧评分加载当成新模型链可运行。
- 预防：检查pyvenv.cfg和真实导入，保留旧环境后按锁文件创建独立环境；启动器预检、记录stderr、给首次冷启动足够时间，进程范围内使用项目临时目录。分别报告网页、预算配置、历史读取和新模型评分的验收状态。

## 2026-10-11 局部覆盖门、确认顺序与原子口径

- 错误模式：局部面积增益门使用全局hypothesis响应；已知无增益三点先占用FoldX预算；显式H输入与无H评分结构做SASA分类比较。
- 预防：门槛绑定量纲、空间范围和对照身份；缺失不能当零。便宜且确定的资格检查放在昂贵确认之前，为双点保留基础机会。源链与规范化ECD统一重原子、探针/点数/半径，再审查真实上下文变化。
- 边界：几何对照不是稳定性证据；统一口径不允许豁免仍阈下的位点。重分区后编号变化须按残基身份解释，历史任务保留。
- 排期与展示：确认计划必须在推荐生成前传入，不能事后只有每区统计却把每个未评分候选写成“排期未知”。直径标签依约束使用全图或诱导图的实际口径，不能凭是否存在region_metrics决定。新长组合标签须验窄屏父容器的最小宽度。

## 2026-10-10 分区可行性、理想大小与画布生命周期

- 错误模式：最少4残基的硬验收被误认为已满足每区约十多个残基；局部合并失败被说成几何必然；同一结果刷新重复创建共享画布实例。
- 触发条件：直径受限的稀疏表面图、固定中心数及首个可行边界修复，或adoptResult/renderAll都负责初始化。
- 预防：大小目标和硬底线独立记录，先核对冲突下界与可行见证，再提出目标；局部失败不证明无解。画布初始化只由一个入口负责，任务切换销毁旧实例并防止旧异步结果覆盖。
- 科学/工程边界：图距离是当前图代理，不能称连续物理定理；缩减首屏外置证据须保留完整文件/hash。浏览器当前能载入不等于用户原故障已定因。

## 2026-10-09 版本标识的语义与时间范围

- 错误模式：Torch版本充当工具版本；DiscoTope包1.0被误当模型系列；没有vendor Git时获取了父工作台提交。
- 触发条件：工具包/方法/运行时版本混合，或事后补查历史来源。
- 预防：软件包、模型系列、Git和权重分别记录；检查vendor .git后查询Git。事后环境识别必须说明时间范围，权重一致不能证明历史代码一致；未来worker应在运行时记录真实版本与源码指纹。

## 2026-10-09 数值证据与规则结论分开

- 错误模式：布尔审查复用数值空值渲染；缺pLDDT被当成实验结构缺少置信度；有先验但未达高支持条件被说成无证据。
- 触发条件：混合模型评分、保护/连通规则及上传实验结构的同一详情表。
- 预防：显式检查类型并展示规则事实；从原始结构方法识别pLDDT适用性，实验局部精度未知仍如实说明；区分评分缺失与未达到支持条件。正式主面板重新核对两模型有限值，缓存状态不能替代证据。
- 2026-10-11补充：浏览器功能PASS不保证审查解释可读。最终检查真实卡片正文，对已知上下文暴露分类变化给出具体原因；新增共享标签回归后同时重建离线与ZIP，核对分析SHA不变并重跑正式/离线桌面手机验收。

## 2026-10-09 真实窄屏数据验收

- 错误模式：短测试数据通过响应式检查，真实候选编号在常驻序列栏撑宽页面。
- 触发条件：组合切换后显示长编号，窄屏容器存在最小内容宽度。
- 预防：内部编号只放证据详情；用真实新评分任务切换首选与备选后检查页面宽度，正式网页与离线报告均验证390像素。

- 分区公平不能只用hypothesis锚定数量推断；保障实际本区独立位点组与可验证机会，跨区候选不替代专属首选。同位点不同替换、重复窗口和跨区展示不能重复计模型成本；候选不足、验证未安排、计算失败、科学门未通过分别记录。
- API内存结果与磁盘JSON应在类型上也一致。把按突变阶数索引的资源审计嵌入推荐前，整数键转为字符串，并验证JSON往返等值；不能只验证两个已序列化文件相等。
- 后端原始代码、候选哈希和多尺度审计指标不能拼成推荐正文。默认卡片说明用途和关键风险，完整数值用可展开表格，主页面与离线页面共用渲染逻辑，未知指标不补零。

- V4设计的旧region_id可能表示生成窗口，不能直接用于当前显示分区的构建体分组或定位。按实际突变成员映射显示区，跨区组合在对应区复用展示但不重复计算预算；浏览器验收须检查逐区建议与实际构建体卡片两处，不能只检查VS编号存在。

- 正式网页交付必须核对生产任务登记、真实API结果及文件和无模拟接口的浏览器渲染。报告目录存在或模拟请求通过不能证明用户入口可见；重选历史评分应发布独立任务，原任务保持不变，并明确无新增模型计算。

- 结构网页或模型报告必须同时列逻辑配对槽位数与唯一结构字节数；“5个槽位”不等于“5个不同WT结构”，更不等于5个独立骨架。候选扩展可能只出现在逐条pair_manifest而未进入旧汇总字段，展示前按逐条清单重算并标出重复内容。

- When a tool finishes but a downstream mesh-QC gate fails, account for the attempted solver call in `finally`; cache a typed unresolved result so resume does not repeat a deterministic failure. Never convert an open/nonconverged mesh to zero error or let one failed slot erase valid sibling slots. Apply this whenever numerical bounds are computed over fixed replicate slots; keep non-mesh identity/tool failures fatal.
- When a quorum allows partial valid replicates, exclude nonfinite values before taking the numeric maximum; plain `max` lets a single unknown poison otherwise-valid slots. Test both the accepted quorum and the below-quorum unknown case. Likewise, iterate a frozen candidate population until the target is met or the pool is exhausted; a heuristic scan cap must not masquerade as an exhausted search.
- A non-closed SES produced from a PQR is a local geometry failure: cache that query as invalid/unknown with its exact reason and continue independent models. Catch only the known topology error; identity, parse, and other solver errors remain fatal. Keep the unresolved reason in both cache metadata and candidate reports.
- SHA-256 text is hexadecimal data, so uppercase and lowercase strings represent the same digest. Compare normalized hex strings, but still hash and verify the exact frozen bytes; do not interpret a case-only mismatch as input corruption.
- Freeze canonical maps by aligning actual polymer sequences while preserving author residue IDs and insertion codes; reconcile current reference sequences with mmCIF `struct_ref_seq_dif` annotations. Retired UniProt accessions and fusion/polyprotein mature domains need explicit mapping provenance. Any unresolved sequence/deposition conflict excludes that position from a candidate freeze until reviewed.
- When a frozen contact-scale cohort has one open/non-closed SES mesh, classify that complex's geodesic scale as unknown and continue other independent structures; catch only the exact topology error and let tool/identity failures stop the run. Keep downloaded CIFs cached and atomically committed so a resumed deterministic rank order does not pay network cost again.
- HTTP transfer failures include `http.client.IncompleteRead`, not only TLS/URL exceptions. Retry transient transfer errors, write to a `.part` file and atomically replace only after content validation; never promote a partial coordinate file into the frozen source set.
- Candidate-freeze maps can contain placeholder author IDs such as `?` for unresolved polymer positions. Do not parse those as coordinates or infer a canonical position; index only numeric author residue IDs, then verify every observed coordinate residue against the frozen map.
- On Windows, a project `.venv\Scripts\python.exe` can appear as a low-CPU launcher whose child is the base interpreter named in `pyvenv.cfg`. Before diagnosing duplicate scientific runs from matching command lines, verify parentage, creation time, venv base executable, and CPU activity; report one logical invocation when the parent is only a launcher.
- RCSB polymer-entity sequences can help prioritize Fab-antigen structures, but do not replace coordinate-derived antigen sequences for this cohort: unresolved residues, truncations, multi-chain antigens, and a failed first structure's QC can change which later homolog is accepted. Preserve the frozen candidate order and evaluate later homologs until one passes the actual geometry QC.

- 对无向图边端点做排序时必须逐边沿末轴规范化；沿节点轴排序会制造自环并使测地距离变为零。结构接触抽样应无放回；路径长度非正应作为拓扑故障，而不是接触尺度。触发条件：表面网格或抗体接触图构图、汇总接触点距离。
- 扩展一个候选时，不要只按候选 ID 判定全局数值合同变化。先逐项比较真实影响输入/网格几何；只有 append-only 候选变更且所有网格和数值协议逐字/数值相同才可延续旧合同，并把旧合同hash、追加ID和理由写入历史。候选移除、重排或真实网格变化必须拒绝。
- 同时报告实际物理阳性面积与膨胀后的扫描机会覆盖。若低阶几何膨胀已让候选池达到100%，新增突变仍可增加raw阳性或旧合同hard覆盖，但不能宣称提高扫描机会；结构尺度是机会模型，不是结合实验概率。
- 研究报告的清单字段可能滞后于底层配对清单。核验应以逐条的pair identity和hash为准，保留原始manifest并在独立完成收据中给出可复算计数；不因修正摘要而覆盖原运行身份。

- 异步结果交付必须将可渲染结果与COMPLETE状态一次提交。进度回调提前发布终态，会使SSE消费者关流并读取空结果。核查时在完成回调与返回之间受控暂停，验证所有可见终态都有对应结果；前端异常应显式显示并保留恢复读取路径。结果数量低于目标与页面渲染失败需分别验证，不能由已生成文件或低数量推断页面故障原因。

- 分区验收必须同时核对配置的最小残基数/面积、连通性、直径以及完整唯一归属。不能因几何通过就称整体正常；配置中的minimum可能仅是贪心合并触发阈值。整区并集都失败也不证明小区不可避免，应检查不损伤供体的边界重分配。新增拆分修复时需保存真实分支表面残片回归，避免只用线状图覆盖最小大小；方法文档与例外语义必须同步。
- 新硬约束导致旧正向流程失败时，先搜索合法边界调整见证，不可直接改测试为预期失败。仅完整吸收可能把受体远端节点挤入无处可放的状态；须检查将受体重建为两个均满足大小、连通、直径的区域。PD-L1和CA-only PD1说明该分支不能省略。局部搜索失败仍不证明全局无解。分区失败的图导出必须显式无归属：空sector列表仍需拒绝未覆盖图，防止保存诊断时覆盖原异常。触发条件：增加分区最终验收及失败证据导出。
- 新增engine诊断文件时，同时核对reporting清单、pipeline最终files白名单和ZIP真实成员；磁盘存在及engine文件列表都不能证明最终包包含。回归应逐字节比较诊断、graph、display和sectors的磁盘/ZIP内容，再宣称导出一致。触发条件：新增审计或诊断导出。
- 生成窗口anchor不是实际突变位置归属。判断显示分区是否有合格组合时，应按真实突变残基求交，并允许一个跨区物理组合服务多个显示区；投影不能重复计预算。触发条件：科学假设窗口重叠或显示分区重建。
- 最近门槛候选必须保留独立门数、缺失证据数和已知超限量，缺失值不能当0风险；模型分歧、缺失FoldX与真实能量超限分别展示。待审查建议不升级正式资格，不生成正式Stage2控制。触发条件：每区要求候选但无全部通过者。
- 覆盖优先的最终选择无法弥补上游候选截断偏置。审查所有site pool、beam、全局配额和昂贵模型队列；同一位点组多种氨基酸排列不能垄断预算。直接残基面积无法度量相邻扰动重叠；邻域代理与真实组合mask必须分别报告。触发条件：用户要求优化相邻多点组合的覆盖收益。
- 邻域包含也可造成零新增面积，即使Jaccard低于100%；双点应比较最佳单点，三点比较最佳精确双点，逐位点leave-one-out不能省略。整区覆盖须与实际成员取交集，不能拿邻区面积顶替；零面积分母记null，缺支持不伪造0。触发条件：空间互补或整区覆盖排序。
- 多窗口物理去重后应统一全局几何分母，不能让第一个被保留的窗口决定全局排名。配置的requested shortlist可能是软容量；保留更多独立组合会增大有效SPURS输入数，不能仅因参数未改就声称模型成本不变。触发条件：科学窗口重叠及配额保留。


## 2026-10-10 显示大小目标与画布生命周期

- 错误模式：把最少4残基硬门等同每区十多个的优化目标；固定26/18尺度下盲目要求全区至少10。触发：分区虽硬验收合格但大量4–9小区。预防：先在固定图核查几何下界，明确偏好与硬门，显式配置尺度，多起点/边界平衡并导出大小例外；局部失败不能解释全局无解。
- 错误模式：结果采用和面板重绘各初始化一次插件，旧实例/回包可污染新任务，超大成对诊断拖慢显示。预防：一次初始化、实例令牌/单次销毁、单结构请求；正式/离线共用轻量投影，完整审计文件引用必须可解引用，不能丢科学详情。
- 错误模式：只检查桌面正常就认定手机可用。预防：实际390像素检查document宽度；grid隐式列须允许缩小，不以隐藏溢出掩盖。
- 错误模式：换盘独立任务失去repair缓存后出现FoldX原生崩溃。预防：绑定输入与binary身份，使用已有短执行别名；崩溃保持任务失败，不能用SPURS-only放行。
# 2026-10-10 FoldX超时与完整证据恢复

- **错误模式**：外层worker允许3600秒，但内部BuildModel批次仍在1800秒截断；SPURS成功不能替代FoldX证据。半成品WT/Mut结构成对也不等于已有有效Dif评分。
- **触发条件**：较慢替换、5次重复及并行批次；超时使末尾组合未完成，统一报“评分不完整”掩盖执行原因。
- **预防**：内外时限一致并留收集余量，记录每次原生执行TIMEOUT/退出码；先核对每组全部5条有限Dif、坐标突变、PdbList与哈希，再复用完整组。部分组保持未验证，独立恢复任务补齐，保留原失败任务与目录。
