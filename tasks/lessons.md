# Lessons

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
