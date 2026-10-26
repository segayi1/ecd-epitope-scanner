# Binder 输入适配与文章来源批次

## 目标与范围

用户确认对 binder 本体作表面扫描与突变设计；11 个输入包含 5 个 AF3 单链和 6 个文章来源复合物。本轮仅将 6 个文章来源复合物投入完整模型流程，5 个 AF3 暂缓。文章下载来源并不意味着结构经过实验解析。

## 输入核查

`workspace_data/binder_input_preflight_20261016/inventory.json` 保存原始文件 SHA256、全部 17 条链序列、作者编号/插入码及旧自动选择。`binder_chain_acceptance.json` 保存 11 个 binder A 链的实际 SASA/独立能量链导出核查；各源文件哈希保持。六份复合物 B-factor 均为零，不能解释成 pLDDT；五份单链的 B-factor 含义也未由独立文件验证。源目录内 `__MACOSX` AppleDouble 文件不作为结构。

旧 Ig 序列启发式将六个复合物的 B 链列为 antibody_chains；这不构成生物身份依据。新 binder 模式要求多链输入显式给出 binder 链，单链自动选择唯一可用链。本批使用 A；B 链只保留复合物原始上下文/接触注释，不进入 binder 单链 SASA 或能量评分结构。

## 实现

### 小蛋白分区尺度（用户后续澄清）

用户要求针对约100 aa的小蛋白，每区目标约8个表面残基，**总分区数4–12仅为软范围**，不是固定8区。工程范围按实际分析链≤150 aa；自动偏好8/6–10，长链仍为12/10–15，显式尺寸优先。硬门仍为每区至少4残基/50Å²、诱导连通、完整唯一覆盖、全图测地≤30Å/Cα≤22Å；物理预览冻结的旧26/18策略不改，构建预算不随分区数改变。

config.display_sectoring_parameters负责长度默认，engine只用未合并调用参数判断显式覆盖；API/PipelineRequest支持显式8–15。网页目标留空才自动，不再每次默传12；i18n/资源版本同步。sectors的4–12偏好只加入软优化与偏离诊断。scripts/replay_v4_selection.py新增--repartition-small-protein，在独立目录重建分区、归属、建议与报告，保留原评分身份和来源，不生成新候选、不调用模型；坐标派生binder不再要求不存在的uniprot.json。

原始AIB234/AIB226计算完成后保留原任务和评分，独立生成AIB234_v2及其余5例重分区交付；每例保留源候选身份/评分，重分区新模型调用为0。具体区域尺寸、节点数、NEAR_THRESHOLD原因和HTML/ZIP验收见`small_protein_delivery_acceptance.json`及`article_batch_summary_20261016.md`。

重分区用冻结候选池时不能声称完整新区候选搜索：六例独立验收均PASS，36区中29区PASSED、7区NEAR_THRESHOLD（AIB234三、AIB226二、AIB239一、AIB201一）。区级缺FoldX、SPURS超门或模型差异超门的细节见`workspace_data/binder_input_preflight_20261016/small_protein_delivery_acceptance.json`。独立verify_small_protein_repartition.py核验全部节点、几何、评分身份、JSON/CSV/FASTA/ZIP哈希一致；六份保存候选池重放的新模型调用均为0。最终总览`article_batch_summary_20261016.md`。

最终相关79项回归19.12秒通过，含大小默认/显式8/API/导出/AF3/旧小区硬门；JS语法、i18n与范围内diff空白检查通过，仅有已有Starlette/httpx弃用和LF/CRLF提示。收据small_protein_verification.json。六例原始任务与派生分区交付已完成；服务在任务空闲后重载以加载新Python/API默认。未跑全仓或实验，不把分区几何和模型代理提升为结合结论。

历史首轮FAIL收据保留；最终检查器运行后六个文章来源原始任务均PASS，批次及区级科学门缺口分别记录，详见`batch_acceptance.json`和总览报告。

- `PipelineRequest.analysis_target` 位于 dataclass 末尾，默认 antigen，避免已有位置参数移位；binder 仅适用于上传结构，移除固定抗原 assay profile，实验格式/膜方向保留未知。
- `pipeline.py` 从选定链坐标派生序列与编号，不虚构 UniProt 映射。binder结果导出膜锚为 null，独立 V4 来源保留 analysis_target，序列及后续验证文案使用 binder 语义。
- `/api/jobs/upload` 与上传对话框传递 analysis_target；CLI 支持 `--analysis-target binder --chain A`。原抗原模式默认不变。
- binder 不使用 DiscoTope 抗体表位先验；SASA、ProCaliper、SPURS single/multi、FoldX 和原稳定性/空间/序列门沿用。`_adapter` 显式 enabled=false 优先于环境 worker 命令。
- 不把 PDB B-factor 自动当作 pLDDT，也不将普通 binder 上传标为官方 AF3。没有配套置信度时，坐标质量标 NOT_ASSESSED；设计分数为计算证据，不等于结合或折叠的实验结论。

## 验证与运行状态

11文件/17链精确映射及11个binder链SASA/能量结构检查通过；AIB203单体和AIB239复合物的禁用模型预览均完整生成表面分区及报告导出，分别4区/64残基和3区/64残基。这两份预览不是完整科学评分。

新增API3项、pipeline5项及既有PDB/assay/AF3回归合计25项通过，耗时695.88秒；最终修改后8项新增再次通过，与25项重叠不相加。有一条已有Starlette/httpx弃用提示。UI语法、i18n、实际FormData与中英文切换通过；真实浏览器确认上传对象选择及255px窄对话框无横向溢出，图见`workspace_data/binder_input_preflight_20261016/binder_upload_zh.png`。首次真实浏览器加载暴露出旧app.js/i18n缓存与新HTML不一致，新增资源版本参数后验证正常；重复的固定抗原说明已移除。

六份文章来源输入与原始上传哈希关联在`article_batch_manifest.json`，strict、binder A链、FoldX最多60组合×5次，未放宽科学门。六个源任务均由`verify_article_batch.py`独立验收PASS；逐例模型证据、产物/ZIP一致性及科学缺口记录在`batch_acceptance.json`。

临时heartbeat `binder` 仅接续审核已提交的六任务；批次和分区交付完成后停止。`verify_article_batch.py`实际审计六例真实模型/300配对和导出，一律以最终`batch_acceptance.json`为准。

早前AIB203小预算真实评分冒烟在用户澄清分批意图后已停止；停止凭据及部分文件保留，不作为完整运行或本轮文章批次的一部分。
