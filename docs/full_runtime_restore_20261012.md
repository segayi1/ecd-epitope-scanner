# 系统重装后完整科学评分链恢复

## 目标与验收

以独立 PD-L1 新任务验证原生 SASA、ProCaliper、DiscoTope、SPURS 单点与真实多点组合、FoldX，以及正式 API、网页和导出。不能以历史评分加载或单独 worker 成功代替完整任务成功。FoldX 默认最多60个组合，每组合5次，原科学门槛不变。

审计目录：`workspace_data/full_runtime_restore_20261012`。冻结输入来自 `G:/ecd_acceptance/pdl1_review_repair_20261011/119b97d86de9`；新任务不导入历史候选评分、不复用历史 paired-model 缓存。输入及历史文件哈希保存在新任务的 `job.json`。

## 已确认故障与修复

首次完整任务 `07bbf0dd36ca` 的 Windows→WSL 桥返回 `Wsl/E_ACCESSDENIED`。`.venv/Scripts/python.exe` 继承 Low Mandatory Level，实际进程也是 Low；直接 Linux worker 已证明模型环境和权重仍可用，因此不是模型损坏。

只对 `.venv/Scripts/python.exe` 和 `pythonw.exe` 设置 Medium 完整性标签，未递归修改项目或 WSL 磁盘权限：

```powershell
icacls .venv\Scripts\python.exe /setintegritylevel M
icacls .venv\Scripts\pythonw.exe /setintegritylevel M
```

随后以实际 `.venv` Python 执行 WSL `/bin/true` 成功，新任务进程为 Medium。已启动的进程不会随文件标签修改而提升身份，因此旧 Low 网站服务已停止并重新启动。启动器增加实际 Python→WSL 执行预检；失败时给出原因并停止启动服务。请从正常用户会话运行启动器。

## 分阶段实测

- ProCaliper 0.2.1 / BioPython 1.88：215条有限 SASA 输出。冷导入约44秒，25秒探测超时不足以判断不可用；生产180秒限制足够。此工具仍仅承担 SASA 交叉检查角色。
- SPURS single：真实 CUDA 推理，215个位点各19个非WT替换均有限。
- SPURS multi checkpoint 冒烟：Q65A有限输出，仅验证接口；真实双点/三点以完整任务为准。
- DiscoTope-3.0：真实 CUDA 推理215条，score、calibrated_score、rSASA均有限。
- FoldX 5.1：真实双点 D73N/L50I（局部编号 D55N/L32I）5次 BuildModel，5组 WT/Mut 身份及哈希通过，无 paired cache；中位数0.560916 kcal/mol。
- GPU三次冒烟均 RTX3060 / torch2.5.1+cu124，无OOM或CPU回退；完整源码与权重哈希见 `gpu/gpu_smoke_summary.json`。CPU证据见 `cpu/cpu_runtime_receipt.json`。
- 运行配置、GPU恢复逻辑、FoldX执行路径与解析、PD1 AF3映射/置信度及完整fixture流水线回归：25项通过，收据 `runtime_regression.xml`。
- 新服务wrapper 32432/listener 26740均为Medium，收据 `gpu/server_integrity_receipt.json`。真实生产Python→WSL及健康API预检通过，收据 `launcher_preflight.json`。

## 完整任务

任务 `484899e022f3` 在修复后以真实新评分运行；收据为 `full_run_medium/job.json`。完成标准由 `run.py` 的断言约束：五个适配器全部成功、每残基215条、无适配器缓存命中、FoldX 1–60条、正式主面板两模型均有限、覆盖与映射一致、源文件未修改。失败保留 FAILED 状态和日志，不交付成 COMPLETE。

该任务已 COMPLETE，完整真实运行与独立验收均通过：

| 算法 | 实际新输出 | 设备 | 适配器耗时 |
|---|---|---|---|
| 原生 SASA | 215条；SASA、rSASA、侧链SASA有限且非负 | CPU | 包含在分析阶段 |
| ProCaliper | 215条 | CPU | 1.33秒 |
| SPURS single | 215个位点的完整替换矩阵 | CUDA | 82.44秒 |
| SPURS multi | 476个真实组合 | CUDA | 148.90秒 |
| DiscoTope | 215条 | CUDA | 38.82秒 |
| FoldX | 60组×5次，共300组新WT/Mut配对 | CPU/8路 | 976.25秒 |

五个适配器均 `success`，缓存命中为false，worker警告为空；FoldX每组new_model_counts=5、reused_model_counts=0。300组独立(batch,group,run)身份、WT/Mut/PdbList/Dif哈希合同与有限评分逐一核验通过。规范化输入和残基映射哈希与冻结来源一致，历史八个来源文件均未修改。

105个表面节点完整唯一分配到10区，大小10/9/9/10/9/8/11/15/13/11。10个专属首选均通过，正式主面板8个构建体，全部同时有有限SPURS multi与FoldX。9/9/9/8的四区偏离10–15软目标，未放宽连通或直径约束。

独立API/导出审核核验65个manifest文件、1022条展示引用及JSON/CSV/ZIP/结构下载一致性，PASS。推荐CSV25行是首选与备选记录，不是25个主面板构建体。真实正式网页显示COMPLETE/100%、10专属首选、两模型数值与展开准入证据；结构定位与序列选择一致，只有1个canvas且无console error。截图 `full_run_medium/real_site.jpg`。

收据：`scientific_acceptance.json`、`native_sasa_acceptance.json`、`foldx_pair_acceptance.json`、`export_acceptance.json`、`browser_acceptance.json`，均位于 `full_run_medium`。正式任务：[484899e022f3](http://127.0.0.1:8765/?job=484899e022f3)。

## 边界

计算成功不证明实验结合或真实扰动范围。ProCaliper 与原生 SASA 的软件版本及计算角色分别记录；SPURS源码有既有修改，不能只用上游commit代替真实源码哈希。10–15个表面残基为分区大小软偏好，不能为凑数放宽连通和直径硬约束。本次保守性查询执行成功但只有1条可用直系同源记录，明确为insufficient_orthologs/未知评分；未导入IUPred证据、未执行可选研究物理预览，不将缺失当作0。
