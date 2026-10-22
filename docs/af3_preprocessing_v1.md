# AF3多模型预处理 v1

算法标识：`af3_preprocessing_v1`。输入清单schema：1。

## 使用与范围

从仓库父目录运行CLI：

```powershell
.\ecd_epitope_scanner\.venv\Scripts\python.exe -m ecd_epitope_scanner.cli analyze --af3-directory "AF3输出文件夹" --analysis-target binder --chain A --output "新结果目录"
```

抗原ECD沿用已有accession/sequence和ECD范围配置。API新增`POST /api/jobs/af3-directory`与`POST /api/jobs/{job_id}/af3-directory`，原单结构/confidence入口保留。

- 冻结原始CIF/JSON，验证原始原子行、身份映射、pLDDT与PAE；不自动读取输入引用的外部URL/路径。
- 先按目标f70（距最佳0.05以内）与平均pLDDT（距最佳5分以内）形成质量池。
- 来源组等权、组内sample等权；共同可靠核心采用固定关系图和标准加权medoid。
- 共识分数距最小值0.05Å以内按可统一比较的合法pTM、目标pLDDT、覆盖度、原始ranking和SHA排序；不足核心时明确回退。
- 选择一个原始构象及其原confidence进入公共流程。额外模型只作诊断，不新增位点排除条件、不融合坐标、不逐模型跑评分。
- 目录最多100模型/总1GiB，单JSON64MiB，token/PAE维度最多2048。核心、全局、局部及域间诊断单独归档。

## 已知来源标签事项

同目录Server输出的sample_id fallback目前使用父目录名，没有按model_index区分不同样本。当前五例均单job/seed、每个模型实际权重0.2，比较结论不变；修复时须保留旧冻结manifest可重放。该事项记录于tasks/todo.md。
