# 2026-10-08 结果展示完成状态竞态审查

## 目标与范围

定位最新作业 `ee5b17b0f7ea` 的计算、状态持久化、API读取和页面展示链条。使用保存结果、临时SQLite库和隔离浏览器；不重新运行科学求解，不修改历史运行或生产代码。

## 检查与证据

1. 最新作业创建/完成时间为上海时间2026-10-08 11:15:54–11:26:23，当前数据库状态COMPLETE、error为空。保存结果含125个残基、11个显示区域、9条设计，其中3条入选主探针。
2. 当前 `/api/jobs/ee5b17b0f7ea` 返回HTTP 200（41,288,897字节）。SQLite result_json与analysis.json解析后完全相同，没有非有限JSON数值。真实Chrome打开 `/?job=ee5b17b0f7ea`，统计、设计列表与结构页面正常载入，无页面或控制台错误。
3. `pipeline.py` 在返回analysis之前调用 `progress("COMPLETE", ...)`。`jobs.py:267` 的进度回调立即发布COMPLETE状态，但 `result_json=result` 在run_pipeline返回后才保存。两次更新之间存在可见的 `COMPLETE + terminal=true + result=null` 状态。
4. 用最新保存的实际result对象替代科学计算，在临时JobStore中暂停于COMPLETE进度回调之后：读取确为COMPLETE、terminal=true、result=null；放行后再次读取，result完整且与保存对象相同。新增科学调用0。
5. `api.py` 的SSE遇terminal即终止；`web/app.js:286` 收到COMPLETE后关流并只读取一次作业，随后直接 `adoptResult(job.result)`。隔离Chrome拦截作业响应为COMPLETE+null，实际前端代码稳定重现：

```text
显示状态：完成
SSE已关闭：true
bad SSE event TypeError: Cannot read properties of null (reading 'residues')
  at adoptResult (app.js:224:34)
  at source.onmessage (app.js:290:11)
```

## 结论与验收

定位到计算结果交付的完成状态竞态，并分别复现后端空结果时间窗口和前端渲染失败。只选3个构建体的预算停止策略不会使该结果页面空白，页面会正常提示未达到最少目标。

原作业运行时的浏览器控制台及瞬时API响应未保存，因此不能证明当时一定命中了该窗口；当前代码缺陷和症状一致，已用相同保存输入受控复现。

修复应在JobRunner统一交付位置，把COMPLETE与result_json放在同一次提交中；前端遇空结果也应继续等待并明确展示读取错误。此次仅定位和复现，尚未应用修复。
