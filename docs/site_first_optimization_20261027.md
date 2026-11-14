# 位点先行优化：实现与验收记录

## 2026-10-27 位点先行 v2 实施及验证中

用户批准方案；新增显式协议ordinary_site_first_v2（request.metadata.joint_v4_ml.execution_protocol.name），继续保留原生产默认。全部合法双/三位点集合先计算几何，机制用至多3个真实witness预排；当前位点集合才展开每位点3→6代表，缺失表位证据记录gap且不得获表位原则推荐。三点对区域最强合法双点做增益门，并用确切AA亲本做对应原则提前比较。SPURS最多4/队列预取/批64，FoldX跨队首批量最多8候选、3/10槽和总300预算，WT修复任务内按输入及二进制哈希验收后复用；模型与配对证据按候选独立、技术失败保留批内成功且停止新派发。

针对性发现并修复：跨原则canonical对象覆盖、跨区同突变单次模型/独立正式门、错开队列复用旧FoldX、SPURS额度耗尽继续消费缓存、物理槽run编号完整性、单候选pair解析失败不吞同批成功；新协议正确接受真实表位状态PREDICTED并区分缺失。卡片测试使用明确标注fixture的侧链暴露数据；无原子侧链证据保持未知。

最终源码对应targeted v3为34 passed/10.66秒，perf comparator最后修改专项3 passed；此前综合225 passed/43.11秒是前一源码快照。最终综合回归仍在运行，stdout会话73309，完成证据预定test_output/site_first_20261027/final_regression.xml，不能据运行中状态宣称通过。六位点×19替换输入、模拟模型：35位点集，准备0.05385秒、代表展开0.00593秒、全流程2.112秒、主进程峰值142,585,856B，记录six_site_19_replacements_final.json；前版同输入模拟99.423秒但搜索空间不同，两者不构成真实模型总耗时比较。

真实性能尚未验收：无真实SPURS/FoldX基准、无13例重跑。site_first_performance.json是打包后独立sidecar，峰值仅当前Python进程、不含子进程/GPU；compare_site_first_performance.py只给初步runtime判定，需独立确认GPU/模型身份、RepairPDB冷热缓存及科学证据，不能单靠其宣称50%通过。新协议无全替换索引DB，显式启用，不自动改变现有默认或启动工作。
