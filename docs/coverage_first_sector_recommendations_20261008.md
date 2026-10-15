# 覆盖优先与逐区建议验收

主面板以合格表面 `geodesic_0` 新增面积优先，之后信息、风险及ID；原 `NO_POSITIVE_NET_INFORMATION_GAIN`不再否决正覆盖增益。每区按实际突变残基匹配建议，原anchor只是生成窗口。门槛未改，待审查者不进入正式主面板及Stage2。主预算独立，跨区组合可共用。

固定CTLA4来源作业 `5068b62450ec`，新输出 `workspace_data/ctla4_coverage_first_regions_20261008_final`：8个通过资格门的主组合，9条逐区建议（7 PASSED、2 NEAR_THRESHOLD），直接残基面积覆盖0.256274。此数值不是扰动阳性面积或实验结合覆盖。

- VS008：M134W+Y140W，SPURS组合值0.648233高于0.5；FoldX缺失，保留 `COMBINATION_REJECT`。
- VS009：E155D+S50I，模型差值1.978935高于1.5，保留 `MODEL_DISAGREEMENT`；不解释为已证实不稳定。

近门槛排序按独立未过门数、缺失证据数、已知归一化超限量、覆盖、风险、ID。缺失FoldX不是零能量。诊断保存模型值/阈值/超限及missing标志。

验证：97项相关回归通过；完整PD-L1流水线补验4项通过（183.32秒），检查逐区数量、原因、JSON及ZIP逐字节一致。初次新增断言误用panel字段，修正为selected_panel后重跑通过，生产输出未因此改变。真实Chrome显示9区/8主组合/2待审查，点击对应结构和序列正常，无page_errors；证据 `workspace_data/ctla4_coverage_first_regions_20261008_browser_final`。逐区JSON/CSV、analysis、selected panel和ZIP的候选/状态/原因一致。

选择重放新增科学调用0，9个来源文件hash保持不变。早期重放因序列化键类型和Stage2格式问题失败的输出保留，没有覆盖历史。旧运行不自动重写。2026-10-08尝试确认无活跃作业/子进程后重启服务，自动审批拒绝，原因为blocked by policy；命令未执行，运行后端尚未加载新源码，新策略待正常服务重启生效。分区内组合内部冗余尚未解决，后续方案见 `intra_sector_design_optimization_20261008.md`。
