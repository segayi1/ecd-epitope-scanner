# Binder / 抗原共用工作台（2026-10-25）

两种分析对象共用已有 V4 主流程和七步工作台。此次将 binder 修复后的目标显示范围、伙伴链开关、定位与表示重建行为推广到普通抗原；保留分析对象字段及其科学语义。未重新评分或覆盖历史运行。

## 分区按分析链长度设置

复用 `config.display_sectoring_parameters`：分析链少于120 aa 时每区目标7个表面残基，120–150 aa 目标8个，更长目标12个。`joint_design_v4/engine.py` 根据实际分析记录数调用，与 antigen / binder 标签无关。显式用户目标优先；连通、直径、最小面积和唯一归属硬门继续使用已有实现，小蛋白4–12区是软范围。

修正 HTML / JS / 中英文说明仍要求8–15的过时限制，改为后端已支持的7–15；同步API旧测试的无效下界为6并验证7能够提交。没有引入按对象特判的新分区算法。

## 展示与算法边界

- `web/app.js`：有伙伴时统一使用完整输入结构，首次载入和重试走相同路径；所有非空结果都配置精确目标残基；第四步提供伙伴开关，离开第四步自动隐藏。
- `web/renderer.js`：使用固定版本PDBe Mol*既有 `getBundle` / `tryCreateComponent` 建立精确目标组件，避免展示同链ECD外残基；伙伴独立组件，表示重建后恢复可见性，复位聚焦精确目标。支持单体及原有整链调用。
- `web/index.html`更新脚本缓存版本；新离线站点由现有 `build_offline_site`打包同一前端。历史离线包保持原样。
- DiscoTope是主要评分差异：抗原使用抗体表位先验，binder禁用；SASA、ProCaliper、SPURS、FoldX、稳定性/空间门与预算共用。binder另外保留多链需指定目标链、可选参考序列映射、未知assay/膜方向与来源置信度规则。伙伴仅提供接触注释和展示。

## 验收

- 父目录项目venv运行8份相关测试：72 passed；收据 `test_output/unified_target_20261025/regression.xml`。包含长度119/120/150/151边界、7的API/JS传递、精确残基与插入码查询、单体/复合物、重载/复位、binder输入、普通V4预算和离线导出。
- 两份真实保存输入重新导出到独立验收目录：抗原 `05ecca7e4cc7`（125残基）和binder `article_repacked/f606506c8e4a`（77残基）。Chrome断网file浏览器均PASS：实际目标组件残基集合等于分析结果；第四步开关、离开隐藏、表示切换、同一viewer、七步及390px无横向溢出；pageerror与HTTP请求均为0。
- 浏览器收据与截图：`test_output/unified_target_20261025/browser_acceptance.json`、`antigen.png`、`binder.png`。JS语法、i18n JSON及本次前端/测试diff空白检查通过；项目记忆原有段落的行尾空白警告保留，未清理无关文本。
- 初轮45 passed / 1 failed来自API测试仍要求拒绝7；用户明确长度适配后将过时测试与UI同步，最终72项全部通过。保留一条已有Starlette/httpx弃用提示。

未跑全仓测试，未重启服务，未重新计算科学模型；真实浏览器验收针对新生成离线站点，未重新实测在线服务。
