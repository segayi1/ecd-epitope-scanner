# 网页内容获取能力核查(2026-10-16)

> **时点说明(2026-10-16 补注)**:本文件核查动作发生在 **2026-10-16**,属于对 2026-10-13 那次调研的通道补核;标题日期指核查当日,不是调研当日。正文引用本文件作为"检索只走 PubMed E-utilities"的依据时,请同时看到 `docs/extension_directions.md` 的调研时点(2026-10-13)与补核时点(2026-10-16)。此处曾被认为与正文时序不符,现予澄清:两个日期各自正确,指向两件事。

## 结论

**不需要装插件。** 本机已有两条可用通道,文献调研不受阻。

| 通道 | 状态 | 实测证据 |
|---|---|---|
| `web_fetch`(内置工具) | **可用** | `GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=36366745&retmode=json` → HTTP 200,返回完整元数据(标题/期刊/年份/DOI/PMCID) |
| Shell 内 `Invoke-WebRequest`(pwsh) | **可用** | `GET .../esearch.fcgi?db=pubmed&term=epitope+prediction&retmode=json&retmax=2` → HTTP 200,`count=18240`,返回 idlist |
| `web_search`(内置工具) | **不可用** | `Error: DeepSeek search has no API key for "DEEPSEEK_API_KEY"; ... the web Models page writes it` |

## 早先「TLS 失败」的重新判读

会话早期 `Invoke-WebRequest` 报过传输层失败,当时据此改走 `web_fetch`。现在重测同一类请求 HTTP 200 通过。**结论:那次失败是暂时的,不是通道能力缺失**,不应作为「网络不通」的证据保留。真实且稳定复现的缺口只有 `web_search` 一条。

## `web_search` 的三条修复路径(工具自己给出的原文)

1. 经 credentials service 存 `DEEPSEEK_API_KEY`(web Models 页面写入);
2. 在启动环境里 export;
3. 在 `web-search-deepseek` 配置里写 literal `apiKey`;
4. 或:用 DeepSeek Account 模型开对话,则改用账号登录态搜索。

**我没有采取任何一条**:凭据属于用户,不该由我代填或代造。此项交人决定。

## 为什么不需要装插件(已排查)

市场快照 `~/.dsh/profiles/desktop/.dsh-market/discovery-compatibility-v1.json` 中确有几个名称相关的插件:

| 插件 | 版本 | 判断 |
|---|---|---|
| `dsh-free-search` | 0.7.2 | 提供搜索工具,但**未安装**;且它是又一入口,不解决 PubMed 检索本身 |
| `@anysearch/anysearch-dsh` | 0.1.6 | 依赖 `@deepseek-ai/dsh-web`,**未安装** |
| `dsh-builtin-browser` | 0.4.1 | 浏览器控制,**未安装** |
| `@a9i5k4/dsh-literature` | 0.3.6 | 文献向,**未安装** |

已安装 bundle 只有:`@agents-anywhere/dsh-bridge-next` 2.0.2、`@furongjun1999/dsh-memory` 0.8.1、`billion-context` 0.1.197、`clearai-dsh` 0.5.5、`dsh-better-sidebar` 0.26.1、`dsh-codex-migrate` 2.0.4、`dsh-deja` 0.22.1、`dsh-univer-office` 0.3.7、`dsh-workbuddy-connect` 0.7.6、`dshmarket` 1.66.20。

浏览器扩展 `~/.dsh/browser-extension/manifest.json` v0.1.4(权限含 `tabs`/`activeTab`/`scripting`、host `http(s)://*/*`)——但这是**扩展侧**的能力,与我的工具面无关,不能假定我能调用它。

**关键判断**:PubMed E-utilities 是结构化检索 API,比通用网页搜索更精确、可引用、可复现(每条命中都有 PMID/DOI)。在已有 `web_fetch` + shell HTTP 的情况下,装一个通用搜索插件的收益远小于引入新依赖与兼容性风险的成本。

## 附:检索接口的两个坑(实测)

- `sort=date` **不被接受**:返回 warning `Unknown sort schema date ignored`,静默退化成 best-match。要按时间排序必须用 `sort=pub_date`。
- 同一条检索式的关键词会被自动改写,例如 `epitope` 被映射到 `"difluprednate"[Supplementary Concept]`,导致召回噪声(如把电化学传感器论文拉进来)。**必须人工筛**,不能直接采信命中。
