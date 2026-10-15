# SPURS与DiscoTope实际版本识别

| 工具 | 软件包版本 | 方法/模型系列 | 当前源码修订 |
|---|---|---|---|
| SPURS单点 | spurs 1.0.0 | SPURS | 70cbb99af5cdf1f6da392389aaf93fd245610efb |
| SPURS组合 | spurs 1.0.0 | SPURS multi，无独立模型语义版本记录 | 同上 |
| DiscoTope | discotope3 1.0 | DiscoTope-3.0 | vendor无Git元数据，未知；新版worker记录源码文件指纹 |

## 依据与边界

- 按冻结任务6aef485cc7c1命令使用Ubuntu-22.04内spurs_v4/discotope_v4两个实际Python环境，读取importlib.metadata版本、editable direct_url和模块路径，核对本地setup.py/README及vendor Git。
- SPURS当前只记录requirements.inference.txt有跟踪修改；包版本不代表仓库完全未修改。DiscoTope没有vendor Git，不能将工作台父仓库提交当成其修订。
- 当前108个独立模型文件与冻结任务model_files_sha256全部匹配；SPURS单点、组合checkpoint不同，分别保留原权重摘要。
- 原运行没有记录软件版本或代码提交，不能追溯证明当时源码完全相同。新版投影明确标为“当前环境补充识别，非历史运行时版本记录”，原worker_runs完整保留；包身份及权重匹配不构成历史源码冻结。
- PyTorch 2.5.1+cu124/CUDA12.4是计算运行时，不是SPURS或DiscoTope的软件版本。

## 更新与交付

- SPURS/DiscoTope worker在推理结果直接记录包版本、Git提交/跟踪修改、源码文件指纹；DiscoTope另记模型系列DiscoTope-3.0。原模型权重和评分逻辑不变。
- engine将实际worker版本写入tool_inventory及每个推荐组合来源，正式网页与离线共用展示。DiscoTope模型3.0和安装包1.0分别标注，哈希/路径仍放折叠详情。
- 保存评分重选支持提供绑定冻结manifest和权重的版本识别审计；不得绑定其他任务的版本审计。JSON/CSV/ZIP均保留补充识别说明。
- 正式新任务f2b213d9b0ad：http://127.0.0.1:8765/?job=f2b213d9b0ad。评分来源6aef，新增模型调用0；仍9区7通过2待审查、8正式主构建体，历史265条任务保留。

## 证据

- workspace_data/model_version_identification_20261009/installed_identity.json：环境身份、108文件哈希逐项比较与原运行限制。
- 同目录spurs_worker_identity.json/discotope3_worker_identity.json：实际环境执行新版身份函数，版本/提交/源码指纹自检通过，不导入推理模型或调用GPU评分。
- 同目录export_verification.json：API/JSON/CSV/ZIP、下载字节、清单及冻结源一致。
- 界面5项聚焦回归通过；正式及离线真实Chrome1500/390像素共4项验收通过，三条工具来源版本/补充身份说明准确，无模拟路由、页面异常或溢出。收据与截图见workspace_data/model_version_browser_20261009_final/browser_verification.json。

没有安装、升级、重算模型或覆盖历史证据；原运行代码提交仍无法追溯证明。
