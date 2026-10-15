# 系统重装后网站恢复与FoldX预算60

## 原因与范围

2026-10-12正式服务8765未监听。旧`.venv`引用`G:/anaconda`并开启system-site-packages；旧基础环境与新系统Python均无法导入FastAPI。临时目录变量使用短路径，当前执行环境无法写入系统临时目录，Python退回当前目录或报无可用临时目录。旧启动器没有依赖预检，冷启动超时阈值仅25秒且缺少固定启动日志。

## 修复

- 原环境保留为`.venv_pre_reinstall_20261012`并忽略Git追踪。使用新系统Python3.12.10创建独立`.venv`，按`requirements.lock`安装53依赖；`uv pip check`通过。具体安装版本见验收目录`restored_environment.txt`。
- 双击启动器保留原8765/loopback及父目录包解析，加入FastAPI/Uvicorn预检，后台隐藏启动并输出`workspace_data/server.stdout.log`与`server.stderr.log`；等待上限60秒。TEMP/TMP仅在启动器进程树中指向`workspace_data/runtime_tmp`，不修改用户或系统全局变量。重建提示改为独立环境/锁文件。
- `config.JOINT_V4_ML_PARAMETERS.foldx_confirmation.max_candidates`及engine的三个缺省回退值为60；每组合5次、FoldX2.0及SPURS0.5/1.0/分歧1.5门保持。共享设置同样用于专用路径的缺省预算，其选择逻辑保留；V3独立81预算未改。
- 当前方法和启动说明同步更新；历史42候选任务、基准性能记录与科学产物不重写。

## 验收

- 实际启动器在TTY下成功启动后台服务，冷启动约23秒、后续约2秒。结束启动器后服务保持。健康API为ok。
- 历史PD-L1任务18c0b8ee73fa保持COMPLETE；展示API返回10个专属首选；结构下载131489 bytes，SHA256见`site_api_receipt.json`。
- 实际正式浏览器载入结构画布、完成状态、10个分区首选与推荐详情；一个canvas，未捕获console error。截图为`workspace_data/system_restore_20261012/restored_site.jpg`。没有提交新科学设计任务。
- 预算专项3项通过；网站/展示API/投影/画布生命周期/V4回归58项通过，两批有重叠不求唯一总数。收据分别为`foldx_budget60.xml`和`site_regression.xml`。

## 未验收的模型运行

本轮恢复网站访问和预算配置，不以历史结果可加载代替新模型可运行。当前工具执行环境下，WSL Ubuntu-22.04已注册、VHD仍存在，但执行探测返回`Wsl/E_ACCESSDENIED`；ProCaliper旧独立解释器的导入探测25秒超时。未重装Linux模型环境、未修改VHD权限、未运行新的SPURS/DiscoTope/FoldX评分。具体探测保留在`model_runtime_probe.json`；这些结果不能单独证明模型文件损坏，新的科学评分链仍需恢复/验证。
