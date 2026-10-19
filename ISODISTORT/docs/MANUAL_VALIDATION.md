# 当前手工与长时验证操作

下载步骤与路径只维护在 [DOWNLOAD_CHECKLIST.md](DOWNLOAD_CHECKLIST.md)。
检验任务见 [validation_checklist_20261019.json](manifests/validation_checklist_20261019.json)。
旧命令、输入签名、运行结果和 checkpoint 已在 [历史索引](history/README.md) 归档；不按旧进度续跑。

## 工程检查

从 `E:\CRIS` 运行；若执行策略禁止脚本，按下面的进程级 Bypass 调用，不修改系统策略。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest tests_repo -q --tb=line
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest ISODISTORT/tests -q --tb=line
.\.venv\Scripts\python.exe -m pytest E:/GD/tests_dev/test_isodistort_to_gd.py -q --tb=line
git -c core.safecrlf=false diff --check
```

生产模块回归按 `tests/backend`、`tests/application`、`tests/frontend` 存放，每模块一个文件。
单模块命令与清单见 [测试入口](../tests/README.md)。准备、下载、诊断和审计脚本直接放在 `tests/`；
原 `tests/manual/` 已移除。`tests/validation/` 是只读比较实现包，调用接口保持不变。
路径、工具内容或签名材料变化后，必须新建运行签名，不能恢复历史 checkpoint。

## 独立输入网格复核

此工具只读取输入清单、母相和通用 k 表，不读取官网结果；结论仅限有限代表输入。
使用新的报告路径，避免覆盖冻结证据。

```powershell
.\.venv\Scripts\python.exe -X utf8 ISODISTORT/tests/verify_parametric_k_grid.py --manifest ISODISTORT/docs/manifests/official_output_download_manifest.json --report ISODISTORT/output/validation/<新批次>/grid-report.json
```

## 文件清单与归档检查

主登记为 `docs/manifests/official_archive_reuse_20261019.json`，包含按内容识别的全部候选与文件指纹；
`official_result_list_stage_20261019.json` 分别记录列表、候选四件套和待检验状态。
一次性重建脚本和机器证据位于 `output/validation/rebuild_20261019/`，含已执行的归档/搬移步骤，
不可把有写操作的重建脚本当作日常审计重复执行。

`tests/register_official_archive_reuse.py` 仅复核三个已知 Method 2 查询的文件合同，
输出 `docs/manifests/official_archive_reuse_known_method2.json`。它不是全部查询发现器或完整清单，
不能用它覆盖主登记；新增官网材料应先按完整查询身份与当前主登记核对。

## 本地重跑与只读对照

`local_web_output_refresh_manifest.json` 当前为 `unfrozen`；历史 runtime 路径保留为无有效签名的
兼容入口，旧下载器拒绝启动是预期行为。新的检验任务清单已经冻结，但缺少官网列表、独立
输入或运行签名的任务仍未达到执行条件。

Method 3/4 旧机器清单已撤销可执行案例；完整旧文件只在历史归档中。先证明输入来源与合同，
再生成新批次，不得为恢复工具运行而把官网结果回填为输入。

本地独立求解、导出并冻结之后，才可以调用公共只读比较入口：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m tests.validation compare --local "<本地冻结CIF>" --reference "<官网CIF>" --json
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m tests.validation batch --local-dir "<本地目录>" --reference-dir "<官网目录>" --json
```

母相、方法、k/nmod、IR/OPD、完整子群身份、单位/坐标/归一化、实际外部读取、配置/源码/依赖
和输出指纹必须登记；每个数值比较须引用已有理论合同和适用容差。不存在输入、结果或有效
配对时记未执行/不支持/inconclusive，不把零案例或进程成功退出记成科学通过。

## 严格可信参考比较

普通模式报告结构/元数据语义一致；严格模式额外要求可信登记指纹及字节一致。
期望指纹从独立冻结的登记读取，不能临时对待比较文件算一个哈希再充作可信来源。
缺少可信哈希且未发现其他不匹配时返回 inconclusive；已有语义或字节差异仍为 fail，均退出码 1。
解析/调用错误退出码为 2。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m tests.validation compare --local "<本地冻结CIF>" --reference "<官网CIF>" --reference-sha256 "<独立登记SHA256>" --strict --json
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m tests.validation batch --local-dir "<本地目录>" --reference-dir "<官网目录>" --hash-manifest "<参考相对路径到可信SHA256的JSON>" --strict --json
```

`tests/validation/backend/provenance/` 提供共享文件快照、哈希、运行环境及原子 JSON 写入服务；
各审计定义自己的签名材料，尚非所有报告均使用原子写入。源码职责迁移后必须新建签名，
不得沿用旧 search_methods 路径、Method 4 v1 provenance 或未冻结 Method 3 checkpoint。
Method 3/4 准备器通过共享注册表覆盖生产 backend/application 源码和运行配置；
实现或配置变化后的旧 checkpoint/输入批次不得继续作为同一次运行复用。
公用科学计算入口为 `application.session`；外部 GD 读取 `displacement_snapshot()`。
浏览器人工验收应记录上传→搜索→模式→导出时 revision、job_id、取消结果及关闭后的资源释放；
失败动作后还要确认原候选/输入/模式未被改写，不以页面打开代替这项验收。
