# 当前手工与长时验证操作

当前任务见 [开发计划](DEVELOPMENT_PLAN.md)，浏览器下载步骤见
[下载清单](DOWNLOAD_CHECKLIST.md)。历史命令、计数与修复过程见
[可恢复归档](history/README.md)。

## 运行前

### Agent 的统一测试与只读比较入口

以下命令从 **CRIS 根目录**执行，其余章节按各自说明切换目录：

```powershell
.\run_cris.ps1 -m pytest ISODISTORT/tests -q --tb=line
.\run_cris.ps1 -m pytest tests_repo -q --tb=line
.\run_cris.ps1 -m tests.validation compare --local "local.cif" --reference "reference.cif" --json
.\run_cris.ps1 -m tests.validation batch --local-dir "local" --reference-dir "reference" --json
```

单对输入可改用默认目录下的相对路径；批量以相同相对路径配对。
共同参数：`--lattice-tol`、`--coord-tol`、`--scalar-tol`、`--ignore-atom-order`、
`--strict`、`--json`。单对支持 `--reference-sha256`；批量支持 `--pattern` 与
`--hash-manifest`（相对路径到参考 SHA-256 的 JSON 映射）。未显式传目录时读取
`tests/validation/config/settings.yaml`。模式分解、超胞生成或修改参数必须在比较前
独立完成；参考答案只能进入只读比较环节。没有实际比较用例必须报告未验证。

### 下载输入网格的精确审计

从 CRIS 根目录运行（不调用 ISO、不读取官网计算结果）：

```powershell
.\.venv\Scripts\python.exe -X utf8 ISODISTORT/tests/manual/verify_parametric_k_grid.py
.\.venv\Scripts\python.exe -m pytest ISODISTORT/tests/test_parametric_k_grid.py -q
```

默认读取 `docs/manifests/official_output_download_manifest.json` 的输入白名单，重新核验
母相 SHA-256、标准 setting、全部标签、参数与坐标、符号一般小余群和同星碰撞。
默认报告写入 `output/validation/reciprocal_scope_20261017/grid-report.json`；后续修订用
`--manifest <输入清单> --report <新报告路径>`，不要覆盖本轮冻结证据。
该通过状态只证明声明的输入网格合同，不是官网结果或物理模式验收。

旧 `measure_method3_closure.py` 已停用，任何调用都返回 2，禁止恢复旧 checkpoint。
必须先完成开发计划中的独立输入枚举与心化/指数证明，再实现新的受检验范围工具。

### 长时科学审计

命令从 `CRIS/ISODISTORT/` 执行，经 `..\run_cris.ps1` 使用实体 .venv。
先核对冻结输入与运行签名；给 report、checkpoint、dest-root 选择新的批次目录。
以下历史日期是路径示例，不得覆盖原有证据。API 审计不能替代本轮真实浏览器下载。
固定 Method 4 daughter 已在 manifest 中按哈希引用，本轮不得重新生成它们。

`tests/manual/download_local_web_outputs.mjs` 的运行签名修订当前**没有生效项**：旧清单
`manifests/local_web_output_refresh_manifest.json` 已于 2026-10-16 置为 `superseded`，
脚本会据此拒绝运行，这是预期行为（见 [重置记录](manifests/validation_reset_20261016.json)）。
旧修订 `..._runtime_revision_20261011.json` 与 `..._20261015.json` 仅作历史记录，不得再用；
新官网集审计通过、新的冻结清单写入后，必须新增一份 runtime revision 才能恢复下载。
源码指纹覆盖 `backend/`、`application/`、`frontend/` 的实际实现，网页合同签名指向实际会话、
工作流与消息实现；修订记录的是该次 Windows 工作区的原始文件字节，源码或换行变化后须核对
并新增修订，不能把旧批次结果当成新代码的验收。

## 可用命令

```text
..\run_cris.ps1 tests/manual/validate_method_outputs.py
..\run_cris.ps1 tests/manual/validate_method_outputs.py --live-method12
..\run_cris.ps1 tests/manual/validate_method_outputs.py --live-method3
..\run_cris.ps1 tests/manual/audit_method3_downloads.py "../output_compare/EuAl4 Parent.cif/官网/Method3" --json-output output/validation/method3_official_download_audit_eual4.json
..\run_cris.ps1 tests/manual/audit_method3_downloads.py "../output_compare/NdNiO2 own.cif/官网/Method3" --json-output output/validation/method3_official_download_audit_ndnio2.json
..\run_cris.ps1 tests/manual/audit_method3_embedding_routes.py --omit-raw-output
..\run_cris.ps1 tests/manual/audit_method1_4310_live.py --checkpoint output/validation/method1_4310_current_checkpoint_20261006.json --report output/validation/method1_4310_current_audit_20261006.json --dest-root output/validation/method1_4310_current_cifs_20261006 --restart
..\run_cris.ps1 tests/manual/audit_method12_numeric_semantics.py --live-current-source --case EU_M2_LD1_C1_NMOD0 --case ND_M2_Y1_C1_NMOD0 --live-root output/validation/method2_numeric_live_current_source --live-report output/validation/method2_numeric_semantic_live_current.json --live-checkpoint output/validation/method2_numeric_semantic_live_current.checkpoint.json --live-restart
..\run_cris.ps1 tests/manual/diagnose_method2_export_preparation.py "../experiment_data/NdNiO2 own.cif" --output output/validation/ndnio2_export_audit_20261006 --k-point Y --k-parameter 1/3 --label Y_full
..\run_cris.ps1 tests/manual/compare_method3_official_local.py --report output/validation/method3_official_local_comparison_20261006_current.json --checkpoint output/validation/method3_official_local_checkpoint_20261006_current.json --restart
..\run_cris.ps1 tests/manual/audit_method3_basis_metamorphism.py
..\run_cris.ps1 tests/manual/validate_method4_local.py --json-output output/validation/method4_local_validation_20261006_current_v2.json
..\run_cris.ps1 tests/manual/validate_method4_local.py --parent "EuAl4 Parent.cif" --context gamma --case-id F02-invalid-lattice --json-output output/validation/strain_audit_f02_end_to_end.json
..\run_cris.ps1 tests/manual/audit_method4_official.py --local-report output/validation/method4_local_validation_20261006_current_v2.json --json-output output/validation/method4_official_audit_20261006_current_v2.json
..\run_cris.ps1 -m pytest tests -q --tb=line
```

## 验收与续跑约束

- NdNiO2 诊断命令加 `--list-only` 只枚举候选；完整诊断保存模式身份原因与导出准备
  状态，不发布浏览器文件。两种模式均不能替代真实网页 nmod=0/1 全候选验收。
- 四格式批次若部分成功或取消，保存实际生成的
  `export_candidate_status.json/.txt`；未生成的报告不可补造为浏览器原件。
- Method 3 官网结果表必须含候选，或明确的合法零候选响应；核对母相、目标 SG、
  basis、centering、Types 与完整候选身份。缺页、重复角色、provisional/skipped
  项不能升级为全量通过。诊断选项 `--allow-candidate-inventory`、
  `--accept-provisional` 不改变这个要求。
- HTML 按内容角色识别，须保留 .html/.htm 后缀与原候选目录归属。
- basis/origin 差异先记录字面结果，再以精确 Seitz/整数幺模 witness 判等价。
- checkpoint 只按各脚本合同续跑；`pending_final_signature` 不可复用，
  数值审计必须 restart；参数 k 生成缓存不能稳定签名时禁用跨进程复用。
  运行结束重验源码和证据签名，漂移则结果失效。
- Method 4 核对 daughter 哈希、拟合系数、As/Ap/normfactor、Å residual、
  原子匹配、秩/条件数及预期拒绝；筛选单例不能代表完整矩阵。
- 4310 `--max-candidates 1` 只作 partial 烟测；全量仍须逐候选核验身份、
  模式数、BUSH 覆盖、秩与 CIF 语义。数值向量与归一化另需独立证据。

## 扩展验证入口

```text
..\run_cris.ps1 tests/manual/run_web.py spotcheck
..\run_cris.ps1 tests/manual/run_web.py m134
..\run_cris.ps1 tests/manual/run_web.py method2_ld
..\run_cris.ps1 tests/manual/run_batch.py cif30
..\run_cris.ps1 tests/manual/run_batch.py external
```

样本入口为 `tests/cifs_30/`、`cifs_external/` 与只读 `../experiment_data/`。
运行配置唯一来源为 `resources/config/settings.yaml`；网页入口为 `main_web.py`。
新增外部样本前按 [EXTERNAL_CIF_SOURCES.md](EXTERNAL_CIF_SOURCES.md) 确定来源与许可。
