# 手工/长时验证脚本

实际脚本位于 `tests/manual/`，不会被 pytest 自动收集；用于生成 CIF 与批量回归。

Method 1–4 的计划与待办见 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)；已完成批次、计数、
源码签名和证据边界的唯一来源是 [BUGFIX_VALIDATION_REPORT.md](BUGFIX_VALIDATION_REPORT.md)。
current-source 报告必须与最终源码签名一致，旧签名不能冒充新结论。2026-10-06 的长时
官网/live 矩阵是重构前的历史证据；2026-10-07 修复后的回归不能代替这些矩阵重算。

旧 4310 live 的 `75/125` 与 `50/125` 只属于当时源码；
`method1_4310_live_postfix_cif_precision_reanalysis_20261003.json` 的精度比较口径错误，
已撤回。任何脚本或人工汇报都不得把这两份旧结果当作当前源码通过证据。

以下命令都从 `CRIS/ISODISTORT/` 执行。`..\run_cris.ps1` 从仓库根实体 `.venv`
加载依赖，并绕过本机 OneDrive 内解释器启动 WSL 时的访问拒绝。现有保存输出的
只读审计，以及可选的真实 Method 3 iso/WSL 烟雾：

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
..\run_cris.ps1 tests/manual/prepare_method4_inputs.py --force
..\run_cris.ps1 tests/manual/validate_method4_local.py --json-output output/validation/method4_local_validation_20261006_current_v2.json
..\run_cris.ps1 tests/manual/validate_method4_local.py --parent "EuAl4 Parent.cif" --context gamma --case-id F02-invalid-lattice --json-output output/validation/strain_audit_f02_end_to_end.json
..\run_cris.ps1 tests/manual/audit_method4_official.py --local-report output/validation/method4_local_validation_20261006_current_v2.json --json-output output/validation/method4_official_audit_20261006_current_v2.json
..\run_cris.ps1 -m pytest tests -q --tb=line
```

命令中的历史日期用于定位已保存批次。重跑前给 report、checkpoint 和结果目录选用新的
时间戳，保留已有认证证据；任何源码、输入或官网证据变化后都要生成新报告，不能只复制
历史计数。报告文件 SHA-256、最终回归计数和适用范围见验证报告。

`diagnose_method2_export_preparation.py` 是 NdNiO2 参数-k 的逐候选只读长时诊断器；只需
候选清单时在命令末尾加 `--list-only`。完整运行会把每条模式的身份状态/原因与单候选
导出准备结果写入指定 `output/validation/` 目录，不发布生产导出。修复后的网页/API 批量
导出另会在部分成功或安全取消时，把 `export_candidate_status.json/.txt` 放在 ZIP/磁盘
批次根部。两类报告用途不同，不能用 25 秒的 `--list-only` 或原有两个代表例替代
NdNiO2 `Y a=1/3` 的 nmod=0 与 nmod≥1 全 48 候选验收。

审计按 CIF 内部候选身份配对，以精确整幺模变换证明等价子格，并将等价表示
归一到官网设置后比较；不会改写 `output_compare/`。Method 3 下载审计器只认
含 `Finish selecting...` 且带 `orderparam` 候选，或明确带
`There are no subgroups...` + `Try again` 的合法零候选结果表，并逐行核对
`(SG,symbol,basis,origin,s,i)` 与四类核心文件。CLI 默认绑定 Method 3 manifest，
同时核对 parent SG、目标 SG、输入 basis、Default/显式 centering 和 Types，
因此不能用另一查询的整套 HTML+导出来获得假通过。只要仍有误存/缺失结果页，
该脚本会写出完整 JSON 后以非零状态结束，这是“批次未齐”的预期信号。
Method 1/2/4 使用的共享 HTML resolver 不看 basename，而按页面结构、角色与内容
SHA-256 识别证据；在原候选目录内改成任意名称仍有效，只要后缀保持 `.html` 或
`.htm`（大小写不敏感）。需要的角色必须恰好一页：0 页或多页都 fail-closed。跨候选
目录移动会改变科学归属，改成非 HTML 后缀会退出清单，二者都不是安全的“仅改名”。
生产搜索、四格式导出与 ZIP 不读取这些归档；该规则只影响只读验证证据。
`audit_method3_embedding_routes.py` 对每条官网 embedding 运行精确的 ISO
`DISPLAY DIRECTION` 探针并保留命令、完整 stdout、版本及输入哈希；它只判定
是否存在 single-IR route，不伪造实际 coupled-IR 组合。当前 ISO 9.6.1 没有
直接枚举 Method 3 首屏 affine embedding 的命令，`DISPLAY DIRECTION` / `DISPLAY
ISOTROPY COUPLED` 只用于已知 embedding 之后的 route/COPL 第二阶段。
本地—官网脚本的原子 checkpoint 签名覆盖比较/审计脚本、manifest、`settings.yaml`、
全部 `backend/**/*.py`、母相 CIF、官网结果及四类核心文件、iso 二进制、
`const.dat`、配置数据目录下的 `data_*.txt`、Python/平台及关键依赖版本；缺失依赖同样会改变签名。
生成的 `i*.iso` 参数-k 缓存不能通过现有 API 稳定做内容签名，所以这类案例明确禁用跨进程 checkpoint 复用。`audit_method12_numeric_semantics.py` 的中断 checkpoint 处于 `pending_final_signature`，也必须用 `--restart` 从头运行；只有正常结束、完成终签名复核并标记为 `certified` 的 checkpoint 才可复用。其他审计是否支持中断续跑，以各脚本的 checkpoint 契约为准。当前 schema-5 比较器先冻结 manifest 与官网案例字节快照，签名和解析消费同一份内容，并在结束时重新签名；源码或证据漂移会使报告和 checkpoint 失效。
比较先保留逐字段 literal 结果，再用 `method3_affine_equivalence.py` 作独立的证明级
晶体学判定：以 `Fraction` 重建 Seitz 操作，在子群 primitive translation lattice
上取模，并在有限母/子平移商中检查母空间群共轭
`(Q,q)(R,t)(Q,q)⁻¹=(QRQ⁻¹,Qt+q-QRQ⁻¹q)`。只有能给出 exact/共轭
witness 的 basis/origin 差异才计为 `affine_equivalent`；验证错误或无法证明的项仍是差异。

默认只使用首屏结果表证明完整性的 authoritative 案例；只要存在 provisional 或 skipped 案例，route/比较 CLI 默认返回非零。`--allow-candidate-inventory` 与 `--accept-provisional` 仅用于将来诊断不完整批次，不能把候选目录升级成完整官网集合。当前 40 个案例和 77 条 embedding 全部 authoritative；route 审计为 `61 single_ir_exact + 16 coupled_ir_required`，0 条 indeterminate、0 个错误/跳过。该分类描述每条官网 embedding 是否可由一个 IR 单独稳定，不因产品实现变化而改写；当前默认空间群产品查询已覆盖 77/77。

`audit_method3_basis_metamorphism.py` 不新增官网金标准；它把 40 个查询各改写成两种精确 GL(3,Z) representative，并先证明 basis+centering 仍生成同一 primitive lattice，再要求稳定 embedding ID 多重集完全相同。报告的 `source_signature` 覆盖审计器、manifest、搜索/仿射/格子核心和两份母相 CIF，源码变化后必须重跑。`validate_method4_local.py` 校验固定 daughter CIF 的 SHA-256、原始拟合系数、`As/Ap/normfactor`、Å residual、匹配、秩/条件数和失败用例；不带筛选参数的当前报告已经完成双母相 24/24，仍只证明其签名对应的冻结矩阵，不外推到 4310 或跨晶系。带 EuAl4 F02 三个筛选参数的命令仍只生成单例报告。`audit_method4_official.py` 只读 `output_compare/`，核对冻结输入哈希、官网完整子群身份、模式/应变幅值、CIF/TOPAS/IsoVIZ 一致性，以及可唯一配对模式的 `normfactor` 和 `As/Ap`；当前机器报告为 `output/validation/method4_official_audit_20261006_current_v2.json`。

CIF 位移模型与 microscopic provenance 已接入生产 core 和四个 writer。模式查询先从目标
embedding 的 exact invariant direction 生成 primitive-integer `VALUE DIRECTION VECTOR,...`；
mapper、公共 I/O 合同及 writer 边界会重复核验直接 ISO 列、可证明的 BUSH 全域延拓、
逐原子父子映射、完整平移陪集、`(parent/child SG,B,q,primary IR/OPD)`、primary selector、
不可变结构快照和原子顺序。非零 strain + displacive 仍按官方 lattice-coordinate 独立约定
共用同一合同。当前合同只接受单物种满占位点；混合/部分占位、轴矢量、磁和占位模式不能
由 polar-vector 位移链的通过结果外推。

空间群产品查询默认启用 Stage A/B 和 coupled 稳定子交集证明；只有 `known_single_ir` 或经 exact fixed-space 可达性与交集见证双重证明的 `exact_fixed_space` 行可选择。Displacive 的逐物种选择会传入 Stage-B 表示；部分物种的共同位移保留为相对自由度，只有覆盖全部结构位点时才扣除三维全晶体刚体平移。`include_affine_only_diagnostics=True` 仍会额外显示不可达/未解析 Stage-A 行，但这些诊断行不得进入 Method 2 或导出。

`validate_method_outputs.py --live-method12` 使用独立的源码/输入签名与原子 checkpoint；
旧 `live_method12_report.json` 只属于历史签名。当前 Method 2 数值语义使用上面的
`audit_method12_numeric_semantics.py --live-current-source` 命令，schema-6 报告为两个
nmod=0 代表例与四种导出各 2/2 pass；需要明确抛弃兼容断点时使用 `--live-restart`。
`audit_method1_4310_live.py` 对 4310 的官网 Method 1 候选全集重新枚举并逐项计算，
核对完整候选身份、位移模式数、BUSH 子胞覆盖、模式基秩和零振幅 CIF 科学语义。
官网目录始终只读；生成 CIF、原子 checkpoint 和报告只写入
`output/validation/`。断点仅在审计脚本、验证器、最终 `backend` 源码、配置、
后端二进制/数据、母相、官网参考和运行时版本的联合签名完全相同时复用；
`--restart` 强制从头计算，`--max-candidates 1` 只用于明确标记为 partial 的单例烟雾。
N1+ 4D1 和整体 `z+0.1` 母相仍是定向根因回归；正式当前全量证据为 schema-3
`method1_4310_current_audit_20261006.json`。该报告状态为 `complete-passed`，运行前后签名均为
`96ff44cd307432ae2d8ca672afed6ff7c9e7dedbfc39ee07cfd71fefeaf94d12`，125 个官网与 live
候选一一配对，模式数、BUSH 覆盖和 CIF 语义均为 125/125，0
missing/extra/duplicate/failure；报告 SHA-256 为
`be7bee90ea47a03bb0aac56debb455a882b0d5d198c593130bd6b11db2f9fa6b`。该有限批次不能
证明任意晶体，也没有独立证明每个数值模式向量、归一化或 primary IR/OPD 分解唯一性。

`method1_4310_live_final6_20261004.json` 在运行期间发生源码变化，只能作为
115/125 的根因诊断；final7 是较早稳定源码快照，两者都不能覆盖当前 schema-3 报告。
以下命令用于同签名复验或源码变化后的重跑：

```powershell
..\run_cris.ps1 tests\manual\audit_method1_4310_live.py `
  --checkpoint output\validation\method1_4310_current_checkpoint_20261006.json `
  --report output\validation\method1_4310_current_audit_20261006.json `
  --dest-root output\validation\method1_4310_current_cifs_20261006 `
  --restart
```

作为历史环境事实，final7 执行时曾把虚拟环境放在 `%LOCALAPPDATA%\CRIS\venv-cris-py312`，并由根 `.venv`
junction 进入，以避开 OneDrive 内 Python 子进程的 `Wsl/E_ACCESSDENIED`；当时 doctor 为
24 pass、0 fail，Python → WSL 返回 `/home/devoutwang`。当前 `.venv` 已改为 CRIS 内的
实体目录，所有需要 WSL 的命令通过 `run_cris.ps1` 用 `.venv/pyvenv.cfg` 记录的外部基础
Python 建立进程，并用 `-S` 禁用其全局 site-packages，只从实体 `.venv` 加载第三方依赖。
启动器已实测 `spglib 2.7.0`、网页模块、
WSL home 和 `IsoDistort` 初始化。final7 使用 Python 3.12.5 运行 3174.65331 秒，报告 SHA-256 为
`383ed579c80d8f6252a2268e298ec1ee956b5f7fadaa2698bbf0053c68af9dc1`，并满足
`signature_unchanged_at_completion=true`。当前报告另见上述新签名与 SHA-256；任何重跑仍须先核对
结束签名稳定字段，再解释 summary。

```text
..\run_cris.ps1 tests/manual/make_cifs_30.py
..\run_cris.ps1 tests/manual/fetch_cod_cifs.py
..\run_cris.ps1 tests/manual/run_web.py spotcheck
..\run_cris.ps1 tests/manual/run_web.py m134
..\run_cris.ps1 tests/manual/run_web.py method2_ld
..\run_cris.ps1 tests/manual/run_batch.py cif30
..\run_cris.ps1 tests/manual/run_batch.py external
```

CIF 数据仍在 `tests/cifs_30/` 与 `cifs_external/`。

## 用户需要配置 / 放置的路径

本目录脚本**沿用**上层 ISODISTORT 配置，一般不必单独改路径：

| 项目 | 说明 |
| --- | --- |
| **母相 / 批量 CIF** | 默认读 `tests/cifs_30/`、`cifs_external/`，以及仓库 `experiment_data/`（只读）。换样本时把 CIF 放进对应测试数据目录，或改脚本参数里的路径 |
| **iso / WSL / output** | 与正式程序相同：`resources/config/settings.yaml`、`resources/isobyu/`，见 [ISODISTORT/README.md](../README.md) §3.3 |
| **Method 3/4 官网批次** | 固定输入与下载要求见 `docs/manifests/`；运行断点和机器结果仍在 `output/validation/` |
| **网页 spotcheck** | 依赖本机已能启动 `scripts/main_web.py`（端口见 `runtime.web_port`） |

外部 CIF 来源见 [EXTERNAL_CIF_SOURCES.md](EXTERNAL_CIF_SOURCES.md)。跨项目总表见仓库根 [README.md](../../README.md)。
