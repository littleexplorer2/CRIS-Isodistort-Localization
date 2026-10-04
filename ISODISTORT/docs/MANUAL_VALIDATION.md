# 手工/长时验证脚本

实际脚本位于 `tests_dev/manual/`，不会被 pytest 自动收集；用于生成 CIF 与批量回归。

Method 1–4 的官网差分、候选归一、四格式导出、网页/终端/API 一致性与科学验收统一见 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)。Method 1/2 的 current-source live 报告必须与最终源码签名一致，旧签名不能冒充新结论。三种母相的 Method 1 官网下载静态门禁已通过；4310 的页头、物理轨道、逐物种轨道作用域和完整母相 setting 规范化已有定向回归，final7 又在签名 `0f8900bb5add71bc1de5e9d4c562142d8595539daf757dee074bf614c919a3b5` 下完成 125/125 候选 live：身份、模式数、BUSH 覆盖和零振幅 CIF 语义均为 125/125，0 失败。Method 3 官网集已经全部核定：EuAl4 20/20、35 个 embedding；NdNiO2 20/20、42 个 embedding。coupled 产品接入后已用对应签名源码从头重跑 40 组比较，结果为 13 exact + 27 affine-equivalent + 0 differences、0 errors；特殊-k与参数-k代表 coupled 行的完整位移模式数分别与官网 `2/2`、`22/22` 一致。Method 4 的官网归档与本地冻结快照 24 例已经闭合；当前源码只完成 EuAl4 F02 单例 WSL live，不能代替 24/24 重跑或跨晶系验证。精确状态见验证报告，扩展样本、性能和并发阶段只在开发计划维护。

旧 4310 live 的 `75/125` 与 `50/125` 只属于当时源码；
`method1_4310_live_postfix_cif_precision_reanalysis_20261003.json` 的精度比较口径错误，
已撤回。任何脚本或人工汇报都不得把这两份旧结果当作当前源码通过证据。

以下命令都从 `CRIS/ISODISTORT/` 执行。`..\run_cris.ps1` 从仓库根实体 `.venv`
加载依赖，并绕过本机 OneDrive 内解释器启动 WSL 时的访问拒绝。现有保存输出的
只读审计，以及可选的真实 Method 3 iso/WSL 烟雾：

```text
..\run_cris.ps1 tests_dev/manual/validate_method_outputs.py
..\run_cris.ps1 tests_dev/manual/validate_method_outputs.py --live-method12
..\run_cris.ps1 tests_dev/manual/validate_method_outputs.py --live-method3
..\run_cris.ps1 tests_dev/manual/audit_method3_downloads.py "../output_compare/EuAl4 Parent.cif/官网/Method3" --json-output output/validation/method3_official_download_audit_eual4.json
..\run_cris.ps1 tests_dev/manual/audit_method3_downloads.py "../output_compare/NdNiO2 own.cif/官网/Method3" --json-output output/validation/method3_official_download_audit_ndnio2.json
..\run_cris.ps1 tests_dev/manual/audit_method3_embedding_routes.py --omit-raw-output
..\run_cris.ps1 tests_dev/manual/compare_method3_official_local.py --restart
..\run_cris.ps1 tests_dev/manual/audit_method3_basis_metamorphism.py
..\run_cris.ps1 tests_dev/manual/prepare_method4_inputs.py --force
..\run_cris.ps1 tests_dev/manual/validate_method4_local.py --update-manifest-status
..\run_cris.ps1 tests_dev/manual/validate_method4_local.py --parent "EuAl4 Parent.cif" --context gamma --case-id F02-invalid-lattice --json-output output/validation/strain_audit_f02_end_to_end.json
..\run_cris.ps1 tests_dev/manual/audit_method4_official.py --parent "EuAl4 Parent.cif"
```

审计按 CIF 内部候选身份配对，以精确整幺模变换证明等价子格，并将等价表示
归一到官网设置后比较；不会改写 `output_compare/`。Method 3 下载审计器只认
含 `Finish selecting...` 且带 `orderparam` 候选，或明确带
`There are no subgroups...` + `Try again` 的合法零候选结果表，并逐行核对
`(SG,symbol,basis,origin,s,i)` 与四类核心文件。CLI 默认绑定 Method 3 manifest，
同时核对 parent SG、目标 SG、输入 basis、Default/显式 centering 和 Types，
因此不能用另一查询的整套 HTML+导出来获得假通过。只要仍有误存/缺失结果页，
该脚本会写出完整 JSON 后以非零状态结束，这是“批次未齐”的预期信号。
`audit_method3_embedding_routes.py` 对每条官网 embedding 运行精确的 ISO
`DISPLAY DIRECTION` 探针并保留命令、完整 stdout、版本及输入哈希；它只判定
是否存在 single-IR route，不伪造实际 coupled-IR 组合。当前 ISO 9.6.1 没有
直接枚举 Method 3 首屏 affine embedding 的命令，`DISPLAY DIRECTION` / `DISPLAY
ISOTROPY COUPLED` 只用于已知 embedding 之后的 route/COPL 第二阶段。
本地—官网脚本的原子 checkpoint 签名覆盖比较/审计脚本、manifest、`settings.yaml`、
全部 `isocore/**/*.py`、母相 CIF、官网结果及四类核心文件、iso 二进制、
`const.dat`、配置数据目录下的 `data_*.txt`、Python/平台及关键依赖版本；缺失依赖同样会改变签名。
生成的 `i*.iso` 参数-k 缓存不能通过现有 API 稳定做内容签名，所以这类案例明确禁用跨进程 checkpoint 复用；其他内容已签名案例中断后可安全续跑。要强制丢弃兼容断点并从头运行时使用 `--restart`。schema-3
比较先保留逐字段 literal 结果，再用 `method3_affine_equivalence.py` 作独立的证明级
晶体学判定：以 `Fraction` 重建 Seitz 操作，在子群 primitive translation lattice
上取模，并在有限母/子平移商中检查母空间群共轭
`(Q,q)(R,t)(Q,q)⁻¹=(QRQ⁻¹,Qt+q-QRQ⁻¹q)`。只有能给出 exact/共轭
witness 的 basis/origin 差异才计为 `affine_equivalent`；验证错误或无法证明的项仍是差异。

默认只使用首屏结果表证明完整性的 authoritative 案例；只要存在 provisional 或 skipped 案例，route/比较 CLI 默认返回非零。`--allow-candidate-inventory` 与 `--accept-provisional` 仅用于将来诊断不完整批次，不能把候选目录升级成完整官网集合。当前 40 个案例和 77 条 embedding 全部 authoritative；route 审计为 `61 single_ir_exact + 16 coupled_ir_required`，0 条 indeterminate、0 个错误/跳过。该分类描述每条官网 embedding 是否可由一个 IR 单独稳定，不因产品实现变化而改写；当前默认空间群产品查询已覆盖 77/77。

`audit_method3_basis_metamorphism.py` 不新增官网金标准；它把 40 个查询各改写成两种精确 GL(3,Z) representative，并先证明 basis+centering 仍生成同一 primitive lattice，再要求稳定 embedding ID 多重集完全相同。报告的 `source_signature` 覆盖审计器、manifest、搜索/仿射/格子核心和两份母相 CIF，源码变化后必须重跑。`validate_method4_local.py` 校验固定 daughter CIF 的 SHA-256、原始拟合系数、`As/Ap/normfactor`、Å residual、匹配、秩/条件数和失败用例；不带筛选参数的报告只证明其签名对应的冻结快照，带 EuAl4 F02 三个筛选参数的命令生成当前源码单例 WSL 报告，仍不能升级为 24/24 或跨晶系结论。`audit_method4_official.py` 只读 `output_compare/`，核对冻结输入哈希、官网完整子群身份、模式/应变幅值、CIF/TOPAS/IsoVIZ 一致性，以及可唯一配对模式的 `normfactor` 和 `As/Ap`；机器报告写入 `output/validation/method4_official_audit.json`。

CIF 位移模型与 microscopic provenance 已接入生产 core 和四个 writer。模式查询先从目标
embedding 的 exact invariant direction 生成 primitive-integer `VALUE DIRECTION VECTOR,...`；
mapper、公共 I/O 合同及 writer 边界会重复核验直接 ISO 列、可证明的 BUSH 全域延拓、
逐原子父子映射、完整平移陪集、`(parent/child SG,B,q,primary IR/OPD)`、primary selector、
不可变结构快照和原子顺序。非零 strain + displacive 仍按官方 lattice-coordinate 独立约定
共用同一合同。当前合同只接受单物种满占位点；混合/部分占位、轴矢量、磁和占位模式不能
由 polar-vector 位移链的通过结果外推。

空间群产品查询默认启用 Stage A/B 和 coupled 稳定子交集证明；只有 `known_single_ir` 或经 exact fixed-space 可达性与交集见证双重证明的 `exact_fixed_space` 行可选择。Displacive 的逐物种选择会传入 Stage-B 表示；部分物种的共同位移保留为相对自由度，只有覆盖全部结构位点时才扣除三维全晶体刚体平移。`include_affine_only_diagnostics=True` 仍会额外显示不可达/未解析 Stage-A 行，但这些诊断行不得进入 Method 2 或导出。

`validate_method_outputs.py --live-method12` 使用独立的源码/输入签名与原子 checkpoint；当前源码长测结束前只引用 `live_method12_report.json` 的实时状态，不预填最终通过数。需要明确抛弃兼容断点时使用 `--live-method12-restart`。
`audit_method1_4310_live.py` 对 4310 的官网 Method 1 候选全集重新枚举并逐项计算，
核对完整候选身份、位移模式数、BUSH 子胞覆盖、模式基秩和零振幅 CIF 科学语义。
官网目录始终只读；生成 CIF、原子 checkpoint 和报告只写入
`output/validation/`。断点仅在审计脚本、验证器、最终 `isocore` 源码、配置、
后端二进制/数据、母相、官网参考和运行时版本的联合签名完全相同时复用；
`--restart` 强制从头计算，`--max-candidates 1` 只用于明确标记为 partial 的单例烟雾。
N1+ 4D1 和整体 `z+0.1` 母相仍是定向根因回归；正式全量证据为
`method1_4310_live_final7_20261004.json`。该报告状态为 `complete-passed`，运行结束时
签名未变化，125 个官网与 live 候选一一配对，模式数、BUSH 覆盖和 CIF 语义均为
125/125，0 missing/extra/duplicate/failure。该有限批次不能证明任意晶体，也没有独立证明
每个数值模式向量和归一化。

`method1_4310_live_final6_20261004.json` 在运行期间发生源码变化，只能作为
115/125 的根因诊断，不能续作或引用为最终源码证据。final7 已在普通用户
PowerShell 中从 `CRIS/ISODISTORT/` 完成；以下命令用于同签名复验或源码变化后的重跑：

```powershell
..\run_cris.ps1 tests_dev\manual\audit_method1_4310_live.py `
  --checkpoint output\validation\method1_4310_live_final7_20261004.checkpoint.json `
  --report output\validation\method1_4310_live_final7_20261004.json `
  --dest-root output\validation\method1_4310_live_final7_20261004_generated `
  --restart
```

final7 执行时曾把虚拟环境放在 `%LOCALAPPDATA%\CRIS\venv-cris-py312`，并由根 `.venv`
junction 进入，以避开 OneDrive 内 Python 子进程的 `Wsl/E_ACCESSDENIED`；当时 doctor 为
24 pass、0 fail，Python → WSL 返回 `/home/devoutwang`。当前 `.venv` 已改为 CRIS 内的
实体目录，所有需要 WSL 的命令通过 `run_cris.ps1` 用 `.venv/pyvenv.cfg` 记录的外部基础
Python 建立进程，并用 `-S` 禁用其全局 site-packages，只从实体 `.venv` 加载第三方依赖。
启动器已实测 `spglib 2.7.0`、网页模块、
WSL home 和 `IsoDistort` 初始化。final7 使用 Python 3.12.5 运行 3174.65331 秒，报告 SHA-256 为
`383ed579c80d8f6252a2268e298ec1ee956b5f7fadaa2698bbf0053c68af9dc1`，并满足
`signature_unchanged_at_completion=true`。任何重跑仍须先核对该字段，再解释 summary。

```text
..\run_cris.ps1 tests_dev/manual/make_cifs_30.py
..\run_cris.ps1 tests_dev/manual/fetch_cod_cifs.py
..\run_cris.ps1 tests_dev/manual/run_web.py spotcheck
..\run_cris.ps1 tests_dev/manual/run_web.py m134
..\run_cris.ps1 tests_dev/manual/run_web.py method2_ld
..\run_cris.ps1 tests_dev/manual/run_batch.py cif30
..\run_cris.ps1 tests_dev/manual/run_batch.py external
..\run_cris.ps1 tests_dev/manual/run_batch.py terminal
```

CIF 数据仍在 `tests_dev/cifs_30/` 与 `cifs_external/`。

## 用户需要配置 / 放置的路径

本目录脚本**沿用**上层 ISODISTORT 配置，一般不必单独改路径：

| 项目 | 说明 |
| --- | --- |
| **母相 / 批量 CIF** | 默认读 `tests_dev/cifs_30/`、`cifs_external/`，以及仓库 `experiment_data/`（只读）。换样本时把 CIF 放进对应测试数据目录，或改脚本参数里的路径 |
| **iso / WSL / output** | 与正式程序相同：`config/settings.yaml`、`isobyu/`，见 [ISODISTORT/README.md](../README.md) §3.4 |
| **Method 3/4 官网批次** | 固定输入与下载要求见 `docs/manifests/`；运行断点和机器结果仍在 `output/validation/` |
| **网页 spotcheck** | 依赖本机已能启动 `main_web.py`（端口见 `runtime.web_port`） |

外部 CIF 来源见 [EXTERNAL_CIF_SOURCES.md](EXTERNAL_CIF_SOURCES.md)。跨项目总表见仓库根 [README.md](../../README.md)。
