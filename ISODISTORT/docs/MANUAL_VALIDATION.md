# 手工/长时验证脚本

实际脚本位于 `tests_dev/manual/`，不会被 pytest 自动收集；用于生成 CIF 与批量回归。

Method 1–4 的官网差分、候选归一、四格式导出、网页/终端/API 一致性与科学验收统一见 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)。Method 1/2 的 current-source live 报告必须与最终源码签名一致，旧签名不能冒充新结论。Method 3 官网集已经全部核定：EuAl4 20/20、35 个 embedding；NdNiO2 20/20、42 个 embedding。coupled 产品接入后已用当前源码从头重跑 40 组比较，结果为 13 exact + 27 affine-equivalent + 0 differences、0 errors；特殊-k与参数-k代表 coupled 行的完整位移模式数分别与官网 `2/2`、`22/22` 一致。Method 4 双母相本地/官网 24 例也已闭合，精确状态见验证报告。扩展样本、性能和并发阶段只在开发计划维护。

现有保存输出的只读审计，以及可选的真实 Method 3 iso/WSL 烟雾：

```text
python tests_dev/manual/validate_method_outputs.py
python tests_dev/manual/validate_method_outputs.py --live-method12
python tests_dev/manual/validate_method_outputs.py --live-method3
python tests_dev/manual/audit_method3_downloads.py "../output_compare/EuAl4 Parent.cif/官网/Method3" --json-output output/validation/method3_official_download_audit_eual4.json
python tests_dev/manual/audit_method3_downloads.py "../output_compare/NdNiO2 own.cif/官网/Method3" --json-output output/validation/method3_official_download_audit_ndnio2.json
python tests_dev/manual/audit_method3_embedding_routes.py --omit-raw-output
python tests_dev/manual/compare_method3_official_local.py --restart
python tests_dev/manual/audit_method3_basis_metamorphism.py
python tests_dev/manual/prepare_method4_inputs.py --force
python tests_dev/manual/validate_method4_local.py --update-manifest-status
python tests_dev/manual/audit_method4_official.py --parent "EuAl4 Parent.cif"
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

`audit_method3_basis_metamorphism.py` 不新增官网金标准；它把 40 个查询各改写成两种精确 GL(3,Z) representative，并先证明 basis+centering 仍生成同一 primitive lattice，再要求稳定 embedding ID 多重集完全相同。报告的 `source_signature` 覆盖审计器、manifest、搜索/仿射/格子核心和两份母相 CIF，源码变化后必须重跑。`validate_method4_local.py` 校验固定 daughter CIF 的 SHA-256、原始拟合系数、`As/Ap/normfactor`、Å residual、匹配、秩/条件数和失败用例；该报告只证明本地逆一致性。`audit_method4_official.py` 只读 `output_compare/`，核对冻结输入哈希、官网完整子群身份、模式/应变幅值、CIF/TOPAS/IsoVIZ 一致性，以及可唯一配对模式的 `normfactor` 和 `As/Ap`；机器报告写入 `output/validation/method4_official_audit.json`。

空间群产品查询默认启用 Stage A/B 和 coupled 稳定子交集证明；只有 `known_single_ir` 或经 exact fixed-space 可达性与交集见证双重证明的 `exact_fixed_space` 行可选择。Displacive 的逐物种选择会传入 Stage-B 表示；部分物种的共同位移保留为相对自由度，只有覆盖全部结构位点时才扣除三维全晶体刚体平移。`include_affine_only_diagnostics=True` 仍会额外显示不可达/未解析 Stage-A 行，但这些诊断行不得进入 Method 2 或导出。

`validate_method_outputs.py --live-method12` 使用独立的源码/输入签名与原子 checkpoint；当前源码长测结束前只引用 `live_method12_report.json` 的实时状态，不预填最终通过数。需要明确抛弃兼容断点时使用 `--live-method12-restart`。

```text
python tests_dev/manual/make_cifs_30.py
python tests_dev/manual/fetch_cod_cifs.py
python tests_dev/manual/run_web.py spotcheck
python tests_dev/manual/run_web.py m134
python tests_dev/manual/run_web.py method2_ld
python tests_dev/manual/run_batch.py cif30
python tests_dev/manual/run_batch.py external
python tests_dev/manual/run_batch.py terminal
```

请使用项目 `.venv` 中的 Python。CIF 数据仍在 `tests_dev/cifs_30/` 与 `cifs_external/`。

## 用户需要配置 / 放置的路径

本目录脚本**沿用**上层 ISODISTORT 配置，一般不必单独改路径：

| 项目 | 说明 |
| --- | --- |
| **母相 / 批量 CIF** | 默认读 `tests_dev/cifs_30/`、`cifs_external/`，以及仓库 `experiment_data/`（只读）。换样本时把 CIF 放进对应测试数据目录，或改脚本参数里的路径 |
| **iso / WSL / output** | 与正式程序相同：`config/settings.yaml`、`isobyu/`，见 [ISODISTORT/README.md](../README.md) §3.4 |
| **Method 3/4 官网批次** | 固定输入与下载要求见 `docs/manifests/`；运行断点和机器结果仍在 `output/validation/` |
| **网页 spotcheck** | 依赖本机已能启动 `main_web.py`（端口见 `runtime.web_port`） |

外部 CIF 来源见 [EXTERNAL_CIF_SOURCES.md](EXTERNAL_CIF_SOURCES.md)。跨项目总表见仓库根 [README.md](../../README.md)。
