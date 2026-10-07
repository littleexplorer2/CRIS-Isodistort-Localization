# ISODISTORT 本地化开发计划

本文件是 `ISODISTORT/` 唯一的开发目标与后续验证计划。开发规则、修改边界和说明性文件职责见 [agent.md](../agent.md)；已经完成的修复、实测结果和机器报告索引见 [BUGFIX_VALIDATION_REPORT.md](BUGFIX_VALIDATION_REPORT.md)。

实现任何目标时都禁止按单个母相、IR 或 OPD 硬编码官网答案。应修复可推广的晶体学算法，并把无法实现的能力明确列为限制。

## 分区定义

本文件固定保留以下两个常驻分区：

1. **待修的程序漏洞**：已确认或高度可疑、且有用户可见症状或明确触发条件的现有代码缺陷。每条包含症状、触发条件、证据、根因链、影响与验收判据。
2. **待执行的优化计划**：尚未实施的能力、覆盖扩展、性能与工程改进。每条包含目标与验收标准。

除此之外的内容（当前顺序与状态、新会话续接执行顺序、已完成目标的历史记录）不属于上述两个分区，保持原样、不做硬分类。

现有内容分类映射：

| 现有内容 | 分区 |
| --- | --- |
| §待修的程序漏洞 | 待修的程序漏洞 |
| §待执行的优化计划 | 待执行的优化计划 |
| §新会话续接执行顺序 | 未归档（跨分区调度顺序，保持原样） |
| §当前顺序与状态、目标 1/2、目标 3 已完成部分、目标 4 已完成矩阵、目标 5 §5.1 与 §5.2.1 | 未归档（进度与证据记录，保持原样） |
| 目标 3 §3.3 第 4 项、目标 4 前置与收尾、目标 5 §5.2 第 2–5 项、§保留的产品限制、§后续扩展验收 | 待执行的优化计划（正文保留在原处，此处登记索引） |

审计与重跑协议：本计划的任何审计结论都必须记录源码/输入/运行时签名、只读边界和机器报告位置；不能由单母相、单 IR 或单候选的通过外推到其他晶体或全矩阵。

## 待修的程序漏洞

条目按严重程度排序。P0 = 当前源码在受支持输入下无法完成交付路径；P1 = 会输出错误或不可复核内容的风险；P2 = 可用性与可诊断性缺陷。

### [P0] BUG-001 Method 2 参数 k 批量导出被整批拒绝（28/48 候选，NdNiO2 Y `a=1/3`）

- **修复状态（2026-10-06）**：通用实现已完成，真实 48 候选验收仍待长时重跑。writer
  门禁现只检查候选实际写出的列；身份未解析候选单独跳过并写入根级
  `export_candidate_status.json/.txt`，已验证候选继续按原子批次发布。其它计算/writer
  异常仍整批 fail-closed。本条保留在待修分区，直到 nmod=0 与 nmod≥1 的 48 候选矩阵
  按下述验收判据完成，不能用单元测试或原有 2/2 代表例提前关闭。
- **症状**：网页 Method 2 四格式批量导出在运行约两小时后抛出
  `ExportPreparationError: Method 2 batch export preparation failed for 28 candidate(s): …`，
  每个候选的错误文本是
  `ValueError: displacement modes lack verified ISO microscopic identities: <mode keys>`；
  由于批量发布是原子的，成功的候选也一并被丢弃，目标目录里 0 个新文件。
- **触发条件**：母相 `experiment_data/NdNiO2 own.cif`（`P4/mmm #123`），Method 2 指定
  k 点 `Y`、参数 `a=1/3`，Types = strain + displacive。本地枚举得到 48 个参数 k 候选
  （Y1/Y2/Y3/Y4 × P1–P4、C1–C5、S1/S2、4D1），用户运行的失败集合是其中的 28 个
  （其余候选的导出结果被整批回滚，无法从该次运行分别判定）。
- **证据**：用户 2026-10-06 网页运行的完整错误文本；本机只读枚举复现
  （`ISODISTORT/output/validation/ndnio2_export_audit_20261006/`）得到 48 个参数 k 候选，
  与用户报错中的候选身份、basis、origin 逐项一致；[output_compare/NdNiO2 own.cif/现有网页版交互/Method2/](<../../output_compare/NdNiO2 own.cif/现有网页版交互/Method2/>)
  下保存的 9/24 产物在同一候选身份上存在 4 个非空文件（`subgroup.cif`、
  `Complete modes details.txt`、`data.isoviz`、`topas.str`），其
  `Complete modes details.txt` 记录同一 `Y1 P2 / 51 Pmma / basis {(0,2,0),(0,0,1),(3,0,0)} /
  origin (0.5,0.5,0) / s=6 i=12` 与 10 个带 ISO 站点 IR 标签的位移模式。
- **复现与定位入口**：只读诊断器
  `tests/manual/diagnose_method2_export_preparation.py`（本轮新增，只写报告，不改程序）。
  枚举与逐候选模式/导出准备：
  `.\run_cris.ps1 ISODISTORT\tests\manual\diagnose_method2_export_preparation.py "experiment_data\NdNiO2 own.cif" --output "ISODISTORT\output\validation\ndnio2_export_audit_20261006" --k-point Y --k-parameter 1/3 --label Y_full`；
  只枚举约 25 秒，同命令加 `--list-only`。报告逐候选记录
  `mode_identity_status_counts`、每条模式的 `status/source/reason`、
  `export_preparation_ok|failed` 与异常原文。注意：该诊断器必须经 `run_cris.ps1`
  或已安装硬链接的 `.venv` 启动，否则 OneDrive 内的 venv stub 会触发
  `Wsl/E_ACCESSDENIED`。
- **本轮实测的规模信号**：`--list-only` 枚举 48 个候选约 25 秒；同一命令去掉
  `--list-only` 后，在本机运行超过 50 分钟仍未完成第一个候选（进程持续占满单核，
  且这 50 分钟内没有新的 ISO 子进程，说明主要成本落在 `_exact_child_operations` /
  `_exact_displacive_child` / `_displacive_export_data_for_subgroup` 的精确有理运算，
  不在后端查询），期间诊断器与网页都不会输出任何逐候选进度。该信号与 BUG-004 的
  成本/可观测性问题一致，但不改变 BUG-001 的根因链。
- **代码路径（根因链）**：
  以下 1–6 描述修复前链路，便于与用户原始报错和旧报告对应：
  1. `features/method2/superspace.py::compute_parametric_modes` 先用 SMODES/Gram–Schmidt
     构造完整固定子空间，并**按设计**把每条列标成
     `ModeIdentity.unresolved(reason="smodes_basis_has_no_iso_microscopic_identity")`
     （约 3634–3704 行）。
  2. 随后尝试用 ISO 微观列整空间替换该基：`iso.list_invariant_directions` →
     `iso.calc_parametric_microscopic_distortion_modes` → `_canonicalize_rootless_supplement`
     → `_replace_rootless_supplement_basis`（约 3741–3826、1474–1681 行）。
  3. 任一步抛出 `ValueError / IsodistortError / LinAlgError` 时，函数只把原因写进
     `ParametricModeResult.note`（`"parameter-k microscopic canonicalization rejected: …"`），
     模式列表**保持全部 unresolved** 并正常返回，不报错、不降级、不标记候选。
  4. `backend/api/core_api.py::_displacive_export_data_for_subgroup`（约 3489–3502 行）
     要求当前会话的**每一个**位移模式都满足
     `mode_identity.status == "verified" and source == "iso_microscopic"` 且
     `microscopic_provenance is not None`，否则抛出上述 `ValueError`。该门槛针对整个
     位移缓存，而不是该候选实际用到的、或该 writer 实际写出的模式子集。
  5. `_collect_export_specs`（约 3965–4034 行）按候选汇总异常，最后一次性抛出
     `ExportPreparationError`（`backend/api/core_api.py` 333–380 行）；`_render_export_batch`
     与 `_publish_export_batch` 都要求全批次通过，因此 20 个成功候选也被丢弃。
  6. 网页边界 `frontend/web/server.py`（1087–1097 行）把该异常原样转成 JSON 500 返回。
- **影响**：受支持母相 + 受支持 Method 2 路径 + 受支持四格式导出，在真实晶体上
  无法交付；已投入的两小时计算全部作废，且失败原因不是用户可操作的输入错误。
- **验收判据**：
  1. NdNiO2 `Y a=1/3`（nmod=0 与 nmod≥1 各一轮）48/48 候选要么写出四格式、要么以
     “该候选的模式身份未解析”这一显式、可操作的原因单独列出；不允许一个候选的
     未解析身份阻止其他候选发布。
  2. 每个未解析候选的错误文本必须同时给出：未解析的物理轨道、该轨道失败的判据
     （`iso_microscopic_block_missing` / `iso_microscopic_output_unavailable` /
     rank 或 principal-angle 不匹配 / `embedding translation lattice is not a parent
     sublattice` / `_direction_symbol_for_invariant` 唯一性失败）以及被查询的
     `(parent SG, k, KVALUE, IR, direction selector)`。
  3. 修复必须落在通用生成算法上，不得按母相、IR 或 OPD 写特例。

### [P1] BUG-002 参数 k 模式身份的失败原因在交付边界被吞掉

- **修复状态（2026-10-06）**：代码与定向回归已完成；全 48 候选的真实失败判据分布仍随
  BUG-001 长时验收确认。canonicalization 的实际异常类别/文本与
  `(parent SG, k, KVALUE, IR, direction selector)` 现写回每一条仍未解析的
  `ModeIdentity.reason`，API、网页警告、导出异常及机器报告共用这一事实源。
- **症状**：用户只看到 “lack verified ISO microscopic identities” 与一长串
  `…__unresolved(a)` 键，看不到真正的拒绝判据。
- **根因**：唯一的诊断载体是 `ParametricModeResult.note`，而
  `search_methods.py`（1204–1245 行）与 `core_api.py::_calc_displacive_modes`
  把它放进元数据（`coupled_mode_note` / `parametric_note`），既不进入
  `distortion_modes`，也不进入导出异常；`_displacive_export_data_for_subgroup`
  只读 `mode_identity.reason`，而这里是统一的
  `smodes_basis_has_no_iso_microscopic_identity`。
- **影响**：无法区分“ISO 查询失败”“方向解析不唯一”“整空间证明被拒”三类根因，
  使 debug 与回归定位成本极高（配合 BUG-004 的两小时单次运行，一轮定位就要数小时）。
- **验收判据**：失败候选的 `mode_identity.reason`（或等价的失败记录）必须携带
  具体判据与查询上下文，并在网页错误文本、`ExportPreparationError.message` 和
  机器报告中保持一致。

### [P1] BUG-003 参数 k 候选在模式身份未解析时仍可被选择与批量导出

- **修复状态（2026-10-06）**：代码与定向回归已完成。候选尚未计算模式时保持
  `not_checked`；一旦同一 writer 门禁判定为 `unresolved`，API 返回
  `export_ready=false`，网页显示具体警告并把该行改为 `selectable=false`。批量路径无需
  预先点选，直接按相同门禁逐候选跳过并写机器报告。
- **症状**：`_filter_subgroups_for_search`（`core_api.py` 992–1032 行）按 IR 活性
  过滤参数 k 子群，不检查该子群的模式身份链是否可解析；这些候选照常出现在网页
  Method 2 结果表里并可被勾选，直到批量导出才失败。
- **根因**：可选择性（`selectable`）与身份可解析性是两条独立逻辑，前者未接入后者。
- **影响**：违反 [agent.md](../agent.md) 第 5 条“诊断候选必须 `selectable=false`”，
  用户承担了整批失败的全部代价。
- **验收判据**：身份未解析的候选必须显式标记为不可交付（或至少在 UI 与批量导出
  前给出同样显式的警告），且该标记必须由同一身份判据驱动，不得复制成第二套规则。

### [P1] BUG-006 非 Γ 母胞模式的 fallback 只生成 cosine，无法表达独立 sine/OPD 分量

- **修复状态（2026-10-07）**：静默 cosine 合成已删除。三条母胞单实列路径共用同一
  提升门禁：一般非自共轭 k 拒绝；自共轭模式在上传坐标系的真实母相平移格上验证
  精确 ±1 字符、中心化等价位点的位移协变性和子胞周期性；`opd_direction` 缺少
  来源契约时明确报错。已有完整子胞列路径继续使用。48 项正反例通过，包括极小
  非零波矢、显式 Γ 的中心化协变、合法子晶格与 `k=None` 缩胞场周期性。以下症状
  与根因记录修复前行为；成对实列/复列及 k-star/OPD 来源合同的能力扩展仍按验收判据待做。
- **症状与最小反例**：`features/method4/distortion_engine.py` 的
  `generate_single_mode`、`generate_modes` 与 `lift_mode_displacements` 都把母胞位移乘以
  `cos(2πk·t)`；`opd_direction` 虽声明为“选择相位符号”，实际未参与计算。对单原子
  `q=(0,0,1/4)`、`1×1×4` 超胞，`opd_direction=[1,0]` 与 `[0,1]` 得到相同的
  `[+1,0,-1,0]` 序列，独立的 sine 序列 `[0,+1,0,-1]` 不可达。
- **触发路径**：产品通常优先使用已经证明的完整子胞列 `mode_displacements_sc`；但其为空时，
  `backend/api/core_api.py::_lifted_mode_displacements` 与 `_method4_reference_modes` 仍会调用
  `lift_mode_displacements`，因此该 fallback 不是死代码。
- **根因与边界**：一般非 Γ 实位移是复 Bloch 振幅（或等价的 cosine/sine 两个实分量）与
  k-star/OPD 系数的组合。当前数据契约只有一个母胞实位移向量、`k_vector` 和未消费的
  `opd_direction`，不足以唯一重建相位原点、共轭臂及各分量系数；不能通过猜测一个 sine
  公式或按 IR 名称写特例修补。
- **影响**：当完整子胞列缺失而 fallback 被调用时，Method 4 参考子空间可能降秩，导致
  合法 sine/相移畸变无法分解，或把不同 OPD 分量当作同一模式。
- **验收判据**：扩展唯一模式事实源，使每个非 Γ 模式携带经来源证明的复列/成对实列、
  k-star 臂、相位原点与 OPD 系数（或直接携带完整子胞列）；生成、分解和四种 writer
  消费同一契约。缺少这些字段时必须 fail-closed，不得从母胞单列合成。至少覆盖
  `q=1/4` 的正交 cosine/sine、区边界自共轭 k、中心化点阵及原点平移协变性。

### [P2] BUG-004 批量导出按候选串行重算模式，缺少成本预估与中断

- **修复状态（2026-10-06，部分完成）**：已加入逐候选进度、完成数和基于已耗时均值的
  剩余时间估计；取消在当前候选完成后的安全边界生效，已完成候选与未运行清单一起原子
  发布；网页只在锁内创建独立导出快照，长时计算不再占用 `_SESSION_LOCK`。逐候选串行
  重算和大型超胞的精确有理运算成本尚未优化，因此本条仍开放。
- **症状**：48 个候选各自触发一次完整的模式计算与精确子胞构造，一次批量导出以小时计；
  期间既没有逐候选进度，也没有成本预估，用户只能等待，失败后全部作废。
- **根因**：`_collect_export_specs` 对缓存不匹配的候选逐个调用 `search_method_2`，
  再逐个进入 `_spec_for_subgroup` → `_displacive_export_data_for_subgroup`
  （精确 `_exact_child_operations`、`_exact_displacive_child`、逐原子父子映射与来源
  核验都对全超胞做有理运算），并且整个过程在 `_SESSION_LOCK` 内执行
  （`frontend/web/server.py` 491–516 行）。
- **本轮实测**：同一母相同一查询下，`--list-only` 枚举 48 个候选约 25 秒；加入逐候选
  模式与导出准备后，单次进程运行 50 分钟仍停在第一个候选上，且该时段无新的 ISO
  子进程，说明瓶颈是 Python 侧的精确有理运算而不是后端查询。
- **影响**：单点失败即浪费全部运行时间；同一母相的其他标签页在整个批量导出期间
  无法提交任何 mutation。
- **验收判据**：给出按候选的数量级成本预估与逐候选进度；支持安全中断并保留已完成
  候选的机器记录；批量导出不再长时间独占会话锁。

### [P2] BUG-005 大型整数超胞缺少覆盖与秩/条件数证据，固定空间证明退化为整批拒绝

- **症状**：失败集中在 `subgroup_index ≥ 36`（3×12、6×12、12×12 型超胞）的候选，
  而 4D1 与低指数子群可以通过。
- **怀疑**：`_replace_rootless_supplement_basis` 要求逐物理轨道
  `rank(reference) == rank(candidate) == 列数` 且 principal-angle 匹配，任一不满足即
  整空间拒绝；大超胞下 SMODES 投影基与 ISO 微观列的顺序、跨轨道归属或稀疏域覆盖
  更容易出现不一致，而当前没有针对“大超胞”这一维度的定向覆盖证据。
- **状态**：待实测确认（按 BUG-001 的验收判据逐候选记录失败判据后再定根因），
  在确认前不得当作已定根因。

## 待执行的优化计划

本分区为优化计划的常驻空间，当前已登记的方向如下；新增条目按同一格式追加。

1. **参数 k 全矩阵复验（P0，实现已完成、验收待执行）**：NdNiO2 `Y a=1/3`
   nmod=0 与 nmod≥1 的全 48 候选逐候选确认“写出四格式”或“显式未解析原因”，并核对
   根级 JSON/TXT 报告与网页状态一致；原有 2/2 代表例不能替代该矩阵。
2. **身份解析覆盖面（P1）**：把精确微观身份链扩展到更多晶系、点阵类型、超胞大小与
   混合/部分占位，并为每个新增覆盖层留下正反例。验收：每层都有定向回归与报告。
3. **批量导出性能（P1，进度/中断/锁隔离已完成）**：研究可证明安全的候选间
   ISO/精确子胞缓存复用，减少仍然存在的串行重算。验收：同规模批量的墙钟时间与峰值
   内存可量化下降，结果与无缓存路径逐候选相同；不得为速度弱化身份或 writer 门禁。
4. **参数 k 子群数据库治理（P2）**：首次生成可能长达数小时的数据库应有集中位置、
   来源与版本记录、进度提示和复用证据；当前它会写进只读 `resources/isobyu/` 数据目录
   （ISO 经 `~/.id/data` 符号链接写入），审计必须用 `_stage_dir` 重定向。验收：生成
   过程可复核、可复用、失败可诊断。
5. **性能与资源账本（P2）**：为模式计算、批量导出、四格式 writer 记录 P50/P95、
   峰值内存、WSL 子进程数与 ZIP 大小；与 §后续扩展验收第 4 项合并执行。
6. **并发与回收测试（P2）**：多标签、多进程、快速重复提交、ZIP/GenDB 中断和缓存
   并发；与 §后续扩展验收第 5 项合并执行。
7. **未归档的优化想法**：本分区预留，后续审计结论先追加到 §待修的程序漏洞，确认
   属于能力扩展或性能改进时再登记到本节。

## 当前顺序与状态

1. **(3+d) superspace 位移模式：已完成。**
2. **ISODISTORT → GD 笔记本输入：已完成。**
3. **Method 1/2/3 官网对照：前两种母相已完成既定集合，4310 Method 1 已完成。**
   三种母相的 Method 1 官网静态下载门禁已核对 330/330 个候选和 1320/1320 个
   核心文件。EuAl4/NdNiO2 的 Method 1/2 曾在签名 `e3347be0…5908` 下完成
   275/275 个候选的 live 对照；该签名早于本轮母相 setting/物理轨道修复，只能作为
   对应源码快照的证据。4310 已修复页头、物理 `orbit_id` 与逐物种允许轨道传递；
   对应源码签名下，整体 `z+0.1` 原点平移输入仍得到 125 个候选，N1+ 4D1 的
   192 个位移模式与官网一致。final7 保留为历史快照；当前 schema-3 报告已在签名
   `96ff44cd307432ae2d8ca672afed6ff7c9e7dedbfc39ee07cfd71fefeaf94d12` 下从头完成：
   125 个官网与 live 候选一一配对，模式数、BUSH 覆盖和零振幅 CIF 语义均为
   125/125，0 失败，结束时签名未变化。该结果闭合 4310 Method 1，但不证明每个数值
   模式向量/归一化或 primary IR/OPD 分解唯一性，也不能外推到任意晶体。
   Method 2 当前 schema-6 数值报告只覆盖 EuAl4 LD1 C1 与 NdNiO2 Y1 C1 两个 nmod=0
   参数-k代表例，四导出各 2/2 通过；不等于全 Method 2 矩阵。该代表例结论只覆盖那两个
   候选：同一 NdNiO2 `Y a=1/3` 查询的其余候选当前无法准备导出，见 §待修的程序漏洞
   BUG-001。Method 3 官网黄金集和
   当前 schema-5 签名报告为 40/40 查询、77/77 authoritative
   embedding：13 组逐字段一致、27 组经精确 Seitz/母群共轭证明 affine 等价、
   0 差异、0 错误。母相仍集中在四方体系，不能据此宣称完整科研级适用范围。
4. **Method 4 分解对照与 debug：当前源码冻结双母相矩阵已完成。** 当前源码本地
   24/24 通过；官网审计为 23 个完全通过、EuAl4 G05 auto-origin 1 个证据警告、
   0 个 inconclusive/fail。F01 物种拒绝、F02 均匀应变成功和 F03 robust 距离阈值拒绝
   路径均已闭合。symmetry-adapted 应变模式现已接入 CIF、IsoVIZ、Complete modes 和
   TOPAS 固定实际晶胞；仍没有 4310 或跨晶系矩阵。
5. **第三晶体 `4310_tetra.cif`：Method 1 已完成，Method 2–4 待建立。**
   当前解析器能从 FullProf 多数据块 CIF 的真实结构块读出 34 原子、约化式
   La4Ni3O10 和 `I4/mmm (#139)`；页头、8 个独立物理 Wyckoff 轨道、逐物种允许
   轨道和非标准 setting/origin 的完整母相规范化已经修复。Method 1 的 125 个官网
   候选已完整归档；N1+ 4D1 的 192 个位移模式及整体 `z+0.1` 原点平移代表例在其
   对应源码签名下通过，当前 schema-3 报告又完成 125/125 全量 live。
   Method 2–4 尚未建立 4310 的正式矩阵，前两种母相结果不得自动外推。
6. **CIF 位移模式与 microscopic provenance：生产 core 与四 writer 已接通。**
   core 现从精确 child Hall frame 构造稳定 atom ID、逐原子父子映射和物理轨道，并把
   未混合 ISO microscopic 来源列、查询的 `(parent/child SG,B,q,primary IR/OPD)` 与
   exact `VECTOR` 方向一起交给 `DisplaciveExportData`。稀疏官网点域只在公共域满秩、
   唯一变基及完整 BUSH 域延拓全部可复核时使用；proof 绑定有序 reference/candidate
   序列，合法 BUSH 基换序保持结果不变，重排或换绑后的 proof 会 fail-closed；直接来源与
   延拓证据分开保存。下一步是把同一合同扩展到更多晶系及混合/部分占位。

不要再把“参数 k 无模式”“部分特殊 k 模式不全”或“GD 输入尚未转换”列为未完成项。

## 新会话续接执行顺序

1. 以 `method1_4310_current_audit_20261006.json` 为 4310 Method 1 当前基线；只有源码、
   输入或运行时联合签名变化时才重跑全量 125 项，不能用 `final7`、`final6` 或旧 75/125 报告
   覆盖它。新会话先读 `output/validation/CRIS_DEBUG_HANDOFF_20261004.md`、本计划与
   验证报告，再检查 `git status`；当前工作树含大量未提交的在途修改，不得 reset。
2. 先处理 §待修的程序漏洞 的 P0 条目：当前源码在 NdNiO2 `Y a=1/3` 上无法完成
   Method 2 四格式批量导出，且该缺陷会掩盖其他覆盖结论。
3. 再建立 4310 Method 2 正式矩阵：覆盖 Gamma、非 Gamma 特殊 k 和一条公度参数 k，
   同时核对 `nmod=0` 与 `nmod>=1`、site-mode 身份、模式维数/归一化/子空间和四种
   writer。每一层先做最小反例与不变量检查，再扩成完整官网差分。
4. Method 2 闭合后生成 4310 Method 3 精确机器清单和查询上下文，覆盖 identity、
   oriented-cell、special-supercell 与 parameter-k。清单冻结前不让用户下载；冻结后
   按 `DOWNLOAD_CHECKLIST.md` 人工保存官网结果并先做完整性审计，再用于算法结论。
5. 只从已核定的 Method 1/2/3 路径生成 4310 Method 4 daughter 矩阵，冻结 CIF、哈希、
   basis/origin/matching 与候选身份；随后再让用户逐项上传和保存官网证据。先核对输入
   签名，再比较幅值、`As/Ap/normfactor`、residual、原子匹配、应变与四类导出。
6. 三个四方母相的 Method 1–4 路径闭合后，继续七晶系、P/A/B/C/I/F/R 点阵、一般/
   特殊 Wyckoff、混合占位、性能和并发扩展验收；在这些证据完成前，不宣称对任意
   晶体均已正确。

## 目标 1：完整的 (3+d) superspace 模式计算（已完成）

完成标准：

1. `backend` 通用计算参数 k 的公度锁定与 (3+d) 模式，不为 LD1 写特例。
2. 参数 k 子群导出的 CIF、IsoVIZ、Complete modes details 和 TOPAS 含真实位移模式。
3. nmod=0 保留折叠进子群晶胞的全部母相 k；Method 2 的 nmod=1/2/3 只保留当前 q 的谐波和 Γ。
4. 输出可稳定提供 GD/机器学习所需的模式基矢、名称、归一化与振幅边界。

状态与证据见 `BUGFIX_VALIDATION_REPORT.md`。

## 目标 2：生成 GD 笔记本输入（已完成）

工作位置是与 CRIS 同级的 `GD/`。入口为该目录的 `isodistort_to_gd.py`，使用 CRIS `.venv`，默认从同级 `CRIS/` 定位项目根（非同级布局用 `CRIS_ROOT` 覆盖），输出到 `generated/{irrep}_{opd}/`。禁止修改 tianren 参考笔记本和 `D:\OneDrive\...` 原始精修数据。

转换结果必须包含笔记本实际读取的模式名称、归一化因子、振幅上界、displacive mode 表和可导入的 `*_alris_functions.py`。实验衍射表 `All Combined.csv` 不属于 ISODISTORT 产物。

## 目标 3：Method 1/2/3 官网对照（前两种母相已完成）

对照输入 `webpage_info/`、`output_compare/`、`experiment_data/` 和 `resources/isobyu/` 均只读。新的官网下载放入用户指定的新目录或现有对照目录的明确新批次中，不覆盖旧证据。

### 3.1 Method 1（EuAl4/NdNiO2 保存产物已复验）

1. 每种母相按官网存档选择 strains + displacive 及相应物种。
2. 对照整张结果表的 `Irrep`、`OPD`、`Dir`、`SG`、`basis`、`origin`、`s`、`i`、`k-active`；本地额外的 `idx` 不参与。
3. 核对全部候选与核心文件；抽查非 GM 特殊 k 的 ZIP、CIF 语义、模式标签/数量及 VESTA/IsoVIZ 可打开性。
4. basis/origin 仅为等价设置时，用精确有理数与整数幺模变换证明同一子格，再归一到官网设置；报告区分直接匹配与等价匹配。

### 3.2 Method 2（当前保存产物已复验）

1. EuAl4 覆盖 LD `(0,0,g)`、`g=1/6`；NdNiO2 覆盖 Y `(a,1/2,0)`、`a=1/3`；每个母相另含一个特殊 k。
2. 参数 k 分别验证 nmod=0 和 nmod=1；nmod=2/3 与 nmod=1 等价，不重复冒充独立 q 测试。
3. 至少对 EuAl4 LD1 C1 及 NdNiO2 一个参数 k 子群逐项比较模式标签、数量、归一化和四类导出。
4. 源码签名 `e3347be0…5908` 下的四批 live 结果分别为 EuAl4 Method 1 `123/123`、Method 2 `22/22`、NdNiO2 Method 1 `82/82`、Method 2 `48/48`；missing、extra、duplicate 和 mode-count mismatch 全为 0。该结果只证明所记录的源码快照；本轮母相 setting/物理轨道代码改变后，不再称为“当前最终源码”结果。精确报告索引见 `BUGFIX_VALIDATION_REPORT.md`。
5. 当前源码的 schema-6 数值语义 live 报告只复验 EuAl4 `LD1 C1, g=1/6, nmod=0`
   与 NdNiO2 `Y1 C1, a=1/3, nmod=0`：两例及其 CIF/IsoVIZ/Complete modes/TOPAS
   四种导出均通过。报告为
   `output/validation/method2_numeric_semantic_live_current.json`，SHA-256 为
   `c4830c1895681f19031c6c41113ea5d4039e14e1af2dbad06fe05661edc65cfb`。它是两个代表例的
   当前源码证据，不把历史 22/48 全候选矩阵、nmod=1 或其他晶体自动升级为当前结论。

### 3.3 Method 3（两种母相均 20/20 官网集合已核定，coupled 产品链路已接入）

每个母相下载并验证 20 组不同查询，即 EuAl4 20 组 + NdNiO2 20 组，共 40 组官网运行。每组是一个确定的空间群、direct 子格基矢和官网 **Default** centering，以及由此产生的整张结果表。不得把旧本地 40/40 路径命中预检写成“Method 3 官网差分完成”；正式判据是官网结果表的完整 `(SG, basis, origin, s, i)` embedding 集合。

40 个精确输入冻结在 `docs/manifests/method3_download_manifest.json`。2026-09-24 的历史 WSL/iso 预检曾为 40/40 路径命中、83 个本地候选，但其算法把 Method 1 单-IR/OPD 路径当作 Method 3 embedding，并把 P 错当 Default，故已标记为历史失效结果；`candidate_count` 只是旧本地原始计数。当前已用精确有理数修复 direct 子格解释：Default 采用目标空间群默认 centering，P/A/B/C/I/F/R 各自转换到 primitive translation lattice；恒等和非恒等 basis 均在母相点群轨道下做 exact GL(3,Z) 同格匹配，并拒绝不是母晶格子格的输入。

已实现的一参数公度路径不再猜测 `1/d`：它保留字符串 Fraction 的精确值，解析一变量仿射 k 坐标，在母相完整 reciprocal k-star 上求 `0 ≤ p < 1` 的所有 `T k(p) ∈ Z³` 解，并按母相中心化 direct primitive translations 定义的 reciprocal-lattice 等价关系排除特殊-k 端点重复。对同一物理 k-star 的不同参数代表元会在调用后端前去重；NdNiO2 `Y(a=1/3)` 与 `Y(a=2/3)` 不再生成两套不同原点的重复行。它仍只是单-IR、一参数线子集，不得写成完整 Method 3。

下载目录已按 manifest 预建。官网结果直接放入 `output_compare/<母相 CIF>/官网/Method3/<短案例号>/`；重新生成的本地网页结果放入 `output_compare/<母相 CIF>/现有网页版交互/Method3/<短案例号>/`。两种母相各有 20 个案例目录；每个目录只保存该次查询的完整结果页与导出 ZIP，避免不同查询的同名候选相互覆盖。逐组输入和目录名见 `docs/DOWNLOAD_CHECKLIST.md`。

2026-09-27 的通用只读审计器 `tests/manual/audit_method3_downloads.py` 确认：EuAl4 20/20 个案例、35 个官方候选和 140/140 个核心文件均通过，0 错误/警告；NdNiO2 20/20、42 个候选和 168/168 个核心文件通过，0 错误、12 条 `noncanonical_*_page_location` 布局提示。审计器从页面内容核对查询上下文和完整候选表，并且只按 details 中完全一致的 `(SG,basis,origin,s,i)` 配对，所以这些提示不表示漏下载或内容冲突。机器报告分别为 `output/validation/method3_official_download_audit_eual4.json` 与 `method3_official_download_audit_ndnio2.json`。

新增条款下的算法审计确认，官网第一阶段每一行是唯一 `(space-group type, supercell basis, supercell origin)`，并自动考虑合适取向与原点。历史 single-IR 产品曾得到 13 组逐字段完全一致、16 组精确仿射等价、11 组欠枚举；缺少的 16 个全部为 coupled-IR-only。当前实现从 Stage-A affine embedding 出发，在 `strain ⊕ displacive` 表示中用 Stage-B exact fixed-space 判据拒绝不可达项，再把枚举到的单-IR 稳定子提升到同一 `N_G(T_s)/T_s` 有限商，要求其精确交集等于目标嵌入；最后按母群仿射共轭保留一个 Method-3 首屏代表元。`tests/manual/compare_method3_official_local.py --restart` 当前 schema-5 报告为 40/40 authoritative 查询、13 组逐字段相同、27 组精确 affine 等价、0 差异、0 本地错误/跳过，逐组候选数均一致。机器报告为 `output/validation/method3_official_local_comparison_20261006_current.json`，SHA-256 为 `796ae2f39235c6d61422cb9d02403da11e795a60815f6a0590ea02b45c6c45f7`。脚本让签名与解析消费同一份冻结 manifest/官网案例字节快照，并在结束时重新签名；漂移会使报告和 checkpoint fail-closed，任何后续源码或证据变化都必须从头重跑。

通用只读 route 审计器 `tests/manual/audit_method3_embedding_routes.py` 对 40 组、77 条 authoritative embedding 执行 ISO `DISPLAY DIRECTION`：61 条 `single_ir_exact`、16 条 `coupled_ir_required`，0 条不确定、0 个错误/跳过。机器报告为 `output/validation/method3_embedding_route_audit.json`。项目所带 ISO 9.6.1 没有直接枚举 Method 3 首屏 affine embedding 的命令；`DISPLAY DIRECTION` 和 `DISPLAY ISOTROPY COUPLED` 都属于已知 embedding 之后的第二阶段探针，不能拿来代替第一阶段。`SHOW DOMAIN` 会把 SG12 的 3 个官网 embedding 错扩成 8 个畴，故已否决。

当前空间群查询的 coupled 产品链路如下；尚未覆盖的搜索域继续按论文与空间群仿射代数推进：

1. **阶段 A 精确核心（空间群产品查询已启用）**：按所选子格的 normalizer 有限商 `N_G(T_s)/T_s` 枚举闭合点子群与 affine lifts，施加 Seitz/factor-set 闭包，识别目标 SG 类型；spglib 只提议标准 setting，最终由精确算子重建复核。全量诊断为 40/40 案例、142 个候选，覆盖 77/77 官网 embedding；65 个表面多余项中 62 个只是母群共轭多重性，3 个为 `ND-12`/`ND-13` 的真正新轨道。
2. **阶段 B exact fixed-space 可达性（已接产品）**：在 `strain ⊕ displacive` 有限商表示上用 exact character/generated-subgroup 判据验证 `Stab_G(Fix(H))=H`。40 组、77 个官网 embedding 全部可达；142 个阶段 A 候选中 139 个可达，3 个新轨道全部判 `embedding_infeasible`，恰好消除 Stage-A 假阳性。报告为 `output/validation/method3_stage_b_feasibility_audit.json`。
3. **coupled 见证与完整模式（已接产品）**：对每个可达 route-less embedding，把可能的单-IR 稳定子连同其平移陪集提升到目标有限商，用动态规划寻找最小精确交集；只有交集严格等于目标 `H` 才标为 `exact_fixed_space` 并允许进入 Method 2。模式不把任一见证组合冒充唯一主 IR，而是以 nmod=0 计算目标子群全部折叠 k 的完整 displacive fixed space。代表实测：EuAl4 M3-EU-04 coupled C2/m 为本地/官网 `2/2` 个模式；NdNiO2 M3-ND-18 参数-k coupled Pmm2 为 `22/22` 个模式。
4. **仍待扩展的 Inverse Landau/COPL 搜索域**：空间群查询当前只从特殊 k 与一参数公度线枚举单-IR 稳定子；任意/多参数 k、point-group-only affine 枚举和 reciprocal-sublattice 尚未实现。显式“唯一 primary coupled IR/OPD”分解也未声明；产品给出的是经精确见证的 embedding 与完整固定子空间模式。
5. **并发安全（当前最低闭环已完成）**：状态 API、搜索、选择和导出已共用 `RLock`；单调 `revision` 会拒绝旧标签的 stale mutation/export，上传文件使用 UUID 名。真正的 per-client/session 独立 `IsoDistort` 实例仅作可选后续扩展。
6. **大超胞成本预算（已完成）**：`runtime.method3_max_parametric_values` 和 `runtime.method3_max_backend_queries` 已进入 `resources/config/settings.yaml`，`0` 表示不限制。所有预算超限均显式失败，禁止截断；阶段 A/B 也具有独立边界。Stage B 使用 exact character/generated-subgroup 判据避免构造不必要的稠密 `3N×3N` 有理矩阵，但物理判据不变。

已修复的 ZIP 状态问题不再保留为计划项：批量导出会把同一请求 `nmod` 显式传给每个候选，只有缓存 `nmod` 相同才复用现有模式/结构；否则统一重算并在结束后恢复原会话，避免不同 superspace 模式基混入一包。

覆盖要求：

1. 不同空间群号、恒等或轴向重排的单位体积代表基矢，以及至少两种非恒等整数超胞；官网批次统一使用 Default centering，P/no-centering 不再冒充已覆盖项。
2. 至少 4 组能推断公度参数 k；每组对照候选数量和除 `idx` 外全部字段。
3. 每组至少导出一个子群 ZIP；参数 k 子群按 Method 2 的完整 superspace 模式标准检查。
4. reciprocal 子格本地不支持，不计入 20 组，保留为明确限制。

原先的两组冒烟已纳入完整清单（分别为 `M3-EU-17`、`M3-ND-17`）：

- **EuAl4 Parent.cif / `M3-EU-17`**：Types=strain+displacive，SG `99 P4mm`，direct/Default，basis `diag(1,1,6)`；真实结果表和候选导出已通过审计，官网 1 个 embedding 与本地逐字段一致。
- **NdNiO2 own.cif / `M3-ND-17`**：Types=strain+displacive，SG `47 Pmmm`，direct/Default，basis `{(-3,0,0),(0,0,1),(0,2,0)}`；官网 2 个 embedding 均已由结果表和精确候选身份核定，本地默认 single-IR 路径逐项覆盖。

## 目标 4：Method 4 模式分解验证与 debug（冻结快照与官网 24/24，官网 1 个证据 warning）

### 4.1 前置条件

1. 目标 3 的 Method 1/2 保存产物已经复验，Method 3 官网黄金集也已完整归档并通过只读审计；本前置条件现已满足。
2. 官网输入必须同时保存母相 CIF、女儿相 CIF、选中的 path/子群模式上下文、Method 4 结果页和下载的 txt/csv，避免只保存幅度表而无法复现归一化。

### 4.1.1 官网下载要求与时机

Method 4 需要官网下载，但**不要下载随机案例**。本地准备脚本 `tests/manual/prepare_method4_inputs.py` 按重冻结源码签名 `23c951b9…b279` 管理 §4.2 的 24 个女儿相 CIF；机器清单为 `docs/manifests/method4_download_manifest.json`，输入根目录为 `output/validation/method4_inputs/`。清单记录每个 CIF 的 SHA-256、完整 `(k,IR,OPD,SG,basis,origin,s,i)` 上下文与全部本地模式标签；验证脚本直接从该上下文重建参数-k候选，不再为一个已冻结案例重新枚举整张 Method 3 表。P02/P03 的 `dmax=0.30/0.35 Å` 输入保持有效，F02 已重分类为均匀应变成功例。完整映射与门禁顺序见 `docs/DOWNLOAD_CHECKLIST.md`。

EuAl4 12 例已全部下载并通过当前数值审计。G01–G06、P01–P03 核心数值与上传哈希通过；G05 auto-origin 仅缺 basis HTML并维持证据警告。F01/F03 分别按物种不匹配和 robust 距离超阈值被官网拒绝。F02 的 30 个位移模式全零，6 个应变模式中 4 个非零；本地从母相 metric、冻结 basis 与女儿相 metric 解出唯一对称正定 `M=I+epsilon`，重建相对残差约 `1.0e-15`，与官网 IsoVIZ 应用张量最大差 `5.08e-6`。NdNiO2 本地与官网 12 例也全部通过：G01–P03 已通过上传哈希、输入控件、P1 结果身份、Γ 与参数-k 幅度、完整分支、显式/自动原点、P1 微噪声分配与导出结构差分；F01 按物种不匹配拒绝，F02 均匀应变差分通过，F03 的正确冻结哈希、robust `dmax=0.1 Å` 截图和匹配失败语义通过。

正式批次固定为 24 组官网运行（EuAl4 12 组、NdNiO2 12 组）。每组保存：母相 CIF 标识、原始女儿相 CIF、所用 path/子群及 Complete modes details、Method 4 完整结果页、幅度表 txt/csv；预期失败案例保存官网错误页/错误文案。若官网没有某种下载格式，保留完整 HTML 和可复制的结果表，不自行补造文件。归档审计按 HTML 内容结构和角色识别页面：原候选目录内可任意改 basename，但须保留 `.html`/`.htm` 后缀；要求角色为 0 页或多页都会 fail-closed，跨候选目录移动不视为改名。另有“未准备模式即调用”这一项只做本地流程错误测试，不要求上传 CIF。

### 4.2 验证数据矩阵

对 EuAl4 与 NdNiO2 分别建立至少 12 个可复现案例，覆盖：

1. 同晶胞 Γ/特殊 k path 与参数 k 公度超胞 path 各至少一组。
2. 零畸变、单模式正/负振幅、两个正交模式、多模式混合和接近振幅上界。
3. 原子行重排、周期边界换像、允许的整体原点平移；验证 nearest-site 匹配不依赖 CIF 行序。
4. 母相同尺寸与超胞女儿相；至少一组含轻微数值噪声。本地受限模式空间应给出非零 residual；官网若因 P1 完整模式空间吸收噪声，则记录额外噪声模式，不能强求两者 residual 相同或把完整 P1 拟合误判为伪模式。
5. 明确失败案例：元素/化学计量不匹配、原子距离超过阈值、所选 basis 与无畸变参考胞不一致、尚未准备模式就调用 Method 4。对固定 basis 的任意非奇异正定晶格 metric，六维均匀应变原则上都能解释，不能再把普通晶格形变伪造为“incompatible lattice”失败例。

优先用本地 `generate_distortion` 从已知幅度生成女儿相，再用 Method 4 回归，形成“生成 → 分解”的闭环；官方差分使用完全相同的母相、女儿相与 path。

### 4.3 对照项目

1. 模式全集、模式标签、排序及幅度的符号/尺度；若官网基矢仅差整体符号，先证明基矢等价，再比较相应反号幅度。
2. `RMS residual`、最大绝对 residual、原子匹配和原点平移；零噪声可表示输入应接近机器精度。只有模式空间相同时才直接比较 residual；官网 P1 完整空间与本地受限子群空间不同时，改为比较共享模式并记录额外 P1 噪声分量。
3. 网页和 Python API 必须调用同一 `IsoDistort.search_method_4` 路径，并给出等价幅度和 residual。
4. Distortion 的 Method 4 txt/csv 必须包含全部筛选命中行；Method 4 不产生子群 ZIP，错误提示需与网页/API 一致。
5. 超胞模式提升必须与生成方向互逆；检查晶格坐标/笛卡尔坐标、单位和归一化，不以调幅度常数掩盖错误。

### 4.4 Debug 原则与完成标准

先以可控的合成女儿相定位算法错误，再做官网差分。发现差异时依次检查原子对应、周期最短位移、原点、超胞变换、模式矩阵秩/条件数和幅度规范；不得为某个模式名称写补丁。

完成标准：所有无噪声闭环在约定容差内恢复输入幅度，带噪声案例 residual 合理，错误输入稳定拒绝，网页/API 和 txt/csv 一致，并完成两种母相的官网差分。修复记录与结果统一追加到 `BUGFIX_VALIDATION_REPORT.md`。

当前状态：当前源码本地冻结矩阵 24/24 通过，schema-3 报告 SHA-256 为
`e19a005dd77bad4b12c7dcc2600889c8e89fe3a689d572c542e6f626a954bb5d`。官网合计
23 个完全通过、EuAl4 G05 auto-origin 1 个证据警告、0 个 inconclusive/fail，schema-5
报告 SHA-256 为 `946dba00162ca813060153c4f38daf3131030ffe352525bcc439a2062c5f9519`，
综合状态为 `pass_with_warnings`。NdNiO2 F03 的正确 daughter SHA-256 为
`02b41e63…332e36`，填写截图确认 frozen basis、automatic origin、robust `dmax=0.1 Å`
和未选 manual mapping，官网按预期返回匹配失败。两母相归档/当前源码冻结矩阵已闭合；
4310、跨晶系及下列未实现能力仍是后续门槛。

## 目标 5：第三晶体 `4310_tetra.cif` 四种 Method 对照（Method 1 已完成）

### 5.1 是否可用及门禁

`experiment_data/4310_tetra.cif` 是 FullProf 输出的多数据块 CIF：首个 `data_global` 不含原子结构，真正结构在后续数据块。当前 `pymatgen` 解析会对空块给出警告，但能继续读取真实结构块；实测得到 34 原子、La4Ni3O10、`a=b=3.8516 Å`、`c=27.967 Å` 和 `I4/mmm (#139)`。文件没有显式 superspace CIF 标签，本身是三维母相；已完成的 `(3+d)` 模块解决的是后续参数-k 子群的模式计算与导出，不是这个 CIF 的基本载入问题。官网页头对应 8 个独立物理轨道，其中 5 个都属于 `4e`；本地现在以稳定物理 `orbit_id` 而非单独 Wyckoff 字母区分它们，并按所选物种把允许轨道传给 BUSH/smodes。

三种晶体的 Method 1 官网和现有网页版保存侧门禁均已满足：两侧各 `330/330` 个候选、
各 `1320/1320` 个核心文件完整且身份正确，所以无需再次下载，审查和 debug 已获有效
证据支撑。4310 子集在两侧各为 125 个候选、500 个核心文件，N1+ 4D1 当前目录内容
正确。此前疑似发生覆盖的具体旧目标目录已无足够证据还原，只能记为 `inconclusive`，
不能猜测到某个子群名下。官网与旧本地产物的 175 项科学差异全部来自 4310，其中
125 条为模式数差异、50 条为 CIF 语义差异；另有 108 条 BUSH 覆盖诊断。这些是修复前
基线，不是下载错误，也不能继续冒充当前源码结论。当前 schema-3 报告已以带源码、
输入与运行时签名的 125 候选 live 重算闭合 Method 1；后续按下列顺序推进 Method 2–4。

### 5.2 后续正式工作清单

1. **Method 1（已完成）**：下载完整性门禁已通过；当前修复了页头标准代表、稳定物理
   `orbit_id`、逐物种允许轨道传递，以及 FINDSYM 要求换 basis/origin 时对完整母相
   做一致规范化的根因。对应源码签名下，原始 4310 与整体 `z+0.1` 输入都返回
   125 个候选，N1+ 4D1 均得到官网的 192 个位移模式。当前 schema-3 报告已完成 125/125 live
   身份、模式数、BUSH 覆盖和 CIF 语义对照，0 失败；C8、P5|4D1、P2|C1、X1-|C1
   以及 C10/C12 的模式向量、标签和来源另有定向 live 证据。`final6` 在运行期间发生
   源码漂移，只保留为 115/125 的根因诊断；final7 只保留为较早稳定快照。当前全量报告
   也不独立证明 125 项中每个数值模式向量/归一化或 primary IR/OPD 分解唯一性，
   不改变更多晶系仍待扩展的限制。
2. **Method 2**：覆盖 Γ、非 Γ 特殊 k 与至少一条公度参数 k；核对 `nmod=0` 和 `nmod>=1` 的保留规则、site-mode 标签、维数、归一化、模式子空间及四格式导出。参数-k 是本晶体必须单列的 superspace 回归，不得只复用 EuAl4 的结论。
3. **Method 3**：建立 identity、oriented-cell、special-supercell、parameter-k 的分层官方清单；逐项证明 `(SG,basis,origin,s,i)` 身份、精确等价子格和 single/coupled-IR route，不按同空间群移植现有答案。
4. **Method 4**：从已核定的 Method 1/2/3 path 生成冻结女儿相，覆盖零、正/负单模式、正交双模式、参数-k、多模式/近上界、换像/重排/原点/噪声和明确失败输入；验证幅值、`As/Ap/normfactor`、residual、原子匹配与导出，再做完全相同输入的官网差分。
5. **接口与科研判据**：网页、Python API 以及 txt/csv/ZIP 必须共享身份与数值；记录源码/依赖/输入签名和资源成本。4310 与前两晶体全通过只能证明这三个四方母相的覆盖范围；科研级交付结论还必须结合后续扩展验收的七晶系、点阵、Wyckoff、组成、性能和并发证据，不能由第三个同晶系案例单独推出。

## 保留的产品限制

- symmetry-adapted 应变模式已实现：CIF 写 max-component-one 的 `q_raw`，IsoVIZ 写单位 Frobenius 范数的 `q_unit`，Complete modes 同时写两者；TOPAS 按官网行为只写 `B M P` 的固定实际晶胞，不生成官网不存在的 strain 精修参数。该链路依赖 ISO rank-[12] 宏观基、目标 embedding invariant direction 和独立 metric fixed-space 整空间核验；证据不全即 fail-closed。当前源码已完成 EuAl4/NdNiO2 冻结矩阵 24/24 本地重跑；4310 与跨晶系泛化仍待扩展验收。
- CIF 位移模型与四 writer 的共享合同已接入生产 core：精确 emitted frame、稳定 atom ID、
  父子映射、完整子群查询上下文和未混合 microscopic source-column provenance 会在 mapper、
  公共 I/O 合同及 writer 边界重复核验。合同暂只接受单物种满占位点；混合/部分占位与更多
  晶系仍待扩展验收。
- Method 2 的 nmod=2/3 不增加第二、第三条独立 q。
- Method 3 reciprocal 子格不支持。
- Method 3 空间群查询默认启用精确 affine 阶段 A、fixed-space 阶段 B 和 coupled 稳定子交集见证；在当前双母相 40 组权威查询中已覆盖 77/77 embedding，并连接完整位移 fixed-space 模式。仍未覆盖 arbitrary/multi-parameter k、point-group-only affine 枚举、reciprocal 子格，以及唯一 primary coupled IR/OPD 分解；这些范围不能由当前通过结果外推。
- rotational-only Types 目前借用 smodes 位移活性近似过滤，未实现独立的刚性转动/轴矢量模式生成。
- polar-vector 位移标签使用 exact `DISPLAY DIRECTION` → primitive-integer
  `VALUE DIRECTION VECTOR,...` → oriented microscopic 列的证据链。官网稀疏点域可在唯一、
  满秩且数值稳定的 BUSH 域延拓后使用；纯数值 SMODES fallback 和证据不完整的列保持
  `unresolved`。非 Γ 非零 transporter 在缺 component-to-star-arm 相位身份时仍明确拒绝；
  轴矢量、磁、占位模式的同等级身份链尚未实现。

## 后续扩展验收（尚未实施）

以下内容由原 `tests/manual/METHOD_VALIDATION_PLAN.md` 合并而来；本文件仍是唯一有效开发计划：

1. 冻结可公开再分发的分层 CIF 样本，覆盖七晶系、P/A/B/C/I/F/R 点阵、一般/特殊 Wyckoff、单/多元素、混合占位和不同 ASU 大小；来源统一维护在 `docs/EXTERNAL_CIF_SOURCES.md`。
2. 建立受限速保护的官网周期性重抓与版本漂移报告；官网不可用、证据缺失或配对歧义只记为 inconclusive，不计入通过率。
3. 扩展 TOPAS 解析/本体测试及 VESTA、IsoVIZ 分层抽检；未安装商业或图形工具时不得把静态检查表述为本体兼容证明。
4. 进行三次冷/热重复，记录 P50/P95、峰值内存、WSL 进程、ZIP 大小和规范化哈希。
5. 增加多标签、多进程、快速重复提交、ZIP/GenDB 中断和缓存并发测试，并验证服务、端口及 WSL 子进程均能回收。
6. 在看到扩展结果之前冻结发布阈值；固定母相全通过只证明已覆盖范围，不构成对未知结构零缺陷的保证。
