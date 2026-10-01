# ISODISTORT 本地化开发计划

本文件是 `ISODISTORT/` 唯一的开发目标与后续验证计划。开发规则、修改边界和说明性文件职责见 [agent.md](../agent.md)；已经完成的修复、实测结果和机器报告索引见 [BUGFIX_VALIDATION_REPORT.md](BUGFIX_VALIDATION_REPORT.md)。

实现任何目标时都禁止按单个母相、IR 或 OPD 硬编码官网答案。应修复可推广的晶体学算法，并把无法实现的能力明确列为限制。

## 当前顺序与状态

1. **(3+d) superspace 位移模式：已完成。**
2. **ISODISTORT → GD 笔记本输入：已完成。**
3. **Method 1/2/3 官网对照：当前双母相范围已完成；跨晶系扩展仍继续。** Method 1/2 的保存产物与最终源码 live 计算均覆盖 275/275 个官方候选，候选集合、模式数和完整性全部一致；代表性的标签/向量/归一化/四格式数值语义也已通过。Method 3 官网黄金集已完整核定：EuAl4 20/20、35 个 embedding、140/140 个核心文件、0 错误/警告；NdNiO2 20/20、42 个 embedding、168/168 个核心文件、0 错误，另有 12 条不影响内容身份的网页包层级提示。当前默认产品重算 40/40 查询、77/77 authoritative embedding：13 组逐字段一致、27 组经精确 Seitz/母群共轭证明 affine 等价、0 差异、0 错误。原先 16 条 coupled-only 已由 Stage A/B、精确稳定子交集见证和完整 fixed-space 模式链路接入。80/80 个等价 basis/centering 变形通过，但母相仍集中在四方体系，不能据此宣称完整科研级适用范围。
4. **Method 4 分解对照与 debug：当前双母相矩阵已完成。** 本地与官网均为 24/24 个有效冻结输入；官网审计为 23 个完全通过、EuAl4 G05 auto-origin 1 个证据警告、0 个失败。F01 物种拒绝、F02 均匀应变成功和 F03 robust 距离阈值拒绝路径均已闭合。该结论只覆盖当前两种四方母相和已列能力，不外推为跨晶系科研级精度。
5. **第三晶体 `4310_tetra.cif`：可用性预检完成，Method 3 coupled 门禁已解除，尚未开始四种 Method 的正式矩阵。** 当前解析器能从 FullProf 多数据块 CIF 的真实结构块读出 34 原子、约化式 La4Ni3O10 和 `I4/mmm (#139)`；可作为结构复杂度对照组。下一步可据独立冻结清单开始该晶体的 Method 1–4 验证，但双母相结果不得自动外推为它已通过。

不要再把“参数 k 无模式”“部分特殊 k 模式不全”或“GD 输入尚未转换”列为未完成项。

## 目标 1：完整的 (3+d) superspace 模式计算（已完成）

完成标准：

1. `isocore` 通用计算参数 k 的公度锁定与 (3+d) 模式，不为 LD1 写特例。
2. 参数 k 子群导出的 CIF、IsoVIZ、Complete modes details 和 TOPAS 含真实位移模式。
3. nmod=0 保留折叠进子群晶胞的全部母相 k；Method 2 的 nmod=1/2/3 只保留当前 q 的谐波和 Γ。
4. 输出可稳定提供 GD/机器学习所需的模式基矢、名称、归一化与振幅边界。

状态与证据见 `BUGFIX_VALIDATION_REPORT.md`。

## 目标 2：生成 GD 笔记本输入（已完成）

工作位置是桌面 `GD（未同步git）`，不是仓库内只读的 `CRIS/GD/`。入口为该目录的 `isodistort_to_gd.py`，使用 CRIS `.venv`，输出到 `generated/{irrep}_{opd}/`。禁止修改 tianren 参考笔记本和 `D:\OneDrive\...` 原始精修数据。

转换结果必须包含笔记本实际读取的模式名称、归一化因子、振幅上界、displacive mode 表和可导入的 `*_alris_functions.py`。实验衍射表 `All Combined.csv` 不属于 ISODISTORT 产物。

## 目标 3：Method 1/2/3 官网对照（已完成）

对照输入 `webpage_info/`、`output_compare/`、`experiment_data/` 和 `isobyu/` 均只读。新的官网下载放入用户指定的新目录或现有对照目录的明确新批次中，不覆盖旧证据。

### 3.1 Method 1（当前保存产物已复验）

1. 每种母相按官网存档选择 strains + displacive 及相应物种。
2. 对照整张结果表的 `Irrep`、`OPD`、`Dir`、`SG`、`basis`、`origin`、`s`、`i`、`k-active`；本地额外的 `idx` 不参与。
3. 核对全部候选与核心文件；抽查非 GM 特殊 k 的 ZIP、CIF 语义、模式标签/数量及 VESTA/IsoVIZ 可打开性。
4. basis/origin 仅为等价设置时，用精确有理数与整数幺模变换证明同一子格，再归一到官网设置；报告区分直接匹配与等价匹配。

### 3.2 Method 2（当前保存产物已复验）

1. EuAl4 覆盖 LD `(0,0,g)`、`g=1/6`；NdNiO2 覆盖 Y `(a,1/2,0)`、`a=1/3`；每个母相另含一个特殊 k。
2. 参数 k 分别验证 nmod=0 和 nmod=1；nmod=2/3 与 nmod=1 等价，不重复冒充独立 q 测试。
3. 至少对 EuAl4 LD1 C1 及 NdNiO2 一个参数 k 子群逐项比较模式标签、数量、归一化和四类导出。
4. 最终源码签名 `e3347be0…5908` 下的四批 live 结果分别为 EuAl4 Method 1 `123/123`、Method 2 `22/22`、NdNiO2 Method 1 `82/82`、Method 2 `48/48`；missing、extra、duplicate 和 mode-count mismatch 全为 0。精确报告索引见 `BUGFIX_VALIDATION_REPORT.md`。

### 3.3 Method 3（两种母相均 20/20 官网集合已核定，coupled 产品链路已接入）

每个母相下载并验证 20 组不同查询，即 EuAl4 20 组 + NdNiO2 20 组，共 40 组官网运行。每组是一个确定的空间群、direct 子格基矢和官网 **Default** centering，以及由此产生的整张结果表。不得把旧本地 40/40 路径命中预检写成“Method 3 官网差分完成”；正式判据是官网结果表的完整 `(SG, basis, origin, s, i)` embedding 集合。

40 个精确输入冻结在 `docs/manifests/method3_download_manifest.json`。2026-09-24 的历史 WSL/iso 预检曾为 40/40 路径命中、83 个本地候选，但其算法把 Method 1 单-IR/OPD 路径当作 Method 3 embedding，并把 P 错当 Default，故已标记为历史失效结果；`candidate_count` 只是旧本地原始计数。当前已用精确有理数修复 direct 子格解释：Default 采用目标空间群默认 centering，P/A/B/C/I/F/R 各自转换到 primitive translation lattice；恒等和非恒等 basis 均在母相点群轨道下做 exact GL(3,Z) 同格匹配，并拒绝不是母晶格子格的输入。

已实现的一参数公度路径不再猜测 `1/d`：它保留字符串 Fraction 的精确值，解析一变量仿射 k 坐标，在母相完整 reciprocal k-star 上求 `0 ≤ p < 1` 的所有 `T k(p) ∈ Z³` 解，并按母相中心化 direct primitive translations 定义的 reciprocal-lattice 等价关系排除特殊-k 端点重复。对同一物理 k-star 的不同参数代表元会在调用后端前去重；NdNiO2 `Y(a=1/3)` 与 `Y(a=2/3)` 不再生成两套不同原点的重复行。它仍只是单-IR、一参数线子集，不得写成完整 Method 3。

下载目录已按 manifest 预建。官网结果直接放入 `output_compare/<母相 CIF>/官网/Method3/<可读案例目录>/`；重新生成的本地网页结果放入 `output_compare/<母相 CIF>/现有网页版交互/Method3/<可读案例目录>/`。两种母相各有 20 个案例目录；每个目录只保存该次查询的完整结果页与导出 ZIP，避免不同查询的同名候选相互覆盖。逐组输入和目录名见 `docs/DOWNLOAD_CHECKLIST.md`。

2026-09-27 的通用只读审计器 `tests_dev/manual/audit_method3_downloads.py` 确认：EuAl4 20/20 个案例、35 个官方候选和 140/140 个核心文件均通过，0 错误/警告；NdNiO2 20/20、42 个候选和 168/168 个核心文件通过，0 错误、12 条 `noncanonical_*_page_location` 布局提示。审计器从页面内容核对查询上下文和完整候选表，并且只按 details 中完全一致的 `(SG,basis,origin,s,i)` 配对，所以这些提示不表示漏下载或内容冲突。机器报告分别为 `output/validation/method3_official_download_audit_eual4.json` 与 `method3_official_download_audit_ndnio2.json`。

新增条款下的算法审计确认，官网第一阶段每一行是唯一 `(space-group type, supercell basis, supercell origin)`，并自动考虑合适取向与原点。历史 single-IR 产品曾得到 13 组逐字段完全一致、16 组精确仿射等价、11 组欠枚举；缺少的 16 个全部为 coupled-IR-only。当前实现从 Stage-A affine embedding 出发，在 `strain ⊕ displacive` 表示中用 Stage-B exact fixed-space 判据拒绝不可达项，再把枚举到的单-IR 稳定子提升到同一 `N_G(T_s)/T_s` 有限商，要求其精确交集等于目标嵌入；最后按母群仿射共轭保留一个 Method-3 首屏代表元。`tests_dev/manual/compare_method3_official_local.py --restart` 的当前源码结果为 40/40 authoritative 查询、13 组逐字段相同、27 组精确 affine 等价、0 差异、0 本地错误/跳过，逐组候选数均一致。机器报告为 `output/validation/method3_official_local_comparison.json`；任何后续源码变化都必须凭签名从头重跑，不能沿用旧结论。

通用只读 route 审计器 `tests_dev/manual/audit_method3_embedding_routes.py` 对 40 组、77 条 authoritative embedding 执行 ISO `DISPLAY DIRECTION`：61 条 `single_ir_exact`、16 条 `coupled_ir_required`，0 条不确定、0 个错误/跳过。机器报告为 `output/validation/method3_embedding_route_audit.json`。项目所带 ISO 9.6.1 没有直接枚举 Method 3 首屏 affine embedding 的命令；`DISPLAY DIRECTION` 和 `DISPLAY ISOTROPY COUPLED` 都属于已知 embedding 之后的第二阶段探针，不能拿来代替第一阶段。`SHOW DOMAIN` 会把 SG12 的 3 个官网 embedding 错扩成 8 个畴，故已否决。

当前空间群查询的 coupled 产品链路如下；尚未覆盖的搜索域继续按论文与空间群仿射代数推进：

1. **阶段 A 精确核心（空间群产品查询已启用）**：按所选子格的 normalizer 有限商 `N_G(T_s)/T_s` 枚举闭合点子群与 affine lifts，施加 Seitz/factor-set 闭包，识别目标 SG 类型；spglib 只提议标准 setting，最终由精确算子重建复核。全量诊断为 40/40 案例、142 个候选，覆盖 77/77 官网 embedding；65 个表面多余项中 62 个只是母群共轭多重性，3 个为 `ND-12`/`ND-13` 的真正新轨道。
2. **阶段 B exact fixed-space 可达性（已接产品）**：在 `strain ⊕ displacive` 有限商表示上用 exact character/generated-subgroup 判据验证 `Stab_G(Fix(H))=H`。40 组、77 个官网 embedding 全部可达；142 个阶段 A 候选中 139 个可达，3 个新轨道全部判 `embedding_infeasible`，恰好消除 Stage-A 假阳性。报告为 `output/validation/method3_stage_b_feasibility_audit.json`。
3. **coupled 见证与完整模式（已接产品）**：对每个可达 route-less embedding，把可能的单-IR 稳定子连同其平移陪集提升到目标有限商，用动态规划寻找最小精确交集；只有交集严格等于目标 `H` 才标为 `exact_fixed_space` 并允许进入 Method 2。模式不把任一见证组合冒充唯一主 IR，而是以 nmod=0 计算目标子群全部折叠 k 的完整 displacive fixed space。代表实测：EuAl4 M3-EU-04 coupled C2/m 为本地/官网 `2/2` 个模式；NdNiO2 M3-ND-18 参数-k coupled Pmm2 为 `22/22` 个模式。
4. **仍待扩展的 Inverse Landau/COPL 搜索域**：空间群查询当前只从特殊 k 与一参数公度线枚举单-IR 稳定子；任意/多参数 k、point-group-only affine 枚举和 reciprocal-sublattice 尚未实现。显式“唯一 primary coupled IR/OPD”分解也未声明；产品给出的是经精确见证的 embedding 与完整固定子空间模式。
5. **并发安全（当前最低闭环已完成）**：状态 API、搜索、选择和导出已共用 `RLock`；单调 `revision` 会拒绝旧标签的 stale mutation/export，上传文件使用 UUID 名。真正的 per-client/session 独立 `IsoDistort` 实例仅作可选后续扩展。
6. **大超胞成本预算（已完成）**：`runtime.method3_max_parametric_values` 和 `runtime.method3_max_backend_queries` 已进入 `config/settings.yaml`，`0` 表示不限制。所有预算超限均显式失败，禁止截断；阶段 A/B 也具有独立边界。Stage B 使用 exact character/generated-subgroup 判据避免构造不必要的稠密 `3N×3N` 有理矩阵，但物理判据不变。

已修复的 ZIP 状态问题不再保留为计划项：批量导出会把同一请求 `nmod` 显式传给每个候选，只有缓存 `nmod` 相同才复用现有模式/结构；否则统一重算并在结束后恢复原会话，避免不同 superspace 模式基混入一包。

覆盖要求：

1. 不同空间群号、恒等或轴向重排的单位体积代表基矢，以及至少两种非恒等整数超胞；官网批次统一使用 Default centering，P/no-centering 不再冒充已覆盖项。
2. 至少 4 组能推断公度参数 k；每组对照候选数量和除 `idx` 外全部字段。
3. 每组至少导出一个子群 ZIP；参数 k 子群按 Method 2 的完整 superspace 模式标准检查。
4. reciprocal 子格本地不支持，不计入 20 组，保留为明确限制。

原先的两组冒烟已纳入完整清单（分别为 `M3-EU-17`、`M3-ND-17`）：

- **EuAl4 Parent.cif / `M3-EU-17`**：Types=strain+displacive，SG `99 P4mm`，direct/Default，basis `diag(1,1,6)`；真实结果表和候选导出已通过审计，官网 1 个 embedding 与本地逐字段一致。
- **NdNiO2 own.cif / `M3-ND-17`**：Types=strain+displacive，SG `47 Pmmm`，direct/Default，basis `{(-3,0,0),(0,0,1),(0,2,0)}`；官网 2 个 embedding 均已由结果表和精确候选身份核定，本地默认 single-IR 路径逐项覆盖。

## 目标 4：Method 4 模式分解验证与 debug（本地与官网 24/24，官网 1 个证据 warning）

### 4.1 前置条件

1. 目标 3 的 Method 1/2 保存产物已经复验，Method 3 官网黄金集也已完整归档并通过只读审计；本前置条件现已满足。
2. 用户已于 2026-09-27 明确指示开始 Method 4；启动门禁已满足。
3. 官网输入必须同时保存母相 CIF、女儿相 CIF、选中的 path/子群模式上下文、Method 4 结果页和下载的 txt/csv，避免只保存幅度表而无法复现归一化。

### 4.1.1 官网下载要求与时机

Method 4 需要官网下载，但**不要下载随机案例**。本地准备脚本 `tests_dev/manual/prepare_method4_inputs.py` 按重冻结源码签名 `23c951b9…b279` 管理 §4.2 的 24 个女儿相 CIF；机器清单为 `docs/manifests/method4_download_manifest.json`，输入根目录为 `output/validation/method4_inputs/`。清单记录每个 CIF 的 SHA-256、完整 `(k,IR,OPD,SG,basis,origin,s,i)` 上下文与全部本地模式标签；验证脚本直接从该上下文重建参数-k候选，不再为一个已冻结案例重新枚举整张 Method 3 表。P02/P03 的 `dmax=0.30/0.35 Å` 输入保持有效，F02 已重分类为均匀应变成功例。完整映射与门禁顺序见 `docs/DOWNLOAD_CHECKLIST.md`。

EuAl4 12 例已全部下载并通过当前数值审计。G01–G06、P01–P03 核心数值与上传哈希通过；G05 auto-origin 仅缺 basis HTML并维持证据警告。F01/F03 分别按物种不匹配和 robust 距离超阈值被官网拒绝。F02 的 30 个位移模式全零，6 个应变模式中 4 个非零；本地从母相 metric、冻结 basis 与女儿相 metric 解出唯一对称正定 `M=I+epsilon`，重建相对残差约 `1.0e-15`，与官网 IsoVIZ 应用张量最大差 `5.08e-6`。NdNiO2 本地与官网 12 例也全部通过：G01–P03 已通过上传哈希、输入控件、P1 结果身份、Γ 与参数-k 幅度、完整分支、显式/自动原点、P1 微噪声分配与导出结构差分；F01 按物种不匹配拒绝，F02 均匀应变差分通过，F03 的正确冻结哈希、robust `dmax=0.1 Å` 截图和匹配失败语义通过。

正式批次固定为 24 组官网运行（EuAl4 12 组、NdNiO2 12 组）。每组保存：母相 CIF 标识、原始女儿相 CIF、所用 path/子群及 Complete modes details、Method 4 完整结果页、幅度表 txt/csv；预期失败案例保存官网错误页/错误文案。若官网没有某种下载格式，保留完整 HTML 和可复制的结果表，不自行补造文件。另有“未准备模式即调用”这一项只做本地流程错误测试，不要求上传 CIF。

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
3. 网页、终端和 Python API 必须调用同一 `IsoDistort.search_method_4` 路径，并给出等价幅度和 residual。
4. Distortion 的 Method 4 txt/csv 必须包含全部筛选命中行；Method 4 不产生子群 ZIP，错误提示需与网页/终端一致。
5. 超胞模式提升必须与生成方向互逆；检查晶格坐标/笛卡尔坐标、单位和归一化，不以调幅度常数掩盖错误。

### 4.4 Debug 原则与完成标准

先以可控的合成女儿相定位算法错误，再做官网差分。发现差异时依次检查原子对应、周期最短位移、原点、超胞变换、模式矩阵秩/条件数和幅度规范；不得为某个模式名称写补丁。

完成标准：所有无噪声闭环在约定容差内恢复输入幅度，带噪声案例 residual 合理，错误输入稳定拒绝，网页/终端/API 和 txt/csv 一致，并完成两种母相的官网差分。修复记录与结果统一追加到 `BUGFIX_VALIDATION_REPORT.md`。

当前状态：重冻结清单下本地和官网均为 24/24 个有效输入。官网合计 23 个完全通过、EuAl4 G05 auto-origin 1 个证据警告、0 个失败，综合报告为 `pass_with_warnings`。NdNiO2 F03 的正确 daughter SHA-256 为 `02b41e63…332e36`，填写截图确认 frozen basis、automatic origin、robust `dmax=0.1 Å` 和未选 manual mapping，官网按预期返回匹配失败。两母相 Method 4 当前门禁已关闭；科研级外推仍受晶系范围及下列未实现能力限制。

## 目标 5：第三晶体 `4310_tetra.cif` 四种 Method 对照（等待 Method 3 门禁）

### 5.1 是否可用及门禁

`experiment_data/4310_tetra.cif` 是 FullProf 输出的多数据块 CIF：首个 `data_global` 不含原子结构，真正结构在后续数据块。当前 `pymatgen` 解析会对空块给出警告，但能继续读取真实结构块；实测得到 34 原子、La4Ni3O10、`a=b=3.8516 Å`、`c=27.967 Å` 和 `I4/mmm (#139)`。文件没有显式 superspace CIF 标签，本身是三维母相；已完成的 `(3+d)` 模块解决的是后续参数-k 子群的模式计算与导出，不是这个 CIF 的基本载入问题。因此该文件现在**具备作为新 debug 母相的输入条件**，但尚未完成任一 Method 的端到端或官网验收。

硬门禁：EuAl4 与 NdNiO2 的 Method 1、2、3、4 必须全部结束 debug，所有计划内官网案例、失败路径、证据警告和机器报告均已闭合，才可为 4310 建立清单、运行后端、生成女儿相或下载官网数据。门禁前只允许只读解析和计划维护，不得把预检写成 Method 通过。

### 5.2 门禁解除后的工作清单

1. **Method 1**：冻结母相 SHA-256 和官方查询设置，核对完整候选集合、字段、多重度及 CIF/IsoVIZ/modes/TOPAS 四格式；对比同为 `I4/mmm` 但具有 8 个不等价位点、34 个常规胞原子的复杂度增量。
2. **Method 2**：覆盖 Γ、非 Γ 特殊 k 与至少一条公度参数 k；核对 `nmod=0` 和 `nmod>=1` 的保留规则、site-mode 标签、维数、归一化、模式子空间及四格式导出。参数-k 是本晶体必须单列的 superspace 回归，不得只复用 EuAl4 的结论。
3. **Method 3**：建立 identity、oriented-cell、special-supercell、parameter-k 的分层官方清单；逐项证明 `(SG,basis,origin,s,i)` 身份、精确等价子格和 single/coupled-IR route，不按同空间群移植现有答案。
4. **Method 4**：从已核定的 Method 1/2/3 path 生成冻结女儿相，覆盖零、正/负单模式、正交双模式、参数-k、多模式/近上界、换像/重排/原点/噪声和明确失败输入；验证幅值、`As/Ap/normfactor`、residual、原子匹配与导出，再做完全相同输入的官网差分。
5. **接口与科研判据**：网页、终端、Python API 以及 txt/csv/ZIP 必须共享身份与数值；记录源码/依赖/输入签名和资源成本。4310 与前两晶体全通过只能证明这三个四方母相的覆盖范围；科研级交付结论还必须结合后续扩展验收的七晶系、点阵、Wyckoff、组成、性能和并发证据，不能由第三个同晶系案例单独推出。

## 保留的产品限制

- Method 4 已实现均匀应变张量分解和报告；symmetry-adapted 应变模式的生成、标签与导出仍未实现。
- Method 2 的 nmod=2/3 不增加第二、第三条独立 q。
- Method 3 reciprocal 子格不支持。
- Method 3 空间群查询默认启用精确 affine 阶段 A、fixed-space 阶段 B 和 coupled 稳定子交集见证；在当前双母相 40 组权威查询中已覆盖 77/77 embedding，并连接完整位移 fixed-space 模式。仍未覆盖 arbitrary/multi-parameter k、point-group-only affine 枚举、reciprocal 子格，以及唯一 primary coupled IR/OPD 分解；这些范围不能由当前通过结果外推。
- rotational-only Types 目前借用 smodes 位移活性近似过滤，未实现独立的刚性转动/轴矢量模式生成。
- polar-vector 位移标签已实现 Wyckoff representative/transporter 分类，但所有 site group 与轴矢量、磁、占位模式的通用 character-table decomposition 尚未完成；超出已验证范围的标签不能冒充官网完整 site-symmetry 分解。

## 后续扩展验收（尚未实施）

以下内容由原 `tests_dev/manual/METHOD_VALIDATION_PLAN.md` 合并而来；本文件仍是唯一有效开发计划：

1. 冻结可公开再分发的分层 CIF 样本，覆盖七晶系、P/A/B/C/I/F/R 点阵、一般/特殊 Wyckoff、单/多元素、混合占位和不同 ASU 大小；来源统一维护在 `docs/EXTERNAL_CIF_SOURCES.md`。
2. 建立受限速保护的官网周期性重抓与版本漂移报告；官网不可用、证据缺失或配对歧义只记为 inconclusive，不计入通过率。
3. 扩展 TOPAS 解析/本体测试及 VESTA、IsoVIZ 分层抽检；未安装商业或图形工具时不得把静态检查表述为本体兼容证明。
4. 进行三次冷/热重复，记录 P50/P95、峰值内存、WSL 进程、ZIP 大小和规范化哈希。
5. 增加多标签、多进程、快速重复提交、ZIP/GenDB 中断和缓存并发测试，并验证服务、端口及 WSL 子进程均能回收。
6. 在看到扩展结果之前冻结发布阈值；固定母相全通过只证明已覆盖范围，不构成对未知结构零缺陷的保证。
