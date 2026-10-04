# 修复与验证报告

本文件只记录已经修复的问题、验证证据、机器报告位置和仍然存在的限制。未完成目标与后续顺序见 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)，人工命令见 [MANUAL_VALIDATION.md](MANUAL_VALIDATION.md)。

## 证据边界

- `experiment_data/`、`webpage_info/`、`output_compare/`、`isobyu/` 始终只读；生产代码不读取官网下载目录来生成答案。
- 官网输出用于差分验证，不作为硬编码表。实现依据空间群仿射作用、直接/倒易格对偶、k-star、小群/表示、Wyckoff 轨道和 fixed-space 等晶体学定义。
- CIF/basis/origin 的字面差异只有在精确有理数、整数幺模变换和 Seitz 共轭给出 witness 后，才记为等价。
- 每份长时报告都带源码、配置、输入、依赖和二进制指纹；签名改变后必须重算，旧报告不能证明当前源码。

## 当前官网与保存输出

### EuAl4/NdNiO2 Method 1/2 保存产物

只读静态审计已再次通过：

- EuAl4 Method 1：123 个候选；Method 2：22 个候选。
- NdNiO2 Method 1：82 个候选；Method 2：48 个候选。
- 合计 275/275 个候选成功配对；235 个 CIF 直接语义匹配，40 个为精确等价 setting。
- 275/275 个保存的 displacive mode 数一致，1100/1100 个核心文件存在且可解析。

源码签名 `e3347be0…5908` 的 live 审计曾从头完成：EuAl4 Method 1 `123/123`、Method 2 `22/22`、NdNiO2 Method 1 `82/82`、Method 2 `48/48`；四批候选和模式数全部一致，missing、extra、duplicate、mismatch 均为 0。一次 OneDrive 瞬时原子 replace 锁只中断报告落盘；同签名 checkpoint 重新枚举候选身份并验证完整指纹后安全续跑，没有复用旧源码结果。报告为 `output/validation/live_method12_report.json`。本轮母相 setting 与物理轨道代码已改变源码签名，因此该报告是对应源码快照的有效证据，不能称为当前源码的最终 live 结论。

### Method 1 官网与现有网页版保存产物

- 只读静态门禁核对官网和现有网页版保存侧各 `330/330` 个候选、各 `1320/1320` 个核心文件；候选身份、目录归属、四格式内部身份、文件类型、非空性和重复身份检查均通过，当前无需重下。
- 其中 4310 两侧各有 `125/125` 个候选、各 `500/500` 个核心文件。官网与本地旧产物之间的 175 项科学差异全部属于 4310：125 项模式数差异、50 项 CIF 语义差异；这是程序修复前的 debug 基线，不是漏下载、误下载或错放。
- `N1+_4D1_SG2` 当前内容正确。由于被覆盖前的旧文件与操作记录未保留，无法从现有完整目录证明此前究竟覆盖到哪个子群位置；历史覆盖目标状态为 `inconclusive`，没有按名称猜测。
- 这项静态结论只证明官网下载材料完整，不能代替 live 计算；当前 live 结论另由下述
  final7 最终源码/输入/运行时签名报告提供。
- 旧 4310 live 报告的 `75/125` 通过、`50/125` CIF 语义失败只描述当时源码；后续生成的 `method1_4310_live_postfix_cif_precision_reanalysis_20261003.json` 使用了错误的精度比较口径，已撤回并失效。两者都不得写成当前代码结果或验收证据。

### Method 3 官网下载

- EuAl4：20/20 个查询、35/35 个 embedding、140/140 个核心文件，0 错误、0 警告。
- NdNiO2：20/20 个查询、42/42 个 embedding、168/168 个核心文件，0 错误；12 条提示仅表示 ND-15～ND-20 的结果页/details 网页包位于非推荐层级。查询上下文、完整首屏表和完整 `(SG,basis,origin,s,i)` 均能唯一配对，证据等级仍为 authoritative。
- 两种母相合计 40 个查询、77 条 authoritative embedding；当前不需要补下载。

对应报告：

- `output/validation/method3_official_download_audit_eual4.json`
- `output/validation/method3_official_download_audit_ndnio2.json`

## 已完成的通用修复

### 晶格、结构与模式基础

1. 周期原子匹配改为真实晶格度量下的 minimum image，并按物种求全局完美匹配；不再用分数坐标逐分量回绕或贪心占位。
2. `load_structure()` 立即固定母相 CIF 的标签、Wyckoff 轨道和原子顺序，避免后续标准化改变导出身份。
3. 超胞与倒格折叠使用 parent/child primitive translation lattice 的精确关系；允许中心化母格的分数 conventional basis，不再对 basis 做 `round()`。
4. 中心化 reciprocal canonicalization 枚举所有满足 centering phase 的邻近倒格平移；修复半整数边界的银行家舍入把同一 P-star 重复计数、继而漏掉 Γ/M 次级模式的问题。EuAl4 P2/P4 的 C1、P1、P3 六项已分别恢复官网模式数 `5/5、2/2、4/4`。
5. 模式独立性使用笛卡尔空间秩而非 `dot > 0.98`；既保留近乎平行但独立的模式，也删除真正线性相关项。

### (3+d) superspace、标签与导出

1. 参数 k 使用精确字符串分数、完整 reciprocal k-star、公度条件 `T k ∈ Z³` 和母相 centered reciprocal 等价；不再猜测 `1/d`，也不重复计算同一物理 star。
2. nmod=0 保留折叠到子胞 Γ 的完整母相模式空间；nmod>0 只保留当前调制谐波和 Γ。ZIP 对每个候选显式传递 nmod，且只复用 nmod 相同的缓存。
3. 自共轭特殊 k 用含平移的 little-group 作用补齐简并分量；参数 k 的正交相位分量保持独立逻辑。
4. 位点模式先通过轨道 transporter 拉回 Wyckoff representative，再按该 representative 的 site point group 分类；分量字母与重复 copy 编号在每个 site-irrep family 内独立分配。
5. Complete modes、IsoVIZ 和 TOPAS 共用归一化：`Ap = As / sqrt(s)`，`dmax = |As| × normfactor × max||uB||`。零振幅与模式基导出已逐格式核对。

### ISO microscopic 位点身份与 fail-closed 标签

1. 位点模式标签不再由单个位移矢量的最大笛卡尔分量或手写 `A1/E/...` 规则猜测。
   程序先把目标 embedding 的 exact invariant direction 化为 primitive-integer
   `VALUE DIRECTION VECTOR,...`，再读取 ISO oriented `DISPLAY DISTORTION` 列；命名 OPD
   token 只校验 primary 子群身份，不代替实际方向。
2. microscopic 列按 `(global irrep, Wyckoff)` 分组，并在公共点域与 BUSH 参考基比较秩、
   独立性和主角度。官网表只打印稀疏代表点时，只有唯一满秩变基
   `R_common T=C_common` 可解、条件数与残差合格，且 `C_full=R_full T` 的完整域摘要可复核，
   才使用延拓列。直接 ISO 列的行/精确向量/规范浮点摘要与延拓证据分别保存，不能把后者
   冒充官网原始输出。
3. 每列同时绑定 parent/k、child SG、精确 `B,q`、primary IR/OPD/selector、物理轨道、
   site irrep/component、frame/atom order 和来源查询顺序。等价位点的整数 transporter
   及已证明相位进入独立摘要；非 Γ 情况缺 component-to-star-arm 相位证据时 fail-closed。
   缺块、重复身份、非唯一/病态变基、不一致空间及纯数值 fallback 均保持 `unresolved`。
4. 稀疏点域 proof 现在还绑定生成它的有序 BUSH 参考基和有序 microscopic 候选列；
   延拓前重算 reference/candidate digest，并逐列复核 source token、`source_order_key`、
   identity、精确行/向量及浮点来源摘要。proof 被复用于重排或换绑后的列会在任何
   same-domain early return 前拒绝；后端同时要求输出索引严格保持 canonical source order。
   参考基换序本身通过 `R_common T=C_common` 的行置换保持结果不变，不再把合法的 BUSH
   基顺序差异误判成 microscopic 身份错位。

### CIF 位移模式模型、精确来源与 writer 合同

1. 新的纯科学模型显式保存 emitted child frame、原子 ID、Seitz 算子、子群轨道、
   transporter、自由坐标、未混合模式列及其正尺度。它以 `R^T G R = G` 检查参考
   lattice metric，对 free-coordinate 列做 max-component-one 规范化，并同步 full atom
   modes、primitive-cell RSS normfactor 与 `As` 语义；符号未决、零列、列混合、重复/
   缺失原子或模式身份都会稳定 fail-closed。
2. ISO microscopic parser 与 mapper 保留 source block/column、精确行与向量摘要、
   查询上下文、frame/atom-order ID 及 signed scalar；可变外层字段、BUSH 行、当前位移列、
   模式来源和目标子群会逐项与冻结身份比较。非对角列混合、错误轨道路由、错误 query
   context 或来源不一致都会拒绝。
3. `DisplaciveExportData` 已成为 CIF、IsoVIZ、Complete modes 和 TOPAS 四个 writer 的
   共享入口；同一已验证模式不会被重复应用，参考坐标和最终坐标保持分离。非零 strain
   与位移共存时，按官网“lattice-coordinate mode 不随 strain 改变”的约定保留位移分数
   坐标，并把最终晶格统一重建为 `B M P`。生产 core 现从目标 Hall setting 构造精确
   emitted frame、稳定 atom ID、物理轨道和父子映射，并只接受与 canonical ISO 列逐项
   相同的 lifted arrays；随后把完整合同交给四个 writer。证据不完整时仍在发布前
   fail-closed。
4. writer 合同现显式保存不可变的 parent/reference/final 结构快照、精确父子嵌入和
   逐原子 `ParentChildSiteMapping`。除参考晶格必须满足 `B P` 外，每个原子还须满足
   `x_parent=x_child B+q`，并通过化学身份、物理轨道、母相 primitive-site class 和完整
   平移陪集核验。IR、OPD、母/子空间群、`B,q`、embedding ID 及由平移格推导的 `s`
    和 primary `VECTOR` 都绑定为不可变身份；调用方不能篡改这些字段、静默改变
    `As→Ap` 或输出错误方向文本。模式 key、显示标签、global IR 与原子顺序只从已验证的
   `ModeIdentity + MicroscopicColumnProvenance` 派生，不能再由平行字符串或数组覆盖。
5. 当前四格式共同合同只接受单物种满占位点；混合/部分占位在所有 writer 都能无损表达
   前保持 fail-closed。批量模式计算或 writer 异常会按候选汇总；磁盘和 ZIP 都在发布前
    完成全批次渲染。新/空根目录直接发布；已有非空根目录发布到内容寻址的
    `.isodistort-batch-v1-<digest>.ready` 目录，并写逐文件大小/SHA-256 manifest；文件锁与
    进程锁保证重叠批次只能整批成功或整批拒绝。成功的空模式表仍与失败区分。
    IsoVIZ 的 atom type 直接来自不可变物理轨道映射；空列、跨类型列或同元素轨道
   错配都会拒绝，不再按标签或元素回退。

### symmetry-adapted 均匀应变模式与导出

1. 应变使用母胞行晶格 `P`、子群 basis `B` 和对称乘子 `M=I+E`，满足 `M(PP^T)M=B^-1(DD^T)B^-T`，并以 `B M P` 重建女儿胞。工程 Voigt 顺序为 `(E11,E22,E33,2E23,2E13,2E12)`，内积权重为 `diag(1,1,1,1/2,1/2,1/2)`。
2. 规范模式来自 ISO rank-[12] 宏观表与目标 embedding 的 `DISPLAY DIRECTION`；程序再用实际嵌入点群和母胞 metric 独立计算 fixed space，并要求整空间维数与 span 一致。缺失或不一致的证据会中止候选导出，不退回坐标轴基或猜 `GM` 标签。
3. CIF 写 max-component-one 的 `q_raw`，IsoVIZ 写 `q_unit=normfactor*q_raw`，Complete modes 同时写 `q_raw`、`q_unit`、`normfactor`、幅度和两种总和。TOPAS 与官网一致，不生成 strain-mode 精修参数，只写由应用应变得到的固定实际晶胞；参考子胞与 `B P` 不一致时明确失败。
4. 定向应变测试覆盖固定空间、官方 F02 数值契约、三种模式导出与 TOPAS 固定晶胞，
   共 89 passed、1 skipped。当前源码另以 WSL 完成 EuAl4 F02 单例 live：1/1 pass，
   metric 重建相对残差约 `1.00e-15`，报告为
   `output/validation/strain_audit_f02_end_to_end.json`。该单例不等于 current-source
   24/24 全量重跑，也不提供跨晶系证据。

### Method 3 精确 direct-lattice 与 coupled fixed-space 产品链路

1. direct basis、Default/P/A/B/C/I/F/R centering 和母相点群轨道均以 `Fraction`/GL(3,Z) 处理；Default 采用目标空间群默认 centering，P 不再冒充 Default。
2. 恒等与非恒等 basis 统一做同格/母群仿射共轭，不再允许恒等 basis 匹配任意更大子格；上传坐标中的母平移格从 identity-rotation 纯平移恢复，不靠 HM 字符猜测。
3. 40 组官网下载结果的 route 审计得到 61 条 `single_ir_exact` 与 16 条 `coupled_ir_required`，0 条 indeterminate、0 错误/跳过；该分类仍作为独立金标准保留。
4. 历史 single-IR 产品比较为 13 组逐字段相同、16 组精确仿射等价、11 组欠枚举。当前产品接入 coupled 链路后，从头重算得到 13 组逐字段相同、27 组精确仿射等价、0 差异、0 本地错误/跳过；40/40 查询的官网/本地候选数逐组相同，77/77 embedding 全部匹配。
5. Stage A 在有限 normalizer 商 `N_G(T_s)/T_s` 中枚举闭合点子群与 affine lifts。全量结果为 40/40 查询、142 个候选，覆盖 77/77 官网 embedding；65 个表面多余项中 62 个是母群共轭多重性，3 个是新轨道。
6. Stage B 在 `strain ⊕ displacive` 表示上用 exact Reynolds/fixed-space 稳定子判据验证可达性。77/77 官网 embedding 可达；142 个 Stage-A 候选中 139 个可达，3 个新轨道均为 `embedding_infeasible`，恰好排除 Stage-A 假阳性。
7. 空间群产品查询现默认运行 Stage A/B。对 route-less 可达 embedding，把枚举到的单-IR 稳定子及其平移陪集提升到目标有限商，用精确交集动态规划寻找见证；只有交集严格等于目标 subgroup 才标为 `exact_fixed_space`。这类行可进入 Method 2，以 nmod=0 生成目标子群全部折叠 k 的完整位移 fixed space；不把某一个见证组合宣称为唯一 primary COPL 分解。无法证明、物理不可达或超出支持表示的诊断行仍不可选择。

### Method 3 签名复验与等价输入扩展

1. Method 4 归一化改动后的历史重跑仍为 13 exact + 16 affine-equivalent + 11 欠枚举。coupled 产品接入后已再次用 `--restart` 从对应签名源码重算正式报告：13 exact + 27 affine-equivalent + 0 differences，官方/本地均为 77 条，0 本地错误/跳过。`affine_equivalent` 只在 exact Fraction/Seitz/母群共轭 witness 成立时记录，不冒充 basis/origin 字面相同；后续源码变化若要更新 current-source 结论必须重跑。
2. `audit_method3_basis_metamorphism.py` 对每个官方查询构造循环换轴和带符号轴旋转两种 GL(3,Z) representative；先以精确有理格选择与原查询相同 primitive lattice 的 Default/显式 centering，再比较稳定 embedding ID 多重集。coupled 产品和逐物种 Stage-B 接入后已在源签名 `ea06c9ca…0234` 下从头完成：40 个案例、80/80 个变体全部通过，候选数逐例一致，missing/extra 均为 0；其中 4 个变体须显式选 `B` centering。
3. 这项变形审计证明已覆盖查询对等价 basis/centering 表示不敏感；它没有新增母相晶系，也没有补齐任意/多参数 k、point-group-only 或 reciprocal-sublattice 能力。

### Method 4 本地分解 debug

1. 原实现以“原子数是否变化”判断是否需要子群胞，导致 determinant=1 的旋转子群胞直接与母相设置比较；参数-k 又忽略已计算的完整 supercell mode vectors。现在始终以所选子群的零幅度 child cell 为参考，并优先使用经形状/有限值校验的 `mode_displacements_sc`。
2. 原 `provided_origin_shift` 只写入 metadata、不参与匹配或 residual；现在会从女儿相分数坐标中显式扣除。网页与终端增加 matching、robust Å 阈值和已知原点输入，三者仍调用同一 `IsoDistort.search_method_4`。
3. 原子匹配从分数坐标贪心改为按物种的全局最小笛卡尔距离分配；robust 阈值明确为 Å。同原子数也必须通过晶格 metric 门禁，错误晶格不再被静默接受。
4. 最小二乘在笛卡尔 Å 空间进行，RMS/max residual 明确以 Å 报告；模式基若秩亏则拒绝不唯一幅度，并在 metadata 记录秩、列数和条件数。
5. 在当时的本地契约下，最终签名清单 24/24 固定案例通过，两个母相的“未准备模式”也稳定拒绝。EuAl4/ NdNiO2 无噪声最大 RMS 分别为 `4.99e-16 Å` / `4.62e-9 Å`；噪声案例为 `1.242e-4 Å` / `3.597e-5 Å`。物种、8% 晶格变化和 0.1 Å 超阈值输入均被本地拒绝；其中 8% 晶格变化后来由官网 F02 证明为均匀应变成功例，故这一项不再计作正确拒绝。
6. G01 官网对照揭示原实现把无量纲/分数模式的最小二乘系数直接命名为 amplitude，未换算官网 `As/Ap`。现在保留 `raw_coefficients`，另由女儿原胞内的笛卡尔模式范数计算 `normfactor`、`As=raw/normfactor` 和 `Ap=As/sqrt(s)`；女儿 CIF 声明的空间群决定 conventional-to-primitive 重数，不能用坐标重新识别对称性后静默提升 P1。
7. EuAl4 G01 已完成官网只读审计：官网结果为 `1 P1`、basis `{(0,1,1),(1,0,0),(0,0,-1)}`、origin `(0,0,0)`、`s=2,i=32`；30 个位移与 6 个应变幅值全零。官网 CIF 与冻结 daughter 的最大物种匹配距离为 `0 Å`，晶格 metric 相对差 `5.0444e-8`；CIF/TOPAS/IsoVIZ 与 details 的模式数一致。
8. G01 中 4 个可唯一配对的本地模式 `normfactor` 与官网五位小数显示值全部一致，最大绝对差 `4.734e-6 Å^-1`。这证明这些模式的方向尺度，但全零输入本身不能验证非零 `As/Ap` 和符号；首次审计还发现缺少上传副本及动态 form 截图/PDF，前者已在 G02 审计前补齐。
9. EuAl4 G02 非零门禁通过：官网唯一非零重叠模式为 `As=0.70432 Å`、`Ap=0.49803 Å`、`dmax=0.35216 Å`；本地 `As=0.7043199999999998 Å`，绝对差 `1.11e-16 Å`，本地 `Ap=0.4980294481253090 Å`，与官网五位小数值差 `5.52e-7 Å`。4/4 个重叠模式的 `As/Ap/normfactor` 均在官网显示精度容差内，上传 daughter SHA-256 也与清单一致。
10. G01 上传副本已经补齐且哈希一致；G01/G02 的普通 HTML 仍不能证明动态 form 值，因未归档填写完成后的截图/PDF而保留非致命提示。结果页、details 和三个导出文件能证明服务器实际接受的 basis 与最终数值，但后续案例必须额外保存控件截图/PDF。
11. EuAl4 G03 负号/周期换像通过：官网 `As=-0.70432 Å`、`Ap=-0.49803 Å`，与本地在官网显示精度内一致；导出结构在周期最短距离下最大差 `1.36e-15 Å`。G04 两模式通过：官网/本地共享模式分别为 `As=(0.52824,-0.78141) Å`、`Ap=(0.37352,-0.55254) Å`。
12. EuAl4 G06 的两个目标主模式通过：官网 `As=(1.58468,2.45578) Å`，与本地差分别为 `4.32e-6`、`9.10e-7 Å`。微噪声在官网完整 30 维 P1 空间与本地受限 4 模式空间中的投影不要求相同；官网五位分数坐标导出的整体原点对齐后最大结构差 `1.25e-4 Å`，低于由五位坐标精度推得的 `1.447e-4 Å` 容差。
13. G05 揭示人工清单的显式-origin 符号错误：daughter 坐标由 `+t_child` 生成，而官网字段要求新超胞原点相对母相的逆平移 `-t_child·B`。EuAl4 正确值是 `(-1/4,-1/8,-1/8)`；自动运行确实找到同一原点并恢复 G04 的两项幅值。通用生成脚本与两种母相清单已修正，24 个 daughter SHA-256 均未改变。
14. G05 已用修正后的显式原点重做并通过：官网完整身份为 `1 P1`、basis `{(0,1,1),(1,0,0),(0,0,-1)}`、origin `(-0.25,-0.125,-0.125)`、`s=2,i=32`；两项非零 `As=(-0.78141,0.52824) Å`、`Ap=(-0.55254,0.37352) Å`。显式与 auto-origin 的 details/CIF/TOPAS/IsoVIZ 数值一致，导出 CIF 对冻结 daughter 的最大匹配距离为 `2.22e-16 Å`。G01–G06 后续补录的 completed-controls 截图均已识别；G01–G04、G06 的 basis HTML 也已归档并完全通过，只有 G05 auto-origin 缺 basis HTML，保留证据警告而非数值失败。
15. EuAl4 P01 参数-k 零超胞通过：官网为 `1 P1`、basis `{(1,0,0),(0,1,0),(0,0,6)}`、origin `(0,0,0)`、`s=12,i=192`，180 个位移模式与 6 个应变模式全零；48 个本地参数-k 模式均能按完整 `(k,IR,site,site-irrep,component)` 身份唯一映射。导出结构最大差 `2.23e-4 Å`，低于官网五位分数坐标精度推得的 `3.98e-4 Å` 容差。参数-k 多维实基的非目标分量可能在同一子空间内发生旋转，审计只对案例目标分量严格比较 normfactor，避免把基约定差异误报成物理差异。
16. P02/P03 首次官网下载明确返回 `Nearest-site method failed.`，没有 details/CIF/TOPAS/IsoVIZ，不能计作成功。诊断发现旧固定原始系数使最大正确原子位移达到 `4.68846/3.90492 Å`；官网 robust 算法会枚举阈值内候选映射，60 原子超胞在 4–5 Å 阈值下组合数不可控。生成器现通过 `10^-4` 小扰动测定线性笛卡尔斜率，再通用缩放到 `dmax=0.30/0.35 Å`；两种母相仅 P02/P03 共 4 个 daughter 哈希改变，其余 20 个不变。新签名 `6e6b9ff0…b039` 曾在旧契约下本地 24/24 通过；F02 官网差分随后使该总状态降级。
17. EuAl4 新 P02/P03 已用精确冻结哈希重新下载并通过：两例均得到 `1 P1 / diag(1,1,6) / s=12,i=192` 与 180 个位移模式，导出结构对 daughter 的最大距离分别为 `2.887e-4/2.963e-4 Å`，低于 `3.979e-4 Å` 坐标舍入容差。P03 暴露审计器曾错误逐分量比较等价二维 LD1 实基；改为仅在目标分量完整覆盖同一模式子空间时比较基不变欧氏范数，并拒绝对部分子空间使用该放宽。P03 本地/官网 `As` 范数为 `0.85744161/0.85744126 Å`，`Ap` 范数为 `0.24752207/0.24752496 Å`，差值 `3.47e-7/2.89e-6 Å` 均低于官网五位显示精度容差 `5.1e-6 Å`；相关审计单测 5/5 通过。
18. EuAl4 F01/F03 失败路径已通过：官网分别返回物种不匹配与 robust `dmax=0.1 Å` 无法匹配，上传 daughter 哈希、basis/origin/matching 截图和错误语义均一致。F02 则被官网成功接受为 `1 P1 / {(0,1,1),(1,0,0),(0,0,-1)} / origin 0 / s=2,i=32`；30 个位移模式全零，6 个应变模式中 `GM1+ strain_1=0.06016`、`GM1+ strain_2=-0.00036`、`GM2+ strain=-0.06016`、`GM5+ strain(a)=0.09683` 非零，导出 CIF 与冻结 daughter 的原子坐标匹配距离为 `0 Å`。这证明 8% 晶格变化是可分解的纯均匀应变，本地 metric 门禁把它当作“无效晶格”拒绝属于能力/测试契约缺陷。
19. Method 4 现按官网开源实现的 `M=I+epsilon` 约定，在母相 conventional axes 中解对称正定方程 `M G_parent M = B^-1 G_daughter B^-T`，用 metric 消除 CIF 任意刚体旋转；Voigt 顺序为 `(xx,yy,zz,2yz,2xz,2xy)`。EuAl4 F02 得到 `(0,0.085083962,-0.000364592,0.136943037,0,0)`，晶格重建相对残差 `1.00e-15`；官网 `data.isoviz` 的幅度乘六维模式向量得到 `(0,0.085079475,-0.00036,0.136937954,0,0)`，最大差 `5.08e-6`，逐分量均落在由两组五位小数推导的舍入界内。NdNiO2 F02 本地恢复 `yy=0.08`、重建残差 `4.81e-16`。冻结参数-k上下文新增 `k_point_label/k_parameters/k_coordinates`，24 例验证可直接复现所选候选而不依赖重新枚举整张 Method 3 表；最终本地 24/24 通过，EuAl4 官网 12 例为 11 pass + G05 1 warning + 0 fail。
20. NdNiO2 G01–G03 官网输入哈希、basis `{(0,1,0),(0,0,1),(1,0,0)}`、automatic origin、nearest-site 与截图证据均通过；三例均返回 `1 P1 / s=1,i=16`，G03 的 origin `(-1,0,0)` 与零原点仅差母格平移。G02/G03 官网目标 `As/Ap=±0.31366 Å`，本地为 `±0.313664 Å`，差 `4.0e-6 Å`；导出 CIF 对冻结 daughter 的匹配距离均为 `0 Å`，六维应用应变差仅为浮点噪声。首次审计因本地氧 `f` 位列名 `a,b,c,d` 与官网同一四维空间的 `B3u(a,b)+B2u(a,b)` 不可逐字符串配对而误报 3 个失败；审计器现只在双方维数完整且相同位点/IR、相同归一化时比较基不变 `As/Ap` 范数，并对不完整子空间保持严格拒绝。回归测试 9/9 通过，重新审计为 G01–G03 3/3 pass。
21. NdNiO2 G04–G06 的冻结哈希、basis、matching、截图和完整导出均通过。G04 两个本地目标 `As=(0.235248,-0.137228) Å` 与官网 `(0.23525,-0.13723) Å` 逐项相差约 `2e-6 Å`；G05 显式 parent-fractional origin `(0,-1/8,-1/4)` 与 automatic-origin 运行找到的同一原点完全一致，details/CIF/TOPAS/IsoVIZ 幅度相同。G06 官网自动原点为 `(0,-1e-5,0)`，两个共享目标的官网/本地差为 `4.37e-6/2.46e-5 Å`；后者超过纯五位显示舍入但小于从原点对齐导出结构推导的 `6.12e-5 Å` 上界。审计器现仅对清单明确标注的带噪输入使用“显示舍入 + aligned-export Cartesian RSS”的 Cauchy–Schwarz 幅度界；无噪声案例仍保持 `5.1e-6 Å` 严格阈值。G06 导出结构最大匹配距离为 `3.82e-5 Å`，低于 `1e-4 Å` 输出精度门限；这次初始 1 个失败是审计模型缺陷，不是下载或分解失败。
22. NdNiO2 P01–P03 的冻结上传哈希、basis `{(-3,0,0),(0,0,1),(0,2,0)}`、automatic origin、nearest-site、填写截图、basis HTML 和全部导出均通过；三例官网均为 `1 P1 / s=6,i=96`、72 个位移模式与 6 个应变模式。P01 全零保持。P02 自动原点 `(1,0,0)` 对 `k=(1/3,1/2,0)` 产生 `k·t=1/3` 相位，本地单一归一化 Y1 分量被官网旋转为 `As=(-0.25981,-0.45000) Å`；官网范数 `0.51961643 Å` 与本地 `0.51961526 Å` 差 `1.17e-6 Å`。P03 原点 `(0,-1,0)` 对 Y1/X1+ 都产生半周期相位，两个分支反号；Y1/X1+ 的 As 范数差分别为 `3.88e-6/4.45e-6 Å`。三例导出 CIF 对冻结 daughter 的最大距离不超过 `4.04e-5 Å`，低于五位坐标推导的 `1.20e-4 Å` 门限。首次审计把非整数 origin phase 后的实分量强行逐项比较而误报 P02/P03；现仅当活动 `k·t mod 1 != 0`、全部本地分支分量都是显式目标、完整官网分支存在且身份无重复时比较归一化分支范数，零相位与不完整分支继续严格逐分量比较。相关审计单测 14/14 通过。
23. NdNiO2 F01 的上传哈希与冻结值一致，官网明确返回 `Types of atoms in subgroup do not match types of atoms in parent.`，basis/automatic-origin/nearest-site 截图与错误页齐全。F02 也完整通过：官网 `1 P1 / {(0,1,0),(0,0,1),(1,0,0)} / origin 0 / s=1,i=16`，原子坐标精确重建，应用应变官网 `yy=0.0800024254`、本地 `yy=0.08`，差 `2.43e-6`，小于 IsoVIZ 五位舍入界。F03 首次截图虽正确选择 Robust 与 `dmax=0.1 Å`，但当时归档 `daughter.cif` 的 SHA-256 为 F02 的 `9d8f28df…e1ba3b`，不是冻结 F03 的 `02b41e63…332e36`；截图也显示 F02 晶格 `a=4.234464 Å`，而正确 F03 应为 `a=3.920800 Å` 且 Nd 坐标 `(0.86,0.81,0.50)`。其 details/CIF/TOPAS/IsoVIZ 与 F02 逐文件相同，故官网成功是重复运行 F02，不是 robust 阈值失效；该运行当时被正确判为 fail 并要求只重跑 F03。
24. NdNiO2 F03 已用正确冻结输入重跑并通过审计：归档 daughter 与机器清单 SHA-256 均为 `02b41e63c770b9c4e5ff0a622429841be5e74478f1ca478c9755bff459332e36`，晶格为 `a=3.920800, b=3.281000, c=3.920800 Å`，Nd 坐标为 `(0.86,0.81,0.50)`。填写截图确认 basis `{(0,1,0),(0,0,1),(1,0,0)}`、automatic origin、robust matching、`dmax=0.1 Å` 且 manual mapping 未选；官网错误页返回 `Failed to find match. Try using a larger value for dmax.`，没有任何成功导出。整批官网报告因此更新为 24/24 已下载且输入有效、23 pass + 1 pass_with_warnings + 0 fail，状态 `pass_with_warnings`；唯一 warning 仍是 EuAl4 G05 auto-origin 缺 basis HTML，不影响其已核定的结果身份与数值结论。

### 4310 FullProf、母相 setting 与物理轨道修复

1. `4310_tetra.cif` 的结构层原本能解析，但完整 `load_structure()` 会卡在源 CIF 原子标签/Wyckoff 映射。根因不是 superspace：`parent_header.parse_cif_atom_site_rows()` 使用跨整文件、带回溯的 loop 正则，在 FullProf 模板的大量非结构 loop 和含下划线值上出现灾难性回溯。
2. 解析器现改为逐行识别 `loop_`、标签和记录，只在找到 `_atom_site_fract_x` 的 loop 后分词；不按 4310 的材料名、标签或坐标写特例。最小 FullProf 多数据块回归及原有 NdNiO2 标签/顺序回归 3/3 通过。
3. 当前完整 API 能识别 34 原子、约化式 La4Ni3O10、`I4/mmm (#139)` 和 8 个独立物理 Wyckoff 轨道；其中 5 个轨道都属于 `4e`。页头现在使用 FINDSYM 的标准代表与参数，按源 CIF 标签顺序显示 8 行，不再把体心等价的非标准坐标直接当作官网代表。
4. 每个物理轨道现有由空间群、Wyckoff 类型、物种和完整等价坐标集合生成的稳定 `orbit_id`。`DistortionMode`、BUSH/smodes 查询、母胞到子胞映射和导出均携带该身份；逐物种 displacive 作用域传入允许的 `orbit_id`，不会再把同字母的 5 个 `4e` 轨道折叠或把未选物种的同字母轨道混入。
5. 若 FINDSYM 检出等价的非标准 basis/origin，程序先用每个物理轨道的临时独立物种把完整母相规范到 conventional standard setting，再恢复真实物种，并复核空间群、轨道多重度和一一对应；页头、搜索、模式和导出随后都使用同一规范结构。无法完成一致变换时明确失败，不再混用两套坐标。
6. 在对应源码签名的定向实测中，原始 4310 的 N1+ 4D1 子群得到官网一致的 192 个位移模式；把母相全部原子整体平移 `z+0.1` 后会触发完整母相规范化，仍得到 125 个 Method 1 候选，N1+ 4D1 仍为 192 个模式。另一个 F-centered primitive 输入会完整恢复为 conventional cell，证明实现不依赖 4310 名称或固定坐标。
7. `method1_4310_live_final6_20261004.json` 完成过 125 个候选的诊断运行，得到
   115 pass + 10 fail、125/125 模式数和 BUSH 覆盖；但运行期间源码签名改变，报告状态为
   `complete-source-drift`，只能用于定位 C8、P5|4D1 和次级模式延拓根因，不能作为最终
   125/125 证据。随后对应定向 live 已得到 C8 `112/112`、P5|4D1 `90/90`、P2|C1
   `18/18`、X1-|C1 `16/16` 个 verified microscopic identity，并成功生成代表 CIF。
8. final7 已在源码/输入/运行时联合签名
   `0f8900bb5add71bc1de5e9d4c562142d8595539daf757dee074bf614c919a3b5` 下从头完成。
   `method1_4310_live_final7_20261004.json` 状态为 `complete-passed`，结束时签名未变化：
   官网/live/配对候选均为 125，125/125 模式数一致，125/125 BUSH 覆盖通过，125/125
   零振幅 CIF 语义通过，0 missing/extra/duplicate/failure。报告运行 3174.65331 秒，
   SHA-256 为 `383ed579c80d8f6252a2268e298ec1ee956b5f7fadaa2698bbf0053c68af9dc1`。
   该结果闭合 4310 Method 1 的当前全量官网矩阵，但不独立证明每个数值模式向量和
   归一化，也不能外推到任意晶体；Method 2–4 仍未建立 4310 的正式官网矩阵。

### 接口、并发与工具链

1. Method 3 embedding 使用稳定内容 ID，不依赖 Python `id()` 或列表下标。
2. 网页状态、搜索、选择和导出共用 `RLock` 与单调 revision；旧标签页的 stale mutation/export 会明确失败，上传暂存名使用 UUID。
3. Method 3 参数值、后端查询、有限商、点子群和 affine lift 均设有“超限即失败”的预算，不做静默截断。
4. ISODISTORT_VALIDATE 按晶格度量、物种占位和可选原子顺序比较 CIF；ISOVIZ_INPUT 对路径、CSV、幅度和启动器进行显式校验。三个项目共用根 `.venv` 与各自 yaml 配置。
5. 取证确认受保护 `isobyu/` 中的 ignored `iso.log` 来自一次旧的裸 WSL 探针：该命令
   显式 `cd` 到二进制/数据目录后运行 `./iso`，时间和命令流与日志完全一致；生产包装器
   始终在 `~/.id/tmp` 运行并用独立 `ISODATA` 链接读取数据。包装器的 shell 链已由分号
   改为 `export ... && cd <stage> && <binary>`，使 stage 不可用时在启动二进制前失败，
   不会退回调用者 cwd 写日志；无 WSL 命令构造回归已通过。现有受保护日志未被清理或改写。

### ZIP 子群目录短名与防覆盖

1. 导出目录名统一由 `isocore/io/distortion_formats.py` 生成：Method 1 为 `<IR>_<OPD>_SG<number>`，Method 2 为 `<IR>_<OPD>`；保留 IR 的 `+`/`-` 和 `4D1` 等 OPD token，并替换 Windows 非法字符、控制字符和保留设备名。
2. Method 3 有外部 `stable_case_id` 时使用其安全化形式；否则从整批候选的完整身份生成与输入顺序无关的 `M3-<10位摘要>`，内部候选使用 `C<序号>_SG<number>`，不含 Eu/Nd 或样例硬编码。
3. 短名碰撞按 Windows 大小写不敏感规则检测，确定性追加候选身份摘要或序号。磁盘导出遇到既有目录会选择新名字并以 `exist_ok=False` 建立，ZIP 也拒绝重复成员路径，因此不会静默覆盖；内部文件名仍为 `subgroup.cif`、`data.isoviz`、`Complete modes details.txt` 和 `topas.str`。
4. 命名、安全字符、碰撞、Method 3 稳定案例号、ZIP/API 与目录防覆盖定向测试为 9 passed；网页 Method 参数透传测试为 1 passed。较宽的格式测试为 40 passed、1 skipped、1 deselected。

### CRIS 0.4.0 部署、文档与清理

1. 根目录新增唯一部署入口 `setup_cris.py`，统一提供 `install`、`download` 和
   `doctor`。它只创建或复用根 `.venv`，支持按子项目安装、开发依赖、离线 wheelhouse
   和 JSON 诊断；三个旧 `main_requirement.py` 仅保留兼容转发，不再复制安装逻辑。
2. 三个包统一为 `0.4.0`，根 `VERSION` 是仓库版本来源。各 `pyproject.toml` 从本项目
   `requirements.txt` 动态读取运行时依赖，避免包元数据和安装脚本维护两份版本约束；
   ISODISTORT 的运行时声明补齐 `scipy>=1.8`，并与代码要求统一为 `spglib>=2.0`。
3. 根 README 已改为新用户入口，覆盖获取代码、系统组件、在线/离线安装、部署诊断、
   三项目快速上手、数据边界、更新和故障排查；三个子项目 README 继续承担各自输入、
   工作流、参数、输出、限制和验证命令。4310 已从输入预检、Method 1 官网静态门禁和
   根因定向回归推进到 final7 125/125 全量 live；该结论仍不能外推到 Method 2–4。
4. 2026-10-02 曾把 26 个手工探针文件、294 个旧顶层输出、历史 checkpoint/中间审计、
   测试临时目录和构建产物迁入 Git 忽略的
   `output/archive/completed-debug-20261002/`。final7 和本次交接证据闭合后，2026-10-04
   已删除其中不再引用的手工探针、旧输出、中间报告和全部测试/项目临时副本，约
   10.99 MB；只保留下条部署证据引用的 `package-builds/` 三个 wheel。
   `output/validation/` 继续只保留最终证据、Method 4 冻结输入、DSH 会话和测试契约
   需要的 `method3_stage_a_all_audit.json`；`output/tmp/` 作为程序所需的空运行目录保留。
5. 原先受旧 `CodexSandboxOffline` ACL 阻挡的 31 个根 `.test-tmp-*` 目录及其他非受保护
   缓存已经在 2026-10-04 修复权限后删除；当前没有根测试临时目录或非受保护的 Python、
   pytest、Ruff、mypy 缓存。`experiment_data/__pycache__/` 位于规则指定的只读黄金输入树，
   为避免改写原始数据范围而未触碰。
6. 三个包均以 `--no-index --no-deps --no-build-isolation` 从该次产品化源码成功构建 0.4.0
   wheel，文件保存在归档的 `package-builds/`。部署诊断结果为：ISODISTORT 14 PASS；
   ISODISTORT_VALIDATE 9 PASS + 1 WARN（尚无用户 CIF）；ISOVIZ_INPUT 8 PASS + 1 WARN
   （尚无用户 CSV，未做 GUI 数据识别）。WSL 探针、`iso`、`smodes`、10 个
   `data_*.txt`、Java 和 IsoVIZ 启动器均已识别。
7. 2026-10-02 部署阶段源码快照的完整回归为 ISODISTORT `370 passed, 4 skipped`
   （在允许 WSL 的会话中运行）、
   ISODISTORT_VALIDATE `19 passed`、ISOVIZ_INPUT `20 passed`、部署器单测 `5 passed`；
   三组改动文件的 Ruff 检查通过。沙箱内曾有 8 个 WSL 用例失败，沙箱外重跑全部通过，
   因此前者只记录为执行环境限制，不记为产品失败。
8. 本机 OneDrive 工作区内的 Python 虚拟环境启动器会被 WSL 拒绝，即使解释器带有效
   Python Software Foundation 签名。final7 执行时曾把环境放在
   `%LOCALAPPDATA%\CRIS\venv-cris-py312`，根 `.venv` 用 junction 保持调用路径；该阶段
   Python 3.12.5 → WSL、标准 `setup_cris.py install --dev` 和 doctor 均通过。
   2026-10-04 按用户要求将环境本体移回 CRIS，根 `.venv` 现为实体目录。需要 WSL 的
   命令由 `run_cris.ps1` 启动 `.venv/pyvenv.cfg` 记录的外部基础 Python，以 `-S` 禁用
   外部解释器的全局 site-packages，并只注入实体 `.venv/Lib/site-packages` 与项目路径；
   实测 `spglib 2.7.0`、网页模块、WSL `/home/devoutwang` 与 `IsoDistort` 初始化均通过。
   安装器的 Windows 子进程环境固定
   启用 UTF-8，避免中文区域设置下 pip 以 GBK 解码 UTF-8 requirements 注释而失败。

## 当前验证矩阵

- 三种母相 Method 1 官网与现有网页版保存产物静态门禁：PASS，两侧各 330/330 候选、各 1320/1320 核心文件；175 项科学差异全部是 4310 旧程序的 125 项模式数差异和 50 项 CIF 语义差异，不是下载错误。4310 N1+ 4D1 当前内容正确，历史覆盖目标为 `inconclusive`。
- EuAl4/NdNiO2 Method 1/2 保存产物：PASS，275/275 候选，1100/1100 核心文件。
- EuAl4/NdNiO2 Method 1/2 live 快照：签名 `e3347be0…5908` 下 PASS，275/275 候选与模式数一致，0 missing/extra/duplicate/mismatch；本轮代码改变后不再称为当前源码结论。
- 4310 Method 1 全量 live：final7 签名为
  `0f8900bb5add71bc1de5e9d4c562142d8595539daf757dee074bf614c919a3b5`，运行结束时签名未变化；
  125/125 候选身份配对、模式数、BUSH 覆盖和零振幅 CIF 语义全部通过，0 失败。页头、
  8 个物理轨道、逐物种轨道作用域、完整母相 setting 规范化，以及原始/整体 `z+0.1`
  输入的 N1+ 4D1 `192/192` 另有定向回归。
- 4310 Method 1 C10/C12 定向 live：报告 schema 为 `method1-4310-c10-c12-live-v2`，
  运行前后源码签名均为
  `ff79ff5e5201c4219379366ecda9d9abc867f9bffcbc20e61e5f7bcdc0bb2cf8`，其中
  `isocore/**/*.py` 摘要为
  `e3c1e1f38863df1c94e245517539c7e9437d9cd53256b357de119ceeb1d32528`。
  N1+ C10 SG15 与 N1+ C12 SG2 各有 96/96 个 raw、verified identity、direct provenance、
  mapped 和 public-contract 模式；各含 48 个 verified domain extension、48 个 unique
  source-order key，以及 96 个 unique identity/amplitude key。C10 的 exact selector/query
  digest 为 `VECTOR,A,B,B,-A` / `09d35014e8f7e0edabc20e49b9ed9fc597621b3ac2c5c2e3d9eef93e93162087`，
  最大归一化延拓残差为 `3.4558256586127476e-08`；C12 为 `VECTOR,A,0,B,0` /
  `9b4f9c83bdac1ab86611a31332f7fdef49d7b5f598815d1e03069b1140a57dc5`，残差为
  `2.663862690288989e-08`。两例 setting 均为 equivalent 且 conclusive，`setting_issues=[]`。
  这是两个代表候选的定向数值证据，与随后完成的 final7 125/125 全量身份、计数、
  BUSH 覆盖和 CIF 语义证据互补。
- ISO microscopic 位点身份：exact `DISPLAY DIRECTION` → primitive-integer `VECTOR` 查询、
  完整子空间验证、稀疏点域唯一延拓、重复/缺失/篡改身份拒绝和 Bloch transporter
  fail-closed 定向回归通过。有序 reference/candidate proof 绑定、BUSH 基换序不变性、
  proof 重排/换绑拒绝和 canonical source-order 后端门禁也已覆盖。非 Γ
  component-to-star-arm 相位证据缺失时仍拒绝；final7 已在冻结签名下完成 125/125
  全量 live，但该批次不把缺失证据时的 fail-closed 限制改写为已实现能力。
- CIF 位移模型/来源：生产 core 已构造精确 frame、轨道/free-row、canonical scale、
  primitive norm、未混合 provenance、逐原子父子映射、完整查询/子群身份和不可变 writer
  快照。批量 writer 异常、IsoVIZ 轨道错配、错误 query context、模式/方向篡改和非满占位点
  均 fail-closed；磁盘/ZIP 发布前全批次渲染，非空目标采用内容寻址 ready 目录和 manifest。
  非零 strain + displacive 四 writer 组合已通过定向回归。
- symmetry-adapted 应变：`q_raw/q_unit/normfactor`、CIF/IsoVIZ/Complete modes 共用契约及 TOPAS 固定实际晶胞定向测试 89 passed、1 skipped；当前源码 EuAl4 F02 单例 WSL live 为 1/1 pass。24/24 全量重跑与跨晶系扩展仍为 `pending`。
- Method 1/2 已记录的签名数值语义抽检：3/3 案例通过；Complete modes、IsoVIZ、TOPAS、CIF 共 12/12 类导出通过。模式数分别为 Eu LD1 `48/48`、Nd Y1 `24/24`、Eu X4− `5/5`，且标签 family、笛卡尔子空间、normfactor、As/Ap/dmax 全部满足声明的判据。
- Method 3 官网下载：PASS，40/40 查询，77/77 embedding，308/308 核心文件。
- Method 3 route：61 single-IR + 16 coupled-only，0 不确定。
- Method 3 默认产品比较：13 exact + 27 affine-equivalent + 0 differences；77/77 embedding 匹配，逐组候选数一致，0 错误/跳过。
- Method 3 Stage A：40/40，77/77 官网 embedding 覆盖。
- Method 3 Stage B：77/77 官网 embedding 可达；3/3 新轨道不可达。
- Method 3 等价 basis/centering 变形：PASS，40 个案例、80/80 个变体，0 missing/extra。
- Method 4 本地闭环：重冻结清单 24/24 通过；官网 24/24 个案例均使用正确冻结输入，审计为 23 个完全通过、EuAl4 G05 auto-origin 1 个证据警告、0 个失败，综合状态 `pass_with_warnings`。F01 物种拒绝、F02 均匀应变成功与 F03 robust 距离阈值拒绝路径均通过。
- 网页/终端对齐：Method 3 均显示并导出 `route_resolution` 与完整 known routes；Method 4 均显示 `As/Ap/raw/normfactor`、residual 与六分量均匀应变，终端原点输入采用与网页相同的空分量补零语义。终端 ZIP 和目录导出均与网页一样显式传递当前 nmod，不再由可变会话状态隐式决定参数 k 模式。网页女儿相上传文件在每次 Method 4 成功或失败后删除，当前母相上传文件在替换或服务退出时删除。
- 实体 `.venv` 迁移后的启动链：Windows 网页启动改由隐藏的 PowerShell helper 调用 Windows Shell；在首选端口被独占时，受支持入口 `.\run_cris.ps1 ISODISTORT\main_web.py` 实测自动打开 `http://127.0.0.1:8001/`，服务收到根页面、静态资源、初始化 API 和持续心跳，`/api/state` 返回 200。直接 `.venv\Scripts\python.exe` 的限制会沿普通子进程和 ShellExecute 进程链保留，不能在 Python 内自愈；网页与终端入口因此在启动服务/交互前 fail-closed，并打印唯一正确的 `run_cris.ps1` 命令。
- 此前生产源码快照的全量回归：654 passed、3 skipped、0 failed，用时 1003.51 秒；
  全套在真实 WSL/ISO 可用的非 sandbox 环境执行，skip 为可选能力门禁。该回归覆盖
  microscopic 查询上下文/完整域延拓、生产 CIF 位移合同、四 writer、批量原子发布及既有
  功能，但早于本轮 proof 绑定加固，不能冒充当前最终源码全量结果。
- 当前源码的目标回归为 160 passed、1 skipped；完整 `ISODISTORT/tests_dev` 为
  669 passed、3 skipped、0 failed、1305 warnings，用时 1108.93 秒。WSL/ISO 用例在普通
  用户上下文中实际执行；此前 630 passed、32 skipped、10 个 `Wsl/E_ACCESSDENIED`
  失败只保留为旧环境诊断，已被本次结果取代。`setup_cris.py` 与 proof 绑定、二进制
  启动链等改动文件的 Ruff 检查通过。
- DSH 只读诊断会话
  `output/validation/dsh_sessions/20261004-122244-method1-identity-review-stdin/` 保存完整 prompt、
  metadata、stdout/stderr；Session ID 为 `session-d129b86d-2b15-43da-9b12-dd5a133966e1`。
  DeepSeek 的 72 passed、2 skipped 及候选原因只作诊断材料；BUSH 排序疑点经矩阵方向和反例
  复核后排除，邻近的 proof 重绑定缺口由 Codex 独立复现、修复并以上述目标回归验证。
- ISODISTORT_VALIDATE：19 passed。
- ISOVIZ_INPUT：20 passed。
- 统一部署器：6 passed；三项目 0.4.0 wheel 均构建成功；本轮全项目 `doctor --dev` 为 24 pass、2 个缺少用户输入的 warn、0 fail。

关键机器报告为 `method3_embedding_route_audit.json`、`method3_official_local_comparison.json`、
`method3_basis_metamorphic_audit.json`、`method3_stage_a_diagnostic_audit.json`、
`method3_stage_b_feasibility_audit.json`、`method4_local_validation.json`、
`method4_official_audit.json`、`strain_audit_f02_end_to_end.json`、`live_method12_report.json`、
`audit_method1_official_downloads_20261003.json`、
`audit_method1_overwrite_forensics_20261003.json`、`method1_4310_live_final6_20261004.json`、
`method1_4310_live_final7_20261004.json` 和
`method1_4310_vector_c10_c12_extension_final_20261004/summary.json`；全部位于
`output/validation/` 且不提交 Git。长时命令见 [MANUAL_VALIDATION.md](MANUAL_VALIDATION.md)。

## 仍保留的限制与证据缺口

- Method 3 空间群查询已连接 coupled embedding、逐物种 displacive Stage-B 表示与完整 fixed-space 位移模式；只选择部分物种时不会把相对共同位移误删为全晶体刚体平移。但不声明稳定子交集见证中的某一组 IR 是唯一 primary COPL 分解；arbitrary/multi-parameter k、point-group-only affine 枚举，以及 rotational/occupational/magnetic 的同等级 Stage-B 表示仍未实现。
- Method 3 reciprocal-sublattice 输入未实现；rotational-only 仍借用位移活性作近似筛选。
- 当前 site-irrep 标签只接受 exact `VECTOR` 查询和完整 BUSH/microscopic 子空间核验的 ISO
  oriented 身份；稀疏表必须经过唯一、满秩且数值稳定的完整域延拓。纯数值 fallback、
  缺少次级 IR 精确 k 身份，以及非 Γ 非零 transporter 缺 component-to-star-arm 证据时均
  保持 `unresolved`。轴矢量、磁和占位模式尚无同等级身份链。
- 生产 core 已接入 CIF 位移合同和四个 writer，并在 mapper、公共 I/O 与发布边界重验来源
  查询上下文。合同暂只支持单物种满占位点；混合/部分占位与更多晶系仍无同等级生产证据。
- 延拓证据是进程内不可变数据合同，不是密码学签名。公共校验器可从 coefficient columns
  重建 `T` 并复核摘要、秩、条件数、列完整性、逐物理轨道分组和 full-basis 摘要；公共域
  原始矩阵未跨层保存，因此 reference/canonical condition、公共域残差和 transport phase
  只能核对构造期摘要。若将来需要抵抗能同步伪造全部字段的外部不可信输入，必须额外保存
  可重放 transcript/矩阵或使用签名。
- 4310 Method 1 已有官网静态全集、对应源码签名的定向根因回归和 final7 125/125
  候选全量 live；Method 2–4 尚无该母相正式矩阵。
- 官网 Method 1–3 的默认子群导出为 `As=0`，因此这些路径的非零振幅数值归一化主要由公式、单元测试和本地生成验证；Method 4 已有非零官网对照，nmod=1 尚无单独归档的官网四格式非零参考对。
- 本轮已识别 Java 与 IsoVIZ 启动器，但按自动验证边界没有启动 GUI；20 个静态测试只能
  证明 `.isoviz` 生成、路径和启动参数逻辑，不能证明窗口实际打开或 IsoVIZ 正确认出数据。
- Method 4 官网 EuAl4 与 NdNiO2 的 G01–G06、P01–P03、F01–F03 已验证零保持、正负幅值、原子重排/周期换像、完整 Γ/参数-k 等价多维模式基、双模式、显式/自动原点、参数-k origin phase、超胞零/非零分解、P1 噪声分配、应用均匀应变张量和拒绝路径。EuAl4 G05 auto-origin 仍仅缺 basis HTML；清晰填写截图、accepted result identity、完整导出和显式-origin 数值对照已足以维持 `pass_with_warnings`。symmetry-adapted strain 的生成、标签和 CIF/IsoVIZ/Complete modes 导出以及 TOPAS 固定实际晶胞现已实现；当前源码只完成 EuAl4 F02 单例 WSL live，尚未重跑 24/24。官网按子群算子自动枚举全部允许原点、occupancy/magnetic/rotational 分解仍未验收。

## 主要科学依据

- [ISODISTORT Method 3 help](https://iso.byu.edu/isodistorthelp.php)
- [ISOTROPY user documentation](https://iso.byu.edu/isotropy_doc.php)
- Campbell, Stokes, Tanner & Hatch, *J. Appl. Cryst.* 39 (2006), [DOI 10.1107/S0021889806014075](https://doi.org/10.1107/S0021889806014075)
- Stokes & Campbell, *Acta Cryst.* A73 (2017), Appendix B, DOI `10.1107/S2053273316017629`
- Hatch & Stokes, *Phys. Rev. B* 65 (2002), DOI `10.1103/PhysRevB.65.014113`
- Stokes, van Orden & Campbell, *J. Appl. Cryst.* 49 (2016), DOI `10.1107/S160057671601311X`

这些文献定义算法与不变量；官网结果只用于检验实现是否复现同一科学语义。
