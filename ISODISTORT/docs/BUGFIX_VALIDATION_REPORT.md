# 修复与验证报告

本文件只记录已经修复的问题、验证证据、机器报告位置和仍然存在的限制。未完成目标与后续顺序见 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)，人工命令见 [MANUAL_VALIDATION.md](MANUAL_VALIDATION.md)。

## 证据边界

- `experiment_data/`、`webpage_info/`、`output_compare/`、`isobyu/` 始终只读；生产代码不读取官网下载目录来生成答案。
- 官网输出用于差分验证，不作为硬编码表。实现依据空间群仿射作用、直接/倒易格对偶、k-star、小群/表示、Wyckoff 轨道和 fixed-space 等晶体学定义。
- CIF/basis/origin 的字面差异只有在精确有理数、整数幺模变换和 Seitz 共轭给出 witness 后，才记为等价。
- 每份长时报告都带源码、配置、输入、依赖和二进制指纹；签名改变后必须重算，旧报告不能证明当前源码。

## 当前官网与保存输出

### Method 1/2 保存产物

只读静态审计已再次通过：

- EuAl4 Method 1：123 个候选；Method 2：22 个候选。
- NdNiO2 Method 1：82 个候选；Method 2：48 个候选。
- 合计 275/275 个候选成功配对；235 个 CIF 直接语义匹配，40 个为精确等价 setting。
- 275/275 个保存的 displacive mode 数一致，1100/1100 个核心文件存在且可解析。

最终源码 live 审计已从头完成：EuAl4 Method 1 `123/123`、Method 2 `22/22`、NdNiO2 Method 1 `82/82`、Method 2 `48/48`；四批候选和模式数全部一致，missing、extra、duplicate、mismatch 均为 0。报告签名为 `e3347be0…5908`。一次 OneDrive 瞬时原子 replace 锁只中断报告落盘；同签名 checkpoint 重新枚举候选身份并验证完整指纹后安全续跑，没有复用旧源码结果。报告为 `output/validation/live_method12_report.json`。

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

### Method 3 精确 direct-lattice 与 inverse-Landau 诊断

1. direct basis、Default/P/A/B/C/I/F/R centering 和母相点群轨道均以 `Fraction`/GL(3,Z) 处理；Default 采用目标空间群默认 centering，P 不再冒充 Default。
2. 恒等与非恒等 basis 统一做同格/母群仿射共轭，不再允许恒等 basis 匹配任意更大子格；上传坐标中的母平移格从 identity-rotation 纯平移恢复，不靠 HM 字符猜测。
3. 40 组官网下载结果的 route 审计得到 61 条 `single_ir_exact` 与 16 条 `coupled_ir_required`，0 条 indeterminate、0 错误/跳过。默认产品 single-IR 路径与 61 条身份一一匹配，未伪造 16 条 coupled route。
4. 最终源码下的本地—官网比较为 13 组逐字段相同、16 组精确仿射等价、11 组欠枚举；61/61 条本地 embedding 均匹配 single-IR 官网身份，0 条多余，缺少的 16 条恰为 coupled-only，0 本地错误/跳过。
5. Stage A 在有限 normalizer 商 `N_G(T_s)/T_s` 中枚举闭合点子群与 affine lifts。全量结果为 40/40 查询、142 个候选，覆盖 77/77 官网 embedding；65 个表面多余项中 62 个是母群共轭多重性，3 个是新轨道。
6. Stage B 在 `strain ⊕ displacive` 表示上用 exact Reynolds/fixed-space 稳定子判据验证可达性。77/77 官网 embedding 可达；142 个 Stage-A 候选中 139 个可达，3 个新轨道均为 `embedding_infeasible`，恰好排除 Stage-A 假阳性。
7. Stage A/B 仍是默认关闭的诊断能力。没有真实 IR/COPL route 和模式基的行保持 `selectable=false`，网页、终端和 Python API 均拒绝继续导出。

### 接口、并发与工具链

1. Method 3 embedding 使用稳定内容 ID，不依赖 Python `id()` 或列表下标。
2. 网页状态、搜索、选择和导出共用 `RLock` 与单调 revision；旧标签页的 stale mutation/export 会明确失败，上传暂存名使用 UUID。
3. Method 3 参数值、后端查询、有限商、点子群和 affine lift 均设有“超限即失败”的预算，不做静默截断。
4. ISODISTORT_VALIDATE 按晶格度量、物种占位和可选原子顺序比较 CIF；ISOVIZ_INPUT 对路径、CSV、幅度和启动器进行显式校验。三个项目共用根 `.venv` 与各自 yaml 配置。

## 当前验证矩阵

- Method 1/2 静态官网差分：PASS，275/275 候选，1100/1100 核心文件。
- Method 1/2 最终源码 live：PASS，275/275 候选与模式数一致，0 missing/extra/duplicate/mismatch；签名 `e3347be0…5908`。
- Method 1/2 当前源码数值语义抽检：3/3 案例通过；Complete modes、IsoVIZ、TOPAS、CIF 共 12/12 类导出通过。模式数分别为 Eu LD1 `48/48`、Nd Y1 `24/24`、Eu X4− `5/5`，且标签 family、笛卡尔子空间、normfactor、As/Ap/dmax 全部满足声明的判据。
- Method 3 官网下载：PASS，40/40 查询，77/77 embedding，308/308 核心文件。
- Method 3 route：61 single-IR + 16 coupled-only，0 不确定。
- Method 3 默认产品比较：13 exact + 16 affine-equivalent + 11 欠枚举；61/61 single-IR 匹配，16 条 coupled-only 未伪造，0 多余/错误/跳过。
- Method 3 Stage A：40/40，77/77 官网 embedding 覆盖。
- Method 3 Stage B：77/77 官网 embedding 可达；3/3 新轨道不可达。
- ISODISTORT：330 passed、4 skipped；skip 为环境/可选能力门禁，0 failed。
- ISODISTORT_VALIDATE：19 passed。
- ISOVIZ_INPUT：12 passed。

关键机器报告为 `method3_embedding_route_audit.json`、`method3_official_local_comparison.json`、`method3_stage_a_diagnostic_audit.json`、`method3_stage_b_feasibility_audit.json` 和 `live_method12_report.json`；全部位于 `output/validation/` 且不提交 Git。长时命令见 [MANUAL_VALIDATION.md](MANUAL_VALIDATION.md)。

## 仍保留的限制与证据缺口

- 完整 Method 3 产品还缺 coupled-IR 的 IR/COPL 分解、route/modes 连接，以及 arbitrary/multi-parameter k 和 point-group-only affine 枚举；默认首屏仍是 single-IR 子集。
- Method 3 reciprocal-sublattice 输入未实现；rotational-only 仍借用位移活性作近似筛选。
- 当前 site-symmetry 标签实现了 polar-vector representative/transporter 分类，但尚未覆盖所有 site group 与全部模式类型的通用 character-table decomposition。
- 官网默认导出均为 `As=0`，因此非零振幅的数值归一化主要由公式、单元测试和本地生成验证；nmod=1 尚无单独归档的官网四格式非零参考对。
- GUI 可打开性在未安装 VESTA/IsoVIZ 的环境中只能做格式/解析静态验证。
- Method 4 只准备并哈希了 24 个 daughter CIF；按用户门禁尚未运行、下载或 debug。

## 主要科学依据

- [ISODISTORT Method 3 help](https://iso.byu.edu/isodistorthelp.php)
- [ISOTROPY user documentation](https://iso.byu.edu/isotropy_doc.php)
- Stokes & Campbell, *Acta Cryst.* A73 (2017), Appendix B, DOI `10.1107/S2053273316017629`
- Hatch & Stokes, *Phys. Rev. B* 65 (2002), DOI `10.1103/PhysRevB.65.014113`
- Stokes, van Orden & Campbell, *J. Appl. Cryst.* 49 (2016), DOI `10.1107/S160057671601311X`

这些文献定义算法与不变量；官网结果只用于检验实现是否复现同一科学语义。
