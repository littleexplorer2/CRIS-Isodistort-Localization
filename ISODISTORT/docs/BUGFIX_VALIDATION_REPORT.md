# 修复与验证报告

本文件只记录已经修复的问题、验证证据、机器报告位置和仍然存在的限制。未完成目标与后续顺序见 [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)，人工命令见 [MANUAL_VALIDATION.md](MANUAL_VALIDATION.md)。

## 证据边界

- `experiment_data/`、`webpage_info/`、`output_compare/`、`resources/isobyu/` 始终只读；生产代码不读取官网下载目录来生成答案。
- 官网输出用于差分验证，不作为硬编码表。实现依据空间群仿射作用、直接/倒易格对偶、k-star、小群/表示、Wyckoff 轨道和 fixed-space 等晶体学定义。
- CIF/basis/origin 的字面差异只有在精确有理数、整数幺模变换和 Seitz 共轭给出 witness 后，才记为等价。
- 每份长时报告都带源码、配置、输入、依赖和二进制指纹；签名改变后必须重算，旧报告不能证明当前源码。

### 2026-10-07 重构修复后最终复验

本轮已恢复迁移时遗失的校验模块和此前实现，修复路径、包边界、安装布局及计算核，
未提交或推送。最终回归在真实 Windows/WSL/ISO 环境执行，不以 sandbox 下的访问失败
冒充程序缺陷或通过结果。

- **完整回归**：`1027 passed, 2 skipped, 0 failed`，0 collection error，
  `1911 warnings`，用时 `1416.76s`。命令从 CRIS 根目录执行：
  `.\run_cris.ps1 -m pytest ISODISTORT\tests -q --tb=line --basetemp ISODISTORT\output\validation\pytest-temp-20261007-final-v6 --junitxml ISODISTORT\output\validation\post_restructure_pytest_20261007_final_v6.xml`。
  JUnit 记录 1029 项、0 errors/failures；两项 skip 是缺少官网/本地 X4− P3 对照 CIF，
  以及 Windows 无法运行原生 POSIX shell 集成测试。原先的合成 CIF 缺失与官网 HTML
  路径假 skip 已处理，并实际运行相关测试。告警包含 spglib API 弃用和 CIF 源数据提示，
  未隐藏，也不视为全部输入数据已无瑕疵。
- **冻结证据**：运行前后 139 个源码/测试/配置/入口文件 SHA-256 均为
  `6dcf96ab945603838b4117c4c0400f5d141d1985dbfef51ee7ce9266d63e5354`；
  30 个合成 CIF 的前后摘要均为
  `17834254c5b8b7d5ce88626e06515a82861099fe2c0da9cc53348538aad54797`。
  文档和输出不属于这个源码摘要。文件范围、算法、依赖、输入哈希、命令与 JUnit 摘要见
  `output/validation/post_restructure_environment_20261007_final.json`。
- **关联回归**：`ISODISTORT_VALIDATE/tests_dev` 19 passed，
  `ISOVIZ_INPUT/tests_dev` 20 passed（静态检查，未打开或验收 GUI），根部署器测试
  10 passed。对应 JUnit 为 `post_restructure_validate_20261007_final.xml`、
  `post_restructure_isoviz_20261007_final.xml`、`post_restructure_setup_20261007_final_v4.xml`。
  Ruff 全树与根部署器通过；tracked `git diff --check` 通过，另对 132 个未跟踪源码/配置
  文件执行 CRLF-aware no-index 空白检查通过，未改写换行符。
- **环境诊断**：最终 `setup_cris.py doctor --project all --dev --json` 为
  24 pass、2 warn、0 fail；警告仅为对比目录缺 CIF 与幅度目录缺 CSV，GUI 未验证。
  报告为 `output/validation/post_restructure_doctor_20261007.json`。
  22 个官网运行资源与迁移前 Git blob 内容完全一致，见
  `output/validation/post_restructure_readonly_resources_20261007.json`；黄金输入未改。
- **安装包**：隔离构建的
  `output/validation/wheel-20261007-final-v6/isodistort-0.4.0-py3-none-any.whl`，SHA-256
  `ce8e9816e266f6dfc10161ebc927abf2096ab4d7053484a1dad518c9f0f568b7`。
  92 个归档条目中，87 个源码/资源与冻结工作树逐字节相同，包含全部 22 个官网资源；
  CRC、白名单和 console entrypoint 检查通过。独立解包可导入网页入口，运行目录在包外，
  新坐标系门禁与极小波矢保留烟测通过。烟测借用共享 `.venv` 依赖，未执行完整 pip
  安装或 sdist 构建，也不等于原生 Linux 全链验收。

本节只证明已执行的有限回归及精确不变量。EuAl4 官网 OPD 的 123 对精确嵌入核验详见
下文第 21 项；它不证明全部模式幅度、所有晶体或所有模式类型。重构前的 Method 1–4
长时矩阵保留为下述历史源码快照；新源码的全矩阵、跨晶系及未支持来源合同仍按开发
计划验收，不能用本轮单元/定向测试计数替代。

### 2026-10-07 E 盘迁移准备复验

- GD 的 CRIS 根改为 `CRIS_ROOT` 环境变量或同级 `CRIS/`；GD 与 ISOVIZ_INPUT 的振幅
  目录改为 `BEST_MODEL_PARAMETERS_DIR`、同级 `Best_Model_Parameters/`、历史 Desktop
  的有序回退。配置、生产者、消费者、单元测试与 README 已同步，不把 `E:` 写死进源码。
- 迁移适配后的真实 Windows/WSL/ISO 全量回归为 `1027 passed, 2 skipped, 0 failed`，
  `1911 warnings`，用时 `1303.83s`。JUnit 为
  `output/validation/migration_e_drive_pytest_20261007.xml`，SHA-256 为
  `94fd15554ef67ceae9fe4fa263bcd27a7172529bf93784b1ffdf80d45553c5e5`；1029 项、
  0 errors/failures。两个 skip 的边界与上节相同。
- 关联门禁：ISOVIZ_INPUT `22 passed`，GD 转换/CSV 路径测试 `7 passed`，
  ISODISTORT_VALIDATE `19 passed`，根部署器 `10 passed`；变更文件的 Ruff 检查通过。
  未安装 pytest 的 GD 独立 TensorFlow `.venv` 没有冒充已执行测试，GD 测试实际使用
  CRIS `.venv` 中已声明的开发依赖。
- `setup_cris.py doctor --project all --dev --json` 为 24 pass、2 warn、0 fail；两个
  warning 仍是 compare 无成对 CIF 与共享振幅目录无 CSV，GUI 未验收。
- 已删除可访问的 Python/pytest/ruff 缓存、构建/egg-info、错误落入 resources 的运行
  产物及 pytest 临时树，共 2891 个文件、约 180 MiB。另有 54 个历史 pytest 目录的
  Windows ACL 同时拒绝当前用户、takeown 与 WSL root；它们仍位于 `ISODISTORT/output/`
  或其 `validation/` 子目录，不能宣称已清除。迁移前需由有权访问其创建者 SID 的账户或
  提升的 Windows 管理员删除；黄金数据、官网资源、科学 JSON/XML 报告与用户 CSV 未动。

### 2026-10-06 历史源码最终复验

本节对应重构前记录的源码签名。2026-10-07 目录迁移、恢复与计算核修改后，以下长时
官网/live 报告仍是历史快照，不能作为本轮最终源码已通过同一完整矩阵的证明。

- **Method 1 / 4310**：`output/validation/method1_4310_current_audit_20261006.json`
  为 schema 3，状态 `complete-passed`；报告 SHA-256 为
  `be7bee90ea47a03bb0aac56debb455a882b0d5d198c593130bd6b11db2f9fa6b`，运行前后
  源码/输入/运行时签名均为
  `96ff44cd307432ae2d8ca672afed6ff7c9e7dedbfc39ee07cfd71fefeaf94d12`。
  官网/live 候选身份、模式数、BUSH 覆盖和零振幅 CIF 语义均为 `125/125`，0 失败。
  这证明该有限候选矩阵及报告所列判据，不独立证明每个数值模式向量、归一化或
  primary IR/OPD 分解的唯一性，也不能外推到其他晶体。
- **Method 2 / 两个参数-k代表例**：
  `output/validation/method2_numeric_semantic_live_current.json` 为 schema 6，报告
  SHA-256 为 `c4830c1895681f19031c6c41113ea5d4039e14e1af2dbad06fe05661edc65cfb`。
  EuAl4 `LD1 C1, g=1/6, nmod=0` 与 NdNiO2 `Y1 C1, a=1/3, nmod=0` 均通过，
  每例的 CIF、IsoVIZ、Complete modes 和 TOPAS 四种导出均为 pass，运行期间源码与
  官网/母相证据稳定。`2/2` 只代表这两个 nmod=0 案例，不代表全部 Method 2 候选、
  nmod=1、4310 或跨晶系矩阵。
- **Method 3 / 双母相官网集合**：
  `output/validation/method3_official_local_comparison_20261006_current.json` 为 schema 5，
  报告 SHA-256 为
  `796ae2f39235c6d61422cb9d02403da11e795a60815f6a0590ea02b45c6c45f7`；40/40 个
  authoritative 查询得到 13 exact + 27 affine-equivalent，0 difference、0 error，
  且运行前后 source signature 一致。审计器把 manifest 与每个官网案例先捕获为不可变
  字节快照，签名与解析消费同一份内容；结束时重新签名，发生漂移就使报告与 checkpoint
  失效，不能把中途变化的归档拼成一次通过。范围仍限于这两种四方母相的 40 个查询。
- **Method 4 / 双母相冻结矩阵**：
  `output/validation/method4_local_validation_20261006_current_v2.json` 为 schema 3，
  当前源码本地闭环 `24/24` 通过，报告 SHA-256 为
  `e19a005dd77bad4b12c7dcc2600889c8e89fe3a689d572c542e6f626a954bb5d`；
  `output/validation/method4_official_audit_20261006_current_v2.json` 为 schema 5，
  报告 SHA-256 为
  `946dba00162ca813060153c4f38daf3131030ffe352525bcc439a2062c5f9519`，结果为
  23 pass + EuAl4 G05 auto-origin 1 warning + 0 inconclusive/fail。warning 仍只表示缺少
  basis HTML；该案例的其余身份和数值证据已通过。这个 `24/24` 只覆盖 EuAl4/NdNiO2
  冻结案例，不构成 4310 或跨晶系验收。
- **完整回归**：最终 `ISODISTORT/tests` 为 `884 passed, 2 skipped, 0 failed`，
  `1343 warnings`，用时 `1263.04s`；命令经根 `run_cris.ps1`，并把 pytest 临时目录
  显式放在 `ISODISTORT/output/validation/`，以允许 WSL 与本机 HTTP 合同测试正常运行。
  全量后仅调整“部分成功 ZIP”的完成提示文案；最终 `test_web.py` 再全量复验为
  `49 passed, 0 failed`、`288 warnings`，用时 `327.35s`。

### 2026-10-06 Method 2 参数-k批量交付边界修复

1. `compute_parametric_modes` 仍以 SMODES 构造数值 fixed-space，但 ISO microscopic
   canonicalization 失败时，不再只写顶层 note；实际异常类别/文本及
   `(parent SG, k, KVALUE, IR, direction selector)` 同时写入每一条仍未解析的
   `ModeIdentity.reason`。API、网页与导出机器报告因此使用同一诊断来源。
2. authoritative writer 门禁改为只核验该候选实际 emitted 的位移列。缓存中没有被该
   writer 使用的模式不能再否决候选；实际写出的每列仍必须具有 verified
   ISO microscopic identity 与完整 provenance，缺失源记录、列混合或数组不一致继续
   fail-closed。
3. Method 2 补算批次把 `UnresolvedModeIdentityError` 作为逐候选不可交付状态处理：已验证
   候选继续发布，跳过候选及完整科学身份/原因写入根级
   `export_candidate_status.json` 和 `export_candidate_status.txt`。ZIP 与磁盘都把成功候选
   和报告先完整渲染，再做一次原子发布；其它模式计算或 writer 异常仍整批拒绝，未放宽
   科学门禁。
4. `search_method_2` 返回同一门禁的 `mode_identity_status`、`export_ready` 与失败明细；
   网页在已计算候选为 unresolved 时显示具体警告并设 `selectable=false`。批量导出提供
   当前候选、完成/成功/不可交付/失败数与剩余时间估计；取消只在候选边界生效，保留已
   完成候选并记录未运行候选。
5. 网页在 `_SESSION_LOCK` 内二次复核 revision、复制结构/候选/模式状态并创建带独立后端
   的 export worker，随后在锁外执行长时计算；因此实时会话可继续使用，且修改不会混入
   已启动批次。同一时刻仍只允许一个批量导出任务。
6. 定向证据覆盖详细身份原因传播、只核验 emitted 列、未知身份逐候选跳过、ZIP/磁盘
   部分成功报告、原子发布、候选边界取消、进度/ETA、真实快照隔离、HTTP 状态/取消接口
   及长任务释放会话锁；末端界面会按最终状态显示 `ready/total`，不会把部分成功写成全部
   完成。完整回归与最终 Web 复验计数见本节上方；本轮没有执行超过 50 分钟仍停在首个
   候选的 NdNiO2 48 候选长时矩阵，因此不能把这些回归写成 48/48 科学验收，也不能据此
   确认 BUG-005 的大型超胞秩/条件数假设。

### 官网 HTML 改名与内容身份

- 当前 Method 1/2/4 归档审计按 HTML 结构、页面角色和内容 SHA-256 识别页面，不把
  basename 当作科学身份；在原候选目录内任意改名，只要后缀仍为不区分大小写的
  `.html` 或 `.htm`，无需重新下载。
- 每个要求的页面角色必须在规定目录内恰好匹配一页；0 页为 `missing_match`，多页为
  `ambiguous_match`，两者都 fail-closed。把文件移到另一个候选目录会改变科学归属，
  改成非 HTML 后缀会使文件不再进入清单，这两种操作不属于安全的“仅改文件名”。
- 这套解析只用于只读验证证据。生产搜索、四格式生成和 ZIP 下载均从当前计算状态生成，
  不读取 `webpage_info/` 或 `output_compare/` 归档来制造结果。

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
- 这项静态结论只证明官网下载材料完整，不能代替 live 计算；2026-10-06 live 结论由
  `method1_4310_current_audit_20261006.json` 的 schema-3 稳定签名报告提供。final7
  继续作为较早源码快照的历史证据保留。
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
   前保持 fail-closed。门禁只检查候选实际 emitted 的列；Method 2 identity 未解析候选
   单独跳过并进入根级 JSON/TXT 状态报告，其它批量模式计算或 writer 异常仍按候选汇总后
   整批拒绝。磁盘和 ZIP 都在发布前完成本次成功候选与报告的全量渲染。新/空根目录直接发布；已有非空根目录发布到内容寻址的
    `.isodistort-batch-v1-<digest>.ready` 目录，并写逐文件大小/SHA-256 manifest；文件锁与
    进程锁保证重叠批次只能整批成功或整批拒绝。成功的空模式表仍与失败区分。
    IsoVIZ 的 atom type 直接来自不可变物理轨道映射；空列、跨类型列或同元素轨道
   错配都会拒绝，不再按标签或元素回退。

### symmetry-adapted 均匀应变模式与导出

1. 应变使用母胞行晶格 `P`、子群 basis `B` 和对称乘子 `M=I+E`，满足 `M(PP^T)M=B^-1(DD^T)B^-T`，并以 `B M P` 重建女儿胞。工程 Voigt 顺序为 `(E11,E22,E33,2E23,2E13,2E12)`，内积权重为 `diag(1,1,1,1/2,1/2,1/2)`。
2. 规范模式来自 ISO rank-[12] 宏观表与目标 embedding 的 `DISPLAY DIRECTION`；程序再用实际嵌入点群和母胞 metric 独立计算 fixed space，并要求整空间维数与 span 一致。缺失或不一致的证据会中止候选导出，不退回坐标轴基或猜 `GM` 标签。
3. CIF 写 max-component-one 的 `q_raw`，IsoVIZ 写 `q_unit=normfactor*q_raw`，Complete modes 同时写 `q_raw`、`q_unit`、`normfactor`、幅度和两种总和。TOPAS 与官网一致，不生成 strain-mode 精修参数，只写由应用应变得到的固定实际晶胞；参考子胞与 `B P` 不一致时明确失败。
4. 定向应变测试覆盖固定空间、官方 F02 数值契约、三种模式导出与 TOPAS 固定晶胞，
   共 89 passed、1 skipped。较早的 EuAl4 F02 单例 WSL live 为 1/1 pass，metric 重建
   相对残差约 `1.00e-15`，报告为 `output/validation/strain_audit_f02_end_to_end.json`；
   当前源码随后已完成双母相冻结矩阵 24/24 全量本地重跑，见
   `method4_local_validation_20261006_current_v2.json`。两者都不提供跨晶系证据。

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
2. 原 `provided_origin_shift` 只写入 metadata、不参与匹配或 residual；现在会从女儿相分数坐标中显式扣除。网页增加 matching、robust Å 阈值和已知原点输入，网页与 API 仍调用同一 `IsoDistort.search_method_4`。
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
2. 网页实时状态、搜索、选择及导出快照创建共用 `RLock` 与单调 revision；旧标签页的 stale mutation/export 会明确失败。长时导出在独立快照上于锁外运行，并通过单任务锁、进度状态和候选边界取消协调；上传暂存名使用 UUID。
3. Method 3 参数值、后端查询、有限商、点子群和 affine lift 均设有“超限即失败”的预算，不做静默截断。
4. ISODISTORT_VALIDATE 按晶格度量、物种占位和可选原子顺序比较 CIF；ISOVIZ_INPUT 对路径、CSV、幅度和启动器进行显式校验。三个项目共用根 `.venv` 与各自 yaml 配置。
5. 取证确认受保护 `resources/isobyu/` 中的 ignored `iso.log` 来自一次旧的裸 WSL 探针：该命令
   显式 `cd` 到二进制/数据目录后运行 `./iso`，时间和命令流与日志完全一致；生产包装器
   始终在 `~/.id/tmp` 运行并用独立 `ISODATA` 链接读取数据。包装器的 shell 链已由分号
   改为 `export ... && cd <stage> && <binary>`，使 stage 不可用时在启动二进制前失败，
   不会退回调用者 cwd 写日志；无 WSL 命令构造回归已通过。现有受保护日志未被清理或改写。

### ZIP 子群目录短名与防覆盖

1. 导出目录名统一由 `features/export/distortion_formats.py` 生成：Method 1 为 `<IR>_<OPD>_SG<number>`，Method 2 为 `<IR>_<OPD>`；保留 IR 的 `+`/`-` 和 `4D1` 等 OPD token，并替换 Windows 非法字符、控制字符和保留设备名。
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
9. 2026-10-04 对同一桌面用户和同一 Ubuntu 做了三路重复探针：普通 PowerShell → WSL
   `10/10` 成功，实体 `.venv\\Scripts\\python.exe` → WSL `0/10` 且均为
   `Wsl/E_ACCESSDENIED`，`run_cris.ps1` → WSL `10/10` 成功。由此排除“WSL 服务随机
   失效”和“需要管理员权限”，把根因收敛为 OneDrive 内 Python 进程映像的稳定访问限制。
   同一 cwd 的补充 A/B 中，`.venv` 启动器加 `-S` 仍失败，外部基础 Python 加或不加
   `-S` 均成功，进一步排除了 `-S` 和常规 site 初始化本身。
   包装层现会分别解码 Linux UTF-8 输出和 `wsl.exe` 的 UTF-16LE 系统诊断，并在该错误
   出现时给出唯一受支持的 runner 命令；doctor 也提供同一结论。受支持路径的完整
   ISODISTORT doctor 为 15 PASS、0 WARN、0 FAIL，真实包装层得到
   `/home/devoutwang/.id/tmp`；定向单测 9 passed，包装层与真实二进制集成套件为
   38 passed、2 skipped，改动文件 Ruff 通过。localhost proxy 警告来自 Windows/Clash
   自动代理与 WSL NAT，不影响本地 ISOTROPY 二进制执行。DSH 只读复核保存在
   `output/validation/dsh_sessions/20261004-215458-wsl-access-review/`；它支持上述归因和
   不自动提权/重启的结论，但其低完整性沙箱会拒绝所有 WSL distro launch，因此只作
   证据审阅，不能替代普通用户上下文的 10/10 实测。
10. 2026-10-05 恢复旧版 `scripts/main_web.py` 的直接启动契约，同时保留实体 `.venv`：进一步
    A/B 发现，不是所有 OneDrive 内的 Python 映像都会失败；复制到该目录的基础 Python
    仍为 `Wsl/E_ACCESSDENIED`，而基础 Python 的 NTFS 硬链接从同一目录调用 WSL 返回 0，
    并保持 `.venv` 的 `sys.prefix`、site-packages 和 `IsoDistort()` 初始化正常。安装器现只在
    “venv stub 失败且基础 Python 成功”这一可检验条件成立时，用经过 `sys.prefix` 与 WSL
    双重验证的硬链接替换 stub，并保留 `python-venv-launcher.exe` 回退副本；无法证明该条件
    时不修改环境，仍可使用 `run_cris.ps1`。修复后的直接 `.venv\Scripts\python.exe`
    启动烟测中，根页面、`/api/state` 与 `/api/shutdown` 均返回 200，服务正常释放端口；
    启动器/Web 定向回归 6 passed，随后用修复后的直接解释器运行安装器与全部 Web
    测试为 50 passed，直接 doctor 为 15 PASS、0 WARN、0 FAIL，相关 Ruff 检查通过。
    首轮完整回归另暴露提交 `483fe07` 删除 `resources/isobyu/smodes_sample.out` 后测试仍依赖该
    只读目录的既有夹具归属问题；历史样本已迁到 `tests/fixtures/` 并加精确 ignore
    例外，未回写 `resources/isobyu/`。最终完整回归为 677 passed、3 skipped、0 failed、1305
    warnings，用时 1101.99 秒。该修复只改部署、入口与测试夹具归属，不改 `backend`
    科学算法、导出或黄金数据。
11. 2026-10-05 修复 Method 2 参数 k 候选的四格式生成与 ZIP 下载。根因是参数 k 的
    SMODES 数值基虽然完整，但其模式身份按设计保持 `unresolved`，而四个 writer 共用的
    发布门禁只接受经 ISO 微观列和来源绑定验证的模式，因此整批候选在渲染前统一被拒绝。
    现从选定 embedding 的 exact `DISPLAY DIRECTION` 保留每个谐波的 KVALUE、k 坐标和
    不变方向，再用同一 ISO 会话获取 `DISPLAY DISTORTION` 微观列；重复 irrep 依靠真实
    查询表边界归属，中间无位移列的查询不会破坏顺序。只有逐轨道完整子空间、秩、来源
    和身份全部匹配时才用 canonical ISO 基替换 SMODES 基；表边界缺失、重复 irrep 的空
    结果无法唯一归属、列数或子空间不一致时仍保持 `unresolved` 并拒绝导出。
    EuAl4 LD `g=1/6` 实际全批次验证覆盖 22/22 候选，生成 340496 字节 ZIP、22 个目录和
    88 个非空文件，CIF/IsoVIZ/Complete modes/TOPAS 各 22 个，四类语义标记均 0 失败；
    目标回归 147 passed，最终完整 `ISODISTORT/tests` 为 681 passed、2 skipped、
    0 failed、1332 warnings，用时 1107.40 秒。DSH 只读诊断材料保存在
    `output/validation/dsh_sessions/20261005-180108-method2-export-identity/`；会话在输出根因
    证据后由 Codex 中断，结论已由上述真实 WSL/ISO 全量导出与回归独立复核。
12. 2026-10-05 删除 ISODISTORT 终端版交互及其专用内容：移除 `main_terminal.py`、
    `runtime_launcher.py`、终端专用测试、手工批处理 `terminal` 子命令、专用英语文案，
    并同步安装 doctor、项目规则与用户/开发文档。网页、Python API、`run_cris.ps1`
    兼容启动器及 `backend` 科学计算/四格式导出路径均保留。残留门禁确认两个入口文件
    与旧字节码均不存在，`run_batch.py terminal` 以 argparse exit 2 拒绝且只列出
    `cif30` / `external`。定向网页/API/安装器回归为 67 passed；doctor 为 15 PASS、
    0 WARN、0 FAIL；最终完整 `ISODISTORT/tests` 为 675 passed、2 skipped、
    0 failed、1332 warnings，用时 1222.64 秒。相对上一轮 681 passed 少 6 项，均为
    被删除的终端专用用例，不是现存功能回归。DSH 只读诊断材料保存在
    `output/validation/dsh_sessions/20261005-200608-remove-terminal-ui/`。
13. 2026-10-07 修复四部分目录迁移后的路径与导入回归：`irreps_cdml` 从模块位置解析
    `resources/isobyu/`，`core_api` 从模块位置解析 CRIS 仓库根，配置中的输出目录继续落在
    `ISODISTORT/output/` 而不是 `resources/output/`。Method 1–4 包入口改为只导出各自拥有的
    公共 API，并延迟加载具体实现，消除 Method 2/3 模块在全新进程中先导入时出现的循环
    导入；调用方改从符号所属功能包导入，不再经 Method 1 反向重导出。另按 Git 中迁移前
    的同一 README 内容和当前目录结构恢复被误删的连字符、命令参数、公式、URL 与 Markdown
    语法，未凭空重写科学结论。新增架构回归 6 项通过，k 点/空间群表头回归 5 项通过，
    定向运行时回归 6 项通过；221 项相关测试可完整收集，改动范围 Ruff 与
    `git diff --check` 通过。这里没有宣称本轮全测试套件已恢复全绿；完整回归仍以本次最终
    测试记录为准。
14. 2026-10-07 补齐重构后 wheel/sdist 的非 Python 资源契约：`resources/` 增加包标记，
    setuptools 只发现有 `__init__.py` 的常规包，并以显式 `package-data` 白名单收录
    `frontend/web/index.html`、`static/`、`resources/config/*.yaml` 以及只读
    `resources/isobyu/` 的运行文件/数据库；`*.log` 与 `resources/output/` 不会进入分发包。
    新增回归既核对白名单的完整展开结果，也在临时 `site-packages` 中重建等价安装布局，
    从该布局独立导入配置、CDML 表和网页服务并读取全部关键资源；另以 AST 防止重复定义。
    实际 wheel 已通过 pip 的隔离构建（未修改实体 `.venv` 依赖），归档包含配置、网页与全部
    冻结运行资源，未含 `.log` 或 `resources/output/`；独立解包后配置与网页可导入，ISO
    路径解析到解包资源。提供 `isodistort-web` 命令，调用同一网页入口。sdist 本轮未构建。
15. 2026-10-07 从本地 Codex 历史会话证据恢复迁移时丢失的
    `tests/manual/official_html_resolver.py` 与 `method4_provenance.py`，并恢复此前已经完成的
    Method 1–4 校验器实现；随后按新目录结构修正源签名、配置相对路径与测试 fixture。
    不再把缺失内容归因于“测试比实现新”。numeric semantics 校验保留 numpy/非有限值
    JSON、TOPAS 别名链、完整模式来源证明、源码/输入/上下文签名及结束复签；中途中断的
    `pending_final_signature` checkpoint 拒绝复用，帮助文本与手工文档同步说明重启要求。
    全部测试可收集，原 3 个 collection error 已消失。
16. 2026-10-07 补强周期与调制计算：分数坐标边界按周期等价归并，CIF 空间群元数据与首个
    实际可解析结构块绑定，坐标表达式解析要求完整消费；SMODES 子胞映射中有未覆盖原子
    时整列拒绝，不把缺失映射伪装成零位移。谐波阶数在母相原胞平移基上求解精确有理
    同余并以广义 CRT 合并；四分之一相位由有限平移商的 Bézout 构造，字符阶数不超过 2
    时不生成独立 quadrature，波矢文本格式不设晶体学分母上限。正反例包含 13/30 谐波、
    1/49 波矢、1/400 长周期、I 心倒格矢等价及真实零位移原子。新增波矢 Γ 判断与已有
    星臂相位函数使用不同名称，保留既有 fail-closed 星臂证据门禁。
17. 2026-10-07 修复安装布局的写入与执行：运行目录由配置层唯一解析；源码默认仍为项目
    `output/`，wheel 默认落入用户状态目录，支持绝对 `ISODISTORT_RUNTIME_ROOT` 覆盖与
    YAML 显式绝对路径。Linux 每个 wrapper 使用私有短暂存目录；非可执行原文件仅复制
    并对副本设置权限，doctor 使用同一路径准备并区分直接/暂存/失败。只读运行资源未变。
    Windows 上已覆盖复制、缓存失效、失败与 WSL 路径兼容；真实 POSIX shell 集成测试
    在本 Windows 环境跳过，不能据此宣称完整原生 Linux 科学计算验收。
18. 2026-10-07 修复 BUG-006 的静默 cosine fallback：三条母胞单实列路径共用
    `_lift_real_column`，精确相位原语下沉至 `backend/utils/lattice.py` 并由 superspace
    复用；非自共轭字符缺少成对实列/复列来源时拒绝，自共轭模式验证真实母相平移格、
    中心化位点的位移协变性与子胞周期性。Method 4 相同原子数路径也执行同一验证。
    完整子胞列保持可用，无法解释的 `opd_direction` 不再静默忽略。15 项新增回归覆盖
    q=1/4、中心化半 k、精确 sign-only、原点协变与同原子数 API 路径；相关四组测试为
    101 passed、3 skipped。成对实列/复列及一般 k-star/OPD 来源合同仍是待完善能力。
    最终复核另补 6 个浮点边界反例：`±1e-13`、`1±1e-13` 与 `0.5±1e-13`。
    波矢恢复若会把明确非零的十进制值舍入为整数/半整数，则保留该十进制有理数；
    波矢格式化与相位计算使用同一原语，避免误认 Γ 点或自共轭字符。普通 `1/3` 的
    有理数恢复保留，晶格矩阵的既有近似合同不变。四组定向复核为 90 passed、1 skipped。
    后续复核再修复显式 Γ/等价倒格矢提前绕过中心化位点协变检查，以及 `k=None` 在
    缩胞时丢失不周期位移场的问题。未声明 k 的列只要求在请求子胞平移下为 +1 字符，
    不强制视为整个母相 Γ；显式 k 则必须证明母相原始平移字符。新增 27 项回归后，
    fallback 与 distortion 两组为 102 passed、3 skipped。
19. 2026-10-07 修复通用分数基矢构胞的平移格与守恒合同：原胞平移格的精确 coset
    重建及 spglib 有理数恢复下沉到 `backend/utils/lattice.py`；Method 1 仿射群和
    `features/input_cif/coordinate_transform.py` 复用同一算法，不按 HM 字母猜上传坐标系。
    构胞首先验证有限、非奇异基矢，分数基矢须满足 `B @ inverse(P)` 为整数矩阵，
    并在结果边界验证原子数、物种/占位与 `abs(det(B))` 体积比一致。整数路径不再
    `allclose` 后直接截断。非母相平移的 P 半轴反例拒绝，I/F 原胞与合法分数超胞保留。
20. 2026-10-07 补齐 FINDSYM 运行失败时的坐标系门禁：正常载入链会将 A 心非标准输入
    和 FCC 原胞整体转换到 ISO 标准惯用胞；但明确 `status=unavailable` 的失败分支
    保留原坐标供诊断，不能继续消费标准 k/B。共享会话门禁在 Method 1/3 搜索、
    Method 2、路径选择、共享模式计算和子群导出的后端调用/状态修改/发布之前拒绝，
    包括跳过模式计算的 CIF-only 路径。收尾最小反例使用正交单 Ni 原子、非标准原点：
    FINDSYM 失败后旧 CIF-only ZIP 路径按标准 P-1 操作展开，读回由 1 原子变成 2 原子；
    门禁在构胞、规格收集和公开导出入口复用同一判据，要求修复环境后重新载入。
    不按材料或 A/C 字母特判。缺失该元数据的旧私有 fixture 不因此获得 setting 认证。
21. 2026-10-07 修复官网回归的假 skip 与字面基矢误判：EuAl4 Method 1 OPD 页面改由
    不可变 HTML 清单的内容、POST form 和完整母相/作用域上下文唯一定位，不按文件名
    猜身份。123 个完整 IR/OPD/方向/SG/HM/s/i/k-active 身份一一对应；其中 90 行字面
    相同，33 对为不同 conventional basis。逐对精确 GL(3,Z)、完整 Seitz 群相等以及
    体积/群指数不变量全部通过，严格禁止用 parent conjugacy 合并不同 orientation/domain。
    见证保存于 `output/validation/post_restructure_eual4_opd_witness_20261007.json`；
    测试不读取这份实测报告或硬编码 123/33/代表基。另执行既有 `make_cifs_30.py`，
    生成并读回核验 30 个合成空间群样本，原来缺这些输入而跳过的 3 个测试实际通过。
    合成样本不是实验黄金数据，也不表示 30 种晶体的 Method 1–4 全链验收。

## 当前验证矩阵

以下官网/live 项目保留其对应源码快照的证据；除明确标为本轮重构后复验的项目外，
不得将“当前”字样理解为 2026-10-07 新源码已经重新跑完同一长时矩阵。

- 三种母相 Method 1 官网与现有网页版保存产物静态门禁：PASS，两侧各 330/330 候选、各 1320/1320 核心文件；175 项科学差异全部是 4310 旧程序的 125 项模式数差异和 50 项 CIF 语义差异，不是下载错误。4310 N1+ 4D1 当前内容正确，历史覆盖目标为 `inconclusive`。
- EuAl4/NdNiO2 Method 1/2 保存产物：PASS，275/275 候选，1100/1100 核心文件。
- EuAl4/NdNiO2 Method 1/2 live 快照：签名 `e3347be0…5908` 下 PASS，275/275 候选与模式数一致，0 missing/extra/duplicate/mismatch；本轮代码改变后不再称为当前源码结论。
- 4310 Method 1 2026-10-06 快照全量 live：schema-3 报告签名为
  `96ff44cd307432ae2d8ca672afed6ff7c9e7dedbfc39ee07cfd71fefeaf94d12`，运行结束时签名未变化；
  125/125 候选身份配对、模式数、BUSH 覆盖和零振幅 CIF 语义全部通过，0 失败。页头、
  8 个物理轨道、逐物种轨道作用域、完整母相 setting 规范化，以及原始/整体 `z+0.1`
  输入的 N1+ 4D1 `192/192` 另有定向回归。报告 SHA-256 为
  `be7bee90ea47a03bb0aac56debb455a882b0d5d198c593130bd6b11db2f9fa6b`。
- 4310 Method 1 C10/C12 定向 live：报告 schema 为 `method1-4310-c10-c12-live-v2`，
  运行前后源码签名均为
  `ff79ff5e5201c4219379366ecda9d9abc867f9bffcbc20e61e5f7bcdc0bb2cf8`，其中
  `backend/**/*.py` 摘要为
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
- symmetry-adapted 应变：`q_raw/q_unit/normfactor`、CIF/IsoVIZ/Complete modes 共用契约及 TOPAS 固定实际晶胞定向测试 89 passed、1 skipped；当前源码双母相冻结矩阵本地 24/24 通过。跨晶系扩展仍为 `pending`。
- Method 1/2 已记录的签名数值语义抽检：3/3 案例通过；Complete modes、IsoVIZ、TOPAS、CIF 共 12/12 类导出通过。模式数分别为 Eu LD1 `48/48`、Nd Y1 `24/24`、Eu X4− `5/5`，且标签 family、笛卡尔子空间、normfactor、As/Ap/dmax 全部满足声明的判据。
- Method 2 2026-10-06 快照参数-k数值语义：schema-6 live 报告中的 EuAl4 LD1 C1 与
  NdNiO2 Y1 C1 两个 nmod=0 案例均通过，四种导出各 `2/2` pass；报告 SHA-256 为
  `c4830c1895681f19031c6c41113ea5d4039e14e1af2dbad06fe05661edc65cfb`。该结果不外推到
  全部 Method 2 候选、nmod=1、4310 或跨晶系。
- Method 2 参数 k 四格式批量发布的较早签名批次：EuAl4 LD `g=1/6` 22/22 候选成功，ZIP 内
  88/88 文件非空，CIF/IsoVIZ/Complete modes/TOPAS 各 22 个且语义标记 0 失败；LD1 与
  重复 LD5 谐波均使用 exact KVALUE 的 ISO microscopic identity/provenance，歧义分区仍
  fail-closed。该证据只证明当时签名的批次与代码路径；2026-10-06 的最终数值结论仅为上项
  schema-6 两个 nmod=0 代表例，不外推为所有候选、母相和晶系的官网一致性。
- Method 3 官网下载：PASS，40/40 查询，77/77 embedding，308/308 核心文件。
- Method 3 route：61 single-IR + 16 coupled-only，0 不确定。
- Method 3 2026-10-06 默认产品比较：schema-5 报告为 13 exact + 27 affine-equivalent +
  0 differences；77/77 embedding 匹配，逐组候选数一致，0 错误/跳过，source signature
  稳定；报告 SHA-256 为
  `796ae2f39235c6d61422cb9d02403da11e795a60815f6a0590ea02b45c6c45f7`。
- Method 3 Stage A：40/40，77/77 官网 embedding 覆盖。
- Method 3 Stage B：77/77 官网 embedding 可达；3/3 新轨道不可达。
- Method 3 等价 basis/centering 变形：PASS，40 个案例、80/80 个变体，0 missing/extra。
- Method 4 2026-10-06 快照本地闭环：24/24 通过，schema-3 报告 SHA-256 为
  `e19a005dd77bad4b12c7dcc2600889c8e89fe3a689d572c542e6f626a954bb5d`；官网 24 个案例
  均使用正确冻结输入，schema-5 审计为 23 个完全通过、EuAl4 G05 auto-origin 1 个
  证据警告、0 个 inconclusive/fail，报告 SHA-256 为
  `946dba00162ca813060153c4f38daf3131030ffe352525bcc439a2062c5f9519`。F01 物种拒绝、
  F02 均匀应变成功与 F03 robust 距离阈值拒绝路径均通过。
- 网页/API 对齐：Method 3 均提供 `route_resolution` 与完整 known routes；Method 4 均提供 `As/Ap/raw/normfactor`、residual 与六分量均匀应变。网页 ZIP 显式传递当前 nmod，不再由可变会话状态隐式决定参数 k 模式。网页女儿相上传文件在每次 Method 4 成功或失败后删除，当前母相上传文件在替换或服务退出时删除。
- 实体 `.venv` 迁移后的旧启动链证据仍成立：标准 venv stub 以及它启动的普通子进程和
  ShellExecute 进程链无法在 Python 内自愈；`run_cris.ps1` 曾是唯一受支持入口。该结论现
  被上文 2026-10-05 的“安装阶段硬链接修复”收窄：修复后的直接
  `.venv\Scripts\python.exe ISODISTORT\scripts\main_web.py` 已完成根页面、初始化 API 和停止 API
  烟测；`run_cris.ps1` 保留为无法建立硬链接时的兼容入口。
- 此前生产源码快照的全量回归：654 passed、3 skipped、0 failed，用时 1003.51 秒；
  全套在真实 WSL/ISO 可用的非 sandbox 环境执行，skip 为可选能力门禁。该回归覆盖
  microscopic 查询上下文/完整域延拓、生产 CIF 位移合同、四 writer、批量原子发布及既有
  功能，但早于本轮 proof 绑定加固，不能冒充当前最终源码全量结果。
- 2026-10-06 完整 `ISODISTORT/tests` 为 `884 passed, 2 skipped, 0 failed`；此前
  669/675/681 等完整回归计数同样只代表对应历史源码快照，不代表重构后当前源码。
- DSH 只读诊断会话
  `output/validation/dsh_sessions/20261004-122244-method1-identity-review-stdin/` 保存完整 prompt、
  metadata、stdout/stderr；Session ID 为 `session-d129b86d-2b15-43da-9b12-dd5a133966e1`。
  DeepSeek 的 72 passed、2 skipped 及候选原因只作诊断材料；BUSH 排序疑点经矩阵方向和反例
  复核后排除，邻近的 proof 重绑定缺口由 Codex 独立复现、修复并以上述目标回归验证。
- ISODISTORT_VALIDATE：19 passed。
- ISOVIZ_INPUT：20 passed。
- 统一部署器本轮回归：10 passed；ISODISTORT 0.4.0 wheel 已重新构建并独立解包烟测，
  另两个项目的 wheel 为历史构建证据。本轮全项目 `doctor --dev` 为 24 pass、
  2 个缺少用户输入的 warn、0 fail。

2026-10-06 源码快照的机器报告为 `method1_4310_current_audit_20261006.json`、
`method2_numeric_semantic_live_current.json`、
`method3_official_local_comparison_20261006_current.json`、
`method4_local_validation_20261006_current_v2.json` 和
`method4_official_audit_20261006_current_v2.json`。历史关键报告包括
`method3_embedding_route_audit.json`、`method3_official_local_comparison.json`、
`method3_basis_metamorphic_audit.json`、`method3_stage_a_diagnostic_audit.json`、
`method3_stage_b_feasibility_audit.json`、`method4_local_validation.json`、
`method4_official_audit.json`、`strain_audit_f02_end_to_end.json`、`live_method12_report.json`、
`audit_method1_official_downloads_20261003.json`、
`audit_method1_overwrite_forensics_20261003.json`、`method1_4310_live_final6_20261004.json`、
`method1_4310_live_final7_20261004.json` 和
`method1_4310_vector_c10_c12_extension_final_20261004/summary.json`；全部位于
`output/validation/` 且不提交 Git。长时命令见 [MANUAL_VALIDATION.md](MANUAL_VALIDATION.md)。

## 2026-10-08 E 盘迁移恢复与环境隔离

- 三个同级目录已恢复为 `E:\CRIS`、`E:\GD` 和
  `E:\Best_Model_Parameters`；原桌面目录不存在。CRIS/GD 均以 Python 3.12.5
  重建实体 `.venv` 并隔离重放迁移前 freeze，两个环境的 `pip check` 均通过；GD 从自身
  环境导入 TensorFlow 2.21.0。全项目 doctor 为 24 pass、2 warn、0 fail；两项 warning
  分别是用户尚未提供 compare CIF 和 Best_Model_Parameters 下尚无振幅 CSV，不是迁移失败。
- 恢复时实际复现了跨环境污染：同一 PowerShell 进程调用 `run_cris.ps1` 后，原实现会遗留
  CRIS 的 `PYTHONPATH`、`VIRTUAL_ENV` 和 `PATH`，使 GD 安装器误把 CRIS 包判为已安装，
  初次恢复的 GD 环境因此出现 22 项依赖破损且 TensorFlow 缺 `typing_extensions`。
  `run_cris.ps1` 现用 `try/finally` 恢复调用者环境；GD `main_requirement.py` 启动子进程前
  删除 `PYTHONHOME`、`PYTHONPATH`、`VIRTUAL_ENV` 和 `__PYVENV_LAUNCHER__`。Windows
  启动器动态回归及 GD 子进程回归均覆盖该根因。
- GD 转换器仍引用重构前已删除的 `isocore.api` 与 `isocore.io`。现改为从
  `backend.api` 和 `features.export.distortion_formats` 载入唯一生产实现，并以模块装载回归
  固定该契约。实际重生成 EuAl4 LD1 C1 成功：母相/子相空间群为 139/99，48 个唯一模式、
  60 个原子，生成和安装的 `LD1_C1_alris_functions.py` SHA-256 相同，零振幅 TensorFlow
  坐标为有限的 `(60, 3)` 数组。48 个正振幅界限满足
  `bound × dmax = 1 Å`，最大绝对残差 `1.46e-7`；metadata 的 CIF 为
  `E:\CRIS\experiment_data\EuAl4 Parent.cif`，并已纠正 `maxamp` 语义说明。该结果只验证
  当前 EuAl4 LD1 C1 生成合同和数值不变量，不外推为全部材料、晶系或完整训练已验证。
- 排除 `.git`、`.venv`、历史迁移记录和测试输出后，对三个 E 盘工作区扫描旧桌面三根路径
  为零命中；受保护实验数据、`resources/isobyu`、天仁原始 notebook 和外部
  `D:\OneDrive\...` 科学数据均未修改。
- 迁移后完整回归：`ISODISTORT/tests` 为 `1027 passed, 2 skipped, 0 failed`
  （1269.52 秒）；根部署器 11 passed；ISOVIZ_INPUT 22 passed；
  ISODISTORT_VALIDATE 19 passed；GD tests 在 GD 环境和 CRIS 环境各 9 passed；GD
  TensorFlow 运行时烟测通过。完整套件的 2 项 skip 为既有可选能力门禁；第三方弃用和 CIF
  宽容解析 warning 保留为 warning，未计作额外科学证据。

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
- 4310 Method 1 已有官网静态全集、对应历史源码签名的定向根因回归和 schema-3 125/125
  候选全量 live；Method 2–4 尚无该母相正式矩阵。
- 官网 Method 1–3 的默认子群导出为 `As=0`，因此这些路径的非零振幅数值归一化主要由公式、单元测试和本地生成验证；Method 4 已有非零官网对照，nmod=1 尚无单独归档的官网四格式非零参考对。
- 本轮已识别 Java 与 IsoVIZ 启动器，但按自动验证边界没有启动 GUI；20 个静态测试只能
  证明 `.isoviz` 生成、路径和启动参数逻辑，不能证明窗口实际打开或 IsoVIZ 正确认出数据。
- Method 4 官网 EuAl4 与 NdNiO2 的 G01–G06、P01–P03、F01–F03 已验证零保持、正负幅值、原子重排/周期换像、完整 Γ/参数-k 等价多维模式基、双模式、显式/自动原点、参数-k origin phase、超胞零/非零分解、P1 噪声分配、应用均匀应变张量和拒绝路径。EuAl4 G05 auto-origin 仍仅缺 basis HTML；清晰填写截图、accepted result identity、完整导出和显式-origin 数值对照已足以维持 `pass_with_warnings`。symmetry-adapted strain 的生成、标签和 CIF/IsoVIZ/Complete modes 导出以及 TOPAS 固定实际晶胞现已实现，当前源码已重跑双母相 24/24；官网按子群算子自动枚举全部允许原点、4310/跨晶系矩阵以及 occupancy/magnetic/rotational 分解仍未验收。

## 主要科学依据

- [ISODISTORT Method 3 help](https://iso.byu.edu/isodistorthelp.php)
- [ISOTROPY user documentation](https://iso.byu.edu/isotropy_doc.php)
- Campbell, Stokes, Tanner & Hatch, *J. Appl. Cryst.* 39 (2006), [DOI 10.1107/S0021889806014075](https://doi.org/10.1107/S0021889806014075)
- Stokes & Campbell, *Acta Cryst.* A73 (2017), Appendix B, DOI `10.1107/S2053273316017629`
- Hatch & Stokes, *Phys. Rev. B* 65 (2002), DOI `10.1103/PhysRevB.65.014113`
- Stokes, van Orden & Campbell, *J. Appl. Cryst.* 49 (2016), DOI `10.1107/S160057671601311X`
- Wagner & Schönleber, *Acta Cryst.* B65 (2009), [superspace introduction](https://journals.iucr.org/b/issues/2009/03/00/bk5084/index.html)：有理调制、公度周期、cosine/sine 分量及相位约定。
- [IUCr CIF specifications](https://www.iucr.org/resources/cif/spec)：data block、loop 与结构数据的语法边界。
- [International Tables, Vol. C, §1.1](https://it.iucr.org/Cb/ch1o1v0001/)：原/惯用胞、中心化平移与倒格矢的定义及体积关系。
- [IUCr Superstructure definition](https://dictionary.iucr.org/Superstructure)：子胞平移群必须为母相平移群的子群。

这些文献定义算法与不变量；官网结果只用于检验实现是否复现同一科学语义。
