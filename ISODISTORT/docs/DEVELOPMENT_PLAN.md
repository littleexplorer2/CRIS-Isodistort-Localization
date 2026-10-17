# ISODISTORT 当前开发与验证计划

本文件只维护当前待办与后续计划。完整历史原文已无损压缩，恢复方式见
[历史索引](history/README.md)。下载步骤见 [DOWNLOAD_CHECKLIST.md](DOWNLOAD_CHECKLIST.md)，
验证命令见 [MANUAL_VALIDATION.md](MANUAL_VALIDATION.md)，证据入口见
[BUGFIX_VALIDATION_REPORT.md](BUGFIX_VALIDATION_REPORT.md)。
2026-10-16 的重置决定、作废范围与归档指纹见
[重置记录](manifests/validation_reset_20261016.json)：旧的 74 项冻结清单与其批次产物
已作废归档，当前**没有**生效的本地刷新冻结清单。

## 当前执行顺序

**当前范围（2026-10-17）**：本轮执行输入网格的理论修复及回归；官网新差分仍等下载完成。
旧 Method 3 闭包因答案污染及晶格判据缺陷而失效，不能仅凭源码稳定或清空缓存恢复。

1. **A 阶段（可立即执行，用户在官网操作）**
   - Method 2 共 **98 个查询**已冻结：16 个固定特殊 k 标签 × nmod{0,1} = 32 个，33 个参数化
     标签 × nmod{0,1} = 66 个（代表值 `a=1/6`、`b=1/4`、`g=1/3`，仿射模板一般位置精确验证 pass；不代表非公度模式验证）。目录已建好。
   - 官网补页**已完成并核对为 `pass`**：EuAl4 Method 4 G05 auto-origin 的 basis 页
     （SHA-256 `9e2bfc51…`）＋提交态截图 `Screenshot 2026-10-17 005445.png`（`230d2182…`）
     共同证明基矢 `(0,1,1),(1,0,0),(0,0,-1)`、`Automatic origin detection`、`nearest-site`、
     robust 未启用；被替换的旧截图已在记录中登记。操作与验收记录见
     [下载清单第 4 节](DOWNLOAD_CHECKLIST.md)。
   - 新官网材料写入
     `ISODISTORT/output/validation/official_full_download/20261016-official-full-v1/`；
     `output_compare/` 与 `webpage_info/` 保持只读。
2. **B 阶段（独立重建范围）**：Method 3/4 暂不下载。旧闭包工具已停用，历史案例分类匹配
   不能证明覆盖。须从用户输入和通用晶格数据定义有限原胞指数/基矢界限，明确母/子心化、
   setting 与等价关系，并证明枚举范围。Method 4 另行冻结独立 daughter 输入及变换上界；
   不继承旧闭包种子、checkpoint 或推算规模。
3. 下载完成后运行只读审计：逐页面/候选登记身份、四格式、非空性与 SHA-256；官网确实未提供的
   项显式记 `not_available`，不得补造。审计未通过前不进入算法差分。
4. 按预先登记的独立用户输入重建本地刷新清单（新 batch id + 新 runtime revision），不得
   从官网候选结果反推输入；本地独立冻结后再作集合对照。旧 74 项清单只作历史参考。
5. 用新官网集重跑 Method 1 → 2 → 3 → 4 本地对照；旧报告的数字与结论只适用于旧官网集，
   不得作为新范围的验收。
6. 遇到实际 blocker，保留原始错误，修通用根因，再续跑相同案例；不得删例、换输入或降低身份
   门槛来获得通过。
7. 最后执行跨晶系、组成、工具兼容、性能和并发扩展验收。

## 待修的程序漏洞

P0 为交付阻断，P1 为结果或可复核性风险，P2 为性能与可诊断性问题。
已有实现但尚缺全矩阵证据的项目，必须完成验收后才关闭。

### [P0] BUG-001 参数 k 批量导出：全量验收待完成

- 输入：NdNiO2 Y `a=1/3`，strain + 全物种 displacive，`nmod=0/1`。
- 待做：通过冻结案例 `M2-ND-Y-1_3-NMOD0/1` 验收各轮 48 个候选。
- 判据：每个候选生成四格式，或明确记录模式身份未解析的原因；一个未解析候选
  不得阻止其他合法候选发布。明确跳过可证明错误隔离有效，但不能称为全文件下载完成。
- 失败必须记录物理轨道、具体判据、异常原文与
  `(parent SG,k,KVALUE,IR,direction selector)`。
- 定位入口：`tests/manual/diagnose_method2_export_preparation.py`。
  既有证据：`output/validation/ndnio2_export_audit_20261006/`。
- 通用修复已实现，真实长时验收待完成；若重跑成功，关闭验收即可。

### [P1] BUG-002 失败原因传播：随 BUG-001 验收

逐候选确认 `ModeIdentity.reason`、网页警告、导出状态与机器报告使用同一具体失败原因
和查询上下文。验收必须能区分 ISO 查询失败、方向解析歧义与空间证明拒绝。

### [P1] BUG-003 不可交付候选：随 BUG-001 验收

核对未计算时的 `not_checked` 与已判定 `unresolved` 的区别；后者必须显式不可交付，
网页 `selectable=false`、API `export_ready=false` 与批量状态报告使用同一门槛。
不要为标记状态重复实现身份算法。

### [P1] BUG-006 非 Γ 模式来源合同扩展

静默 cosine fallback 已移除。剩余工作是建立经证明的复列/成对实列、k-star 臂、
相位原点与 OPD 系数合同，或提供完整子胞列；生成、分解与四 writer 共用。
缺少来源字段继续 fail-closed。验收覆盖 `q=1/4` 的正交 cosine/sine、自共轭区边界、
中心化点阵和原点平移协变性。

### [P1] THEORY-001 完整 setting/Wyckoff 元数据审定

依据理论资料库与标准表，建立 Hall setting、原点、坐标和 Wyckoff 字母的完整对应，
覆盖一般位置、特殊参数退化和相同重数的不同位置。writer 的部分表和 spglib 标准化
结果须用实际仿射变换校验，不能只凭空间群号或重数赋字母；未知字段保持未知。
完整认证前，不把当前有限回归外推至全部 230 空间群。

### [P1] THEORY-002 缺失物理表示与数值边界

occupational 需建立标量站点表示、母相占位基线、成分约束及物理可行范围；rotational
需轴矢量与刚性单元约束；magnetic 需时间奇轴矢量、磁群及 setting 合同。表示完成前
继续明确拒绝正式计算，不能复用位移活性或坐标排序二分模式。
对版本特定应变规范、极端病态晶格的容差、非公度谐波空间封闭和任意 k 搜索完备性
补充独立证据；理论义务与源码对应见参考资料库的实现审查清单。

### [P2] BUG-004 大批量导出成本

逐候选串行重算和大型超胞精确有理运算仍需优化。先记录候选级耗时与峰值资源，
定位 `_exact_child_operations`、`_exact_displacive_child` 和导出准备的成本，
再验证可安全复用的缓存。实测进度、取消后保留产物、状态报告与会话锁隔离；
性能优化不得改变身份或 writer 门槛。

### [P2] BUG-005 大型超胞证明覆盖不足（待证实）

对高指数超胞逐候选记录物理轨道分组、reference/candidate 秩、条件数、principal angle、
稀疏域覆盖和失败原因。与低指数正例比较后才判根因；现有规模相关性不是根因证明。

## 待执行的优化计划

### P0：官网已冻结有限输入集下载 → 本地清单重建 → 本地刷新

1. **官网下载（当前 P0，由用户执行）**：按 [下载清单](DOWNLOAD_CHECKLIST.md) 中三母相
   Method 2 已冻结查询逐项保存，Method 3/4 等独立范围冻结后再执行；不重下已登记的旧官网件，不动 `output_compare/` 只读归档。
2. **下载后只读审计**：新增支持新查询目录的 Method 2 只读审计，为 Method 4 审计补充
   `--official-root`；Method 3 使用既有审计入口但需显式匹配新输入清单。具体缺口见下载清单第 9 节；
   候选数、四格式、页面角色、身份与 SHA-256 逐项核对，官网确未提供者记 `not_available`。
3. **清单重建**：从独立用户输入生成新的 `local_web_output_refresh_manifest.json`
   （新 batch id）与新的 runtime revision；旧清单保持 `superseded` 状态可查。
4. **本地刷新**：按 Method 1 → 2 → 3 → 4 经真实本地网页执行；Method 1–3 保存完整
   TXT/CSV、四格式原始 ZIP 及实际产生的状态报告，Method 4 保存 TXT/CSV 或预期拒绝的
   原始响应。查询、候选和文件分别计数，完成状态绑定最终运行签名，未完成案例可追溯。

### Method 2/3/4 扩展矩阵（三母相通用）

1. Method 2：按已冻结代表输入与 `nmod=0/1` 下载；变更参数后重新运行精确网格审计；
   核对 site-mode 身份、维数、归一化、子空间与四格式。已归档的 X 点等官网件只登记，
   不重下，也不能决定其它 k 点的预期答案。
2. Method 3：从独立输入枚举并证明有限晶格/setting 范围；四类历史标签只作描述，
   核对 `(SG,basis,origin,s,i)`、等价子格、single/coupled route。先审计官网结果表和
   候选全集，再作算法差分。
3. Method 4：从核定路径生成并冻结 daughter CIF、哈希、basis/origin/matching 和候选
   身份，覆盖零、正负单模式、正交双模式、参数 k、多模式/近上界、重排/换像/原点/
   小噪声、均匀应变、物种不符和距离阈值拒绝。核对幅度、As/Ap/normfactor、
   residual、原子匹配、应变及实际可用导出。

### P1：泛用性与能力缺口

- 扩展精确微观身份链到更多晶系、点阵、超胞、混合/部分占位，逐层保留正反例。
- Method 3：arbitrary/multi-parameter k、point-group-only affine、reciprocal 子格、
  唯一 primary coupled IR/OPD 分解，以及 rotational/occupational/magnetic
  同等级表示尚未完成；这些物理类型的正式计算须维持明确拒绝，不能以位移活性替代。
- 位移来源合同暂限单物种满占位；补齐非 Γ transporter 的 component-to-star-arm
  相位证据，扩展轴矢量、磁和占位模式身份链。
- Method 4：核验全部允许原点自动枚举，并扩展 occupancy/magnetic/rotational 分解。
- nmod=1 非零四格式官网参考对仍是证据缺口：新下载清单把每个 k 标签的 `nmod=1` 收进范围，
  旧归档不重下。
- 如需抵抗外部伪造 provenance，须补可重放 transcript/矩阵或签名；现有进程内合同
  不能当作密码学认证。

### P2：运行治理与扩展验收

- 参数 k 生成缓存：统一可写位置、版本/来源、进度与失败记录；不得写入只读
  `resources/isobyu/`，运行前核验 WSL staging/符号链接去向。
- 冻结可再分发的七晶系、P/A/B/C/I/F/R 点阵、一般/特殊 Wyckoff、单/多元素、
  混合占位、不同 ASU 大小的样本；来源见 [EXTERNAL_CIF_SOURCES.md](EXTERNAL_CIF_SOURCES.md)。
- 后续单独建立限速的官网版本漂移检查；官网不可用、证据不足或身份歧义记 inconclusive。
  只下载当轮冻结清单内的项；已登记归档只补已定义的缺页，不做全量重下。
- 扩展 TOPAS 本体/解析测试和 VESTA、IsoVIZ GUI 分层抽检；静态解析不代表工具正确识别。
- 三次冷/热重复，记录 P50/P95、峰值内存、WSL 进程数、ZIP 大小与规范化摘要。
- 多标签、多进程、快速提交、ZIP/GenDB 中断与缓存并发；核验服务、端口及子进程回收。
- 扩展结果产生前冻结发布阈值；三个四方母相的通过不能外推到任意晶体。
