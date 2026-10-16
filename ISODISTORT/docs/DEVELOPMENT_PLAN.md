# ISODISTORT 当前开发与验证计划

本文件只维护当前待办与后续计划。完整历史原文已无损压缩，恢复方式见
[历史索引](history/README.md)。下载步骤见 [DOWNLOAD_CHECKLIST.md](DOWNLOAD_CHECKLIST.md)，
验证命令见 [MANUAL_VALIDATION.md](MANUAL_VALIDATION.md)，证据入口见
[BUGFIX_VALIDATION_REPORT.md](BUGFIX_VALIDATION_REPORT.md)。

## 当前执行顺序

1. 使用已冻结的 [本地重下载清单](manifests/local_web_output_refresh_manifest.json)，
   开始前复核来源 manifest、母相与 daughter 哈希并记录实际工作树、配置、依赖和后端签名。
   当前批次 74 项（Method 1/2/3/4 分别 3/7/40/24），下载尚未启动。
2. 经真实本地网页按 Method 1 → 2 → 3 → 4 执行；每项保存完整输入、输出、状态、
   候选身份与文件清单。旧 API/live 报告不能代替本轮浏览器下载。
3. 遇到实际 blocker，保留原始错误，修通用根因，再续跑相同案例；对修改影响的已完成
   案例重新验收。不得删除失败项、换输入或降低身份门槛来获得通过。
4. 仅尝试补 EuAl4 Method 4 G05 auto-origin 的官网 basis HTML；
   无法取得则保留有依据的 not_available/pass_with_warnings。
5. 依次建立并执行 4310 Method 2 扩展、Method 3、Method 4 正式矩阵。
6. 执行跨晶系、组成、工具兼容、性能和并发扩展验收。

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

### P0：本地网页版输出刷新

按统一 manifest 固定的 74 项执行，不穷举任意输入或筛选组合。
Method 1–3 保存完整 TXT/CSV、四格式原始 ZIP及实际产生的状态报告；
Method 4 保存 TXT/CSV 或预期拒绝的原始响应。查询、候选和文件分别计数。
完成状态必须绑定最终运行签名；所有未完成案例保持可追溯。
官网既有核心归档只读，当前只补 G05 basis 页。操作详见下载清单。

### 4310 新增矩阵

1. Method 2：先冻结 Γ、非 Γ 特殊 k、公度参数 k 的输入清单；参数 k 覆盖
   `nmod=0/1`，核对 site-mode 身份、维数、归一化、子空间与四格式。
   已有 X 查询在当前批次；扩展批次不重复官网下载它。
2. Method 3：Method 2 闭合后，冻结 identity、oriented-cell、special-supercell、
   parameter-k 查询，核对 `(SG,basis,origin,s,i)`、等价子格、single/coupled route。
   先审计官网结果表和候选全集，再作算法差分。
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
- nmod=1 非零四格式官网参考对仍是证据缺口，后续另定义输入；当前不据此要求重下旧官网。
- 如需抵抗外部伪造 provenance，须补可重放 transcript/矩阵或签名；现有进程内合同
  不能当作密码学认证。

### P2：运行治理与扩展验收

- 参数 k 生成缓存：统一可写位置、版本/来源、进度与失败记录；不得写入只读
  `resources/isobyu/`，运行前核验 WSL staging/符号链接去向。
- 冻结可再分发的七晶系、P/A/B/C/I/F/R 点阵、一般/特殊 Wyckoff、单/多元素、
  混合占位、不同 ASU 大小的样本；来源见 [EXTERNAL_CIF_SOURCES.md](EXTERNAL_CIF_SOURCES.md)。
- 后续单独建立限速的官网版本漂移检查；官网不可用、证据不足或身份歧义记 inconclusive。
  当前批次仍只补已定义缺页。
- 扩展 TOPAS 本体/解析测试和 VESTA、IsoVIZ GUI 分层抽检；静态解析不代表工具正确识别。
- 三次冷/热重复，记录 P50/P95、峰值内存、WSL 进程数、ZIP 大小与规范化摘要。
- 多标签、多进程、快速提交、ZIP/GenDB 中断与缓存并发；核验服务、端口及子进程回收。
- 扩展结果产生前冻结发布阈值；三个四方母相的通过不能外推到任意晶体。
