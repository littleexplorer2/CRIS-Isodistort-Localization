# 当前验证证据与缺口

## 2026-10-20 安装入口、启动器与版本验证

仅修复上轮只读审查中的第 2、3 项：`ISODISTORT/main_requirement.py` 根据自身位置找到
CRIS 根部署器，不依赖当前工作目录；ISOVIZ_INPUT 的 YAML 默认目录、缺省回退、运行提示、
部署提示和 Git 忽略规则统一到 `CRIS/launchers/`。`ISOVIZ` / `ISOVIZ_JAR` 的优先级保持不变。
根 `.gitkeep` 保留用户指定位置，旧安装路径和启动器路径不恢复兼容分支。

科学影响判断：改动只涉及部署器寻址、外部可视化启动器发现和版本控制排除；没有修改
CSV/模式输入选择合同、`amp` 的 Best Model Parameter 约定、晶胞、单位、符号、归一化或容差。
无需新晶体学推导；验证不变量为工作目录无关、仅默认搜索根启动器目录、显式环境变量优先、
本机启动器不进入提交。没有用官网数据作为本地计算答案。

最小反例先执行：异地目录运行安装入口报 `ModuleNotFoundError: setup_cris`，配置及默认回退
仍发现旧目录，登记为 4 failed、1 passed；修后同范围及启动器分发测试 9 passed。
首次沙箱测试因默认临时目录权限出现 4 个 setup error，改用本轮专属 basetemp 后复现上述真实缺陷，
未改断言。沙箱部署诊断中的 Java 探测超时；主机只读诊断随后 8 pass、1 warn、0 fail，
自动发现 `E:/CRIS/launchers/ISOViz.lnk`，未启动 GUI。warn 为缺少振幅 CSV。

本轮证据目录为 `ISODISTORT/output/validation/launcher_release_20261019/`：

- `repo.log`：59 passed，13.01 秒，含架构、任意 cwd 安装入口、ISODISTORT 实际 wheel 门禁。
- `isoviz.log`：42 passed，3.01 秒，含配置、环境覆盖、旧目录拒绝及 ISOVIZ_INPUT 实际 wheel 门禁。
- `targeted.log`：9 passed，0.70 秒，与上述结果重叠，不累加通过数量。
- 首次联合回归 `full.log`：1 failed、1392 passed、4 skipped，1186.97 秒。
  GD 安装器测试使用通用模块名，因 ISODISTORT 根目录新增同名入口而导入错误模块；
  现按 GD 自身文件路径加载被测安装器，不修改 `_run` 环境隔离断言，也不修改生产算法。
  修后架构与 GD 联合定向检查 `collision-after.log` 为 53 passed，14.54 秒。
  `collision-before.log` 仅作缩小诊断：该进程结束前文件已修订，回溯显示的源码行可能不同；
  失败依据以源码稳定的首次完整日志为准。
- 最终联合回归 `full-final.log` / `full-final.exit`：**1396 passed、4 skipped、2026 warnings，
  1192.69 秒，退出码 0**。跳过项为 LD1 g=1/6 本地子群数据库未生成、X4− P3 官网/本地 CIF
  配对缺失、原生 POSIX 专用用例和已撤销的 Method 3 输入冻结。前一项不是已下载官网文件漏读；
  本轮未启动长时 GenDB。警告来自依赖弃用及 CIF/元素解析，不作为科学验收通过依据。
- 完整回归后，仅清除 6 个文件的行尾空白或末尾多余空行；`format-equivalence.json`
  逐文件证明含位置的 AST 与序列化编译结果完全相同。原始字节保留在 `format-before/`，
  发布源码指纹另存 `publication-source.json`，不声称其原始字节等于完整回归快照。
  清理后仓库门禁、workflows、两项审计工具及导出/快照定向检查
  `post-format.log` 为 **110 passed、3 warnings，17.11 秒**，与其他测试覆盖重叠。
- 新改动 Python 文件 `ruff check --select F,E9` 通过；尚未重新安装环境，也未做真实 GUI/训练验收。

完整联合回归命令：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest ISODISTORT/tests E:/GD/tests_dev -q -ra -x --tb=short -o faulthandler_timeout=120 --basetemp=ISODISTORT/output/validation/launcher_release_20261019/full-final-temp
```

`download-checklist-audit.json` 将人工清单与冻结查询逐项对应：EuAl4 30、NdNiO2 42、4310 30，
合计 102 项，无缺项或重复，参数与目的路径一致。99 个页面集未齐、2 个部分齐、1 个齐；
已知缺 5 张结果列表。清单增加执行顺序和范围说明，没有扩大查询输入或宣布下载完成。
未知 IR/候选需以新列表登记，Method 3/4 仍待独立输入；这不是所有方法的最终逐文件完整清单。

第 1 项前端迟到响应、第 4 项动作/事务策略双登记按用户要求暂缓；WSL 混合编码诊断、
独立科学验收和真实浏览器/GD 接收的既有缺口仍见 [开发计划](DEVELOPMENT_PLAN.md)。
本轮未改前端状态接收和动作登记逻辑，相关两文件与上轮交付快照逐字节一致。
`final-verification.json` 在格式清理前确认本轮完整回归的 236 个源码/配置/测试文件和 7 个依赖
无漂移；5085 个受保护文件路径、字节数和 SHA-256 均保持一致，无新增或删除。
14 份修前归档逐条解压核验及 ZIP CRC 通过。这是文件完整性核查，不是操作系统逐次读取追踪。

冻结记录后继第 4 版仅更新 4 份支撑文档，35 个文件指纹逐项核对；下载范围、独立输入和
科学验收状态不变，旧运行签名不可复用。发布前的索引范围、源码与冻结指纹检查见
`publication-check.json`；实际提交及远程分支核对结果见本次交付消息和本机 `delivery.json`。
清理前核对绝对路径和 reparse point，仅移除本轮 9 个临时测试目录，计 5340 文件、
972,132,674 字节；清理后剩余 0 个目标目录。日志、指纹、检查脚本及格式清理前字节均保留，
见 `cleanup-plan.json` / `cleanup-result.json`。

历史文件恢复见 [归档索引](history/README.md)。本机 765,015,040 字节的旧派生产物包保留并排除
提交；版本推送不包含官网原件、实验输入或本机启动器。前后版本文件指纹与依赖固定在
`run-environment.json` / `run-environment-final.json`；两者仅有 GD 安装器测试加载路径的差异，
不会把工程检查冒充新的逐候选科学验收。

## 2026-10-19 测试整合与使用说明更新

### 修改范围与覆盖保留

生产层级维持 `frontend → application → backend`。原 52 个 pytest 文件整合为
14 个生产模块检查文件和 18 个根目录工具/工程检查文件；每个生产模块只有一个
`test_<模块>.py`，按所属层存放。26 个手工工具（25 个 Python、1 个 JavaScript）
移到 `tests/` 根，删除旧 `tests/manual/`。8 个根目录 `support_*.py` 保存共享夹具和路径，
已有 `tests/validation/` 仍是唯一比较实现包，不复制其算法。
模块与一键命令见 [测试入口](../tests/README.md)。

迁移是函数体、测试和调用路径的真实合并，没有动态加载旧测试的转发层。
`nodeid_mapping.json` 对照迁移前收集结果，原 **1360 个参数化节点**一一对应；
原 **954 个测试定义、3350 条 assert**保留。旧测试增加一条非空工具清单断言，防止空遍历误通过；
另增加诊断脚本 `--help` 不初始化 WSL 的回归，因此最终收集 **1361 个节点**。
`definition_ast_preservation.json` / `structural_preservation_audit.json` 核查 1184 条定义/赋值，
在明确登记的路径、名称、导入及夹具调整后没有 AST 差异；5 份重复助手改用 3 个既有同体实现。
原 5 个 fixture 定义保留，新增应用层隔离 fixture；未出现顶层重名覆盖或测试助手循环导入。
静态对应证明覆盖保留，不替代实际回归或科学正确性证明。

迁移同时修复了根路径层数、审计源码签名中的多行 `tests/manual` 路径和子进程环境继承；
保留来源漂移拒绝反例，没有放宽签名核验。诊断工具原先导入已删除的 `_parent_map`，
连 `--help` 也失败，现调用公开的 `backend.input_cif.map_supercell_sites`。
模块归属按实际被测合同复核，输入上下文属于 session，k 格过滤属于 Method 2，
应变度量与模式数值语义属于 Method 4。
根架构门禁新增模块检查文件唯一性，并禁止比较实现经任意其他 `tests.*` 工具或助手间接启动计算。

根和子项目 README 已重写部署与使用步骤；`main_requirement.py`、根 `.gitkeep`、
`launchers/` 的位置记录与实际文件一致。按用户范围，迁移后的入口定位、启动器默认配置
和忽略规则只记录现状及可用替代命令，没有修改其运行代码。

### 科学影响、依据与适用条件

本轮没有修改生产科学算法、运行默认值或数值容差；73 个已登记生产源码/配置与上轮最终
快照逐字节一致，范围见 `production-unchanged.json`。测试迁移不改变数学对象、单位、
坐标、模式符号或归一化约定；原数学正反例保留，理论入口仍是
[架构数学合同](reference_library/notes/ARCHITECTURE_CONTRACTS_20261019.md) 及能力映射。
未增加物理近似、材料分支或预计算答案。守恒关系与极限检验没有因文件组织而产生新适用项。

诊断工具的站点映射依赖更新须单独核查：未畸变子胞采用行基 `A_child=B A_parent`，
同一 Cartesian 点满足 `x_parent=x_child B`；对应母位点需满足
`||(x_parent-x_i-n) A_parent|| ≤ ε`，其中 `n∈Z³`，元素/占据相同，且容差内只能有一个像。
`ε` 使用既有 `defaults.symmetry_cartesian_tolerance_angstrom` 配置（当前 0.001 Å），不使用分数坐标欧氏距离。
实现复用 `input_cif/coordinate_transform.py`，未引入新的映射公式；现有映射回归包含斜晶格、
歧义周期像及缺失匹配拒绝。新 `--help` 测试仅证明依赖可导入且不启动 WSL；
诊断工具完整 `--origin-scores` 路线本轮未进行端到端科学验收，不能从 CLI 成功推导数值正确。

来源签名、输入固定、独立计算后冻结再比较的顺序保持不变。官网结果未作为本地求解输入。
`coverage.json` 和当前实现审计只更新入口位置；历史 `code_context.json` 与复用记录 generator
仍保留原快照路径，不能将历史报告改名为新源码的验收证据。

### 验证与可恢复证据

本轮所有日志及核查清单位于 `ISODISTORT/output/validation/test_layout_20261019/`。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest ISODISTORT/tests -q -ra -x --tb=short -o faulthandler_timeout=120 --basetemp=ISODISTORT/output/validation/test_layout_20261019/full-final-temp
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest tests_repo -q -ra --tb=short --basetemp=ISODISTORT/output/validation/test_layout_20261019/repo-complete-temp
```

- 最终完整回归：**1357 passed、4 skipped、2026 warnings，1196.88 秒**，退出码 0，
  见 `full-final.log` / `full-final.exit`。跳过项分别为本地 LD1 g=1/6 子群数据库未生成、
  缺少 X4− P3 官网/本地 CIF 配对、原生 POSIX 专用用例、Method 3 旧输入冻结撤销。
  LD1 跳过不代表已下载官网文件漏读；本轮未自动生成长时子群数据库。
  警告为依赖弃用及 CIF/元素解析提示，未出现未捕获线程异常或 finalizer 异常警告。
- 仓库门禁：**57 passed，11.57 秒**，见 `repo-complete.log`，包含真实 wheel 打包/安装。
  最后模块归属复核后的单独层级门禁 **43 passed，0.57 秒**，见 `layer_gate.log`。
- 定向回归：**414 passed、1 skipped、75 warnings，64.57 秒**，见 `targeted_final.log`；
  跳过项为旧 Method 3 输入冻结撤销，与完整回归重叠，不能累加通过数。
- 25 个 Python 工具编译、19 个 CLI 帮助入口及 JavaScript 语法检查通过。
  全测试与修改的仓库门禁运行 `ruff check --select F,E9` 通过。
- 11 个当前文档的 78 个本地链接、3 个 Python 示例及代码围栏检查通过；
  能力/实现映射的 37 个模块路径、20 个证明标签和 6 个现行反例测试路径均可定位。
  这只检查引用完整性，不是重新审读文献或证明。
- 首次全量因迁移后的 session 子进程路径错误停止：1 failed、44 passed；
  首次定向因环境继承和来源签名旧路径记录 13 failed、387 passed、1 skipped。
  原始 `full.log`、`targeted_after.log`、`tool-help.json` 保留，未改写为成功证据。

最终运行前 `run-environment-final.json` 固定 182 个文件和依赖版本；实际测试经 `run_cris.ps1`
使用 Python 3.12.5 基础解释器、`-S` 和根 `.venv/Lib/site-packages`，不同于采集清单时的 venv 启动器。
`actual-runner-runtime.json` 记录本轮实际调用；运行后 `source-stability.json` 核对全部 182 个文件及
7 个依赖版本无漂移。`protected-after-check.json` 核对 5085 个受保护文件的路径、字节数及 SHA-256
全部一致，无新增或删除；范围包括 `webpage_info/`、`output_compare/`、`experiment_data/`、
`resources/isobyu/`。这不是操作系统逐次文件读取追踪。

调整前 145 个文件已无损归档并验证解压字节 SHA-256，见
[历史索引](history/README.md#2026-10-19测试整合前的源码说明和冻结记录)。旧测试可恢复，但不会留作有效计算入口。
结束后再次核对归档 CRC、压缩包 SHA-256 及 145 个原文件的解压字节，见 `archive-recovery-check.json`。
仅清理本轮列明的 7 个临时目录（7187 文件、1,188,324,131 字节），已验证路径均位于本轮证据目录内、
不含 reparse point；根层日志和审计记录保留。范围及清理后检查见 `cleanup-plan.json` / `cleanup-result.json`。
后继冻结 `freeze_revision=3` 仅更新工具/支撑文档路径及指纹；下载范围、独立输入和科学验收状态
保持原值，旧运行签名不复用。
本轮不新增官网下载、逐候选科学对照、长时 GenDB、真实用户浏览器操作或 GD TensorFlow 接收验收。
进一步架构检查只读；入口根定位、启动器默认目录、前端迟到响应和动作策略双登记见
[开发计划](DEVELOPMENT_PLAN.md)，未实施修复。

## 2026-10-19 架构调整（测试整合前）

本轮按计算状态安全、职责与调用方向、公开合同、验证可信度的优先级实施。
改动保留已有官网下载任务、原件及独立实验输入；没有以官网输出构造本地计算输入。
调整前的下载重建报告等六份原文已逐条校验归档，见
[历史索引](history/README.md) 和 [归档](history/architecture-before-20261019.zip)。

### 已完成调整

1. **动作事务和资源所有权。** `frontend/web` 经 `application/workflows` 执行有状态计算/选择时操作独立会话快照，
   输入校验、计算、响应构造及 UTF-8 JSON 编码全部成功才单次发布状态并增加 revision。
   异常保留原会话。ISO/SMODES/FINDSYM 各实例拥有私有运行目录；快照复制生成数据库及
   临时结果，避免对象已隔离而底层文件仍共享。构造失败、关闭和克隆失败均回收资源。
   真实回归的线程栈定位到 Windows 超时后无界等管道 EOF；改用文件捕获、有界 wait，
   并使超时实例失效，避免未知 Linux 后代继续写入的缓存被复用。
   新母相上传从新计算实例开始，已失效引擎也可恢复；成功才替换，失败保留旧输入和状态。
   另以 20 秒自然空闲反例定位到 WSL 的 `/tmp` 随发行版退出清空；对象仍存活且清理器未执行。
   WSL 暂存迁到持久文件系统的 `/var/tmp`，保留独占所有权和关闭回收，避免会话中途丢目录。
2. **导出任务唯一状态。** `workflows/jobs.py` 管理 job_id、来源 revision、进度、取消与终态；
   删除旧进度/取消/锁的平行全局实现。过期任务不能取消或覆盖新任务；HTTP/UI 传递同一身份。
   下载响应携带本任务终态，避免下一任务启动后把部分 ZIP 显示成全部完成。
   关闭服务先结束在途请求，再释放会话。坏 JSON、空/重复参数及不支持请求得到明确错误。
3. **求解职责归位并删除旧代码。** 删除 `method1/search_methods.py`，将原 Method 1/2/3
   计算体分别放入所属模块的 `search.py`；会话仅协调各 solver。Method 4 继续由
   `method4/decomposition.py` 负责。共有仿射格算法迁入 `utils/affine_embeddings.py`，
   相路径合同迁入 `models/phase_path.py`，删除旧位置文件；k 格助手、Hall/晶系查询统一来源。
   没有给旧单体增加兼容覆盖层。33 个原计算体的 AST 迁移比对保存在本轮证据目录。
4. **会话、导出和身份合同。** 导出批次发布移入 `session/publication.py`，精确子胞 Seitz
   变换移入 `export/child_frame.py`。ModeIdentityFailure/Status 与 Method 3 route/embedding
   进入正式数据模型，去掉动态私有字段和重复身份。nmod 共用整数/范围验证，不再截断或隐式转换。
   `select_path()` 收敛到 `search_method_2()`，删除另一套模式计算分支；保留当前 nmod/物种作用域，
   统一重复候选与诊断路线拒绝，修复显式类型选择后快照上下文过期的缺陷。
5. **GD 使用公开结果快照。** `DisplacementSnapshot` 固定候选、计算上下文、模式顺序/符号、
   晶格和未移原点超胞；GD 不再读取会话私有字段。输出身份只来自计算结果；含糊候选、
   不完整列、失效上下文及无效模式被拒绝。修复非对称 B 的倒易变换 `h_c=h_p Bᵀ`，
   格指数按精确行列式判断，不再四舍五入；未支持非整数指数显式拒绝。
   float32 转换后再次校验所有输出及倒易阵，超域结果在返回/写文件前拒绝。
6. **共享归一化根因修复。** 删除极小非零列绕过归一化的阈值分支和伪单位因子。
   同一算法验证有限非零模式、有效晶格及正整数计数；可表示域内保留模式方向和符号。
   补充必要分量下溢的明确拒绝，避免“整列非零”掩盖丢失方向后错误归一化；
   病态晶格的误差保证未扩大，正常极端尺度的既有容差未放宽。
   GD 明确使用所表示母胞的均方约定，不把它自动称作 ISO 原胞 As/Ap。
7. **比较与证据统一。** CIF 结构、字节及哈希来自一次不可变读取；严格模式要求
   语义、可信参考哈希和字节均匹配；缺哈希且无其他差异为 inconclusive，已知差异仍为 fail。
   手工工具共享快照、哈希、运行环境及原子 JSON 写入服务，各审计定义签名材料，尚非全部报告均原子写入。
   Method 4 及 Stage A/B 的解析和来源签名绑定同一文件快照，防止旧解析配新指纹；
   Method 4 起止来源漂移及 Stage A/B 恢复签名不匹配时拒绝旧结果。
   Method 3/4 准备器签名改为读取共享生产源码注册表及运行配置，删除易漏项的小列表；
   改变仿射求解实现后，旧通过 checkpoint 与旧输入批次必须失效，不能继续跳过计算。
   独立数学比较器保留独立实现，避免测试直接复用被验证算法。架构门禁阻止比较层
   启动生产计算及 GD 绕过公开会话合同，打包检查验证实际可安装入口。

### 科学影响、依据与适用条件

完整定义、逐步推导、反例、单位和容差见
[架构数学合同](reference_library/notes/ARCHITECTURE_CONTRACTS_20261019.md)。
`coverage.json` 和 `audit/implementation_audit.json` 的模块路径已同步，既有未证明/未支持状态保留。
实际审读范围限于该证明文档列出的项目推导段落；没有宣称通读全部外部文献。

- Seitz 共轭采用可逆有理 B、q 与 Fraction 算术；保持群乘法、恒等元、逆和闭包，
  非对角 B/非零 q 及完整 Hall 操作集测试检查坐标合同。精确路径无浮点容差。
- 行基 A_c=B A_p 对应倒易行指数 h_c=h_p Bᵀ，由相位配对直接推导；
  非对称剪切矩阵反例能区分旧错误算法。精确整数格指数不使用容差。
- 模式采用分数位移 U、以 Å 表示的晶格 A、基准长度 ℓ₀=1 Å，
  V=U ℓ₀ sqrt(N/Σ||UA||²)。检查均方、正负尺度、斜晶格、线性叠加及上界×dmax=1 Å。
  GD float32 回代使用 rtol=2e−7、atol=2e−8；双精度坐标复制 atol=1e−14、rtol=0。
- nmod 0..3 是实现支持范围，不是数学空间维数上限；非法类型/负值和更高维分别拒绝。
- 会话回滚、文件隔离、任务所有权与 SHA-256 同源性是工程不变量；守恒定律或能量最小化
  不适用于这些状态管理改动。它们也不能证明数值算法或参考来源真实。

### 验证命令与结果

在最终源码冻结后执行，两个进程均退出 0：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest ISODISTORT/tests E:/GD/tests_dev -q -ra -x --tb=short -o faulthandler_timeout=120 --basetemp=ISODISTORT/output/validation/architecture_20261019/full-final-temp
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1 -m pytest tests_repo -q -ra --tb=short --basetemp=ISODISTORT/output/validation/architecture_20261019/repo-final-complete-temp
```

- ISODISTORT + GD 全量：**1395 passed、4 skipped、2027 warnings，1186.18 秒**。
  日志 `full-final.log`。跳过项分别为缺少 X4− P3 官网/本地 CIF 配对、Method 3 旧输入冻结撤销、
  原生 POSIX 用例、LD1 g=1/6 的本地 ISO 子群数据库未生成。
  最后一项是本地计算数据库前置条件，不是漏读已下载的官网输出；本轮没有自动启动长时 GenDB。
  警告为依赖弃用及 CIF/元素解析提示，未发现未捕获线程异常或 finalizer 异常警告。
- 仓库架构、启动部署器和实际 wheel 安装门禁：**54 passed，12.11 秒**。
  日志 `repo-final-complete.log`；包含安装后规范入口与已删除旧入口的核查。
- 改动的 104 个 Python 文件运行 `ruff check --select F,E9`，通过。
  CRIS 与 GD 的 `git -c core.safecrlf=false diff --check` 通过；未提交或推送。
- 文档核查覆盖 11 个当前文档、2 个能力/证明映射，检查本地链接、模块路径和证明标签。
- 文件保护核查：4953 个登记官网文件、100 个保留的历史副本/上传输入、3 份独立母相输入、
  25 个通用套件文件指纹一致；调整前六份文档归档 CRC 和逐文件哈希核验通过。

源码快照 `run-environment.json` 登记 195 个源码/配置/测试文件及依赖版本；
`actual-runner-runtime.json` 记录实际 run_cris 使用的 Python 3.12.5 基础解释器、`-S` 与
`E:/CRIS/.venv/Lib/site-packages` 依赖来源，区别于采集清单使用的 venv 启动器。
`source-stability.json` 记录全量前后的源码、通用资源和依赖一致性。
这些指纹是运行材料清单，不声称进行了每一次操作系统文件读取的全链追踪。

定向最小反例、修前失败和修后结果亦保留：归一化下溢、Windows 后代持有输出句柄、
WSL 自然空闲目录丢失、select_path 双入口不一致、准备器漏签名与文件快照漂移。
它们与全量测试有重叠，不累加通过数量。
本轮结束后仅清理明确列出的 57 个临时目录（9186 文件，约 1.48 GB），
保留根层日志、来源清单、反例探针与核查报告；范围及结果见 `cleanup-plan.json`、`cleanup-result.json`。

原始日志、AST 迁移记录与本轮文件核查位于
`ISODISTORT/output/validation/architecture_20261019/`。
第一次全量回归因审查发现归一化反例而主动中止。
第二次记录两项失败（4310 实际模式导出、nmod 错误文案），在末尾 WSL 超时处理发生
无界输出管道等待；取得线程栈后精确停止该 pytest 进程。两次都不是全绿，原始日志保留。
4310 案例单独复查曾通过，但第三次联合回归仍在该项失败（748 passed、3 skipped 后停止），
取得了暂存目录不存在的明确异常；随后用自然空闲反例定位到 tmpfs 生命周期。
单项通过不能替代整套验证；最终整套复查结果见上。

### 当前官网下载与科学验收状态

下载/检验任务仍为批次 `20261019-official-rebuild-v1`：Method 2 共 102 查询，
99 个页面集未齐、2 个部分齐、1 个齐，已知缺 5 张结果列表，已有 94 个四件套复用。
精确路径和动作以 [下载清单](DOWNLOAD_CHECKLIST.md) 为准；本轮没有执行新官网下载。
下载整理的历史测试及文件保护数字保留在上述归档，不能替代当前代码回归。

独立输入和有限界限未核定的 Method 3/4 任务仍 blocked；新逐候选科学验收仍 not_run。
真实用户浏览器全操作链、GD 的 TensorFlow 计算接收、长时 GenDB 中途取消未作本轮人工验收。
真实 WSL 与自动化 HTTP 测试的结果仅按上述检查各自范围解释。
WSL 混合编码诊断文本仍可能乱码，已在开发计划登记；这项文本解码不属于本轮已完成修复。
不宣称已证明全部 230 空间群、无限 k 范围、外部二进制或所有未支持物理表示。

支撑文档/工具发生架构变化后，冻结记录登记为后继修订 `freeze_revision=2`，
保留原冻结文件的可恢复归档及逐项支撑修订；
下载任务、参考原件及输入指纹不因支撑修订变化，旧运行签名不复用。
后续工作只维护在 [开发计划](DEVELOPMENT_PLAN.md)。
