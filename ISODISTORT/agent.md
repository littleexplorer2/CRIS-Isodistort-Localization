# ISODISTORT agent rules

适用于 `ISODISTORT/`。先读仓库根 `AGENTS.md`。

## 项目目标与范围

- 本项目是在现有本地网页版功能范围内对
  [ISODISTORT 官网](https://landau3.byu.edu/isodistort.php)进行本地化的逆向工程。
  当前产品范围以本地网页版已经提供的功能为准；官网存在而本地网页版没有的
  功能暂不新增，包括官网内置在线晶体结构查看和 Distortion Generate / Domains。
  存档 HTML 只作验证证据，不扩大产品范围。
- 算法必须适用于声明的输入域与上述功能范围。有限母相或候选通过只证明已覆盖
  范围，不能宣称所有晶体或跨晶系科研级精度。

## 文档路由

根 `AGENTS.md` 的文档职责继续适用；本子项目另指定：

- 官网批次输入：`docs/manifests/`
- 机器报告/checkpoint：`output/validation/`（不提交）
- 测试比较参数：`tests/validation/config/settings.yaml`；操作命令见
  `docs/MANUAL_VALIDATION.md`，使用与限制见本项目 `README.md`。

## 实现边界

- 新输出必须写到约定的本地结果/验证目录；官网数据隔离执行根规则，不设子项目例外。
  官网结果列表页只写入 `webpage_info/<母相>/` 页面存档（命名见 `webpage_info/README.md`，
  改名/哈希登记在 `_PAGE_INDEX.json`）；`output_compare/` 的查询/候选目录只按真实列表行
  建立，未知列表保持空白，禁止用本地网页版输出补，也禁止覆盖既有归档或复制进计算输入。
- 网页和 Python API 共用 `application/session/` 的计算与导出协调，交互工作流统一由
  `application/workflows/` 提供；网页只处理传输、呈现和浏览器生命周期。
  共享英文文案统一由 `backend/messages/` 提供，前端通过应用层取得。
- 项目仅提供英文界面；不实现语言切换、语言参数或多语言目录。
  `/api/messages` 只提供英文消息模板，不参与计算参数或科学标签转换。
- 运行时默认值先写 `resources/config/settings.yaml`，代码通过 `get_config()` 读取。

## 三层模块职责

依赖方向、目录清单、兼容入口与架构门禁以根 `AGENTS.md` 和
`tests_repo/architecture_manifest.json` 为准。本子项目职责细分如下：

- `backend/wrappers/`：iso/findsym/smodes 调用、WSL 暂存与 ISO 缓存；
  `models/`：模式、来源与共享物理类型/相路径合同；`tables/`：CDML、Kovalev、k 点与空间群/Hall 表；
  `utils/`：精确晶格、仿射嵌入、k/格相容、解析、配置、异常、自检；`messages/`：共享英文消息。
- `backend/input_cif/`：CIF 读取、坐标变换与对称性验证；
  `method1/`：特殊 k 搜索、畴与 occupational 边界；
  `method2/`：指定候选求解、superspace 与母相到子胞映射；`method3/`：任意 k/格搜索、Inverse Landau/COPL 与
  coupled 见证；`method4/`：模式分解（`decomposition.py`）、应变和生成回代；
  `export/`：格式合同与 writer。
- `application/session/`：唯一 `IsoDistort` 会话 API、计算调用与批量发布协调；
  `application/workflows/`：交互动作、会话 revision/锁、上传和导出任务的状态/错误合同。
- `frontend/web/`：HTTP、JSON 与文件响应、静态页面、浏览器启动及标签页生命周期。

另设“根目录保留文件”这一位置优先类别：不论文件在语义上属于入口、配置、依赖、文档
还是规则，都固定直接放在 `ISODISTORT/`，不得按功能移入上述目录。当前成员为
`.gitignore`、`agent.md`、`main_requirement.py`、`main_web.py`、`pyproject.toml`、`README.md`、
`requirements-dev.txt`、`requirements.txt`。

## 科研与验收

根规则中的输入隔离、每次修改的理论核查及历史假设标注为强制前置条件；另须遵守：

1. 官网比较前必须核验结果表、独立输入记录、候选全集及每个候选的完整身份；
   发现漏下、误下、错放或覆盖时，必须列明补齐事项，并暂停依赖该证据的结论。
   核验过程只能用于比较/审计，不得把候选全集或结果字段回填计算；不依赖这些
   结果的理论推导和修复可以继续。
2. 各 Method 的推导必须核对所涉及的空间群仿射作用、
   直接/倒易格对偶、k-star/小群、表示与 IR、OPD、Wyckoff 轨道、超胞、
   setting/origin、调制与 superspace 等。优先使用精确整数/有理数，并检查
   群闭包、行列式/指数、轨道重数、维数和等价类等不变量。
3. 必须按适用范围查阅 Stokes、Campbell、Hatch 的 ISOTROPY/ISODISTORT 论文，
   同时使用 International Tables、权威晶体学/固体物理文献和官方帮助交叉验证；
   官方帮助、论文 DOI、能力映射和补充证明统一见 `docs/reference_library/README.md`；
   使用前核对书目中的实际审读层级及证明的输入假设。第三方库的 Hall/HM、原/惯用胞、
   主动/被动变换、行/列约定必须显式核对。
4. 正式验证批次除根规则要求的来源/指纹外，必须记录各 Method 的查询上下文、
   完整候选身份、逐例状态及失败原因；不得把缺失项合并为通过。
5. 诊断候选必须 `selectable=false`。只有理论依据和输入假设通过核查，且真实
   route/modes 已解析并验证后，才能转为可选择的计算结果。

## 架构与验收

- 跨模块使用绝对导入，同模块内部可相对导入；模块目录及入口变更须同步
  `tests/` 根目录审计脚本的源码指纹，不能让报告只签名兼容壳而遗漏实际实现。
- 状态必须显式传递；候选身份不能只靠目录名、列表下标或对象 `id()`；对照与
  导出须使用各 Method 适用的完整科学身份和查询上下文。并发变更须遵守
  session lock/revision，不复用不兼容的模式缓存。
- 默认比较科学语义，不要求整文件逐字节相同。CIF 应可被 VESTA 解析，
  `.isoviz` 应满足 IsoVIZ 结构；未安装 GUI 时只能报告静态验证。
- `tests/validation/` 的单对与批量比较共用同一比较器；只读取已冻结的本地输出和
  参考输入，不创建或重写输入目录。严格模式必须同时满足语义、参考哈希和字节条件；
  空批次、解析失败与非有限容差不得报告科学验证通过。比较结果只覆盖声明的比较域。
- 修改代码或计算配置后必须先运行相关 pytest；交付前从 CRIS 根目录运行完整
  `.\run_cris.ps1 -m pytest ISODISTORT\tests -q --tb=line`。
  官网差分与长时命令以 `docs/MANUAL_VALIDATION.md` 为准。
