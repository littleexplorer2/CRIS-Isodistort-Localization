# CRIS agent rules

进入子项目后还要读取该目录的 `agent.md`；使用方法读各自 `README.md`。

## 边界

- 永远只读：`experiment_data/`、`webpage_info/`、`output_compare/`、
  `ISODISTORT/resources/isobyu/`、仓库内 `GD/`。
- 与 CRIS 同级的 `GD/` 服从其自身规则；禁止修改 tianren 参考笔记本及
  `D:\OneDrive\...` 原始数据。
- 可改主体：`ISODISTORT/`（除 `resources/isobyu/`）、`ISOVIZ_INPUT/`。
- 两个运行项目及集成测试共用根目录内的实体 `.venv`，除非用户明确要求，不重建环境。普通
  Python 命令使用 `.\.venv\Scripts\python.exe`；本机从 OneDrive 内的解释器调用
  WSL 会触发 `Wsl/E_ACCESSDENIED`，需要 WSL 的 ISODISTORT 命令改用根目录
  `.\run_cris.ps1 <Python 参数>`，包括网页入口与相关测试；其依赖仍只从 `.venv`
  加载。项目仅使用该实体 `.venv`，不要另行创建环境。
- 根部署器与仓库编排测试统一放在 `tests_repo/`；科学测试和比较工具统一放在
  `ISODISTORT/tests/`，ISOVIZ_INPUT 使用自身 `tests_dev/`。
  根门禁运行 `.\run_cris.ps1 -m pytest tests_repo -q --tb=line`。
- 本机 GUI 快捷方式与可执行启动器统一放在 `local_tools/launchers/`，该目录内容不提交；
  `ISOVIZ` / `ISOVIZ_JAR` 环境变量仍可显式覆盖默认查找位置。

## 科学计算硬性规则（所有子项目及验证脚本均适用）

### 官网结果与计算输入必须隔离

- 官网计算结果文件只能由独立的只读比较/审计环节读取，禁止进入被验证的本地
  计算链。计算链包括预处理、搜索、筛选、模式构造、归一化、拟合、缓存和导出。
  禁止通过测试夹具、mock、环境变量、旁路文件、复制字段或预计算缓存间接注入答案。
  生产运行与验证运行必须遵守同一隔离规则，不得设置仅验证时启用的答案通道。
- 官网材料中仅可使用用户在官网也必须主动指定的输入参数，并须有独立的输入记录：
  例如用户提供的母相结构、k 向量、所选 IR/OPD、主动指定的基变换或原点。
  必须按当前 Method 的任务定义判定输入/输出角色，不能因为某字段在另一 Method
  中可作输入，就把本任务应求解的量改称参数。禁止从结果文件反推、挑选或补齐输入。
  官网生成的候选全集、子群/setting 答案、模式基、振幅、归一化因子、畴或映射均
  不得用于替代本地求解；上传或改名不改变这些数据的来源与角色。
- 用户原始实验结构等独立任务输入，以及有出处和版本的通用空间群、表示、Wyckoff
  等正式数据库可以参与计算；它们不得夹带按具体案例存储的官网答案。
  `ISODISTORT/resources/isobyu/` 的通用程序/数据库与按允许输入执行它们所得的本地
  结果，不属于官网案例答案。ISOVIZ_INPUT 对显式输入文件的解析/振幅写入测试，
  只验证该转换合同，不得据此宣称上游结构或模式已经独立求解正确。
- 每次科学对照验证必须先固定输入及参数来源，再独立运行本地计算并冻结输出，
  最后才读取官网结果进行比较。验证记录必须列出输入、参考答案、实际读取的外部
  资源、配置/源码/依赖版本及输出的路径与指纹；核查计算路径不存在答案读取或回填。
  发现污染必须将该次验证标为无效，隔离受污染的派生产物并使相关缓存失效，再按隔离
  流程重跑；必须保留审计证据，禁止改动原始数据、黄金文件或用污染结果宣称通过。

### 每次修改必须完成理论影响核查

- 修改代码、配置、接口、测试或科学说明前，必须确定受影响的计算合同及其理论依据；
  修改后必须核对实现、调用方、默认值和验证是否仍满足该合同。纯展示、工程或规则
  修改也必须检查科学语义影响；只有确认不涉及科学语义，才能记录理由并只执行相应
  工程/文档检查。输入选择、单位标签、默认参数、坐标约定和数据传递不能按纯界面
  修改豁免。验证记录进入对应报告或交付说明，不写入规则文件。
- 每次交付必须记录修改点、科学影响判断、依据及适用条件、推导与实现的对应关系、
  采用的容差、验证命令/结果和剩余缺口；无科学影响时必须写明理由及实际执行的检查。
  科学验证还必须记录公式/单位约定、不变量和正反例证据，禁止仅写“已验证”或通过数。
- 涉及科学语义时，必须明确数学对象、物理/晶体学适用条件、符号、单位、坐标及
  归一化约定，给出可追溯的文献推导或本项目完整推导，再逐步核对代码的输入、变换
  与输出。新推导及适用边界写入 `ISODISTORT/docs/reference_library/` 的证明文档，
  同步能力映射；仅有引用标题、官网相符或测试通过，不构成理论证明。
- 必须按改动选择数学、物理和晶体学检验：例如群作用/闭包、晶格对偶、指数与
  重数、维数、对称性、坐标协变性、单位一致性、适用的守恒关系、正逆变换与极限。
  必须说明不适用项及原因。修复必须包含原缺陷的最小反例；新增算法须覆盖边界、
  不支持输入及能区分正确算法与错误猜测的例子，不能只重复实现或只对照单个答案。
- 禁止将无理论依据的猜测、按材料/标签/期望结果分支、为提高匹配率拼接的经验公式
  或调参补丁纳入生产计算。数值近似必须有推导、误差依据与适用域；合法的参数拟合
  必须先定义物理模型、待估参数和目标函数，不能用官网答案拟合出替代算法。
  容差必须依据数值误差或物理尺度，禁止为掩盖失败放宽容差、删例或改写期望值。
- 理论链缺失、假设不成立或数据不足时，必须明确返回不支持、诊断或 `inconclusive`，
  不得静默采用猜测、伪造成功或允许诊断候选作为正式结果导出。理论核查和实现测试
  均完成前不得宣称科学验收通过；静态架构、解析成功和 GUI 启动均不能替代它们。

### 历史假设与权威结论必须区分

- 科学断言必须区分“官方文档陈述”“文献结论”“本项目推导”“历史假设/待验证”
  或“有适用域的近似”。引用须给出可定位的来源、章节/公式与版本，并标明实际审读
  范围；没有读到全文或对应证明，禁止声称已核验。官方行为说明不自动构成数学证明。
- 发现或继续使用历史假设时，必须在相邻代码注释/科学说明中标明“历史假设，未获
  证明”，并链接其来源记录、适用条件、验证状态和缺口。来源不明须明确写“来源不明”；
  已被反例推翻的假设须标为“已否定”，并从有效计算路径移除，不能只改说明后继续使用。
- 历史实现、旧报告、注释、模型生成的推导、官网单例一致或重复引用均不得作为官方
  证明。只有补齐可审查的推导或精确文献依据并完成适用性核查，才能更新其证据状态；
  必须保留来源链。未证明假设只能用于明确隔离的探索/诊断，不得进入生产计算。
- 新代理必须先检查所改模块在 `ISODISTORT/docs/reference_library/coverage.json`
  中的来源与证明、对应审计记录及邻近假设标注；没有记录不代表已获证明。
  历史归档仅供追溯，不能覆盖当前科学合同。

## 文档职责

- 根 `README.md`：跨项目总览。
- 子项目 `README.md`：用户安装、使用、能力与限制。
- `ISODISTORT/resources/config/settings.yaml`：运行时默认值的事实来源。
- 三层模块目录、依赖方向与兼容边界见“三层架构基线”；可执行清单为
  `tests_repo/architecture_manifest.json`。
- `ISODISTORT/docs/DEVELOPMENT_PLAN.md`：唯一开发计划和待办，固定包含两个常驻分区
  “待修的程序漏洞”和“待执行的优化计划”；正文只维护当前与后续待办。
- `ISODISTORT/docs/BUGFIX_VALIDATION_REPORT.md`：当前验证缺口、在用证据入口与科学依据。
- `ISODISTORT/docs/reference_library/`：全项目理论资料与证明；`sources.json` 是书目
  主数据，`coverage.json` 关联能力、来源、证明和模块；`audit/` 保存理论与实现的
  审计依据。书目衍生文件必须由目录内脚本生成。
- 已完成的 debug、下载和修改历史无损压缩到 `ISODISTORT/docs/history/`；索引记录
  原路径、SHA-256 和恢复方法，归档须验证可恢复后再精简正文。有效操作说明和在用
  机器 manifest 保持可直接读取，不重复维护历史叙述。
- `ISODISTORT/docs/DOWNLOAD_CHECKLIST.md`：人工下载动作与目录。
- `ISODISTORT/docs/MANUAL_VALIDATION.md`：长时/手工验证命令。
- 根规则只维护全仓硬性边界，子项目 `agent.md` 只补充专属职责与门禁；规则文件
  禁止追加修复流水账、批次结果、下载状态或临时进度。删除重复规则时必须保留唯一
  权威入口，不得同时删除仍适用的约束。

## 三层架构基线（后续修改必须遵守）

- 调用方向固定为 `frontend → application → backend`。底层不依赖应用层或前端；
  应用层不依赖前端；前端只经应用层调用计算、文件工作流与运行环境查询。
- 每个模块必须有独立文件夹和 `__init__.py` 入口；模块内部可按职责拆文件。
  跨模块复用现有公开入口或其明确公开的子模块，不复制实现；状态通过调用合同传递。
  ISODISTORT 原有科学子模块的文件级接口登记在清单 `public_submodules` 中；
  少量既有下划线命名合同显式登记在 `legacy_scientific_interfaces`，仅保留原调用。
  未登记的实现文件和包内私有符号不得成为其他模块的新依赖；新功能从包公开入口提供。
- `backend` 保存科学算法、数据模型、格式读写和运行基础设施；`application` 保存
  会话、用例协调、事务/锁与结果合同；`frontend` 保存 HTTP/CLI/GUI 输入输出适配。
  `resources`、配置、文档、测试、输出不是额外的业务层。
- 新增、迁移或重命名模块时，同步可执行清单、调用方、打包配置、兼容入口、
  来源签名和对应文档。所有后续改动以以下已实现目录为基线。

```text
ISODISTORT/
  backend/
    wrappers/ models/ tables/ utils/ messages/
    input_cif/ method1/ method2/ method3/ method4/ export/
  application/
    session/                 # IsoDistort Python 会话与计算/导出协调
    workflows/               # 交互用例、revision、上传、导出任务与进度
  frontend/
    web/                     # HTTP、页面/静态资源、浏览器与标签页生命周期
ISODISTORT/tests/validation/   # 测试支持包，不是产品、独立安装包或新业务层
  backend/
    config/ paths/ cif/ comparison/ results/
  application/
    single/ batch/ workspace/
  frontend/
    cli/
ISOVIZ_INPUT/isoviz_input/
  backend/
    amplitudes/ inputs/ validation/ paths/ launcher/ config_loader/
  application/
    pipeline/ input_selection/ environment/
  frontend/
    cli/ gui/ entrypoint/
```

- ISODISTORT 只提供规范三层入口，不恢复 `features`、`backend.api` 或 `frontend.i18n`
  旧架构转发。Python/GD 会话调用 `application.session`，导出计算调用
  `backend.export.distortion_formats`；变更必须同步外部消费者及其测试。
- ISOVIZ_INPUT 包根的旧平铺模块仅作既有兼容转发，生产层不得经它们绕过分层。
- 根 `setup_cris.py`、`run_cris.ps1` 与子项目 `main*.py` 保留安装/启动入口职责；
  不在入口复制领域算法或应用工作流。
- 验证唯一实现是 `ISODISTORT/tests/validation/`，pytest 与 `tests/manual/` 的审计
  脚本必须调用其公共入口；生产三层禁止导入任何测试/比较代码。禁止重建独立验证
  产品、安装包或复制比较器；只读比较不能启动生产计算或回填参考答案。
- 所有管线变更必须核对 HTTP → workflows → session → backend 的实际调用。
  Method 4 分解由 `backend/method4/decomposition.py` 负责；心化矩阵及指数统一来自
  `backend/utils/lattice.py`，超胞站点映射统一来自 `backend/input_cif/`。
  删除无用代码前必须核对入口、动态导出、测试、打包、手工签名及
  外部 GD 调用；仅名字以 `_` 开头或静态搜索没有命中不能单独证明可以删除。
- 架构门禁：`.\run_cris.ps1 -m pytest tests_repo/test_layer_boundaries.py -q --tb=line`。
  代码/模块/打包变更必须执行该门禁，并按子项目规则运行相应回归与打包验证。
  仅文档或规则变更必须检查引用、规则一致性和 `git diff --check`；若改变科学合同，
  还须完成上述理论核查及受影响的实现验证。
