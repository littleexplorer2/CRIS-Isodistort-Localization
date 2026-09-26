# ISODISTORT agent rules

适用于 `ISODISTORT/`。先读仓库根 `AGENTS.md`。本文件只保存稳定规则；开发
状态见 `docs/DEVELOPMENT_PLAN.md`，已完成修复与实测数字见
`docs/BUGFIX_VALIDATION_REPORT.md`。

## 文档路由

- 用户能力、安装、交互、已知限制：`README.md`
- 唯一开发计划、阶段和完成标准：`docs/DEVELOPMENT_PLAN.md`
- 修复与验证证据：`docs/BUGFIX_VALIDATION_REPORT.md`
- 官网人工下载：`docs/DOWNLOAD_CHECKLIST.md`
- 手工/长时命令：`docs/MANUAL_VALIDATION.md`
- 官网批次输入：`docs/manifests/`
- 端口、路径、容差、预算：`config/settings.yaml`
- 机器报告/checkpoint：`output/validation/`（不提交）

不要在本文件记录进度、候选数、批次结论、bug 历史、下载状态或使用教程。

## 不可违反的边界

- 只读：`../experiment_data/`、`../webpage_info/`、`../output_compare/`、
  `isobyu/`、仓库内 `../GD/`。
- 官网黄金集只能读和审计；新输出写到约定的本地结果/验证目录。
- 网页、终端和 Python API 只能共用 `isocore` 计算与导出逻辑，交互壳不得
  复制算法。英语交互文案统一由 `isocore/i18n/messages.py` 提供。
- 运行时默认值先写 `config/settings.yaml`，代码通过 `get_config()` 读取。
- Method 4 服从开发计划的启动门禁；未收到明确开始指令时只维护计划，不运行、
  下载或修改其实现。

## 科研与实现门禁

1. 不死背官网答案，不得按材料、空间群、IR、OPD 或案例目录写特判。官网输出
   是验证证据，不是算法来源。
2. 每项计算须从适用的晶体学/固体物理定义与公式推导：空间群仿射作用、
   直接/倒易格对偶、k-star/小群、表示与 IR、OPD、Wyckoff 轨道、超胞、
   setting/origin、调制与 superspace 等。优先使用精确整数/有理数，并检查
   群闭包、行列式/指数、轨道重数、维数和等价类等不变量。
3. 必须按适用范围查阅 Stokes、Campbell、Hatch 的 ISOTROPY/ISODISTORT 论文，
   同时使用 International Tables、权威晶体学/固体物理文献和官方帮助交叉验证；
   第三方库的 Hall/HM、原/惯用胞、主动/被动变换、行/列约定必须显式核对。
4. 浮点容差要有数值或物理尺度依据。每个正式验证批次记录问题、输入来源、
   公式/约定、单位、源码与依赖指纹、容差、正反例、逐例状态及失败原因。
5. 无法建立推导链或充分证据时标为假设/诊断/`inconclusive`；不得进入默认
   可交付路径或宣称“与官网一致”。诊断候选必须 `selectable=false`，直到真实
   route/modes 已解析。

## 架构与验收

- 唯一职责：后端包装在 `isocore/backend/`，搜索与模式在
  `isocore/distortion/`，会话 API 在 `isocore/api/`，导出在 `isocore/io/`。
  新建模块前先用 `rg` 查现有实现。
- 状态必须显式传递；候选身份不能只靠列表下标或对象 `id()`；并发变更须遵守
  session lock/revision，不复用不兼容的模式缓存。
- 默认比较科学语义，不要求整文件逐字节相同。CIF 应可被 VESTA 解析，
  `.isoviz` 应满足 IsoVIZ 结构；未安装 GUI 时只能报告静态验证。
- 修改算法/导出后至少运行相关 pytest；交付前运行完整
  `.\.venv\Scripts\python.exe -m pytest ISODISTORT\tests_dev -q --tb=line`。
  官网差分与长时命令以 `docs/MANUAL_VALIDATION.md` 为准。
- 同步规则：用户行为/限制改 `README.md`；默认值改 yaml；计划改
  `DEVELOPMENT_PLAN.md`；已修问题与结果改 `BUGFIX_VALIDATION_REPORT.md`。
