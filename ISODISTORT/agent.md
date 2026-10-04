# ISODISTORT agent rules

适用于 `ISODISTORT/`。先读仓库根 `AGENTS.md`。

## 项目目标与范围

- 本项目是在现有本地网页版功能范围内对
  [ISODISTORT 官网](https://landau3.byu.edu/isodistort.php)进行本地化的逆向工程。
  当前产品范围以本地网页版已经提供的功能为准；官网存在而本地网页版没有的
  功能暂不新增，包括官网内置在线晶体结构查看和 Distortion Generate / Domains。
  存档 HTML 只作验证证据，不扩大产品范围。
- 当前重点是全面审查、debug、结果准确性和算法泛用性；目标是在受支持输入与
  上述功能范围内对任意适用晶体给出正确结果。有限母相或候选通过只证明已覆盖
  范围，不能宣称所有晶体或跨晶系科研级精度。

## 文档路由

根 `AGENTS.md` 的文档职责继续适用；本子项目另指定：

- 官网批次输入：`docs/manifests/`
- 机器报告/checkpoint：`output/validation/`（不提交）

## 实现边界

- 官网黄金集只能读和审计；新输出写到约定的本地结果/验证目录。
- 网页、终端和 Python API 只能共用 `isocore` 计算与导出逻辑，交互壳不得
  复制算法。英语交互文案统一由 `isocore/i18n/messages.py` 提供。
- 运行时默认值先写 `config/settings.yaml`，代码通过 `get_config()` 读取。

## 科研与验收

1. 官网输出是独立验证证据，不是算法来源。开始依赖官网材料的审查或修改前，
   先核验结果表、查询条件、候选全集及每个候选的完整身份；发现漏下、误下、
   错放或覆盖时，先列明需重下或归位部分，并暂停依赖该证据的结论与修改。
   计算逻辑必须来自适用的定义、公式与约定。
2. 每项计算须从适用的晶体学/固体物理定义与公式推导：空间群仿射作用、
   直接/倒易格对偶、k-star/小群、表示与 IR、OPD、Wyckoff 轨道、超胞、
   setting/origin、调制与 superspace 等。优先使用精确整数/有理数，并检查
   群闭包、行列式/指数、轨道重数、维数和等价类等不变量。
3. 必须按适用范围查阅 Stokes、Campbell、Hatch 的 ISOTROPY/ISODISTORT 论文，
   同时使用 International Tables、权威晶体学/固体物理文献和官方帮助交叉验证；
   官方帮助与论文 DOI 的引用入口见 `docs/BUGFIX_VALIDATION_REPORT.md` 的“主要科学
   依据”。第三方库的 Hall/HM、原/惯用胞、主动/被动变换、行/列约定必须显式核对。
4. 浮点容差要有数值或物理尺度依据。每个正式验证批次记录问题、输入来源、
   公式/约定、单位、源码与依赖指纹、容差、正反例、逐例状态及失败原因。
5. 无法建立推导链或充分证据的候选不得进入默认可交付路径或宣称“与官网一致”；
   诊断候选必须 `selectable=false`，直到真实 route/modes 已解析。

## 架构与验收

- 唯一职责：后端包装在 `isocore/backend/`，搜索与模式在
  `isocore/distortion/`，会话 API 在 `isocore/api/`，导出在 `isocore/io/`。
- 状态必须显式传递；候选身份不能只靠目录名、列表下标或对象 `id()`；对照与
  导出须使用各 Method 适用的完整科学身份和查询上下文。并发变更须遵守
  session lock/revision，不复用不兼容的模式缓存。
- 默认比较科学语义，不要求整文件逐字节相同。CIF 应可被 VESTA 解析，
  `.isoviz` 应满足 IsoVIZ 结构；未安装 GUI 时只能报告静态验证。
- 修改算法/导出后至少运行相关 pytest；交付前运行完整
  `.\.venv\Scripts\python.exe -m pytest ISODISTORT\tests_dev -q --tb=line`。
  官网差分与长时命令以 `docs/MANUAL_VALIDATION.md` 为准。
