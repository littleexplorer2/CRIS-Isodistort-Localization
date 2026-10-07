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
- 网页和 Python API 只能共用 `backend/` 与 `features/` 的计算与导出逻辑，网页交互壳
  不得复制算法。英语网页文案统一由 `frontend/i18n/messages.py` 提供。
- 运行时默认值先写 `resources/config/settings.yaml`，代码通过 `get_config()` 读取。

## 四部分结构与职责

本项目源码按四个部分组织；修改某一功能时应只改对应部分，跨部分共用的底层进入 `backend/`：

1. `backend/` —— 通用底层计算与共用设施（所有 Method 通用）：
   - `backend/wrappers/`：iso / findsym / smodes 二进制封装、WSL 进程与暂存、ISO 缓存；
   - `backend/models/`：模式与来源领域模型；
   - `backend/tables/`：CDML / Kovalev / 官网 k 点与空间群查表；
   - `backend/utils/`：精确有理数晶格与群论工具、文本解析、OPD 文本、母相页头、
     配置加载、异常、自检；
   - `backend/api/`：`IsoDistort` 会话 API 入口。
2. `features/` —— 按程序功能划分的中端包，每个包只含该功能的实现：
   - `features/input_cif/`：CIF 解析、坐标变换、对称性校验；
   - `features/method1/`：Method 1 搜索与共用畸变引擎（affine embedding、相路径、
     畴、occupational、`search_methods`）；
   - `features/method2/`：Method 2 的 (3+d)/公度锁定模式计算（`superspace`）与母相→
     子胞模式映射（`distortion_mapper`）；
   - `features/method3/`：Inverse Landau/COPL 可行性与 coupled 见证（`inverse_landau*`、
     `coupled_routes`）；
   - `features/method4/`：分解、均匀应变与生成回代（`strain*`、`distortion_engine`）；
   - `features/export/`：四种下载格式 writer 与 CIF 位移/来源共享合同。
3. `frontend/` —— 网页支撑：`frontend/web/`（HTTP 服务、`index.html`、`static/`）与
   `frontend/i18n/`（英文文案）。
4. 其他全局共享/非代码：`resources/config/`、`resources/isobyu/`、`docs/`、`output/`、
   `tests/`、`scripts/`，以及留在根目录的 `.gitignore`、`pyproject.toml`、
   `requirements*.txt`、`README.md`、`agent.md`。

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

- 唯一职责：通用底层与二进制封装在 `backend/`，各 Method 的搜索与模式在
  `features/method1…4/`，会话 API 在 `backend/api/`，导出在 `features/export/`，
  网页在 `frontend/web/`。新增共用算法一律下沉 `backend/`，不得在 `features/` 之间
  交叉复制。
- `features/` 各包只允许通过 `backend/` 与彼此公开的模块接口协作；包间导入使用绝对
  路径（`from features.method2.superspace import …`），同包内才用相对导入。
- 状态必须显式传递；候选身份不能只靠目录名、列表下标或对象 `id()`；对照与
  导出须使用各 Method 适用的完整科学身份和查询上下文。并发变更须遵守
  session lock/revision，不复用不兼容的模式缓存。
- 默认比较科学语义，不要求整文件逐字节相同。CIF 应可被 VESTA 解析，
  `.isoviz` 应满足 IsoVIZ 结构；未安装 GUI 时只能报告静态验证。
- 修改算法/导出后至少运行相关 pytest；交付前从 CRIS 根目录运行完整
  `.\run_cris.ps1 -m pytest ISODISTORT\tests -q --tb=line`。
  官网差分与长时命令以 `docs/MANUAL_VALIDATION.md` 为准。
