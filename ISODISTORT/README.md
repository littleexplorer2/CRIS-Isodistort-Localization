# ISODISTORT 0.4.0（本地版）

本目录在现有本地网页版的功能范围内，对
[ISODISTORT](https://landau3.byu.edu/isodistort.php)（BYU 的晶体畸变搜索工具）
进行可离线运行的本地化实现。官网存在但本地网页版没有的功能暂不属于逆向范围。
当前开发聚焦于 debug、结果准确性和对任意适用晶体的泛用性。计算按空间群作用、
直接/倒易格、k-star、小群与表示、Wyckoff 轨道、setting/origin、超胞及 superspace
等定义逐步推导，并以 Stokes、Campbell、Hatch 等人的 ISOTROPY/ISODISTORT 论文、
官方帮助和晶体学不变量交叉核验；官网文件用于验证结果，不作为按样例拼公式的来源。

你不需要先懂晶体学术语。用白话理解工作流即可：

1. 准备一份**母相**晶体结构文件（CIF）：描述「还没畸变」时的对称性与原子位置。  
2. 告诉程序要考虑哪些**畸变类型**（应变、原子位移、占位、磁、旋转等）。  
3. 用 Method 1 / 2 / 3 **搜索**可能的子群（对称性降低后的候选结构列表）；或用 Method 4 把已经畸变的结构**分解**成模式幅度。  
4. 在 **Distortion** 区下载筛选后的结果表，以及子群结构文件（CIF、IsoVIZ、模式详情、TOPAS 等）。

计算后端是 ISOTROPY Suite 的 Linux 程序 `iso` / `findsym` / `smodes` 与 `data_*.txt` 数据库（放在本目录的 `resources/isobyu/`，只读）。同一套引擎支持两种用法：**网页**与 **Python API**。网页界面为**英语**（没有语言切换）。界面上看到的英文标签，下文一律用引号标出。

本项目与仓库根目录的实体 `CRIS/.venv` 共用一份 Python 虚拟环境。本机的标准
Windows venv 启动器从 OneDrive 调用 WSL 会被拒绝；统一安装器会先做实际 A/B
探针，只在确认“venv stub 失败而基础 Python 成功”后换成指向基础 Python 的硬链接，
因此网页可继续用 `.venv\Scripts\python.exe` 直接启动和计算。原 stub 会保留为回退
副本，`run_cris.ps1` 仍作为无法建立硬链接时的兼容入口。子项目关系见根目录
[README.md](../README.md)；开发规则见 [agent.md](agent.md)，开发计划见
[docs/DEVELOPMENT_PLAN.md](docs/DEVELOPMENT_PLAN.md)，已完成修复和验证结果见
[docs/BUGFIX_VALIDATION_REPORT.md](docs/BUGFIX_VALIDATION_REPORT.md)。

---

## 1. 这个工具具体做什么

| 步骤 | 你做什么 | 得到什么 |
| --- | --- | --- |
| 加载母相 | 在 "Parent CIF" 选文件并点 "Load" | 页头显示空间群、晶格、Wyckoff 位点 |
| 设畸变类型 | 勾选后点 "Change" | 后续 Method 按你勾选的类型过滤 |
| Method 1–3 | 设过滤条件，点 "OK" | 一张可筛选、排序的子群结果表；点一行可看模式 |
| Method 4 | 上传女儿相 CIF，点 "OK" | 官网口径 `As/Ap`、原始拟合系数、`normfactor` + RMS residual（不是子群列表） |
| Distortion | 选 Method、下表 / 勾格式、下 ZIP | 筛选后的 txt/csv；Method 1–3 的结构文件 ZIP |

本地网页不提供官网内置在线晶体结构查看、Distortion **"Generate"**
（按模式幅度生成畸变结构）或 **"Domains"**（畴列表）；这些功能不在当前逆向范围。
结构可通过导出的 CIF / IsoVIZ 文件交给 VESTA / IsoVIZ 查看。Python API 中保留的
`generate_distortion` / `generate_domains` 仅供脚本和验证使用，不表示当前网页产品
范围或与官网对应功能一致的承诺。

Method 1 点 OK 后的结果表字段与官网序参量方向页（`webpage_info/` 中 `a.` 开头的存档）一致：官网是一条条 radio，本地仍用可筛选、排序的表格展示同样的 token（`Irrep` / `OPD` / `Dir` / `SG` / `basis` / `origin` / `s` / `i` / `k-active`），并多一列 `idx` 供点选计算模式。

请勿修改 `ISODISTORT/resources/isobyu/` 内的文件。

---

## 2. 安装与环境

### 2.1 需要什么

| 项目 | 要求 |
| --- | --- |
| Python | **3.10 或更高**（当前产品化验证使用 3.12） |
| Windows | **必须安装 WSL**（`resources/isobyu/iso` 是 Linux ELF，通过 WSL 调用） |
| Linux 本机 | 可直接跑 `iso`，无需 WSL |
| ISOTROPY 套件 | 源码与本项目 wheel/sdist 中冻结的 Linux 版位于 `ISODISTORT/resources/isobyu/`；安装后用 doctor 核对 |

检查 Python：

```powershell
py -3.10 --version
```

检查 WSL（仅 Windows）：

```powershell
wsl --status
wsl --list --verbose
wsl -e uname -a
```

需要有一个**默认 Linux 发行版**且能启动。若未安装，请按 Microsoft 文档启用「适用于 Linux 的 Windows 子系统」，安装 Ubuntu 等发行版后重启。

### 2.2 克隆仓库

```powershell
git clone https://github.com/littleexplorer2/CRIS-Isodistort-Localization.git
cd CRIS-Isodistort-Localization
```

### 2.3 核对 ISOTROPY Suite（必需）

源码 checkout 与本项目构建的 wheel/sdist 会携带当前冻结的 **Linux** 版运行文件；安装后核对：

```text
ISODISTORT/resources/isobyu/
```

- 必需可执行文件：`iso`、`findsym`、`smodes`（`comsubs` 仍为可选）
- 数据库：全部 `data_*.txt`

若发行副本缺失或 `doctor` 报文件不完整，再从 [ISOTROPY Suite](https://iso.byu.edu/isotropy.php)
下载 Linux 版补齐。**不要**放入 Windows 可执行文件。根安装器不会联网下载、覆盖或修改
该只读目录；打包清单也明确排除其中的运行日志和误写输出。

### 2.4 安装 Python 依赖

在**仓库根目录**（不是只在 ISODISTORT 里）执行：

```powershell
cd <CRIS 根目录>
py -3.10 setup_cris.py install --project isodistort
```

若还要使用另外两个子项目，把 `--project isodistort` 改为 `--project all` 或直接省略。安装器会校验 Python 版本、创建或复用根 `.venv`、按声明的版本约束安装依赖、补齐 `output/tmp/`，再运行只读诊断。它不会自动下载第三方 ISOTROPY 文件。

开发与测试环境：

```powershell
py -3.10 setup_cris.py install --project isodistort --dev
```

安装后可随时单独复查环境：

```powershell
.\run_cris.ps1 setup_cris.py doctor --project isodistort
```

`doctor` 会分别检查 Python 依赖、配置导入、WSL 默认发行版、`iso`/`findsym`/`smodes`、执行权限、`data_*.txt`、`ISODATA` 与运行目录。Linux 下若原文件缺少执行位，封装会在私有暂存目录复制并设置副本权限，doctor 以 `WARN` 明示该执行方式；它只证明路径与权限就绪，不证明科学计算结果。有 `FAIL` 时不要开始正式计算；完整安装参数和离线 wheelhouse 流程见仓库根 [README.md](../README.md)。`--recreate` 只在 `.venv` 损坏或 Python 主/次版本改变时使用。

成功后可像旧版本一样直接启动网页：

```powershell
.\.venv\Scripts\python.exe ISODISTORT\scripts\main_web.py
```

这个入口会启动本地服务、自动打开浏览器，并可直接执行 Method 计算。安装器若检测到
OneDrive/WSL 的特定启动器限制，会自动保留原 stub 并安装经过 WSL 与 `sys.prefix`
验证的硬链接。若安装器明确报告硬链接修复不可用，再关闭服务并改用
`.\run_cris.ps1`；它仍加载同一实体 `.venv` 中的依赖。完整测试继续使用该兼容
启动器。若 PowerShell 阻止本地脚本，可对当前进程执行
`Set-ExecutionPolicy -Scope Process Bypass` 后重试。

### 2.5 安装后至少应看到

- 根目录 `setup_cris.py`、`VERSION` 与 `.venv/`；
- `ISODISTORT/scripts/main_web.py`；
- `ISODISTORT/resources/config/settings.yaml`；
- `ISODISTORT/frontend/web/index.html`；
- `ISODISTORT/resources/isobyu/iso`、`findsym`、`smodes` 以及若干 `data_*.txt`。

wheel 安装布局提供 `isodistort-web` 命令，直接调用同一网页入口；Python API 从
`backend.api` 导入。源码安装仍按上述 `scripts/main_web.py` 或 `run_cris.ps1` 启动。
wheel 不包含 CRIS 根部署器、实验数据、官网归档或 manual 审计脚本。

---

## 3. 如何启动

### 3.1 网页

```powershell
cd <CRIS 根目录>
.\.venv\Scripts\python.exe ISODISTORT\scripts\main_web.py
```

也可在已激活的虚拟环境中执行 `python ISODISTORT\scripts\main_web.py`。默认打开
`http://127.0.0.1:8000/`。端口见 `resources/config/settings.yaml` 的 `runtime.web_port`；被占用时
会自动顺延约 20 个端口，再不行则让系统分配。控制台会打印最终 URL；Windows Shell
已收到打开请求时还会打印 `Browser launch requested`，所有打开方式都失败时则明确
提示手动访问该 URL。

- 右上角 **"Stop"**：停止服务并释放端口。
- 每个标签页都有独立心跳；关闭**最后一个**本地标签后通常约 2 秒内自动停服、服务进程退出并释放端口。刷新页面有短暂宽限，多标签页中关闭其中一个不会误停服务。
- 若浏览器崩溃或关闭信标丢失，则由 `runtime.web_idle_timeout`（默认 60 秒）心跳超时兜底；关闭发生在 ZIP 计算期间时，会先完成/结束该请求再停服。要常驻就保留至少一个标签页。
- 网页是一个本机单用户共享会话，不是每标签独立会话。有状态 API 与结果选择在同一重入锁中原子执行；批量导出只在该锁内复核 `revision` 并建立不可变工作快照，耗时的逐候选计算在锁外进行，因此其他标签仍可读写实时会话，但不会把新旧状态混进已启动的 ZIP。每次修改会递增 `revision`，旧标签再提交或导出时会显式拒绝 stale context 并要求从最新状态重试。
- 顶栏 **"ISODISTORT"**：回到本页；**"SUITE"** / **"HELP"**：打开官网（需联网，计算本身不依赖）。

### 3.2 配置文件：`resources/config/settings.yaml`

源码运行时，相对路径相对 `resources/config/` 目录解析。安装 wheel 后，
`runtime.output_dir` / `temp_dir` 的项目内相对路径重基到当前用户的运行目录，
默认分别为 `<运行目录>/output` 与 `<运行目录>/output/tmp`；不向安装包目录写文件。
Windows 运行目录为 `%LOCALAPPDATA%/ISODISTORT`（未设置或不是绝对路径时使用
`~/AppData/Local/ISODISTORT`），Linux 为 `$XDG_STATE_HOME/isodistort`
（未设置或不是绝对路径时使用 `~/.local/state/isodistort`），macOS 为
`~/Library/Application Support/ISODISTORT`。

可在启动前设置绝对路径环境变量 `ISODISTORT_RUNTIME_ROOT`，让源码与 wheel 均将
相对运行路径重基到指定目录。例如 PowerShell 中
`$env:ISODISTORT_RUNTIME_ROOT = 'D:\ISODISTORT-state'` 对应
`D:\ISODISTORT-state\output` 与 `output\tmp`。YAML 中显式配置的绝对运行路径优先于
此覆盖；二进制与数据库的相对路径仍以 `resources/config/` 为基准。

| 键 | 含义 | 默认 |
| --- | --- | --- |
| `isobyu.bin_dir` / `data_dir` | `iso` 与数据库目录 | `../isobyu` |
| `isobyu.iso_bin` 等 | 可执行文件名 | `iso` / `findsym` / `smodes` |
| `defaults.symmetry_cartesian_tolerance_angstrom` | 上传结构空间群识别的 spglib/pymatgen `symprec`；笛卡尔距离，单位 Å，不是分数坐标容差 | `0.001` |
| `defaults.symmetry_angle_tolerance_degrees` | 空间群识别的晶格角容差（度） | `5.0` |
| `defaults.affine_exact_cartesian_tolerance_angstrom` | Method 3 精确有理仿射数据在 spglib 分类边界的笛卡尔容差下限（Å）；实际值还覆盖 `256·ε_machine·||lattice||₂` | `0.0000001` |
| `defaults.fractional_coordinate_tolerance` | spglib 标准化浮点结果恢复为有理分数时的无量纲残差；恢复后仍做精确 Seitz 重建 | `0.0000001` |
| `defaults.lattice_tolerance` | 晶格/整数矩阵比较的无量纲残差，不得传作 `symprec` | `0.00001` |
| `defaults.default_amplitude` | API 生成畸变时的默认幅度 | `1.0` |
| `defaults.eps` | 全局浮点容差 EPS（波矢约化等） | `0.00001`（与 `lattice_tolerance` 相同） |
| `runtime.web_port` | 网页首选端口 | `8000` |
| `runtime.web_idle_timeout` | 浏览器异常退出或关页信标丢失时的自动停服兜底秒数 | `60` |
| `runtime.temp_dir` | 网页上传等 Windows 侧暂存 | `../../output/tmp` |
| `runtime.output_dir` | Python API 与兼容下载端点的文件输出目录 | `../../output` |
| `runtime.timeout` | 普通 `iso` 调用超时（秒） | `60` |
| `runtime.generation_timeout` | Method 2 生成缺失 isotropy 子群库的超时 | `3600` |
| `runtime.method3_max_parametric_values` | Method 3 精确公度参数搜索值上限；超限时在 `list_irreps` 前报错，不截断 | `0`（不限制） |
| `runtime.method3_max_backend_queries` | Method 3 参数路径的 `list_irreps + list_subgroups` 查询总上限；超限时在任何 `list_subgroups` 前报错 | `0`（不限制） |
| `runtime.method3_max_quotient_order` | Method 3 exact affine/fixed-space 有限商的最大阶数；超限明确失败 | `512` |
| `runtime.method3_max_coupled_states` | coupled 稳定子交集动态规划的最大不同状态数；超限明确失败，`0` 为不限制 | `100000` |

容差采用两层误差模型：来自 CIF/实验或弛豫结构的对称性识别用 Å 制
`symmetry_cartesian_tolerance_angstrom`；Method 3 Stage-A 内部由整数/有理数精确
构造、只在第三方库边界转成浮点的数据用更严格的
`affine_exact_cartesian_tolerance_angstrom`。对分数位移 `df`，实际距离始终按
`||df @ lattice||₂` 计算，因此同一个分数误差会随晶胞长度和形状得到不同的
笛卡尔误差；`lattice_tolerance` / `eps` 则只用于无量纲矩阵或波矢残差。
这里的 `symprec` 语义遵循 spglib 官方定义（“Cartesian distance tolerance”）：
<https://spglib.readthedocs.io/en/stable/variable.html#symprec>；不能用晶格矩阵元素的
无量纲舍入阈值替代。

说明：网页 ZIP **不写入、不扫描** `output/`，也不读取 `webpage_info/` 或
`output_compare/` 的验证归档；网页上的筛选表 txt/csv 在浏览器里生成。Python API 可显式选择写入目录。
Method 2 勾选生成缺失子群库时，`iso` 把可复用的 `i*.iso` 写在 **WSL 短路径暂存目录**（`~/.id/tmp`，由封装自动创建），**不是** `resources/isobyu/`；网页可用 **Manage cached subgroup databases** 列出并批量删除。

### 3.3 用户需要配置 / 放置的路径

换机器或自定义目录时，按下面核对。跨项目关系见仓库根 [README.md](../README.md)。

| 项目 | 默认 | 你要做什么 |
| --- | --- | --- |
| **ISOTROPY 套件** | `ISODISTORT/resources/isobyu/` | 安装后核对 Linux 版 `iso` / `findsym` / `smodes` / 全部 `data_*.txt`（见 §2.3），仅在缺失时从官网下载补齐。**不要**为贴合测试改写二进制或数据库；目录本身只读约定见根 README |
| **改套件位置** | `resources/config/settings.yaml` → `isobyu.bin_dir` / `data_dir` | 仅当二进制不在默认 `../isobyu` 时修改（路径相对 `resources/config/`） |
| **母相 CIF** | 运行时选择 | 网页上传，或由 Python API 传入路径。示例：`experiment_data/EuAl4 Parent.cif`、`NdNiO2 own.cif`（只读，请复制后再改）。页头显示从该 CIF 读取 |
| **临时文件（Windows）** | `runtime.temp_dir` → `../../output/tmp` | 网页上传、写入 WSL 前的 Windows 侧临时目录；空间不够时可改到其它盘（仍相对 `resources/config/`） |
| **isotropy 子群库缓存** | WSL `~/.id/tmp/i*.iso` | 参数 k 点「Generate… if missing」生成；网页按钮可管理删除。**勿**手改 `resources/isobyu/` |
| **API 输出目录** | `runtime.output_dir` → `../../output` | Python API 与兼容下载端点的默认文件目录；网页 Distortion ZIP 走浏览器下载，不依赖此项 |
| **WSL / ISODATA** | 自动（`~/.id/data` → `isobyu`） | Windows 下封装会建短路径 `~/.id` 与 ISODATA 符号链接。一般**不必**在系统里手改 `ISODATA`；只需保证 `wsl` 可用 |
| **网页端口** | `runtime.web_port: 8000` | 端口冲突时改 YAML（非文件路径） |
| **抽检用 IsoVIZ** | 仓库根 `ISOViz.lnk` 或 `ISOVIZ`/`ISOVIZ_JAR` | 打开导出的 `data.isoviz`；安装与查找方式见 [ISOVIZ_INPUT/README.md](../ISOVIZ_INPUT/README.md) |
| **抽检用 VESTA** | 可选根目录 `VESTA.lnk` | 手工打开 `subgroup.cif`；程序不读取该路径做计算（下载与作用见下） |

**不需要用户改路径的部分：** `backend/`、`features/`、`frontend/` 源码内的导入、WSL 侧自动 staging、网页 ZIP 内存打包。

#### VESTA（打开 CIF）

[VESTA](https://jp-minerals.org/vesta/en/) 是常用的**三维晶体结构可视化**软件（Windows / macOS / Linux）。在本项目中用于：

- 打开 Distortion ZIP 里的 **`subgroup.cif`**，查看子群结构、原子坐标与晶胞；
- 验收时确认「CIF 能用 VESTA 正常打开」（见本目录 `agent.md` / 下文已知限制）。

下载：[VESTA 下载页](https://jp-minerals.org/vesta/en/download.html)（64 位 Windows 选 `VESTA-win64.zip`，解压后运行其中的 `VESTA.exe`）。  
可选：在 **CRIS 仓库根目录** 放快捷方式 `VESTA.lnk` 指向该 `VESTA.exe`（已 gitignore）。计算与导出走 Python/`iso`，**不调用** VESTA。

#### IsoVIZ（打开 `.isoviz`）

见 [ISOVIZ_INPUT/README.md](../ISOVIZ_INPUT/README.md)（Java + ISOTROPY Suite 中的 IsoVIZ）。

---

## 4. Parent CIF（母相）加载后显示什么

界面区块标题：**"Parent CIF"**。选择 `.cif` / `.CIF`，点 **"Load"**。必须先 Load 成功，再跑任何 Method。

加载成功后，状态区会显示类似官网的页头（**从当前 CIF 读取**，不是写死某几个结构）。格式与官网 Search 页一致，例如 EuAl4：

```text
Space Group: 139 I4/mmm D4h-17
Lattice parameters: a= …, b= …, c= …, alpha= …, beta= …, gamma= …
Default space-group preferences: …
Eu1 2a (0,0,0),
Al1 4d (0,1/2,1/4),
Al2 4e (0,0,z), z= 0.38000
```

NdNiO2 等其它 CIF 会显示该文件自己的空间群、晶胞与 `_atom_site_label` / 位点顺序（如 `O 2f`、`ND 1d`、`NI 1a`）。若输入使用同一空间群的等价非标准 basis/origin，程序会把晶格、全部原子和每个独立物理 Wyckoff 轨道一起规范到 conventional standard setting，再用这一个内部结构生成页头、搜索结果和模式；不会把标准代表坐标单独套到未变换的输入结构。无法保持空间群、轨道身份或多重度时会明确拒绝载入。

若 FINDSYM 运行失败、会话标为 `parent_orbit_standardization.status=unavailable`，
载入结果仅可用于诊断；Method 1/3 搜索、Method 2、路径选择、共享模式计算及子群
导出会拒绝继续，包括不补算模式的 CIF-only 下载。
修复 FINDSYM/WSL 环境后重新载入母相，不要把尚未核验的上传坐标与 ISO 标准 k/B 混用。

含义（白话）：

| 行 | 意思 |
| --- | --- |
| **Space Group** | 空间群编号 + 国际符号 + Schoenflies 符号（如 `D4h-17`） |
| **Lattice parameters** | 晶胞边长与夹角，**五位小数** |
| **Default space-group preferences** | 本地固定的国际标准取位约定（见文末 Preferences 表，只读） |
| **Wyckoff 行** | 来自 CIF 不对称位点：标签（或 type_symbol）、多重度+字母、坐标；自由坐标写成 `z= 0.38000` |

实验示例（只读）：`experiment_data/EuAl4 Parent.cif`、`experiment_data/NdNiO2 own.cif`。官网页面对照见 `webpage_info/<母相名>/`（序号 1–6、2.5、3.5、a.）。Method 2 的 k 点随母相而变（EuAl4：`LD … g=1/6`；NdNiO2：`Y … a=1/3`），验证批次见 [开发计划 §3.2](docs/DEVELOPMENT_PLAN.md#32-method-2当前保存产物已复验)。

---

## 5. Types of distortions to be considered

区块标题：**"Types of distortions to be considered"**。勾选后**必须点 "Change"** 才会生效（页上有英文提示 *"Important: You must click on Change…"*）。只勾选不点 Change，后面的 Method 仍用旧设置。

| 界面标签 | 白话 |
| --- | --- |
| **"Strain"** | 是否考虑晶格应变（单个复选框） |
| **"Displacive"** | 原子位移；每行有 **all / none / 各物种** 复选框 |
| **"Occupational"** | 占位（有序）畸变；同样 all / none / 物种 |
| **"Magnetic"** | 磁相关；同样 all / none / 物种 |
| **"Rotational"** | 旋转模式；同样 all / none / 物种 |

- 复选框互不联动；点 Change 时按 **all > none > 具体物种** 解释。  
- 默认与官网第 2 页一致：Strain + Displacive 各物种勾选。
- Displacive 的逐物种作用域按独立物理 Wyckoff 轨道传给模式后端；同一 Wyckoff 字母下的多个独立轨道不会合并成一个位点。
- 点 Change 后，顶部蓝色状态栏直接使用该次提交返回的最新状态重绘；Method 1 下拉若仍在后台加载，结束时也只会重绘最新状态，不会把旧 Types 文本覆盖回来。

---

## 6. Method 1–4：每个选项与 "OK" 做什么

四个 Method 面板只负责**调参数并计算结果表**。下载一律到页面底部 **Distortion**。

结果表通用能力（Method 1–4）：

- **"Filter"**：对相关列做不区分大小写的子串匹配。  
- **"Show filtered rows only"**：只显示命中行（下载时仍导出**全部命中行**，不受此勾选限制）。  
- 表头 **"▲" / "▼"**：按该列升序/降序。  
- 点击一行（Method 1–3）：按该子群计算模式基矢（显示在 Method 2 下方 modes 区）。按子群 **`idx`** 计算，不要按屏幕上的「第几行」理解。

### 6.1 Method 1: Search over all special k points

在**所有特殊 k 点**上枚举各向同性子群。

| 选项 | 作用 |
| --- | --- |
| **"Crystal system(s):"** | 多选晶系（triclinic / monoclinic / orthorhombic / tetragonal / trigonal / hexagonal / cubic）；逻辑 OR；全不选 = 不过滤 |
| **"Space-group symmetry:"** | 可达子群空间群下拉；空 = 不过滤。选项随母相与 Types 变化 |
| **"Conventional lattice:"** / **"Primitive lattice:"** | 互斥过滤格点类（选一个会清空另一个），选中后只保留该类候选。按母相点群轨道与 GL(3,Z) 等价关系分组；Primitive 使用子群原胞格。对带心母相用子群对称性 + spglib 标准胞归约消除任意 unimodular basis 漂移。若已存官网 Search 页与加载的 CIF 完全匹配，且所有选项的候选数和格点类均验证一致，使用存档里的显示基矢；其他结构使用实时计算的代表元。分类与搜索结果始终实时计算，不按母相/IR 硬编码官网快照。 |
| **"Maximal subgroups only"** | 只保留极大子群 |
| **"OK"** | 开始枚举；首次可能数秒到数十秒 |

结果表列与官网 Method 1 序参量方向页（radio 一行里的字段）对齐，但以**表格**展示，而不是一条一条的 radio：`idx` / `Irrep` / `OPD` / `Dir` / `SG` / `basis` / `origin` / `s` / `i` / `k-active`。`idx` 是本地点选/导出用的子群序号，官网可见文本里没有这一列。下载的 txt/csv 使用相同列。勾选 lattice strain 时会保留 smodes 没有的纯应变 irrep（例如 I4/mmm 的 GM4+ Fmmm）。`basis` / `origin` 直接来自 iso；`k-active` 由通用算法从方向与 k 坐标生成（如 I4/mmm 的 M 点为 `(1,1,1)`）。不查按 (irrep, OPD) 硬编码的官网 token 表。

### 6.2 Method 2: General method - search over specific k points

针对你指定的单个 **k 点**枚举单 IR 子群。多 IR 叠加的当前限制见下方注意事项。

| 选项 | 作用 |
| --- | --- |
| **"Specify k point:"** | 从母相允许的 k 点列表选择 |
| 参数 k 点数值框 | 如 LD 的 a/b/g，必须填齐再 OK |
| **"Generate isotropy subgroups database if missing"** | 若当前参数 k 点在本地 ISOTROPY `data_*.txt` **没有**现成子群列表，则由本机 `iso` **现场生成**并写入 WSL 暂存目录中的 `i*.iso` 缓存（首次可能数分钟～数小时，之后 OK 直接复用；超时见 `generation_timeout`）。这与已删除的 Distortion **"Generate"**（按幅度生成畸变结构）**不是同一功能** |
| **"Manage cached subgroup databases"** | 列出已生成的 `i*.iso` 缓存（母相 SG / IR / kparam），可勾选后批量删除以释放磁盘或强制下次重新生成 |
| **"Change number of superposed IRs:"** + **"Change"** | 当前本地 Method 2 只支持 **1** 个 IR 组分。输入大于 1 时会明确拒绝，不会把各组单 IR 子群表错误拼接为叠加结果 |
| **"OK"** | 枚举当前单个 k 点下全部 IR 的单 IR 子群 |

结果表列大致含：`idx` / `SG` / `k` / `Irrep` / `OPD` / `s` / `i`。

注意：

- **多 IR 叠加当前 fail-closed**：真正的叠加结果需要枚举所选 IR/OPD 稳定子的精确交集，不等于合并多个 k 点的单 IR 列表。本地网页尚未接入这条 coupled-IR 生产链，因此 `nsup > 1`（包括两个组分选择同一 k 点）会在计算前拒绝，不会生成候选、模式或导出文件。需要该功能时请暂用官方 ISODISTORT 网页。
- **`# of independent incommensurate modulations`（nmod）**：网页可编辑 0–3，须点 **Change** 才生效。**0** = 三维锁定：凡折叠进子群晶胞的母相 k（含 GM/M 等次级 IR）都保留。**1、2、3** 在 Method 2 里结果相同：只保留**当前这一个**主波矢 q 的谐波，外加 Γ。Method 2 一次只选一个 k，不会因为把 n 调到 2 或 3 就再引入第二条独立调制。
- **参数 k 点**（如 LD/DT）：可枚举子群（并可勾选 **Generate isotropy subgroups database if missing**）；位移模式由本地 **smodes + 子群对称性**（公度锁定 / (3+d)）计算，不再默认空表。选中该 k 点会提示走该引擎。
- 子群数据库缺失时：网页提供「本地生成 / 去官网」。

### 6.3 Method 3: Search over arbitrary k points for specified space group and lattice

按指定空间群/点群与格子搜索。

| 选项 | 作用 |
| --- | --- |
| **"Select either space group symmetry:"** | 选 230 个空间群之一 |
| **"or point group (crystal class):"** | 或选 32 个点群之一 |
| **"Specify a real-space sublattice… Default / P / A / B / C / I / F / R centering"** | **direct** 实空间子格。Default 在选空间群时使用目标 SG 的默认 Bravais centering、只选点群时使用 P；显式 P/A/B/C/I/F/R 分别把 conventional basis 精确换算为 primitive translation lattice。输入必须确为母晶格子格 |
| **"Specify a primitive reciprocal-space superlattice"** | **reciprocal**：**本地不支持**，请用 direct |
| **"Choose a representative basis:"** | 3×3 基矢：`a'` / `b'` / `c'` 相对母相 a,b,c 的系数（可填分数）。官网会同时检查该 basis 与目标空间群默认 centering 是否构成母晶格的子格；轴序、负号和带心都重要，不能因为体积相同就随意换成单位矩阵。非恒等整数超胞会推断公度参数 k（如 `(0,0,6)` → LD `g=1/6`）并枚举；若子群库缺失，勾选 Method 2 的 **Generate isotropy subgroups database if missing** |
| **"OK"** | 空间群查询得到 single-IR 与经 exact fixed-space/稳定子交集证明的 coupled affine embeddings；未覆盖范围见下方限制 |

网页/API 传入的分数字符串全程保留为精确有理数；恒等与非恒等 basis 都只匹配母相点群轨道下的同一 lattice class，不会再把任意更大超胞当作命中。对一参数线，引擎可解析 `-a+1` / `1-a` / `a+1/2` 等仿射坐标，在母相完整 reciprocal k-star 上精确求出 `0 ≤ a < 1` 的全部公度解，并用母相中心化格子的倒格矢等价关系排除已由特殊 k 枚举的端点；同一物理 k-star 的参数路线会按精确倒格矢等价关系去重。这不等于支持 Method 3 的 reciprocal-sublattice 输入。

官网首屏候选的科学身份是 `(SG,basis,origin,s,i)` embedding；表中 `basis` / `origin` 优先显示 iso 的官网式精确分数原文，`known routes` 只是本地已知 k/IR/OPD 来源的诊断信息，不能把 route 数误当作官网候选数。single-IR embedding 按各 route 的 Types 活性过滤；coupled embedding 则在所选 `strain ⊕ displacive` 表示中直接做 exact fixed-space 判定。Displacive 的逐物种作用域会进入该表示；只移动部分物种时，其共同位移是相对光学自由度，不会被误删为全晶体刚体平移。参数 k 在活性检查与选中后的位移模式中先做官网参数尺度 → iso 内部尺度转换，随后走 smodes/(3+d) 完整模式（nmod）。

官网对照不能只比较 basis/origin 字符串。生产代码和验证工具都会在母相分数坐标中用精确有理数重建 Seitz 操作，按子群 primitive translation lattice 取模。空间群查询的阶段 A 按有限平移商 `N_G(T_s)/T_s`、Seitz 闭包和 cocycle 条件枚举 affine lifts；阶段 B 在所选 `strain ⊕ displacive` 表示中用 exact character/fixed-space 点态稳定子判据排除物理不可达项。对没有单-IR route 的可达 embedding，所有候选单-IR 稳定子连同平移陪集被提升到同一有限商，只有其精确交集等于目标 `H` 才成为可选择的 `exact_fixed_space` 行；母群仿射共轭只在 Method 3 首屏保留一个代表元，物理畴枚举不属于当前本地网页版范围。选中后以 nmod=0 计算该子群全部折叠 k 的完整位移 fixed space，不把任选的一对 IR 冒充唯一 primary route。当前 40 组权威官网查询的 77 条 embedding 全部匹配：13 组逐字段一致、27 组精确 affine 等价、0 差异、0 错误；此前缺少的 16 条 coupled-only 已接入。只选 point group 时仍使用 single-IR 子集；项目附带的 ISO 9.6.1 没有直接枚举首屏 embedding 的命令，`DISPLAY DIRECTION` / `DISPLAY ISOTROPY COUPLED` 是已知 embedding 后的 route 工具，不能用它们或 `SHOW DOMAIN` 代替完整第一阶段。

Method 3 的空间群下拉已对 `webpage_info/EuAl4 Parent.cif/2. ISODISTORT_ search.html` 做 230 项全量回归：编号、旧版 IT Hermann–Mauguin 简写和 Schoenflies 标号必须逐项一致。

### 6.4 Method 4: Mode decomposition of a distorted structure

| 选项 | 作用 |
| --- | --- |
| **"Upload distorted structure from CIF file:"** | 上传**已经畸变**的女儿相 CIF |
| **"Atom-matching method:"** | `nearest-site` 使用按物种的全局最小笛卡尔距离分配；`robust` 还会拒绝超过阈值的匹配 |
| **"Robust distance threshold (angstrom):"** | robust 匹配的物理距离上限，单位 Å；不是随晶胞变化的分数坐标距离 |
| **"Known origin shift..."** | 可选的女儿相分数坐标原点平移 `(x,y,z)`；留空表示 0/未指定 |
| **"OK"** | 输出**全部**位移模式的 `As`、`Ap`、原始拟合系数和 `normfactor`（可 Filter / 排序）；勾选 strain 且已有精确目标子群时，还输出官网顺序的 symmetry-adapted strain mode 幅度、CIF raw-coordinate sum、应用工程应变 `q=(E11,E22,E33,2E23,2E13,2E12)`，以及单位为 Å 的 RMS/max residual |

这是**分解**，不是子群列表。**不能**作为 Distortion ZIP 的结构来源，但可以在 Distortion 下载该幅度表的 txt/csv。

当前源码的 displacive 生成→分解闭环能恢复无噪声系数并报告带噪声 residual；原子匹配使用真实晶格下的全局一一映射，`As` 按女儿原胞内笛卡尔模式范数归一化，`Ap=As/sqrt(s)`。均匀应变使用官网的母胞基坐标约定：母胞行晶格为 `P`、子群 basis 为 `B`，程序由女儿相 metric 解对称正定 `M=I+E`，满足 `M (P P^T) M = B^-1 (D D^T) B^-T`，并以 `B M P` 重建女儿胞。这里的工程 Voigt `q` 对非正交母胞不是 Cartesian 张量分量。

官网 CIF 的 `_iso_strain_value` 和 strain matrix 使用 `q_raw`，IsoVIZ 与 Method 4 实际晶格使用 `q_unit=normfactor*q_raw`；模式幅度 `a_i` 因而分别给出 `Σa_i q_raw_i` 与应用应变 `Σa_i q_unit_i`，两个六维量不会混用。规范模式来自本地 ISO 的 rank-[12] `DISPLAY DISTORTION`，再由目标 embedding 的 `DISPLAY DIRECTION` 约束；程序还会用实际嵌入点群操作和实际母胞 metric 独立计算固定子空间并核对完整 span。缺少或不通过这条证据链时会明确失败，不会按晶系硬编码模式或猜 `GM` 标签。

Method 4 的冻结输入、官网对照、机器报告及对应源码签名见
[验证报告](docs/BUGFIX_VALIDATION_REPORT.md)。双母相测试不能外推为 4310 或跨晶系科研级
验收；历史报告中的“current”只指其记录的源码快照，不能用来认证修改后的实现。

---

## 7. Distortion（统一下载）

区块标题：**"Distortion"**。没有 Generate / Domains，也没有模式幅度输入框。

| 控件 | 作用 |
| --- | --- |
| **"Export source (exactly one Method):"** | 下拉选 **一个** Method 1/2/3/4（默认 Method 2）。不能多选 |
| **"Download filtered (txt)"** / **"Download filtered (csv)"** | 按所选 Method **当前** Filter 与排序，下载全部命中行。须先对该 Method 点过 OK |
| **"Download formats:"** | 可多选，仅用于下面的 ZIP： |
| → **"CIF file"** | 子群 CIF（ISODISTORT 6.12 布局：子群设置、不对称单元、`iso_*` 循环；一般为零振幅） |
| → **"Save interactive distortion"** | `.isoviz`（给 IsoVIZ 用） |
| → **"Complete modes details"** | 完整模式详情 `.txt` |
| → **"TOPAS.STR"** | TOPAS 结构文件 |
| **"Download all (ZIP)"** | 打包 Method **1 / 2 / 3** 当前筛选命中的子群文件（无筛选则全部），**只含勾选格式**，不扫描 `output/`、`webpage_info/` 或 `output_compare/`；ZIP 从启动时的会话快照生成。选 Method 4 再点 ZIP 会提示改用表格下载。勾选 isoviz / modes / topas 时，会对子群补跑 Method 2（可能很慢）；界面显示当前候选、完成数和估算剩余时间，可请求在当前候选结束后安全取消。已完成候选仍发布，未运行候选写入状态报告。**参数 k 点**（如 LD）走 smodes/(3+d) 完整模式（nmod，默认 0）。仅 CIF 或 URL 带 `compute_modes=0` 时可跳过补算 |

生产 core 会把已验证的 ISO microscopic 列映射到精确 emitted child frame，并把稳定原子
ID、父子原子映射、查询的 `(parent/child SG, B, q, primary IR/OPD)`、精确 `VECTOR`
方向和未混合来源列一起交给 `DisplaciveExportData`。官方 microscopic 表只覆盖部分代表点
时，程序只有在公共点域满秩、变基唯一且 BUSH 全域延拓通过秩、条件数、残差和摘要复核后
才发布。模式身份未解析的 Method 2 候选不会进入四格式 writer，也不会阻止同批其他已验证候选发布；ZIP/磁盘批次根部写出 `export_candidate_status.json` 与
`export_candidate_status.txt`，记录候选完整身份、物理轨道、实际拒绝判据和查询上下文。身份以外的计算或 writer 异常仍使整批 fail-closed。

勾选 strain 时，CIF、IsoVIZ 和 Complete modes 共用同一个应变数据契约和模式顺序。CIF 写完整的 mode number/label/value/norm/raw matrix，IsoVIZ 写 `q_unit` 且 strain `maxamp=0.1`，Complete modes 同时说明 `q_raw`、`q_unit`、`normfactor`、幅度和两种六维和。官网 TOPAS 文件不导出 strain-mode 精修参数；本地同样只写由 `B M(Σa_i q_unit_i) P` 得到的固定实际晶胞。若参考子胞不能与 `B P` 严格对应，TOPAS 导出会报错而不编造晶格关系。

共享合同已覆盖非零 strain 与已验证位移同时存在的情况：依据
[ISODISTORT 官方说明](https://iso.byu.edu/isodistorthelp.php)，位移模式的幅度、
归一化因子和 lattice-coordinate 模式向量不随 strain 改变；程序因此保留已验证的
位移分数坐标，并只把最终晶格重建为 `B M P`。生产 core 和四个 writer 使用同一合同；
跨晶系、混合/部分占位及尚未建立同等级身份链的模式类型仍需另行验收。

symmetry-adapted strain 的定向回归及 Method 4 冻结矩阵证据统一见
[验证报告](docs/BUGFIX_VALIDATION_REPORT.md)，以其中的源码签名和范围为准。

压缩包名形如 `isodistort_methodN.zip`。Method 1/2 解压后直接是短子群目录；若有身份未解析或用户取消后尚未运行的候选，根部还会有上述两份状态报告。
Method 3 先有一个稳定案例目录：API 可用 `stable_case_id` 指定，网页则根据
当次导出的候选集合生成通用的 `M3-<摘要>`，不包含特定晶体或样例名称。

```text
GM1+_P1_SG139/                   # Method 1：IR_OPD_SG<number>
  subgroup.cif
  data.isoviz
  Complete modes details.txt
  topas.str
LD1_C1/                          # Method 2：IR_OPD
  subgroup.cif
  data.isoviz
  Complete modes details.txt     # 官网为 HTML 页；本地用 .txt
  topas.str
export_candidate_status.json     # 仅部分成功/安全取消时出现
export_candidate_status.txt
M3-a1b2c3d4e5/                   # Method 3：自动生成的稳定案例目录（示意）
  C01_SG139/
    subgroup.cif
    data.isoviz
    Complete modes details.txt
    topas.str
```

`+` / `-` 与 `4D1` 等 IR/OPD token 会保留，Windows 非法字符会替换为下划线。
若不同候选得到同一短名，程序会追加候选完整身份的短摘要；目录导出遇到已有同名
目录也会另建摘要/序号目录，不覆盖原文件。四种内部文件名保持不变。

单行点选后出现的模式表在 Method 2 区域，仅供查看。

### Space-Group Preferences（只读）

页底 **"Space-Group Preferences"** 固定为国际标准取位（Monoclinic axes a(b)c、cell choice 1、Orthorhombic abc、Trigonal hexagonal、Origin choice 2、Superspace standard）。本地 `iso` 不能改这些选项；Method 1 表里的 `basis`/`origin` 对已收录母相（如 I4/mmm #139）按官网序参量页显示。

---

## 8. Python API 快速示例

在已安装依赖、且 `PYTHONPATH` 含 `ISODISTORT/`（或从该目录运行）时：

```python
from backend.api import IsoDistort

iso = IsoDistort()
iso.load_structure(r"C:\path\to\parent.cif")
iso.set_distortion_scope({
    "displacive": ["*"],
    "occupational": [],
    "strain": [],
    "magnetic": [],
    "rotational": [],
})
iso.set_distortion_types(["strain", "displacive"])

m1 = iso.search_method_1(crystal_system="tetragonal")
m1_candidates = [item.subgroup for item in m1]
selected = iso.search_method_2(
    m1_candidates[0].index,
    candidates=m1_candidates,
)
iso.export_subgroups(
    "out_batch",
    formats=["cif", "isoviz", "modes", "topas"],
    subgroups=m1_candidates,
)
```

源码按四个部分组织：`backend/`（通用底层）、`features/`（按功能划分的中端）、`frontend/`（网页支撑）、`resources/` 等共享文件。Python 侧只有这一套包名，公开导入路径是 `from backend.api import IsoDistort`；历史上并存的 `isocore` / `isodistort` 兼容命名空间已随结构重构移除，不再提供别名。`features/` 内部跨包导入使用绝对路径，同包内才用相对导入。

当程序同时保留多个 Method 的结果表时，必须把候选池显式传给 `search_method_2(..., candidates=...)`，并把导出池显式传给 `export_subgroups(..., subgroups=...)`。各 Method 的显示索引都会从 0 开始，不能把索引当成跨表唯一 ID，也不要从界面层直接改 `iso.subgroups`、`mode_displacements` 或私有 `iso._iso`。单一路径脚本仍可省略 `candidates`，沿用最近一次 API 产生的默认候选池。

`generate_distortion` / `generate_mixed_distortion` / `generate_domains` 仍在 API 中，
仅供脚本和内部验证；网页不调用，也不属于当前网页逆向范围。

---

## 9. 已知限制（影响日常使用）

1. **Windows 必须经 WSL** 调用 Linux 版 `iso`。  
   Linux 原生运行会在 `/tmp/id-*` 建立每个封装实例独立的私有短暂存目录。缺少执行位的二进制会复制到此目录并仅对副本设置 `0700`；原套件文件保持只读。原文件已经可执行时直接使用原文件。
2. **symmetry-adapted 应变模式采用 fail-closed**：勾选 strain 后，程序必须取得 ISO rank-[12] 宏观基、目标 embedding 的 invariant directions，并通过实际 metric fixed-space 的整空间核验，才会写 CIF / IsoVIZ / Complete modes，并据此计算 TOPAS 的固定实际晶胞。ISO 输出缺失、解析不完整、模式数不符或 span 核验失败都会中止该候选导出；不会退回无标签坐标轴基或猜 `GM` 标签。TOPAS 按官网行为不提供 strain-mode 精修参数。
3. **参数 k 点（LD/DT 等）**：可枚举子群（+ Generate DB）；位移模式由 **smodes + 子群恒等表示** 计算（三维锁定 / (3+d) 谐波，由 nmod 控制）。网页不再把这类结果默认降级为仅 CIF。
4. **nmod / (3+d) superspace**：本地可编辑 nmod（0–3）。0 = 公度锁定，保留全部折叠 k。n≥1 = 只保留 Method 2 所选那一个 q 的谐波，再加上 Γ；1、2、3 不会增加第二条独立调制。标签格式对齐官网 `Parent[k]IR(opd)[Site:letter:dsp]siteIR(comp)`。IsoVIZ 对非 Γ 的 k 保留 `[kx,ky,kz]IR[...]` 前缀。二维 IR 的第二个实分量用母相平移（含心平移）做相位正交；旋转星臂分别保存展示坐标与实际相位坐标，自共轭特殊 k 的简并基由含平移的母相 little group 补齐。折叠商按母相中心化倒格矢计算（I 心整数奇偶余类不会误合并），搜索范围跟着子群基矢走，不限于 5 个母相单胞。最终模式空间在笛卡尔坐标中用修正 Gram–Schmidt 求精确数值秩，不再用固定点积阈值误删或重复计数。
5. **Method 3**：reciprocal-sublattice 输入不支持。空间群 direct 查询已接入 closed affine lifts、`strain ⊕ displacive` exact fixed-space 可达性、single-IR 稳定子精确交集见证和完整 fixed-space 位移模式，因此当前双母相 40 组/77 条权威 embedding 无欠枚举。搜索单-IR 稳定子的 k 域仍限于特殊 k 与一参数公度线；任意/多参数 k、point-group-only affine 枚举，以及 rotational/occupational/magnetic 的同等级 Stage-B 表示尚未完成。`exact_fixed_space` 行不声明某一组 IR 是唯一 primary COPL 分解，但可计算目标子群的完整折叠-k位移模式；不会用 `SHOW DOMAIN` 畴数或样例硬编码填满。
6. **Method 4**：只有勾选 strain 且会话中保留精确目标子群 embedding 时，才报告官网顺序的 symmetry-adapted strain mode 幅度；低层 API 未提供 canonical basis 时会把模式状态标为 unresolved，同时仍可报告由 metric 解出的应用 `q`。一般非自共轭 k 必须使用已证明的完整子胞位移列；只有一个母胞实列时明确拒绝，不能由 cosine 猜测 sine/OPD。显式 Γ/自共轭实列均须验证母相平移字符、中心化位点协变性与子胞周期性；未声明 k 不强制为 Γ，但位移场仍须在请求子胞中周期。occupancy / magnetic / rotational 分解尚未验收。当前实测结果以 [验证报告](docs/BUGFIX_VALIDATION_REPORT.md) 为准。
7. **magnetic**：带 `m` 前缀的 IR 默认不进入流程。
8. **occupational**：本地为 ±1 占据近似，校验失败会标明。
9. **rotational-only Types**：当前用 smodes 的位移活性作为 rotational route 的近似筛选，尚无独立的刚性转动/轴矢量模式生成器；不应把 rotational-only 结果称为与官网严格等价。
10. **官网内置在线结构查看、Distortion Generate / Domains**：不属于当前本地
    网页逆向范围；使用导出的 CIF / IsoVIZ 与外部 VESTA / IsoVIZ 查看结构。
11. **界面仅英语**。
12. **导出验收（见本目录 `agent.md`）**：网页委托 `IsoDistort.export_subgroups_zip` / `_collect_export_specs` 完成导出。ZIP 内每子群文件夹含 `subgroup.cif` / `data.isoviz` / `topas.str` / `Complete modes details.txt`。
    - **modes `.txt`**：完整写入本地计算结果即可，**不要求**与官网 HTML 逐字节一致。  
    - **CIF / isoviz / TOPAS**：内容与格式尽量靠官网；`.cif` 须能用 **[VESTA](https://jp-minerals.org/vesta/en/)**（[下载](https://jp-minerals.org/vesta/en/download.html)）打开，`.isoviz` 须能用 **ISOViz**（ISOTROPY IsoVIZ）打开。根目录可放 `VESTA.lnk` / `ISOViz.lnk` 便于抽检。  
    - TOPAS / IsoVIZ 位移向量按**原胞笛卡尔 Σ‖Δr‖²=1** 归一化（惯用胞求和除以 centering 重数）。官网 `As` 的物理最大位移为 `dmax=|As|·normfactor·max_i‖u_i B‖`，因此滑条/TOPAS 对称界为 `maxamp=1/dmax(As=1)`；不能按点阵类型硬编码 `√2`/`2`。
    - CIF、TOPAS 和 ISOVIZ 共用目标子群的原点与位点轨道；不能再用零振幅母相的对称性合并子群位点。ISOVIZ 会按原始 CIF 位点顺序写类型和子位点，并为边界周期像写对应模式向量。其他允许的表示差异包括 symop 顺序、周期像及子位点排序、basis 点群等价代表元。CIF、ISOVIZ 和 Complete modes 的应变模式标签、顺序、`q_raw`、`q_unit`、`normfactor` 和幅度由同一契约生成；TOPAS 只写固定实际晶胞，不包含 strain-mode 参数。VALIDATE 默认语义比较，`--strict` 仅排版调试。
13. **位点模式标号采用 fail-closed**：程序先用目标 embedding 的精确
    `DISPLAY DIRECTION` 构造 primitive-integer `VALUE DIRECTION VECTOR,...` 查询，再从
    oriented `DISPLAY DISTORTION` 表读取 `(global irrep, Wyckoff, site irrep, component)`
    身份。命名 OPD token 只作为 primary 子群身份门禁，不能替代 exact invariant direction。
    只有同一 `(global irrep, Wyckoff)` 的列空间与 BUSH 空间在秩、独立性和主角度上吻合，
    且每列来源完整、唯一时，标签才标为 `verified`。若官网只打印稀疏点域，程序会在已证明
    的公共点域求唯一 `R_common T=C_common`，再以 `C_full=R_full T` 延拓到完整 BUSH 点域；
    该延拓有独立证据，不会冒充官网直接打印的行。缺块、非唯一/病态变基、残差超限、
    未证明的非 Γ 平移相位和纯数值 fallback 均保持 `unresolved`。
14. **位移导出采用已验证合同并 fail-closed**：CIF、IsoVIZ、Complete modes 和 TOPAS
    四个 writer 已统一接入 `DisplaciveExportData`。共享模型验证 Seitz frame、子群轨道、
    自由坐标、max-component-one 模式尺度、primitive norm 与未混合 microscopic 来源；
    每列模式还必须绑定精确有理 k、物理 `orbit_id`、已解析 frame/atom-order 和唯一来源列。
    父相、参考子胞与最终结构使用不可变快照；显式父子原子映射逐项验证
    `x_parent=x_child B+q`、化学身份、物理轨道和完整平移陪集。IR、OPD、母/子空间群、
    `B,q`、embedding ID、primary `VECTOR` 与由平移格推导的 `s` 也会在每次 writer 调用前
    重新核对，不能由可变子群字段改写 `As→Ap`、模式标签或查询 setting。生产 core 只接受
    与 canonical ISO 列逐项一致的 lifted arrays，匿名或混合来源列会拒绝。合同当前只接受
    单物种满占位点；混合/部分占位在四种格式都有无损表示之前同样拒绝。门禁只检查该候选
    实际写出的位移列，不能由同一缓存中未写出的其它列误伤。Method 2 的身份未解析候选被
    单独排除并进入根级状态报告；其它准备或 writer 异常仍整批拒绝。磁盘与 ZIP 都先把
    本次允许发布的全部候选和状态报告渲染完；新/空目标直接发布，已有非空目标写入带内容摘要的
    `.isodistort-batch-v1-<digest>.ready` 目录及逐文件大小/SHA-256 manifest，并以进程锁和
    文件锁保证并发批次只能整批成功或整批拒绝。成功但确实没有模式仍与失败区分。
    IsoVIZ 的 parent-atom 类型只来自
    权威物理轨道映射，同元素的不同轨道不会合并或按显示标签猜测。共享合同已定向验证
    非零 strain + displacive 的组合；已覆盖的真实非 P1 结果与最终全量状态见验证报告。
15. **Method 1 候选 setting**：本地 iso 与官网可能给出同一嵌入子群的不同晶胞基矢表示。正式审计先按完整 IR/OPD/方向/SG/s/i/k-active 身份配对；同一候选的代表胞差异须用精确整数幺模变换证明相同子格，并核验完整 Seitz 群相同，再统一表示比较结构。仅晶格相同或属于母群共轭不能证明同一候选；不同 orientation/domain 不会因此合并。短目录名和字面 basis 不是科学身份。
16. **Method 2 短目录名**：当前导出和对照目录使用 `<IR>_<OPD>`（如 `LD5_C4`）。短名只用于定位；若目录名与 CIF 内 OPD 行不一致，按目录名硬配对会得到假差异，审计仍须以 CIF 内完整身份和 `# … k-active=` 行为准。
17. **子群 CIF 的 ASU 行数 / 部分晶胞边长**：子群设定、原点与不对称单元取位与官网不完全相同时，ASU 行数或 a/b/c 可能不同；应检查展开后的原子数、化学计量和结构匹配。分数基矢必须属于实际上传母相的平移格，扩胞需包含新晶胞内全部母胞平移，并保持原子数/物种与体积比守恒；非奇异但不满足子晶格关系的矩阵会拒绝。导出时原点须使目标子群对称操作映射到同种原子。勿按文件文本差异硬编码答案。

更细的差异用 `output_compare/<母相>/{官网,现有网页版交互}/Method1|2` 做 diff；上表是用户最常撞到的几条。

---

## 10. 测试如何运行

从仓库根目录：

```powershell
cd <CRIS 根目录>
py -3.10 setup_cris.py install --project isodistort --dev
.\run_cris.ps1 -m pytest ISODISTORT\tests -q --tb=line --basetemp .test-tmp-isodistort
```

若已经进入 `ISODISTORT/` 目录：

```powershell
..\run_cris.ps1 -m pytest tests -q
```

可选：从根目录执行 `.\run_cris.ps1 -m ruff check ISODISTORT`。依赖 WSL 的用例在 WSL 不可用时会跳过。

手工/长时验证说明见 [docs/MANUAL_VALIDATION.md](docs/MANUAL_VALIDATION.md)，实际脚本位于 `tests/manual/`，例如：

```powershell
..\run_cris.ps1 tests\manual\run_web.py spotcheck
..\run_cris.ps1 tests\manual\run_batch.py cif30
```

Method 1–4 的开发阶段、后续顺序和完成标准维护在
[docs/DEVELOPMENT_PLAN.md](docs/DEVELOPMENT_PLAN.md)；已修 bug、人工可读的精确统计和机器报告索引维护在
[docs/BUGFIX_VALIDATION_REPORT.md](docs/BUGFIX_VALIDATION_REPORT.md)。当前人工逐组下载步骤与可读目录名见
[docs/DOWNLOAD_CHECKLIST.md](docs/DOWNLOAD_CHECKLIST.md)。运行
`tests/manual/validate_method_outputs.py` 可检查候选/CIF、四个核心文件和模式数；
`--live-method12` 使用带源码/输入签名的原子检查点，可安全续跑。

---

## 11. 常见问题（简表）

| 现象 | 处理 |
| --- | --- |
| `wsl` 找不到 / WSL 报错 | 安装并设默认发行版 |
| `iso` 找不到 / Permission denied | 确认 `resources/isobyu/iso` 为 Linux 二进制 |
| Method 1 很慢 | 首次枚举全部特殊 k 点，属正常 |
| Method 1 表与官网看起来不同 | 字段对齐官网序参量页，本地用表格而非 radio 长行；`basis`/`origin` 来自 iso，`k-active` 由通用算法生成（不再用按 irrep 硬编码的官网 token 表），个别项可能与官网设定差一个点群等价代表 |
| Types 改了结果不变 | 忘记点 **"Change"** |
| ZIP 报没有子群 | 先对下拉框里选中的 Method 1/2/3 点 OK |
| Method 2 参数 k 点无子群 | 勾选 / 回答 **"Generate isotropy subgroups database if missing"** 后重试（可能极慢）；可用 **"Manage cached subgroup databases"** 查看或删除已生成的 `i*.iso` |
| Method 2 缓存占磁盘 | 在网页使用 **Manage cached subgroup databases** 选择并删除缓存 |
| Method 2 短目录名对上了但 SG/HM 不同 | `<IR>_<OPD>` 只用于定位；先核对 CIF 内完整候选身份和 OPD / `k-active` 行。本地 ZIP 按真实 IR/OPD 生成短名，审计不以目录名作为科学身份 |
| 手动修改了官网下载 HTML 的文件名 | 在原候选目录内可任意修改 basename，无需重下，但必须保留 `.html` 或 `.htm` 后缀。验证器按页面内容结构/角色和 SHA-256 识别；所需角色为 0 页或多页都会拒绝。跨候选目录移动或改成非 HTML 后缀不属于安全改名 |
| 模式数少于官网 / `nstrain=0` | 先核对候选完整 identity、Types 与 nmod。未勾选 strain 时 `nstrain=0` 是预期；已勾选时必须取得并核验 ISO rank-[12] canonical basis，否则导出应明确报错。若位移模式数仍少，运行对应的 current-source live 审计并按 bug 处理 |
| Method 1 候选 basis 与官网不同 | 不按短目录名判定；运行统一审计器，以 CIF 候选身份配对并用精确整幺模变换证明同一子格，再归一为官网表示比较 |
| 端口被占用 | 改 `web_port` 或关掉旧的 `scripts/main_web.py` |
| PowerShell 窗口打印 URL 但没有打开浏览器 | 无论使用直接入口还是 `.\run_cris.ps1`，若没有 `Browser launch requested`，按窗口提示手动打开最终 URL |
| PowerShell 中 `wsl` 正常，但 Python 报 `Wsl/E_ACCESSDENIED` | 这是标准 venv stub 在本机 OneDrive 路径中的进程限制，不是缺少管理员权限。重新运行 `py -3.10 setup_cris.py install --project isodistort` 让安装器做 A/B 探针和硬链接修复；若它明确报告无法修复，再用 `.\run_cris.ps1 <原 Python 参数>` |
| 每次 WSL 启动都提示 localhost proxy 未镜像 | Clash 的 `127.0.0.1` 代理与 WSL NAT 的提示不等于计算失败；本地 ISOTROPY 二进制不需要联网。只有 WSL 内确需联网时才配置 mirrored networking 或可从 WSL 访问的 Windows 地址 |
| OneDrive 路径偶发文件锁 | 可拷到非同步本地盘再试 |

---

## 12. 目录结构（本项目）

源码按四个部分组织；每部分只承担自己的职责，改某个功能只需改对应部分。

```text
ISODISTORT/
├── agent.md                开发规则、修改边界、四部分职责、验收规程
├── README.md            本文件（用户使用说明）
├── pyproject.toml          pytest / ruff / 打包配置（留在根目录）
├── .gitignore              忽略规则（留在根目录）
├── requirements*.txt       依赖清单（留在根目录）
├── backend/                ① 通用底层：所有 Method 共用的计算与设施
│   ├── wrappers/           iso / findsym / smodes 二进制封装、WSL 暂存、ISO 缓存
│   ├── models/             模式与来源领域模型
│   ├── tables/             CDML / Kovalev / 官网 k 点与空间群查表
│   ├── utils/              精确有理数晶格与群论工具、文本解析、配置、异常、自检
│   └── api/                IsoDistort 会话 API 入口
├── features/               ② 中端：按功能划分，每包只含该功能实现
│   ├── input_cif/          CIF 解析、坐标变换、对称性校验
│   ├── method1/            Method 1 搜索与共用畸变引擎
│   ├── method2/            Method 2 参数 k 模式计算与母相→子胞映射
│   ├── method3/            Inverse Landau/COPL 与 coupled 见证
│   ├── method4/            分解、均匀应变与生成回代
│   └── export/             四种下载格式 writer + CIF 位移/来源共享合同
├── frontend/               ③ 网页支撑
│   ├── web/                server.py + index.html + static/
│   └── i18n/               网页英文文案
├── scripts/                入口脚本：main_web.py（网页）、main_requirement.py（安装转发）
├── tests/                  ④ 其他：pytest 测试 + manual/ 手工与长时脚本 + fixtures/
├── resources/
│   ├── config/settings.yaml  运行时默认值
│   └── isobyu/             【只读】ISOTROPY 套件二进制与数据库
├── docs/                   说明文档、验证指南、来源记录、manifests/ 下载清单
└── output/              运行产物（不入库）
```

官网概念帮助：[isodistorthelp.php](https://iso.byu.edu/isodistorthelp.php)。

把 Method 2 完整模式写成 GD 笔记本输入：在与 CRIS 同级的 `GD/` 用 **CRIS `.venv`** 跑 `isodistort_to_gd.py`（默认 EuAl4 LD1 C1）。转换器默认发现同级 `CRIS/`，非同级布局可设置 `CRIS_ROOT`。产物在 GD 的 `generated/`，并安装 `LD1_C1_alris_functions.py`；不要覆盖 `D:\OneDrive\...` 原始精修文件，也不要改 tianren 参考笔记本。详情见 GD 的 `README.md`。
