# ISODISTORT 本地化项目（Isodistort_back）

本项目将 [ISODISTORT 网站](https://iso.byu.edu/isodistort.php)（BYU 的晶体畸变在线计算服务）的核心功能本地化，
基于 ISOTROPY Suite 的 `iso` / `findsym` Linux 二进制（位于 `isobyu/`，只读）在本地复现
“Search Page → Distortion Page”的完整工作流，避免官网服务器崩溃导致的科研进度延误。

## 一、它能做什么（与官网功能对应）

| 官网页面 | 官网功能 | 本地实现 |
| --- | --- | --- |
| Search Page | 上传母相 CIF、设置 Distortion Types | 终端/网页/API 均可 |
| Search Page | Method 1: Search over all special k points | 菜单 3（真实枚举全部特殊 k 点子群 + 过滤） |
| Search Page | Method 2: General method - specific k points | 菜单 4（选择子群并计算畸变模式） |
| Search Page | Method 3: Arbitrary k + point/space group + supercell | 菜单 5（本地近似实现，见“已知差异”） |
| Search Page | Method 4: Mode decomposition | 菜单 6（最小二乘模式分解） |
| Distortion Page | 单模式 / 多模式畸变生成 | 菜单 7-1 / 7-2 |
| Distortion Page | 导出 CIF / POSCAR | 菜单 7-3 |
| Distortion Page | Domains（畴列表） | 菜单 7-4 |

底层计算全部由 `isobyu/iso`（子群枚举、模式基矢、畴）与 `isobyu/findsym`（空间群识别）完成，
与官网使用同一套数据库（`data_*.txt`），因此**子群列表与官网一致**。

## 二、环境要求

- Python >= 3.10
- Windows + WSL（`isobyu` 为 Linux ELF 二进制；WSL 必须可用，`wsl --status` 能通过）
  - 程序首次运行会自动在 WSL 用户主目录下创建短路径暂存目录（`~/.id/`）与
    ISODATA 符号链接（`~/.id/data` → isobyu 数据目录），无需手动配置
- Linux（原生运行，无需 WSL）
- 依赖：`numpy`、`pymatgen`、`pyyaml`、`spglib`

## 三、快速开始

### 1. 安装依赖

```powershell
cd "C:\Users\devou\OneDrive\Desktop\CRIS\Isodistort_back"
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 2. 三种启动方式

本项目提供三种等效的启动/使用方式（底层共用同一套真实 iso/findsym 计算）：

**方式 A：网页交互（推荐，图形界面 + 中/EN/中+EN 切换按钮）**

```powershell
python main_web.py
```

（等价于 `python web\server.py`）启动后自动打开浏览器访问
`http://127.0.0.1:8000/`（端口可在 `config/settings.yaml` 的 `runtime.web_port`
修改；被占用时自动顺延）。页面右上角按钮可在三种界面模式间循环切换：
`中`（全中文）→ `EN`（全英文）→ `中+EN`（核心专有名词英文、其余中文）。
网页端仅依赖 Python 标准库，无需额外安装。

**方式 B：终端交互菜单**

```powershell
python main_terminal.py
```

按菜单操作（主菜单第 9 项可随时在 zh / en / mixed 三种模式间切换）：

```text
Search Page
  1. 重新加载 Parent CIF        （选择母相 CIF 文件）
  2. 设置 Distortion Types      （默认 displacement + strain）
  3. Method 1 ...               （枚举全部特殊 k 点子群，可加过滤条件）
  4. Method 2 ...               （选择子群，计算畸变模式）
  5. Method 3 ...               （点群/空间群 + 超胞搜索）
  6. Method 4 ...               （畸变结构模式分解）
  7. 进入 Distortion Page       （模式生成 / 导出 / 畴）
  8. 查看当前状态
  9. 切换语言 / Switch Language
  0. 退出
```

**方式 C：Python API**

```python
from isocore.api import IsoDistort

iso = IsoDistort(language="en")     # 默认英语；可 "zh"（中文）/ "mixed"（中英混杂）
# iso.set_language("mixed")         # 运行中可随时切换
iso.load_structure("parent.cif")
# ... 后续调用见“Python API 使用”一节
```

> 说明：项目**不提供**命令行子命令式调用（如 `isodistort run --xxx`），
> 也不存在需要删除的 CLI 入口；`python main_terminal.py` / `python main_web.py`
> 即全部启动方式。
> 仓库中 `Isodistort_validate/` 的 `compare_cif.py`/`batch_compare.py` 属于独立的
> 比对检测工具（有意保留其命令行接口，与本项目的启动方式无关）。

> 首次执行 Method 1 需要枚举全部特殊 k 点（数十次 iso 调用），
> 可能耗时 10~60 秒，与官网“数据库查询”等待一致，属正常现象。

### 3. 典型流程（对应官网）

1. 加载母相 CIF（可用 `实验数据与GD代码/EuAl4 Springer (parent).cif` 试运行）
2. 执行 Method 1，得到子群候选列表
3. 记下目标子群的 `idx`
4. 执行 Method 2，输入 `idx` 得到该路径的畸变模式
   - Method 2 还支持“直接 k 点搜索”：选择 k 点 → IR → 子群（对齐官网 Method 2 流程；
     参数 k 点需输入参数，首次查询会询问是否在线生成子群数据库）
5. Distortion Page 生成单模式畸变（输入幅度），程序自动按子群基矢扩胞并导出 CIF
6. 查看畴列表（畴数 = 子群指数，与官网一致）

## 四、Python API 使用

```python
from isocore.api import IsoDistort

iso = IsoDistort(language="en")              # 默认英语（配置 runtime.language）；
                                             # 可 "zh"（中文）/ "mixed"（中英混杂）
iso.load_structure("parent.cif")             # 1. 加载母相
iso.set_language("mixed")                    # 运行中随时切换控制台输出语言

m1 = iso.search_method_1(                    # 2. Method 1
    distortion_types=["displacement", "strain"],
    crystal_system="tetragonal",
    maximal_subgroup_only=True,
)

m2 = iso.search_method_2(                    # 3. Method 2（子群序号来自 Method 1）
    subgroup_idx=m1[0].subgroup.index,
    distortion_type="displacement",
)

iso.generate_distortion(amplitude=0.1)       # 4. 生成畸变（默认按子群基矢扩胞）
iso.export("distorted", formats=["cif", "poscar"])   # 5. 导出
iso.generate_domains()                       # 6. 畴列表
```

### 语言模式切换（三种方式统一）

`isocore/i18n` 统一管理三种界面模式（进程级全局、线程安全）：

| 模式 | 说明 |
| --- | --- |
| `zh` | 全中文 |
| `en` | 全英文（**默认**，可在 `config/settings.yaml` 的 `runtime.language` 修改） |
| `mixed` | 中英混杂：核心科学专有名词用英文（space group、Wyckoff position、order parameter、supercell、isotropy subgroup…），其余输入提示与输出提示用中文衔接 |

| 使用方式 | 切换方法 |
| --- | --- |
| 网页端 | 页面右上角按钮循环切换 `中` → `EN` → `中+EN`（即时重渲染，无需刷新） |
| 终端 | 主菜单第 9 项“切换语言”；或修改 `runtime.language` |
| Python API | `IsoDistort(language="en")` 或运行中 `iso.set_language("mixed")` |

- 界面文案目录：`isocore/i18n/messages.py`（zh/en 两套，键一一对应；mixed 由 zh + 术语表合成）
- 科学术语对照表：`isocore/i18n/terms.py`（约 100 条，来源标注见下节）
- 网页端术语经 `/api/i18n` 下发，与后端共用同一份术语表（单一数据源）

## 五、项目结构

```text
Isodistort_back/
├── main_terminal.py         # 终端交互入口（方式 B）
├── main_web.py              # 网页交互入口（方式 A：启动网页并自动打开浏览器）
├── web/                     # 网页服务实现（纯标准库）
│   ├── server.py            #   本地 HTTP 服务 + JSON API（自动打开浏览器）
│   └── index.html           #   单页界面（中/EN/中+EN 切换按钮）
├── pyproject.toml           # 包配置
├── requirements.txt         # 运行时依赖
├── requirements-dev.txt     # 开发/测试依赖
├── README.md
├── config/
│   └── settings.yaml        # 全局配置（二进制路径、容差、语言、端口、输出目录）
├── isobyu/                  # 官方二进制与数据库（只读，禁止修改）
├── isocore/                 # 核心实现
│   ├── api/                 #   对外 Python API（IsoDistort，方式 C）
│   ├── backend/             #   iso / findsym 二进制封装（WSL 桥接 + 输出解析）
│   ├── structure/           #   CIF 读写、对称分析、坐标变换
│   ├── distortion/          #   子群搜索、模式映射、畸变生成、畴
│   ├── io/                  #   导出（CIF/POSCAR/XYZ/JSON）
│   ├── i18n/                #   中英双语 + 中英混杂：界面文案 + 科学术语对照表
│   └── utils/               #   配置、异常、文本解析
├── isodistort/              # 兼容包名（isodistort.* 重导出 isocore.*）
├── examples/                # 可运行的示例脚本
├── tests/                   # 测试（含真实二进制冒烟测试，WSL 不可用时自动跳过）
└── output/                  # 运行输出（CIF/POSCAR/JSON），已被 .gitignore 忽略
```

## 六、测试与代码检查

```powershell
pip install -r requirements-dev.txt
python -m pytest -q          # 单元测试 + 真实二进制冒烟测试（WSL 不可用时自动跳过）
ruff check .                 # 代码风格检查（配置见 pyproject.toml）
```

- 单元测试不依赖二进制（使用桩对象）
- `tests/test_real_binaries.py` 为真实二进制冒烟测试：WSL 不可用时自动跳过；
  需要真实计算环境时运行（`pytest tests/test_real_binaries.py -v`）

## 七、与官网的已知差异（重要）

以下差异是本地化过程中为“可运行、可维护”做出的取舍，后续版本会逐步对齐：

1. **Distortion Types 过滤时机**：官网在 Search 阶段就按类型过滤子群；
   本地在搜索阶段保留全部子群，类型过滤在模式计算阶段（BUSH）由 iso 自动完成——
   若某子群在你结构的 Wyckoff 位置上没有对应类型的模式，Method 2 会返回空模式列表。
2. **模式振幅语义**：官网使用 As/Ap（超胞归一化振幅 + normfactor）；
   本地将“位移向量（最大分量为 1）× 用户幅度”直接叠加到原子坐标上，
   方向模式与官网一致，数值换算待与官网导出 CIF 批量比对后校准。
3. **模式列表范围**：官网 Distortion Page 列出某 IR 的全部模式（主模式 +
   可共存的次级模式，如 `[Eu0:a:dsp]A2u(a)`、`[Al1:e:dsp]A1(a)`）；
   本地（iso DISPLAY BUSH）只给出主（root）模式对应 OPD 的位移模式，
   次级模式暂不列出。
4. **晶格应变模式未实现**：本地引擎只施加原子位移，不改变晶格参数。
   对**纯位移驱动**的相变（如 I4/mmm→I4mm 极化模式），生成的畸变结构
   对称性正确（已用 spglib 与 findsym 双验证）。但对**铁弹应变**相变
   （如 I4/mmm→Immm 的 GM2+ 应变模式，需 a≠b 的晶格畸变），
   仅施加位移分量不会降低晶系对称性，产物空间群与目标子群不一致。
   此类场景请结合官网输出或后续版本（应变模式开发中）。
5. **参数 k 点（非特殊 k 点）**：
   - 子群枚举：支持（需在线生成子群数据库，对应官网 “Generate isotropy
     subgroups”，可能耗时数分钟到数小时；生成的数据库会缓存，之后秒回）。
   - 模式/畸变生成：**暂不支持**。iso 的 DISPLAY BUSH 仅支持对称 k 点；
     官网对参数 k 点使用 (3+d) 维超空间机制，本地二进制无法完成。
     遇到时会给出明确错误提示。
   - 参数约定：iso 的 k 点坐标用 `2a`/`a+b` 等形式，参数值与官网可能差整数倍
     （如官网 g=1/6 对应 iso 参数 a=1/12），请按 k 点坐标形式自行换算。
6. **Method 3**：官网支持任意超胞基矢的子群在线生成；
   本地枚举仅覆盖特殊 k 点（无参数 k 点），任意基矢场景会提示需要在线生成。
7. **Method 4**：官网支持超胞/基矢/原点选择；本地要求母相与子相原子数一致，
   且需先通过 Method 2 获得模式基矢，再做最小二乘分解。
8. **magnetic 类型**：官网支持磁畸变；本地枚举中 magnetic 相关的不可约表示
   （带 `m` 前缀）不参与默认流程，如需磁畸变请自行扩展。

## 八、与官网输出比对（Isodistort_validate）

仓库中的 `Isodistort_validate/` 是独立的 CIF 比对工具，用于把本地生成的畸变 CIF
与官网导出的参考 CIF 做语义级比对（晶格、坐标、占据率、磁矩、空间群推断），
帮助定位回归。典型用法：

```powershell
cd "C:\Users\devou\OneDrive\Desktop\CRIS\Isodistort_validate"
python compare_cif.py `
  "C:\...\Isodistort_back\output\distorted_GM3-_a0p2.cif" `
  "C:\...\官网导出的参考.cif"
```

更详细的批量比对与可信参考哈希校验见 `Isodistort_validate/README.md`。
注意：只有**同一条相变路径**的本地输出与官网输出才有可比性；
由于本地模式振幅语义（见“已知差异”第 2 条）与官网 As/Ap 存在换算差异，
坐标数值可能不逐位一致，应以“空间群与位移模式一致”作为主要判据。

## 九、常见问题

- **运行时报 `wsl` 相关错误**：确认 WSL 已安装且有默认发行版（`wsl --status`）。
- **findsym 崩溃（Fortran runtime error: End of file）**：通常是路径过长被定长缓冲区
  截断。本程序已自动改用 WSL 侧短路径暂存，正常情况下不会出现；
  若修改过 `settings.yaml` 的 `temp_dir`，请确认其不是深层路径。
- **Method 1 耗时较长**：属正常（枚举全部特殊 k 点）；候选缓存可后续版本优化。
- **导出文件名含 `+`/`-`**：如 `distorted_GM2+_a0p1.cif`，系模式标签所致，
  多数工具可正常读取。
- **网页打不开**：确认 `python web\server.py` 已启动且端口未被占用；
  浏览器需能访问 `127.0.0.1`。

## 十、科学术语规范说明（重要）

项目内全部界面文案、文档与代码注释中的科学术语按以下优先级规范：

- **英文术语**：① ISODISTORT 官网（https://iso.byu.edu/isodistorthelp.php）用语；
  ② 《International Tables for Crystallography》（ITC）英文原版。
- **中文术语**：① 《晶体学名词》（全国科学技术名词审定委员会审定）；
  ② 《国际晶体学表》中文译本。

完整对照表见 `isocore/i18n/terms.py`（每条附来源标注：`[web]`/`[ITC]`/`[CT]`/`[com]`）。
以下为常用词条摘录：

| 英文（官网/ITC） | 中文（《晶体学名词》优先） | 说明 |
| --- | --- | --- |
| space group | 空间群 | CT |
| point group | 点群 | CT |
| crystal system | 晶系 | CT |
| lattice | 点阵 | CT（“晶格”为通行同义词） |
| primitive cell | 原胞 | CT |
| conventional cell | 惯用晶胞 | CT |
| supercell | 超胞 | [web] 官网超胞输入框用语 |
| fractional coordinates | 分数坐标 | CT |
| Wyckoff position | Wyckoff 位置 | CT 附录表 |
| site symmetry | 位置对称性 | CT |
| maximal subgroup | 极大子群 | CT |
| irreducible representation (IR) | 不可约表示 (IR) | CT |
| order parameter | 序参量 | 群论物理通行 |
| isotropy subgroup | 各向同性子群 | [web]（Stokes-Hatch） |
| parent structure | 母体结构 | [web]（“母相”为相变文献通行） |
| distortion | 畸变 | [web] |
| displacement | 位移 | [web] |
| strain | 应变 | [web] |
| occupancy | 占据率 | CT |
| magnetic moment | 磁矩 | CT/[web] |
| domain | 畴 | CT |
| incommensurate / commensurate | 无公度 / 公度 | [web]/[ITC] |
| superspace group | 超空间群 | [web]/[ITC] |
| Hermann-Mauguin symbol | 赫尔曼-莫甘符号 | CT |

代码内部使用的英文标识符（如 `wyckoff_letter`、`supercell`、`basis_vectors`）
是 API 名称而非界面用语，与术语表不冲突。
