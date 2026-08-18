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

**方式 A：网页交互（推荐，图形界面 + 语言下拉菜单）**

```powershell
python main_web.py
```

（等价于 `python web\server.py`）启动后**自动打开默认浏览器**访问实际绑定的地址
（端口默认 8000，可在 `config/settings.yaml` 的 `runtime.web_port` 修改；
被占用时自动顺延，极端情况下由系统分配空闲端口——无需手动输入端口网址）。

界面布局**对齐官网搜索页**（`webpage_info/` 中保存的官网首页/搜索页资源）：
Bootstrap 导航栏 + 若干 `.well` 面板，从上到下依次为：

- **Parent CIF**：上传母相 CIF；页头显示官网同款信息
  （空间群、点阵参数、Wyckoff 位置，以及
  `Default space-group preferences: monoclinic axes a(b)c, ...` 行）
- **Types of distortions to be considered**：`Change` 按钮生效。面板与官网一致：
  `Strain` 单复选框 + `Displacive` / `Occupational` / `Magnetic` / `Rotational`
  四行，每行带 **all / none / 各物种（如 Eu、Al）** 复选框——
  选中的物种范围即该类型模式的作用域（底层计算按物种过滤 Wyckoff 位置；
  occupational 类型由本地占据率模式生成器处理，见“已知差异”）
- **Method 1**（`OK` 按钮）：晶系复选框（可多选，OR 语义）→
  **可达子群空间群下拉**（官网行为：只列出与母相结构相容的子群空间群，
  而非全部 230 个；数据来自真实枚举并按会话缓存）→
  **Conventional lattice / Primitive lattice 下拉**（官网同款：
  惯用格与原胞格两种坐标系的超胞基矢选项，选中即按“所选子格过滤”）
  → Maximal subgroups only 复选框
- **Method 2**（`OK` 按钮）：k 点下拉 → 官网同款 **a= / b= / g= 参数输入框**
  （固定显示；参数 k 点必须填写对应参数）→ IR 下拉 →
  **OPD 下拉**（列出子群结果的序参量方向，可再过滤）→
  `# of independent incommensurate modulations`（默认为 0）→
  `List subgroups` 按钮；面板下方为官网同款
  **“Change number of superposed IRs:”**（默认为 1）输入 + `Change` 按钮
  - 子群/候选结果表格**点击一行**即按该子群计算畸变模式（Method 2）
- **Method 3**（`OK` 按钮）：230 空间群下拉 **或** 32 点群（晶类）下拉（二选一，
  与官网规则一致：空间群优先）→ 官网同款 radio
  **direct（实空间子格，带心 Default/P/A/B/C/I/F/R）/ reciprocal（倒易超格）**
  → a'=… b'=… c'=… 3×3 基矢输入（reciprocal 本地暂不支持，会给出明确提示）
- **Method 4**（`OK` 按钮）：上传畸变 CIF → 幅度表 + RMS 残差
  （匹配方法与阈值参数按官网默认 nearest-site/0.25，API 层仍可自定义）
- **Distortion Page**：模式幅度输入（occupational 模式标注“近似”提示）→
  单模式/混合生成 → 导出（CIF/POSCAR）→ 畴
- **Space-Group Preferences**：完整显示官网 settings 面板的全部选项
  （Monoclinic axes / cell choice、Orthorhombic axes、Trigonal axes、
  Origin choice、Superspace group setting、parent-like basis vectors），
  本地引擎自动采用官网默认值（已按默认预选，界面暂不可修改，仅参考）

页面右上角为**语言下拉菜单**（English / 中文），选中即切换，无需刷新。
网页端除 Python 标准库外仅额外用到 pymatgen（用于生成空间群下拉与格子换算，
项目本就依赖）。

**关闭页面自动停止服务并释放端口**：页面打开后每 5 秒发送一次心跳（`/api/ping`）；
关闭页面（或点击右上角 **Stop** 按钮）会发送 `shutdown` 信标（`/api/shutdown`），
守护线程随即关闭 HTTP 服务、退出进程并**释放端口**。
若心跳因故停止（如页面异常关闭未发出信标），也会在
`runtime.web_idle_timeout`（默认 60 秒）后自动停服。
若浏览器从未打开过页面，则服务常驻，不会误停。

**方式 B：终端交互菜单**

```powershell
python main_terminal.py
```

按菜单操作（主菜单第 9 项可随时在 zh / en 两种语言间切换）：

```text
Search Page
  1. 重新加载 Parent CIF        （选择母相 CIF 文件）
  2. 设置 Distortion Types      （默认 displacive + strain；含每类的作用域物种）
  3. Method 1 ...               （枚举全部特殊 k 点子群，可加过滤条件）
  4. Method 2 ...               （选择子群，计算畸变模式）
  5. Method 3 ...               （点群/空间群 + 超胞搜索）
  6. Method 4 ...               （畸变结构模式分解）
  7. 进入 Distortion Page       （模式生成 / 导出 / 畴）
  8. 查看当前状态
  9. 切换语言
  0. 退出
```

**方式 C：Python API**

```python
from isocore.api import IsoDistort

iso = IsoDistort(language="en")     # 默认英语；可 "zh"（中文）
# iso.set_language("zh")            # 运行中可随时切换
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

1. 加载母相 CIF（可用 `实验数据与GD代码/EuAl4 Parent.cif` 试运行）
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
                                             # 可 "zh"（中文）
iso.load_structure("parent.cif")             # 1. 加载母相
iso.set_language("zh")                       # 运行中随时切换控制台输出语言

# 1b. 畸变类型作用域（对齐官网 all/none/Eu/Al 复选框；"*"=全部物种）
iso.set_distortion_scope({
    "displacive": ["*"],                     # 位移型：全部原子
    "occupational": ["Al"],                  # 占据率型：仅 Al
    "strain": [],                            # 应变：无物种概念
})

m1 = iso.search_method_1(                    # 2. Method 1
    distortion_types=["displacive", "strain"],   # 旧名 displacement/order 亦可用
    crystal_system="tetragonal",
    maximal_subgroup_only=True,
    lattice=[[1, 1, 0], [-1, 1, 0], [0, 0, 1]],  # 官网 Conventional lattice 过滤
)

m2 = iso.search_method_2(                    # 3. Method 2（子群序号来自 Method 1）
    subgroup_idx=m1[0].subgroup.index,
    distortion_type=["displacive", "occupational"],
)
# m2 之后：iso.mode_displacements（位移模式）与 iso.mode_occupancies
# （occupational 占据率模式，键形如 "occ-Al-d"）

iso.generate_distortion(amplitude=0.1)       # 4. 生成畸变（默认按子群基矢扩胞）
iso.export("distorted", formats=["cif", "poscar"])   # 5. 导出
iso.generate_domains()                       # 6. 畴列表
```

> Method 1 下拉数据（可达子群空间群 + Conventional/Primitive lattice 选项，
> 与网页版一致）可通过 `iso.method1_options()` 获取；
> `iso.space_group_preferences()` 返回官网 “Default space-group preferences:” 行。

### 语言切换（三种方式统一）

`isocore/i18n` 统一管理两种语言（进程级全局、线程安全）：

| 语言 | 说明 |
| --- | --- |
| `zh` | 全中文（界面文案为纯中文，仅保留官方专名/缩写如 Method 1、CIF、OPD） |
| `en` | 全英文（**默认**，可在 `config/settings.yaml` 的 `runtime.language` 修改） |

| 使用方式 | 切换方法 |
| --- | --- |
| 网页端 | 页面右上角**下拉菜单**（English / 中文），选中即显示所选语言（无需刷新） |
| 终端 | 主菜单第 9 项“切换语言”；或修改 `runtime.language` |
| Python API | `IsoDistort(language="en")` 或运行中 `iso.set_language("zh")` |

- 界面文案目录：`isocore/i18n/messages.py`（zh/en 两套，键一一对应，均为纯中文/纯英文）
- 科学术语对照表：`isocore/i18n/terms.py`（约 100 条，来源标注见下节）
- 网页端术语经 `/api/i18n` 下发，与后端共用同一份术语表（单一数据源）

## 五、项目结构

```text
Isodistort_back/
├── main_terminal.py         # 终端交互入口（方式 B）
├── main_web.py              # 网页交互入口（方式 A：启动网页并自动打开浏览器）
├── web/                     # 网页服务实现（标准库 + pymatgen 生成空间群下拉）
│   ├── server.py            #   本地 HTTP 服务 + JSON API（端口容错 + 自动开浏览器 + 关页自动停服）
│   ├── index.html           #   单页界面（对齐官网搜索页布局；语言下拉菜单；心跳/关闭信标）
│   └── static/              #   官网样式资源（bootstrap.css / docs.css / help.jpg）
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
│   ├── i18n/                #   中英双语：界面文案 + 科学术语对照表
│   └── utils/               #   配置、异常、文本解析
├── isodistort/              # 兼容包名（isodistort.* 重导出 isocore.*）
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

### 30-CIF 科学验证（网页/API 与终端）

`tests/` 下提供一组覆盖全部 7 个晶系、30 个母相空间群的测试 CIF
（`tests/cifs_30/`，由 `tests/make_cifs_30.py` 生成并自检）以及两个真实
EuAl4 CIF，用于对网页/API 与终端做科研级正确性验证：

```powershell
python tests/run_30cif_validation.py     # API/网页同底层：32 个母相全流程
python tests/run_terminal_validation.py  # 终端（脚本化 stdin）与 API 报告逐项比对
```

验证判据（= 官网保证的不变量）：
- 加载后识别的空间群与 CIF 一致；
- Method 1 可达子群下拉按序号升序、lattice 选项按格点等价去重；
- 对可产生位移模式的子群生成畸变结构，**spglib 校验畸变结构空间群 == 目标子群**；
- 畴数 == 子群指数；
- 官网参考比对：EuAl4 → LD1（参数 k 点）→ P4mm #99 的子群枚举
  （index=24、size=12）与官网导出的 `LD1_C1_subgroup.cif` 一致
  （参数 k 点的模式计算为已知限制，见差异 5）。

已知边界：1 原子原胞的退化结构（如单原子 bcc/fcc）在 Γ 点的极性模式是
刚性平移，纯位移畸变无法降低对称性（需应变耦合），验证中单独归类为
“degenerate-rigid-translation”，不属于算法错误。

## 七、与官网的已知差异（重要）

以下差异是本地化过程中为“可运行、可维护”做出的取舍，后续版本会逐步对齐：

1. **Distortion Types 过滤时机**：官网在 Search 阶段就按类型过滤子群；
   本地在搜索阶段保留全部子群，类型与物种作用域过滤在模式计算阶段生效——
   若某子群在你结构的 Wyckoff 位置上没有对应类型/物种的模式，Method 2
   会返回空模式列表。
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
   - **k 点下拉不显示 Kovalev 编号**（如官网 “GM, k14 (0,0,0)”）：
     iso 本地输出只有 Miller-Love 记号（GM/DT/LD...），官网的 Kovalev 编号
     来自其站点数据库，本地无法复刻，因此显示为 “GM (0,0,0)” 形式。
6. **Method 3**：官网支持任意超胞基矢的子群在线生成；
   本地枚举仅覆盖特殊 k 点（无参数 k 点），任意基矢场景会提示需要在线生成。
   官网的 **reciprocal（倒易空间超格）** 选项本地暂不支持（radio 仍保留，
   选择后会给出明确错误提示）。
7. **Method 4**：官网支持超胞/基矢/原点选择；本地要求母相与子相原子数一致，
   且需先通过 Method 2 获得模式基矢，再做最小二乘分解。
   网页端按官网默认使用 nearest-site/0.25；如需 robust/自定义阈值，请用
   Python API（`search_method_4(atom_matching_method=..., robust_distance_threshold=...)`）。
8. **magnetic 类型**：官网支持磁畸变；本地枚举中 magnetic 相关的不可约表示
   （带 `m` 前缀）不参与默认流程，如需磁畸变请自行扩展。
9. **occupational（占据率）畸变（v1 近似算法）**：本地按子群超胞对选定物种的
   Wyckoff 位点做 +1/-1 二分类占据率调制（+1 类保持全占据，-1 类占据率
   1-amplitude），并用 spglib 校验调制后超胞的对称群是否等于目标子群：
   - 校验通过（validated=True）：模式与子群对称性一致；
   - 校验失败（validated=False）：界面会标注“近似模式”，请谨慎使用。
   官网按 (k, IR, OPD) 精确计算每个轨道占据率的完整算法（含多轨道字符模式
   与富占据型有序）尚未实现；t 子群（无超胞）与奇数分裂不产生占据率模式。
10. **Method 1 的 Conventional/Primitive lattice 选项**：官网选项来自其站点
    数据库；本地选项由真实枚举（会话缓存）得到的子群超胞基矢生成，并按
    **格点等价**（GL(3,Z) 幺模变换）去重、按行列式升序排列——
    Conventional 为惯用坐标表达、Primitive 为原胞坐标表达（同一组超胞的
    两种坐标），标签格式与官网一致（如 `(1,0,0),(0,1,0),(0,0,1)`）。
    去重后选项数与官网同量级（实测 I4/mmm：本地 13 个 vs 官网 12 个）；
    具体选项集合因本地 iso 子群数据库与官网站点数据库的差异，可能与官网
    略有出入。Method 1 的**空间群下拉**只列出真实枚举得到的可达子群
    （与官网“排除不相容对称性”的行为一致），且按序号升序、符号起始列对齐。

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
- **关闭网页后服务自动退出**：属预期行为（关页发 `shutdown` 信标；心跳停止
  `web_idle_timeout` 秒后自动停服并释放端口）。如需服务常驻，保持页面打开即可。
- **Method 3 的带心下拉不生效**：本地 Method 3 为近似实现，目前仅
  空间群/点群过滤与对角超胞匹配生效，`direct_sublattice_centering` 仅为
  对齐官网表单保留（见“已知差异”第 6 条）。

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
