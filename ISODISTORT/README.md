# ISODISTORT 本地化项目

本项目将 [ISODISTORT 网站](https://iso.byu.edu/isodistort.php)（BYU 晶体畸变在线计算服务）的核心功能本地化：基于 ISOTROPY Suite 的 `iso` Linux 二进制（位于 `isobyu/`，只读）在本地复现「Search Page → Distortion Page」的完整工作流，子群枚举、模式基矢与畴的计算结果与官网一致，可在无网络或官网不可用时完成科研级畸变分析。（`isobyu/findsym` 已不参与本地生产流程，仅为测试/参考保留。）

同一底层计算引擎对外提供**三种等效的使用方式**：终端交互菜单、网页图形界面、Python API。

---

## 一、功能特性（与官网功能对应）

| 官网页面 | 官网功能 | 本地实现 |
| --- | --- | --- |
| Search Page | 上传母相 CIF、设置 Distortion Types | 终端 / 网页 / API 均可 |
| Search Page | Method 1: Search over all special k points | 真实枚举全部特殊 k 点子群 + 多条件过滤 |
| Search Page | Method 2: General method - specific k points | 选择子群并计算畸变模式基矢 |
| Search Page | Method 3: Arbitrary k + point/space group + supercell | 本地近似实现（见「与官网的已知差异」） |
| Search Page | Method 4: Mode decomposition | 最小二乘模式分解（支持超胞畸变结构） |
| Distortion Page | 单模式 / 多模式畸变生成 | 支持位移与占据率（occupational）模式 |
| Distortion Page | 导出 CIF / POSCAR | 自动按子群基矢扩胞并导出 |
| Distortion Page | Domains（畴列表） | 畴数 = 子群指数，与官网一致 |

底层计算全部由 `isobyu/iso`（子群枚举、模式基矢、畴）完成，与官网使用同一套数据库（`data_*.txt`），因此**子群列表与官网一致**；结构的空间群/对称性识别由 pymatgen 的 `SpacegroupAnalyzer`（`isocore/structure/symmetry_validator.py`）完成。`isobyu/findsym` 在本地生产流程中**已不参与**（保留为可选/测试组件）。

补充能力：

- **中英双语**：三种使用方式均可随时在 zh / en 间切换；界面文案与科学术语由 `isocore/i18n` 统一管理（术语对照 120+ 条，来源标注规范）。
- **物理自洽性自检**：零振幅回退、子群规则、模式正交性、振幅线性、对称性守恒五项内置检查，任一失败即提示计算异常，避免静默输出错误结构。
- **科研级验证体系**：30-CIF 全覆盖验证、COD 外部真实结构验证、金标准回归（SrTiO₃ 官网三条路径）、三接口一致性、终端交互验证，详见「测试与验证」。

---

## 二、环境要求

- Python ≥ 3.10（开发与测试在 3.14 上验证）
- **Windows + WSL**：`isobyu/` 为 Linux ELF 二进制，WSL 必须可用（`wsl --status` 通过）。
  程序首次运行会自动在 WSL 用户主目录创建短路径暂存目录与 ISODATA 符号链接，无需手动配置。
- **Linux**：原生运行，无需 WSL。
- 依赖：`numpy`、`pymatgen`、`pyyaml`、`spglib`（见 `requirements.txt`）。

---

## 三、安装

统一虚拟环境与依赖安装改为使用仓库内的脚本（只创建一次 `CRIS/.venv`）：

```powershell
cd "C:\Users\devou\OneDrive\Desktop\CRIS"
python ISODISTORT\main_requirement.py
```

开发者端（额外安装 `requirements-dev.txt`）：

```powershell
python ISODISTORT\main_requirement.py --dev
```

### 部署要求（ISOTROPY 套件）

`isobyu/` 目录存放 ISOTROPY 套件的二进制与数据库文件。**必需**的是 `iso` 二进制与 `data_*.txt` 数据库（子群枚举、模式、畴依赖）；`findsym`、`smodes`、`comsubs` 等为可选（本地生产流程仅用到 `iso`，`findsym` 仅为测试/参考保留）。体积较大且不在版本库中，部署时需将套件文件放入 `ISODISTORT/isobyu/`；`iso` 缺失时程序会在调用相关功能时报错。该目录为**只读**，程序不会修改其中的任何文件。

---

## 四、快速开始（三种启动方式）

三种方式共用同一底层引擎（真实 iso 计算；结构对称性识别由 pymatgen 完成）：

### 方式 A：网页交互（推荐）

```powershell
python main_web.py
```

启动后**自动打开默认浏览器**访问实际绑定地址（端口默认 8000，可在 `config/settings.yaml` 的 `runtime.web_port` 修改；被占用时自动顺延）。页面布局对齐官网搜索页，从上到下为：

- **Parent CIF**：上传母相 CIF；页头显示官网同款信息（空间群、点阵参数、Wyckoff 位置、`Default space-group preferences:` 行）。
- **Types of distortions to be considered**：`Strain:` 单复选框 + `Displacive:` / `Occupational:` / `Magnetic:` / `Rotational:` 四行，每行带 **all / none / 各物种** 复选框（标签在复选框前，与官网一致），选中的物种即该类型模式的作用域。**各复选框完全独立、无任何自动联动**（对齐官网：勾选/取消物种不会联动 all/none，全部取消也不会自动勾 none）；点击 **Change** 提交时才按 **all > none > 物种列表** 的顺序解释（all 勾选 = 全部物种；none 勾选或全部未勾 = 不启用该类型；否则按勾选物种）。**默认勾选 Strain + Displacive 全物种**（对齐官网：加载 CIF 后 Displacive 行的各物种复选框逐个勾选，见 `webpage_info` 第 2 页 HTML）。
- **Method 1**：晶系复选框（多选 OR，triclinic/monoclinic/orthorhombic/tetragonal/trigonal/hexagonal/cubic）→ 可达子群空间群下拉（只列出与母相结构相容的子群，来自真实枚举并按会话缓存）→ Conventional/Primitive lattice 下拉（由真实枚举的子群基矢按格点等价去重生成；因本地 iso 9.6.1 与官网站点数据库存在版本差异，选项可能与官网略有不同，界面有明确告示）→ Maximal subgroups only。
- **Method 2**：`Specify k point:` → **superposed IRs** 高级功能（对齐官网）：修改 `Change number of superposed IRs:` 的数值并点击 **Change** 后，显示对应数量的 **`k vector N:` 行**（每行：k 点下拉 + a= / b= / g= 参数），可为多个 primary order parameter 分别选择不同的 k 点；点 **OK** 枚举所选 k 点的全部不可约表示子群（合并显示，含 k/IR/OPD 列）。**子群表支持列筛选与筛选结果下载**：表上方提供 Irrep / OPD / s / i / SG / k 列的文本筛选框（大小写不敏感包含匹配）、「只显示筛选后的行」开关与「下载筛选结果 (txt/csv)」链接。点击子群行显示该路径的模式基矢；**若所选子群属于参数 k 点（如 LD/DT，非公度），iso 二进制无法计算位移模式，页面以提示（非报错）说明并保留子群表**（对应官网 (3+d) 维超空间机制，本地暂不支持，见「已知差异」第 5 条）。`# of independent incommensurate modulations` 依赖官网的 (3+d) 维超空间工作流，本地已从网页/终端交互中移除，请在官网使用该功能。若本地枚举为空（常见于带参数 k 点如 LD/DT），界面**并列给出两个选项、且「调用本地资源计算」优先**：**① 调用本地资源计算** —— 用本机算力生成该 k 点的子群数据库（对应官网 Generate isotropy subgroups，本地 iso 二进制对 "Should the data base be added?" 以空行应答触发生成，可能耗时数分钟到数小时，生成后缓存）—— 若生成失败会显示带 iso 输出末尾的诊断信息；**② 去 ISODISTORT 官网重试** —— 用相同母相 CIF 与 (k 点, 参数) 在官网执行 Method 2（该 k 点在官网可正常生成子群）。服务器已改用多线程 HTTPServer，长耗时生成期间页面心跳照常、不会触发“关页自动停服”。
- **Method 3**：空间群下拉或点群下拉（空间群优先）→ direct/reciprocal radio（reciprocal 本地暂不支持，会给出明确提示）→ 3×3 基矢输入。
- **Method 4**：上传畸变 CIF → 幅度表 + RMS 残差（默认 nearest-site / 0.25，API 层可自定义）。
- **Space-Group Preferences（只读）**：本地 iso 二进制固定采用国际标准取位（即官网默认值），自定义取位（axes / cell choice / origin / SSG 等）会导致 `Syntax error`，因此本地**不提供可交互修改面板**，改为只读表格展示项目固定采用的默认值（Monoclinic axes a(b)c、cell choice 1、Orthorhombic axes abc、Trigonal axes hexagonal、Origin choice 2、SSG standard），并注明无法修改的原因。页头同时显示官网同款 `Default space-group preferences:` 行（`iso.space_group_preferences()`）。

> 网页除官网搜索页选项外，还包含 **Distortion Page 区段**：Method 2 计算模式后，可为每个模式输入幅度（含 occupational 占据率模式）并点击 **Generate** 生成混合畸变、下载导出 CIF，点击 **Domains** 查看畴列表，并可点击 **Download all (ZIP)** 一键打包下载输出目录中的全部生成文件。点击子群行后会显示该子群的完整信息（空间群、k 点、IR、OPD、超胞基矢、原点、s、i，对应官网 order parameter direction 页）。该区段与**终端菜单第 7 项**、**Python API**（`generate_mixed_distortion` / `export` / `generate_domains`）共用同一底层引擎。

页面右上角为**语言下拉菜单**（English / 中文），选中即切换。**关闭页面自动停服并释放端口**：页面每 5 秒发送心跳；关闭页面（或点击 Stop）发送 `shutdown` 信标；心跳停止超过 `runtime.web_idle_timeout`（默认 60 秒）也自动停服。浏览器从未打开过页面则服务常驻。

### 方式 B：终端交互菜单

```powershell
python main_terminal.py
```

主菜单：

```text
Search Page
  1. 重新加载 Parent CIF        （选择母相 CIF 文件）
  2. 设置 Distortion Types      （默认 strain + displacive，对齐官网；含每类的作用域物种）
  3. Method 1 ...               （枚举全部特殊 k 点子群，可加过滤条件）
  4. Method 2 ...               （选择子群，计算畸变模式；支持直接 k 点搜索）
  5. Method 3 ...               （点群/空间群 + 超胞搜索）
  6. Method 4 ...               （畸变结构模式分解）
  7. 进入 Distortion Page       （模式生成 / 导出 / 畴）
  8. 查看当前状态
  9. 切换语言
  0. 退出
```

#### 终端版：文件与路径说明

终端涉及三类路径——**输入 CIF**、**输出产物**、**只读的 ISOTROPY 套件**。理解它们的解析规则，可避免「找不到文件 / 文件写到了别处」的问题。

**① 输入文件（母相 / 子相 CIF）**

- 启动程序（或执行菜单 1「重新加载 Parent CIF」、菜单 6「Method 4」）时，程序会在**项目目录 `ISODISTORT/` 内递归搜索** `.cif` 文件，列出至多 30 个，按序号选择即可。
- **搜索范围只到项目目录**。若你的 CIF 在项目目录之外（例如 `实验数据与GD代码/EuAl4 Parent.cif`），它不会出现在列表里 → 请选 `0. 手动输入路径`，然后输入该文件的路径。
- 手动输入的相对路径是相对于**你运行 `python main_terminal.py` 时的当前工作目录（当前命令行所在目录）**，不是相对于项目目录。为免混淆，建议直接输入**绝对路径**，例如：

  ```powershell
  C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码\EuAl4 Parent.cif
  ```

- `ISODISTORT/output/`（结果目录）与 `ISODISTORT/output/tmp/`（上传/中间暂存）中的 CIF **已被自动排除**，不会混入候选列表（它们多为程序生成的成品）。若要加载某个输出文件，用上述「手动输入路径」即可。
- 加载支持 `.cif` / `.vasp` / POSCAR / `.xyz` 等（按扩展名识别），非法输入会明确报错。

**② 输出文件（导出 CIF / POSCAR、分析结果）**

- 生成的结果（单/多模式畸变 CIF、导出文件、畴列表等）默认写入**输出目录**：`config/settings.yaml` 的 `runtime.output_dir`，即 `ISODISTORT/output/`。
- 菜单 7 → 导出，会提示输入「导出文件名前缀」与格式（如 `cif,poscar`）；导出后程序会打印每个文件的**完整路径**。
- 修改输出目录：编辑 `config/settings.yaml` 的 `runtime.output_dir`（相对 `config/` 目录，如 `../output`；也可写绝对路径）。中间/上传暂存目录 `runtime.temp_dir`（默认 `../output/tmp`）同理。

**③ ISOTROPY 套件（只读，勿移动）**

- `ISODISTORT/isobyu/` 存放 ISOTROPY 套件二进制与数据库（`iso`、`findsym`、`smodes` 及 `data_*.txt` 等），**只读**。部署时把套件文件放入该目录即可（缺失时调用相关功能会报错）。
- 二进制为 Linux ELF 格式，Windows 下**自动经 WSL 调用**；封装层会自动在 WSL 用户主目录建立短路径暂存目录与 `ISODATA` 符号链接，**无需手动配置任何路径**。Linux 下原生运行、无需 WSL。

**路径速查表**

| 项目 | 默认值 | 说明 | 在哪里改 |
| --- | --- | --- | --- |
| 输入 CIF 搜索目录 | `ISODISTORT/`（项目根） | 手动画名单；范围仅项目目录 | 属于「菜单 1/6 的路径输入」，无配置文件 |
| 输出目录 | `ISODISTORT/output/` | 导出 CIF/POSCAR、分析结果 | `settings.yaml → runtime.output_dir` |
| 临时/上传暂存 | `ISODISTORT/output/tmp/` | 中间文件、网页上传 CIF | `settings.yaml → runtime.temp_dir` |
| 二进制与数据库 | `ISODISTORT/isobyu/` | 只读，WSL 经短路径+`ISODATA` 自动链接 | 仅部署时放入套件文件 |

> `settings.yaml` 中的相对路径（如 `../output`）均相对于它所在的 `config/` 目录解析；`config_loader.py` 会自动把它换算为绝对路径。无需（也不建议）手动改动 `isobyu` 相关路径。

### 方式 C：Python API

```python
from isocore.api import IsoDistort

iso = IsoDistort(language="en")     # 默认英语；可 "zh"
iso.load_structure("parent.cif")
# ... 后续调用见「Python API 使用」
```

> 说明：项目不提供命令行子命令式调用；`python main_terminal.py` / `python main_web.py` 即全部启动方式。

> 首次执行 Method 1 需要枚举全部特殊 k 点（数十次 iso 调用），可能耗时 10~60 秒，与官网「数据库查询」等待一致，属正常现象。

### 典型流程（对应官网）

1. 加载母相 CIF（可用 `实验数据与GD代码/EuAl4 Parent.cif` 试运行）。
2. 执行 Method 1，得到子群候选列表。
3. 记下目标子群的 `idx`。
4. 执行 Method 2，输入 `idx` 得到该路径的畸变模式。
5. Distortion Page 生成单模式畸变（输入幅度），自动按子群基矢扩胞并导出 CIF。
6. 查看畴列表（畴数 = 子群指数）。

---

## 五、Python API 使用

```python
from isocore.api import IsoDistort

iso = IsoDistort(language="en")              # 默认英语（配置 runtime.language）；可 "zh"
iso.load_structure("parent.cif")             # 1. 加载母相
iso.set_language("zh")                       # 运行中随时切换控制台输出语言

# 1b. 畸变类型作用域（对齐官网 all/none/物种 复选框；"*" = 全部物种）
iso.set_distortion_scope({
    "displacive": ["*"],                     # 位移型：全部原子
    "occupational": ["Al"],                  # 占据率型：仅 Al
    "strain": [],                            # 应变：无物种概念
})

m1 = iso.search_method_1(                    # 2. Method 1
    distortion_types=["displacive", "strain"],
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

- Method 1 下拉数据（可达子群空间群 + Conventional/Primitive lattice 选项）可通过 `iso.method1_options()` 获取；`iso.space_group_preferences()` 返回官网「Default space-group preferences:」行。
- 加载结构支持 `.cif` / `.vasp` / POSCAR / `.xyz` 等（`read_structure` 按扩展名识别），非法输入会给出明确错误。

### 语言切换（三种方式统一）

| 使用方式 | 切换方法 |
| --- | --- |
| 网页端 | 页面右上角下拉菜单（English / 中文），无需刷新 |
| 终端 | 主菜单第 9 项「切换语言」 |
| Python API | `IsoDistort(language="en")` 或运行中 `iso.set_language("zh")` |

界面文案：`isocore/i18n/messages.py`（zh/en 两套，键一一对应）；科学术语对照表：`isocore/i18n/terms.py`；网页端术语经 `/api/i18n` 下发，与后端共用同一数据源。

> 说明：**终端菜单的提示文案**采用「`t()` 国际化 + 少量硬编码」混排。网页端与 Python API 的界面/控制台输出已全量走 `t()`；终端仍有部分交互提示（如 `_prompt_distortion_types`、`_prompt_supercell`、Method 3/4 的若干说明行）为固定中文/英文，不受语言下拉影响。这是「可运行优先」的取舍，**不影响计算正确性**，后续可按需并入 `messages.py` 统一管理。

---

## 六、项目结构（开发 / 生产环境文件划分）

```text
ISODISTORT/
├── main_terminal.py         # 【生产】终端交互入口
├── main_web.py              # 【生产】网页交互入口（启动网页并自动打开浏览器）
├── web/                     # 【生产】网页服务实现（标准库 + pymatgen）
│   ├── server.py            #   本地 HTTP 服务 + JSON API（端口容错/自动开浏览器/关页自动停服）
│   ├── index.html           #   单页界面（对齐官网搜索页布局；语言下拉；心跳/关闭信标）
│   └── static/              #   官网样式资源
├── isocore/                 # 【生产】核心实现
│   ├── api/                 #   对外 Python API（IsoDistort）
│   ├── backend/             #   iso / findsym / smodes 二进制封装（WSL 桥接 + 输出解析）
│   ├── structure/           #   CIF/POSCAR 读写、对称分析、坐标变换
│   ├── distortion/          #   子群搜索、模式映射、畸变生成、畴
│   ├── io/                  #   导出（CIF/POSCAR/XYZ/JSON）
│   ├── i18n/                #   中英双语：界面文案 + 科学术语对照表
│   └── utils/               #   配置、异常、文本解析
├── isodistort/              # 【生产】兼容包名（isodistort.* 重导出 isocore.*）
├── config/
│   └── settings.yaml        # 【生产】全局配置（二进制路径、容差、语言、端口、输出目录）
├── isobyu/                  # 【生产·部署时放入】ISOTROPY 套件二进制与数据库（只读）
├── pyproject.toml           # 【生产+开发】包配置 / ruff 检查配置
├── requirements.txt         # 【生产】运行时依赖
├── requirements-dev.txt     # 【开发】测试/检查依赖（pytest、ruff）
├── README.md                # 【生产+开发】本文档
├── tests_dev/               # 【开发】测试与验证（见「测试与验证」；数据目录不入版本库）
└── output/                  # 运行输出（CIF/POSCAR/JSON），已忽略，不入版本库
```

**文件归属说明**：

- **客户端/生产环境部署必需**：`main_terminal.py`、`main_web.py`、`web/`、`isocore/`、`isodistort/`、`config/settings.yaml`、`requirements.txt`、`pyproject.toml`（安装用），以及部署时放入的 `isobyu/`。
- **仅开发者端需要（文件名带 `_dev` 后缀）**：`tests_dev/`（全部测试与验证脚本、30-CIF 测试集生成器、COD 外部结构下载器）。生产环境可整体移除该目录。
- `requirements-dev.txt` 中的 `-dev` 后缀即开发者端标记。
- `output/`、`.venv/`、各类缓存目录为运行时产物，不随版本库分发。

---

## 七、测试与验证

测试与验证全部位于 `tests_dev/`（开发者端），依赖 WSL（真实 iso/findsym）的用例在 WSL 不可用时自动跳过。

### 运行方式

```powershell
python ISODISTORT\main_requirement.py --dev
python -m pytest tests_dev -q        # 单元测试 + 金标准 + 三接口 + 鲁棒性（WSL 门控）
ruff check .                         # 代码风格检查（配置见 pyproject.toml）
```

### 验证体系（四层 + 两类外部数据源）

1. **金标准回归**（`test_golden_standard.py`）：以官网示例位点设置的 SrTiO₃（Pm-3m #221）为金标准母相，验证三条文献公认相变路径：**R₄⁺ → I4/mcm #140**、**M₃⁺ → P4/mbm #127**、**Γ₄⁻ → P4mm #99**；`StructureMatcher` 语义比对（坐标容差 ≤1e-5、振幅相对误差 ≤1e-4）、零振幅回退、振幅线性、官网 LD1 参考 CIF 比对。
2. **三接口一致性**（`test_three_interface.py`）：同源输入分别走 Python API 与网页 HTTP 端点，method1_options / 候选数 / 模式数 / 生成畸变空间群完全一致（数值偏差 ≤1e-6）。
3. **30-CIF 全覆盖验证**（`tests_dev/cifs_30/` 由 `make_cifs_30.py` 生成，覆盖全部 7 个晶系 30 个母相空间群 + 2 个真实 EuAl4 CIF；驱动脚本 `run_30cif_validation.py`）：加载识别、Method 1 枚举、Method 2 模式、生成畸变 spglib 对称性 == 目标子群、畴数 == 子群指数、五项物理自检、官网 LD1 → P4mm #99 路径比对。
4. **COD 外部真实结构验证**（`cifs_external/` 由 `fetch_cod_cifs.py` 按元素+空间群从 [Crystallography Open Database](https://www.crystallography.net/cod/) 下载，来源记录在 `SOURCES.md`；驱动脚本 `run_external_validation.py` 同时做 API 全流程 + 网页 HTTP 端点抽查）：用真实晶体结构（钙钛矿、岩盐、纤锌矿、赤铁矿/方解石/刚玉、金红石/锐钛矿、萤石、尖晶石、石英等）验证算法在真实数据上的正确性与网页交互无异常。
5. **终端交互验证**（`run_terminal_validation.py`）：脚本化驱动 `main_terminal.py`，对 30-CIF 与 COD 外部结构逐一比对 Method 1 候选数、Method 2 模式数、生成畸变空间群与 API 结果一致。

### 当前实测结果

| 验证项 | 结果 |
| --- | --- |
| pytest 全套 | 90 passed |
| 30-CIF（32 个母相） | 28 精确匹配 + 2 物理退化（单原子原胞刚性平移）+ 2 共享 Wyckoff 字母限制，**0 失败 / 0 静默错误** |
| 官网参考路径（EuAl4 → LD1 → P4mm #99） | path_match = True（IR LD1、OPD C1、s=12、i=24、c 轴 6 倍超胞） |
| COD 外部结构（23 个，22 个有合格 CIF） | 19 精确匹配 + 4 明确归因（2 共享字母 + 1 物理退化 + 1 取位约定差异），网页 HTTP 抽查全部正常 |
| 金标准三条路径（SrTiO₃） | R₄⁺→#140、M₃⁺→#127、Γ₄⁻→#99 全部匹配 |
| 终端交互验证（`run_terminal_validation.py`，30-CIF + COD 共 55 个） | 候选数 / 模式数 55/55 与 API 一致；畸变对称性 47 精确匹配 + 8 已知边界，**不一致 0** |

所有非通过项均有明确归因（见「已知边界」），**不存在静默输出错误结果**的情形。

---

## 八、与官网的已知差异

以下差异是本地化过程中为「可运行、可维护」做出的取舍：

1. **Distortion Types 过滤时机**：官网在 Search 阶段按类型过滤子群；本地对**参数 k 点**（LD/DT 等）已通过对 `smodes` 的 IR 活性探测 + 物种作用域在枚举后过滤（例如 EuAl4 LD g=1/6 + displacive(Eu,Al) 仅保留 LD1/LD2/LD5，与官网一致）；**特殊 k 点**在 Method 2 子群枚举时亦会按 BUSH 探测过滤无位移模式的子群。Method 1 全特殊 k 点遍历仍保留全部子群再在模式阶段过滤（下拉空间群可能比官网多几个条目，见下）。子群符号与 Schoenflies 符号已对齐官网（简短 HM + Schoenflies）。
2. **模式振幅语义**：官网使用 As/Ap（超胞归一化振幅 + normfactor）；本地以「位移向量（按最大矢量模长归一化为 1）× 用户幅度」叠加到原子坐标，方向模式与官网一致，数值换算待校准。
3. **模式列表范围**：官网 Distortion Page 列出某 IR 的全部模式（主模式 + 可共存次级模式）；本地（iso DISPLAY BUSH）只给出主（root）模式对应 OPD 的位移模式。
4. **晶格应变模式未实现**：本地引擎只施加原子位移、不改变晶格参数。对纯位移驱动的相变（如 I4/mmm→I4mm 极化模式）结果正确；对铁弹应变相变（需 a≠b 晶格畸变）仅位移分量无法降低晶系对称性，此类场景请结合官网输出。
5. **参数 k 点（非特殊 k 点）**：子群枚举支持（需在线生成子群数据库并缓存到 WSL 暂存目录，不写入只读的 `isobyu/`）；iso 的 DISPLAY ISOTROPY 流程要求**先选 IR 再设 KVALUE**（与 DISPLAY IRREP 流程顺序相反），生成提示以**空行（回车）应答触发生成**（“Enter RETURN to continue”，实测 iso 9.6.1）。模式/畸变生成暂不支持（官网使用 (3+d) 维超空间机制，本地二进制无法完成），遇到时会给出明确错误提示。`# of independent incommensurate modulations`（nmod）功能依赖官网 (3+d) 维超空间计算链路，本地网页/终端交互已移除该选项，需到官网使用。**官网 a/b/g 参数会自动换算为 iso KVALUE**（例如 I4/mmm LD 点官网 `g=1/6` → iso 内部 `2a=1/6` 即 KVALUE `1/12`，见 `isocore/data/kpoints_official.py` 的 `official_kparams_to_iso`）；未收录空间群仍按原值传入。**k 点下拉显示**：对已收录的母相空间群（目前 I4/mmm/EuAl4，`isocore/data/kpoints_official.py`），k 点下拉显示官网同款「Miller-Love 记号 + Kovalev 编号 + 官网坐标」（如 `GM, k14 (0,0,0)`）；未收录空间群回退 iso 原始坐标（无 Kovalev 编号）。Kovalev 编号与官网坐标参数化来自官网站点数据库（CDML k 点表），本地 iso 二进制不提供，需逐空间群收录。
6. **Method 3**：本地枚举仅覆盖特殊 k 点；官网的 reciprocal（倒易空间超格）选项本地暂不支持（选择后给出明确错误提示）。supercell_basis（3×3 子格基矢）按格点等价（GL(3,Z)）过滤枚举出的特殊 k 点子群；`direct_sublattice_centering` 仅支持官网默认 `d`（P/A/B/C/I/F/R 会明确报错，不再静默忽略）。
7. **Method 4**：需先通过 Method 2 获得模式基矢再做最小二乘分解；**支持畸变结构为母相的超胞**（官网 Method 4 的常规情形）：原子数不一致时自动从当前相变路径（或两晶格矩阵）确定超胞基矢，把母相与模式位移提升到超胞坐标系（含非零 k 点 Bloch 相位调制）再分解。网页端按官网默认 nearest-site/0.25，API 层可自定义匹配方法与阈值。
8. **magnetic 类型**：本地枚举中 magnetic 相关不可约表示（带 `m` 前缀）不参与默认流程，磁畸变需自行扩展。
9. **occupational（占据率）畸变（v1 近似）**：本地按子群超胞对选定物种的 Wyckoff 位点做 +1/-1 二分类占据率调制，并用 spglib 校验调制后超胞对称群是否等于目标子群；校验失败时界面标注「近似模式」。官网按 (k, IR, OPD) 精确计算每个轨道占据率的完整算法尚未实现。
10. **Method 1 的 Conventional/Primitive lattice 选项**：官网选项来自其站点数据库（官网 iso 版本预定义）；本地选项由真实枚举的子群超胞基矢生成：Conventional 按惯用格点等价（GL(3,Z) 幺模变换）分类、Primitive 在原胞坐标下分类后转回惯用坐标显示（与官网 isoplattice 下拉的显示语义一致），均保持子群枚举顺序。因本地 iso 9.6.1 与官网站点数据库的子群基矢存在版本差异，选项数量与官网略有出入（实测 EuAl4/I4/mmm：Conventional 13 vs 官网 12、Primitive 13 vs 官网 9），这是**数据库版本差异**而非算法错误，网页 Method 1 面板会显示明确告示说明；子群枚举本身与官网一致。
11. **多物种共享 Wyckoff 字母时的模式映射**：当多个物种落在同一个 Wyckoff 位置（如 Pnma 钙钛矿中 A 位与 O1 同处 4c）时，iso 的 DISPLAY BUSH 只输出一份符号化代表点，本地映射器无法按位置区分物种。此类结构（P4/nmm、R-3c、Pnma 钙钛矿的 LaMnO₃/CaTiO₃ 等）被归类为共享字母限制，**不会静默输出错误结果**。规避：改用不同 Wyckoff 位置的母相模型，或结合官网输出。
12. **多维模式的生成方向（v1 近似）**：官网按 OPD 方向确定多维 IR 各分量的权重；本地生成时取各分量等权求和再归一化（OPD 通用方向近似），对全对称（GM1+）模式与 1D 模式结果正确，特定 OPD 方向的多维模式方向可能与官网略有差异。
13. **取位/原点约定差异（真实数据场景）**：当数据库记录的结构位点设置与 iso 的轨道符号化约定不一致时（如 COD 尖晶石 MgAl₂O₄ 的 O 32e），BUSH 的符号化代表点（`x,-x+1/4,...`）无法与结构实际位置对齐，模式映射结果不可用（畸变结构退化为 P1）。此类结构被归类为取位约定限制并明确标记，不属于静默错误。

---

## 九、科学术语规范

界面文案、文档与代码注释中的科学术语按以下优先级规范：

- **英文术语**：① ISODISTORT 官网用语；② 《International Tables for Crystallography》（ITC）英文原版。
- **中文术语**：① 《晶体学名词》（全国科学技术名词审定委员会）；② 《国际晶体学表》中文译本。

完整对照表见 `isocore/i18n/terms.py`（每条附来源标注：`[web]`/`[ITC]`/`[CT]`/`[com]`）。常用词条：space group 空间群、point group 点群、crystal system 晶系、lattice 点阵、primitive cell 原胞、conventional cell 惯用晶胞、supercell 超胞、fractional coordinates 分数坐标、Wyckoff position Wyckoff 位置、irreducible representation 不可约表示、order parameter direction 序参量方向、subgroup 子群、domain 畴等。

---

## 十、常见问题

- **运行时报 `wsl` 相关错误**：确认 WSL 已安装且有默认发行版（`wsl --status`）。
- **findsym 崩溃（Fortran runtime error: End of file）**：仅在直接运行 `findsym`（如 `tests_dev/` 的二进制测试）时出现，通常因路径过长被定长缓冲区截断；程序已自动改用 WSL 侧短路径暂存，若修改过 `settings.yaml` 的 `temp_dir` 请确认其不是深层路径。本地生产流程不调用 `findsym`，故不影响正常使用。
- **Method 1 耗时较长**：属正常（枚举全部特殊 k 点）；候选按会话缓存，重复调用秒回。
- **导出文件名含 `+`/`-`**：如 `distorted_GM2+_a0p1.cif`，系模式标签所致，多数工具可正常读取。
- **网页打不开**：确认 `python main_web.py` 已启动且端口未被占用；浏览器需能访问 `127.0.0.1`。
- **关闭网页后服务自动退出**：属预期行为（关页发 `shutdown` 信标；心跳停止 `web_idle_timeout` 秒后自动停服并释放端口）。如需服务常驻，保持页面打开即可。
- **Method 3 的带心下拉**：本地 Method 3 仅支持默认带心 `d`；选择 P/A/B/C/I/F/R 会得到明确报错（本地 iso 无法按任意带心再生成子群数据库，见「已知差异」第 6 条）。
- **isobyu 缺失**：部署时需将 ISOTROPY 套件文件放入 `ISODISTORT/isobyu/`（见「部署要求」）。
