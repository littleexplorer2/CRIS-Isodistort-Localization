# ISODISTORT（本地版）

本目录把 [ISODISTORT](https://iso.byu.edu/isodistort.php)（BYU 的晶体畸变搜索工具）做成**可离线运行**的程序。

你不需要先懂晶体学术语。用白话理解工作流即可：

1. 准备一份**母相**晶体结构文件（CIF）：描述「还没畸变」时的对称性与原子位置。  
2. 告诉程序要考虑哪些**畸变类型**（应变、原子位移、占位、磁、旋转等）。  
3. 用 Method 1 / 2 / 3 **搜索**可能的子群（对称性降低后的候选结构列表）；或用 Method 4 把已经畸变的结构**分解**成模式幅度。  
4. 在 **Distortion** 区下载筛选后的结果表，以及子群结构文件（CIF、IsoVIZ、模式详情、TOPAS 等）。

计算后端是 ISOTROPY Suite 的 Linux 程序 `iso` / `smodes` 与 `data_*.txt` 数据库（放在本目录的 `isobyu/`，只读）。同一套引擎支持三种用法：**网页**、**终端菜单**、**Python API**。网页与终端界面为**英语**（没有语言切换）。界面上看到的英文标签，下文一律用引号标出。

本项目与仓库根目录的 `CRIS/.venv` 共用一份 Python 虚拟环境。子项目之间的关系见仓库根目录 [README.md](../README.md)。

---

## 1. 这个工具具体做什么

| 步骤 | 你做什么 | 得到什么 |
| --- | --- | --- |
| 加载母相 | 在 "Parent CIF" 选文件并点 "Load" | 页头显示空间群、晶格、Wyckoff 位点 |
| 设畸变类型 | 勾选后点 "Change" | 后续 Method 按你勾选的类型过滤 |
| Method 1–3 | 设过滤条件，点 "OK" | 一张可筛选、排序的子群结果表；点一行可看模式 |
| Method 4 | 上传女儿相 CIF，点 "OK" | 模式幅度表 + RMS residual（不是子群列表） |
| Distortion | 选 Method、下表 / 勾格式、下 ZIP | 筛选后的 txt/csv；Method 1–3 的结构文件 ZIP |

**官网 Distortion 页上的 "Generate"（按模式幅度生成畸变结构）和 "Domains"（畴列表）已从本地网页/终端删除**——本项目目标是「得到子群结构信息或文件」，Method OK + ZIP 已给出零振幅超胞 CIF 与模式文件。Python API 仍保留 `generate_distortion` / `generate_domains` 供脚本使用。

Method 1 点 OK 后的结果表字段与官网序参量方向页（`webpage_info/` 中 `a.` 开头的存档）一致：官网是一条条 radio，本地仍用可筛选、排序的表格展示同样的 token（`Irrep` / `OPD` / `Dir` / `SG` / `basis` / `origin` / `s` / `i` / `k-active`），并多一列 `idx` 供点选计算模式。

请勿修改 `ISODISTORT/isobyu/` 内的文件。

---

## 2. 安装与环境

### 2.1 需要什么

| 项目 | 要求 |
| --- | --- |
| Python | **3.10 或更高**（开发在 3.12.5 上验证过） |
| Windows | **必须安装 WSL**（`isobyu/iso` 是 Linux ELF，通过 WSL 调用） |
| Linux 本机 | 可直接跑 `iso`，无需 WSL |
| ISOTROPY 套件 | 自行下载 Linux 版，放入 `ISODISTORT/isobyu/`（仓库**不附带**二进制） |

检查 Python：

```powershell
python --version
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

### 2.3 放入 ISOTROPY Suite（必需）

从 [ISOTROPY Suite](https://iso.byu.edu/isotropy.php) 下载 **Linux** 版，把至少这些文件放进：

```text
ISODISTORT/isobyu/
```

- 可执行文件：`iso`、`smodes`（可选还有 `findsym`、`comsubs`）  
- 数据库：全部 `data_*.txt`

**不要**把 Windows 可执行文件放进去。`main_requirement.py` **不会**自动下载套件；若缺少 `isobyu/`，只会新建空目录并打印下载地址。

### 2.4 安装 Python 依赖

在**仓库根目录**（不是只在 ISODISTORT 里）执行：

```powershell
cd <CRIS 根目录>
python ISODISTORT\main_requirement.py
```

该脚本会：

1. 确认 Python ≥ 3.10  
2. 在 Windows 上检查 WSL  
3. 确保存在 `ISODISTORT/output/`（及 `output/tmp/`）  
4. 检查 `isobyu/`（没有则建空目录并提醒下载）  
5. 若缺少 `ISODISTORT_VALIDATE/compare/{item,true}` 则自动创建  
6. 若缺少 `ISOVIZ_INPUT/input_content/{data.csv,subgroup.isoviz}` 则自动创建（ISOVIZ 子项目不使用 `output/`）  
7. 创建或复用根目录 `CRIS/.venv`，只安装尚未存在的依赖  
8. 若 `isobyu` 里已有数据库，会检查并配置 **ISODATA**（`iso` 读数据库用的环境变量；运行时由配置自动设置，WSL 侧会建短路径符号链接，一般**不必**在 Windows 系统属性里永久 `setx ISODATA`）

可选参数：

| 参数 | 作用 |
| --- | --- |
| `--dev` | 额外安装 `requirements-dev.txt`（pytest / ruff） |
| `--recreate` | 删除后重建 `.venv` |

成功后请用提示的解释器，例如：

```powershell
.\.venv\Scripts\python.exe ISODISTORT\main_web.py
```

激活虚拟环境（可选）：

```powershell
.\.venv\Scripts\Activate.ps1
```

**若运行很久像卡死，且终端提示符没有 `(.venv)`：** 多半是虚拟环境没有正确加载。请关闭当前进程，用上面的 `.\.venv\Scripts\python.exe …` 重新启动（或先 `Activate.ps1` 再运行）。不要用未激活的系统 Python 去跑长计算。

若 PowerShell 禁止执行脚本，直接用上面的全路径即可。

### 2.5 安装后至少应看到

- `ISODISTORT/main_web.py`、`main_terminal.py`、`main_requirement.py`  
- `ISODISTORT/config/settings.yaml`  
- `ISODISTORT/web/index.html`  
- `ISODISTORT/isobyu/iso` 以及若干 `data_*.txt`

---

## 3. 如何启动

### 3.1 网页

```powershell
cd <CRIS 根目录>
.\.venv\Scripts\python.exe ISODISTORT\main_web.py
```

默认打开 `http://127.0.0.1:8000/`。端口见 `config/settings.yaml` 的 `runtime.web_port`；被占用时会自动顺延约 20 个端口，再不行则让系统分配。控制台会打印最终 URL。

- 右上角 **"Stop"**：停止服务并释放端口。  
- 关闭浏览器标签后约 **60 秒**（`runtime.web_idle_timeout`）自动停服；要常驻就不要关页。  
- 顶栏 **"ISODISTORT"**：回到本页；**"SUITE"** / **"HELP"**：打开官网（需联网，计算本身不依赖）。

### 3.2 终端

```powershell
.\.venv\Scripts\python.exe ISODISTORT\main_terminal.py
```

先选母相 CIF（在 `ISODISTORT/` 下最多列出 30 个 `.cif`，排除 `output/`；仓库外文件请选手动输入路径），再进入 Search Page 菜单。方括号里的值是默认值，直接回车即采用。

非交互 (3+d) 超空间内核（不进菜单；`d` 即 nmod）：

```powershell
.\.venv\Scripts\python.exe ISODISTORT\main_terminal.py --superspace-d 1 --space-group 139 --q-vectors 0,0,1/6 --k-label LD --export ISODISTORT\output\superspace_nmod1.json
```

### 3.3 配置文件：`config/settings.yaml`

相对路径均相对 `config/` 目录解析。

| 键 | 含义 | 默认 |
| --- | --- | --- |
| `isobyu.bin_dir` / `data_dir` | `iso` 与数据库目录 | `../isobyu` |
| `isobyu.iso_bin` 等 | 可执行文件名 | `iso` / `findsym` / `smodes` |
| `defaults.position_tolerance` | 分数坐标容差 | `0.001` |
| `defaults.lattice_tolerance` | 晶格容差 | `0.00001` |
| `defaults.default_amplitude` | API 生成畸变时的默认幅度 | `1.0` |
| `defaults.eps` | 全局浮点容差 EPS（超空间对称判定 / 波矢约化） | `0.00001`（与 `lattice_tolerance` 相同） |
| `defaults.max_nmod` | (3+d) 附加维度上限（官网 nmod） | `3` |
| `runtime.web_port` | 网页首选端口 | `8000` |
| `runtime.web_idle_timeout` | 关页后自动停服秒数 | `60` |
| `runtime.temp_dir` | 上传暂存 | `../output/tmp` |
| `runtime.output_dir` | 终端导出等 | `../output` |
| `runtime.timeout` | 普通 `iso` 调用超时（秒） | `60` |
| `runtime.generation_timeout` | Method 2 勾选生成缺失子群库时的超时 | `3600` |

说明：网页 ZIP **不写入、不扫描** `output/`；网页上的筛选表 txt/csv 在浏览器里生成。终端菜单 7 才会把表和结构文件写到 `output/`。

---

## 4. Parent CIF（母相）加载后显示什么

界面区块标题：**"Parent CIF"**。选择 `.cif` / `.CIF`，点 **"Load"**。必须先 Load 成功，再跑任何 Method。

加载成功后，状态区会显示类似官网的页头，例如（EuAl4）：

```text
Space Group: 139 I4/mmm D4h-17
Lattice parameters: a= …, b= …, c= …, alpha= …, beta= …, gamma= …
Default space-group preferences: …
Eu1 2a (0,0,0),
Al1 4d (0,1/2,1/4),
Al2 4e (0,0,z), z= 0.38000
```

含义（白话）：

| 行 | 意思 |
| --- | --- |
| **Space Group** | 空间群编号 + 国际符号 + Schoenflies 符号（如 `D4h-17`） |
| **Lattice parameters** | 晶胞边长与夹角，**五位小数** |
| **Default space-group preferences** | 本地固定的国际标准取位约定（见文末 Preferences 表，只读） |
| **Wyckoff 行** | 每种独立原子位点：物种序号、多重度+字母、坐标；自由坐标写成 `z= 0.38000` 这种形式 |

实验示例（若仓库中仍有、且你未改只读目录）：`experiment_data/EuAl4 Parent.cif`。

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
| **"Conventional lattice:"** / **"Primitive lattice:"** | 互斥过滤超胞基矢（选一个会清空另一个）。标签形如 `(1,0,0),(0,1,0),(0,0,1)`。本地 iso 数据库与官网不完全相同，下拉可能略有差别；对常见母相（如 I4/mmm #139）会尽量对齐官网标签与顺序 |
| **"Maximal subgroups only"** | 只保留极大子群 |
| **"OK"** | 开始枚举；首次可能数秒到数十秒 |

结果表列与官网 Method 1 序参量方向页（radio 一行里的字段）对齐，但以**表格**展示，而不是一条一条的 radio：`idx` / `Irrep` / `OPD` / `Dir` / `SG` / `basis` / `origin` / `s` / `i` / `k-active`。`idx` 是本地点选/导出用的子群序号，官网可见文本里没有这一列。下载的 txt/csv 使用相同列。勾选 lattice strain 时会保留 smodes 没有的纯应变 irrep（例如 I4/mmm 的 GM4+ Fmmm）。`k-active` 使用官网列出的特殊 k 坐标（I4/mmm 的 M 点为 `(1,1,1)`，不会收成 `(0,0,1)`）。对 I4/mmm 等已收录母相，`basis` / `origin` / `k-active` 与官网序参量页同一套 token。

### 6.2 Method 2: General method - search over specific k points

针对你指定的 **k 点**（可多组叠加 IR）枚举子群。

| 选项 | 作用 |
| --- | --- |
| **"Specify k point:"** | 从母相允许的 k 点列表选择 |
| 参数 k 点数值框 | 如 LD 的 a/b/g，必须填齐再 OK |
| **"Generate isotropy subgroups database if missing"** | 若当前 k 点在本地 `data_*.txt` **没有**现成子群列表，则由本机 `iso` **现场生成**（可能很慢，超时见 `generation_timeout`）。这与已删除的 Distortion **"Generate"**（按幅度生成畸变结构）**不是同一功能** |
| **"Change number of superposed IRs:"** + **"Change"** | 改叠加 IR 组数；**必须点 Change** 才会出现多组 k vector |
| **"OK"** | 枚举该（组）k 点全部 IR 的子群 |

结果表列大致含：`idx` / `SG` / `k` / `Irrep` / `OPD` / `s` / `i`。

注意：

- 本地 **nmod**（`# of independent incommensurate modulations`）即超空间附加维度 **d**。`nmod=0` 为公度三维（走 `iso`）；`nmod=1..3` 走 `isocore` 的 **(3+d) 超空间内核**（取位固定为 **standard (IT-C)**，与页底 Space-Group Preferences 一致）。参数 k 点（LD/DT）在 `nmod=0` 时仍不能用 `iso` 算位移模式；设 `nmod≥1` 后可用本地内核生成超空间模式并投影回三维。  
- **参数 k 点**（如 LD/DT）可以列出子群；位移模式在 **nmod=0** 时本地 `iso` 无法计算，选中该 k 点会提示。设 **nmod≥1** 后走 (3+d) 超空间内核（IT-C），模式表与 ZIP 中的模式类格式可被填充。表与 CIF 结构在 nmod=0 时仍可下载。  
- 子群数据库缺失时：网页提供「本地生成 / 去官网」；终端同样有对应选项。

### 6.3 Method 3: Search over arbitrary k points for specified space group and lattice

按指定空间群/点群与格子搜索。

| 选项 | 作用 |
| --- | --- |
| **"Select either space group symmetry:"** | 选 230 个空间群之一 |
| **"or point group (crystal class):"** | 或选 32 个点群之一 |
| **"Specify a real-space sublattice… Default / P / A / B / C / I / F / R centering"** | **direct** 实空间子格。本地 Method 3 **只真正支持 Default (`d`)**；选其它带心时后端会按默认处理或提示 |
| **"Specify a primitive reciprocal-space superlattice"** | **reciprocal**：**本地不支持**，请用 direct |
| **"Choose a representative basis:"** | 3×3 基矢：`a'` / `b'` / `c'` 相对母相 a,b,c 的系数（可填分数） |
| **"OK"** | 得到候选子群表 |

结果表列大致含：`idx` / `SG` / `k point` / `Irrep` / `point group`。

### 6.4 Method 4: Mode decomposition of a distorted structure

| 选项 | 作用 |
| --- | --- |
| **"Upload distorted structure from CIF file:"** | 上传**已经畸变**的女儿相 CIF |
| **"OK"** | 输出**全部**模式幅度（可 Filter / 排序）和 RMS residual |

这是**分解**，不是子群列表。**不能**作为 Distortion ZIP 的结构来源，但可以在 Distortion 下载该幅度表的 txt/csv。

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
| **"Download all (ZIP)"** | 打包 Method **1 / 2 / 3** 当前筛选命中的子群文件（无筛选则全部），**只含勾选格式**，**不扫描** `output/`。选 Method 4 再点 ZIP 会提示改用表格下载。勾选了 isoviz / modes / topas 时，会对**非参数 k 点**子群补跑 Method 2 以填充模式（可能较慢，界面有进度条）；参数 k 点（如 LD）在 **nmod=0** 时本地 iso 无法算位移模式（对应文件中模式为空），**nmod≥1** 时由 (3+d) 超空间内核填充。仅 CIF 或 URL 带 `compute_modes=0` 时可跳过补算 |

压缩包名形如 `isodistort_methodN.zip`，解压后**直接是各子群文件夹**（与官网下载结构一致）：

```text
LD1 C1/                          # Method 2/3：IR + OPD
  subgroup.cif
  data.isoviz
  Complete modes details.txt     # 官网为 HTML 页；本地用 .txt
  topas.str
GM1+ P1 (a) 139 I4/mmm, .../     # Method 1：完整 OPD 行作文件夹名
  subgroup.cif
  ...
```

单行点选后出现的模式表在 Method 2 区域，仅供查看。

### Space-Group Preferences（只读）

页底 **"Space-Group Preferences"** 固定为国际标准取位（Monoclinic axes a(b)c、cell choice 1、Orthorhombic abc、Trigonal hexagonal、Origin choice 2、Superspace standard）。本地 `iso` 不能改这些选项；Method 1 表里的 `basis`/`origin` 对已收录母相（如 I4/mmm #139）按官网序参量页显示。

---

## 8. 网页 vs 终端

网页与终端调用**同一套** `isocore` API（Method 1–4 参数、子群枚举、`export_subgroups_zip` 的 `wrapping=None` / Method 1 OPD 文件夹名 / 模式补算策略一致）。差异主要在交互壳：

| | 网页 | 终端 |
| --- | --- | --- |
| 启动 | `main_web.py` | `main_terminal.py` |
| 交互 | 一页上全部面板 | 分步菜单（1–8 / 0） |
| 超空间 / nmod / 参数 k | Method 2 的 nmod 输入（d=0..3）+ (3+d) 面板；参数 k 在 nmod=0 时提示 iso 无位移模式 | Method 2 询问 nmod；菜单 9 跑超空间；`--superspace-d` |
| 子群库缺失 | 按钮：本地生成 / 官网链接 | 数字选项：本地生成 / 打印官网 URL |
| 下载 | 浏览器 ZIP / txt/csv | 菜单 **7. Distortion** → ZIP 或目录写到 `output/` |
| 模式补算 | 默认开启；URL `compute_modes=0` 可关 | 导出时询问（默认 yes，对应网页默认） |

终端主菜单摘要：

1. Reload parent CIF  
2. Set distortion types  
3–6. Method 1–4（表交互：`f` 筛选、`s` 排序、`c` 清除、`only`、输入 `idx` 算模式、`q` 结束）  
7. Distortion（导出结构文件或筛选表；参数 k 点在 nmod=0 时会提示模式为空）  
8. Show current state（含固定 Space-Group Preferences / nmod）  
9. (3+d) superspace（nmod / k_s / q-vectors / 导出 JSON）  
0. Exit  

结构文件格式提示默认：`cif,isoviz,modes,topas`。**没有** Generate / Domains。

---

## 9. Python API 快速示例

在已安装依赖、且 `PYTHONPATH` 含 `ISODISTORT/`（或从该目录运行）时：

```python
from isocore.api import IsoDistort

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
iso.export_subgroups("out_batch", formats=["cif", "isoviz", "modes", "topas"])

# (3+d) 超空间内核（nmod = d；不依赖 iso DISPLAY BUSH）
ss = iso.run_superspace(nmod=1, q_vectors=[[0, 0, 1 / 6]], k_point_label="LD")
```

`export_subgroups` 使用当前 `iso.subgroups`（最近一次写入的列表）。网页 ZIP 会显式传入所选 Method 的列表，避免三种 Method 混在一起。

`generate_distortion` / `generate_mixed_distortion` / `generate_domains` 仍在 API 中，网页和终端不再调用。

---

## 10. 已知限制（影响日常使用）

1. **Windows 必须经 WSL** 调用 Linux 版 `iso`。  
2. **应变模式未实现**：不改晶格参数。  
3. **参数 k 点**：可枚举子群；`nmod=0` 时位移模式仍不能由本地 `iso` 计算。设 `nmod≥1` 后由 `isocore.superspace` 生成 (3+d) 模式并投影到三维。  
4. **nmod**（独立非公度调制数 = d）：网页 Method 2 与 (3+d) 面板、终端菜单 9、CLI `--superspace-d`、API `run_superspace` / `search_method_2(..., number_of_independent_modulations=nmod)` 均已接通。上限 `defaults.max_nmod`（默认 3）。取位固定 **standard (IT-C)**。  
5. **Method 3**：reciprocal 不支持；带心仅 Default (`d`)。  
6. **magnetic**：带 `m` 前缀的 IR 默认不进入流程。  
7. **occupational**：本地为 ±1 占据近似，校验失败会标明。  
8. **Distortion Generate / Domains**：官网有，本地网页/终端已去掉。  
9. **界面仅英语**。  
10. **ISOVIZ / modes / TOPAS** 按官网选项语义在本地生成：文件名与文件夹布局对齐官网（`subgroup.cif` / `data.isoviz` / `topas.str` / `Complete modes details.txt`）；`.isoviz` 使用官网 `!tag` 布局，可供 IsoVIZ / ISOVIZ_INPUT 使用，但原子轨道展开与模式向量未必与网站逐字节一致。CIF 按官网 Distortion 页格式写出（子群空间群、常规原点、不对称单元、Hall 设置、`k-active` 参数代入）。**参数 k 点**（LD 等）在 `nmod=0` 时位移模式循环可能为 0；`nmod≥1` 时由本地 (3+d) 内核填充。特殊 k 点子群在勾选模式类格式时会补跑 Method 2 填充模式。

更细的差异列表见历史开发记录；上表是用户最常撞到的几条。

---

## 11. 测试如何运行

从仓库根目录：

```powershell
cd <CRIS 根目录>
python ISODISTORT\main_requirement.py --dev
.\.venv\Scripts\python.exe -m pytest ISODISTORT\tests_dev -q
```

(3+d) 超空间内核（不依赖 WSL / `iso`）：

```powershell
python -m pytest ISODISTORT/tests_dev/test_3pd.py -q
```

若已经进入 `ISODISTORT/` 目录：

```powershell
..\.venv\Scripts\python.exe -m pytest tests_dev -q
```

可选：在 `ISODISTORT/` 下执行 `ruff check .`。依赖 WSL 的用例在 WSL 不可用时会跳过。

手工/长时验证（生成 CIF、网页抽查、批量回归）见 `tests_dev/manual/`，例如：

```powershell
..\.venv\Scripts\python.exe tests_dev\manual\run_web.py spotcheck
..\.venv\Scripts\python.exe tests_dev\manual\run_batch.py cif30
```

说明见 `tests_dev/manual/README.md`。

---

## 12. 常见问题（简表）

| 现象 | 处理 |
| --- | --- |
| `wsl` 找不到 / WSL 报错 | 安装并设默认发行版 |
| `iso` 找不到 / Permission denied | 确认 `isobyu/iso` 为 Linux 二进制 |
| Method 1 很慢 | 首次枚举全部特殊 k 点，属正常 |
| Method 1 表与官网看起来不同 | 字段与官网序参量页同一套 token，只是本地用表格而不是 radio 长行；I4/mmm 的 GM4+ / GM5+ basis / M 点 k-active 已按官网对齐 |
| Types 改了结果不变 | 忘记点 **"Change"** |
| ZIP 报没有子群 | 先对下拉框里选中的 Method 1/2/3 点 OK |
| Method 2 参数 k 点无子群 | 勾选 **"Generate isotropy subgroups database if missing"** 后重试（可能极慢） |
| 端口被占用 | 改 `web_port` 或关掉旧的 `main_web.py` |
| OneDrive 路径偶发文件锁 | 可拷到非同步本地盘再试 |

---

## 13. 目录结构（本项目）

```text
ISODISTORT/
├── main_web.py / main_terminal.py / main_requirement.py
├── web/                 网页（server.py + index.html + static/）
├── isocore/             计算核心（api / backend / structure / distortion / io / superspace）
├── config/settings.yaml
├── isobyu/              【只读】ISOTROPY 套件
├── tests_dev/           开发测试（pytest）+ manual/ 手工脚本
└── output/              运行产物（不入库）
```

官网概念帮助：[isodistorthelp.php](https://iso.byu.edu/isodistorthelp.php)。
