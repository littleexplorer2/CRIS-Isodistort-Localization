# ISODISTORT 本地化

本目录把 [ISODISTORT](https://iso.byu.edu/isodistort.php)（BYU 晶体畸变在线计算）做到离线可用：从**母相结构文件**出发，在本地得到**子群结构信息或结构文件**。计算后端是 ISOTROPY Suite 的 `iso` Linux 二进制与 `data_*.txt` 数据库。

同一套引擎提供三种用法：**网页**、**终端菜单**、**Python API**。网页与终端界面为**英语**，没有语言切换。

本项目的目标工作流是：上传母相 CIF → 设定 distortion types → 在 Method 1 / 2 / 3 / 4 点 OK 得到结果表（可筛选、排序）→ 在 Distortion 区下载**筛选后的结果表**（1–4）以及按**一个 Method** 导出的**筛选后**子群结构文件 ZIP（仅 1 / 2 / 3：CIF / ISOVIZ / Complete modes details / TOPAS.STR）。官网 Distortion 页上的 **Generate**（按模式幅度生成畸变结构）和 **Domains**（畴列表）**不参与**“得到子群结构信息或文件”，已从网页和终端中删除。Python API 仍保留 `generate_distortion` / `generate_domains`，供脚本使用。

---

## 1. 仓库边界（哪些可以改）

| 路径 | 可否修改 | 说明 |
| --- | --- | --- |
| `ISODISTORT/`（不含下表例外） | 可以 | 本地化主体：计算、网页、终端、API、文档、测试 |
| `ISODISTORT/isobyu/` | **不可以** | 从 [iso.byu.edu](https://iso.byu.edu/isotropy.php) 下载的 Linux 二进制与数据库，只读 |
| `ISODISTORT_VALIDATE/` | 可以 | 比较本地 CIF 与官网 CIF |
| `ISOVIZ_INPUT/` | 可以 | 把振幅 CSV 写入 `.isoviz` 并打开 IsoVIZ；`data.csv/` 与 `subgroup.isoviz/` 不入库 |
| `webpage_info/` | **不可以** | 官网各步 HTML 存档 |
| `experiment_data/` | **不可以** | 实验母相 CIF 等原始数据 |
| `GD/` | **不可以** | 梯度下降拟合代码与笔记本 |

官网交互对照（`webpage_info/`，勿改其中文件）：在 1 首页上传 `experiment_data/EuAl4 Parent.cif` → 2 勾选 strains 与 displacive（Eu, Al），Method 2 选 LD、K10 (0,0,g)、g=1/6 → 连续 OK 进入 3、4、5 → 6 为导出页（CIF / Save interactive distortion / Complete modes details / TOPAS.STR）。

---

## 2. 从 GitHub 下载与部署

仓库地址：

```text
https://github.com/littleexplorer2/CRIS-Isodistort-Localization.git
```

### 2.1 克隆

需要已安装 [Git](https://git-scm.com/)。在要放置项目的目录打开终端：

```powershell
git clone https://github.com/littleexplorer2/CRIS-Isodistort-Localization.git
cd CRIS-Isodistort-Localization
```

若 GitHub 要求登录，按提示使用 HTTPS 凭据或 SSH 密钥。克隆失败时检查：网络能否访问 github.com、仓库是否为私有、磁盘空间是否足够。

### 2.2 放入 ISOTROPY 套件（必需）

克隆下来的仓库**不包含**可运行的 `iso` 二进制。请从 [ISOTROPY Suite](https://iso.byu.edu/isotropy.php) 下载 Linux 版，把 `iso`、`smodes` 与全部 `data_*.txt` 放进：

```text
ISODISTORT/isobyu/
```

`python ISODISTORT/main_requirement.py` 若发现没有 `isobyu/`，会新建空目录并打印上述下载地址；**不会**自动下载套件。没有 `ISODISTORT/output/` 时会自动创建。

不要改该目录里的文件内容。缺 `iso` 或数据库时，网页/终端在真正计算时会报错。

### 2.3 目录是否完整

克隆并放入 `isobyu` 后，至少应能看到：

- `ISODISTORT/main_web.py`、`main_terminal.py`、`main_requirement.py`
- `ISODISTORT/config/settings.yaml`
- `ISODISTORT/web/index.html`
- `ISODISTORT_VALIDATE/main.py`
- `ISODISTORT/isobyu/iso` 以及若干 `data_*.txt`

---

## 3. 运行前检查：依赖与环境

### 3.1 操作系统

- **Windows 10/11 + WSL**（推荐开发环境）：`isobyu/iso` 是 Linux ELF，必须通过 WSL 调用。
- **原生 Linux**：可直接运行 `iso`，无需 WSL。

在 Windows 上检查 WSL：

```powershell
wsl --status
wsl --list --verbose
```

需要有一个**默认 Linux 发行版**且能启动。若提示未安装 WSL：按 Microsoft 文档启用“适用于 Linux 的 Windows 子系统”，安装 Ubuntu 等发行版，重启后再试。`wsl -e uname -a` 应打印 Linux 内核信息。

### 3.2 Python

需要 **Python 3.10 或更高**（开发在 3.12.5 上验证过）。检查：

```powershell
python --version
```

若命令找不到，安装 Python 并勾选 “Add python.exe to PATH”，然后新开一个终端。不要用 Windows 商店里过旧的 Python 2。

### 3.3 本机还缺什么时不要直接点网页

| 检查项 | 怎样确认 | 缺了会怎样 |
| --- | --- | --- |
| Python ≥ 3.10 | `python --version` | 建虚拟环境失败 |
| WSL（仅 Windows） | `wsl --status` | 调用 `iso` 失败 |
| `ISODISTORT/isobyu/iso` | 资源管理器中能看到该文件 | Method OK 时报底层错误 |
| `data_*.txt` 在 `isobyu/` | 同目录有多个 data 文件 | 子群枚举为空或报数据库错误 |
| 磁盘可写 `ISODISTORT/output/` | 程序会自动建目录 | 上传 CIF / 终端导出失败 |

---

## 4. 安装 Python 依赖

在**仓库根目录**只建一份虚拟环境 `CRIS/.venv`（覆盖 ISODISTORT、ISODISTORT_VALIDATE 与 ISOVIZ_INPUT）：

```powershell
cd <你克隆下来的仓库根目录>
python ISODISTORT\main_requirement.py
```

开发（pytest / ruff）再加：

```powershell
python ISODISTORT\main_requirement.py --dev
```

已有 `.venv` 则复用；`--recreate` 会删掉后重建。脚本会检查 Python、WSL、`ISODISTORT/output/`（没有则新建）、`isobyu/`（没有则新建空目录并提醒去官网下载），并只安装尚未存在的 `requirements.txt` 依赖。**不会**自动下载 ISOTROPY 套件。

成功后控制台会打印建议命令，例如用 `.venv\Scripts\python.exe` 启动网页或终端。以后运行请始终用该解释器，避免装到了系统 Python 里却用另一个 Python 启动。

激活虚拟环境（可选）：

```powershell
.\.venv\Scripts\Activate.ps1
```

若 PowerShell 禁止执行脚本，可改用上面的 `.\.venv\Scripts\python.exe` 全路径，不必激活。

---

## 5. 网页版：启动与整页功能说明

### 5.1 启动

```powershell
cd <仓库根目录>
.\.venv\Scripts\python.exe ISODISTORT\main_web.py
```

默认打开 `http://127.0.0.1:8000/`。端口在 `ISODISTORT/config/settings.yaml` 的 `runtime.web_port`；被占用时自动顺延约 20 个端口，再不行则让系统分配。控制台会打印最终 URL。浏览器未自动打开时，把该 URL 粘贴到地址栏。

- 右上角 **Stop**：停止服务并释放端口。
- 关闭标签页后约 **60 秒**（`runtime.web_idle_timeout`）自动停服。要常驻就不要关页。
- 界面语言固定为英语，没有 Language 下拉菜单。

### 5.2 顶栏

- **ISODISTORT**：回到本页（本地搜索页）。
- **SUITE / HELP**：打开官网套件首页与帮助（新标签）。需要联网；计算本身不依赖这两项。

### 5.3 Parent CIF（母相）

- 选择 `.cif` / `.CIF` 文件，点 **Load**。
- 页头会显示空间群、晶格、Wyckoff、Default space-group preferences。
- 必须先 Load 再跑任何 Method；否则 OK 会提示先加载母相。
- 实验示例（若仓库中仍有该文件）：`experiment_data/EuAl4 Parent.cif`。不要修改该目录。

### 5.4 Types of distortions to be considered

- **Strain** 复选框：是否考虑应变。
- **Displacive / Occupational / Magnetic / Rotational**：每行 **all / none / 各物种**。复选框互不联动。
- 默认 Strain + Displacive 全物种（与官网第 2 页一致）。
- **必须点 Change** 才会提交；只勾选不点 Change，后面的 Method 仍用旧设置。

### 5.5 Method 1: Search over all special k points

- **Crystal system(s)**：多选，逻辑为 OR；全不选表示不过滤晶系。
- **Space-group symmetry**：可达子群空间群下拉；空为不过滤。
- **Conventional lattice / Primitive lattice**：互斥过滤超胞基矢。本地 iso 数据库（9.6.1）与官网不完全相同，下拉选项可能略有差别；页上有英文说明。
- **Maximal subgroups only**：只保留极大子群。
- **OK**：枚举全部特殊 k 点上的各向同性子群。首次可能数秒到数十秒，属正常。
- 结果与 Method 2 / 3 一样是**表格**（`idx` / `SG` / `k` / `Irrep` / `OPD` / `crystal system` / `maximal`），可 Filter 与 ▲/▼，而不是官网 Order parameter 页那种一行一条的 radio 文本。
- **筛选**：Filter 输入框对 `SG` / `k point` / `Irrep` / `OPD` / `crystal system` / `maximal subgroup` 做不区分大小写的子串匹配。可勾选 **Show filtered rows only** 只显示命中行。`idx` 不参与筛选（那是子群编号）。
- **排序**：`SG`、`k point`、`Irrep`、`OPD`、`crystal system`、`maximal subgroup` 表头右侧各有 **▲（升序）** 和 **▼（降序）**。每次只按当前点的那一列、那一个方向排；点击行计算模式时按子群 `idx`，不会因排序错位。
  - `SG`：按空间群序号数字排。
  - `k point` / `Irrep` / `OPD`：按标签字符串排（含数字时按自然序）。
  - `crystal system`：按三斜→单斜→正交→四方→三方→六方→立方。
  - `maximal subgroup`：升序为非极大在前，降序为极大在前。
- **点击一行**：按该子群计算模式基矢，结果显示在 Method 2 下方的 modes 区。
- Method 区**没有**下载按钮。筛选后的表、以及全部子群结构文件，都在页面底部 **Distortion** 下载。

### 5.6 Method 2: General method — specific k points

- **Specify k point**：从母相允许的 k 点列表选择。
- 参数 k 点（如 LD 的 a/b/g）会出现数值框，必须填齐再 OK。
- **Generate isotropy subgroups database if missing**：勾选后，若当前 k 点在本地 `data_*.txt` **没有**现成的各向同性子群列表，就现场生成这份列表（官网 “Generate isotropy subgroups”）。首次可能较慢，之后缓存。勾选框后有一行英文说明。
- **OK**：枚举该（组）k 点全部 IR 的子群。
- 若枚举为空：可勾选上面的数据库生成后重试，或按提示去官网。
- 结果表：Filter（SG / k / Irrep / OPD / s / i）、**▲/▼ 排序**、**Show filtered rows only**。点击行查看模式基矢。表旁**不再**提供 txt/csv 下载。
- **Change number of superposed IRs**：改叠加 IR 组数后必须点 **Change**，才会出现多组 k vector。
- 本地**没有**官网的 “# of independent incommensurate modulations”（nmod）；页上有英文说明。
- 参数 k 点子群可以列出，但本地不能算位移模式；点行会看到英文警告。筛选/排序后的表仍可在 Distortion 下载。

### 5.7 Method 3: space/point group + lattice

- 选 **space group** 或 **point group**（晶类）。
- **direct**：实空间子格 + Default/P/A/B/C/I/F/R 带心。本地 Method 3 **只支持 Default (`d`)**；选其它带心时后端会按默认处理或提示。
- **reciprocal**：倒易超格。**本地不支持**，请用 direct。
- 3×3 **basis** 输入 a'、b'、c'。
- **OK** 后得到候选表（`idx`、`SG`、`k point`、`Irrep`、`point group`）。与 Method 1 一样可 **Filter**、**▲/▼ 排序**、**Show filtered rows only**；点击行计算该子群模式。下载在 Distortion。

### 5.8 Method 4: Mode decomposition

- 上传**已经畸变**的女儿相 CIF，点 **OK**。
- 输出**全部**模式幅度（可 Filter / ▲▼ 排序）和 RMS residual。表可滚动，不再只显示绝对值前 20。
- 这是分解，不是子群列表；**不能**作为 Distortion **ZIP** 的 Method 来源，但可以在 Distortion 下载该幅度表的 txt/csv。

### 5.9 Distortion（统一下载）

四个 Method 面板只负责**调参数和计算**。本区负责**下载**：筛选后的结果表（Method 1–4）以及子群结构文件 ZIP（Method 1–3）。没有 Generate / Domains，也没有模式幅度输入框。

- **Export source**：下拉选 **一个** Method（1 / 2 / 3 / 4），默认 Method 2。不能多选。
- **Filtered result table**：按所选 Method **当前**的 Filter 与排序，下载全部命中行的 txt 或 csv（不受 “Show filtered rows only” 限制：只要筛中就会写入文件）。须先对该 Method 点过 OK。
- **Download formats**（可多选，对应官网第 6 页；仅用于下面的 ZIP）：
  - CIF file
  - Save interactive distortion（`.isoviz`）
  - Complete modes details（`.txt`）
  - TOPAS.STR
- **Download all (ZIP)**：打包 **Method 1 / 2 / 3** 当前筛选命中的子群结构文件（无筛选则全部），**只含勾选的格式**，**不扫描** `output/`。Method 4 没有子群 CIF，选 Method 4 再点 ZIP 会提示改用表格下载。
- 请先对该 Method 点过 OK。否则 ZIP 会提示没有子群，表格下载会提示没有结果表。
- 压缩包名 `isodistort_methodN.zip`，内部：

```text
isodistort_method2/
  LD1 C1/
    LD1 C1 CIF.cif
    LD1 C1 Save interactive distortion.isoviz
    LD1 C1 Complete modes details.txt
    LD1 C1 TOPAS.STR
```

- CIF 一般为该子群**零振幅超胞**。参数 k 点上本地算不出位移模式时，模式类文件会注明限制。
- 单行点选后出现的模式表在 Method 2 区域，仅用于查看，不再用于 Generate。

### 5.10 Space-Group Preferences

只读。本地 `iso` 固定国际标准取位（单斜 a(b)c、cell choice 1、正交 abc、三方 hexagonal、原点 2、超空间 standard）。不能在本地改这些选项。

### 5.11 为何删除 Generate 与 Domains

- **Generate**：把用户填写的模式幅度作用到结构上，得到**有限振幅畸变 CIF**。本项目要的是子群列表及其结构文件；Method OK + ZIP 已给出零振幅超胞 CIF 与模式文件，不需要再点 Generate。
- **Domains**：列出畴变体（指数、生成元、空间群），**不写出结构文件**，对“得到子群结构信息或文件”没有帮助。
- 不要把 Method 2 的 **Generate isotropy subgroups database if missing** 与已删除的 Distortion **Generate** 搞混：前者是补参数 k 点子群库，后者才是按幅度生成畸变结构。

---

## 6. 终端版：启动与菜单说明

```powershell
cd <仓库根目录>
.\.venv\Scripts\python.exe ISODISTORT\main_terminal.py
```

界面为英语。方括号里的值是默认值，直接回车即采用。没有“切换语言”菜单项。

启动后会先要求选择**母相 CIF**：

- 在 `ISODISTORT/` 下最多列出 30 个 `.cif`（不含 `output/` 与 `output/tmp/`）。
- 仓库外的文件（如 `experiment_data/EuAl4 Parent.cif`）请选手动输入路径，建议绝对路径。

### 6.1 主菜单（Search Page）

- **1. Reload parent CIF**：重新加载母相；会清空上次 Method 1/2/3/4 结果。
- **2. Set distortion types**：输入类型编号或名称（逗号分隔，如 `1,3`）。再为非 strain 类型指定物种作用域：`all` / `none` / `Eu,Al`。
- **3. Method 1**：可选晶系（逗号分隔，如 `tetragonal`）、可选子群空间群号、是否只要极大子群、Conventional/Primitive lattice（输入如 `C3`、`P2`，`0` 表示不选）。计算后进入与网页相同的表交互：
  - `f <col>=text` 筛选（列：`sg` / `k` / `irrep` / `opd` / `cs` / `max`）
  - `s <col> a|d` 升序或降序
  - `c` 清除筛选；`only` 切换是否只显示命中行
  - 输入子群 `idx` 计算模式；`q` 结束（下载到菜单 7）
- **4. Method 2**：先打印 “Generate isotropy subgroups database if missing” 的说明，再询问是否生成缺失库。选择 superposed IR 组数与各组 k 点；参数 k 点会按顺序要 a/b/g。子群表命令与 Method 1 相同（列：`sg` / `k` / `irrep` / `opd` / `s` / `i`）。参数 k 点只枚举子群、不算位移模式。
- **5. Method 3**：点群或空间群、direct/reciprocal（reciprocal 会改回 direct）、带心（本地仅 `d`）、可选 3×3 基矢。结果表命令同上（列含 `pg`）。
- **6. Method 4**：再选一个女儿 CIF，列出**全部**模式幅度及 RMS，可用 `f` / `s` / `c` / `only`（无 idx）。
- **7. Distortion**：见下一小节。
- **8. Show current state**：母相空间群、types、scope、各 Method 条数、已映射模式数。
- **0. Exit**。

### 6.2 菜单 7：Distortion（下载结果）

对应网页 Distortion：先在 Search Page 算完 Method，再到这里下载。结构文件写入目录而不是浏览器 ZIP。

- **1 / 2 / 3**：分别导出 Method 1、2、3 的子群结构文件（必须先跑过对应 Method）。若该 Method 表上有 Filter，则只导出命中的子群。
- **4**：按所选 Method（1–4）**当前**筛选与排序，写出结果表 `isodistort_methodN_filtered.txt` 或 `.csv`。
- 结构文件格式提示：默认 `cif,isoviz,modes,topas`。
- 结构文件输出目录：默认 `ISODISTORT/output/isodistort_methodN/`。
- 每个子群一个文件夹，文件名规则与网页 ZIP 相同。
- **没有**单模式 Generate、混合 Generate、导出当前畸变结构、Domains。

---

## 7. Python API（简要）

```python
from isocore.api import IsoDistort

iso = IsoDistort()  # language= 参数仍可写，但会被忽略（仅英语）
iso.load_structure(r"C:\path\to\parent.cif")
iso.set_distortion_scope({
    "displacive": ["*"], "occupational": [], "strain": [],
    "magnetic": [], "rotational": [],
})
iso.set_distortion_types(["strain", "displacive"])

m1 = iso.search_method_1(crystal_system="tetragonal")
iso.export_subgroups("out_batch", formats=["cif", "isoviz", "modes", "topas"])
```

`export_subgroups` 使用当前 `iso.subgroups`（最近一次写入的列表）。网页 ZIP 会显式传入所选 Method 的列表，避免三种 Method 混在一起。

`generate_distortion` / `generate_mixed_distortion` / `generate_domains` 仍在 API 中，网页和终端不再调用。

---

## 8. 路径与配置文件

文件：`ISODISTORT/config/settings.yaml`。其中的相对路径相对 `config/` 解析。

| 键 | 含义 | 默认 |
| --- | --- | --- |
| `isobyu.bin_dir` / `data_dir` | `iso` 与数据库目录 | `../isobyu` |
| `isobyu.iso_bin` 等 | 可执行文件名 | `iso` / `findsym` / `smodes` |
| `defaults.position_tolerance` | 分数坐标容差 | `0.001` |
| `defaults.lattice_tolerance` | 晶格容差 | `0.00001` |
| `defaults.default_amplitude` | API 生成畸变时的默认幅度 | `1.0` |
| `defaults.default_supercell` | 未用子群基矢时的超胞 | `[1,1,1]` |
| `runtime.web_port` | 网页首选端口 | `8000` |
| `runtime.web_idle_timeout` | 关页后自动停服秒数 | `60` |
| `runtime.temp_dir` | 上传暂存 | `../output/tmp` |
| `runtime.output_dir` | 终端导出等 | `../output` |
| `runtime.timeout` | 普通 `iso` 调用超时 | `60` |
| `runtime.generation_timeout` | Method 2 勾选生成缺失子群库时的超时 | `3600` |

已删除 `runtime.language`：界面固定英语。网页 ZIP **不写入、不读取** `output/`。网页上的筛选表 txt/csv 在浏览器里生成；终端菜单 7 选项 4 才把表写到 `output/`。

---

## 9. 目录结构

```text
ISODISTORT/
├── main_web.py / main_terminal.py / main_requirement.py
├── web/                 网页（server.py + index.html + static/）
├── isocore/             计算核心
│   ├── api/             IsoDistort
│   ├── backend/         iso / smodes / findsym 封装
│   ├── structure/       CIF 与超胞
│   ├── distortion/      Method 1–4、模式、畴（畴仅 API）
│   ├── io/              CIF / ISOVIZ / modes details / TOPAS.STR / POSCAR
│   ├── i18n/            英语界面字符串（messages.py）
│   └── utils/           配置、解析
├── isodistort/          包名别名（重导出 isocore）
├── config/settings.yaml
├── isobyu/              【只读】ISOTROPY 套件
├── tests_dev/           开发测试（生产可删）
└── output/              运行产物（不入库）
```

---

## 10. 测试与 CIF 验证

```powershell
python ISODISTORT\main_requirement.py --dev
cd ISODISTORT
.\.venv\Scripts\python.exe -m pytest tests_dev -q
ruff check .
```

若当前目录已是 `ISODISTORT`，pytest 用该目录下的解释器即可。依赖 WSL 的用例在 WSL 不可用时会跳过。

本地 CIF 与官网 CIF 的语义比较见仓库内 `ISODISTORT_VALIDATE/`（该目录 README）。网页 ZIP 解压后的 `LD1 C1/LD1 C1 CIF.cif` 可与官网第 6 页 CIF 配对。

---

## 11. 与官网的已知差异

1. **Types 过滤时机**：参数 k 点按 smodes + 物种过滤；Method 1 在枚举后按 smodes 活性 IR + 默认排除 magnetic（`m*`）收缩列表，使 Order parameter 行与官网一致。
2. **振幅**：官网 As/Ap + normfactor；本地 API 生成畸变为「最大模长归一化位移 × 用户幅度」。网页已不再提供幅度输入。
3. **模式列表**：官网含次级模式；本地 DISPLAY BUSH 以主 OPD 位移模式为主。
4. **应变模式未实现**：不改晶格参数。
5. **参数 k 点**：可枚举子群；位移模式依赖官网 (3+d) 超空间，本地会提示。nmod 已从网页/终端移除。
6. **Method 3**：只覆盖特殊 k 点；reciprocal 不支持；带心仅默认 `d`。
7. **Method 4**：需先有模式基矢；女儿相为超胞时会提升坐标再分解。
8. **magnetic**：带 `m` 前缀的 IR 不进入默认流程。
9. **occupational**：本地为 ±1 占据近似，校验失败会标明。
10. **Method 1 lattice 下拉**：由本地枚举基矢去重生成，数量可能与官网略有差别。
11. **多物种共享 Wyckoff 字母**：BUSH 无法分物种时不静默给错结构。
12. **多维模式**：各分量等权求和再归一化。
13. **取位约定**：CIF 位点与 iso 符号点对不齐时会标明。
14. **界面语言**：仅英语；官网为英语，本地不再提供中文 UI。
15. **Distortion Generate / Domains**：官网有；本地网页/终端已去掉。子群结构文件由 Distortion ZIP/目录导出；结果表 txt/csv 也在 Distortion 下载。

ISOVIZ / Complete modes details / TOPAS.STR 按官网选项语义在本地生成，不是官网站点 PHP 的逐字节拷贝。

---

## 12. 运行过程中可能遇到的问题

- **`wsl` 不是内部或外部命令 / WSL 报错**：安装并设置默认发行版；`wsl --status` 应成功。
- **`iso` 找不到或 Permission denied**：确认 `ISODISTORT/isobyu/iso` 存在且为 Linux 二进制；不要把 Windows 可执行文件放进去。
- **Method 1 很慢**：首次枚举全部特殊 k 点，属正常；同一次会话内再跑会快一些。
- **网页打不开**：看控制台打印的 URL；检查防火墙是否拦了本机 8000；换 `web_port`。
- **关页即停服**：预期行为。点 Stop 或关标签后等待超时。
- **Load 后 Method OK 仍说没有 CIF**：确认 Load 成功（绿色提示），刷新状态栏应显示空间群。
- **Types 改了但结果没变**：忘记点 **Change**。
- **Method 1 表头排序后点行错了**：`idx` 列是子群编号，不是当前行号；请点要算的那一行的 `idx`，不要按屏幕上的第几行数。
- **ZIP 报没有子群**：先跑下拉框里选中的那一个 Method（1/2/3），再 Download all。
- **表格下载报没有结果**：先对该 Method 点 OK；Method 4 只能下 txt/csv，不能下 ZIP。
- **ZIP 里文件特别多、名字对不上子群**：旧版曾打包整个 `output/`；当前版本只含所选 Method 的子群文件夹（且不受表上 Filter 限制）。
- **想只要筛过的子群列表**：用 Distortion 的 Download filtered (txt/csv)，不要用 ZIP。
- **Method 2 参数 k 点没有子群**：勾选 “Generate isotropy subgroups database if missing” 后重试，或去官网；生成可能极慢，超时见 `generation_timeout`。勾选框下方的英文会说明这是在补本地子群库。
- **参数 k 点点行没有模式**：本地限制，不是排序或筛选错误。
- **Method 3 reciprocal / 非 Default 带心**：本地不支持或仅 `d`；改用 direct + Default。
- **终端找不到实验 CIF**：不在 `ISODISTORT/` 内时用手动路径。
- **端口被占用**：改 `web_port` 或关掉占用 8000 的旧 `main_web.py`。
- **OneDrive 路径**：仓库若在 OneDrive 下，偶发文件锁；计算失败时可把仓库拷到本地非同步盘再试。
- **中文路径 / 控制台乱码**：终端启动时会尝试把 stdout 设为 UTF-8；仍乱码可在 Windows 终端使用 UTF-8 代码页。界面本身是英文。
- **pytest 收集失败**：用 `--dev` 装过依赖，并在 `ISODISTORT/` 下用同一 `.venv` 运行。

---

## 13. 术语

英文以官网与 ITC 为准。界面不再翻译成中文。官网帮助：[isodistorthelp.php](https://iso.byu.edu/isodistorthelp.php)。
