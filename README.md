# CRIS

CRIS 是一个**本地化晶体学工具集合**（monorepo）。你不需要先学晶体学：可以把这里理解成「从母相晶体结构出发，在本地搜索可能的畸变子结构、核对结果是否与官网一致、再把拟合得到的振幅写进可视化文件」的一整套流水线。

仓库地址示例：

```text
https://github.com/littleexplorer2/CRIS-Isodistort-Localization.git
```

底层计算依赖 [ISOTROPY Suite](https://iso.byu.edu/isotropy.php)（尤其是其中的 ISODISTORT / IsoVIZ 相关工具）。本仓库把常用流程做成可离线使用的 Python 程序；网页与终端界面为**英语**。

---

## 三个子项目分别做什么

| 子项目 | 一句话 | 详细说明 |
| --- | --- | --- |
| **ISODISTORT/** | 本地「子群搜索 + 导出结构」 | 上传母相 CIF → 勾选畸变类型 → Method 1–4 搜索/分解 → 在 Distortion 下载结果表与子群结构 ZIP。Method 2 的 nmod 为 (3+d) 超空间附加维度 d（0=公度，1–3=本地内核，IT-C）；网页页底 (3+d) 板块是本地检查台，官网没有。见 [ISODISTORT/README.md](ISODISTORT/README.md) |
| **ISODISTORT_VALIDATE/** | 核对本地 CIF 是否算对 | 把本地 CIF 放入 `compare/item/`，把官网参考 CIF 放入 `compare/true/`（批量比较须改名一一对应），用 `main.py` 比较并输出 PASS/FAIL。见 [ISODISTORT_VALIDATE/README.md](ISODISTORT_VALIDATE/README.md) |
| **ISOVIZ_INPUT/** | 把振幅 CSV 写入 `.isoviz` 并启动 IsoVIZ | 从 `input_content/` 读取 CSV 与子群 `.isoviz`，写入 `amp` 后自动打开 Java 版 IsoVIZ（不使用 `output/`）。见 [ISOVIZ_INPUT/README.md](ISOVIZ_INPUT/README.md) |

各子项目 README **只讲该项目本身**。跨项目怎么串起来，只在本文件说明。

---

## 典型工作流（三者如何配合）

下面是一条常见路径（你也可以只用其中一段）：

1. **准备母相结构**  
   例如实验母相 CIF（仓库里若有 `experiment_data/`，那是只读原始数据，不要改）。

2. **用 ISODISTORT 搜索子群并导出**  
   - 启动网页（`ISODISTORT/main_web.py`）或终端（`ISODISTORT/main_terminal.py`）。  
   - 加载 Parent CIF → 设置 Types of distortions → 对 Method 1 / 2 / 3 点 OK → 在 Distortion 下载筛选后的表，以及 CIF / `.isoviz` / modes / TOPAS 等 ZIP。  
   - Method 4 是把「已经畸变」的女儿相 CIF 分解成模式幅度，不生成子群 ZIP。

3. **（可选）用 ISODISTORT_VALIDATE 核对**  
   把本地 ZIP 里某个 `subgroup.cif`（或旧版 `… CIF.cif`）拷到 `ISODISTORT_VALIDATE/compare/item/`，把官网第 6 页导出的同名子群 CIF（官网常为 `subgroup.cif`）拷到 `compare/true/`。**批量比较时必须把 `true/` 里官网下载的文件改成与 `item/` 完全相同的相对路径和文件名。** 然后运行 `python ISODISTORT_VALIDATE/main.py`（菜单或 `compare` / `batch` 子命令）。不要再传入自定义路径。`compare/` 整目录不入库，由 `main_requirement.py` 在缺失时自动创建。

4. **（可选）用 ISOVIZ_INPUT 看拟合振幅**  
   把振幅 CSV 放入 `ISOVIZ_INPUT/input_content/data.csv/`，把对应子群 `.isoviz` 放入 `input_content/subgroup.isoviz/`（整个 `input_content/` 不入库，缺失时由安装脚本自动创建）。运行 `ISOVIZ_INPUT/main.py` 读取输入后会**直接启动 IsoVIZ**，本子项目不使用 `output/` 写出结果。

```text
母相 CIF
   │
   ▼
ISODISTORT  ──►  子群表 / CIF / .isoviz / modes / TOPAS
   │                    │
   │                    ├──► ISODISTORT_VALIDATE（与官网 CIF 比对）
   │                    │
   │                    └──► ISOVIZ_INPUT（CSV 振幅 → 启动 IsoVIZ）
```

---

## 绝对不要修改的目录

这些目录是原始数据、官网存档或第三方二进制，**请只读使用**：

| 路径 | 原因 |
| --- | --- |
| `experiment_data/` | 实验母相 CIF 等原始数据 |
| `GD/` | 梯度下降拟合代码与笔记本（本流水线的上游数据来源之一，但不是本仓库要改的部分） |
| `webpage_info/` | 官网各步 HTML 存档（对照交互顺序用） |
| `ISODISTORT/isobyu/` | 从 iso.byu.edu 下载的 Linux 二进制（`iso`、`smodes` 等）与 `data_*.txt` 数据库 |

可以改的主体是：`ISODISTORT/`（除 `isobyu/`）、`ISODISTORT_VALIDATE/`、`ISOVIZ_INPUT/`。

---

## 共享虚拟环境

三个子项目共用仓库根目录下的一份虚拟环境：

```text
CRIS/.venv
```

在根目录任选其一安装即可（都会创建/复用 `.venv`）：

```powershell
cd <CRIS 根目录>
python ISODISTORT\main_requirement.py
```

或：

```powershell
python ISODISTORT_VALIDATE\main_requirement.py
```

或：

```powershell
python ISOVIZ_INPUT\main_requirement.py
```

之后请始终用 `.venv\Scripts\python.exe` 运行各子项目入口，避免装到系统 Python。各子项目 README 里写有本项目的启动命令。

**若运行卡死很久且终端没有出现 `(.venv)`：** 通常是虚拟环境未加载成功。关闭当前进程后，用 `.\.venv\Scripts\python.exe …`（或先 `.\.venv\Scripts\Activate.ps1`）重新运行即可。

---

## 用户需要配置 / 放置的路径（总表）

多数路径有合理默认值；**换机器或换安装位置时**请按下面核对。细节以各子项目 README 的同名小节为准。

| 用途 | 默认位置 / 做法 | 何时需要改 |
| --- | --- | --- |
| **ISOTROPY Linux 二进制与 `data_*.txt`** | 放入 `ISODISTORT/isobyu/` | 首次安装必做；若放别处则改 `ISODISTORT/config/settings.yaml` 的 `isobyu.bin_dir` / `data_dir` |
| **母相 CIF（计算输入）** | 网页上传，或终端选择/粘贴路径；示例只读目录 `experiment_data/` | 每次计算指定你的 CIF，**不要改** `experiment_data/` 里的原始文件 |
| **临时目录 / 终端导出目录** | `ISODISTORT/output/tmp`、`ISODISTORT/output`（相对 `config/` 写在 `settings.yaml`） | 仅当磁盘空间或策略要求换盘时改 `runtime.temp_dir` / `output_dir` |
| **VALIDATE 成对 CIF** | 固定目录 `ISODISTORT_VALIDATE/compare/item/` 与 `compare/true/` | **不要改程序路径**；把文件放进这两处，并使相对路径/文件名一一对应 |
| **ISOVIZ 振幅 CSV + 子群 `.isoviz`** | `ISOVIZ_INPUT/input_content/data.csv/` 与 `…/subgroup.isoviz/` | 把日常输入放进这里，或运行时用 `--data` / `--structure` 指定任意路径 |
| **IsoVIZ 程序** | 根目录 `ISOViz.lnk`（已 gitignore），或根目录 `.jar`/`.exe`，或环境变量 `ISOVIZ` / `ISOVIZ_JAR` | 本机安装 IsoVIZ 后任选一种方式让程序能找到它；还需 **Java** 在 `PATH` 中 |
| **VESTA（人工打开 CIF）** | 可选：根目录放 `VESTA.lnk`（已 gitignore）指向本机 `VESTA.exe` | 仅手工抽检用；程序不强制读该快捷方式 |

### 外部可视化工具：VESTA 与 IsoVIZ

本仓库**不附带**下列程序，需自行安装；用于打开 ISODISTORT 导出的结构文件做人工抽检。

| 工具 | 作用（在本项目中） | 下载 |
| --- | --- | --- |
| **[VESTA](https://jp-minerals.org/vesta/en/)** | 三维晶体结构可视化。用它打开 Distortion 导出的 **`subgroup.cif`**，检查晶胞、原子位点与对称是否合理。免费（学术/非商业等，以官网许可为准）。 | [下载页](https://jp-minerals.org/vesta/en/download.html)（Windows 常用 `VESTA-win64.zip`，解压后运行 `VESTA.exe`） |
| **IsoVIZ**（ISOTROPY Suite） | 交互查看畸变模式振幅。用它打开 **`data.isoviz`**。 | 随 [ISOTROPY Suite](https://iso.byu.edu/isotropy.php) 安装；启动方式见 [ISOVIZ_INPUT/README.md](ISOVIZ_INPUT/README.md) |

建议在仓库根目录放置本机快捷方式（已在 `.gitignore`，不会入库）：

- `VESTA.lnk` → 指向解压后的 `VESTA.exe`
- `ISOViz.lnk` → 指向 IsoVIZ 启动器 / `.jar`

双击快捷方式或「打开方式」即可检查导出文件；**计算流程不依赖**这两个快捷方式是否存在。
| **WSL / `ISODATA`** | 由封装自动处理（短路径暂存 + 符号链接） | 一般**不必**改系统环境变量；只需装好默认 WSL 发行版 |
| **网页端口** | `settings.yaml` → `runtime.web_port`（默认 `8000`） | 端口被占用时修改（不是文件路径，但属本机配置） |

跨项目流水线仍见上文「典型工作流」。

---

## 仓库地图（还可能看到什么）

```text
CRIS/
├── README.md                 ← 本文件：总览与跨项目关系
├── agent.md                  ← AI Agent 工作指南（修改边界、思考方式、验证清单）
├── .venv/                    ← 共享虚拟环境（不入库）
├── ISODISTORT/               ← 本地 ISODISTORT（网页 / 终端 / API）
├── ISODISTORT_VALIDATE/      ← CIF 语义比较（compare/item vs compare/true）
├── ISOVIZ_INPUT/             ← CSV → 启动 IsoVIZ（input_content/ 不入库）
├── experiment_data/          ← 【勿改】实验数据
├── webpage_info/             ← 【勿改】官网 HTML 存档
├── GD/                       ← 【勿改】梯度下降相关
├── ISOViz.lnk                ← （本机）IsoVIZ 快捷方式，已在 .gitignore
└── VESTA.lnk                 ← （可选，本机）VESTA 快捷方式，已在 .gitignore
```

---

## 从哪里开始读

| 你想做的事 | 打开 |
| --- | --- |
| AI / Agent 修改与验收本仓库 | [agent.md](agent.md) |
| 安装环境、跑网页/终端搜索子群 | [ISODISTORT/README.md](ISODISTORT/README.md) |
| 比较本地 CIF 与官网 CIF | 放入 [ISODISTORT_VALIDATE/compare/](ISODISTORT_VALIDATE/compare/) 后运行 `ISODISTORT_VALIDATE/main.py`，见 [ISODISTORT_VALIDATE/README.md](ISODISTORT_VALIDATE/README.md) |
| 把振幅写入 IsoVIZ 并打开 | 放入 [ISOVIZ_INPUT/input_content/](ISOVIZ_INPUT/input_content/) 后见 [ISOVIZ_INPUT/README.md](ISOVIZ_INPUT/README.md) |

官网帮助（概念背景，非本仓库文档）：[ISODISTORT Help](https://iso.byu.edu/isodistorthelp.php)、[ISOTROPY Suite](https://iso.byu.edu/isotropy.php)。
