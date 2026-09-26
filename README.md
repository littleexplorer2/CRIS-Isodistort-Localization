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
| **ISODISTORT/** | 本地「子群搜索 + 导出结构」 | 上传母相 CIF → 勾选畸变类型 → Method 1–4 搜索/分解 → Distortion 导出结果表与子群结构 ZIP。Method 2 可生成缺失子群库（WSL `~/.id/tmp/i*.iso`），网页与终端均可 Manage 缓存；参数 k 点的位移模式走本地 smodes/(3+d)（nmod）。见 [ISODISTORT/README.md](ISODISTORT/README.md) |
| **ISODISTORT_VALIDATE/** | 核对本地 CIF 是否算对 | 把本地 CIF 放入 `compare/item/`，把官网参考 CIF 放入 `compare/true/`（批量比较须改名一一对应），用 `main.py` 比较并输出 PASS/FAIL。见 [ISODISTORT_VALIDATE/README.md](ISODISTORT_VALIDATE/README.md) |
| **ISOVIZ_INPUT/** | 把振幅 CSV 写入 `.isoviz` 并启动 IsoVIZ | 从桌面 `Best_Model_Parameters/` 读取 GD CSV，再读入绝对路径的子群 `.isoviz`，写入 `amp` 后打开 Java 版 IsoVIZ。见 [ISOVIZ_INPUT/README.md](ISOVIZ_INPUT/README.md) |

仓库级长期规则见 [`AGENTS.md`](AGENTS.md)；各子项目 README **只讲该项目本身**，子项目特有规则与配置分别放在其 `agent.md` 和 `config/settings.yaml`。跨项目怎么串起来，只在本文件说明。

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
   GD 会把振幅 CSV 写到桌面 `Best_Model_Parameters/<irrep>_<structure_type>/`。运行 `ISOVIZ_INPUT/main.py`，输入该子文件夹名、CSV 文件名，以及晶体 `.isoviz` 的绝对路径（可带引号），程序会写入 `amp` 并启动 IsoVIZ。

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
| `GD/` | 仓库内占位（gitignore）。实际拟合代码在桌面 `GD（未同步git）`；目标 2 转换器写在该桌面目录，**不要改本仓库 `GD/`** |
| `webpage_info/` | 官网各步 HTML 存档（按母相分子目录：`EuAl4 Parent.cif/`、`NdNiO2 own.cif/` 等；对照交互顺序用，**勿改**） |
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

## 用户需要配置 / 放置的路径

换机器时按**各子项目自己的 README** 核对，不要把细节写回本文件：

| 子项目 | 配置文件 | 说明 |
| --- | --- | --- |
| ISODISTORT | [ISODISTORT/config/settings.yaml](ISODISTORT/config/settings.yaml) | `iso` 路径、端口、临时/导出目录。见 [ISODISTORT/README.md](ISODISTORT/README.md) |
| ISODISTORT_VALIDATE | [ISODISTORT_VALIDATE/config/settings.yaml](ISODISTORT_VALIDATE/config/settings.yaml) | `compare/item` 与 `compare/true`、默认容差。见 [ISODISTORT_VALIDATE/README.md](ISODISTORT_VALIDATE/README.md) |
| ISOVIZ_INPUT | [ISOVIZ_INPUT/config/settings.yaml](ISOVIZ_INPUT/config/settings.yaml) | 桌面 CSV 目录、IsoVIZ 启动器查找。见 [ISOVIZ_INPUT/README.md](ISOVIZ_INPUT/README.md) |

仓库根目录可以放本机快捷方式（已 gitignore）：`VESTA.lnk`、`ISOViz.lnk`。计算流程不依赖它们是否存在；IsoVIZ 查找顺序以 ISOVIZ_INPUT 的 yaml 为准。

跨项目流水线仍见上文「典型工作流」。

---

## 仓库地图（还可能看到什么）

```text
CRIS/
├── AGENTS.md                      ← 仓库级边界、工作流与文档职责
├── README.md                      ← 本文件：总览与跨项目关系
├── .venv/                         ← 三个子项目共用的虚拟环境（不入库）
├── ISODISTORT/                    ← 本地 ISODISTORT
│   ├── README.md / agent.md
│   └── config/settings.yaml
├── ISODISTORT_VALIDATE/           ← CIF 语义比较
│   ├── README.md / agent.md
│   └── config/settings.yaml
├── ISOVIZ_INPUT/                  ← CSV → 启动 IsoVIZ
│   ├── README.md / agent.md
│   └── config/settings.yaml
├── experiment_data/               ← 【勿改】实验数据
├── webpage_info/                  ← 【勿改】官网 HTML 存档
├── GD/                            ← 【勿改】仓库占位（gitignore）
├── ISOViz.lnk / VESTA.lnk         ← （本机）快捷方式，已 gitignore
```

根目录只保留总览、仓库级精简规则、gitignore 与共享 `.venv`。各子项目的使用说明、特有 Agent 规则和运行配置都放在各自目录里；开发进度和验证报告不得写入 `AGENTS.md` / `agent.md`。

---

## 从哪里开始读

| 你想做的事 | 打开 |
| --- | --- |
| 三个子项目如何串起来 | 本文件 |
| 改 ISODISTORT / 开发目标 | 规则见 [ISODISTORT/agent.md](ISODISTORT/agent.md)，计划见 [ISODISTORT/docs/DEVELOPMENT_PLAN.md](ISODISTORT/docs/DEVELOPMENT_PLAN.md)，使用说明见 [ISODISTORT/README.md](ISODISTORT/README.md) |
| 比较本地 CIF 与官网 CIF | [ISODISTORT_VALIDATE/agent.md](ISODISTORT_VALIDATE/agent.md)、[ISODISTORT_VALIDATE/README.md](ISODISTORT_VALIDATE/README.md) |
| 把振幅写入 IsoVIZ 并打开 | [ISOVIZ_INPUT/agent.md](ISOVIZ_INPUT/agent.md)、[ISOVIZ_INPUT/README.md](ISOVIZ_INPUT/README.md) |

官网帮助（概念背景，非本仓库文档）：[ISODISTORT Help](https://iso.byu.edu/isodistorthelp.php)、[ISOTROPY Suite](https://iso.byu.edu/isotropy.php)。
