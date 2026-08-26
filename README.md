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
| **ISODISTORT/** | 本地「子群搜索 + 导出结构」 | 上传母相 CIF → 勾选畸变类型 → Method 1–4 搜索/分解 → 在 Distortion 下载结果表与子群结构 ZIP。见 [ISODISTORT/README.md](ISODISTORT/README.md) |
| **ISODISTORT_VALIDATE/** | 核对本地 CIF 是否算对 | 把本地导出的 CIF 与官网（或其它可信）参考 CIF 做语义比较，输出 PASS/FAIL。见 [ISODISTORT_VALIDATE/README.md](ISODISTORT_VALIDATE/README.md) |
| **ISOVIZ_INPUT/** | 把振幅 CSV 写入 `.isoviz` 并打开 IsoVIZ | 读取梯度下降得到的 Best Model Parameter，写入子群 `.isoviz` 的 `amp`，再启动 Java 版 IsoVIZ。见 [ISOVIZ_INPUT/README.md](ISOVIZ_INPUT/README.md) |

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
   把本地 ZIP 里某个 `… CIF.cif` 与官网第 6 页导出的同名子群 CIF 配对比较，确认本地化没有算错结构。

4. **（可选）用 ISOVIZ_INPUT 看拟合振幅**  
   若你已有振幅 CSV（例如梯度下降给出的 Best Model Parameter）和对应子群 `.isoviz`，用本工具写入 `amp` 并打开 IsoVIZ 查看畸变结构。

```text
母相 CIF
   │
   ▼
ISODISTORT  ──►  子群表 / CIF / .isoviz / modes / TOPAS
   │                    │
   │                    ├──► ISODISTORT_VALIDATE（与官网 CIF 比对）
   │                    │
   │                    └──► ISOVIZ_INPUT（CSV 振幅 → 写入 .isoviz → 打开 IsoVIZ）
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
python ISOVIZ_INPUT\main_requirement.py
```

之后请始终用 `.venv\Scripts\python.exe` 运行各子项目入口，避免装到系统 Python。各子项目 README 里写有本项目的启动命令。

---

## 仓库地图（还可能看到什么）

```text
CRIS/
├── README.md                 ← 本文件：总览与跨项目关系
├── .venv/                    ← 共享虚拟环境（不入库）
├── ISODISTORT/               ← 本地 ISODISTORT（网页 / 终端 / API）
├── ISODISTORT_VALIDATE/      ← CIF 语义比较
├── ISOVIZ_INPUT/             ← CSV → .isoviz → IsoVIZ
├── experiment_data/          ← 【勿改】实验数据
├── webpage_info/             ← 【勿改】官网 HTML 存档
├── GD/                       ← 【勿改】梯度下降相关
└── ISOViz.lnk                ← （可选）IsoVIZ 快捷方式，已在 .gitignore
```

---

## 从哪里开始读

| 你想做的事 | 打开 |
| --- | --- |
| 安装环境、跑网页/终端搜索子群 | [ISODISTORT/README.md](ISODISTORT/README.md) |
| 比较本地 CIF 与官网 CIF | [ISODISTORT_VALIDATE/README.md](ISODISTORT_VALIDATE/README.md) |
| 把振幅写入 IsoVIZ 并打开 | [ISOVIZ_INPUT/README.md](ISOVIZ_INPUT/README.md) |

官网帮助（概念背景，非本仓库文档）：[ISODISTORT Help](https://iso.byu.edu/isodistorthelp.php)、[ISOTROPY Suite](https://iso.byu.edu/isotropy.php)。
