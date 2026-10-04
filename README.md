# CRIS 0.4.0

CRIS 是一套本地晶体学工作流。它把三个原本需要分别操作的环节放在同一个仓库中：

1. 从母相 CIF 搜索可能的畸变子群并导出结构；
2. 把本地导出的 CIF 与官网参考 CIF 做晶体学语义比较；
3. 把拟合得到的模式振幅写入 `.isoviz`，再交给 IsoVIZ 查看。

其中 `ISODISTORT/` 是对官网现有本地网页版功能范围的本地化逆向工程。当前工作重点是
修复缺陷、提高数值与晶体学结果准确性，并把算法推广到任意适用晶体；官网存在但本地
网页没有的在线结构查看、Distortion Generate 和 Domains 暂不实现。计算逻辑以晶体学
定义、可检验不变量及 Stokes、Campbell、Hatch 等人的 ISOTROPY/ISODISTORT 论文为
依据，官网下载结果只作为独立差分证据。

第一次接触本项目时，先按下文“第一次安装”完成环境检查，再只阅读你要使用的子项目 README。程序界面为英语，安装和使用说明为中文。

> 验证结论只适用于已覆盖的输入和功能；有限母相或候选通过不能外推为任意晶体
> 或跨晶系科研级精度。最新实测证据见
> [BUGFIX_VALIDATION_REPORT.md](ISODISTORT/docs/BUGFIX_VALIDATION_REPORT.md)，
> 当前工作顺序见 [DEVELOPMENT_PLAN.md](ISODISTORT/docs/DEVELOPMENT_PLAN.md)。

## 三个子项目

| 子项目 | 用途 | 输入 | 输出 / 动作 |
| --- | --- | --- | --- |
| [ISODISTORT](ISODISTORT/README.md) | 本地子群搜索、Method 4 模式分解 | 母相 CIF；Method 4 还需女儿相 CIF | 结果表、CIF、`.isoviz`、modes、TOPAS ZIP |
| [ISODISTORT_VALIDATE](ISODISTORT_VALIDATE/README.md) | 比较本地与官网 CIF | `compare/item/` 与 `compare/true/` 中同相对路径的 CIF | PASS/FAIL、差异说明、可选 JSON |
| [ISOVIZ_INPUT](ISOVIZ_INPUT/README.md) | 把 GD 振幅写入 `.isoviz` | `Best_Model_Parameters` CSV + 对应 `.isoviz` | 系统临时副本，并可启动 IsoVIZ |

三个子项目共用仓库根目录内的实体 `.venv`。仓库级版本只由 [`VERSION`](VERSION) 定义；三个 Python 包与文档在发布时保持同一版本。本机的 OneDrive 安全限制会阻止该目录内的 Python 子进程启动 WSL，因此需要 WSL 的命令统一由根目录 `run_cris.ps1` 启动；第三方 Python 包仍全部从 `.venv` 加载。

## 第一次安装（Windows）

### 1. 获取代码

安装 Git 后，在 PowerShell 中执行：

```powershell
git clone https://github.com/littleexplorer2/CRIS-Isodistort-Localization.git CRIS
cd CRIS
```

也可以从 GitHub 下载 ZIP 并解压；以下命令均假定当前目录是包含
`setup_cris.py`、`ISODISTORT/`、`ISODISTORT_VALIDATE/` 和 `ISOVIZ_INPUT/`
的 CRIS 根目录。

### 2. 安装系统组件

完整使用三个子项目需要：

| 组件 | 最低要求 | 哪个项目需要 |
| --- | --- | --- |
| Python | 3.10 或更高；当前开发环境为 3.12 | 全部 |
| WSL + 默认 Linux 发行版 | `wsl -e sh -c "echo ok"` 能执行 | ISODISTORT（Windows） |
| Linux 版 ISOTROPY Suite | `iso`、`findsym`、`smodes`、完整 `data_*.txt` | ISODISTORT |
| Java + IsoVIZ | `java -version` 可执行，且能定位 `.lnk` / `.jar` / `.exe` | ISOVIZ_INPUT |

只使用 `ISODISTORT_VALIDATE` 时不需要 WSL、ISOTROPY、Java 或 IsoVIZ。

先在 PowerShell 检查 Python 与 WSL：

```powershell
py -3.10 --version
wsl --status
wsl --list --verbose
wsl -e sh -c "echo CRIS_WSL_OK"
```

若 `py -3.10` 不存在，可把后续命令中的它替换为本机可用的 Python 3.10+ 启动命令，例如 `python`。

### 3. 放置不能由 pip 安装的程序

ISOTROPY Suite 与 IsoVIZ 需要从 [ISOTROPY Suite 官网](https://iso.byu.edu/isotropy.php)人工获取，本仓库不会伪造、改写或自动下载这些第三方文件。

把 Linux 版 ISOTROPY 文件放到只读目录：

```text
ISODISTORT/isobyu/
  iso
  smodes
  findsym          （Wyckoff 标准代表与位点参数）
  comsubs          （可选功能会使用）
  data_*.txt       （完整数据库）
```

Windows 下不要放 Windows 版 `iso.exe`；本项目通过 WSL 运行 Linux ELF。

IsoVIZ 任选一种配置方式：

- 把 `ISOViz.lnk`、`ISOViz.jar` 或受支持的 `.exe` 放到仓库根目录；
- 设置环境变量 `ISOVIZ` 或 `ISOVIZ_JAR` 为完整路径；
- 详细查找顺序见 [ISOVIZ_INPUT/config/settings.yaml](ISOVIZ_INPUT/config/settings.yaml)。

### 4. 创建共享 Python 环境

在 CRIS 根目录执行唯一安装入口：

```powershell
py -3.10 setup_cris.py install
```

它会：

- 在 `CRIS/.venv` 创建或复用实体虚拟环境；
- 按三个项目的 `requirements.txt` 校验版本并安装不满足的依赖；
- 创建本地运行目录，但不改原始数据或第三方二进制；
- 检查包导入、WSL、ISOTROPY 数据、Java、IsoVIZ 启动器和输入目录；
- 以 `PASS`、`WARN`、`FAIL` 明确区分成功、缺少用户输入和部署错误。

常用安装选项：

```powershell
# 包含 pytest / ruff，供开发和完整测试使用
py -3.10 setup_cris.py install --dev

# 只安装一个子项目及其运行目录
py -3.10 setup_cris.py install --project isodistort
py -3.10 setup_cris.py install --project validate
py -3.10 setup_cris.py install --project isoviz

# 仅在 .venv 已损坏或 Python 大版本改变时使用
py -3.10 setup_cris.py install --recreate
```

旧的三个 `main_requirement.py` 仍可调用，但只作为兼容转发；新说明和新部署统一使用根目录 `setup_cris.py`。

### 5. 运行部署诊断

安装后或换机器、改配置后运行：

```powershell
.\run_cris.ps1 setup_cris.py doctor --project all
```

只诊断一个子项目：

```powershell
.\run_cris.ps1 setup_cris.py doctor --project isodistort
.\.venv\Scripts\python.exe setup_cris.py doctor --project validate
.\.venv\Scripts\python.exe setup_cris.py doctor --project isoviz
```

`run_cris.ps1` 从 `.venv/pyvenv.cfg` 读取基础 Python 路径，并把
`.venv/Lib/site-packages` 置于导入路径最前。它不会另建环境或从系统环境取代项目依赖。

自动化环境可加 `--json`。有 `FAIL` 时退出码为 1；只有 `PASS/WARN` 时为 0。`compare/` 尚无成对 CIF、`Best_Model_Parameters/` 尚无 CSV 会记为 `WARN`，表示程序已部署但还没有用户输入，不等同于计算失败。

## 离线安装 Python 依赖

在与目标机器相同的操作系统、CPU 架构和 Python 版本上联网下载：

```powershell
py -3.10 setup_cris.py download --wheelhouse C:\cris-wheels
```

把整个 `C:\cris-wheels` 复制到离线机器，再在 CRIS 根目录执行：

```powershell
py -3.10 setup_cris.py install --wheelhouse C:\cris-wheels
```

开发依赖两边都加 `--dev`。下载命令默认只接受与当前平台兼容的 wheel，从而避免离线机器临时编译失败；只有明确能够处理源码构建时才使用 `--allow-source`。该流程只下载 Python 包，ISOTROPY、IsoVIZ、Java、WSL 仍按上一节人工准备。

## 5 分钟快速上手

### 启动本地 ISODISTORT 网页

```powershell
cd <CRIS 根目录>
.\run_cris.ps1
```

不传参数时启动器默认运行 `ISODISTORT\main_web.py`。终端会打印实际 URL，通常为
`http://127.0.0.1:8000/`，并通过 Windows Shell 打开默认浏览器。浏览器中依次完成：

1. 在 `Parent CIF` 选择母相文件并点 `Load`；
2. 在 `Types of distortions` 选择畸变类型并点 `Change`；
3. 在 Method 1、2 或 3 填参数并点 `OK`；Method 4 另需女儿相 CIF；
4. 在 `Distortion` 过滤结果并下载表格或结构 ZIP。

完整字段、四种 Method、配置和限制见 [ISODISTORT 使用说明](ISODISTORT/README.md)。

### 比较一对 CIF

把两边文件放成相同相对路径：

```text
ISODISTORT_VALIDATE/compare/item/LD1_C1/subgroup.cif
ISODISTORT_VALIDATE/compare/true/LD1_C1/subgroup.cif
```

运行：

```powershell
.\.venv\Scripts\python.exe ISODISTORT_VALIDATE\main.py compare "LD1_C1/subgroup.cif"
```

默认比较晶格、周期分数坐标、物种、占据率、磁矩和空间群语义；字节排版不同不会单独导致失败。详见 [ISODISTORT_VALIDATE 使用说明](ISODISTORT_VALIDATE/README.md)。

### 把振幅写入 IsoVIZ

先确认 GD 已生成：

```text
<桌面>/Best_Model_Parameters/<IR_OPD>/<IR_OPD>_best_model_parameters.csv
```

再运行：

```powershell
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py
```

按提示选择 CSV 并粘贴对应 `.isoviz` 的绝对路径。程序读取 `Best Model Parameter`，只修改系统临时副本中的 `amp`，不会覆盖 CSV 或原始 `.isoviz`。详见 [ISOVIZ_INPUT 使用说明](ISOVIZ_INPUT/README.md)。

## 典型跨项目流程

```text
母相 CIF
   │
   ▼
ISODISTORT ──► 结果表 / CIF / .isoviz / modes / TOPAS
   │                         │
   │                         ├──► ISODISTORT_VALIDATE（与官网 CIF 做语义比较）
   │                         │
   │                         └──► ISOVIZ_INPUT（GD CSV 振幅 → IsoVIZ 临时副本）
   │
   └── Method 4：母相 + 女儿相 CIF → 模式幅度与 residual
```

这些步骤可以独立使用。验证工具不生成结构，IsoVIZ 输入工具不搜索子群，三个项目不会复制彼此的核心算法。

## 配置与本机路径

| 子项目 | 配置事实来源 | 主要内容 |
| --- | --- | --- |
| ISODISTORT | [config/settings.yaml](ISODISTORT/config/settings.yaml) | ISOTROPY 路径、容差、端口、超时、输出目录 |
| ISODISTORT_VALIDATE | [config/settings.yaml](ISODISTORT_VALIDATE/config/settings.yaml) | 固定比较目录、默认容差、严格模式 |
| ISOVIZ_INPUT | [config/settings.yaml](ISOVIZ_INPUT/config/settings.yaml) | CSV 根目录、启动器环境变量和查找顺序 |

相对路径均按各配置文件所在的 `config/` 目录解析。不要把本机绝对路径写入源码；需要覆盖时优先改 yaml 或使用说明中列出的环境变量/命令参数。

## 数据安全与目录边界

以下内容始终只读：

- `experiment_data/`：实验母相和原始输入；
- `webpage_info/`：官网交互与页面存档；
- `output_compare/`：官网/本地对照证据；
- `ISODISTORT/isobyu/`：第三方二进制和数据库；
- 仓库内 `GD/`；桌面 `GD（未同步git）` 另有自己的规则。

运行输出只写入各项目约定的 `output/`、`compare/`、系统临时目录或用户明确指定的位置。安装脚本不会修改上述只读数据。

## 更新已有安装

拉取新版本后通常不需要重建 `.venv`：

```powershell
git pull
py -3.10 setup_cris.py install
.\run_cris.ps1 setup_cris.py doctor --project all
```

`pip install -r` 会重新解析版本约束；已经满足的包不会重复安装。只有虚拟环境损坏或切换 Python 主/次版本时才加 `--recreate`。

## 常见问题

| 现象 | 处理 |
| --- | --- |
| `py` 或 `python` 找不到 | 安装 Python 3.10+，安装时勾选加入 PATH，重新打开 PowerShell |
| PowerShell 不允许 `Activate.ps1` | 不必激活；ISODISTORT 使用 `.\run_cris.ps1`，另两个子项目可直接使用 `.\.venv\Scripts\python.exe` |
| WSL 检查失败 | 先手工运行 `wsl -e sh -c "echo ok"`；确认存在默认发行版且已完成首次初始化 |
| 直接运行 `main_web.py` / `main_terminal.py` 提示不能调用 WSL | 从 CRIS 根目录运行 `.\run_cris.ps1`（网页）或 `.\run_cris.ps1 ISODISTORT\main_terminal.py`；本机需要由 OneDrive 外的基础 Python 建立首个进程，包仍来自实体 `.venv` |
| 服务打印 URL 但浏览器没有出现 | 先确认使用 `.\run_cris.ps1`；成功分发时终端会打印 `Browser launch requested`，失败时会明确提示手动访问最终 URL |
| 找不到 `iso` / `findsym` / `smodes` / `data_*.txt` | 重新核对 `ISODISTORT/isobyu/` 与 `ISODISTORT/config/settings.yaml` |
| Python 包缺失或版本不满足 | 重新运行 `py -3.10 setup_cris.py install`；不要手工逐个猜包名 |
| 离线安装提示缺 wheel | 用相同平台和 Python 版本重新执行 `download`，并完整复制 wheelhouse |
| VALIDATE 没有可比较文件 | 两侧 CIF 必须使用完全相同的相对路径和文件名 |
| 找不到 Java / IsoVIZ | 运行 `doctor --project isoviz`，按其 `Fix` 提示配置 Java 与启动器 |
| IsoVIZ 进程已分发但结果未确认 | 进程启动不等于 GUI 正确识别；按 ISOVIZ_INPUT README 做独立人工验收 |

## 文档入口

- 新用户安装、跨项目关系和部署诊断：本文件；
- ISODISTORT 使用、能力和限制：[ISODISTORT/README.md](ISODISTORT/README.md)；
- CIF 比较使用说明：[ISODISTORT_VALIDATE/README.md](ISODISTORT_VALIDATE/README.md)；
- CSV → IsoVIZ 使用说明：[ISOVIZ_INPUT/README.md](ISOVIZ_INPUT/README.md)；
- 开发计划：[DEVELOPMENT_PLAN.md](ISODISTORT/docs/DEVELOPMENT_PLAN.md)；
- 已完成修复与验证证据：[BUGFIX_VALIDATION_REPORT.md](ISODISTORT/docs/BUGFIX_VALIDATION_REPORT.md)；
- 官网下载操作：[DOWNLOAD_CHECKLIST.md](ISODISTORT/docs/DOWNLOAD_CHECKLIST.md)；
- 手工和长时验证命令：[MANUAL_VALIDATION.md](ISODISTORT/docs/MANUAL_VALIDATION.md)。

仓库规则见 [`AGENTS.md`](AGENTS.md)，子项目规则见各目录的 `agent.md`。这些规则文件只定义长期边界，不代替用户说明。
