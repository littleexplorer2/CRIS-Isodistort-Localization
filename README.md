# CRIS

CRIS 将晶体畸变分析与模式可视化整理为两个本地工具，共用根目录的 Python 环境：

- **[ISODISTORT](ISODISTORT/README.md)**：载入母相 CIF，搜索畸变子群，计算模式并导出 CIF、IsoVIZ、模式详情和 TOPAS 文件；也可用母相与女儿相结构进行 Method 4 模式分解。
- **[ISOVIZ_INPUT](ISOVIZ_INPUT/README.md)**：读取 GD 拟合产生的振幅 CSV 和对应 `.isoviz`，把振幅写入临时副本并交给 IsoVIZ 查看。

ISODISTORT 的独立比较工具位于 `ISODISTORT/tests/validation/`，用于比较已经冻结的本地结果与官网参考。它是测试支持包，不是第三个运行项目，也不参与生产计算。仓库版本以 [`VERSION`](VERSION) 为准。

本地 ISODISTORT 的功能范围以现有网页为准。官网内置在线结构查看、Distortion Generate 和 Domains 不在当前网页范围内。有限测试通过只证明已覆盖的输入域；当前能力、限制及验证证据分别见子项目 README 和[验证报告](ISODISTORT/docs/BUGFIX_VALIDATION_REPORT.md)。

## 已有环境：直接开始

本机工作区是 `E:\CRIS`。已经安装依赖时，在 PowerShell 执行：

```powershell
cd E:\CRIS
.\run_cris.ps1 setup_cris.py doctor --project isodistort
.\run_cris.ps1
```

最后一条命令默认运行 `ISODISTORT\main_web.py`。终端显示实际访问地址并尝试打开浏览器；在网页中按 **Parent CIF → Load → Types → Change → Method → OK → Distortion** 操作。详细字段及下载位置见 [ISODISTORT 使用说明](ISODISTORT/README.md)。

需要把已有振幅送到 IsoVIZ 时：

```powershell
cd E:\CRIS
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py
```

程序默认在 `CRIS/launchers/` 查找启动器，快捷方式的目标程序须可用；也可用 `ISOVIZ` / `ISOVIZ_JAR` 指定其他位置。首次安装请从下一节开始。

## 第一次部署（Windows）

### 1. 准备源码和系统环境

新机器可克隆仓库，也可解压源码 ZIP：

```powershell
git clone https://github.com/littleexplorer2/CRIS-Isodistort-Localization.git CRIS
cd CRIS
```

后文命令均在包含 `setup_cris.py`、`run_cris.ps1` 和两个子项目的 **CRIS 根目录**运行。已有 `E:\CRIS` 工作区时使用原目录，无需再克隆一份。

需要的组件：

- Python **3.10 或更高**。命令示例使用 `py -3.10`；如果已安装的是 3.12，请统一替换为 `py -3.12`，或确认版本后使用 `python`。
- Windows 上的 ISODISTORT 需要 **WSL 和已完成首次初始化的默认 Linux 发行版**，用于运行 Linux 版 ISOTROPY。
- ISOVIZ_INPUT 的可视化需要 **Java 和 IsoVIZ**。只使用 ISODISTORT 计算时无需安装它们。
- VESTA 是查看导出 CIF 的可选外部程序。

先检查：

```powershell
py -3.10 --version
wsl --status
wsl --list --verbose
wsl -e sh -c "echo CRIS_WSL_OK"
```

WSL 未安装或没有默认发行版时，先完成 Windows 的 WSL 安装与发行版初始化。Python 包安装不能代替这一步。

### 2. 核对外部程序

ISODISTORT 使用随项目冻结的 Linux ISOTROPY 文件：

```text
ISODISTORT/resources/isobyu/
  iso
  findsym
  smodes
  data_*.txt
  comsubs           可选
```

该目录是只读计算资源，不存放具体晶体的官网答案。安装器不会下载、更新或改写其中的文件；如果分发副本缺失，先核对分发内容和 [ISOTROPY Suite 官网](https://iso.byu.edu/isotropy.php)。Windows 下也需要 Linux 版程序，不能用 `iso.exe` 替换。

本机快捷方式放在 `CRIS/launchers/`，目前有 `ISOViz.lnk` 和 `VESTA.lnk`。ISOVIZ_INPUT 默认搜索此目录；使用其他位置或名称时可显式覆盖：

```powershell
$env:ISOVIZ = (Resolve-Path .\launchers\ISOViz.lnk).Path
java -version
```

环境变量设置只对当前 PowerShell 及其子进程生效。快捷方式不包含目标软件本身；换机器后需检查目标。`VESTA.lnk` 供手工打开 CIF 使用。根 `.gitignore` 排除整个 `launchers/`，本机程序与快捷方式不提交。

### 3. 安装 Python 依赖

安装两个运行项目：

```powershell
py -3.10 setup_cris.py install
```

只需要一个项目时：

```powershell
py -3.10 setup_cris.py install --project isodistort
# 或
py -3.10 setup_cris.py install --project isoviz
```

安装器在 `CRIS/.venv` 创建或复用**同一个实体环境**，安装所选项目声明的依赖，创建必要运行目录，然后执行部署诊断。`--dev` 额外安装测试和检查依赖。已有环境通常直接复用；`--recreate` 会删除并重建 `.venv`，不是常规安装或更新步骤。

**安装入口以根目录 `setup_cris.py` 为准。** `ISODISTORT/main_requirement.py` 与 `main_web.py` 同在子项目根，作为兼容入口转发到同一个部署器，自动选择 `install --project isodistort`；通过绝对路径也可从其他工作目录调用。

### 4. 检查部署结果

```powershell
.\run_cris.ps1 setup_cris.py doctor --project all
```

只用 ISODISTORT 时将 `all` 改为 `isodistort`；只用 ISOVIZ_INPUT 时也可运行：

```powershell
.\.venv\Scripts\python.exe setup_cris.py doctor --project isoviz
```

诊断的 `PASS` 表示该项部署检查通过，`WARN` 需按提示判断，`FAIL` 表示运行条件未满足。有 `FAIL` 时退出码为 1；仅 `PASS/WARN` 时为 0。尚无 `Best_Model_Parameters` CSV 属于缺少用户输入，会显示 `WARN`。可加 `--json` 输出机器报告。部署检查通过不等于科学结果已通过验证。

`run_cris.ps1` 从 `.venv/pyvenv.cfg` 找到基础 Python，再加载实体 `.venv` 内的依赖，以避开本机曾出现的 OneDrive/WSL 启动限制。它不创建另一套环境，也不要求先激活 `.venv`。需要 WSL 的网页、脚本和测试统一使用这个入口。

若 PowerShell 的执行策略阻止启动本地脚本，可运行：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1
```

### 离线安装

在与目标机器相同的操作系统、CPU 架构和 Python 版本上先下载依赖：

```powershell
py -3.10 setup_cris.py download --wheelhouse C:\cris-wheels
```

把完整目录复制到离线机器，再在 CRIS 根目录执行：

```powershell
py -3.10 setup_cris.py install --wheelhouse C:\cris-wheels
```

开发环境在两条命令中都加 `--dev`；单项目安装在两边使用相同的 `--project`。默认仅下载兼容 wheel；`--allow-source` 留给能自行处理源码构建的部署。这个流程只解决 Python 依赖，WSL、Java、IsoVIZ 仍需单独准备。

## 输入、结果和目录放在哪里

```text
CRIS/
  setup_cris.py / run_cris.ps1 / VERSION   仓库部署、启动和版本入口
  .venv/                                两个运行项目共用的 Python 环境
  .gitkeep                              根目录中的占位文件，不是运行入口
  launchers/                            本机 IsoVIZ / VESTA 快捷方式
  ISODISTORT/                           科学计算与本地网页
    main_web.py                         网页入口
    main_requirement.py                 转发到统一安装器
    backend/ application/ frontend/     三层产品代码
    resources/                          配置、只读通用程序与数据库
    tests/                              自动测试、手工审计与 validation 支持包
    docs/                               使用之外的计划、证据、理论与清单
    output/                             本地运行和验证产物
  ISOVIZ_INPUT/                         振幅写入及 IsoVIZ 启动工具
  tests_repo/                           仓库部署、架构与打包检查
  experiment_data/                      只读实验输入
  webpage_info/                         官网页面存档
  output_compare/                       官网输出归档
```

源码 ZIP/Git 不保证包含本机只读数据、实验文件、快捷方式或拟合 CSV；这些材料需保留或按其来源取得，安装器不会生成它们。

- **网页下载**：保存到浏览器的下载位置或你在“另存为”中选择的位置。网页 ZIP 当次计算生成，不会自动落到 `ISODISTORT/output/`，也不会扫描旧输出目录。
- **Python API 输出**：按调用中指定的目录保存；配置中的运行输出默认位于 `ISODISTORT/output/`。
- **机器验证报告**：放在 `ISODISTORT/output/validation/`；正式操作见[手工验证说明](ISODISTORT/docs/MANUAL_VALIDATION.md)。
- **GD 振幅 CSV**：默认从 CRIS 同级的 `Best_Model_Parameters/` 查找；ISOVIZ_INPUT 写入系统临时副本，保留原 CSV 和 `.isoviz`。

建议保留本机同级布局：

```text
E:\CRIS
E:\GD
E:\Best_Model_Parameters
```

GD 使用其自身环境与规则，不与 CRIS 的环境混装。不同布局可按子项目说明设置 `CRIS_ROOT`、`BEST_MODEL_PARAMETERS_DIR`。迁移后先核对这些输入位置、快捷方式目标和 `.venv/pyvenv.cfg` 中的基础解释器；虚拟环境不能假定可直接跨路径搬移。

## 配置与科学数据边界

配置各有唯一入口：

- [ISODISTORT 配置](ISODISTORT/resources/config/settings.yaml)：ISOTROPY 路径、科学容差、端口、超时与本地输出目录。
- [ISOVIZ_INPUT 配置](ISOVIZ_INPUT/config/settings.yaml)：CSV 路径、启动器环境变量与搜索目录。
- [比较器配置](ISODISTORT/tests/validation/config/settings.yaml)：已有比较文件目录、容差和严格模式。

配置内的相对路径按相应配置目录解析，具体覆盖方式见子项目 README。`resources/isobyu/`、实验原始输入、既有官网页面与输出保持只读。官网下载与归档按[下载清单](ISODISTORT/docs/DOWNLOAD_CHECKLIST.md)及其登记规则执行，不在本 README 重复维护下载状态。

**官网结果只用于独立审计，不能进入本地计算。** 搜索、模式构造、归一化、拟合、缓存和导出都不得读取或回填官网答案。科学比较先固定独立输入，再运行并冻结本地结果，最后才读取官网参考；比较器不会生成或修补结构。

## 维护与故障排查

更新代码前先检查并保留本地修改。依赖变化后重新执行 `setup_cris.py install`，通常无需重建环境，再运行对应 `doctor`。

- **终端打印地址但没有浏览器**：手动打开终端显示的最终 URL。打印启动请求只证明已分发，不证明浏览器页面已正确显示。
- **`Wsl/E_ACCESSDENIED`**：从 CRIS 根目录使用 `.\run_cris.ps1 <Python 参数>`；核对基础 Python 仍存在。不要另建一套项目环境绕过。
- **缺少二进制、数据库或 Python 包**：按 `doctor` 的具体失败项核对分发资源、配置和依赖。
- **找不到 IsoVIZ**：确认 `launchers/` 内启动器名称、快捷方式目标及所需 Java；使用其他目录时设置 `ISOVIZ` / `ISOVIZ_JAR` 绝对路径。
- **计算超时后不能继续**：网页重新上传母相，Python 重新创建会话；不要继续复用已经失效的运行实例。
- **没有参考材料或比较项**：缺失不算通过，也不能用本地产物填官网目录。按下载清单处理。

开发检查从 CRIS 根目录运行：

```powershell
.\run_cris.ps1 -m pytest tests_repo -q --tb=line
.\run_cris.ps1 -m pytest ISODISTORT\tests -q --tb=line
```

这两条分别检查仓库集成与 ISODISTORT；需要先安装 `--dev` 依赖。按模块单独检查的路径和命令见 [ISODISTORT 功能检查](ISODISTORT/tests/README.md)。真实 WSL 计算、官网差分及 GUI 验收有各自条件，不能用静态检查或进程启动代替。

## 继续阅读

- [ISODISTORT README](ISODISTORT/README.md)：网页字段、四种 Method、导出、Python API 和能力限制。
- [ISOVIZ_INPUT README](ISOVIZ_INPUT/README.md)：CSV/文本输入、模式匹配、临时副本和 GUI 验收。
- [DOWNLOAD_CHECKLIST](ISODISTORT/docs/DOWNLOAD_CHECKLIST.md)：唯一人工下载操作清单。
- [MANUAL_VALIDATION](ISODISTORT/docs/MANUAL_VALIDATION.md)：长时及手工验证命令。
- [BUGFIX_VALIDATION_REPORT](ISODISTORT/docs/BUGFIX_VALIDATION_REPORT.md)：当前验证证据和缺口。
- [DEVELOPMENT_PLAN](ISODISTORT/docs/DEVELOPMENT_PLAN.md)：待修问题及后续计划。
- [理论参考资料库](ISODISTORT/docs/reference_library/README.md)：科学来源、完整推导、适用域和能力映射。
- [AGENTS.md](AGENTS.md) 与各子项目 `agent.md`：开发边界和模块职责。
