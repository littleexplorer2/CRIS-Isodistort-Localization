# ISOVIZ_INPUT 0.4.0

把振幅表（CSV）里的数值写进官方 IsoVIZ 的子群结构文件（`.isoviz`），并**自动启动** Java 版 IsoVIZ，查看畸变后的晶体结构。本子项目**不使用** `output/`：读完输入文件和参数后直接打开 IsoVIZ，而不是把结果写成仓库里的产物文件。

你不需要先懂晶体学。可以把它理解成：

1. 你已经有一份「每个原子模式该拧到多少」的表格（CSV，常见列名 **Best Model Parameter**）。  
2. 你已经有一份对应子群的 `.isoviz`（里面有各个模式的进度条字段 `amp`）。  
3. 本工具按模式名字（或顺序别名）把 CSV 里的数写进 `amp`，再启动 IsoVIZ。

IsoVIZ 属于 [ISOTROPY Suite](https://iso.byu.edu/isotropy.php)。手动拖动每个模式进度条既慢又容易出错；本子项目用 Python 自动填写。

本项目与仓库根目录的 `CRIS/.venv` 共用一份 Python 虚拟环境。说明与配置都在本目录：[`README.md`](README.md)、[`agent.md`](agent.md)、[`config/settings.yaml`](config/settings.yaml)。子项目之间的关系见仓库根目录 [README.md](../README.md)。若要先核对 ISODISTORT 导出的 CIF 是否与官网一致，见 [ISODISTORT_VALIDATE/README.md](../ISODISTORT_VALIDATE/README.md)。

---

## 目录结构

```text
ISOVIZ_INPUT/
  agent.md                Agent 工作指南（修改边界、验收）
  README.md               本文件（使用说明）
  config/settings.yaml    桌面 CSV 目录、input_content、IsoVIZ 查找规则
  main.py                 读取输入并启动 IsoVIZ
  main_requirement.py     旧命令兼容转发；实际安装逻辑在根 setup_cris.py
  requirements.txt        运行时依赖（PyYAML）
  requirements-dev.txt    开发依赖（pytest）
  pyproject.toml
  isoviz_input/           包：配置加载、CSV 解析、写 .isoviz、启动 IsoVIZ
  tests_dev/              开发测试（用 tests_dev/fixtures/，不依赖你本机样本）
  input_content/          本地输入（gitignore，不入库）
    data.csv/             振幅 CSV 文件夹（目录名就是 data.csv）
    subgroup.isoviz/      官方子群 .isoviz 文件夹（目录名就是 subgroup.isoviz）
```

日常振幅 CSV 由 GD 笔记本自动写到桌面：

```text
<桌面>/Best_Model_Parameters/
  <irrep>_<structure_type>/
    <irrep>_<structure_type>_best_model_parameters.csv
```

例如 `Best_Model_Parameters/LD1_C1/LD1_C1_best_model_parameters.csv`。列出文件夹和 CSV 时按名字不区分大小写排序；CSV 行按模式序号 `a1, a2, …, a10` 排序。子群 `.isoviz` 在运行时用绝对路径提供。`input_content/` 仍会自动创建（历史兼容，gitignore），但默认不再从那里读振幅。

---

## 用户需要配置 / 放置的路径

| 项目 | 默认 / 做法 | 你要做什么 |
| --- | --- | --- |
| **振幅 CSV** | 桌面 `Best_Model_Parameters/<irrep>_<structure_type>/`（目录名见 yaml） | GD 自动写出；运行 `main.py` 时输入该子文件夹名和 CSV 文件名 |
| **子群 `.isoviz`** | 运行时输入**绝对路径** | 可带英文或中文引号；也可用 `--structure` |
| **IsoVIZ 可执行体** | `config/settings.yaml` → `isoviz` | **本机必配其一**：本目录或仓库根的 `ISOViz.lnk` / `.jar` / `.exe`，或环境变量 `ISOVIZ` / `ISOVIZ_JAR` |
| **Java** | 系统 `PATH` 中的 `java` | 安装 JRE/JDK；用 `.jar` 启动时必需 |
| **临时启动文件** | 系统临时目录（`tempfile`） | 程序自动写入再交给 IsoVIZ；**不要**也不需要配置本子项目的 `output/` |
| **Python / venv** | 仓库根 `CRIS/.venv` | 与其它子项目共用 |

**一般不必改：** 换机器时核对本目录 `config/settings.yaml`。只有桌面不在默认位置、或 IsoVIZ 不在查找目录里时才改 yaml / 环境变量。

Agent 约定见 [agent.md](agent.md)。跨项目总览见仓库根 [README.md](../README.md)。

---

## 安装

### 1. Python 环境

需要 **Python ≥ 3.10**。在 **CRIS 根目录**执行统一安装入口：

```powershell
cd <CRIS 根目录>
py -3.10 setup_cris.py install --project isoviz
```

可选：

| 参数 | 作用 |
| --- | --- |
| `--dev` | 额外安装 `requirements-dev.txt`（pytest） |
| `--recreate` | 强制重建 `CRIS/.venv` |

脚本会：

1. 确认 Python 版本；
2. 创建或复用 `CRIS/.venv`，并按版本约束安装依赖；
3. 补齐仓库内历史兼容的 `input_content/` 目录；
4. 检查 Java、IsoVIZ 启动器以及桌面 `Best_Model_Parameters` 中是否已有 CSV。缺少用户 CSV 记为 `WARN`，不会伪装成 GUI 已验证。

只做只读诊断：

```powershell
.\.venv\Scripts\python.exe setup_cris.py doctor --project isoviz
```

### 2. Java 与 IsoVIZ（必需）

- 安装 JRE/JDK，保证终端里能运行 `java -version`。  
- 自行安装 IsoVIZ（ISOTROPY Suite 的一部分）。  
- 任选一种方式让本工具找到它：

| 方式 | 说明 |
| --- | --- |
| 快捷方式 | 把快捷方式放到 **ISOVIZ_INPUT/** 或仓库**根目录**，命名为 `ISOViz.lnk`（或 `IsoVIZ.lnk` / `ISOVIZ.lnk`）。根目录快捷方式已在根 `.gitignore` 中；查找顺序见 `config/settings.yaml` |
| 可执行文件 / JAR | 上述查找目录中放置 `IsoViz.exe` / `ISOViz.exe`，或 `IsoViz.jar` / `ISOViz.jar` / `isoviz.jar` |
| 环境变量 | 设置 `ISOVIZ` 或 `ISOVIZ_JAR` 指向 `.exe` / `.jar` 的完整路径 |
| Windows 文件关联 | 若 `.isoviz` 已关联到 IsoVIZ，脚本也可直接 `startfile` 打开 |

本工具会把填好振幅的内容写成临时文件再交给 IsoVIZ，**不会**在本子项目里维护 `output/`。

---

## 使用：`main.py`

从 CRIS 根目录（推荐交互方式）：

```powershell
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py
```

程序会：

1. 列出桌面 `Best_Model_Parameters` 下的子文件夹，请**输入文件夹名**（或列表编号）
2. 列出该文件夹中的 CSV，请**输入文件名**（可省略 `.csv`，或输入列表编号）
3. 请输入晶体 `.isoviz` 的**绝对路径**。从资源管理器复制时可能带 `"..."` 或 `“...”`，程序会去掉引号

```text
Folders in C:\Users\...\Desktop\Best_Model_Parameters:
  1. LD1_C1
Best_Model_Parameters folder name: LD1_C1
CSV files in ...\LD1_C1:
  1. LD1_C1_best_model_parameters.csv
CSV file name: LD1_C1_best_model_parameters.csv
Absolute path of the crystal .isoviz file: "D:\data\LD1_C1.isoviz"
```

两个输入齐了之后，程序把 CSV 振幅写入该 `.isoviz` 的 `amp`，再通过 CRIS 根目录的 IsoVIZ 快捷方式/JAR/EXE 打开填好的临时文件。成功时打印匹配到的模式数、CSV 未用名字、以及保持原振幅的 IsoVIZ 模式。

### 参数一览

| 参数 | 是否必填 | 含义 |
| --- | --- | --- |
| （无参数） | 默认 | 询问文件夹名、CSV 文件名、`.isoviz` 绝对路径 |
| `--folder` | 可选 | `Best_Model_Parameters` 下的子文件夹名 |
| `--csv-name` | 可选 | 该文件夹中的 CSV 文件名 |
| `--data` | 可选 | 振幅 CSV 的完整路径（跳过文件夹/文件名询问） |
| `--structure` | 可选 | `.isoviz` 路径；可带引号 |

### CSV 需要什么列

CSV **必须有表头**。识别列名时不区分大小写。常用（与梯度下降脚本写出的表兼容）：

| 列名（示例） | 用途 |
| --- | --- |
| **Mode Name**（或 `modelabel` / `label` / `name`） | 与 `.isoviz` 中的模式标签匹配，例如 `[0,0,1/6]LD1[Eu1:a:dsp]A2u(a)` |
| **Best Model Parameter**（或 `amplitude` / `amp` / `value`） | 写入 IsoVIZ 进度条使用的 **`amp`**（与 `maxamp` 同一单位）。**不要**误用 Normalized Amplitude |
| **Mode** | 备选：`a1`, `a2`, … 按 `.isoviz` 文件中 **strain 模式再 displacive 模式**的出现顺序对应 |
| Maximum Mode Amplitude / maxamp | 可选；解析会读入，匹配主要仍靠名字或 Mode 别名 |

模式名里可能含逗号；推荐用带引号的 CSV（例如 pandas 默认写出格式）。未匹配到的 IsoVIZ 模式保持原振幅。

### 完整示例

```powershell
cd <CRIS 根目录>

# 1) 准备环境（若尚未做过）
py -3.10 setup_cris.py install --project isoviz

# 2) 运行 GD 笔记本写出 Desktop\Best_Model_Parameters\LD1_C1\...csv
# 3) 交互输入文件夹名、CSV 文件名、.isoviz 绝对路径
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py

# 或非交互：
# .\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py --folder LD1_C1 --csv-name LD1_C1_best_model_parameters.csv --structure "D:\data\LD1_C1.isoviz"
```

### 建议你第一次这样试用

1. 先运行 GD 笔记本的保存单元，确认
   `Desktop\Best_Model_Parameters\LD1_C1\LD1_C1_best_model_parameters.csv`
   实际存在。只有 `LD1_C1` 空目录时，主程序无法继续。
2. 使用与这份 CSV **同一个 IR/OPD/path** 导出的官方 `.isoviz`；不能把
   `LD1_C1` 振幅写进另一个子群文件后再据此判断结构是否正确。
3. 从 CRIS 根目录运行 `.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py`，依次输入
   `LD1_C1`、CSV 文件名和 `.isoviz` 的绝对路径。
4. 先看终端报告：`[matched]` 应大于 0；`[csv leftover]` 表示 CSV 中有模式未写入；
   `isoviz modes without CSV values` 表示这些滑条保留原值。这两类不一定都是错误，
   但在判定结构前必须逐项解释。
5. IsoVIZ 打开后，核对母相/子群、模式总数和标签，并抽查终端列出的 `amp`；
   然后观察结构是否随这些振幅变化。`Best Model Parameter` 是写入值，
   `Normalized Amplitude` 不是。

如果默认查找位置没有 IsoVIZ 启动器，仅有 Java 还不够；请把 `ISOViz.lnk`
放到 CRIS 根目录或 `ISOVIZ_INPUT/`，或者设置 `ISOVIZ` / `ISOVIZ_JAR`。源
`.isoviz` 和 CSV 都不会被修改；被打开的是系统临时目录中的副本。

---

## 测试

```powershell
cd <CRIS 根目录>
.\.venv\Scripts\python.exe -m pytest ISOVIZ_INPUT\tests_dev -q --tb=line --basetemp .test-tmp-isoviz
```

开发依赖：

```powershell
py -3.10 setup_cris.py install --project isoviz --dev
```

测试使用 `tests_dev/fixtures/` 内的样本，不依赖你本机的 `input_content/`。

### GUI 打开与数据识别验收

普通 pytest 现在自动覆盖三层：CSV/标签匹配与 `amp` 写入、补丁后模式身份和结构
元数据保持、JAR/EXE/Windows 快捷方式是否把临时 `.isoviz` 作为明确参数传给
IsoVIZ。这些测试不弹出窗口。

此前的测试只把启动函数替换成 mock，并检查它收到的文本，**不能证明窗口真的打开，
也不能证明 Java IsoVIZ 已识别结构和模式**。现在新增了显式手工验收器：

```powershell
# 只生成临时文件、静态报告和 GUI 中应看到的值；不启动 GUI
.\.venv\Scripts\python.exe ISOVIZ_INPUT\tests_dev\manual\validate_isoviz_gui.py `
  --data "C:\path\to\Best_Model_Parameters\LD1_C1\LD1_C1_best_model_parameters.csv" `
  --structure "D:\data\LD1_C1.isoviz"

# 之后需要补跑的真实 GUI 验收；会启动并询问两个独立问题
.\.venv\Scripts\python.exe ISOVIZ_INPUT\tests_dev\manual\validate_isoviz_gui.py `
  --data "C:\path\to\Best_Model_Parameters\LD1_C1\LD1_C1_best_model_parameters.csv" `
  --structure "D:\data\LD1_C1.isoviz" `
  --launch
```

`--launch` 后必须分别确认：

1. IsoVIZ 窗口确实打开了脚本生成的临时文件；
2. 窗口中的 parent/child 信息、模式标签和 `amp` 与终端期望清单一致。

脚本把临时 `.isoviz` 和 `validation_report.json` 放在系统临时目录。没有
`--launch` 时 GUI 状态为 `not_run`；只完成进程分发、没有人工确认时也不得写成
“IsoVIZ 已正确识别”。这是 GUI 没有机器可读加载回执时的必要边界。

---

## 常见问题

| 现象 | 处理 |
| --- | --- |
| 提示找不到 Java | 安装 JRE/JDK，并把 `java` 加入 PATH |
| 提示找不到 IsoVIZ | 在根目录放 `ISOViz.lnk`，或设置 `ISOVIZ` / `ISOVIZ_JAR`，或关联 `.isoviz` 扩展名 |
| `[matched] 0 mode(s)` | 核对 CSV 的 Mode Name 是否与 `.isoviz` 标签一致；或改用 `Mode=a1,a2,…` 按文件顺序 |
| 改完 CSV 再跑仍像旧的 | 确认 IsoVIZ 打开的是本次启动生成的临时文件，而不是原始输入 `.isoviz` |

官网套件入口：[isotropy.php](https://iso.byu.edu/isotropy.php)。
