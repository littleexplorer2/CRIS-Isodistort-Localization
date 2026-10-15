# ISOVIZ_INPUT 0.4.0

ISOVIZ_INPUT 是位于 GD 振幅结果与外部 IsoVIZ 软件之间的**本地输入适配器**。用户在桌面界面中上传文件或粘贴文本；程序识别振幅表与官方子群 `.isoviz`，只修改临时副本中的模式 `amp`，通过静态门禁后再请求 IsoVIZ 打开。

它不搜索子群、不比较 CIF、不拟合振幅，也不把“进程已分发”表述为“IsoVIZ 已正确识别数据”。写入值来自 CSV 的 **Best Model Parameter**，不是 Normalized Amplitude。

## 程序结构与职责

```text
ISOVIZ_INPUT/
  main.py                       唯一运行入口；默认 GUI，--cli 为命令行
  main_requirement.py           统一安装器的兼容入口
  .gitignore                    本地输入、缓存和历史 output 排除规则
  config/settings.yaml          路径和 IsoVIZ 启动器查找规则
  isoviz_input/
    frontend/                   前端：仅调用应用层公开接口
      gui/                      上传、粘贴、报告和用户动作
      cli/                      脚本化/兼容命令行适配
      entrypoint/               GUI / CLI 统一入口选择
    application/                应用层：用例编排及输入/环境查询
      pipeline/                 唯一用例管线与启动门禁
      input_selection/          文件/文本选择及输入路径查询接口
      environment/              IsoVIZ 启动器可用性查询
    backend/                    底层：独立能力模块
      inputs/                   文件/粘贴文本统一为 TextDocument
      amplitudes/               CSV 识别、标签匹配和 amp 写入
      validation/               身份、非振幅内容和写入值不变量
      launcher/                 .lnk/.exe/.jar/文件关联启动适配
      paths/                    Best_Model_Parameters 与用户路径解析
      config_loader/            settings.yaml 读取
    amplitudes.py 等             旧导入路径的薄兼容层；不保存第二套算法
  tests_dev/                    自动测试与显式手工 GUI 验收器
```

每个职责目录通过 `__init__.py` 暴露公开接口，实现保存在目录内独立文件中。
前端只调用应用层，应用层调用底层；底层不能反向依赖。包根旧接口只转发，不保存第二套实现。

主要处理管线：

```text
GUI / CLI
   │  PipelineRequest(TextDocument, TextDocument)
   ▼
application.pipeline
   ├─► backend.amplitudes  ─► 匹配报告 + patched text
   ├─► backend.validation  ─► pass / pass_with_warnings / fail
   ├─► 系统临时 .isoviz（fail 时不生成可启动文件）
   └─► backend.launcher    ─► 外部 IsoVIZ 进程分发
```

`application.input_selection` 统一前端文件/文本选择，`application.environment` 提供启动器查询；
因此上传与粘贴共享同一解析器，GUI 与 CLI 共享同一管线；任何新输入方式都不能绕过静态验证或直接调用启动器。

日常 CSV 通常来自同级目录：

```text
<CRIS 与 GD 的同级目录>/Best_Model_Parameters/
  <irrep>_<structure_type>/<irrep>_<structure_type>_best_model_parameters.csv
```

默认 GUI 启动器统一放在仓库根 `local_tools/launchers/`，实际快捷方式、JAR 和 EXE 由根 `.gitignore` 排除。源 CSV、源 `.isoviz` 和同级 `Best_Model_Parameters` 都不会被本程序改写。

## 用户需要配置 / 放置的路径

| 项目 | 默认 / 做法 | 用户动作 |
| --- | --- | --- |
| 振幅 CSV | GUI 选择文件或粘贴完整 CSV；文件对话框优先打开同级 `Best_Model_Parameters` | 使用 GD 导出的 Best Model Parameter |
| 子群 `.isoviz` | GUI 选择文件或粘贴完整文本 | 使用与 CSV 相同 IR / OPD / path 的官方文件 |
| IsoVIZ | `config/settings.yaml` 的 `isoviz` 段 | 放到 `local_tools/launchers/`，或设置 `ISOVIZ` / `ISOVIZ_JAR` |
| Java | 系统 `PATH` | 仅从 `.jar` 启动时必需 |
| 临时文件 | 系统临时目录 | 自动生成；不需要仓库 `output/` |

项目与仓库根目录的 `.venv` 共用 Python 环境。更换目录布局时可设置 `BEST_MODEL_PARAMETERS_DIR`；运行时默认值以 [config/settings.yaml](config/settings.yaml) 为准。

源码运行入口仍为 `main.py`；安装 wheel 后可执行 `isoviz-input`，两者调用相同前端入口。
构建过程从 `config/settings.yaml` 自动生成包内配置资源，源码中只维护这一份配置。
独立安装包不包含本机启动器，使用 `ISOVIZ` / `ISOVIZ_JAR` 指定启动器，
使用 `BEST_MODEL_PARAMETERS_DIR` 指定外部振幅目录。

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
4. 检查 Java、IsoVIZ 启动器以及解析到的 `Best_Model_Parameters` 中是否已有 CSV。缺少用户 CSV 记为 `WARN`，不会伪装成 GUI 已验证。

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
| 快捷方式 | 把快捷方式放到仓库根的 **`local_tools/launchers/`**，命名为 `ISOViz.lnk`（或 `IsoVIZ.lnk` / `ISOVIZ.lnk`）。目录内容已在根 `.gitignore` 中排除；查找顺序见 `config/settings.yaml` |
| 可执行文件 / JAR | 在 `local_tools/launchers/` 中放置 `IsoViz.exe` / `ISOViz.exe`，或 `IsoViz.jar` / `ISOViz.jar` / `isoviz.jar` |
| 环境变量 | 设置 `ISOVIZ` 或 `ISOVIZ_JAR` 指向 `.exe` / `.jar` 的完整路径 |
| Windows 文件关联 | 若 `.isoviz` 已关联到 IsoVIZ，脚本也可直接 `startfile` 打开 |

本工具会把填好振幅的内容写成临时文件再交给 IsoVIZ，**不会**在本子项目里维护 `output/`。

---

## 使用：桌面界面（默认）

从 CRIS 根目录运行：

```powershell
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py
```

界面左右各有一个输入区：

1. 在“振幅 CSV”中点击“选择文件…”上传 CSV，或把完整 CSV 文本复制后点“粘贴剪贴板”。
2. 在“子群 `.isoviz`”中用同样方式上传或粘贴官方结构内容。
3. 先点“检查并生成临时文件”。报告会列出匹配模式、未使用 CSV 行、保留原值的 IsoVIZ 模式和静态校验状态。
4. 状态不是 `fail` 时，再点“生成并启动 IsoVIZ”。
5. 在 IsoVIZ 中人工核对 parent/child、模式标签和 `amp`。界面显示“已请求启动”只证明进程分发，不证明窗口已打开或数据已识别。

输入框内容是本次运行的实际输入：选择文件后仍可手动修订，管线会处理修订后的文本。每次准备生成新的系统临时副本；源文件不变。

## 使用：命令行兼容方式

脚本化调用显式加 `--cli`：

```powershell
# 生成、校验并启动
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py --cli `
  --data "E:\Best_Model_Parameters\LD1_C1\LD1_C1_best_model_parameters.csv" `
  --structure "D:\data\LD1_C1.isoviz"

# 只准备临时文件，不启动外部软件
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py --cli `
  --data "E:\Best_Model_Parameters\LD1_C1\LD1_C1_best_model_parameters.csv" `
  --structure "D:\data\LD1_C1.isoviz" `
  --prepare-only
```

为兼容原有脚本，直接传 `--data` / `--structure` 也会进入 CLI。没有文件参数但指定 `--cli` 时，才使用原来的文件夹/文件名终端询问流程。

| 参数 | 含义 |
| --- | --- |
| `--gui` | 显式启动桌面界面 |
| `--cli` | 使用命令行适配器 |
| `--folder` + `--csv-name` | 从 `Best_Model_Parameters` 定位 CSV；必须成对提供 |
| `--data` | 振幅 CSV 完整路径 |
| `--structure` | `.isoviz` 完整路径 |
| `--prepare-only` | 只做解析、静态校验和临时文件生成 |

### CSV 需要什么列

CSV **必须有表头**。识别列名时不区分大小写。常用（与梯度下降脚本写出的表兼容）：

| 列名（示例） | 用途 |
| --- | --- |
| **Mode Name**（或 `modelabel` / `label` / `name`） | 与 `.isoviz` 中的模式标签匹配，例如 `[0,0,1/6]LD1[Eu1:a:dsp]A2u(a)` |
| **Best Model Parameter**（或 `amplitude` / `amp` / `value`） | 写入 IsoVIZ 进度条使用的 **`amp`**（与 `maxamp` 同一单位）。**不要**误用 Normalized Amplitude |
| **Mode** | 备选：`a1`, `a2`, … 按 `.isoviz` 文件中 **strain 模式再 displacive 模式**的出现顺序对应 |
| Maximum Mode Amplitude / maxamp | 可选；解析会读入，匹配主要仍靠名字或 Mode 别名 |

模式名里可能含逗号；推荐用带引号的 CSV（例如 pandas 默认写出格式）。未匹配到的 IsoVIZ 模式保持原振幅。

### 建议你第一次这样试用

1. 先运行 GD 笔记本的保存单元，确认
   `Best_Model_Parameters\LD1_C1\LD1_C1_best_model_parameters.csv`
   实际存在。
2. 使用与这份 CSV **同一个 IR/OPD/path** 导出的官方 `.isoviz`；不能把
   `LD1_C1` 振幅写进另一个子群文件后再据此判断结构是否正确。
3. 启动 GUI，分别上传两个文件；也可打开文件后复制完整文本粘贴。
4. 先看界面报告：匹配数必须大于 0；“CSV 中未使用的模式”和“保留原振幅的
   IsoVIZ 模式”不一定都是错误，但在判定结构前必须逐项解释。
5. IsoVIZ 打开后，核对母相/子群、模式总数和标签，并抽查界面报告列出的 `amp`；
   然后观察结构是否随这些振幅变化。`Best Model Parameter` 是写入值，
   `Normalized Amplitude` 不是。

如果默认查找位置没有 IsoVIZ 启动器，仅有 Java 还不够；请把 `ISOViz.lnk`
放到 CRIS 根目录的 `local_tools/launchers/`，或者设置 `ISOVIZ` / `ISOVIZ_JAR`。源
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
入口回归会验证真实 CLI 的 `--prepare-only` 流程、源输入保持、临时 wheel 构建与
离开源码目录后的安装包入口；构建依赖 `setuptools` 和 `wheel` 已列入开发依赖。
这些检查不启动外部 IsoVIZ，也不证明 GUI 已识别结构。

### GUI 打开与数据识别验收

普通 pytest 自动覆盖四层：上传与粘贴输入的统一读取、CSV/标签匹配与 `amp` 写入、
补丁后模式身份及非振幅内容保持、JAR/EXE/Windows 快捷方式是否把临时 `.isoviz`
作为明确参数传给 IsoVIZ。零匹配会阻止启动。这些测试不弹出窗口。

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
| 提示找不到 IsoVIZ | 在 `local_tools/launchers/` 放 `ISOViz.lnk`，或设置 `ISOVIZ` / `ISOVIZ_JAR`，或关联 `.isoviz` 扩展名 |
| `[matched] 0 mode(s)` | 核对 CSV 的 Mode Name 是否与 `.isoviz` 标签一致；或改用 `Mode=a1,a2,…` 按文件顺序 |
| 改完 CSV 再跑仍像旧的 | 确认 IsoVIZ 打开的是本次启动生成的临时文件，而不是原始输入 `.isoviz` |

官网套件入口：[isotropy.php](https://iso.byu.edu/isotropy.php)。
