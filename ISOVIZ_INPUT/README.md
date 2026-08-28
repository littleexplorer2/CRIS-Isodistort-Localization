# ISOVIZ_INPUT

把振幅表（CSV）里的数值写进官方 IsoVIZ 的子群结构文件（`.isoviz`），并**自动启动** Java 版 IsoVIZ，查看畸变后的晶体结构。本子项目**不使用** `output/`：读完输入文件和参数后直接打开 IsoVIZ，而不是把结果写成仓库里的产物文件。

你不需要先懂晶体学。可以把它理解成：

1. 你已经有一份「每个原子模式该拧到多少」的表格（CSV，常见列名 **Best Model Parameter**）。  
2. 你已经有一份对应子群的 `.isoviz`（里面有各个模式的进度条字段 `amp`）。  
3. 本工具按模式名字（或顺序别名）把 CSV 里的数写进 `amp`，再启动 IsoVIZ。

IsoVIZ 属于 [ISOTROPY Suite](https://iso.byu.edu/isotropy.php)。手动拖动每个模式进度条既慢又容易出错；本子项目用 Python 自动填写。

本项目与仓库根目录的 `CRIS/.venv` 共用一份 Python 虚拟环境。子项目之间的关系见仓库根目录 [README.md](../README.md)。若要先核对 ISODISTORT 导出的 CIF 是否与官网一致，请把官网答案放入 `ISODISTORT_VALIDATE/compare/true/`、把本地 CIF 放入 `compare/item/`，再运行 `ISODISTORT_VALIDATE/main.py`（不要再传入自定义路径）。`compare/` 不入库，缺失时由 `ISODISTORT_VALIDATE/main_requirement.py` 自动创建。批量比较时须把 `compare/true/` 中官网下载的文件改成与 `compare/item/` 相同的相对路径和文件名。

---

## 目录结构

```text
ISOVIZ_INPUT/
  main.py                 读取输入并启动 IsoVIZ
  main_requirement.py     检查/创建 CRIS/.venv，只安装缺失依赖；补齐 input_content/
  requirements.txt        运行时依赖（当前多为标准库即可）
  requirements-dev.txt    开发依赖（pytest）
  pyproject.toml
  README.md
  isoviz_input/           包：CSV 解析、写 .isoviz、启动 IsoVIZ
  tests_dev/              开发测试（用 tests_dev/fixtures/，不依赖你本机样本）
  input_content/          本地输入（gitignore，不入库）
    data.csv/             振幅 CSV 文件夹（目录名就是 data.csv）
    subgroup.isoviz/      官方子群 .isoviz 文件夹（目录名就是 subgroup.isoviz）
```

请把日常用的振幅 CSV 放进 `input_content/data.csv/`，把对应的子群 `.isoviz` 放进 `input_content/subgroup.isoviz/`。整个 `input_content/` **不会上传远程仓库**。缺失时由 `main_requirement.py`（以及 `main.py` 启动时）自动创建这两个子文件夹。

---

## 用户需要配置 / 放置的路径

| 项目 | 默认 / 做法 | 你要做什么 |
| --- | --- | --- |
| **振幅 CSV** | `ISOVIZ_INPUT/input_content/data.csv/` | 把 GD/拟合导出的 CSV 放进该文件夹，或启动时用 `--data <任意路径>` |
| **子群 `.isoviz`** | `ISOVIZ_INPUT/input_content/subgroup.isoviz/` | 放入官网或本地 ISODISTORT 导出的 `data.isoviz`（可改名），或用 `--structure <任意路径>` |
| **IsoVIZ 可执行体** | 见下表「Java 与 IsoVIZ」 | **本机必配其一**：根目录 `ISOViz.lnk`、根目录 `.jar`/`.exe`、或环境变量 `ISOVIZ` / `ISOVIZ_JAR` |
| **Java** | 系统 `PATH` 中的 `java` | 安装 JRE/JDK；用 `.jar` 启动时必需 |
| **临时启动文件** | 系统临时目录（`tempfile`） | 程序自动写入再交给 IsoVIZ；**不要**也不需要配置本子项目的 `output/` |
| **Python / venv** | 仓库根 `CRIS/.venv` | 与其它子项目共用 |

**一般不必改：** 包内 `isoviz_input/paths.py` 的常量（已指向上述 `input_content/`）。若只想换输入文件，优先用命令行参数或往默认文件夹放文件。

跨项目总表见仓库根 [README.md](../README.md)。

---

## 安装

### 1. Python 环境

需要 **Python ≥ 3.10**。在 **CRIS 根目录**执行：

```powershell
cd <CRIS 根目录>
python ISOVIZ_INPUT\main_requirement.py
```

可选：

| 参数 | 作用 |
| --- | --- |
| `--dev` | 额外安装 `requirements-dev.txt`（pytest） |
| `--recreate` | 强制重建 `CRIS/.venv` |

脚本会：

1. 确认 Python 版本  
2. 若缺少 `input_content/`、`input_content/data.csv/`、`input_content/subgroup.isoviz/` 则自动新建  
3. 检查 **Java**（`java` / `javaw` 是否在 PATH；IsoVIZ 是 Java 程序）  
4. 创建或复用 `CRIS/.venv`，只 `pip install` 尚未安装的包  
5. 查找 IsoVIZ 启动方式（见下一小节）

也可以用上游统一安装脚本（同样使用 `CRIS/.venv`，并会补齐上述输入文件夹）：

```powershell
python ISODISTORT\main_requirement.py
```

### 2. Java 与 IsoVIZ（必需）

- 安装 JRE/JDK，保证终端里能运行 `java -version`。  
- 自行安装 IsoVIZ（ISOTROPY Suite 的一部分）。  
- 任选一种方式让本工具找到它：

| 方式 | 说明 |
| --- | --- |
| 快捷方式 | 把快捷方式放到仓库**根目录**，命名为 `ISOViz.lnk`（或 `IsoVIZ.lnk` / `ISOVIZ.lnk`）。该文件已在根 `.gitignore` 中 |
| 可执行文件 / JAR | 根目录放置 `IsoViz.exe` / `ISOViz.exe`，或 `IsoViz.jar` / `ISOViz.jar` / `isoviz.jar` |
| 环境变量 | 设置 `ISOVIZ` 或 `ISOVIZ_JAR` 指向 `.exe` / `.jar` 的完整路径 |
| Windows 文件关联 | 若 `.isoviz` 已关联到 IsoVIZ，脚本也可直接 `startfile` 打开 |

本工具会把填好振幅的内容写成临时文件再交给 IsoVIZ，**不会**在本子项目里维护 `output/`。

---

## 使用：`main.py`

从 CRIS 根目录：

```powershell
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py --data <振幅.csv> --structure <子群.isoviz>
```

### 参数一览

| 参数 | 是否必填 | 含义 |
| --- | --- | --- |
| `--data` | 建议填写 | 振幅 CSV 路径。省略时：若 `input_content/data.csv/` 中有 `.csv` 会先列出供选择，否则提示输入路径 |
| `--structure` | 建议填写 | 子群 `.isoviz` 路径。省略时：若 `input_content/subgroup.isoviz/` 中有 `.isoviz` 会先列出，否则提示输入路径 |

交互示例（不传参数）：

```powershell
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py
```

程序可能打印：

```text
Files in ...\input_content\data.csv:
  1. my_amplitudes.csv
Choose a number or paste a path:
```

运行成功时会显示匹配到的模式数、CSV 里未用到的名字、以及 `.isoviz` 里没有对应 CSV 值而保持原振幅（多为 0）的模式，然后自动启动 IsoVIZ。若一个都没匹配上，会打印该 `.isoviz` 里前若干条模式标签作为提示。

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
cd C:\Users\devou\OneDrive\Desktop\CRIS

# 1) 准备环境（若尚未做过）
python ISOVIZ_INPUT\main_requirement.py

# 2) 把 CSV 和 .isoviz 放进 input_content（也可直接用 --data / --structure 指定任意路径）
#    ISOVIZ_INPUT\input_content\data.csv\best.csv
#    ISOVIZ_INPUT\input_content\subgroup.isoviz\LD1_C1.isoviz

# 3) 读取输入并启动 IsoVIZ
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py `
  --data ISOVIZ_INPUT\input_content\data.csv\best.csv `
  --structure ISOVIZ_INPUT\input_content\subgroup.isoviz\LD1_C1.isoviz
```

---

## 测试

```powershell
cd <CRIS 根目录>
.\.venv\Scripts\python.exe -m pytest ISOVIZ_INPUT\tests_dev
```

开发依赖：

```powershell
python ISOVIZ_INPUT\main_requirement.py --dev
```

测试使用 `tests_dev/fixtures/` 内的样本，不依赖你本机的 `input_content/`。

---

## 常见问题

| 现象 | 处理 |
| --- | --- |
| 提示找不到 Java | 安装 JRE/JDK，并把 `java` 加入 PATH |
| 提示找不到 IsoVIZ | 在根目录放 `ISOViz.lnk`，或设置 `ISOVIZ` / `ISOVIZ_JAR`，或关联 `.isoviz` 扩展名 |
| `[matched] 0 mode(s)` | 核对 CSV 的 Mode Name 是否与 `.isoviz` 标签一致；或改用 `Mode=a1,a2,…` 按文件顺序 |
| 改完 CSV 再跑仍像旧的 | 确认 IsoVIZ 打开的是本次启动生成的临时文件，而不是原始输入 `.isoviz` |

官网套件入口：[isotropy.php](https://iso.byu.edu/isotropy.php)。
