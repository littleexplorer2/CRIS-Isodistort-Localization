# ISODISTORT_VALIDATE

把**本地算出的 CIF**与**官网导出的参考 CIF**放在固定目录里比较，判断两边是不是「同一种晶体结构」。你不需要先学晶体学：工具会检查晶格、原子、坐标、占据率、磁矩、空间群声明等，并给出 **PASS** 或 **FAIL**。

比较目录是固定的，**不要再传入自定义路径**：

```text
ISODISTORT_VALIDATE/compare/
  true/     官网下载的标准答案 CIF
  item/     需要验证的本地 CIF
```

两侧按**相对路径**配对（含子目录）。例如：

```text
compare/item/LD1 C1/subgroup.cif
compare/true/LD1 C1/subgroup.cif
```

**批量比较时请务必改名：** 把官网下载到 `compare/true/` 的 CIF 改成与 `compare/item/` 中本地文件**完全相同**的相对路径和文件名，才能一一对应。程序启动后只要发现两侧有对不上的文件名，会立刻列出这些文件并提醒你改名。

默认看的是**晶体学语义**是否一致，而不是文件排版是否一模一样。若要求字节级完全相同，请加 `--strict`。

本项目与仓库根目录的 `CRIS/.venv` 共用一份 Python 虚拟环境。单对比较和批量比较都从同一个入口启动：`main.py`。

---

## 它比较什么

对每一对「`compare/item` 中的 CIF ↔ `compare/true` 中的同名 CIF」，会报告：

| 检查项 | 说明 |
| --- | --- |
| 字节是否完全一致 | 文件内容逐字节相同 |
| UTF-8 文本是否一致 | 文本内容相同（可能与字节略有差别，如换行） |
| 结构与关键元数据 | 晶格参数、原子数、元素、周期分数坐标、原子标签、占据率、磁矩、CIF 里声明的空间群 |
| spglib 推断 | 程序根据坐标独立推断的空间群，以及是否与 CIF 声明一致 |

**适合拿来比的文件**：本地 Distortion 导出的子群 CIF，例如解压 ZIP 后：

```text
LD1 C1/subgroup.cif
```

把该文件拷到 `compare/item/` 下（可保留子目录），把官网第 6 页对该子群导出的 CIF（通常也叫 `subgroup.cif`）拷到 `compare/true/` 下同一相对路径。

**不要拿来当「结构比对」输入的**：

- 结果表的 `.txt` / `.csv`（那是列表，不是结构）  
- 历史遗留的、路径对不上的零散 CIF  

`compare/`（含 `item/` 与 `true/`）**不会上传远程仓库**（已 gitignore）。一对文件 PASS ≠ 整个上游工具无 bug；应对多种空间群、模式、超胞做批量回归。

---

## 安装

需要 **Python ≥ 3.10**。本项目与仓库根目录的 `CRIS/.venv` 共用一份虚拟环境。在仓库根目录执行本子项目的准备脚本：

```powershell
cd "C:\Users\devou\OneDrive\Desktop\CRIS"
python ISODISTORT_VALIDATE\main_requirement.py
```

该脚本会：

1. 确认 Python ≥ 3.10  
2. 创建或复用 `CRIS/.venv`，只安装尚未存在的依赖（已下载的包不会重新下载）  
3. 检查运行本工具**不需要**额外环境变量（无 ISODATA / WSL）  
4. 若缺少 `compare/`、`compare/item/`、`compare/true/` 则自动创建  

开发测试额外依赖：

```powershell
python ISODISTORT_VALIDATE\main_requirement.py --dev
```

也可以继续用上游统一安装脚本（同样会创建/复用 `.venv`，并补上 VALIDATE 的 compare 目录）：

```powershell
python ISODISTORT\main_requirement.py
python ISODISTORT\main_requirement.py --dev
```

之后请用虚拟环境里的解释器运行本目录脚本，例如：

```powershell
.\.venv\Scripts\python.exe ISODISTORT_VALIDATE\main.py
```

或先 `cd ISODISTORT_VALIDATE` 再运行 `python main.py`（需已激活该 venv）。

---

## 启动入口：`main.py`

单对比较和批量比较都走这一个文件。

```powershell
cd ISODISTORT_VALIDATE
python main.py
```

无参数时进入菜单：

| 选项 | 作用 |
| --- | --- |
| **1. 比较一对 CIF** | 从 `compare/item` 与 `compare/true` 列出相对路径并选择一对；再询问容差、是否忽略原子顺序、是否严格字节比较、可选参考 SHA-256 |
| **2. 批量回归验证** | 比较两个固定文件夹中的全部配对；询问匹配模式（默认 `*.cif`）、同样的容差选项、可选 hash manifest、是否输出 JSON。**使用前请把 `true/` 中官网文件改名，使其与 `item/` 一一对应** |
| **3. 查看验证说明** | 打印本工具摘要（细节以本 README 为准） |
| **0. 退出** | 结束 |

菜单**不再询问文件或目录路径**。晶格 / 分数坐标 / 占据率与磁矩容差默认都是 `1e-5`。

进入菜单、单对比较或批量比较时，若 `item/` 与 `true/` 里有对不上的相对路径，会立刻打印这些文件名，并请你把 `true/` 中官网 CIF 改成与 `item/` 相同。没有成对文件时，单对/批量比较不会继续问容差。

### 「是否忽略原子排列顺序」是什么意思

CIF 里每个原子占一行。这个选项决定两边原子怎么对上号：

| 回答 | 含义 |
| --- | --- |
| **n（默认）** | 按文件里的行号一一对应：第 1 行对第 1 行，第 2 行对第 2 行。元素种类或分数坐标对不上，或者只是行顺序不同，都会判为不一致。两边导出格式接近时用这个。 |
| **y** | 不按行号，按「同种元素 + 分数坐标足够接近」配对。同一套原子只是写成不同顺序时（例如 `Eu, Al, Al` 对 `Al, Eu, Al`）可以选 y。这**不会**忽略坐标、元素种类、占据率或磁矩的真实差异。 |

命令行对应开关是 `--ignore-atom-order`（加上即相当于回答 y）。

### 命令行（同一入口）

```powershell
python main.py compare "LD1 C1/subgroup.cif"
python main.py compare
python main.py batch
python main.py batch --json > report.json
python main.py batch --pattern "subgroup.cif"
```

`compare` 在 `compare/item` 里只有一个 CIF 时可以省略相对路径。`batch` 始终比较固定目录，按相对路径配对。

### `compare` 参数

| 参数 | 默认 | 含义 |
| --- | --- | --- |
| `relative_path` | （可省略；仅当 item 中只有一个 CIF 时） | 相对于 `compare/item` 与 `compare/true` 的路径 |
| `--lattice-tol` | `1e-5` | 晶格参数容差 |
| `--coord-tol` | `1e-5` | 分数坐标（周期）容差 |
| `--scalar-tol` | `1e-5` | 占据率、磁矩等容差 |
| `--ignore-atom-order` | 关 | 忽略 CIF 原子**行顺序**，按「元素 + 分数坐标」配对；默认仍按行号对应。不会忽略真实的坐标或元素差异 |
| `--reference-sha256` | 无 | 参考文件的期望 SHA-256；不匹配则失败（防止拿错参考文件） |
| `--strict` | 关 | 除语义一致外，还要求**字节完全一致**才 PASS |
| `--structure-only` | — | **已废弃别名**；默认本来就是语义比较 |
| `--json` | 关 | 打印完整 JSON 结果 |

容差示例：

```powershell
python main.py compare "LD1 C1/subgroup.cif" `
  --lattice-tol 1e-5 --coord-tol 1e-5 --scalar-tol 1e-5
```

可信参考哈希：

```powershell
Get-FileHash "ISODISTORT_VALIDATE\compare\true\sample.cif" -Algorithm SHA256
python main.py compare sample.cif --reference-sha256 "这里填写64位十六进制"
```

### `batch` 参数

| 参数 | 默认 | 含义 |
| --- | --- | --- |
| `--pattern` | `*.cif` | 递归匹配模式 |
| `--lattice-tol` / `--coord-tol` / `--scalar-tol` | `1e-5` | 同单对比较 |
| `--ignore-atom-order` | 关 | 同单对比较 |
| `--hash-manifest` | 无 | JSON：键为相对路径，值为参考文件 SHA-256 |
| `--strict` | 关 | 每对还要求字节一致 |
| `--json` | 关 | 机器可读汇总 |

manifest 示例：

```json
{
  "LD1 C1/subgroup.cif": "参考文件的64位SHA256"
}
```

文本模式下每行形如 `PASS …` / `FAIL …`，最后有：

```text
Summary: total=…, passed=…, failed=…
```

某相对路径只在一侧存在时，记为失败，原因含 `missing matching CIF`。这通常就是 `true/` 里官网文件名还没改成与 `item/` 一致。

---

## PASS / FAIL 与退出码

### 默认（不加 `--strict`）

- **PASS**：解析后的结构与关键元数据一致，且（若提供了）参考哈希校验通过。  
- **FAIL**：语义不一致，或参考哈希不匹配。

字节不同但结构相同 → **仍可 PASS**（常见于排版、字段顺序、浮点写法差异）。

### 加 `--strict`

- **PASS**：语义一致 **并且** 字节完全一致。  
- 否则 **FAIL**。

### 进程退出码

| 码 | 含义 |
| --- | --- |
| **0** | 通过（默认语义；`--strict` 时还要求字节一致） |
| **1** | 解析成功但比较失败，或批量中有失败用例 |
| **2** | 路径、参数、依赖或 CIF 解析错误 |

### 读结果时的提示

| 现象 | 通常意味着 |
| --- | --- |
| `byte_exact=否` 且结构一致 | 排版差异，一般不当作算错 |
| 晶格 / 分数坐标差异 | 查超胞、坐标变换、原子匹配 |
| 占据率 / 磁矩差异 | 非几何量可能有问题 |
| spglib 推断不一致 | 优先查坐标；CIF 声明 `P1` 但 spglib 给出更高对称性时，可能是导出未对称化，需结合参考文件判断 |

---

## 测试

```powershell
cd ISODISTORT_VALIDATE
python -m pytest tests_dev -q
```

或从仓库根目录：

```powershell
.\.venv\Scripts\python.exe -m pytest ISODISTORT_VALIDATE\tests_dev -q
```

---

## 本目录文件

```text
ISODISTORT_VALIDATE/
  main.py                 唯一启动入口（交互菜单 / compare / batch）
  main_requirement.py     依赖 / 环境 / compare 目录准备
  isodistort_validate/    计算核心（比较算法、固定路径）
  compare/                本地比对目录（gitignore，不入库）
    true/                 官网标准答案 CIF（批量比较前须改名以匹配 item/）
    item/                 待验证 CIF
  requirements.txt
  requirements-dev.txt
  pyproject.toml
  tests_dev/              pytest
  README.md
```
