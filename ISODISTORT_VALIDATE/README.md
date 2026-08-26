# ISODISTORT_VALIDATE

把**本地算出的 CIF**与**可信参考 CIF**（通常来自官网导出）放在一起比较，判断两边是不是「同一种晶体结构」。你不需要先学晶体学：工具会检查晶格、原子、坐标、占据率、磁矩、空间群声明等，并给出 **PASS** 或 **FAIL**。

默认看的是**晶体学语义**是否一致，而不是文件排版是否一模一样。若要求字节级完全相同，请加 `--strict`。

本项目与仓库根目录的 `CRIS/.venv` 共用一份 Python 虚拟环境。

---

## 它比较什么

对每一对「本地 CIF ↔ 参考 CIF」，会报告：

| 检查项 | 说明 |
| --- | --- |
| 字节是否完全一致 | 文件内容逐字节相同 |
| UTF-8 文本是否一致 | 文本内容相同（可能与字节略有差别，如换行） |
| 结构与关键元数据 | 晶格参数、原子数、元素、周期分数坐标、原子标签、占据率、磁矩、CIF 里声明的空间群 |
| spglib 推断 | 程序根据坐标独立推断的空间群，以及是否与 CIF 声明一致 |

**适合拿来比的文件**：本地 Distortion 导出的子群 CIF，例如解压 ZIP 后：

```text
isodistort_method2/LD1 C1/LD1 C1 CIF.cif
```

与官网对该子群导出的 CIF 按相同相对路径配对。

**不要拿来当「结构比对」输入的**：

- 结果表的 `.txt` / `.csv`（那是列表，不是结构）  
- 历史遗留的、路径对不上的零散 CIF  

一对文件 PASS ≠ 整个上游工具无 bug；应对多种空间群、模式、超胞做批量回归。

---

## 安装

需要 **Python ≥ 3.10**。在仓库根目录执行（会创建/复用 `CRIS/.venv`，并安装本目录 `requirements.txt`）：

```powershell
cd "C:\Users\devou\OneDrive\Desktop\CRIS"
python ISODISTORT\main_requirement.py
```

开发测试额外依赖（与上游共用 `--dev`）：

```powershell
python ISODISTORT\main_requirement.py --dev
```

之后请用虚拟环境里的解释器运行本目录脚本，例如：

```powershell
.\.venv\Scripts\python.exe ISODISTORT_VALIDATE\main.py
```

或先 `cd ISODISTORT_VALIDATE` 再运行 `python main.py`（需已激活该 venv）。

---

## 终端菜单：`main.py`

```powershell
cd ISODISTORT_VALIDATE
python main.py
```

菜单：

| 选项 | 作用 |
| --- | --- |
| **1. 比较一对 CIF** | 输入本地路径、参考路径、容差、是否忽略原子顺序、是否严格字节比较、可选参考 SHA-256 |
| **2. 批量回归验证** | 输入两个目录、文件匹配模式（默认 `*.cif`）、同样的容差选项、可选 hash manifest、是否输出 JSON |
| **3. 查看验证说明** | 打印本工具摘要（细节以本 README 为准） |
| **0. 退出** | 结束 |

菜单里询问的「晶格容差 / 分数坐标容差 / 占据率与磁矩容差」默认都是 `1e-5`。

---

## 命令行：单对比较 `compare_cif.py`

```powershell
cd ISODISTORT_VALIDATE
python compare_cif.py "本地.cif" "参考.cif"
```

路径建议用本机**绝对路径**。示例：

```powershell
python compare_cif.py `
  "C:\path\to\isodistort_method2\LD1 C1\LD1 C1 CIF.cif" `
  "C:\path\to\official.cif"
```

### 全部参数

| 参数 | 默认 | 含义 |
| --- | --- | --- |
| `local_cif` | （必填位置参数） | 本地 CIF |
| `reference_cif` | （必填位置参数） | 参考 / 官网 CIF |
| `--lattice-tol` | `1e-5` | 晶格参数容差 |
| `--coord-tol` | `1e-5` | 分数坐标（周期）容差 |
| `--scalar-tol` | `1e-5` | 占据率、磁矩等容差 |
| `--ignore-atom-order` | 关 | 按「元素 + 周期坐标」匹配原子，而不是按文件里的行顺序 |
| `--reference-sha256` | 无 | 参考文件的期望 SHA-256；不匹配则失败（防止拿错参考文件） |
| `--strict` | 关 | 除语义一致外，还要求**字节完全一致**才 PASS |
| `--structure-only` | — | **已废弃别名**；默认本来就是语义比较 |
| `--json` | 关 | 打印完整 JSON 结果 |

容差示例：

```powershell
python compare_cif.py local.cif official.cif `
  --lattice-tol 1e-5 --coord-tol 1e-5 --scalar-tol 1e-5
```

可信参考哈希：

```powershell
Get-FileHash "C:\path\to\official.cif" -Algorithm SHA256
python compare_cif.py local.cif official.cif --reference-sha256 "这里填写64位十六进制"
```

---

## 命令行：批量比较 `batch_compare.py`

两个目录按**相对路径**配对。例如：

```text
local_cases/LD1 C1/LD1 C1 CIF.cif
reference_cases/LD1 C1/LD1 C1 CIF.cif
```

```powershell
python batch_compare.py "C:\path\to\local_cases" "C:\path\to\reference_cases"
python batch_compare.py local_cases reference_cases --json > report.json
python batch_compare.py local_cases reference_cases --pattern "*CIF.cif"
```

### 全部参数

| 参数 | 默认 | 含义 |
| --- | --- | --- |
| `local_dir` | （必填） | 本地 CIF 根目录 |
| `reference_dir` | （必填） | 参考 CIF 根目录 |
| `--pattern` | `*.cif` | 递归匹配模式 |
| `--lattice-tol` / `--coord-tol` / `--scalar-tol` | `1e-5` | 同单对比较 |
| `--ignore-atom-order` | 关 | 同单对比较 |
| `--hash-manifest` | 无 | JSON：键为相对路径，值为参考文件 SHA-256 |
| `--strict` | 关 | 每对还要求字节一致 |
| `--json` | 关 | 机器可读汇总 |

manifest 示例：

```json
{
  "LD1 C1/LD1 C1 CIF.cif": "参考文件的64位SHA256"
}
```

文本模式下每行形如 `PASS …` / `FAIL …`，最后有：

```text
Summary: total=…, passed=…, failed=…
```

某相对路径只在一侧存在时，记为失败，原因含 `missing matching CIF`。

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
  main.py              交互菜单
  compare_cif.py       单对比较 CLI
  batch_compare.py     批量比较 CLI
  requirements.txt
  tests_dev/           pytest
  README.md
```
