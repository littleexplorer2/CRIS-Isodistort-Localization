# ISODISTORT_VALIDATE

比较 **ISODISTORT 本地导出的 CIF** 与 **官网导出的参考 CIF** 是否为同一结构，用来检查本地化有没有算错。本目录可以改；`webpage_info/`、`experiment_data/` 与 `GD/` 不能改。

默认 **PASS** 看的是晶体学语义（晶格、原子、坐标、占据率、磁矩、声明空间群），不是文件排版。字节级回归用 `--strict`。

与 `ISODISTORT` 共用仓库根目录的 `CRIS/.venv`。

---

## 和本地批量导出如何对应

网页与终端界面为英语。四个 Method 只计算；下载在 Distortion：

- **子群结构文件**（Method 1/2/3）：每次只含 **一个 Method** 的子群；若结果表有 Filter，ZIP 只含命中行，CIF 形如 `isodistort_method2/LD1 C1/LD1 C1 CIF.cif`。
- **结果表 txt/csv**（Method 1–4）：是筛选后的列表，不是结构 CIF，不要拿来和官网第 6 页 CIF 做结构比较。

网页已不再提供 Generate（有限振幅畸变 CIF）；不要用 `output/` 里旧的 `mixed_*.cif` 去对官网第 6 页参考集。

```text
isodistort_method2.zip
  isodistort_method2/
    LD1 C1/LD1 C1 CIF.cif
    LD5 P6/LD5 P6 CIF.cif
```

与官网第 6 页对该子群导出的 CIF 按相对路径配对比较即可。不要拿 `output/` 里旧的 `mixed_*.cif`、其它 Method 的结果、或终端单文件导出去对整批官网参考集。

Method 4 是畸变结构分解，没有子群列表可打包。Python API 若仍调用 `generate_distortion` 写出的 CIF，可用 `compare_cif.py` 与官网对应文件一对一比较，不要拿来充当 Method ZIP 里的子群 CIF。

---

## 安装

```powershell
cd "C:\Users\devou\OneDrive\Desktop\CRIS"
python ISODISTORT\main_requirement.py
```

开发依赖（本目录 `tests_dev/`）与 ISODISTORT 共用 `--dev`。下面命令在 `ISODISTORT_VALIDATE/` 下运行，并已激活该 venv。

---

## 终端菜单

```powershell
python main.py
```

- 比较一对本地 / 参考 CIF
- 按相对路径批量比较两个目录（可改文件匹配模式，默认 `*.cif`）
- 设置晶格、坐标、占据率/磁矩容差
- 是否忽略原子顺序、是否严格字节比较
- 参考文件 SHA-256 或批量 hash manifest
- 文本或 JSON 报告

菜单第 3 项「查看验证说明」是本工具的摘要，详细约定以本 README 为准。

---

## 单对比较

两个输入都支持本机绝对路径。批量 ZIP 解压后的 CIF 示例：

```powershell
python compare_cif.py `
  "C:\path\to\isodistort_method2\LD1 C1\LD1 C1 CIF.cif" `
  "C:\path\to\official.cif"
```

容差：

```powershell
python compare_cif.py local.cif official.cif `
  --lattice-tol 1e-5 --coord-tol 1e-5 --scalar-tol 1e-5
```

原子顺序不同但元素+周期坐标应能对上时加 `--ignore-atom-order`。要求文件逐字节相同加 `--strict`。`--json` 打印完整结果。`--structure-only` 是已废弃别名，默认本来就是语义比较。

### 可信参考文件哈希

先算官网参考文件的 SHA-256，再传给比较器，避免拿错参考文件：

```powershell
Get-FileHash "C:\path\to\official.cif" -Algorithm SHA256
python compare_cif.py local.cif official.cif --reference-sha256 "这里填写官网文件的SHA256"
```

哈希不匹配会失败。

---

## 批量比较

两个目录按**相对路径**配对。把 ZIP 里的 `isodistort_methodN/` 当作本地根，官网参考目录做成同样的子群文件夹结构：

```text
local_cases/LD1 C1/LD1 C1 CIF.cif
reference_cases/LD1 C1/LD1 C1 CIF.cif
```

```powershell
python batch_compare.py "C:\path\to\local_cases" "C:\path\to\reference_cases"
python batch_compare.py local_cases reference_cases --json > report.json
python batch_compare.py local_cases reference_cases --pattern "*CIF.cif"
```

会报告总用例数、通过/失败数、缺失配对、以及每对的字节/结构/哈希结果。可用 `--hash-manifest manifest.json`（键为相对路径，值为参考文件 SHA-256）：

```json
{
  "LD1 C1/LD1 C1 CIF.cif": "参考文件的64位SHA256"
}
```

---

## 比较层次与返回码

`compare_cif.py` 会同时给出：

- 字节是否完全一致
- UTF-8 文本是否一致
- 结构与关键元数据是否一致（晶格、原子数、元素、周期分数坐标、原子标签、占据率、磁矩、CIF 声明空间群）
- spglib 独立推断的空间群（以及是否与 CIF 声明一致）

| 码 | 含义 |
| --- | --- |
| 0 | 语义比较通过（`--strict` 时还要求字节一致） |
| 1 | 解析成功但语义/哈希失败，或批量有失败用例 |
| 2 | 路径、参数、依赖或 CIF 解析错误 |

解读：

- `byte_exact=否` 且 `structure_equal=是`：多半是排版、字段顺序、浮点写法，不单独当 bug。
- `lattice differs` / `fractional coordinates differ`：查超胞、幅度、原子匹配。
- `occupancies` / `magnetic moments`：非几何量可能有问题。
- `spglib-inferred space group differs`：两边对称性判断不同，优先查坐标。
- CIF 声明 `P1` 但 spglib 给出更高对称性：可能是导出未对称化，要结合官网文件判断，不能单凭这一项断定本地算错。

一对文件通过 ≠ 整个 ISODISTORT 无 bug。应对不同空间群、模式、超胞和幅度做批量参考集。

---

## 测试

```powershell
cd ISODISTORT_VALIDATE
python -m pytest tests_dev -q
```
