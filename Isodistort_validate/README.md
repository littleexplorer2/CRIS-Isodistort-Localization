# Isodistort Validate

用于比较 `Isodistort_back` 生成的 CIF 与官网导出的参考 CIF，帮助发现本地化 ISODISTORT 的输出回归。

## 检查层次

`compare_cif.py` 同时给出以下结果：

- **字节完全一致**：二进制内容逐字节相同，仅作为严格格式检查。
- **文本内容一致**：按 UTF-8 读取后内容相同。
- **结构与关键元数据一致**：比较晶格、原子数、元素、周期分数坐标、原子标签、占据率、磁矩和 CIF 声明空间群。
- **spglib 独立诊断**：分别从两个解析结构推断空间群，并检查两边推断是否一致；同时报告 CIF 声明空间群是否与独立推断相符。

默认 `PASS` 的标准是语义结构与关键元数据一致，且参考文件哈希（如果提供）正确。CIF 排版、行尾、字段顺序等差异不会再导致默认误报。需要逐字节回归时使用 `--strict`。

## 安装

```powershell
cd "C:\Users\devou\OneDrive\Desktop\CRIS\Isodistort_validate"
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

也可以使用已经安装 `numpy`、`pymatgen`、`spglib` 的 `Isodistort_back\.venv`。

## 单对 CIF 比较

两个输入都支持本机绝对路径：

```powershell
python compare_cif.py `
  "C:\path\to\Isodistort_back\output\local.cif" `
  "C:\path\to\official\official.cif"
```

默认容差分别用于不同物理量：

```powershell
python compare_cif.py "C:\path\to\local.cif" "C:\path\to\official.cif" `
  --lattice-tol 1e-5 `
  --coord-tol 1e-5 `
  --scalar-tol 1e-5
```

如果两个文件的原子排序不同，但希望按“元素 + 周期坐标”匹配：

```powershell
python compare_cif.py "C:\path\to\local.cif" "C:\path\to\official.cif" --ignore-atom-order
```

严格要求文件字节完全一致：

```powershell
python compare_cif.py "C:\path\to\local.cif" "C:\path\to\official.cif" --strict
```

### 可信参考文件哈希

先计算官网参考文件的 SHA-256：

```powershell
Get-FileHash "C:\path\to\official.cif" -Algorithm SHA256
```

再把哈希传给比较器：

```powershell
python compare_cif.py "C:\path\to\local.cif" "C:\path\to\official.cif" `
  --reference-sha256 "这里填写官网文件的SHA256"
```

哈希不匹配会返回失败，避免误把错误的参考文件当作权威基准。

## 批量回归

批量模式按照相对路径配对两个目录中的 CIF。例如：

```text
local_cases/
  tetragonal/case01.cif
reference_cases/
  tetragonal/case01.cif
```

运行：

```powershell
python batch_compare.py `
  "C:\path\to\local_cases" `
  "C:\path\to\reference_cases"
```

输出 JSON 汇总：

```powershell
python batch_compare.py `
  "C:\path\to\local_cases" `
  "C:\path\to\reference_cases" `
  --json > report.json
```

批量模式会报告：

- 总用例数、通过数、失败数
- 缺失的配对文件
- 每对文件的字节结果、结构结果、哈希结果和差异原因

可用 `--hash-manifest manifest.json` 校验可信参考集。manifest 格式：

```json
{
  "tetragonal/case01.cif": "参考文件的64位SHA256"
}
```

## 返回码

- `0`：比较通过。
- `1`：文件成功解析但语义比较失败、哈希失败或批量存在失败用例。
- `2`：路径、参数、依赖或 CIF 解析错误。

## 如何解读结果

- `byte_exact=否`、`structure_equal=是`：通常只是 CIF 排版、字段顺序、浮点格式或非关键文本不同，不自动视为 bug。
- `lattice differs`：检查晶胞参数、超胞设置和晶胞表示。
- `fractional coordinates differ`：检查位移映射、畸变幅度、原子匹配和周期坐标。
- `occupancies` 或 `magnetic moments` 差异：说明非几何物理量可能有 bug。
- `spglib-inferred space group differs`：两份结构的独立对称性判断不同，应优先检查坐标和结构设置。
- `declared space group does not match spglib inference`：这是独立诊断告警；某些 CIF 会显式声明 `P1`，但坐标实际具有更高对称性，因此需要结合官网文件和预期空间群判断，不能仅凭这一项断言本地程序错误。

一个参考文件通过不等于整个 ISODISTORT 无 bug。应建立覆盖不同空间群、畸变模式、超胞、磁结构和幅度的批量参考集。
