# ISODISTORT 官网与本地输出下载清单

本清单把 `docs/manifests/` 中的机器 JSON 转写为人工操作步骤。JSON 是脚本读取的精确数据源；本文件用于在官网或本地网页逐项操作、核对并存放下载文件。本文与机器清单及实际下载目录已于 2026-09-27 交叉审计。

## 1. 统一目录结构

Method 1–4 都按以下顶层结构管理：

```text
output_compare/
  <母相 CIF>/
    官网/
      Method1/
      Method2/
      Method3/<可读案例目录>/
      Method4/<可读案例目录>/
    现有网页版交互/
      Method1/
      Method2/
      Method3/<可读案例目录>/
      Method4/<可读案例目录>/
```

Method 1/2 一次搜索即可得到整批候选，因此候选目录直接放在 Method 根目录。Method 3/4 需要多次独立运行，且不同运行可能出现同名候选，所以增加一层可读案例目录，防止下载互相覆盖。

## 2. 当前下载状态与顺序

1. **Method 1/2：现在不需要下载。** 两种母相的官网与本地网页输出均已存在；审计结果为 EuAl4 `123+22`、NdNiO2 `82+48`，合计 275 个候选、1100 个本地核心文件，275/275 个位移模式数一致。
2. **Method 3：现在不需要下载。** EuAl4 20/20 张结果表、35/35 个候选通过审计；NdNiO2 20/20 张结果表、42/42 个候选也通过审计。两种母相合计 77 份 CIF、77 份 isoviz、77 份 TOPAS 和 77 份 details 均存在且非空。NdNiO2 `ND-15`～`ND-20` 仍有 12 条“网页包层级非标准”提示，但结果页查询上下文、完整候选表以及 details 中的精确 `(SG,basis,origin,s,i)` 均可无歧义配对，故它们是已核定的 authoritative 证据，不是漏下载或内容错误。
3. **Method 4：现在不要运行或下载。** 目录和 24 个固定 daughter CIF 已准备好，等 Method 3 官网差分完成并收到明确开始指令后再做。

本地差分不会扩大上述下载范围：40 个案例、77 条 embedding 已进入本地算法验证；basis/origin 的字面差异会另用精确仿射群运算判定，coupled-IR 欠枚举属于本地算法整改项，不应通过重复下载已经核准的官网文件来处理。

以后若重新生成 Method 1/2，仍直接替换对应 `现有网页版交互/Method1`、`Method2`，不要改动 `官网/` 参考。

## 3. Method 3 下载规则

每组都选择：

- Types：`strain` + 所有物种的 `displacive`
- Method：`Method 3`
- Sublattice type：`direct`
- Centering：**40 组全部选择 Default**，即单词 `Default` 后面的圆圈（官网表单值为 `d`）
- **不要选择 P**。官网的 P 表示 primitive / no centering；它不是 Default，在 SG 71 Immm 等带心目标空间群上会得到 “There are no subgroups...”
- 输入表中给出的目标空间群和三行 basis；点群保持 `Not selected`
- basis 必须逐项照抄，包括行顺序和负号；即使某组 basis 的体积为 1，也不能自行改成单位矩阵，否则官网可能报 “Basis vectors must be lattice vectors of the parent”

每个案例目录保存：

1. 完整 Method 3 **结果表页面** HTML（不是 Complete modes details 页面）；
2. 使用浏览器“网页，完整”保存时，HTML 同名 `_files` 资源目录可原样保留以便离线打开。候选全集证据仍是案例根目录的结果 HTML；`_files` 不能代替结果页，也不算候选的四类核心文件；
3. 结果表返回的**每个候选**的完整导出，保存 `subgroup.cif`、`data.isoviz`、`topas.str` 和 Complete modes details HTML；只有一个候选时可以直接放在案例根目录；若官网明确返回 0 个候选，只保存该结果页，不创建或伪造候选文件；
4. 多候选案例按官网结果行显示的 SG、basis、origin、s、i 建独立候选子目录（若下一页能无歧义确认 IR/OPD，也可用 `IR OPD` 命名），不得把文件直接解压到案例根目录，否则同名文件会相互覆盖；
5. “预期候选”仅用于快速确认输入没有填错，不能只下载该候选而漏掉同表其他候选。

判断保存正确的结果表：页面正文必须包含 `Finish selecting the distortion mode`，并且满足以下二者之一：含一个或多个 `name="orderparam"` 候选；或明确同时显示 `There are no subgroups...` 与 `Try again`，表示合法的零候选结果。只有 Method 3 输入控件的 `ISODISTORT_ search.html` 不是结果表；每个候选自己的 Complete modes details 也不能代替整张结果表。

多候选目录名中的分数不要再直接删除 `/`，否则 `1/2` 会与整数 `12` 混淆。建议把 `/` 写成 `over`、负号写成 `m`，例如 `1/2 → 1over2`、`-1/4 → m1over4`；目录名只是可读标签，最终身份始终以 `subgroup.cif` 的 `# Subgroup:` 行为准。

推荐结构：

```text
M3-... - SG... - IR OPD/
  method3-result.html
  <多候选时：官网结果行身份>/
    subgroup.cif
    data.isoviz
    topas.str
    ISODISTORT_ complete modes details.html
```

### 3.1 EuAl4 Parent.cif：20组（20/20 已核完，无需下载）

官网保存根目录：`output_compare/EuAl4 Parent.cif/官网/Method3/`

| 案例目录 | 目标 SG | centering | basis：a′；b′；c′ | 官网表行数 / 当前动作 |
| --- | ---: | --- | --- | --- |
| `M3-EU-01 - SG139 - GM1+ P1` | 139 | Default | `[1,0,0]; [0,1,0]; [0,0,1]` | 1，已核完 |
| `M3-EU-02 - SG71 - GM2+ P1` | 71 | Default | `[1,0,0]; [0,1,0]; [0,0,1]` | 1，已核完 |
| `M3-EU-03 - SG107 - GM3- P1` | 107 | Default | `[1,0,0]; [0,1,0]; [0,0,1]` | 1，已核完 |
| `M3-EU-04 - SG12 - GM5+ P1` | 12 | Default | `[0,1,1]; [-1,0,0]; [0,-1,0]` | 3，已核完 |
| `M3-EU-05 - SG8 - GM5- C1` | 8 | Default | `[-1,1,0]; [0,0,1]; [1,0,0]` | 3，已核完 |
| `M3-EU-06 - SG44 - GM5- P1` | 44 | Default | `[0,1,0]; [0,0,1]; [1,0,0]` | 2，已核完 |
| `M3-EU-07 - SG123 - M1+ P1` | 123 | Default | `[1,0,0]; [0,1,0]; [0,0,1]` | 1，已核完 |
| `M3-EU-08 - SG129 - M3- P1` | 129 | Default | `[1,0,0]; [0,1,0]; [0,0,1]` | 1，已核完 |
| `M3-EU-09 - SG47 - X1+ C1` | 47 | Default | `[1,1,0]; [-1,1,0]; [0,0,1]` | 1，已核完 |
| `M3-EU-10 - SG125 - X1- P1` | 125 | Default | `[1,-1,0]; [1,1,0]; [0,0,1]` | 1，已核完 |
| `M3-EU-11 - SG119 - P1 C1` | 119 | Default | `[1,1,0]; [-1,1,0]; [0,0,2]` | 2，已核完 |
| `M3-EU-12 - SG140 - P2 P3` | 140 | Default | `[-1,1,0]; [-1,-1,0]; [0,0,2]` | 2，已核完 |
| `M3-EU-13 - SG2 - N1+ 4D1` | 2 | Default | `[2,0,0]; [0,2,0]; [-1,-1,1]` | 4，已核完 |
| `M3-EU-14 - SG65 - N1+ P1` | 65 | Default | `[2,0,0]; [0,0,-2]; [0,1,0]` | 2，已核完 |
| `M3-EU-15 - SG122 - P5 P1` | 122 | Default | `[-1,1,0]; [-1,-1,0]; [0,0,2]` | 3，已核完 |
| `M3-EU-16 - SG24 - P5 C3` | 24 | Default | `[0,0,-2]; [-1,-1,0]; [-1,1,0]` | 2，已核完 |
| `M3-EU-17 - SG99 - LD1 C1` | 99 | Default | `[1,0,0]; [0,1,0]; [0,0,6]` | 1，已核完 |
| `M3-EU-18 - SG123 - LD1 P1` | 123 | Default | `[1,0,0]; [0,1,0]; [0,0,6]` | 1，已核完 |
| `M3-EU-19 - SG105 - LD2 C1` | 105 | Default | `[1,0,0]; [0,1,0]; [0,0,6]` | 1，已核完 |
| `M3-EU-20 - SG20 - LD5 C4` | 20 | Default | `[-1,1,0]; [-1,-1,0]; [0,0,6]` | 2，已核完 |

EU-15/16/17 的新结果表已分别核实为 3/2/1 行，且与原有候选逐项一致，因此核心文件无需重下。

`M3-EU-07` 与 `M3-EU-18` 都是 SG 123，属于正常且有意保留的两组不同验证条件，并非清单重复：前者验证特殊 k 的 `M1+ P1`（basis 为单位矩阵，官网导出 `s=2, i=2`），后者验证参数 k 的 `LD1 P1`（basis 为 `diag(1,1,6)`，官网导出 `s=12, i=12`）。空间群号只描述所得结构的对称类型，不能唯一确定活动 k 点、IR/OPD、超胞或畸变路径。

### 3.2 NdNiO2 own.cif：20组（20/20 已核完，无需下载）

官网保存根目录：`output_compare/NdNiO2 own.cif/官网/Method3/`

| 案例目录 | 目标 SG | centering | basis：a′；b′；c′ | 官网表行数 / 当前动作 |
| --- | ---: | --- | --- | --- |
| `M3-ND-01 - SG123 - GM1+ P1` | 123 | Default | `[1,0,0]; [0,1,0]; [0,0,1]` | 1，已核完 |
| `M3-ND-02 - SG47 - GM2+ P1` | 47 | Default | `[1,0,0]; [0,1,0]; [0,0,1]` | 1，已核完 |
| `M3-ND-03 - SG99 - GM3- P1` | 99 | Default | `[1,0,0]; [0,1,0]; [0,0,1]` | 1，已核完 |
| `M3-ND-04 - SG115 - GM4- P1` | 115 | Default | `[1,0,0]; [0,1,0]; [0,0,1]` | 1，已核完 |
| `M3-ND-05 - SG2 - GM5+ C1` | 2 | Default | `[0,0,1]; [1,0,0]; [0,1,0]` | 1，已核完 |
| `M3-ND-06 - SG10 - GM5+ P1` | 10 | Default | `[0,0,-1]; [-1,0,0]; [0,1,0]` | 2，已核完 |
| `M3-ND-07 - SG6 - GM5- C1` | 6 | Default | `[0,1,0]; [0,0,1]; [1,0,0]` | 2，已核完 |
| `M3-ND-08 - SG25 - GM5- P1` | 25 | Default | `[0,1,0]; [0,0,1]; [1,0,0]` | 2，已核完 |
| `M3-ND-09 - SG123 - M1+ P1` | 123 | Default | `[1,1,0]; [-1,1,0]; [0,0,1]` | 2，已核完 |
| `M3-ND-10 - SG127 - M2+ P1` | 127 | Default | `[1,1,0]; [-1,1,0]; [0,0,1]` | 2，已核完 |
| `M3-ND-11 - SG129 - M2- P1` | 129 | Default | `[1,1,0]; [-1,1,0]; [0,0,1]` | 2，已核完 |
| `M3-ND-12 - SG139 - A1+ P1` | 139 | Default | `[1,1,0]; [-1,1,0]; [0,0,2]` | 3，已核完 |
| `M3-ND-13 - SG140 - A2+ P1` | 140 | Default | `[1,1,0]; [-1,1,0]; [0,0,2]` | 2，已核完 |
| `M3-ND-14 - SG47 - X1+ C1` | 47 | Default | `[2,0,0]; [0,2,0]; [0,0,1]` | 3，已核完 |
| `M3-ND-15 - SG123 - Z1+ P1` | 123 | Default | `[1,0,0]; [0,1,0]; [0,0,2]` | 2，已核完 |
| `M3-ND-16 - SG11 - Z5+ C1` | 11 | Default | `[0,1,0]; [0,0,2]; [1,0,0]` | 2，已核完 |
| `M3-ND-17 - SG47 - Y1 P1` | 47 | Default | `[-3,0,0]; [0,0,1]; [0,2,0]` | 2，已核完 |
| `M3-ND-18 - SG25 - Y1 C1` | 25 | Default | `[0,2,0]; [0,0,1]; [3,0,0]` | 5，已核完 |
| `M3-ND-19 - SG49 - Y2 P1` | 49 | Default | `[-3,0,0]; [0,0,1]; [0,2,0]` | 2，已核完 |
| `M3-ND-20 - SG51 - Y2 P2` | 51 | Default | `[0,2,0]; [-3,0,0]; [0,0,1]` | 4，已核完 |

`ND-15`～`ND-20` 的结果页与一个 details 网页包仍位于非推荐层级，因此审计报告保留 12 条布局提示；审计器会递归识别唯一结果页，并且只在 details 打印的完整子群身份与候选完全一致时配对。该提示不要求重新下载或移动受保护的 `output_compare/`，也不降低 42 条候选内容的 authoritative 等级。旧 `preflight.candidate_count` 只是单-IR 子集预检数量，不能判断官网下载是否完整。

## 4. Method 4 预建目录与未来下载清单

**当前不要执行或下载 Method 4。** 下列目录只是提前建立，待 Method 3 验证结束后使用。

未来每组保存：

1. 实际上传的固定 `daughter.cif`（必须与机器清单 SHA-256 一致）；
2. 母相 CIF 标识，以及所用 path/子群的 Complete modes details；
3. 完整 Method 4 结果页 HTML；
4. 官网提供的 amplitude table TXT/CSV；若没有下载按钮，保存完整 HTML 表格；
5. 三个预期失败案例保存官网错误页或完整错误文字；
6. Method 4 不产生子群 ZIP，不要为它下载 Method 1–3 的结构包。

官网结果放 `官网/Method4/<案例目录>/`；对应本地结果放 `现有网页版交互/Method4/<案例目录>/`。

### 4.1 EuAl4 Parent.cif：12组

- Γ 上下文：SG 12，`GM5+ P1`
- 参数-k 上下文：SG 99，`LD1 C1`，basis `[1,0,0]; [0,1,0]; [0,0,6]`

| 案例目录 | 上传的 daughter CIF | 目的与预期 |
| --- | --- | --- |
| `G01 - gamma - zero` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G01-zero/daughter.cif` | Γ 零振幅；预期：成功 |
| `G02 - gamma - positive - reordered` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G02-positive-reordered/daughter.cif` | 正振幅、原子倒序；预期：成功 |
| `G03 - gamma - negative - wrapped` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G03-negative-wrapped/daughter.cif` | 负振幅、周期换像；预期：成功 |
| `G04 - gamma - two-mode` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G04-two-mode/daughter.cif` | 两模式混合；预期：成功 |
| `G05 - gamma - origin-shift` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G05-origin-shift/daughter.cif` | 整体原点平移；预期：成功 |
| `G06 - gamma - near-bound-noise` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G06-near-bound-noise/daughter.cif` | 接近振幅上界并加入微噪声；预期：成功且 residual 非零 |
| `P01 - parameter-k - zero-supercell` | `ISODISTORT/output/validation/method4_inputs/EuAl4/P01-zero-supercell/daughter.cif` | LD 参数-k 超胞零振幅；预期：成功 |
| `P02 - parameter-k - positive-supercell` | `ISODISTORT/output/validation/method4_inputs/EuAl4/P02-positive-supercell/daughter.cif` | LD 参数-k 正振幅；预期：成功 |
| `P03 - parameter-k - mixed-reordered-wrapped` | `ISODISTORT/output/validation/method4_inputs/EuAl4/P03-mixed-reordered-wrapped/daughter.cif` | 两模式、倒序、周期换像；预期：成功 |
| `F01 - expected-failure - species-mismatch` | `ISODISTORT/output/validation/method4_inputs/EuAl4/F01-species-mismatch/daughter.cif` | 元素/化学计量不符；应拒绝 |
| `F02 - expected-failure - invalid-lattice` | `ISODISTORT/output/validation/method4_inputs/EuAl4/F02-invalid-lattice/daughter.cif` | 晶格不兼容；应拒绝 |
| `F03 - expected-failure - distance-threshold` | `ISODISTORT/output/validation/method4_inputs/EuAl4/F03-distance-threshold/daughter.cif` | 原子匹配距离超阈值；应拒绝 |

### 4.2 NdNiO2 own.cif：12组

- Γ 上下文：SG 6，`GM5- C1`
- 参数-k 上下文：SG 47，`Y1 P1`，basis `[3,0,0]; [0,2,0]; [0,0,1]`

| 案例目录 | 上传的 daughter CIF | 目的与预期 |
| --- | --- | --- |
| `G01 - gamma - zero` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G01-zero/daughter.cif` | Γ 零振幅；预期：成功 |
| `G02 - gamma - positive - reordered` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G02-positive-reordered/daughter.cif` | 正振幅、原子倒序；预期：成功 |
| `G03 - gamma - negative - wrapped` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G03-negative-wrapped/daughter.cif` | 负振幅、周期换像；预期：成功 |
| `G04 - gamma - two-mode` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G04-two-mode/daughter.cif` | 两模式混合；预期：成功 |
| `G05 - gamma - origin-shift` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G05-origin-shift/daughter.cif` | 整体原点平移；预期：成功 |
| `G06 - gamma - near-bound-noise` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G06-near-bound-noise/daughter.cif` | 接近振幅上界并加入微噪声；预期：成功且 residual 非零 |
| `P01 - parameter-k - zero-supercell` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/P01-zero-supercell/daughter.cif` | Y 参数-k 超胞零振幅；预期：成功 |
| `P02 - parameter-k - positive-supercell` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/P02-positive-supercell/daughter.cif` | Y 参数-k 正振幅；预期：成功 |
| `P03 - parameter-k - mixed-reordered-wrapped` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/P03-mixed-reordered-wrapped/daughter.cif` | 两模式、倒序、周期换像；预期：成功 |
| `F01 - expected-failure - species-mismatch` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/F01-species-mismatch/daughter.cif` | 元素/化学计量不符；应拒绝 |
| `F02 - expected-failure - invalid-lattice` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/F02-invalid-lattice/daughter.cif` | 晶格不兼容；应拒绝 |
| `F03 - expected-failure - distance-threshold` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/F03-distance-threshold/daughter.cif` | 原子匹配距离超阈值；应拒绝 |

另有一项只在本地测试：未选择/计算模式就调用 Method 4 必须被拒绝；该项没有 daughter CIF，也不需要官网下载。
