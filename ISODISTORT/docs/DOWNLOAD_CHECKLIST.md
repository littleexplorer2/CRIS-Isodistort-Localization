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
3. **Method 4：现在不需要下载。** 本地重冻结矩阵与官网均为 24/24 个有效输入；官网审计为 23 pass + EuAl4 G05 auto-origin 1 warning + 0 fail。NdNiO2 F01 已按物种不符拒绝、F02 均匀应变成功、F03 按 robust `dmax=0.1 Å` 匹配失败，三者均通过。G05 warning 仅因 auto-origin 运行缺 basis HTML；已有填写截图与完整数值证据，不要求重跑。

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

## 4. Method 4 官网下载清单（已完成，供后续复核）

只上传 `docs/manifests/method4_download_manifest.json` 中重冻结签名 `23c951b9…b279` 所列文件；上传前核对逐例 SHA-256。不要重新生成、手改或另找相似 CIF。

执行顺序：

1. EuAl4 `G01`–`G06` 已完成核心审计；实际上传副本均与冻结 SHA-256 一致。G01–G04、G06 的 basis/matching HTML 与填写完成截图齐全并完全通过；G05 显式-origin 证据齐全，仅 `auto-origin/` 没有 basis HTML，凭有效截图与 accepted result identity 维持 `pass_with_warnings`。旧数值文件均无需重跑或重下。
2. `P01 - parameter-k - zero-supercell` 已完全通过审计：官网为 `1 P1`、basis `{(1,0,0),(0,1,0),(0,0,6)}`、origin `(0,0,0)`、`s=12,i=192`，180 个位移模式和 6 个应变模式全零；basis HTML、填写截图与完整结果均可正常解析。
3. P02/P03 旧输入的最大原子位移为 `4.68846 Å` / `3.90492 Å`，官网结果页明确为 `Nearest-site method failed.`；两个错误页已保存在 `large-displacement-nearest-failed/`。新冻结输入把实际最大位移降到 `0.30 Å` / `0.35 Å`，同时保留参数-k、非零双分量、倒序和周期换像覆盖，均已用 automatic origin + nearest-site 成功分解并通过审计。P03 的二维 LD1 分量在官网与本地等价基之间发生旋转，完整子空间的 `As/Ap` 欧氏范数在官网显示精度内一致。
4. EuAl4 F01/F03 已通过预期拒绝审计。F02 官网成功输出 30 个全零位移模式与 6 个应变模式，其中 4 个非零；本地已按官网 `M=I+epsilon` 约定恢复六个应用 Voigt 分量，最大官网差 `5.08e-6`，小于由 `data.isoviz` 五位幅值/模式向量推得的逐分量舍入界。旧 F02 文件无需重下。
5. NdNiO2 G01–G06 已完全通过。官网把氧 `f` 位的同一四维 Γ 子空间写成 `B3u(a,b)+B2u(a,b)`，本地连续列名为 `a,b,c,d`；审计在确认完整维数和相同归一化后比较基不变 `As/Ap` 范数，不再按不可比的单列字符串强配对。六例均得到 `1 P1`、相同冻结 basis、`s=1,i=16`；G05 的显式与自动原点运行一致。G06 的带噪共享幅度使用归档结构推导的严格误差界通过，不需要重跑或重下。
6. NdNiO2 P01–P03 已完全通过：三例官网均为 `1 P1`、basis `{(-3,0,0),(0,0,1),(0,2,0)}`、`s=6,i=96`、72 个位移模式与 6 个应变模式；P01 全零，P02/P03 的参数-k分支在自动原点相位旋转后与本地归一化范数一致。三例导出结构最大距离不超过 `4.04e-5 Å`，不需要重跑或重下。
7. NdNiO2 F01/F02 已通过：F01 错误页明确显示母女相原子类型不匹配；F02 得到 `1 P1`、basis `{(0,1,0),(0,0,1),(1,0,0)}`、`s=1,i=16`，导出结构原子/物种完全一致，晶格与应变差分均在舍入误差界内。
8. NdNiO2 F03 首次运行曾误上传 F02 daughter，官网成功输出因此只是 F02 的重复运行，不是 robust 阈值失效；该错误证据未计入正式矩阵。
9. NdNiO2 F03 已用正确冻结输入重跑并通过。归档与冻结 SHA-256 均为 `02b41e63c770b9c4e5ff0a622429841be5e74478f1ca478c9755bff459332e36`；截图显示 `a=3.92080, b=3.28100, c=3.92080 Å`、Nd `(0.86,0.81,0.50)`、basis `{(0,1,0),(0,0,1),(1,0,0)}`、automatic origin、robust matching、`dmax=0.1 Å` 且 manual mapping 未选。官网按预期返回 `Failed to find match. Try using a larger value for dmax.`，无成功导出；当前无人工下载动作。
10. 官网 Method 4 的直接输入是母相、daughter CIF、basis、origin 和 matching；通常不会先让用户点选 `IR/OPD/SG/s/i`。下方及清单中的 `resolved_contexts.*.candidate_identity` 是本地模式基的审计身份：官网若给出多个近似晶格 basis，选择与其中 `basis` 完全相同的一项，但不要在不存在相应控件时强行寻找或填写 IR/OPD。`resolved_contexts.*.mode_labels` 是本地实际拟合矩阵的完整列顺序；官网 Complete modes details 的完整标签/顺序也必须保存，不能只抄案例中非零的前一两个模式。
11. 冻结的 daughter CIF 使用 P1 完整原子表保存，以保留行重排、周期换像和微噪声，不是待上传的“子群导出文件”。不要修改其 P1 标记、合并原子或另存成高对称 ASU；官网可能因此列出比本地受限子群模式基更多的零幅度/噪声模式，这正是差分需要记录的内容。EuAl4 G01/G02 实际为 P1、`s=2,i=32`；NdNiO2 G01–G06 实际为 P1、`s=1,i=16`，P01–P03 为 P1、`s=6,i=96`。它们与本地受限子群模式上下文不同是预期现象。
12. `G05` 的 CIF 坐标是用女儿胞分数平移 `+(1/8,1/4,0)` 生成的；官网字段要求“新超胞原点相对母相原点”，因此必须使用其**逆平移**并换算到母胞：EuAl4 为 `(-1/4,-1/8,-1/8)`，NdNiO2 为 `(0,-1/8,-1/4)`。两者冻结在 `official_origin_shift_parent_fractional`。G05 案例根目录保存显式 origin 运行的六个标准文件及上传副本；`auto-origin/` 保存自动枚举运行的六个同名标准文件和截图/PDF。禁止互相覆盖。
13. `F03` 选择 robust matching、阈值 `0.1 Å`；其余案例先用 nearest-site。若官网只能显示默认单位/阈值，原样记录界面值，不自行换算。

每组保存：

1. 实际上传的固定 `daughter.cif` 副本（必须与机器清单 SHA-256 一致，建议命名 `uploaded-daughter.cif`）；
2. basis/origin/matching 页面 HTML，另加一张**填写完成后**的截图或打印 PDF；仅“保存网页”不能保留 JavaScript/动态输入状态。若已通过案例只剩清晰截图与 accepted result identity，不必重跑数值文件，但会保留 HTML 缺失的证据警告；
3. 完整 Method 4 结果页和 Complete modes details HTML；
4. 官网逐项提供的 `subgroup.cif`、`topas.str`、`data.isoviz`，以及 amplitude TXT/CSV（若有）；
5. 预期失败案例保存官网错误页或完整错误文字，并保留填写完成后的 matching/阈值截图；F02 已由官网实证为 strain 成功例，必须保存完整成功产物；
6. Method 4 不产生 Method 1–3 那种整批子群 ZIP；不要用其他 Method 的 ZIP 代替上述逐项文件。

官网结果放 `官网/Method4/<案例目录>/`；对应本地结果放 `现有网页版交互/Method4/<案例目录>/`。

### 4.1 EuAl4 Parent.cif：12组

- Γ 上下文：`GM5+ P1 (a,0)`，SG `12 C2/m`，basis `{(0,1,1),(1,0,0),(0,0,-1)}`，origin `(0,0,0)`，`s=1, i=4`
- 参数-k 上下文：`LD1 C1 (a,b)`，SG `99 P4mm`，basis `{(1,0,0),(0,1,0),(0,0,6)}`，origin `(0,0,0)`，`s=12, i=24`，k-active `(0,0,1/6)`
- `G05` 显式 origin：官网母胞分数坐标 `(-1/4,-1/8,-1/8)`；这是生成 CIF 时女儿胞正平移 `(1/8,1/4,0)` 经 basis 换算后的逆平移，符号不可混用

| 案例目录 | 上传的 daughter CIF | 目的与预期 |
| --- | --- | --- |
| `G01 - gamma - zero` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G01-zero/daughter.cif` | Γ 零振幅；核心、basis HTML 与截图完全通过 |
| `G02 - gamma - positive - reordered` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G02-positive-reordered/daughter.cif` | 正振幅、原子倒序；核心、basis HTML 与截图完全通过 |
| `G03 - gamma - negative - wrapped` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G03-negative-wrapped/daughter.cif` | 负振幅、周期换像；核心、basis HTML 与截图完全通过 |
| `G04 - gamma - two-mode` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G04-two-mode/daughter.cif` | 两模式混合；核心、basis HTML 与截图完全通过 |
| `G05 - gamma - origin-shift` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G05-origin-shift/daughter.cif` | 显式与 auto-origin 均通过；仅 auto-origin 缺 basis HTML，保留证据警告 |
| `G06 - gamma - near-bound-noise` | `ISODISTORT/output/validation/method4_inputs/EuAl4/G06-near-bound-noise/daughter.cif` | 主模式、P1 噪声分配、basis HTML 与截图完全通过；不要求 residual 相同 |
| `P01 - parameter-k - zero-supercell` | `ISODISTORT/output/validation/method4_inputs/EuAl4/P01-zero-supercell/daughter.cif` | 180 个位移模式零保持、basis HTML 与完整导出完全通过 |
| `P02 - parameter-k - positive-supercell` | `ISODISTORT/output/validation/method4_inputs/EuAl4/P02-positive-supercell/daughter.cif` | `dmax=0.30 Å`，SHA-256 `330d0250…b73b`；非零参数-k、180 模式与完整导出通过 |
| `P03 - parameter-k - mixed-reordered-wrapped` | `ISODISTORT/output/validation/method4_inputs/EuAl4/P03-mixed-reordered-wrapped/daughter.cif` | `dmax=0.35 Å`，SHA-256 `68eb78ec…41e`；二维 LD1 子空间范数、倒序/周期换像与完整导出通过 |
| `F01 - expected-failure - species-mismatch` | `ISODISTORT/output/validation/method4_inputs/EuAl4/F01-species-mismatch/daughter.cif` | 官网已按元素/化学计量不符拒绝；审计通过 |
| `F02 - homogeneous-strain（目录保留旧名）` | `ISODISTORT/output/validation/method4_inputs/EuAl4/F02-invalid-lattice/daughter.cif` | 官网与本地应用应变张量差分通过；旧目录名仅为保持归档稳定，不重下 |
| `F03 - expected-failure - distance-threshold` | `ISODISTORT/output/validation/method4_inputs/EuAl4/F03-distance-threshold/daughter.cif` | 官网已按 robust `dmax=0.1 Å` 匹配失败拒绝；审计通过 |

### 4.2 NdNiO2 own.cif：12组

- Γ 上下文：`GM5- C1 (a,b)`，SG `6 Pm`，basis `{(0,1,0),(0,0,1),(1,0,0)}`，origin `(0,0,0)`，`s=1, i=8`
- 参数-k 上下文：`Y1 P1 (a,0;0,0)`，SG `47 Pmmm`，basis `{(-3,0,0),(0,0,1),(0,2,0)}`，origin `(0,0,0)`，`s=6, i=12`，k-active `(1/3,1/2,0)`
- `G05` 显式 origin：官网母胞分数坐标 `(0,-1/8,-1/4)`；这是生成 CIF 时女儿胞正平移 `(1/8,1/4,0)` 经 basis 换算后的逆平移，符号不可混用
- F03 正确归档 SHA-256 为 `02b41e63c770b9c4e5ff0a622429841be5e74478f1ca478c9755bff459332e36`；已通过审计，无需再次上传。

| 案例目录 | 上传的 daughter CIF | 目的与预期 |
| --- | --- | --- |
| `G01 - gamma - zero` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G01-zero/daughter.cif` | Γ 零振幅；官网审计通过 |
| `G02 - gamma - positive - reordered` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G02-positive-reordered/daughter.cif` | 正振幅、原子倒序；官网审计通过 |
| `G03 - gamma - negative - wrapped` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G03-negative-wrapped/daughter.cif` | 负振幅、周期换像；官网审计通过 |
| `G04 - gamma - two-mode` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G04-two-mode/daughter.cif` | 两模式混合；官网审计通过 |
| `G05 - gamma - origin-shift` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G05-origin-shift/daughter.cif` | 显式与自动原点一致；官网审计通过 |
| `G06 - gamma - near-bound-noise` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/G06-near-bound-noise/daughter.cif` | 主模式与额外 P1 噪声模式按结构推导误差界通过；不要求 residual 相同 |
| `P01 - parameter-k - zero-supercell` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/P01-zero-supercell/daughter.cif` | Y 参数-k 超胞零振幅；官网审计通过 |
| `P02 - parameter-k - positive-supercell` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/P02-positive-supercell/daughter.cif` | Y 参数-k 正振幅与 origin phase；官网审计通过 |
| `P03 - parameter-k - mixed-reordered-wrapped` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/P03-mixed-reordered-wrapped/daughter.cif` | 两模式、倒序、周期换像与半周期反号；官网审计通过 |
| `F01 - expected-failure - species-mismatch` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/F01-species-mismatch/daughter.cif` | 官网已按元素/化学计量不符拒绝；审计通过，不重跑 |
| `F02 - homogeneous-strain（目录保留旧名）` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/F02-invalid-lattice/daughter.cif` | 8% 晶格变化的均匀应变正例；官网导出与差分审计通过，不重跑 |
| `F03 - expected-failure - distance-threshold` | `ISODISTORT/output/validation/method4_inputs/NdNiO2/F03-distance-threshold/daughter.cif` | 官网已按 robust `dmax=0.1 Å` 匹配失败拒绝；哈希、截图和错误语义通过 |

另有一项只在本地测试：未选择/计算模式就调用 Method 4 必须被拒绝；该项没有 daughter CIF，也不需要官网下载。
