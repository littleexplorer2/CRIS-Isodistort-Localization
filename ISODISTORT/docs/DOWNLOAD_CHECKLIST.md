# ISODISTORT 官网有限代表输入下载清单（人工动作与目录）

本文件是人工下载动作与目录的唯一入口。2026-10-16 重置后，当前阶段先把**官网三母相已冻结的有限输入集的输出**下载完成，再重建本地网页版刷新清单并重新开始验证与 debug。作废范围、
归档指纹与开放上界见[重置记录](manifests/validation_reset_20261016.json)；范围枚举的证据
与不确定点见 `../output/tmp/official_scope/`（临时研究产物）。机器可读范围见
[official_output_download_manifest.json](manifests/official_output_download_manifest.json)。
完成后的逐例结论写入 [BUGFIX_VALIDATION_REPORT.md](BUGFIX_VALIDATION_REPORT.md)。

## 0. 硬性边界

- `output_compare/` 与 `webpage_info/` 永远只读：已登记的 2026-08…10 官网归档只登记、
  不重下，也不把新下载混入其中（其文件清单指纹已登记，不能失效）。
- 新官网材料写入批次根
  `ISODISTORT/output/validation/official_full_download/20261016-official-full-v1/`。
  该目录是**官网参考只读区**：只由用户的下载动作写入；之后只允许只读审计读取，
  禁止进入预处理、搜索、筛选、模式构造、归一化、拟合、缓存、导出的任何一环。
- 官网页面里的候选全集、子群/setting 答案、模式基、振幅、归一化因子不得用于替代本地
  求解，也不得从结果页反推输入；上传或改名不改变数据来源与角色。
- 有限性：只下载本清单冻结的项。参数化 k 必须先给参数值；Method 3/4 必须先给上界。
- 下载即登记：每项记录查询/候选身份、文件清单、字节数与 SHA-256；官网确实未提供者记
  `not_available`，禁止补造或用别的案例顶替。

## 1. 范围总表（有限可判定集）

| 母相 | Method 1 | Method 2 固定特殊 k | Method 2 参数化 k | Method 3 | Method 4 |
|---|---|---|---|---|---|
| `EuAl4 Parent.cif` (#139) | 123 候选已归档 → 只登记 | 5 标签 × nmod{0,1} = 10 查询（可立即开始） | 10 标签代表值已冻结：LD 已登记，其余 9 新下载 | 待独立输入枚举与上界 | 待独立 daughter 输入与上界 |
| `NdNiO2 own.cif` (#123) | 82 已归档 → 只登记 | 6 标签 × 2 = 12 查询（可立即开始） | 15 标签代表值已冻结：Y 已登记，其余 14 新下载 | 同上 | 同上 |
| `4310_tetra.cif` (#139) | 125 已归档 → 只登记 | 5 标签 × 2 = 10 查询（可立即开始） | 10 标签代表值已冻结，全部新下载（固定 X 旧件单独登记） | 官网目录为空，需新建 | 官网目录为空，需新建 |

固定特殊 k 标签（无自由参数；来源 `ISODISTORT/backend/tables/kpoints_official.py`，
已与 `webpage_info/<母相>/2*.html` 的 k 下拉逐条核对一致）：

- #139（EuAl4、4310）：`GM (0,0,0)`、`X (1/2,1/2,0)`、`M (1,1,1)`、`N (1/2,0,1/2)`、
  `P (1/2,1/2,1/2)`
- #123（NdNiO2）：`GM (0,0,0)`、`M (1/2,1/2,0)`、`X (0,1/2,0)`、`A (1/2,1/2,1/2)`、
  `R (0,1/2,1/2)`、`Z (0,0,1/2)`

## 2. 目录、命名与文件门禁

批次根：`ISODISTORT/output/validation/official_full_download/20261016-official-full-v1/`

```text
<根>/<母相>/官网/Method2/<K>_nmod0/          # K = GM/X/M/N/P/A/R/Z
<根>/<母相>/官网/Method2/<K>_nmod1/
<根>/<母相>/官网/Method2/<查询目录>/<IR>_<OPD>/  # 候选必须位于其查询内，不能跨 k/nmod 覆盖
<根>/<母相>/官网/Method3/<稳定案例号>/
<根>/<母相>/官网/Method4/<稳定案例号>/
<根>/<母相>/官网补页/Method4/               # 仅 EuAl4 G05 auto-origin basis
```

每个查询/案例目录内的必备文件：

| 文件 | 要求 |
|---|---|
| `result_table.html` | Method 2/3 的完整结果表原页；必须含 `Finish selecting the distortion mode` 与候选行，或明确显示 `There are no subgroups` 加 `Try again`。只含输入控件的搜索页不能代替 |
| `<IR>_<OPD>/subgroup.cif` | 结果表中**每一个**候选一份，非空 |
| `<IR>_<OPD>/data.isoviz` | 同上 |
| `<IR>_<OPD>/topas.str` | 同上 |
| `<IR>_<OPD>/complete_modes_details.html` | 同上，含完整模式标签与顺序 |
| `*_files` 资源目录 | 浏览器“网页，完整”生成的资源目录保持原名，不计入候选核心文件 |

每次提交前保留 k 参数、nmod、Types、母相哈希与查询 ID 的独立输入记录及控件截图；
保存的 HTML 默认值不能代替提交态。

命名：候选目录用 `<IR>_<OPD>`（如 `X1+_P1`、`LD1_C1`）；Method 3/4 多候选用
`C01_SG<n>`、`C02_SG<n>`…，按结果行顺序，避免同名覆盖。Method 2 若 IR/OPD 重名，也加稳定行号并登记完整身份。HTML 保留 `.html`/`.htm` 后缀与
原名，审计按内容角色识别，不按文件名。
官网不提供本地网页版的 TXT/CSV 合计导出，这不算官网缺件。

## 3. A 阶段（可立即开始）：Method 2 固定特殊 k

每个查询：打开官网 Method 2 → 选 k 标签（固定点无需参数）→
`# of independent incommensurate modulations` 取 `0`（nmod0）或 `1`（nmod1）→ 提交 →
保存完整结果表 → 对结果表里每个候选保存四件套。

| 母相 | k 标签 | Kovalev | 坐标 | nmod=0 候选数（审计参考） | nmod=1 |
|---|---|---|---|---|---|
| `EuAl4 Parent.cif` | GM | k14 | (0,0,0) | 10 | 未知，下载后登记 |
| `EuAl4 Parent.cif` | X | k13 | (1/2,1/2,0) | 24 | 未知 |
| `EuAl4 Parent.cif` | M | k15 | (1,1,1) | 9 | 未知 |
| `EuAl4 Parent.cif` | N | k11 | (1/2,0,1/2) | 52 | 未知 |
| `EuAl4 Parent.cif` | P | k12 | (1/2,1/2,1/2) | 28 | 未知 |
| `NdNiO2 own.cif` | GM | k17 | (0,0,0) | 11 | 未知 |
| `NdNiO2 own.cif` | M | k18 | (1/2,1/2,0) | 12 | 未知 |
| `NdNiO2 own.cif` | X | k15 | (0,1/2,0) | 18 | 未知 |
| `NdNiO2 own.cif` | A | k20 | (1/2,1/2,1/2) | 11 | 未知 |
| `NdNiO2 own.cif` | R | k16 | (0,1/2,1/2) | 21 | 未知 |
| `NdNiO2 own.cif` | Z | k19 | (0,0,1/2) | 9 | 未知 |
| `4310_tetra.cif` | GM | k14 | (0,0,0) | 11 | 未知 |
| `4310_tetra.cif` | X | k13 | (1/2,1/2,0) | 24 | 未知 |
| `4310_tetra.cif` | M | k15 | (1,1,1) | 10 | 未知 |
| `4310_tetra.cif` | N | k11 | (1/2,0,1/2) | 52 | 未知 |
| `4310_tetra.cif` | P | k12 | (1/2,1/2,1/2) | 28 | 未知 |

计数口径：`output_compare/<母相>/官网/Method1/` 的候选目录按 IR 前缀分区，三母相合计
123/82/125，与目录总数、冻结 manifest、`method_outputs_static_current_20261005.json`
自洽。它是**审计参考**，不是官网 Method 2 的预定答案，也不能决定其它 k 点的期望值。
“Method 1 分区 ⊆ Method 2 候选集”目前只有一个点被独立证实：4310 `X`
（24 个 (IR, 子群类) 身份键相同、`subgroup.cif` 除注释空格逐字节相同、
`!numberOfModulations` 全为 0）；其余标签只是约定，不能当结论。

`nmod` 只取 0 与 1，这是本批输入范围的选择。**没有证明 `nmod≥2` 与 `nmod=1` 等价**；
更高维调制或其它参数值必须另立输入清单，本批不能证明其正确性或完整性。

4310 `X` nmod=0 的四件套与旧归档重复。两种做法都可，需在批次记录里写明选了哪种：

- 字面完整：照下载，与旧归档重复也保留；
- 去重：只存 `result_table.html`，并在审计里登记“候选身份与已归档 Method 1 相同”。

## 4. A 阶段：官网补页 —— EuAl4 Method 4 G05 auto-origin 的 basis/matching 页

> **状态：2026-10-17 已补齐并核对为 `pass`。** 两项证据：
>
> 1. **basis 页**：`output_compare/EuAl4 Parent.cif/官网/Method4/G05/auto-origin/ISODISTORT_ distorted structure (basis).html`
>    （SHA-256 `9e2bfc515d6c13ef7e322041fe16d162976d08936d92dc51193f2e289d020c9a`，9055 B）＋同名
>    `_files/`（5 个文件，与显式那份资源集合相同）。目录、文件名、页面角色（title/h1 均为
>    `ISODISTORT: distorted structure (basis)`）与隐藏字段（`origintype=method4`、`spacegroup=139 I4/mmm`）
>    都正确，且与显式那份**不是同一份**（后者 SHA-256 `ce01fea7…`）；两者只差一个隐藏字段（服务端
>    会话文件名），说明保存下来的 HTML 是这一步表单的**默认态**。
> 2. **提交态截图**：`…/auto-origin/Screenshot 2026-10-17 005445.png`
>    （SHA-256 `230d21826028ca67f4fae55bad27236f6deb00c097a04c4e0c415ffe8abd197c`，290,929 B，1917×1077）
>    证明提交配置：`Specify basis as` 选中且三行系数 `a'=0a+1b+1c`、`b'=1a+0b+0c`、`c'=0a+0b-1c`
>    （即行矢 `(0,1,1)`、`(1,0,0)`、`(0,0,-1)`）；`Automatic origin detection` 选中（`Specify origin as`
>    未选）；`Nearest-site method` 选中（`More robust but slower method` 未选、`d_max=1 Å` 未用）；
>    页头 `Include strain, displacive ALL distortions`；上传的女儿相为 `P1`、10 原子；母相头为
>    `139 I4/mmm`、`Eu1 2a / Al1 4d / Al2 4e z=0.38000`。默认态那页无法证明的部分由此关闭。
>
> **证据替换已登记**：原先被审计的 `Screenshot 2026-09-27 204517.png`（SHA-256 `cc388949…`）
> 已不在 `output_compare` 任何位置；2026-10-06 的审计报告仍引用它，作为其自身快照的历史记录
> **保持不改写**。新旧截图的路径与哈希差异已写入本清单与机器 manifest。
>
> 为它准备的落盘目录 `…/official_full_download/20261016-official-full-v1/EuAl4 Parent.cif/官网补页/`
> 已在核对后按用户要求删除（删除前确认 0 个文件）。以下 4.1–4.6 保留为操作与验收记录。

### 4.1 缺的是什么

`auto-origin/` 目录里已有 `subgroup.cif`、`topas.str`、`data.isoviz`、
`ISODISTORT_ distortion.html`（结果页）、`ISODISTORT_ complete modes details.html` 和填写完成
的截图 `Screenshot 2026-09-27 204517.png`，**缺**该次运行的
`ISODISTORT_ distorted structure (basis).html`（basis/matching 输入页）。同案的显式原点那次
（`G05/`）有这页（SHA-256 `ce01fea7520cce2fe92a9075c54ce0356d51df67f5e49cc481967e27b721b26c`），
不能拿它冒充 auto-origin 那次。两次运行的子群身份相同（见 4.5），缺的是 auto-origin 那次
“超胞基矢 + 原点处理 + 匹配方式”的提交态记录，属**证据完整性**缺口，不是数值缺口。

### 4.2 上传哪个结构（冻结输入，不得改动）

| 项 | 值 |
|---|---|
| 母相 | `experiment_data/EuAl4 Parent.cif`，SHA-256 `2a01942bff7aa74d6156e4907de207ff864b134ccdc85f38f88d53e2d6ebf5e7` |
| 上传的 daughter CIF | `ISODISTORT/output/validation/method4_inputs/EuAl4/G05-origin-shift/daughter.cif` |
| daughter SHA-256 | `a503025301b508690dcf3941f2080d39f1abd7d440ca0cef63a010a875a32dda`（1214 B，Eu2Al8，10 原子） |
| Types | `strain` + 全物种 `displacive`（改动后要按页面提示点 `Change` 才生效） |

不要重新生成 daughter CIF；上传的就是这份哈希固定的文件。

### 4.3 官网操作（两步）

**第一步：Method 4 入口页**

1. 打开官网载入母相 CIF，Types 设为 strain + 全物种 displacive。
2. `Upload distorted structure from CIF file:` 选 4.2 的 `daughter.cif`。
3. `Atom-matching method:` 选 `nearest-site`；`Robust distance threshold` 不启用（留默认）。
4. `Known origin shift (daughter fractional x,y,z; optional):` **留空**。这是 auto-origin 那次
   与显式那次唯一的输入差别：显式那次填 daughter 分数坐标 `0.125, 0.25, 0.0`，即母相分数
   坐标 `(-0.25, -0.125, -0.125)`。
5. 点 `OK` 进入第二步。

**第二步：`ISODISTORT: distorted structure (basis)` 页（就是缺的那一页）**

| 控件（官方字段名） | 应选/应填 | 说明 |
|---|---|---|
| 基矢选择（`inputbasis` 的 `list` + `basisselect`，或 `specify` + `basis11..basis33`） | 选到目标超胞，或选 `specify` 后逐格填 | 目标超胞基矢按行 `(0,1,1)`、`(1,0,0)`、`(0,0,-1)`；`specify` 时九个整数按 `0 1 1 / 1 0 0 / 0 0 -1` 填，保序保号 |
| `Automatic origin detection`（`chooseorigin=true`） | **选它** | auto-origin 那次的定义；不要选 `Specify origin as` |
| `Wyckoff-site matching`（`trynearest=true`） | 选它 | 对应 `nearest-site` |
| `More robust but slower method`（`dmax`） | 不启用，留默认 | 与两份已归档运行一致（robust 未启用） |

**提交前先截屏或打印 PDF**，把上面四项的填写状态拍下来：保存下来的 HTML 只能证明表单结构，
不能证明提交时的勾选状态。

### 4.4 保存到哪、怎么命名

```text
ISODISTORT/output/validation/official_full_download/20261016-official-full-v1/EuAl4 Parent.cif/官网补页/Method4/
```

- 用浏览器“网页，完整”保存第二步页面，保持原名 `ISODISTORT_ distorted structure (basis).html`，
  连同 `ISODISTORT_ distorted structure (basis)_files/` 资源目录一起放入；HTML 按原名引用资源，
  不要改名或只存单文件。
- 官网若同时给出结果页或幅度表，按官网原文件名一并放入本目录。
- 截图/PDF 建议命名 `controls-auto-origin-<YYYYMMDD-HHMMSS>.png`，并注明是提交前状态。

### 4.5 验收（新页面必须满足）

1. 标题（`<title>`/`<h1>`）为 `ISODISTORT: distorted structure (basis)`，审计角色等价于
   `distorted_structure_basis`。
2. 该次分解的子群身份与 `auto-origin/` 已有产物一致：`SG 1 (P1)`、
   basis `{(0,1,1),(1,0,0),(0,0,-1)}`、origin `(-0.25,-0.125,-0.125)`、`s=2`、`i=32`。
3. 幅度表与 `auto-origin/` 已有件一致：30 个 displacive + 6 个 strain 模式；例如
   `[0,0,0]GM2+(a)[Al1:d:dsp]B2(a)` 的 `As=-0.78141 / Ap=-0.55254`、
   `[0,0,0]GM5+(a,b)[Al1:d:dsp]E(a)` 的 `As=0.52824 / Ap=0.37352`。
4. 新页面不得覆盖 `G05/` 里的显式原点那份；若与它逐字节相同，必须在审计记录里写明
   “与显式原点页相同”，不能当成独立证据。

**已知的旧证据矛盾（写进审计，不要掩盖）**：官网审计把本案例记为
`automatic_origin_checked_in_saved_html: true`，但已归档的显式原点那次 basis 页里
`chooseorigin=false`（`Specify origin`，`origin1..3 = 0,0,0`）。所以提交态只能靠本次截图/PDF
确认；若本次复现出的是显式那次的结果而不是 `auto-origin/` 的产物，保留 `pass_with_warnings`
并在审计里记录差异，不得改写已有数值证据。

### 4.6 取不到时

按第 8 节写 `not_available`：母相、Method 4、案例 `G05-origin-shift`、上传 daughter 的哈希、
四项控件设置、尝试时间、官网原始响应（错误文本或页面截图）、结论；保留原有
`pass_with_warnings`，不改动已接受的数值证据。

## 5. B 阶段：Method 3（尚未冻结，暂不下载）

旧闭包脚本已撤销，运行 `tests/manual/measure_method3_closure.py` 将明确退出（代码 2）。
它从官网 Method 1 输出提取基矢/子群再传回本地搜索，违反输入隔离；“40 例分类一致”、
114 种目标群、HNF 代表类数及据此估算的规模均不能作为范围完备性的证据。
旧 checkpoint 无论是否签名都不得恢复，源代码与报告原文保存在历史归档中。

新的有限范围须从独立输入定义：母相与目标 setting、原胞平移、允许的平移指数和基矢界限。
令 Pp/Pc 为母/子相原胞行基矢，B 为子惯用胞相对母惯用胞的行变换：
先验证 `Pc B Pp^-1` 为整数矩阵，再以 `abs(det(Pc B))/abs(det(Pp))` 计算平移指数。
单看 `det(B)` 或伴随矩阵奇偶会漏掉心化；不能预先断言 I 母相只允许 P/I 子胞。
“二阶平移特征”也不等价于所有晶系的特殊 k 点，四类历史标签不能作穷尽性证明。

独立枚举、等价关系和有限覆盖证明完成后，才能冻结具体查询并恢复该阶段。
这不阻止当前已核对的 Method 2 查询下载。

## 6. Method 4（尚未冻结，暂不下载新案例）

先独立核定路径、模式身份、daughter CIF 哈希、基矢/原点/matching 参数，以及每路径的
模式数、幅度、应变、噪声和样本数上界。旧闭包不能给出这一范围，也不保留其推算的总数。
已有 G05 补页及提交态截图属于独立历史证据，状态见第 4 节，不需要重下。

冻结后的每例保存上传参数记录、完整结果页与官网实际提供的四件套；预期拒绝须保留原始
错误页。旧案例只是已有覆盖，不代表新范围完备。

## 7. 参数化 k：代表值已冻结

统一代表值仍为 `a=1/6`、`b=1/4`、`g=1/3`。已归档的 EuAl4 `LD g=1/6` 与
NdNiO2 `Y a=1/3` 保留原输入，只登记不重下。**参数互不相等不是一般位置的充分条件。**

2026-10-17 精确复核通过：`tests/manual/verify_parametric_k_grid.py` 直接读取机器清单中的
输入参数，采用倒空间作用 `R^-T` 和原胞平移判据 `P delta ∈ Z³`，并用仿射模板的符号
系数计算一般位置小余群，逐项检查模倒格平移和同 k-star 碰撞。没有使用官网计算结果。

当前报告：`output/validation/reciprocal_scope_20261017/grid-report.json`：51 个标签
（16 固定、35 参数化）、0 issue、无模倒格/同星碰撞。旧报告明确失效：I 心必须满足
`h+k+l` 为偶数，故 GM 与 M=(1,1,1) **不等价**；EuAl4/4310 的 Q 小余群阶数为 2，W 为 4。
代表值无需改动。这里的“一般位置”仅指给定仿射模板的小余群；有理代表值不能证明
非公度物理模式、整个布里渊区覆盖或全部参数域正确。推导见资料库补充证明
`sec:grid-duality`，报告与源码签名登记于机器清单。

| 母相 | 标签 | 模板 | 代表值 | 坐标 | 动作 |
|---|---|---|---|---|---|
| EuAl4 | LD k10 | (0,0,g) | g=1/6 | (0 , 0 , 1/6) | 已归档，只登记 |
| NdNiO2 | Y k6 | (a,1/2,0) | a=1/3 | (1/3 , 1/2 , 0) | 已归档，只登记 |
| EuAl4 | A k3 | (a,a,g) | a=1/6, g=1/3 | (1/6 , 1/6 , 1/3) | 新查询 ×2（nmod 0/1） |
| EuAl4 | B k1 | (a,0,g) | a=1/6, g=1/3 | (1/6 , 0 , 1/3) | 新查询 ×2（nmod 0/1） |
| EuAl4 | C k2 | (a,b,0) | a=1/6, b=1/4 | (1/6 , 1/4 , 0) | 新查询 ×2（nmod 0/1） |
| EuAl4 | DT k7 | (a,a,0) | a=1/6 | (1/6 , 1/6 , 0) | 新查询 ×2（nmod 0/1） |
| EuAl4 | GP k0 | (a,b,g) | a=1/6, b=1/4, g=1/3 | (1/6 , 1/4 , 1/3) | 新查询 ×2（nmod 0/1） |
| EuAl4 | Q k5 | (1/2,b,1/2) | b=1/4 | (1/2 , 1/4 , 1/2) | 新查询 ×2（nmod 0/1） |
| EuAl4 | SM k6 | (a,0,0) | a=1/6 | (1/6 , 0 , 0) | 新查询 ×2（nmod 0/1） |
| EuAl4 | W k9 | (1/2,1/2,g) | g=1/3 | (1/2 , 1/2 , 1/3) | 新查询 ×2（nmod 0/1） |
| EuAl4 | Y k8 | (a,-a+1,0) | a=1/6 | (1/6 , 5/6 , 0) | 新查询 ×2（nmod 0/1） |
| 4310 | A k3 | (a,a,g) | a=1/6, g=1/3 | (1/6 , 1/6 , 1/3) | 新查询 ×2（nmod 0/1） |
| 4310 | B k1 | (a,0,g) | a=1/6, g=1/3 | (1/6 , 0 , 1/3) | 新查询 ×2（nmod 0/1） |
| 4310 | C k2 | (a,b,0) | a=1/6, b=1/4 | (1/6 , 1/4 , 0) | 新查询 ×2（nmod 0/1） |
| 4310 | DT k7 | (a,a,0) | a=1/6 | (1/6 , 1/6 , 0) | 新查询 ×2（nmod 0/1） |
| 4310 | GP k0 | (a,b,g) | a=1/6, b=1/4, g=1/3 | (1/6 , 1/4 , 1/3) | 新查询 ×2（nmod 0/1） |
| 4310 | LD k10 | (0,0,g) | g=1/3 | (0 , 0 , 1/3) | 新查询 ×2（nmod 0/1） |
| 4310 | Q k5 | (1/2,b,1/2) | b=1/4 | (1/2 , 1/4 , 1/2) | 新查询 ×2（nmod 0/1） |
| 4310 | SM k6 | (a,0,0) | a=1/6 | (1/6 , 0 , 0) | 新查询 ×2（nmod 0/1） |
| 4310 | W k9 | (1/2,1/2,g) | g=1/3 | (1/2 , 1/2 , 1/3) | 新查询 ×2（nmod 0/1） |
| 4310 | Y k8 | (a,-a+1,0) | a=1/6 | (1/6 , 5/6 , 0) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | B k3 | (0,b,g) | b=1/4, g=1/3 | (0 , 1/4 , 1/3) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | C k5 | (a,a,g) | a=1/6, g=1/3 | (1/6 , 1/6 , 1/3) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | D k1 | (a,b,0) | a=1/6, b=1/4 | (1/6 , 1/4 , 0) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | DT k8 | (0,b,0) | b=1/4 | (0 , 1/4 , 0) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | E k2 | (a,b,1/2) | a=1/6, b=1/4 | (1/6 , 1/4 , 1/2) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | F k4 | (a,1/2,g) | a=1/6, g=1/3 | (1/6 , 1/2 , 1/3) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | GP k0 | (a,b,g) | a=1/6, b=1/4, g=1/3 | (1/6 , 1/4 , 1/3) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | LD k13 | (0,0,g) | g=1/3 | (0 , 0 , 1/3) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | S k11 | (a,a,1/2) | a=1/6 | (1/6 , 1/6 , 1/2) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | SM k10 | (a,a,0) | a=1/6 | (1/6 , 1/6 , 0) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | T k7 | (a,1/2,1/2) | a=1/6 | (1/6 , 1/2 , 1/2) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | U k9 | (0,b,1/2) | b=1/4 | (0 , 1/4 , 1/2) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | V k14 | (1/2,1/2,g) | g=1/3 | (1/2 , 1/2 , 1/3) | 新查询 ×2（nmod 0/1） |
| NdNiO2 | W k12 | (0,1/2,g) | g=1/3 | (0 , 1/2 , 1/3) | 新查询 ×2（nmod 0/1） |

计数：参数化标签 35 个（其中 33 个是新查询），
查询 70 个（新查询 66 个）。
连同固定特殊 k 的 32 个查询，Method 2 本轮共需下载 98 个查询。

98 指查询次数，不是候选数或文件数；4310 X nmod0 可按第 3 节登记旧件。
本地运行规模预估属于可选性能工作，不是下载前置条件。保存每个查询实际返回的全部候选，
不能用 Method 1 的历史计数截断结果。参数目录名称见每条查询的 `download_directories`（相对批次根）；98 个既有目录均已核对。

## 8. 未提供项与失败记录

每个 `not_available` 记录：母相、Method、查询/案例身份、参数、尝试时间、官网原始响应
（错误页或文本）、结论与后续动作。禁止用其它案例的输出顶替，禁止改名冒充。

## 9. 下载完成后的只读审计

1. **Method 2**（新根目录）：目前没有可复用的维护版审计工具，需要新增一个只读审计：
   解析 `result_table.html` 得到候选全集，核对每个候选的四件套非空、身份
   `(SG, basis, origin, s, i)`、`!numberOfModulations`，以及结果表候选数与目录数一致。
2. **Method 3**：`tests/manual/audit_method3_downloads.py <root> --json-output <path>`
   可直接指向新根目录；所属 manifest 由根路径或案例名推断，必要时加 `--parent-cif`。
3. **Method 4**：`tests/manual/audit_method4_official.py` 的官网根目前硬编码
   `output_compare`，需要新增 `--official-root`（后续实现项）才能审计新根目录。
4. **Method 1**：只登记；`output/validation/audit_method1_official_downloads_20261003.json`
   的 330/330 结论继续作为旧归档记录，不用于新范围验收。

从 CRIS 根目录通过 `powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_cris.ps1`
执行项目命令；不把旧 DS 的 pwsh/沙箱错误当成当前下载阻断。需要 WSL 的测试须在允许 WSL
的环境执行，缺失依赖时明确记录未运行。

## 10. 之后：本地清单重建与重跑

1. 审计通过后，从预先登记的独立用户输入重建 `local_web_output_refresh_manifest.json`
   （新 batch id）与 runtime revision；不得将官网候选集、基矢或模式回填本地输入。
   旧清单保持 `superseded`；本地独立求解并冻结后，再比较完整候选集合。
2. `download_local_web_outputs.mjs` 只有在 `status == "frozen"` 且
   `parent_manifest_sha256` 与新 revision 一致时才会运行；这是恢复本地刷新的前置条件。
3. 本地产物写入
   `ISODISTORT/output/validation/web_output_refresh/<新批次>/<母相>/现有网页版/<Method>/`；
   官网材料与本地产物分树存放，审计只读官网树。
4. 用新官网集重跑 Method 1 → 2 → 3 → 4 对照；旧报告的数字与结论不得作为新范围验收。
