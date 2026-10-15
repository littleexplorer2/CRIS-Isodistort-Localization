# ISODISTORT 下载与上传清单

本文件列当前本地网页版产物重下和仍可能需要用户执行的官网补页操作。已完成批次的逐例
结果、验证结论和机器报告索引写入 [BUGFIX_VALIDATION_REPORT.md](BUGFIX_VALIDATION_REPORT.md)，
精确机器输入写入 `docs/manifests/`；这里不重复维护完成案例表。

## 当前动作

**下一优先动作是把当前已有权威查询的本地网页版输出全部重跑并保存；官网不重跑，只补
已经定义但确实缺失的页面。**

### 执行入口

使用 [统一冻结清单](manifests/local_web_output_refresh_manifest.json) 中的 74 项
（Method 1/2/3/4：3/7/40/24）。当前下载尚未启动；从 `M1-EU-ALL` 开始，
启动前核对来源文件、输入哈希和实际运行签名。冻结历史见 [归档索引](history/README.md)。

### “全部”的有限范围

- 母相只包括已有正式网页/对照基线的 `EuAl4 Parent.cif`、`NdNiO2 own.cif` 和
  `4310_tetra.cif`。`EuAl4 Springer.cif` 是辅助输入，`LD1_C1_subgroup.cif` 是子群结构，
  都不是第四、第五种母相。
- Method 1：三种母相按各自已核定的标准 Types 查询，保存整张候选表和全部候选输出。
- Method 2：按统一冻结清单保存 EuAl4 LD `g=1/6`、NdNiO2 Y `a=1/3`、两者各一个
  X 特殊-k查询（EuAl4 `(1/2,1/2,0)`、NdNiO2 `(0,1/2,0)`），以及 4310 已归档的
  X `(1/2,1/2,0)` 查询。参数-k四格式 ZIP 分别保存 `nmod=0` 与
  `nmod=1`；`nmod=2/3` 与 1 语义相同，不重复下载。尚未建立的 4310 Gamma/参数-k扩展
  不是当前可下载集合，须先按下文冻结 manifest。
- Method 3：只使用 `method3_download_manifest.json` 已冻结的 EuAl4/NdNiO2 共 40 个查询。
- Method 4：只使用 `method4_download_manifest.json` 已冻结的 EuAl4/NdNiO2 共 24 个
  daughter 案例。

“全部”指上述每个查询的全部结果行，不指穷举 Method 1 的筛选组合、所有 Method 2 k 点、
所有 Method 3 空间群/超胞或任意 daughter CIF。没有冻结输入的无限参数空间不能用几个
临时案例冒充全量。

### 本地网页版重下内容

- Method 1–3：清空筛选条件后分别下载完整 TXT、CSV 和勾选 CIF、IsoVIZ、Complete modes、
  TOPAS 的原始 ZIP；保存 ZIP 根部 `export_candidate_status.json/.txt`，并核对全部候选的
  四种文件、大小和 SHA-256。筛选只改变同一结果集的视图，不逐一下载所有筛选组合。
- Method 4：24 个冻结案例逐例运行；成功案例下载完整 TXT/CSV，失败案例保存原始网页/API
  错误与输入哈希。Method 4 没有 ZIP，不补造四格式文件。
- 每次运行记录母相和输入哈希、源码/依赖/运行时签名、Types、Method 参数、Method 2
  `nmod`、完整查询/候选身份、开始/结束时间和状态。旧机器报告只作预检，不能代替本轮
  浏览器产物。

### 官网当前补页

仅尝试补 EuAl4 Method 4 G05 auto-origin 的 basis HTML，写入新补页目录。
无法取得时记录 not_available 并保留 pass_with_warnings；已有官网核心文件不重下载。
官网未提供的本地 TXT/CSV 导出格式不算官网缺件。
4310 尚未定义的扩展查询须先冻结，再作为新增批次执行。

本地重跑和官网补页分别写入：

```text
ISODISTORT/output/validation/web_output_refresh/<YYYYMMDD-HHMMSS>/<母相>/现有网页版/<Method>/
ISODISTORT/output/validation/web_output_refresh/<YYYYMMDD-HHMMSS>/<母相>/官网补页/<Method>/
```

`output_compare/` 永远只读；审计完成后是否由用户另行归档，不在下载脚本中自动决定。

## 4310 Method 2/3/4 的触发条件

### Method 2

Codex 先冻结覆盖 Gamma、至少一个非 Gamma 特殊 k 和一条公度参数 k 的正式 manifest，
并分别列出 `nmod=0` 与 `nmod>=1`、完整 Types、KVALUE、IR/OPD、预期候选身份和四格式
门禁。现有 X 点 24 个候选只作为官网已归档项登记，不进入官网重下载清单，也不能决定
其它 k 点的预期答案。

现有 X 点查询在本轮本地重跑中重新下载 TXT/CSV 与四格式 ZIP；其官网24个候选只登记已有
文件与哈希，不重下。扩展 manifest 通过自检后，本地运行全部新增查询，官网只执行这些
此前未定义的新查询。若参数 k 子群生成或本地身份门禁失败，保留进度、状态报告和原始
错误为 `blocked`，修复通用算法后续跑，不得删减 manifest 来获得“全通过”。

### Method 3

Codex 必须先生成并自检机器 manifest，至少覆盖 identity、oriented-cell、
special-supercell 和 parameter-k。manifest 要冻结母相 SHA-256、查询参数、稳定案例号、
目标目录和期望候选身份。清单未冻结前不要在官网操作，也不要凭名称手写案例。

manifest 完成后，每个案例按其中的值选择：

- Types：`strain` 与所有物种的 `displacive`；
- Method：`Method 3`；
- Sublattice type：`direct`；
- Centering：`Default`，除非 manifest 明确给出其他值；
- 目标空间群与三行 basis：逐项照抄，保留顺序与负号。

每个案例保存完整结果表 HTML。结果表必须包含
`Finish selecting the distortion mode`，并包含一个或多个 `name="orderparam"` 候选，
或同时明确显示 `There are no subgroups...` 与 `Try again`。只含输入控件的搜索页和
单个候选的 Complete modes details 都不能代替结果表。

结果表中的每个候选都保存：

1. `subgroup.cif`；
2. `data.isoviz`；
3. `topas.str`；
4. Complete modes details HTML。

多候选按结果行顺序放入 `C01_SG<number>`、`C02_SG<number>` 等短目录，避免同名文件
覆盖。下载后先核对结果表候选数、完整 `(SG,basis,origin,s,i)`、文件类型、非空性和
重复身份；门禁通过前不开始算法差分。

目标根目录为：

```text
ISODISTORT/output/validation/web_output_refresh/<批次>/4310_tetra.cif/官网新增/Method3/<稳定案例号>/
```

### Method 4

4310 的 Method 1、2、3 代表路径通过后，Codex 才能从这些路径生成 daughter CIF，并在
`ISODISTORT/output/validation/method4_inputs/4310/` 冻结输入哈希、basis、origin、
matching、候选身份和预期结果。覆盖范围至少包括：

- 零、正/负单模式和正交双模式；
- 参数-k、多模式及接近允许上界的案例；
- 原子重排、周期换像、显式/自动原点和小噪声；
- 均匀应变；
- 物种不符与距离阈值失败路径。

manifest 冻结后逐项上传完全相同的 daughter CIF。成功案例保存填写完成后的
basis/origin/matching 截图或 PDF、结果页，以及官网实际提供的 Complete modes details、
CIF、TOPAS、IsoVIZ；若官网另提供幅度表才保存该原始文件。预期失败案例保存错误页和完整填写状态。
下载后先核对上传文件 SHA-256、完整身份、模式数和产物集合，再做数值差分。

目标根目录为：

```text
ISODISTORT/output/validation/web_output_refresh/<批次>/4310_tetra.cif/官网新增/Method4/<稳定案例号>/
```

## 文件与目录约定

以 manifest 的稳定案例 ID 隔离查询；短目录名不代表科学身份。
浏览器“网页，完整”生成的 `*_files` 资源保持原名，不计入候选核心文件。
新产物写入 `ISODISTORT/output/validation/`，`output_compare/` 保持只读。
历史布局与下载叙述可从 [归档](history/README.md) 完整恢复。
