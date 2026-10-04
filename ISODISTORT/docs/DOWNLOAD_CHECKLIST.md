# ISODISTORT 人工下载与上传清单

本文件只列仍可能需要用户执行的官网操作。已完成批次的逐例结果、验证结论和机器
报告索引写入 [BUGFIX_VALIDATION_REPORT.md](BUGFIX_VALIDATION_REPORT.md)，精确机器输入
写入 `docs/manifests/`；这里不重复维护完成案例表。

## 当前动作

**截至 2026-10-04，没有需要重新下载或上传的文件。**

- Method 1：4310、EuAl4、NdNiO2 的官网与现有网页版保存集均已通过完整性门禁。
- Method 2：EuAl4、NdNiO2 的保存集已复验；4310 尚未形成需要下载的正式批次。
- Method 3：EuAl4、NdNiO2 的 40 个查询、77 个候选已核定，不重下。
- Method 4：EuAl4、NdNiO2 的 24 个案例已核定；EuAl4 G05 的 auto-origin 仅有证据
  警告，现有截图与结果身份足够，不要求重跑。

上述完成状态不得推导为 4310 Method 2–4 已通过，也不得推导为任意晶体均正确。

## 4310 Method 3/4 的触发条件

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
output_compare/4310_tetra.cif/官网/Method3/<manifest 中的稳定案例号>/
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
basis/origin/matching 截图或 PDF、结果页、Complete modes details、CIF、TOPAS、
IsoVIZ 和 amplitude TXT/CSV（若官网提供）；预期失败案例保存错误页和完整填写状态。
下载后先核对上传文件 SHA-256、完整身份、模式数和产物集合，再做数值差分。

目标根目录为：

```text
output_compare/4310_tetra.cif/官网/Method4/<manifest 中的稳定案例号>/
```

## 文件与目录约定

```text
output_compare/
  <母相 CIF>/
    官网/
      Method1/
      Method2/
      Method3/<短案例号>/
      Method4/<短案例号>/
    现有网页版交互/
      Method1/
      Method2/
      Method3/<短案例号>/
      Method4/<短案例号>/
```

Method 1/2 的候选目录使用 `<IR>_<OPD>_SG<number>` 与 `<IR>_<OPD>`；Method 3/4
使用 manifest 中的稳定短案例号隔离独立运行。目录短名只用于定位，候选身份始终由
完整科学身份和查询上下文确定。完整旧名与新路径见只读
`output_compare/_FOLDER_NAME_MAP.json`。

浏览器“网页，完整”生成的 `*_files` 是 HTML 内部资源目录，保留浏览器原名。它不算
候选核心文件，也不能替代结果表。新本地结果先写入 `ISODISTORT/output/validation/`
并完成审计，程序不得直接覆盖只读 `output_compare/`；是否归档由用户决定。

4310 的在线结构查看页因模式数超过 200 而省略幅度输入框，相关 HTML 已保存在只读
`webpage_info/`。本地网页版没有在线结构查看功能，该文件只作官网行为存档；高模式数
仍要纳入模式生成与导出压力回归。
