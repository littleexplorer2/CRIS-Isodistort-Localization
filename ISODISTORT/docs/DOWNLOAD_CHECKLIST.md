# 官网下载清单（2026-10-19 重新核对并冻结）

本清单按当前文件内容重新建立；旧冻结、下载进度和 debug 结论全部作废归档。三种母相全部 k 标签各做 nmod=0/1，合计 102 个查询。冻结的是有限下载与检验任务，不是科学验收通过。

这里的“全部 k 标签”指各标签及本清单指定的参数代表点，不是连续参数空间的全部数值。
当前可执行的 Method 2 查询均已列入，但未取得的 IR/结果列表仍决定最终候选数和逐文件名称；
Method 3/4 的新增任务尚待独立输入和有限范围确定，不能把本清单称作所有 Method 的最终完整文件表。

冻结范围及文件 SHA-256 见 [本批次冻结记录](manifests/checklist_freeze_20261019.json)；
逐项检验任务见 [检验清单](manifests/validation_checklist_20261019.json)。

## 当前应做什么

可按下列顺序执行：先补“先补齐的页面”中的 5 张列表，然后按正文晶体顺序、各晶体的 nmod=0/1
顺序处理查询。每个查询先保存 IR 页，再保存每个 IR 的结果列表，最后按列表真实候选逐行核对四件套；
已有且身份一致的文件不重下。新列表如出现尚未登记的 IR/候选，先补登记并建立其真实目录，不能沿用
旧候选数判断下载完成。标为待定义的 Method 3/4 暂停，不自行猜参数。

- 先补 5 张已知缺失的结果列表页；已有 94 个 Method 2 候选四件套不重下。
- 其余 99 个查询先保存 IR 列表和每个 IR 的结果列表，再按真实候选行下载缺失四件套。
- 4310 X nmod0：列表和 24 个候选均已有，直接复用。
- Method 1 三晶体已有 330 个候选；Method 3 40 个既有案例含 77 个候选，四件套都找到。Method 4 有 24 个案例根和额外响应/auto-origin 记录；已有文件全部保留，未重新宣称数值验证通过。
- Method 3/4 的旧冻结及 4310 空列表冻结已撤销；新案例需要独立输入，当前列为待定义，不能从官网答案反推输入。

## 统一文件管理

- 原始母相：`E:\CRIS\experiment_data\<母相文件名>`，只读。
- Method 2 页面：`E:\CRIS\webpage_info\<母相文件名>\`，平铺保存。
- IR 列表：`3. <K>_nmod<N> - ISODISTORT_ irreducible representation.html`。
- 每个 IR 的结果列表：`4.<i> <K>_nmod<N>_<IR> - ISODISTORT_ order parameter direction.html`；i 按该查询 IR 下拉顺序从 1 开始，每换 K/nmod 重置。保留 IR 的正负号。
- 四件套：`E:\CRIS\output_compare\<母相文件名>\官网\Method2\<K>\nmod<N>\<IR>_<OPD>\`。
- 文件为 `subgroup.cif`、`data.isoviz`、`topas.str`、模式详情 HTML。现存官网原名 `ISODISTORT_ complete modes details.html` 继续使用；新下载也用此名称，不加 6. 编号。保留 HTML 原引用的 `_files` 资源目录。
- Method 1 既有 `<IR>_<OPD>_SG<n>` 目录保留；Method 3 统一为 `<案例>/` 放结果列表、`<案例>/C<序号>_SG<n>/` 放四件套；Method 4 保留 `<案例>/`、`auto-origin/` 和明确标记的历史失败响应。
- 历史网页修复副本归 `webpage_info/<母相>/_variants/repaired_20261019/`，不当作新下载。所有搬移都有旧路径、新路径、SHA-256，不覆盖或删除官网输出。
- 已有候选可以依据其真实文件身份归位；列表缺失仍单独标记，不能把归位当作候选全集已证实。未知候选不预建占位目录。

## 官网操作

上传对应母相；Types 选 strain 和全物种 displacive，点 Change。Method 2 选择下列 K/参数，设置 # of independent incommensurate modulations 为 0 或 1。IR 页 OPD 留空，Use only one arm of the star 不勾选。保存网页选“网页，完整”，HTML 和资源目录一并保留。已有文件按当前索引复用，不用 nmod0 顶替 nmod1。

## 先补齐的页面

- `E:\CRIS\webpage_info\EuAl4 Parent.cif\4.2 LD_g1-6_nmod0_LD2 - ISODISTORT_ order parameter direction.html`
- `E:\CRIS\webpage_info\EuAl4 Parent.cif\4.3 LD_g1-6_nmod0_LD5 - ISODISTORT_ order parameter direction.html`
- `E:\CRIS\webpage_info\NdNiO2 own.cif\4.2 Y_a1-3_nmod0_Y2 - ISODISTORT_ order parameter direction.html`
- `E:\CRIS\webpage_info\NdNiO2 own.cif\4.3 Y_a1-3_nmod0_Y3 - ISODISTORT_ order parameter direction.html`
- `E:\CRIS\webpage_info\NdNiO2 own.cif\4.4 Y_a1-3_nmod0_Y4 - ISODISTORT_ order parameter direction.html`

以下完整查询清单中，参数化 K 后缀就是文件名中的 K，例如 a1-6 表示 a=1/6。未取得 IR 页时，IR 数量和全部 4.i 文件名暂不能确定，须按官网实际下拉框。

## EuAl4 Parent.cif

页面保存位置：`E:\CRIS\webpage_info\EuAl4 Parent.cif\`。

### nmod=0

- **GM**（k14；无自由参数）：保存 `3. GM_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GM_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **X**（k13；无自由参数）：保存 `3. X_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> X_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **M**（k15；无自由参数）：保存 `3. M_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> M_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **N**（k11；无自由参数）：保存 `3. N_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> N_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **P**（k12；无自由参数）：保存 `3. P_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> P_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **A_a1-6-g1-3**（k3；a=1/6, g=1/3）：保存 `3. A_a1-6-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> A_a1-6-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **B_a1-6-g1-3**（k1；a=1/6, g=1/3）：保存 `3. B_a1-6-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> B_a1-6-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **C_a1-6-b1-4**（k2；a=1/6, b=1/4）：保存 `3. C_a1-6-b1-4_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> C_a1-6-b1-4_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **DT_a1-6**（k7；a=1/6）：保存 `3. DT_a1-6_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> DT_a1-6_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **GP_a1-6-b1-4-g1-3**（k0；a=1/6, b=1/4, g=1/3）：保存 `3. GP_a1-6-b1-4-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GP_a1-6-b1-4-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **Q_b1-4**（k5；b=1/4）：保存 `3. Q_b1-4_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Q_b1-4_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **SM_a1-6**（k6；a=1/6）：保存 `3. SM_a1-6_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> SM_a1-6_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **W_g1-3**（k9；g=1/3）：保存 `3. W_g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> W_g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **Y_a1-6**（k8；a=1/6）：保存 `3. Y_a1-6_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Y_a1-6_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **LD_g1-6**（k10；g=1/6）：已有 22 个候选四件套；只补 4.2 LD2、4.3 LD5。

### nmod=1

- **GM**（k14；无自由参数）：保存 `3. GM_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GM_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **X**（k13；无自由参数）：保存 `3. X_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> X_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **M**（k15；无自由参数）：保存 `3. M_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> M_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **N**（k11；无自由参数）：保存 `3. N_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> N_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **P**（k12；无自由参数）：保存 `3. P_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> P_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **A_a1-6-g1-3**（k3；a=1/6, g=1/3）：保存 `3. A_a1-6-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> A_a1-6-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **B_a1-6-g1-3**（k1；a=1/6, g=1/3）：保存 `3. B_a1-6-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> B_a1-6-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **C_a1-6-b1-4**（k2；a=1/6, b=1/4）：保存 `3. C_a1-6-b1-4_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> C_a1-6-b1-4_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **DT_a1-6**（k7；a=1/6）：保存 `3. DT_a1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> DT_a1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **GP_a1-6-b1-4-g1-3**（k0；a=1/6, b=1/4, g=1/3）：保存 `3. GP_a1-6-b1-4-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GP_a1-6-b1-4-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **Q_b1-4**（k5；b=1/4）：保存 `3. Q_b1-4_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Q_b1-4_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **SM_a1-6**（k6；a=1/6）：保存 `3. SM_a1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> SM_a1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **W_g1-3**（k9；g=1/3）：保存 `3. W_g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> W_g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **Y_a1-6**（k8；a=1/6）：保存 `3. Y_a1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Y_a1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **LD_g1-6**（k10；g=1/6）：保存 `3. LD_g1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> LD_g1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。

## NdNiO2 own.cif

页面保存位置：`E:\CRIS\webpage_info\NdNiO2 own.cif\`。

### nmod=0

- **GM**（k17；无自由参数）：保存 `3. GM_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GM_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **M**（k18；无自由参数）：保存 `3. M_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> M_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **X**（k15；无自由参数）：保存 `3. X_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> X_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **A**（k20；无自由参数）：保存 `3. A_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> A_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **R**（k16；无自由参数）：保存 `3. R_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> R_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **Z**（k19；无自由参数）：保存 `3. Z_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Z_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **B_b1-4-g1-3**（k3；b=1/4, g=1/3）：保存 `3. B_b1-4-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> B_b1-4-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **C_a1-6-g1-3**（k5；a=1/6, g=1/3）：保存 `3. C_a1-6-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> C_a1-6-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **D_a1-6-b1-4**（k1；a=1/6, b=1/4）：保存 `3. D_a1-6-b1-4_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> D_a1-6-b1-4_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **DT_b1-4**（k8；b=1/4）：保存 `3. DT_b1-4_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> DT_b1-4_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **E_a1-6-b1-4**（k2；a=1/6, b=1/4）：保存 `3. E_a1-6-b1-4_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> E_a1-6-b1-4_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **F_a1-6-g1-3**（k4；a=1/6, g=1/3）：保存 `3. F_a1-6-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> F_a1-6-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **GP_a1-6-b1-4-g1-3**（k0；a=1/6, b=1/4, g=1/3）：保存 `3. GP_a1-6-b1-4-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GP_a1-6-b1-4-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **LD_g1-3**（k13；g=1/3）：保存 `3. LD_g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> LD_g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **S_a1-6**（k11；a=1/6）：保存 `3. S_a1-6_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> S_a1-6_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **SM_a1-6**（k10；a=1/6）：保存 `3. SM_a1-6_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> SM_a1-6_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **T_a1-6**（k7；a=1/6）：保存 `3. T_a1-6_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> T_a1-6_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **U_b1-4**（k9；b=1/4）：保存 `3. U_b1-4_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> U_b1-4_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **V_g1-3**（k14；g=1/3）：保存 `3. V_g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> V_g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **W_g1-3**（k12；g=1/3）：保存 `3. W_g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> W_g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **Y_a1-3**（k6；a=1/3）：已有 48 个候选四件套；只补 4.2 Y2、4.3 Y3、4.4 Y4。

### nmod=1

- **GM**（k17；无自由参数）：保存 `3. GM_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GM_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **M**（k18；无自由参数）：保存 `3. M_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> M_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **X**（k15；无自由参数）：保存 `3. X_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> X_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **A**（k20；无自由参数）：保存 `3. A_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> A_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **R**（k16；无自由参数）：保存 `3. R_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> R_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **Z**（k19；无自由参数）：保存 `3. Z_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Z_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **B_b1-4-g1-3**（k3；b=1/4, g=1/3）：保存 `3. B_b1-4-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> B_b1-4-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **C_a1-6-g1-3**（k5；a=1/6, g=1/3）：保存 `3. C_a1-6-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> C_a1-6-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **D_a1-6-b1-4**（k1；a=1/6, b=1/4）：保存 `3. D_a1-6-b1-4_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> D_a1-6-b1-4_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **DT_b1-4**（k8；b=1/4）：保存 `3. DT_b1-4_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> DT_b1-4_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **E_a1-6-b1-4**（k2；a=1/6, b=1/4）：保存 `3. E_a1-6-b1-4_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> E_a1-6-b1-4_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **F_a1-6-g1-3**（k4；a=1/6, g=1/3）：保存 `3. F_a1-6-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> F_a1-6-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **GP_a1-6-b1-4-g1-3**（k0；a=1/6, b=1/4, g=1/3）：保存 `3. GP_a1-6-b1-4-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GP_a1-6-b1-4-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **LD_g1-3**（k13；g=1/3）：保存 `3. LD_g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> LD_g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **S_a1-6**（k11；a=1/6）：保存 `3. S_a1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> S_a1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **SM_a1-6**（k10；a=1/6）：保存 `3. SM_a1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> SM_a1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **T_a1-6**（k7；a=1/6）：保存 `3. T_a1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> T_a1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **U_b1-4**（k9；b=1/4）：保存 `3. U_b1-4_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> U_b1-4_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **V_g1-3**（k14；g=1/3）：保存 `3. V_g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> V_g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **W_g1-3**（k12；g=1/3）：保存 `3. W_g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> W_g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **Y_a1-3**（k6；a=1/3）：保存 `3. Y_a1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Y_a1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。

## 4310_tetra.cif

页面保存位置：`E:\CRIS\webpage_info\4310_tetra.cif\`。

### nmod=0

- **GM**（k14；无自由参数）：保存 `3. GM_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GM_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **X**（k13；无自由参数）：已有完整列表和 24 个候选；不下载。
- **M**（k15；无自由参数）：保存 `3. M_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> M_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **N**（k11；无自由参数）：保存 `3. N_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> N_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **P**（k12；无自由参数）：保存 `3. P_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> P_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **A_a1-6-g1-3**（k3；a=1/6, g=1/3）：保存 `3. A_a1-6-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> A_a1-6-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **B_a1-6-g1-3**（k1；a=1/6, g=1/3）：保存 `3. B_a1-6-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> B_a1-6-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **C_a1-6-b1-4**（k2；a=1/6, b=1/4）：保存 `3. C_a1-6-b1-4_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> C_a1-6-b1-4_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **DT_a1-6**（k7；a=1/6）：保存 `3. DT_a1-6_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> DT_a1-6_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **GP_a1-6-b1-4-g1-3**（k0；a=1/6, b=1/4, g=1/3）：保存 `3. GP_a1-6-b1-4-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GP_a1-6-b1-4-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **LD_g1-3**（k10；g=1/3）：保存 `3. LD_g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> LD_g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **Q_b1-4**（k5；b=1/4）：保存 `3. Q_b1-4_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Q_b1-4_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **SM_a1-6**（k6；a=1/6）：保存 `3. SM_a1-6_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> SM_a1-6_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **W_g1-3**（k9；g=1/3）：保存 `3. W_g1-3_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> W_g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- **Y_a1-6**（k8；a=1/6）：保存 `3. Y_a1-6_nmod0 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Y_a1-6_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。

### nmod=1

- **GM**（k14；无自由参数）：保存 `3. GM_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GM_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **X**（k13；无自由参数）：保存 `3. X_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> X_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **M**（k15；无自由参数）：保存 `3. M_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> M_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **N**（k11；无自由参数）：保存 `3. N_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> N_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **P**（k12；无自由参数）：保存 `3. P_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> P_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **A_a1-6-g1-3**（k3；a=1/6, g=1/3）：保存 `3. A_a1-6-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> A_a1-6-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **B_a1-6-g1-3**（k1；a=1/6, g=1/3）：保存 `3. B_a1-6-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> B_a1-6-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **C_a1-6-b1-4**（k2；a=1/6, b=1/4）：保存 `3. C_a1-6-b1-4_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> C_a1-6-b1-4_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **DT_a1-6**（k7；a=1/6）：保存 `3. DT_a1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> DT_a1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **GP_a1-6-b1-4-g1-3**（k0；a=1/6, b=1/4, g=1/3）：保存 `3. GP_a1-6-b1-4-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> GP_a1-6-b1-4-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **LD_g1-3**（k10；g=1/3）：保存 `3. LD_g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> LD_g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **Q_b1-4**（k5；b=1/4）：保存 `3. Q_b1-4_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Q_b1-4_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **SM_a1-6**（k6；a=1/6）：保存 `3. SM_a1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> SM_a1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **W_g1-3**（k9；g=1/3）：保存 `3. W_g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> W_g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- **Y_a1-6**（k8；a=1/6）：保存 `3. Y_a1-6_nmod1 - ISODISTORT_ irreducible representation.html`，再逐 IR 保存 `4.<i> Y_a1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。

## 下载后的检验

每次记录真实路径、字节数、SHA-256、母相、k、nmod、IR、OPD，以及完整 (SG,basis,origin,s,i) 身份；列表与候选逐行对应。文件已有、列表齐全、身份一致、科学验收是不同状态。
本地独立运行并冻结输出后才读取官网答案比较；本轮已删除旧本地产物，所有新数值/浏览器验收从未执行状态开始。
机器清单：[下载范围](manifests/official_output_download_manifest.json)、[逐查询页面与候选现状](manifests/official_result_list_stage_20261019.json)、[官网文件登记](manifests/official_archive_reuse_20261019.json)、[检验清单](manifests/validation_checklist_20261019.json)。
本轮操作与限制见 [验证报告](BUGFIX_VALIDATION_REPORT.md)，归档恢复见 [历史索引](history/README.md)。
