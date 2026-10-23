# 官网下载清单（2026-10-20 按 Method 2 算法缺口重新筛选）

**当前先补第 3、4 项的列表；第 1、2、5 项列表已下载，第 6 项暂停。**
本次新增页面已核对并登记：第 1 项 3 个 IR、12 行候选，第 2 项 4 个 IR、8 行，
第 5 项 1 个 IR、232 行。列表齐全不表示候选四件套已经下载或科学验证通过。
第 6 项保存的是官网错误响应“参数使 k 点落在 C，而非 GP”，不是有效 IR 列表；
保留原件，暂不重试或据此下载四件套，不自行改参数凑结果。

原有 3 个查询的列表与 94 个候选四件套继续复用，不重下。原来的全部 k 标签、
每点 nmod=0/1 共 102 个查询不再是必下载范围；退出当前任务的 93 个查询也不再计作缺失下载。
冻结的是 9 个选定查询：6 组列表齐全、2 组列表待补、1 个错误响应暂停；
科学验收仍未执行。下文保留原第 1–6 项编号，避免混淆已下载文件。

已补齐的五张列表现已登记：EuAl4 的 `4.2 LD2`、`4.3 LD5`，
NdNiO2 的 `4.2 Y2`、`4.3 Y3`、`4.4 Y4`；全部是原 nmod=0 查询，
完整文件名和指纹见 [页面索引](../../webpage_info/_PAGE_INDEX.json)。
这五张列表及其资源不需要再次下载。

## 筛选原则与现有材料

- EuAl4：已有 LD/k10、g=1/6、nmod=0，1 张 IR 页、3 张 OPD 列表、22 个候选。
- NdNiO2：已有 Y/k6、a=1/3、nmod=0，1 张 IR 页、4 张 OPD 列表、48 个候选。
- 4310：已有 X/k13、nmod=0，1 张 IR 页、8 张 OPD 列表、24 个候选。
- 三晶体 Method 1 的 330 个候选已有特殊 k 的模式/结构材料，因此不再为了重复检查共享
  特殊 k 算法而逐个补 Method 2 的特殊 k 查询。Method 1 原件只按来源引用，
  **不能改称 Method 2 查询结果，也不能证明 Method 2 各特殊 k 列表完整**。

最初按正 nmod、二参数/三参数波矢和更大 star 的证据缺口选取下面的查询。
现在已取得 LD/Y 的正 nmod 列表和 GP 的 nmod0 列表，二参数 B 的列表仍待补，
GP 的正 nmod 请求则须先核定实际输入和官网参数处理。
已有文件可用于哪些检查、为何新增、哪些仍无法据此证明，见
[筛选与理论边界](reference_library/notes/METHOD2_TARGETED_EVIDENCE_20261020.md)。
上述三个原有查询的列表行数与候选目录数相符；新查询仅补齐列表，
逐候选完整身份、四件套和数值比较仍须另行完成。

## 共用官网设置

1. 从 `E:\CRIS\experiment_data\` 上传任务指定的原始母相。
   **跨晶体切换任务时重新上传母相，并核对页面顶部 Space Group、晶格参数及元素。**
   EuAl4/4310 都是 `139 I4/mmm`，但晶格及元素不同；NdNiO2 是 `123 P4/mmm`。
2. Types of distortions 勾选 **Strain + 全物种 Displacive**，点击 Change；
   不选 Occupational、Magnetic、Rotational。
3. 使用 Method 2，保持单个 IR（number of superposed IRs = 1），按任务填 k、参数和 nmod。
   未使用的参数留空。参数用下文的精确分数，不自行改值。
4. 保存 IR 列表；随后依官网实际下拉顺序逐个选择 IR，OPD 留空，
   `Use only one arm of the star` 不勾选，获取并保存每个 IR 的全部 OPD 列表。
5. 按真实列表逐行核对候选，下载尚缺的四件套。没有实际列表行，不预建候选目录，
   不用本地程序输出补官网结果。

官网的 nmod 表示独立非公度调制数，**不是谐波最高阶，也不是 a/b/g 参数个数**。
本清单保留独立输入合同中的有理代表值；nmod=1 查询用于取得官网实际响应和比较合同，
不能把它称作真实无理波矢的完整验证。若官网拒绝该设置、返回空列表或不提供某种导出，
保存完整响应与输入上下文，记为待判断/不适用；不要改参数凑结果或拿 nmod=0 文件填入。
这些异常响应也采用对应的 `3.` 或 `4.<i>` 页面路径，记录实际角色，不冒充正常结果列表。

## 统一文件与 HTML 编号

- 页面平铺在 `E:\CRIS\webpage_info\<母相文件名>\`。
- IR 列表的流程序号始终是 **3.**。
- OPD 列表序号为 **4.1、4.2、4.3……**，严格按本次 IR 下拉顺序。
  **每换 k/参数/nmod 都从 4.1 重新开始**；不同查询靠文件名后缀区分。
  新查询的 IR 数量和标签尚未知，不能按旧查询猜序号。
- HTML 均保存为“网页，完整”；保留对应 `_files` 目录及原资源引用。
- 每个候选四件套放在任务给出的 `<IR>_<OPD>` 目录，文件名固定：
  `subgroup.cif`、`data.isoviz`、`topas.str`、
  `ISODISTORT_ complete modes details.html`。
  模式详情 HTML 在四件套目录中**不加 6. 前缀**；保留其资源目录。
  `Distortion file` 不是 `data.isoviz`。
- IR 的 `+/-` 符号保留；若真实列表出现同名 IR/OPD 的不同完整身份，
  先按 SG、basis、origin、s、i 消歧登记再建目录，不能覆盖。

## 本次新增下载顺序
### 1. EuAl4 Parent.cif：LD / k10，nmod=1

**列表已齐全，不重下。** LD1/LD2/LD5 分别为 2/2/8 行。正 nmod 页面显示
`g=0.16667`，与原输入 `1/6` 的五位显示一致，但不能从舍入显示反证提交值精确相等。
结果为 4 维超空间身份，后续不能套用 nmod0 的三维身份比较或候选目录。

- 设置：**g=1/6**；官网坐标 `(0,0,1/6)`。
- 补足的证据：与已有同参数 nmod=0 构成对照，补 I 心格下正 nmod 查询、模式筛选和导出上下文的官网响应证据。
- IR HTML：`E:\CRIS\webpage_info\EuAl4 Parent.cif\3. LD_g1-6_nmod1 - ISODISTORT_ irreducible representation.html`。
- OPD HTML：`E:\CRIS\webpage_info\EuAl4 Parent.cif\4.<i> LD_g1-6_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- 四件套目录：`E:\CRIS\output_compare\EuAl4 Parent.cif\官网\Method2\LD_g1-6\nmod1\<IR>_<OPD>\`。

### 2. NdNiO2 own.cif：Y / k6，nmod=1

**列表已齐全，不重下。** Y1/Y2/Y3/Y4 各 2 行。正 nmod 页面显示 `a=0.33333`，
按五位显示记录，四维身份和四件套仍须单独核定。

- 设置：**a=1/3**；官网坐标 `(1/3,1/2,0)`。
- 补足的证据：与已有同参数 nmod=0 构成对照，补 P 格、带固定 1/2 分量的正 nmod 查询证据。
- IR HTML：`E:\CRIS\webpage_info\NdNiO2 own.cif\3. Y_a1-3_nmod1 - ISODISTORT_ irreducible representation.html`。
- OPD HTML：`E:\CRIS\webpage_info\NdNiO2 own.cif\4.<i> Y_a1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- 四件套目录：`E:\CRIS\output_compare\NdNiO2 own.cif\官网\Method2\Y_a1-3\nmod1\<IR>_<OPD>\`。

### 3. 4310_tetra.cif：B / k1，nmod=0

**待补列表。先重新上传 `E:\CRIS\experiment_data\4310_tetra.cif`。**
确认页首为 `139 I4/mmm`，晶格 `a=b=3.85160 Å，c=27.96700 Å`，元素包含 La/Ni/O。
下拉框必须是 **`B, k1 (a,0,g)`**；填写 **a=1/6、b 留空、g=1/3、nmod=0**。

若看到截图中的 **`B, k3 (0,b,g)`**，其定义与已存 NdNiO2 的 `123 P4/mmm`
搜索页一致，尚未处于本项预期上下文。先检查母相/空间群，不要把 a 的值改填 b
后继续保存到 4310 目录。仅凭裁剪截图不能确认当前实际上传的文件。

- 设置：**a=1/6, g=1/3**；官网坐标 `(1/6,0,1/3)`。
- 补足的证据：补二参数波矢、8 臂 star，以及多个独立物理轨道在参数 k 模式构造中的处理证据。
- IR HTML：`E:\CRIS\webpage_info\4310_tetra.cif\3. B_a1-6-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`。
- OPD HTML：`E:\CRIS\webpage_info\4310_tetra.cif\4.<i> B_a1-6-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- 四件套目录：`E:\CRIS\output_compare\4310_tetra.cif\官网\Method2\B_a1-6-g1-3\nmod0\<IR>_<OPD>\`。

### 4. 4310_tetra.cif：B / k1，nmod=1

**待补列表。沿用第 3 项核对后的 4310 母相与 B/k1 定义，仅将 nmod 改为 1。**
填写 **a=1/6、b 留空、g=1/3、nmod=1**。若官网返回错误，保留响应，
暂停本项的候选导出；本地有理小群检查不保证官网正 nmod 请求被接受。

- 设置：**a=1/6, g=1/3**；官网坐标 `(1/6,0,1/3)`。
- 补足的证据：补二参数波矢、8 臂 star，以及多个独立物理轨道在参数 k 模式构造中的处理证据。
- IR HTML：`E:\CRIS\webpage_info\4310_tetra.cif\3. B_a1-6-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`。
- OPD HTML：`E:\CRIS\webpage_info\4310_tetra.cif\4.<i> B_a1-6-g1-3_nmod1_<IR> - ISODISTORT_ order parameter direction.html`。
- 四件套目录：`E:\CRIS\output_compare\4310_tetra.cif\官网\Method2\B_a1-6-g1-3\nmod1\<IR>_<OPD>\`。

### 5. NdNiO2 own.cif：GP / k0，nmod=0

**列表已齐全，不重下。** 只有 GP1，一个 OPD 列表含 232 行。
这些是候选行数，不能记成已取得 232 个四件套。

- 设置：**a=1/6, b=1/4, g=1/3**；官网坐标 `(1/6,1/4,1/3)`。
- 补足的证据：补三参数一般波矢、16 臂 star 和一般方向的相位/共轭臂处理证据；选择原子较少的 NdNiO2 控制下载与后续计算成本。
- IR HTML：`E:\CRIS\webpage_info\NdNiO2 own.cif\3. GP_a1-6-b1-4-g1-3_nmod0 - ISODISTORT_ irreducible representation.html`。
- OPD HTML：`E:\CRIS\webpage_info\NdNiO2 own.cif\4.<i> GP_a1-6-b1-4-g1-3_nmod0_<IR> - ISODISTORT_ order parameter direction.html`。
- 四件套目录：`E:\CRIS\output_compare\NdNiO2 own.cif\官网\Method2\GP_a1-6-b1-4-g1-3\nmod0\<IR>_<OPD>\`。

### 6. NdNiO2 own.cif：GP / k0，nmod=1

**暂停，不作为已完成列表，也不继续下载。**
原计划输入为 a=1/6、b=1/4、g=1/3、nmod=1；
已保存的同名 HTML 正文是 `Your choice of parameters puts the k point at C not GP`，
没有 IR 选项或 OPD 列表。

该响应没有保留所提交的 k/参数/nmod 字段，任务归属目前依据文件名和用户报告；
不能据此断言是某个参数被填错，也不能确认原计划参数必然会触发这一响应。
需要先核定提交上下文及正 nmod 参数处理，再决定是否重试或调整独立输入合同。
暂不提供未经核实的新分数。

原件保留在
`E:\CRIS\webpage_info\NdNiO2 own.cif\3. GP_a1-6-b1-4-g1-3_nmod1 - ISODISTORT_ irreducible representation.html`，
索引角色为 `error_response`；以后即使重试成功，也不能覆盖这份错误证据。

## 暂不新增的任务

本次不再要求上述六项以外的 Method 2 下载。特别是已有的 EuAl4 LD nmod0、
NdNiO2 Y nmod0、4310 X nmod0 不重下；其他特殊 k 及同类参数 k 不为凑齐标签而下载。
这表示它们退出当前必下载范围，不表示这些查询通过了验证。

继续核对已有真实列表，并补第 3、4 项，再检查是否存在所需的模式/轨道/相位见证。若没有，
必须明确指出缺少的可观测量，才能增补另一点；不自动恢复旧的 102 查询清单。
不保证原六个新增查询构成数学上最小或完备的测试集；第 6 项已遇到错误响应，
其待检验缺口仍保留。也不以相同 star 大小证明算法等价。
nmod=2/3、真实非公度完备性、多 IR 耦合以及不支持的物理类型不在本批次验收范围。

Method 1 的现有材料继续引用。Method 3/4 新下载仍待独立输入、有限范围和来源核定，
旧冻结不恢复，现存官网材料全部保留。ISOVIZ_INPUT 不纳入当前计划。

## 下载后的登记与检验

记录真实路径、SHA-256、母相、k/参数、nmod、IR/OPD 和完整候选身份。
页面齐全、四件套齐全、身份一致、科学数值通过分别记录，不能相互替代。
独立固定输入并计算本地结果，冻结后才只读比较官网输出；官网结果不得进入计算链。
未取得页面时，最终候选数与文件总数仍未知。

机器依据：[下载范围](manifests/official_output_download_manifest.json)、
[页面与候选](manifests/official_result_list_stage_20261019.json)、
[筛选理由](manifests/method2_algorithm_coverage_20261020.json)、
[官网原件登记](manifests/official_archive_reuse_20261019.json)、
[检验任务](manifests/validation_checklist_20261019.json)、
[冻结记录](manifests/checklist_freeze_20261019.json)。
旧清单完整保存在 [历史归档](history/README.md)，恢复只能解压到新的目录，不能覆盖当前记录。
