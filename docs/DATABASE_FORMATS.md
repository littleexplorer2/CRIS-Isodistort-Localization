# ISOTROPY 数据文件的实际格式与能力边界

审计日期：2026-10-20（Asia/Hong_Kong）。目标只读目录：
`E:\CRIS\ISODISTORT\resources\isobyu`。
本报告与证据均写入独立的 `CRIS_REA`；没有执行源目录中的程序或修改数据。

## 结论与证据等级

本地 `const.dat` 和十个 `data_*.txt` 共 77,898,142 字节，**主要是可读的、带字段名的文本数据库**，
不是十个不透明二进制数据库。九个 `data_*.txt` 和 `const.dat` 的所有字节均属于 ASCII 可打印字符或
TAB/LF/CR；`data_ssg.txt` 除此以外包含 140 个 NUL 字节。所有这些文件均使用 LF 换行，
没有 CRLF 或裸 CR。该结论由逐字节扫描支持，不能从 `.txt` 扩展名单独得出。

证据入口：

- [完整文件清单、SHA-256、大小及未变更核查](../evidence/database/inventory.json)。
- `evidence/database/<原文件名>.inventory.json`：每个候选字段的行号、零基字节偏移、
  数据区间、区间 SHA-256、原始十六进制样本、词法 token 数量及首批 token。
- [可复现的只读扫描器](../tools/inspect_databases.py)。
- [官方文档快照身份与下载来源](../evidence/database/official/sources.json)。

本报告使用三个等级：

1. **直接观察**：文件字节、标识符、计数、数值范围、官方页面明确说出的行为，以及
   带原始记录的反汇编/反编译输出。反编译输出本身是观察，不代表其类型恢复必然正确。
2. **语义推断**：字段名与官方概念相符，足以定位后续逆向目标，尚不等于完整 schema。
3. **未恢复**：下标、维数、压缩编码、正负号、矩阵行列约定、单位和使用条件尚须读取例程、
   反汇编/伪代码及受控运行互相印证。不能把推断直接放入科学计算。

## 字节布局

十个 `data_*.txt` 都呈现重复的结构：单独一行的小写字段标识符，随后若干行空白分隔的整数、
逻辑字面量或双引号字符串；末尾有单独的 `end`。扫描器把符合
`[a-z][a-z0-9_]*` 的整行当作**候选字段**，不假设它们已经构成经验证的数据库 schema。
报告中的字段数包括末尾 `end`。

字符串经常有尾随空格，例如 `data_irreps.txt` 从偏移 0 开始的 `irrep_label` 段，
第一个样本是 `"GM1     "`。这些空格、引号、换行及 NUL 在扫描和指纹中原样保留。
`T`/`F` 行被当作逻辑字面量数据，不能误认成字段名。

`const.dat` 不采用命名字段结构。88 行的可见排列是整数标识符、十进制常数，
部分行再给可读表达式。例如原文件前几行含 `1 0.000000000000000 0`、
`2 1.000000000000000 1`、`4 0.500000000000000 1/2`。
这些是文件直接观察。常数表的正负配对现已获得下述原生反汇编证据；
“所有矩阵代码如何索引它”仍未恢复。

### NUL 的准确位置

`data_ssg.txt` 中唯一的非文本字节种类是 `0x00`，共 140 个，组成 28 段，每段 5 字节。
第一段从零基偏移 **17,499,217** 开始，最后一段从 **17,499,433** 开始，
均位于 `ssgn_lattice_nlabel_janner` 字段的数据区内。
字段标识符自身从偏移 **17,496,613**、第 **227879** 行开始。
这些 NUL 串被双引号包围，前面可见 `"3.217"`，后面进入 `issgn_lattice_nmod` 字段。

这更像固定长度字符条目的未填充区域，但这只是推断。没有来源程序或读取路径的进一步证据，
不能判定为合法占位、下载损坏或有意义的标签，也不能自动删除或替换 NUL。
`rg` 的默认文本搜索会在此处停止并报告 binary；调查整个文件需字节扫描或明确的二进制搜索模式。

## 各文件保存了什么

下述“用途”是由真实字段名加官方命令概念得到的语义定位。各字段的完整清单、
词法类型/数量及偏移保存在相应 JSON 中，不在文档重复维护全部数据。

### `data_space.txt`

584,128 字节，68 个候选字段。直接可见点操作标签和操作表（`point_op_label`、
`ipoint_op`、`ipoint_op_mlt`、`ipoint_op_inverse`），点群、晶格、k 点、
空间群符号，以及 setting、origin、轴选择和 Wyckoff 关联字段。
`space_label` 从偏移 47,642 开始，包含 **230 个字符串**。
`ipoint_op` 从偏移 3,927 开始，有 **648 个整数，范围 −1 到 1**；
不能仅凭 648=72×9 就把每九个整数解释成已验证的矩阵。

用途定位：通用三维空间群操作、命名和坐标/setting 转换基础数据。
**未恢复**：操作编码对应哪些平移分母、矩阵按行还是按列、setting 转换的主动/被动方向。

### `data_wyckoff.txt`

728,638 字节，23 个候选字段。可见 `iwyckoff_count`、`iwyckoff_pointer`、
`wyckoff_label`、`iwyckoff_fract`、位点点群及其元素/IR/特征标，
还有 `iwyckoff_pg_vector_basis`、`iwyckoff_pg_ellipsoid_basis` 和 subduction 字段。
开头的 `iwyckoff_count` 有 **230 个整数**；偏移 1,560 的 `wyckoff_label`
有 **1,731 个字符串**。这两个存储数量不能混称“1731 个晶体结构”。

用途定位：通用 Wyckoff 轨道、位点对称性、位移/椭球基以及表示关联。
**未恢复**：参数化坐标的编码，轨道与位点点群索引映射，基矢规范和归一化。

### `data_irreps.txt`

1,104,194 字节，24 个候选字段。可见不同体系的 IR 标签
（`irrep_label_bc`、`irrep_label_kov`、`irrep_label_zak` 等）、所属空间群、
image、Lifshitz/Felix、生成元、矩阵及 faintness 字段。
首段 `irrep_label` 有 **4,777 个字符串**。
偏移 404,559 的 `irrep_matrices` 有 **170,000 个整数，范围 0 到 15**。

用途定位：表示标签对照、表示和群的关联及表示矩阵编码。
**未恢复**：整数代码怎样还原为实际矩阵元素。把 `2` 直接当数值 2，或者未经读取例程
验证就映射到 `const.dat`，均不可靠。

### `data_little.txt`

10,612,956 字节，56 个候选字段。可见 little-group/k 标签与参数、
小群操作、IR、复共轭关系、k-star、Miller–Love/惯用胞的转换索引、
完整表示矩阵、real 变换及 subduction 数据。
`little_irr_full_count` 在偏移 1,998,846，值为 **10294**。
`little_irr_full_matrices` 在偏移 4,907,025，有 **2,220,000 个整数，范围 0 到 89**。
`little_k_count` 只有 **14 个整数**，不能误解为“每个空间群一项”。

用途定位：一般 k 向量、小群表示、星及诱导到完整空间群表示所需的数据。
**未恢复**：k 参数系数/分母的打包规则、星臂顺序、相位与复共轭约定、实表示的基变换。

### `data_images.txt`

146,986 字节，28 个候选字段。可见 image 标签、维数、阶、类型、Landau/active、
生成元、对角化和 subgroup/order-parameter-direction 相关数组。
`image_label` 包含 **141 个字符串**；`image_dir_label` 在偏移 39,178，
有 **700 个字符串**，样本包括 `P1`、`C1`、`S1`。

用途定位：表示的 image group 及其不变方向/子群的预制结构。
**未恢复**：具体 image 矩阵和空间群表示的对应、参数方向的约束及子群链的索引。

### `data_isotropy.txt`

6,112,902 字节，37 个候选字段。可见父群、子群、IR、OPD 标签/维数/自由参数、
basis、origin、star arms、ferroic、domain 和 subduction 字段，并保留多组 `_old` 字段。
`isotropy_parent` 在偏移 88,096，有 **15,239 个整数，范围 1 到 230**；
`isotropy_orderparam_label` 同样有 **15,239 个字符串**。
`isotropy_basis` 有 **137,151=15,239×9 个整数**，
`isotropy_origin` 有 **60,956=15,239×4 个整数**。

用途定位：以 IR/OPD 组织的通用 isotropy-subgroup 表，以及相关嵌入、畴和二级表示信息。
数量比是直接观察，可作为读取例程维数核查目标；它们**并未证明**
origin 的第四项一定是公分母、basis 的排列顺序或 `_old` 当前仍被使用。

### `data_magnetic.txt`

2,553,677 字节，57 个候选字段。开头 `mag_count` 值 **1651**，
`mag_bns_label` 与 `mag_unified_label` 均各有 **1651 个字符串**。
还有 OG 标号、磁性点操作/点群、格、type-IV supergroup 嵌入，
以及 `mag_iso_*` 的父/子群、IR、OPD、basis 和 origin 字段。

用途定位：通用磁空间群的命名、操作、嵌入和磁性 isotropy-subgroup 数据。
**未恢复**：反时间操作编码、BNS/OG/UNI 设置转换的完整合同、磁矩和伪矢量约定。

### `data_diperiodic.txt`

10,994 字节，14 个候选字段。`dg_label` 有 **80 个字符串**，
另有 point group、lattice、三维空间群关联、操作、生成元及惯用晶格/分母字段。

用途定位：二周期群（层群）的通用数据。仅有文件不能证明本次任一命令入口能够完整使用它，
也不能把其 80 项和三维空间群 230 项混用。

### `data_ssg.txt`

25,412,801 字节，56 个候选字段。`issgn_count=16697`，存在 label、
`nmod`/`nstars`、k vectors、矩阵、操作、空间群格、Yamamoto/Janner 标签、
超胞心化、k 参数依赖、IR 变换及 block diagonalization 字段。
`issgn_nmod` 有 **17,000 个整数，范围 0 到 3**；**17,000 是存储 token 数，
不能取代文件显式的 16,697 计数**。

用途定位：通用 superspace-group 相关数据；数据可容纳多种调制维数是字段证据，
不是本地 CLI 或 ISODISTORT 全部超空间能力已经可用的证据。
本文件的 NUL 情形见上文。

### `data_ssgmag.txt`

30,627,646 字节，12 个候选字段。`issgmag_count=325127`，另有非磁 SSG
和磁群关联、标签、操作数量/索引、点操作、分数平移与 `transmag` 字段。

用途定位：磁性 superspace-group 数据及其关联。**325127 是该文件声明的 count，
本次没有独立证明其去重、等价类或完备性含义**。

### `const.dat`

3,220 字节，无命名字段，88 行；整数代码、十进制实常数及部分表达式的文本列表。
用途定位：供程序使用的常数表。

本地 `iso` 的 GNU `objdump` 静态证据进一步恢复了常数访问合同：

- [init_constants_ 反汇编](../evidence/native/static/init_constants_.asm.txt)
  在 `0x457014–0x45703b` 调用 Fortran runtime 读取 4 字节整数索引和 8 字节 real，
  目标地址为 `const_common_ + index×8`。`iso.strings.txt` 中的 `const.dat`
  与这一例程中的字符串地址相符。
- 同一例程 `0x45708e–0x4570be` 从 index=2 开始，每次前进两个 double，
  将当前值通过固定 XOR 掩码写到紧邻的下一项。
  [掩码原始字节与 ELF 地址映射](../evidence/database/constant_sign_mask.json)
  证实地址 `0x5a6f10` 的 8 字节是 `0x8000000000000000`，只翻转 binary64 符号位。
  因而当前样本中，偶索引常数对应的下一奇索引确实是其相反数。
- [constant_ 反汇编](../evidence/native/static/constant_.asm.txt)
  的 `0x45717a–0x4571a3` 检查索引范围并按 `index×8` 返回 double；
  `0x4571aa–0x4571b4` 是首次初始化调用。

这属于**GNU 工具的静态反汇编观察**，不是 REA/Ghidra 反编译，也不是实际执行轨迹。
样本 SHA-256 为 `242b2c6ef93b00007d580eef27a4540018d94467e52e9a0cb81c4c835b37ddc8`。
它没有证明 `irrep_matrices` 或 `little_irr_full_matrices` 的每一整数都可直接传入 `constant_`；
某些值也可能是索引、哨兵或别的编码，完整调用与解码链仍待恢复。

## REA/Ghidra 对读取与使用路径的补充

下述五条记录来自对同一 SHA-256 的 `vendor/isobyu/iso` 副本执行的真实 REA MCP，
提供者为 **Ghidra 12.1.3**；原始响应均含结构化证据，且未置 `isError: true`。
成功响应没有显式 `isError` 字段，不能把缺省状态描述成原始 JSON 写了 `false`。
它们与上一节 GNU 静态证据
分别保留，不能把先前失败的请求也算作成功。原始响应、提取伪代码及 SHA-256 的关联见
[伪代码清单](../evidence/native/pseudocode/manifest.json)；
[本次数据库路径复核记录](../evidence/database/rea_readpath_review.json)保存实际审读范围、
源响应指纹、定位片段和未恢复项。伪代码是反编译器输出，不是官方 Fortran 源码，也未执行为新内核。

### 常数文件：文本记录到内存值

- `init_constants_`：
  [原始 REA 响应](../evidence/native/rea-linux-iso-recovered/07-init-constants.json)，
  Evidence ID `ev_c3934e6205016da351a9a92ccc97319261a336f2661c83e95d6a8d0a5440cfe6`。
  [可读伪代码](../evidence/native/pseudocode/init-constants.c)显示先读取 `ISODATA`，
  去掉路径尾随空格，再直接拼接 `const.dat`；没有在此处补 `/`。每轮 Fortran
  读取调用只传入一个 4 字节整数索引和一个 8 字节实数目标，然后结束该次读取。
  文本中的可读表达式没有对应的第三个取值调用。开文件失败进入 `bomb2_`。
- `constant_`：
  [原始 REA 响应](../evidence/native/rea-linux-iso-recovered/08-constant.json)，
  Evidence ID `ev_f52b50abf8a82a8fdff67c6627669397e75935686ac1fa0b0ba3e89d81c6b0e0`。
  [可读伪代码](../evidence/native/pseudocode/constant.c)显示通过初始化标志首次加载，
  检查传入索引至少为 1、至多为保存的上界，随后取一个 8 字节值。

反编译的填充语句为 `puVar3[1] = *puVar3 ^ 0x8000000000000000`，与 GNU 指令和
已读掩码一致。对有效已加载的偶数索引 `2j`，紧邻 `2j+1` 的 binary64 值取相反数。
伪代码把最大索引变量赋为循环保留索引加 1；这不是对任意损坏、乱序或越界常数文件的
健壮性证明。输入中原有 `1=0`、`2=1`、`4=1/2` 的意义来自文件记录；
相反数项的生成来自程序证据，二者不再混作字段名猜测。

`&const_common_ + index*2` 是 Ghidra 对已恢复基指针类型的算术展示；
结合 GNU 地址步幅与 `_gfortran_transfer_real(...,8)`，实际槽位步幅为 **8 字节**，
不能误读成每个代码跨两个 double。`constant_` 被恢复为 `undefined8` 返回类型，
仅靠这个类型名也不能否认其 double 语义或宣布完整 Fortran 声明已恢复。

### 表示计算：已找到使用链，矩阵代码仍有缺口

- `get_irrep4_`：
  [原始 REA 响应](../evidence/native/rea-linux-iso-recovered/10-get-irrep4.json)，
  Evidence ID `ev_14bd60fa1fca5963680c7f1cfd61949e8da31466cb5b6fab89b10b6f87edab71`。
  [可读伪代码](../evidence/native/pseudocode/get-irrep4.c)以第二个参数索引
  `common_little_` 的多段内存，取 2 字节和 1 字节条目，建立操作/参数相关缓存；
  可见 `generate_space_group_`、`vsub_`、`get_kvec_`、`ivdot_`、`vaddn_` 调用。
  之后调用 `get_irrep4_ssgn_`。首个参数为 2 时复制矩阵缓冲区，为 1 时按对角步幅
  累加到标量输出。这是矩阵/特征标路径的代码证据。
- `get_irreps_`：
  [原始 REA 响应](../evidence/native/rea-linux-iso-recovered/11-get-irreps.json)，
  Evidence ID `ev_876f3330452ee37940cf889e2523bf929f49ca0cfd8c7a84ee870ada181f9d86`。
  [可读伪代码](../evidence/native/pseudocode/get-irreps.c)支持直接索引与另一个编号映射：
  直接索引为零时，第三个参数须在 **1–4777** 内，然后从 `common_little_` 读取
  2 字节映射项。这与 `irrep_label` 的 4777 项数量相符，但数量一致仍不能独自证明
  该内存段就是哪一个文本字段。矩阵路径会核对调用方容量与内存中的维数；
  其首参数 3/4 分支另读取 `common_magnetic_` 的映射/符号相关值，再调用
  `get_irrep4_` 并将输出乘以选定的 `+1` 或 `−1`。

这两段函数没有直接调用 `constant_`；核心矩阵构造继续委托给本次未反编译的
`get_irrep4_ssgn_`。所以不能据此把 `irrep_matrices` 或
`little_irr_full_matrices` 的全部整数直接解释成常数表下标，特别是文件中可见的 `0`。
同样，`common_little_ + 偏移` 的取值没有自动恢复“文本字段 → 内存地址”的加载映射。
伪代码中的 `char`/`short`、局部数组大小和 `param_N` 是反编译器恢复结果；
Fortran 字符串隐藏长度、按引用参数、数组维数/列优先顺序及复数表示约定仍须独立核查。

### isotropy 文件：命名与按需生成协调，不是完整表读取 schema

`get_isotropy_file_` 的
[原始 REA 响应](../evidence/native/rea-linux-iso-recovered/14-get-isotropy-file.json)
对应 Evidence ID `ev_2310ee83c419aac9e2bae0bac220527bd066f4f3a5ac69c79e2eca5861c873b7`；
[可读伪代码](../evidence/native/pseudocode/get-isotropy-file.c)直接显示：

- 在 `superspace` 标志为零的分支，依据第二参数的 1–5 模式处理命名与生成；
  **这里的模式数字不能直接当作 ISODISTORT 网页 Method 1–4**。
  可见 `i     xx.iso`、`v     xx.iso` 模板及五位补零整数格式 `(i5.5)`；
  其他分支的常量字节形成 `a     xx.iso`，kernel 分支把前缀改成 `k`。
- 调用 `get_isotropy_filename_` 后检查两个返回位置的值。
  与其后“找不到数据库”“正在生成”的文本相结合，可把它们解释为存在/生成状态的控制量；
  这属于有邻近文本支持的语义推断，尚未审读该下层函数的实际文件查询或锁实现。
  找不到时若条件满足便调用 `get_isotropy_`，而非从固定文本表直接返回全部答案。
- 第一模式的一个缺失分支在 `www` 标志为零时提示用户继续，随后调用生成函数；
  `www` 非零会绕过这段交互。这个符号只证明本地样本有相应分支，
  **不能证明官网服务器就调用本二进制或使用相同配置**。
- 第五模式核对给定序参量维数与读出的表示维数；不相等时置状态 2。
  取消或显示“正在生成”的分支置状态 1。该函数中 `superspace` 非零会把文件名输出
  填为空格；这只是此函数的边界，不能外推为所有程序均不支持超空间。

本段函数没有直接打开 `data_isotropy.txt`、逐字段读取其 37 段、或解析 `.iso` 文件正文。
文件正文 schema、缓存位置、完整键构造、并发锁与生成算法仍在下层函数或其他入口中。
本次 `MAIN__` 伪代码请求
[15-main.json](../evidence/native/rea-linux-iso-recovered/15-main.json)
返回 `isError: true`；没有用空输出补造读取 schema，也没有重复失败请求。
因此前文 basis/origin 数量比、`_old` 字段和有效计数/容量的限制仍然成立。

这次复核只增加有来源的静态解释，没有改变生产解析器、数据库字节、科学公式或容差。
核查对象是响应身份、提取指纹、调用语句和局部控制流；未运行表示乘法、群闭包或
任意结构模式计算，不能由这些证据宣称科学内核已被完整重建。

## 官方文档支持哪些边界

[Linux 套件下载页](https://iso.byu.edu/isolinux.php)明确把上述系列列为 Data files，
把 `iso`、`findsym`、`smodes`、`comsubs` 列为可执行程序，并要求 `ISODATA`
指向资源目录且以 `/` 结束。当前源目录缺少该页列出的 `match_wyckoff.dat` 和
`symmetry.dat`，却多了页面未列出的 `data_ssgmag.txt`。
这是**文档列表与本地库存的差异**，不能单凭差异判定安装损坏、所缺文件必需或程序完整。
该页本身也不是服务端 ISODISTORT 的部署清单。

[ISOTROPY 命令手册](https://iso.byu.edu/isotropy_doc.php)的
`DISPLAY ISOTROPY` 和 `DISPLAY ISOTROPY COUPLED` 段说明：
参数化 k 的部分子群和耦合序参量子群从专门文件读取，不存在时可按用户选择计算生成。
因此数据库不是所有任意输入的预计算答案全集。
手册的 setting、irrep、basis 与 Wyckoff 命令也表明命名和坐标约定是操作合同的一部分。
本次实际审读的是命令目录、上述 DISPLAY 段及引言/setting 说明，
**没有把手册引用书目当作已阅读全文的数学证明**。

[官方版本历史](https://iso.byu.edu/isotropyversions.php)在 2021-09、v9.5 项说明：
ISOTROPY 移除了当时过时且不完整的 (3+1) superspace-group 功能，并指向 ISOSUBGROUP；
2022-01 的 v9.6/v9.6.1 项继续描述磁群标签与 IR 标签改动。
所以不能由 `data_ssg*.txt` 存在推断 `iso` 的 superspace 命令还可使用。
这个历史说明同样**不能证明眼前 `data_ssg.txt` 就是当年那份过时表**；
数据的来源版本还需要配套二进制、下载包和版本指纹建立关联。

以上三页的原始 HTML、URL、抓取时间和 SHA-256 已保存在 `evidence/database/official/`。
最初搜索返回的 `/iso/isolinux.php` 路径实际 GET 为 404；正文使用从官方首页
“ISOTROPY Software Suite for Linux” 链接确认的 `/isolinux.php`。

## 对本地计算与重建的意义

这些文件是通用晶体学/表示论数据，能够作为同版本本地程序的计算资源；它们不含本课题组
特定材料的官网计算答案，也不能替代从独立结构输入到结果的计算链。
“可读取文件”仅证明存储可访问，尚不证明逐字段的科学含义已恢复。

稳妥的第一条本地路径是保持原始数据和官方二进制配套，以显式输入、隔离工作目录和输出指纹
建立受控调用合同。直接重写科学内核之前，仍需恢复并验证：

- 各读库例程实际打开哪些文件、怎样分配数组、怎样处理容量与有效计数。
- 常数和矩阵代码的解码、下标基数、负号、复数/实表示转换。
- 点操作与平移的组合、setting/origin、直接格/倒易格、k-star 与矩阵基的约定。
- 通用数据查询、按需子群生成与结构专属模式计算之间的调用边界。
- 本地套件缺少而官网需要的入口、服务端程序、缓存/会话与结构导出层。

解码后的科学验收应检查操作闭包与逆元、表示乘法、维数、轨道重数、
坐标变换协变性及模式归一化。它们是后续验收项目，**本次词法扫描没有执行这些科学检验**。
不得为了补齐官网页面功能而按材料名称、OPD 标签或单次官网答案填表。

## 验证记录与复现

本次修改只增加只读证据扫描器及报告。扫描器没有晶体学变换、参数筛选、默认容差、
拟合或模式构造，因此没有改变科学计算语义；其不变量是原始字节、偏移、指纹和 source hash
前后相等。字节比较采用完全一致，不使用数值容差。物理单位/归一化尚未解码，故不适用于
本次词法检查，也不得据此宣称科学验证通过。

从 `CRIS_REA` 执行：

```powershell
& 'E:\CRIS\.venv\Scripts\python.exe' .\tools\inspect_databases.py --self-test
& 'E:\CRIS\.venv\Scripts\python.exe' .\tools\inspect_databases.py `
  --source 'E:\CRIS\ISODISTORT\resources\isobyu' `
  --output '.\evidence\database'
```

第一次命令通过：含空格的引号串、CRLF 字节偏移、NUL、逻辑值行及空文件的反例均满足断言。
第二次命令返回 0：11 个文件各自扫描前后 SHA-256 一致；每个候选字段的数据区间和样本可复查。
另一次独立字节复核逐一检查了 **375 个候选字段**的标识符区间、数据区间 SHA-256
和十六进制样本，结果见 [verification.json](../evidence/database/verification.json)。
扫描器拒绝把输出目录设在源目录内部。以上是工程与格式证据，不是全量科学回归或官网一致性验收。
未执行原项目任何回归、未写入原项目、未用官网结果文件作本地计算输入。
