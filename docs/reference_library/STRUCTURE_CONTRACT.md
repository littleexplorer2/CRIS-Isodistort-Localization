# 母相 CIF、晶格与对称轨道合同

适用实现：`cris_rea/backend/structure/`；定向测试：`tests/test_structure.py`。
这是新项目独立实现的普通三维结构层，不复制原 CRIS 算法，不读取官网案例答案。
本模块没有文件读写、网络、子进程、空间群查表或自动标准化行为。

## 来源及实际审读范围

本次在 2026-10-20 查阅以下官方定义；其作用是确定输入语义，不是给科学结论套用官方标签。

- [IUCr CIF 1.1 File Syntax](https://www.iucr.org/what-we-do/digital-standards/cif/cif1/file-syntax)：
  审读 §2.4、§7、§14–17、§22–26 及 Appendix A 的 Numeric/loop 语法。
  据此处理大小写不敏感的 tag、单层 loop、引号、分号文本、ASCII 和带标准不确定度的数字。
  本实现支持下面明确列出的子集；不自称完整 CIF 验证器。
- coreCIF 2.4.5 的
  [分数坐标定义](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_fract_.html)和
  [xyz 操作定义](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ispace_group_symop_operation_xyz.html)：
  实际审读 Definition、Type、关联旧 tag。位置用晶格基表示，操作为仿射作用。
- coreCIF 2.4.5 的
  [occupancy](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_occupancy.html)、
  [disorder group](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_disorder_group.html)和
  [Wyckoff symbol](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_Wyckoff_symbol.html)：
  实际审读各 Definition/Type/默认值。占位率缺省为 1；非空 disorder 标签拒绝；
  Wyckoff 字母只保存声明，不据标签生成坐标。
- [旧 multiplicity tag 的弃用说明](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_symmetry_multiplicity.html)
  明确指出历史软件曾混用轨道重数与位点群阶。本实现保存实际 tag，旧 tag 加 notice；
  与精确轨道不一致时诊断并提示此歧义，绝不自动把重数取倒数或改成位点群阶。
- 只读本地发行包 `vendor/isobyu/findsym.txt`，开头及 “Input using key words” 的
  `latticeParameters`、`latticeBasisVectors`、`atomType`、`atomPosition`、
  `atomOccupation` 段：长度 Å、角度 degree、分数坐标和原子类型区别的工具合同。
  该手册自称 Version 6 / January 2018，不能冒充当前二进制源代码或最新完整说明。
  既有 [独立合成输入](../../evidence/native/cif-converter-positive/synthetic_p1.cif)
  只作为格式回归输入；未读取或反推官网结构。

下列数学说明是**本项目推导**，与上列官方字段定义区分。没有使用历史猜测或个案答案。

## 支持范围与公开入口

从 `cris_rea.backend.structure` 导入：

- `parse_cif(text_or_bytes, metric_relative_tolerance=1e-8) → CrystalStructure`。
  只接收内容，不接收路径；源字节 SHA-256 随结构保存。返回前验证操作与所有轨道。
- `parse_cif_number(token) → CifNumber`，包括精确 `value`、`standard_uncertainty` 和原 token。
- `parse_symmetry_operation(text) → SymmetryOperation`。
- `validate_symmetry_operations(ops, lattice=None, metric_relative_tolerance=1e-8)`。
- `expand_structure(structure, metric_relative_tolerance=1e-8) → tuple[ExpandedSite,...]`。
- 数据对象 `Lattice`、`AtomicSite`、`CrystalStructure`、`SymmetryOperation`、
  `ExpandedSite` 为冻结 dataclass；失败统一为 `StructureError(code, message, **details)`。

接受一个 named `data_` block、ASCII CIF1 内容、普通 3D 六个晶格参数、一个 atom loop，
以及显式 xyz operation loop。loop 行可以跨物理行。quoted string、注释和分号长文本用于
普通元数据；数值字段必须是未加引号的 CIF 数字。数字支持尾部 ESD，例如 `1.234e-2(5)`。
`1/3` 不是 CIF 数字 token；它可以出现在 xyz 操作表达式中。缺失坐标 `?`/`.` 拒绝。

原子必须显式提供 label、type_symbol 和三个 fract 坐标。类型作为区分物种的字符串保存，
不从 label 推测元素；例如 `Fe2+` 与 `Fe3+` 不合并。occupancy 缺省 1，保留单组分
`0 < occupancy < 1` 并加 notice；后续工作流可进一步拒绝此情形。本模块明确不支持
混合占位、dummy/零占位或利用 ESD 容忍超过 1 的占位，因此范围比完整 coreCIF 更窄。

新旧 xyz tag、可选 IT number/H–M 名称、Wyckoff 字母和 multiplicity 见实现中的 alias。
没有操作时返回 `missing_symmetry_operations`，即使声明编号 1 也不自动补恒等。
声明 SG number/name、Wyckoff 字母不是已验证的群识别；调用方须以独立工具/通用数据库核对。
本模块验证的是**输入操作集的数学一致性**，不是证明其等于该声明编号的完整空间群。

CIF2、多 block、save frame/global/stop、磁性或调制字段、显式 disorder、Cartesian-only
原子输入、缺类型、重复 tag/语义 alias、不完整 loop 行都明确拒绝。位移椭球等无关字段
可以存在但不产生 ADP 计算能力。资源上限为输入 1 MiB、384 个操作、数字指数绝对值 300、
有效数字 token 的实现长度限制；这些是工程限制，不是晶体学完备性声明。

## 数据对象与约定

`Lattice(a,b,c,alpha,beta,gamma)` 的长度是 Å、角是 degree。
`matrix` 是三行三列 tuple，**列**为直接晶格向量；`metric`、`volume`、
`fractional_to_cartesian`、`cartesian_to_fractional` 均采用这一约定。
三角函数/Cartesian 使用 Python binary64；CIF 数值中心值、操作平移和原子分数坐标使用
`Fraction`，不以 `limit_denominator` 或容差做自动有理化。

`AtomicSite` 字段包括 `label`、`type_symbol`（只读别名 `symbol`）、`fractional`、
`occupancy`、`wyckoff_symbol`、`site_symmetry_multiplicity`、
`coordinate_uncertainties`、`occupancy_uncertainty`、`multiplicity_source_tag`。
输入坐标可以超出 [0,1)，原始中心值保留；展开结果才取模。

`ExpandedSite` 保留 `source_label`、零基 `source_index`（别名 `asymmetric_index`）、
全部生成操作的零基 `operation_indices`，以及类型、坐标、占位率。
它的确定性显示 label 为 `源label:轨道序号`；显示 label 不参与原子同一性判定。
轨道顺序为源站点输入顺序，再按显式操作首次生成顺序。

## 晶格推导及实现对应

取分数坐标列向量 \(x\)，Cartesian 列向量 \(r=A x\)，且
\[
a_1=(a,0,0),\quad
a_2=(b\cos\gamma,b\sin\gamma,0),\quad
a_3=(c\cos\beta,c q,c\sqrt{1-\cos^2\beta-q^2}),
\qquad q=\frac{\cos\alpha-\cos\beta\cos\gamma}{\sin\gamma}.
\]
这是选取右手 Cartesian 基后的构造：前两向量分别满足长度 a/b 与夹角 γ；
由 \(a_1\cdot a_3=ac\cos\beta\) 确定第一分量，
由 \(a_2\cdot a_3=bc\cos\alpha\) 确定第二分量，长度 c 再确定正第三分量。
因此 \(A^T A=G\)，其中对角为 \(a^2,b^2,c^2\)，非对角为相应长度乘积与夹角余弦。

行列式
\[
V=\det A=abc\sqrt{1+2\cos\alpha\cos\beta\cos\gamma
 -\cos^2\alpha-\cos^2\beta-\cos^2\gamma}>0
\]
给出 Å³ 的胞体积。`Lattice.matrix` 对应上述构造；逆变换按上三角 A 从 z、y、x
逐项回代，不转置晶格、不使用倒易格替代直接格。负长/非法角、非正 Gram、数值奇异
以及体积上溢/下溢明确报错。

数值边界 `abs(sin(gamma)) <= 1e-12` 或
`1-cos(beta)^2-q^2 <= 1e-12` 返回 `singular_lattice`；它限制接近退化时的相消/除法误差，
是保守计算域，不是物理上所有正体积胞都非法。

## 对称操作、群与轨道的推导

操作为 \(g=(R,t)\)，按主动分数坐标作用 \(g(x)=Rx+t\)。
R 为整数且 det R=±1，故将整数格双射到自身；t 是精确有理向量。
`parse_symmetry_operation` 用受限 AST 展开仿射表达式，不执行 `eval`；
只允许坐标一次式、数值加减/乘除，最终检查 R 的整数性和可逆性。

以给定胞的整数平移为等价关系，t 每分量取 \(t-\lfloor t\rfloor\)。
因为 R 整数，取模不会破坏群作用；复合为
\[
(R_g,t_g)\circ(R_h,t_h)=(R_gR_h,\ R_gt_h+t_g)\pmod{\mathbb Z^3}.
\]
`SymmetryOperation.compose` 按此顺序实现。有限集合须无重复、包含恒等且逐对闭合；
所有元素又可逆，因此每个元素的有限幂中存在逆元，这些条件足以保证输入集合为群。
未提供的操作不会通过闭包生成后悄悄加入输入。

Cartesian 长度不变要求
\[
R^T G R=G .
\]
`validate_symmetry_operations` 令 \(D=\operatorname{diag}(a,b,c)\)，计算
\(\max_{ij}|[D^{-1}(R^TGR-G)D^{-1}]_{ij}|\)，即逐项残差除以相应轴长度的乘积，
与显式 `metric_relative_tolerance` 比较。这相当于在单位长度的基方向比较内积，
避免很长的 c 轴掩盖 a/b 的失配。默认 \(10^{-8}\) 是 dimensionless 的输入一致性界限，
并非 FINDSYM 的 Å 容差，也不声称它由实验 ESD 推出。其保守量级远大于普通良态胞的
binary64 舍入误差；拒绝近退化胞降低数值条件恶化。调用方若显式改阈值，须记录值和理由；
本模块限定 [0,10⁻⁴]，不因某个结果失败自动放宽。检查不修改晶格或操作。

对源位置 x，轨道 \(O_x=\{g(x)\bmod 1:g\in G\}\) 通过精确 Fraction tuple 去重。
\(g_1x=g_2x\) 当且仅当 \(g_2^{-1}g_1\) 在稳定子 \(G_x\) 内，因此
\(|O_x|=|G|/|G_x|\)。生成同一位置的所有操作索引被保存；特殊位置不会因多条操作重复加原子。
已知 multiplicity 若不等于轨道大小则拒绝。每个轨道再施加任意群元仍是同一集合，
这是 `expand_structure` 的可检验协变不变量。

不同独立站点生成同一位置时，不做坐标平均、占位率相加或类型择一：
类型不同报 `mixed_occupancy_unsupported`，同类型报 `overlapping_independent_sites`。
这是保守的输入合同，不能接收同一非平凡轨道全部展开后又声明为独立代表的 CIF。
全胞输入可使用与其实际语义一致的显式 P1 操作；不得为消除诊断随意删掉真实对称信息。

### 十进制近似的边界

`0.333333` 与 `1/3` 是不同有理数；例如三个以小数近似写出的 1/3 平移不会精确闭合，
返回 `nonclosed_symmetry_operations` 并指出失败的操作对。没有自动吸附、猜分母或扩大容差。
如果真实输入需要重建操作，应由调用方依据来源、明确容差和正式通用数据作独立恢复并记录，
再传入新输入；不能从官网结果反推操作。

原子坐标也不吸附到特殊位置。`.00001(10)` 的区间含零并不授权把中心值改为零；
在反演下它仍得到两个精确位置。ESD 被保留但本模块不传播协方差、不拟合、不声称轨道
重数与带测量误差的理想空间群已完全一致。声明重数可以暴露这种冲突，缺声明时需上层核查。

模式导出也须尊重此边界。例如 P3 操作 `-y,x-y,z` 固定 `(1/3,2/3,z)`，但有限小数
近似一般不再精确固定；只输出一个代表并声称重数 1 会被本解析器正确拒绝。上层
`backend/modes/export.py:checked_cif` 因此先验证 H 不对称单元的实际重新解析/展开，
若失败则尝试 P1 全胞原子序列化并再次核对。舍入只在序列化核查中有明确界限，既不
改变本模块精确群作用，也不授权坐标吸附。P1 是文件表示，所求 H 及精确参考位置仍在
`modes.json`，表示原因和残差在 `report.json.cif_serialization`。把该 CIF 重新作为
母相输入时，本模块遵从其 P1 声明，不从自定义 H 注释重新建立操作。

## 验证与剩余边界

定向命令（在新项目目录）：

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
& 'E:\CRIS\.venv\Scripts\python.exe' -m unittest discover -s tests -p test_structure.py -v
```

使用已存在解释器只读加载标准库；不执行原项目回归、不修改原项目。
测试包含独立合成 P1、ESD 与指数、legacy/modern tag、CRLF/quote/text field、
无操作/无坐标/错误 loop/重复 tag、磁/调制/混占拒绝、晶格 Gram/体积/正逆变换、
glide 的平方、仿射复合、精确 1/3 与错误小数反例、三重六方旋转与错误直角晶格、
中心化、反演特殊/一般位置、轨道协变和声明重数反例。Cartesian 数值检查采用
11–13 位小数级断言，精确算术部分使用严格相等；这些不是所有输入的数值误差上界。

本模块改变了结构解析与坐标计算的科学语义，故不能归类为“纯界面修改”。
数学对象、约定、实现对应和反例已在上文给出。测试结果与源码指纹见
[STRUCTURE_VALIDATION.json](STRUCTURE_VALIDATION.json)。
它们验证此有限结构合同，不证明任意材料、所有 CIF 方言、群识别、IR/OPD、
模式完备性、应变、磁性或 ISODISTORT 网站的端到端等价性。
