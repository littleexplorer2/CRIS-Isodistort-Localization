# ISO 命令与精确协议边界

本报告记录当前 `vendor/isobyu/iso`（ISOTROPY 9.6.1，2022-01）的有限实测，
以及 `cris_rea/backend/iso_protocol.py` 的输入输出合同。它不宣称复现官网
ISODISTORT 的归一化、模式命名或所有 Method。实验输入是独立构造的空间群、
k、OPD、Wyckoff 位置与晶格参数；未读取任何官网案例答案作为计算输入。

官方命令文档的唯一归档是
`evidence/database/official/isotropy_doc.php.html`。本报告实际审读的相关段落为
DISPLAY DIRECTION（464–468）、DISTORTION（470–474）、LABEL ELEMENT
（558–574）、LABEL VECTOR（615–626）、SCREEN（638–646）、SHOW BASIS
（770–785）、CARTESIAN（787–806）、ELEMENTS（962–976）、FREQUENCY
（约1006–1030）、MICROSCOPIC（1169–1178）、UNITCELL（1368–1370）、
VALUE BASIS/CELL（1428–1453）、DIRECTION（1502–1526）、DOMAIN（1529–1550）、
IRREP（1652–1679）、KVALUE（1715–1728）、LATTICE PARAMETER（1771–1780）、
ORIGIN（1821–1831）、WYCKOFF/XYZ（1953–1971）。
官方行为说明不是这些科学计算的独立数学证明。

## 可运行的最小命令

每次运行启动一个新 ISO 进程。当前二进制启动时声明
`International (new ed.) with conventional basis vectors`；协议解析器核对这条
声明。显式命令 `SETTING INTERNATIONAL NEW` 虽出现在归档手册中，本二进制
实测却报 `isotropy_setting_commands 84` 并以零退出。因此 builder 保持已确认
的初始 setting，不发送这条命令。失败记录完整保留在
`builder-*-` 初批目录和 `screen-1000-primary/`，修正后记录带
`-default-setting` 后缀。禁止把初批零退出当作成功。
下方早期调查命令用 SCREEN 240；生产 builder 在确认 SCREEN 1000 可用后
改为1000以减少折行，解析器仍支持折行并拒绝带截断星号的输出。

母群的本原基矢与点余类代表：

```text
PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
LABEL ELEMENT INTERNATIONAL
VALUE PARENT 225
SHOW PARENT
SHOW BASIS
SHOW ELEMENTS
DISPLAY PARENT
QUIT
```

主序参量到子群：

```text
PAGE NOBREAK
SCREEN 240
LABEL VECTOR CONVENTIONAL
LABEL ELEMENT INTERNATIONAL
VALUE PARENT 221
VALUE IRREP GM4-
VALUE DIRECTION P1
SHOW PARENT
SHOW IRREP
SHOW DIRECTION VECTOR
SHOW SUBGROUP
SHOW INDEX
SHOW SIZE
SHOW BASIS
SHOW ORIGIN
SHOW ELEMENTS
DISPLAY ISOTROPY

QUIT
```

空行的作用是同意官方程序在当前隔离运行目录创建缺失的参数 k 子群表。
缺表提示明确要求 RETURN；`Y` 会取消显示而不是确认。若无需建表，空行只是
空命令。`allow_generate=False` 使用连续两个 QUIT：缺表时第一个作为非空拒绝，
第二个退出；无缺表时第一个直接退出。创建出的 `.iso` 表是官方本次计算结果，
保存在当前运行及证据副本内，不导入其他运行缓存。

对于有理参数 k，在 VALUE IRREP 后加 `VALUE KVALUE n,p1,...,pn`。
例如 SG1 的 `GP1GQ1` 使用 `VALUE KVALUE 3,1/3,1/4,0`。
特殊点省略 KVALUE；线、面和一般位置分别需要1、2、3个参数。这里的参数是
所选 k 族的参数，不能未经识别就把任意三维 k 坐标当作三参数输入。
VALUE KPOINT 会取消先前的 IRREP 选择，不能放在已选 IRREP 后随意补上。
生产流程在设置参数前先通过 `build_irrep_query` 查询 SHOW KPOINT、SHOW
DIMENSION 和 DISPLAY IRREP，`parse_irrep_output` 根据返回的有理仿射表达式
识别 a、b、g 参数；`IrrepMetadata.k_vector` 严格拒绝缺失或额外参数，并精确
求值为 Fraction。实测覆盖0/1/2/3参数，包括 `(a,a,a)` 仅需一个参数、
`(a,a,b)` 需要两个参数。参数个数不由 IR 名称猜测。

返回向量是名义 k 星代表，不保证是所选 OPD 的实际激活臂，不能拿这个单一
代表与 H 的平移直接做一般周期验收。当前 DISPLAY ISOTROPY 表格也不再次
回显用户赋值后的参数，特殊化后的 k 分支及星臂选择仍有独立核验边界。

目标 H 的全部相容 IR 可以通过另一新进程查询：

```text
PAGE NOBREAK
SCREEN 240
VALUE PARENT 221
VALUE SUBGROUP 123
VALUE BASIS 1,0,0 0,0,1 0,-4,0
VALUE ORIGIN 0,0,0
SHOW IRREP
SHOW KPOINT
SHOW DIRECTION VECTOR
DISPLAY DIRECTION
QUIT
```

该例返回 GM1+、GM3+、DT1(k参数1/4)、X1+。每行的 Subgroup 是该表示方向
自身的稳定子，可以严格包含目标 H，不应据此覆盖先前选定的 H。
DISPLAY DIRECTION 的数值 OPD 系数有时仅3位小数（如1.732），保留为说明
文本，不能作为高精度模式计算输入。

`SHOW FREQUENCY DIR` 在特殊 k 的 DISPLAY ISOTROPY 中提供更方便的
IR/方向/畴编号，例如221 GM4- P1得到 GM1+ P1(1)、GM3+ P1(3)、GM4- P1(1)。
但本二进制对于实测的 DT1(1/4) 与 GP1GQ1(1/3,1/4,0)，Frequency 栏为空。
空栏并不证明没有次级模式；通用有限超胞应交叉检查 DISPLAY DIRECTION。

## 精确坐标、基矢和心化

所有协议矩阵作用于列向量。若 ISO 打印三条基矢 b1,b2,b3，
解析结果 `basis` 是以它们为列的矩阵 B，坐标关系为
`x_parent = B x_subgroup + origin`。`primitive_basis` 采用同一列约定。
Seitz 操作是 `x -> R x + t`；R、t、B、origin 都用 Fraction 保存。
绝不先把 t 对母胞模1：例如 `(x,-y-1,-z)` 的整数平移在四倍超胞中有意义。

SHOW ELEMENTS 返回相对于本原平移子群的余类代表，而非相对于任意常规胞的
全部有限操作。225 Fm-3m 实测只有48个代表，本原基矢行列式绝对值1/4；
加上四个心化平移，常规胞有192个代表。160 R3m 是6个点余类、三个心化
平移，常规胞有18个代表。对应精确平移闭包核对在 `audit.json`。

221 GM4- P3 的 H=160，B 的行列式绝对值3，而 ISO Size=1。
将 H 自己的本原基矢 P 乘入 B 得到整数幺模矩阵 BP，行列式绝对值1。
因此 Size 是母子**本原胞**体积比，不是 `abs(det(B))`。实际作用群需要
结合 H 的本原平移（映射为 BP）和打印的 Seitz 代表，并在目标胞坐标中闭包。
本探针没有提供独立 CIF symmetry loop，所以未声称完成 CIF 操作等价验证。

## 位移投影及归一化边界

选定 IR/方向后，使用 `VALUE WYCKOFF <letter>`、`VALUE WYCKOFF XYZ x y z`、
`SHOW MICROSCOPIC VECTOR`、`DISPLAY DISTORTION` 获取原始极性位移投影。
同一 Wyckoff 字母的多个独立 site 必须分别设置 xyz 并分别显示；xyz 是当前
全局选项，不是任意多 site 的一张结构输入表。DISPLAY DISTORTION 只使用
VALUE IRREP 列表的第一个 IR。指定次级模式时应保留正确 VALUE DOMAIN，
不能把每个 IR 的默认畴混用。

`VALUE CELL b1 b2 b3` 指定显示的超胞窗口。它不检验周期相容性：
SG1、k=(1/3,1/4,0) 的周期12模式，在 CELL=I 下仍返回一母胞投影且不报错。
所以非空输出不能证明已经得到完整周期模式；必须独立核验周期与原子数。
`SHOW UNITCELL` 会主动只显示母胞，完整超胞工作流不应使用。

原始投影没有单位归一保证。SG1、晶格(4,5,6,80°,85°,75°)的 GM1 原始
分数位移为三条坐标轴；SHOW CARTESIAN 后变为晶格向量 a,b,c，长度为4、5、6，
且互不正交。221 GM4- P1 的多site例甚至有原始长度4与8的不同副本。
需要明确选定笛卡尔内积及归一化规范，不能把这些列称为官网归一模式。

SG1 的有理 k 双site例实测每site三块，每块12个位置、两列 cos/sin。
本项目独立核对这一**平移群例子**的公式为
`cos(2π k·t)` 和 `-sin(2π k·t)`，三个笛卡尔/分数轴分别重复。
数字0.866等只有3位；与三角函数差异约2.54e-5，落在输出舍入界5e-4内。
这不是放宽科学求解容差，也不是所有空间群模式的公式。

同一12胞的 DISPLAY DIRECTION 返回 GM1、Y1 和五个不同参数的 GP1GQ1，
实表示维数1+1+5×2=12。双site完整位移空间有72维，单一 primary IR 投影
只有12维，直接显示主 IR 无法代替完整 H 允许的畸变空间。

SHOW CARTESIAN 的 Point 输出也不能盲读：自由 xyz 实测正确转换，
但221的固定 b 点在晶格a=4时仍打印 `(1/2,1/2,1/2)`，同时位移已乘4。
此显示差异单独保留，不能据此从 Point 列反构 CIF 坐标。
极性模式完全禁戒时（221 GM1+ at1a）程序只打印提示符，没有表头；
仍需要独立维数判据区分物理零空间和查询失败。

## 实现、验证与复现入口

`iso_protocol.py` 只构造命令和解析精确空间群数据，不执行外部程序，不打开
evidence，不归一化或生成模式。执行仍走现有 RunService/runtime 隔离链。
这次模块增加改变了输入/坐标合同，已按上文的矩阵、Seitz 平移与心化不变量
核查；没有新增科学近似或更改科学默认容差。拒绝浮点/小数坐标作为精确基矢，
拒绝截断、显式错误、多个子群结果、丢失字段或不匹配的选择。

探针位于 `evidence/phase2/iso-contract/<case>/`，每例都有运行前冻结的 input、
来源记录、原始 stdout/stderr、完整运行证据；`generated-work/` 保存小型本次
生成表和日志。`index.json` 索引其哈希，`audit.json` 保存独立核查结果。
工程测试 `tests/test_iso_protocol.py` 读取这些冻结的本地二进制输出，只测试协议；
它不读取网站案例、不构成整个结构到模式工作流的科学验收。
最终19项协议定向测试全部通过；真实命令、退出码、日志和源码/测试SHA记录在
`evidence/phase2/iso-contract/protocol-tests-with-k-metadata/record.json`。
较早的14项记录保留在 `protocol-tests/`，不代表新增 metadata 的验证。
通过项包括0/1/2/3参数的真实表达式、缺/多参数反例、最终宽屏 builder 的
真实母群输出，以及 DT1(1/4) P2 的非零原点和整数 Seitz 平移。
独立 Codex 代理已只读复核 basis 转置、Seitz 平移保留、唯一结果及截断拒绝；
未运行全量回归，该复核也不替代最终结构工作流的科学验收。

```powershell
python -m unittest discover -s tests -p test_iso_protocol.py -v
python -m cris_rea run iso --input evidence/phase2/iso-contract/elements-221-gm4minus-p1/input.txt
```

后一个命令将创建新的 runs UUID，不覆盖已有证据。
`tools/probe_iso_contract.py` 是调查脚本，已有 case 目录时主动拒绝覆盖；
按留存 input 复跑更直接。`tools/audit_iso_contract.py` 只读取冻结输出进行
有限数学检查及小型产物归档，不能被生产计算链导入。
