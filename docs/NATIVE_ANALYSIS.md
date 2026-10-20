# Linux 原生程序：证据与计算路径

本报告针对 `evidence/environment/source_manifest.json` 冻结的二进制，分析日期
2026-10-20（Asia/Hong_Kong）。源目录仅只读；静态分析和执行都针对桌面副本。
官网下载套件不等于官网 PHP 服务端源码。这里没有取得其源代码。

## 证据等级与限制

REA MCP `open_binary` 识别 `iso` 为 x86-64 ELF，
SHA-256 `242b2c6ef93b00007d580eef27a4540018d94467e52e9a0cb81c4c835b37ddc8`。
首次 Windows MCP 的 `binary_overview` 和 `search_procedures` 返回
`provider_unavailable`；Ghidra 未配置且 Windows P0 不接受 ELF，其他提供器不可用。
完整错误在 `evidence/native/rea-provider-0.json` 与 `rea-provider-1.json`。
初次 REA 导出包有 0 个 Evidence record，**不把它当作反编译成功**。

随后在桌面项目内安装固定版本 REA 6.7.0、Ghidra 12.1.3 和 Linux JDK，
通过真实 stdio MCP 完成 `initialize`、`tools/list` 和 `tools/call`。
首次 Linux 分析因 Windows 挂载盘的 Unix socket 失败；调整 IPC 路径后，
`evidence/native/rea-linux-iso-recovered/04-overview.json` 已成功返回
`ev_10fa485dcea3acd26f8ce756770157bd685a4c3ffd9c7789c8b18315152dafd2`，
以及下述函数级 Evidence。失败记录与恢复记录分别保存，不覆盖。
环境及通信根因见 [NATIVE_RUNTIME](NATIVE_RUNTIME.md) 和
[GHIDRA_TRANSPORT_RECOVERY](GHIDRA_TRANSPORT_RECOVERY.md)。

`evidence/native/pseudocode/` 是从 REA 原始响应提取的可读伪代码，
其 `manifest.json` 关联 Evidence ID、原响应指纹和派生文件指纹。
这些文件**不是原始 Fortran 源码，也不能直接编译**。REA 明确报告此桥验证版本
为 Ghidra 12.1.4，本次观察版本为其接受的 12.1.3；`project_` 另有
`Type propagation algorithm not settling` 警告。类型和复杂数组解释仍需逐项核对。

本次成功导出 `evidence/native/rea-iso-bundle.json`（12 条 Evidence）
及关闭时的 `rea-iso-snapshot.json`，关闭响应保存在第 17 条 MCP 记录。
恢复后的 overview 观察到 2,015 个 procedure、6,409 个 string 和 32 个 segment；
这些是分析器识别数量，不代表逐一审查的算法数量。本次取得十个目标函数的可读伪代码。
对巨大 `MAIN__` 的伪代码请求仍返回 `execution_failure`（第 15 条记录），
没有给出成功伪代码；其 GNU 反汇编仍可读，但主调度器的完整高级结构尚未恢复。
批次脚本正常结束不抹去这一工具级失败。

`evidence/native/static/*.elf.txt` 是 GNU readelf 输出；`*.symbols.txt` 是 nm
输出；`*.asm.txt` 是 objdump 反汇编；`*.strings.txt` 是字符串及文件偏移。
它们是独立工具证据，不冒充 REA Evidence ID，也不是原 Fortran 源代码。
符号地址只适用于此 SHA 的 ELF。符号存在不代表某次运行调用了该函数。

## 已观察到的程序和调用关系

五个 ELF 是 `iso`、`findsym`、`findsym_cifinput`、`smodes`、`comsubs`。
各自完整 ELF 元数据保存在 static 目录。`iso` 保留 Fortran 函数与 COMMON/module
符号，包含 `_gfortran_*` 运行时；不能据此恢复源代码原始类型声明或所有数组边界。
运行横幅确认 ISO 9.6.1（Jan 2022）、SMODES 1.2.5（Sep 2021）；
随附 smodes.txt 的 1.2.3 是文档版本，不能替代实际程序版本。

已定位静态调用：

- `isotropy_display_distortion_` (`0x4260a0`) 调用 `get_irreps_`
  (`0x426d37`)、`project_` (`0x4293a2`)、`project_local_` (`0x42af42`)、
  `project_vector_local_` (`0x42b049`)、`project_vector_` (`0x42c476`)。
  这支持“畸变输出通过表示与投影程序构造”的调用层结论，不单独证明全部投影公式。
- `orderparam_to_subgroup_` (`0x49e2bd`) 调用 `generate_space_group_`
  (`0x49e321`)、`newlat_order_` (`0x49e392`)、`get_new_fractionals_`
  (`0x49e3ee`)、`get_irreps_` (`0x49e505`) 和 `orderparam_check_`
  (`0x49e49c`)。可见晶格、群操作、表示和序参量约束共同参与筛选。
- `orderparam_check_` 的 `0x4b4861..0x4b487c` 是乘加循环，
  `0x4b4813` 减去序参量分量，`0x4b4828` 调用 `ncmp_`；
  非零比较在 `0x4b482f` 跳至失败返回。与逐分量检验 `D(g)η−η=0` 相符。
  这是由汇编推导的局部解释；该例程的通用参数类型和内部容差尚待完整核验。
- `project_` 调用 `transform_irrep2complex_` (`0x49f68a`)、`get_irreps_`
  (`0x49faed`/`0x49fb43`) 和复数行约简 `zrowop2_` (`0x49fd3a` 等)。
  本次未证明其每个归一化步骤与网页模式向量约定相同。

### 新取得的函数级证据

- `orderparam_check_`：
  `ev_532c54c72d30b144fda62477b3ffc34e3b611d5054f43e5e06a71ad360dd6eaa`。
  REA 同时返回伪代码、反汇编、基本块和引用。对参数与容差的独立复核集中在
  [ORDER_PARAMETER_AUDIT](ORDER_PARAMETER_AUDIT.md)，不把单个例程扩展成全空间群证明。
- `get_irreps_`：
  `ev_876f3330452ee37940cf889e2523bf929f49ca0cfd8c7a84ee870ada181f9d86`。
  实际可见 IR/小群索引检查、`common_little_` 和 `common_magnetic_` 数据访问、
  `get_irrep4_` 调用，以及返回标量或矩阵的分支。标量/矩阵模式和额外符号处理可见，
  但未将所有参数编号都恢复成官方命名；尤其不据此宣称所有磁性表示约定已核验。
- `project_`：
  `ev_1131e1fb43cc1bbde64db04c22132ae5b6ee430b157e4ce0589a9f3a02ca8cc5`。
  提取文件 `project.c` 第 214、219 行调用两次表示计算；第 226–240 行的内层循环
  累加两个表示输出的乘积。随后第 246–280 行求该候选的最大绝对分量，
  判非零后第 288 行用该值缩放；第 315、446 行调用 `zrowop2_` 做复数行操作，
  根据 `ncmpz_` 的非零结果保留独立候选。
  这是局部数据流解释，**最大分量缩放不是整个晶体模式的物理归一化证明**，
  也不是“每次投影均采用简单欧氏单位长度”的依据。首行类型传播警告适用。
- `project_vector_`：
  `ev_995c4e8a0963d7a321775c56b99b779feec99547d7d3799555092da781337c30`。
  可见对 Wyckoff 相关共享表的索引、`constant_` 解码及浮点乘加；
  这补足了常数表确实参与向量生成的读取链，尚未恢复所有字段布局。

## 常数表如何解码

`init_constants_` 在 `0x456ed3` 读取 ISODATA，拼接 const.dat，随后进行
Fortran 文本读取。`0x45701f` 读 4-byte 整数索引，`0x45703b` 读 8-byte real
到 `const_common_ + index*8`。在 `0x4570aa..0x4570be`，程序对偶数位置的值
变号并写至紧邻的奇数位置。`constant_` 的 `0x457199..0x4571a3`
按索引取得 double，之前有索引范围检查。

因此 const.dat 中 `2 → 1`、`4 → 1/2` 对应内存中的 `3 → −1`、
`5 → −1/2`；这是机器码与表数据相互支持的解释，不是拟合材料答案。
还不能据此断定所有 data_*.txt 的整数矩阵都使用这套编码：
每个字段的读取和解码路径必须分别追踪。特别不能把 0..15 一律直接当作实数矩阵。

## 隔离执行及反例

`tools/probe_native.py` 的输入先按官方命令语法和通用立方例子固定，然后运行。
`evidence/native/runtime/*/record.json` 保存输入/程序/所有 ISODATA 文件指纹、
命令、平台、时间、退出码和变化清单；stdout/stderr 原文可直接审查。
每次都创建独立 cwd 与数据库副本，不读官网结果。

`iso-gm4minus-p1`：输入 221、GM4−、P1、Wyckoff a，
实际输出子群 99 P4mm，指数 6、胞体积比 1，
基矢 `(0,1,0),(0,0,1),(1,0,0)`，原点 `(0,0,0)`，
点 `(0,0,0)` 的投影矢量 `(1,0,0)`。
这里只确认这个命令实例，不等于完成任意 CIF → 全模式 → 归一化导出的 Method 2。

`iso-invalid-keyword`：实际输出 `Syntax error`，但退出码是 0。
这是“进程成功不等于计算成功”的最小反例，执行器必须保留领域诊断。

`findsym-shipped-input`：首次传长绝对路径，Fortran 实际打开被截断的文件名，
随后 EOF，退出码 2。保留失败 stderr，不把它记为 FINDSYM 验收通过。
调用侧应将输入放在工作目录中，仅传短文件名，并重跑验证。
随后 `findsym-short-filename` 使用短相对路径实际成功，随附 NaCl 输入识别为
225 Fm-3m。新执行器采用该短路径合同。输出主标题 7.1.6 与内嵌 CIF 注释
7.1.3 的差异原样保留，不能偷偷统一版本文本。

`findsym-cif-short-filename` 中未显式提供 symmetry-operator loop 的合成 CIF
报告读取错误，尽管退出码为零。另一个先行固定、含恒等操作的合成 P1 CIF，
经冻结应用成功转换并将未编辑输出交给 FINDSYM：
[正例证据](../evidence/native/cif-converter-positive/README.md)。
它验证该输入的转换/格式合同，不是任意 CIF 兼容性证明。

`comsubs-shipped-input` 成功完成随附的有界公共子群样例，结果保存在原始 stdout。
用户输入的搜索边界是结论的一部分；不能把有限搜索结果理解为所有超胞的穷尽证明。

`smodes-shipped-input`：实际处理随附输入中的 GM、Z、DT 1/4，
输出对称适配模式。GM 的维数合计 `3 + 1 + 8 = 12 = 3×4`，
与常规 P 晶胞的四个原子位移自由度相符。这里只是维数核查，
并未验证全部模式的正交归一性、声子频率或力常数。

## 数学依据与独立核查

以下为本项目推导，不冒称官方源码：
对已指定表示 `D` 和非零序参量 `η`，稳定子定义为
`Hη={g∈G | D(g)η=η}`。因为 `D(gh)η=D(g)D(h)η=η` 且逆也固定 η，
它构成子群。非零极性向量的正负方向不可混同；把 `±η` 一律当同一不变向量会错误扩大子群。
对有限正交表示，`P_H=|H|⁻¹Σ_{h∈H}D(h)` 满足 `P_H²=P_H`，
并且每个 H 不变向量由 P_H 保持。此公式描述不变空间投影；
一般 IR 的矩阵元投影还需要其表示矩阵和完整的适用条件，不能仅由本式替代。

本次最小例子是 Γ 点的极性三维立方表示：群用全部 48 个 3×3
带符号置换矩阵精确枚举。η=(1,0,0) 的稳定子有 8 个操作，
指数 48/8=6；投影为 diag(1,0,0)。`tools/verify_cubic_stabilizer.py`
用整数/有理数核对群闭包、稳定子闭包、投影幂等和向量不变性，
另覆盖 (1,1,0)、(1,1,1)、(1,2,3)、零向量与反演反例。
不使用浮点容差，也不靠官网期望值构造群。
这是对本例指数/方向的独立依据，不推导空间群编号，也不证明任意 k/IR。

位移方向为无量纲基矢系数，未指定振幅和晶格长度就不解释为 Å。
Γ 点本例无非平凡平移相位；不能外推到非零 k、磁性/反幺正、超空间、
多 k 耦合或网站振幅归一化。没有模拟动力学或能量，能量守恒不适用于本检查。

公开命令合同来自 [ISOTROPY 用户手册](https://iso.byu.edu/isotropy_doc.php)
的 Display Isotropy、Display Distortion、Value/Show 命令章节。
本报告审读上述章节，没有声称审读其书目中全部专著或已核验全部理论。
