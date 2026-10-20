# 本地调用合同与科学影响核查

## 问题、责任和理论边界

原始程序入口解决的是：把用户明确提供的原始输入交给已获许可的官方 Linux 程序，
在独立目录保存输入、输出、进程状态与资源指纹。该原始执行器没有实现或推导空间群求解、
表示构造、序参量方向分类、模式归一化、拟合或子群枚举算法。

原始包装层不存在新增数值近似、拟合参数、科学容差或经验公式；不以网站输出
反推输入，不按材料或答案选择算法。进程完成和科学正确性是两个独立命题，
`completed_uninterpreted` 明确只记录前者，`scientific_validation` 始终为
`not_performed`。科学结论必须由独立的后续研究与验证提供。

结构化 `modes` 入口在这一执行器之上增加了独立科学合同，见下节。
本文件其余“原始调用”叙述不覆盖结构化求解器的坐标、归一化或验收状态。

## 结构化求解合同

调用链为 CLI/HTTP → `application.modes.ModeWorkflow` →
`structure`/`iso_protocol`/`runtime`/`modes`。与原始执行器共享项目级串行锁。
用户输入是原始母相 CIF、单个 IR/OPD、该 IR 符号 k 的独立有理参数，
以及可选本地基振幅。CIF 必须显式提供普通三维对称操作和空间群号。
不自动猜测 setting、不从官网候选结果补参数。

官方 ISO 提供通用 IR 的 k 参数维数、母/子群点余类与原胞基、OPD 子群操作、
基变换与原点、群指数和原胞体积比。所有调用均在独立数据库副本内执行。
参数个数由实际官方符号 k 表达式检查，缺省和多余参数均拒绝。
名义 k 代表不自动等于 OPD 选中的活跃星臂；当前不声称独立复核了完整 IR 矩阵及相位。
`DISPLAY DIRECTION` 相容表示文本作为注释保存，三位小数系数不参与本地模式计算。

生产模式由 H 在所声明有限胞原子上的极矢量作用构造，完整性、位移/应变度量、
坐标及心化约定、数值界和反例见
`reference_library/LOCAL_MODE_DERIVATION.md`。基为本地 H 不变笛卡尔正交基，
可混合母相 IR，不等同官网逐 IR 模式列。状态 `validated_local_subspace`
只表示这些局部合同检查通过，不代表所有空间群、官网逐列或动力学验收。

每次保留 parent.cif、请求、所有 ISO 原始命令/输出与运行 ID、代码/资源/输出 SHA、
Python/NumPy 版本、矩阵检验残差、导出文件。失败保存
`rejected_or_inconclusive`；失败流程禁止通过下载接口提供正式导出。
JSON/CSV 包含完整模式，CIF 在写出后重新解析检查；无法精确保留特殊位点的十进制
H 非对称单元时使用明确标记的 P1 全胞序列化，求得的 H 仍保存在 JSON 及自定义标签中。
这不是把求得的子群改成 P1。详见导出报告的 `cif_serialization`。

## 原始调用的输入到输出对应关系

1. CLI 文件或网页上传得到字节序列；网页直接输入/编辑的文本编码为 UTF-8。
2. 校验程序白名单、非空、1 MiB 上限、不含 NUL 和有限正超时，原样写入 `input.dat`。
   不解析或改写原子、坐标、k、IR、OPD、晶格设定或容差。
3. 从 `vendor/isobyu/` 复制所选可执行文件、`const.dat` 和 `data_*.txt`，逐项核对
   复制前后 SHA-256。既有 `.log`、参考输出和样例输入不自动复制到计算资源。
4. 在新 `work/` 运行，并把 `ISODATA` 指向当前可写副本，末尾带 `/`；所有命令
   使用参数数组，不经 shell。输入 argv 使用 `../input.dat`，`ISODATA` 使用
   `../isodata/`，避免旧 Fortran 字符缓冲区截断长绝对路径。
5. 收集 stdout、stderr、程序自行产生的文件和数据库运行后指纹，保存真实退出状态。
   不把输出内容补写到输入，也不自动宣称计算接受了输入。

FINDSYM 手册第 “Input using key words” 节要求 `findsym inputfilename`；
同手册开头给出 `findsym_cifinput ciffilename` 生成 FINDSYM 输入的合同。
SMODES 手册开头和 “input” 节、COMSUBS 手册开头明确采用标准输入。
ISO 的受限命令生成入口使用 `VALUE PARENT`、`VALUE IRREP`、`VALUE DIRECTION`、
`SHOW`、`DISPLAY ISOTROPY`，可选 `VALUE WYCKOFF`、`SHOW MICROSCOPIC VECTOR`、
`DISPLAY DISTORTION`。这些是官方命令序列的序列化，参数直接来自用户，原始命令
在运行前可检查，实际适用性由官方程序与独立科学核查确定。实现不保存组合答案，
不自动推断 k 参数、Wyckoff 坐标、晶格参数或完整 CIF 模式导出。
上述陈述仅核对本地发行包手册，不声称完成底层算法证明。

## 原始调用的坐标、单位、规范化与适用条件

不增加统一坐标约定，沿用每个工具及其输入格式的官方约定：

- FINDSYM 关键词格式：晶格长度及位置容差按手册为 Å，晶格角度为度；原子位置是
  以 conventional lattice 基为单位的无量纲坐标；磁矩单位为 Bohr magneton。
  关键词格式和传统格式的默认容差存在差异，包装层不覆盖默认值。
- SMODES：输入晶格角度为度，Wyckoff 参数按照手册指定设定；输出 k 坐标按 conventional
  reciprocal lattice，模式位移使用手册规定的笛卡尔坐标。程序并不保证所有 Wyckoff
  参数约定都与表格一致，手册要求查看实际输出位置。
- COMSUBS：长度单位任意但必须自洽，角度为度；单斜 unique axis b / cell choice 1、
  菱方用六方轴、origin choice 2。边界参数属于用户的原始模型输入，包装层不调整。
- ISO：由用户的官方命令及对应手册决定；本层不自动解释或转换 k/IR/OPD。

包装层的 byte-preservation 不变量可以用 SHA-256 检查。空间群闭包、模式正交性、
维数、指数重数、正逆变换等不变量不属于当前工程包装层的计算结果，本次没有执行
这些科学验收，不能以工程测试替代。

## 运行限制与失败语义

Linux 原生调用；Windows 经 WSL `python3` 启动同一 Linux supervisor。缺少 WSL、
Python、动态库或权限会保存错误，不使用原生 Windows 模拟计算回退。
Linux 子进程用独立进程组启动，超时和输出上限触发进程组终止；若外层 supervisor
本身失联，标记 `supervisor_timeout`，明确无法确认子进程终止。

`program_diagnostic`、`process_failed`、`timed_out`、`output_limit`、`runtime_error` 和 supervisor 失败都
不能被解释为科学搜索无解。`completed_uninterpreted` 也不能被解释为科学正确；
Fortran 程序可在退出码零时打印错误。stdout/stderr 超过 20 MiB 或超时会截断计算，
因此这些运行不可用于宣称搜索完备。

输出诊断仅识别显式 `Syntax error`、行首 `Error`、`This program has bombed`、
Fortran runtime error，不试图推断任意文本的科学语义。未命中这些有限模式不是成功证明。

独立 cwd 和 ISODATA 副本只隔离默认工作文件。并没有 chroot、容器或系统级文件权限
沙箱，用户必须信任原始命令输入。禁止把该服务器开放给不可信用户。

## 验证范围

工程测试覆盖原始输入传递、拒绝未知程序/路径、大小与超时边界、结果状态传递、
独立运行目录、CLI 和 HTTP 的会话边界。测试中使用替身进程/假资源仅验证控制流，
明确不是官方科学样本，也不能写成真实二进制验收证据。

真实官方二进制运行和网站交互证据另存 `evidence/`，比较过程只能读取冻结的本地结果。
结构化求解的实际科学验证入口为 `docs/PHASE2_VALIDATION.md`。
