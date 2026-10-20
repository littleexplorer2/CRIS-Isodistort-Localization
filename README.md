# CRIS REA 本地晶体学项目

本项目提供 **母相 CIF + IR / OPD / 必要的 k 参数 → 子群 → 位移与应变子空间 → CIF
导出** 的本地流程，并保留官方 ISOTROPY 程序的原始输入控制台。计算使用已获许可的
官方 Linux 程序及其通用数据库，再由本项目独立实现结构处理和子群不变子空间算法。
没有复制旧 CRIS 的科学实现，官网案例答案不进入计算链。

当前模式输出是给定子群 **H 允许的完整位移与均匀对称应变空间**，采用本项目的正交
归一基。基向量可能混合母相的不同不可约表示，尚未恢复官网逐 IR 的列分类、列顺序和
归一化约定。因此它提供有明确适用域的本地结构参数化，不能视为全部 ISODISTORT
Method 的等价复现。具体理论和验收边界分别见
[实际算法推导](docs/reference_library/LOCAL_MODE_DERIVATION.md)、
[模式理论合同](docs/reference_library/MODE_CONTRACT.md)和
[本轮验证报告](docs/PHASE2_VALIDATION.md)。

逆向结论与原始证据从 [调查总览](docs/INVESTIGATION.md) 进入；后续工作见
[开发计划](docs/DEVELOPMENT_PLAN.md)。

## 安装与启动

需要 Python 3.10 或以上及 `requirements.txt` 中的 NumPy。Windows 另需可运行的
WSL Linux 环境，其中有 `python3` 及官方程序所需的 Linux 动态库；Linux 可直接运行。
Python 科学计算运行在启动项目所用的解释器中，Windows 主机运行时不要求再在 WSL
安装 NumPy。

优先使用已有、已具备依赖的 Python。本轮实际使用
`C:\Users\devou\AppData\Local\Programs\Python\Python312\python.exe`
及 NumPy 2.5.3；WSL 只负责官方进程。**不要修改或重建旧 CRIS 环境**。
以下命令中的 `python` 均指你选定的解释器：

```powershell
Set-Location 'C:\Users\devou\OneDrive\Desktop\CRIS_REA'
python -c "import sys, numpy; print(sys.executable); print(numpy.__version__)"
```

如果所选环境缺少依赖，在你管理的项目环境中安装：

```powershell
python -m pip install -r requirements.txt
```

可以为新项目使用自己的环境；不需要修改旧项目源码、数据库或重新创建其环境。
运行前确认 `vendor/isobyu/` 已包含许可允许使用的程序和数据库。版本库不携带官方
发行包及每次运行的数据库副本；外部文献 PDF 仅保留本地缓存，版本库跟踪来源 URL、
指纹和本项目分析。新克隆可从已获许可的发行包恢复：

```powershell
python tools/import_vendor.py 'E:\CRIS\ISODISTORT\resources\isobyu'
```

将参数换成实际发行包目录。导入工具依据
`evidence/environment/source_manifest.json` 核验每个文件的 SHA-256；缺失、版本不同
或目标已有不同内容时拒绝继续，不修改来源，不覆盖不同版本。更换官方版本须另行审查
并更新来源记录。

启动本地网页：

```powershell
python -m cris_rea doctor
python -m cris_rea serve --port 8765
```

浏览器打开 <http://127.0.0.1:8765>，停止服务按 `Ctrl+C`。仅监听 `127.0.0.1`，不作为
局域网或公网服务直接发布。`doctor` 检查资源与 Linux/WSL Python 启动能力，不执行
科学计算，也不保证所有二进制动态库都已齐备。环境缺失或权限不足会返回真实错误。
若 Codex 沙箱拒绝启动 WSL，请在允许执行 WSL 的本机终端运行同一命令。

本目录有独立 Git 历史，远程为 `littleexplorer2/CRIS-Isodistort-Localization`，
专用分支为 `codex/cris-rea`。其他机器可获取：

```powershell
git clone --branch codex/cris-rea --single-branch https://github.com/littleexplorer2/CRIS-Isodistort-Localization.git CRIS_REA
```

只在新目录内同步成果，不把整个目录覆盖回旧 CRIS。

## 母相 CIF 到本地模式

网页上方的 **母相 CIF → 本地模式** 是主入口：

1. 粘贴独立母相 CIF，或上传文件。上传保留原始字节，编辑预览后按 UTF-8 文本提交。
2. 填写该母相对应的 IR 和 OPD 标签。需要自由 k 参数的 IR，在 **ISO 独立 k 参数**
   中按官方顺序以空格分隔填写，例如 `1/3`；这里接收独立参数，不是任意完整 k 三分量。
3. 首次将位移和应变振幅留空，点击 **计算本地模式**，生成零振幅参考结构及模式基。
4. 查看 `modes.json` 的 `mode_names`、基列和归一化。需要畸变时，保留相同输入参数，
   按导出列顺序填写整组振幅再运行。位移振幅单位为 Å，应变振幅无量纲；数量须与对应
   基维数一致。位移列中包含单独命名的刚体平移列，不要额外拼接或删掉三个自由度。
5. 结果显示母相、目标子群、原子数与基维数；下载五个导出文件并阅读 `report.json`。

主流程支持：

- 一个 ASCII CIF 1.1 数据块，普通三维、非磁性、公度结构。
- 六个晶格参数、分数坐标、原子标签和类型；允许 `5.4321(8)` 这类含 ESD 的数值。
- 显式 `_space_group_symop_operation_xyz` 或兼容旧标签的对称操作 loop，以及母相
  国际空间群编号。操作须完整、闭合并保持晶格度量，且与官方 ISO
  `INTERNATIONAL NEW` 的 conventional 标准设定一致。
- 固定组成、完整占位；没有给占位率时按 CIF 默认值 1。提供的 Wyckoff 字母会保留，
  提供的站点重数须与操作展开一致；当前不独立识别或证明 Wyckoff 字母本身。

不提供对称操作的 CIF 不会仅凭空间群编号猜操作；本流程也不自动将非标准 setting、
近似晶格或特殊位置理想化。分数坐标和操作按输入值精确解析，小数 `0.333333` 不会
静默变成 `1/3`，含测量误差的位置也不自动吸附到特殊位置。这可能导致群闭包、
重数或原子映射诊断；需要有来源的独立输入修正或另行恢复，而不是读取官网答案补齐。
解析器本身可保留单组分部分占位，但模式流程明确拒绝它。磁性、调制、非公度、混合
占位/无序、多主 IR 和任意刚体约束不在本流程支持域内。详见
[结构合同](docs/reference_library/STRUCTURE_CONTRACT.md)。

当前有限胞实现限制导出胞最多 256 个原子、群闭包最多 384 个操作，超胞平移枚举
最多 1,000,000 次候选检查；超限明确拒绝，不截断后当作完整结果。

CLI 使用同一应用流程。以下 IR/OPD 只是语法示例，须替换为你的母相所支持的选择：

```powershell
python -m cris_rea modes --cif .\my-parent.cif --irrep GM4- --opd P1 --timeout 60
```

附带的自构三斜输入可直接运行，也可上传到网页；它不是实验数据或预计算答案：

```powershell
python -m cris_rea modes --cif examples/triclinic_parent.cif --irrep GP1GQ1 --opd C1 --k-values 1/3 1/4 0
```

有独立 k 参数时添加 `--k-values 1/3`；多个参数可以空格分隔，含负分数时使用
`--k-values=1/3,-1/4,0` 这种带等号的逗号格式，单个负分数也可用 `--k-values=-1/3`。
顺序和分数原样传递给同一协议校验；空的逗号项会拒绝。可选 `--amplitudes` 和
`--strain-amplitudes` 后各接按列排列的有限数值列表。留空振幅表示全零，不从文件名、
标签或其他运行推断。`--timeout` 是每次官方调用的上限，工作流最多调用 ISO 五次，
总耗时可以超过此值。

模式工作流只在检查完成时返回 `validated_local_subspace` 和 CLI 退出码 0。
拒绝或无法判定会返回非零，并在已建立的运行目录保存诊断；这与下方原始程序的
`completed_uninterpreted` 含义不同。

## 模式约定与输出

`basis_convention` 为 `CRIS_REA_H_invariant_Cartesian_orthonormal_v1`。
晶格矩阵以 Cartesian 列向量表示，长度单位 Å；原子分数坐标位于导出的子群
conventional 晶胞中。位移基在该胞内全部原子上的 Cartesian 平方和为 1，不按质量
加权。正交列系数的欧氏范数等于未应变胞中的总位移 RSS，不能直接套用官网模式
振幅。内部位移与允许的刚体平移分开命名，但均属于导出位移基。

应变采用无量纲 Kelvin/Frobenius 正交归一基，顺序为
`xx, yy, zz, sqrt(2)*yz, sqrt(2)*xz, sqrt(2)*xy`；施加 `F = I + E`，要求 F 正定。
它不是声子本征模，不计算频率、力常数或相变能量。所选振幅保持目标子群 H；零振幅
或特殊组合可能有更高对称，当前不识别其最大空间群。官方兼容 IR 列表只作注释，
不作为本地基系数，也不证明本地列已按母相 IR 分类。

每次模式工作流在 `runs/<workflow-UUID>/` 保存原始 `parent.cif`、`workflow.json`、
每次 ISO 查询及原始响应，并在 `exports/` 保存：

- `modes.json`：精确参考坐标、H 操作、基变换/原点、位移与应变基、列名称、选定振幅、
  归一化、检查结果和来源。
- `modes.csv`：每个模式、每个原子的 Cartesian 位移基系数及单位约定。
- `reference.cif`：未施加振幅和应变的参考结构。
- `distorted.cif`：施加所选位移和应变的结构。
- `report.json`：验证残差、适用范围、限制与 CIF 序列化记录。

CIF 导出先尝试 H 的不对称单元表示，再实际重新解析、展开，核对原子数、物种、
全胞位置和晶格度量。如果十进制输出破坏特殊位置的精确约束（例如三重轴上的
`1/3, 2/3`），则改为 **P1 全胞原子表示** 并再次核查，不放宽解析容差或修改坐标。
这是序列化表示选择，**不是把求得的子群改成 P1**。`modes.json` 保留真实 H 和精确
参考坐标，CIF 的 `_cris_rea_preserved_subgroup_number` 标明所保持的 H；
`report.json.cif_serialization` 记录表示、原因和舍入残差。把这种 P1 文件作为新的
母相输入时，主流程会按它声明的 P1 处理，不自动从自定义标签恢复原母相。

## 官方原始工具控制台

网页下方 **Raw official-tool console** 和 CLI `run` 保留官方原始输入：

- `iso`：官方 ISOTROPY 命令语言，通过标准输入传递。
- `findsym`：按 `vendor/isobyu/findsym.txt`，以文件参数传入关键词或传统格式，
  识别结构空间群或磁空间群。其原始工具能力不等于上方模式流程支持磁结构。
- `findsym_cifinput`：把 CIF 文件转换为 FINDSYM 输入；原始控制台不自动续跑 FINDSYM，
  可将标准输出作为下一次 FINDSYM 输入。
- `smodes`：按 `vendor/isobyu/smodes.txt`，通过标准输入计算位移对称模式。
- `comsubs`：按 `vendor/isobyu/comsubs.txt` 搜索公共子群与结构路径；范围受用户输入和
  官方默认值限制，超时不能解释为没有子群或路径。

原始工具包装层保留输入字节，不改变晶格、坐标、容差、单位、归一化或 k 向量；网页
编辑后按 UTF-8 发送。不同程序的科学约定以对应手册为准。例如运行发行包示例：

```powershell
python -m cris_rea run findsym --input vendor/isobyu/findsym_sample.in --timeout 30
python -m cris_rea run smodes --input vendor/isobyu/smodes_sample.in --timeout 60
```

先预览受限 ISO 查询，再按需加 `--run`：

```powershell
python -m cris_rea query --parent 221 --irrep GM4- --opd P1 --wyckoff a
python -m cris_rea query --parent 221 --irrep GM4- --opd P1 --wyckoff a --run
```

参数由用户指定，没有该组合的预存答案。网页 **Prepare a restricted ISO query**
提供相同功能；此原始命令辅助表单不处理 k 自由参数或母相 CIF，不是上方模式流程。

原始工具的 `completed_uninterpreted` 和 CLI 退出码 0 只表示进程正常退出，仍须
检查完整输出。有限的 `Syntax error`、行首 `Error`、`This program has bombed`、
Fortran runtime error 等诊断会标为 `program_diagnostic` 并返回非零；其余输出并不
因此自动获得科学验收。

## 运行证据与边界

每次底层官方调用都另有 `runs/<run-UUID>/`，保存 `input.dat`、`stdout.txt`、
`stderr.txt`、`execution.json`、`process.json`、`evidence.json`，以及独立的 `work/`
和可写 `isodata/` 副本。指纹覆盖输入、程序、通用数据库、源码和产物，并保留数据库
运行前后状态。模式工作流的 `binary_runs` 链接这些底层记录。

计算只读取 `vendor/isobyu/` 的官方通用资源和独立输入，不读取 `evidence/` 中官网
结果、旧 CRIS 输出或参考答案。数据库副本保留以便审计，每次运行会占用磁盘；确认
不再需要后，可手动删除对应的完整 UUID 目录。不要只删模式记录引用的底层目录而
仍宣称证据链完整。

目录隔离不是操作系统沙箱。官方命令语言可能包含文件操作，仅执行可信输入。服务
使用本机绑定、Host/Origin 检查、会话令牌和同一进程内按项目共享的串行锁；不提供多用户权限或
文件系统隔离。下载只允许工作流导出的五个文件并核对已保存的指纹。

单个输入上限 1 MiB，每次官方调用超时上限 300 秒，标准输出与标准错误合计上限
20 MiB；这些是工程保护，不是科学筛选参数。官方程序自行创建的其他文件不受该输出
字节上限约束。

## 工程结构与检查

调用方向为 `frontend → application → backend`：

- `cris_rea/frontend/`：CLI、HTTP、网页输入与结果展示。
- `cris_rea/application/`：模式工作流、共享串行锁、运行记录和下载协调。
- `cris_rea/backend/structure/`：独立 CIF 解析、晶格、精确群作用和轨道展开。
- `cris_rea/backend/modes/`：公度超胞、H 不变空间、畸变和导出。
- `cris_rea/backend/iso_protocol.py`：官方命令与输出的有限、明确解析合同。
- `tests/`：结构/数学反例、协议与前端工程检查；测试替身不用于真实科学验收。

新项目检查命令：

```powershell
python -m unittest discover -s tests -v
```

这些检查不调用旧 CRIS 的全量回归。具体已执行命令、真实官方程序运行和有限验收范围
以 [本地验证报告](docs/LOCAL_VALIDATION.md) 为准；测试通过、原始进程退出或页面展示
均不能替代 ISODISTORT 全功能的独立科学验收。
