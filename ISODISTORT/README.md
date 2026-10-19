# ISODISTORT 本地版

ISODISTORT 本地版用于从母相晶体结构搜索畸变子群、计算对称适配模式，并把结果导出给后续分析与可视化程序。它在当前本地网页的功能范围内实现 [BYU ISODISTORT](https://landau3.byu.edu/isodistort.php) 的工作流，提供同一计算引擎的**英文网页**和 **Python 会话 API**。

计算使用 `resources/isobyu/` 中冻结的 Linux ISOTROPY 程序 `iso`、`findsym`、`smodes` 及通用数据库。官网下载的具体晶体结果只作独立审计参考，不能被搜索、模式、拟合、缓存或导出读取来补答案。

先按下面的快速开始运行一份自己的 CIF，再根据任务阅读对应 Method。科学能力受输入条件和实现范围限制，不能由有限测试外推到任意晶体。最新证据见[验证报告](docs/BUGFIX_VALIDATION_REPORT.md)，公式、推导和适用域见[理论参考资料库](docs/reference_library/README.md)。

## 快速开始

以下为本机 Windows 工作区命令。没有部署环境时，先按根 [README 的安装步骤](../README.md#第一次部署windows)准备 Python 3.10+、WSL 默认 Linux 发行版和项目资源。

```powershell
cd E:\CRIS
py -3.10 setup_cris.py install --project isodistort
.\run_cris.ps1 setup_cris.py doctor --project isodistort
.\run_cris.ps1
```

示例中的 `py -3.10` 可统一替换为本机已安装的 Python 3.10+ 命令。已经部署时跳过 `install`。两个运行项目共用根 `.venv`；不要在 `ISODISTORT/` 再建环境。本目录的 `main_requirement.py` 是兼容安装入口，转发到同一个根部署器并选择 ISODISTORT，不另建环境。

`run_cris.ps1` 不传参数时运行本目录的 `main_web.py`；也可显式运行：

```powershell
.\run_cris.ps1 ISODISTORT\main_web.py
```

它从根 `.venv` 读取依赖，并用记录的基础 Python 启动 WSL，避免本机曾出现的 OneDrive/venv 启动限制。无需激活虚拟环境。终端会显示实际 URL 并尝试打开浏览器；通常从 `http://127.0.0.1:8000/` 开始，端口占用时寻找空闲端口，以终端地址为准。

在网页中：

1. **Parent CIF**：选择母相 `.cif`，点 **Load**，确认空间群、晶胞和原子位点。
2. **Types of distortions to be considered**：选择畸变类型与物种，点 **Change**。
3. **Method 1、2 或 3**：填写搜索条件，点 **OK**，查看候选表。Method 4 则另外上传已经畸变的女儿相 CIF 做模式分解。
4. **Distortion**：选择刚才计算的 Method，下载筛选表；Method 1–3 还可勾选格式下载结构 ZIP。

下载文件保存到**浏览器的下载位置或“另存为”所选目录**，不会自动写入本项目的 `output/`。第一次可先下载少量候选，核对状态报告和实际文件，再进行大批量导出。

## 输入与界面操作

### 母相结构

必须先成功加载母相。母相 CIF 描述参考晶格、空间群和原子位置；Method 4 的女儿相 CIF 描述已经畸变的结构。`experiment_data/` 中的实验文件是只读输入，本机可按实际文件选用，源码克隆不保证包含它们。

加载后显示的空间群、晶胞、Wyckoff 位点和位点参数来自当前输入。等价的非标准 basis/origin 会连同晶格、原子及独立物理轨道一起转到 conventional standard setting；无法保持空间群、轨道身份或多重度时明确拒绝。

如果 FINDSYM 失败并显示 `parent_orbit_standardization.status=unavailable`，本次加载只能用于诊断，后续搜索、模式和子群导出（包括 CIF-only）会拒绝。先修复环境，再重新加载母相。

### 畸变类型与物种

当前正式计算范围是 **Strain**（晶格应变）和 **Displacive**（原子位移）。默认选择这两类；Displacive 可选择全部或部分物种。每次改动须点 **Change**，只勾选不会生效。`all / none / 具体物种` 同时出现时按该优先级解释。

**Occupational、Magnetic、Rotational 尚不支持正式计算**，与已支持类型混合请求也会明确拒绝。界面有这些标签不代表相关物理模式已经实现。

物种选择按独立物理 Wyckoff 轨道传给计算；同元素或同一 Wyckoff 字母下的不同轨道不会合并。

### 结果表

Method 1–4 都提供 **Filter**、排序及筛选表下载。`Show filtered rows only` 控制显示；下载仍包含全部筛选命中行。Method 1–3 点选候选后计算该候选模式，模式表显示在 Method 2 区域。

`idx` 只标识所属候选表中的行，不是跨 Method 的全局身份。IR/OPD、空间群、basis、origin、方向及查询上下文共同决定科学身份；目录短名或屏幕第几行都不能代替它。

## 四种 Method 如何选择

### Method 1：搜索全部特殊 k 点

适合还没有指定传播波矢、希望先查看特殊 k 点下可达子群的任务。

- **Crystal system(s)**：晶系之间取 OR，全不选表示不按晶系过滤。
- **Space-group symmetry**：指定子群空间群，留空不限制。
- **Conventional lattice / Primitive lattice**：互斥的格子类筛选，分别使用惯用胞或子群原胞。它们按点群轨道和精确格等价关系分组，不是简单比较矩阵文本。
- **Maximal subgroups only**：只保留极大子群。
- **OK**：执行搜索，显示 `idx / Irrep / OPD / Dir / SG / basis / origin / s / i / k-active`。

同一格子可能选用不同的等价代表基矢。正式对照须检查完整候选身份及精确仿射等价性，不能仅凭 basis 字符串相同或不同作结论。

### Method 2：指定一个 k 点

适合已知主波矢及其参数的任务。先在 **Specify k point** 选择母相允许的 k 点；参数线需填齐相应数值，再点 **OK** 搜索其单 IR 子群。

**Change number of superposed IRs** 当前只支持 `1`。多 IR 叠加需要稳定子的精确交集，不能靠拼接几个单 IR 列表实现；网页请求 `nsup > 1` 会在计算前拒绝。

**# of independent incommensurate modulations** 使用整数 `nmod=0..3`，改动后点旁边 **Change**：

- `0`：三维公度锁定，保留折叠到子胞内的 k，包括允许的次级 IR。
- `1、2、3`：当前 Method 2 只选一个主 q，因此都保留该 q 的谐波及 Γ；提高数值不会额外引入第二条独立调制。

**Generate isotropy subgroups database if missing** 在通用数据库缺少当前参数 k 子群库时，由本机 `iso` 生成。首次可能很慢，超时由配置 `runtime.generation_timeout` 控制；生成的 `i*.iso` 属于当前运行会话，关闭会话后不保留。**Manage cached subgroup databases** 可查看、选择和删除该会话缓存。

这个勾选项只生成搜索所需的通用子群库，与网页范围之外的 Distortion Generate（按幅度生成畸变结构）是不同功能。参数 k 的位移模式使用本地 smodes 与子群对称性/超空间管线；不是默认只导出 CIF。

### Method 3：指定目标空间群和格子

适合已知目标空间群及其与母相晶格关系、需要搜索对应 embedding 的任务。

1. 选择 **space group symmetry** 或 **point group**。
2. 使用 **real-space sublattice（direct）**。Default 在空间群查询中采用目标空间群的默认 centering；仅选点群时采用 P，也可显式选择 P/A/B/C/I/F/R。
3. 填写相对母相的三条基矢 `a' / b' / c'`，可输入分数。基矢连同 centering 必须定义母平移格的真实子格，不能只按体积判断。
4. 点 **OK**。缺少参数 k 子群库时，可使用 Method 2 区域的生成选项。

**primitive reciprocal-space superlattice 输入当前不支持。** direct 查询内部使用精确倒格关系，并不意味着支持这个界面输入模式。

空间群查询枚举 affine embeddings，并在 `strain ⊕ displacive` 表示中检查可达性。具有精确 fixed-space 和稳定子交集见证的 coupled 候选可以选择；它们不声称某一组 IR 是唯一 primary 分解。`known routes` 是已知来源路线，不能把路线条数当作 embedding 数。

当前单 IR 路线搜索域是特殊 k 与一参数公度线；任意/多参数 k 尚未完整覆盖。仅按 point group 查询仍是 single-IR 子集，没有空间群查询的完整 affine 枚举。Occupational/Magnetic/Rotational 也没有同等级可达性表示。搜索规模超过配置上限会明确失败，不以截断结果冒充完整结果；Method 3 的一般完备性仍须按[能力映射](docs/reference_library/coverage.json)中的条件判断。

### Method 4：分解已有畸变结构

先加载母相，再在本面板上传女儿相 CIF：

- **Atom-matching method**：`nearest-site` 按物种做全局最小笛卡尔距离的一一匹配；`robust` 额外拒绝超过阈值的匹配。
- **Robust distance threshold (angstrom)**：物理距离阈值，单位 Å。
- **Known origin shift**：可选的女儿相分数坐标原点平移 `(x,y,z)`。
- **OK**：输出模式幅度、原始拟合系数、`normfactor` 及 residual。

结果是分解幅度表，不是子群列表。可下载 txt/csv，**不能用作 Distortion 结构 ZIP 的来源**。

位移幅度 `As` 按女儿原胞内笛卡尔模式范数定义，`Ap=As/sqrt(s)`，其中 `s` 是对应原胞数比。RMS/max residual 使用 Å。报告小 residual 仍不足以证明模式身份和输入映射正确，必须满足完整模式列、平移字符及子胞周期性合同。

勾选 strain 且有精确目标子群 embedding 时，才报告规范 symmetry-adapted strain 模式。母胞行晶格记为 `P`、子群 basis 为 `B`、女儿胞为 `D`；应变由 metric 方程求得，并以 `B M P` 重建实际晶胞。`q=(E11,E22,E33,2E23,2E13,2E12)` 在非正交母胞中不是 Cartesian 张量分量。`q_raw` 与 `q_unit=normfactor*q_raw` 分别服务不同格式，不能混用；完整定义见[理论资料库](docs/reference_library/README.md)。

缺少规范应变基、非自共轭 k 的完整子胞列或必要映射证据时，程序返回未解析状态或明确拒绝，不从单个 cosine 列猜 sine/OPD。Occupancy、Magnetic、Rotational 分解尚未验收。一般输入域的科学验收状态见[验证报告](docs/BUGFIX_VALIDATION_REPORT.md)。

## Distortion：下载哪些文件

在 **Export source (exactly one Method)** 选择一个已计算的 Method：

- **Download filtered (txt/csv)**：下载该 Method 当前筛选和排序后的表格，Method 1–4 均可用。
- **Download formats**：为结构 ZIP 勾选 CIF、Save interactive distortion、Complete modes details、TOPAS.STR。
- **Download all (ZIP)**：导出 Method 1/2/3 当前筛选命中的全部候选。没有筛选就导出全部，只包含勾选的格式。

每个成功候选最多含四个文件：

```text
subgroup.cif                子群结构和模式元数据，可交给 VESTA
data.isoviz                 IsoVIZ 交互输入
Complete modes details.txt  本地完整模式文本；官网对应材料通常是 HTML
topas.str                   TOPAS 结构输入
```

ZIP 名形如 `isodistort_method2.zip`。Method 1 候选目录形如 `GM1+_P1_SG139/`，Method 2 形如 `LD1_C1/`；这些只是命名示例，不代表本次必有该候选。Method 3 先有 `M3-<摘要>/` 案例目录，再含候选目录。短名冲突时追加完整身份摘要，不能用目录名代替科学身份。

ZIP 从启动时的独立会话快照计算，不扫描 `output/`、`webpage_info/` 或 `output_compare/`。下载 isoviz、modes 或 topas 时可能需要逐候选补算模式；仅 CIF 可跳过模式补算，但仍须通过母相和子群结构门禁。

计算中可请求在当前候选结束后安全取消。部分候选模式身份未解析，或用户取消后尚有候选未执行时，ZIP 会附 `export_candidate_status.json` 和 `.txt`，列明完整身份、状态和原因；必须检查这些报告，不能把 ZIP 下载成功等同于全部候选通过。其他准备或 writer 异常会拒绝整批发布。

原子位移和应变的四种格式共用已验证的数据合同。规范应变缺失时不会退回猜测标签或按晶系拼基；TOPAS 只写重建后的固定晶胞，不提供 strain-mode 精修参数。当前位移导出要求单物种满占位点；混合或部分占位会拒绝。

网页不提供在线晶体查看、Distortion Generate 或 Domains。查看结构可手工打开 `CRIS/launchers/VESTA.lnk`，把 `subgroup.cif` 交给 VESTA；IsoVIZ 的振幅写入和启动见 [ISOVIZ_INPUT README](../ISOVIZ_INPUT/README.md)。默认搜索 `CRIS/launchers/`，也可用 `ISOVIZ` 显式指定。程序分发成功与 GUI 正确识别文件是两个独立验收步骤。

## 配置、缓存和会话生命周期

运行默认值唯一维护在 [`resources/config/settings.yaml`](resources/config/settings.yaml)，修改后重启服务。主要配置包括：

- `isobyu.bin_dir / data_dir`：通用程序与数据库路径；相对路径以 `resources/config/` 为基准。
- `runtime.web_port / web_idle_timeout`：首选端口与浏览器心跳失联后的停服时间。
- `runtime.output_dir / temp_dir`：源码默认分别对应 `ISODISTORT/output/` 和 `output/tmp/`；浏览器 ZIP 不按此目录保存。
- `runtime.timeout / generation_timeout`：普通外部调用与生成子群库的超时。
- `runtime.method3_max_*`：精确参数搜索、查询数、有限商及 coupled 状态规模限制；各项零值含义见配置，不能一概当作无限制。
- `defaults.*`：科学容差及默认参数。实验/CIF 对称性识别使用 Å 制笛卡尔距离；精确仿射数据在库边界使用另一容差；分数坐标、矩阵和波矢残差是无量纲量，不能混用。

Windows 各 wrapper 在 WSL **`/var/tmp/id-XXXXXXXX`** 建立独占短目录。只读通用 ISODATA 与生成时的私有副本分离；事务克隆独立复制已生成资源，互不共享可变缓存。`/var/tmp` 用于跨 WSL 自然空闲退出保留当前 Windows 会话的暂存数据；正常关闭实例会删除自己的目录。无需手工更改系统 `ISODATA`。

外部调用超时后，该运行实例失效。**网页重新上传母相，Python 新建 `IsoDistort`**，不能继续复用旧缓存。运行器对本地子进程等待设有上限，不会无限等待输出管道关闭；它不重启 WSL，也不承诺所有 Linux 后代都已终止或整个会话清理只耗一秒。

网页是本机单用户共享会话，多标签页共享同一状态。变更按 revision 提交，旧标签的计算/选择/导出可能提示 stale context，此时刷新状态后重试。上传新母相建立新运行实例，加载成功后才替换旧状态；失败不覆盖当前结果。只读查询和缓存管理不递增 revision，但外部超时仍可能使实例失效。

导出使用独立快照及任务 ID；取消只作用于对应任务。已开始的 ZIP 不混入后来上传的新结构。Python 脚本必须在结束时调用 `close()` 释放资源。

网页右上角 **Stop** 停服；关闭最后一个标签页后会自动停服，浏览器崩溃或关闭信标丢失时由心跳超时兜底。计算期间关闭页面可能需等待当前请求结束。要保持服务可用，应保留至少一个标签页。

源码以外的 wheel 安装提供 `isodistort-web` 命令，输出使用用户运行目录而非安装包目录：Windows 为 `%LOCALAPPDATA%/ISODISTORT`，Linux 为 `$XDG_STATE_HOME/isodistort`（缺省 `~/.local/state/isodistort`）。可用绝对路径环境变量 `ISODISTORT_RUNTIME_ROOT` 覆盖相对运行目录；YAML 中显式绝对路径优先。wheel 不包含根部署器、实验输入、官网归档或手工审计脚本。自带程序是 Linux 版，不能据此宣称支持任意操作系统直接运行。

## Python API

公共会话入口是 `application.session.IsoDistort`。在 CRIS 根保存自己的脚本，例如 `search_parent.py`，用下面的方式运行，确保使用同一依赖环境及 WSL 入口：

```powershell
.\run_cris.ps1 .\search_parent.py
```

脚本示例只搜索并打印候选，不替用户决定具体候选。将输入路径改为实际母相 CIF：

```python
from pathlib import Path
from application.session import IsoDistort

parent = Path(r"E:\CRIS\experiment_data\EuAl4 Parent.cif")
iso = IsoDistort()
try:
    iso.load_structure(parent)
    iso.set_distortion_types(["strain", "displacive"])
    results = iso.search_method_1()
    for item in results:
        print(item.subgroup)
finally:
    iso.close()
```

需要从候选继续计算时，使用 `search_method_2(chosen_index, candidates=...)` 显式传入所属候选池；导出使用 `export_subgroups(destination, subgroups=..., formats=[...])`。多候选导出模式格式时，可显式设 `compute_missing_modes=True` 补算缺失模式。多个 Method 的索引可重复，不能把索引跨池使用。磁盘导出在完成渲染后整批发布；已有非空目标会使用带摘要的 `.isodistort-batch-v1-<digest>.ready` 目录和文件 manifest，避免覆盖既有结果。

外部程序通过 `iso.displacement_snapshot()` 获取只读结果，包括完整候选、计算上下文、原子/物种顺序、Å 制行晶格、分数坐标及有序原始模式列。快照是**应用子原点之前的基超胞**，不能当作已移动原点的 CIF 坐标。不要读取或修改会话私有缓存。

GD 使用此公开合同，以及 `backend.export.distortion_formats.cart_normalized_mode_matrix`。后者要求有限非零模式列、有限可逆的 3×3 晶格与正整数 `N`（不接受布尔值），并明确使用数值 Å 约定：对返回列 `V` 和行晶格 `A`，满足 `Σ|V A|²/N=1`，最大位移与对称幅度上界的乘积为 1。任意给定 `N` 并不自动具有官网 `As` 的含义。必要分量在缩放、点积或输出中下溢时会拒绝；统一极大/极小尺度的良态输入可计算，不保证任意病态晶格的前向精度。推导见[架构合同](docs/reference_library/notes/ARCHITECTURE_CONTRACTS_20261019.md)。

`generate_distortion / generate_mixed_distortion / generate_domains` 保留为脚本和验证 API；它们不扩大当前网页产品范围。旧 `features`、`backend.api` 和旧 Method 1 总引擎入口不提供转发，消费者应使用当前公开包。

## 科学使用边界

请把下面的限制与各 Method 的输入说明一起使用：

- **诊断不等于可用结果**：`selectable=false`、模式身份 `unresolved`、来源缺失或规范化失败的候选，不能作为已验证结果导出。缺失不以空表或猜测成功替代。
- **位移身份必须完整**：目标 embedding 的精确方向、物理轨道、k、原子顺序和唯一来源列需一致。稀疏 microscopic 表只有在公共点域变基唯一且完整 BUSH 延拓通过秩、条件数和残差核验后才能使用。不能用同名 IR 或近邻复制补齐未知列。
- **结构与模式要在同一坐标系**：导出核验父子原子映射 `x_parent=x_child B+q`、物种、物理轨道及平移陪集，保持体积比、原子数和化学计量关系。未知 Wyckoff 字母写 `?`，不由重数猜测。
- **格式不同不必然是结构不同**：原子顺序、周期像、等价 basis/origin 和对称操作顺序可能不同；也不能只因晶格一样就认定同一候选。Method 1 审计还需区分不同 orientation/domain。
- **模式类型与数值域有限**：当前正式类型、nmod、Method 2 单 IR、Method 3 搜索域、Method 4 映射及部分占位限制见前文。数值近似的适用域不能靠放宽容差扩大。
- **理论与实测分开**：通用算法的推导、实现门禁、有限回归、官网对照和外部 GUI 验收提供不同证据。任何单例、历史报告或静态解析都不能替代当前源码和输入域的科学验证。

## 集成验证与测试

`tests/validation/` 是 pytest 和手工审计共用的只读比较实现，生产代码不能导入它。比较前应固定独立输入与参数来源，先运行并冻结本地输出，再读取官网参考；不得从官网结果反推、挑选或补齐计算输入。

从 CRIS 根目录运行自动测试：

```powershell
.\run_cris.ps1 -m pytest ISODISTORT\tests -q --tb=line
.\run_cris.ps1 -m pytest tests_repo -q --tb=line
```

需先用根安装器安装 `--dev` 依赖；真实外部计算还需 WSL 和完整 ISOTROPY。测试按产品层/模块组织，完整路径和单模块命令见 [tests/README.md](tests/README.md)；手工审计工具直接位于 `tests/` 根，其命令和长时条件统一见 [MANUAL_VALIDATION.md](docs/MANUAL_VALIDATION.md)。测试运行成功也只覆盖实际执行的项目；跳过或缺材料的验证不算通过。

比较已有的一对 CIF 或两个目录：

```powershell
.\run_cris.ps1 -m tests.validation compare --local "本地输出.cif" --reference "参考原件.cif" --json
.\run_cris.ps1 -m tests.validation batch --local-dir "本地输出目录" --reference-dir "参考目录" --json
```

路径请替换为实际文件/目录。显式指定路径不会移动参考原件。省略目录时使用 [比较器配置](tests/validation/config/settings.yaml)中 `output/validation/compare/{item,true}` 的默认目录；批量按相同相对路径配对，不生成、重命名或改写输入。

默认比较既有设置下的晶格、周期分数坐标、物种、占据率、磁矩和空间群语义。`--ignore-atom-order` 采用完整二分匹配；任意 setting/origin 等价性需要另做精确仿射审计，不能只用这个比较器的失败否定结构等价。

`--strict` 还要求可信参考 SHA-256 及字节一致。单对用 `--reference-sha256`，批量用 `--hash-manifest` 提供独立登记指纹。语义相同但缺可信指纹是 `inconclusive`；已有差异仍是 `fail`。退出码 0 为已执行比较通过，1 为差异或 inconclusive，2 为输入/解析错误或空批次。普通比较通过不能声明参考来源已验证。

官网文件保存、HTML 序号、已下载项目与补缺顺序仅在 [DOWNLOAD_CHECKLIST.md](docs/DOWNLOAD_CHECKLIST.md)及其 manifest 维护，本 README 不复制清单。已有官网归档、原始实验输入和 `resources/isobyu/` 保持只读。

## 目录与后续阅读

```text
ISODISTORT/
  main_web.py              源码网页入口
  main_requirement.py      转发到统一安装器
  pyproject.toml           包元数据、打包与工具配置
  requirements*.txt        运行及开发依赖
  backend/                 科学算法、共享模型、格式、表和运行设施
    method1/ ... method4/  各 Method 自身的职责
    wrappers/ models/ tables/ utils/ messages/ input_cif/ export/
  application/
    session/               IsoDistort、结果合同与发布
    workflows/             用户动作、事务、上传与导出任务
  frontend/web/            HTTP、英文页面、浏览器与标签页生命周期
  tests/                   按层/模块归组的测试及根级手工审计脚本
    backend/               后端各模块的 test_<模块>.py
    application/           会话与工作流检查
    frontend/              网页检查
    validation/            共用只读比较与来源核查支持包
    fixtures/              测试输入
  resources/
    config/settings.yaml   运行默认值
    isobyu/                只读通用二进制与数据库
  docs/                    计划、报告、理论、清单及归档
  output/                  本地运行与验证产物
```

调用方向固定为 **frontend → application → backend**。网页经 workflows 调用 session，Python 消费者使用同一 session；Method 4 分解由 `backend/method4/decomposition.py` 负责。共享格、仿射运算及模型各有唯一实现，不通过旧入口绕过分层。

- 安装、离线依赖、跨项目目录：[CRIS README](../README.md)。
- 运行能力和科学缺口：[BUGFIX_VALIDATION_REPORT](docs/BUGFIX_VALIDATION_REPORT.md)。
- 当前与后续待办：[DEVELOPMENT_PLAN](docs/DEVELOPMENT_PLAN.md)。
- 长时审计：[MANUAL_VALIDATION](docs/MANUAL_VALIDATION.md)。
- 文献、证明、能力映射：[理论参考资料库](docs/reference_library/README.md)。
- 模块边界和开发验收：[agent.md](agent.md)、[根 AGENTS.md](../AGENTS.md)。
