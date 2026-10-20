# CRIS REA 本地晶体学项目

这是从空目录建立的本地项目，用标准库 Python 调用已获许可的官方 Linux
ISOTROPY 程序，并保存可追溯运行记录。没有复制旧 CRIS 的科学实现，也没有新增猜测性
晶体学算法。

课题组当前优先目标是 **母相 + k / IR / OPD → 子群与畸变模式**。当前版本提供
官方工具的原始输入执行入口；母相结构处理、模式构造及网站各 Method 的完整本地
编排仍需进一步研究和验证，不能把这个控制台视为已经复现 ISODISTORT 网站。

研究结论与原始证据从 [调查总览](docs/INVESTIGATION.md) 进入；
后续实现顺序见 [开发计划](docs/DEVELOPMENT_PLAN.md)。

## 启动

需要 Python 3.10 或以上。Windows 另需可运行的 WSL Linux 环境，其中有 `python3`
及官方程序所需的 Linux 动态库；Linux 可直接运行。项目无 pip 依赖，不创建环境。
运行前确认 `vendor/isobyu/` 已包含许可允许使用的程序和数据库。

版本库保存源码、说明与证据清单，不携带官方发行包及每次运行的数据库副本。
在新的克隆目录中，可以从你已获许可的官方发行包恢复 `vendor/isobyu/`：

```powershell
python tools/import_vendor.py 'E:\CRIS\ISODISTORT\resources\isobyu'
```

将参数换成实际的发行包目录。导入工具先按
`evidence/environment/source_manifest.json` 核验每个文件的 SHA-256，遇到缺失、
版本不同或目标已有不同内容时拒绝继续；不修改来源目录，也不覆盖不同版本的目标。
更换官方版本需要另行审查和更新来源记录，不能跳过校验。

本目录使用独立 Git 历史和专用分支跟踪，便于从空项目审查成果；它不是旧 CRIS
工作树的覆盖式更新。请在本目录的分支内同步，不把整个新目录复制回旧 CRIS。
远程为 `littleexplorer2/CRIS-Isodistort-Localization`，专用分支为
`codex/cris-rea`。在其他机器获取这份独立项目时使用：

```powershell
git clone --branch codex/cris-rea --single-branch https://github.com/littleexplorer2/CRIS-Isodistort-Localization.git CRIS_REA
```

在 PowerShell 中：

```powershell
Set-Location 'C:\Users\devou\OneDrive\Desktop\CRIS_REA'
python -m cris_rea doctor
python -m cris_rea serve --port 8765
```

浏览器打开 <http://127.0.0.1:8765>。选择程序，粘贴原始输入或上传文件，设置超时，
点击 **Run official program**。界面为英文。停止服务按 `Ctrl+C`。
仅监听 `127.0.0.1`，不能作为局域网或公网服务直接发布。

`doctor` 检查文件是否存在、Linux/WSL Python 是否可以启动；它不执行科学计算，
也不等于确认每个二进制所需的动态库均已安装。WSL 不可用、权限不足或依赖缺失会如实
返回错误，不会用模拟结果替代。若 Codex 沙箱拒绝启动 WSL，需要在允许执行 WSL 的
本机终端运行同一条命令。

## 可调用程序与输入

- `iso`：官方 ISOTROPY 命令语言，经标准输入传递；界面可按用户指定的空间群编号、
  IR/OPD 标签和可选 Wyckoff 字母生成受限查询。它不推断 k 参数、原子坐标或标签合法性。
- `findsym`：按 `vendor/isobyu/findsym.txt`，以文件参数传入关键词或传统格式，
  识别结构空间群或磁空间群。
- `findsym_cifinput`：以 CIF 文件参数生成 FINDSYM 输入。它只转换，当前不会自动再
  运行 FINDSYM；其标准输出可作为下一次 FINDSYM 的输入。
- `smodes`：按 `vendor/isobyu/smodes.txt`，通过标准输入计算位移对称模式。
- `comsubs`：按 `vendor/isobyu/comsubs.txt`，通过标准输入搜索公共子群与结构路径。
  搜索范围受用户输入和官方程序默认值限制；超时不意味着没有结果或不存在子群。

文件上传保留原始字节；在网页中编辑内容后改为 UTF-8 文本。包装层不改变晶格、
坐标、单位、空间群设定、容差、归一化或 k 向量。请注意不同官方程序的约定可能不同，
以对应手册为准。

例如，运行官方发行包自带的 FINDSYM 输入示例：

```powershell
python -m cris_rea run findsym --input vendor/isobyu/findsym_sample.in --timeout 30
```

运行 SMODES 示例：

```powershell
python -m cris_rea run smodes --input vendor/isobyu/smodes_sample.in --timeout 60
```

输出 JSON 中包含原始输出预览、执行状态和本次目录。CLI 退出码 `0` 只对应
`completed_uninterpreted`（进程退出码为零）；官方程序可能在退出码为零时报告输入错误，
因此仍须读完整输出。进程非零、超时或运行环境错误返回非零 CLI 退出码。
显式 `Syntax error`、行首 `Error`、`This program has bombed` 和 Fortran runtime error
会标为 `program_diagnostic`，即使官方进程退出码为零也返回非零 CLI 退出码；这只是
有限的错误识别，不能证明其余输出均无问题。

准备一条受限 ISO 查询（先只显示命令，不执行）：

```powershell
python -m cris_rea query --parent 221 --irrep GM4- --opd P1 --wyckoff a
```

核对命令后，添加 `--run` 通过同一执行器运行：

```powershell
python -m cris_rea query --parent 221 --irrep GM4- --opd P1 --wyckoff a --run
```

这只是一个公开输入示例，参数均由用户指定，程序没有预存该组合的答案。
网页 **Prepare a restricted ISO query** 区域提供相同功能，点击
**Prepare ISO commands** 后可查看原始命令，再点击运行。当前表单不处理带自由参数
的 k 向量、Wyckoff 坐标、完整母相 CIF 或归一化模式导出；需要这些内容时不能把
默认/缺省输入的输出当作已满足课题目标。

## 运行产物与边界

每次计算在 `runs/<UUID>/` 内产生：

- `input.dat`：提交给程序的原始字节。
- `work/`：独立工作目录，接收官方程序生成的日志等文件。
- `isodata/`：本次允许使用的可执行文件、`const.dat` 和 `data_*.txt` 的可写副本。
- `stdout.txt`、`stderr.txt`：原始标准输出和错误。
- `execution.json`、`process.json`：真实执行命令、环境覆盖项、进程状态和耗时。
- `evidence.json`：输入、代码、官方资源和产物的 SHA-256；执行状态、运行环境、
  修改后的数据库指纹、科学验证状态。

计算入口只从 `vendor/isobyu/` 选择官方通用资源，不能选择任意程序路径；不会读取
`evidence/` 中的网站结果、旧 CRIS 输出或参考答案参与计算。运行结束后仍保留
数据库副本，便于核查程序是否改写了它们。每次副本会占用磁盘；手工确认不再需要的
运行后，可自行删除对应的完整 UUID 目录。

**目录隔离不是操作系统沙箱。** 官方命令语言可能包含文件操作，仅执行可信输入。
服务采取本机绑定、Host/Origin 检查和会话令牌，拒绝跨站执行；它不提供多用户权限、
文件系统隔离或不可信用户代码执行服务。单个输入上限 1 MiB，超时上限 300 秒，
标准输出与标准错误合计上限 20 MiB；上限只属于工程保护，不是科学筛选参数。
程序自行创建的其他文件不受输出字节上限约束。

## 工程结构与检查

调用方向为 `frontend → application → backend`：

- `cris_rea/frontend/`：CLI、HTTP、网页。
- `cris_rea/application/`：串行运行协调、运行记录读取。
- `cris_rea/backend/`：输入边界、资源副本、官方进程监管与证据指纹。
- `tests/`：工程合同测试。测试替身不用于科学结果验收。
- `docs/LOCAL_CONTRACT.md`：科学语义与执行合同、适用范围。

```powershell
python -m unittest discover -s tests -v
```

工程测试不能替代真实 Linux 程序运行、公式推导、独立数值不变量或网站对照。当前不
宣布全套 ISODISTORT 科学验收通过；具体实测证据以 `evidence/` 中对应报告为准。
