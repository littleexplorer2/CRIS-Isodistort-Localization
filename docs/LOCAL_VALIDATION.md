# 本地执行器工程验证记录

本记录只覆盖从空目录建立的 `cris_rea` 包及其执行、输入与证据合同。
官方程序的科学计算、网站交互与逆向结论见各自证据；本记录不替代它们。

## 本次修改与科学影响

新增三层标准库 Python 实现，执行五个白名单官方 Linux 程序；Windows 经 WSL 调用。
增加有限 ISO 命令生成、原始输入文件/文本入口、运行隔离目录、资源复制校验、原始
输出与证据记录，提供 CLI 和仅监听回环地址的英文网页。

包装层没有新增科学公式、拟合或数值近似。它的科学相关职责是确保输入选择、
字节传递、命令顺序和文件资源正确，不静默改变坐标、单位、晶格设定、IR/OPD、
容差或归一化。受限查询只序列化用户参数；标签是否适用于所选空间群，及未指定的
k 参数和 Wyckoff 坐标是否足够，必须单独核对。详细依据和适用范围见
[LOCAL_CONTRACT.md](LOCAL_CONTRACT.md)。

本层采用的工程界限为输入 1 MiB、stdout/stderr 合计 20 MiB、运行上限 300 秒；
这些不是科学容差或搜索完备性边界。无新数值容差。

## 最终可追溯证据

以下两个目录均包含原始 `stdout.log`、`stderr.log` 和 `record.json`；记录包括真实
命令、开始/结束 UTC 时间、退出码、平台/Python 信息、源码与测试 SHA-256。
最后保留日志的重跑是为了补齐证据文件，同时覆盖审查后新增的定向反例。

- Windows：`evidence/validation/20261020T135545223268Z-win32/`
  - 启动：`python tools/run_engineering_tests.py`
  - 记录内实际测试命令：该 Python 解释器执行
    `-m unittest discover -s tests -v`。
  - 28 项测试，23 项通过，5 项仅在 Linux 执行的符号链接/进程组测试跳过。
  - 退出码 0。Node v24.22.0 用于网页输入竞态的 DOM 适配器测试。
- WSL Linux：`evidence/validation/20261020T135602600670Z-linux/`
  - 启动：
    `wsl.exe --exec python3 /mnt/c/Users/devou/OneDrive/Desktop/CRIS_REA/tools/run_engineering_tests.py`
  - 记录内实际测试命令：`/usr/bin/python3 -m unittest discover -s tests -v`。
  - 28 项测试，27 项通过，1 项因 WSL 无 Node 的网页竞态测试跳过。
  - 退出码 0，执行了真实 Linux Python 子进程的非零退出、超时进程组终止与输出上限。

两个平台的覆盖互补；不把跳过项写成对应平台已经执行。网页竞态测试使用精简 DOM
适配器，验证输入选择与异步顺序，不代表完整真实浏览器验收。

## 真实本地浏览器到官方 ISO 的执行证据

2026-10-20 14:03–14:04 UTC，REA 启动受其管理的 Chromium 浏览器，实际访问
`http://127.0.0.1:8765/`，填写/点击本地网页，通过 HTTP 提交，再由 Windows
Python 调用 WSL Python 和每次运行目录内的官方 `iso` 副本。结束后 REA 终止
自己启动的浏览器。以下记录证明这两个有限场景的真实执行与数据传递；
它们不是对全部浏览器行为、全部 ISO 命令或晶体学正确性的验收。

- 正常查询：浏览器填写 parent `221`、irrep `GM4-`、OPD `P1`、Wyckoff `a`，
  点击生成输入，再点击运行。`POST /api/iso-query` 和 `POST /api/run` 均返回
  HTTP 200；这只是 HTTP 成功，不等于科学正确。
  [脱敏摘要](../evidence/validation/live-http/success/summary.json)、
  [可访问性快照](../evidence/validation/live-http/success/accessibility.txt)、
  [未编辑截图](../evidence/validation/live-http/success/final.png)和原始输入输出
  保存在同一目录。run ID 为 `8a69195a-4f22-4574-b980-451275c41fb5`，
  REA evidence ID 为
  `ev_5308c0d89e7d80cb537e045f28c06e98e49df50bc000e63b8e3064cffbda39d4`。
  官方程序退出码为 0，网页与保存结果均为 `completed_uninterpreted`，
  `scientific_validation` 为 `not_performed`，未命中文本诊断。
  原始 stdout 报告 ISO 9.6.1、`99 P4mm`、index 6、size 1、P1、
  basis `(0,1,0),(0,0,1),(1,0,0)`、origin `(0,0,0)`，并输出该 Wyckoff
  位置的投影向量 `(1,0,0)`。这里只转述原始输出，不把单例输出升级为理论证明。
- 故意错误命令：浏览器提交含 `THIS_IS_NOT_A_COMMAND` 的原始 ISO 输入，
  `POST /api/run` 返回 HTTP 200，官方程序仍退出 0，但 stdout 包含
  `***Syntax error: THIS_IS_NOT_A_COMMAND`。应用层和网页均给出
  `program_diagnostic`，诊断 code 为 `syntax_error`。
  [脱敏摘要](../evidence/validation/live-http/negative/summary.json)与
  [可访问性快照](../evidence/validation/live-http/negative/accessibility.txt)
  对应 run `6a34fad6-bf5d-4157-87ca-086bc24adc13`、REA evidence
  `ev_c2a71bb109f6829f44b92d418433920119e504aca8c4519deebfbcd4404a784c`。
  该反例直接覆盖“官方退出 0 不得自动视为成功”的工程合同。

每个场景分别复制了 `input.dat`、`stdout.txt`、`stderr.txt`、`evidence.json`、
`execution.json`、`process.json` 和 supervisor 的两个原始流；没有复制可执行文件、
数据库或工作目录。逐字节核对确认 HTTP 提交输入等于实际 `input.dat`，
输入 SHA-256 与运行 manifest 一致，HTTP 返回的 evidence 与保存 manifest 一致，
本例完整 HTTP stdout/stderr 与保存的原始 UTF-8 流一致。两次运行记录中的
生产源码指纹均与提取时磁盘源码一致；磁盘指纹本身不证明某个长期运行进程已加载
同一版本。输入 SHA-256 分别为
`05ff2d09697e74ff9f0c61de4620088ad3dd7e25818ed1a6faa4abaf508d69df` 和
`2524bdb4ada5e289c6ce048a206293591077ba3aa5c4cc1cdea67ba8449c70ec`。

两次官方 stderr 均为空；supervisor stderr 各保留 260 字节的 WSL localhost
代理/NAT 环境提示，按原始 UTF-16LE 字节保存，不当作官方 ISO 错误，也不丢弃。
成功场景的截图是 REA 捕获 PNG 的原始字节，视口位于输入框和运行按钮附近，
不能单靠该截图判断上方的结果内容；结果状态另由可访问性快照、HTTP 响应及本地
运行文件交叉佐证。负例未请求截图，未补造图片。

### 脱敏、来源与 REA 完整性边界

`tools/capture_web_local_safe.py` 从原始 REA bundle 导出上述派生材料。原 bundle
含本次回环 HTTP 会话令牌，保留本地且由 `.gitignore` 精确排除这两个路径；
完整 DOM、页面脚本和 header values 不进入可提交材料。提取器在内存识别令牌、
对可导出文本作替换，并检查导出字节不含已知令牌；截图只做 base64 解码与 hash
核对，没有编辑图片。JSON 摘要明确标记为派生证据，不冒充原始 REA record。

- `rea-local-browser-bundle.json` 原始 SHA-256：
  `a2c4347d4c318a9fa882950421273cf0e8f0923c8a00424337d6bea112eff4e3`。
- `rea-local-negative-bundle.json` 原始 SHA-256：
  `9f371be60567c062a37f602a4ef9dc950b128d95119dacbd9f4f033c7d8a61df`。
  后者也含前一个场景，提取时按 evidence ID 选取对应负例。

REA 对两条记录均给出 `status: incomplete`、`missing_sections: ["events"]`、
`equality_eligible: false`。具体 `selected network content unavailable` 记录是
`request-1` 的 `GET /` 和 `request-2` 的 `GET /api/programs` 请求体状态
`not_exposed`；不能把这个状态解释为已经证明请求无 body。两条记录均保留
`events_dropped: 0`，没有 budget/truncation 事件遗漏，所有被选中的 POST
请求体及其响应体均已捕获。保留 REA 的不完整结论，不改成 complete。

网络响应字节来自 Playwright `Response.body()`，其提供者明确不保证等于原始
网络线上的字节，底层可能为缺失 body 另行加载；事件顺序仅是接收顺序，且无请求
initiator stack，不能据此单独证明因果关系。本记录通过 HTTP payload、实际本地
运行文件和源码指纹补充关联证据，不宣称消除了这些提供者限制。

执行 `python tools/capture_web_local_safe.py` 完成提取，随后独立核对两个原始
bundle SHA-256、19 个导出文件的长度/指纹及复制文件与原始运行 manifest 的指纹。
提取脚本 AST 检查通过；两个原 bundle 的 `git check-ignore -v` 均命中上述精确
排除项。按 Git 可提交候选文件清单扫描当时的 304 个文件，未发现这两次捕获中
识别到的会话令牌；`git diff --check` 无报错。这些是证据与脱敏检查，没有重新
执行计算或扩大科学验收范围。

## 可检验不变量及反例

- 输入文件字节（含 CRLF）不变，证据 SHA-256 等于提交字节；每次运行使用新的 UUID。
- 修改本次 `isodata/` 副本不改变 `vendor/`；已有样例日志和网站证据不进入复制资源。
- 拒绝未知程序、路径形式的程序名、NUL、空/超长输入、非有限/越界超时和目录穿越。
- 对 Windows junction/符号链接使用解析后路径边界；Linux 实际符号链接反例确认
  在 `runs` 指向外部目录时写入前拒绝，并拒绝读取链接形式的 `evidence.json`。
- 复制后文件参数使用 `../input.dat`，`ISODATA` 使用 `../isodata/`，避免已发现的
  旧 Fortran 长路径缓冲区截断；测试核对实际组装的参数。
- 退出码非零、超时、输出超限不会成为完成状态；退出码零但含 `Syntax error`、
  连续 `***` 提示符后的错误、CIF 转换器显式错误/终止文本均标记诊断。
- HTTP 校验 Host、Origin 与会话令牌，拒绝超长请求和非法可选字段；未送完请求体
  有有限读取超时，不无限等待。
- 文件正在读取时不允许提交旧文本；先选文件的迟到结果不能覆盖后选文件、手工
  编辑或新生成的查询。测试对比真正交给 API 的内容。
- ISO 查询生成允许合法格式的用户参数，拒绝换行命令注入；预览操作不执行程序。

工程测试使用的假 ELF 字节、假数据文件、假 stdout 和 Python 小进程仅是控制流替身。
它们不表示官方科学样本，不参与任何科学正确性结论。

## 环境问题与剩余边界

最初 Windows 测试使用沙箱默认临时目录时遇到 ACL 拒绝；改为在本项目
`runs/engineering-test-temp/` 内创建临时目录后执行正常。记录脚本显式设置其子进程
的 TEMP/TMP/TMPDIR，不改系统设置。另一次早期测试因模拟平台后标准库探测被测试
替身误拦截而失败，已把替身范围固定；这些不是官方程序科学失败。

没有对任意 ISO 命令提供操作系统级文件沙箱；只接受可信输入。有限文本诊断未覆盖
所有可能错误，未命中不等于科学正确。官方程序自己生成的其他文件不受 stdout/stderr
字节上限约束。外层 supervisor 失联时明确无法确认子进程是否已退出。

本层没有验证群闭包、晶格指数、模式正交性、物理归一化或完整网站方法等价性，
也没有实现完整母相 CIF 到子群/畸变模式的工作流。这些仍属于单独科学研究与验证任务。

## 合成 CIF 转换器正例

另以预先冻结的独立合成 CIF 执行真实 `findsym_cifinput → findsym` 两步调用：
一般三斜晶格、显式 `P1`、恒等操作 `x,y,z`、Si/O 两种原子。
[输入来源与操作说明](../evidence/native/cif-converter-positive/README.md)、
[运行前输入清单](../evidence/native/cif-converter-positive/input_manifest.json)和
[14 项格式核对](../evidence/native/cif-converter-positive/verification.json)均已保存，
同目录保留两个程序的原始 stdout/stderr、运行 evidence 与执行记录。

转换器输出逐字节作为 FINDSYM 输入；晶格、类型、分数坐标、占位率和恒等操作保持一致，
小数采用精确 Decimal 比较，未放宽容差。两程序退出码均为 0、未命中已识别的文本诊断；
FINDSYM 接受转换后的格式并输出 P1 及含两原子的 CIF。两个输入清单均先于各自运行写入，
应用源码指纹未变化，数据库副本只保留在忽略的 `runs/`，没有复制进证据目录。

这是一个明确标为合成的格式/下游接受正例，未读取官网结果来反推输入。
输入本来已声明 P1 和操作，故此结果不证明独立发现空间群，也不覆盖任意 CIF 语法、
磁性、其他 setting、数值鲁棒性或全面科学正确性。
