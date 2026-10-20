# ISODISTORT 官网协议与能力证据

记录日期：2026-10-20（香港）；网页自报版本：6.12.2，Jul 2025。
所有本报告产物均在 CRIS_REA。没有向官网上传用户实验文件，没有写入旧 CRIS。

## 结论与证据层次

已经观察到官网采用 HTML 表单驱动的多步计算流程：母相导入后，按
`input` 字段选择工作阶段，使用 `POST isodistortform.php` 返回下一张
表单或结果。Method 2 的母相 → k → IR → OPD → 模式结果链已走通一个
公开教学样例。它证明这个样例的协议与官网行为，不证明本地实现已经
等价，也不证明全部方法、全部晶体或全部参数已覆盖。

浏览器和 HTTP 证据不能看到 PHP 源码、服务器进程、实际调用的可执行
文件或数据库访问。本报告不把 URL 后缀、字段名或临时文件名视为服务器
内部实现的证明。与本地 Linux 二进制的对应关系，必须由本地静态和运行
证据单独建立。

证据类别：

- **运行观察**：REA 实际启动 Edge 无头浏览器、执行声明动作、记录请求，
  以及随后使用相同真实表单字段进行的 HTTPS POST。
- **官网陈述**：帮助文档描述的功能、算法策略及服务器限制；没有在本次
  实验中全部验证。
- **静态观察**：公开 JavaScript 和 HTML 的代码与字段。
- **分析建议**：本地化架构建议，不是对官网内部架构的断言。

## 证据入口

`evidence/web/public_pages_manifest.json` 记录每份公开资源的 URL、UTC
时间、响应状态、字节数和 SHA-256。HTML 原文、派生文本及表单结构与
manifest 同目录。正式引用应以原文的哈希为准。

`evidence/web/rea_web_bundle.json` 是完整 REA 导出，包含四份浏览器观察、
一次 JavaScript 分析和一次分析视图。重点 Evidence：

- 首页：`ev_5d4b5d781cc47fbdf4a0fe315f19a23543f7e51470b5d58b3f69c00a6940030b`。
- 自动加载母相前的表单：
  `ev_9250e87bd2044fe7a551ac2f685eb5b88f06b77676095b6150eaa7cb14f4b876`。
- 母相搜索页：
  `ev_6f70e9243bf29d87473062a0d6f4b56e0d19d1cbee6237698337876e728f4695`。
- 浏览器实际点击 Method 2、取得 GM 的 IR 选择响应：
  `ev_67da9fc70410b1fac55997e21067150dcdd527239ad9f170f057e632e2bb7abe`。
- 公开 `isodistort.js` 的静态分析：
  `ev_c7cf66e69ed4b3790677a5db91a49aa2e9e77cd9bfa33134eac50da36c3d64ae`；
  摘要视图 `ev_d944470e948ce6ed65b29836e7be8d7ac18b35af14d5b4daf721718a2df58f84`。

`evidence/web/rea_extraction_manifest.json` 将捕获响应关联到 Evidence ID、
事务编号、文件与哈希；提取过程对 21 个正文逐个验证 SHA-256。
同名 `ev_*` 目录保存该记录及可读 DOM、网络请求和响应。

REA 浏览器局限原样保留在 bundle。部分 SwingJS 响应超过选定正文预算，
所以浏览器记录不是完整网站资源归档。Playwright 的 `Response.body()`
结果不是原始线上字节的保证；该 API 也可能执行后续读取。它没有请求
发起者堆栈，事件先后不能单独证明因果。所有自有浏览器均报告
`terminated-owned-process`。

## 已观察的交互合同

首页 [isodistort.php](https://landau3.byu.edu/isodistort.php)：

- 母相 CIF：`POST isodistortuploadfile.php`，编码
  `multipart/form-data`，字段 `input=uploadparentcif`、
  `toProcess=<file>`。
- 保存畸变导入：同一上传端点，`input=distort`、
  `origintype=dfile`、`toProcess=<file>`。
- 公开钙钛矿入口 `isodistortparent.html` 在页面加载后自动提交
  `input=autoloadparent` 到 `isodistortform.php`。

搜索页包含畸变类型选择及四种 Method 的独立表单。畸变类型修改使用
`input=changesearch`；站点选择和应变选择经 Change 提交才作用于后续
页面。搜索表单携带大量隐藏状态：

- 母相：`spacegroup`、`lattparam`、`dlattparam`、`wycount`；
  每站点 `wyatomNNN`、`wyatomtypeNNN`、`wyckoffNNN`、
  `wyparamNNN`、`wyoccNNN` 等。
- 选择：`includedisplaciveNNN`、`includestrain` 等。
- setting：`settingaxesm`、`settingcell`、`settingorigin`、
  `settingaxesh`、`settingaxeso`、`settingssg`。
- 内部数值引用：`wypointerNNN`、`parentsetting` 等。观察到它们被
  往返传递，不能据此保证其跨版本稳定性或视为可公开复用的标准编号。

Method 2 的实际阶段：

1. `input=kvector`：`irrepcount`、`kvec1`、
   `kparam11/kparam21/kparam31`、`nmodstar1`，同时携带母相状态。
   本次默认 GM 对应选项 `1 *GM, k12 (0,0,0)`。参数索引属于网页
   协议；不是通用晶体学 ID。
2. 返回 IR 页面，下一次提交 `input=irrep`，选择 `irrep1`，
   可选填写 `irrepopd1`。返回的规范化状态包括 `kvecnumber1`、
   `kparam1` 等。
3. 返回 OPD 单选列表，每项显示 OPD、子群、基、原点、`s`、`i`；
   提交 `input=distort`、`origintype=method2`、
   `orderparam=<所选返回值>`。
4. 返回 distortion 页，包含模式及应变振幅输入、导出选项。
   下一阶段 `input=displaydistort`，用 `origintype` 选择导出。
   本次观察到 `modesfilename=isodistort_64169.iso`，
   `isosubgroup=12414`、`subgroupsym=99`、`subgroupsetting=172`。
   这些都是本次官网返回的状态，不得预置进本地求解器。
5. `origintype=completemodesdetails` 返回完整模式详情。
   页面还有 `isovizdistortion`、`isovizdiffraction`、
   `structurefile`、`distortionfile`、`domains`、`primary`、
   `modesdetails`、`topas`、`fullprof`、`irreps`、`tree` 等选项。
   本次没有逐个请求或验收这些导出格式。

其他方法的搜索入口也直接出现在母相搜索页：Method 1 使用
`input=isosubgroup`；Method 3 使用 `input=parentbasis`；
Method 4 以 `input=uploadsubgroupcif`、`toProcess` 通过 multipart
提交至 `isodistortuploadfile.php`。完整控件原文可在
`evidence/web/method2/01_search.html` 和 `.forms.json` 中审查。
本次没有完成 Method 1、3、4 的运行结果流程。

## 本次真实官网样例

输入来源是官网公开教学母相，不是实验数据：空间群 221，立方晶胞
`a=4.2 Å`，Sr 在 1b、Ti 在 1a、O 在 3d。选择位移及应变，
一个 IR，GM，`nmodstar1=0`，GM4-，P1 `(a,0,0)`，振幅保留零。

完整逐步证据位于 `evidence/web/method2/`：

- `01_search`：公开母相提交 → 搜索页。
- `02_irrep`：GM 提交 → IR 菜单。
- `03_opd`：GM4- 提交 → 六个 OPD 候选。
- `04_distortion`：P1 提交 → 子群及振幅表单。
- `05_complete_modes`：完整模式详情。

每一步保存 `.request.txt`、`.fields.json`、`.html`、
`.forms.json`、`.metadata.json`。元数据固定请求和响应指纹及
UTC 时间。`01`—`03` 来自脚本主流程，`04`—`05` 在检查上一步后
调用相同 `post()` 完成；脚本现已包含全部五步。

结果页实际显示 P4mm（99），基
`{(0,1,0),(0,0,1),(1,0,0)}`，零原点，`s=1`、`i=6`；
4 个位移模式和 2 个应变模式。零振幅使本次畸变后坐标等于未畸变坐标。
这不是非零振幅、非公度、磁结构、模式分解或泛化正确性的验收。
没有把这些官网结果输入本地计算，也没有据此宣称本地匹配通过。

## 客户端 JavaScript

公开文件 [isodistort.js](https://landau3.byu.edu/isodistort.js) 的 SHA-256：
`7cbff80027215bb4c0df6645df64c97a1ff1c298cd57280946319b7d8ae89f5d`。

直接审读可见：

- `fetchIsoDistort()` 收集表单输入、已选单选/复选框，按当前 form
  的 action/method 调用 `$.ajax`。
- `clickIsoDistort()` 请求 `isovizdistortion` 并将返回内容交给
  `createIsoDistortApp()`。
- 后者配置 SwingJS，入口名为 `org.byu.isodistort.IsoDistortApp`
  或 `IsoDiffractApp`，资源目录为 `swingjs/j2s`，并把响应数据、
  formData 与 document 作为参数。

这些是代码结构证据。本次没有点击 interactive viewer，也没有验证
SwingJS 程序内部的计算。代码中的外部 `serverURL` 是配置字符串，
本次不把它当作已发生的网络调用。

REA 只分析指定 `isodistort.js`（1967 字节、282 个 AST 节点），
未分析依赖；静态图报告 partial/unknown，并保留动态调用的未知项。
其 “0 findings” 不能解释为无功能或无请求。

## 官方能力陈述与本次审读范围

权威入口：[官方 Help](https://landau3.byu.edu/isodistorthelp.php)，
原文 SHA-256 为
`93bb70d6eafd98a5d657229456889077a645a00cae9ab35b4f4816de35b79d2c`。
本次实际审读介绍、母相/setting、IR 矩阵、畸变类型、四种方法、实时
子群、原子匹配和模式详情章节；论文仅识别引用，没有宣称读过其证明。

- [Method 1](https://landau3.byu.edu/isodistorthelp.php#isosub)：
  单 IR、特殊 k 点，使用预计算子群并按晶系、空间群、晶格等过滤。
- [Method 2](https://landau3.byu.edu/isodistorthelp.php#kpoint)：
  指定 k、IR、OPD，支持耦合及非公度；文档目前上限为 3 个独立调制。
- [Method 3](https://landau3.byu.edu/isodistorthelp.php#method3)：
  指定点群或空间群与代表晶格，搜索相容的取向、原点及畸变对称性。
- [Method 4](https://landau3.byu.edu/isodistorthelp.php#method4)：
  母相与畸变 CIF 的模式分解；要求相容的晶胞与一一原子映射。
- [实时子群](https://landau3.byu.edu/isodistorthelp.php#genisosub)：
  特殊点单 IR 等有表；其余情况可能按需生成并临时保存。
  官网声明此类任务一小时终止、临时文件每周清理。
- [原子匹配](https://landau3.byu.edu/isodistorthelp.php#trynearest)：
  有最近站点和阈值内组合搜索；分解任务声明 30 分钟终止。

总体功能包括位移、占位、磁矩、转动、应变，结构/衍射可视化、模式、
畴、IR 矩阵、子群树及精修格式导出。本次仅验证上文的有限流程。

## 科学合同核查及缺口

本工作只增加协议采集脚本和证据文档，没有实现或修改科学算法。
采集不得把官网下发的 IR/OPD 索引、子群、基、模式或归一化数值
当作独立本地求解输入；完整页面仅作为隔离参考。

Help 的 [Modes Details](https://landau3.byu.edu/isodistorthelp.php#modesdetails)
给出的合同是：位移模式使用无量纲超胞坐标，归一化用未应变晶格；
实际位移用应变后晶格。微观振幅区分 `As` 与
`Ap = As sqrt(Vp/Vs)`，这里体积为原胞体积；应变另行处理。
这属于官方定义，尚未逐一追到本地二进制的实现。

帮助中的数值算例存在待核查不一致：同一 `(1,-1,0)` 向量在前式使用
7.72 Å，随后位移代入显示 5.46 Å 和 0.0533 Å。不能把该印刷数字
直接作为正确答案或修改本地容差的依据。这里仅标记文本矛盾，不据此
断言官网运行算法错误。

本次精确字节校验没有数值容差；没有执行科学误差比较。未证明的部分：
PHP → 原生程序的调用链、命令参数、进程隔离、数据库版本配对、
所有高级方法、非零振幅导出、浏览器可视化及全量科学正确性。

已执行的定向检查：三个采集脚本通过 Python AST 语法解析；7 个公开
GET 文件与 manifest 的 SHA-256 相符；5 个 POST 返回与各自元数据
的 SHA-256 相符；21 个 REA 提取正文与工具记录中的 SHA-256 相符。
这些是采集完整性检查，不是科学验收。

## 复现与本地化建议

采集脚本在 `tools/capture_web_*.py`，只依赖 Python 标准库，
不应由生产计算导入。执行解释器可以是本项目解释器；本次使用现有
`E:\CRIS\.venv\Scripts\python.exe`，工作目录始终为 CRIS_REA。
使用这个解释器不代表写入旧项目。

1. `capture_web_public.py` 获取明确列出的公开 GET 页面。
2. REA 的完整浏览器参数和动作保存在 Evidence 的 `parameters` 中。
3. `capture_web_extract_rea.py` 从指定 REA bundle 提取并校验正文。
4. `capture_web_method2.py` 使用已存公开示例母相走五步 POST；
   遇到已存响应即停止，避免覆盖证据。重新采集应选择新的输出目录。

本地化应把显式输入、规范化会话、科学引擎、结果导出分开。
网页隐藏字段可以启发用例状态，但生产端应自行验证输入并生成本地
结果，不能信任客户端提供的数据库指针。官网的临时结果文件名可以
启发本地 job/result 标识，但不得把其部署细节写成已经逆向确认的事实。
