# CRIS — Agent Guide

本文件供 Cursor / AI Agent 在本仓库内工作时阅读。优先遵守此处的修改边界、思考方式与验证清单；细节以各子项目 README 与 `ISODISTORT/config/settings.yaml` 为准。

---

## 1. 项目介绍

**CRIS** 是一个针对开源晶体学工具 **ISODISTORT** 套件的本地化项目，目前包含：

- **ISODISTORT 本地化**（子群搜索、畸变模式、网页/终端/API、结果导出）
- **IsoVIZ 数据导入**（振幅 CSV → `.isoviz` → 启动 IsoVIZ）

仓库目标：在本地复现官网常用交互与计算；导出文件**内容与格式尽可能接近官网**，并以工具可打开与晶体学正确性为硬验收。网页与终端共用 `isocore` 计算/导出路径，**仅交互方式不同**。

**导出验收（与 `ISODISTORT/README.md`、`config/settings.yaml` 一致）**：

- `Complete modes details.txt`：完整写入本地计算结果即可；**不要求**与官网 HTML 逐字节一致。
- `subgroup.cif` / `data.isoviz` / `topas.str`：字段与布局尽量靠官网；不以全文件逐字节为 DoD。
- **可用性**：`.cif` 须能用 **[VESTA](https://jp-minerals.org/vesta/en/)**（[下载](https://jp-minerals.org/vesta/en/download.html)；可选根目录 `VESTA.lnk`）打开；`.isoviz` 须能用 **ISOViz**（ISOTROPY IsoVIZ）打开。VESTA 用于三维查看子群 CIF，本仓库不附带、计算流程不调用。
- VALIDATE 默认做语义比较；`--strict` 仅用于排版调试，不是默认完成门槛。

### 1.1 目录与修改权限

| 路径 | 作用 | 可否修改 |
| --- | --- | --- |
| `experiment_data/` | 实验母相 CIF 等原始数据（如 `EuAl4 Parent.cif`） | **禁止** |
| `GD/` | 梯度下降拟合相关代码与笔记本 | **禁止** |
| `webpage_info/` | 官网各页 HTML 存档（按序号标注交互与跳转顺序） | **禁止** |
| `ISODISTORT/isobyu/` | 从 [iso.byu.edu](https://iso.byu.edu/isotropy.php) 下载的 Linux 二进制与 `data_*.txt` 数据库 | **禁止** |
| `ISODISTORT/`（除 `isobyu/`） | ISODISTORT 本地化主体（`isocore/`、`web/`、`config/`、`tests_dev/` 等） | **允许** |
| `ISODISTORT_VALIDATE/` | 比较本地 CIF 与官网 CIF，检查本地化是否有 bug | **允许** |
| `ISOVIZ_INPUT/` | 振幅 CSV + 子群 `.isoviz` → 写入并打开 IsoVIZ | **允许** |
| `output_compare/`（若存在） | 官网 / 本地导出对照黄金集（通常 gitignore） | 只读对照，勿当源码改 |

### 1.2 官网交互对照（`webpage_info/`）

存档页以序号表示先后顺序，典型验收路径：

1. **1 首页**：上传 `experiment_data/EuAl4 Parent.cif`
2. **2 Search**：畸变类型勾选 **strains**、**displacive（Eu, Al）**；如果用Method 2 的话，就选 **LD, k10 (0,0,g)，g=1/6**
3. 多次 **OK** → 依次进入 **3、4、5**
4. **6 Distortion**：选择导出格式并下载结果
5. **a order parameter direction files**：只勾选畸变类型为 **strains**、**displacive（Eu, Al）**，并使用Method 1 得到的结果导出页
6. **2.5 Search**, **3.5**： 畸变类型勾选 **strains**、**displacive（Eu, Al）**；Method 2 选 **LD, k10 (0,0,g)，g=1/6**，并调整 **Change number of superposed IRs:** 的参数得到的页面

本地实现应对齐该流程的语义（Types、Method、导出选项），UI 为**英语**。

### 1.3 使用方式

ISODISTORT 支持三种入口（共用 `isocore`）：

1. **网页**：`ISODISTORT/main_web.py`（`web/server.py` + `index.html`）
2. **终端**：`ISODISTORT/main_terminal.py`
3. **API**：`isocore.api`（如 `IsoDistort`）

ISOVIZ_INPUT：把振幅写入子群 `.isoviz` 并启动 IsoVIZ；开发用 `input_content/data.csv/` 与 `input_content/subgroup.isoviz/` **不上传远程仓库**。

### 1.4 共享环境

```text
CRIS/.venv
```

用根目录任一门户安装/复用：

```powershell
python ISODISTORT\main_requirement.py
# 或 ISODISTORT_VALIDATE\main_requirement.py / ISOVIZ_INPUT\main_requirement.py
```

一律用 `.\.venv\Scripts\python.exe` 跑测试与入口。Windows 上 `iso` 经 **WSL** 调用。

---

## 2. 思考要求（强制）

**不要为了生成正确结果而死记硬背正确答案，而是思考如何优化底层的计算代码，从而让所有这种类型的结果都能正确运行。**

展开含义：

1. **禁止**把某个 `output_compare` / 官网 CIF 的位点、HM 字符串、模式个数「抄进」特例 `if irrep == "X4-"` 式硬编码，只为通过单测或单例对比。
2. **应当**在 `isocore` 中修通用算法：原点选择、Hall/OC2、Wyckoff 判定、OPD/`k-active`、超胞变换、BUSH/smodes 作用域、导出命名等，使 **同一类** 子群/k 点/空间群都受益。
3. 对照官网时：先区分 **实现 bug**（可修）与 **引擎/算法限制**（本地 `iso` 无 strain / secondary OP、参数 k 在 nmod=0 无位移模式等）；限制写进 README「已知限制」，不要用假数据伪装已实现。
4. **不以全文件逐字节相同为默认完成标准**。优先修通用渲染与算法，使 CIF/isoviz/TOPAS 字段与布局尽量靠官网，modes `.txt` 写全本地结果；并用 VESTA / IsoVIZ 抽检可打开。需要排版对照时再用 VALIDATE `--strict` 或二进制 diff 定位差异。

---

## 3. 每次修改完成后的必做事项

### 3.1 文档与配置同步

凡改动行为、导出布局、Method 语义、已知限制或端口/超时等，必须同步更新：

- `ISODISTORT/README.md`（及根目录 `README.md` 若跨项目说明变化）
- `ISODISTORT/config/settings.yaml`（及相关注释）
- 若改 VALIDATE / ISOVIZ_INPUT 行为：对应子项目 README
- 本文件 `agent.md` 仅在流程/边界变化时更新

### 3.2 完整验证测试科目

在声称「本轮修改完成」之前，按范围执行下列科目（能跑尽则跑尽；WSL 不可用时注明跳过项）。

#### A. 静态 / 单元与集成（默认必跑）

在仓库根目录：

```powershell
cd <CRIS 根目录>
.\.venv\Scripts\python.exe -m pytest ISODISTORT\tests_dev -q --tb=line
```

覆盖科目（`tests_dev/`）：

| 科目 | 文件 | 检查意图 |
| --- | --- | --- |
| 基础与 Wyckoff 头 | `test_basics.py` | 母相识别、显示格式 |
| 结构 I/O | `test_structure.py` | CIF/结构读写 |
| 搜索 Method 1–4 | `test_search_methods.py` | 查询对象与过滤逻辑 |
| 畸变映射 | `test_distortion.py` | 模式→原子位移等 |
| 导出格式 | `test_formats.py` | OPD 行、ZIP 布局、CIF/isoviz/TOPAS 标记 |
| (3+d) 超空间 | `test_3pd.py` | 不依赖 WSL 的 superspace 内核 |
| 网页 API | `test_web.py` | 上传、Method、导出相关 HTTP |
| WSL / iso | `test_wsl.py` | 二进制可达性（无 WSL 则 skip） |
| 鲁棒性 | `test_robustness.py` | 异常输入与边界 |

最小烟雾（改动很小时至少跑）：

```powershell
.\.venv\Scripts\python.exe -m pytest ISODISTORT\tests_dev\test_formats.py ISODISTORT\tests_dev\test_basics.py -q --tb=line
```

可选：

```powershell
cd ISODISTORT
..\.venv\Scripts\python.exe -m ruff check .
```

#### B. 与官网结果对照（改 CIF / Method / 导出时）

对照目标是**尽量接近官网** + **工具可打开**，不是默认要求每一个字节相同。

1. 用本地网页或终端对 `EuAl4 Parent.cif` 走与 `webpage_info` 一致的 Types，导出 Method 1 和/或 Method 2 ZIP。
2. 与 `output_compare/官网/`（若有）或自备官网下载做对照，覆盖：
   - 文件夹命名（Method 1 完整 OPD 行且删除 `/`；Method 2 `IR OPD`）
   - 文件名：`subgroup.cif` / `data.isoviz` / `topas.str` / `Complete modes details.txt`
   - 关键字段：HM/Hall、Wyckoff、模式标签与数量、isoviz `!tag` 布局、TOPAS 模式块
   - 已知引擎限制（strain / secondary / nmod=0 等）允许模式更少；不要用假数据填满
3. 发现差异时：区分 bug vs 引擎限制；能通用修则修 `isocore` 渲染/算法，不要改对照黄金文件迁就本地输出。
4. CIF 语义对照（默认）：

```powershell
# 将本地 CIF → ISODISTORT_VALIDATE/compare/item/
# 官网 CIF → compare/true/（相对路径与文件名必须一致）
.\.venv\Scripts\python.exe ISODISTORT_VALIDATE\main.py compare
.\.venv\Scripts\python.exe ISODISTORT_VALIDATE\main.py batch
```

排版调试时再加 `--strict`。非 CIF 可用文本/二进制 diff 辅助，但 `byte_exact=否`  alone 不否决本轮完成。

5. **可用性抽检**（改 CIF / isoviz 时必做）：多选几个子群（勿只测 GM1+），用 **[VESTA](https://jp-minerals.org/vesta/en/download.html)** 打开若干 `subgroup.cif`，用 **ISOViz** 打开若干 `data.isoviz`（根目录可放 `VESTA.lnk` / `ISOViz.lnk`）。

#### C. 手工 / 长时科目（大改搜索、网页或批量导出时）

见 `ISODISTORT/tests_dev/manual/README.md`：

```powershell
cd ISODISTORT
..\.venv\Scripts\python.exe tests_dev\manual\run_web.py spotcheck
..\.venv\Scripts\python.exe tests_dev\manual\run_web.py method2_ld
..\.venv\Scripts\python.exe tests_dev\manual\run_batch.py cif30
```

#### D. ISOVIZ_INPUT（仅当改动该子项目时）

- 确认 `input_content/` 下 CSV 与 `.isoviz` 成对可读
- 跑通 `ISOVIZ_INPUT/main.py`（或该项目 README 中的测试说明），振幅写入后能启动 IsoVIZ

#### E. 完成标准（DoD）

- [ ] 未修改禁止目录
- [ ] 相关 pytest 已跑；失败已修复或说明为环境 skip
- [ ] 涉及导出时：内容/格式尽量靠官网；modes 写全；CIF/isoviz 已用 VESTA/IsoVIZ 多子群抽检可打开
- [ ] 行为变化已写入 README / `settings.yaml` / 本文件（若流程变化）
- [ ] 未引入「单例硬编码」冒充通用正确性
- [ ] 已知引擎限制与真实 bug 已在说明中区分
- [ ] 网页与终端仅交互不同，计算/导出走同一 `IsoDistort` API

---

## 4. 实现约定（Agent 速查）

### 4.1 代码范围

- 计算与导出逻辑放在 `ISODISTORT/isocore/`（`api` / `backend` / `distortion` / `io` / `superspace` / `utils`）。
- 网页与终端只做交互壳，避免复制一套计算逻辑。
- 不要提交 `output/`、`.pytest_cache`、`__pycache__`、大型对照 HTML 资源树（除非用户明确要求）。

### 4.2 导出命名（与官网 Windows 下载对齐）

- Method 1 文件夹：完整 OPD 行；非法字符处理时 **删除** `/`（`I4/mmm`→`I4mmm`，`1/2`→`12`），不要换成空格。
- Method 2/3 文件夹：`IR OPD`（如 `LD1 C1`）。
- ZIP 根下直接是各子群文件夹（无多余 `isodistort_methodN/` 包裹层，除非产品文档另有约定）。
- 单文件：`subgroup.cif`、`data.isoviz`、`topas.str`、`Complete modes details.txt`（官网 modes 为 HTML；本地用 txt）。

### 4.3 对照数据

- `webpage_info/`：交互与文案对照。
- `output_compare/`（本地）：官网 vs 现有导出；用于找格式/计算差，**不要**把黄金 CIF 整目录提交进 git（若已被 ignore）。
- `ISODISTORT_VALIDATE/compare/`：成对 CIF 比较；默认语义；`--strict` 仅排版调试。根目录 `VESTA` / `ISOViz` 快捷方式用于可用性抽检。

### 4.4 Git

- 仅在用户明确要求时 commit / push。
- 不强制 push；不改 git config；不用破坏性 rebase。

### 4.5 界面语言

网页与终端面向用户的文案保持 **英语**（`isocore/i18n`）；本仓库说明文档可为中文。

---

## 5. 常见坑

| 现象 | 处理 |
| --- | --- |
| 终端无 `(.venv)` 提示却「卡住」 | 用 `.\.venv\Scripts\python.exe` 显式启动；结束残留进程 |
| Method 1 很慢 | 首次枚举全部特殊 k 点，属预期 |
| Method 2 参数 k 无子群 | 需生成 isotropy 子群库（可能极慢） |
| ZIP 无子群 | 先对 Method 1/2/3 点 OK 再导出 |
| 模式数少于官网 | 查 README 已知限制（strain / secondary / no root mode / LD nmod=0）后再判 bug |
| VALIDATE 语义 PASS 但 `byte_exact=否` | 默认**不算**失败；优先语义与 VESTA/IsoVIZ 可打开。仅调试排版时用 `--strict` |
| IsoVIZ 打开崩溃 | 检查 `!displacivemodelist` 向量数是否等于该 `parentatom` 在 `!atomcoordlist` 中的行数 |
| 网页 vs 终端导出不一致 | 二者应都走 `export_subgroups_zip` / `_collect_export_specs`；差异只应在 UX |

---

## 6. 相关入口一览

```text
ISODISTORT/main_web.py          网页
ISODISTORT/main_terminal.py     终端
ISODISTORT/main_requirement.py  依赖与 .venv
ISODISTORT/config/settings.yaml 运行时配置
ISODISTORT_VALIDATE/main.py     CIF 比较
ISOVIZ_INPUT/main.py            振幅 → IsoVIZ
```

官网帮助：[ISODISTORT help](https://iso.byu.edu/isodistorthelp.php) · [ISOTROPY Suite](https://iso.byu.edu/isotropy.php)
