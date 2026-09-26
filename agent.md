# GD — Agent Guide

## 0. 项目提示词（原文）

> 本文件夹（桌面 `GD（未同步git）`）存放 EuAl4 畸变模式振幅的梯度下降拟合代码。`GD_modified.ipynb` 是可修改的工作副本；`LD1_C1_gradient_descent_tianren.ipynb` 是原始参考笔记本，**不允许修改**。`.venv/` 是本目录专用虚拟环境（Python 3.12 + TensorFlow），不要手改其中第三方包文件。实验衍射表、模式名、归一化因子、振幅上界等数据在笔记本里以外部路径读入，那些原始数据文件不允许改。`save_best_model_parameters.py` 负责把最佳模式振幅写成桌面 `Best_Model_Parameters/` 下的 CSV，允许修改。`main_requirement.py` 负责创建/复用 `.venv` 并安装依赖，允许修改。说明文档 `README.md`、本文件 `agent.md` 允许修改。

> 以上为项目提示词原文。以下章节在其基础上展开可执行的目录权限、思考方式、验证科目、实现约定与常见坑；冲突时以本提示词的边界（哪些文件不允许修改）为准。

本文件供 Cursor / AI Agent 在本 GD 目录工作时阅读。优先遵守此处的修改边界、思考方式与验证清单；使用细节以 `README.md` 为准。

---

## 1. 项目介绍

本目录用 **TensorFlow / Keras 梯度下降**，把实验衍射强度拟合到一组 ISODISTORT 位移模式振幅上。物理上：母相 EuAl4（四方，空间群 I4/mmm）沿 LD1 / C1 子群发生畸变；每个模式振幅改变原子坐标，进而改变结构因子与强度。机器学习上：HKL 是特征，归一化实验强度是标签，可学习权重是 `mode_num` 个模式振幅，损失是晶体学 R 因子。

与 CRIS 的关系：拟合得到的 Best Model Parameters 由 `save_best_model_parameters.py` 写到桌面 CSV；CRIS 里的 `ISOVIZ_INPUT` 再把该 CSV 写入 `.isoviz` 并打开 IsoVIZ。本目录**不**实现 ISODISTORT 搜索或 CIF 导出。

### 1.1 目录与修改权限

| 路径 | 作用 | 可否修改 |
| --- | --- | --- |
| `GD_modified.ipynb` | 可迭代的拟合工作副本（含自动导出 CSV） | **允许** |
| `LD1_C1_gradient_descent_tianren.ipynb` | 原始参考笔记本 | **禁止** |
| `save_best_model_parameters.py` | 把最佳振幅写成桌面 CSV | **允许** |
| `isodistort_to_gd.py` | 把本地 ISODISTORT Method 2 结果写成笔记本输入 | **允许** |
| `generated/` | 转换器输出（不覆盖 `D:\` 原始精修文件） | **允许**（可整目录重生） |
| `tests_dev/` | 转换器格式短测 | **允许** |
| `main_requirement.py` | 创建 GD `.venv`、安装依赖 | **允许** |
| `requirements.txt` / `requirements-dev.txt` | 运行时 / 开发依赖 | **允许** |
| `README.md` / `agent.md` | 使用说明与 Agent 指南 | **允许** |
| `.venv/` | 本目录虚拟环境 | **禁止**手改 site-packages；只用安装脚本维护 |
| 笔记本读入的外部数据（`D:\OneDrive\...` 下的 CSV / txt 等） | 实验强度与模式定义 | **禁止** |
| `LD1_C1_alris_functions.py` | 由 `isodistort_to_gd.py` 生成的结构因子模块 | **允许重生**；不要手改或把某次拟合振幅写死进去 |

### 1.2 使用方式

1. （可选）用 CRIS `.venv` 跑 `isodistort_to_gd.py`，生成 `generated/LD1_C1/` 下的 txt 与同目录 `LD1_C1_alris_functions.py`。不要覆盖 `D:\OneDrive\...` 原始文件。
2. 在本目录运行 `python main_requirement.py`（需要本机 Python 3.12.x）。
3. 用 `GD/.venv/Scripts/python.exe` 作为 Jupyter 内核。
4. **改代码时只打开 `GD_modified.ipynb`。** 要对照原始流程，只读 `LD1_C1_gradient_descent_tianren.ipynb`。要把本地生成的模式表接上工作副本，把第一个单元格的 `path` 指到 `generated\LD1_C1`（`sub_name` 仍拼 `{path}\{irrep}_{structure_type}`）。实验 CSV 仍须自备，转换器不生成 `All Combined.csv`。
5. 按单元格顺序运行：配置 → 读数据 → 定义 `fun_tf` / `FunAsLayer` → 训练 → 打印最佳振幅。`GD_modified.ipynb` 会在分析振幅后自动写出桌面 CSV。

### 1.3 环境

```text
GD（未同步git）/.venv
```

本目录**不与** CRIS `.venv` 共用（GD 需要 Python 3.12 + TensorFlow）。一律用：

```powershell
cd <本 GD 目录>
python main_requirement.py
.\.venv\Scripts\python.exe
```

**解释器寻址约定**（与 `main_requirement.py` 一致）：`py` 启动器只按 major.minor 登记（`py --list` 只有 `-V:3.12`）。`py -3.12.5` 会报 `No suitable Python runtime found`。要指定 3.12 写 `py -3.12`；精确解释器写 `.\.venv\Scripts\python.exe`；脚本不要用补丁号拼启动命令，按 **major.minor** 校验。

---

## 2. 思考要求（强制）

**不要为了得到某一次好看的振幅而死记硬背数值，而是思考如何改进通用的物理模型、损失与优化过程，使同一类衍射拟合都能正确运行。**

以下两条为每轮 Bug 修复的强制要求：

I. 修 bug 时，不要只把打印结果改成「看起来对」的振幅表，而是思考此处应有的结构因子 / R 因子 / 模式振幅约束；不允许把某次拟合的 48 个振幅抄进 `if irrep == "LD1"` 式硬编码。

II. 每次修改完成后，用可重复的检查确认行为（导入、CSV 导出、必要的短测）；若结果仍有出入则继续修改，不要靠改外部实验 CSV 来迁就代码。

展开含义：

1. **禁止**把某次 `best_pars_overall` 的数值写死进笔记本或 `save_best_model_parameters.py`，只为「导出一份指定 CSV」。
2. **应当**改通用部分：结构因子前向（`fun_tf` / `LD1_C1_alris_functions`）、振幅约束、损失（R 因子）、HKL 对齐、归一化、CSV 命名与排序。同一套逻辑应对任意 `irrep` / `structure_type` 在改路径后仍能跑。
3. 区分 **实现 bug**（索引对错、归一化用错列、CSV 列名与 IsoVIZ 不兼容）与 **优化/数据限制**（局部极小、实验噪声、模式数不足）。限制写进 README，不要用假强度伪装已经拟合成功。
4. **不以「再跑满 1000 epoch 且 R 因子必须低于某阈值」为默认完成标准**（训练很慢、依赖本机 GPU/数据）。默认完成门槛见第 3 节：导入、导出契约、短测、文档。用户明确要求复现拟合时，再跑训练单元格。

---

## 3. 每次修改完成后的必做事项

### 3.0 文档、配置同步与清理

**强制规则：**

I. 每次行为变化后，检查并按实际行为修改所有受影响的说明文件（`README.md`、本文件若流程/边界变化、CSV 导出约定）。

II. 工作完成以后，整理或清理与项目代码无关的文件，包括本轮产生的临时截图、缓存、日志、测试写出的假 CSV。确需保留的验收产物必须放入约定目录（桌面 `Best_Model_Parameters/<irrep>_<structure_type>/`）并在交付说明中注明。

III. 项目版本与历史变更由 Git 管理（若本目录之后纳入版本库）。仅在用户明确要求时 commit / push；不改 git config；不用破坏性 rebase。

### 3.1 文档与配置同步

凡改动训练入口、损失、输入路径约定、CSV 命名/排序或已知限制，必须同步更新：

- 本目录 `README.md`
- 若改导出布局：`save_best_model_parameters.py` 的文档字符串
- 本文件 `agent.md` 仅在流程/边界变化时更新

### 3.2 完整验证测试科目

在声称「本轮修改完成」之前，按范围执行下列科目（能跑尽则跑尽；缺数据或未装 TensorFlow 时注明跳过项）。

#### A. 静态 / 短测（改代码后默认必跑）

在本 GD 目录，用 **GD `.venv`**：

```powershell
cd <本 GD 目录>
.\.venv\Scripts\python.exe -c "import tensorflow as tf; import numpy; import pandas; print(tf.__version__)"
.\.venv\Scripts\python.exe -c "from save_best_model_parameters import save_best_model_parameters, default_output_path; print(default_output_path(irrep='LD1', structure_type='C1'))"
```

若已安装开发依赖（`python main_requirement.py --dev`）：

```powershell
.\.venv\Scripts\python.exe -m pytest -q --tb=line
```

转换器格式短测（不需要 TensorFlow / WSL）：

```powershell
C:\Users\devou\OneDrive\Desktop\CRIS\.venv\Scripts\python.exe -m pytest tests_dev -q --tb=line
```

可选：

```powershell
.\.venv\Scripts\python.exe -m ruff check save_best_model_parameters.py main_requirement.py
```

改 `save_best_model_parameters.py` 时额外检查：

| 科目 | 检查意图 |
| --- | --- |
| 桌面根目录 | 自动创建 `Desktop/Best_Model_Parameters/`（OneDrive Desktop 优先） |
| 子文件夹 / 文件名 | `{irrep}_{structure_type}` / `{irrep}_{structure_type}_best_model_parameters.csv` |
| 列 | `Mode, Mode Name, Best Model Parameter, Maximum Mode Amplitude, Normalized Amplitude` |
| 行序 | `a1, a2, …, a10` 按整数序号，不是字符串字典序 |
| 不污染用户数据 | 测试写出的假 CSV 用完删除，不要留在用户桌面当真结果 |

#### B. 笔记本与拟合（仅当改 `GD_modified.ipynb` 的训练/前向时）

1. 内核必须是 `GD/.venv`（Python 3.12），不要用 CRIS `.venv`。
2. 自上而下能执行到「定义 `fun_tf` / `FunAsLayer`」；缺 `D:\` 数据或 `LD1_C1_alris_functions` 时在交付说明里写明，不要改禁止修改的原始数据来让单元格变绿。
3. **不要**把满 1000 epoch 当作默验收尾。用户要求复现拟合时，再跑训练循环，并核对：
   - 打印的 Best model parameters 条数 = `mode_num`
   - CSV 已出现在桌面对应子文件夹（仅 `GD_modified.ipynb`）
   - R 因子相对改前没有无解释地变差
4. **禁止**改 `LD1_C1_gradient_descent_tianren.ipynb` 来验证；对照时只读。

#### C. 与 IsoVIZ 衔接（仅当改 CSV 列名、模式名或导出路径时）

CRIS `ISOVIZ_INPUT` 依赖本导出的列名与 `Mode=aN` 顺序。改导出后，在 CRIS 根目录跑：

```powershell
.\.venv\Scripts\python.exe -m pytest ISOVIZ_INPUT\tests_dev -q --tb=line
```

注意：这是 **CRIS `.venv`**，不是 GD `.venv`。

#### D. 完成标准（DoD）

- [ ] 未修改禁止文件（尤其 `LD1_C1_gradient_descent_tianren.ipynb`、外部实验数据、`.venv` 第三方树）
- [ ] 相关短测 / pytest 已跑；失败已修复或说明为环境 skip（缺 TF、缺 D:\ 数据等）
- [ ] 涉及导出时：桌面路径、文件名、列、排序与 README 一致；测试假文件已清理
- [ ] 行为变化已写入 README / 本文件（若流程变化）
- [ ] 未引入「单次拟合数值硬编码」冒充通用正确性
- [ ] 优化限制与真实 bug 已在说明中区分
- [ ] 所有受影响的说明文件已更新
- [ ] 已整理本轮非项目代码文件；保留的验收 CSV 位于约定桌面目录且已说明

---

## 4. 实现约定（Agent 速查）

### 4.1 代码范围

- 可学习振幅、训练循环、可视化：只改 `GD_modified.ipynb`。
- CSV 落地、命名、排序：只改 `save_best_model_parameters.py`。
- 环境：只改 `main_requirement.py` 与 requirements 文件。
- ISODISTORT → 笔记本输入：只改 `isodistort_to_gd.py`，输出到 `generated/`；用 CRIS `.venv` 跑。
- 结构因子内核是生成的 `LD1_C1_alris_functions.py`：需要改算法时改转换器再重生，不要手改生成文件或把某次拟合振幅写死进去。
- 不要提交 `.venv/`、`__pycache__/`、`.pytest_cache/`、大型笔记本输出副本（除非用户明确要求）。

### 4.2 CSV 导出命名

- 根：用户桌面 `Best_Model_Parameters/`（自动创建）。
- 子文件夹：`{irrep}_{structure_type}`，非法字符压成 `_`。
- 文件：`{irrep}_{structure_type}_best_model_parameters.csv`（同名覆盖为最新一次拟合）。
- `Mode` 列为 `a1, a2, …`（与 IsoVIZ 顺序别名兼容）；写出前按模式整数序号排序。

### 4.3 对照数据

- `LD1_C1_gradient_descent_tianren.ipynb`：原始流程与超参对照，只读。
- 外部 `*_displacive_modes.txt` / `*_max_bound_vectors.txt` / `All Combined.csv`：实验与模式定义，只读。
- 不要把某次最佳振幅表抄进仓库当「标准答案」。

### 4.4 Git

- 仅在用户明确要求时 commit / push。
- 不强制 push；不改 git config；不用破坏性 rebase。

### 4.5 防止模块漂移、接口失控与上下文债务（强制）

1. **先找唯一责任模块，再改代码**：导出归 `save_best_model_parameters.py`，环境归 `main_requirement.py`，拟合归 `GD_modified.ipynb`，结构因子归 `alris` 模块。不要在笔记本里复制一份第二套 CSV 写出逻辑。
2. **接口契约先于实现**：变更 `save_best_model_parameters()` 的参数、CSV 列名或默认路径前，先列出生产者（笔记本）、消费者（ISOVIZ_INPUT）、默认值与错误语义；同步 README 与短测。禁止调用不存在的「猜测 API」。
3. **禁止靠改只读笔记本过关**：需要新行为就改 `GD_modified.ipynb` 或 `.py`，不要改 `LD1_C1_gradient_descent_tianren.ipynb`。
4. **状态必须显式传递**：`best_pars_overall`、`mode_names`、`max_mode_amps`、`irrep`、`structure_type` 作为函数参数传入导出函数，不要在导出模块里读 Jupyter `globals()`。
5. **小步闭环**：一次只解决一个可描述的问题；先能复现（短脚本或单元格），再改，再跑第 3 节科目。不要为通过单次拟合添加结构/IR 特例。
6. **上下文恢复清单**：开始工作先读本文件、`README.md`、目标笔记本/模块；记录本轮不变量（只读文件、CSV 契约、验收命令）。上下文不足时重新读取事实。
7. **完成前审计**：检查是否误改了 tianren 笔记本、是否把振幅写死、是否把测试 CSV 留在用户桌面；确认文档与代码一致。

---

## 5. 常见坑

| 现象 | 处理 |
| --- | --- |
| Jupyter 内核不是 3.12 / 没有 TensorFlow | 选 `GD/.venv`；先跑 `python main_requirement.py` |
| 终端无 `(.venv)` 却像卡住 | 用 `.\.venv\Scripts\python.exe` 显式启动 |
| `py -3.12.5` 报 No suitable Python runtime | 改用 `py -3.12` 或 `.\.venv\Scripts\python.exe` |
| `import LD1_C1_alris_functions` 失败 | 在本目录跑一次 `isodistort_to_gd.py`（CRIS `.venv`）；不要从 tianren 笔记本里删掉这行来「修好」导入 |
| 读不到 `D:\OneDrive\...` | 本机路径写在第一个代码单元格；改工作副本里的 `path` / `data_path`，不要改外部原始 CSV |
| `best_pars_overall is not defined` | 必须先跑完训练循环；不要为了导出去编造振幅 |
| CSV 有了但 IsoVIZ 匹配 0 条 | 列名须含 Best Model Parameter；模式名与 `.isoviz` 标签规范化后应对得上，或依赖 `Mode=a1,a2,…` |
| 训练极慢 / tf.function retracing 警告 | 前向应尽量固定 HKL 网格并缓存；不要在循环里反复建新 `tf.function` |
| 把 CRIS `.venv` 拿来跑本笔记本 | GD 需要 3.12+TF；CRIS 环境通常不够 |

---

## 6. 相关入口一览

```text
GD_modified.ipynb                      可修改的拟合笔记本（自动导出 CSV）
LD1_C1_gradient_descent_tianren.ipynb  只读原始笔记本（用法见 README.md）
isodistort_to_gd.py                    ISODISTORT → 笔记本输入（用 CRIS .venv）
generated/LD1_C1/                      转换器输出（txt + alris 源）
LD1_C1_alris_functions.py              笔记本同目录导入的结构因子模块
save_best_model_parameters.py          桌面 CSV 导出
main_requirement.py                    GD/.venv 与依赖
README.md                              tianren 笔记本原理与用法
```
