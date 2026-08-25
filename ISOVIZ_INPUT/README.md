# ISOVIZ_INPUT

把梯度下降得到的振幅 CSV 写入官方 IsoVIZ 的 `.isoviz` 子群结构文件，并一键打开 Java 版 IsoVIZ 查看畸变后的晶体结构。

IsoVIZ 属于 [ISOTROPY Suite](https://iso.byu.edu/isotropy.php)，用进度条手动调每个原子模式既慢又容易出错。本子项目用 Python 把 CSV 中的 **Best Model Parameter** 写进 `.isoviz` 里对应模式的 `amp` 字段，再启动 IsoVIZ。

与 `ISODISTORT`、`ISODISTORT_VALIDATE` 共用仓库根目录的 `CRIS/.venv`。

---

## 仓库边界（给后续 AI / 协作者）

CRIS 文件夹是一个针对开源晶体学工具 ISOTROPY 套件的本地化项目，目前包含 ISODISTORT 本地化与 IsoVIZ 数据导入。

- `experiment_data/` 存放实验母相 CIF 等原始数据，**不允许修改**。
- `GD/` 存放梯度下降拟合相关代码与笔记本，**不允许修改**。
- `webpage_info/` 记录了 ISODISTORT 官网各个页面的 HTML 文件，分别以序号标注交互与跳转的先后顺序：在 1 首页上传 `experiment_data` 中的 `EuAl4 Parent.cif`；在 2 上勾选畸变类型为 strains, displacive（Eu, Al）；在 2 上选择 method2（LD, K10 (0, 0, G)，g=1/6）；然后陆续点击几次 OK 依次进入 3、4、5、6；6 是选择导出结果的文件格式的导出页。`webpage_info/` 的内容**不允许修改**。
- `ISODISTORT_VALIDATE/` 用于比较 CIF 文件内容是否一致，通过比较 ISODISTORT 生成的结果与官网生成的结果是否一致，从而检查本地化是否有 bug，**允许修改**。
- `ISODISTORT/` 是 ISODISTORT 本地化主体，**允许修改**；其中 `isobyu/` 是从官网 https://iso.byu.edu/isotropy.php 下载的 Linux 二进制程序与数据库，**不允许修改**。ISODISTORT 有终端交互、网页交互、API 调用三种使用方法。网页与终端为英语界面。
- `ISOVIZ_INPUT/` 用于把振幅 CSV 与子群 `.isoviz` 一键写入 IsoVIZ 并打开查看畸变结构，**允许修改**。开发用的 `data.csv/` 与 `subgroup.isoviz/` **不上传远程仓库**。

---

## 目录结构

```text
ISOVIZ_INPUT/
  main.py                 一键写入并打开 IsoVIZ
  main_requirement.py     检查/创建 CRIS/.venv，只安装缺失的依赖
  requirements.txt        运行时依赖（当前仅标准库）
  requirements-dev.txt    开发依赖（pytest）
  pyproject.toml
  README.md
  isoviz_input/           Python 包：CSV 解析、写 .isoviz、启动 IsoVIZ
  tests/                  单元测试（使用 tests/fixtures/，不依赖本地样本）
  output/                 写入后的 .isoviz（不入库）
  data.csv/               本地开发用 CSV（git 忽略）
  subgroup.isoviz/        本地开发用官方 .isoviz（git 忽略）
```

---

## 安装

在 **CRIS 根目录**：

```powershell
python ISOVIZ_INPUT\main_requirement.py
```

也可以一次性准备三个子项目（同样使用 `CRIS/.venv`，已安装的包不会再下载）：

```powershell
python ISODISTORT\main_requirement.py
```

准备脚本会：

1. 确认 Python ≥ 3.10
2. 复用或创建 `CRIS/.venv`
3. 按包名检查，只 `pip install` 尚未安装的依赖
4. 检查 Java（IsoVIZ 是 Java 程序）以及根目录的 `ISOViz.lnk`
5. 若缺少 `ISOVIZ_INPUT/output/` 则自动新建

IsoVIZ 本体需要自行安装。推荐把快捷方式放到仓库根目录并命名为 `ISOViz.lnk`（该文件已在根 `.gitignore` 中），或设置环境变量 `ISOVIZ` / `ISOVIZ_JAR` 指向 `.exe` / `.jar`。Windows 上若 `.isoviz` 已关联 IsoVIZ，脚本会直接打开生成的文件。

---

## 使用

```powershell
cd <CRIS 根目录>
.\.venv\Scripts\python.exe ISOVIZ_INPUT\main.py --data <振幅.csv> --structure <子群.isoviz>
```

不传参数时会提示输入路径；若本地存在 `data.csv/` 或 `subgroup.isoviz/`，会先列出其中的文件供选择。

常用选项：

- `--output`：指定写出路径（默认 `ISOVIZ_INPUT/output/<原名>_patched.isoviz`，不覆盖输入文件）
- `--no-open`：只写文件，不启动 IsoVIZ

CSV 需含表头。`GD/save_best_model_parameters.py` 写出的列可直接用：

| 列 | 用途 |
| --- | --- |
| Mode Name | 与 `.isoviz` 中的 `modelabel` 匹配，例如 `[0,0,1/6]LD1[Eu1:a:dsp]A2u(a)` |
| Best Model Parameter | 写入 IsoVIZ 进度条使用的 `amp`（与 `maxamp` 同一单位，**不要**用 Normalized Amplitude） |
| Mode | 备选：`a1`, `a2`, ... 按 `.isoviz` 文件中 strain 模式再 displacive 模式的出现顺序 |

未匹配到的 IsoVIZ 模式保持原振幅（一般为 0）。

---

## 测试

```powershell
.\.venv\Scripts\python.exe -m pytest ISOVIZ_INPUT\tests
```
