# ISODISTORT 本地化项目

本目录把 [ISODISTORT](https://iso.byu.edu/isodistort.php)（BYU 晶体畸变在线计算）做到离线可用：用 ISOTROPY Suite 的 `iso` Linux 二进制与 `data_*.txt` 数据库，在本地复现 Search Page → Distortion Page 工作流。子群枚举、模式基矢、畴的计算结果与官网同一套数据库一致。

同一套计算引擎对外提供三种用法：**网页**、**终端菜单**、**Python API**。

---

## 仓库边界（哪些可以改）

CRIS 仓库按目录分工。改代码前请遵守：

| 路径 | 可否修改 | 说明 |
| --- | --- | --- |
| `ISODISTORT/`（本目录，不含下表例外） | 可以 | 本地化主体：计算、网页、终端、API、文档、测试 |
| `ISODISTORT/isobyu/` | **不可以** | 从 [iso.byu.edu](https://iso.byu.edu/isotropy.php) 下载的 Linux 二进制与数据库，只读 |
| `ISODISTORT_VALIDATE/` | 可以 | 比较本地 CIF 与官网 CIF，用来查回归 |
| `webpage_info/` | **不可以** | 官网各步 HTML 存档，文件名序号即交互顺序 |
| `实验数据与GD代码/` | **不可以** | 实验母相 CIF 等原始数据 |

官网交互对照（`webpage_info/`，勿改其中文件）：在 1 首页上传 `实验数据与GD代码/EuAl4 Parent.cif` → 2 勾选 strains 与 displacive（Eu, Al），Method 2 选 LD、K10 (0,0,g)、g=1/6 → 连续 OK 进入 3、4、5 → 6 为导出页（CIF / Save interactive distortion / Complete modes details / TOPAS.STR 等）。

---

## 一、能做什么

| 官网 | 本地 |
| --- | --- |
| 上传母相 CIF、设置 Distortion Types | 网页 / 终端 / API |
| Method 1：全部特殊 k 点子群 + 过滤 | 真实 `iso` 枚举 + 会话缓存 |
| Method 2：指定 k 点（可叠加 IR） | 枚举子群；点击行计算模式 |
| Method 3：点群/空间群 + 超胞 | 本地近似（见「已知差异」） |
| Method 4：畸变结构模式分解 | 最小二乘；支持超胞女儿相 |
| Distortion：幅度、Generate、Domains | 三种接口共用引擎 |
| Distortion：批量导出 | 网页 **只能选一个 Method** 的子群；勾选 CIF / ISOVIZ / Complete modes details / TOPAS.STR 打 ZIP |

底层子群/模式/畴由 `isobyu/iso` 完成。空间群识别用 pymatgen `SpacegroupAnalyzer`。`findsym` 不参与生产流程。

另外提供：中英双语（`isocore/i18n`）、五项物理自检、`tests_dev/` 与 `ISODISTORT_VALIDATE` 验证。

---

## 二、环境与安装

- Python ≥ 3.10（开发在 3.14 上验证）
- **Windows + WSL**（`isobyu/` 为 Linux ELF），或原生 **Linux**
- 依赖见 `requirements.txt`：`numpy`、`pymatgen`、`pyyaml`、`spglib`

在仓库根目录只建一份虚拟环境 `CRIS/.venv`：

```powershell
cd "C:\Users\devou\OneDrive\Desktop\CRIS"
python ISODISTORT\main_requirement.py
python ISODISTORT\main_requirement.py --dev    # 开发：pytest / ruff
```

脚本会检查 Python、WSL、`isobyu/iso` 与 `data_*.txt`；已有 `.venv` 则复用（`--recreate` 重建）。**不会**自动下载 ISOTROPY 套件。

把套件放进 `ISODISTORT/isobyu/`（必需 `iso` + `data_*.txt`）。该目录只读。Windows 下程序会在 WSL 家目录建短路径暂存和 `ISODATA` 符号链接，不必手改路径。

---

## 三、三种启动方式

都从 `ISODISTORT/` 运行，共用 `isocore`。

### A. 网页（推荐）

```powershell
python main_web.py
```

默认 `http://127.0.0.1:8000/`（`config/settings.yaml` → `runtime.web_port`；占用则顺延）。自动打开浏览器。右上角切换 English / 中文。关页后心跳超时（默认 60 s）自动停服。

页面自上而下对齐官网搜索页：

1. **Parent CIF**：上传母相。页头显示空间群、点阵、Wyckoff、Default space-group preferences。
2. **Types**：Strain 复选框 + Displacive / Occupational / Magnetic / Rotational 各行 all / none / 物种。复选框互不联动；点 **Change** 才提交。默认 Strain + Displacive 全物种（与 `webpage_info` 第 2 页一致）。
3. **Method 1**：晶系（多选 OR）→ 可达子群空间群 → Conventional/Primitive lattice → Maximal only。
4. **Method 2**：k 点；可改 superposed IR 数量后出现多组 k vector。OK 后列出子群（可筛选/下载表格）。点行计算模式。参数 k 点（LD/DT）本地往往只能枚举子群、不能算位移模式（见差异第 5 条）。枚举为空时可「生成本地子群数据库」或去官网。
5. **Method 3**：空间群或点群 + 实空间基矢（reciprocal 本地不支持）。
6. **Method 4**：上传女儿 CIF，得到模式幅度与 RMS。
7. **Distortion**：幅度、Generate、Domains、**按 Method 批量下载**（下一节）。
8. **Space-Group Preferences**：只读，固定国际标准取位。

### Distortion 批量下载（ZIP）

先在上面跑完 **某一个** Method（1 / 2 / 3），再在 Distortion 区下载。Method 4 是分解，没有「子群列表」可打包。

- **导出数据来源**：下拉菜单，**只能选一个** Method，不能多选。
- **格式**（可多选，对应官网第 6 页，手册 [modeparams](https://landau3.byu.edu/isodistorthelp.php#modeparams)）：
  - CIF file
  - Save interactive distortion
  - Complete modes details
  - TOPAS.STR
- 点 **Download all (ZIP)**。压缩包**只含该 Method 算出的子群**，不会把 `output/` 里历史文件打进去。

解压得到 `isodistort_methodN/`，其下按子群分文件夹，文件名「子群 + 格式」：

```text
isodistort_method2.zip
  isodistort_method2/
    LD1 C1/
      LD1 C1 CIF.cif
      LD1 C1 Save interactive distortion.isoviz
      LD1 C1 Complete modes details.txt
      LD1 C1 TOPAS.STR
    LD5 P6/
      ...
```

CIF 为该子群零振幅超胞（若当前子群已 Generate，则该子群的 CIF 用畸变结构）。参数 k 点上本地算不出位移模式时，模式类文件会注明限制。Generate 旁的单文件 CIF 仍走 `/api/download`，与 ZIP 无关。

### B. 终端

```powershell
python main_terminal.py
```

菜单：加载 CIF、Types、Method 1–4、Distortion（生成/导出/畴）、状态、语言、退出。

输入 CIF：在 `ISODISTORT/` 内搜 `.cif`（最多 30 个）；项目外的文件（如 `实验数据与GD代码/EuAl4 Parent.cif`）用「手动输入路径」，建议绝对路径。`output/` 与 `output/tmp/` 中的 CIF 不出现在列表里。

终端「导出」写到 `runtime.output_dir`（默认 `ISODISTORT/output/`），格式如 `cif,poscar`。按 Method 分文件夹的 ZIP 批量导出在**网页 Distortion**（或 API `export_subgroups`）。

### C. Python API

```python
from isocore.api import IsoDistort

iso = IsoDistort(language="en")
iso.load_structure(r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码\EuAl4 Parent.cif")
iso.set_distortion_scope({
    "displacive": ["*"], "occupational": [], "strain": [],
    "magnetic": [], "rotational": [],
})
iso.set_distortion_types(["strain", "displacive"])

m1 = iso.search_method_1(crystal_system="tetragonal")
iso.search_method_2(subgroup_idx=m1[0].subgroup.index)

iso.generate_distortion(amplitude=0.1)
iso.export("distorted", formats=["cif", "poscar"])          # 当前畸变结构
iso.export_subgroups("out_batch", formats=["cif", "isoviz", "modes", "topas"])
iso.generate_domains()
```

`export_subgroups` 使用当前 `iso.subgroups`（即最近一次 Method 1 / 2 k 点枚举 / Method 3 写入的列表）。网页 ZIP 则显式传入所选 Method 的列表，避免三种 Method 的结果混在一起。

语言：网页右上角；终端菜单第 9 项；API `iso.set_language("zh")`。

---

## 四、路径

| 用途 | 默认 | 配置 |
| --- | --- | --- |
| 终端搜输入 CIF | `ISODISTORT/` | 无；项目外请用绝对路径 |
| Generate / 终端导出 | `ISODISTORT/output/` | `settings.yaml` → `runtime.output_dir` |
| 上传暂存 | `ISODISTORT/output/tmp/` | `runtime.temp_dir` |
| 二进制与库 | `ISODISTORT/isobyu/` | 部署时放入，只读 |

`settings.yaml` 里的相对路径相对 `config/` 解析。网页 ZIP **不写入、不读取** `output/`。

---

## 五、目录结构

```text
ISODISTORT/
├── main_web.py / main_terminal.py / main_requirement.py
├── web/                 网页（server.py + index.html + static/）
├── isocore/             计算核心
│   ├── api/             IsoDistort
│   ├── backend/         iso / smodes / findsym 封装
│   ├── structure/       CIF 与超胞
│   ├── distortion/      Method 1–4、模式、畴
│   ├── io/              CIF / ISOVIZ / modes details / TOPAS.STR / POSCAR
│   ├── i18n/            中英文案与术语
│   └── utils/           配置、解析
├── isodistort/          包名别名（重导出 isocore）
├── config/settings.yaml
├── isobyu/              【只读】ISOTROPY 套件
├── tests_dev/           开发测试（生产可删）
└── output/              运行产物（不入库）
```

---

## 六、测试与验证

```powershell
python ISODISTORT\main_requirement.py --dev
cd ISODISTORT
python -m pytest tests_dev -q
ruff check .
```

依赖 WSL 的用例在 WSL 不可用时跳过。本地 CIF 与官网 CIF 的语义比较用仓库内 `ISODISTORT_VALIDATE/`（见该目录 README）。

主要层：金标准 SrTiO₃ 三条路径；API 与网页一致性；30-CIF；COD 外部结构；终端脚本驱动。批量导出的命名与「只含一个 Method」由 `tests_dev/test_distortion_export.py`、`test_web_server.py` 覆盖。

---

## 七、与官网的已知差异

1. **Types 过滤时机**：参数 k 点按 smodes + 物种过滤（EuAl4 LD g=1/6 + displacive 仅 LD1/LD2/LD5，与官网一致）；Method 1 全特殊 k 点仍先枚举再在模式阶段过滤。
2. **振幅**：官网 As/Ap + normfactor；本地为「最大模长归一化位移 × 用户幅度」，方向一致，数值待校准。
3. **模式列表**：官网含次级模式；本地 DISPLAY BUSH 以主 OPD 位移模式为主。
4. **应变模式未实现**：不改晶格参数；铁弹 a≠b 场景请对照官网。
5. **参数 k 点**：可枚举子群（首次可能要生成本地库，数分钟到数小时）；位移模式依赖官网 (3+d) 超空间，本地会明确提示。nmod 已从网页/终端移除。官网 a/b/g 会换算为 iso KVALUE（见 `kpoints_official.py`）。
6. **Method 3**：只覆盖特殊 k 点；reciprocal 不支持；带心仅默认 `d`。
7. **Method 4**：需先有模式基矢；女儿相为超胞时会提升坐标再分解。
8. **magnetic**：带 `m` 前缀的 IR 不进入默认流程。
9. **occupational**：本地为 ±1 占据近似，校验失败会标明。
10. **Method 1 lattice 下拉**：由本地枚举基矢去重生成，数量可能与官网数据库略有差别（版本差异，界面有告示）。
11. **多物种共享 Wyckoff 字母**：BUSH 无法分物种时不静默给错结构。
12. **多维模式**：各分量等权求和再归一化（OPD 通用方向近似）。
13. **取位约定**：CIF 位点与 iso 符号点对不齐时会标明，不静默当正确结果。

ISOVIZ / Complete modes details / TOPAS.STR 按官网选项语义在本地生成，不是官网站点 PHP 的逐字节拷贝。

---

## 八、术语与常见问题

术语：英文以官网与 ITC 为准；中文以《晶体学名词》为准。对照表 `isocore/i18n/terms.py`。

- **WSL 报错**：`wsl --status`，需有默认发行版。
- **Method 1 慢**：首次枚举全部特殊 k 点，属正常；之后会话缓存秒回。
- **网页打不开**：确认 `main_web.py` 已启动，访问控制台打印的 URL。
- **关页即停服**：预期行为。要常驻就别关标签页。
- **isobyu 缺失**：把套件放入 `ISODISTORT/isobyu/`。
- **ZIP 报没有子群**：先跑下拉框里选中的那一个 Method，再下载。
- **ZIP 里文件特别多、名字对不上子群**：旧版曾打包整个 `output/`；当前版本只含所选 Method 的子群文件夹。
