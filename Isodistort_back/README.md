# isodistort - 本地 ISODISTORT 工作流

本项目基于 ISOTROPY Suite 的 `iso` / `findsym` 二进制，提供一个可本地运行的 Python 工作流，用于晶体结构畸变分析。

当前版本提供：
- 统一入口 `main.py`（终端交互式主程序）
- Python API `isodistort.api.IsoDistort`
- 对应 ISODISTORT Help 的四种搜索方法（Method 1-4）

## 交互流程总览

`main.py` 中的交互分为两页：

- Search Page
  - 设置 Distortion Types
  - Method 1: Search over all special k points
  - Method 2: General method - specific k points
  - Method 3: Arbitrary k + point-group/space-group + supercell
  - Method 4: Mode decomposition of distorted structure
- Distortion Page
  - 单模式生成
  - 多模式混合
  - 导出
  - 畴生成

对应 API 与示例脚本映射：

| main.py 菜单项 | API 方法 | 示例脚本 |
| --- | --- | --- |
| Method 1 | `search_method_1` | `examples/01_basic_workflow.py` |
| Method 2 | `search_method_2` | `examples/01_basic_workflow.py` |
| Method 3 | `search_method_3` | `examples/usage_examples.py` |
| Method 4 | `search_method_4` | `examples/usage_examples.py` |
| 单模式生成 | `generate_distortion` | `examples/01_basic_workflow.py` |
| 多模式混合 | `generate_mixed_distortion` | `examples/02_mixed_modes.py` |
| 导出 | `export` | `examples/01_basic_workflow.py` |
| 畴生成 | `generate_domains` | `main.py` 交互菜单 |

## 项目结构

```text
Isodistort_back/
├── main.py
├── pyproject.toml
├── requirements.txt
├── requirements-dev.txt
├── README.md
├── config/
│   └── settings.yaml
├── isobyu/                 # 官方二进制与数据库（只读）
├── isocore/                # 核心实现
├── isodistort/             # 对外包名兼容层
├── examples/
└── tests/
```

## 环境要求

- Python >= 3.9
- Windows 下建议安装 WSL（`isobyu` 为 Linux ELF 二进制）
- 依赖：`numpy`、`pymatgen`、`pyyaml`、`spglib`

## 快速开始

### 1. 安装依赖

```powershell
cd "c:\Users\devou\OneDrive\Desktop\CRIS\Isodistort_back"
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 2. 运行统一入口

```powershell
python main.py
```

随后按菜单逐步执行搜索与畸变生成。

## 示例脚本

- `examples/01_basic_workflow.py`
  - 演示 Method 1 -> Method 2 -> 单模式生成 -> 导出。
- `examples/02_mixed_modes.py`
  - 演示 Method 1 -> Method 2 -> 多模式混合。
- `examples/03_low_level_api.py`
  - 演示底层 `FindsymWrapper` 与 `IsoWrapper` 直接调用。
- `examples/usage_examples.py`
  - 汇总 Search/Distortion 语义一致的 API 示例函数。

注意：示例中的 `parent.cif` / `daughter.cif` 是占位路径，运行前请改为你的真实文件路径。

## Python API 示例

```python
from isodistort.api import IsoDistort

iso = IsoDistort()
iso.load_structure("parent.cif")

# Method 1
m1 = iso.search_method_1(
    distortion_types=["displacement", "strain"],
    crystal_system="tetragonal",
)

# Method 2
if m1:
    m2 = iso.search_method_2(
        subgroup_idx=m1[0].subgroup.index,
        distortion_type="displacement",
        k_point_label="X",
        k_point_coordinates=["1/2", "0", "0"],
        number_of_superposed_irs=1,
    )

# Distortion generation
iso.generate_distortion(amplitude=1.0, supercell=[1, 1, 1])
iso.export("distorted_output", formats=["cif", "poscar"])
```

## 配置文件

全局配置在 `config/settings.yaml`：
- `isobyu.bin_dir` / `isobyu.data_dir`
- `defaults.default_amplitude` / `defaults.default_supercell`
- `runtime.temp_dir` / `runtime.output_dir` / `runtime.timeout`

## 测试

```powershell
pip install -r requirements-dev.txt
pytest -q
```

## 说明

- `isobyu/` 下文件属于底层工具与数据库，不建议修改。
- 本地流程已提供四种 Method 的统一入口与可执行链路。
- 由于网页 GUI 和本地终端交互形式不同，界面表现不会完全一致，但功能路径与参数语义已对齐当前实现。
