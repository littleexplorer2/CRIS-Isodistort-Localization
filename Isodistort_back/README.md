# isodistort - 晶体畸变分析工具

基于 ISOTROPY Software Suite（iso + findsym）的离线封装，实现 ISODISTORT 核心的结构畸变分析能力。

## 项目架构

### 程序完整框架
```
isodistort_back/
├── isobyu/              ISOTROPY Linux 二进制程序文件与数据库（只读，不可修改）
├── isodistort/          对外公开包名（兼容入口）
├── isocore/             项目核心代码框架
├── config/              配置文件（settings.yaml）
├── examples/            示例脚本
├── tests/               单元测试
├── output/              默认输出目录
├── pyproject.toml       项目配置
├── output/              默认输出目录
├── main.py/             项目主入口
├── requirement.txt/     客户端依赖清单
├── requirement-dev.txt/ 开发阶段依赖清单
└── README.md/           用户使用说明
```

### 项目核心代码框架
```
isocore/
├── backend/        【底层封装层】封装 iso / findsym 二进制（纯群论计算）
├── structure/      【晶体结构层】CIF 读写、对称处理、坐标变换、位点匹配
├── distortion/     【畸变计算层】相变路径、模式映射、畸变引擎、畴生成（核心）
├── io/             【结果输出层】结构文件导出、结果序列化
├── api/            【接口层】Python API（不提供独立命令行脚本）
└── utils/          【工具层】配置、异常、文本解析
```

### 对应 ISODISTORT 12 步骤的模块分布

| 阶段 | 步骤 | 模块 | 实现方式 |

一、结构输入与对称识别     
1. 读取 CIF             | structure.cif_io            | ❌ 自研（pymatgen）
2. 识别空间群与 Wyckoff | backend.findsym_wrapper      | ✅ 封装 findsym 
3. 对称性校验           | structure.symmetry_validator | ⚖️ 混合

二、子群枚举              
4. 枚举各向同性子群      | backend.iso_wrapper         | ✅ 封装 iso
5. 相变路径参数组装      | distortion.phase_path       | ❌ 自研

三、畸变模式分解 
6. 计算畸变模式基矢      | backend.iso_wrapper         | ✅ 封装 iso
7. Wyckoff 位点分裂     | structure.site_mapping      | ⚖️ 混合

四、畸变映射 
8. 基矢→原子坐标映射     | distortion.distortion_mapper| ❌ 自研（核心难点）

五、畸变生成 
9. 幅度缩放与多模式混合  | distortion.distortion_engine | ❌ 自研
10. 畴变体生成          | distortion.domain_generator  | ⚖️ 混合

六、结果输出
11. 结构文件导出        | io.structure_exporter        | ❌ 自研（pymatgen）

七、可视化  
12. 3D 结构可视化       | 外部软件VESTA                 | ✅ 已有



## 环境要求

- Python >= 3.9
- **WSL2**：isobyu 中的二进制为 Linux ELF 格式
  - Windows 用户需启用 WSL，代码自动检测并通过 WSL 调用
- 依赖：numpy, pymatgen, pyyaml
  - 可视化：导出 CIF 并使用外部软件（如 VESTA）查看



## 使用方法

### 1. 在Windows Powershell中创建虚拟环境并安装依赖

```powershell
cd "c:\Users\devou\OneDrive\Desktop\CRIS\Isodistort_back"（修改成自己电脑上的路径）
python -m venv .venv  # 创建虚拟环境
.\.venv\Scripts\Activate.ps1  # 激活虚拟环境

# 在虚拟环境中安装运行依赖
pip install -r requirements.txt

# 运行示例脚本（完整工作流：加载 CIF -> 枚举子群 -> 生成畸变 -> 导出 CIF）
python examples\01_basic_workflow.py

# 可选做：运行单元测试
pip install pytest
pytest -q
```

### 2. Python API 使用

> 注意：本项目**不提供独立的命令行可执行脚本**；请通过 Python API 或运行示例脚本来使用。

```python
from isodistort.api import IsoDistort

iso = IsoDistort()

# 加载结构，修改为你的cif文件的绝对路径
iso.load_structure("your_cif_path.cif")

# 枚举子群
iso.list_subgroups()

# 选择相变路径（子群序号默认为3 + 畸变类型默认为displacement）
iso.select_path(subgroup_idx=3, distortion_type="displacement")

# 生成畸变结构（畸变幅度默认为0.05）
iso.generate_distortion(amplitude=0.05)

# 导出为cif或poscar格式文件
iso.export("distorted_output", formats=["cif", "poscar"])
```

### 3. 查看导出 CIF（使用 VESTA）

导出后的 CIF 文件位于配置的输出目录（默认 `output/`），请使用VESTA等第三方软件打开



## 核心复用 vs 自研

### ✅ 直接复用 isobyu 的能力（100% 对标官方）
- 空间群与 Wyckoff 位点识别（findsym）
- 各向同性子群全量枚举与相变路径（iso）
- 四类畸变模式的群论分解与基矢计算（iso）
- 位点对称性分裂的理论计算（iso）
- 畴变体对称操作矩阵（iso）
- 全量晶体对称数据库（data_*.txt）

### ❌ 需要自研的部分
- CIF 结构读写与位点匹配
- 畸变基矢到原子坐标的映射（核心难点，决定结果一致性）
- 交互界面与参数校验
- 结构文件导出（CIF/POSCAR 等）



## 注意事项

1. **isobyu 文件夹为只读**：官方二进制与数据库，请勿修改
2. **一致性对标**：畸变映射逻辑需与在线版 ISODISTORT 反复对标调试
   - 坐标变换约定、原点选择、基矢归一化是主要误差来源
3. **Windows 环境**：自动通过 WSL 调用 Linux 二进制，需确保已安装 WSL
4. **高阶功能**：smodes、comsubs 暂未封装，架构预留扩展位



## 开发路线

- [x] 项目骨架与分层架构
- [x] backend 层：iso / findsym 封装
- [x] structure 层：CIF 读写、对称校验、坐标变换、位点匹配
- [x] distortion 层：相变路径、畸变映射、畸变引擎、畴生成
- [x] io 层：结构导出、结果序列化
- [x] vis 层：基础可视化（已移除，使用外部 VESTA 或其他工具查看导出 CIF）
- [x] api 层：Python API + CLI
- [x] api 层：Python API（无独立 CLI）
- [ ] 畸变映射与在线版一致性对标（核心待完善）
- [ ] smodes / comsubs 扩展封装