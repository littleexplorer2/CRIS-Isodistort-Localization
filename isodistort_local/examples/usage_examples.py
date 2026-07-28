"""
ISODISTORT Local 使用示例

演示如何使用Python API进行批量畸变模式计算
"""

from isodistort_local import (
    ISOTROPYConfig,
    DistortionConfig,
    BatchDistortionCalculator,
    ResultExporter,
)


def example_basic():
    """基本用法：计算单个CIF的所有畸变模式"""
    
    print("示例1: 基本批量计算")
    print("=" * 50)
    
    # 1. 配置ISOTROPY路径
    iso_config = ISOTROPYConfig(
        isotropy_path="/path/to/isobyu",  # 替换为你的ISOTROPY路径
        # use_wsl=True,  # Windows下使用WSL时开启
    )
    
    # 2. 配置畸变计算参数
    dist_config = DistortionConfig(
        k_points=["GM", "X", "M", "R"],  # 计算这些k点
        distortion_types=["displacive", "strain"],  # 位移+应变
        method=1,  # 方法1：枚举子群
        mode_amplitude=1.0,
    )
    
    # 3. 创建计算器
    calculator = BatchDistortionCalculator(
        iso_config=iso_config,
        dist_config=dist_config,
        output_dir="./output",
    )
    
    # 4. 执行计算
    result = calculator.calculate_from_cif("SrTiO3.cif")
    
    # 5. 查看结果
    print(f"\n计算完成！")
    print(f"  总模式数: {result.total_modes()}")
    print(f"  k点数量: {len(result.results_by_kpoint)}")
    
    for k, v in result.results_by_kpoint.items():
        print(f"  k={k}: {len(v.modes)} 个模式")
    
    # 6. 导出结果
    exporter = ResultExporter(output_dir="./output")
    exported = exporter.export_all(result)
    
    print(f"\n导出文件:")
    for fmt, path in exported.items():
        print(f"  {fmt}: {path}")
    
    return result


def example_smodes():
    """使用SMODES计算位移模式"""
    
    print("\n\n示例2: SMODES位移模式计算")
    print("=" * 50)
    
    from isodistort_local.input import CrystalStructure, InputGenerator
    
    # 创建结构（或从CIF加载）
    structure = CrystalStructure(
        title="SrTiO3",
        space_group_number=221,
        lattice_params=(3.905, 3.905, 3.905, 90, 90, 90),
        wyckoff_sites=[
            {"symbol": "Sr", "letter": "a", "x": 0, "y": 0, "z": 0},
            {"symbol": "Ti", "letter": "b", "x": 0.5, "y": 0.5, "z": 0.5},
            {"symbol": "O", "letter": "c", "x": 0.5, "y": 0.5, "z": 0},
        ],
    )
    
    iso_config = ISOTROPYConfig(isotropy_path="/path/to/isobyu")
    dist_config = DistortionConfig(k_points=["GM", "R", "M"])
    
    calculator = BatchDistortionCalculator(
        iso_config=iso_config,
        dist_config=dist_config,
    )
    
    # 使用SMODES计算
    smodes_result = calculator.calculate_with_smodes(structure)
    
    print(f"找到 {len(smodes_result.get('modes', []))} 个模式")
    
    return smodes_result


def example_from_pymatgen():
    """从pymatgen结构对象开始"""
    
    print("\n\n示例3: 从pymatgen结构对象开始")
    print("=" * 50)
    
    try:
        from pymatgen.core.structure import Structure
        from pymatgen.io.cif import CifParser
    except ImportError:
        print("需要安装pymatgen: pip install pymatgen")
        return
    
    # 读取CIF
    parser = CifParser("structure.cif")
    struct = parser.get_structures()[0]
    
    # 配置和计算
    iso_config = ISOTROPYConfig(isotropy_path="/path/to/isobyu")
    dist_config = DistortionConfig(k_points=["GM", "X"])
    
    calculator = BatchDistortionCalculator(
        iso_config=iso_config,
        dist_config=dist_config,
    )
    
    # 从pymatgen结构计算
    from isodistort_local.input import InputGenerator
    crystal_struct = InputGenerator.structure_from_pymatgen(struct)
    result = calculator.calculate(crystal_struct)
    
    print(f"完成，共 {result.total_modes()} 个模式")
    
    return result


if __name__ == "__main__":
    print("ISODISTORT Local - 使用示例")
    print("=" * 60)
    print()
    print("注意：运行示例前请修改ISOTROPY路径")
    print("下载ISOTROPY: https://iso.byu.edu/isolinux.php")
    print()
    
    # 取消注释运行对应示例
    # example_basic()
    # example_smodes()
    # example_from_pymatgen()
