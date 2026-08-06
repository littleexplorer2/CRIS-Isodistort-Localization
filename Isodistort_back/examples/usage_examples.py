"""
isodistort 使用示例

演示如何使用当前公开的 Python API 进行畸变计算。
"""

from isocore.api import IsoDistort


def example_basic():
    """基本用法：从 CIF 加载结构并生成单模式畸变。
    
    Relative path: examples/usage_examples.py"""

    print("示例1: 基本工作流")
    print("=" * 50)

    iso = IsoDistort()
    iso.load_structure("SrTiO3.cif")
    iso.list_subgroups()
    iso.select_path(subgroup_idx=3, distortion_type="displacement")
    distorted = iso.generate_distortion(amplitude=0.05)
    iso.export("SrTiO3_distorted", formats=["cif"])

    print(f"完成，共 {len(distorted)} 个原子")
    return distorted


def example_mixed_modes():
    """基本多模式叠加示例。
    
    Relative path: examples/usage_examples.py"""

    print("\n\n示例2: 多模式混合")
    print("=" * 50)

    iso = IsoDistort()
    iso.load_structure("SrTiO3.cif")
    iso.select_path(subgroup_idx=3, distortion_type="displacement")

    mixed = iso.generate_mixed_distortion(
        contributions={"GM1+": 0.02, "GM4-": 0.03},
        supercell=[2, 2, 2],
    )
    iso.export("SrTiO3_mixed_distortion", formats=["cif"])

    print(f"畸变后体积: {mixed.volume:.3f} Å³")
    return mixed


if __name__ == "__main__":
    print("isodistort - 使用示例")
    print("=" * 60)
    print()
    print("注意：运行前请准备好 CIF 文件路径")
    print()

    # 取消注释运行对应示例
    # example_basic()
    # example_mixed_modes()
