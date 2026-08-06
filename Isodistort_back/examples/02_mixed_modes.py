"""
示例 2：多模式混合畸变

演示如何叠加多个不可约表示的畸变模式，
生成更复杂的畸变结构。
"""
from isocore.api import IsoDistort


def main():
    """Relative path: examples/02_mixed_modes.py"""
    
    iso = IsoDistort()

    # 加载结构
    iso.load_structure("../examples/sample_structures/NaCl.cif")

    # 选择路径
    iso.select_path(subgroup_idx=3, distortion_type="displacement")

    # 多模式混合：各模式指定不同幅度
    # {模式标号: 幅度系数}
    contributions = {
        "GM1+": 0.02,
        "GM4-": 0.03,
    }

    print("生成多模式混合畸变...")
    print(f"模式贡献: {contributions}")

    mixed = iso.generate_mixed_distortion(
        contributions=contributions,
        supercell=[2, 2, 2],
    )

    # 导出
    iso.export("NaCl_mixed_distortion", formats=["cif"])
    print(f"畸变后体积: {mixed.volume:.3f} Å³")


if __name__ == "__main__":
    main()
