"""示例 2：与 Distortion Page 对齐的多模式混合。

流程：
1) Method 1 + Method 2 获取可用模式
2) 从返回模式中自动选择前两个标签
3) 生成 mixed distortion 并导出
"""

from isocore.api import IsoDistort


def main():
    """Relative path: examples/02_mixed_modes.py"""

    iso = IsoDistort()
    iso.load_structure("parent.cif")

    method1 = iso.search_method_1(distortion_types=["displacement", "strain"])
    if not method1:
        raise RuntimeError("Method 1 未返回候选子群")

    method2 = iso.search_method_2(
        subgroup_idx=method1[0].subgroup.index,
        distortion_type="displacement",
    )
    mode_labels = [m.irrep_label for m in method2.modes]
    if len(mode_labels) < 2:
        raise RuntimeError("可用模式少于 2 个，无法演示混合模式")

    contributions = {
        mode_labels[0]: 0.02,
        mode_labels[1]: 0.03,
    }

    print("生成多模式混合畸变...")
    print(f"模式贡献: {contributions}")

    mixed = iso.generate_mixed_distortion(
        contributions=contributions,
        supercell=[1, 1, 1],
    )

    iso.export("mixed_distorted_output", formats=["cif"])
    print(f"畸变后体积: {mixed.volume:.3f} Å^3")


if __name__ == "__main__":
    main()
