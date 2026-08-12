"""与 main.py 菜单语义一致的 API 示例集合。"""

from isocore.api import IsoDistort


def example_search_and_single_mode(parent_cif: str = "parent.cif"):
    """示例 A：Search Page -> Distortion Page（单模式）。

    Relative path: examples/usage_examples.py
    """

    print("示例 A: Method 1 + Method 2 + single-mode")
    print("=" * 60)

    iso = IsoDistort()
    iso.load_structure(parent_cif)

    m1 = iso.search_method_1(distortion_types=["displacement", "strain"])
    if not m1:
        raise RuntimeError("Method 1 没有候选结果")

    m2 = iso.search_method_2(
        subgroup_idx=m1[0].subgroup.index,
        distortion_type="displacement",
    )
    if not m2.modes:
        raise RuntimeError("Method 2 没有返回模式")

    distorted = iso.generate_distortion(
        irrep_label=m2.modes[0].irrep_label,
        amplitude=1.0,
        supercell=[1, 1, 1],
    )
    iso.export("example_single_mode", formats=["cif"])
    print(f"完成，原子数: {len(distorted)}")
    return distorted


def example_search_and_mixed_mode(parent_cif: str = "parent.cif"):
    """示例 B：Search Page -> Distortion Page（多模式混合）。

    Relative path: examples/usage_examples.py
    """

    print("\n示例 B: Method 1 + Method 2 + mixed-mode")
    print("=" * 60)

    iso = IsoDistort()
    iso.load_structure(parent_cif)

    m1 = iso.search_method_1(distortion_types=["displacement", "strain"])
    if not m1:
        raise RuntimeError("Method 1 没有候选结果")

    m2 = iso.search_method_2(
        subgroup_idx=m1[0].subgroup.index,
        distortion_type="displacement",
    )
    labels = [m.irrep_label for m in m2.modes]
    if len(labels) < 2:
        raise RuntimeError("模式不足 2 个，无法演示 mixed-mode")

    mixed = iso.generate_mixed_distortion(
        contributions={labels[0]: 0.02, labels[1]: 0.03},
        supercell=[1, 1, 1],
    )
    iso.export("example_mixed_mode", formats=["cif"])
    print(f"完成，体积: {mixed.volume:.3f} Å^3")
    return mixed


def example_method_4_decomposition(parent_cif: str = "parent.cif", daughter_cif: str = "daughter.cif"):
    """示例 C：Method 4 模式分解。

    Relative path: examples/usage_examples.py
    """

    print("\n示例 C: Method 4 decomposition")
    print("=" * 60)

    iso = IsoDistort()
    iso.load_structure(parent_cif)

    m1 = iso.search_method_1(distortion_types=["displacement", "strain"])
    if not m1:
        raise RuntimeError("Method 1 没有候选结果")
    iso.search_method_2(subgroup_idx=m1[0].subgroup.index, distortion_type="displacement")

    result = iso.search_method_4(
        distorted_cif_path=daughter_cif,
        atom_matching_method="nearest-site",
        robust_distance_threshold=0.25,
    )
    print(f"RMS residual: {result.rms_residual:.8e}")
    return result


if __name__ == "__main__":
    print("isodistort - 使用示例（与 main.py 流程对齐）")
    print("=" * 60)
    print("请先准备 parent.cif / daughter.cif，或在函数参数中传入实际路径。")

    # 按需取消注释：
    # example_search_and_single_mode("parent.cif")
    # example_search_and_mixed_mode("parent.cif")
    # example_method_4_decomposition("parent.cif", "daughter.cif")
