"""示例 1：与 Search Page 对齐的最小完整流程。

对应主入口 `main.py` 的菜单路径：
1) Method 1 过滤候选子群
2) Method 2 选择特定子群并计算模式
3) Distortion Page 生成单模式并导出
"""

from isocore.api import IsoDistort


def main():
    """Relative path: examples/01_basic_workflow.py"""

    iso = IsoDistort()

    # 请替换为你自己的 parent CIF 文件路径
    cif_path = "parent.cif"
    print("=== Step 1: Load parent structure ===")
    iso.load_structure(cif_path)

    print("\n=== Step 2: Method 1 search ===")
    method1 = iso.search_method_1(
        distortion_types=["displacement", "strain"],
        crystal_system=None,
        maximal_subgroup_only=False,
    )
    if not method1:
        raise RuntimeError("Method 1 未返回候选子群")

    subgroup_idx = method1[0].subgroup.index
    print(f"Method 1 首个候选 subgroup_idx = {subgroup_idx}")

    print("\n=== Step 3: Method 2 search ===")
    method2 = iso.search_method_2(
        subgroup_idx=subgroup_idx,
        distortion_type="displacement",
        k_point_label=None,
        k_point_coordinates=None,
        number_of_superposed_irs=1,
    )
    if not method2.modes:
        raise RuntimeError("Method 2 未返回可用模式")

    print("\n=== Step 4: Generate distortion ===")
    distorted = iso.generate_distortion(
        irrep_label=method2.modes[0].irrep_label,
        amplitude=1.0,
        supercell=[1, 1, 1],
    )
    print(f"生成完成，原子数 = {len(distorted)}")

    print("\n=== Step 5: Export ===")
    iso.export("distorted_output", formats=["cif", "poscar"])
    print("导出完成。")


if __name__ == "__main__":
    main()
