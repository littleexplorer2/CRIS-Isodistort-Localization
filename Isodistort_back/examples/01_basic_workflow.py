"""示例 1：与 Search/Distortion 页对齐的最小完整流程。

流程：
1) Method 1 过滤候选子群
2) Method 2 选择特定子群并计算模式
3) Distortion Page 生成单模式并导出

运行前请把 CIF_PATH 改为你的母相 CIF 文件路径。
"""

from isocore.api import IsoDistort

# 请替换为你自己的 parent CIF 文件路径
CIF_PATH = r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码\EuAl4 Springer (parent).cif"


def main() -> None:
    iso = IsoDistort()

    print("=== Step 1: Load parent structure ===")
    iso.load_structure(CIF_PATH)

    print("\n=== Step 2: Method 1 search ===")
    method1 = iso.search_method_1(
        distortion_types=["displacement", "strain"],
        crystal_system=None,
        maximal_subgroup_only=False,
    )
    if not method1:
        raise RuntimeError("Method 1 未返回候选子群")

    subgroup_idx = method1[0].subgroup.index
    print(f"Method 1 首个候选 subgroup_idx = {subgroup_idx} "
          f"({method1[0].subgroup.space_group_symbol})")

    print("\n=== Step 3: Method 2 search ===")
    method2 = iso.search_method_2(
        subgroup_idx=subgroup_idx,
        distortion_type="displacement",
    )
    if not method2.modes:
        raise RuntimeError("Method 2 未返回可用模式（可尝试 Method 1 的其他候选）")

    print("\n=== Step 4: Generate distortion ===")
    distorted = iso.generate_distortion(
        irrep_label=method2.modes[0].irrep_label,
        amplitude=0.1,
    )
    print(f"生成完成，原子数 = {len(distorted)}")

    print("\n=== Step 5: Export ===")
    iso.export("distorted_output", formats=["cif", "poscar"])
    print("导出完成。")


if __name__ == "__main__":
    main()
