"""示例 1：与 Search/Distortion 页对齐的最小完整流程。

流程：
1) Method 1 过滤候选子群
2) Method 2 选择特定子群并计算模式
3) Distortion Page 生成单模式并导出

运行前请把 CIF_PATH 改为你的母相 CIF 文件路径。
"""

import sys
from pathlib import Path

# 允许直接以脚本方式运行（python examples/01_basic_workflow.py）
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from isocore.api import IsoDistort
from isocore.utils import IsodistortError

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

    # 依次尝试候选，找到第一个能产生位移模式的子群
    # （部分候选在当前结构的 Wyckoff 位置上没有位移模式，跳过即可）
    method2 = None
    for item in method1:
        try:
            res = iso.search_method_2(
                subgroup_idx=item.subgroup.index,
                distortion_type="displacement",
            )
        except (IsodistortError, ValueError, RuntimeError):
            continue
        if res.modes:
            method2 = res
            subgroup_idx = item.subgroup.index
            print(f"选中子群 subgroup_idx = {subgroup_idx} "
                  f"({item.subgroup.space_group_symbol})")
            break
    if method2 is None:
        raise RuntimeError("Method 1 候选均无可用位移模式")

    print("\n=== Step 3: Method 2 search ===")
    print(f"模式数: {len(method2.modes)}")

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
