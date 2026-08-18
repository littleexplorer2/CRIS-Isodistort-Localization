"""与 main_terminal.py 菜单语义一致的 API 示例集合。

运行前请把 CIF_PATH（与可选的 DAUGHTER_CIF）改为你的真实文件路径。
"""

import sys
from pathlib import Path

# 允许直接以脚本方式运行（python examples/usage_examples.py）
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from isocore.api import IsoDistort
from isocore.utils import IsodistortError

CIF_PATH = r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码\EuAl4 Springer (parent).cif"
DAUGHTER_CIF = r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码\LD1_C1_subgroup.cif"


def _first_subgroup_with_modes(iso: IsoDistort, m1) -> "IsoDistort":
    """返回第一个能产生位移模式的 Method 2 结果。

    部分候选在当前结构的 Wyckoff 位置上没有位移模式（官网同样只显示
    “允许位移的 IR”），逐个尝试并跳过即可。
    """
    for item in m1:
        try:
            res = iso.search_method_2(
                subgroup_idx=item.subgroup.index,
                distortion_type="displacement",
            )
        except (IsodistortError, ValueError, RuntimeError):
            continue
        if res.modes:
            return res
    raise RuntimeError("Method 1 候选均无法产生位移模式")


def example_search_and_single_mode() -> None:
    """示例 A：Search Page -> Distortion Page（单模式）。"""
    print("示例 A: Method 1 + Method 2 + single-mode")
    print("=" * 60)

    iso = IsoDistort()
    iso.load_structure(CIF_PATH)

    m1 = iso.search_method_1(distortion_types=["displacement", "strain"])
    if not m1:
        raise RuntimeError("Method 1 没有候选结果")

    m2 = _first_subgroup_with_modes(iso, m1)
    distorted = iso.generate_distortion(
        irrep_label=m2.modes[0].irrep_label,
        amplitude=0.1,
    )
    iso.export("example_single_mode", formats=["cif"])
    print(f"完成，原子数: {len(distorted)}")


def example_search_and_mixed_mode() -> None:
    """示例 B：Search Page -> Distortion Page（多模式混合）。"""
    print("\n示例 B: Method 1 + Method 2 + mixed-mode")
    print("=" * 60)

    iso = IsoDistort()
    iso.load_structure(CIF_PATH)

    m1 = iso.search_method_1(distortion_types=["displacement", "strain"])
    if not m1:
        raise RuntimeError("Method 1 没有候选结果")

    m2 = _first_subgroup_with_modes(iso, m1)
    labels = [m.irrep_label for m in m2.modes]
    if len(labels) < 2:
        raise RuntimeError("模式不足 2 个，无法演示 mixed-mode")

    mixed = iso.generate_mixed_distortion(
        contributions={labels[0]: 0.02, labels[1]: 0.03},
    )
    iso.export("example_mixed_mode", formats=["cif"])
    print(f"完成，体积: {mixed.volume:.3f} A^3")


def example_method_4_decomposition() -> None:
    """示例 C：Method 4 模式分解（需先获得模式基矢）。"""
    print("\n示例 C: Method 4 decomposition")
    print("=" * 60)

    iso = IsoDistort()
    iso.load_structure(CIF_PATH)

    m1 = iso.search_method_1(distortion_types=["displacement", "strain"])
    if not m1:
        raise RuntimeError("Method 1 没有候选结果")
    _first_subgroup_with_modes(iso, m1)

    result = iso.search_method_4(
        distorted_cif_path=DAUGHTER_CIF,
        atom_matching_method="nearest-site",
        robust_distance_threshold=0.25,
    )
    print(f"RMS residual: {result.rms_residual:.8e}")
    ranked = sorted(result.amplitudes.items(), key=lambda kv: abs(kv[1]), reverse=True)
    for label, amp in ranked[:10]:
        print(f"  {label:<12s} {amp:+.8f}")


if __name__ == "__main__":
    print("isodistort - 使用示例（与 main_terminal.py 流程对齐）")
    print("=" * 60)

    # 按需取消注释：
    # example_search_and_single_mode()
    # example_search_and_mixed_mode()
    # example_method_4_decomposition()
