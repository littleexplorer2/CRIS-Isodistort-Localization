"""示例 3：底层 API 直接调用（backend 层）。

跳过 `IsoDistort` 高层编排，直接调用：
- `FindsymWrapper`（结构识别）
- `IsoWrapper`（k 点枚举、子群枚举——Method 1 的底层能力）
"""
import sys
from pathlib import Path

# 允许直接以脚本方式运行（python examples/03_low_level_api.py）
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from isocore.backend import FindsymWrapper, IsoWrapper


def demo_find_sym() -> None:
    """演示 findsym 空间群识别"""
    print("=== Findsym 空间群识别示例 ===")
    findsym = FindsymWrapper()

    # NaCl 结构参数（F 心，常规胞 2 个原子）
    result = findsym.identify(
        lattice_params=[5.63, 5.63, 5.63, 90, 90, 90],
        atom_types=["Na", "Cl"],
        atom_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        centering="F",
    )

    print(f"空间群: #{result.space_group_number} ({result.space_group_symbol})")
    print(f"Wyckoff 位置数: {len(result.wyckoff_sites)}")
    for site in result.wyckoff_sites:
        print(f"  {site['wyckoff_letter']}: {len(site['atoms'])} 个等效点")


def demo_iso_subgroups() -> None:
    """演示 iso 子群枚举（Method 1 的底层数据源）"""
    print("\n=== ISO 子群枚举示例（SG 225, k=GM, IR=GM5-） ===")
    iso = IsoWrapper()

    kpoints = iso.list_k_points(225)
    print(f"k 点列表（{len(kpoints)} 个）：")
    print("  " + ", ".join(kp.label for kp in kpoints))

    irreps = iso.list_irreps(225, "GM")
    print(f"GM 点的不可约表示（{len(irreps)} 个）：")
    print("  " + ", ".join(f"{ir.label}({ir.dimension}D)" for ir in irreps))

    subgroups = iso.list_subgroups(225, "GM", "GM5-")
    print(f"\nGM5- 的各向同性子群（{len(subgroups)} 个）：")
    for sg in subgroups:
        print(
            f"  {sg.describe()} | index={sg.subgroup_index} "
            f"max={sg.is_maximal}"
        )


if __name__ == "__main__":
    demo_find_sym()
    demo_iso_subgroups()
