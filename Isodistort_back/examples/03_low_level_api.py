"""示例 3：底层 API 直接调用（backend 层）。

该示例跳过 `IsoDistort` 高层编排，直接调用：
- `FindsymWrapper`（结构识别）
- `IsoWrapper`（子群搜索，接近 Method 1 的底层能力）
"""
from isocore.backend import IsoWrapper, FindsymWrapper


def demo_find_sym():
    """演示 findsym 空间群识别

    Relative path: examples/03_low_level_api.py"""

    print("=== Findsym 空间群识别示例 ===")
    findsym = FindsymWrapper()

    # NaCl 结构参数
    result = findsym.identify(
        lattice_params=[5.63, 5.63, 5.63, 90, 90, 90],
        atom_types=["Na", "Cl"],
        atom_positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
        centering="F",
    )

    print(f"空间群: #{result.space_group_number} ({result.space_group_symbol})")
    print(f"Wyckoff 位点数: {len(result.wyckoff_sites)}")
    for site in result.wyckoff_sites:
        print(f"  {site['wyckoff_letter']}: {len(site['atoms'])} 个等效点")


def demo_list_subgroups():
    """演示 iso 子群枚举

    Relative path: examples/03_low_level_api.py"""

    print("\n=== ISO 子群枚举示例（对应 Method 1 的底层子群候选） ===")
    iso = IsoWrapper()

    subgroups = iso.list_subgroups(
        parent_sg=225,
        distortion_type="displacement",
    )

    print(f"共 {len(subgroups)} 个各向同性子群:")
    for sg in subgroups:
        print(f"  #{sg.index:2d}  SG #{sg.space_group_number}  {sg.space_group_symbol}")


if __name__ == "__main__":
    demo_find_sym()
    demo_list_subgroups()
