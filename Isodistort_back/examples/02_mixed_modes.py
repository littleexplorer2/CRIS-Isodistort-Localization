"""示例 2：与 Distortion Page 对齐的多模式混合。

流程：
1) Method 1 + Method 2 获取可用模式
2) 从返回模式中自动选择前两个标签
3) 生成 mixed distortion 并导出

运行前请把 CIF_PATH 改为你的母相 CIF 文件路径。
"""

from isocore.api import IsoDistort
from isocore.utils import IsodistortError

CIF_PATH = r"C:\Users\devou\OneDrive\Desktop\CRIS\实验数据与GD代码\EuAl4 Springer (parent).cif"


def main() -> None:
    iso = IsoDistort()
    iso.load_structure(CIF_PATH)

    method1 = iso.search_method_1(distortion_types=["displacement", "strain"])
    if not method1:
        raise RuntimeError("Method 1 未返回候选子群")

    # 依次尝试候选，找到第一个能产生 ≥2 个模式的子群
    # （部分候选在当前结构的 Wyckoff 位点上没有位移模式，跳过即可）
    method2 = None
    for item in method1:
        try:
            res = iso.search_method_2(
                subgroup_idx=item.subgroup.index,
                distortion_type="displacement",
            )
        except (IsodistortError, ValueError, RuntimeError):
            continue
        if len(res.modes) >= 2:
            method2 = res
            break
    if method2 is None:
        raise RuntimeError("未找到可产生两个以上模式的子群")

    mode_labels = [m.irrep_label for m in method2.modes]
    contributions = {
        mode_labels[0]: 0.02,
        mode_labels[1]: 0.03,
    }

    print("生成多模式混合畸变...")
    print(f"模式贡献: {contributions}")

    mixed = iso.generate_mixed_distortion(
        contributions=contributions,
    )

    iso.export("mixed_distorted_output", formats=["cif"])
    print(f"畸变后体积: {mixed.volume:.3f} A^3")


if __name__ == "__main__":
    main()
