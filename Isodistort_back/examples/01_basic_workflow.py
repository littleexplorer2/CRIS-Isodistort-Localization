"""
示例 1：完整工作流 - 从 CIF 到畸变结构导出

演示 IsoDistort API 的标准使用流程：
1. 加载 CIF 结构
2. 枚举所有各向同性子群
3. 选择相变路径
4. 生成畸变结构
5. 导出 CIF / POSCAR
"""
from isocore.api import IsoDistort


def main():
    # 初始化
    iso = IsoDistort()

    # ---- 步骤 1：加载结构 ----
    # 替换为你的 CIF 文件路径
    cif_path = "../examples/sample_structures/NaCl.cif"
    print("=== 步骤 1：加载结构 ===")
    iso.load_structure(cif_path)

    # ---- 步骤 2：枚举子群 ----
    print("\n=== 步骤 2：枚举各向同性子群 ===")
    subgroups = iso.list_subgroups(distortion_type="displacement")

    # ---- 步骤 3：选择相变路径 ----
    # 选择第 3 个子群（示例，实际根据需要选择）
    print("\n=== 步骤 3：选择相变路径 ===")
    iso.select_path(subgroup_idx=3, distortion_type="displacement")

    # ---- 步骤 4：生成畸变结构 ----
    print("\n=== 步骤 4：生成畸变结构 ===")
    distorted = iso.generate_distortion(
        irrep_label=None,  # None 表示使用第一个模式
        amplitude=0.05,    # 畸变幅度
        supercell=[1, 1, 1],
    )

    # ---- 步骤 5：导出 ----
    print("\n=== 步骤 5：导出结构 ===")
    iso.export("NaCl_distorted", formats=["cif", "poscar"])

    # ---- 步骤 6：可视化对比 (内置可视化已移除) ----
    print("\n=== 步骤 6：结构对比 ===")
    print("已在输出目录导出 CIF 文件，使用 VESTA 或其他工具打开以查看畸变结构。")


if __name__ == "__main__":
    main()
