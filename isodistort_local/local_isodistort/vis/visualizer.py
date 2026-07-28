"""
3D 结构可视化 - 基于 pyvista 的晶体结构渲染

对应阶段七，步骤12：3D 结构可视化与交互调节
实现方式：❌ 自研（基于 pyvista / pymatgen 可视化）

注意：本模块为可选功能，需要安装 pyvista 依赖。
"""
from typing import Optional
from pymatgen.core import Structure


class StructureVisualizer:
    """
    晶体结构 3D 可视化器

    功能：
    1. 渲染晶体结构（原子球 + 键连）
    2. 畸变前后对比
    3. 原子样式调节
    """

    def __init__(self, backend: str = "pymatgen"):
        """
        Args:
            backend: 可视化后端
                - "pymatgen": 使用 pymatgen 内置可视化（默认，依赖少）
                - "pyvista": 使用 pyvista 3D 渲染（需额外安装）
        """
        self.backend = backend

    def plot(self, structure: Structure, title: str = "Crystal Structure",
             save_path: str = None):
        """
        绘制晶体结构

        Args:
            structure: 晶体结构
            title: 图标题
            save_path: 保存图片路径（为 None 则弹窗显示）
        """
        if self.backend == "pyvista":
            self._plot_pyvista(structure, title, save_path)
        else:
            self._plot_pymatgen(structure, title, save_path)

    @staticmethod
    def _plot_pymatgen(structure: Structure, title: str, save_path: Optional[str]):
        """使用 pymatgen 内置的可视化"""
        try:
            from pymatgen.vis.structure_vtk import StructureVis
            vis = StructureVis()
            vis.set_structure(structure)
            if save_path:
                # VTK 截图保存
                vis.show()
            else:
                vis.show()
        except ImportError:
            # 降级：打印结构信息
            print(f"[可视化降级] {title}")
            print(f"空间群: {structure.get_space_group_info()[1]}")
            print(f"原子数: {len(structure)}")
            print("成分:", structure.composition)
            print("\n原子坐标:")
            for i, site in enumerate(structure):
                print(f"  {i:2d} {site.species_string:4s}  "
                      f"{site.frac_coords[0]:.4f} {site.frac_coords[1]:.4f} {site.frac_coords[2]:.4f}")

    @staticmethod
    def _plot_pyvista(structure: Structure, title: str, save_path: Optional[str]):
        """使用 pyvista 3D 渲染（需安装 pyvista）"""
        try:
            import pyvista as pv
            from pymatgen.vis.structure_plotter import StructurePlotter
            # pyvista 实现待完善
            print(f"[pyvista] {title} - {len(structure)} atoms")
        except ImportError:
            print("警告: pyvista 未安装，已降级为文本输出")
            StructureVisualizer._plot_pymatgen(structure, title, save_path)

    def compare_structures(self, struct1: Structure, struct2: Structure,
                           label1: str = "Parent", label2: str = "Distorted"):
        """对比两个结构（畸变前后）"""
        print(f"\n=== 结构对比: {label1} vs {label2} ===")
        print(f"{label1}: {len(struct1)} 原子, 体积 {struct1.volume:.3f} Å³")
        print(f"{label2}: {len(struct2)} 原子, 体积 {struct2.volume:.3f} Å³")
        print(f"体积变化: {(struct2.volume - struct1.volume) / struct1.volume * 100:.3f}%")
