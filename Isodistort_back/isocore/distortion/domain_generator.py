"""
畴变体生成器 - 基于对称操作生成各畴结构

对应阶段五，步骤10：畴变体生成与切换
实现方式：⚖️ 混合实现（对称操作复用 iso，坐标变换自研）
"""
import numpy as np
from pymatgen.core import Structure
from typing import List

from ..backend import IsoWrapper
from ..structure import apply_affine_transform, wrap_to_unit_cell
from .phase_path import PhasePath


class DomainGenerator:
    """
    畴变体生成器

    功能：
    1. 获取相变对应的所有畴变体对称操作
    2. 对畸变结构应用对称操作，生成各畴变体
    3. 畴切换支持
    """

    def __init__(self, iso_wrapper: IsoWrapper = None):
        """Relative path: isocore/distortion/domain_generator.py"""
        
        self.iso = iso_wrapper or IsoWrapper()

    def get_domain_count(self, path: PhasePath) -> int:
        """获取畴变体总数

        Relative path: isocore/distortion/domain_generator.py"""

        ops = self.iso.get_domain_operations(
            path.parent_sg_number, path.subgroup_index
        )
        return len(ops)

    def generate_domains(self, distorted_structure: Structure,
                    path: PhasePath) -> List[Structure]:
        """
        生成所有畴变体结构

        Args:
            distorted_structure: 基准畸变结构（畴 1）
            path: 相变路径

        Returns:
            List[Structure]: 所有畴变体结构
        
        Relative path: isocore/distortion/domain_generator.py"""

        domain_ops = self.iso.get_domain_operations(
            path.parent_sg_number, path.subgroup_index
        )

        domains = []
        for op_matrix in domain_ops:
            domain = self._apply_domain_operation(distorted_structure, op_matrix)
            domains.append(domain)

        return domains

    @staticmethod
    def _apply_domain_operation(structure: Structure,
                                op_matrix: np.ndarray) -> Structure:
        """
        对结构应用畴对称操作（4x4 仿射变换）

        Args:
            structure: 原始结构
            op_matrix: (4, 4) 仿射变换矩阵 [R | t]

        Returns:
            变换后的结构
        
        Relative path: isocore/distortion/domain_generator.py"""

        # 分解：旋转部分 R (3x3)，平移部分 t (3,)
        R = op_matrix[:3, :3]
        t = op_matrix[:3, 3]

        # 分数坐标变换：x' = R x + t
        new_coords = structure.frac_coords @ R.T + t
        new_coords = wrap_to_unit_cell(new_coords)

        # 晶格也同步变换
        new_lattice_matrix = structure.lattice.matrix @ R.T
        from pymatgen.core import Lattice
        new_lattice = Lattice(new_lattice_matrix)

        return Structure(
            lattice=new_lattice,
            species=structure.species,
            coords=new_coords,
            coords_are_cartesian=False,
        )
