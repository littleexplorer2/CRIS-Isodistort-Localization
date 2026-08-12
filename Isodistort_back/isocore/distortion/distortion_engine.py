"""
畸变引擎 - 幅度缩放、多模式混合、生成畸变结构

对应阶段五，步骤9：畸变幅度缩放与多模式混合
实现方式：❌ 自研
"""
import numpy as np
from pymatgen.core import Structure
from typing import Dict, List

from .distortion_mapper import DistortionMapper
from ..structure import build_supercell, wrap_to_unit_cell
from ..utils import get_config


class DistortionEngine:
    """
    畸变计算引擎

    功能：
    1. 单模式畸变：按指定幅度生成畸变结构
    2. 多模式混合：线性叠加多个不可约表示的畸变
    3. 晶格应变畸变：同步更新晶格参数
    4. 超胞构建：生成对应子群的超胞结构
    """

    def __init__(self, mapper: DistortionMapper = None):
        """Relative path: isocore/distortion/distortion_engine.py"""
        self.mapper = mapper or DistortionMapper()
        cfg = get_config()
        self.default_amplitude = cfg.defaults["default_amplitude"]

    def generate_single_mode(self, parent_structure: Structure,
                            mode_displacements: np.ndarray,
                            amplitude: float = None,
                            supercell: List[int] = None) -> Structure:
        """
        生成单模式畸变后的结构

        Args:
            parent_structure: 母相结构
            mode_displacements: 模式位移向量 (N, 3)，分数坐标
            amplitude: 畸变幅度
            supercell: 超胞大小 [a, b, c]

        Returns:
            Structure: 畸变后的结构
        
        Relative path: isocore/distortion/distortion_engine.py"""

        amplitude = amplitude or self.default_amplitude
        supercell = supercell or [1, 1, 1]

        # 1. 构建超胞
        structure = build_supercell(parent_structure, supercell)

        # 2. 扩展位移向量到超胞
        n_prim = len(parent_structure)
        n_super = len(structure)
        super_disp = np.tile(mode_displacements, (n_super // n_prim, 1))

        # 3. 应用畸变：新坐标 = 原坐标 + 幅度 × 位移
        scaled_disp = amplitude * super_disp
        new_frac_coords = structure.frac_coords + scaled_disp
        new_frac_coords = wrap_to_unit_cell(new_frac_coords)

        # 4. 创建新结构
        distorted = Structure(
            lattice=structure.lattice,
            species=structure.species,
            coords=new_frac_coords,
            coords_are_cartesian=False,
        )

        return distorted

    def generate_mixed_mode(self, parent_structure: Structure,
                            mode_contributions: Dict[str, float],
                            all_displacements: Dict[str, np.ndarray],
                            supercell: List[int] = None) -> Structure:
        """
        生成多模式混合畸变结构

        Args:
            parent_structure: 母相结构
            mode_contributions: {irrep_label: amplitude} 各模式的幅度
            all_displacements: {irrep_label: displacement_array} 所有模式的位移向量
            supercell: 超胞大小

        Returns:
            Structure: 混合畸变后的结构
        
        Relative path: isocore/distortion/distortion_engine.py"""

        supercell = supercell or [1, 1, 1]

        # 叠加所有模式的位移
        total_disp = None
        for irrep, amp in mode_contributions.items():
            if irrep not in all_displacements:
                continue
            if total_disp is None:
                total_disp = amp * all_displacements[irrep]
            else:
                total_disp += amp * all_displacements[irrep]

        if total_disp is None:
            return parent_structure.copy()

        return self.generate_single_mode(
            parent_structure, total_disp, amplitude=1.0, supercell=supercell
        )

    def apply_strain(self, structure: Structure,
                    strain_tensor: np.ndarray) -> Structure:
        """
        应用晶格应变畸变

        Args:
            structure: 原始结构
            strain_tensor: (3, 3) 应变张量

        Returns:
            Structure: 应变后的结构
        
        Relative path: isocore/distortion/distortion_engine.py"""
        
        F = np.eye(3) + strain_tensor
        new_lattice_matrix = structure.lattice.matrix @ F

        from pymatgen.core import Lattice
        new_lattice = Lattice(new_lattice_matrix)

        return Structure(
            lattice=new_lattice,
            species=structure.species,
            coords=structure.frac_coords,
            coords_are_cartesian=False,
        )
