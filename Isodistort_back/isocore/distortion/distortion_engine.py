"""
畸变引擎 - 幅度缩放、多模式混合、生成畸变结构

对应阶段五，步骤9：畸变幅度缩放与多模式混合

生成流程（与官网 Distortion Page 语义对齐）：
1. 在母相原胞上应用模式位移：新坐标 = 原坐标 + 幅度 × 位移向量
2. 再按子群超胞基矢（3x3 矩阵或 [a,b,c]）扩胞

已知差异（见 README）：
- 位移向量直接以母相分数坐标单位应用（方向模式正确）；
  官网的 As/Ap 振幅与 normfactor 归一化换算尚未实现，
  待与官网导出的 CIF 批量比对后校准。
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from pymatgen.core import Lattice, Structure

from ..structure import build_supercell, wrap_to_unit_cell
from ..utils import get_config
from .distortion_mapper import DistortionMapper

SupercellSpec = Sequence[int] | Sequence[Sequence[float]] | None


class DistortionEngine:
    """
    畸变计算引擎

    功能：
    1. 单模式畸变：按指定幅度生成畸变结构
    2. 多模式混合：线性叠加多个不可约表示的畸变
    3. 晶格应变畸变：同步更新晶格参数
    4. 超胞构建：生成对应子群的超胞结构
    """

    def __init__(self, mapper: DistortionMapper | None = None) -> None:
        self.mapper = mapper or DistortionMapper()
        cfg = get_config()
        self.default_amplitude = cfg.defaults["default_amplitude"]

    def generate_single_mode(self, parent_structure: Structure,
                             mode_displacements: np.ndarray,
                             amplitude: float | None = None,
                             supercell: SupercellSpec = None) -> Structure:
        """
        生成单模式畸变后的结构

        Args:
            parent_structure: 母相结构
            mode_displacements: 模式位移向量 (N, 3)，母相分数坐标单位
            amplitude: 畸变幅度（默认取配置 default_amplitude）
            supercell: 超胞规格，可为 [a,b,c] 整数列表或 3x3 矩阵
                （子群基矢，母相格单位）

        Returns:
            Structure: 畸变后的结构
        """
        amplitude = self.default_amplitude if amplitude is None else amplitude
        if mode_displacements.shape[0] != len(parent_structure):
            raise ValueError(
                f"位移向量长度 {mode_displacements.shape[0]} 与母相原子数 "
                f"{len(parent_structure)} 不一致"
            )

        # 1. 在母相原胞上应用位移
        scaled_disp = amplitude * np.asarray(mode_displacements, dtype=float)
        new_frac_coords = parent_structure.frac_coords + scaled_disp
        new_frac_coords = wrap_to_unit_cell(new_frac_coords)

        distorted = Structure(
            lattice=parent_structure.lattice,
            species=parent_structure.species,
            coords=new_frac_coords,
            coords_are_cartesian=False,
        )

        # 2. 扩胞（None / [1,1,1] 时保持原胞）
        if supercell is not None:
            distorted = build_supercell(distorted, supercell)

        return distorted

    def generate_mixed_mode(self, parent_structure: Structure,
                            mode_contributions: dict[str, float],
                            all_displacements: dict[str, np.ndarray],
                            supercell: SupercellSpec = None) -> Structure:
        """
        生成多模式混合畸变结构

        Args:
            parent_structure: 母相结构
            mode_contributions: {irrep_label: amplitude} 各模式的幅度
            all_displacements: {irrep_label: displacement_array} 所有模式的位移向量
            supercell: 超胞规格

        Returns:
            Structure: 混合畸变后的结构
        """
        total_disp: np.ndarray | None = None
        for irrep, amp in mode_contributions.items():
            if irrep not in all_displacements:
                continue
            contribution = amp * np.asarray(all_displacements[irrep], dtype=float)
            total_disp = contribution if total_disp is None else total_disp + contribution

        if total_disp is None:
            raise ValueError("未提供任何有效的模式贡献")

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
        """
        F = np.eye(3) + np.asarray(strain_tensor, dtype=float)
        new_lattice_matrix = structure.lattice.matrix @ F

        new_lattice = Lattice(new_lattice_matrix)

        return Structure(
            lattice=new_lattice,
            species=structure.species,
            coords=structure.frac_coords,
            coords_are_cartesian=False,
        )
