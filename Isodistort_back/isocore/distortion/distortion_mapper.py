"""
畸变模式映射器 - 将 iso 输出的抽象基矢映射到每个原子的位移向量

对应阶段四，步骤8：畸变基矢到原子坐标的映射
实现方式：❌ 自研（项目核心难点）

核心逻辑：
iso 输出的畸变基矢是对应 Wyckoff 位点的抽象向量，
需要将其分配到该位点的每个等效原子上，并完成母子群坐标变换，
才能得到实际晶体中每个原子的位移方向。
"""
import numpy as np
from pymatgen.core import Structure
from typing import List, Dict

from ..backend import DistortionMode
from ..structure import SiteMapper
from ..utils import DistortionModeError


class DistortionMapper:
    """
    畸变模式 -> 原子位移 映射器

    这是 ISODISTORT 最核心的自研模块，
    负责将 iso 计算出的抽象群论模式基矢，转化为每个原子的实际位移向量。
    """

    def __init__(self, site_mapper: SiteMapper = None):
        """Relative path: isocore/distortion/distortion_mapper.py"""
        self.site_mapper = site_mapper or SiteMapper()

    def map_modes_to_atoms(self, structure: Structure,
                        wyckoff_sites: List[Dict],
                        modes: List[DistortionMode]) -> Dict:
        """
        将畸变模式映射到结构中每个原子

        Args:
            structure: 母相晶体结构
            wyckoff_sites: 结构的 Wyckoff 位点分组信息
            modes: iso 计算出的畸变模式列表

        Returns:
            dict: 包含每个模式对应的原子位移向量
                {
                    "irrep_label": {
                        "mode_index": int,
                        "displacements": np.ndarray (N_atoms x 3),
                    }
                }
        
        Relative path: isocore/distortion/distortion_mapper.py"""

        n_atoms = len(structure)
        result = {}

        for mode in modes:
            displacements = np.zeros((n_atoms, 3))

            # 按 Wyckoff 位点分配基矢
            for site_info in wyckoff_sites:
                letter = site_info["wyckoff_letter"]
                indices = site_info["equivalent_indices"]
                multiplicity = site_info["multiplicity"]

                # 获取该位点对应的模式基矢分量
                site_basis = self._extract_site_basis(mode, letter, multiplicity)
                if site_basis is None:
                    continue

                # 将基矢分配给每个等效原子
                for i, atom_idx in enumerate(indices):
                    if i < len(site_basis):
                        displacements[atom_idx] = np.array(site_basis[i])

            result[mode.irrep_label] = {
                "mode": mode,
                "displacements": displacements,
            }

        return result

    @staticmethod
    def _extract_site_basis(mode: DistortionMode, wyckoff_letter: str,
                            multiplicity: int) -> List[List[float]]:
        """
        从模式基矢中提取对应 Wyckoff 位点的分量

        简化实现：假设基矢按位点顺序排列，每个位点占 multiplicity * 3 个分量。
        实际使用时需根据 iso 输出的具体格式调整。

        TODO: 需与 iso 输出格式严格对齐，这是结果一致性的关键。
        
        Relative path: isocore/distortion/distortion_mapper.py"""

        if not mode.basis_vectors:
            return None

        # 简化版：直接返回基矢，假设每个基矢对应一个原子的 xyz 位移
        # 实际 iso 输出格式需要根据具体输出解析
        return mode.basis_vectors

    @staticmethod
    def normalize_displacements(displacements: np.ndarray) -> np.ndarray:
        """归一化位移向量（最大位移为 1）

        Relative path: isocore/distortion/distortion_mapper.py"""

        max_norm = np.max(np.linalg.norm(displacements, axis=1))
        if max_norm > 0:
            return displacements / max_norm
        return displacements
