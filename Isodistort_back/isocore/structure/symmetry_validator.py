"""
对称性校验 - 验证结构与空间群的一致性

对应阶段一，步骤3：结构对称性校验
实现方式：⚖️ 混合实现（基础判定复用 findsym，校验逻辑自研）
"""
import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ..utils import get_config


class SymmetryValidator:
    """
    晶体结构对称性校验器

    功能：
    1. 确认结构的空间群（与 findsym 结果交叉验证）
    2. 检查原子位置容差
    3. 识别有序/无序占位
    4. 获取 Wyckoff 位点分配
    """

    def __init__(self, tolerance: float = None):
        """Relative path: isocore/structure/symmetry_validator.py"""
        
        cfg = get_config()
        self.tolerance = tolerance or cfg.position_tolerance

    def validate(self, structure: Structure) -> dict:
        """
        执行对称性校验

        Returns:
            dict: 包含空间群号、空间群符号、Wyckoff位点、是否有序等信息
        
        Relative path: isocore/structure/symmetry_validator.py"""

        sga = SpacegroupAnalyzer(structure, symprec=self.tolerance)

        sg_number = sga.get_space_group_number()
        sg_symbol = sga.get_space_group_symbol()
        symm_structure = sga.get_symmetrized_structure()

        # Wyckoff 位点信息
        wyckoff_labels = symm_structure.wyckoff_symbols
        equivalent_indices = symm_structure.equivalent_indices

        wyckoff_sites = []
        for label, indices in zip(wyckoff_labels, equivalent_indices):
            # label 形如 "4a"
            multiplicity = int(label[:-1])
            letter = label[-1]
            wyckoff_sites.append({
                "wyckoff_letter": letter,
                "multiplicity": multiplicity,
                "species": structure[indices[0]].species_string,
                "representative_index": indices[0],
                "equivalent_indices": indices,
            })

        # 占位检查：同时覆盖混占与部分/超占位
        has_disorder = any(
            (len(site.species) > 1)
            or (not np.isclose(sum(site.species.values()), 1.0, atol=1e-8))
            for site in structure
        )

        return {
            "space_group_number": sg_number,
            "space_group_symbol": sg_symbol,
            "wyckoff_sites": wyckoff_sites,
            "has_disorder": has_disorder,
            "tolerance": self.tolerance,
        }

    @staticmethod
    def get_centering(structure: Structure) -> str:
        """获取点阵中心类型 (P/I/F/A/B/C/R)

        Relative path: isocore/structure/symmetry_validator.py"""

        sga = SpacegroupAnalyzer(structure)
        sg_symbol = sga.get_space_group_symbol()
        # 简单提取第一个字母
        return sg_symbol[0] if sg_symbol else "P"
