"""
位点匹配 - 将 CIF 中真实原子与 iso 输出的理论子位点对应

对应阶段三，步骤7：Wyckoff 位点分裂分析（实际匹配部分）
实现方式：⚖️ 混合实现（理论分裂复用 iso，匹配逻辑自研）
"""
import numpy as np
from pymatgen.core import Structure
from typing import List, Dict

from ..utils import SiteMappingError, get_config
from .coordinate_transform import coordinates_are_equal


class SiteMapper:
    """
    原子位点匹配器

    职责：
    1. 将母相结构中的原子按 Wyckoff 位点分组
    2. 结合 iso 输出的位点分裂结果，建立母子位点对应关系
    3. 为畸变模式基矢分配到具体原子做准备
    """

    def __init__(self, tolerance: float = None):
        """Relative path: isocore/structure/site_mapping.py"""
        cfg = get_config()
        self.tolerance = tolerance or cfg.position_tolerance

    def map_sites(self, parent_structure: Structure,
                parent_wyckoff_sites: List[Dict],
                child_site_splitting: Dict[str, List[Dict]]) -> Dict:
        """
        建立母相位点 -> 子相位点的映射关系

        Args:
            parent_structure: 母相结构
            parent_wyckoff_sites: 母相 Wyckoff 位点信息（来自 symmetry_validator）
            child_site_splitting: 每个母相位点在子群下的分裂结果
                key: 母相 Wyckoff 字母, value: 子位点列表

        Returns:
            dict: 完整的位点映射，包含每个原子对应的子位点信息
        
        Relative path: isocore/structure/site_mapping.py"""

        mapping = {
            "parent_sites": parent_wyckoff_sites,
            "site_correspondence": {},  # 母位点 -> 子位点列表
            "atom_assignments": [],     # 每个原子的位点分配
        }

        for parent_site in parent_wyckoff_sites:
            letter = parent_site["wyckoff_letter"]
            child_sites = child_site_splitting.get(letter, [])

            mapping["site_correspondence"][letter] = {
                "parent_multiplicity": parent_site["multiplicity"],
                "child_sites": child_sites,
                "species": parent_site["species"],
            }

            # 为每个等效原子分配子位点
            for atom_idx in parent_site["equivalent_indices"]:
                atom = parent_structure[atom_idx]
                mapping["atom_assignments"].append({
                    "atom_index": atom_idx,
                    "species": atom.species_string,
                    "frac_coords": atom.frac_coords.tolist(),
                    "parent_wyckoff": letter,
                    "child_wyckoff_candidates": [s["wyckoff_letter"] for s in child_sites],
                })

        return mapping

    def match_atoms_to_child_sites(self, child_structure: Structure,
                                mapping: Dict) -> List[Dict]:
        """
        将子结构中的原子与理论子位点一一匹配

        Args:
            child_structure: 子相结构（畸变后）
            mapping: 位点映射信息

        Returns:
            每个子位点对应的原子索引
        
        Relative path: isocore/structure/site_mapping.py"""

        assignments = []
        used_atoms = set()

        for atom_assign in mapping["atom_assignments"]:
            parent_coords = np.array(atom_assign["frac_coords"])

            # 在子结构中寻找最匹配的原子
            best_idx = None
            best_dist = float("inf")

            for i, site in enumerate(child_structure):
                if i in used_atoms:
                    continue
                if site.species_string != atom_assign["species"]:
                    continue

                if coordinates_are_equal(parent_coords, site.frac_coords,
                                        self.tolerance):
                    dist = np.linalg.norm(parent_coords - site.frac_coords)
                    if dist < best_dist:
                        best_dist = dist
                        best_idx = i

            if best_idx is not None:
                used_atoms.add(best_idx)
                assignments.append({
                    "parent_atom_index": atom_assign["atom_index"],
                    "child_atom_index": best_idx,
                    "parent_wyckoff": atom_assign["parent_wyckoff"],
                })
            else:
                raise SiteMappingError(
                    f"无法为原子 {atom_assign['atom_index']} "
                    f"({atom_assign['species']}) 找到匹配的子位点"
                )

        return assignments
