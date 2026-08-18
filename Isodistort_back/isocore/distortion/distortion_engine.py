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

    def generate_modes(self, parent_structure: Structure,
                       supercell: SupercellSpec = None,
                       parent_displacements: np.ndarray | None = None,
                       occupancy_patterns: list[tuple[np.ndarray, float]] | None = None,
                       ) -> Structure:
        """
        组合生成畸变结构：可同时施加原子位移与占据率（occupational）调制。

        Args:
            parent_structure: 母相结构
            supercell: 超胞规格（子群基矢 3x3 矩阵或 [a,b,c]）
            parent_displacements: (N_parent, 3) 位移向量（母相分数坐标单位），
                每个超胞副本按对应母相原子施加；None 表示无位移
            occupancy_patterns: [(pattern, amplitude), ...]，pattern 为
                长度 = 超胞原子数的 ±1 数组（+1 类保持全占据，-1 类占据率
                1-amplitude）；None 表示无占据率调制

        Returns:
            Structure: 畸变后的结构（占据率调制时位点为部分占据）
        """
        if parent_displacements is not None \
                and parent_displacements.shape[0] != len(parent_structure):
            raise ValueError(
                f"位移向量长度 {parent_displacements.shape[0]} 与母相原子数 "
                f"{len(parent_structure)} 不一致"
            )
        sc = build_supercell(parent_structure, supercell) if supercell is not None \
            else parent_structure.copy()

        coords = np.asarray(sc.frac_coords, dtype=float)
        if parent_displacements is not None:
            idx = self._map_parent_indices(parent_structure, sc, supercell)
            coords = coords + np.asarray(parent_displacements, dtype=float)[idx]
            coords = wrap_to_unit_cell(coords)

        def _elem_symbol(spec) -> str:
            """位点成分 -> 元素符号（兼容 Element / Composition / str）。"""
            if isinstance(spec, str):
                return spec
            if hasattr(spec, "symbol") and not hasattr(spec, "elements"):
                return str(spec.symbol)  # Element
            els = spec.elements  # Composition
            return str(els[0].symbol) if len(els) == 1 else str(spec)

        species = [_elem_symbol(site.species) for site in sc]
        if occupancy_patterns:
            n_sc = len(sc)
            for pattern, amp in occupancy_patterns:
                if len(pattern) != n_sc:
                    raise ValueError(
                        f"占据率模式长度 {len(pattern)} 与超胞原子数 {n_sc} 不一致"
                    )
                amp_capped = min(max(float(amp), 0.0), 1.0)
                new_species = []
                for j, sym in enumerate(species):
                    if pattern[j] > 0:
                        new_species.append(sym)
                    else:
                        occ = max(1.0 - amp_capped, 0.01)
                        if occ >= 0.999:
                            new_species.append(sym)
                        else:
                            new_species.append({sym: occ})
                species = new_species

        return Structure(sc.lattice, species, coords, coords_are_cartesian=False)

    @staticmethod
    def _map_parent_indices(parent: Structure, sc: Structure,
                            supercell) -> np.ndarray:
        """把超胞原子映射回母相原子索引（几何最近邻 + 同物种，鲁棒于排序）。"""
        basis = np.asarray(supercell, dtype=float)
        n, m = len(parent), len(sc)
        pc = np.asarray(parent.frac_coords, dtype=float)
        pj = np.asarray(sc.frac_coords, dtype=float) @ basis
        idx = np.zeros(m, dtype=int)
        for j in range(m):
            best_dist, best_i = float("inf"), -1
            for i in range(n):
                if sc[j].species_string != parent[i].species_string:
                    continue
                d = pj[j] - pc[i]
                d -= np.round(d)
                dist = float(np.linalg.norm(d))
                if dist < best_dist:
                    best_dist, best_i = dist, i
            idx[j] = best_i
        return idx

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
