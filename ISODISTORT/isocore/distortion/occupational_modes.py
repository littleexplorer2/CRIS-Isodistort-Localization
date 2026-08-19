"""
Occupational 畸变模式生成器（v1，本地近似实现）

对应官网 Distortion Types 中的 “Occupational”（占据率/原子有序）畸变：
官网把某 Wyckoff 位点在子群作用下分裂出的不同轨道赋予不同的占据率
（如棋盘式 +1/-1 的交替占据），形成有序化模式。

本地 v1 算法（已在 README「已知差异」中说明为近似实现）：
1. 按子群超胞基矢构建超胞，把超胞原子映射回母相原子（几何最近邻匹配，
   不依赖 pymatgen 内部原子排序）。
2. 对选定物种的每个母相 Wyckoff 位点，取其全部超胞副本；
   按“超胞平移奇偶”分为 +1 / -1 两类（对应 k 子群最常见的棋盘式有序，
   如 2x2x1 反相畴有序），两类原子数相等。
3. 用 spglib 校验：对施加了 +1/-1 占据率调制的超胞做对称分析，
   若识别出的空间群号与目标子群一致，标记 validated=True；
   否则仍返回模式但标记 validated=False（用户界面与文档会给出提示）。
4. 幅度语义：+1 类位点保持全占据（占据率 1.0），-1 类位点占据率 = 1 - amplitude
   （amplitude ∈ [0, 1)，即部分占据/空位型有序）。

注意：
- 仅支持把某位点分裂为两类（+1/-1）的情况；同一母相位点在超胞中
  只有 1 个副本（无分裂）时不产生占据率模式；
- t 子群（无超胞，基矢为单位阵）下位点不分裂，通常不产生占据率模式；
- 官网按 (k, IR, OPD) 精确计算每个轨道占据率的完整算法
  （含多轨道字符模式与占据率 > 1 的富占据型有序）本地尚未实现，
  见 README 已知差异。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ..structure import build_supercell
from ..utils import get_config


@dataclass
class OccupationalMode:
    """占据率畸变模式（作用于子群超胞中的具体原子）。"""

    label: str                            # 如 occ-Al-d（物种-位点）
    species: str                          # 作用的元素符号
    wyckoff_letter: str                   # 母相 Wyckoff 位置字母
    basis_vectors: list[list[float]] = field(default_factory=list)  # 子群超胞基矢
    pattern: np.ndarray = None            # 每个超胞原子 +1/-1（对齐超胞原子顺序）
    parent_index: np.ndarray = None       # 每个超胞原子对应的母相原子索引
    validated: bool = False               # spglib 校验是否命中子群空间群
    note: str = ""                        # 附加说明（如“近似模式”）
    k_point_label: str = ""
    irrep_label: str = ""


class OccupationalModeGenerator:
    """按子群超胞生成 +1/-1 交替占据率模式（v1 近似）。"""

    def __init__(self, tolerance: float | None = None) -> None:
        cfg = get_config()
        self.tolerance = tolerance or cfg.position_tolerance

    # ----------------------------------------------------------------
    # 主入口
    # ----------------------------------------------------------------

    def generate(self, parent: Structure, wyckoff_sites: list[dict],
                 subgroup, species_scope: set[str]) -> list[OccupationalMode]:
        """
        生成子群路径下指定物种的占据率模式。

        Args:
            parent: 母相结构
            wyckoff_sites: SymmetryValidator 输出的 Wyckoff 位点信息
            subgroup: SubgroupInfo（须含 space_group_number / basis_vectors /
                k_point_label / irrep_label）
            species_scope: 允许产生占据率模式的元素集合（如 {"Al"}）

        Returns:
            List[OccupationalMode]
        """
        if not subgroup.basis_vectors or len(subgroup.basis_vectors) != 3:
            return []
        basis = np.asarray(subgroup.basis_vectors, dtype=float)
        det = abs(round(float(np.linalg.det(basis))))
        if det < 2:
            # 无超胞（t 子群 / 原胞型）：位点不分裂，无占据率模式（v1 限制）
            return []

        sc = build_supercell(parent, basis)
        parent_idx, _trans = self._map_supercell_to_parent(parent, sc, basis)

        # 母相 Wyckoff 字母 -> 母相原子索引
        letter_to_indices: dict[str, list[int]] = {}
        species_by_letter: dict[str, str] = {}
        for site in wyckoff_sites:
            letter_to_indices.setdefault(
                site["wyckoff_letter"], []
            ).extend(site["equivalent_indices"])
            species_by_letter[site["wyckoff_letter"]] = site["species"]

        modes: list[OccupationalMode] = []
        for letter, indices in letter_to_indices.items():
            species = species_by_letter[letter]
            if species not in species_scope:
                continue
            # 该位点全部超胞副本（按母相原子索引分组）
            sc_indices = [j for j in range(len(sc)) if parent_idx[j] in set(indices)]
            if len(sc_indices) < 2:
                continue
            pattern = self._binary_split(sc, sc_indices, basis, parent_idx)
            if not np.any(pattern):
                continue
            occ = self._occupancy_pattern(pattern, sc_indices)
            if occ is None:
                continue
            validated = self._validate(parent, basis, occ, subgroup)
            modes.append(OccupationalMode(
                label=f"occ-{species}-{letter}",
                species=species,
                wyckoff_letter=letter,
                basis_vectors=subgroup.basis_vectors,
                pattern=occ,
                parent_index=parent_idx,
                validated=validated,
                note="" if validated else "v1 近似：占据率模式未通过子群对称性校验",
                k_point_label=subgroup.k_point_label,
                irrep_label=subgroup.irrep_label,
            ))
        return modes

    # ----------------------------------------------------------------
    # 实现细节
    # ----------------------------------------------------------------

    @staticmethod
    def _map_supercell_to_parent(parent: Structure, sc: Structure,
                                 basis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """把超胞原子映射回母相原子（几何最近邻 + 同物种）。"""
        n, m = len(parent), len(sc)
        pc = np.asarray(parent.frac_coords, dtype=float)
        rc = np.asarray(sc.frac_coords, dtype=float)
        pj = rc @ basis  # 超胞分数坐标 -> 母相分数坐标
        parent_idx = np.zeros(m, dtype=int)
        trans = np.zeros((m, 3), dtype=float)
        for j in range(m):
            best_dist, best_i, best_d = float("inf"), -1, None
            for i in range(n):
                if sc[j].species_string != parent[i].species_string:
                    continue
                d = pj[j] - pc[i]
                d -= np.round(d)
                dist = float(np.linalg.norm(d))
                if dist < best_dist:
                    best_dist, best_i, best_d = dist, i, d
            parent_idx[j] = best_i
            trans[j] = best_d
        return parent_idx, trans

    def _binary_split(self, sc: Structure, sc_indices: list[int],
                      basis: np.ndarray, parent_idx: np.ndarray) -> np.ndarray:
        """把某位点的超胞副本分为 +1 / -1 两类（数量相等）。"""
        rc = np.asarray(sc.frac_coords, dtype=float)
        # 以“母相分数坐标”排序（对常见棋盘式有序，等价于按超胞平移奇偶分半）
        coords = [tuple(np.round(rc[j] @ basis, 6)) for j in sc_indices]
        order = sorted(range(len(sc_indices)), key=lambda k: coords[k])
        half = len(order) // 2
        pattern = np.zeros(len(sc), dtype=float)
        for rank, k in enumerate(order):
            pattern[sc_indices[k]] = 1.0 if rank < half else -1.0
        return pattern

    @staticmethod
    def _occupancy_pattern(pattern: np.ndarray,
                           sc_indices: list[int]) -> np.ndarray | None:
        """从 ±1 二分模式构造占据率模式；两类数量不等则放弃。"""
        n_plus = int(np.sum(pattern[sc_indices] > 0))
        n_minus = len(sc_indices) - n_plus
        if n_plus == 0 or n_minus == 0 or n_plus != n_minus:
            return None
        return pattern

    def _validate(self, parent: Structure, basis: np.ndarray,
                  occ: np.ndarray, subgroup) -> bool:
        """spglib 校验：+1/-1 交替占据调制的超胞对称群是否 == 目标子群。"""
        try:
            sc = build_supercell(parent, basis)
            new_species = []
            for j, site in enumerate(sc):
                if occ[j] > 0:
                    new_species.append(site.species_string)          # 全占据（+1 类）
                else:
                    new_species.append({site.species_string: 0.5})   # 半占据（-1 类）
            ordered = Structure(sc.lattice, new_species, sc.frac_coords,
                                coords_are_cartesian=False)
            sg = SpacegroupAnalyzer(ordered, symprec=self.tolerance)
            return sg.get_space_group_number() == subgroup.space_group_number
        except Exception:  # noqa: BLE001 - 校验失败不阻断，返回 False
            return False
