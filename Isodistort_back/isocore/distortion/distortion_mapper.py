"""
畸变模式映射器 - 将 iso（DISPLAY BUSH）输出的模式基矢映射到每个原子的位移向量

对应阶段四，步骤8：畸变基矢到原子坐标的映射

映射规则（基于实测 DISPLAY BUSH 输出）：
- BUSH 每行给出某 Wyckoff 位置一个“代表原子”的位移向量（可能多个，
  对应模式的不同分量）；本实现取第一个向量作为该位点的位移模式。
- 若某位点只有一个代表原子行（常见于 Gamma 点均匀模式），则将该位移
  均匀作用于该位点的全部等效原子。
- 若某位点有多个代表原子行，则每个原子按“周期最近邻”匹配到最近的
  代表原子，取其位移。

已知差异（见 README）：
- 位移向量为超胞分数坐标单位、最大分量为 1 的约定（与官网一致），
  但尚未实现官网 As/Ap 振幅的 normfactor 归一化换算。
"""

import numpy as np
from pymatgen.core import Structure

from ..backend import BushMode, DistortionMode


class DistortionMapper:
    """
    畸变模式 -> 原子位移 映射器

    负责将 iso 计算出的抽象模式基矢，转化为每个原子的实际位移向量。
    """

    def __init__(self) -> None:
        pass

    def map_modes_to_atoms(self, structure: Structure,
                           wyckoff_sites: list[dict],
                           modes: list[DistortionMode]) -> dict:
        """
        将畸变模式映射到结构中每个原子

        Args:
            structure: 母相晶体结构
            wyckoff_sites: 结构的 Wyckoff 位置分组信息（来自 SymmetryValidator）
            modes: iso 计算出的畸变模式列表（含 bush_modes）

        Returns:
            dict: 每个模式对应的原子位移向量
                {
                    "irrep_label": {
                        "mode": DistortionMode,
                        "displacements": np.ndarray (N_atoms x 3),
                    }
                }
        """
        n_atoms = len(structure)

        # Wyckoff 字母 -> 原子索引列表
        letter_to_indices: dict[str, list[int]] = {}
        for site_info in wyckoff_sites:
            letter_to_indices.setdefault(
                site_info["wyckoff_letter"], []
            ).extend(site_info["equivalent_indices"])

        result: dict[str, dict] = {}
        for mode in modes:
            displacements = np.zeros((n_atoms, 3))

            # 按位点字母分组 BUSH 行
            by_letter: dict[str, list[BushMode]] = {}
            for bush in mode.bush_modes:
                by_letter.setdefault(bush.wyckoff_letter, []).append(bush)

            for letter, bushes in by_letter.items():
                indices = letter_to_indices.get(letter, [])
                if not indices:
                    continue

                # 每个代表原子的位移（取第一个向量分量）
                reps: list[tuple] = []
                for bush in bushes:
                    if bush.displacements:
                        reps.append((
                            np.asarray(bush.point, dtype=float),
                            np.asarray(bush.displacements[0], dtype=float),
                        ))
                if not reps:
                    continue

                if len(reps) == 1:
                    # 均匀模式：同一位移作用于该位点的全部等效原子
                    disp = reps[0][1]
                    for idx in indices:
                        displacements[idx] = disp
                else:
                    # 多代表原子：按周期最近邻匹配
                    for idx in indices:
                        coord = np.asarray(structure[idx].frac_coords)
                        best_disp = reps[0][1]
                        best_dist = float("inf")
                        for rep_point, rep_disp in reps:
                            d = self._periodic_distance(coord, rep_point)
                            if d < best_dist:
                                best_dist = d
                                best_disp = rep_disp
                        displacements[idx] = best_disp

            result[mode.irrep_label] = {
                "mode": mode,
                "displacements": displacements,
            }

        return result

    @staticmethod
    def _periodic_distance(a: np.ndarray, b: np.ndarray) -> float:
        """分数坐标的最小镜像距离。"""
        delta = a - b
        delta -= np.round(delta)
        return float(np.linalg.norm(delta))

    @staticmethod
    def normalize_displacements(displacements: np.ndarray) -> np.ndarray:
        """归一化位移向量（最大位移为 1）"""
        max_norm = np.max(np.linalg.norm(displacements, axis=1))
        if max_norm > 0:
            return displacements / max_norm
        return displacements
