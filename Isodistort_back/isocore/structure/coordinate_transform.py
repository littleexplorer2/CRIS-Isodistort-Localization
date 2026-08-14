"""
坐标变换工具 - 分数坐标周期化、超胞构建、坐标等价判断

对应阶段四、五：畸变生成与位点匹配的坐标基础
"""
import numpy as np
from pymatgen.core import Structure


def wrap_to_unit_cell(coords: np.ndarray) -> np.ndarray:
    """将分数坐标归一化到 [0, 1) 单胞范围内"""
    return coords % 1.0


def build_supercell(structure: Structure,
                    scaling_matrix: list | np.ndarray) -> Structure:
    """
    构建超胞

    Args:
        structure: 原始结构
        scaling_matrix: 扩胞矩阵，可以是 [a, b, c] 整数列表或 3x3 矩阵
            （行向量为母相格单位下的新基矢，与 pymatgen make_supercell 约定一致）

    Returns:
        Structure: 超胞结构
    """
    return structure * scaling_matrix


def coordinates_are_equal(c1: np.ndarray, c2: np.ndarray,
                          tolerance: float = 1e-4) -> bool:
    """
    判断两个分数坐标是否等价（考虑周期性边界）

    计算最小镜像距离，小于容差则认为等价。
    """
    diff = np.abs(c1 - c2)
    diff = np.minimum(diff, 1.0 - diff)
    return np.all(diff < tolerance)
