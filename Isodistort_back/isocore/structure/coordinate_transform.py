"""
坐标变换工具 - 分数/笛卡尔、基矢变换、对称操作

对应阶段四、五，步骤8、10：坐标变换基础
实现方式：❌ 自研
"""
import numpy as np
from pymatgen.core import Structure, Lattice


def frac_to_cart(coords: np.ndarray, lattice: Lattice) -> np.ndarray:
    """分数坐标 -> 笛卡尔坐标"""
    return coords @ lattice.matrix


def cart_to_frac(coords: np.ndarray, lattice: Lattice) -> np.ndarray:
    """笛卡尔坐标 -> 分数坐标"""
    return coords @ lattice.inv_matrix


def apply_affine_transform(coords: np.ndarray,
                        transform_matrix: np.ndarray) -> np.ndarray:
    """
    对坐标应用 4x4 仿射变换矩阵

    Args:
        coords: (N, 3) 坐标数组
        transform_matrix: (4, 4) 仿射变换矩阵

    Returns:
        (N, 3) 变换后的坐标
    """
    n = len(coords)
    # 转为齐次坐标
    homogeneous = np.hstack([coords, np.ones((n, 1))])
    transformed = homogeneous @ transform_matrix.T
    return transformed[:, :3]


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

    Returns:
        Structure: 超胞结构
    """
    return structure * scaling_matrix


def transform_lattice(lattice: Lattice,
                    transform_matrix: np.ndarray) -> Lattice:
    """
    对晶格基矢应用变换矩阵

    Args:
        lattice: 原始晶格
        transform_matrix: (3, 3) 变换矩阵（右乘）

    Returns:
        新的 Lattice 对象
    """
    new_matrix = lattice.matrix @ transform_matrix
    return Lattice(new_matrix)


def coordinates_are_equal(c1: np.ndarray, c2: np.ndarray,
                        tolerance: float = 1e-4) -> bool:
    """
    判断两个分数坐标是否等价（考虑周期性边界）

    计算最小镜像距离，小于容差则认为等价。
    """
    diff = np.abs(c1 - c2)
    diff = np.minimum(diff, 1.0 - diff)
    return np.all(diff < tolerance)
