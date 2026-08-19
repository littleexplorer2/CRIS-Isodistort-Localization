"""
坐标变换工具 - 分数坐标周期化、超胞构建、坐标等价判断

对应阶段四、五：畸变生成与位点匹配的坐标基础
"""
import numpy as np
from pymatgen.core import Lattice, Structure

_MERGE_TOL = 1e-4


def wrap_to_unit_cell(coords: np.ndarray) -> np.ndarray:
    """将分数坐标归一化到 [0, 1) 单胞范围内"""
    return coords % 1.0


def build_supercell(structure: Structure,
                    scaling_matrix: list | np.ndarray) -> Structure:
    """
    构建超胞（或按子群基矢构建新格子）。

    Args:
        structure: 原始结构
        scaling_matrix: 扩胞矩阵，可以是 [a, b, c] 整数列表或 3x3 矩阵
            （行向量为母相格单位下的新基矢，与 pymatgen make_supercell 约定一致）

    Returns:
        Structure: 超胞/新格子结构

    说明：
    - 整数 3x3 矩阵走 pymatgen 快速路径（自动枚举全部格点拷贝）；
    - 分数 3x3 矩阵（如带心母相的子群基矢，行列式 < 1 的“亚胞”）由本地
      实现：新格子 = basis @ 旧格子，旧原子映射进新格子后按周期合并重复位点。
      pymatgen 的 ``Structure * matrix`` 只接受整数矩阵，分数矩阵会崩溃
      （LinAlgError: Singular matrix）。
    """
    sm = np.asarray(scaling_matrix, dtype=float)
    if sm.shape == (3, 3):
        if np.allclose(sm, np.round(sm)):
            return structure * sm.astype(int)
        return _build_fractional_cell(structure, sm)
    # [a, b, c] 整数列表：交给 pymatgen
    return structure * sm


def _build_fractional_cell(structure: Structure,
                           basis: np.ndarray) -> Structure:
    """按分数基矢构建新格子：新格子 = basis @ 旧格子；原子映射 + 周期合并。

    用于子群基矢行列式 < 1 的场景（带心母相）：子群格子是母相格子的
    亚胞，原子映射后可能出现重合位点（如 Fm-3m 4 原子 -> I4/mmm 2 原子），
    按周期等价合并。位移畸变会自然分离重合位点（降对称），不会误并。
    """
    new_lattice = Lattice(basis @ structure.lattice.matrix)
    frac = np.asarray(structure.frac_coords, dtype=float)
    # frac_new @ new_lattice == frac_old @ old_lattice
    # => frac_new @ (basis @ old) == frac_old @ old => frac_new = frac_old @ inv(basis)
    frac_new = frac @ np.linalg.inv(basis)
    frac_new = wrap_to_unit_cell(frac_new)

    species: list = []
    coords: list[np.ndarray] = []
    for i, f in enumerate(frac_new):
        dup = False
        for j, g in enumerate(coords):
            d = np.abs(f - g)
            d = np.minimum(d, 1.0 - d)
            # 注意：species[j] 是「已收集」坐标的物种——coords 列表可能因
            # 前面跳过重复位点而与 structure 下标错位，不能直接用
            # structure[j]（多物种亚胞如 NaCl 4 Na + 4 Cl -> 2 Na + 2 Cl
            # 时，错位会把 Cl 与 Na 比较，导致重合 Cl 无法合并、畸变结构
            # 出现原子重叠，spglib 无法确定对称性）。
            if structure[i].species == species[j] and np.all(d < _MERGE_TOL):
                dup = True
                break
        if not dup:
            species.append(structure[i].species)
            coords.append(f)
    return Structure(new_lattice, species, coords, coords_are_cartesian=False)


def coordinates_are_equal(c1: np.ndarray, c2: np.ndarray,
                          tolerance: float = 1e-4) -> bool:
    """
    判断两个分数坐标是否等价（考虑周期性边界）

    计算最小镜像距离，小于容差则认为等价。
    """
    diff = np.abs(c1 - c2)
    diff = np.minimum(diff, 1.0 - diff)
    return np.all(diff < tolerance)
