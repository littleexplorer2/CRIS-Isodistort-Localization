"""
坐标变换工具 - 分数坐标周期化、超胞构建、坐标等价判断

对应阶段四、五：畸变生成与位点匹配的坐标基础
"""
from fractions import Fraction
from numbers import Integral, Real

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from backend.utils.config_loader import get_config
from backend.utils.lattice import (
    RationalMatrix,
    as_fraction,
    crystallographic_fraction,
    determinant,
    is_integral,
    is_sublattice,
    rational_matrix,
    to_float_rows,
    translation_lattice_from_cosets,
)

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
    - 分数 3x3 矩阵必须属于实际母相平移格（不能只检查行列式）。
      对合法基矢（如带心母相的子群基矢）由本地实现：枚举落入新单胞
      的全部母胞平移，再按新单胞周期合并重复位点。这同时适用于亚胞和超胞。
      pymatgen 的 ``Structure * matrix`` 只接受整数矩阵，分数矩阵会崩溃
      （LinAlgError: Singular matrix）。
    """
    raw = np.asarray(scaling_matrix, dtype=object)
    if raw.shape == (3,):
        raw = np.diag(raw)
    if raw.shape != (3, 3):
        raise ValueError("Cell basis must be a three-component scaling vector or a 3x3 matrix")
    basis = rational_matrix([
        [_cell_basis_fraction(value) for value in row] for row in raw
    ])
    if is_integral(basis):
        # Do not classify with allclose and then truncate: 1.999999999 is
        # neither 1 nor 2. Only proved integer translations use this branch.
        integer_basis = [[int(value) for value in row] for row in basis]
        return structure * integer_basis
    parent_lattice = _parent_translation_lattice(structure)
    if not is_sublattice(basis, parent_lattice):
        raise ValueError("Child cell basis is not a sublattice of the actual parent translation lattice")
    expected_count = len(structure) * abs(determinant(basis))
    if expected_count.denominator != 1:
        raise ValueError("Child cell basis does not conserve an integer atom count")
    child = _build_fractional_cell(structure, np.asarray(to_float_rows(basis)))
    if len(child) != expected_count.numerator:
        raise ValueError("Fractional cell construction failed atom-count conservation")
    expected_composition = structure.composition * float(abs(determinant(basis)))
    if not child.composition.almost_equals(expected_composition, rtol=0, atol=1e-8):
        raise ValueError("Fractional cell construction failed species/occupancy conservation")
    return child


def _cell_basis_fraction(value: object) -> Fraction:
    """Recover float rationals without rounding a stated noninteger to an integer.

    Reconstructing 1/3 from its binary float is allowed when the recovered
    rational converts back to that same float. If reconstruction changes the
    float (e.g. 1-1e-13 to 1), retain the stated decimal instead. Exact text and
    Fraction inputs keep their original values.
    """
    recovered = as_fraction(value)
    if isinstance(value, Real) and not isinstance(value, Integral):
        number = float(value)
        if float(recovered) != number:
            return Fraction(str(number))
    return recovered


def _parent_translation_lattice(structure: Structure) -> RationalMatrix:
    """Detect actual translations in uploaded axes without depending on Method 1."""
    cfg = get_config()
    dataset = SpacegroupAnalyzer(
        structure,
        symprec=cfg.symmetry_cartesian_tolerance_angstrom,
        angle_tolerance=cfg.symmetry_angle_tolerance_degrees,
    ).get_symmetry_dataset()
    if dataset is None:
        raise ValueError("spglib could not determine the parent translation lattice")
    return translation_lattice_from_cosets([
        [crystallographic_fraction(value) for value in translation]
        for rotation, translation in zip(dataset.rotations, dataset.translations, strict=True)
        if np.array_equal(rotation, np.eye(3, dtype=int))
    ])


def _build_fractional_cell(structure: Structure,
                           basis: np.ndarray) -> Structure:
    """Build a cell for any nonsingular rational basis in parent coordinates.

    Mapping only the atoms in one parent cell loses copies when ``det(basis)``
    exceeds one. Enumerate every integer parent-cell translation intersecting
    the new parallelepiped, then merge equivalent child-cell positions.
    """
    new_lattice = Lattice(basis @ structure.lattice.matrix)
    inverse = np.linalg.inv(basis)
    corners = np.array([
        np.array([x, y, z], dtype=float) @ basis
        for x in (0, 1) for y in (0, 1) for z in (0, 1)
    ])
    species: list = []
    coords: list[np.ndarray] = []
    labels: list[str | None] = []
    original_indices: list[int] = []
    for i, site in enumerate(structure):
        parent_frac = np.asarray(site.frac_coords, dtype=float)
        low = np.floor(corners.min(axis=0) - parent_frac).astype(int) - 1
        high = np.ceil(corners.max(axis=0) - parent_frac).astype(int) + 1
        for nx in range(low[0], high[0] + 1):
            for ny in range(low[1], high[1] + 1):
                for nz in range(low[2], high[2] + 1):
                    translation = np.asarray([nx, ny, nz], dtype=float)
                    child_frac = (parent_frac + translation) @ inverse
                    if np.any(child_frac < -_MERGE_TOL) or np.any(child_frac > 1 + _MERGE_TOL):
                        continue
                    child_frac = wrap_to_unit_cell(child_frac)
                    child_frac[np.isclose(child_frac, 1.0, atol=_MERGE_TOL)] = 0.0
                    if any(
                        site.species == species[j]
                        and np.all(np.abs((child_frac - old) - np.round(child_frac - old)) < _MERGE_TOL)
                        for j, old in enumerate(coords)
                    ):
                        continue
                    species.append(site.species)
                    coords.append(child_frac)
                    labels.append(site.label)
                    original_indices.append(i)
    property_names = {key for site in structure for key in site.properties}
    properties = {
        key: [structure[i].properties.get(key) for i in original_indices]
        for key in property_names
    }
    return Structure(
        new_lattice, species, coords, coords_are_cartesian=False,
        site_properties=properties, labels=labels,
    )


def coordinates_are_equal(c1: np.ndarray, c2: np.ndarray,
                          tolerance: float = 1e-4) -> bool:
    """
    判断两个分数坐标是否等价（考虑周期性边界）

    计算最小镜像距离，小于容差则认为等价。
    """
    # Fractional coordinates are defined modulo every integer lattice
    # translation, not only modulo one wrap of coordinates already in
    # ``[0, 1)``.  Reducing the signed difference first keeps the distance
    # non-negative for arbitrary equivalent images such as ``x`` and
    # ``x + (2, -3, 1)``.
    diff = np.asarray(c1, dtype=float) - np.asarray(c2, dtype=float)
    diff -= np.round(diff)
    return bool(np.all(np.abs(diff) < tolerance))
