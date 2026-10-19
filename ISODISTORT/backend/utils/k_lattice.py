"""Exact reciprocal-lattice compatibility and finite translation quotient grid."""

from __future__ import annotations

import math
from collections.abc import Sequence
from fractions import Fraction
import numpy as np
from backend.utils.lattice import as_fraction, determinant, rational_matrix


def param_value_candidates(M: np.ndarray) -> list[Fraction]:
    """Return the finite translation quotient's complete scalar grid.

    This shared quotient-lattice contract is used by superspace calculations.  It derives the
    exponent (largest Smith invariant) from exact determinants/minors instead
    of scanning every divisor up to ``det(M)`` or computing a float
    determinant.  Method 3 line searches use
    :func:`backend.method3.search._commensurate_line_parameters`, which is more selective and
    solves the actual affine line rather than guessing ``1/d`` values.
    """
    exact = rational_matrix(M)
    if any(value.denominator != 1 for row in exact for value in row):
        raise ValueError("Supercell matrix must contain exact integers")
    det = abs(int(determinant(exact)))
    if det <= 1:
        return []

    entries = [[int(value) for value in row] for row in exact]
    minors: list[int] = []
    for row_a in range(3):
        for row_b in range(row_a + 1, 3):
            for col_a in range(3):
                for col_b in range(col_a + 1, 3):
                    minors.append(
                        abs(
                            entries[row_a][col_a] * entries[row_b][col_b]
                            - entries[row_a][col_b] * entries[row_b][col_a]
                        )
                    )
    gcd_minors = 0
    for value in minors:
        gcd_minors = math.gcd(gcd_minors, value)
    exponent = det // gcd_minors if gcd_minors else det
    return sorted({Fraction(n, exponent) for n in range(1, exponent)})


def k_compatible_with_supercell(
    kvec: Sequence[str | int | float | Fraction],
    translations: Sequence[Sequence[str | int | float | Fraction]],
) -> bool:
    """k 是否为超胞倒格点（M k ∈ ℤ³）。

    ``basis_vectors`` 与 :func:`build_supercell` 一致，每一行是一条用
    母相基矢表示的子格基矢。因此子格平移 ``m_i`` 与母相
    倒格坐标 ``k`` 的相位是第 ``i`` 行与 ``k`` 的点积，即
    ``M @ k``。转置只会在对角超胞上偶然给出相同结果，轴置换或
    剪切子格会错误丢弃本应折叠到 Gamma 的波矢。
    """
    if len(kvec) != 3:
        return False
    try:
        matrix = rational_matrix(translations)
        vector = tuple(as_fraction(value) for value in kvec)
    except (TypeError, ValueError, ZeroDivisionError):
        return False
    return all(sum(row[col] * vector[col] for col in range(3)).denominator == 1 for row in matrix)
