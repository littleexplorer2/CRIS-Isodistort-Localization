"""(3+d) 参数校验：nmod / d、波矢长度、对称操作闭合。"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ..utils import (
    DimensionMismatchError,
    InputError,
    NumericalSingularError,
)
from .constants import PHYSICAL_DIM, eps, max_nmod


def validate_nmod(value: object) -> int:
    """校验官网 nmod（= 超空间附加维度 d），返回 int。

    合法：0..max_nmod（含 0 = 普通三维）。非法类型 / 负数 / 过大抛 ``InputError``。
    """
    if value is None:
        raise InputError("nmod / superspace d is required; got None")
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise InputError(
            f"nmod / superspace d must be an integer 0..{max_nmod()}; got {value!r}"
        )
    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise InputError("nmod / superspace d is empty")
        try:
            if "." in text or "e" in text.lower():
                as_float = float(text)
                if abs(as_float - round(as_float)) > eps():
                    raise InputError(
                        f"nmod / superspace d must be an integer; got {value!r}"
                    )
                value = round(as_float)
            else:
                value = int(text)
        except ValueError as exc:
            raise InputError(
                f"nmod / superspace d must be an integer 0..{max_nmod()}; got {value!r}"
            ) from exc
    elif isinstance(value, float):
        if abs(value - round(value)) > eps():
            raise InputError(f"nmod / superspace d must be an integer; got {value!r}")
        value = round(value)
    nmod = int(value)
    if nmod < 0:
        raise InputError(f"nmod / superspace d must be >= 0; got {nmod}")
    cap = max_nmod()
    if nmod > cap:
        raise InputError(
            f"nmod / superspace d={nmod} exceeds max_nmod={cap} "
            f"(official independent incommensurate modulations are d=1,2,3)"
        )
    return nmod


def validate_space_group_number(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise InputError(f"space_group_number must be an integer 1..230; got {value!r}")
    number = int(value)
    if not 1 <= number <= 230:
        raise InputError(f"space_group_number must be in 1..230; got {number}")
    return number


def dim_3pd(nmod: int) -> int:
    return PHYSICAL_DIM + int(nmod)


def validate_vector_dim(vec: Sequence[float], expected: int, *, name: str = "vector") -> np.ndarray:
    arr = np.asarray(vec, dtype=float).reshape(-1)
    if arr.size != expected:
        raise DimensionMismatchError(
            f"{name} length {arr.size} does not match expected dimension {expected}"
        )
    if not np.all(np.isfinite(arr)):
        raise InputError(f"{name} contains non-finite values")
    return arr


def validate_matrix_shape(
    matrix: Sequence[Sequence[float]],
    expected: tuple[int, int],
    *,
    name: str = "matrix",
) -> np.ndarray:
    arr = np.asarray(matrix, dtype=float)
    if arr.shape != expected:
        raise DimensionMismatchError(
            f"{name} shape {arr.shape} does not match expected {expected}"
        )
    if not np.all(np.isfinite(arr)):
        raise InputError(f"{name} contains non-finite values")
    return arr


def require_invertible(matrix: np.ndarray, *, name: str = "matrix") -> np.ndarray:
    square = np.asarray(matrix, dtype=float)
    if square.ndim != 2 or square.shape[0] != square.shape[1]:
        raise DimensionMismatchError(f"{name} must be square; got shape {square.shape}")
    try:
        det = float(np.linalg.det(square))
    except np.linalg.LinAlgError as exc:
        raise NumericalSingularError(f"{name} is singular") from exc
    if abs(det) < eps():
        raise NumericalSingularError(f"{name} is singular (det={det})")
    return square
