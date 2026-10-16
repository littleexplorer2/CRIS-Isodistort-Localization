"""CDML / Kovalev k 点显示约定与官网↔iso 参数化换算。

这是**通用空间群 k 点元数据**（Miller-Love ↔ Kovalev、参数名与比例），
不是按 irrep/OPD/算例硬编码的「正确答案」表。

官网 ISODISTORT 的 Method 2 k 点下拉显示 “GM, k14 (0,0,0)” 等，其中
- Kovalev 编号（k14 等）来自 CDML 约定；
- 坐标采用官网的参数化约定（如 SM 用 (a,0,0)，而 iso 二进制输出 (2a,0,0)）。

本地 iso 二进制只提供 Miller-Love 记号 + iso 自身参数化，故对已知母相空间群
在此维护显示/换算表，使网页下拉与官网一致。未覆盖的空间群回退 iso 原始输出。

用户界面输入的 a/b/g 参数值为**官网约定**；调用 iso 的 ``VALUE KVALUE``
前须换算为 iso 内部参数（例如官网 LD ``g=1/6`` 对应 iso ``2a=1/6`` 即
KVALUE ``1/12``）。

格式：{空间群号: {k 点 label: (kovalev 编号, [官网坐标分量], [官网参数名])}}
"""
from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction
from typing import TYPE_CHECKING

from backend.utils.text_parser import k_coordinate_parameters, parse_affine_coordinate

if TYPE_CHECKING:
    from backend.wrappers import KPointInfo

# I4/mmm (#139) CDML / Kovalev 与官网参数化
KPOINT_OFFICIAL: dict[int, dict[str, tuple[str, list[str], list[str]]]] = {
    139: {
        "GM": ("k14", ["0", "0", "0"], []),
        "SM": ("k6", ["a", "0", "0"], ["a"]),
        "LD": ("k10", ["0", "0", "g"], ["g"]),
        "DT": ("k7", ["a", "a", "0"], ["a"]),
        "X": ("k13", ["1/2", "1/2", "0"], []),
        "M": ("k15", ["1", "1", "1"], []),
        "N": ("k11", ["1/2", "0", "1/2"], []),
        "P": ("k12", ["1/2", "1/2", "1/2"], []),
        "Q": ("k5", ["1/2", "b", "1/2"], ["b"]),
        "W": ("k9", ["1/2", "1/2", "g"], ["g"]),
        "Y": ("k8", ["a", "-a+1", "0"], ["a"]),
        "A": ("k3", ["a", "a", "g"], ["a", "g"]),
        "B": ("k1", ["a", "0", "g"], ["a", "g"]),
        "C": ("k2", ["a", "b", "0"], ["a", "b"]),
        "GP": ("k0", ["a", "b", "g"], ["a", "b", "g"]),
    },
    # P4/mmm (#123) — NdNiO2 parent; labels from Distortion-page k-point list
    123: {
        "GM": ("k17", ["0", "0", "0"], []),
        "SM": ("k10", ["a", "a", "0"], ["a"]),
        "DT": ("k8", ["0", "b", "0"], ["b"]),
        "LD": ("k13", ["0", "0", "g"], ["g"]),
        "M": ("k18", ["1/2", "1/2", "0"], []),
        "X": ("k15", ["0", "1/2", "0"], []),
        "A": ("k20", ["1/2", "1/2", "1/2"], []),
        "R": ("k16", ["0", "1/2", "1/2"], []),
        "Z": ("k19", ["0", "0", "1/2"], []),
        "Y": ("k6", ["a", "1/2", "0"], ["a"]),
        "V": ("k14", ["1/2", "1/2", "g"], ["g"]),
        "W": ("k12", ["0", "1/2", "g"], ["g"]),
        "S": ("k11", ["a", "a", "1/2"], ["a"]),
        "T": ("k7", ["a", "1/2", "1/2"], ["a"]),
        "U": ("k9", ["0", "b", "1/2"], ["b"]),
        "B": ("k3", ["0", "b", "g"], ["b", "g"]),
        "C": ("k5", ["a", "a", "g"], ["a", "g"]),
        "D": ("k1", ["a", "b", "0"], ["a", "b"]),
        "E": ("k2", ["a", "b", "1/2"], ["a", "b"]),
        "F": ("k4", ["a", "1/2", "g"], ["a", "g"]),
        "GP": ("k0", ["a", "b", "g"], ["a", "b", "g"]),
    },
}


def official_kparams_to_iso(parent_sg: int,
                            k_label: str,
                            official_values: Sequence[str],
                            iso_kp: KPointInfo) -> list[str]:
    """Solve affine coordinate equations in exact arithmetic.

    Unlisted settings use the displayed ISO parameters directly. A listed
    convention must match all three reciprocal coordinates; an unproved
    conversion raises instead of silently passing through unrelated values.
    """
    if not official_values:
        return []
    override = KPOINT_OFFICIAL.get(parent_sg, {}).get(k_label.strip())
    if override is None:
        return [str(v).strip() for v in official_values]
    _kovalev, official_coords, official_names = override
    iso_coords = list(iso_kp.coordinates or [])
    iso_names = k_coordinate_parameters(iso_coords)
    if len(official_values) != len(official_names) or len(iso_names) != len(official_names):
        raise ValueError("k-parameter count does not match the declared coordinate convention")
    if len(iso_coords) != 3 or len(set(iso_names)) != len(iso_names):
        raise ValueError("ISO k coordinates require three components and distinct parameters")

    targets = [_substitute_k_component(c, official_names, official_values)
               for c in official_coords]
    rows = []
    for expression, target in zip(iso_coords, targets, strict=True):
        parsed = parse_affine_coordinate(expression, iso_names)
        if parsed is None:
            raise ValueError("Unsupported symbolic k-coordinate conversion")
        coefficients, offset = parsed
        rows.append([*coefficients, Fraction(target) - offset])
    # Exact elimination solves A p = k_requested - b, including coupled
    # parameters. Each row operation is invertible; check all residual rows
    # and full column rank before accepting the unique coordinate solution.
    pivot_rows = []
    for column in range(len(iso_names)):
        pivot = next((i for i in range(len(pivot_rows), len(rows)) if rows[i][column]), None)
        if pivot is None:
            raise ValueError("Non-unique affine k-coordinate parameter conversion")
        row_index = len(pivot_rows)
        rows[row_index], rows[pivot] = rows[pivot], rows[row_index]
        divisor = rows[row_index][column]
        rows[row_index] = [value / divisor for value in rows[row_index]]
        for i in range(len(rows)):
            if i == row_index:
                continue
            factor = rows[i][column]
            rows[i] = [left - factor * right
                       for left, right in zip(rows[i], rows[row_index], strict=True)]
        pivot_rows.append(row_index)
    if any(not any(row[:-1]) and row[-1] for row in rows):
        raise ValueError("Inconsistent affine k-coordinate parameter conversion")
    converted = [str(rows[i][-1]) for i in pivot_rows]
    expected = [_substitute_k_component(c, official_names, official_values) for c in official_coords]
    obtained = [_substitute_k_component(c, iso_names, converted) for c in iso_coords]
    try:
        agrees = [Fraction(c) for c in expected] == [Fraction(c) for c in obtained]
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("Unresolved symbolic k-coordinate conversion") from exc
    if not agrees:
        raise ValueError("Converted k parameters do not preserve the reciprocal vector")
    return converted


def official_special_k_coords(
    parent_sg: int,
    k_label: str,
    fallback: Sequence[str] | None = None,
    k_parameters: Sequence[str] | None = None,
) -> list[str]:
    """Listed special-k coordinates (CDML / website), with parameters substituted.

    For parametric points (e.g. LD ``(0,0,g)`` with ``g=1/6``) returns
    ``["0", "0", "1/6"]`` so ``k-active`` matches the Distortion-page CIF.
    """
    override = KPOINT_OFFICIAL.get(int(parent_sg) if parent_sg else 0, {}).get(
        (k_label or "").strip()
    )
    if override is not None:
        _kovalev, coords, params = override
        if not params:
            return [str(c).strip() for c in coords]
        values = [str(v).strip() for v in (k_parameters or [])]
        return [_substitute_k_component(c, params, values) for c in coords]
    return [str(c) for c in (fallback or [])]


def _substitute_k_component(
    component: str,
    param_names: Sequence[str],
    param_values: Sequence[str],
) -> str:
    """Evaluate a rational affine component without discarding its offset.

    Unbound or unsupported symbolic expressions stay symbolic; callers that
    require numerical coordinates must reject them explicitly.
    """
    text = str(component).strip()
    parsed = parse_affine_coordinate(text, param_names)
    if parsed is None:
        return text
    coefficients, offset = parsed
    if len(param_values) != len(coefficients):
        return text
    return str(offset + sum((coefficient * Fraction(str(value).strip())
                             for coefficient, value in zip(coefficients, param_values, strict=True)),
                            Fraction(0)))


def format_k_point_display(
    parent_sg: int,
    k_label: str,
    k_parameters: Sequence[str] | None = None,
    k_coordinates: Sequence[str] | None = None,
) -> str:
    """``LD, k10 (0,0,g), g=1/6`` for CIF / modes headers."""
    label = (k_label or "").strip() or "?"
    entry = KPOINT_OFFICIAL.get(int(parent_sg) if parent_sg else 0, {}).get(label)
    if entry is None:
        coords = ",".join(str(c) for c in (k_coordinates or []))
        return f"{label} ({coords})" if coords else label
    kovalev, coords, params = entry
    coord_text = "(" + ",".join(str(c) for c in coords) + ")"
    parts = [f"{label}, {kovalev} {coord_text}"]
    values = [str(v).strip() for v in (k_parameters or [])]
    for name, raw in zip(params, values, strict=False):
        parts.append(f"{name}={raw}")
    return ", ".join(parts)
