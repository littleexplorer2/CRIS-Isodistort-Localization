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

import re
from collections.abc import Sequence
from fractions import Fraction
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..backend import KPointInfo

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


def _param_coeff_in_component(component: str, param_name: str) -> Fraction | None:
    """官网坐标分量中对某参数名的线性系数（如 ``g``→1，``2a`` 对 ``a``→2）。"""
    comp = component.strip()
    if not re.search(r"[a-z]", comp):
        return None
    if comp == param_name:
        return Fraction(1)
    m = re.fullmatch(rf"(-?\d+){re.escape(param_name)}", comp)
    if m:
        return Fraction(int(m.group(1)))
    if param_name in comp:
        return Fraction(1)
    return None


def _iso_kvalue_scale(iso_component: str, iso_param_token: str) -> Fraction | None:
    """iso 坐标分量对 KVALUE 所设参数的线性系数（``2a`` + token ``2a`` → 2）。"""
    token = iso_param_token.strip()
    m = re.match(r"^(-?\d+)([a-z])", token)
    if m:
        return Fraction(abs(int(m.group(1))))
    if re.match(r"^[a-z]", token):
        if iso_component.strip() == token:
            return Fraction(1)
        letter = token[0]
        cm = re.fullmatch(rf"(-?\d+){re.escape(letter)}", iso_component.strip())
        if cm:
            return Fraction(abs(int(cm.group(1))))
        return Fraction(1)
    return None


def official_kparams_to_iso(parent_sg: int,
                            k_label: str,
                            official_values: Sequence[str],
                            iso_kp: "KPointInfo") -> list[str]:
    """把官网 UI 输入的 k 点参数值换算为 iso ``VALUE KVALUE`` 所需数值。

    典型：I4/mmm #139 的 LD 点官网 ``(0,0,g)``、``g=1/6`` → iso 坐标
    ``(0,0,2a)``、KVALUE ``1/12``（使倒格矢坐标仍为 ``1/6``）。

    未收录覆盖表或无法推断换算时原样返回（兼容直接传入 iso 值的脚本）。
    """
    if not official_values:
        return []
    override = KPOINT_OFFICIAL.get(parent_sg, {}).get(k_label.strip())
    if override is None:
        return [str(v).strip() for v in official_values]
    _kovalev, official_coords, official_param_names = override
    if len(official_values) != len(official_param_names):
        return [str(v).strip() for v in official_values]
    iso_params = list(iso_kp.parameters or [])
    iso_coords = list(iso_kp.coordinates or [])
    if not iso_params:
        return [str(v).strip() for v in official_values]

    converted: list[str] = []
    for i, (name, raw_val) in enumerate(
            zip(official_param_names, official_values, strict=False)):
        off_val = Fraction(str(raw_val).strip())
        off_idx = next(
            (j for j, c in enumerate(official_coords)
             if _param_coeff_in_component(c, name) is not None),
            None,
        )
        if off_idx is None or off_idx >= len(iso_coords):
            converted.append(str(raw_val).strip())
            continue
        off_scale = _param_coeff_in_component(official_coords[off_idx], name)
        iso_param = iso_params[i] if i < len(iso_params) else iso_params[-1]
        iso_scale = _iso_kvalue_scale(iso_coords[off_idx], iso_param)
        if off_scale is None or iso_scale is None or iso_scale == 0:
            converted.append(str(raw_val).strip())
        else:
            converted.append(str(off_val * off_scale / iso_scale))
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
    """Replace ``g`` / ``2a`` tokens using website parameter values."""
    text = str(component).strip()
    if not re.search(r"[a-z]", text):
        return text
    for name, raw in zip(param_names, param_values, strict=False):
        if not name or raw == "":
            continue
        try:
            val = Fraction(str(raw).strip())
        except (ValueError, ZeroDivisionError):
            continue
        scale = _param_coeff_in_component(text, name)
        if scale is None:
            continue
        return str(val * scale)
    return text


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

