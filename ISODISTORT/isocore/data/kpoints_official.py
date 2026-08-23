"""官网 k 点显示覆盖表（Kovalev 编号 + 官网坐标参数化）。

官网 ISODISTORT 的 Method 2 k 点下拉显示 “GM, k14 (0,0,0)” 等，其中
- Kovalev 编号（k14 等）来自官网站点数据库（CDML k 点编号）；
- 坐标采用官网的参数化约定（如 SM 用 (a,0,0)，而 iso 二进制输出 (2a,0,0)）。

本地 iso 二进制只提供 Miller-Love 记号 + iso 自身参数化，故对已知母相空间群
在此维护覆盖表，使网页下拉与官网一致。未覆盖的空间群回退 iso 原始输出
（见 README「与官网的已知差异」第 5 条）。

用户界面输入的 a/b/g 参数值为**官网约定**；调用 iso 的 ``VALUE KVALUE``
前须换算为 iso 内部参数（例如官网 LD ``g=1/6`` 对应 iso ``2a=1/6`` 即
KVALUE ``1/12``）。

数据来源：官网选项页 HTML（webpage_info/）。
格式：{空间群号: {k 点 label: (kovalev 编号, [官网坐标分量], [官网参数名])}}
"""
from __future__ import annotations

import re
from collections.abc import Sequence
from fractions import Fraction
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..backend import KPointInfo

# EuAl4 母相：I4/mmm #139（官网 search_method2.html）
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
