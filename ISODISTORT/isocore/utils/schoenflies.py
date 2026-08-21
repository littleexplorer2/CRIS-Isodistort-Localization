"""空间群/点群 Schoenflies 符号（对齐官网下拉显示，如 “Ci-1”、“D4h-17”）。

来源：spglib 内置数据库（`get_spacegroup_type` 的 `schoenflies` 字段，
形如 “Ci^1”），转换为官网格式（上标 “^” 改为 “-”）。
"""
from __future__ import annotations

import spglib

# hall_number -> number 的映射缓存（spglib 以 hall number 索引）
_CACHE: dict[int, str] = {}


def _schoenflies(number: int) -> str:
    """返回空间群 Schoenflies 符号（官网格式，如 Ci-1 / D4h-17）。"""
    if number in _CACHE:
        return _CACHE[number]
    sym = None
    for hall in range(1, 531):
        t = spglib.get_spacegroup_type(hall)
        if t is not None and t["number"] == number:
            sym = t.get("schoenflies") or ""
            break
    if sym:
        sym = sym.replace("^", "-")
    _CACHE[number] = sym or ""
    return _CACHE[number]


def schoenflies_symbol(number: int) -> str:
    """空间群 Schoenflies 符号（官网格式）。"""
    return _schoenflies(number)


def hm_symbol(number: int) -> str:
    """空间群简短 Hermann–Mauguin 符号（官网下拉显示，如 C2 / Cm / P2/m）。

    pymatgen 的 ``SpaceGroup.symbol`` 对单斜等晶系返回带轴设置的完整符号
    （如 “C121”），与官网下拉的简短符号（“C2”）不一致；spglib 的
    ``international_short`` 即官网同款简短符号。
    """
    for hall in range(1, 531):
        t = spglib.get_spacegroup_type(hall)
        if t is not None and t["number"] == number:
            return t.get("international_short") or str(number)
    return str(number)


# 32 个晶类（点群）的 Schoenflies 简写（官网 pointgroupsym 下拉显示，
# 如 “triclinic: 1     C1”）。按 (晶系, HM 点群符号) 组织。
POINT_GROUP_SCHOENFLIES: dict[str, str] = {
    "1": "C1", "-1": "Ci",
    "2": "C2", "m": "Cs", "2/m": "C2h",
    "222": "D2", "mm2": "C2v", "mmm": "D2h",
    "4": "C4", "-4": "S4", "4/m": "C4h",
    "422": "D4", "4mm": "C4v", "-42m": "D2d", "4/mmm": "D4h",
    "3": "C3", "-3": "C3i", "32": "D3", "3m": "C3v", "-3m": "D3d",
    "6": "C6", "-6": "C3h", "6/m": "C6h",
    "622": "D6", "6mm": "C6v", "-6m2": "D3h", "6/mmm": "D6h",
    "23": "T", "m-3": "Th", "432": "O", "-43m": "Td", "m-3m": "Oh",
}

# 晶类所属晶系（官网 pointgroupsym 下拉显示晶系前缀）
POINT_GROUP_SYSTEM: dict[str, str] = {
    "1": "triclinic", "-1": "triclinic",
    "2": "monoclinic", "m": "monoclinic", "2/m": "monoclinic",
    "222": "orthorhombic", "mm2": "orthorhombic", "mmm": "orthorhombic",
    "4": "tetragonal", "-4": "tetragonal", "4/m": "tetragonal",
    "422": "tetragonal", "4mm": "tetragonal", "-42m": "tetragonal",
    "4/mmm": "tetragonal",
    "3": "trigonal", "-3": "trigonal", "32": "trigonal", "3m": "trigonal",
    "-3m": "trigonal",
    "6": "hexagonal", "-6": "hexagonal", "6/m": "hexagonal",
    "622": "hexagonal", "6mm": "hexagonal", "-6m2": "hexagonal",
    "6/mmm": "hexagonal",
    "23": "cubic", "m-3": "cubic", "432": "cubic", "-43m": "cubic",
    "m-3m": "cubic",
}
