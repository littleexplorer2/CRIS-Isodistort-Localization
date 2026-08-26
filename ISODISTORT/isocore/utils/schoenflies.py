"""空间群/点群 Schoenflies 符号（对齐官网下拉显示，如 “Ci-1”、“D4h-17”）。

数据来自预生成静态表（`isocore/data/spacegroup_meta.py`，源自 spglib
数据库），查询 O(1)、零运行时开销。
"""
from __future__ import annotations

from ..data.spacegroup_meta import SPACE_GROUP_META


def schoenflies_symbol(number: int) -> str:
    """空间群 Schoenflies 符号（官网格式，如 Ci-1 / D4h-17）。"""
    entry = SPACE_GROUP_META.get(int(number))
    return entry[1] if entry else ""


def hm_symbol(number: int) -> str:
    """空间群简短 Hermann–Mauguin 符号（官网下拉显示，如 C2 / Cm / P2/m）。

    与 ISODISTORT 一致，采用 International Tables 旧简写（如 #64 Cmca，而非
    2005 年后的 Cmce）。pymatgen 的 ``SpaceGroup.symbol`` 对单斜等晶系返回
    带轴设置的完整符号（如 “C121”），与官网下拉的简短符号不一致。
    """
    entry = SPACE_GROUP_META.get(int(number))
    return entry[0] if entry else str(number)


# 32 个晶类（点群）的 Schoenflies 简写（官网 pointgroupsym 下拉显示，
# 如 “triclinic: 1     C1”）。按 HM 点群符号组织。
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
    "3": "hexagonal", "-3": "hexagonal", "32": "hexagonal",
    "3m": "hexagonal", "-3m": "hexagonal",
    "6": "hexagonal", "-6": "hexagonal", "6/m": "hexagonal",
    "622": "hexagonal", "6mm": "hexagonal", "-6m2": "hexagonal",
    "6/mmm": "hexagonal",
    "23": "cubic", "m-3": "cubic", "432": "cubic", "-43m": "cubic",
    "m-3m": "cubic",
}
