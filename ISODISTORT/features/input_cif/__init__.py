"""features.input_cif - 中端「解析输入 CIF」部分。

只负责母相/女儿相 CIF 的读取、坐标与超胞变换、对称性校验；不包含任何 Method
搜索或导出逻辑。
"""
from .cif_io import read_cif, read_cif_space_group_number, read_structure
from .coordinate_transform import (
    build_supercell,
    coordinates_are_equal,
    wrap_to_unit_cell,
)
from .symmetry_validator import SymmetryValidator

__all__ = [
    "SymmetryValidator",
    "build_supercell",
    "coordinates_are_equal",
    "read_cif",
    "read_cif_space_group_number",
    "read_structure",
    "wrap_to_unit_cell",
]
