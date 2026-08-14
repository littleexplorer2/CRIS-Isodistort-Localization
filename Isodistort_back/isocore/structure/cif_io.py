"""
CIF 结构读写 - 基于 pymatgen 的晶体结构输入输出

对应阶段一，步骤1：读取晶体结构文件
"""
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser


def read_cif(file_path: str | Path, primitive: bool = False) -> Structure:
    """
    读取 CIF 文件，返回 pymatgen Structure 对象

    Args:
        file_path: CIF 文件路径
        primitive: 是否转为原胞

    Returns:
        Structure: pymatgen 晶体结构对象
    """
    parser = CifParser(str(file_path))
    structure = parser.parse_structures(primitive=primitive)[0]
    return structure
