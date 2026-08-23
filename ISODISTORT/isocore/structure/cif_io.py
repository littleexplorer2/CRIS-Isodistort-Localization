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


def read_structure(file_path: str | Path) -> Structure:
    """按扩展名读取常见晶体结构文件（CIF / VASP POSCAR / xyz）。

    格式兼容层（对应“格式兼容用例”）：科研用户常用的结构文件格式均可作为
    母相输入。pymatgen 的 ``Structure.from_file`` 自动识别扩展名，未知
    格式抛 ``ValueError``（明确报错，不静默）。

    Args:
        file_path: 结构文件路径（.cif / .vasp / POSCAR / CONTCAR / .xyz）

    Returns:
        Structure: pymatgen 晶体结构对象
    """
    path = Path(file_path)
    if path.suffix.lower() in (".cif", ".vasp", ".xyz") or path.name.upper() in (
        "POSCAR", "CONTCAR",
    ):
        return Structure.from_file(str(path))
    raise ValueError(
        f"不支持的结构文件格式: {path.name}（支持 CIF / VASP POSCAR / xyz）"
    )
