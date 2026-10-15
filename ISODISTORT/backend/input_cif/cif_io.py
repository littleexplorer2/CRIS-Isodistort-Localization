"""
CIF 结构读写 - 基于 pymatgen 的晶体结构输入输出

对应阶段一，步骤1：读取晶体结构文件
"""
from pathlib import Path
from typing import Any

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser


def _first_parsed_structure_block(parser: CifParser) -> dict[str, Any] | None:
    """Return the CIF block that produces ``parse_structures()[0]``.

    FullProf files commonly start with a metadata-only ``data_global`` block.
    Iterating the parser's sanitized blocks through the same structure builder
    used by :meth:`CifParser.parse_structures` keeps metadata lookup attached to
    the actual first structure instead of to the first block in the file.
    """
    # Pymatgen exposes no public mapping from parsed structures back to their
    # source blocks, so this mirrors its own ordered block iteration.
    for block in parser._cif.data.values():
        try:
            structure = parser._get_structure(
                block,
                primitive=False,
                symmetrized=False,
                check_occu=True,
            )
        except (KeyError, ValueError):
            # These are exactly the per-block failures skipped by
            # CifParser.parse_structures(on_error="warn").
            continue
        if structure is not None:
            return dict(block.data)
    return None


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


def read_cif_space_group_number(file_path: str | Path) -> int | None:
    """Return the space-group number declared by a CIF, if present.

    Method 4 must respect the daughter setting declared in the uploaded CIF.
    Detecting symmetry from coordinates instead can silently promote a P1
    full-atom input to a higher-symmetry group and change the primitive-cell
    normalization used for official As/Ap amplitudes.
    """
    parser = CifParser(str(file_path))
    block = _first_parsed_structure_block(parser)
    if block is None:
        return None
    # CIF data names are case-insensitive.  Normalizing here also collapses
    # the historical underscore/dot spellings handled below.
    normalized = {str(key).lower(): value for key, value in block.items()}
    for key in (
        "_space_group_it_number",
        "_space_group.it_number",
        "_symmetry_int_tables_number",
    ):
        value = normalized.get(key)
        if isinstance(value, list):
            value = value[0] if value else None
        if value in (None, "", "?", "."):
            continue
        try:
            number = int(float(str(value).strip()))
        except ValueError:
            continue
        if 1 <= number <= 230:
            return number
    return None


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
