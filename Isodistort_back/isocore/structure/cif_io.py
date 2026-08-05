"""
CIF 结构读写 - 基于 pymatgen 的晶体结构输入输出

对应阶段一，步骤1：读取晶体结构文件
实现方式：❌ 自研（基于 pymatgen 库封装）
"""
from pathlib import Path
from typing import Optional

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser, CifWriter


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
    structure = parser.get_structures(primitive=primitive)[0]
    return structure


def structure_from_params(lattice_params: list, species: list,
                        coords: list, coords_are_cartesian: bool = False
                        ) -> Structure:
    """
    从参数直接构建 Structure 对象

    Args:
        lattice_params: [a, b, c, alpha, beta, gamma]
        species: 元素符号列表
        coords: 坐标列表
        coords_are_cartesian: 是否笛卡尔坐标
    """
    from pymatgen.core import Lattice
    lattice = Lattice.from_parameters(*lattice_params)
    return Structure(
        lattice=lattice,
        species=species,
        coords=coords,
        coords_are_cartesian=coords_are_cartesian,
    )


def get_lattice_params(structure: Structure) -> tuple:
    """获取结构的晶格参数 (a, b, c, alpha, beta, gamma)"""
    lat = structure.lattice
    return (lat.a, lat.b, lat.c, lat.alpha, lat.beta, lat.gamma)


def get_atom_info(structure: Structure) -> tuple:
    """
    从 Structure 提取 findsym 所需的原子信息

    Returns:
        (atom_types, atom_positions): 去重后的种类列表，以及每个位点的坐标
    """
    # 按等效位点合并（取每个 Wyckoff 的代表原子）
    # 简化实现：直接取所有不等价原子
    species = []
    positions = []
    seen = set()

    for site in structure:
        spec = site.species_string
        # 简单去重：同一元素 + 近似坐标视为等价（实际应由对称性判断）
        key = (spec, tuple(round(x, 4) for x in site.frac_coords))
        if key not in seen:
            seen.add(key)
            species.append(spec)
            positions.append(site.frac_coords.tolist())

    return species, positions
