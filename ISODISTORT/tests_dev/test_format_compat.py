"""
格式兼容测试（第一层 §1 格式兼容用例 + 第二层接口一致性）。

覆盖：
- 不同写法的 CIF（小数/分数坐标、不同空间群设置、惯用/原胞）；
- VASP POSCAR / CONTCAR 输入（read_structure + API load_structure）；
- 同一结构不同格式读入后识别空间群一致（三接口同源一致性基础）。
"""
from __future__ import annotations

import pytest
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.api import IsoDistort
from isocore.io import StructureExporter
from isocore.structure import read_cif, read_structure


def _sg(path) -> int:
    return SpacegroupAnalyzer(read_cif(str(path)), symprec=1e-3) \
        .get_space_group_number()


def test_same_structure_different_cif_writings(tmp_path):
    """同一 SrTiO₃ 结构：不同写法 CIF 识别一致（Pm-3m #221）。"""
    lattice = [[3.905, 0, 0], [0, 3.905, 0], [0, 0, 3.905]]
    base = Structure(
        lattice,
        ["Sr", "Ti", "O", "O", "O"],
        [[0, 0, 0], [0.5, 0.5, 0.5],
         [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
    )
    a = tmp_path / "srtio3_default.cif"
    b = tmp_path / "srtio3_symprec.cif"
    c = tmp_path / "srtio3_primitive.cif"
    CifWriter(base).write_file(str(a))
    CifWriter(base, symprec=1e-3).write_file(str(b))
    CifWriter(base.get_primitive_structure()).write_file(str(c))
    assert _sg(a) == 221
    assert _sg(b) == 221
    # 原胞写法：P 心原胞读回仍识别为 Pm-3m
    assert SpacegroupAnalyzer(
        read_cif(str(c)), symprec=1e-3).get_space_group_number() == 221


def test_poscar_roundtrip_read_structure(tmp_path):
    """VASP POSCAR 输入：read_structure 读取并识别空间群。"""
    lattice = [[3.6, 0, 0], [0, 3.6, 0], [0, 0, 3.6]]
    base = Structure(lattice, ["Fe", "Fe"],
                     [[0, 0, 0], [0.5, 0.5, 0.5]])
    poscar_path = tmp_path / "POSCAR"
    Poscar(base).write_file(str(poscar_path))
    s = read_structure(str(poscar_path))
    assert len(s) == 2
    sg = SpacegroupAnalyzer(s, symprec=1e-3).get_space_group_number()
    assert sg == 229  # 体心立方 Im-3m


def test_api_load_structure_poscar(tmp_path):
    """API load_structure 直接接受 POSCAR（格式兼容 → 与 CIF 识别一致）。"""
    lattice = [[3.6, 0, 0], [0, 3.6, 0], [0, 0, 3.6]]
    base = Structure(lattice, ["Fe", "Fe"],
                     [[0, 0, 0], [0.5, 0.5, 0.5]])
    poscar_path = tmp_path / "POSCAR"
    Poscar(base).write_file(str(poscar_path))
    cif_path = tmp_path / "fe.cif"
    CifWriter(base).write_file(str(cif_path))

    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(poscar_path)
    sg_poscar = iso.symmetry_info["space_group_number"]

    iso2 = IsoDistort(language="en")
    iso2.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso2.load_structure(cif_path)
    sg_cif = iso2.symmetry_info["space_group_number"]

    assert sg_poscar == sg_cif == 229
    assert len(iso.structure) == len(iso2.structure)


def test_unknown_format_raises(tmp_path):
    """未知格式（如 .txt）必须明确报错，不静默。"""
    bad = tmp_path / "input.txt"
    bad.write_text("some text", encoding="utf-8")
    with pytest.raises(ValueError):
        read_structure(str(bad))


def test_export_poscar_reread_same_sg(tmp_path):
    """导出 POSCAR -> 重读 -> 空间群一致（导出/导入闭环）。"""
    lattice = [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]]
    base = Structure(lattice, ["Fe"], [[0, 0, 0]])
    exporter = StructureExporter(tmp_path)
    out = exporter.to_poscar(base, "fe_export")
    s = read_structure(str(out))
    assert len(s) == 1
    sg = SpacegroupAnalyzer(s, symprec=1e-3).get_space_group_number()
    assert sg == 221  # 简单立方 Pm-3m


def test_primitive_vs_conventional_equivalent():
    """原胞/惯用胞写法：同一晶体的格点等价性（供 Method1 lattice 去重依据）。"""
    from pymatgen.analysis.structure_matcher import StructureMatcher  # noqa: PLC0415 - 测试内局部导入

    # 体心立方 Fe（2 原子惯用胞；原胞 1 原子）
    cubic = Structure([[3.6, 0, 0], [0, 3.6, 0], [0, 0, 3.6]],
                      ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    conv = SpacegroupAnalyzer(cubic).get_conventional_standard_structure()
    prim = cubic.get_primitive_structure()
    assert SpacegroupAnalyzer(conv, symprec=1e-3).get_space_group_number() == 229
    assert SpacegroupAnalyzer(prim, symprec=1e-3).get_space_group_number() == 229
    # 两种写法的格点应等价（原胞归约后由 StructureMatcher 判定）
    m = StructureMatcher(primitive_cell=True)
    assert m.fit(conv, prim), "惯用胞与原胞应格点等价"
