"""Structure utilities and input-format compatibility."""
from __future__ import annotations

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from isocore.api import IsoDistort
from isocore.io import StructureExporter
from isocore.structure import (
    SymmetryValidator,
    build_supercell,
    coordinates_are_equal,
    read_cif,
    read_structure,
    wrap_to_unit_cell,
)


# --- from test_structure.py ---

def test_coordinate_utils():
    # wrap_to_unit_cell

    coords = np.array([[1.2, -0.3, 0.5]])
    wrapped = wrap_to_unit_cell(coords)
    assert abs(wrapped[0, 0] - 0.2) < 1e-6
    assert abs(wrapped[0, 1] - 0.7) < 1e-6

    # coordinates_are_equal
    assert coordinates_are_equal(np.array([0.1, 0.2, 0.3]),
                                np.array([0.1, 0.2, 0.3]))
    assert coordinates_are_equal(np.array([0.0, 0.0, 0.0]),
                                np.array([1.0, 1.0, 1.0]))
    assert not coordinates_are_equal(np.array([0.0, 0.0, 0.0]),
                                    np.array([0.5, 0.0, 0.0]))

    print("✅ 坐标工具测试通过")


def test_supercell():
    # 简单立方原胞

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["Na"], [[0, 0, 0]])
    supercell = build_supercell(struct, [2, 2, 2])
    assert len(supercell) == 8
    print("✅ 超胞构建测试通过")


def test_fractional_supercell_subcell():
    """回归：子群基矢行列式 < 1（带心母相的亚胞）不能走 pymatgen 整数路径。

    对应真实 bug：Fm-3m（4 原子惯用胞）的 I4/mmm 子群基矢
    det=0.5，pymatgen ``Structure * matrix`` 只接受整数矩阵，
    曾抛 LinAlgError: Singular matrix。
    """
    lattice = Lattice.cubic(3.6)
    struct = Structure(
        lattice,
        ["Fe", "Fe", "Fe", "Fe"],
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
    )
    basis = [[-0.5, 0.5, 0.0], [-0.5, -0.5, 0.0], [0.0, 0.0, 1.0]]
    sub = build_supercell(struct, basis)
    # 亚胞：4 原子 -> 2 原子（重合位点按周期合并）
    assert len(sub) == 2
    sg = SpacegroupAnalyzer(sub, symprec=1e-3).get_space_group_number()
    # 零畸变的 I4/mmm 亚胞仍为立方 Fm-3m（c/a=√2 的 bct 就是 fcc 晶格）
    assert sg == 225

    # 整数矩阵仍走原路径（回归保护）
    assert len(build_supercell(struct, np.eye(3, dtype=int))) == 4
    print("✅ 分数基矢亚胞构建测试通过")


def test_fractional_supercell_subcell_multi_species_merge():
    """回归：多物种亚胞的周期合并不得因索引错位漏并（原子重叠）。

    对应真实 bug：NaCl（Fm-3m，4 Na + 4 Cl）的 Imm2 子群基矢 det=0.5，
    4 个 Cl 映射到亚胞后两两重合。旧合并逻辑用 ``structure[j]`` 取
    「已收集坐标」对应的物种——当 Na 先被合并跳过（coords 索引与
    structure 索引错位）后，Cl 的物种检查会拿 Cl 与 Na 比较，重合的
    Cl 无法合并，畸变结构出现原子重叠，spglib 抛
    SymmetryUndeterminedError（NaCl/MgO 外部真实结构复现）。
    """
    lattice = Lattice.cubic(5.64)
    struct = Structure(
        lattice,
        ["Na", "Na", "Na", "Na", "Cl", "Cl", "Cl", "Cl"],
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0],
         [0.5, 0.5, 0.5], [0, 0, 0.5], [0, 0.5, 0], [0.5, 0, 0]],
    )
    # Imm2 子群基矢（det=0.5）：8 原子 -> 4 原子（2 Na + 2 Cl）
    basis = [[0.0, 0.0, -1.0], [-0.5, 0.5, 0.0], [0.5, 0.5, 0.0]]
    sub = build_supercell(struct, basis)
    assert len(sub) == 4, f"亚胞应合并为 4 原子，实际 {len(sub)}"
    # 无原子重叠：合并后恰为 2 Na + 2 Cl（str 为 "Na1"/"Cl1" 形式）
    species_list = [str(s.species) for s in sub]
    assert species_list.count("Cl1") == 2 and species_list.count("Na1") == 2
    print("✅ 多物种亚胞合并回归测试通过")


def test_symmetry_validator():
    # NaCl 结构（P 原胞：Na 与 Cl 各一个，晶系为 Pm-3m #221）

    lattice = Lattice.cubic(5.63)
    struct = Structure(
        lattice,
        ["Na", "Cl"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    validator = SymmetryValidator()
    result = validator.validate(struct)

    assert result["space_group_number"] == 221
    assert len(result["wyckoff_sites"]) == 2
    assert not result["has_disorder"]

    print("✅ 对称性校验测试通过")
    print(f"   空间群: #{result['space_group_number']} ({result['space_group_symbol']})")
    print(f"   Wyckoff 位置: {[s['wyckoff_letter'] for s in result['wyckoff_sites']]}")




# --- from test_format_compat.py ---

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

