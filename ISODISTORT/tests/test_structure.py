"""Structure utilities and input-format compatibility."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import backend.api.core_api as core_api_module
from backend.api import IsoDistort
from backend.utils.lattice import (
    centering_primitive_matrix,
    determinant,
    is_sublattice,
    rational_matrix,
    to_float_rows,
    translation_lattice_from_cosets,
)
from features.export import StructureExporter
from features.input_cif import (
    SymmetryValidator,
    build_supercell,
    coordinates_are_equal,
    read_cif,
    read_cif_space_group_number,
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
    assert coordinates_are_equal(np.array([0.2, -1.8, 3.2]),
                                np.array([2.2, 0.2, 0.2]))
    assert not coordinates_are_equal(np.array([0.0, 0.0, 0.0]),
                                    np.array([0.5, 0.0, 0.0]))
    # An unwrapped difference larger than one must never become a negative
    # "distance" and pass the tolerance check.
    assert not coordinates_are_equal(np.zeros(3), np.full(3, 2.2))

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


def test_fractional_basis_expands_all_parent_translations():
    """A centered-lattice basis can be fractional and enclose several atoms."""
    parent = Structure(
        Lattice.cubic(3.0), ["Na", "Na"], [[0, 0, 0], [0.5, 0.5, 0.5]]
    )
    basis = [[1, 0, 0], [0, 1, 0], [0.5, 0.5, 1.5]]
    child = build_supercell(parent, basis)
    assert len(child) == 3
    assert child.composition == parent.composition * 1.5
    assert np.isclose(child.volume, parent.volume * 1.5)


def test_fractional_basis_rejects_nonparent_translation_even_with_integral_volume():
    """det=1 does not make a half-axis of a primitive crystal a translation."""
    parent = Structure(Lattice.cubic(3.0), ["Na"], [[0, 0, 0]])
    with pytest.raises(ValueError, match="actual parent translation lattice"):
        build_supercell(parent, np.diag([0.5, 1.0, 2.0]))


@pytest.mark.parametrize("centering", ["I", "F"])
@pytest.mark.parametrize("primitive_index", [1, 2])
@pytest.mark.parametrize("orientation", [1, -1])
def test_fractional_centered_cells_conserve_volume_atoms_and_multispecies(
    centering: str, primitive_index: int, orientation: int,
):
    translations = (
        [[0, 0, 0], [0.5, 0.5, 0.5]] if centering == "I"
        else [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
    )
    offset = np.array([0.25, 0.125, 0.2])
    count = len(translations)
    parent = Structure(
        Lattice.cubic(4.0), ["Na"] * count + ["Cl"] * count,
        translations + [np.asarray(translation) + offset for translation in translations],
        site_properties={"species_origin": ["Na"] * count + ["Cl"] * count},
    )
    basis = np.asarray(to_float_rows(centering_primitive_matrix(centering)))
    basis[0] *= primitive_index * orientation

    child = build_supercell(parent, basis)

    assert len(child) == 2 * primitive_index
    assert child.composition["Na"] == primitive_index
    assert child.composition["Cl"] == primitive_index
    assert child.volume == pytest.approx(parent.volume * primitive_index / count)
    assert all(site.properties["species_origin"] == site.specie.symbol for site in child)


@pytest.mark.parametrize("near_integer", [2.0 - 1e-9, 1.0 - 1e-13, 1.0 + 1e-13])
@pytest.mark.parametrize("vector_input", [False, True])
def test_near_integer_cell_basis_does_not_silently_round_or_truncate(
    near_integer: float, vector_input: bool,
):
    parent = Structure(Lattice.cubic(3.0), ["Na"], [[0, 0, 0]])
    factors = [near_integer, 1, 1]
    basis = factors if vector_input else np.diag(factors)
    with pytest.raises(ValueError, match="actual parent translation lattice"):
        build_supercell(parent, basis)
    assert len(build_supercell(parent, [2.0, 1.0, 1.0])) == 2


@pytest.mark.parametrize(
    "basis", [[float("nan"), 1, 1], [float("inf"), 1, 1], [0, 1, 1], [[1, 0], [0, 1]]],
)
def test_cell_basis_requires_finite_nonsingular_three_dimensional_input(basis):
    parent = Structure(Lattice.cubic(3.0), ["Na"], [[0, 0, 0]])
    with pytest.raises(ValueError):
        build_supercell(parent, basis)


def test_shared_translation_lattice_recovers_nonconventional_input_cell():
    """Pure translations identify an uploaded doubled cell, not its HM letter."""
    cosets = [[0, 0, 0], [0.5, 0, 0]]
    lattice = translation_lattice_from_cosets(cosets)
    expected = rational_matrix([[0.5, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert abs(determinant(lattice)) == abs(determinant(expected))
    assert is_sublattice(lattice, expected)
    assert is_sublattice(expected, lattice)


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


def test_symmetry_validator_tolerance_is_cartesian_and_cell_scaled():
    small = Structure(Lattice.cubic(4.0), ["Na"], [[0, 0, 0]])
    large = Structure(Lattice.cubic(10.0), ["Na"], [[0, 0, 0]])
    validator = SymmetryValidator(tolerance=1e-3, angle_tolerance_degrees=2.0)

    small_result = validator.validate(small)
    large_result = validator.validate(large)
    assert small_result["symmetry_tolerance_angstrom"] == 1e-3
    assert small_result["symmetry_angle_tolerance_degrees"] == 2.0
    assert small_result["fractional_norm_tolerance_upper_bound"] == pytest.approx(
        1e-3 / 4.0
    )
    assert large_result["fractional_norm_tolerance_upper_bound"] == pytest.approx(
        1e-3 / 10.0
    )


def test_symmetry_validator_rejects_nonphysical_tolerances():
    with pytest.raises(ValueError, match="positive Å"):
        SymmetryValidator(tolerance=0.0)
    with pytest.raises(ValueError, match="positive degree"):
        SymmetryValidator(angle_tolerance_degrees=float("nan"))




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


def test_cif_space_group_number_comes_from_first_structure_block(tmp_path):
    """FullProf metadata blocks must not override the loaded structure block."""
    path = tmp_path / "fullprof_multiblock.cif"
    path.write_text(
        """data_global
_space_group_IT_number 1
_publ_section_title 'metadata only'

data_structure
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 6.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_IT_number 139
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1 Fe 0 0 0
""",
        encoding="utf-8",
    )

    assert len(read_cif(path)) > 0
    assert read_cif_space_group_number(path) == 139


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


@pytest.mark.parametrize("entrypoint", ["load", "set"])
def test_failed_parent_install_rolls_back_complete_method_state(
    tmp_path: Path,
    monkeypatch,
    entrypoint: str,
) -> None:
    """A failed new parent must not mix with the prior Method 1--4 caches."""

    # Exercise the state transaction without requiring any external wrapper.
    iso = object.__new__(IsoDistort)
    iso._sym_val = SimpleNamespace(validate=lambda _structure: {})
    old_structure = Structure(Lattice.cubic(3.0), ["Fe"], [[0, 0, 0]])
    new_structure = Structure(Lattice.cubic(4.0), ["Ni"], [[0, 0, 0]])
    old_path = (tmp_path / "old.cif").resolve()
    old_symmetry = {"space_group_number": 221, "space_group_symbol": "Pm-3m"}
    old_findsym = object()
    old_candidate = object()
    old_mode = object()
    old_distorted = old_structure.copy()
    old_displacements = {"old-mode": np.array([[0.1, 0.0, 0.0]])}
    iso.structure = old_structure
    iso.structure_path = old_path
    iso.symmetry_info = old_symmetry
    iso._findsym = old_findsym
    iso.subgroups = [old_candidate]
    iso._selected_subgroup = old_candidate
    iso.distortion_modes = [old_mode]
    iso.mode_displacements = old_displacements
    iso.distorted_structure = old_distorted
    iso._mode_cache_key = ("old",)
    iso._generated_structure_cache_key = ("old-generated",)

    new_symmetry = {"space_group_number": 221, "space_group_symbol": "Pm-3m"}
    monkeypatch.setattr(iso._sym_val, "validate", lambda _structure: new_symmetry)
    monkeypatch.setattr(iso, "parent_wyckoff_display", lambda: [])

    def reject_new_parent() -> None:
        assert iso.structure is new_structure
        assert iso.symmetry_info is new_symmetry
        # Model a failed canonicalization after it has temporarily replaced
        # both primary fields, not merely a failure before the first write.
        iso.structure = Structure(Lattice.cubic(5.0), ["Co"], [[0, 0, 0]])
        iso.symmetry_info = {"space_group_number": 1, "space_group_symbol": "P1"}
        iso._findsym = object()
        raise ValueError("inconsistent parent setting")

    monkeypatch.setattr(iso, "_attach_standard_parent_orbits", reject_new_parent)
    if entrypoint == "load":
        new_path = tmp_path / "new.cif"
        new_path.write_text("fixture intercepted by read_cif", encoding="utf-8")
        monkeypatch.setattr(core_api_module, "read_cif", lambda _path: new_structure)
        with pytest.raises(ValueError, match="inconsistent parent setting"):
            iso.load_structure(new_path)
    else:
        with pytest.raises(ValueError, match="inconsistent parent setting"):
            iso.set_structure(new_structure)

    assert iso.structure is old_structure
    assert iso.structure_path == old_path
    assert iso.symmetry_info is old_symmetry
    assert iso._findsym is old_findsym
    assert iso.subgroups == [old_candidate]
    assert iso._selected_subgroup is old_candidate
    assert iso.distortion_modes == [old_mode]
    assert iso.mode_displacements is old_displacements
    assert iso.distorted_structure is old_distorted
    assert iso._mode_cache_key == ("old",)
    assert iso._generated_structure_cache_key == ("old-generated",)


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


def test_export_filename_cannot_escape_output_directory(tmp_path):
    """分数型 k 标签含路径分隔符时仍必须导出到指定目录。"""
    base = Structure([[4, 0, 0], [0, 4, 0], [0, 0, 4]], ["Fe"], [[0, 0, 0]])
    exporter = StructureExporter(tmp_path)
    out = exporter.to_cif(base, r"distorted_R4+[1/2,1/2,1/2]\\mode")

    assert out.parent == tmp_path
    assert out.exists()
    assert "/" not in out.name and "\\" not in out.name
    assert out.name == "distorted_R4+[12,12,12]mode.cif"


def test_primitive_vs_conventional_equivalent():
    """原胞/惯用胞写法：同一晶体的格点等价性（供 Method1 lattice 去重依据）。"""
    from pymatgen.analysis.structure_matcher import StructureMatcher

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
