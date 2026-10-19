"""Tests for backend.input_cif; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.input_cif import SymmetryValidator
from backend.input_cif import build_supercell
from backend.input_cif import coordinates_are_equal
from backend.input_cif import map_supercell_sites
from backend.input_cif import read_cif
from backend.input_cif import read_cif_space_group_number
from backend.input_cif import read_structure
from backend.input_cif import wrap_to_unit_cell
from backend.method4 import DistortionEngine
from backend.utils.lattice import centering_primitive_matrix
from backend.utils.lattice import determinant
from backend.utils.lattice import is_sublattice
from backend.utils.lattice import rational_matrix
from backend.utils.lattice import to_float_rows
from backend.utils.lattice import translation_lattice_from_cosets
from pymatgen.core import Lattice
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from tests.support_application_session import _robustness_CIFS_DIR
import numpy as np
import pytest


# Original coverage: test_basics.py.

def test_orbit_identity_is_periodic_at_fractional_cell_boundary():
    lattice = Lattice.cubic(4.0)
    exact = Structure(lattice, ["He"], [[0.0, 0.25, 0.5]])
    wrapped_noise = Structure(lattice, ["He"], [[-1e-12, 0.25, 0.5]])
    validator = SymmetryValidator()

    exact_id = validator.validate(exact)["wyckoff_sites"][0]["orbit_id"]
    wrapped_id = validator.validate(wrapped_noise)["wyckoff_sites"][0]["orbit_id"]
    assert exact_id == wrapped_id

# Original coverage: test_pipeline_contracts.py.

def test_supercell_mapping_conserves_sites_translations_and_species():
    parent = Structure(Lattice.from_parameters(3, 4, 5, 75, 85, 65), ["Si", "O"], [[0, 0, 0], [.2, .3, .4]])
    basis = np.array([[2, 1, 0], [0, 1, 0], [0, 0, 1]])
    child = build_supercell(parent, basis)
    indices, images = map_supercell_sites(parent, child, basis)
    assert np.bincount(indices).tolist() == [2, 2]
    assert np.allclose(child.frac_coords @ basis, parent.frac_coords[indices] + images)
    assert all(site.species == parent[index].species for site, index in zip(child, indices))
    reordered = Structure.from_sites(list(reversed(child)))
    other, translations = map_supercell_sites(parent, reordered, basis)
    assert other.tolist() == indices[::-1].tolist()
    assert np.array_equal(translations, images[::-1])

@pytest.mark.parametrize("species,coords", [(["O"], [[0, 0, 0]]), (["Si"], [[.2, .2, .2]])])
def test_mapping_rejects_missing_parent_images(species, coords):
    parent = Structure(Lattice.cubic(3), ["Si"], [[0, 0, 0]])
    with pytest.raises(ValueError, match="compatible parent sites"):
        map_supercell_sites(parent, Structure(parent.lattice, species, coords), np.eye(3))

def test_mapping_rejects_two_translations_of_one_parent_site():
    parent = Structure(Lattice.cubic(1), ["Si"], [[0, 0, 0]])
    child = Structure(parent.lattice, ["Si"], [[.5, 0, 0]])
    with pytest.raises(ValueError, match="2 compatible parent sites"):
        map_supercell_sites(parent, child, np.eye(3), tolerance=.51)

# Original coverage: test_robustness.py.

def test_invalid_basis_matrix_raises():
    """非法基矢（非 3x3）必须报错。"""
    struct = Structure(Lattice.cubic(5.0), ["A"], [[0, 0, 0]])
    engine = DistortionEngine()
    with pytest.raises(ValueError):
        engine.generate_single_mode(struct, np.zeros((1, 3)),
                                    amplitude=0.1, supercell=[[1, 0], [0, 1]])

def test_corrupt_cif_raises_clear_error(tmp_path):
    """损坏的 CIF 必须给出清晰错误而非崩溃。"""
    bad = tmp_path / "bad.cif"
    bad.write_text("this is not a cif file at all\n", encoding="utf-8")
    with pytest.raises(ValueError):
        read_cif(str(bad))

def test_occupancy_is_not_silently_clamped():
    """A missing scalar mode model cannot be repaired by clipping amplitudes."""
    parent = Structure(Lattice.cubic(4), [{"Al": 0.6}], [[0, 0, 0]])
    with pytest.raises(NotImplementedError, match="occupational"):
        DistortionEngine().generate_modes(
            parent, occupancy_patterns=[(np.array([-1.0]), 5.0)],
        )
    assert parent[0].species["Al"] == 0.6

# --- from test_robustness_extended.py ---

def test_read_cif_missing_file_raises(tmp_path):
    with pytest.raises(OSError):
        read_cif(str(tmp_path / "nope.cif"))

def test_read_cif_empty_file_raises(tmp_path):
    bad = tmp_path / "empty.cif"
    bad.write_text("", encoding="utf-8")
    with pytest.raises((ValueError, KeyError, IndexError)):
        read_cif(str(bad))

def test_read_cif_garbage_raises(tmp_path):
    bad = tmp_path / "garbage.cif"
    bad.write_text("not a cif at all\n" * 10, encoding="utf-8")
    with pytest.raises((ValueError, KeyError, IndexError)):
        read_cif(str(bad))

def test_huge_supercell_atom_count():
    """超大超胞 (10,10,10)：原子数精确 = 1000×。"""
    struct = Structure(Lattice.cubic(3.0), ["Fe"], [[0, 0, 0]])
    d = DistortionEngine().generate_single_mode(
        struct, np.zeros((1, 3)), amplitude=0.0, supercell=[10, 10, 10]
    )
    assert len(d) == 1000

def test_read_cif_roundtrip_stable(tmp_path):
    """CIF 读回两次：识别空间群一致（可重复性）。"""
    cif = _robustness_CIFS_DIR / "sg123.cif"
    if not cif.exists():
        pytest.skip("测试 CIF 不存在")
    s1 = read_cif(str(cif))
    s2 = read_cif(str(cif))
    g1 = SpacegroupAnalyzer(s1, symprec=1e-3).get_space_group_number()
    g2 = SpacegroupAnalyzer(s2, symprec=1e-3).get_space_group_number()
    assert g1 == g2 == 123
    assert np.allclose(s1.frac_coords, s2.frac_coords, atol=1e-12)

# Original coverage: test_structure.py.

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

def test_unknown_format_raises(tmp_path):
    """未知格式（如 .txt）必须明确报错，不静默。"""
    bad = tmp_path / "input.txt"
    bad.write_text("some text", encoding="utf-8")
    with pytest.raises(ValueError):
        read_structure(str(bad))

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
