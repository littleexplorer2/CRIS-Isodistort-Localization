"""Tests for backend.utils; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.utils import affine_embeddings as affine_embeddings_module
from backend.utils import config_loader
from backend.utils import format_opd_line
from backend.utils import get_config
from backend.utils.affine_embeddings import AffineEmbedding
from backend.utils.affine_embeddings import AffineOperation
from backend.utils.affine_embeddings import ParentAffineGroup
from backend.utils.affine_embeddings import _hall_database_operations
from backend.utils.affine_embeddings import _hall_setting_transforms
from backend.utils.affine_embeddings import affine_equivalence
from backend.utils.affine_embeddings import embedding_from_identity
from backend.utils.affine_embeddings import enumerate_lifted_point_subgroups
from backend.utils.affine_embeddings import enumerate_target_affine_embeddings
from backend.utils.affine_embeddings import parent_affine_group
from backend.utils.affine_embeddings import standardized_embedding_identity
from backend.utils.config_loader import Config
from backend.utils.lattice import centering_multiplicity
from backend.utils.lattice import centering_primitive_matrix
from backend.utils.lattice import determinant
from backend.utils.lattice import identity_matrix
from backend.utils.lattice import inverse
from backend.utils.lattice import is_sublattice
from backend.utils.lattice import multiply
from backend.utils.lattice import rational_matrix
from backend.utils.lattice import real_translation_character_sign
from backend.utils.lattice import same_lattice
from backend.utils.lattice import wavevector_translation_phases
from backend.utils.opd_format import active_k_star_coordinates
from backend.utils.opd_format import format_k_active
from backend.utils.opd_format import format_number
from backend.utils.opd_format import k_star_tuples
from backend.utils.opd_format import resolve_active_k_star_coordinates
from backend.utils.parent_header import format_fixed_coord
from backend.utils.parent_header import format_wyckoff_site
from backend.utils.parent_header import format_wyckoff_sites
from backend.utils.parent_header import format_wyckoff_sites_from_cif
from backend.utils.schoenflies import hm_symbol
from backend.utils.schoenflies import schoenflies_symbol
from backend.utils.self_check import check_linearity
from backend.utils.self_check import check_mode_orthogonality
from backend.utils.self_check import check_subgroup_rule
from backend.utils.self_check import check_symmetry_conservation
from backend.utils.self_check import check_zero_amplitude
from backend.utils.text_parser import parse_invariant_direction_table
from backend.utils.text_parser import parse_macroscopic_tensor_blocks
from backend.utils.text_parser import parse_microscopic_vector_blocks
from backend.utils.text_parser import parse_symbolic_tensor_expression
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from pymatgen.core import Lattice
from pymatgen.core import Structure
from tests.support_backend_method3 import _i4mmm_parent
from tests.support_paths import TESTS_ROOT
from types import SimpleNamespace
import numpy as np
import os
import pytest
import spglib


# Original coverage: test_basics.py.

# --- from test_config.py ---

def test_config_load():
    cfg = get_config()
    # isobyu 二进制路径可解析（存在性由部署保证；WSL 下为 Windows 侧路径）
    assert cfg.iso_bin.name == "iso"
    assert cfg.findsym_bin.name == "findsym"
    assert cfg.temp_dir.exists(), "临时目录已创建"
    assert cfg.output_dir.exists(), "输出目录已创建"
    assert "ISODATA" in os.environ, "ISODATA 环境变量已设置"
    print("配置加载测试通过")
    print(f"   ISODATA = {os.environ['ISODATA']}")

def test_crystallographic_tolerances_have_explicit_units_and_roles():
    cfg = get_config()

    assert cfg.symmetry_cartesian_tolerance_angstrom == 1e-3
    assert cfg.position_tolerance == cfg.symmetry_cartesian_tolerance_angstrom
    assert cfg.symmetry_angle_tolerance_degrees == 5.0
    assert cfg.affine_exact_cartesian_tolerance_angstrom == 1e-7
    assert cfg.fractional_coordinate_tolerance == 1e-7
    assert cfg.lattice_tolerance == 1e-5
    # Exact synthetic affine data have a tighter Cartesian error budget than
    # uploaded/refined structures; matrix and fractional residuals are
    # dimensionless and are intentionally separate values.
    assert (
        cfg.affine_exact_cartesian_tolerance_angstrom
        < cfg.symmetry_cartesian_tolerance_angstrom
    )

# --- from test_parent_header.py ---

def test_format_fixed_coord_fractions():
    assert format_fixed_coord(0.0) == "0"
    assert format_fixed_coord(0.5) == "1/2"
    assert format_fixed_coord(0.25) == "1/4"

def test_format_wyckoff_site_free_z():
    line = format_wyckoff_site("Al", 2, 4, "e", [0.0, 0.0, 0.38])
    assert line == "Al2 4e (0,0,z), z= 0.38000"

def test_format_wyckoff_sites_eual4_like():
    lattice = Lattice.tetragonal(4.402, 11.163)
    structure = Structure(
        lattice,
        ["Eu", "Eu", "Al", "Al", "Al", "Al", "Al", "Al", "Al", "Al"],
        [
            [0, 0, 0], [0.5, 0.5, 0.5],
            [0, 0.5, 0.25], [0.5, 0, 0.25], [0.5, 0, 0.75], [0, 0.5, 0.75],
            [0, 0, 0.38], [0, 0, 0.62], [0.5, 0.5, 0.12], [0.5, 0.5, 0.88],
        ],
    )
    sites = [
        {"species": "Eu", "multiplicity": 2, "wyckoff_letter": "a",
         "representative_index": 0},
        {"species": "Al", "multiplicity": 4, "wyckoff_letter": "d",
         "representative_index": 2},
        {"species": "Al", "multiplicity": 4, "wyckoff_letter": "e",
         "representative_index": 6},
    ]
    lines = format_wyckoff_sites(structure, sites)
    assert lines == [
        "Eu1 2a (0,0,0)",
        "Al1 4d (0,1/2,1/4)",
        "Al2 4e (0,0,z), z= 0.38000",
    ]

def test_cif_label_matches_nonrepresentative_member_of_orbit(tmp_path):
    structure = Structure(
        Lattice.cubic(4.0),
        ["H", "H"],
        [[0.1, 0, 0], [0.9, 0, 0]],
    )
    sites = [{
        "species": "H",
        "multiplicity": 2,
        "wyckoff_letter": "a",
        "representative_index": 0,
        "equivalent_indices": [0, 1],
    }]
    cif = tmp_path / "equivalent_member.cif"
    cif.write_text(
        """data_test
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
H_second H 0.9 0 0
""",
        encoding="utf-8",
    )

    lines = format_wyckoff_sites_from_cif(cif, structure, sites, tol=0.01)

    assert lines is not None
    assert lines[0].startswith("H_second 2a ")
    assert sites[0]["display_label"] == "H_second"

def test_cif_label_matching_checks_species_and_repeated_orbit_distance(tmp_path):
    structure = Structure(
        Lattice.cubic(4.0),
        ["Na", "Na", "K"],
        [[0.1, 0, 0], [0.2, 0, 0], [0.201, 0, 0]],
    )
    sites = [
        {
            "species": "Na", "multiplicity": 1, "wyckoff_letter": "e",
            "representative_index": 0, "equivalent_indices": [0],
        },
        {
            "species": "Na", "multiplicity": 1, "wyckoff_letter": "e",
            "representative_index": 1, "equivalent_indices": [1],
        },
        {
            "species": "K", "multiplicity": 1, "wyckoff_letter": "e",
            "representative_index": 2, "equivalent_indices": [2],
        },
    ]
    cif = tmp_path / "nearby_orbits.cif"
    cif.write_text(
        """data_test
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
NaB Na 0.2009 0 0
K1 K 0.2010 0 0
NaA Na 0.1000 0 0
""",
        encoding="utf-8",
    )

    lines = format_wyckoff_sites_from_cif(cif, structure, sites, tol=0.01)

    assert lines is not None
    assert [site["display_label"] for site in sites] == ["NaA", "NaB", "K1"]
    assert [line.split()[0] for line in lines] == ["NaB", "K1", "NaA"]

def test_parent_wyckoff_from_cif_ndnio2():
    """CIF atom_site labels/order drive the web header (not memorized EuAl4)."""

    from application.session import IsoDistort

    cris = TESTS_ROOT.parent.parent
    cif = cris / "experiment_data" / "NdNiO2 own.cif"
    if not cif.is_file():
        import pytest
        pytest.skip("NdNiO2 own.cif not in experiment_data")
    iso = IsoDistort()
    iso.load_structure(cif)
    labels_after_load = {
        site.get("display_label")
        for site in iso.symmetry_info["wyckoff_sites"]
    }
    assert {"O", "ND", "NI"} <= labels_after_load
    lines = iso.parent_wyckoff_display()
    assert lines[0].startswith("O 2f")
    assert any(ln.startswith("ND 1d") for ln in lines)
    assert any(ln.startswith("NI 1a") for ln in lines)

def test_parent_wyckoff_from_cif_4310_uses_standard_orbit_representatives():
    """Body-centering-equivalent CIF rows must use FINDSYM standard points."""

    from application.session import IsoDistort

    cris = TESTS_ROOT.parent.parent
    cif = cris / "experiment_data" / "4310_tetra.cif"
    if not cif.is_file():
        import pytest
        pytest.skip("4310_tetra.cif not in experiment_data")

    iso = IsoDistort()
    iso.load_structure(cif)
    sites = iso.symmetry_info["wyckoff_sites"]
    assert len({site["orbit_id"] for site in sites}) == 8
    assert len({site["orbit_id"] for site in sites if site["wyckoff_letter"] == "e"}) == 5
    assert iso.parent_wyckoff_display() == [
        "La1 4e (0,0,z), z= 0.43204",
        "La2 4e (0,0,z), z= 0.30148",
        "Ni1 2a (0,0,0)",
        "Ni2 4e (0,0,z), z=-0.13885",
        "O1 8g (0,1/2,z), z= 0.36070",
        "O3 4c (0,1/2,0)",
        "O4 4e (0,0,z), z=-0.21680",
        "O2 4e (0,0,z), z=-0.06780",
    ]

# Original coverage: test_distortion.py.

# --- from test_self_check.py ---

class _StubIso:
    """最小桩：仅提供自检所需的 structure / generate_distortion。"""

    def __init__(self, struct: Structure, disp: np.ndarray):
        self.structure = struct
        self._disp = np.asarray(disp, dtype=float)
        self.mode_displacements = {}

    def generate_distortion(self, irrep_label: str, amplitude: float):
        coords = self.structure.frac_coords + amplitude * self._disp
        return Structure(self.structure.lattice, self.structure.species,
                         coords % 1.0, coords_are_cartesian=False)

def _self_parent() -> Structure:
    return Structure(Lattice.cubic(4.0), ["Fe", "Fe"],
                     [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])

def test_symmetry_conservation_same_structure():
    p = _self_parent()
    ok, detail = check_symmetry_conservation(p, p)
    assert ok, detail

def test_symmetry_conservation_extra_ops_breaks():
    """畸变结构出现母相没有的操作（如单斜母相被误建成更高对称）-> 不守恒。"""
    tetragonal = Structure(Lattice.tetragonal(4.0, 8.0), ["Fe", "Fe"],
                           [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    cubic = Structure(Lattice.cubic(4.0), ["Fe", "Fe"],
                      [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    ok, _ = check_symmetry_conservation(tetragonal, cubic)
    assert not ok  # 立方对称操作 ⊄ 四方对称操作

def test_subgroup_rule_same_sg():
    p = _self_parent()
    ok, detail = check_subgroup_rule(229, 229, p, p)
    assert ok, detail

def test_zero_amplitude_returns_parent():
    stub = _StubIso(_self_parent(), np.array([[0.2, 0, 0], [0, 0, 0]]))
    ok, detail = check_zero_amplitude(stub, "X")
    assert ok, detail

def test_linearity_doubling():
    stub = _StubIso(_self_parent(), np.array([[0.2, 0, 0], [0, 0, 0]]))
    ok, detail = check_linearity(stub, "X")
    assert ok, detail

def test_mode_orthogonality():
    ok, _ = check_mode_orthogonality({
        "A": {"displacements": np.array([[1, 0, 0], [0, 0, 0]])},
        "B": {"displacements": np.array([[0, 1, 0], [0, 0, 0]])},
    })
    assert ok
    bad, _ = check_mode_orthogonality({
        "A": {"displacements": np.array([[1, 0, 0], [0, 0, 0]])},
        "B": {"displacements": np.array([[2, 0, 0], [0, 0, 0]])},
    })
    assert not bad

# Original coverage: test_formats.py.

# --- from test_opd_format.py ---


def test_format_number_preserves_exact_fractions_without_rounding_to_zero():
    assert format_number(Fraction(1, 97)) == "1/97"
    assert format_number(Fraction(2_000_001, 2_000_000)) == "2000001/2000000"
    assert format_number("1/2000001") == "1/2000001"
    assert format_number(0.25) == "1/4"
    approximate = format_number(1 / 97)
    assert approximate == repr(1 / 97)
    assert approximate != "0"
    assert format_number(0.5000000001) == repr(0.5000000001)

def test_format_opd_line_gm1_matches_official_visible_text():
    line = format_opd_line(
        irrep_label="GM1+",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        space_group_number=139,
        space_group_symbol="I4/mmm",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        size=1,
        subgroup_index=1,
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    assert line == (
        "GM1+     P1   (a) 139 I4/mmm, basis={(1,0,0),(0,1,0),(0,0,1)}, "
        "origin=(0,0,0), s=1, i=1, k-active= (0,0,0)"
    )

def test_format_opd_line_pads_sg_number_and_x_point():
    line = format_opd_line(
        irrep_label="GM2+",
        opd_symbol="P1",
        opd_dir_raw="(a)",
        space_group_number=71,
        space_group_symbol="Immm",
        basis_raw="(1,0,0),(0,1,0),(0,0,1)",
        origin_raw="(0,0,0)",
        size=1,
        subgroup_index=2,
        k_coordinates=["0", "0", "0"],
        parent_sg=139,
    )
    assert line.startswith("GM2+     P1   (a)  71 Immm")

def test_x_star_and_k_active_two_arms():
    star = k_star_tuples(["1/2", "1/2", "0"], 139)
    assert star[0] == "(1/2,1/2,0)"
    assert len(star) == 2
    assert "(1/2,1/2,1)" in star
    both = format_k_active("(a;a)", ["1/2", "1/2", "0"], 139)
    one = format_k_active("(a;0)", ["1/2", "1/2", "0"], 139)
    assert both == " (1/2,1/2,0),(1/2,1/2,1)"
    assert one == " (1/2,1/2,0)"

def test_n_star_four_arms_partial_active():
    star = k_star_tuples(["1/2", "0", "1/2"], 139)
    assert star == [
        "(1/2,0,1/2)",
        "(1/2,1,1/2)",
        "(0,1/2,1/2)",
        "(1,1/2,1/2)",
    ]
    assert format_k_active("(a;a;0;0)", ["1/2", "0", "1/2"], 139) == (
        " (1/2,0,1/2),(1/2,1,1/2)"
    )
    assert format_k_active("(a;0;a;0)", ["1/2", "0", "1/2"], 139) == (
        " (1/2,0,1/2),(0,1/2,1/2)"
    )

def test_m_point_k_active_keeps_cdml_coords():
    star = k_star_tuples(["1", "1", "1"], 139)
    assert star == ["(1,1,1)"]
    assert format_k_active("(a)", ["1", "1", "1"], 139) == " (1,1,1)"

def test_y_star_keeps_exact_thirds_after_reciprocal_rotation():
    star = k_star_tuples(["1/3", "1/2", "0"], 123)

    assert star[:2] == ["(1/3,1/2,0)", "(1/2,1/3,0)"]
    assert len(star) == 4
    assert all("0.333333" not in arm for arm in star)
    assert format_k_active("(a;a)", ["1/3", "1/2", "0"], 123) == (
        " (1/3,1/2,0),(1/2,1/3,0)"
    )

def test_exact_vector_column_selects_conjugate_dt_and_self_conjugate_x_arms():
    # VECTOR tells us how many conjugate classes are active; the exact child
    # lattice, not local k-star list order, uniquely identifies which classes.
    assert active_k_star_coordinates(
        ("0", "1/3", "0"), 123, "VECTOR,0,0,A,0", 0,
        subgroup_translation_lattice=((3, 0, 0), (0, 2, 0), (0, 0, 1)),
    ) == (("1/3", "0", "0"), ("2/3", "0", "0"))
    assert active_k_star_coordinates(
        ("0", "1/2", "0"), 123, "VECTOR,A,B", 0,
        subgroup_translation_lattice=((3, 0, 0), (0, 2, 0), (0, 0, 1)),
    ) == (("0", "1/2", "0"),)
    assert active_k_star_coordinates(
        ("0", "1/2", "0"), 123, "VECTOR,A,B", 1,
        subgroup_translation_lattice=((2, 0, 0), (0, 3, 0), (0, 0, 1)),
    ) == (("1/2", "0", "0"),)
    exact = resolve_active_k_star_coordinates(
        ("0", "1/2", "0"), 123, "VECTOR,A,B", 0,
        subgroup_translation_lattice=((3, 0, 0), (0, 2, 0), (0, 0, 1)),
    )
    assert exact.resolution_kind == "exact_column_support"
    assert exact.active_class_count == exact.compatible_class_count == 1

def test_active_k_star_resolution_rejects_inconsistent_vector_dimensions():
    with pytest.raises(ValueError, match="cannot be partitioned"):
        active_k_star_coordinates(
            ("0", "1/3", "0"), 123, "VECTOR,A,0,0", 0,
            subgroup_translation_lattice=((3, 0, 0), (0, 2, 0), (0, 0, 1)),
        )

def test_active_k_star_resolution_marks_ambiguous_support_as_envelope():
    resolution = resolve_active_k_star_coordinates(
        ("0", "1/2", "0"), 123, "VECTOR,A,0", 0,
        subgroup_translation_lattice=((2, 0, 0), (0, 2, 0), (0, 0, 1)),
    )

    assert resolution.coordinates == (
        ("0", "1/2", "0"),
        ("1/2", "0", "0"),
    )
    assert resolution.resolution_kind == "compatible_envelope"
    assert resolution.active_class_count == 1
    assert resolution.compatible_class_count == 2

def test_active_k_star_resolution_rejects_too_few_compatible_classes():
    with pytest.raises(ValueError, match="fewer compatible k-star classes"):
        resolve_active_k_star_coordinates(
            ("0", "1/2", "0"), 123, "VECTOR,A,A", 0,
            subgroup_translation_lattice=((3, 0, 0), (0, 2, 0), (0, 0, 1)),
        )

def test_parent_header_reads_atom_sites_after_fullprof_global_block(tmp_path):
    """FullProf metadata loops must not trap atom-site discovery in an earlier block."""
    from backend.utils.parent_header import parse_cif_atom_site_rows

    path = tmp_path / "fullprof_multiblock.cif"
    path.write_text(
        """data_global
loop_
_publ_manuscript_incl_extra_item
_publ_manuscript_incl_extra_info
_publ_manuscript_incl_extra_defn
'_pd_proc_ls_prof_cR_factor' 'corrected background' no
data_structure
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
La1 1.0 1.0 0.43204 La
Ni1 0.5 0.5 0.5 Ni
""",
        encoding="utf-8",
    )

    assert parse_cif_atom_site_rows(path) == [
        {"label": "La1", "type_symbol": "La", "frac": (1.0, 1.0, 0.43204)},
        {"label": "Ni1", "type_symbol": "Ni", "frac": (0.5, 0.5, 0.5)},
    ]

def test_parent_header_fractional_distance_is_integer_image_invariant():
    from backend.utils.parent_header import _frac_dist

    assert _frac_dist([0.2, -1.8, 3.2], [2.2, 0.2, 0.2]) == pytest.approx(0.0)
    assert _frac_dist([0.0, 0.0, 0.0], [2.2, 2.2, 2.2]) == pytest.approx(0.2)

# Original coverage: test_method3_affine_standardization.py.

@pytest.mark.parametrize(
    ("source_hall", "preferred_hall"),
    ((408, 409), (402, 403)),
)
def test_standardization_converts_spglib_origin_choice_to_preferred_hall(
    source_hall: int,
    preferred_hall: int,
) -> None:
    """SG 129/125 choice 1 must become exact ISODISTORT choice 2."""

    source_type = spglib.get_spacegroup_type(source_hall)
    assert source_type is not None
    lattice = identity_matrix()
    operations = _hall_database_operations(source_hall)
    parent = ParentAffineGroup(
        lattice=lattice,
        operations=operations,
        space_group_number=int(source_type.number),
        hall_number=source_hall,
        symbol=str(source_type.international_short),
        standard_basis=lattice,
        standard_origin=(Fraction(0), Fraction(0), Fraction(0)),
    )
    embedding = AffineEmbedding(lattice, operations, source_hall)

    identity = standardized_embedding_identity(
        embedding,
        parent,
        ((4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 8.0)),
        int(source_type.number),
    )

    assert identity is not None
    assert "hall_number" not in identity
    reconstructed = embedding_from_identity(identity, parent)
    assert reconstructed.hall_number == preferred_hall
    assert (
        affine_equivalence(
            reconstructed,
            embedding,
            parent,
            allow_parent_conjugacy=False,
        )
        == "equal"
    )

def test_hall_setting_transform_supports_conventional_axis_relabelling() -> None:
    """The setting bridge is not restricted to origin-choice-only examples."""

    transforms = _hall_setting_transforms(4, 3)  # P 2 (c) -> P 2y (b)

    assert transforms
    assert any(matrix != identity_matrix() for matrix, _shift in transforms)

def test_explicit_hall_number_must_belong_to_requested_space_group() -> None:
    lattice = identity_matrix()
    operations = _hall_database_operations(408)
    parent = ParentAffineGroup(
        lattice=lattice,
        operations=operations,
        space_group_number=129,
        hall_number=408,
        symbol="P4/nmm",
    )

    with pytest.raises(ValueError, match="belongs to SG 125, not SG 129"):
        embedding_from_identity(
            {
                "space_group_number": 129,
                "symbol": "P4/nmm",
                "basis": lattice,
                "origin": (0, 0, 0),
                "hall_number": 403,
            },
            parent,
        )

# Original coverage: test_pipeline_contracts.py.

@pytest.mark.parametrize("letter,index", [("P", 1), ("A", 2), ("B", 2), ("C", 2), ("I", 2), ("F", 4), ("R", 3)])
def test_centering_index_matches_exact_volume_ratio(letter, index):
    # Independent crystallographic conventional-cell coset counts.
    assert centering_multiplicity(letter) == index
    assert abs(determinant(centering_primitive_matrix(letter))) * index == 1

# Original coverage: test_real_mode_fallback.py.

def test_centered_reciprocal_character_is_checked_on_primitive_translations():
    phases = wavevector_translation_phases(
        [Fraction(1, 2), 0, 0], centering_primitive_matrix("I"),
    )
    assert set(phases) == {Fraction(-1, 4), Fraction(1, 4)}
    with pytest.raises(ValueError, match="paired real/complex"):
        real_translation_character_sign([Fraction(1, 2), 0, 0], [Fraction(1, 2)] * 3)

# Original coverage: test_runtime_paths.py.

@pytest.fixture
def layout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    package_root = tmp_path / "arbitrarily-named-installation"
    config_path = package_root / "resources" / "config" / "settings.yaml"
    config_path.parent.mkdir(parents=True)
    monkeypatch.setattr(config_loader, "PROJECT_ROOT", package_root)
    monkeypatch.setattr(config_loader, "CONFIG_PATH", config_path)
    monkeypatch.delenv("ISODISTORT_RUNTIME_ROOT", raising=False)
    monkeypatch.setattr(config_loader.sys, "platform", "win32")
    user_data = tmp_path / "user-data"
    monkeypatch.setenv("LOCALAPPDATA", str(user_data))
    cfg = object.__new__(config_loader.Config)
    cfg._cfg = {
        "runtime": {"output_dir": "../../output", "temp_dir": "../../output/tmp"},
        "isobyu": {"bin_dir": "../isobyu", "iso_bin": "iso"},
    }
    return cfg, package_root, user_data

def _mark_source_tree(package_root: Path) -> None:
    (package_root / "pyproject.toml").touch()
    (package_root / "main_web.py").touch()

def test_source_checkout_keeps_project_output_default(layout) -> None:
    cfg, package_root, _ = layout
    _mark_source_tree(package_root)

    assert cfg.runtime_root == package_root
    assert cfg.output_dir == package_root / "output"
    assert cfg.temp_dir == package_root / "output" / "tmp"

def test_installed_read_only_layout_uses_user_state(layout, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg, package_root, user_data = layout
    original_mkdir = Path.mkdir

    def read_only_package_mkdir(path: Path, *args, **kwargs) -> None:
        if path == package_root or package_root in path.parents:
            raise PermissionError("the installed package is read-only")
        original_mkdir(path, *args, **kwargs)

    monkeypatch.setattr(Path, "mkdir", read_only_package_mkdir)

    assert cfg.output_dir == user_data / "ISODISTORT" / "output"
    assert cfg.temp_dir == user_data / "ISODISTORT" / "output" / "tmp"
    assert cfg.output_dir.is_dir()
    assert cfg.temp_dir.is_dir()
    assert cfg.iso_bin == package_root / "resources" / "isobyu" / "iso"
    assert not (package_root / "output").exists()

@pytest.mark.parametrize("source_checkout", [False, True])
def test_runtime_root_override_rebases_relative_paths_only(
    layout, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source_checkout: bool
) -> None:
    cfg, package_root, _ = layout
    if source_checkout:
        _mark_source_tree(package_root)
    override = tmp_path / "external-state"
    monkeypatch.setenv("ISODISTORT_RUNTIME_ROOT", str(override))

    assert cfg.output_dir == override / "output"
    assert cfg.temp_dir == override / "output" / "tmp"
    assert cfg.iso_bin == package_root / "resources" / "isobyu" / "iso"
    explicit_output = tmp_path / "explicit-output"
    cfg._cfg["runtime"]["output_dir"] = str(explicit_output)
    assert cfg.output_dir == explicit_output

def test_relative_override_is_rejected_before_directory_creation(
    layout, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg, package_root, _ = layout
    monkeypatch.setenv("ISODISTORT_RUNTIME_ROOT", "relative-state")
    with pytest.raises(ValueError, match="must be an absolute path"):
        _ = cfg.output_dir
    assert not (package_root / "output").exists()

def test_installed_relative_runtime_path_cannot_escape_root(layout, tmp_path: Path) -> None:
    cfg, _, _ = layout
    with pytest.raises(ValueError, match="use an absolute path"):
        cfg.resolve_runtime_path("../../../outside")
    assert cfg.resolve_runtime_path(str(tmp_path / "outside")) == tmp_path / "outside"

@pytest.mark.parametrize(
    ("platform", "expected_suffix"),
    [
        ("win32", "AppData/Local/ISODISTORT"),
        ("darwin", "Library/Application Support/ISODISTORT"),
        ("linux", ".local/state/isodistort"),
    ],
)
def test_user_state_platform_fallbacks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, platform: str, expected_suffix: str
) -> None:
    monkeypatch.setattr(config_loader.sys, "platform", platform)
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.delenv("LOCALAPPDATA", raising=False)
    monkeypatch.delenv("XDG_STATE_HOME", raising=False)
    assert config_loader._user_runtime_root() == tmp_path / expected_suffix

def test_xdg_state_home_is_used_only_when_absolute(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(config_loader.sys, "platform", "linux")
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path / "custom-state"))
    assert config_loader._user_runtime_root() == tmp_path / "custom-state" / "isodistort"
    monkeypatch.setenv("XDG_STATE_HOME", "relative-state")
    assert config_loader._user_runtime_root() == tmp_path / ".local" / "state" / "isodistort"

# Original coverage: test_search_methods.py.

def test_exact_affine_spglib_boundary_uses_configured_angstrom_tolerance(
    monkeypatch,
) -> None:
    captured = {}

    def _identify(rotations, translations, lattice=None, symprec=None):
        captured["symprec"] = symprec
        return SimpleNamespace(number=1)

    monkeypatch.setattr(
        affine_embeddings_module.spglib,
        "get_spacegroup_type_from_symmetry",
        _identify,
    )
    identity = identity_matrix()
    embedding = AffineEmbedding(
        lattice=identity,
        operations=(
            AffineOperation(
                identity,
                (Fraction(0), Fraction(0), Fraction(0)),
            ),
        ),
    )

    assert (
        affine_embeddings_module.identify_embedding_space_group(
            embedding, np.eye(3)
        )
        == 1
    )
    assert captured["symprec"] == (
        Config().affine_exact_cartesian_tolerance_angstrom
    )

    assert (
        affine_embeddings_module.identify_embedding_space_group(
            embedding, np.eye(3) * 1e10
        )
        == 1
    )
    expected_roundoff_guard = 256 * np.finfo(float).eps * 1e10
    assert captured["symprec"] == pytest.approx(expected_roundoff_guard)

def test_method3_exact_affine_identity_accepts_alternative_setting() -> None:
    parent = parent_affine_group(_i4mmm_parent(), 1e-5)
    first = {
        "space_group_number": 71,
        "symbol": "Immm",
        "basis": [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
        "origin": [0, 0, 0],
    }
    second = {
        "space_group_number": 71,
        "symbol": "Immm",
        "basis": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "origin": [0, 0, 0],
    }

    assert affine_equivalence(
        embedding_from_identity(first, parent),
        embedding_from_identity(second, parent),
        parent,
        allow_parent_conjugacy=False,
    ) == "equal"

def test_parent_affine_group_uses_uploaded_cell_translation_lattice() -> None:
    conventional = _i4mmm_parent()
    conventional_group = parent_affine_group(conventional, 1e-5)
    assert same_lattice(
        conventional_group.lattice,
        centering_primitive_matrix("I"),
    )

    primitive = conventional.get_primitive_structure()
    primitive_group = parent_affine_group(primitive, 1e-5)
    assert primitive_group.space_group_number == 139
    assert same_lattice(primitive_group.lattice, identity_matrix())

    nonstandard = conventional.copy()
    shear = rational_matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    nonstandard.make_supercell([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    nonstandard_group = parent_affine_group(nonstandard, 1e-5)
    expected = multiply(centering_primitive_matrix("I"), inverse(shear))
    assert nonstandard_group.space_group_number == 139
    assert same_lattice(nonstandard_group.lattice, expected)

def test_parent_conjugacy_does_not_erase_embedding_multiplicity() -> None:
    parent = parent_affine_group(_i4mmm_parent(), 1e-5)
    first = embedding_from_identity(
        {
            "space_group_number": 140,
            "symbol": "I4/mcm",
            "basis": [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]],
            "origin": ["-1/2", "-1/2", 0],
        },
        parent,
    )
    second = embedding_from_identity(
        {
            "space_group_number": 140,
            "symbol": "I4/mcm",
            "basis": [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]],
            "origin": [0, 0, "1/2"],
        },
        parent,
    )

    assert affine_equivalence(
        first, second, parent, allow_parent_conjugacy=False
    ) is None
    assert affine_equivalence(first, second, parent) == "parent_conjugate"
    embeddings = [first, second]
    assert len(embeddings) == 2

def test_method3_exact_affine_identity_finds_parent_conjugacy() -> None:
    parent = parent_affine_group(_i4mmm_parent(), 1e-5)
    first = {
        "space_group_number": 140,
        "symbol": "I4/mcm",
        "basis": [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]],
        "origin": ["-1/2", "-1/2", 0],
    }
    second = {
        "space_group_number": 140,
        "symbol": "I4/mcm",
        "basis": [[-1, 1, 0], [-1, -1, 0], [0, 0, 2]],
        "origin": [0, 0, "1/2"],
    }

    assert affine_equivalence(
        embedding_from_identity(first, parent),
        embedding_from_identity(second, parent),
        parent,
    ) == "parent_conjugate"

def test_method3_finite_quotient_finds_embeddings_without_ir_routes() -> None:
    """Stage A is group-theoretic: no single-IR route is needed to find lifts."""

    identity = identity_matrix()
    twofold = rational_matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=(
            AffineOperation(identity, zero),
            AffineOperation(twofold, zero),
        ),
        space_group_number=3,
        hall_number=3,
        symbol="P2",
    )
    child_lattice = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])

    embeddings = enumerate_lifted_point_subgroups(
        parent, child_lattice, [[identity, twofold]]
    )

    assert len(embeddings) == 2
    assert all(len(embedding.operations) == 2 for embedding in embeddings)
    assert all(isinstance(embedding, AffineEmbedding) for embedding in embeddings)

def test_method3_affine_lifts_are_solved_from_point_group_generators() -> None:
    """A cyclic point group has one lift variable, not one per rotation."""

    identity = identity_matrix()
    quarter_turn = rational_matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    half_turn = multiply(quarter_turn, quarter_turn)
    inverse_quarter = multiply(half_turn, quarter_turn)
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=tuple(
            AffineOperation(rotation, zero)
            for rotation in (identity, quarter_turn, half_turn, inverse_quarter)
        ),
        space_group_number=75,
        hall_number=349,
        symbol="P4",
    )
    child_lattice = rational_matrix([[4, 0, 0], [0, 4, 0], [0, 0, 1]])

    # There are 16 lifts of the single C4 generator.  The former per-rotation
    # Cartesian estimate was 16**3 and incorrectly rejected this exact search
    # at the same boundary before solving the factor-set constraints.
    embeddings = enumerate_lifted_point_subgroups(
        parent,
        child_lattice,
        [[identity, quarter_turn, half_turn, inverse_quarter]],
        max_lift_assignments=16,
    )

    assert len(embeddings) == 16
    assert all(len(embedding.operations) == 4 for embedding in embeddings)

def test_method3_finite_quotient_restricts_to_lattice_normalizer() -> None:
    identity = identity_matrix()
    quarter_turn = rational_matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    half_turn = rational_matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    inverse_quarter = rational_matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=tuple(
            AffineOperation(rotation, zero)
            for rotation in (identity, quarter_turn, half_turn, inverse_quarter)
        ),
        space_group_number=75,
        hall_number=349,
        symbol="P4",
    )
    # The 2a x b x c sublattice is not invariant under a 90-degree rotation.
    child_lattice = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])

    with pytest.raises(ValueError, match="outside the parent group"):
        enumerate_lifted_point_subgroups(
            parent,
            child_lattice,
            [[identity, quarter_turn, half_turn, inverse_quarter]],
        )

def test_method3_affine_quotient_limit_precedes_translation_coset_allocation(
    monkeypatch,
) -> None:
    identity = identity_matrix()
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=(AffineOperation(identity, zero),),
        space_group_number=1,
        hall_number=1,
        symbol="P1",
    )
    child_lattice = rational_matrix(
        [[1_000_000, 0, 0], [0, 1, 0], [0, 0, 1]]
    )

    def unexpected_coset_allocation(*_args, **_kwargs):
        raise AssertionError("translation cosets must not be allocated after preflight")

    monkeypatch.setattr(
        affine_embeddings_module,
        "translation_cosets",
        unexpected_coset_allocation,
    )
    with pytest.raises(
        RuntimeError,
        match=r"quotient order 1000000.*no translation cosets",
    ):
        enumerate_lifted_point_subgroups(
            parent,
            child_lattice,
            [[identity]],
            max_quotient_order=512,
        )

def test_method3_affine_lift_limit_precedes_translation_coset_allocation(
    monkeypatch,
) -> None:
    identity = identity_matrix()
    twofold = rational_matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=(
            AffineOperation(identity, zero),
            AffineOperation(twofold, zero),
        ),
        space_group_number=3,
        hall_number=3,
        symbol="P2",
    )
    child_lattice = rational_matrix(
        [[1_000_000, 0, 0], [0, 1, 0], [0, 0, 1]]
    )

    def unexpected_coset_allocation(*_args, **_kwargs):
        raise AssertionError("translation cosets must not be allocated after preflight")

    monkeypatch.setattr(
        affine_embeddings_module,
        "translation_cosets",
        unexpected_coset_allocation,
    )
    with pytest.raises(
        RuntimeError,
        match=r"1000000 generator-lift assignments.*no translation cosets",
    ):
        enumerate_lifted_point_subgroups(
            parent,
            child_lattice,
            [[identity, twofold]],
            max_quotient_order=2_000_000,
            max_lift_assignments=10,
        )

def test_method3_target_type_enumeration_identifies_synthetic_p2_lifts() -> None:
    identity = identity_matrix()
    twofold = rational_matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    zero = (Fraction(0), Fraction(0), Fraction(0))
    parent = ParentAffineGroup(
        lattice=identity,
        operations=(
            AffineOperation(identity, zero),
            AffineOperation(twofold, zero),
        ),
        space_group_number=3,
        hall_number=3,
        symbol="P2",
    )
    child_lattice = rational_matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])

    embeddings = enumerate_target_affine_embeddings(
        parent,
        child_lattice,
        3,
        [[3.0, 0, 0], [0, 4.0, 0], [0, 0, 5.0]],
    )

    assert len(embeddings) == 2

def test_sublattice_check_general_matrix():
    """官网 Conventional/Primitive lattice 过滤的子格判定（一般 3x3 矩阵）。"""
    ident = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    doubled = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    fmmm = [[1, 1, 0], [-1, 1, 0], [0, 0, 1]]
    assert is_sublattice(rational_matrix(doubled), rational_matrix(doubled)) is True
    assert is_sublattice(rational_matrix(ident), rational_matrix(ident)) is True
    # 倍格是单位格的子格；单位格不是倍格的子格
    assert is_sublattice(rational_matrix(doubled), rational_matrix(ident)) is True
    assert is_sublattice(rational_matrix(ident), rational_matrix(doubled)) is False
    # Fmmm 型非对角超胞（L_fmmm = {x+y 为偶数}）：倍格(全偶坐标)是其子格
    assert is_sublattice(rational_matrix(fmmm), rational_matrix(fmmm)) is True
    assert is_sublattice(rational_matrix(doubled), rational_matrix(fmmm)) is True
    assert is_sublattice(rational_matrix(fmmm), rational_matrix(doubled)) is False

# Original coverage: test_site_irreps.py.

def test_parse_exact_embedding_directions():
    rows = parse_invariant_direction_table(
        """
        *********Irrep (ML) Dir       Subgroup     Size
        GM1+       (a)       139 I4/mmm   1
        GM5+       (a,b)     2 P-1        1
        N1+        (a,b,c,d) 2 P-1        8
        *
        """
    )
    assert rows == [
        {
            "irrep_label": "GM1+",
            "direction_raw": "(a)",
            "subgroup_number": 139,
            "subgroup_symbol": "I4/mmm",
            "size": 1,
        },
        {
            "irrep_label": "GM5+",
            "direction_raw": "(a,b)",
            "subgroup_number": 2,
            "subgroup_symbol": "P-1",
            "size": 1,
        },
        {
            "irrep_label": "N1+",
            "direction_raw": "(a,b,c,d)",
            "subgroup_number": 2,
            "subgroup_symbol": "P-1",
            "size": 8,
        },
    ]

def test_parse_parametric_embedding_directions_retains_exact_k_context():
    rows = parse_invariant_direction_table(
        """
        *************Irrep (ML) k params k vector  Dir   Subgroup   Size
        GM1+                (0,0,0)   (a)   139 I4/mmm 1
        LD1        1/12     (0,0,1/6) (a,b) 99 P4mm    12
        LD1        1/6      (0,0,1/3) (a,b) 99 P4mm    6
        *
        """
    )

    assert rows[0]["k_coordinates"] == ("0", "0", "0")
    assert rows[0]["k_parameters"] == ()
    assert rows[1]["k_parameters"] == ("1/12",)
    assert rows[1]["k_coordinates"] == ("0", "0", "1/6")
    assert rows[2]["k_parameters"] == ("1/6",)
    assert rows[2]["k_coordinates"] == ("0", "0", "1/3")

def test_microscopic_parser_preserves_fractions_rows_and_repeated_blocks():
    rows = parse_microscopic_vector_blocks(
        """
        Irrep (ML) Wyckoff Irrep Point        Projected Vectors
        N1+        g       A1    (0,1/2,z)    (0,0,1/3), (1,0,0)
                                 (1,1/2,z+1)  (0,0,-1/3), (-1,0,0)
        N1+        g       A1    (0,1/2,z)    (0,1/2,0), (0,-1/2,0)
                                 (1,1/2,z+1)  (0,-1/2,0), (0,1/2,0)
        *
        """
    )
    assert len(rows) == 2
    assert [row["source_order"] for row in rows] == [0, 1]
    assert all(row["site_irrep"] == "A1" for row in rows)
    assert len(rows[0]["rows"]) == 2
    assert rows[0]["rows"][0]["point_raw"] == ("0", "1/2", "z")
    assert rows[0]["rows"][0]["displacements"][0] == (
        Fraction(0), Fraction(0), Fraction(1, 3)
    )
    assert rows[1]["rows"][1]["displacements"][1] == (
        Fraction(0), Fraction(1, 2), Fraction(0)
    )

def test_microscopic_parser_preserves_decimal_quantization_intervals():
    rows = parse_microscopic_vector_blocks(
        "GM1+ a A1 (0,0,0) (0.577,-1.155,1/2)"
    )

    assert rows[0]["rows"][0]["displacements"][0] == (
        Fraction(577, 1000),
        Fraction(-231, 200),
        Fraction(1, 2),
    )
    assert rows[0]["rows"][0]["displacement_half_steps"][0] == (
        Fraction(1, 2000),
        Fraction(1, 2000),
        Fraction(0),
    )

def test_macroscopic_parser_preserves_copies_components_and_decimal_evidence():
    tetragonal = """
    *******Irrep (ML) Basis Functions
    GM1+       xx+yy
               zz
    GM2+       xx-yy
    GM4+       xy
    GM5+       yz,-xz
    *
    """
    rows = parse_macroscopic_tensor_blocks(tetragonal)
    assert [(row["global_irrep"], len(row["components"])) for row in rows] == [
        ("GM1+", 1), ("GM1+", 1), ("GM2+", 1), ("GM4+", 1), ("GM5+", 2)
    ]
    assert rows[0]["components"][0]["coefficients"] == (
        Fraction(1), Fraction(1), Fraction(0),
        Fraction(0), Fraction(0), Fraction(0),
    )
    assert rows[-1]["components"][1]["coefficients"] == (
        Fraction(0), Fraction(0), Fraction(0),
        Fraction(0), Fraction(-1), Fraction(0),
    )

    cubic = parse_symbolic_tensor_expression("1.732xx-1.732yy")
    assert cubic["quality"] == "printed_decimal"
    assert cubic["coefficients"][:2] == (Decimal("1.732"), Decimal("-1.732"))
    radical = parse_symbolic_tensor_expression("sqrt(3)xx")
    assert radical["quality"] == "unresolved"
    assert radical["coefficients"] is None

# Original coverage: test_theory_contracts.py.

def test_gl3z_equivalence_has_no_relative_float_tolerance():
    basis = [[1, Fraction(1000000001, 10), 0], [0, 1, 0], [0, 0, 1]]
    assert not same_lattice(rational_matrix(basis), rational_matrix(np.eye(3, dtype=int)))
    basis[0][1] = 100000000
    assert same_lattice(rational_matrix(basis), rational_matrix(np.eye(3, dtype=int)))

# Original coverage: test_web.py.

def test_parent_header_schoenflies():
    assert schoenflies_symbol(139) == "D4h-17"
    # Method 3 / Method 1 下拉与官网一致的旧 IT 简写
    assert hm_symbol(39) == "Abm2"
    assert hm_symbol(41) == "Aba2"
    assert hm_symbol(64) == "Cmca"
    assert hm_symbol(67) == "Cmma"
    assert hm_symbol(68) == "Ccca"
