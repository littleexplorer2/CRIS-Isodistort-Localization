"""Distortion engine, mapper, occupational modes, self-check, multi-mode."""
from __future__ import annotations

import shutil
import subprocess
from itertools import product
from pathlib import Path

import numpy as np
import numpy.linalg as la
import pytest
from data_dir import experiment_data_dir
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from backend.api import IsoDistort
from backend.models.iso_mode_models import ModeIdentity
from backend.utils.self_check import (
    check_linearity,
    check_mode_orthogonality,
    check_subgroup_rule,
    check_symmetry_conservation,
    check_zero_amplitude,
)
from backend.wrappers import BushMode, DistortionMode, SubgroupInfo
from features.input_cif import build_supercell
from features.method1 import OccupationalModeGenerator
from features.method1.phase_path import (
    DEFAULT_DISTORTION_TYPES,
    normalize_distortion_types,
)
from features.method1.search_methods import IsoSearchEngine, Method4Query
from features.method2 import DistortionMapper
from features.method4 import DistortionEngine

CIFS_DIR = Path(__file__).resolve().parent / "cifs_30"
DATA_DIR = experiment_data_dir()


@pytest.mark.parametrize("gamma", [10.0, 90.0])
@pytest.mark.parametrize("matching", ["nearest-site", "robust"])
@pytest.mark.parametrize("strained", [False, True])
@pytest.mark.parametrize("shift", [(0.0, 0.0, 0.0), (1.125, -2.25, 0.375)])
def test_method4_fits_parent_metric_minimum_image(gamma, matching, strained, shift):
    """A zero residual must recover the shortest physical displacement branch.

    Fractional daughter coordinates already use the undistorted comparison
    cell's basis: homogeneous strain changes its metric, not those axes.
    """
    lattice = Lattice.from_parameters(1.0, 1.0, 1.0, 90.0, 90.0, gamma)
    parent = Structure(lattice, ["Na"], [[0.0, 0.0, 0.0]])
    delta = np.array([0.49, 0.48, 0.0])
    distorted = Structure(
        Lattice.cubic(1.0) if strained else lattice,
        ["Na"],
        [delta + np.asarray(shift)],
        to_unit_cell=False,
    )
    # Independent finite enumeration is sufficient for this one-cell fixture.
    images = np.asarray(list(product(range(-2, 3), repeat=3)))
    candidates = delta + images
    cartesian_candidates = candidates @ lattice.matrix
    distances = np.linalg.norm(cartesian_candidates, axis=1)
    shortest = candidates[np.argmin(distances)]
    if gamma == 10.0:
        assert np.allclose(shortest, [-0.51, 0.48, 0.0], atol=1e-12)
        assert np.min(distances) < np.linalg.norm(delta @ lattice.matrix) / 10.0
    else:
        assert np.allclose(shortest, delta, atol=1e-12)

    modes = {axis: np.eye(3)[index].reshape((1, 3))
             for index, axis in enumerate(("x", "y", "z"))}
    engine = IsoSearchEngine(iso_wrapper=None)
    result = engine.method_4_decompose(
        parent, distorted, modes,
        Method4Query(
            atom_matching_method=matching,
            robust_distance_threshold=0.2 if gamma == 10.0 else 0.75,
            provided_origin_shift=shift,
        ),
    )
    coefficients = np.array([result.raw_coefficients[axis] for axis in modes])
    reconstructed_cartesian = coefficients @ lattice.matrix
    assert result.assignments == [0]
    assert np.allclose(coefficients, shortest, atol=1e-12)
    assert np.linalg.norm(reconstructed_cartesian) == pytest.approx(
        np.min(distances), abs=1e-12,
    )
    assert result.rms_residual == pytest.approx(0.0, abs=1e-12)
    # Generation also reaches the observed site modulo lattice translations.
    regenerated = DistortionEngine().generate_modes(
        parent, parent_displacements=coefficients.reshape((1, 3)),
    )
    endpoint_difference = regenerated.frac_coords - delta
    assert np.allclose(
        endpoint_difference - np.round(endpoint_difference), 0.0, atol=1e-12,
    )


def test_method4_does_not_exactly_fit_wrong_long_branch_in_skew_cell():
    lattice = Lattice.from_parameters(1.0, 1.0, 1.0, 90.0, 90.0, 10.0)
    parent = Structure(lattice, ["Na"], [[0.0, 0.0, 0.0]])
    long_delta = np.array([[0.49, 0.48, 0.0]])
    distorted = Structure(lattice, ["Na"], long_delta)
    result = IsoSearchEngine(iso_wrapper=None).method_4_decompose(
        parent, distorted, {"long_branch": long_delta}, Method4Query(),
    )
    # The wrapped endpoint alone cannot justify an exact fit to a large jump.
    shortest_cartesian = np.array([-0.51, 0.48, 0.0]) @ lattice.matrix
    mode_cartesian = long_delta[0] @ lattice.matrix
    expected_coefficient = np.dot(shortest_cartesian, mode_cartesian) / np.dot(
        mode_cartesian, mode_cartesian,
    )
    expected_residual = shortest_cartesian - expected_coefficient * mode_cartesian
    assert result.raw_coefficients["long_branch"] == pytest.approx(
        expected_coefficient, abs=1e-12,
    )
    assert result.rms_residual == pytest.approx(
        np.linalg.norm(expected_residual) / np.sqrt(3.0), abs=1e-12,
    )
    assert result.rms_residual > 0.04


@pytest.mark.parametrize("coordinate", [0.5, -0.5, 1.5, -1.5])
def test_method4_preserves_orthogonal_half_cell_tie_convention(coordinate):
    lattice = Lattice.orthorhombic(1.0, 2.0, 3.0)
    parent = Structure(lattice, ["Na"], [[0.0, 0.0, 0.0]])
    distorted = Structure(lattice, ["Na"], [[coordinate, 0.0, 0.0]])
    result = IsoSearchEngine(iso_wrapper=None).method_4_decompose(
        parent, distorted, {"x": np.array([[1.0, 0.0, 0.0]])}, Method4Query(),
    )
    assert result.raw_coefficients["x"] == pytest.approx(
        coordinate - np.round(coordinate), abs=1e-12,
    )
    assert abs(result.amplitudes["x"]) == pytest.approx(0.5, abs=1e-12)
    assert result.rms_residual == pytest.approx(0.0, abs=1e-12)


def _wsl_available() -> bool:
    if shutil.which("wsl.exe") is None:
        return False
    try:
        result = subprocess.run(  # noqa: PLW1510
            ["wsl.exe", "--status"],  # noqa: S607
            capture_output=True, text=True, timeout=30,
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, OSError):
        return False


# --- from test_distortion.py ---

def test_default_distortion_modes():
    """默认畸变类型应为 strain 单种（对齐官网：Types 面板默认只勾选 Strain）。"""
    # 默认对齐官网 Types 面板：Strain 勾选 + Displacive 全物种勾选
    assert DEFAULT_DISTORTION_TYPES == ["strain", "displacive"]
    assert normalize_distortion_types(None) == ["strain", "displacive"]
    assert normalize_distortion_types(["strain", "displacement", "strain"]) == [
        "strain",
        "displacive",  # 旧名 displacement 自动映射为 displacive
    ]
    assert normalize_distortion_types("order") == ["occupational"]


@pytest.mark.parametrize(
    "values",
    [
        np.zeros((1, 3)),
        np.zeros((3, 3)),
        np.zeros((2, 2)),
        np.array([[0.0, 0.0, np.nan], [0.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, np.inf], [0.0, 0.0, 0.0]]),
    ],
)
def test_supercell_modes_fail_fast_instead_of_padding_or_truncating(values) -> None:
    with pytest.raises(ValueError, match="Mode 'bad'"):
        IsoDistort._validated_supercell_displacement("bad", values, 2)


def test_supercell_mode_validator_accepts_exact_finite_shape() -> None:
    values = np.arange(6, dtype=float).reshape((2, 3))
    checked = IsoDistort._validated_supercell_displacement("ok", values, 2)
    assert np.array_equal(checked, values)


def test_distortion_engine_single_mode():
    # 简单立方

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["A", "B"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    # 人工构造位移：A 原子沿 x 方向位移
    displacements = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ])

    engine = DistortionEngine()
    distorted = engine.generate_single_mode(
        struct, displacements, amplitude=0.05
    )

    # 验证位移
    dx = distorted[0].frac_coords[0] - struct[0].frac_coords[0]
    assert abs(dx - 0.05) < 1e-6, f"期望位移 0.05，实际 {dx}"

    # B 原子不动
    dy = distorted[1].frac_coords[1] - struct[1].frac_coords[1]
    assert abs(dy) < 1e-6

    print("✅ 单模式畸变测试通过")


def test_distortion_engine_supercell():

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["A"], [[0, 0, 0]])

    displacements = np.array([[0.01, 0.0, 0.0]])
    engine = DistortionEngine()
    super_distorted = engine.generate_single_mode(
        struct, displacements, amplitude=1.0, supercell=[2, 2, 2]
    )

    assert len(super_distorted) == 8
    print("✅ 超胞畸变测试通过")


def test_mixed_mode():

    lattice = Lattice.cubic(5.0)
    struct = Structure(lattice, ["A"], [[0, 0, 0]])

    disp_x = np.array([[1.0, 0.0, 0.0]])
    disp_y = np.array([[0.0, 1.0, 0.0]])

    engine = DistortionEngine()
    total_disp = 0.1 * disp_x + 0.2 * disp_y
    mixed = engine.generate_modes(struct, parent_displacements=total_disp)

    dx = mixed[0].frac_coords[0]
    dy = mixed[0].frac_coords[1]
    assert abs(dx - 0.1) < 1e-6
    assert abs(dy - 0.2) < 1e-6

    print("✅ 多模式混合畸变测试通过")




# --- from test_distortion_mapper.py ---

def _mapper_parent() -> Structure:
    """四方母相：Eu 2a + Al 2e（(0,0,0.38)/(0,0,0.62)）。"""
    return Structure(
        Lattice.tetragonal(4.4, 11.2),
        ["Eu", "Eu", "Al", "Al"],
        [[0, 0, 0], [0.5, 0.5, 0.5], [0, 0, 0.38], [0, 0, 0.62]],
    )


def _mapper_wyckoff() -> list[dict]:
    return [
        {"wyckoff_letter": "a", "multiplicity": 2, "species": "Eu",
         "representative_index": 0, "equivalent_indices": [0, 1]},
        {"wyckoff_letter": "e", "multiplicity": 2, "species": "Al",
         "representative_index": 2, "equivalent_indices": [2, 3]},
    ]


def _mapper_mode() -> DistortionMode:
    mode = DistortionMode(irrep_label="GM1+", opd_symbol="P1",
                          k_point_label="GM")
    mode.bush_modes = [
        BushMode("GM1+", "P1", "a", [0, 0, 0], ["0", "0", "0"], [[0, 0, 0]]),
        # e 位点：两个相反位移的代表（(0,0,z) 与 (0,0,-z)）
        BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "z"], [[0, 0, 1]]),
        BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "-z"], [[0, 0, -1]]),
    ]
    return mode


def test_mapper_parameterized_z_points():
    """(0,0,z) 与 (0,0,-z) 必须解析为不同位置：z=0.38 -> +z，z=0.62 -> -z。"""
    mapper = DistortionMapper()
    result = mapper.map_modes_to_atoms(_mapper_parent(), _mapper_wyckoff(), [_mapper_mode()])
    disp = result["GM1+"]["displacements"]
    assert np.allclose(disp[0], 0) and np.allclose(disp[1], 0)  # Eu 不动
    assert np.allclose(disp[2], [0, 0, 1])   # (0,0,0.38) -> +z
    assert np.allclose(disp[3], [0, 0, -1])  # (0,0,0.62) -> -z


def test_mapper_coordinate_token_parser_consumes_the_complete_expression():
    evaluate = DistortionMapper._eval_token
    assert evaluate("-2a + 1", {"a": 0.25}) == pytest.approx(0.5)
    assert evaluate("x-y+1/2", {"x": 0.2, "y": 0.3}) == pytest.approx(0.4)

    for invalid in ("x*y", "x/2", "sin(x)", "x?1"):
        with pytest.raises(ValueError, match="BUSH coordinate"):
            evaluate(invalid, {"x": 0.2, "y": 0.3})
    with pytest.raises(ValueError, match="unresolved BUSH coordinate parameter"):
        evaluate("q", {})


def test_mapper_uniform_single_rep():
    """单代表行（均匀模式）：同一位移作用于该位点全部原子。"""
    mode = _mapper_mode()
    mode.bush_modes = [
        BushMode("GM1+", "P1", "e", [0, 0, 0], ["0", "0", "z"], [[0, 0, 1]]),
    ]
    mapper = DistortionMapper()
    result = mapper.map_modes_to_atoms(_mapper_parent(), _mapper_wyckoff(), [mode])
    disp = result["GM1+"]["displacements"]
    assert np.allclose(disp[2], [0, 0, 1])
    assert np.allclose(disp[3], [0, 0, 1])


def test_mapper_periodic_distance_uses_true_skew_cell_minimum_image():
    """Component-wise fractional wrapping can miss the closest lattice image."""
    structure = Structure(
        Lattice.from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 10.0),
        ["Fe"],
        [[0.0, 0.0, 0.0]],
    )
    mapper = DistortionMapper()
    mapper._lattice_setup(structure)

    first = np.asarray([0.0, 0.0, 0.0])
    second = np.asarray([0.49, 0.49, 0.0])
    component_wrapped = second - first
    component_wrapped -= np.round(component_wrapped)
    component_wrapped_distance = float(
        np.linalg.norm(component_wrapped @ structure.lattice.matrix)
    )

    distance = mapper._periodic_distance(first, second)
    expected = structure.lattice.get_distance_and_image(first, second)[0]
    assert distance == pytest.approx(expected)
    assert distance < 1.0
    assert component_wrapped_distance > 9.0


def test_bush_supercell_mapping_preserves_explicit_phase_pattern():
    """BUSH representatives, not a guessed primary-k cosine, define phases."""
    parent = Structure(Lattice.cubic(4.0), ["H"], [[0, 0, 0]])
    sites = [{
        "wyckoff_letter": "a",
        "multiplicity": 1,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    mode = DistortionMode(
        irrep_label="X1+",
        wyckoff_site="a",
        amplitude_key="X1+__a__P1(1)__a",
        bush_modes=[
            BushMode("X1+", "P1(1)", "a", [0, 0, 0],
                     ["0", "0", "0"], [[1, 0, 0]]),
            BushMode("X1+", "P1(1)", "a", [1, 0, 0],
                     ["1", "0", "0"], [[-1, 0, 0]]),
        ],
    )

    mapped = DistortionMapper().map_bush_modes_to_supercell(
        parent, sites, [mode], [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
    )[mode.amplitude_key]

    # Parent-fractional arrows +/-1 become child-fractional +/-1/2.
    assert mapped.shape == (2, 3)
    assert np.allclose(sorted(mapped[:, 0]), [-0.5, 0.5])
    assert np.allclose(mapped[:, 1:], 0.0)


def test_bush_supercell_mapping_is_invariant_to_unreduced_translation_basis():
    """A GL(3,Z) basis change must not change nearest-image coverage."""
    parent = Structure(Lattice.cubic(1.0), ["H"], [[0, 0, 0]])
    sites = [{
        "wyckoff_letter": "a",
        "multiplicity": 1,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    mode = DistortionMode(
        irrep_label="GM1+",
        wyckoff_site="a",
        amplitude_key="GM1+__a__P1__a",
        bush_modes=[BushMode(
            "GM1+", "P1", "a", [-0.98, -0.049, 0],
            ["-0.98", "-0.049", "0"], [[1, 0, 0]],
        )],
    )
    # This unimodular, strongly sheared basis generates exactly Z^3.  In its
    # fractional coordinates the representative residual is
    # (0.49, 0.049, 0): component-wise rounding gives 0.981 Å, while the true
    # nearest lattice image is only sqrt(0.02^2 + 0.049^2) Å away.
    translation_basis = [[1, 0, 0], [10, 1, 0], [0, 0, 1]]

    mapped = DistortionMapper().map_bush_modes_to_supercell(
        parent,
        sites,
        [mode],
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        cartesian_tolerance=0.1,
        subgroup_translation_lattice=translation_basis,
    )[mode.amplitude_key]

    assert np.allclose(mapped, [[1, 0, 0]])


def test_bush_supercell_mapping_rejects_incomplete_representatives():
    parent = Structure(Lattice.cubic(4.0), ["H"], [[0, 0, 0]])
    sites = [{
        "wyckoff_letter": "a",
        "multiplicity": 1,
        "species": "H",
        "representative_index": 0,
        "equivalent_indices": [0],
    }]
    mode = DistortionMode(
        irrep_label="X1+",
        wyckoff_site="a",
        amplitude_key="X1+__a__P1(1)__a",
        bush_modes=[BushMode(
            "X1+", "P1(1)", "a", [0, 0, 0],
            ["0", "0", "0"], [[1, 0, 0]],
        )],
    )

    with pytest.raises(ValueError, match="does not cover child atom"):
        DistortionMapper().map_bush_modes_to_supercell(
            parent, sites, [mode], [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
        )


def test_zone_center_bush_rows_repeat_over_centered_parent_translations():
    """A Gamma column printed for one primitive parent cell covers a larger child."""
    parent = Structure(
        Lattice.cubic(4.0),
        ["H"] * 4,
        [[0, 0, 0.2], [0, 0, 0.8], [0.5, 0.5, 0.7], [0.5, 0.5, 0.3]],
    )
    orbit_id = "centered-orbit"
    sites = [{
        "wyckoff_letter": "e",
        "orbit_id": orbit_id,
        "equivalent_indices": [0, 1, 2, 3],
        "standard_representative_parameters": {"z": 0.2},
    }]
    mode = DistortionMode(
        irrep_label="GM",
        wyckoff_site="e",
        wyckoff_orbit_id=orbit_id,
        amplitude_key="gamma-column",
        mode_identity=ModeIdentity(
            # The fixture explicitly contains the I-centering translation
            # (1/2,1/2,1/2); its source identity must declare that lattice.
            parent_sg=229,
            global_irrep="GM",
            k_coordinates=("0", "0", "0"),
            wyckoff_letter="e",
            orbit_id=orbit_id,
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[
            BushMode("GM", "P1", "e", [0, 0, 0],
                     ["0", "0", "z"], [[0, 0, 1]]),
            BushMode("GM", "P1", "e", [0, 0, 0],
                     ["0", "0", "-z"], [[0, 0, -1]]),
        ],
    )
    basis = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mapped = DistortionMapper().map_bush_modes_to_supercell(
        parent, sites, [mode], basis,
    )[mode.amplitude_key]
    child = build_supercell(parent, basis)
    assert mapped.shape == (len(child), 3)
    for index, site in enumerate(child):
        parent_z = (site.frac_coords @ np.asarray(basis))[2] % 1
        expected = 0.5 if np.isclose(parent_z, (0.2, 0.7)).any() else -0.5
        assert np.allclose(mapped[index], [0, 0, expected])


def test_bush_rows_repeat_only_along_translations_invariant_under_k_star():
    """An in-plane boundary k permits c repetition but changes a/b phase."""
    parent = Structure(Lattice.tetragonal(4.0, 8.0), ["H"], [[0, 0, 0]])
    sites = [{"wyckoff_letter": "a", "equivalent_indices": [0]}]
    mode = DistortionMode(
        irrep_label="X",
        wyckoff_site="a",
        amplitude_key="boundary-column",
        mode_identity=ModeIdentity(
            parent_sg=123,
            global_irrep="X",
            k_coordinates=("1/2", "1/2", "0"),
            wyckoff_letter="a",
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[
            BushMode("X", "P1", "a", [0, 0, 0],
                     ["0", "0", "0"], [[1, 0, 0]]),
            BushMode("X", "P1", "a", [1, 0, 0],
                     ["1", "0", "0"], [[-1, 0, 0]]),
            BushMode("X", "P1", "a", [0, 1, 0],
                     ["0", "1", "0"], [[-1, 0, 0]]),
            BushMode("X", "P1", "a", [1, 1, 0],
                     ["1", "1", "0"], [[1, 0, 0]]),
        ],
    )
    basis = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mapped = DistortionMapper().map_bush_modes_to_supercell(
        parent, sites, [mode], basis,
    )[mode.amplitude_key]
    child = build_supercell(parent, basis)
    for index, site in enumerate(child):
        x, y, _z = site.frac_coords @ np.asarray(basis)
        expected = 0.5 * (-1) ** (round(x) + round(y))
        assert np.allclose(mapped[index], [expected, 0, 0])


def test_explicit_multiarm_bush_rows_can_select_one_star_arm():
    """Complete BUSH rows can prove that one arm of a two-arm star is active."""
    parent = Structure(Lattice.tetragonal(4.0, 8.0), ["H"], [[0, 0, 0]])
    sites = [{"wyckoff_letter": "a", "equivalent_indices": [0]}]
    mode = DistortionMode(
        irrep_label="X1+",
        wyckoff_site="a",
        amplitude_key="single-arm-column",
        opd_symbol="C1",
        opd_component="a",
        mode_identity=ModeIdentity(
            parent_sg=123,
            global_irrep="X1+",
            k_coordinates=("0", "1/2", "0"),
            wyckoff_letter="a",
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[
            BushMode("X1+", "C1", "a", [0, 0, 0],
                     ["0", "0", "0"], [[1, 0, 0]]),
            BushMode("X1+", "C1", "a", [1, 0, 0],
                     ["1", "0", "0"], [[-1, 0, 0]]),
            BushMode("X1+", "C1", "a", [0, 1, 0],
                     ["0", "1", "0"], [[1, 0, 0]]),
            BushMode("X1+", "C1", "a", [1, 1, 0],
                     ["1", "1", "0"], [[-1, 0, 0]]),
        ],
    )
    basis = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]
    mapped = DistortionMapper().map_bush_modes_to_supercell(
        parent, sites, [mode], basis,
    )[mode.amplitude_key]
    child = build_supercell(parent, basis)
    for index, site in enumerate(child):
        x, _y, _z = site.frac_coords @ np.asarray(basis)
        assert np.allclose(mapped[index], [0.5 * (-1) ** round(x), 0, 0])


def test_sparse_multiarm_bush_without_column_arm_weights_fails_closed():
    """A star identity and local column label do not identify its active arm.

    At the one printed translation, both legitimate X-star completions
    ``u=(-1)^x e_x`` and ``u=(-1)^y e_x`` have the same arrow.  The C1 OPD
    name and local ``a`` component do not bind this microscopic column to
    either arm, so repeating it would invent missing source information.
    """
    parent = Structure(Lattice.tetragonal(4.0, 8.0), ["H"], [[0, 0, 0]])
    mode = DistortionMode(
        irrep_label="X1+",
        wyckoff_site="a",
        amplitude_key="ambiguous-single-arm-column",
        opd_symbol="C1",
        opd_component="a",
        mode_identity=ModeIdentity(
            parent_sg=123,
            global_irrep="X1+",
            k_coordinates=("0", "1/2", "0"),
            wyckoff_letter="a",
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[BushMode(
            "X1+", "C1", "a", [0, 0, 0],
            ["0", "0", "0"], [[1, 0, 0]],
        )],
    )
    with pytest.raises(ValueError, match="does not cover child atom"):
        DistortionMapper().map_bush_modes_to_supercell(
            parent,
            [{"wyckoff_letter": "a", "equivalent_indices": [0]}],
            [mode],
            [[2, 0, 0], [0, 2, 0], [0, 0, 1]],
        )


def test_integer_conventional_k_does_not_hide_centering_phase():
    """An integer conventional k can still be nonzero on a centered lattice."""
    parent = Structure(
        Lattice.cubic(4.0), ["H", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    mode = DistortionMode(
        irrep_label="boundary",
        wyckoff_site="a",
        amplitude_key="centered-boundary",
        mode_identity=ModeIdentity(
            parent_sg=229,
            global_irrep="boundary",
            k_coordinates=("1", "0", "0"),
            wyckoff_letter="a",
            source="iso_microscopic",
            status="verified",
        ),
        bush_modes=[BushMode(
            "boundary", "P1", "a", [0, 0, 0],
            ["0", "0", "0"], [[1, 0, 0]],
        )],
    )
    with pytest.raises(ValueError, match="does not cover child atom"):
        DistortionMapper().map_bush_modes_to_supercell(
            parent,
            [{"wyckoff_letter": "a", "equivalent_indices": [0, 1]}],
            [mode],
            [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
        )


def test_bush_mapping_uses_standard_representative_for_one_physical_orbit():
    """A centered equivalent input point must use the standard free parameter."""
    parent = Structure(
        Lattice.cubic(4.0),
        ["H", "He"],
        [[0.5, 0.5, 0.25], [0.0, 0.0, 0.1]],
    )
    sites = [
        {
            "orbit_id": "hydrogen-e",
            "wyckoff_letter": "e",
            "multiplicity": 1,
            "species": "H",
            "representative_index": 0,
            "equivalent_indices": [0],
            "standard_representative_frac_coords": [0.0, 0.0, -0.25],
            "standard_representative_symmform": "0,0,z",
            "standard_representative_parameters": {"z": -0.25},
        },
        {
            "orbit_id": "helium-e",
            "wyckoff_letter": "e",
            "multiplicity": 1,
            "species": "He",
            "representative_index": 1,
            "equivalent_indices": [1],
        },
    ]
    mode = DistortionMode(
        irrep_label="N1+",
        wyckoff_site="e",
        wyckoff_orbit_id="hydrogen-e",
        amplitude_key="N1+__e@hydrogen-e__P1__a",
        bush_modes=[BushMode(
            "N1+", "P1", "e", [0, 0, 0],
            ["0", "0", "z"], [[1, 0, 0]],
        )],
    )

    mapped = DistortionMapper().map_bush_modes_to_supercell(
        parent,
        sites,
        [mode],
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        subgroup_translation_lattice=[
            [1, 0, 0], [0, 1, 0], [0.5, 0.5, 0.5],
        ],
    )[mode.amplitude_key]

    assert np.allclose(mapped[0], [1, 0, 0])
    assert np.allclose(mapped[1], [0, 0, 0])


# --- from test_occupational_modes.py ---

def _occ_parent() -> Structure:
    """简单四方晶胞：Eu(0,0,0) + Al(0.25,0.25,0.25)。"""
    return Structure(Lattice.tetragonal(4.0, 8.0),
                     ["Eu", "Al"], [[0, 0, 0], [0.25, 0.25, 0.25]])


def _occ_wyckoff() -> list[dict]:
    return [
        {"wyckoff_letter": "a", "multiplicity": 1, "species": "Eu",
         "representative_index": 0, "equivalent_indices": [0]},
        {"wyckoff_letter": "b", "multiplicity": 1, "species": "Al",
         "representative_index": 1, "equivalent_indices": [1]},
    ]


def _occ_subgroup(basis) -> SubgroupInfo:
    return SubgroupInfo(
        index=0, space_group_number=2, space_group_symbol="P-1",
        basis_vectors=basis, k_point_label="X", irrep_label="X1+",
    )


def test_occupational_mode_doubled_cell():
    """2x1x1 超胞：Al 位点分裂为两类 -> 产生 +1/-1 占据率模式。"""
    parent = _occ_parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _occ_wyckoff(), _occ_subgroup([[2, 0, 0], [0, 1, 0], [0, 0, 1]]),
                         {"Al"})
    assert modes, "倍胞下 Al 位点应产生占据率模式"
    m = modes[0]
    assert m.species == "Al"
    assert m.label == "occ-Al-b"
    # 模式值只含 +1/-1/0，且 +1 与 -1 数量相等（两类各半）
    assert set(np.unique(m.pattern)) <= {1.0, -1.0, 0.0}
    assert np.count_nonzero(m.pattern == 1) == np.count_nonzero(m.pattern == -1)


def test_occupational_mode_no_split_identity_cell():
    """单位基矢（t 子群，无超胞）：位点不分裂 -> 无占据率模式。"""
    parent = _occ_parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _occ_wyckoff(),
                         _occ_subgroup([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), {"Al"})
    assert modes == []


def test_occupational_scope_filters_species():
    """作用域只含 Eu 时，Al 位点不产生模式。"""
    parent = _occ_parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _occ_wyckoff(),
                         _occ_subgroup([[2, 0, 0], [0, 1, 0], [0, 0, 1]]), {"Eu"})
    assert all(m.species == "Eu" for m in modes)


def test_occupancy_exportable_structure():
    """生成的占据率结构可被 pymatgen 构建（部分占据位点）。"""
    parent = _occ_parent()
    gen = OccupationalModeGenerator(tolerance=1e-4)
    modes = gen.generate(parent, _occ_wyckoff(),
                         _occ_subgroup([[2, 0, 0], [0, 1, 0], [0, 0, 1]]), {"Al"})
    engine = DistortionEngine()
    m = modes[0]
    sc = engine.generate_modes(
        parent, m.basis_vectors,
        parent_displacements=None,
        occupancy_patterns=[(m.pattern, 0.2)],
    )
    assert len(sc) == 4  # 2 原子 x 2 倍胞
    occs = {s.species_string for s in sc if "Al" in s.species_string}
    assert any(o != "Al" for o in occs)  # 至少一个 Al 位点变为部分占据


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


# --- from test_multi_mode.py ---

def _load_iso(path: Path) -> IsoDistort:
    iso = IsoDistort(language="en")
    iso.set_distortion_scope({
        "displacive": ["*"], "occupational": [], "strain": [],
        "magnetic": [], "rotational": [],
    })
    iso.load_structure(path)
    return iso


def _subgroup_with_modes(iso: IsoDistort, min_modes: int = 2) -> tuple:
    """找首个可产生 >= min_modes 个位移模式的子群。"""
    candidates = iso.search_method_1(
        distortion_types=["displacive", "strain"],
    )
    for cand in candidates:
        try:
            iso.search_method_2(subgroup_idx=cand.subgroup.index,
                                distortion_type=["displacive"])
        except Exception:  # noqa: BLE001,S112 - 无模式子群跳过
            continue
        if len(iso.mode_displacements) >= min_modes:
            return cand.subgroup, list(iso.mode_displacements.keys())
    raise AssertionError("未找到含 >= 2 个位移模式的子群")



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试")
def test_mixed_mode_linear_superposition(tmp_path):
    """混合畸变 = 各单模式位移的线性叠加（StructureMatcher 原子对应后比对）。

    位移向量以母相原胞为单位定义（mode_displacements）；混合生成的超胞
    原子位移（相对母相，周期约化）应等于各模式位移的相位调制线性叠加。
    子群基矢可能是旋转幺模胞，原子顺序会重排，故用 StructureMatcher
    建立“混合畸变原子 -> 母相原子”的对应后再比对位移。
    """
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = CIFS_DIR / "sg139.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")
    iso = _load_iso(cif)
    subgroup, labels = _subgroup_with_modes(iso, min_modes=2)
    label_a, label_b = labels[:2]

    base_disp = {
        lab: np.asarray(iso.mode_displacements[lab]["displacements"], dtype=float)
        for lab in (label_a, label_b)
    }
    # 若子群超胞因子 > 1，则改用超胞因子 1 的 Γ 点子群验证线性叠加
    if abs(round(float(la.det(subgroup.basis_vectors)), 3)) != 1.0:
        iso2 = _load_iso(cif)
        subs = iso2.list_subgroups_at("GM", "GM1+")
        for sg in subs:
            try:
                iso2.search_method_2(subgroup_idx=sg.index,
                                     distortion_type=["displacive"])
            except Exception:  # noqa: BLE001,S112 - 无模式子群跳过
                continue
            if len(iso2.mode_displacements) >= 2:
                labels2 = list(iso2.mode_displacements)[:2]
                base_disp = {
                    lab: np.asarray(iso2.mode_displacements[lab]["displacements"],
                                    dtype=float)
                    for lab in labels2
                }
                amp_a, amp_b = 0.05, -0.03
                mixed = iso2.generate_mixed_distortion(
                    {labels2[0]: amp_a, labels2[1]: amp_b})
                _assert_linear_superposition(iso2, mixed, base_disp,
                                             labels2[0], labels2[1],
                                             amp_a, amp_b)
                return
        pytest.skip("无超胞因子 1 的多模式子群，线性叠加改由对称性校验覆盖")

    amp_a, amp_b = 0.05, -0.03
    mixed = iso.generate_mixed_distortion(
        {label_a: amp_a, label_b: amp_b}
    )
    _assert_linear_superposition(iso, mixed, base_disp,
                                 label_a, label_b, amp_a, amp_b)


def _assert_linear_superposition(iso, mixed, base_disp, label_a, label_b,
                                 amp_a, amp_b) -> None:
    """按“同物种 + 周期最小镜像”原子对应比对混合畸变位移与线性叠加。

    混合结构以子群超胞坐标表达（超胞 = basis @ 母相格子）；位移以母相
    分数坐标定义，故先把超胞坐标变换回母相坐标（parent = sc @ basis）
    再与期望线性叠加比对。
    """
    expected = amp_a * base_disp[label_a] + amp_b * base_disp[label_b]
    n_parent = len(iso.structure)
    assert len(mixed) == n_parent, \
        f"该子群超胞因子应为 1（{n_parent} -> {len(mixed)}）"

    # 子群基矢（母相格单位）；无超胞时为单位阵
    basis = (np.asarray(iso.phase_path.supercell_basis(), dtype=float)
             if iso.phase_path is not None else np.eye(3))
    if basis.shape != (3, 3):
        basis = np.eye(3)

    pc = np.asarray(iso.structure.frac_coords, dtype=float)
    mc_parent = np.asarray(mixed.frac_coords, dtype=float) @ basis
    max_diff = 0.0
    for i in range(n_parent):
        # 同物种中取最小镜像距离对应的混合原子
        best_j, best_dist = -1, float("inf")
        for j in range(n_parent):
            if mixed[j].species_string != iso.structure[i].species_string:
                continue
            delta = mc_parent[j] - pc[i]
            delta -= np.round(delta)
            dist = float(np.linalg.norm(delta))
            if dist < best_dist:
                best_dist, best_j = dist, j
        delta = mc_parent[best_j] - pc[i]
        delta -= np.round(delta)
        diff = delta - expected[i]
        diff -= np.round(diff)
        max_diff = max(max_diff, float(np.max(np.abs(diff))))
    assert max_diff < 1e-6, \
        f"混合位移与线性叠加偏差 {max_diff:.2e}"



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试")
def test_mixed_mode_symmetry_conserved(tmp_path):
    """多模式混合畸变：spglib 对称性 == 目标子群（叠加不破坏对称性）。"""
    cif = DATA_DIR / "EuAl4 Parent.cif"
    if not cif.exists():
        cif = CIFS_DIR / "sg139.cif"
    if not cif.exists():
        pytest.skip("无可用母相 CIF")
    iso = _load_iso(cif)
    subgroup, labels = _subgroup_with_modes(iso, min_modes=2)

    mixed = iso.generate_mixed_distortion({
        labels[0]: 0.05,
        labels[1]: -0.03,
    })
    sg = SpacegroupAnalyzer(mixed, symprec=1e-3).get_space_group_number()
    assert sg == subgroup.space_group_number, \
        f"混合畸变对称性 #{sg}，应为 #{subgroup.space_group_number}"
    # 原子数 = 母相 × 超胞倍数（基矢行列式）
    det = abs(round(float(la.det(subgroup.basis_vectors)), 3))
    if det >= 1:
        assert len(mixed) == round(det * len(iso.structure)), \
            f"混合畸变原子数 {len(mixed)}，期望 {det * len(iso.structure)}"



@pytest.mark.skipif(not _wsl_available(), reason="WSL 不可用，跳过多模式真实计算测试")
def test_mixed_mode_invalid_label_raises(tmp_path):
    """混合生成传入未知模式标号：明确报错，不静默。"""
    cif = CIFS_DIR / "sg001.cif"
    if not cif.exists():
        pytest.skip("测试 CIF 不存在")
    iso = _load_iso(cif)
    subs = iso.list_subgroups_at("GM", "GM1")
    for sg in subs:
        try:
            iso.search_method_2(subgroup_idx=sg.index,
                                distortion_type=["displacive"])
        except Exception:  # noqa: BLE001,S112 - 无模式子群跳过
            continue
        if iso.mode_displacements:
            break
    else:
        pytest.skip("该母相无可计算位移模式")
    with pytest.raises(ValueError):
        iso.generate_mixed_distortion({"NOT_A_MODE": 0.1})
