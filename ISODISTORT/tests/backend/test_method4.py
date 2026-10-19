"""Tests for backend.method4; migrated without changing scientific expectations."""
from __future__ import annotations

from backend.export import StructureExporter
from backend.input_cif import build_supercell
from backend.method4 import DistortionEngine
from backend.method4 import Method4Query
from backend.method4 import ModeDecomposer
from backend.method4.strain import decompose_homogeneous_strain
from backend.method4.strain_modes import ENGINEERING_VOIGT_METRIC
from backend.method4.strain_modes import apply_canonical_strain_basis
from backend.method4.strain_modes import canonical_strain_definitions_from_iso
from backend.method4.strain_modes import compute_homogeneous_strain_modes
from backend.method4.strain_modes import engineering_voigt_to_tensor
from backend.method4.strain_modes import parent_basis_metric_perturbation
from backend.method4.strain_modes import parent_basis_strain_action
from backend.method4.strain_modes import tensor_to_engineering_voigt
from backend.models.iso_mode_models import InvariantDirection
from backend.models.iso_mode_models import MacroscopicTensorBasis
from backend.models.iso_mode_models import MacroscopicTensorBlock
from backend.utils.lattice import centering_primitive_matrix
from backend.utils.lattice import rational_vector
from fractions import Fraction
from itertools import product
from pymatgen.core import Lattice
from pymatgen.core import Structure
from tests.support_backend_export import _macro_component
from tests.support_backend_method4 import _parent
from tests.support_paths import TESTS_ROOT
import itertools
import math
import numpy as np
import pytest
import re


# Original coverage: test_distortion.py.

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
    engine = ModeDecomposer()
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
    result = ModeDecomposer().method_4_decompose(
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
    result = ModeDecomposer().method_4_decompose(
        parent, distorted, {"x": np.array([[1.0, 0.0, 0.0]])}, Method4Query(),
    )
    assert result.raw_coefficients["x"] == pytest.approx(
        coordinate - np.round(coordinate), abs=1e-12,
    )
    assert abs(result.amplitudes["x"]) == pytest.approx(0.5, abs=1e-12)
    assert result.rms_residual == pytest.approx(0.0, abs=1e-12)

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

# Original coverage: test_real_mode_fallback.py.

def _run_column_operation(operation, parent, values, k, basis=None):
    engine = DistortionEngine()
    basis = [1, 1, 1] if basis is None else basis
    if operation == "lift":
        reference, columns = engine.lift_mode_displacements(
            parent, basis, {"mode": values}, k_vector=k,
        )
        return reference, columns["mode"]
    if operation == "single":
        result = engine.generate_single_mode(parent, values, 1.0, basis, k_vector=k)
    else:
        result = engine.generate_modes(parent, basis, parent_displacements=values, k_vector=k)
    reference = build_supercell(parent, basis)
    observed = result.frac_coords - reference.frac_coords
    observed -= np.round(observed)
    return reference, observed

@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("k", [[0, 0, 0], [2, 0, 0]])
def test_explicit_primitive_gamma_requires_centering_covariance(operation, k):
    parent = _parent(centered=True)
    opposite = np.asarray([[0.01, 0, 0], [-0.01, 0, 0]])
    with pytest.raises(ValueError, match=r"does not prove.*translation character"):
        _run_column_operation(operation, parent, opposite, k)

    same = np.abs(opposite)
    _reference, observed = _run_column_operation(operation, parent, same, k)
    assert np.allclose(observed, same, atol=1e-14, rtol=0.0)

@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("basis", [[1, 1, 1], [2, 1, 1]])
def test_unspecified_wavevector_keeps_the_undeclared_parent_column_contract(operation, basis):
    parent = _parent(centered=True)
    values = np.asarray([[0.01, 0, 0], [-0.01, 0, 0]])
    reference, observed = _run_column_operation(operation, parent, values, None, basis)
    parent_positions = reference.frac_coords @ np.diag(basis)
    expected = np.zeros_like(observed)
    expected[:, 0] = np.where(np.isclose(parent_positions[:, 0] % 1.0, 0.5), -0.01, 0.01)
    assert np.allclose(observed @ np.diag(basis), expected, atol=1e-14, rtol=0.0)

@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
def test_unspecified_column_cannot_lose_a_distinct_vector_on_primitive_reduction(operation):
    parent = _parent(centered=True)
    basis = np.asarray(centering_primitive_matrix("I"), dtype=float)
    opposite = np.asarray([[0.01, 0, 0], [-0.01, 0, 0]])
    with pytest.raises(ValueError, match="displacement column is not periodic"):
        _run_column_operation(operation, parent, opposite, None, basis)

@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("k", [None, [0, 0, Fraction(1, 2)]])
def test_child_translation_lattice_must_be_a_parent_sublattice(operation, k):
    # k·B is integral and det(B)=1, but a/2 is not a translation of this P parent.
    basis = np.diag([0.5, 1, 2])
    with pytest.raises(ValueError, match="not a sublattice"):
        _run_column_operation(operation, _parent(), np.asarray([[0.01, 0, 0]]), k, basis)

@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("k", [None, [0, 0, 0]])
def test_valid_centered_primitive_child_preserves_atom_density_and_vectors(operation, k):
    parent = _parent(centered=True)
    basis = np.asarray(centering_primitive_matrix("I"), dtype=float)
    values = np.asarray([[0.01, 0, 0], [0.01, 0, 0]])
    reference, observed = _run_column_operation(operation, parent, values, k, basis)
    assert len(reference) == 1
    assert reference.volume == pytest.approx(parent.volume / 2)
    assert np.allclose(
        observed @ reference.lattice.matrix,
        values[:1] @ parent.lattice.matrix,
        atol=1e-14, rtol=0.0,
    )

@pytest.mark.parametrize("component", [
    1e-13, -1e-13, 1.0 + 1e-13, 1.0 - 1e-13,
    0.5 + 1e-13, 0.5 - 1e-13,
])
def test_float_wavevector_is_not_rounded_across_a_real_phase_boundary(component):
    k = rational_vector([0, 0, component])
    assert k[2] == Fraction(str(component))
    assert (2 * k[2]).denominator != 1
    with pytest.raises(ValueError, match="Non-self-conjugate wavevector"):
        DistortionEngine().generate_single_mode(
            _parent(), np.asarray([[0.01, 0.0, 0.0]]), 1.0, [1, 1, 2],
            k_vector=[0, 0, component],
        )

@pytest.mark.parametrize("operation", ["single", "mixed", "lift"])
@pytest.mark.parametrize("centered,k,basis", [
    (False, [0, 0, Fraction(1, 4)], [1, 1, 4]),
    (True, [Fraction(1, 2), 0, 0], [2, 1, 1]),
])
def test_missing_independent_complex_phase_is_rejected(operation, centered, k, basis):
    parent = _parent(centered)
    values = np.tile([0.01, 0.0, 0.0], (len(parent), 1))
    engine = DistortionEngine()
    with pytest.raises(ValueError, match="Non-self-conjugate wavevector"):
        if operation == "single":
            engine.generate_single_mode(parent, values, 1.0, basis, k_vector=k)
        elif operation == "mixed":
            engine.generate_modes(parent, basis, parent_displacements=values, k_vector=k)
        else:
            engine.lift_mode_displacements(parent, basis, {"mode": values}, k_vector=k)

@pytest.mark.parametrize("opd", [[1, 0], [0, 1]])
def test_opd_components_cannot_be_silently_ignored(opd):
    with pytest.raises(ValueError, match="opd_direction cannot be interpreted"):
        DistortionEngine().generate_single_mode(
            _parent(), np.asarray([[0.01, 0.0, 0.0]]), 1.0, [1, 1, 4],
            k_vector=[0, 0, Fraction(1, 4)], opd_direction=opd,
        )

@pytest.mark.parametrize("origin", [(0, 0, 0), (0.23, -0.17, 0.4)])
def test_self_conjugate_real_mode_has_exact_signs_and_origin_covariance(origin):
    parent = _parent(origin=origin)
    engine = DistortionEngine()
    values = np.asarray([[0.01, 0.0, 0.0]])
    reference, lifted = engine.lift_mode_displacements(
        parent, [1, 1, 2], {"mode": values}, k_vector=[0, 0, Fraction(1, 2)],
    )
    parent_positions = reference.frac_coords @ np.diag([1, 1, 2])
    images = np.round(parent_positions - parent.frac_coords[0]).astype(int)
    expected = np.asarray([(-1) ** int(image[2]) for image in images]) * 0.01
    assert np.array_equal(lifted["mode"][:, 0], expected)
    generated = engine.generate_single_mode(
        parent, values, 1.0, [1, 1, 2], k_vector=[0, 0, Fraction(1, 2)],
    )
    observed = generated.frac_coords - reference.frac_coords
    observed -= np.round(observed)
    assert np.allclose(observed, lifted["mode"], atol=1e-14)

def test_centered_self_conjugate_real_mode_requires_internal_site_character():
    parent = _parent(centered=True)
    engine = DistortionEngine()
    k = [1, 0, 0]  # Real -1 on I centering, +1 on conventional translations.
    good = np.asarray([[0.01, 0, 0], [-0.01, 0, 0]])
    reference, lifted = engine.lift_mode_displacements(
        parent, [1, 1, 1], {"mode": good}, k_vector=k,
    )
    assert len(reference) == 2
    assert np.array_equal(lifted["mode"], good)
    with pytest.raises(ValueError, match=r"does not prove.*translation character"):
        engine.generate_modes(
            parent, parent_displacements=np.abs(good), k_vector=k,
        )

def test_sign_only_wavevector_requires_a_periodic_child_cell():
    with pytest.raises(ValueError, match="not periodic"):
        DistortionEngine().generate_single_mode(
            _parent(), np.asarray([[0.01, 0, 0]]), k_vector=[0, 0, Fraction(1, 2)],
        )

def test_nonzero_small_wavevector_is_not_fuzzed_into_gamma():
    with pytest.raises(ValueError, match="Non-self-conjugate wavevector"):
        DistortionEngine().generate_single_mode(
            _parent(), np.asarray([[0.01, 0, 0]]), k_vector=[0, 0, "1/10000000000"],
        )

# Original coverage: test_robustness.py.

# --- from test_robustness.py ---

def test_negative_amplitude_is_allowed():
    """负振幅：位移方向翻转，不崩溃（合法数学输入）。"""
    struct = Structure(Lattice.cubic(5.0), ["A", "B"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    disp = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    engine = DistortionEngine()
    d = engine.generate_single_mode(struct, disp, amplitude=-0.1)
    assert len(d) == 2
    # 原子 0 沿 -x 位移（≠ 原坐标）
    assert abs(d[0].frac_coords[0]) > 1e-6

def test_wrong_sized_displacements_raise():
    """位移向量与原子数不一致必须报错，不能静默。"""
    struct = Structure(Lattice.cubic(5.0), ["A", "B"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    engine = DistortionEngine()
    with pytest.raises(ValueError):
        engine.generate_single_mode(struct, np.zeros((3, 3)), amplitude=0.1)

def test_extreme_amplitude_no_crash():
    """极大振幅（1e6）：数值上允许，但必须不崩溃、不产生 NaN。"""
    struct = Structure(Lattice.cubic(5.0), ["Fe", "Fe"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    disp = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    d = DistortionEngine().generate_single_mode(struct, disp, amplitude=1e6)
    assert np.all(np.isfinite(d.frac_coords))
    assert len(d) == 2

def test_tiny_amplitude_returns_parent_within_tolerance():
    """极小振幅（1e-12）：位移在容差内，结构 ≈ 母相。"""
    struct = Structure(Lattice.cubic(5.0), ["Fe", "Fe"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    disp = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    d = DistortionEngine().generate_single_mode(struct, disp, amplitude=1e-12)
    diff = d.frac_coords - struct.frac_coords
    diff -= np.round(diff)
    assert np.max(np.abs(diff)) < 1e-6

def test_many_atoms_structure_no_crash():
    """原子数极多的结构（200 原子）：引擎不崩溃。"""
    coords = []
    for i in range(200):
        coords.append([(i % 10) / 10.0, ((i // 10) % 10) / 10.0,
                       (i // 100) / 2.0])
    struct = Structure(Lattice.cubic(20.0), ["Fe"] * 200, coords)
    disp = np.tile([1.0, 0.0, 0.0], (200, 1))
    d = DistortionEngine().generate_single_mode(struct, disp, amplitude=0.1)
    assert len(d) == 200
    assert np.all(np.isfinite(d.frac_coords))

def test_nan_displacement_rejected():
    """含 NaN 的位移向量必须报错（不能静默产生非法结构）。"""
    struct = Structure(Lattice.cubic(5.0), ["Fe"], [[0, 0, 0]])
    disp = np.array([[np.nan, 0.0, 0.0]])
    with pytest.raises((ValueError, FloatingPointError)):
        DistortionEngine().generate_single_mode(struct, disp, amplitude=0.1)

def test_zero_length_displacement_ok():
    """全零位移 + 非零振幅 = 母相（合法，不报错）。"""
    struct = Structure(Lattice.cubic(5.0), ["Fe"], [[0, 0, 0]])
    d = DistortionEngine().generate_single_mode(
        struct, np.zeros((1, 3)), amplitude=0.5
    )
    diff = d.frac_coords - struct.frac_coords
    diff -= np.round(diff)
    assert np.max(np.abs(diff)) < 1e-9

# ------------------------------------------------------------------
# 第三层 §3：可重复性（确定性群论计算，结果必须逐次一致）
# ------------------------------------------------------------------

def test_engine_deterministic_repeat():
    """同一输入两次生成：逐字节一致（确定性）。"""
    struct = Structure(Lattice.cubic(4.0), ["Fe", "O"],
                       [[0, 0, 0], [0.5, 0.5, 0.5]])
    disp = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    engine = DistortionEngine()
    d1 = engine.generate_single_mode(struct, disp, amplitude=0.1)
    d2 = engine.generate_single_mode(struct, disp, amplitude=0.1)
    assert d1.frac_coords.shape == d2.frac_coords.shape
    assert np.allclose(d1.frac_coords, d2.frac_coords, atol=1e-12)
    # 导出 CIF 逐字节一致
    out1 = StructureExporter().to_cif(d1, "rep_test_1")
    out2 = StructureExporter().to_cif(d2, "rep_test_2")
    assert out1.read_bytes() == out2.read_bytes(), \
        "同一输入两次导出的 CIF 必须逐字节一致"

# Original coverage: test_search_methods.py.

def test_method_4_decomposition_recovery():
    engine = ModeDecomposer()

    parent = Structure(Lattice.cubic(5.0), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    distorted = parent.copy()
    distorted.translate_sites([0], [0.03, 0.0, 0.0], frac_coords=True)

    mode_displacements = {
        "GM1+": np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        "R5-": np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
    }

    result = engine.method_4_decompose(
        parent,
        distorted,
        mode_displacements,
        Method4Query(atom_matching_method="nearest-site"),
    )

    assert abs(result.raw_coefficients["GM1+"] - 0.03) < 1e-8
    assert abs(result.amplitudes["GM1+"] - 0.15) < 1e-8
    assert abs(result.parent_cell_amplitudes["GM1+"] - 0.15) < 1e-8
    assert abs(result.mode_normfactors["GM1+"] - 0.2) < 1e-8
    assert abs(result.amplitudes["R5-"]) < 1e-8
    assert result.rms_residual < 1e-10

def test_method_4_applies_origin_shift_and_reports_cartesian_residuals():
    engine = ModeDecomposer()
    parent = Structure(
        Lattice.orthorhombic(4.0, 5.0, 6.0),
        ["Na", "Cl"],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )
    mode = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    shift = np.array([0.125, 0.25, 0.0])
    distorted = Structure(
        parent.lattice,
        list(parent.species),
        np.asarray(parent.frac_coords) + 0.03 * mode + shift,
        coords_are_cartesian=False,
        to_unit_cell=False,
    )

    result = engine.method_4_decompose(
        parent,
        distorted,
        {"GM1+": mode},
        Method4Query(provided_origin_shift=shift),
    )

    assert result.raw_coefficients["GM1+"] == pytest.approx(0.03, abs=1e-10)
    assert result.amplitudes["GM1+"] == pytest.approx(0.12, abs=1e-10)
    assert result.rms_residual == pytest.approx(0.0, abs=1e-10)
    assert result.metadata["residual_unit"] == "angstrom"

def test_method_4_uses_global_cartesian_atom_assignment():
    engine = ModeDecomposer()
    parent = Structure(
        Lattice.orthorhombic(10.0, 5.0, 5.0),
        ["Na", "Na"],
        [[0.10, 0.0, 0.0], [0.20, 0.0, 0.0]],
    )
    distorted = Structure(
        parent.lattice,
        ["Na", "Na"],
        [[0.19, 0.0, 0.0], [0.00, 0.0, 0.0]],
    )
    mode = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])

    result = engine.method_4_decompose(
        parent,
        distorted,
        {"mode": mode},
        Method4Query(),
    )

    assert result.assignments == [1, 0]

def test_method_4_decomposes_homogeneous_lattice_strain():
    engine = ModeDecomposer()
    parent = Structure(Lattice.cubic(5.0), ["Na"], [[0.0, 0.0, 0.0]])
    distorted = Structure(Lattice.cubic(5.4), ["Na"], [[0.0, 0.0, 0.0]])

    result = engine.method_4_decompose(
        parent,
        distorted,
        {"mode": np.ones((1, 3))},
        Method4Query(),
    )

    assert result.strain_applied_engineering_q_parent_basis == pytest.approx(
        {"11": 0.08, "22": 0.08, "33": 0.08, "2*23": 0, "2*13": 0, "2*12": 0}
    )
    assert result.strain_mode_amplitudes == {}
    assert result.strain_raw_coordinate_sum_parent_basis is None
    assert result.metadata["strain_mode_status"] == "unresolved_missing_canonical_basis"
    assert result.metadata["strain_reconstruction_relative_metric_residual"] < 1e-14

def test_method_4_matches_eual4_f02_canonical_strain_contract():
    engine = ModeDecomposer()
    parent_lattice = np.diag([4.402, 4.402, 11.163])
    basis = np.array(
        [[0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]
    )
    labels = [
        "GM1+strain_1(a)",
        "GM1+strain_2(a)",
        "GM2+strain(a)",
        "GM4+strain(a)",
        "GM5+strain(a)",
        "GM5+strain(b)",
    ]
    irreps = ["GM1+", "GM1+", "GM2+", "GM4+", "GM5+", "GM5+"]
    q_raw = np.array(
        [
            [1, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [1, -1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
        ],
        dtype=float,
    )
    normfactors = np.array(
        [1 / np.sqrt(2), 1, 1 / np.sqrt(2), np.sqrt(2), np.sqrt(2), np.sqrt(2)]
    )
    q_unit = q_raw * normfactors[:, None]
    expected_amplitudes = np.array([0.06016, -0.00036, -0.06016, 0, 0.09683, 0])
    expected_raw_sum = expected_amplitudes @ q_raw
    expected_applied_q = expected_amplitudes @ q_unit
    multiplier = np.eye(3) + np.array(
        [
            [expected_applied_q[0], expected_applied_q[5] / 2, expected_applied_q[4] / 2],
            [expected_applied_q[5] / 2, expected_applied_q[1], expected_applied_q[3] / 2],
            [expected_applied_q[4] / 2, expected_applied_q[3] / 2, expected_applied_q[2]],
        ]
    )
    reference = Structure(
        Lattice(basis @ parent_lattice), ["Na"], [[0.0, 0.0, 0.0]]
    )
    distorted = Structure(
        Lattice(basis @ multiplier @ parent_lattice),
        ["Na"],
        [[0.0, 0.0, 0.0]],
    )

    result = engine.method_4_decompose(
        reference,
        distorted,
        {"displacive-test-mode": np.array([[1.0, 0.0, 0.0]])},
        Method4Query(
            reference_parent_lattice=parent_lattice.tolist(),
            parent_to_child_basis=basis.tolist(),
            strain_mode_labels=labels,
            strain_mode_irrep_labels=irreps,
            strain_mode_q_raw=q_raw.tolist(),
            strain_mode_q_unit=q_unit.tolist(),
            strain_mode_normfactors=normfactors.tolist(),
        ),
    )

    assert list(result.strain_mode_amplitudes) == labels
    assert list(result.strain_mode_amplitudes.values()) == pytest.approx(
        expected_amplitudes, abs=2.0e-10
    )
    assert result.strain_raw_coordinate_sum_parent_basis == pytest.approx(
        dict(
            zip(
                ("11", "22", "33", "2*23", "2*13", "2*12"),
                expected_raw_sum,
                strict=True,
            )
        ),
        abs=2.0e-10,
    )
    assert result.strain_applied_engineering_q_parent_basis == pytest.approx(
        dict(
            zip(
                ("11", "22", "33", "2*23", "2*13", "2*12"),
                expected_applied_q,
                strict=True,
            )
        ),
        abs=2.0e-10,
    )
    assert result.metadata["strain_mode_status"] == "canonical_iso_rank_12"
    assert [mode["label"] for mode in result.strain_modes] == labels
    assert [mode["q_raw"] for mode in result.strain_modes] == pytest.approx(q_raw)
    assert [mode["q_unit"] for mode in result.strain_modes] == pytest.approx(q_unit)

def test_method_4_rejects_basis_incompatible_with_reference_cell():
    engine = ModeDecomposer()
    parent = Structure(Lattice.cubic(5.0), ["Na"], [[0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match="basis is incompatible"):
        engine.method_4_decompose(
            parent,
            parent.copy(),
            {"mode": np.ones((1, 3))},
            Method4Query(
                reference_parent_lattice=np.eye(3).tolist(),
                parent_to_child_basis=np.eye(3).tolist(),
            ),
        )

def test_method_4_robust_threshold_is_cartesian_angstrom():
    engine = ModeDecomposer()
    parent = Structure(
        Lattice.orthorhombic(10.0, 2.0, 2.0),
        ["Na"],
        [[0.0, 0.0, 0.0]],
    )
    distorted = Structure(parent.lattice, ["Na"], [[0.06, 0.0, 0.0]])

    with pytest.raises(ValueError, match="Cannot match"):
        engine.method_4_decompose(
            parent,
            distorted,
            {"mode": np.ones((1, 3))},
            Method4Query(
                atom_matching_method="robust",
                robust_distance_threshold=0.25,
            ),
        )

def test_method_4_rejects_rank_deficient_mode_basis():
    engine = ModeDecomposer()
    parent = Structure(Lattice.cubic(5.0), ["Na"], [[0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match="linearly dependent"):
        engine.method_4_decompose(
            parent,
            parent.copy(),
            {
                "mode-a": np.array([[1.0, 0.0, 0.0]]),
                "mode-b": np.array([[2.0, 0.0, 0.0]]),
            },
            Method4Query(),
        )

def test_method_4_reports_official_As_and_Ap_normalization():
    engine = ModeDecomposer()
    parent = Structure(Lattice.cubic(4.0), ["Na"], [[0.0, 0.0, 0.0]])
    mode = np.array([[1.0, 0.0, 0.0]])
    distorted = Structure(parent.lattice, ["Na"], [[0.03, 0.0, 0.0]])

    result = engine.method_4_decompose(
        parent,
        distorted,
        {"mode": mode},
        Method4Query(primitive_cell_multiplicity=2, supercell_size=4),
    )

    assert result.raw_coefficients["mode"] == pytest.approx(0.03)
    assert result.mode_normfactors["mode"] == pytest.approx(np.sqrt(2) / 4)
    assert result.amplitudes["mode"] == pytest.approx(0.12 / np.sqrt(2))
    assert result.parent_cell_amplitudes["mode"] == pytest.approx(0.06 / np.sqrt(2))

def test_method_4_separates_internal_and_max_component_one_coordinates():
    engine = ModeDecomposer()
    parent = Structure(Lattice.cubic(4.0), ["Na"], [[0.0, 0.0, 0.0]])
    internal_mode = np.array([[0.0, 0.0, 1.0 / 6.0]])
    internal_coefficient = 0.18
    distorted = Structure(
        parent.lattice,
        ["Na"],
        [[0.0, 0.0, internal_coefficient / 6.0]],
    )

    result = engine.method_4_decompose(
        parent,
        distorted,
        {"mode": internal_mode},
        Method4Query(primitive_cell_multiplicity=1, supercell_size=4),
    )

    scale = 1.0 / 6.0
    canonical_mode = internal_mode / scale
    canonical_coefficient = internal_coefficient * scale
    normfactor = 1.0 / 4.0
    amplitude = canonical_coefficient / normfactor
    assert result.raw_coefficients["mode"] == pytest.approx(internal_coefficient)
    assert result.metadata["canonical_mode_fractional_scale"]["mode"] == pytest.approx(
        scale
    )
    assert result.metadata["canonical_raw_coefficients"]["mode"] == pytest.approx(
        canonical_coefficient
    )
    assert result.mode_normfactors["mode"] == pytest.approx(normfactor)
    assert result.amplitudes["mode"] == pytest.approx(amplitude)
    assert result.parent_cell_amplitudes["mode"] == pytest.approx(amplitude / 2.0)
    assert np.allclose(
        internal_coefficient * internal_mode,
        canonical_coefficient * canonical_mode,
    )
    assert np.allclose(
        canonical_coefficient * canonical_mode,
        amplitude * normfactor * canonical_mode,
    )

# Original coverage: test_strain.py.

_REPO = TESTS_ROOT.parent.parent

_F02 = _REPO / "output_compare" / "EuAl4 Parent.cif" / "官网" / "Method4" / "F02"

def test_homogeneous_strain_identity_is_zero_and_contract_is_explicit() -> None:
    parent = np.diag([4.0, 5.0, 6.0])
    basis = np.array(
        [[0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]
    )

    result = decompose_homogeneous_strain(parent, basis, basis @ parent)

    assert result.applied_engineering_q_parent_basis == pytest.approx(
        np.zeros(6), abs=1e-14
    )
    assert result.coordinate_frame == "parent_lattice_basis"
    assert result.voigt_order == ("11", "22", "33", "2*23", "2*13", "2*12")
    assert result.column_lattice_action == "C_deformed = C_parent @ M"
    assert result.row_lattice_action == "P_deformed = M @ P_parent"
    assert result.child_lattice_action == "L_child = B @ M @ P_parent"
    assert result.metric_equation == "M @ G_parent @ M = H_parent_basis"
    assert result.relative_metric_residual < 1e-14

def _read_f02_official_strain() -> tuple[np.ndarray, np.ndarray]:
    """Return official raw CIF values and normalized IsoVIZ vector sum."""

    cif = (_F02 / "subgroup.cif").read_text(encoding="utf-8")
    raw = np.array(
        [
            float(
                re.search(
                    rf"^\s*{index}\s+E_{index}\s+([-+0-9.]+)\s*$",
                    cif,
                    re.MULTILINE,
                ).group(1)
            )
            for index in range(1, 7)
        ]
    )
    isoviz = (_F02 / "data.isoviz").read_text(encoding="utf-8")
    block = isoviz.split("!strainmodelist", 1)[1].split(
        "!displacivemodelist", 1
    )[0]
    applied = np.zeros(6)
    pairs = re.findall(
        r"^\s*\d+\s+([-+0-9.]+)\s+[-+0-9.]+\s+\d+\s+\S+\s*$"
        r"\n\s*((?:[-+0-9.]+\s+){5}[-+0-9.]+)\s*$",
        block,
        re.MULTILINE,
    )
    assert len(pairs) == 6
    for amplitude_text, vector_text in pairs:
        amplitude = float(amplitude_text)
        vector = np.array([float(value) for value in vector_text.split()])
        applied += amplitude * vector
    return raw, applied

def test_homogeneous_strain_matches_authoritative_eual4_f02_metric_contract() -> None:
    """F02 separates raw CIF coordinates from the applied normalized q."""

    parent = np.diag([4.402, 4.402, 11.163])
    basis = np.array(
        [[0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]
    )
    distorted = basis @ parent
    distorted[0] *= 1.08
    raw_cif, applied_isoviz = _read_f02_official_strain()

    result = decompose_homogeneous_strain(parent, basis, distorted)

    assert raw_cif == pytest.approx(
        [0.0, 0.12032, -0.00036, 0.09683, 0.0, 0.0], abs=5e-8
    )
    assert applied_isoviz == pytest.approx(
        [0.0, 0.0850794752, -0.00036, 0.1369379543, 0.0, 0.0],
        abs=6e-6,
    )
    assert not np.allclose(raw_cif, applied_isoviz, atol=1e-3)
    assert result.applied_engineering_q_parent_basis == pytest.approx(
        applied_isoviz, abs=6e-6
    )
    assert result.strain_tensor_parent_basis[1, 2] == pytest.approx(
        result.applied_engineering_q_parent_basis[3] / 2.0, abs=1e-14
    )
    assert result.relative_metric_residual < 1e-14

def test_oblique_nonidentity_basis_recovers_parent_basis_multiplier() -> None:
    parent = Lattice.from_parameters(4.2, 5.1, 6.4, 77.0, 84.0, 68.0).matrix
    basis = np.array(
        [[2.0, 0.0, 0.0], [0.0, 1.0, 1.0], [-1.0, 0.0, 1.0]]
    )
    q = np.array([0.02, -0.02, 0.03, 0.04, -0.012, 0.02])
    multiplier = np.eye(3) + engineering_voigt_to_tensor(q)
    angle = 0.37
    rigid_rotation = np.array(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    distorted = basis @ multiplier @ parent @ rigid_rotation

    result = decompose_homogeneous_strain(parent, basis, distorted)

    assert result.multiplier_parent_basis == pytest.approx(multiplier, abs=2e-9)
    assert result.strain_tensor_parent_basis == pytest.approx(
        multiplier - np.eye(3), abs=2e-9
    )
    assert result.applied_engineering_q_parent_basis == pytest.approx(
        q, abs=2e-9
    )
    assert result.reconstructed_child_metric == pytest.approx(
        distorted @ distorted.T, abs=2e-7
    )
    parent_metric = parent @ parent.T
    child_parent_metric = (
        np.linalg.inv(basis)
        @ (distorted @ distorted.T)
        @ np.linalg.inv(basis).T
    )
    assert (
        result.multiplier_parent_basis
        @ parent_metric
        @ result.multiplier_parent_basis
    ) == pytest.approx(child_parent_metric, abs=2e-7)
    assert result.relative_metric_residual < 2e-9

def test_deprecated_generic_aliases_share_parent_basis_objects() -> None:
    result = decompose_homogeneous_strain(
        np.eye(3), np.eye(3), 1.02 * np.eye(3)
    )

    with pytest.deprecated_call(match="multiplier_parent_basis"):
        assert result.strain_plus_identity is result.multiplier_parent_basis
    with pytest.deprecated_call(match="strain_tensor_parent_basis"):
        assert result.tensor is result.strain_tensor_parent_basis
    with pytest.deprecated_call(match="applied_engineering_q_parent_basis"):
        assert (
            result.voigt_engineering
            is result.applied_engineering_q_parent_basis
        )

def test_homogeneous_strain_rejects_singular_basis() -> None:
    with pytest.raises(ValueError, match="basis must be nonsingular"):
        decompose_homogeneous_strain(np.eye(3), np.zeros((3, 3)), np.eye(3))

# Original coverage: test_strain_modes.py.

def _rotation_z(angle_degrees: float) -> np.ndarray:
    angle = math.radians(angle_degrees)
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return np.array(
        [[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]]
    )

def _sign_group() -> list[np.ndarray]:
    return [
        np.diag(signs).astype(float)
        for signs in itertools.product((-1.0, 1.0), repeat=3)
    ]

def _proper_cubic_group() -> list[np.ndarray]:
    rotations = []
    for permutation in itertools.permutations(range(3)):
        permuted = np.zeros((3, 3), dtype=float)
        for row, column in enumerate(permutation):
            permuted[row, column] = 1.0
        for signs in itertools.product((-1.0, 1.0), repeat=3):
            candidate = np.diag(signs) @ permuted
            if float(np.linalg.det(candidate)) > 0.0:
                rotations.append(candidate)
    return rotations

def _fractional_from_cartesian(
    rotations_cartesian: list[np.ndarray], lattice: np.ndarray
) -> list[np.ndarray]:
    columns = lattice.T
    return [
        np.linalg.solve(columns, rotation @ columns)
        for rotation in rotations_cartesian
    ]

@pytest.mark.parametrize(
    ("name", "lattice", "rotations_cartesian", "expected_dimension"),
    (
        (
            "triclinic",
            Lattice.from_parameters(4.1, 5.2, 6.3, 76.0, 83.0, 69.0).matrix,
            [np.eye(3), -np.eye(3)],
            6,
        ),
        (
            "monoclinic",
            Lattice.monoclinic(4.1, 5.2, 6.3, 107.0).matrix,
            [
                np.eye(3),
                np.diag([-1.0, 1.0, -1.0]),
                -np.eye(3),
                np.diag([1.0, -1.0, 1.0]),
            ],
            4,
        ),
        (
            "orthorhombic",
            Lattice.orthorhombic(4.1, 5.2, 6.3).matrix,
            _sign_group(),
            3,
        ),
        (
            "tetragonal",
            Lattice.tetragonal(4.1, 6.3).matrix,
            [_rotation_z(value) for value in (0.0, 90.0, 180.0, 270.0)],
            2,
        ),
        (
            "hexagonal",
            Lattice.hexagonal(4.1, 6.3).matrix,
            [_rotation_z(60.0 * value) for value in range(6)],
            2,
        ),
        (
            "cubic",
            Lattice.cubic(4.1).matrix,
            _proper_cubic_group(),
            1,
        ),
    ),
)
def test_fixed_dimensions_for_general_crystal_systems(
    name: str,
    lattice: np.ndarray,
    rotations_cartesian: list[np.ndarray],
    expected_dimension: int,
) -> None:
    rotations = _fractional_from_cartesian(rotations_cartesian, lattice)
    result = compute_homogeneous_strain_modes(rotations, lattice)

    assert result.diagnostics.fixed_dimension == expected_dimension, name
    assert len(result.modes) == expected_dimension
    assert result.diagnostics.invariance_max_relative_metric_error < 1.0e-12
    assert result.diagnostics.orthonormality_max_abs_error < 1.0e-12
    assert result.diagnostics.fixed_projector_idempotence_max_error < 1.0e-12
    assert result.diagnostics.fixed_projector_constraint_max_error < 1.0e-10

def test_engineering_voigt_round_trip_and_official_normalization() -> None:
    tensor = np.array(
        [[0.3, -0.2, 0.4], [-0.2, 0.1, 0.5], [0.4, 0.5, -0.7]]
    )
    assert engineering_voigt_to_tensor(
        tensor_to_engineering_voigt(tensor)
    ) == pytest.approx(tensor)

    result = compute_homogeneous_strain_modes(
        [np.eye(3), -np.eye(3)], np.eye(3)
    )
    for mode in result.modes:
        assert np.max(np.abs(mode.q_raw)) == pytest.approx(1.0, abs=1.0e-14)
        assert mode.q_unit == pytest.approx(
            mode.normfactor * mode.q_raw, abs=1.0e-14
        )
        assert (
            mode.q_unit @ ENGINEERING_VOIGT_METRIC @ mode.q_unit
        ) == pytest.approx(1.0, abs=1.0e-14)
        assert mode.normfactor == pytest.approx(
            1.0
            / math.sqrt(
                float(mode.q_raw @ ENGINEERING_VOIGT_METRIC @ mode.q_raw)
            ),
            abs=1.0e-14,
        )
        assert mode.irrep_label is None
        assert mode.label_status == "unlabeled"

    # Official Complete modes gives sqrt(2) for a raw unit engineering shear.
    assert result.modes[3].q_raw == pytest.approx([0, 0, 0, 1, 0, 0])
    assert result.modes[3].normfactor == pytest.approx(math.sqrt(2.0))
    assert engineering_voigt_to_tensor(result.modes[3].q_raw)[1, 2] == 0.5

def test_tetragonal_d4h_engineering_directions_match_official_vectors() -> None:
    lattice = Lattice.tetragonal(4.402, 11.163).matrix
    rotations = _fractional_from_cartesian(
        [_rotation_z(value) for value in (0.0, 90.0, 180.0, 270.0)],
        lattice,
    )
    fixed = compute_homogeneous_strain_modes(rotations, lattice)

    # The D4h-invariant GM1+ directions use the same normalization as the
    # official F02 Complete modes and IsoVIZ files.
    assert fixed.modes[0].q_raw == pytest.approx([1, 1, 0, 0, 0, 0])
    assert fixed.modes[0].q_unit == pytest.approx(
        [1 / math.sqrt(2), 1 / math.sqrt(2), 0, 0, 0, 0]
    )
    assert fixed.modes[1].q_unit == pytest.approx([0, 0, 1, 0, 0, 0])

    quarter_turn = parent_basis_strain_action(rotations[1], lattice)
    official_directions = {
        "GM2": np.array([1, -1, 0, 0, 0, 0]) / math.sqrt(2),
        "GM4": np.array([0, 0, 0, 0, 0, math.sqrt(2)]),
        "GM5a": np.array([0, 0, 0, math.sqrt(2), 0, 0]),
        "GM5b": np.array([0, 0, 0, 0, math.sqrt(2), 0]),
    }
    assert quarter_turn @ official_directions["GM2"] == pytest.approx(
        -official_directions["GM2"], abs=1.0e-12
    )
    assert quarter_turn @ official_directions["GM4"] == pytest.approx(
        -official_directions["GM4"], abs=1.0e-12
    )
    assert quarter_turn @ official_directions["GM5a"] == pytest.approx(
        official_directions["GM5b"], abs=1.0e-12
    )
    assert quarter_turn @ official_directions["GM5b"] == pytest.approx(
        -official_directions["GM5a"], abs=1.0e-12
    )

def test_eual4_f02_canonical_macro_modes_close_raw_unit_amplitude_contract() -> None:
    blocks = (
        MacroscopicTensorBlock(
            "GM1+", (_macro_component("xx+yy", [1, 1, 0, 0, 0, 0]),), 1, 1
        ),
        MacroscopicTensorBlock(
            "GM1+", (_macro_component("zz", [0, 0, 1, 0, 0, 0]),), 2, 2
        ),
        MacroscopicTensorBlock(
            "GM2+", (_macro_component("xx-yy", [1, -1, 0, 0, 0, 0]),), 3
        ),
        MacroscopicTensorBlock(
            "GM4+", (_macro_component("xy", [0, 0, 0, 0, 0, 1]),), 4
        ),
        MacroscopicTensorBlock(
            "GM5+",
            (
                _macro_component("yz", [0, 0, 0, 1, 0, 0]),
                _macro_component("-xz", [0, 0, 0, 0, -1, 0]),
            ),
            5,
        ),
    )
    directions = tuple(
        InvariantDirection(label, direction, 1, "P1", 2)
        for label, direction in (
            ("GM1+", "(a)"),
            ("GM2+", "(a)"),
            ("GM4+", "(a)"),
            ("GM5+", "(a,b)"),
        )
    )
    macro = MacroscopicTensorBasis(
        rank_signature="[12]",
        coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
        blocks=blocks,
        status="verified",
        invariant_directions=directions,
    )
    definitions = canonical_strain_definitions_from_iso(macro)
    fixed = compute_homogeneous_strain_modes(
        [np.eye(3)], Lattice.tetragonal(4.402, 11.163).matrix
    )
    result = apply_canonical_strain_basis(fixed, definitions)

    assert [mode.canonical_label for mode in result.modes] == [
        "GM1+strain_1(a)",
        "GM1+strain_2(a)",
        "GM2+strain(a)",
        "GM4+strain(a)",
        "GM5+strain(a)",
        "GM5+strain(b)",
    ]
    expected_raw = np.array(
        [
            [1, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [1, -1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
        ],
        dtype=float,
    )
    assert np.vstack([mode.q_raw for mode in result.modes]) == pytest.approx(
        expected_raw
    )
    amplitudes = np.array([0.06016, -0.00036, -0.06016, 0.0, 0.09683, 0.0])
    raw_sum = expected_raw.T @ amplitudes
    unit_matrix = np.column_stack([mode.q_unit for mode in result.modes])
    applied_q = unit_matrix @ amplitudes
    recovered, *_ = np.linalg.lstsq(unit_matrix, applied_q, rcond=None)
    assert recovered == pytest.approx(amplitudes, abs=1.0e-14)
    assert raw_sum == pytest.approx(
        [0.0, 0.12032, -0.00036, 0.09683, 0.0, 0.0]
    )
    assert applied_q == pytest.approx(
        [0.0, 0.0850791, -0.00036, 0.1369383, 0.0, 0.0], abs=6e-7
    )

    # The ISO parent macro prints the sixth basis vector as -xz, while the
    # final website mode flips it to the canonical first-nonzero-positive +xz.
    # A nonzero synthetic amplitude proves the sign is applied before fitting.
    signed_amplitudes = amplitudes.copy()
    signed_amplitudes[5] = -0.01234
    signed_applied = unit_matrix @ signed_amplitudes
    signed_recovered, *_ = np.linalg.lstsq(
        unit_matrix, signed_applied, rcond=None
    )
    assert result.modes[5].q_raw == pytest.approx([0, 0, 0, 0, 1, 0])
    assert signed_recovered[5] == pytest.approx(-0.01234, abs=1.0e-14)
    assert signed_applied[4] == pytest.approx(
        -0.01234 * math.sqrt(2.0), abs=1.0e-14
    )

def test_cubic_tetragonal_mode_is_scaled_to_website_raw_coordinates() -> None:
    lattice = Lattice.cubic(4.1).matrix
    fixed = compute_homogeneous_strain_modes(
        _fractional_from_cartesian(
            [_rotation_z(value) for value in (0.0, 90.0, 180.0, 270.0)],
            lattice,
        ),
        lattice,
    )
    macro = MacroscopicTensorBasis(
        rank_signature="[12]",
        coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
        blocks=(
            MacroscopicTensorBlock(
                "GM1+",
                (_macro_component("xx+yy+zz", [1, 1, 1, 0, 0, 0]),),
                1,
            ),
            MacroscopicTensorBlock(
                "GM3+",
                (_macro_component("xx+yy-2zz", [1, 1, -2, 0, 0, 0]),),
                2,
            ),
        ),
        status="verified",
        invariant_directions=(
            InvariantDirection("GM1+", "(a)", 123, "P4/mmm", 1),
            InvariantDirection("GM3+", "(a)", 123, "P4/mmm", 1),
        ),
    )

    definitions = canonical_strain_definitions_from_iso(macro)
    result = apply_canonical_strain_basis(fixed, definitions)

    assert result.modes[1].q_raw == pytest.approx(
        [0.5, 0.5, -1.0, 0.0, 0.0, 0.0], abs=1.0e-14
    )
    assert result.modes[1].normfactor == pytest.approx(
        1.0 / math.sqrt(1.5), abs=1.0e-14
    )
    assert result.modes[1].q_unit == pytest.approx(
        np.array([0.5, 0.5, -1.0, 0.0, 0.0, 0.0]) / math.sqrt(1.5),
        abs=1.0e-14,
    )

    unscaled = list(definitions)
    unscaled[1] = type(definitions[1])(
        label=definitions[1].label,
        q_raw=np.array([1.0, 1.0, -2.0, 0.0, 0.0, 0.0]),
        irrep_label=definitions[1].irrep_label,
    )
    with pytest.raises(ValueError, match="max-component-one"):
        apply_canonical_strain_basis(fixed, unscaled)

def test_decimal_invariant_direction_is_preserved_then_validated() -> None:
    lattice = Lattice.tetragonal(4.1, 6.3).matrix
    fixed = compute_homogeneous_strain_modes(
        _fractional_from_cartesian(
            [_rotation_z(value) for value in (0.0, 90.0, 180.0, 270.0)],
            lattice,
        ),
        lattice,
    )
    macro = MacroscopicTensorBasis(
        rank_signature="[12]",
        coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
        blocks=(
            MacroscopicTensorBlock(
                "GMDEC",
                (
                    _macro_component("xx", [1, 0, 0, 0, 0, 0]),
                    _macro_component("yy", [0, 1, 0, 0, 0, 0]),
                ),
                1,
            ),
            MacroscopicTensorBlock(
                "GMZ",
                (_macro_component("zz", [0, 0, 1, 0, 0, 0]),),
                2,
            ),
        ),
        status="partial",
        invariant_directions=(
            InvariantDirection("GMDEC", "(1.732a,1.732a)", 123, "P4/mmm", 1),
            InvariantDirection("GMZ", "(b)", 123, "P4/mmm", 1),
        ),
    )

    definitions = canonical_strain_definitions_from_iso(macro)
    result = apply_canonical_strain_basis(fixed, definitions)

    assert result.modes[0].q_raw == pytest.approx([1, 1, 0, 0, 0, 0])
    assert result.fixed_space_validation_max_error < 1.0e-12

def test_oblique_parent_uses_metric_action_not_fractional_euclidean_action() -> None:
    lattice = Lattice.from_parameters(4.2, 5.1, 6.4, 77.0, 84.0, 68.0).matrix
    inversion = -np.eye(3)
    result = compute_homogeneous_strain_modes([np.eye(3), inversion], lattice)

    assert len(result.modes) == 6
    metric = lattice @ lattice.T
    for mode in result.modes:
        delta_metric = parent_basis_metric_perturbation(mode.q_unit, metric)
        assert inversion.T @ delta_metric @ inversion == pytest.approx(
            delta_metric, abs=1.0e-12
        )

def test_basis_is_deterministic_under_operation_order_and_duplicates() -> None:
    lattice = Lattice.tetragonal(4.1, 6.3).matrix
    rotations = _fractional_from_cartesian(
        [_rotation_z(value) for value in (0.0, 90.0, 180.0, 270.0)],
        lattice,
    )
    expected = compute_homogeneous_strain_modes(rotations, lattice)
    reordered = compute_homogeneous_strain_modes(
        [*reversed(rotations), rotations[0]], lattice
    )

    assert reordered.diagnostics.point_group_order == 4
    assert np.stack([mode.q_unit for mode in reordered.modes]) == pytest.approx(
        np.stack([mode.q_unit for mode in expected.modes]), abs=1.0e-14
    )

def test_rejects_incomplete_rotation_set() -> None:
    lattice = Lattice.tetragonal(4.1, 6.3).matrix
    rotations = _fractional_from_cartesian(
        [np.eye(3), _rotation_z(90.0)], lattice
    )
    with pytest.raises(ValueError, match="not closed"):
        compute_homogeneous_strain_modes(rotations, lattice)

# Original coverage: test_theory_contracts.py.

@pytest.mark.parametrize("amplitude", [0.0, 0.2, 2.0])
def test_legacy_binary_occupancy_contract_cannot_rewrite_parent(amplitude):
    parent = Structure(Lattice.cubic(4), [{"Si": 0.6}], [[0, 0, 0]])
    # The old path changed 0.6 to 1 even at zero amplitude and clipped at 0.01.
    with pytest.raises(NotImplementedError, match="occupational"):
        DistortionEngine().generate_modes(parent, occupancy_patterns=[(np.ones(1), amplitude)])
