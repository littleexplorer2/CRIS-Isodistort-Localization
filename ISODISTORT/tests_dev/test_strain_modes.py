from __future__ import annotations

import itertools
import math
from fractions import Fraction

import numpy as np
import pytest
from pymatgen.core import Lattice

from isocore.backend.iso_mode_models import (
    InvariantDirection,
    MacroscopicTensorBasis,
    MacroscopicTensorBlock,
    SymbolicTensorComponent,
)
from isocore.distortion.strain_modes import (
    ENGINEERING_VOIGT_METRIC,
    apply_canonical_strain_basis,
    canonical_strain_definitions_from_iso,
    compute_homogeneous_strain_modes,
    engineering_voigt_to_tensor,
    parent_basis_metric_perturbation,
    parent_basis_strain_action,
    tensor_to_engineering_voigt,
)


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


def _macro_component(expression: str, values: list[int]) -> SymbolicTensorComponent:
    return SymbolicTensorComponent(
        expression_raw=expression,
        coefficient_basis=("xx", "yy", "zz", "yz", "xz", "xy"),
        coefficients=tuple(Fraction(value) for value in values),
        quality="exact_rational",
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
