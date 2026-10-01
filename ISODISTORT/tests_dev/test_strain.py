from __future__ import annotations

import numpy as np
import pytest

from isocore.distortion.strain import decompose_homogeneous_strain


def test_homogeneous_strain_identity_is_zero() -> None:
    parent = np.diag([4.0, 5.0, 6.0])
    basis = np.array([[0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])

    result = decompose_homogeneous_strain(parent, basis, basis @ parent)

    assert result.voigt_engineering == pytest.approx(np.zeros(6), abs=1e-14)
    assert result.relative_metric_residual < 1e-14


def test_homogeneous_strain_matches_official_eual4_f02() -> None:
    """Official 6.12.2 ISOVIZ mode vectors independently fix this tensor.

    F02 amplitudes/mode vectors give applied Voigt values approximately
    (0, 0.08508, -0.00036, 0.13694, 0, 0).  The uploaded daughter was made by
    scaling child vector a by 1.08; CIF rounding accounts for the last digits.
    """

    parent = np.diag([4.402, 4.402, 11.163])
    basis = np.array([[0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])
    distorted = basis @ parent
    distorted[0] *= 1.08

    result = decompose_homogeneous_strain(parent, basis, distorted)
    official_from_isoviz = np.array(
        [
            0.0,
            0.06016 * 0.70711 + (-0.06016) * (-0.70711),
            -0.00036,
            0.09683 * 1.41421,
            0.0,
            0.0,
        ]
    )

    assert result.voigt_engineering == pytest.approx(
        official_from_isoviz, abs=6e-6
    )
    assert result.relative_metric_residual < 1e-14


def test_homogeneous_strain_uses_metric_not_cif_cartesian_orientation() -> None:
    parent = np.diag([4.0, 5.0, 6.0])
    basis = np.eye(3)
    strain_plus_identity = np.array(
        [[1.02, 0.01, 0.0], [0.01, 0.98, 0.02], [0.0, 0.02, 1.03]]
    )
    distorted = basis @ strain_plus_identity @ parent
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])

    result = decompose_homogeneous_strain(parent, basis, distorted @ rotation)

    assert result.strain_plus_identity == pytest.approx(
        strain_plus_identity, abs=1e-12
    )


def test_homogeneous_strain_rejects_singular_basis() -> None:
    with pytest.raises(ValueError, match="basis must be nonsingular"):
        decompose_homogeneous_strain(np.eye(3), np.zeros((3, 3)), np.eye(3))
