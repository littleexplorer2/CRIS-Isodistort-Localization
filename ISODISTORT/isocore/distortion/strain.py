"""Homogeneous-strain decomposition in the ISODISTORT lattice convention.

ISODISTORT stores the three direct-lattice vectors as the columns of a
Cartesian basis matrix ``C`` and applies homogeneous strain as ``C @ M``,
where ``M = I + epsilon`` is symmetric.  For the row-vector lattice
convention used by pymatgen, a child basis ``B`` therefore has lattice
``B @ M @ P`` before an arbitrary rigid Cartesian rotation, where ``P`` is
the conventional parent lattice.

Only lattice metrics are observable in a CIF.  The symmetric positive
solution of ``M Gp M = H`` removes that arbitrary rotation and gives the
same strain tensor used by ISODISTORT.  Engineering shear order is
``(xx, yy, zz, 2yz, 2xz, 2xy)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HomogeneousStrainResult:
    """Rotation-free homogeneous strain recovered from two lattice metrics."""

    strain_plus_identity: np.ndarray
    tensor: np.ndarray
    voigt_engineering: np.ndarray
    reconstructed_child_metric: np.ndarray
    relative_metric_residual: float


def _matrix3(values, description: str) -> np.ndarray:
    matrix = np.asarray(values, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError(f"{description} must be a 3x3 matrix")
    if not bool(np.all(np.isfinite(matrix))):
        raise ValueError(f"{description} must contain only finite values")
    return matrix


def _positive_symmetric_sqrt(matrix: np.ndarray, description: str) -> np.ndarray:
    symmetric = (matrix + matrix.T) / 2.0
    values, vectors = np.linalg.eigh(symmetric)
    scale = max(float(np.max(np.abs(values))), 1.0)
    cutoff = np.finfo(float).eps * scale * 128.0
    if float(np.min(values)) <= cutoff:
        raise ValueError(f"{description} must be symmetric positive definite")
    return (vectors * np.sqrt(values)) @ vectors.T


def decompose_homogeneous_strain(
    parent_lattice,
    parent_to_child_basis,
    distorted_child_lattice,
) -> HomogeneousStrainResult:
    """Recover the ISODISTORT homogeneous strain from a distorted child cell.

    Args:
        parent_lattice: Conventional parent direct-lattice vectors as rows.
        parent_to_child_basis: Child direct-lattice vectors in parent
            fractional coordinates (the basis entered on Method 4).
        distorted_child_lattice: Distorted child direct-lattice vectors as
            rows.  Its arbitrary Cartesian orientation is discarded.
    """

    parent = _matrix3(parent_lattice, "parent lattice")
    basis = _matrix3(parent_to_child_basis, "parent-to-child basis")
    distorted = _matrix3(distorted_child_lattice, "distorted child lattice")
    if abs(float(np.linalg.det(parent))) <= np.finfo(float).eps:
        raise ValueError("parent lattice must be nonsingular")
    if abs(float(np.linalg.det(basis))) <= np.finfo(float).eps:
        raise ValueError("parent-to-child basis must be nonsingular")
    if abs(float(np.linalg.det(distorted))) <= np.finfo(float).eps:
        raise ValueError("distorted child lattice must be nonsingular")

    parent_metric = parent @ parent.T
    distorted_metric = distorted @ distorted.T
    basis_inverse = np.linalg.inv(basis)
    strained_parent_metric = basis_inverse @ distorted_metric @ basis_inverse.T

    parent_sqrt = _positive_symmetric_sqrt(parent_metric, "parent metric")
    parent_inverse_sqrt = np.linalg.inv(parent_sqrt)
    middle = parent_sqrt @ strained_parent_metric @ parent_sqrt
    strain_plus_identity = (
        parent_inverse_sqrt
        @ _positive_symmetric_sqrt(middle, "transformed distorted metric")
        @ parent_inverse_sqrt
    )
    strain_plus_identity = (
        strain_plus_identity + strain_plus_identity.T
    ) / 2.0
    tensor = strain_plus_identity - np.eye(3)
    voigt = np.array(
        [
            tensor[0, 0],
            tensor[1, 1],
            tensor[2, 2],
            2.0 * tensor[1, 2],
            2.0 * tensor[0, 2],
            2.0 * tensor[0, 1],
        ],
        dtype=float,
    )

    reconstructed_lattice = basis @ strain_plus_identity @ parent
    reconstructed_metric = reconstructed_lattice @ reconstructed_lattice.T
    metric_scale = max(float(np.linalg.norm(distorted_metric)), 1.0)
    residual = float(
        np.linalg.norm(reconstructed_metric - distorted_metric) / metric_scale
    )
    return HomogeneousStrainResult(
        strain_plus_identity=strain_plus_identity,
        tensor=tensor,
        voigt_engineering=voigt,
        reconstructed_child_metric=reconstructed_metric,
        relative_metric_residual=residual,
    )


__all__ = ["HomogeneousStrainResult", "decompose_homogeneous_strain"]
