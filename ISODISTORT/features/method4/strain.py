"""Recover ISODISTORT homogeneous strain from parent and daughter metrics.

Let ``P`` contain parent direct-lattice vectors as rows and ``B`` be a child
basis in parent fractional coordinates. ISODISTORT's symmetric parent-basis
multiplier is ``M = I + E``. With column basis ``C=P.T`` it acts on the right,
``C -> C @ M``; in the row convention used by pymatgen it acts on the left,
``P -> M @ P``. The child lattice is therefore ``B @ M @ P`` up to an
unobservable rigid Cartesian rotation.

Writing ``G_parent=P@P.T`` and
``H=B^-1@(D@D.T)@B^-T``, the unique symmetric-positive solution obeys
``M @ G_parent @ M = H``. The applied engineering coordinates are
``q=(E11,E22,E33,2E23,2E13,2E12)``. They are parent-lattice-basis coordinates,
not Cartesian tensor components for an oblique parent cell.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from features.method4.strain_modes import (
    CHILD_LATTICE_ACTION,
    COLUMN_LATTICE_ACTION,
    ENGINEERING_VOIGT_ORDER,
    METRIC_EQUATION,
    PARENT_BASIS_COORDINATE_FRAME,
    ROW_LATTICE_ACTION,
    tensor_to_engineering_voigt,
)


@dataclass(frozen=True, slots=True)
class HomogeneousStrainResult:
    """Rotation-free applied strain in the official parent-basis convention."""

    multiplier_parent_basis: np.ndarray
    strain_tensor_parent_basis: np.ndarray
    applied_engineering_q_parent_basis: np.ndarray
    reconstructed_child_metric: np.ndarray
    relative_metric_residual: float
    coordinate_frame: str = PARENT_BASIS_COORDINATE_FRAME
    voigt_order: tuple[str, ...] = ENGINEERING_VOIGT_ORDER
    column_lattice_action: str = COLUMN_LATTICE_ACTION
    row_lattice_action: str = ROW_LATTICE_ACTION
    child_lattice_action: str = CHILD_LATTICE_ACTION
    metric_equation: str = METRIC_EQUATION

    @property
    def strain_plus_identity(self) -> np.ndarray:
        """Deprecated alias for :attr:`multiplier_parent_basis`."""

        warnings.warn(
            "strain_plus_identity is deprecated; use multiplier_parent_basis",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.multiplier_parent_basis

    @property
    def tensor(self) -> np.ndarray:
        """Deprecated alias for :attr:`strain_tensor_parent_basis`."""

        warnings.warn(
            "tensor is deprecated; use strain_tensor_parent_basis",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.strain_tensor_parent_basis

    @property
    def voigt_engineering(self) -> np.ndarray:
        """Deprecated alias for :attr:`applied_engineering_q_parent_basis`."""

        warnings.warn(
            "voigt_engineering is deprecated; use "
            "applied_engineering_q_parent_basis",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.applied_engineering_q_parent_basis


def _matrix3(values: object, description: str) -> np.ndarray:
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
    parent_lattice: object,
    parent_to_child_basis: object,
    distorted_child_lattice: object,
) -> HomogeneousStrainResult:
    """Recover the official applied parent-basis strain from a daughter cell.

    ``distorted_child_lattice`` may have any rigid Cartesian orientation. Only
    its row metric is used. The returned engineering ``q`` is the physical
    normalized-vector sum ``sum(amplitude_i * q_unit_i)``. Official CIF
    ``_iso_strain_value`` rows instead use the raw-coordinate sum
    ``sum(amplitude_i * q_raw_i)``; those are intentionally separate fields in
    the export contract.
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
    strained_parent_metric = (
        basis_inverse @ distorted_metric @ basis_inverse.T
    )

    parent_sqrt = _positive_symmetric_sqrt(parent_metric, "parent metric")
    parent_inverse_sqrt = np.linalg.inv(parent_sqrt)
    middle = parent_sqrt @ strained_parent_metric @ parent_sqrt
    multiplier = (
        parent_inverse_sqrt
        @ _positive_symmetric_sqrt(middle, "transformed distorted metric")
        @ parent_inverse_sqrt
    )
    multiplier = (multiplier + multiplier.T) / 2.0
    strain_tensor = multiplier - np.eye(3)
    applied_q = tensor_to_engineering_voigt(strain_tensor)

    reconstructed_lattice = basis @ multiplier @ parent
    reconstructed_metric = reconstructed_lattice @ reconstructed_lattice.T
    metric_scale = max(float(np.linalg.norm(distorted_metric)), 1.0)
    residual = float(
        np.linalg.norm(reconstructed_metric - distorted_metric) / metric_scale
    )
    return HomogeneousStrainResult(
        multiplier_parent_basis=multiplier,
        strain_tensor_parent_basis=strain_tensor,
        applied_engineering_q_parent_basis=applied_q,
        reconstructed_child_metric=reconstructed_metric,
        relative_metric_residual=residual,
    )


__all__ = ["HomogeneousStrainResult", "decompose_homogeneous_strain"]
