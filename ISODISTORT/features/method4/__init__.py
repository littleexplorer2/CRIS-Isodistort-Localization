"""Public API for Method 4 distortion generation, strain, and decomposition."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "DistortionEngine": ("features.method4.distortion_engine", "DistortionEngine"),
    "HomogeneousStrainResult": ("features.method4.strain", "HomogeneousStrainResult"),
    "decompose_homogeneous_strain": (
        "features.method4.strain",
        "decompose_homogeneous_strain",
    ),
    "CHILD_LATTICE_ACTION": ("features.method4.strain_modes", "CHILD_LATTICE_ACTION"),
    "COLUMN_LATTICE_ACTION": ("features.method4.strain_modes", "COLUMN_LATTICE_ACTION"),
    "ENGINEERING_VOIGT_METRIC": (
        "features.method4.strain_modes",
        "ENGINEERING_VOIGT_METRIC",
    ),
    "ENGINEERING_VOIGT_ORDER": (
        "features.method4.strain_modes",
        "ENGINEERING_VOIGT_ORDER",
    ),
    "METRIC_EQUATION": ("features.method4.strain_modes", "METRIC_EQUATION"),
    "PARENT_BASIS_COORDINATE_FRAME": (
        "features.method4.strain_modes",
        "PARENT_BASIS_COORDINATE_FRAME",
    ),
    "ROW_LATTICE_ACTION": ("features.method4.strain_modes", "ROW_LATTICE_ACTION"),
    "CanonicalStrainModeDefinition": (
        "features.method4.strain_modes",
        "CanonicalStrainModeDefinition",
    ),
    "HomogeneousStrainMode": ("features.method4.strain_modes", "HomogeneousStrainMode"),
    "HomogeneousStrainModeDiagnostics": (
        "features.method4.strain_modes",
        "HomogeneousStrainModeDiagnostics",
    ),
    "HomogeneousStrainModeResult": (
        "features.method4.strain_modes",
        "HomogeneousStrainModeResult",
    ),
    "apply_canonical_strain_basis": (
        "features.method4.strain_modes",
        "apply_canonical_strain_basis",
    ),
    "canonical_strain_definitions_from_iso": (
        "features.method4.strain_modes",
        "canonical_strain_definitions_from_iso",
    ),
    "compute_homogeneous_strain_modes": (
        "features.method4.strain_modes",
        "compute_homogeneous_strain_modes",
    ),
    "engineering_voigt_to_tensor": (
        "features.method4.strain_modes",
        "engineering_voigt_to_tensor",
    ),
    "parent_basis_metric_perturbation": (
        "features.method4.strain_modes",
        "parent_basis_metric_perturbation",
    ),
    "parent_basis_strain_action": (
        "features.method4.strain_modes",
        "parent_basis_strain_action",
    ),
    "tensor_to_engineering_voigt": (
        "features.method4.strain_modes",
        "tensor_to_engineering_voigt",
    ),
}

__all__ = sorted(_EXPORTS)  # noqa: PLE0605 - lazy export map is the source of truth


def __getattr__(name: str) -> Any:
    try:
        module_name, attribute = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
    value = getattr(import_module(module_name), attribute)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted((*globals(), *__all__))
