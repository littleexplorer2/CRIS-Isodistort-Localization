"""Method 1 search and the reusable distortion-search primitives it owns.

The public symbols are loaded lazily. Several Method 2/3 modules depend on
Method 1's affine and search primitives, so eager imports here would make a
module's importability depend on whichever feature happened to be imported
first. Symbols owned by other Method packages are deliberately not re-exported
from this namespace.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "AffineEmbedding": ("backend.method1.affine_embeddings", "AffineEmbedding"),
    "AffineOperation": ("backend.method1.affine_embeddings", "AffineOperation"),
    "AffineQuotient": ("backend.method1.affine_embeddings", "AffineQuotient"),
    "ParentAffineGroup": ("backend.method1.affine_embeddings", "ParentAffineGroup"),
    "affine_equivalence": ("backend.method1.affine_embeddings", "affine_equivalence"),
    "build_affine_quotient": ("backend.method1.affine_embeddings", "build_affine_quotient"),
    "embedding_from_identity": ("backend.method1.affine_embeddings", "embedding_from_identity"),
    "enumerate_lifted_point_subgroups": (
        "backend.method1.affine_embeddings",
        "enumerate_lifted_point_subgroups",
    ),
    "enumerate_target_affine_embeddings": (
        "backend.method1.affine_embeddings",
        "enumerate_target_affine_embeddings",
    ),
    "locate_embedding_subgroup": (
        "backend.method1.affine_embeddings",
        "locate_embedding_subgroup",
    ),
    "parent_affine_group": ("backend.method1.affine_embeddings", "parent_affine_group"),
    "DomainGenerator": ("backend.method1.domain_generator", "DomainGenerator"),
    "OccupationalMode": ("backend.method1.occupational_modes", "OccupationalMode"),
    "OccupationalModeGenerator": (
        "backend.method1.occupational_modes",
        "OccupationalModeGenerator",
    ),
    "DEFAULT_DISTORTION_TYPES": (
        "backend.method1.phase_path",
        "DEFAULT_DISTORTION_TYPES",
    ),
    "DISTORTION_TYPES": ("backend.method1.phase_path", "DISTORTION_TYPES"),
    "TYPE_ALIASES": ("backend.method1.phase_path", "TYPE_ALIASES"),
    "PhasePath": ("backend.method1.phase_path", "PhasePath"),
    "normalize_distortion_types": (
        "backend.method1.phase_path",
        "normalize_distortion_types",
    ),
    "require_supported_physical_types": (
        "backend.method1.phase_path", "require_supported_physical_types",
    ),
    "IsoSearchEngine": ("backend.method1.search_methods", "IsoSearchEngine"),
    "Method1Query": ("backend.method1.search_methods", "Method1Query"),
    "Method1ResultItem": ("backend.method1.search_methods", "Method1ResultItem"),
    "Method2Query": ("backend.method1.search_methods", "Method2Query"),
    "Method2Result": ("backend.method1.search_methods", "Method2Result"),
    "Method3Query": ("backend.method1.search_methods", "Method3Query"),
    "Method3ResultItem": ("backend.method1.search_methods", "Method3ResultItem"),
    "Method4Query": ("backend.method1.search_methods", "Method4Query"),
    "Method4Result": ("backend.method1.search_methods", "Method4Result"),
}

__all__ = sorted(_EXPORTS)  # noqa: PLE0605 - lazy export map is the source of truth


def __getattr__(name: str) -> Any:
    """Resolve one documented public symbol without eager feature imports."""
    try:
        module_name, attribute = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
    value = getattr(import_module(module_name), attribute)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted((*globals(), *__all__))
