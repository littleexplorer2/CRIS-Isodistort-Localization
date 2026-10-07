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
    "AffineEmbedding": ("features.method1.affine_embeddings", "AffineEmbedding"),
    "AffineOperation": ("features.method1.affine_embeddings", "AffineOperation"),
    "AffineQuotient": ("features.method1.affine_embeddings", "AffineQuotient"),
    "ParentAffineGroup": ("features.method1.affine_embeddings", "ParentAffineGroup"),
    "affine_equivalence": ("features.method1.affine_embeddings", "affine_equivalence"),
    "build_affine_quotient": ("features.method1.affine_embeddings", "build_affine_quotient"),
    "embedding_from_identity": ("features.method1.affine_embeddings", "embedding_from_identity"),
    "enumerate_lifted_point_subgroups": (
        "features.method1.affine_embeddings",
        "enumerate_lifted_point_subgroups",
    ),
    "enumerate_target_affine_embeddings": (
        "features.method1.affine_embeddings",
        "enumerate_target_affine_embeddings",
    ),
    "locate_embedding_subgroup": (
        "features.method1.affine_embeddings",
        "locate_embedding_subgroup",
    ),
    "parent_affine_group": ("features.method1.affine_embeddings", "parent_affine_group"),
    "DomainGenerator": ("features.method1.domain_generator", "DomainGenerator"),
    "OccupationalMode": ("features.method1.occupational_modes", "OccupationalMode"),
    "OccupationalModeGenerator": (
        "features.method1.occupational_modes",
        "OccupationalModeGenerator",
    ),
    "DEFAULT_DISTORTION_TYPES": (
        "features.method1.phase_path",
        "DEFAULT_DISTORTION_TYPES",
    ),
    "DISTORTION_TYPES": ("features.method1.phase_path", "DISTORTION_TYPES"),
    "TYPE_ALIASES": ("features.method1.phase_path", "TYPE_ALIASES"),
    "PhasePath": ("features.method1.phase_path", "PhasePath"),
    "normalize_distortion_types": (
        "features.method1.phase_path",
        "normalize_distortion_types",
    ),
    "IsoSearchEngine": ("features.method1.search_methods", "IsoSearchEngine"),
    "Method1Query": ("features.method1.search_methods", "Method1Query"),
    "Method1ResultItem": ("features.method1.search_methods", "Method1ResultItem"),
    "Method2Query": ("features.method1.search_methods", "Method2Query"),
    "Method2Result": ("features.method1.search_methods", "Method2Result"),
    "Method3Query": ("features.method1.search_methods", "Method3Query"),
    "Method3ResultItem": ("features.method1.search_methods", "Method3ResultItem"),
    "Method4Query": ("features.method1.search_methods", "Method4Query"),
    "Method4Result": ("features.method1.search_methods", "Method4Result"),
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
