"""Public API for Method 2 mapping and parametric-k mode calculations."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "DistortionMapper": ("backend.method2.distortion_mapper", "DistortionMapper"),
    "verify_mapped_microscopic_columns": ("backend.method2.distortion_mapper", "verify_mapped_microscopic_columns"),
    "ParametricModeResult": ("backend.method2.superspace", "ParametricModeResult"),
    "compute_parametric_modes": ("backend.method2.superspace", "compute_parametric_modes"),
    "Method2Query": ("backend.method2.search", "Method2Query"),
    "Method2Result": ("backend.method2.search", "Method2Result"),
    "Method2Solver": ("backend.method2.search", "Method2Solver"),
    "validate_single_irrep_count": ("backend.method2.search", "validate_single_irrep_count"),
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
