"""Public API for Method 2 mapping and parametric-k mode calculations."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "DistortionMapper": ("features.method2.distortion_mapper", "DistortionMapper"),
    "verify_mapped_microscopic_columns": (
        "features.method2.distortion_mapper",
        "verify_mapped_microscopic_columns",
    ),
    "ParametricModeResult": ("features.method2.superspace", "ParametricModeResult"),
    "compute_parametric_modes": ("features.method2.superspace", "compute_parametric_modes"),
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
