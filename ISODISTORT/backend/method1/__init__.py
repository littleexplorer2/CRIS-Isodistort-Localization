"""Public Method 1 search, occupational modes and domains."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "DomainGenerator": ("backend.method1.domain_generator", "DomainGenerator"),
    "OccupationalMode": ("backend.method1.occupational_modes", "OccupationalMode"),
    "OccupationalModeGenerator": ("backend.method1.occupational_modes", "OccupationalModeGenerator"),
    "Method1Query": ("backend.method1.search", "Method1Query"),
    "Method1ResultItem": ("backend.method1.search", "Method1ResultItem"),
    "Method1Solver": ("backend.method1.search", "Method1Solver"),
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
