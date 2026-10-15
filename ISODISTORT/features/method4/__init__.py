"""Compatibility exports; new callers use :mod:`backend.method4`."""
from backend import method4 as _implementation

__all__ = _implementation.__all__

def __getattr__(name):
    return getattr(_implementation, name)

def __dir__():
    return sorted(set(globals()) | set(dir(_implementation)))
