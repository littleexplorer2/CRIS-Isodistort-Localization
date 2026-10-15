"""Compatibility exports; new callers use :mod:`backend.method3`."""
from backend import method3 as _implementation

__all__ = _implementation.__all__

def __getattr__(name):
    return getattr(_implementation, name)

def __dir__():
    return sorted(set(globals()) | set(dir(_implementation)))
