"""Compatibility alias sharing the canonical launcher implementation state."""
import sys

from .backend.launcher import runtime as _implementation

sys.modules[__name__] = _implementation
