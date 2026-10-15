"""Compatibility alias sharing the canonical paths implementation state."""
import sys

from .backend.paths import resolver as _implementation

sys.modules[__name__] = _implementation
