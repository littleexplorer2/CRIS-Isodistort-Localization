"""Compatibility alias sharing the canonical validation implementation state."""
import sys

from .backend.validation import checks as _implementation

sys.modules[__name__] = _implementation
