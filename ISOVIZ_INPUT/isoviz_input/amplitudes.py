"""Compatibility alias sharing the canonical amplitudes implementation state."""
import sys

from .backend.amplitudes import core as _implementation

sys.modules[__name__] = _implementation
