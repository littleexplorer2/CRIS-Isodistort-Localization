"""Compatibility alias sharing the canonical config_loader implementation state."""
import sys

from .backend.config_loader import settings as _implementation

sys.modules[__name__] = _implementation
