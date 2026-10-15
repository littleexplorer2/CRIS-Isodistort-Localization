"""Compatibility alias; implementation lives in backend.paths."""

import sys

from .backend.paths import service as _implementation

sys.modules[__name__] = _implementation
