"""Compatibility alias; implementation lives in backend.config."""

import sys

from .backend.config import service as _implementation

sys.modules[__name__] = _implementation
