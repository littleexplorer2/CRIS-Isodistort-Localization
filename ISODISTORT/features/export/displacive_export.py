"""Compatibility import; implementation lives in :mod:`backend.export.displacive_export`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.export.displacive_export")
