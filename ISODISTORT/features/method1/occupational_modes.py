"""Compatibility import; implementation lives in :mod:`backend.method1.occupational_modes`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.method1.occupational_modes")
