"""Compatibility import; implementation lives in :mod:`backend.method3.inverse_landau`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.method3.inverse_landau")
