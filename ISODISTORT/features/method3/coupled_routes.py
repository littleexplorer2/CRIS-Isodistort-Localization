"""Compatibility import; implementation lives in :mod:`backend.method3.coupled_routes`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.method3.coupled_routes")
