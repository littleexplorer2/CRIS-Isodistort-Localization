"""Compatibility import; implementation lives in :mod:`backend.input_cif.cif_io`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.input_cif.cif_io")
