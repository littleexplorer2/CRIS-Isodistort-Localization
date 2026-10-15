"""Compatibility import; implementation lives in :mod:`backend.export.cif_displacive_model`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.export.cif_displacive_model")
