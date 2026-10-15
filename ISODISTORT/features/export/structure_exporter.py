"""Compatibility import; implementation lives in :mod:`backend.export.structure_exporter`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.export.structure_exporter")
