"""Compatibility import; implementation lives in :mod:`backend.export.isodistort_isoviz`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.export.isodistort_isoviz")
