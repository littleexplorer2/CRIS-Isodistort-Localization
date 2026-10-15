"""Compatibility import; implementation lives in :mod:`backend.method4.distortion_engine`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.method4.distortion_engine")
