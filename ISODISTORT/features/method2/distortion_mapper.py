"""Compatibility import; implementation lives in :mod:`backend.method2.distortion_mapper`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.method2.distortion_mapper")
