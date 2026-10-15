"""Compatibility import; implementation lives in :mod:`backend.method1.domain_generator`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.method1.domain_generator")
