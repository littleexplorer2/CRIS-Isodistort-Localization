"""Compatibility import; implementation lives in :mod:`backend.export.result_serializer`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.export.result_serializer")
