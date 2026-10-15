"""Compatibility import; implementation lives in :mod:`application.session.core_api`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("application.session.core_api")
