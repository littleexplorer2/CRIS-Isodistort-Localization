"""Compatibility import; implementation lives in :mod:`backend.messages.messages`."""
import sys
from importlib import import_module

sys.modules[__name__] = import_module("backend.messages.messages")
