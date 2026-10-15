"""Compatibility module alias; new callers use :mod:`application.session`."""
import sys

from application import session as _implementation

sys.modules[__name__ + ".core_api"] = _implementation.core_api
sys.modules[__name__] = _implementation
