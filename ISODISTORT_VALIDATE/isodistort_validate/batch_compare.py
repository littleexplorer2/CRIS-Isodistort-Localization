"""Compatibility imports; use application.batch or frontend.cli for new code."""

from .application.batch import run_batch
from .frontend.cli.batch import _build_parser as _build_parser
from .frontend.cli.batch import main

__all__ = ["run_batch", "main"]

if __name__ == "__main__":
    raise SystemExit(main())
