"""Compatibility imports; use application.single or frontend.cli for new code."""

from .application.single import ComparisonResult, compare_cif
from .frontend.cli.compare import _build_parser as _build_parser
from .frontend.cli.compare import _print_result as _print_result
from .frontend.cli.compare import main

__all__ = ["ComparisonResult", "compare_cif", "main"]

if __name__ == "__main__":
    raise SystemExit(main())
