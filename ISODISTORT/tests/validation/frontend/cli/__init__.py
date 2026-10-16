"""Noninteractive entry point for agent validation commands."""

from __future__ import annotations

import sys
from .batch import main as batch_main
from .compare import main as compare_main

__all__ = ["main"]


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"-h", "--help", "help"}:
        print("python -m tests.validation {compare|batch} [options]")
        return 0
    if args[0] == "compare":
        return compare_main(args[1:])
    if args[0] == "batch":
        return batch_main(args[1:])
    print(f"Unknown validation command: {args[0]}", file=sys.stderr)
    return 2
