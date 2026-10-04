"""Compatibility wrapper for the repository-level CRIS installer."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

setup_main = importlib.import_module("setup_cris").main


if __name__ == "__main__":
    print(
        "[compat] Use `py -3.10 setup_cris.py install --project validate`; "
        "forwarding to the unified installer."
    )
    raise SystemExit(setup_main(["install", "--project", "validate", *sys.argv[1:]]))
