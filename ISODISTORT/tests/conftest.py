"""Make ``data_dir.py`` importable when pytest uses ISODISTORT as rootdir."""
from __future__ import annotations

import sys
from pathlib import Path

_DIR = str(Path(__file__).resolve().parent)
if _DIR not in sys.path:
    sys.path.append(_DIR)
