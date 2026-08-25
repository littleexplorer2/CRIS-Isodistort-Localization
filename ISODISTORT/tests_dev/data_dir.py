"""Locate CRIS/experiment_data (falls back to the old Chinese folder name)."""
from __future__ import annotations

from pathlib import Path

_CRIS_ROOT = Path(__file__).resolve().parents[2]


def experiment_data_dir() -> Path:
    for name in ("experiment_data", "实验数据与GD代码"):
        path = _CRIS_ROOT / name
        if path.is_dir():
            return path
    return _CRIS_ROOT / "experiment_data"
