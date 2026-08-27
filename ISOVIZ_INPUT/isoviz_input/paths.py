"""Fixed local input folders for ISOVIZ_INPUT (not uploaded to git)."""
from __future__ import annotations

from pathlib import Path

ISOVIZ_ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = ISOVIZ_ROOT / "input_content"
DATA_DIR = INPUT_ROOT / "data.csv"
STRUCTURE_DIR = INPUT_ROOT / "subgroup.isoviz"


def ensure_input_content() -> tuple[Path, Path]:
    """Create ``input_content/data.csv`` and ``input_content/subgroup.isoviz`` if missing."""
    for folder in (INPUT_ROOT, DATA_DIR, STRUCTURE_DIR):
        folder.mkdir(parents=True, exist_ok=True)
    return DATA_DIR, STRUCTURE_DIR
