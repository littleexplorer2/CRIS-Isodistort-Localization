"""Read-only CIF pairing for the integrated validation suite.

Official reference CIFs go in ``compare/true/``. Local files to check go in
``compare/item/``. Pairing is by relative path under those two folders.
Paths come from ``tests/validation/config/settings.yaml``.
"""

from __future__ import annotations

import sys
from pathlib import Path

from ..config import get_config

_cfg = get_config()
VALIDATE_ROOT = _cfg.project_root
COMPARE_ROOT = _cfg.compare_root
ITEM_DIR = _cfg.item_dir
TRUE_DIR = _cfg.true_dir
DEFAULT_PATTERN = _cfg.pattern

BATCH_PAIRING_HINT = (
    "批量比较按相对路径配对：compare/true/ 与 compare/item/ 中对应文件的相对路径须一致。"
    "原始参考文件保持只读；如名称不同，请使用 compare --local / --reference 显式配对。"
)

def list_relative_cifs(root: Path, pattern: str | None = None) -> list[str]:
    if pattern is None:
        pattern = DEFAULT_PATTERN
    if not root.is_dir():
        return []
    return sorted(path.relative_to(root).as_posix() for path in root.rglob(pattern) if path.is_file())


def pairing_status(pattern: str | None = None) -> tuple[list[str], list[str], list[str]]:
    """Return ``(paired, item_only, true_only)`` relative CIF paths."""
    if pattern is None:
        pattern = DEFAULT_PATTERN
    item = set(list_relative_cifs(ITEM_DIR, pattern))
    true = set(list_relative_cifs(TRUE_DIR, pattern))
    return sorted(item & true), sorted(item - true), sorted(true - item)


def format_unpaired_warning(item_only: list[str], true_only: list[str]) -> str | None:
    """Describe unmatched paths without modifying reference files."""
    if not item_only and not true_only:
        return None
    lines = [
        "发现文件名对不上。批量比较按相对路径一一对应；请检查配对，不要修改参考原件。",
        "compare/true/ 与 compare/item/ 相对路径不同的文件，可用 compare --local / --reference 显式配对。",
    ]
    if item_only:
        lines.append("只在 compare/item/ 中：")
        lines.extend(f"  {name}" for name in item_only)
    if true_only:
        lines.append("只在 compare/true/ 中：")
        lines.extend(f"  {name}" for name in true_only)
    return "\n".join(lines)


def warn_unpaired_filenames(*, file=None, pattern: str | None = None) -> bool:
    """Print a pairing reminder as soon as unpaired CIF names are found.

    Returns True when a mismatch exists. ``file`` defaults to stderr so JSON
    reports on stdout stay intact.
    """
    _paired, item_only, true_only = pairing_status(pattern)
    text = format_unpaired_warning(item_only, true_only)
    if text is None:
        return False
    print(text, file=sys.stderr if file is None else file)
    return True


def resolve_pair(relative_path: str) -> tuple[Path, Path]:
    """Map a path relative to ``compare/item`` and ``compare/true``."""
    text = str(relative_path).strip().replace("\\", "/")
    if not text:
        raise ValueError("relative CIF path is empty")
    rel = Path(text)
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError(f"CIF path must be relative to compare/item and compare/true (got {relative_path!r})")
    return ITEM_DIR / rel, TRUE_DIR / rel
