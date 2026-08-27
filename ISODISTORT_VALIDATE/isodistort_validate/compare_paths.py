"""Fixed CIF compare directories for ISODISTORT_VALIDATE.

Official reference CIFs go in ``compare/true/``. Local files to check go in
``compare/item/``. Pairing is by relative path under those two folders.
Both folders live next to the package (``ISODISTORT_VALIDATE/compare/``),
not inside ``isodistort_validate/``.
"""
from __future__ import annotations

import sys
from pathlib import Path

VALIDATE_ROOT = Path(__file__).resolve().parents[1]
COMPARE_ROOT = VALIDATE_ROOT / "compare"
ITEM_DIR = COMPARE_ROOT / "item"
TRUE_DIR = COMPARE_ROOT / "true"

BATCH_PAIRING_HINT = (
    "批量比较按相对路径配对：请把官网下载到 compare/true/ 的 CIF "
    "改成与 compare/item/ 中本地文件完全相同的相对路径和文件名。"
)

IGNORE_ATOM_ORDER_HELP = (
    "是否忽略原子排列顺序：\n"
    "  n（默认）= 按 CIF 文件里的原子行号一一对应（第 1 行对第 1 行）。\n"
    "    两边导出格式接近时用这个。元素种类或坐标对不上，或者只是行顺序不同，都会判为不一致。\n"
    "  y = 不按行号，按「同种元素 + 分数坐标足够接近」配对。\n"
    "    同一套原子只是写成不同顺序时（例如 Eu,Al,Al 对 Al,Eu,Al）可选 y。\n"
    "    这不会忽略坐标、元素种类、占据率或磁矩的真实差异。"
)


def ensure_compare_dirs() -> tuple[Path, Path]:
    """Create ``compare/item`` and ``compare/true`` when they are missing."""
    ITEM_DIR.mkdir(parents=True, exist_ok=True)
    TRUE_DIR.mkdir(parents=True, exist_ok=True)
    return ITEM_DIR, TRUE_DIR


def list_relative_cifs(root: Path, pattern: str = "*.cif") -> list[str]:
    if not root.is_dir():
        return []
    return sorted(
        path.relative_to(root).as_posix()
        for path in root.rglob(pattern)
        if path.is_file()
    )


def pairing_status(pattern: str = "*.cif") -> tuple[list[str], list[str], list[str]]:
    """Return ``(paired, item_only, true_only)`` relative CIF paths."""
    ensure_compare_dirs()
    item = set(list_relative_cifs(ITEM_DIR, pattern))
    true = set(list_relative_cifs(TRUE_DIR, pattern))
    return sorted(item & true), sorted(item - true), sorted(true - item)


def format_unpaired_warning(item_only: list[str], true_only: list[str]) -> str | None:
    """Human-readable reminder to rename files so relative paths match."""
    if not item_only and not true_only:
        return None
    lines = [
        "发现文件名对不上。比较按相对路径一一对应，请先修改文件名后再继续。",
        "请把 compare/true/ 中官网下载的 CIF 改成与 compare/item/ 完全相同的相对路径和文件名。",
    ]
    if item_only:
        lines.append("只在 compare/item/ 中：")
        lines.extend(f"  {name}" for name in item_only)
    if true_only:
        lines.append("只在 compare/true/ 中：")
        lines.extend(f"  {name}" for name in true_only)
    return "\n".join(lines)


def warn_unpaired_filenames(*, file=None, pattern: str = "*.cif") -> bool:
    """Print a rename reminder as soon as unpaired CIF names are found.

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
    ensure_compare_dirs()
    text = str(relative_path).strip().replace("\\", "/")
    if not text:
        raise ValueError("relative CIF path is empty")
    rel = Path(text)
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError(
            "CIF path must be relative to compare/item and compare/true "
            f"(got {relative_path!r})"
        )
    return ITEM_DIR / rel, TRUE_DIR / rel
