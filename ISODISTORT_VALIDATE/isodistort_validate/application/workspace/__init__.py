"""Public workspace settings and file-pair selection used by interfaces."""

from __future__ import annotations

from ...backend import paths as cpaths
from ...backend.config import Config, get_config

__all__ = [
    "VALIDATE_ROOT",
    "COMPARE_ROOT",
    "ITEM_DIR",
    "TRUE_DIR",
    "DEFAULT_PATTERN",
    "BATCH_PAIRING_HINT",
    "IGNORE_ATOM_ORDER_HELP",
    "ensure_compare_dirs",
    "list_relative_cifs",
    "pairing_status",
    "format_unpaired_warning",
    "warn_unpaired_filenames",
    "resolve_pair",
    "get_config",
    "Config",
    "default_relative_path",
]


def __getattr__(name):
    if name in cpaths.__all__:
        return getattr(cpaths, name)
    raise AttributeError(name)


def default_relative_path() -> str:
    paired, item_only, true_only = cpaths.pairing_status()
    warning = cpaths.format_unpaired_warning(item_only, true_only)
    if len(paired) == 1:
        return paired[0]
    if not paired:
        raise FileNotFoundError(
            warning
            or (
                "compare/item 与 compare/true 中没有 CIF。"
                "请把官网参考 CIF 放入 compare/true/，把待验证 CIF 放入 compare/item/，"
                "并保持相对路径相同。"
            )
        )
    preview = ", ".join(paired[:8])
    extra = "" if len(paired) <= 8 else f" ... ({len(paired)} files)"
    raise ValueError(
        "compare/item 与 compare/true 中有多对 CIF，请给出相对路径，"
        "或从 CRIS 根目录改用 .\\.venv\\Scripts\\python.exe "
        "ISODISTORT_VALIDATE\\main.py batch。"
        f"可用文件: {preview}{extra}"
    )
