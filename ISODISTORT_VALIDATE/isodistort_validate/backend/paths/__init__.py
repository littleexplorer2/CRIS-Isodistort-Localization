"""Public fixed-directory pairing API; attributes reflect the shared path state."""

from . import service

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
]


def __getattr__(name):
    if name in __all__:
        return getattr(service, name)
    raise AttributeError(name)
