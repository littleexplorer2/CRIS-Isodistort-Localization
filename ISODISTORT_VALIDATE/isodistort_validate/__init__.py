"""CIF semantic comparison of local ISODISTORT output vs official reference files."""

from .batch_compare import run_batch
from .compare_cif import ComparisonResult, compare_cif
from .compare_paths import (
    BATCH_PAIRING_HINT,
    COMPARE_ROOT,
    IGNORE_ATOM_ORDER_HELP,
    ITEM_DIR,
    TRUE_DIR,
    VALIDATE_ROOT,
    ensure_compare_dirs,
    pairing_status,
    resolve_pair,
    warn_unpaired_filenames,
)

__all__ = [
    "BATCH_PAIRING_HINT",
    "COMPARE_ROOT",
    "ComparisonResult",
    "IGNORE_ATOM_ORDER_HELP",
    "ITEM_DIR",
    "TRUE_DIR",
    "VALIDATE_ROOT",
    "compare_cif",
    "ensure_compare_dirs",
    "pairing_status",
    "resolve_pair",
    "run_batch",
    "warn_unpaired_filenames",
]
