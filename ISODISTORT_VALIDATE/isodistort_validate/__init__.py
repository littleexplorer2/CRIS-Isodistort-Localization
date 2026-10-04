"""CIF semantic comparison of local ISODISTORT output vs official reference files."""

__version__ = "0.4.0"

from .batch_compare import run_batch
from .compare_cif import ComparisonResult, compare_cif
from .compare_paths import (
    BATCH_PAIRING_HINT,
    COMPARE_ROOT,
    DEFAULT_PATTERN,
    IGNORE_ATOM_ORDER_HELP,
    ITEM_DIR,
    TRUE_DIR,
    VALIDATE_ROOT,
    ensure_compare_dirs,
    pairing_status,
    resolve_pair,
    warn_unpaired_filenames,
)
from .config_loader import Config, get_config

__all__ = [
    "BATCH_PAIRING_HINT",
    "COMPARE_ROOT",
    "DEFAULT_PATTERN",
    "IGNORE_ATOM_ORDER_HELP",
    "ITEM_DIR",
    "TRUE_DIR",
    "VALIDATE_ROOT",
    "ComparisonResult",
    "Config",
    "__version__",
    "compare_cif",
    "ensure_compare_dirs",
    "get_config",
    "pairing_status",
    "resolve_pair",
    "run_batch",
    "warn_unpaired_filenames",
]
