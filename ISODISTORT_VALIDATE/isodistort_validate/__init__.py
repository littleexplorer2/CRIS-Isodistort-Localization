"""CIF semantic comparison of local ISODISTORT output vs official reference files."""

__version__ = "0.4.0"

from .application.batch import run_batch
from .application.single import ComparisonResult, compare_cif
from .application.workspace import (
    BATCH_PAIRING_HINT,
    COMPARE_ROOT,
    DEFAULT_PATTERN,
    IGNORE_ATOM_ORDER_HELP,
    ITEM_DIR,
    TRUE_DIR,
    VALIDATE_ROOT,
    Config,
    ensure_compare_dirs,
    get_config,
    pairing_status,
    resolve_pair,
    warn_unpaired_filenames,
)

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
