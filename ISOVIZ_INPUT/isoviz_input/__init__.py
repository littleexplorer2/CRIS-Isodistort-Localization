"""Feed GD amplitude CSV files into official IsoVIZ ``.isoviz`` structures."""

__version__ = "0.4.0"

from .amplitudes import (
    ModeAmplitude,
    PatchReport,
    apply_amplitudes,
    list_mode_headers,
    patch_isoviz_file,
    read_amplitude_csv,
)
from .config_loader import Config, get_config
from .launcher import find_isoviz_launcher, open_isoviz
from .paths import (
    DATA_DIR,
    INPUT_ROOT,
    STRUCTURE_DIR,
    best_model_root,
    ensure_best_model_root,
    ensure_input_content,
    resolve_best_model_csv,
    resolve_isoviz_path,
    strip_user_path,
)
from .validation import section_scalar, validate_patched_text

__all__ = [
    "DATA_DIR",
    "INPUT_ROOT",
    "STRUCTURE_DIR",
    "Config",
    "ModeAmplitude",
    "PatchReport",
    "__version__",
    "apply_amplitudes",
    "best_model_root",
    "ensure_best_model_root",
    "ensure_input_content",
    "find_isoviz_launcher",
    "get_config",
    "list_mode_headers",
    "open_isoviz",
    "patch_isoviz_file",
    "read_amplitude_csv",
    "resolve_best_model_csv",
    "resolve_isoviz_path",
    "section_scalar",
    "strip_user_path",
    "validate_patched_text",
]
