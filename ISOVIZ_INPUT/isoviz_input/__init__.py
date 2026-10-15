"""Prepare amplitude inputs and dispatch validated structures to IsoVIZ."""

__version__ = "0.4.0"

from .application import (
    LaunchReceipt,
    PipelineRequest,
    PreparedIsoVIZ,
    dispatch_to_isoviz,
    prepare_isoviz,
)
from .backend.amplitudes import (
    ModeAmplitude,
    PatchReport,
    apply_amplitudes,
    list_mode_headers,
    parse_amplitude_csv,
    patch_isoviz_file,
    read_amplitude_csv,
)
from .backend.config_loader import Config, get_config
from .backend.inputs import TextDocument, document_from_text, load_text_document
from .backend.launcher import find_isoviz_launcher, open_isoviz
from .backend.paths import (
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
from .backend.validation import section_scalar, validate_patched_text

__all__ = [
    "DATA_DIR",
    "INPUT_ROOT",
    "STRUCTURE_DIR",
    "Config",
    "LaunchReceipt",
    "ModeAmplitude",
    "PatchReport",
    "PipelineRequest",
    "PreparedIsoVIZ",
    "TextDocument",
    "__version__",
    "apply_amplitudes",
    "best_model_root",
    "dispatch_to_isoviz",
    "document_from_text",
    "ensure_best_model_root",
    "ensure_input_content",
    "find_isoviz_launcher",
    "get_config",
    "list_mode_headers",
    "load_text_document",
    "open_isoviz",
    "parse_amplitude_csv",
    "patch_isoviz_file",
    "prepare_isoviz",
    "read_amplitude_csv",
    "resolve_best_model_csv",
    "resolve_isoviz_path",
    "section_scalar",
    "strip_user_path",
    "validate_patched_text",
]
