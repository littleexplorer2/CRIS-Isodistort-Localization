"""Local input folders for ISOVIZ_INPUT.

Amplitude CSVs come from the desktop ``Best_Model_Parameters`` tree written by
the GD notebooks (folder name and optional parent path: ``config/settings.yaml``).
Subgroup ``.isoviz`` files are supplied as an absolute path.
"""
from __future__ import annotations

import os
from pathlib import Path

from .config_loader import get_config

_cfg = get_config()
ISOVIZ_ROOT = _cfg.project_root
INPUT_ROOT = _cfg.input_root
DATA_DIR = _cfg.data_dir
STRUCTURE_DIR = _cfg.structure_dir

DESKTOP_FOLDER_NAME = _cfg.best_model_folder_name
_QUOTES = "\"'`“”‘’「」『』"


def user_desktop() -> Path:
    """Return the current user's Desktop, preferring OneDrive Desktop on Windows."""
    env_home = Path(os.environ.get("USERPROFILE") or Path.home())
    candidates = [
        env_home / "OneDrive" / "Desktop",
        env_home / "Desktop",
        Path.home() / "OneDrive" / "Desktop",
        Path.home() / "Desktop",
    ]
    for path in candidates:
        if path.is_dir():
            return path
    fallback = env_home / "Desktop"
    fallback.mkdir(parents=True, exist_ok=True)
    return fallback


def best_model_root() -> Path:
    cfg = get_config()
    parent = cfg.best_model_parent
    if parent is not None:
        return parent / cfg.best_model_folder_name
    return user_desktop() / cfg.best_model_folder_name


def ensure_best_model_root() -> Path:
    root = best_model_root()
    root.mkdir(parents=True, exist_ok=True)
    return root


def ensure_input_content() -> tuple[Path, Path]:
    """Create the historical ``input_content`` folders if missing."""
    for folder in (INPUT_ROOT, DATA_DIR, STRUCTURE_DIR):
        folder.mkdir(parents=True, exist_ok=True)
    return DATA_DIR, STRUCTURE_DIR


def strip_user_path(raw: str) -> str:
    """Strip wrapping quotes (ASCII, curly, or CJK) from a pasted path or name."""
    text = (raw or "").strip()
    while len(text) >= 2 and text[0] in _QUOTES and text[-1] in _QUOTES:
        text = text[1:-1].strip()
    if len(text) >= 2 and text[0] in _QUOTES:
        text = text[1:].strip()
    if text.endswith("\\") and not text.endswith("\\\\"):
        # A trailing backslash before a closing quote is a common paste artifact.
        text = text.rstrip("\\")
    return text


def list_named_subfolders(root: Path | None = None) -> list[Path]:
    folder = root if root is not None else ensure_best_model_root()
    if not folder.is_dir():
        return []
    return sorted(
        (path for path in folder.iterdir() if path.is_dir()),
        key=lambda path: path.name.lower(),
    )


def list_csv_files(folder: Path) -> list[Path]:
    if not folder.is_dir():
        return []
    return sorted(folder.glob("*.csv"), key=lambda path: path.name.lower())


def resolve_best_model_folder(folder_name: str, *, root: Path | None = None) -> Path:
    """Resolve a subfolder of ``Best_Model_Parameters`` from a typed name or list index."""
    base = root if root is not None else ensure_best_model_root()
    folder_token = strip_user_path(folder_name)
    if not folder_token:
        raise ValueError("A Best_Model_Parameters folder name is required.")
    folders = list_named_subfolders(base)
    if folder_token.isdigit() and folders:
        index = int(folder_token)
        if 1 <= index <= len(folders):
            return folders[index - 1]
    folder = base / folder_token
    if not folder.is_dir():
        raise FileNotFoundError(f"Best_Model_Parameters folder not found: {folder}")
    return folder


def resolve_best_model_csv(folder_name: str, filename: str, *, root: Path | None = None) -> Path:
    """Resolve ``Best_Model_Parameters/<folder>/<file>.csv`` from typed names."""
    folder = resolve_best_model_folder(folder_name, root=root)
    file_token = strip_user_path(filename)
    if not file_token:
        raise ValueError("A CSV file name is required.")

    files = list_csv_files(folder)
    if file_token.isdigit() and files:
        index = int(file_token)
        if 1 <= index <= len(files):
            return files[index - 1]

    name = file_token
    if not name.lower().endswith(".csv"):
        name = f"{name}.csv"
    path = folder / name
    if not path.is_file():
        raise FileNotFoundError(f"CSV file not found: {path}")
    return path


def resolve_isoviz_path(raw: str) -> Path:
    text = strip_user_path(raw)
    if not text:
        raise ValueError("An absolute .isoviz path is required.")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"IsoVIZ structure file not found: {path}")
    if path.suffix.lower() != ".isoviz":
        raise ValueError(f"Expected a .isoviz file, got: {path}")
    return path
