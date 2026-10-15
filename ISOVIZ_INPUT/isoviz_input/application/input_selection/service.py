"""Resolve an explicit CSV selection without creating or changing user inputs."""
from pathlib import Path

from ...backend.paths import strip_user_path


def resolve_data_path(raw: str) -> Path:
    path = Path(strip_user_path(raw)).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Data file not found: {path}")
    if path.suffix.lower() != ".csv":
        raise ValueError(f"Expected a .csv file, got: {path}")
    return path
