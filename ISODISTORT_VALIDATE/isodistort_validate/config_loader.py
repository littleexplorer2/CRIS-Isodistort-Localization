"""Load ISODISTORT_VALIDATE/config/settings.yaml."""
from __future__ import annotations

from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class Config:
    """Subproject settings (singleton). Paths are resolved from ``config/``."""

    _instance: Config | None = None
    _loaded = False

    def __new__(cls) -> Config:
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self) -> None:
        if self._loaded:
            return
        self._load_config()
        self._loaded = True

    def _load_config(self) -> None:
        if not CONFIG_PATH.is_file():
            raise FileNotFoundError(f"Missing config file: {CONFIG_PATH}")
        with CONFIG_PATH.open(encoding="utf-8") as handle:
            payload = yaml.safe_load(handle) or {}
        if not isinstance(payload, dict):
            raise ValueError(f"Invalid YAML mapping in {CONFIG_PATH}")
        self._cfg = payload

    def resolve_path(self, rel_path: str) -> Path:
        path = Path(rel_path)
        if path.is_absolute():
            return path
        return (CONFIG_PATH.parent / rel_path).resolve()

    @property
    def project_root(self) -> Path:
        return PROJECT_ROOT

    @property
    def compare_root(self) -> Path:
        return self.resolve_path(str(self._cfg["compare"]["root"]))

    @property
    def item_dir(self) -> Path:
        return self.compare_root / str(self._cfg["compare"]["item_dir"])

    @property
    def true_dir(self) -> Path:
        return self.compare_root / str(self._cfg["compare"]["true_dir"])

    @property
    def pattern(self) -> str:
        return str(self._cfg["compare"].get("pattern", "*.cif"))

    @property
    def lattice_tolerance(self) -> float:
        return float(self._cfg["defaults"]["lattice_tolerance"])

    @property
    def coordinate_tolerance(self) -> float:
        return float(self._cfg["defaults"]["coordinate_tolerance"])

    @property
    def scalar_tolerance(self) -> float:
        return float(self._cfg["defaults"]["scalar_tolerance"])

    @property
    def ignore_atom_order(self) -> bool:
        return bool(self._cfg["defaults"].get("ignore_atom_order", False))

    @property
    def strict(self) -> bool:
        return bool(self._cfg["defaults"].get("strict", False))


def get_config() -> Config:
    return Config()
