"""Load ISODISTORT_VALIDATE/config/settings.yaml."""

from __future__ import annotations

from pathlib import Path

import yaml

PACKAGE_ROOT = next(
    parent
    for parent in Path(__file__).resolve().parents
    if parent.name == "isodistort_validate" and (parent / "__init__.py").is_file()
)
PROJECT_ROOT = PACKAGE_ROOT.parent
_SOURCE_CONFIG_PATH = PROJECT_ROOT / "config" / "settings.yaml"
_PACKAGED_CONFIG_PATH = PACKAGE_ROOT / "_resources" / "settings.yaml"
CONFIG_PATH = (
    _SOURCE_CONFIG_PATH
    if (PROJECT_ROOT / "pyproject.toml").is_file() and _SOURCE_CONFIG_PATH.is_file()
    else _PACKAGED_CONFIG_PATH
)


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
        # Preserve source-relative defaults after packaging, while explicit
        # CONFIG_PATH overrides continue to resolve beside the selected file.
        base = _SOURCE_CONFIG_PATH.parent if CONFIG_PATH == _PACKAGED_CONFIG_PATH else CONFIG_PATH.parent
        return (base / rel_path).resolve()

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
