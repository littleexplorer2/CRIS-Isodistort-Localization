"""Load ISOVIZ_INPUT/config/settings.yaml."""
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
    def input_root(self) -> Path:
        return self.resolve_path(str(self._cfg["paths"]["input_root"]))

    @property
    def data_dir(self) -> Path:
        return self.input_root / str(self._cfg["paths"]["data_dirname"])

    @property
    def structure_dir(self) -> Path:
        return self.input_root / str(self._cfg["paths"]["structure_dirname"])

    @property
    def best_model_folder_name(self) -> str:
        return str(self._cfg["paths"].get("best_model_folder_name", "Best_Model_Parameters"))

    @property
    def best_model_parent(self) -> Path | None:
        raw = str(self._cfg["paths"].get("best_model_parent") or "").strip()
        if not raw:
            return None
        return Path(raw).expanduser()

    @property
    def isoviz_env_vars(self) -> list[str]:
        values = self._cfg.get("isoviz", {}).get("env_vars") or ["ISOVIZ", "ISOVIZ_JAR"]
        return [str(name) for name in values]

    @property
    def launcher_names(self) -> list[str]:
        values = self._cfg.get("isoviz", {}).get("launcher_names") or []
        return [str(name) for name in values]

    @property
    def launcher_search_dirs(self) -> list[Path]:
        values = self._cfg.get("isoviz", {}).get("launcher_search_dirs") or ["..", "../.."]
        return [self.resolve_path(str(item)) for item in values]


def get_config() -> Config:
    return Config()
