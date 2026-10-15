"""Public configuration API; settings.yaml remains the default source."""

from .service import (
    CONFIG_PATH,
    PROJECT_ROOT,
    Config,
    get_config,
)

__all__ = ["Config", "get_config", "PROJECT_ROOT", "CONFIG_PATH"]
