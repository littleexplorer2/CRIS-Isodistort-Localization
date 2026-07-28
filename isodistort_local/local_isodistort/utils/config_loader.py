"""
配置加载器 - 读取 settings.yaml 并设置 ISODATA 环境变量
"""
import os
import yaml
from pathlib import Path


# 项目根目录（local_isodistort 包的上一级）
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
CONFIG_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class Config:
    """全局配置单例"""

    _instance = None
    _loaded = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if self._loaded:
            return
        self._load_config()
        self._setup_environment()
        self._loaded = True

    def _load_config(self):
        with open(CONFIG_PATH, "r", encoding="utf-8") as f:
            self._cfg = yaml.safe_load(f)

    def _setup_environment(self):
        """设置 ISODATA 环境变量，供 iso/findsym 读取数据库"""
        data_dir = self.resolve_path(self._cfg["isobyu"]["data_dir"])
        os.environ["ISODATA"] = str(data_dir)

    def resolve_path(self, rel_path: str) -> Path:
        """将配置中的相对路径解析为绝对路径"""
        p = Path(rel_path)
        if p.is_absolute():
            return p
        return (CONFIG_PATH.parent / rel_path).resolve()

    # ---- 快捷属性 ----

    @property
    def iso_bin(self) -> Path:
        return self.resolve_path(self._cfg["isobyu"]["bin_dir"]) / self._cfg["isobyu"]["iso_bin"]

    @property
    def findsym_bin(self) -> Path:
        return self.resolve_path(self._cfg["isobyu"]["bin_dir"]) / self._cfg["isobyu"]["findsym_bin"]

    @property
    def temp_dir(self) -> Path:
        d = self.resolve_path(self._cfg["runtime"]["temp_dir"])
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def output_dir(self) -> Path:
        d = self.resolve_path(self._cfg["runtime"]["output_dir"])
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def timeout(self) -> int:
        return self._cfg["runtime"]["timeout"]

    @property
    def position_tolerance(self) -> float:
        return self._cfg["defaults"]["position_tolerance"]

    @property
    def defaults(self) -> dict:
        return self._cfg["defaults"]


def get_config() -> Config:
    """获取全局配置实例"""
    return Config()
