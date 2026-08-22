"""
配置加载器 - 读取 settings.yaml 并设置 ISODATA 环境变量
"""
import os
from pathlib import Path

import yaml

# 项目根目录（isodistort 包的上一级）
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
CONFIG_PATH = PROJECT_ROOT / "config" / "settings.yaml"


class Config:
    """全局配置单例
    （创建实例）
    加载yaml文件中的配置
    设置 ISODATA 环境变量"""

    _instance = None
    _loaded = False

    def __new__(cls):
        """重写 __new__ 方法以实现单例模式，这个类只有一个实例：检查instance是否为空，若为空则直接创建一个新实例

        """
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
        """从 settings.yaml 加载配置

        """
        with open(CONFIG_PATH, encoding="utf-8") as f:
            self._cfg = yaml.safe_load(f)

    def _setup_environment(self):
        """读取 iso/findsym 所属的数据库路径信息，并将其设置为 ISODATA 环境变量

        """
        data_dir = self.resolve_path(self._cfg["isobyu"]["data_dir"])
        os.environ["ISODATA"] = str(data_dir)

    def resolve_path(self, rel_path: str) -> Path:
        """将配置中的相对路径解析为绝对路径

        """
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
        """普通子进程调用超时（秒）。"""
        return self._cfg["runtime"]["timeout"]

    @property
    def generation_timeout(self) -> int:
        """子群数据库在线生成超时（秒）。

        对应官网 “Generate isotropy subgroups”：非特殊 k 点的子群数据库
        需在线生成，可能耗时数分钟到数小时，故单独配置更长超时。
        """
        return self._cfg["runtime"].get("generation_timeout", 3600)

    @property
    def language(self) -> str:
        """默认界面语言（"zh" 中文 / "en" English）。"""
        return self._cfg["runtime"].get("language", "zh")

    @property
    def web_port(self) -> int:
        """网页端监听端口（web/server.py 使用）。"""
        return int(self._cfg["runtime"].get("web_port", 8000))

    @property
    def web_idle_timeout(self) -> int:
        """网页端空闲自动停服时间（秒）。

        页面关闭后心跳停止，超过该时长守护线程自动关闭服务并释放端口。
        """
        return int(self._cfg["runtime"].get("web_idle_timeout", 60))

    @property
    def position_tolerance(self) -> float:
        """原子位置容差（分数坐标）。"""
        return self._cfg["defaults"]["position_tolerance"]

    @property
    def defaults(self) -> dict:
        """计算默认参数。"""
        return self._cfg["defaults"]


def get_config() -> Config:
    """获取全局配置实例

    """
    return Config()
