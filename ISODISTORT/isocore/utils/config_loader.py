"""
配置加载器 - 读取 settings.yaml 并设置 ISODATA 环境变量
"""
import math
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

    def _nonnegative_runtime_int(self, key: str, default: int = 0) -> int:
        """Return a non-negative integer runtime limit.

        Method 3 cost limits use ``0`` as the explicit unlimited sentinel.
        Reject invalid or negative YAML values at the configuration boundary so
        callers cannot accidentally reinterpret a bad limit as "unlimited".
        """
        value = self._cfg["runtime"].get(key, default)
        try:
            parsed = int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid configuration runtime.{key}={value!r}: "
                "expected a non-negative integer (0 = unlimited)"
            ) from exc
        if isinstance(value, float) and not value.is_integer():
            raise ValueError(
                f"Invalid configuration runtime.{key}={value!r}: "
                "expected a non-negative integer (0 = unlimited)"
            )
        if parsed < 0:
            raise ValueError(
                f"Invalid configuration runtime.{key}={parsed}: "
                "expected a non-negative integer (0 = unlimited)"
            )
        return parsed

    @property
    def method3_max_parametric_values(self) -> int:
        """Maximum exact parametric-k searches in Method 3; 0 is unlimited."""
        return self._nonnegative_runtime_int("method3_max_parametric_values")

    @property
    def method3_max_backend_queries(self) -> int:
        """Maximum Method 3 ``list_irreps`` + ``list_subgroups`` calls."""
        return self._nonnegative_runtime_int("method3_max_backend_queries")

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

    def _positive_default_float(self, key: str, fallback: str | None = None) -> float:
        """Read one finite positive numerical tolerance from ``defaults``.

        Coordinate tolerances are scientific input, not boolean switches.  A
        zero, negative or non-finite value would silently change a symmetry
        problem into an ill-defined exact/always-failing comparison, so reject
        it at the configuration boundary.
        """
        defaults = self._cfg["defaults"]
        if key in defaults:
            raw = defaults[key]
        elif fallback is not None and fallback in defaults:
            raw = defaults[fallback]
        else:
            raise KeyError(f"missing defaults.{key}")
        try:
            value = float(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid configuration defaults.{key}={raw!r}: "
                "expected a finite positive number"
            ) from exc
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(
                f"Invalid configuration defaults.{key}={raw!r}: "
                "expected a finite positive number"
            )
        return value

    @property
    def symmetry_cartesian_tolerance_angstrom(self) -> float:
        """spglib/pymatgen ``symprec`` for uploaded structures, in Å.

        The fallback keeps old user configuration files readable.  Despite
        its historical name, ``position_tolerance`` was always passed to
        ``symprec`` and therefore also had Cartesian-length semantics.
        """
        return self._positive_default_float(
            "symmetry_cartesian_tolerance_angstrom", "position_tolerance"
        )

    @property
    def symmetry_angle_tolerance_degrees(self) -> float:
        """pymatgen/spglib lattice-angle tolerance, in degrees."""
        return self._positive_default_float("symmetry_angle_tolerance_degrees")

    @property
    def affine_exact_cartesian_tolerance_angstrom(self) -> float:
        """spglib length tolerance for exact synthetic affine data, in Å."""
        return self._positive_default_float(
            "affine_exact_cartesian_tolerance_angstrom"
        )

    @property
    def fractional_coordinate_tolerance(self) -> float:
        """Dimensionless residual for float-to-rational coordinate recovery."""
        return self._positive_default_float("fractional_coordinate_tolerance")

    @property
    def position_tolerance(self) -> float:
        """Backward-compatible alias for the Cartesian Å symmetry tolerance."""
        return self.symmetry_cartesian_tolerance_angstrom

    @property
    def lattice_tolerance(self) -> float:
        """Dimensionless lattice/integer-matrix residual tolerance."""
        return self._positive_default_float("lattice_tolerance")

    @property
    def eps(self) -> float:
        """全局浮点比较容差 EPS。

        与 ``defaults.eps`` 对齐；缺省时回退到 ``lattice_tolerance``（0.00001）。
        """
        defaults = self._cfg["defaults"]
        if "eps" not in defaults:
            return self.lattice_tolerance
        return self._positive_default_float("eps")

    @property
    def defaults(self) -> dict:
        """计算默认参数。"""
        return self._cfg["defaults"]


def get_config() -> Config:
    """获取全局配置实例

    """
    return Config()
