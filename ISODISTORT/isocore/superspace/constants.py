"""(3+d) 超空间内核常量。

命名与官网 / 本项目对齐：
- ``nmod`` = 独立非公度调制数 = 附加维度 d
- ``EPS`` 与 ``config/settings.yaml`` 的 ``defaults.eps`` / ``lattice_tolerance`` 相同
- 空间群符号用 Hermann–Mauguin + Schoenflies（IT 旧简写）
- 超空间群取位固定为 standard (IT-C)
"""
from __future__ import annotations

# 与 defaults.lattice_tolerance / defaults.eps 对齐（0.00001）
EPS = 1e-5
# 官网 nmod 常用范围 (3+1)/(3+2)/(3+3)
MAX_NMOD = 3
PHYSICAL_DIM = 3
# JSON/YAML 对象类型标记（兼容 nmod=1 的 (3+1) 文件）
SERIAL_KIND = "superspace"
SERIAL_KIND_3P1 = "superspace_3p1"


def eps() -> float:
    """运行时 EPS：优先读配置，失败则回退模块常量。"""
    try:
        from ..utils import get_config  # noqa: PLC0415

        return float(get_config().eps)
    except (OSError, KeyError, TypeError, ValueError):
        return EPS


def max_nmod() -> int:
    try:
        from ..utils import get_config  # noqa: PLC0415

        return int(get_config().max_nmod)
    except (OSError, KeyError, TypeError, ValueError):
        return MAX_NMOD
