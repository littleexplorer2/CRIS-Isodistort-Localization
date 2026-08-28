"""空间群 / k 点 / irrep 元数据（CDML 显示约定等）。

禁止在此放置按 (irrep, OPD) 或具体算例硬编码的“官网正确答案”表；
对齐应通过底层算法修正，而非特例数据库。
"""
from .irreps_cdml import lookup_irrep_kovalev
from .kpoints_official import KPOINT_OFFICIAL

__all__ = [
    "KPOINT_OFFICIAL",
    "lookup_irrep_kovalev",
]
