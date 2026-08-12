"""distortion 包 - 畸变业务层（项目核心）"""
from .phase_path import (
    PhasePath,
    DISTORTION_TYPES,
    DEFAULT_DISTORTION_TYPES,
    normalize_distortion_types,
)
from .distortion_mapper import DistortionMapper
from .distortion_engine import DistortionEngine
from .domain_generator import DomainGenerator
from .search_methods import (
    IsoSearchEngine,
    Method1Query,
    Method1ResultItem,
    Method2Query,
    Method2Result,
    Method3Query,
    Method3ResultItem,
    Method4Query,
    Method4Result,
)

__all__ = [
    "PhasePath",
    "DISTORTION_TYPES",
    "DEFAULT_DISTORTION_TYPES",
    "normalize_distortion_types",
    "DistortionMapper",
    "DistortionEngine",
    "DomainGenerator",
    "IsoSearchEngine",
    "Method1Query",
    "Method1ResultItem",
    "Method2Query",
    "Method2Result",
    "Method3Query",
    "Method3ResultItem",
    "Method4Query",
    "Method4Result",
]
