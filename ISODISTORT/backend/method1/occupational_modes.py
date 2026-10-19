"""Compatibility interface for occupational modes, currently fail-closed.

Sorting sites into two equal halves does not construct a scalar irrep, fix
composition, or establish the requested embedded subgroup. A space-group
number match cannot establish k/IR/OPD identity. Until those contracts exist,
this interface must not generate a candidate pattern.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from pymatgen.core import Structure

from backend.utils import get_config
from backend.models.phase_path import require_supported_physical_types


@dataclass
class OccupationalMode:
    """占据率畸变模式（作用于子群超胞中的具体原子）。"""

    label: str                            # 如 occ-Al-d（物种-位点）
    species: str                          # 作用的元素符号
    wyckoff_letter: str                   # 母相 Wyckoff 位置字母
    basis_vectors: list[list[float]] = field(default_factory=list)  # 子群超胞基矢
    pattern: np.ndarray = None            # 每个超胞原子 +1/-1（对齐超胞原子顺序）
    parent_index: np.ndarray = None       # 每个超胞原子对应的母相原子索引
    validated: bool = False               # spglib 校验是否命中子群空间群
    note: str = ""                        # 附加说明（如“近似模式”）
    k_point_label: str = ""
    irrep_label: str = ""


class OccupationalModeGenerator:
    """Retain the legacy API without producing unproved scalar modes."""

    def __init__(
        self,
        tolerance: float | None = None,
        *,
        angle_tolerance_degrees: float | None = None,
    ) -> None:
        cfg = get_config()
        self.tolerance = (
            cfg.symmetry_cartesian_tolerance_angstrom
            if tolerance is None
            else float(tolerance)
        )
        self.angle_tolerance_degrees = (
            cfg.symmetry_angle_tolerance_degrees
            if angle_tolerance_degrees is None
            else float(angle_tolerance_degrees)
        )
        if not np.isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError("symmetry tolerance must be a finite positive Å value")
        if (
            not np.isfinite(self.angle_tolerance_degrees)
            or self.angle_tolerance_degrees <= 0
        ):
            raise ValueError("angle tolerance must be a finite positive degree value")

    def generate(self, parent: Structure, wyckoff_sites: list[dict],
                 subgroup, species_scope: set[str]) -> list[OccupationalMode]:
        """Require a proved scalar/composition representation before generation."""
        require_supported_physical_types(["occupational"])
        return []  # unreachable until the representation contract is implemented
