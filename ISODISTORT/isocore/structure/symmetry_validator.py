"""
对称性校验 - 验证结构与空间群的一致性

对应阶段一，步骤3：结构对称性校验
实现方式：⚖️ 基于 pymatgen 的 SpacegroupAnalyzer 判定空间群与对称性
"""
import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ..utils import get_config


class SymmetryValidator:
    """
    晶体结构对称性校验器

    功能：
    1. 确认结构的空间群（与 findsym 结果交叉验证）
    2. 检查原子位置容差
    3. 识别有序/无序占位
    4. 获取 Wyckoff 位置分配
    """

    def __init__(
        self,
        tolerance: float | None = None,
        *,
        angle_tolerance_degrees: float | None = None,
    ):
        cfg = get_config()
        # pymatgen forwards ``symprec`` to spglib as an absolute Cartesian
        # distance in the structure's lattice unit (Å here), not as a
        # fractional-coordinate residual.  Keep ``tolerance`` as a compatible
        # public argument, but make its physical meaning explicit internally.
        self.symprec_angstrom = (
            cfg.symmetry_cartesian_tolerance_angstrom
            if tolerance is None
            else float(tolerance)
        )
        self.angle_tolerance_degrees = (
            cfg.symmetry_angle_tolerance_degrees
            if angle_tolerance_degrees is None
            else float(angle_tolerance_degrees)
        )
        if not np.isfinite(self.symprec_angstrom) or self.symprec_angstrom <= 0:
            raise ValueError("symmetry tolerance must be a finite positive Å value")
        if (
            not np.isfinite(self.angle_tolerance_degrees)
            or self.angle_tolerance_degrees <= 0
        ):
            raise ValueError("angle tolerance must be a finite positive degree value")
        # Compatibility for callers/reporters that used ``validator.tolerance``.
        self.tolerance = self.symprec_angstrom

    @staticmethod
    def fractional_norm_upper_bound(
        structure: Structure, cartesian_tolerance_angstrom: float
    ) -> float:
        """Convert a Cartesian distance threshold to a safe fractional bound.

        For row-vector coordinates ``dr_cart = dr_frac @ A``, singular-value
        analysis gives ``||dr_frac|| <= ||dr_cart|| / sigma_min(A)``.  The
        result is dimensionless and depends on the full cell metric, so it
        remains meaningful for skewed as well as differently scaled cells.
        It is diagnostic only; spglib receives the Cartesian tolerance itself.
        """
        singular_values = np.linalg.svd(
            np.asarray(structure.lattice.matrix, dtype=float), compute_uv=False
        )
        smallest = float(np.min(singular_values))
        if not np.isfinite(smallest) or smallest <= 0:
            raise ValueError("structure lattice must be finite and nonsingular")
        return float(cartesian_tolerance_angstrom) / smallest

    def validate(self, structure: Structure) -> dict:
        """
        执行对称性校验

        Returns:
            dict: 包含空间群号、空间群符号、Wyckoff位点、是否有序等信息

        """
        sga = SpacegroupAnalyzer(
            structure,
            symprec=self.symprec_angstrom,
            angle_tolerance=self.angle_tolerance_degrees,
        )

        sg_number = sga.get_space_group_number()
        sg_symbol = sga.get_space_group_symbol()
        symm_structure = sga.get_symmetrized_structure()

        # Wyckoff 位置信息
        wyckoff_labels = symm_structure.wyckoff_symbols
        equivalent_indices = symm_structure.equivalent_indices

        wyckoff_sites = []
        for label, indices in zip(wyckoff_labels, equivalent_indices, strict=False):
            # label 形如 "4a"
            multiplicity = int(label[:-1])
            letter = label[-1]
            wyckoff_sites.append({
                "wyckoff_letter": letter,
                "multiplicity": multiplicity,
                "species": structure[indices[0]].species_string,
                "representative_index": indices[0],
                "equivalent_indices": indices,
            })

        # 占位检查：同时覆盖混占与部分/超占位
        has_disorder = any(
            (len(site.species) > 1)
            or (not np.isclose(sum(site.species.values()), 1.0, atol=1e-8))
            for site in structure
        )

        return {
            "space_group_number": sg_number,
            "space_group_symbol": sg_symbol,
            "wyckoff_sites": wyckoff_sites,
            "has_disorder": has_disorder,
            # Keep the old field for API compatibility, but add unambiguous
            # units and the lattice-dependent fractional-coordinate bound.
            "tolerance": self.symprec_angstrom,
            "symmetry_tolerance_angstrom": self.symprec_angstrom,
            "symmetry_angle_tolerance_degrees": self.angle_tolerance_degrees,
            "fractional_norm_tolerance_upper_bound": (
                self.fractional_norm_upper_bound(
                    structure, self.symprec_angstrom
                )
            ),
        }

    @staticmethod
    def get_centering(structure: Structure) -> str:
        """获取点阵中心类型 (P/I/F/A/B/C/R)

        """
        cfg = get_config()
        sga = SpacegroupAnalyzer(
            structure,
            symprec=cfg.symmetry_cartesian_tolerance_angstrom,
            angle_tolerance=cfg.symmetry_angle_tolerance_degrees,
        )
        sg_symbol = sga.get_space_group_symbol()
        # 简单提取第一个字母
        return sg_symbol[0] if sg_symbol else "P"
