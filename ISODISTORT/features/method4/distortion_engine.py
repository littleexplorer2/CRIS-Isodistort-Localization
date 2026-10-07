"""
畸变引擎 - 幅度缩放、多模式混合、生成畸变结构

对应阶段五，步骤9：畸变幅度缩放与多模式混合

生成流程（与官网 Distortion Page 语义对齐）：
1. 在母相原胞上应用模式位移：新坐标 = 原坐标 + 幅度 × 位移向量
2. 再按子群超胞基矢（3x3 矩阵或 [a,b,c]）扩胞

输入约定（见 README）：
- 本引擎接收原始分数位移列；As/Ap 与 normfactor 的换算由上层共享模式合同负责。
- 母胞单实列只能证明 Γ 或自共轭实平移字符；一般非 Γ 模式须消费完整子胞列。
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from pymatgen.core import Structure

from backend.utils import get_config
from backend.utils.lattice import (
    RationalVector,
    is_sublattice,
    rational_matrix,
    rational_vector,
    real_translation_character_sign,
    wavevector_translation_phases,
)
from features.input_cif import build_supercell, wrap_to_unit_cell
from features.method1.affine_embeddings import parent_affine_group
from features.method2.distortion_mapper import DistortionMapper

SupercellSpec = Sequence[int] | Sequence[Sequence[float]] | None


def _as_basis_matrix(supercell: SupercellSpec) -> np.ndarray:
    """把超胞规格规范化为 3x3 基矢矩阵。

    兼容两种输入（与 build_supercell 一致）：
    - [a, b, c] 整数列表：按对角矩阵补成 3x3；
    - 3x3 矩阵（行向量为母相格单位下的新基矢）：原样返回。
    这样下方坐标变换（coords @ basis、inv(basis)）始终作用于 3x3 矩阵，
    避免 1D 数组触发 np.linalg.inv 崩溃。
    """
    arr = np.asarray(supercell, dtype=float)
    if arr.ndim == 1:
        if arr.size != 3:
            raise ValueError(f"超胞必须是 [a,b,c] 或 3x3 矩阵，收到: {supercell}")
        return np.diag(arr)
    if arr.shape != (3, 3):
        raise ValueError(f"超胞必须是 [a,b,c] 或 3x3 矩阵，收到: {supercell}")
    return arr


class DistortionEngine:
    """
    畸变计算引擎

    功能：
    1. 单模式畸变：按指定幅度生成畸变结构
    2. 多模式混合：线性叠加多个不可约表示的畸变（generate_modes）
    3. 超胞构建：生成对应子群的超胞结构
    """

    def __init__(self, mapper: DistortionMapper | None = None) -> None:
        self.mapper = mapper or DistortionMapper()
        cfg = get_config()
        self.default_amplitude = cfg.defaults["default_amplitude"]

    def generate_single_mode(self, parent_structure: Structure,
                             mode_displacements: np.ndarray,
                             amplitude: float | None = None,
                             supercell: SupercellSpec = None,
                             k_vector: Sequence[float] | None = None,
                             opd_direction: Sequence[float] | None = None,
                             ) -> Structure:
        """
        生成单模式畸变后的结构

        Args:
            parent_structure: 母相结构
            mode_displacements: 模式位移向量 (N, 3)，母相分数坐标单位
            amplitude: 畸变幅度（默认取配置 default_amplitude）
            supercell: 超胞规格，可为 [a,b,c] 整数列表或 3x3 矩阵
                （子群基矢，母相格单位）
            k_vector: 母相倒格分数坐标。仅支持 Γ 或母相平移群上
                自共轭的实字符（±1）；一般非 Γ 模式必须提供完整子胞列。
            opd_direction: 旧接口保留；单个母胞实列不足以解释 OPD，
                提供此参数时明确报错，避免静默忽略。

        Returns:
            Structure: 畸变后的结构
        """
        amplitude = self.default_amplitude if amplitude is None else amplitude
        if mode_displacements.shape[0] != len(parent_structure):
            raise ValueError(
                f"位移向量长度 {mode_displacements.shape[0]} 与母相原子数 "
                f"{len(parent_structure)} 不一致"
            )
        if not np.all(np.isfinite(mode_displacements)):
            raise ValueError(
                "位移向量含 NaN/Inf，无法生成畸变结构（请检查模式计算是否出错）"
            )
        if not np.isfinite(amplitude):
            raise ValueError(f"振幅必须为有限数值，收到 {amplitude!r}")

        if opd_direction is not None:
            raise ValueError(
                "opd_direction cannot be interpreted from a single real parent-cell "
                "column; provide proven complete-supercell displacement columns"
            )
        return self.generate_modes(
            parent_structure, supercell,
            parent_displacements=amplitude * np.asarray(mode_displacements, dtype=float),
            k_vector=k_vector,
        )

    def generate_modes(self, parent_structure: Structure,
                       supercell: SupercellSpec = None,
                       parent_displacements: np.ndarray | None = None,
                       occupancy_patterns: list[tuple[np.ndarray, float]] | None = None,
                       k_vector: Sequence[float] | None = None,
                       ) -> Structure:
        """
        组合生成畸变结构：可同时施加原子位移与占据率（occupational）调制。

        Args:
            parent_structure: 母相结构
            supercell: 超胞规格（子群基矢 3x3 矩阵或 [a,b,c]）
            parent_displacements: (N_parent, 3) 位移向量（母相分数坐标单位），
                每个超胞副本按对应母相原子施加；None 表示无位移
            occupancy_patterns: [(pattern, amplitude), ...]，pattern 为
                长度 = 超胞原子数的 ±1 数组（+1 类保持全占据，-1 类占据率
                1-amplitude）；None 表示无占据率调制
            k_vector: 母相倒格分数坐标；仅支持 Γ/自共轭实平移字符，
                缺少完整复列/成对实列的一般非 Γ 模式明确报错。

        Returns:
            Structure: 畸变后的结构（占据率调制时位点为部分占据）
        """
        if parent_displacements is not None \
                and parent_displacements.shape[0] != len(parent_structure):
            raise ValueError(
                f"位移向量长度 {parent_displacements.shape[0]} 与母相原子数 "
                f"{len(parent_structure)} 不一致"
            )
        sc = build_supercell(parent_structure, supercell) if supercell is not None \
            else parent_structure.copy()

        coords = np.asarray(sc.frac_coords, dtype=float)
        if parent_displacements is not None:
            basis = _as_basis_matrix(supercell) if supercell is not None else np.eye(3)
            coords = coords + self._lift_real_column(
                parent_structure, sc, basis, parent_displacements, k_vector,
            )
            coords = wrap_to_unit_cell(coords)

        def _elem_symbol(spec) -> str:
            """位点成分 -> 元素符号（兼容 Element / Composition / str）。"""
            if isinstance(spec, str):
                return spec
            if hasattr(spec, "symbol") and not hasattr(spec, "elements"):
                return str(spec.symbol)  # Element
            els = spec.elements  # Composition
            return str(els[0].symbol) if len(els) == 1 else str(spec)

        species = [site.species for site in sc]
        if occupancy_patterns:
            n_sc = len(sc)
            for pattern, amp in occupancy_patterns:
                if len(pattern) != n_sc:
                    raise ValueError(
                        f"占据率模式长度 {len(pattern)} 与超胞原子数 {n_sc} 不一致"
                    )
                amp_capped = min(max(float(amp), 0.0), 1.0)
                new_species = []
                for j, site_species in enumerate(species):
                    sym = _elem_symbol(site_species)
                    # 占据率语义：+1 类全占据（1.0）；-1 类占据率 = 1 - amplitude；
                    # pattern == 0 的位点未被该模式调制，应保持全占据，不可降低。
                    if pattern[j] < 0:
                        occ = max(1.0 - amp_capped, 0.01)
                        if occ >= 0.999:
                            new_species.append(sym)
                        else:
                            new_species.append({sym: occ})
                    else:
                        new_species.append(sym)
                species = new_species

        return Structure(sc.lattice, species, coords, coords_are_cartesian=False)

    @staticmethod
    def _map_parent_indices(parent: Structure, sc: Structure,
                            supercell) -> np.ndarray:
        """把超胞原子映射回母相原子索引（几何最近邻 + 同物种，鲁棒于排序）。"""
        basis = _as_basis_matrix(supercell)
        n, m = len(parent), len(sc)
        pc = np.asarray(parent.frac_coords, dtype=float)
        pj = np.asarray(sc.frac_coords, dtype=float) @ basis
        idx = np.zeros(m, dtype=int)
        for j in range(m):
            best_dist, best_i = float("inf"), -1
            for i in range(n):
                if sc[j].species_string != parent[i].species_string:
                    continue
                d = pj[j] - pc[i]
                d -= np.round(d)
                dist = float(np.linalg.norm(d))
                if dist < best_dist:
                    best_dist, best_i = dist, i
            idx[j] = best_i
        return idx

    @staticmethod
    def _real_parent_wavevector(
        parent: Structure,
        basis: np.ndarray,
        displacement: np.ndarray,
        k_vector: Sequence[float] | None,
    ) -> RationalVector | None:
        """Prove that a single real parent-cell column has enough phase data.

        ``k=-k+G`` means a real translation character only when G belongs to
        the *primitive* reciprocal lattice. Integer conventional components
        alone are insufficient for centered parents. The affine-group owner
        supplies the detected translation lattice in the actual input axes.
        """
        parent_group = parent_affine_group(parent)
        child_basis = rational_matrix(basis)
        if not is_sublattice(child_basis, parent_group.lattice):
            raise ValueError("requested child cell is not a sublattice of the parent translation lattice")
        if k_vector is None:
            # An undeclared parent-cell field need not be Gamma on the full
            # primitive parent lattice. It must still be periodic under the
            # child generators, particularly when a centered cell is reduced.
            k = rational_vector((0, 0, 0))
            phases = ()
            covariance_generators = child_basis
        else:
            k = rational_vector(k_vector)
            phases = wavevector_translation_phases(k, parent_group.lattice)
            if any((2 * value).denominator != 1 for value in phases):
                raise ValueError(
                    "Non-self-conjugate wavevector cannot be lifted from a single real "
                    "parent-cell column: paired real/complex Bloch columns, k-star/OPD "
                    "and phase-origin evidence are missing; use complete-supercell columns"
                )
            if any(value.denominator != 1 for value in
                   wavevector_translation_phases(k, child_basis)):
                raise ValueError("wavevector is not periodic in the requested child cell")
            covariance_generators = parent_group.lattice
        # A centered conventional cell already contains several primitive
        # translates. Their supplied vectors must realize the same character
        # as the newly generated cell copies, rather than assigning a hidden
        # complex phase or zeroing a missing component.
        cartesian = displacement @ np.asarray(parent.lattice.matrix, dtype=float)
        scale = max(1.0, float(np.max(np.abs(cartesian))))
        vector_tolerance = 256.0 * np.finfo(float).eps * scale
        position_tolerance = get_config().symmetry_cartesian_tolerance_angstrom
        parent_coords = np.asarray(parent.frac_coords, dtype=float)
        for translation in covariance_generators:
            translated = parent_coords + np.asarray(translation, dtype=float)
            distances = parent.lattice.get_all_distances(translated, parent_coords)
            for source, site in enumerate(parent):
                candidates = [
                    target for target, candidate in enumerate(parent)
                    if candidate.species == site.species
                    and distances[source, target] <= position_tolerance
                ]
                if len(candidates) != 1:
                    raise ValueError("parent translation cannot be mapped uniquely to equivalent sites")
                target = candidates[0]
                image = np.round(translated[source] - parent_coords[target]).astype(int)
                effective_translation = tuple(
                    component - int(shift)
                    for component, shift in zip(translation, image, strict=True)
                )
                sign = real_translation_character_sign(k, effective_translation)
                if not np.allclose(
                    cartesian[target], sign * cartesian[source],
                    atol=vector_tolerance, rtol=0.0,
                ):
                    if k_vector is None:
                        raise ValueError(
                            "parent-cell displacement column is not periodic "
                            "in the requested child cell"
                        )
                    raise ValueError(
                        "parent-cell displacement column does not prove the "
                        "self-conjugate translation character on equivalent sites"
                    )
        return None if all(value.denominator == 1 for value in phases) else k

    def _lift_real_column(
        self,
        parent: Structure,
        child: Structure,
        basis: np.ndarray,
        values: np.ndarray,
        k_vector: Sequence[float] | None,
    ) -> np.ndarray:
        displacement = np.asarray(values, dtype=float)
        if displacement.shape != (len(parent), 3) or not np.all(np.isfinite(displacement)):
            raise ValueError("parent-cell displacement column must be finite with shape (N, 3)")
        k = self._real_parent_wavevector(parent, basis, displacement, k_vector)
        indices = self._map_parent_indices(parent, child, basis)
        if np.any(indices < 0):
            raise ValueError("child atoms cannot be mapped to the parent-cell column")
        lifted = displacement[indices] @ np.linalg.inv(basis)
        if k is not None:
            parent_frac = np.asarray(child.frac_coords, dtype=float) @ basis
            parent_coords = np.asarray(parent.frac_coords, dtype=float)
            translations = np.round(parent_frac - parent_coords[indices]).astype(int)
            for row, translation in enumerate(translations):
                lifted[row] *= real_translation_character_sign(k, translation)
        return lifted

    def lift_mode_displacements(self, parent_structure: Structure,
                                supercell: SupercellSpec,
                                mode_displacements: dict,
                                k_vector: Sequence[float] | None = None,
                                ) -> tuple[Structure, dict]:
        """把母相坐标系的模式位移提升到超胞坐标系（Method 4 超胞分解用）。

        与 generate_modes 的位移施加逻辑严格互逆：
        - 超胞原子 j 对应母相原子 idx[j]（几何最近邻 + 同物种）；
        - 位移由母相分数坐标换算到超胞分数坐标：Δf_sc = Δf_parent @ inv(B)；
        - 自共轭 k 的副本使用经母相原始平移群验证的实字符 ±1；
          一般非 Γ 模式缺少完整复列/成对实列时 fail-closed。

        Args:
            parent_structure: 母相结构
            supercell: 超胞规格（子群基矢 3x3 矩阵或 [a,b,c]）
            mode_displacements: {label: (N_parent, 3)} 母相分数坐标位移
            k_vector: k 点坐标（母相倒格分数单位）；None/Γ 点不调制

        Returns:
            (超胞母相结构, {label: (N_sc, 3) 超胞分数坐标位移})
        """
        basis = _as_basis_matrix(supercell)
        sc = build_supercell(parent_structure, basis)
        lifted: dict = {}
        for label, disp in mode_displacements.items():
            lifted[label] = self._lift_real_column(
                parent_structure, sc, basis, disp, k_vector,
            )
        return sc, lifted
