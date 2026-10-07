"""
畸变模式映射器 - 将 iso（DISPLAY BUSH）输出的模式基矢映射到每个原子的位移向量

对应阶段四，步骤8：畸变基矢到原子坐标的映射

映射规则（基于实测 DISPLAY BUSH 输出）：
- BUSH 每个位移列是子群固定子空间的一条独立基矢。
  ``IsoWrapper`` 已按列拆分，因此这里的每个 ``BushMode`` 恰含一个向量。
- 若某位点只有一个代表原子行（常见于 Gamma 点均匀模式），则将该位移
  均匀作用于该位点的全部等效原子。
- 若某位点有多个代表原子行，每个代表对应模式分裂出的一个子轨道：
  - 先按**结构真实坐标**解析代表点中的自由参数（如 ``(0,0,z)`` 中的 z
    取该位点原子的实际坐标，否则两个相反代表点会被解析成同一点）；
  - 再按“模母相格点等价”（格点平移副本取同一位移，如带心平移）把位点
    原子分配到各代表；无法匹配的原子按周期最近邻回退。

下游 Method 4 会把这里保留的内部列尺度与官网 complete-mode 的
child-fractional max-component-one 展示尺度显式区分，并按 primitive-cell
Cartesian RSS 定义换算 normfactor、As 与 Ap。
"""

import re
from collections.abc import Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
from typing import Literal

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.groups import SpaceGroup

from backend.models.iso_mode_models import (
    MicroscopicColumnProvenance,
    ModeIdentity,
    microscopic_atom_order_id,
    validate_microscopic_mode_source,
    validate_signed_cartesian_columns,
)
from backend.utils.lattice import as_fraction, centering_primitive_matrix
from backend.utils.opd_format import (
    k_star_fraction_vectors,
    resolve_active_k_star_coordinates,
)
from backend.wrappers import (
    BushMode,
    DistortionMode,
    MicroscopicDomainExtensionEvidence,
    SubgroupInfo,
)
from backend.wrappers.iso_wrapper import (
    microscopic_bush_column_digest,
    microscopic_bush_row_digest,
    microscopic_bush_source_float_column_digest,
    microscopic_subgroup_translation_lattice,
    validate_microscopic_distortion_mode_source,
    validate_microscopic_domain_extension_bindings,
)
from features.input_cif import build_supercell

MappedColumnRelationStatus = Literal[
    "canonical",
    "unverified",
    "scalar_verified",
]


@dataclass(frozen=True)
class MappedMicroscopicColumn:
    """One mapped column plus immutable ISO source-column evidence."""

    mode_id: str
    amplitude_key: str
    frame_id: str
    atom_order_id: str
    displacements: np.ndarray
    provenance: MicroscopicColumnProvenance
    relation_status: MappedColumnRelationStatus = "unverified"
    current_over_canonical_signed_scale: float | None = None
    cartesian_residual_angstrom: float | None = None
    atom_ids: tuple[str, ...] = ()
    mode_identity: ModeIdentity | None = None
    domain_extension_evidence: MicroscopicDomainExtensionEvidence | None = None

    def __post_init__(self) -> None:
        if not self.mode_id.strip() or not self.amplitude_key.strip():
            raise ValueError("mapped microscopic columns require stable mode and amplitude IDs")
        if not self.frame_id.strip():
            raise ValueError("mapped microscopic column frame ID cannot be empty")
        if not self.atom_order_id.strip():
            raise ValueError("mapped microscopic column atom-order ID cannot be empty")
        atom_ids = tuple(str(value).strip() for value in self.atom_ids)
        if not atom_ids or any(not value for value in atom_ids):
            raise ValueError("mapped microscopic column atom IDs cannot be empty")
        if len(set(atom_ids)) != len(atom_ids):
            raise ValueError("mapped microscopic column atom IDs must be unique")
        if self.atom_order_id != microscopic_atom_order_id(atom_ids):
            raise ValueError("mapped microscopic atom-order ID does not match its atom IDs")
        identity = self.mode_identity
        if identity is None:
            raise ValueError("mapped microscopic column requires a ModeIdentity")
        validate_microscopic_mode_source(
            identity,
            self.provenance,
            mode_id=self.mode_id,
            require_resolved_provenance=True,
        )
        if self.provenance.source_frame_id != self.frame_id:
            raise ValueError("mapped microscopic provenance frame mismatch")
        if self.provenance.source_atom_order_id != self.atom_order_id:
            raise ValueError("mapped microscopic provenance atom-order mismatch")
        extension = self.domain_extension_evidence
        if extension is not None and (
            extension.canonical_source_column_token
            != self.provenance.exact_source_token
        ):
            raise ValueError("mapped domain extension is bound to another source column")
        if extension is not None and (
            extension.canonical_source_order_key != self.provenance.source_order_key
        ):
            raise ValueError("mapped domain extension is bound to another source index")
        array = np.asarray(self.displacements, dtype=float)
        if array.ndim != 2 or array.shape[1] != 3 or not np.all(np.isfinite(array)):
            raise ValueError("mapped microscopic displacement must be finite (n_atoms, 3)")
        if array.shape[0] != len(atom_ids):
            raise ValueError("mapped microscopic atom IDs must match displacement rows")
        immutable = np.array(array, dtype=float, copy=True)
        immutable.setflags(write=False)
        object.__setattr__(self, "atom_ids", atom_ids)
        object.__setattr__(self, "displacements", immutable)
        scale = self.current_over_canonical_signed_scale
        residual = self.cartesian_residual_angstrom
        if self.relation_status == "unverified":
            if scale is not None or residual is not None:
                raise ValueError("unverified mapped columns cannot expose scalar diagnostics")
        elif self.relation_status == "canonical":
            if scale != 1.0 or residual != 0.0:
                raise ValueError("canonical mapped columns require unit scale and zero residual")
        elif self.relation_status == "scalar_verified":
            if scale is None or not np.isfinite(scale) or scale == 0.0:
                raise ValueError("verified mapped columns require one finite nonzero scale")
            if residual is None or not np.isfinite(residual) or residual < 0:
                raise ValueError("verified mapped columns require a finite residual")
        else:
            raise ValueError("unknown mapped microscopic column relation status")


@dataclass(frozen=True)
class _ExactTranslationLattice:
    basis: tuple[tuple[Fraction, Fraction, Fraction], ...]
    inverse: np.ndarray
    metric: Lattice


def _make_exact_translation_lattice(
    basis: Sequence[Sequence[Fraction | int]],
    parent_lattice: np.ndarray,
) -> _ExactTranslationLattice:
    exact_basis = tuple(
        tuple(Fraction(value) for value in row) for row in basis
    )
    float_basis = np.asarray(exact_basis, dtype=float)
    if float_basis.shape != (3, 3) or abs(float(np.linalg.det(float_basis))) < 1e-12:
        raise ValueError("exact translation lattice must be nonsingular 3x3")
    lattice = np.asarray(parent_lattice, dtype=float)
    if lattice.shape != (3, 3) or not np.all(np.isfinite(lattice)):
        raise ValueError("parent lattice must be finite 3x3")
    return _ExactTranslationLattice(
        basis=exact_basis,
        inverse=np.linalg.inv(float_basis),
        metric=Lattice(float_basis @ lattice),
    )


def _snap_to_exact_translation(
    parent_delta: np.ndarray,
    lattice: _ExactTranslationLattice,
    cartesian_tolerance: float,
) -> tuple[Fraction, Fraction, Fraction] | None:
    """Return the nearest proven lattice translation, or ``None``.

    Floating coordinates are used only to choose and tolerance-check the
    nearest lattice image in the physical metric.  The returned translation
    is rebuilt from the integer image and exact rational basis, so downstream
    phase authorization never rounds a floating ``k dot t`` value.
    """

    delta = np.asarray(parent_delta, dtype=float)
    if delta.shape != (3,) or not np.all(np.isfinite(delta)):
        return None
    fractional = delta @ lattice.inverse
    distance, image = lattice.metric.get_distance_and_image(
        fractional,
        np.zeros(3),
    )
    if float(distance) > cartesian_tolerance:
        return None
    coefficients = tuple(int(value) for value in image)
    return tuple(
        sum(
            Fraction(coefficients[row]) * lattice.basis[row][column]
            for row in range(3)
        )
        for column in range(3)
    )


def _has_exact_unit_phase(
    k_arms: Sequence[Sequence[Fraction]],
    translation: Sequence[Fraction],
) -> bool:
    return all(
        sum(
            Fraction(arm[index]) * Fraction(translation[index])
            for index in range(3)
        ).denominator == 1
        for arm in k_arms
    )


def verify_mapped_microscopic_columns(
    canonical: Sequence[MappedMicroscopicColumn],
    current: Sequence[MappedMicroscopicColumn],
    reference_lattice: Sequence[Sequence[float]],
    *,
    rtol: float = 1e-8,
    atol_angstrom: float = 1e-10,
) -> tuple[MappedMicroscopicColumn, ...]:
    """Bind current columns to canonical sources and reject column mixing."""

    if len(canonical) == 0 or len(current) == 0:
        raise ValueError("mapped microscopic column sets cannot be empty")
    canonical_by_id = {column.mode_id: column for column in canonical}
    current_by_id = {column.mode_id: column for column in current}
    if len(canonical_by_id) != len(canonical):
        raise ValueError("canonical mapped microscopic mode IDs are not unique")
    if len(current_by_id) != len(current):
        raise ValueError("current mapped microscopic mode IDs are not unique")
    if set(canonical_by_id) != set(current_by_id):
        raise ValueError("mapped microscopic source columns are missing or unexpected")
    if len({column.frame_id for column in canonical}) != 1:
        raise ValueError("canonical mapped microscopic columns use multiple frames")
    if len({column.atom_order_id for column in canonical}) != 1:
        raise ValueError("canonical mapped microscopic columns use multiple atom orders")
    if len({column.atom_ids for column in canonical}) != 1:
        raise ValueError("canonical mapped microscopic columns use inconsistent atom ID orders")

    ordered_canonical = list(canonical)
    ordered_current = [current_by_id[column.mode_id] for column in ordered_canonical]
    for source, candidate in zip(ordered_canonical, ordered_current, strict=True):
        if source.relation_status != "canonical":
            raise ValueError("canonical mapped source column is not marked canonical")
        if source.mode_identity != candidate.mode_identity:
            raise ValueError("mapped microscopic mode identity mismatch")
        if source.amplitude_key != candidate.amplitude_key:
            raise ValueError("mapped microscopic amplitude-key mismatch")
        if source.frame_id != candidate.frame_id:
            raise ValueError("mapped microscopic mode frame mismatch")
        if source.atom_order_id != candidate.atom_order_id:
            raise ValueError("mapped microscopic atom-order mismatch")
        if source.atom_ids != candidate.atom_ids:
            raise ValueError("mapped microscopic atom ID order mismatch")
        if source.provenance != candidate.provenance:
            raise ValueError("mapped microscopic exact source provenance mismatch")
        if source.domain_extension_evidence != candidate.domain_extension_evidence:
            raise ValueError("mapped microscopic domain-extension evidence mismatch")

    alignment = validate_signed_cartesian_columns(
        [column.displacements for column in ordered_canonical],
        [column.displacements for column in ordered_current],
        np.asarray(reference_lattice, dtype=float),
        rtol=rtol,
        atol=atol_angstrom,
    )
    if not alignment.matched:
        raise ValueError(
            "mapped microscopic columns are not uniquely signed-scaled sources: "
            + alignment.reason
        )
    if alignment.canonical_to_current != tuple(range(len(ordered_canonical))):
        raise ValueError("mapped microscopic column data are attached to the wrong source IDs")
    return tuple(
        replace(
            candidate,
            relation_status="scalar_verified",
            current_over_canonical_signed_scale=scale,
            cartesian_residual_angstrom=residual,
        )
        for candidate, scale, residual in zip(
            ordered_current,
            alignment.current_over_canonical_scales,
            alignment.cartesian_residuals,
            strict=True,
        )
    )


class DistortionMapper:
    """
    畸变模式 -> 原子位移 映射器

    负责将 iso 计算出的抽象模式基矢，转化为每个原子的实际位移向量。
    """

    def __init__(self) -> None:
        self._prim_matrix: np.ndarray | None = None
        self._conv_matrix: np.ndarray | None = None

    def map_modes_to_atoms(self, structure: Structure,
                           wyckoff_sites: list[dict],
                           modes: list[DistortionMode]) -> dict:
        """
        将畸变模式映射到结构中每个原子

        Args:
            structure: 母相晶体结构
            wyckoff_sites: 结构的 Wyckoff 位置分组信息（来自 SymmetryValidator）
            modes: iso 计算出的畸变模式列表（含 bush_modes）

        Returns:
            dict: 每个模式对应的原子位移向量
                {
                    "irrep_label": {
                        "mode": DistortionMode,
                        "displacements": np.ndarray (N_atoms x 3),
                    }
                }
        """
        n_atoms = len(structure)
        self._lattice_setup(structure)

        # 按 (字母, 物种) 分组：同一字母可能被多个物种共用（如 Fe/O 同在 k）
        groups: list[tuple[str, list[int], dict]] = []
        for site_info in wyckoff_sites:
            groups.append((site_info["wyckoff_letter"],
                           list(site_info["equivalent_indices"]),
                           site_info))

        result: dict[str, dict] = {}
        for mode in modes:
            # 按位点字母分组 BUSH 行；多位点时拆成独立幅度键（对齐官网
            # 每个 [Site:letter:dsp]… 各算一个 displacive mode）。
            by_letter: dict[str, list[BushMode]] = {}
            for bush in mode.bush_modes:
                by_letter.setdefault(bush.wyckoff_letter, []).append(bush)

            letter_disps: dict[str, np.ndarray] = {}
            target_orbit = str(getattr(mode, "wyckoff_orbit_id", "") or "")
            for letter, indices, site_info in groups:
                if target_orbit and str(site_info.get("orbit_id") or "") != target_orbit:
                    continue
                bushes = by_letter.get(letter, [])
                if not bushes or not indices:
                    continue
                displacements = np.zeros((n_atoms, 3))

                # 解析每个代表点（含自由参数 -> 结构真实坐标）。
                # 列在 wrapper 层已拆分；如果再收到多列数据，说明上游
                # 违反模式身份约定，必须显式报错而不是静默相加。
                reps: list[tuple[np.ndarray, np.ndarray]] = []
                for bush in bushes:
                    if not bush.displacements:
                        continue
                    if len(bush.displacements) != 1:
                        raise ValueError(
                            "DISPLAY BUSH displacement columns must be split "
                            "into independent DistortionMode objects before mapping"
                        )
                    point = self._resolve_rep_point(
                        bush, structure, indices, site_info=site_info,
                    )
                    vec = np.asarray(bush.displacements[0], dtype=float)
                    reps.append((point, vec))
                if not reps:
                    continue

                # 只保留能“模格点”匹配到本组原子的代表（同一字母不同物种时，
                # 彼此的符号点取位约定可能不同，无法匹配的物种保持不动，
                # 避免错误回退破坏对称性）
                matched_reps = [
                    (p, v) for p, v in reps
                    if any(self._lattice_equiv(
                        np.asarray(structure[i].frac_coords), p) for i in indices)
                ]
                if not matched_reps:
                    continue  # 本组原子保持不动
                if len(matched_reps) == 1:
                    # 均匀模式：同一位移作用于该位点的全部等效原子
                    disp = matched_reps[0][1]
                    for idx in indices:
                        displacements[idx] = disp
                else:
                    # 多代表：优先按“周期等价”（模整数平移）把代表点分配给
                    # 原子。非 Γ k 点模式中，BUSH 会为带心副本给出相位反号的
                    # 独立代表（如 Im-3m H4-：`(0,0,0)→(1,0,1)` 与
                    # `(-1/2,1/2,1/2)→(-1,0,-1)`），这两个代表相差原始格点
                    # 平移、互为“模格点等价”，若按旧逻辑“首个代表胜出”会把
                    # 全部原子分给同一代表、丢失 k 点相位（畸变退化为刚性平移）。
                    # 代表点与原子分数坐标 mod 1 重合时即视为同一原子，直接
                    # 采用该代表的位移；无法周期匹配时回退“模母相格点等价”
                    # （子轨道/带心副本同位移，Γ 点均匀语义），最后退周期最近邻。
                    for idx in indices:
                        coord = np.asarray(structure[idx].frac_coords)
                        assigned = False
                        for point, disp in matched_reps:
                            if self._periodic_equiv(coord, point):
                                displacements[idx] = disp
                                assigned = True
                                break
                        if not assigned:
                            for point, disp in matched_reps:
                                if self._lattice_equiv(coord, point):
                                    displacements[idx] = disp
                                    assigned = True
                                    break
                        if not assigned:
                            # 回退：周期最近邻（仅本组代表内）
                            best_disp = matched_reps[0][1]
                            best_dist = float("inf")
                            for point, disp in matched_reps:
                                d = self._periodic_distance(coord, point)
                                if d < best_dist:
                                    best_dist = d
                                    best_disp = disp
                            displacements[idx] = best_disp

                max_norm = float(np.max(np.linalg.norm(displacements, axis=1)))
                if max_norm > 1e-12:
                    letter_disps[letter] = displacements / max_norm

            if not letter_disps:
                continue
            unique = str(getattr(mode, "amplitude_key", "") or "")
            if unique:
                letter, displacements = next(iter(letter_disps.items()))
                result[unique] = {
                    "mode": mode,
                    "displacements": displacements,
                    "wyckoff_letter": letter,
                }
                continue
            if len(letter_disps) == 1:
                letter, displacements = next(iter(letter_disps.items()))
                result[mode.irrep_label] = {
                    "mode": mode,
                    "displacements": displacements,
                    "wyckoff_letter": letter,
                }
            else:
                for letter, displacements in letter_disps.items():
                    key = f"{mode.irrep_label}__{letter}"
                    result[key] = {
                        "mode": mode,
                        "displacements": displacements,
                        "wyckoff_letter": letter,
                    }

        return result

    def map_bush_modes_to_supercell(
        self,
        structure: Structure,
        wyckoff_sites: list[dict],
        modes: list[DistortionMode],
        basis_vectors: list[list[float]],
        *,
        cartesian_tolerance: float = 1e-4,
        subgroup_operations: Sequence[object] | None = None,
        subgroup_translation_lattice: Sequence[Sequence[float]] | None = None,
        subgroup_context: SubgroupInfo | None = None,
    ) -> dict[str, np.ndarray]:
        """Map complete DISPLAY BUSH basis vectors directly to the child cell.

        BUSH continuation rows already contain the relative signs/phases on
        every representative atom of the selected isotropy subgroup.  For a
        child basis ``B`` (rows in parent fractional coordinates), a child
        atom at parent coordinate ``x`` belongs to a BUSH representative
        ``p`` precisely when ``x-p`` belongs to the subgroup's *primitive*
        translation lattice ``T_H``.  Its arrow is transformed to the printed
        conventional child basis as ``du_child = du_parent @ inv(B)``.  The
        distinction matters for A/B/C/I/F/R centered child groups.

        This is deliberately different from reconstructing every mode with a
        cosine of the *primary* k vector: a single BUSH table can also contain
        secondary irreps at GM/X/N/etc., each with its own phase pattern.
        Missing or conflicting representatives are treated as invalid output
        rather than silently replaced by a guessed Bloch phase.
        """
        if not modes:
            return {}
        basis = np.asarray(basis_vectors, dtype=float)
        if basis.shape != (3, 3) or not np.all(np.isfinite(basis)):
            raise ValueError("subgroup basis must be a finite 3x3 matrix")
        det = float(np.linalg.det(basis))
        if abs(det) < 1e-12:
            raise ValueError("subgroup basis vectors must be linearly independent")
        if cartesian_tolerance <= 0 or not np.isfinite(cartesian_tolerance):
            raise ValueError("cartesian_tolerance must be finite and positive")

        child = build_supercell(structure, basis)
        inverse = np.linalg.inv(basis)
        translation_basis = np.asarray(
            subgroup_translation_lattice
            if subgroup_translation_lattice is not None
            else basis,
            dtype=float,
        )
        if (
            translation_basis.shape != (3, 3)
            or not np.all(np.isfinite(translation_basis))
            or abs(float(np.linalg.det(translation_basis))) < 1e-12
        ):
            raise ValueError(
                "subgroup primitive translation lattice must be a finite "
                "nonsingular 3x3 matrix"
            )
        translation_inverse = np.linalg.inv(translation_basis)
        parent_lattice = np.asarray(structure.lattice.matrix, dtype=float)
        translation_metric = Lattice(translation_basis @ parent_lattice)
        parent_translation_lattices: dict[int, _ExactTranslationLattice] = {}
        parent_coords = np.asarray(structure.frac_coords, dtype=float)
        child_as_parent = np.asarray(child.frac_coords, dtype=float) @ basis

        # A canonical microscopic column may have been extended from ISO's
        # sparse row domain onto the complete BUSH domain.  In that case the
        # working rows are intentionally different from the direct-query rows;
        # the immutable extension evidence, rather than the direct row digest,
        # is the scientific binding.  Validate complete extension groups once
        # up front and bind every current full-domain column to its recorded
        # digest.  Direct columns continue to use the stricter source-row and
        # source-vector checks below.
        extension_bindings = []
        for mode in modes:
            extension = mode.microscopic_domain_extension
            if extension is None:
                continue
            identity = mode.mode_identity
            provenance = mode.microscopic_provenance
            if identity is None or provenance is None:
                raise ValueError(
                    "domain-extended microscopic mode lacks source evidence"
                )
            if (
                microscopic_bush_column_digest(mode.bush_modes)
                != extension.extended_column_digest
            ):
                raise ValueError(
                    "domain-extended microscopic mode differs from its "
                    "verified full-domain column"
                )
            extension_bindings.append((identity, provenance, extension))
        validate_microscopic_domain_extension_bindings(extension_bindings)

        # Map each generated child atom back to one actual parent-cell site.
        # The species check prevents coincident special positions belonging to
        # different chemical orbits from being conflated.
        child_parent_index: list[int] = []
        for child_index, coord in enumerate(child_as_parent):
            best_index = -1
            best_distance = float("inf")
            for parent_index, parent_coord in enumerate(parent_coords):
                if child[child_index].species != structure[parent_index].species:
                    continue
                # Component-wise wrapping of fractional coordinates is not a
                # nearest-image construction in a skew cell.  The direct-
                # lattice metric must choose the integer image that minimizes
                # ``||(coord-parent_coord-n) B||``.
                distance = float(
                    structure.lattice.get_distance_and_image(
                        parent_coord, coord,
                    )[0]
                )
                if distance < best_distance:
                    best_index = parent_index
                    best_distance = distance
            if best_index < 0 or best_distance > cartesian_tolerance:
                raise ValueError(
                    "cannot map child atom to a chemically matching parent site: "
                    f"child index {child_index}, residual {best_distance:.6g} Å"
                )
            child_parent_index.append(best_index)

        parent_groups: list[tuple[str, str, set[int], list[int], dict]] = []
        for site_info in wyckoff_sites:
            indices = [int(value) for value in site_info["equivalent_indices"]]
            parent_groups.append((
                str(site_info["wyckoff_letter"]),
                str(site_info.get("orbit_id") or ""),
                set(indices),
                indices,
                site_info,
            ))

        operation_data: list[tuple[np.ndarray, np.ndarray]] = []
        for operation in subgroup_operations or ():
            rotation = np.asarray(operation.rotation, dtype=float)
            translation = np.asarray(
                operation.translation, dtype=float,
            )
            if rotation.shape != (3, 3) or translation.shape != (3,):
                raise ValueError("subgroup Seitz operations must contain R(3x3), t(3)")
            operation_data.append((rotation, translation))
        if not operation_data:
            operation_data.append((np.eye(3), np.zeros(3)))

        result: dict[str, np.ndarray] = {}
        for mode in modes:
            key = str(mode.amplitude_key or mode.irrep_label)
            if not key:
                raise ValueError("DISPLAY BUSH mode has no stable amplitude key")
            parent_arrows = np.zeros((len(child), 3), dtype=float)
            assigned = np.zeros(len(child), dtype=bool)
            target_letter = str(mode.wyckoff_site or "")
            target_orbit = str(getattr(mode, "wyckoff_orbit_id", "") or "")

            # ISO can print a microscopic column over fewer cells than the
            # selected subgroup basis. It may be repeated only by parent
            # translations with unit phase on every *active* arm proven by
            # that exact VECTOR source column. Without immutable active-arm
            # provenance, retain the conservative full-star rule.
            identity = mode.mode_identity
            provenance = mode.microscopic_provenance
            exact_k_coordinates: tuple[str, ...] = ()
            source_evidence_validated = False
            if (
                identity is not None
                and identity.status == "verified"
                and identity.source == "iso_microscopic"
                and len(identity.k_coordinates) == 3
            ):
                exact_k_coordinates = identity.k_coordinates
            elif identity is not None and identity.status == "partial":
                # A direction-resolved ISO column is intentionally ``partial``
                # until its whole Cartesian space has been compared with the
                # child-fixed reference space.  That comparison can itself
                # require repeating a sparse microscopic table over parent
                # translations, so requiring final ``verified`` status here
                # creates a circular dependency.  Authorize only the exact k
                # fact needed for that repetition, after validating every
                # immutable query/context field and the untouched source rows.
                if provenance is None or subgroup_context is None:
                    raise ValueError(
                        "partial microscopic mode requires source provenance "
                        "and subgroup context before child-cell mapping"
                    )
                provisional_identity = replace(
                    identity,
                    status="verified",
                    reason=None,
                )
                provisional_mode = replace(
                    mode,
                    mode_identity=provisional_identity,
                )
                try:
                    validate_microscopic_distortion_mode_source(
                        provisional_mode,
                        subgroup_context,
                    )
                except ValueError as exc:
                    raise ValueError(
                        "partial microscopic mode source evidence is inconsistent"
                    ) from exc
                if provenance.query_order is None:
                    raise ValueError(
                        "partial microscopic mode query assignment is unresolved"
                    )
                if microscopic_bush_row_digest(mode.bush_modes) != (
                    provenance.exact_row_digest
                ):
                    raise ValueError(
                        "partial microscopic mode rows differ from exact source provenance"
                    )
                if microscopic_bush_source_float_column_digest(mode.bush_modes) != (
                    provenance.source_float_column_digest
                ):
                    raise ValueError(
                        "partial microscopic mode vectors differ from exact source provenance"
                    )
                exact_k_coordinates = identity.k_coordinates
                source_evidence_validated = True

            active_k_coordinates: tuple[tuple[str, ...], ...] = ()
            if provenance is not None and provenance.source_active_k_coordinates:
                if subgroup_context is None:
                    raise ValueError(
                        "active-arm microscopic provenance requires subgroup context"
                    )
                if not source_evidence_validated:
                    try:
                        validate_microscopic_distortion_mode_source(
                            mode,
                            subgroup_context,
                        )
                    except ValueError as exc:
                        raise ValueError(
                            "active-arm microscopic source evidence is inconsistent"
                        ) from exc
                    if provenance.query_order is None:
                        raise ValueError(
                            "active-arm microscopic query assignment is unresolved"
                        )
                    if mode.microscopic_domain_extension is None:
                        if microscopic_bush_row_digest(mode.bush_modes) != (
                            provenance.exact_row_digest
                        ):
                            raise ValueError(
                                "active-arm microscopic rows differ from exact "
                                "source provenance"
                            )
                        if (
                            microscopic_bush_source_float_column_digest(
                                mode.bush_modes
                            )
                            != provenance.source_float_column_digest
                        ):
                            raise ValueError(
                                "active-arm microscopic vectors differ from exact "
                                "source provenance"
                            )
                recomputed_active = resolve_active_k_star_coordinates(
                    provenance.source_k_coordinates,
                    provenance.source_parent_sg,
                    provenance.source_direction_symbol,
                    provenance.source_column_index,
                    subgroup_translation_lattice=(
                        microscopic_subgroup_translation_lattice(subgroup_context)
                    ),
                )
                if (
                    recomputed_active.coordinates
                    != provenance.source_active_k_coordinates
                    or recomputed_active.resolution_kind
                    != provenance.source_active_k_resolution_kind
                ):
                    raise ValueError(
                        "active-arm microscopic provenance disagrees with its exact VECTOR"
                    )
                # For a compatible envelope this intentionally makes the
                # downstream unit-phase repetition gate stricter: repetition
                # is permitted only when every compatible arm has phase one.
                active_k_coordinates = provenance.source_active_k_coordinates

            phase_arms: tuple[tuple[Fraction, Fraction, Fraction], ...] = ()
            parent_translation_lattice: _ExactTranslationLattice | None = None
            if exact_k_coordinates:
                k_coords = exact_k_coordinates
                if active_k_coordinates:
                    phase_arms = tuple(
                        tuple(Fraction(value) for value in arm)
                        for arm in active_k_coordinates
                    )
                else:
                    # Legacy or hand-built sparse modes have no per-column
                    # arm proof and therefore retain the full exact star.
                    phase_arms = tuple(
                        k_star_fraction_vectors(
                            k_coords,
                            identity.parent_sg if identity is not None else 0,
                        )
                    )
                source_parent_sg = (
                    provenance.source_parent_sg
                    if provenance is not None
                    else identity.parent_sg
                )
                parent_translation_lattice = parent_translation_lattices.get(
                    source_parent_sg
                )
                if parent_translation_lattice is None:
                    centering = SpaceGroup.from_int_number(
                        source_parent_sg
                    ).symbol[0]
                    parent_translation_lattice = _make_exact_translation_lattice(
                        centering_primitive_matrix(centering),
                        parent_lattice,
                    )
                    parent_translation_lattices[source_parent_sg] = (
                        parent_translation_lattice
                    )

            for (
                letter,
                orbit_id,
                parent_index_set,
                parent_indices,
                site_info,
            ) in parent_groups:
                if target_letter and letter != target_letter:
                    continue
                if target_orbit and orbit_id != target_orbit:
                    continue
                bushes = [
                    bush for bush in mode.bush_modes
                    if str(bush.wyckoff_letter) == letter
                ]
                if not bushes:
                    continue
                reps: list[tuple[np.ndarray, np.ndarray]] = []
                for bush in bushes:
                    if len(bush.displacements) != 1:
                        raise ValueError(
                            "DISPLAY BUSH displacement columns must be split "
                            "before child-cell mapping"
                        )
                    point = self._resolve_rep_point(
                        bush, structure, parent_indices, site_info=site_info,
                    )
                    arrow = np.asarray(bush.displacements[0], dtype=float)
                    if point.shape != (3,) or arrow.shape != (3,):
                        raise ValueError("BUSH point and displacement must be 3-vectors")
                    for rotation, translation in operation_data:
                        # Seitz convention is x' = R x + t for column
                        # coordinates.  Points/arrows here are row vectors.
                        reps.append((
                            point @ rotation.T + translation,
                            arrow @ rotation.T,
                        ))

                group_children = [
                    index for index, parent_index in enumerate(child_parent_index)
                    if parent_index in parent_index_set
                ]
                for child_index in group_children:
                    coord = child_as_parent[child_index]
                    matches: list[np.ndarray] = []
                    for point, arrow in reps:
                        translation_delta = (coord - point) @ translation_inverse
                        # ``translation_delta`` is fractional in an arbitrary
                        # primitive basis of the subgroup translation lattice.
                        # Rounding its components is not a nearest-image
                        # construction when that basis is skew or unreduced;
                        # it would make BUSH coverage depend on the chosen
                        # (GL(3,Z)-equivalent) lattice basis.  Minimize in the
                        # actual direct-lattice metric instead.
                        distance = float(
                            translation_metric.get_distance_and_image(
                                translation_delta, np.zeros(3),
                            )[0]
                        )
                        if distance <= cartesian_tolerance:
                            matches.append(arrow)
                            continue
                        if not phase_arms or parent_translation_lattice is None:
                            continue
                        parent_delta = coord - point
                        exact_parent_translation = _snap_to_exact_translation(
                            parent_delta,
                            parent_translation_lattice,
                            cartesian_tolerance,
                        )
                        if exact_parent_translation is None:
                            continue
                        if _has_exact_unit_phase(
                            phase_arms,
                            exact_parent_translation,
                        ):
                            matches.append(arrow)
                    if not matches:
                        raise ValueError(
                            "DISPLAY BUSH table does not cover child atom "
                            f"{child_index} in Wyckoff orbit {letter!r}"
                        )
                    reference = matches[0]
                    if any(
                        float(np.linalg.norm((other - reference) @ parent_lattice))
                        > cartesian_tolerance
                        for other in matches[1:]
                    ):
                        raise ValueError(
                            "DISPLAY BUSH representatives assign conflicting arrows "
                            f"to child atom {child_index}"
                        )
                    parent_arrows[child_index] = reference
                    assigned[child_index] = True

            if not np.any(assigned):
                raise ValueError(
                    f"DISPLAY BUSH mode {key!r} did not map to any child atom"
                )
            # Preserve the established UI amplitude convention: the largest
            # parent-fractional arrow norm is one.  The physically essential
            # direction/phase is exact before the child-coordinate transform.
            scale = float(np.max(np.linalg.norm(parent_arrows, axis=1)))
            if scale <= 1e-12:
                raise ValueError(f"DISPLAY BUSH mode {key!r} is identically zero")
            result[key] = (parent_arrows / scale) @ inverse
        return result

    def map_microscopic_columns_to_supercell(
        self,
        structure: Structure,
        wyckoff_sites: list[dict],
        modes: Sequence[DistortionMode],
        basis_vectors: list[list[float]],
        *,
        subgroup_context: SubgroupInfo,
        frame_id: str,
        atom_ids: Sequence[str],
        atom_order_id: str | None = None,
        cartesian_tolerance: float = 1e-4,
        subgroup_operations: Sequence[object] | None = None,
        subgroup_translation_lattice: Sequence[Sequence[float]] | None = None,
    ) -> tuple[MappedMicroscopicColumn, ...]:
        """Map verified ISO columns while retaining their exact source evidence.

        The target frame and atom-order IDs are mandatory integration inputs.
        This method does not derive either identity from array position.
        """

        if not str(frame_id).strip():
            raise ValueError("mapped microscopic target frame ID is required")
        target_atom_ids = tuple(str(value).strip() for value in atom_ids)
        if not target_atom_ids or any(not value for value in target_atom_ids):
            raise ValueError("mapped microscopic target atom IDs are required")
        if len(set(target_atom_ids)) != len(target_atom_ids):
            raise ValueError("mapped microscopic target atom IDs must be unique")
        target_atom_order_id = microscopic_atom_order_id(target_atom_ids)
        if atom_order_id is not None and str(atom_order_id) != target_atom_order_id:
            raise ValueError("mapped microscopic target atom-order ID does not match its atom IDs")
        ordered_modes = list(modes)
        if not ordered_modes:
            return ()
        keys = [str(mode.amplitude_key or mode.irrep_label) for mode in ordered_modes]
        if any(not key for key in keys) or len(keys) != len(set(keys)):
            raise ValueError("mapped microscopic amplitude keys must be nonempty and unique")
        for mode in ordered_modes:
            identity = mode.mode_identity
            provenance = mode.microscopic_provenance
            if identity is None or provenance is None:
                raise ValueError("microscopic mode lacks verified source-column identity")
            validate_microscopic_distortion_mode_source(mode, subgroup_context)
            mapped_basis = tuple(
                tuple(str(as_fraction(value)) for value in row)
                for row in basis_vectors
            )
            if mapped_basis != provenance.source_subgroup_basis:
                raise ValueError(
                    "mapped child basis differs from microscopic query context"
                )
        extension_bindings = []
        for mode in ordered_modes:
            extension = mode.microscopic_domain_extension
            identity = mode.mode_identity
            provenance = mode.microscopic_provenance
            if identity is None or provenance is None:  # guarded above; narrows types
                raise ValueError("microscopic source evidence unexpectedly disappeared")
            if extension is None:
                if microscopic_bush_row_digest(mode.bush_modes) != provenance.exact_row_digest:
                    raise ValueError(
                        "direct microscopic mode rows differ from exact source provenance"
                    )
                if (
                    microscopic_bush_source_float_column_digest(mode.bush_modes)
                    != provenance.source_float_column_digest
                ):
                    raise ValueError(
                        "direct microscopic displacements differ from exact source provenance"
                    )
            elif (
                extension.extended_column_digest
                != microscopic_bush_column_digest(mode.bush_modes)
            ):
                raise ValueError("domain extension does not match the current full column")
            extension_bindings.append((identity, provenance, extension))
        validate_microscopic_domain_extension_bindings(extension_bindings)

        mapped = self.map_bush_modes_to_supercell(
            structure,
            wyckoff_sites,
            ordered_modes,
            basis_vectors,
            cartesian_tolerance=cartesian_tolerance,
            subgroup_operations=subgroup_operations,
            subgroup_translation_lattice=subgroup_translation_lattice,
            subgroup_context=subgroup_context,
        )
        if set(mapped) != set(keys):
            raise ValueError("mapped microscopic output does not cover every source column")
        result: list[MappedMicroscopicColumn] = []
        for mode, key in zip(ordered_modes, keys, strict=True):
            identity = mode.mode_identity
            provenance = mode.microscopic_provenance
            if identity is None or provenance is None:  # guarded above; narrows types
                raise ValueError("microscopic source evidence unexpectedly disappeared")
            bound_provenance = provenance.bind_mapping(
                frame_id=str(frame_id),
                atom_ids=target_atom_ids,
            )
            result.append(MappedMicroscopicColumn(
                mode_id=identity.stable_token,
                amplitude_key=key,
                frame_id=str(frame_id),
                atom_order_id=target_atom_order_id,
                atom_ids=target_atom_ids,
                displacements=mapped[key],
                mode_identity=identity,
                provenance=bound_provenance,
                domain_extension_evidence=mode.microscopic_domain_extension,
                relation_status="canonical",
                current_over_canonical_signed_scale=1.0,
                cartesian_residual_angstrom=0.0,
            ))
        if len({column.mode_id for column in result}) != len(result):
            raise ValueError("mapped microscopic mode identities are not unique")
        return tuple(result)

    # ----------------------------------------------------------------
    # 代表点解析 / 格点等价
    # ----------------------------------------------------------------

    def _lattice_setup(self, structure: Structure) -> None:
        """缓存母相惯用格矩阵与原胞格矩阵（模格点等价判定用）。

        用 structure.get_primitive_structure()（纯胞约化，与母相格点同一格）
        而非 get_primitive_standard_structure()（可能给出旋转后的不同格点，
        实测对 R3m 等带心结构会把带心平移判成非格点）。
        """
        prim = structure.get_primitive_structure()
        self._prim_matrix = np.asarray(prim.lattice.matrix, dtype=float)
        self._conv_matrix = np.asarray(structure.lattice.matrix, dtype=float)

    def _resolve_rep_point(
        self,
        bush: BushMode,
        structure: Structure,
        indices: list[int],
        *,
        site_info: dict | None = None,
    ) -> np.ndarray:
        """把 BUSH 代表点解析为具体分数坐标。

        带自由参数的代表点（如 ``(0,0,z)``、``(x,2x,z)``、``(x+1/2,-y+1,-z)``）
        中，字母参数按其表示的坐标分量取值：x/y/z（或 a/b/g）分别取位点
        第一个原子的第 0/1/2 个分数坐标分量。否则相反符号的代表点
        （``(0,0,z)`` 与 ``(0,0,-z)``）会解析成同一点，导致相反位移的
        代表无法区分。
        """
        _LETTER_AXIS = {"x": 0, "y": 1, "z": 2, "a": 0, "b": 1, "g": 2}
        tokens = [str(t).strip() for t in (bush.point_raw or [])]
        if not tokens:
            # Legacy BUSH objects predate the exact printed coordinate tokens.
            # Verified microscopic columns always carry ``point_raw`` because
            # their immutable source digest is defined on those rows.
            return np.asarray(bush.point, dtype=float)
        if len(tokens) != 3 or any(not token for token in tokens):
            raise ValueError("BUSH point_raw must be one complete coordinate 3-vector")
        if not any(_has_letter(t) for t in tokens):
            exact_point = self._eval_tokens(tokens, {})
            printed_point = np.asarray(bush.point, dtype=float)
            if (
                printed_point.shape != (3,)
                or not np.all(np.isfinite(printed_point))
                or not np.allclose(
                    printed_point,
                    exact_point,
                    rtol=0.0,
                    atol=1e-12,
                )
            ):
                raise ValueError(
                    "numeric constant BUSH point differs from authoritative point_raw"
                )
            return exact_point
        atom0 = np.asarray(structure[indices[0]].frac_coords, dtype=float)
        standard = None
        parameters: dict[str, float] = {}
        if site_info is not None:
            raw_standard = site_info.get("standard_representative_frac_coords")
            if raw_standard is not None:
                standard = np.asarray(raw_standard, dtype=float)
                if standard.shape != (3,) or not np.all(np.isfinite(standard)):
                    raise ValueError(
                        "standard Wyckoff representative must be a finite "
                        "fractional 3-vector"
                    )
            raw_parameters = site_info.get("standard_representative_parameters")
            if raw_parameters is not None:
                if not isinstance(raw_parameters, dict):
                    raise ValueError(
                        "standard Wyckoff representative parameters must be a mapping"
                    )
                try:
                    parameters = {
                        str(name): float(value)
                        for name, value in raw_parameters.items()
                    }
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        "standard Wyckoff representative parameters must be numeric"
                    ) from exc
                if not all(np.isfinite(value) for value in parameters.values()):
                    raise ValueError(
                        "standard Wyckoff representative parameters must be finite"
                    )
        reference = standard if standard is not None else atom0
        letters: dict[str, float] = dict(parameters)
        for tok in tokens:
            for name in re.findall(r"[a-zA-Z]+", tok):
                if name in letters:
                    continue
                axis = _LETTER_AXIS.get(name)
                if axis is None:
                    raise ValueError(
                        f"unsupported BUSH coordinate parameter {name!r}"
                    )
                letters[name] = float(reference[axis])
        point = self._eval_tokens(tokens, letters)
        if point.shape != (3,) or not np.all(np.isfinite(point)):
            raise ValueError("resolved BUSH point must be one finite coordinate 3-vector")
        return point

    @staticmethod
    def _eval_token(tok: str, letters: dict[str, float]) -> float:
        """求值单个坐标 token（可含常数项与线性组合）。

        支持 ``1/2``、``x``、``-x``、``2x``、``-y+1``、``x-1/2``、
        ``-x+3/2``、``y+1/2``、``-2a+1`` 等形式。
        """
        s = tok.replace(" ", "")
        if not s:
            raise ValueError("empty BUSH coordinate token")
        term_pattern = re.compile(
            r"(?P<sign>[+-]?)"
            r"(?:"
            r"(?P<coeff>(?:\d+(?:\.\d*)?|\.\d+)(?:/[1-9]\d*)?)"
            r"(?P<alpha>[a-zA-Z]+)?"
            r"|(?P<bare_alpha>[a-zA-Z]+)"
            r")"
        )
        total = 0.0
        position = 0
        first = True
        while position < len(s):
            match = term_pattern.match(s, position)
            if match is None or match.start() != position:
                raise ValueError(f"invalid BUSH coordinate token {tok!r}")
            sign = match.group("sign")
            if not first and not sign:
                raise ValueError(f"invalid BUSH coordinate token {tok!r}")
            coeff = match.group("coeff")
            alpha = match.group("alpha") or match.group("bare_alpha")
            if alpha:
                if alpha not in letters:
                    raise ValueError(
                        f"unresolved BUSH coordinate parameter {alpha!r}"
                    )
                base = letters[alpha] * (
                    float(Fraction(coeff)) if coeff else 1.0)
            else:
                if coeff is None:  # defensive: the grammar always provides one branch
                    raise ValueError(f"invalid BUSH coordinate token {tok!r}")
                base = float(Fraction(coeff))
            total += base if sign != "-" else -base
            position = match.end()
            first = False
        if position != len(s):
            raise ValueError(f"invalid BUSH coordinate token {tok!r}")
        return total

    @staticmethod
    def _eval_tokens(tokens: list[str], letters: dict[str, float]) -> np.ndarray:
        """按字母参数表求值坐标 token 列表。"""
        return np.asarray(
            [DistortionMapper._eval_token(t, letters) for t in tokens],
            dtype=float,
        )

    def _lattice_equiv(self, a: np.ndarray, b: np.ndarray,
                       atol: float = 1e-3) -> bool:
        """两个分数坐标是否相差一个母相格点（a - b ∈ L）。"""
        delta_cart = (a - b) @ self._conv_matrix
        n = delta_cart @ np.linalg.inv(self._prim_matrix)
        return bool(np.allclose(n, np.round(n), atol=atol))

    @staticmethod
    def _periodic_equiv(a: np.ndarray, b: np.ndarray,
                        atol: float = 1e-3) -> bool:
        """两个分数坐标是否周期等价（相差整数平移，mod 1 相等）。

        用于把 BUSH 代表点匹配到具体原子：代表点与原子分数坐标在 [0,1) 内
        重合即视为同一原子（如带心副本相位反号代表 `(-1/2,1/2,1/2)` 与
        体心原子 `(1/2,1/2,1/2)` mod 1 相等）。
        """
        delta = a - b
        return bool(np.allclose(delta - np.round(delta), 0.0, atol=atol))

    def _periodic_distance(self, a: np.ndarray, b: np.ndarray) -> float:
        """Return the shortest periodic distance in the direct-lattice metric.

        Fractional Euclidean distance has no crystallographic meaning unless
        the basis is orthonormal.  ``get_distance_and_image`` minimizes over
        lattice translations using the actual parent lattice, including skew
        and anisotropic cells.
        """
        if self._conv_matrix is None:
            raise RuntimeError("parent lattice is not initialized")
        return float(Lattice(self._conv_matrix).get_distance_and_image(a, b)[0])


def _has_letter(token: str) -> bool:
    """token 是否含字母（自由参数标记，如 x/y/z）。"""
    return bool(re.search(r"[a-zA-Z]", str(token)))
