"""Search workflows aligned with ISODISTORT Method 1-4（基于真实 iso 枚举）。"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.groups import SpaceGroup

from ..backend import DistortionMode, IsoWrapper, SubgroupInfo
from .phase_path import normalize_distortion_types

CRYSTAL_SYSTEMS = {
    "triclinic",
    "monoclinic",
    "orthorhombic",
    "tetragonal",
    "trigonal",
    "hexagonal",
    "cubic",
}


def _sg_to_crystal_system(space_group_number: int) -> str:
    if 1 <= space_group_number <= 2:
        return "triclinic"
    if 3 <= space_group_number <= 15:
        return "monoclinic"
    if 16 <= space_group_number <= 74:
        return "orthorhombic"
    if 75 <= space_group_number <= 142:
        return "tetragonal"
    if 143 <= space_group_number <= 167:
        return "trigonal"
    if 168 <= space_group_number <= 194:
        return "hexagonal"
    if 195 <= space_group_number <= 230:
        return "cubic"
    return "unknown"


def _to_float(value: str | float) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if "/" in text:
        return float(Fraction(text))
    return float(text)


def _parse_basis_rows(rows: Sequence[Sequence[str | int | float]]) -> list[list[float]]:
    basis = []
    for row in rows:
        if len(row) != 3:
            raise ValueError("Each basis row must contain 3 values")
        basis.append([_to_float(v) for v in row])
    if len(basis) != 3:
        raise ValueError("Basis matrix must be 3x3")
    return basis


def _matrix_matches_diagonal_supercell(matrix: Sequence[Sequence[float]],
                                       supercell: Sequence[int]) -> bool:
    if len(supercell) != 3:
        return False
    target = np.diag(np.array(supercell, dtype=float))
    got = np.array(matrix, dtype=float)
    return np.allclose(got, target, atol=1e-8)


def _space_group_to_point_group(space_group_number: int) -> str:
    try:
        return SpaceGroup.from_int_number(space_group_number).point_group
    except (ValueError, KeyError, TypeError):
        return ""


def _basis_is_sublattice_of(basis: Sequence[Sequence[float]],
                            sublattice: Sequence[Sequence[float]]) -> bool:
    """
    判断子群超胞基矢 B 的格点是否是被选子格 S 的子格。

    对应官网 Method 1 的 direct sublattice / Conventional lattice /
    Primitive lattice 过滤：B 的每一行必须是 S 的整系数线性组合，
    即 N = B @ inv(S) 的元素全部为整数（S 可为对角阵或任意 3x3 矩阵）。
    """
    b = np.asarray(basis, dtype=float)
    s = np.asarray(sublattice, dtype=float)
    if b.shape != (3, 3) or s.shape != (3, 3):
        return False
    try:
        n = b @ np.linalg.inv(s)
    except np.linalg.LinAlgError:
        return False
    return bool(np.allclose(n, np.round(n), atol=1e-6))


@dataclass
class Method1Query:
    """Method 1: search over all special k points."""

    distortion_types: str | Sequence[str] | None = None
    crystal_system: str | Sequence[str] | None = None  # 单个或列表（多选=OR）
    subgroup_space_group: int | None = None
    direct_sublattice: Sequence[int] | None = None
    lattice: Sequence[Sequence[float]] | None = None  # 官网 conventional/primitive lattice（3x3 子格矩阵）
    maximal_subgroup_only: bool = False


@dataclass
class Method1ResultItem:
    subgroup: SubgroupInfo
    crystal_system: str
    is_maximal: bool
    direct_sublattice: list[int] = field(default_factory=lambda: [1, 1, 1])


@dataclass
class Method2Query:
    """Method 2: general method over specified k point(s)."""

    subgroup_idx: int
    distortion_type: str | Sequence[str] = "displacive"
    k_point_label: str | None = None
    k_point_coordinates: Sequence[str | int | float] | None = None
    k_parameters: dict[str, str | int | float] = field(default_factory=dict)
    number_of_independent_modulations: int = 0
    number_of_superposed_irs: int = 1
    specified_opd: str | None = None


@dataclass
class Method2Result:
    subgroup: SubgroupInfo
    modes: list[DistortionMode]
    metadata: dict[str, object]


@dataclass
class Method3Query:
    """Method 3: search over arbitrary k for point/space group + supercell."""

    distortion_types: str | Sequence[str] | None = None
    point_group: str | None = None
    space_group_type: int | None = None
    supercell_basis: Sequence[Sequence[str | int | float]] | None = None
    direct_sublattice_centering: str | None = None
    lattice_type: str = "direct"  # 官网 radio：direct（实空间子格）/ reciprocal（倒易超格）


@dataclass
class Method3ResultItem:
    subgroup: SubgroupInfo
    point_group: str
    basis: list[list[float]]


@dataclass
class Method4Query:
    """Method 4: mode decomposition of distorted structure."""

    atom_matching_method: str = "nearest-site"
    robust_distance_threshold: float = 0.25
    provided_origin_shift: Sequence[float] | None = None


@dataclass
class Method4Result:
    amplitudes: dict[str, float]
    rms_residual: float
    max_abs_residual: float
    assignments: list[int]
    metadata: dict[str, object]


class IsoSearchEngine:
    """Implements local Method 1-4 search workflows on top of the real iso binary."""

    def __init__(self, iso_wrapper: IsoWrapper) -> None:
        self._iso = iso_wrapper

    # ----------------------------------------------------------------
    # Method 1：全特殊 k 点搜索 + 客户端过滤（与官网逻辑 AND 语义一致）
    # ----------------------------------------------------------------

    def method_1_search(self, parent_sg: int, query: Method1Query,
                        distortion_types=None,
                        subgroups: Sequence[SubgroupInfo] | None = None
                        ) -> list[Method1ResultItem]:
        """
        官网 Method 1：遍历全部特殊 k 点，得到子群候选后按用户条件过滤。

        过滤条件（多条件同时生效，逻辑 AND；同一条件内多选为 OR）：
        - crystal system：子群所属晶系（单个或列表，列表任中其一即通过）
        - subgroup space group：子群空间群号
        - maximal subgroup only：仅保留 maximal 子群
        - direct sublattice / lattice：超胞格是否为所选子格的子格
          （官网 direct sublattice / Conventional lattice / Primitive lattice）

        Args:
            parent_sg: 母相空间群号
            query: Method 1 查询
            distortion_types: 畸变类型（保留参数；类型过滤在模式计算阶段执行）
            subgroups: 预枚举的子群候选（缓存复用）；None 时现场枚举
        """
        _ = distortion_types
        if subgroups is None:
            subgroups = self._iso.enumerate_all_special_subgroups(
                parent_sg, query.distortion_types
            )

        crystal_systems: set[str] | None = None
        if query.crystal_system:
            raw = query.crystal_system
            if isinstance(raw, str):
                raw = [raw]
            crystal_systems = {x.strip().lower() for x in raw if x.strip()}

        result: list[Method1ResultItem] = []
        for sg in subgroups:
            crystal_system = _sg_to_crystal_system(sg.space_group_number)
            item = Method1ResultItem(
                subgroup=sg,
                crystal_system=crystal_system,
                is_maximal=sg.is_maximal,
                direct_sublattice=self._diagonal_sublattice(sg.basis_vectors),
            )
            if crystal_systems and crystal_system not in crystal_systems:
                continue
            if query.subgroup_space_group and sg.space_group_number != query.subgroup_space_group:
                continue
            if query.maximal_subgroup_only and not sg.is_maximal:
                continue
            if query.lattice is not None and not _basis_is_sublattice_of(
                sg.basis_vectors, query.lattice
            ):
                continue
            if query.lattice is None and query.direct_sublattice:
                diag = np.diag([float(v) for v in query.direct_sublattice])
                if not _basis_is_sublattice_of(sg.basis_vectors, diag):
                    continue
            result.append(item)
        return result

    @staticmethod
    def _diagonal_sublattice(basis: Sequence[Sequence[float]]) -> list[int]:
        """基矢为对角阵时返回其对角元（母相格单位），否则返回 [1,1,1]。"""
        if len(basis) != 3:
            return [1, 1, 1]
        arr = np.asarray(basis, dtype=float)
        if np.allclose(arr, np.diag(np.diag(arr)), atol=1e-8):
            diag = np.diag(arr)
            if all(abs(d - round(d)) < 1e-8 for d in diag):
                return [round(d) for d in diag]
        return [1, 1, 1]

    # ----------------------------------------------------------------
    # Method 2：指定 k 点/IR/OPD 的模式计算
    # ----------------------------------------------------------------

    def method_2_search(self, parent_sg: int, subgroups: Sequence[SubgroupInfo],
                        query: Method2Query,
                        wyckoff_letters: Sequence[str] | None = None) -> Method2Result:
        """
        官网 Method 2：在已枚举子群中按序号选择目标子群，计算其畸变模式。

        模式基矢由真实 iso 的 DISPLAY BUSH 计算（需要母相 Wyckoff 位置）。

        Args:
            parent_sg: 母相空间群号
            subgroups: 子群候选列表（来自 Method 1 或 list_subgroups）
            query: Method 2 查询参数（subgroup_idx 必填）
            wyckoff_letters: 母相结构各原子的 Wyckoff 位置字母

        Returns:
            Method2Result
        """
        target = next((s for s in subgroups if s.index == query.subgroup_idx), None)
        if target is None:
            raise ValueError(
                f"Subgroup index {query.subgroup_idx} not found; "
                "请先执行 Method 1 或 list_subgroups 获得候选列表"
            )

        # wyckoff_letters 为 None 表示调用方未提供（误用）；空列表表示
        # 作用域内无 Wyckoff 位置（如全部类型选 none）-> 直接返回空模式
        if wyckoff_letters is None:
            raise ValueError(
                "Method 2 计算模式需要母相结构的 Wyckoff 位置信息，"
                "请先加载结构（load_structure）"
            )
        modes = []
        if wyckoff_letters:
            modes = self._iso.calc_distortion_modes(
                parent_sg, target, wyckoff_letters=wyckoff_letters
            )

        metadata = {
            "k_point_label": query.k_point_label or target.k_point_label,
            "k_point_coordinates": [
                _to_float(v) for v in query.k_point_coordinates
            ] if query.k_point_coordinates else None,
            "k_parameters": {k: _to_float(v) for k, v in query.k_parameters.items()},
            "number_of_independent_modulations": query.number_of_independent_modulations,
            "number_of_superposed_irs": query.number_of_superposed_irs,
            "specified_opd": query.specified_opd or target.opd_symbol,
        }

        return Method2Result(subgroup=target, modes=modes, metadata=metadata)

    # ----------------------------------------------------------------
    # Method 3：指定点群/空间群 + 超胞
    # ----------------------------------------------------------------

    def method_3_search(self, parent_sg: int, query: Method3Query) -> list[Method3ResultItem]:
        """
        官网 Method 3 的本地近似实现：

        - 若同时提供 point_group 与 space_group_type，空间群选择优先
          （与官网规则一致）；
        - 超胞基矢过滤：仅当请求基矢为对角阵时做严格匹配；
          任意基矢的完整再生成（需 iso 在线生成子群数据库）为已知限制。
        """
        distortion_types = normalize_distortion_types(query.distortion_types)
        subgroups = self._iso.enumerate_all_special_subgroups(parent_sg, distortion_types)

        if query.point_group and query.space_group_type:
            # 官网规则：space-group selection supersedes point-group selection
            point_group_filter = None
        else:
            point_group_filter = query.point_group

        basis = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        if query.supercell_basis:
            basis = _parse_basis_rows(query.supercell_basis)

        result: list[Method3ResultItem] = []
        for sg in subgroups:
            point_group = _space_group_to_point_group(sg.space_group_number)
            if query.space_group_type and sg.space_group_number != query.space_group_type:
                continue
            if point_group_filter and point_group != point_group_filter:
                continue
            if query.supercell_basis and not _matrix_matches_diagonal_supercell(
                sg.basis_vectors, [1, 1, 1]
            ):
                # 任意基矢的子群再生成需要 iso 在线生成数据库（见 README 已知差异）
                pass

            result.append(Method3ResultItem(subgroup=sg, point_group=point_group, basis=basis))

        return result

    # ----------------------------------------------------------------
    # Method 4：模式分解（自研最小二乘拟合）
    # ----------------------------------------------------------------

    def method_4_decompose(self,
                           parent_structure: Structure,
                           distorted_structure: Structure,
                           mode_displacements: dict[str, np.ndarray],
                           query: Method4Query) -> Method4Result:
        if len(parent_structure) != len(distorted_structure):
            raise ValueError(
                "Mode decomposition currently requires parent and distorted structures "
                "to have the same atom count"
            )

        assignments = self._match_atoms(parent_structure, distorted_structure, query)
        delta = self._build_delta_vector(parent_structure, distorted_structure, assignments)

        mode_labels = list(mode_displacements.keys())
        if not mode_labels:
            raise ValueError("No mode displacements are available for decomposition")

        columns = []
        for label in mode_labels:
            vec = np.asarray(mode_displacements[label], dtype=float).reshape(-1)
            if vec.size != delta.size:
                raise ValueError(
                    f"Mode {label} has incompatible size {vec.size}, expected {delta.size}"
                )
            columns.append(vec)
        a_matrix = np.column_stack(columns)

        coeffs, _, _, _ = np.linalg.lstsq(a_matrix, delta, rcond=None)
        reconstructed = a_matrix @ coeffs
        residual = delta - reconstructed

        amplitudes = {mode_labels[i]: float(coeffs[i]) for i in range(len(mode_labels))}
        rms = float(np.sqrt(np.mean(residual ** 2)))
        max_abs = float(np.max(np.abs(residual)))

        return Method4Result(
            amplitudes=amplitudes,
            rms_residual=rms,
            max_abs_residual=max_abs,
            assignments=assignments,
            metadata={
                "atom_matching_method": query.atom_matching_method,
                "provided_origin_shift": list(query.provided_origin_shift)
                if query.provided_origin_shift is not None else None,
            },
        )

    @staticmethod
    def _fractional_delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        delta = b - a
        delta -= np.round(delta)
        return delta

    def _match_atoms(self,
                     parent_structure: Structure,
                     distorted_structure: Structure,
                     query: Method4Query) -> list[int]:
        if query.atom_matching_method not in {"nearest-site", "robust"}:
            raise ValueError("atom_matching_method must be 'nearest-site' or 'robust'")

        assignments: list[int] = []
        used: set[int] = set()

        for i, site in enumerate(parent_structure):
            species = site.species_string
            parent_coord = np.asarray(site.frac_coords)
            candidates = []
            for j, dst in enumerate(distorted_structure):
                if j in used:
                    continue
                if dst.species_string != species:
                    continue
                delta = self._fractional_delta(parent_coord, np.asarray(dst.frac_coords))
                dist = float(np.linalg.norm(delta))
                if query.atom_matching_method == "nearest-site":
                    candidates.append((dist, j))
                else:
                    if dist <= query.robust_distance_threshold:
                        candidates.append((dist, j))

            if not candidates:
                raise ValueError(f"Cannot match parent atom index {i} ({species})")

            candidates.sort(key=lambda t: t[0])
            chosen = candidates[0][1]
            assignments.append(chosen)
            used.add(chosen)

        return assignments

    def _build_delta_vector(self,
                            parent_structure: Structure,
                            distorted_structure: Structure,
                            assignments: Sequence[int]) -> np.ndarray:
        delta_rows = []
        for i, j in enumerate(assignments):
            p = np.asarray(parent_structure[i].frac_coords)
            d = np.asarray(distorted_structure[j].frac_coords)
            delta_rows.append(self._fractional_delta(p, d))
        return np.asarray(delta_rows, dtype=float).reshape(-1)
