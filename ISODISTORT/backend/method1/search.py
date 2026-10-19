"""Method 1 query, result and solver contract."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import numpy as np
from backend.utils.lattice import rational_matrix, same_lattice_in_point_group_orbit
from backend.wrappers import IsoWrapper, SubgroupInfo
from backend.tables.spacegroup_meta import crystal_system_for_number


def _basis_is_same_lattice(
    basis: Sequence[Sequence[float]],
    selected: Sequence[Sequence[float]],
    parent_rotations: Sequence[np.ndarray] | None = None,
) -> bool:
    """Exact GL(3,Z) class under fractional parent point-group rotations."""
    try:
        rotations = () if parent_rotations is None else parent_rotations
        return same_lattice_in_point_group_orbit(
            rational_matrix(basis),
            rational_matrix(selected),
            [rational_matrix(rotation) for rotation in rotations],
        )
    except (TypeError, ValueError, ZeroDivisionError):
        return False


@dataclass
class Method1Query:
    """Method 1: search over all special k points."""

    distortion_types: str | Sequence[str] | None = None
    crystal_system: str | Sequence[str] | None = None  # 单个或列表（多选=OR）
    subgroup_space_group: int | None = None
    lattice: Sequence[Sequence[float]] | None = None  # conventional 格点类（3x3 矩阵）
    maximal_subgroup_only: bool = False
    # 母相点群旋转（分数坐标）；lattice 过滤时与官网一样合并点群轨道
    parent_rotations: Sequence[Sequence[Sequence[float]]] | None = None


@dataclass
class Method1ResultItem:
    subgroup: SubgroupInfo
    crystal_system: str
    is_maximal: bool


class Method1Solver:
    """Own Method 1 computation; application sessions compose Method services."""

    def __init__(self, iso_wrapper: IsoWrapper) -> None:
        self._iso = iso_wrapper

    def search(
        self,
        parent_sg: int,
        query: Method1Query,
        distortion_types=None,
        subgroups: Sequence[SubgroupInfo] | None = None,
    ) -> list[Method1ResultItem]:
        """
        官网 Method 1：遍历全部特殊 k 点，得到子群候选后按用户条件过滤。

        过滤条件（多条件同时生效，逻辑 AND；同一条件内多选为 OR）：
        - crystal system：子群所属晶系（单个或列表，列表任中其一即通过）
        - subgroup space group：子群空间群号
        - maximal subgroup only：仅保留 maximal 子群
        - lattice：超胞格是否与所选格子相同（含母相点群旋转轨道）
          （官网 Conventional lattice / Primitive lattice）

        Args:
            parent_sg: 母相空间群号
            query: Method 1 查询
            distortion_types: 畸变类型（保留参数；类型过滤在模式计算阶段执行）
            subgroups: 预枚举的子群候选（缓存复用）；None 时现场枚举
        """
        _ = distortion_types
        if subgroups is None:
            subgroups = self._iso.enumerate_all_special_subgroups(parent_sg, query.distortion_types)

        crystal_systems: set[str] | None = None
        if query.crystal_system:
            raw = query.crystal_system
            if isinstance(raw, str):
                raw = [raw]
            crystal_systems = {x.strip().lower() for x in raw if x.strip()}

        result: list[Method1ResultItem] = []
        for sg in subgroups:
            crystal_system = crystal_system_for_number(sg.space_group_number)
            item = Method1ResultItem(
                subgroup=sg,
                crystal_system=crystal_system,
                is_maximal=sg.is_maximal,
            )
            if crystal_systems and crystal_system not in crystal_systems:
                continue
            if query.subgroup_space_group and sg.space_group_number != query.subgroup_space_group:
                continue
            if query.maximal_subgroup_only and not sg.is_maximal:
                continue
            if query.lattice is not None:
                rots = None
                if query.parent_rotations:
                    rots = [np.asarray(r, dtype=float) for r in query.parent_rotations]
                if not _basis_is_same_lattice(sg.basis_vectors, query.lattice, parent_rotations=rots):
                    continue
            result.append(item)
        return result
