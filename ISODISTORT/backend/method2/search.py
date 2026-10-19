"""Method 2 query, result and solver contract."""

from __future__ import annotations
from backend.models.phase_path import validate_modulation_count

from collections.abc import Sequence
from dataclasses import dataclass
from numbers import Integral
from pymatgen.core import Structure
from backend.wrappers import DistortionMode, IsoWrapper, SubgroupInfo


@dataclass
class Method2Query:
    """Method 2: general method over a selected subgroup (k/IR/OPD 由子群自身携带)。"""

    subgroup_idx: int
    distortion_type: str | Sequence[str] = "displacive"
    number_of_independent_modulations: int = 0  # 0 = 3D lock-in; n = (3+n)D harmonics
    # This scalar path represents exactly one selected IR/OPD.  Coupled
    # superposed IRs require a component-aware identity and ISO's coupled
    # isotropy enumeration; treating their single-IR tables as a union is
    # scientifically wrong, so the engine rejects any value other than one.
    number_of_superposed_irreps: int = 1


def validate_single_irrep_count(value: object) -> int:
    """Return one for the only currently supported Method 2 component count."""

    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError("number_of_superposed_irreps must be a positive integer")
    count = int(value)
    if count < 1:
        raise ValueError("number_of_superposed_irreps must be a positive integer")
    if count != 1:
        raise NotImplementedError(
            "Coupled Method 2 superposed-IR calculation is not implemented; "
            "the single-IR subgroup tables cannot be combined as a union"
        )
    return count


@dataclass
class Method2Result:
    subgroup: SubgroupInfo
    modes: list[DistortionMode]
    metadata: dict[str, object]


class Method2Solver:
    """Own Method 2 computation; application sessions compose Method services."""

    def __init__(self, iso_wrapper: IsoWrapper) -> None:
        self._iso = iso_wrapper

    def search(
        self,
        parent_sg: int,
        subgroups: Sequence[SubgroupInfo],
        query: Method2Query,
        wyckoff_letters: Sequence[str] | None = None,
        *,
        wyckoff_orbit_ids: Sequence[str] | None = None,
        structure: Structure | None = None,
        wyckoff_sites: Sequence[dict] | None = None,
        smodes: object | None = None,
        kpoints: Sequence | None = None,
        symmetry_info: dict | None = None,
    ) -> Method2Result:
        """
        官网 Method 2：在已枚举子群中按序号选择目标子群，计算其畸变模式。

        特殊 k：iso DISPLAY BUSH 根模 + rootless Wyckoff 轨道的子群固定
        子空间。参数 k：smodes + 子群对称性（(3+d)/公度锁定）。

        Args:
            parent_sg: 母相空间群号
            subgroups: 子群候选列表（来自 Method 1 或 list_subgroups）
            query: Method 2 查询参数（subgroup_idx 必填）
            wyckoff_letters: 母相结构各原子的 Wyckoff 位置字母
            wyckoff_orbit_ids: 作用域内物理 Wyckoff 轨道的稳定身份

        Returns:
            Method2Result
        """
        target = next((s for s in subgroups if s.index == query.subgroup_idx), None)
        if target is None:
            raise ValueError(
                f"Subgroup index {query.subgroup_idx} not found; 请先执行 Method 1 或 list_subgroups 获得候选列表"
            )

        validate_single_irrep_count(query.number_of_superposed_irreps)

        requested_nmod = validate_modulation_count(query.number_of_independent_modulations)
        # nmod is meaningful only for parameterized k vectors.  Preserve the
        # caller-visible request in metadata for API compatibility, while the
        # special-k computation itself always uses the ordinary 3D lock-in.
        effective_nmod = requested_nmod if target.k_parameters else 0

        # wyckoff_letters 为 None 表示调用方未提供（误用）；空列表表示
        # 作用域内无 Wyckoff 位置（如全部类型选 none）-> 直接返回空模式
        if wyckoff_letters is None:
            raise ValueError("Method 2 计算模式需要母相结构的 Wyckoff 位置信息，请先加载结构（load_structure）")
        modes = []
        extra_meta: dict = {}
        if wyckoff_letters:
            method3_resolution = str(getattr(target, "method3_route_resolution", "") or "")
            if method3_resolution == "exact_fixed_space" and structure is not None and smodes is not None:
                from backend.method2.superspace import compute_parametric_modes  # noqa: PLC0415

                info = dict(symmetry_info or {})
                info.setdefault("space_group_number", parent_sg)
                if wyckoff_sites is not None:
                    info["wyckoff_sites"] = list(wyckoff_sites)
                complete_result = compute_parametric_modes(
                    structure,
                    info,
                    target,
                    wyckoff_letters,
                    smodes,  # type: ignore[arg-type]
                    iso=self._iso,
                    kpoints=kpoints,
                    # A coupled embedding has no single distinguished q.
                    # Its complete 3D fixed space must include every k class
                    # folded by the selected child translation lattice.
                    nmod=0,
                    wyckoff_orbit_ids=wyckoff_orbit_ids,
                )
                modes = complete_result.modes
                extra_meta["supercell_displacements"] = complete_result.supercell_displacements
                extra_meta["mode_labels"] = complete_result.labels
                extra_meta["coupled_mode_note"] = complete_result.note
            elif target.k_parameters and structure is not None and smodes is not None:
                from backend.method2.superspace import compute_parametric_modes  # noqa: PLC0415

                info = dict(symmetry_info or {})
                info.setdefault("space_group_number", parent_sg)
                if wyckoff_sites is not None:
                    info["wyckoff_sites"] = list(wyckoff_sites)
                complete_result = compute_parametric_modes(
                    structure,
                    info,
                    target,
                    wyckoff_letters,
                    smodes,  # type: ignore[arg-type]
                    iso=self._iso,
                    kpoints=kpoints,
                    nmod=effective_nmod,
                    wyckoff_orbit_ids=wyckoff_orbit_ids,
                )
                modes = complete_result.modes
                extra_meta["supercell_displacements"] = complete_result.supercell_displacements
                extra_meta["mode_labels"] = complete_result.labels
                extra_meta["parametric_note"] = complete_result.note
            elif structure is not None and smodes is not None:
                from backend.method2.superspace import (  # noqa: PLC0415
                    compute_special_modes_with_rootless_supplement,
                )

                info = dict(symmetry_info or {})
                info.setdefault("space_group_number", parent_sg)
                if wyckoff_sites is not None:
                    info["wyckoff_sites"] = list(wyckoff_sites)
                complete_result = compute_special_modes_with_rootless_supplement(
                    structure,
                    info,
                    target,
                    wyckoff_letters,
                    self._iso,
                    smodes,  # type: ignore[arg-type]
                    kpoints=kpoints,
                    wyckoff_orbit_ids=wyckoff_orbit_ids,
                )
                modes = complete_result.modes
                extra_meta["supercell_displacements"] = complete_result.supercell_displacements
                extra_meta["mode_labels"] = complete_result.labels
                extra_meta["special_mode_note"] = complete_result.note
            else:
                # Backend-only callers have no structure on which to build the
                # child displacement representation.  They still receive the
                # exact BUSH root modes; the high-level API always supplies the
                # structure and therefore also completes rootless orbits.
                modes = self._iso.calc_distortion_modes(parent_sg, target, wyckoff_letters=wyckoff_letters)

        metadata = {
            # k 点 / IR / OPD 来自所选子群自身（查询时已确定，无需重复传参）
            "k_point_label": target.k_point_label,
            "irrep_label": target.irrep_label,
            "opd_symbol": target.opd_symbol,
            "k_parameters": list(target.k_parameters),
            "number_of_independent_modulations": requested_nmod,
            **extra_meta,
        }

        return Method2Result(subgroup=target, modes=modes, metadata=metadata)
